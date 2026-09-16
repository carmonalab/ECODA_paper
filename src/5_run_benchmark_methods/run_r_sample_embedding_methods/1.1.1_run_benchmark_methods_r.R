# ==============================================================================
# 1.1.1_run_benchmark_methods_r.R — Run one R benchmark method (gloscope,
# mofa, pseudobulk, scitd or composition) for one dataset (Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --method {gloscope,mofa,
#   pseudobulk,scitd,composition} --input_dir --results_dir --pseudobulk_dir
#   --gloscope_cache_dir --log_file [--force]
# Canonical ordinary and batch pseudobulk paths use the raw H5AD CSR Sample
# aggregate plus the direct matrix DESeq2 API.  They never materialize a
# sample-level Seurat object.  Count-backed Seurat remains reserved for
# genuine cell-level methods such as scITD; CT pseudobulk calls
# process_pseudobulk_ct_h5ad_fig() directly.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
# Artifact persistence is provided by benchmark_hpc_utils.R.
# Its checked-load and publication functions are the single R artifact seam.


raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_flags(raw_args)

for (req in c("config_path", "ds_name", "view", "method", "input_dir",
              "results_dir", "log_file")) {
  if (is.null(args[[req]]) || identical(args[[req]], TRUE)) {
    stop("Missing required --", req, " argument")
  }
}
force <- isTRUE(args[["force"]]) || identical(args[["force"]], "TRUE")

method <- args$method
script_identity_arg <- commandArgs(trailingOnly = FALSE)
script_identity <- sub(
  "^--file=",
  "",
  script_identity_arg[grepl("^--file=", script_identity_arg)][1L]
)
message(
  "ECODA_R_DISPATCH_PARSED method=", method,
  " source_root=", Sys.getenv("ECODA_SOURCE_ROOT", unset = ""),
  " script=", script_identity
)
if (!method %in% c("gloscope", "mofa", "pseudobulk", "scitd",
                   "composition")) {
  stop("Unknown method '", method,
       "' (expected gloscope, mofa, pseudobulk, scitd or composition)")
}
combo_supplied <- !is.null(args[["combo"]])
combo_token <- args[["combo"]]
if (combo_supplied && method != "gloscope") {
  stop("--combo is only supported for method gloscope")
}

# A combo shard is intentionally ordinary-only. Batch-effect GloScope keeps
# its existing single-combo behavior and pass-qualified artifact names.

# Method-specific attaches: MOFA2/scITD are needed only by their methods
# (bare create_mofa / initialize_params + make_new_container); gloscope needs
# only the installed namespace (GloScope::gloscope is called qualified).
# composition calls the EPIC::EPIC + GloScope::gloscopeProp drivers BARE
# (in benchmark_methods_r.R; EPIC/GloScope are not attached by any loader),
# so both must be attached for it.
if (method == "mofa") library(MOFA2)
if (method == "scitd") library(scITD)
if (method == "composition") {
  library(EPIC)
  library(GloScope)
}

config <- read_datasets_json(args$config_path, view = args$view)
ds <- args$ds_name
analysis_pass <- args[["analysis_pass"]]
analysis_variant <- Sys.getenv("ANALYSIS_VARIANT", unset = "")
if (!analysis_variant %in% c("", "final", "corrected_final")) {
  stop("Unknown analysis variant: ", analysis_variant)
}
if (identical(analysis_variant, "final") &&
    !identical(analysis_pass, "uncorrected")) {
  stop("final analysis variant requires the uncorrected batch-effect pass")
}
if (identical(analysis_variant, "corrected_final") &&
    !identical(analysis_pass, "corrected")) {
  stop("corrected_final analysis variant requires the corrected batch-effect pass")
}
if (!is.null(analysis_pass) && !analysis_pass %in% c("uncorrected", "corrected")) {
  stop("Unknown analysis pass: ", analysis_pass)
}
if (combo_supplied && !is.null(analysis_pass)) {
  stop("--combo is only supported for ordinary GloScope runs")
}
cache_stem <- if (is.null(analysis_pass)) {
  ds
} else {
  stem <- paste0(ds, "_batch_effect_", analysis_pass)
  if (analysis_variant %in% c("final", "corrected_final")) {
    stem <- paste0(stem, "_final")
  }
  stem
}
entry <- config[[ds]]
if (is.null(entry)) {
  stop("Dataset '", ds, "' not found in ", args$config_path)
}

h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}
# Source identity is run-owned when available; the CT store manifest records
# this identity alongside its H5AD path and scheduler owner.
source_identity <- Sys.getenv("ECODA_SOURCE_IDENTITY", unset = "")
if (!nzchar(source_identity)) {
  run_root_for_identity <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
  if (nzchar(run_root_for_identity)) {
    candidate_identity <- file.path(
      run_root_for_identity, "manifests", "source_identity.json"
    )
    if (file.exists(candidate_identity)) source_identity <- candidate_identity
  }
}
if (!nzchar(source_identity)) source_identity <- NULL
sample_col <- "Sample"
correct_batch_mode <- identical(analysis_pass, "corrected")
python_batch_metadata <- NULL
corrected_final_mode <- identical(analysis_variant, "corrected_final") &&
  correct_batch_mode
if (is.null(analysis_pass)) {
  Sys.unsetenv("ANALYSIS_PASS")
} else {
  Sys.setenv(ANALYSIS_PASS = as.character(analysis_pass))
}
if (corrected_final_mode) {
  # The final corrected lane is explicitly summary-free; the pipeline's
  # compatibility call reads this flag from the worker environment.
  Sys.setenv(ECODA_ALLOW_MISSING_CORRECTED_SUMMARY = "1")
} else {
  # Never let a stale worker environment relax ordinary or uncorrected lanes.
  Sys.unsetenv("ECODA_ALLOW_MISSING_CORRECTED_SUMMARY")
}
batch_keys <- NULL
batch_context <- NULL
effective_batch_keys <- NULL
batch_col <- NULL
h5ad_expected_batch_contract <- NULL
pseudobulk_batch_contract <- NULL
method_batch_contract <- NULL
if (correct_batch_mode) {
  if (is.null(entry$batch_col)) {
    stop("corrected batch-effect view requires a confirmed columns.batch")
  }
  batch_keys <- ecoda_hpc_normalize_keys(
    entry$batch_col,
    sample_col = sample_col,
    biological_label = entry$label_col
  )
  # Preserve the original configured columns.  Only a one-key design uses
  # batch_col as a direct compatibility argument; multi-key correction is
  # supplied through batch_keys and batch_context.
  if (length(batch_keys) == 1L) batch_col <- batch_keys[[1L]]
  h5ad_expected_batch_contract <- ecoda_hpc_batch_contract_identity(
    batch_keys,
    sample_col = sample_col,
    method_id = "preprocess",
    model_id = "hvg_composite_v1"
  )
  pseudobulk_batch_contract <- ecoda_hpc_batch_contract_identity(
    batch_keys,
    sample_col = sample_col,
    method_id = "Pseudobulk",
    model_id = "pseudobulk_limma_fixed_effects_v1"
  )
  method_batch_contract <- switch(
    method,
    composition = ecoda_hpc_batch_contract_identity(
      batch_keys,
      sample_col = sample_col,
      method_id = "ECODA_authors_HR",
      model_id = "limma_fixed_effects_v1"
    ),
    gloscope = ecoda_hpc_batch_contract_identity(
      batch_keys,
      sample_col = sample_col,
      method_id = "GloScope",
      model_id = "embedding_consumer_harmony_v1"
    ),
    pseudobulk = pseudobulk_batch_contract,
    NULL
  )
}
dir.create(args$results_dir, showWarnings = FALSE, recursive = TRUE)

method_rds_stem <- if (is.null(analysis_pass)) ds else cache_stem
method_rds <- file.path(
  args$results_dir,
  paste0(method_rds_stem, "_", method, ".rds")
)

# GloScope and composition use the embedding/obs-only loader.  MOFA and
# pseudobulk use a stricter metadata/HVG-only reader so complete cache paths
# never open or read count values.  scITD remains the sole canonical path
# below that materializes a count-backed Seurat object.
counts_free_method <- method %in% c("gloscope", "composition")
pseudobulk_metadata_method <- method %in% c("mofa", "pseudobulk")
dispatch_branch <- if (counts_free_method) {
  "counts_free"
} else if (pseudobulk_metadata_method) {
  "metadata"
} else if (method == "scitd") {
  "count_backed_seurat"
} else {
  "unknown"
}
message(
  "ECODA_R_DISPATCH_PRE method=", method,
  " branch=", dispatch_branch,
  " analysis_pass=", ifelse(is.null(analysis_pass), "", analysis_pass),
  " force=", force
)
embedding_key <- if (args$view == "batch_effect_corrected") {
  "X_pca_harmony_batch_effect_corrected_hvg2000"
} else if (args$view == "batch_effect_uncorrected") {
  "X_pca_batch_effect_uncorrected_hvg2000"
} else {
  "X_pca_benchmark_analysis_hvg2000"
}
# Validate every cell before any metadata reducer selects a first row per
# Sample. This is retained for ordinary corrected Stage 5 and intentionally
# skipped only by corrected_final, which consumes the exported Feather.
if (correct_batch_mode && !corrected_final_mode) {
  # The Python check validates the selected H5AD's source/preprocessing
  # metadata. Its identity is intentionally the persisted H5AD contract; the
  # Stage 5 consumer contract uses the separate limma model identity.
  validation_method_id <- "preprocess"
  validation_model_id <- "hvg_composite_v1"
  python_batch_metadata <- validate_h5ad_corrected_batch_metadata(
    h5ad_path = h5ad_path,
    batch_keys = as.list(unname(batch_keys)),
    sample_col = sample_col,
    biological_label = entry$label_col,
    method_id = validation_method_id,
    model_id = validation_model_id
  )
}
required_hvg <- if (identical(args$view, "benchmark_analysis")) 3000L else 2000L
hvg_rank_genes <- NULL
embedding_matrices <- NULL
embedding_sample_ids <- NULL
if (pseudobulk_metadata_method) {
  ct_columns <- if (method == "pseudobulk" && is.null(analysis_pass)) {
    c(entry$cell_type_low_res, entry$cell_type_high_res)
  } else {
    character()
  }
  metadata_info <- load_h5ad_pseudobulk_metadata(
    h5ad_path,
    sample_col = sample_col,
    metadata_columns = unique(c(
      entry$label_col,
      if (correct_batch_mode) batch_keys else batch_col,
      ct_columns
    )),
    n_hvg = required_hvg,
    required_nonmissing_columns = unique(c(
      sample_col,
      entry$label_col,
      if (correct_batch_mode) batch_keys else batch_col
    )),
    expected_batch_contract = h5ad_expected_batch_contract,
    view = args$view,
    method = method,
    allow_missing_summary = corrected_final_mode
  )
  obs <- metadata_info$obs
  hvg_rank_genes <- metadata_info$hvg_rank_genes
} else if (method == "gloscope") {
  # GloScope is explicitly counts-free: validate the persisted counts-layer
  # contract through load_h5ad_counts_free(), but materialize only the
  # configured metadata and stored semantic PCA embeddings.  Keep seurat NULL
  # so this branch cannot fall through to the count-backed Seurat path.
  embedding_names <- if (is.null(analysis_pass)) {
    c(
      hvg1000 = "X_pca_benchmark_analysis_hvg1000",
      hvg2000 = "X_pca_benchmark_analysis_hvg2000",
      hvg3000 = "X_pca_benchmark_analysis_hvg3000"
    )
  } else {
    c(hvg2000 = embedding_key)
  }
  adata <- load_h5ad_counts_free(
    h5ad_path,
    unique(c(
      sample_col,
      entry$label_col,
      if (correct_batch_mode) batch_keys else NULL
    )),
    unname(embedding_names),
    obs_prefixes = character(),
    view = args$view,
    method = method,
    expected_batch_contract = h5ad_expected_batch_contract,
    allow_missing_summary = corrected_final_mode
  )
  obs <- py_to_r(adata$obs)
  hvg_rank_genes <- get_hvg_rank_genes(adata)
  embedding_matrices <- lapply(unname(embedding_names), function(key) {
    py_to_r(adata$obsm[[key]])
  })
  names(embedding_matrices) <- names(embedding_names)
  embedding_sample_ids <- as.character(obs[[sample_col]])
} else if (method == "composition") {
  # Composition keeps its existing obs-only path. It consumes the stored
  # hvg2000 PCA embedding, the configured high-resolution cell-type column for
  # batch-effect views, and cached pseudobulk variants.
  composition_obs_columns <- if (
    is.null(analysis_pass) &&
    length(entry$not_suitable_for_auto_annotation) == 0
  ) {
    c("layer2", "scATOMIC_pred")
  } else {
    character()
  }
  composition_cell_type_columns <- if (!is.null(analysis_pass)) {
    entry$cell_type_high_res
  } else {
    c(entry$cell_type_low_res, entry$cell_type_high_res)
  }
  obs_columns <- c(
    sample_col,
    entry$label_col,
    composition_cell_type_columns,
    if (correct_batch_mode) batch_keys else batch_col,
    composition_obs_columns
  )
  adata <- load_h5ad_counts_free(
    h5ad_path,
    unique(obs_columns[!is.na(obs_columns) & nzchar(obs_columns)]),
    embedding_key,
    obs_prefixes = "leiden_res_",
    view = args$view,
    method = method,
    expected_batch_contract = h5ad_expected_batch_contract,
    allow_missing_summary = corrected_final_mode
  )
  obs <- py_to_r(adata$obs)
  hvg_rank_genes <- get_hvg_rank_genes(adata)
} else {
  ad <- import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(h5ad_path, backed = "r")
  obs <- py_to_r(adata$obs)
  validate_benchmark_h5ad_contract(
    adata,
    obs = obs,
    view = args$view,
    method = method,
    expected_batch_contract = h5ad_expected_batch_contract,
    allow_missing_summary = corrected_final_mode
  )
  hvg_rank_genes <- get_hvg_rank_genes(adata)
}

if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}
blind_mode <- is.null(analysis_pass) || analysis_pass == "uncorrected"
if (correct_batch_mode) {
  if (corrected_final_mode) {
    missing_batch_keys <- setdiff(batch_keys, colnames(obs))
    if (length(missing_batch_keys) > 0L) {
      stop(
        "Confirmed batch column(s) missing from obs of ", h5ad_path, ": ",
        paste(missing_batch_keys, collapse = ", ")
      )
    }
    expected_sample_ids <- unique(as.character(obs[[sample_col]]))
    allow_unknown_keys <- if (identical(ds, "Breast_cancer")) {
      intersect(batch_keys, "suspension_dissociation_time")
    } else {
      character()
    }
    batch_context <- ecoda_hpc_load_sample_metadata_contract(
      path = ecoda_hpc_sample_metadata_path(ds),
      expected_sample_ids = expected_sample_ids,
      batch_keys = batch_keys,
      sample_col = sample_col,
      biological_label = entry$label_col,
      required_columns = entry$label_col,
      allow_unknown_keys = allow_unknown_keys
    )
  } else {
    batch_context <- ecoda_hpc_batch_context(
      metadata = obs,
      batch_keys = as.list(unname(batch_keys)),
      sample_col = sample_col,
      biological_label = entry$label_col,
      python_metadata = python_batch_metadata
    )
  }
  effective_batch_keys <- batch_context[["effective_batch_keys"]]
  if (is.null(effective_batch_keys)) {
    effective_batch_keys <- batch_context$validation[["effective_batch_keys"]]
  }
  if (is.null(effective_batch_keys)) effective_batch_keys <- character()
  effective_batch_keys <- unname(as.character(effective_batch_keys))
  if (length(effective_batch_keys) &&
      any(!effective_batch_keys %in% batch_keys)) {
    stop("Corrected worker effective keys are not configured technical columns")
  }
  # Keep non-estimable configured keys in metadata, but pass only an original
  # varying key through the scalar compatibility slot. Multi-key correction
  # uses the separate columns carried by batch_context.
  batch_col <- if (length(effective_batch_keys) == 1L) {
    effective_batch_keys[[1L]]
  } else {
    NULL
  }
}
if (correct_batch_mode) {
  pseudobulk_batch_contract <- ecoda_hpc_augment_batch_contract(
    identity = pseudobulk_batch_contract,
    validation = batch_context$validation,
    method_id = "Pseudobulk",
    batch_keys = batch_keys
  )
  method_batch_contract <- switch(
    method,
    composition = ecoda_hpc_augment_batch_contract(
      identity = method_batch_contract,
      validation = batch_context$validation,
      method_id = "ECODA_authors_HR",
      batch_keys = batch_keys
    ),
    gloscope = ecoda_hpc_augment_batch_contract(
      identity = method_batch_contract,
      validation = batch_context$validation,
      method_id = "GloScope",
      batch_keys = batch_keys
    ),
    pseudobulk = pseudobulk_batch_contract,
    NULL
  )
}

if (!combo_supplied && artifact_checksum_ok(method_rds) && !force) {
  message("Method results already exist and passed checksum/record validation: ", method_rds)
  cached <- read_rds_checked(method_rds)
  if (!is.list(cached)) {
    stop("Method results artifact is not a list: ", method_rds)
  }
  if (correct_batch_mode) {
    ecoda_hpc_validate_batch_contract(
      cached[["batch_contract"]],
      method_batch_contract,
      label = paste0("Method results ", ds, "/", method)
    )
    if (method == "composition") {
      required_composition_methods <- c(
        "ECODA_authors_HR",
        "ECODA_authors_HR_NULL",
        "ECODA_seuratres_2"
      )
      missing_composition_methods <- setdiff(
        required_composition_methods,
        names(cached)
      )
      if (length(missing_composition_methods) > 0L) {
        stop(
          "Method results ", ds, "/", method,
          " is missing corrected composition results: ",
          paste(missing_composition_methods, collapse = ", ")
        )
      }
      for (composition_method in required_composition_methods) {
        composition_contract <- ecoda_hpc_augment_batch_contract(
          identity = ecoda_hpc_batch_contract_identity(
            batch_keys,
            sample_col = sample_col,
            method_id = composition_method,
            model_id = "limma_fixed_effects_v1"
          ),
          validation = batch_context$validation,
          method_id = composition_method,
          batch_keys = batch_keys
        )
        ecoda_hpc_validate_batch_contract(
          cached[[composition_method]][["batch_contract"]],
          composition_contract,
          label = paste0(
            "Method result ", ds, "/", method, "/", composition_method
          )
        )
      }
    }
  }
  shared_replayed <- character()
  shared_rows <- list()
  for (nm in setdiff(
    names(cached),
    c("batch_contract", "batch_contract_identity")
  )) {
    value <- cached[[nm]]
    validate_hpc_timing_bundle(
      value, label = paste0("Method result ", ds, "/", method, "/", nm)
    )
    if ("timing_schema" %in% names(value)) {
      timing_id <- as.character(value[["timing_id"]])
      shared_method <- if ("shared_timing_method" %in% names(value)) {
        as.character(value[["shared_timing_method"]])
      } else if (grepl("^Pseudobulk_CT_", nm)) {
        ct_shared_timing_method_from_timing_id(
          timing_id, fallback_method = nm
        )
      } else {
        "prepare_pseudobulk_shared"
      }
      shared_key <- paste(shared_method, timing_id, sep = "\r")
      if (!shared_key %in% shared_replayed) {
        shared_rows[[shared_key]] <- list(
          method = shared_method,
          time = as.numeric(value[["shared_time_secs"]]),
          mem = value[["shared_mem_GB"]]
        )
        shared_replayed <- c(shared_replayed, shared_key)
      }
    }
    if (!is.null(value$exec_time)) {
      log_exec_row(ds, nm, value$exec_time, args$log_file,
                   mem_gb = value$mem_GB)
    }
  }
  for (shared in shared_rows) {
    log_exec_row(
      ds, shared[["method"]], shared[["time"]], args$log_file,
      mem_gb = shared[["mem"]]
    )
  }
  quit(save = "no", status = 0)
}
pb_variants <- NULL
seurat <- NULL
metadata <- NULL
labels <- NULL
if (method == "gloscope") {
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
} else if (method %in% c("mofa", "pseudobulk")) {
  # Both methods consume direct matrix pseudobulks.  Cache validation occurs
  # inside load_pb_variants before the H5AD raw counts pass; a complete cache
  # set therefore performs no aggregation and creates no Seurat object.
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method ", method)
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  requested_pb_variants <- if (method == "mofa" || is.null(analysis_pass)) {
    PB_VARIANT_NAMES
  } else {
    "hvg2000"
  }
  pb_variants <- load_pb_variants(
    seurat = NULL,
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir,
    ds = ds,
    force = force,
    log_file = args$log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = requested_pb_variants,
    h5ad_path = h5ad_path,
    view = args$view,
    analysis_pass = analysis_pass,
    run_id = Sys.getenv("ECODA_RUN_ID", unset = ""),
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_contract = if (correct_batch_mode) pseudobulk_batch_contract else NULL,
    expected_h5ad_batch_contract = h5ad_expected_batch_contract
  )
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  if (correct_batch_mode) {
    metadata <- ecoda_hpc_apply_batch_context(
      metadata,
      batch_context,
      sample_col = sample_col
    )
  }
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
} else if (method == "composition") {
  # Obs-only path: no Seurat materialization. Consumes the backed h5ad obs
  # (cell-level metadata), the hvg2000 obsm PCA embedding (Avg_PCA_embedding)
  # and the precomputed hvg2000 pseudobulk variant (ECODA_deconv; the submit
  # script auto-prepends prepare_pseudobulk for composition).
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method composition")
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  # Keep a configured Leiden annotation source column alongside its legacy
  # RNA_snn_res.* alias. Parkinson declares leiden_res_5_* as its author
  # high-resolution column; renaming it in place would invalidate that config.
  obs <- rename_leiden_cols(
    obs,
    view = args$view,
    preserve_source = TRUE
  )
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  if (correct_batch_mode) {
    metadata <- ecoda_hpc_apply_batch_context(
      metadata,
      batch_context,
      sample_col = sample_col
    )
  }
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
  obsm_keys <- py_to_r(import_builtins(convert = FALSE)$list(
    adata$obsm$keys()
  ))
  emb_key <- embedding_key
  if (!emb_key %in% obsm_keys) {
    stop("Embedding '", emb_key, "' not found in adata.obsm of ", h5ad_path,
         ". Re-run preprocessing (1.1.1_preprocess.py) for this dataset.")
  }
  pca_emb <- py_to_r(adata$obsm[[emb_key]])
  if (is.null(rownames(pca_emb))) rownames(pca_emb) <- rownames(obs)
  colnames(pca_emb) <- paste0("PC_", seq_len(ncol(pca_emb)))
  pb_variants <- load_composition_pb_variants(
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir,
    ds = ds,
    log_file = args$log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = if (!is.null(analysis_pass)) "hvg2000" else PB_VARIANT_NAMES,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_contract = if (correct_batch_mode) pseudobulk_batch_contract else NULL,
    expected_h5ad_batch_contract = h5ad_expected_batch_contract
  )
} else {
  # scITD is the genuine cell-level count consumer.  Keep its existing
  # count-backed Seurat boundary; ordinary pseudobulk and CT never enter this
  # branch.
  seurat <- load_benchmark_seurat(
    adata, obs, sample_col = sample_col,
    fetch_embedding = NULL,
    counts_layer = "counts"
  )
  # Sample names are already standardized in the preprocessed obs
  # (1.1.1_preprocess.py); no standardize_sample_names() re-application.
  if (length(hvg_rank_genes) > 0) {
    VariableFeatures(seurat) <- hvg_rank_genes[
      seq_len(min(2000, length(hvg_rank_genes)))
    ]
  }
  seurat@misc$label_col <- entry$label_col
  seurat@misc$cell_type_low_res <- entry$cell_type_low_res
  seurat@misc$cell_type_high_res <- entry$cell_type_high_res
  metadata <- get_metadata(seurat)
  labels <- get_labels(seurat, entry$label_col)
}
message(
  "ECODA_R_DISPATCH_FINAL method=", method,
  " branch=", dispatch_branch,
  " seurat_is_null=", is.null(seurat),
  " embedding_matrices_is_null=", is.null(embedding_matrices)
)

results <- switch(
  method,
  gloscope = run_gloscope_hpc(
    seurat, metadata, label_col = entry$label_col,
    sample_col = sample_col,
    gloscope_cache_dir = args$gloscope_cache_dir,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file,
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    combo_token = combo_token,
    embedding_name = if (!is.null(analysis_pass)) {
      sub("^X_", "", embedding_key)
    } else {
      NULL
    },
    embedding_matrices = embedding_matrices,
    embedding_sample_ids = embedding_sample_ids,
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL
  ),
  mofa = run_mofa_hpc(
    metadata, labels, pb_variants,
    results_dir = args$results_dir,
    ds = ds,
    force = force,
    log_file = args$log_file
  ),
  pseudobulk = run_pseudobulk_hpc(
    labels = labels,
    pb_variants = pb_variants,
    sample_col = sample_col,
    results_dir = args$results_dir,
    ds = ds,
    force = force,
    log_file = args$log_file,
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    h5ad_path = h5ad_path,
    ct_col_low_res = entry$cell_type_low_res,
    ct_col_high_res = entry$cell_type_high_res,
    view = args$view,
    analysis_pass = analysis_pass,
    run_id = Sys.getenv("ECODA_RUN_ID", unset = ""),
    source_identity = source_identity,
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL
  ),
  scitd = run_scitd_hpc(
    seurat, label_col = entry$label_col,
    hvg_sets = make_hvg_sets(hvg_rank_genes),
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  composition = run_composition_methods_hpc(
    labels, metadata, pca_emb, pb_variants[["hvg2000"]], obs,
    label_col = entry$label_col,
    ct_col_low_res = entry$cell_type_low_res,
    ct_col_high_res = entry$cell_type_high_res,
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file,
    seurat_res = if (!is.null(analysis_pass)) 2 else c(0.1, 0.4, 2, 5, 20),
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    batch_col = batch_col,
    corrected = correct_batch_mode,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    metadata_validation = if (correct_batch_mode) {
      batch_context$validation
    } else {
      NULL
    },
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    not_suitable_for_auto_annotation = if (
      is.null(entry$not_suitable_for_auto_annotation)
    ) {
      character(0)
    } else {
      entry$not_suitable_for_auto_annotation
    }
  )
)

if (combo_supplied) {
  message(
    "Saved GloScope combo ", combo_token,
    " per-combo bundle/cache; method-level RDS deferred to consolidation"
  )
} else {
  save_rds_atomic(
    results, method_rds, producer = paste0("stage5_", method)
  )
  message("Saved method results: ", method_rds, " (", length(results), " combos)")
}
message("--- ", method, " for ", ds, " complete ---")
