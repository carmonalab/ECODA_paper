# ==============================================================================
# 1.1.1_run_benchmark_methods_r.R — Run one R benchmark method (gloscope,
# mofa, pseudobulk, scitd or composition) for one dataset (Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --method {gloscope,mofa,
#   pseudobulk,scitd,composition} --input_dir --results_dir --pseudobulk_dir
#   --gloscope_cache_dir --log_file [--force]
# Counts-dependent methods load the raw counts layer through reticulate.
# GloScope and composition instead use the h5py/minimal-AnnData loader with
# only required obs columns and precomputed embeddings, so anndata's backed
# open never materializes layers["counts"] for those methods. The scripts set
# seurat@misc$cell_type_low_res / label_col from datasets.json, dispatch on
# method to the T2 driver (benchmark_pipeline.R), write the per-combo bundle
# files (<ds>_<combo>.rds; combo names are method-prefixed) + the method-level
# RDS (<ds>_<method>.rds, a named list of result bundles) + per-combo exec-log
# rows. An optional --combo hvg{n}_pcadims{d} argument selects one ordinary
# GloScope combo; that shard writes only its per-combo bundle/cache and never
# reads or writes the shared method-level RDS. With no --combo, skip-if-exists
# behavior remains unchanged.
#
# Memory: MOFA and batch-mode pseudobulk consume sample metadata and
# precomputed pseudobulks. Missing variants use bounded H5AD CSR aggregation
# into a sample-level object; ordinary pseudobulk CT/scITD retain cell counts.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
# Cache records are run-owned and optional for legacy artifacts.  A record
# permits a size/sidecar comparison before the read boundary; the full MD5
# check below is still mandatory immediately before every readRDS().
ecoda_local_validate_run_id <- function(run_id, label = "artifact producer run ID") {
  if (!is.character(run_id) || length(run_id) != 1L ||
      is.na(run_id) || !nzchar(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id, perl = TRUE)) {
    stop(label, " is invalid")
  }
  run_id
}

ecoda_local_current_run_id <- function() {
  run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  if (!nzchar(run_id)) {
    run_root <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
    if (nzchar(run_root)) {
      run_id <- basename(normalizePath(path.expand(run_root), mustWork = FALSE))
    }
  }
  if (!nzchar(run_id)) return(NULL)
  ecoda_local_validate_run_id(run_id, "ECODA_RUN_ID")
}

ecoda_local_record_run_id <- function(run_id = NULL) {
  if (!is.null(run_id)) {
    return(ecoda_local_validate_run_id(run_id))
  }
  explicit <- Sys.getenv("ECODA_ARTIFACT_PRODUCER_RUN_ID", unset = "")
  if (nzchar(explicit)) {
    return(ecoda_local_validate_run_id(
      explicit, "ECODA_ARTIFACT_PRODUCER_RUN_ID"
    ))
  }
  ecoda_local_current_run_id()
}

ecoda_local_runs_root <- function() {
  runs_root <- Sys.getenv("ECODA_RUNS_ROOT", unset = "")
  if (!nzchar(runs_root)) {
    run_root <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
    if (nzchar(run_root)) {
      runs_root <- dirname(path.expand(run_root))
    }
  }
  if (!nzchar(runs_root)) {
    scratch_root <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
    if (nzchar(scratch_root)) {
      runs_root <- file.path(path.expand(scratch_root), "_ecoda_runs")
    }
  }
  if (!nzchar(runs_root) || !grepl("^/", path.expand(runs_root))) {
    return(NULL)
  }
  path.expand(runs_root)
}

ecoda_local_canonical_path <- function(path) {
  if (!is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(path)) {
    stop("artifact path must be one non-empty string")
  }
  normalizePath(path.expand(path), mustWork = FALSE)
}

ecoda_local_record_path <- function(path, run_id) {
  run_id <- ecoda_local_validate_run_id(run_id)
  runs_root <- ecoda_local_runs_root()
  if (is.null(runs_root)) return(NULL)
  canonical <- ecoda_local_canonical_path(path)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for artifact-record paths")
  }
  key <- digest::digest(canonical, algo = "sha256", serialize = FALSE)
  if (!grepl("^[0-9a-f]{64}$", key, perl = TRUE)) {
    stop("could not derive artifact-record key for ", canonical)
  }
  file.path(
    runs_root, run_id, "manifests", "artifacts",
    paste0(substr(key, 1L, 32L), ".record")
  )
}

ecoda_local_checksum_fields <- function(path) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(path) || !file.exists(sidecar)) return(NULL)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) return(NULL)
  lines <- tryCatch(readLines(sidecar, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) != 3L ||
      any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    return(NULL)
  }
  fields <- setNames(sub("^[^=]*=", "", lines), c("MD5", "SIZE", "PATH"))
  if (!grepl("^[0-9a-f]{32}$", fields[["MD5"]], perl = TRUE) ||
      !grepl("^[0-9]+$", fields[["SIZE"]], perl = TRUE) ||
      !nzchar(fields[["PATH"]])) {
    return(NULL)
  }
  fields
}

ecoda_local_strict_checksum <- function(path, description = "RDS artifact") {
  fields <- ecoda_local_checksum_fields(path)
  if (is.null(fields) ||
      !identical(fields[["PATH"]], path) ||
      !identical(fields[["SIZE"]], as.character(file.info(path)$size))) {
    return(FALSE)
  }
  actual <- tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  isTRUE(length(actual) == 1L && !is.na(actual) &&
         identical(fields[["MD5"]], actual))
}

ecoda_local_read_record <- function(path, run_id = NULL, producer = NULL) {
  run_id <- ecoda_local_record_run_id(run_id)
  if (is.null(run_id)) {
    return(list(state = "unavailable", record_path = NULL))
  }
  record_path <- ecoda_local_record_path(path, run_id)
  if (is.null(record_path) || !file.exists(record_path)) {
    return(list(state = "absent", record_path = record_path, run_id = run_id))
  }
  lines <- tryCatch(readLines(record_path, warn = FALSE), error = function(e) NULL)
  expected_keys <- c("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
  if (is.null(lines) || length(lines) != length(expected_keys) ||
      any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), expected_keys)) {
    return(list(state = "invalid", record_path = record_path, run_id = run_id))
  }
  record <- setNames(sub("^[^=]*=", "", lines), expected_keys)
  canonical <- ecoda_local_canonical_path(path)
  valid <- identical(record[["PATH"]], canonical) &&
    grepl("^[0-9]+$", record[["SIZE"]], perl = TRUE) &&
    grepl("^[0-9a-f]{32}$", record[["MD5"]], perl = TRUE) &&
    identical(record[["RUN_ID"]], run_id) &&
    nzchar(record[["PRODUCER"]]) &&
    identical(record[["STATE"]], "PUBLISHED")
  if (!is.null(producer) &&
      (!is.character(producer) || length(producer) != 1L ||
       is.na(producer) || !identical(record[["PRODUCER"]], producer))) {
    valid <- FALSE
  }
  if (!valid) {
    return(list(state = "invalid", record_path = record_path, run_id = run_id))
  }
  list(
    state = "valid", record_path = record_path, run_id = run_id,
    canonical_path = canonical, record = record
  )
}

ecoda_local_record_only_valid <- function(path, record_info) {
  if (!identical(record_info[["state"]], "valid")) return(FALSE)
  fields <- ecoda_local_checksum_fields(path)
  if (is.null(fields)) return(FALSE)
  info <- file.info(path)
  record <- record_info[["record"]]
  isTRUE(
    !is.na(info$size) && info$size > 0 &&
      identical(record[["SIZE"]], as.character(info$size)) &&
      identical(record[["MD5"]], fields[["MD5"]]) &&
      identical(record[["SIZE"]], fields[["SIZE"]]) &&
      identical(fields[["PATH"]], path)
  )
}

ecoda_local_cache_valid <- function(path, producer = NULL) {
  record_info <- ecoda_local_read_record(path, producer = producer)
  if (identical(record_info[["state"]], "invalid") ||
      identical(record_info[["state"]], "absent")) {
    # A run-bound cache without its expected record is not reusable.  Legacy
    # artifacts are handled only when no run/producer record context exists.
    return(FALSE)
  }
  if (identical(record_info[["state"]], "valid")) {
    return(isTRUE(ecoda_local_record_only_valid(path, record_info)))
  }
  # No record context is the explicit legacy path: strict sidecar/content
  # checking, never existence-only reuse.
  isTRUE(ecoda_local_strict_checksum(path))
}

ecoda_local_read_rds <- function(path, description = "RDS artifact",
                                 producer = NULL) {
  record_info <- ecoda_local_read_record(path, producer = producer)
  if (identical(record_info[["state"]], "invalid") ||
      identical(record_info[["state"]], "absent")) {
    stop(description, " artifact record is missing or malformed: ", path)
  }
  if (identical(record_info[["state"]], "valid") &&
      !ecoda_local_record_only_valid(path, record_info)) {
    stop(description, " artifact record does not match current file: ", path)
  }
  # This full checksum is intentionally the final check before deserialization.
  if (!ecoda_local_strict_checksum(path, description)) {
    stop(description, " checksum/content validation failed: ", path)
  }
  readRDS(path)
}


ecoda_local_write_record <- function(
  path, producer, run_id, checksum = NULL
) {
  run_id <- ecoda_local_validate_run_id(run_id)
  if (!is.character(producer) || length(producer) != 1L ||
      is.na(producer) || !nzchar(producer) ||
      grepl("[\r\n]", producer, perl = TRUE)) {
    stop("artifact producer must be one non-empty line")
  }
  if (is.null(checksum)) {
    if (!ecoda_local_strict_checksum(path)) {
      stop("cannot publish an artifact without a strict checksum: ", path)
    }
    checksum <- ecoda_local_checksum_fields(path)
  }
  if (is.null(checksum)) stop("artifact checksum fields are unavailable: ", path)
  record_path <- ecoda_local_record_path(path, run_id)
  if (is.null(record_path)) {
    stop("artifact-record root is unavailable for ", path)
  }
  dir.create(dirname(record_path), showWarnings = FALSE, recursive = TRUE)
  temporary <- paste0(record_path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(c(
    paste0("PATH=", ecoda_local_canonical_path(path)),
    paste0("SIZE=", checksum[["SIZE"]]),
    paste0("MD5=", checksum[["MD5"]]),
    paste0("RUN_ID=", run_id),
    paste0("PRODUCER=", producer),
    "STATE=PUBLISHED"
  ), temporary, useBytes = TRUE)
  if (!file.rename(temporary, record_path)) {
    stop("could not atomically publish artifact record: ", record_path)
  }
  invisible(record_path)
}

ecoda_local_publish_rds <- function(object, path, producer) {
  save_rds_atomic(object, path)
  run_id <- ecoda_local_current_run_id()
  if (!is.null(run_id)) {
    if (!ecoda_local_strict_checksum(path)) {
      stop("new RDS artifact failed strict publication checksum: ", path)
    }
    checksum <- ecoda_local_checksum_fields(path)
    ecoda_local_write_record(path, producer, run_id, checksum = checksum)
  }
  invisible(NULL)
}


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
if (!is.null(analysis_pass) && !analysis_pass %in% c("uncorrected", "corrected")) {
  stop("Unknown analysis pass: ", analysis_pass)
}
if (combo_supplied && !is.null(analysis_pass)) {
  stop("--combo is only supported for ordinary GloScope runs")
}
cache_stem <- if (is.null(analysis_pass)) {
  ds
} else {
  paste0(ds, "_batch_effect_", analysis_pass)
}
entry <- config[[ds]]
if (is.null(entry)) {
  stop("Dataset '", ds, "' not found in ", args$config_path)
}

h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}
dir.create(args$results_dir, showWarnings = FALSE, recursive = TRUE)

method_rds_stem <- if (is.null(analysis_pass)) ds else cache_stem
method_rds <- file.path(
  args$results_dir,
  paste0(method_rds_stem, "_", method, ".rds")
)
if (!combo_supplied && ecoda_local_cache_valid(method_rds) && !force) {
  message("Method results already exist and passed checksum/record validation: ", method_rds)
  cached <- ecoda_local_read_rds(method_rds, "Method results")
  for (nm in names(cached)) {
    if (!is.null(cached[[nm]]$exec_time)) {
      log_exec_row(ds, nm, cached[[nm]]$exec_time, args$log_file,
                   mem_gb = cached[[nm]]$mem_GB)
    }
  }
  quit(save = "no", status = 0)
}

counts_free_method <- method %in% c("gloscope", "composition", "mofa") ||
  (method == "pseudobulk" && !is.null(analysis_pass))
embedding_key <- if (args$view == "batch_effect_corrected") {
  "X_pca_harmony_batch_effect_corrected_hvg2000"
} else if (args$view == "batch_effect_uncorrected") {
  "X_pca_batch_effect_uncorrected_hvg2000"
} else {
  "X_pca_benchmark_analysis_hvg2000"
}
batch_col <- if (!is.null(analysis_pass) && analysis_pass == "corrected") {
  entry$batch_col
} else {
  NULL
}
sample_col <- "Sample"

if (counts_free_method) {
  embedding_keys <- if (method == "gloscope" && is.null(analysis_pass)) {
    c(
      "X_pca_benchmark_analysis_hvg1000",
      "X_pca_benchmark_analysis_hvg2000",
      "X_pca_benchmark_analysis_hvg3000"
    )
  } else {
    embedding_key
  }
  composition_obs_columns <- if (
    method == "composition" &&
    is.null(analysis_pass) &&
    length(entry$not_suitable_for_auto_annotation) == 0
  ) {
    c("layer2", "scATOMIC_pred")
  } else {
    character()
  }
  obs_columns <- if (method %in% c("gloscope", "mofa", "pseudobulk")) {
    c("Sample", entry$label_col, batch_col)
  } else {
    c(
      "Sample",
      entry$label_col,
      entry$cell_type_low_res,
      entry$cell_type_high_res,
      batch_col,
      composition_obs_columns
    )
  }
  adata <- load_h5ad_counts_free(
    h5ad_path,
    unique(obs_columns[!is.na(obs_columns) & nzchar(obs_columns)]),
    embedding_keys,
    obs_prefixes = if (method == "composition") "leiden_res_" else character(),
    view = args$view,
    method = method
  )
  if (method %in% c("mofa", "pseudobulk")) {
    obs <- load_h5ad_sample_metadata(
      h5ad_path,
      sample_col = sample_col,
      metadata_columns = c(entry$label_col, batch_col)
    )
  } else {
    obs <- py_to_r(adata$obs)
  }
} else {
  ad <- import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(h5ad_path, backed = "r")
  obs <- py_to_r(adata$obs)
  validate_benchmark_h5ad_contract(
    adata,
    obs = obs,
    view = args$view,
    method = method
  )
}

if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}
blind_mode <- is.null(analysis_pass) || analysis_pass == "uncorrected"
correct_batch_mode <- identical(analysis_pass, "corrected")
hvg_rank_genes <- get_hvg_rank_genes(adata)
if (correct_batch_mode) {
  if (is.null(batch_col)) {
    stop("corrected batch-effect view requires a confirmed columns.batch")
  }
  if (!batch_col %in% colnames(obs)) {
    stop("Confirmed batch column '", batch_col, "' not found in obs of ", h5ad_path)
  }
}

pb_variants <- NULL
seurat <- NULL
metadata <- NULL
labels <- NULL

if (method == "mofa") {
  # MOFA consumes only the precomputed pseudobulks. Missing variants are
  # rebuilt through bounded H5AD CSR aggregation, never a full cell Seurat.
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method mofa")
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(pb_variants_missing(
    args$pseudobulk_dir, ds, force, cache_stem = cache_stem
  )) > 0) {
    message("Building sample-level Seurat for on-the-fly pseudobulk variants...")
    seurat <- load_h5ad_pseudobulk_seurat(
      h5ad_path,
      sample_col = sample_col,
      batch_col = batch_col
    )
  }
  pb_variants <- load_pb_variants(
    seurat, sample_col, hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir, ds = ds,
    force = force, log_file = args$log_file, cache_stem = cache_stem,
    batch_col = batch_col, blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = if (!is.null(analysis_pass)) "hvg2000" else PB_VARIANT_NAMES
  )
  # Sample names are already standardized in the preprocessed obs
  # (1.1.1_preprocess.py): no standardize_sample_names() re-application here
  # (it would diverge the labels from the obs names for h5ads that predate
  # the python change, e.g. Adams).
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
} else if (method == "pseudobulk" && !is.null(analysis_pass)) {
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method pseudobulk")
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  if (length(pb_variants_missing(
    args$pseudobulk_dir, ds, force, cache_stem = cache_stem
  )) > 0) {
    message("Building sample-level Seurat for missing pseudobulk variants...")
    seurat <- load_h5ad_pseudobulk_seurat(
      h5ad_path,
      sample_col = sample_col,
      batch_col = batch_col
    )
  }
  pb_variants <- load_pb_variants(
    seurat, sample_col, hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir, ds = ds,
    force = force, log_file = args$log_file, cache_stem = cache_stem,
    batch_col = batch_col, blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = "hvg2000"
  )
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
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
    variants = if (!is.null(analysis_pass)) "hvg2000" else PB_VARIANT_NAMES
  )
} else {
  seurat <- load_benchmark_seurat(
    adata, obs, sample_col = sample_col,
    fetch_embedding = if (method == "gloscope") {
      if (!is.null(analysis_pass)) {
        embedding_key
      } else {
        c("X_pca_benchmark_analysis_hvg1000",
          "X_pca_benchmark_analysis_hvg2000",
          "X_pca_benchmark_analysis_hvg3000")
      }
    } else {
      NULL
    },
    counts_layer = if (method == "gloscope") NULL else "counts"
  )
  # Sample names are already standardized in the preprocessed obs
  # (1.1.1_preprocess.py): no standardize_sample_names() re-application
  # (kept only in the legacy Seurat path of run_benchmark_analysis).
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

  if (method == "pseudobulk") {
    if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
      stop("Missing required --pseudobulk_dir argument for method pseudobulk")
    }
    dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
    pb_variants <- load_pb_variants(
      seurat, sample_col, hvg_rank_genes,
      pseudobulk_dir = args$pseudobulk_dir, ds = ds,
      force = force, log_file = args$log_file, cache_stem = cache_stem,
      batch_col = batch_col, blind = blind_mode,
      correct_batch = correct_batch_mode
    )
  }
  if (method == "gloscope") {
    if (is.null(args$gloscope_cache_dir) ||
        identical(args$gloscope_cache_dir, TRUE)) {
      stop("Missing required --gloscope_cache_dir argument for method gloscope")
    }
    dir.create(args$gloscope_cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

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
    }
  ),
  mofa = run_mofa_hpc(
    metadata, labels, pb_variants,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  pseudobulk = run_pseudobulk_hpc(
    seurat, labels, pb_variants,
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file,
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem
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
    not_suitable_for_auto_annotation = if (is.null(entry$not_suitable_for_auto_annotation)) {
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
  ecoda_local_publish_rds(
    results, method_rds, producer = paste0("stage5_", method)
  )
  message("Saved method results: ", method_rds, " (", length(results), " combos)")
}
message("--- ", method, " for ", ds, " complete ---")
