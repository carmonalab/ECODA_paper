# ==============================================================================
# 1.1.1_prepare_pseudobulk.R — Precompute the shared DESeq2 pseudobulks for
# one dataset (prepare_pseudobulk array of Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --input_dir
#   --pseudobulk_dir --log_file [--force]
# Loads H5AD metadata/HVG ranks without count values, validates the pending
# cache set, then performs one raw Sample aggregation and one shared full-gene
# DESeq2 fit where possible.  Each
# pseudobulks/<ds>_pseudobulk_<variant>.rds record is published atomically
# with schema-2 shared/variant timing fields and one shared timing log row.
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

for (req in c("config_path", "ds_name", "view", "input_dir",
              "pseudobulk_dir", "log_file")) {
  if (is.null(args[[req]]) || identical(args[[req]], TRUE)) {
    stop("Missing required --", req, " argument")
  }
}
force <- isTRUE(args[["force"]]) || identical(args[["force"]], "TRUE")

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
sample_col <- "Sample"
batch_keys <- NULL
batch_context <- NULL
effective_batch_keys <- NULL
pseudobulk_batch_contract <- NULL
h5ad_expected_batch_contract <- NULL
python_batch_metadata <- NULL
batch_col <- NULL
if (!is.null(analysis_pass) && analysis_pass == "corrected") {
  if (is.null(entry$batch_col)) {
    stop("corrected batch-effect view requires a confirmed columns.batch")
  }
  # Keep the configured technical columns separate.  The corrected core
  # derives estimable keys from batch_context; batch_col is only a direct
  # original column for the one-key compatibility boundary.
  batch_keys <- ecoda_hpc_normalize_keys(
    entry$batch_col,
    sample_col = sample_col,
    biological_label = entry$label_col
  )
  if (length(batch_keys) == 1L) batch_col <- batch_keys[[1L]]
}
blind_mode <- is.null(analysis_pass) || analysis_pass == "uncorrected"
correct_batch_mode <- identical(analysis_pass, "corrected")
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
if (correct_batch_mode) {
  h5ad_expected_batch_contract <- ecoda_hpc_batch_contract_identity(
    batch_keys = batch_keys,
    sample_col = sample_col,
    method_id = "preprocess",
    model_id = "hvg_composite_v1"
  )
}
h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}
if (correct_batch_mode && !corrected_final_mode) {
  # Validate every cell before any metadata reducer selects a first row per
  # Sample. This is the historical ordinary corrected contract.
  python_batch_metadata <- validate_h5ad_corrected_batch_metadata(
    h5ad_path = h5ad_path,
    batch_keys = as.list(unname(batch_keys)),
    sample_col = sample_col,
    biological_label = entry$label_col,
    method_id = "preprocess",
    model_id = "hvg_composite_v1"
  )
}
dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)

# Read metadata and ranked HVGs through h5py-only helpers.  These readers
# validate persisted H5AD structure while never materializing count values;
# cache completeness is checked below before the raw CSR pass is requested.
required_hvg <- if (identical(args$view, "benchmark_analysis")) 3000L else 2000L
h5ad_metadata <- load_h5ad_pseudobulk_metadata(
  h5ad_path,
  sample_col = sample_col,
  metadata_columns = if (correct_batch_mode) {
    unique(c(as.list(unname(batch_keys)), entry$label_col))
  } else {
    batch_col
  },
  n_hvg = required_hvg,
  required_nonmissing_columns = unique(c(
    sample_col,
    if (correct_batch_mode) c(batch_keys, entry$label_col) else batch_col
  )),
  expected_batch_contract = h5ad_expected_batch_contract,
  view = args$view,
  method = "preprocessing",
  allow_missing_summary = corrected_final_mode
)
obs <- h5ad_metadata$obs
hvg_rank_genes <- h5ad_metadata$hvg_rank_genes
if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}
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
    stop("Corrected pseudobulk effective keys are not configured technical columns")
  }
  # A constant key remains in metadata but never enters the correction design.
  # Multiple effective keys are supplied through batch_context as separate
  # original columns; no synthetic/composite column is created.
  batch_col <- if (length(effective_batch_keys) == 1L) {
    effective_batch_keys[[1L]]
  } else {
    NULL
  }
  pseudobulk_batch_contract <- ecoda_hpc_augment_batch_contract(
    identity = ecoda_hpc_batch_contract_identity(
      batch_keys = batch_keys,
      sample_col = sample_col,
      method_id = "Pseudobulk",
      model_id = "pseudobulk_limma_fixed_effects_v1"
    ),
    validation = batch_context$validation,
    method_id = "Pseudobulk",
    batch_keys = batch_keys
  )
}

# Sample names are already standardized in the preprocessed obs
# (1.1.1_preprocess.py): do NOT re-apply standardize_sample_names() here —
# it would diverge (hyphen -> underscore) from the obs names for h5ads that
# predate the python change (e.g. Adams), breaking the bundle label match.

requested_variants <- if (is.null(analysis_pass)) {
  PB_VARIANT_NAMES
} else {
  "hvg2000"
}
cache_paths <- file.path(
  args$pseudobulk_dir,
  paste0(cache_stem, "_pseudobulk_", requested_variants, ".rds")
)
cache_valid <- vapply(seq_along(requested_variants), function(index) {
  .pb_variant_cache_valid(
    cache_paths[[index]],
    requested_variants[[index]],
    expected_batch_contract = if (correct_batch_mode) {
      pseudobulk_batch_contract
    } else {
      NULL
    }
  )
}, logical(1))
cached_variants <- list()
reused <- requested_variants[cache_valid & !force]
for (v in reused) {
  cache_index <- match(v, requested_variants)
  value <- read_rds_checked(
    cache_paths[[cache_index]],
    producer = PB_VARIANT_PRODUCERS[[v]]
  )
  validate_pseudobulk_timing_record(value, variant = v)
  if (correct_batch_mode) {
    ecoda_hpc_validate_batch_contract(
      value[["batch_contract"]],
      pseudobulk_batch_contract,
      label = paste0("Pseudobulk cache ", v, " (", cache_paths[[cache_index]], ")")
    )
  }
  cached_variants[[v]] <- value
}
pending <- requested_variants[!cache_valid | force]

if (length(pending) > 0) {
  message("Computing pseudobulk variants: ", paste(pending, collapse = ", "))
  variants <- prepare_pseudobulks_hpc(
    h5ad_path = h5ad_path,
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    variants = pending,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode,
    cache_stem = cache_stem,
    view = args$view,
    analysis_pass = analysis_pass,
    run_id = Sys.getenv("ECODA_RUN_ID", unset = ""),
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_contract = if (correct_batch_mode) pseudobulk_batch_contract else NULL,
    expected_h5ad_batch_contract = h5ad_expected_batch_contract
  )
  for (v in names(variants)) {
    f <- file.path(
      args$pseudobulk_dir,
      paste0(cache_stem, "_pseudobulk_", v, ".rds")
    )
    save_rds_atomic(
      variants[[v]], f, producer = paste0("stage5_prepare_pseudobulk_", v)
    )
    message("  Saved: ", f, " (",
            round(variants[[v]]$time_secs, 1), "s)")
  }
  all_variants <- cached_variants
  for (v in names(variants)) {
    all_variants[[v]] <- variants[[v]]
  }
  all_variants <- all_variants[requested_variants]
  emit_pseudobulk_timing_rows(
    all_variants, ds, log_file = args$log_file
  )
} else {
  # Everything requested is cached: re-emit stored timings on resume without
  # opening the counts layer, aggregating, or constructing a Seurat object.
  for (v in requested_variants) {
    message(
      "Pseudobulk variant already exists: ",
      cache_paths[[match(v, requested_variants)]]
    )
  }
  emit_pseudobulk_timing_rows(
    cached_variants, ds, log_file = args$log_file
  )
}
