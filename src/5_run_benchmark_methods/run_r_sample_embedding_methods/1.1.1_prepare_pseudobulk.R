# ==============================================================================
# 1.1.1_prepare_pseudobulk.R — Precompute the shared DESeq2 pseudobulks for
# one dataset (prepare_pseudobulk array of Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --input_dir
#   --pseudobulk_dir --log_file [--force]
# Loads the preprocessed benchmark view h5ad (raw counts + var["hvg_rank"]
# only; no embeddings), runs prepare_pseudobulks_hpc() and writes
# pseudobulks/<ds>_pseudobulk_<variant>.rds atomically (list(pb, time_secs)),
# with one exec-log row per variant (method "prepare_pseudobulk_<variant>").
# Skip-if-exists per variant unless --force.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

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
if (!is.null(analysis_pass) && !analysis_pass %in% c("uncorrected", "corrected")) {
  stop("Unknown analysis pass: ", analysis_pass)
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
batch_col <- if (!is.null(analysis_pass) && analysis_pass == "corrected") {
  if (is.null(entry$batch_col)) {
    stop("corrected batch-effect view requires a confirmed columns.batch")
  }
  entry$batch_col
} else {
  NULL
}
blind_mode <- is.null(analysis_pass) || analysis_pass == "uncorrected"
correct_batch_mode <- identical(analysis_pass, "corrected")

h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}
dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)

ad <- import("anndata", convert = FALSE)
adata <- ad$read_h5ad(h5ad_path, backed = "r")
obs <- py_to_r(adata$obs)

sample_col <- "Sample"
if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}

seurat <- load_benchmark_seurat(adata, obs, sample_col = sample_col,
                                fetch_embedding = NULL)
# Sample names are already standardized in the preprocessed obs
# (1.1.1_preprocess.py): do NOT re-apply standardize_sample_names() here —
# it would diverge (hyphen -> underscore) from the obs names for h5ads that
# predate the python change (e.g. Adams), breaking the bundle label match.
hvg_rank_genes <- get_hvg_rank_genes(adata)

requested_variants <- if (is.null(analysis_pass)) {
  PB_VARIANT_NAMES
} else {
  "hvg2000"
}
pending <- requested_variants[
  !file.exists(
    file.path(
      args$pseudobulk_dir,
      paste0(cache_stem, "_pseudobulk_", requested_variants, ".rds")
    )
  ) | force
]

if (length(pending) > 0) {
  message("Computing pseudobulk variants: ", paste(pending, collapse = ", "))
  variants <- prepare_pseudobulks_hpc(
    seurat,
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    variants = pending,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode
  )
  for (v in names(variants)) {
    f <- file.path(args$pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds"))
    save_rds_atomic(variants[[v]], f)
    log_exec_row(ds, paste0("prepare_pseudobulk_", v),
                 variants[[v]]$time_secs, args$log_file,
                 mem_gb = variants[[v]]$mem_GB)
    message("  Saved: ", f, " (", round(variants[[v]]$time_secs, 1), "s)")
  }
} else {
  # Everything requested is cached: re-emit stored timings on resume.
  for (v in requested_variants) {
    f <- file.path(args$pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds"))
    cached <- readRDS(f)
    log_exec_row(ds, paste0("prepare_pseudobulk_", v),
                 cached$time_secs, args$log_file,
                 mem_gb = cached$mem_GB)
    message("Pseudobulk variant already exists: ", f)
  }
}

message("--- prepare_pseudobulk for ", ds, " complete ---")
