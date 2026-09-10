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
sample_col <- "Sample"
seurat <- load_h5ad_pseudobulk_seurat(
  h5ad_path,
  sample_col = sample_col,
  batch_col = batch_col
)
obs <- seurat@meta.data
validate_benchmark_h5ad_contract(
  adata,
  obs = obs,
  view = args$view,
  method = "prepare_pseudobulk"
)
if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}

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
  vapply(
    file.path(
      args$pseudobulk_dir,
      paste0(cache_stem, "_pseudobulk_", requested_variants, ".rds")
    ),
    function(path) !ecoda_local_cache_valid(path),
    logical(1)
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
    ecoda_local_publish_rds(
      variants[[v]], f, producer = paste0("stage5_prepare_pseudobulk_", v)
    )
    log_exec_row(ds, paste0("prepare_pseudobulk_", v),
                 variants[[v]]$time_secs, args$log_file,
                 mem_gb = variants[[v]]$mem_GB)
    message("  Saved: ", f, " (", round(variants[[v]]$time_secs, 1), "s)")
  }
} else {
  # Everything requested is cached: re-emit stored timings on resume.
  for (v in requested_variants) {
    f <- file.path(args$pseudobulk_dir, paste0(cache_stem, "_pseudobulk_", v, ".rds"))
    cached <- ecoda_local_read_rds(f, "Pseudobulk variant")
    log_exec_row(ds, paste0("prepare_pseudobulk_", v),
                 cached$time_secs, args$log_file,
                 mem_gb = cached$mem_GB)
    message("Pseudobulk variant already exists: ", f)
  }
}

message("--- prepare_pseudobulk for ", ds, " complete ---")
