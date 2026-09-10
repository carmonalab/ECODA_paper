#!/usr/bin/env Rscript
# Consolidate run-owned GloScope parameter shards into the canonical method RDS.
# The shard workers write one validated bundle per dataset/parameter; this
# serialized tail runs only after the matrix watchdog has gated every shard.

project_root <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT is required")
}
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


args <- parse_flags(commandArgs(trailingOnly = TRUE))
for (required in c("manifest", "results_dir")) {
  if (is.null(args[[required]]) || identical(args[[required]], TRUE) || !nzchar(args[[required]])) {
    stop("Missing required --", required, " argument")
  }
}

manifest <- normalizePath(args$manifest, mustWork = TRUE)
# Preserve the configured path spelling: normalizePath() resolves symlinks,
# while save_rds_atomic() records that spelling in each checksum sidecar.
results_dir <- path.expand(args$results_dir)
if (!dir.exists(results_dir)) {
  stop("Results directory does not exist: ", results_dir)
}
rows <- readLines(manifest, warn = FALSE)
if (!length(rows) || any(!nzchar(rows))) {
  stop("GloScope shard manifest is empty or contains blank rows: ", manifest)
}
fields <- strsplit(rows, "\t", fixed = TRUE)
if (any(lengths(fields) != 4L)) {
  stop("GloScope shard manifest must have four columns: ", manifest)
}

records <- do.call(rbind, lapply(fields, function(values) {
  data.frame(
    dataset = values[[1L]],
    view = values[[2L]],
    method = values[[3L]],
    combo = values[[4L]],
    stringsAsFactors = FALSE
  )
}))
if (any(!nzchar(records$dataset)) || any(records$view != "benchmark_analysis") ||
    any(records$method != "gloscope") || any(!nzchar(records$combo))) {
  stop("GloScope shard manifest contains an invalid row: ", manifest)
}
if (anyDuplicated(paste(records$dataset, records$combo, sep = "\t"))) {
  stop("GloScope shard manifest contains duplicate dataset/combo rows: ", manifest)
}

expected_combos <- c(
  "hvg2000_pcadims10",
  "hvg2000_pcadims30",
  "hvg2000_pcadims50",
  "hvg1000_pcadims30",
  "hvg3000_pcadims30"
)
if (any(!records$combo %in% expected_combos)) {
  stop("GloScope shard manifest contains an unknown parameter combo: ", manifest)
}

for (dataset in unique(records$dataset)) {
  selected <- records[records$dataset == dataset, , drop = FALSE]
  if (!setequal(selected$combo, expected_combos) || nrow(selected) != length(expected_combos)) {
    stop("GloScope shard coverage is incomplete for ", dataset)
  }
  selected <- selected[match(expected_combos, selected$combo), , drop = FALSE]
  combo_names <- paste0("GloScope_", selected$combo)
  shard_paths <- file.path(results_dir, paste0(dataset, "_", combo_names, ".rds"))
  shard_valid <- vapply(shard_paths, ecoda_local_cache_valid, logical(1))
  if (!all(shard_valid)) {
    bad <- shard_paths[!shard_valid]
    stop("GloScope shard checksum/record validation failed: ",
         paste(bad, collapse = ", "))
  }
  bundles <- lapply(
    shard_paths,
    function(path) ecoda_local_read_rds(path, "GloScope shard")
  )
  names(bundles) <- combo_names
  method_path <- file.path(results_dir, paste0(dataset, "_gloscope.rds"))
  ecoda_local_publish_rds(
    bundles, method_path, producer = "stage5_gloscope_consolidate"
  )
  message("Consolidated ", length(bundles), " GloScope shards for ", dataset,
          " -> ", method_path)
}
