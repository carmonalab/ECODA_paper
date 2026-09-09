suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
})

root <- normalizePath(".", mustWork = TRUE)
embedding_dir <- file.path(root, "data", "benchmark", "embeddings")
result_dir <- file.path(root, "data", "benchmark", "results")
baseline_path <- file.path(embedding_dir, "execution_times_OLD.feather")
candidate_backup <- file.path(
  embedding_dir,
  "execution_times.candidate_before_mem_repair_1788900485.feather"
)
staged_path <- file.path(
  embedding_dir,
  "execution_times.mem_repaired_baseline_1788900485.feather"
)
evidence_path <- file.path(
  embedding_dir,
  "execution_times.mem_repair_1788900485.json"
)

expected_baseline_md5 <- "4df267b7492cce9937db26a717d3aa24"
expected_candidate_md5 <- "3d9cdb84633e02e52606feb3b0ffccf4"

md5 <- function(path) unname(tools::md5sum(path))[[1L]]
require_file <- function(path, description) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop(description, " is missing or empty: ", path)
  }
}
require_file(baseline_path, "authoritative baseline")
require_file(candidate_backup, "candidate backup")
if (!identical(md5(baseline_path), expected_baseline_md5)) {
  stop("Authoritative baseline MD5 changed: ", md5(baseline_path))
}
if (!identical(md5(candidate_backup), expected_candidate_md5)) {
  stop("Candidate backup MD5 changed: ", md5(candidate_backup))
}
if (file.exists(staged_path) || file.exists(evidence_path)) {
  stop("Repair target already exists; refusing to overwrite: ", staged_path)
}

baseline <- as.data.frame(read_feather(baseline_path))
required_columns <- c("dataset", "method", "time_secs", "mem_GB")
if (!identical(names(baseline), required_columns)) {
  stop("Baseline schema mismatch: ", paste(names(baseline), collapse = ", "))
}
if (nrow(baseline) != 1376L) stop("Unexpected baseline row count: ", nrow(baseline))
if (anyNA(baseline$dataset) || anyNA(baseline$method) ||
    any(!nzchar(trimws(as.character(baseline$dataset)))) ||
    any(!nzchar(trimws(as.character(baseline$method)))) ||
    anyDuplicated(paste(baseline$dataset, baseline$method, sep = "\r"))) {
  stop("Baseline identifiers are invalid or duplicated")
}
if (any(!is.finite(baseline$time_secs) | baseline$time_secs < 0)) {
  stop("Baseline contains invalid time_secs")
}
if (any(!is.na(baseline$mem_GB) &
        (!is.finite(baseline$mem_GB) | baseline$mem_GB < 0))) {
  stop("Baseline contains invalid non-NA mem_GB")
}

# The validated source is the checksum-protected method-level RDS. Each file
# contains named per-combo result bundles; only finite mem_GB fields inside a
# validated container can be used. Standalone per-combo files are not used.
methods <- c("mofa", "gloscope")
source_files <- unlist(lapply(methods, function(method) {
  list.files(
    result_dir,
    pattern = paste0("_", method, "[.]rds$"),
    full.names = TRUE
  )
}))
source_files <- source_files[!startsWith(basename(source_files), "_")]
if (length(source_files) == 0L) stop("No production method bundles found")

source_rows <- list()
source_file_records <- list()
for (path in sort(source_files)) {
  method <- sub(".*_(mofa|gloscope)[.]rds$", "\\1", basename(path))
  dataset <- sub(paste0("_", method, "[.]rds$"), "", basename(path))
  require_file(path, "source method bundle")
  sidecar_path <- paste0(path, ".md5")
  require_file(sidecar_path, "source bundle checksum sidecar")
  sidecar_lines <- readLines(sidecar_path, warn = FALSE)
  sidecar_keys <- sub("=.*$", "", sidecar_lines)
  sidecar_values <- sub("^[^=]*=", "", sidecar_lines)
  sidecar <- setNames(sidecar_values, sidecar_keys)
  if (!identical(sidecar[["MD5"]], md5(path)) ||
      !identical(sidecar[["SIZE"]], as.character(file.info(path)$size))) {
    stop("Source bundle checksum mismatch: ", path)
  }

  bundles <- readRDS(path)
  if (!is.list(bundles) || is.null(names(bundles)) || any(!nzchar(names(bundles)))) {
    stop("Source bundle has invalid named-list schema: ", path)
  }
  source_file_records[[length(source_file_records) + 1L]] <- list(
    path = path,
    dataset = dataset,
    method = method,
    md5 = md5(path),
    size = unname(file.info(path)$size),
    bundle_count = length(bundles)
  )
  for (combo in names(bundles)) {
    bundle <- bundles[[combo]]
    if (!is.list(bundle) || !"exec_time" %in% names(bundle)) {
      stop("Source bundle combo lacks exec_time: ", path, " / ", combo)
    }
    bundle_time <- as.numeric(bundle[["exec_time"]])
    if (length(bundle_time) != 1L || is.na(bundle_time) ||
        !is.finite(bundle_time) || bundle_time < 0) {
      stop("Source bundle combo has invalid exec_time: ", path, " / ", combo)
    }
    if (!"mem_GB" %in% names(bundle) || is.null(bundle[["mem_GB"]])) next
    bundle_mem <- as.numeric(bundle[["mem_GB"]])
    if (length(bundle_mem) != 1L || is.na(bundle_mem)) next
    if (!is.finite(bundle_mem) || bundle_mem < 0) {
      stop("Source bundle combo has invalid mem_GB: ", path, " / ", combo)
    }
    source_rows[[length(source_rows) + 1L]] <- data.frame(
      dataset = dataset,
      method = combo,
      bundle_time_secs = bundle_time,
      bundle_mem_GB = bundle_mem,
      source_file = path,
      source_md5 = md5(path),
      stringsAsFactors = FALSE
    )
  }
}
source <- do.call(rbind, source_rows)
source_key <- paste(source$dataset, source$method, sep = "\r")
if (anyDuplicated(source_key)) stop("Duplicate source bundle keys")

baseline_key <- paste(baseline$dataset, baseline$method, sep = "\r")
source_index <- match(baseline_key, source_key)
fillable <- is.na(baseline$mem_GB) & !is.na(source_index)
if (!any(fillable)) stop("No missing baseline mem_GB values have validated sources")

repaired <- baseline
repaired$mem_GB[fillable] <- source$bundle_mem_GB[source_index[fillable]]
changed <- data.frame(
  dataset = repaired$dataset[fillable],
  method = repaired$method[fillable],
  baseline_time_secs = baseline$time_secs[fillable],
  source_time_secs = source$bundle_time_secs[source_index[fillable]],
  mem_GB = repaired$mem_GB[fillable],
  source_file = source$source_file[source_index[fillable]],
  source_md5 = source$source_md5[source_index[fillable]],
  stringsAsFactors = FALSE
)

if (!identical(names(repaired), required_columns) || nrow(repaired) != nrow(baseline)) {
  stop("Repaired schema or row count changed")
}
if (!identical(repaired$dataset, baseline$dataset) ||
    !identical(repaired$method, baseline$method) ||
    !identical(repaired$time_secs, baseline$time_secs)) {
  stop("Repair changed baseline identifiers or runtimes")
}
if (any(!is.na(repaired$mem_GB) &
        (!is.finite(repaired$mem_GB) | repaired$mem_GB < 0))) {
  stop("Repaired table contains invalid mem_GB")
}

# Stage the output atomically on the same filesystem. Installation into the
# canonical path is a separate, auditable step after read-back validation.
tmp_path <- file.path(
  embedding_dir,
  paste0(".execution_times.mem_repaired_baseline.tmp.", Sys.getpid())
)
on.exit(unlink(tmp_path), add = TRUE)
arrow::write_feather(repaired, tmp_path)
require_file(tmp_path, "staged repaired Feather")
if (!file.rename(tmp_path, staged_path)) stop("Could not install staged repair: ", staged_path)
read_back <- as.data.frame(read_feather(staged_path))
if (!identical(names(read_back), required_columns) ||
    nrow(read_back) != nrow(repaired) ||
    !identical(read_back$dataset, repaired$dataset) ||
    !identical(read_back$method, repaired$method) ||
    !identical(read_back$time_secs, repaired$time_secs) ||
    !identical(read_back$mem_GB, repaired$mem_GB)) {
  stop("Staged repaired Feather failed read-back validation")
}

# Emit evidence before canonical installation; the installer will update the
# status and canonical output hash after the staged file is accepted.
evidence <- list(
  schema_version = 1L,
  status = "STAGED_FROM_AUTHORITATIVE_BASELINE",
  baseline = list(
    path = baseline_path,
    md5 = md5(baseline_path),
    size = unname(file.info(baseline_path)$size),
    rows = nrow(baseline)
  ),
  candidate_backup = list(
    path = candidate_backup,
    md5 = md5(candidate_backup),
    size = unname(file.info(candidate_backup)$size),
    preserved = TRUE
  ),
  staged_output = list(
    path = staged_path,
    md5 = md5(staged_path),
    size = unname(file.info(staged_path)$size),
    rows = nrow(repaired)
  ),
  source_bundle_files = source_file_records,
  changed_rows = unname(split(changed, seq_len(nrow(changed)))),
  changed_row_count = nrow(changed),
  source_time_mismatch_note = "Bundle memory is copied only into existing baseline rows; baseline time_secs are unchanged. Bundle and baseline times may represent different artifact generations and are recorded above.",
  canonical_installation = "PENDING"
)
json_tmp <- paste0(evidence_path, ".tmp.", Sys.getpid())
write_json(evidence, json_tmp, auto_unbox = TRUE, pretty = TRUE, na = "null")
if (!file.rename(json_tmp, evidence_path)) stop("Could not install repair evidence: ", evidence_path)

cat("STAGED", staged_path, "\n")
cat("EVIDENCE", evidence_path, "\n")
cat("CHANGED_ROWS", nrow(changed), "\n")
print(as.data.frame(table(changed$method)), row.names = FALSE)
