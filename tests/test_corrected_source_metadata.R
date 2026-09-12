#!/usr/bin/env Rscript

# Standalone, deterministic regression for the corrected-mode RDS metadata
# trust boundary.  Every input, identity file, and report is confined to a
# temporary tree; the auditor itself is invoked through the repository Pixi
# environment so this exercises the same command boundary used by callers.

raw_args <- commandArgs(trailingOnly = FALSE)
script_args <- raw_args[grepl("^--file=", raw_args)]
if (length(script_args) != 1L) stop("could not determine test script path", call. = FALSE)
script_path <- sub("^--file=", "", script_args[[1L]])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required", call. = FALSE)
}

main <- function() {
oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(root)

tmp <- tempfile("corrected-source-metadata-")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

snapshot <- file.path(tmp, "snapshot")
tree <- file.path(snapshot, "tree")
identity_dir <- file.path(snapshot, "identity")
runtime_dir <- file.path(snapshot, "runtime")
run_root <- file.path(tmp, "run-root")
for (directory in c(
  file.path(tree, "src", "utils"),
  file.path(snapshot, "aux"),
  identity_dir,
  runtime_dir,
  file.path(run_root, "manifests"),
  file.path(run_root, "reports"),
  file.path(tmp, "input"),
  file.path(tmp, "negative-input")
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

source_root <- normalizePath(tree, mustWork = TRUE)
run_root <- normalizePath(run_root, mustWork = TRUE)

# The auditor resolves its contract from the immutable source tree, rather
# than from the mutable checkout.  Copy only that source dependency into the
# synthetic snapshot.
contract_source <- file.path(root, "src", "utils", "batch_contract.R")
contract_snapshot <- file.path(tree, "src", "utils", "batch_contract.R")
if (!isTRUE(file.copy(contract_source, contract_snapshot, overwrite = TRUE))) {
  stop("could not copy batch contract into fixture snapshot", call. = FALSE)
}

config_path <- file.path(tree, "datasets.json")
config <- list(
  Fixture = list(
    display_name = "Corrected metadata fixture",
    use_for_batch_effect = TRUE,
    columns = list(
      sample = "sample_id",
      label = "label",
      batch = "batch"
    ),
    views = list(
      batch_effect_corrected = list(
        input_file_name = "Fixture.rds",
        output_file_name = "Fixture_batch_effect_corrected_ECODAprocessed.h5ad",
        subset_vars = list(
          subset_group = list(values = "keep", op = "in")
        )
      )
    )
  )
)
jsonlite::write_json(config, config_path, auto_unbox = TRUE, pretty = TRUE)

# A repeated-sample data.frame is the supported RDS source.  Four samples
# are retained by the membership subset and one is dropped; the retained rows
# still contain both batch levels and have an estimable corrected design.
metadata <- data.frame(
  sample_id = rep(c("S1", "S2", "S3", "S4", "S5"), each = 2L),
  label = rep(c("case", "case", "control", "control", "case"), each = 2L),
  batch = rep(c("batch_a", "batch_a", "batch_a", "batch_a", "batch_b", "batch_b", "batch_b", "batch_b", "batch_b", "batch_b"), 1L),
  subset_group = rep(c("keep", "keep", "keep", "keep", "drop"), each = 2L),
  stringsAsFactors = FALSE,
  row.names = paste0("cell", seq_len(10L))
)
source_rds <- file.path(tmp, "input", "Fixture.rds")
saveRDS(metadata, source_rds)
source_rds <- normalizePath(source_rds, mustWork = TRUE)

# Bind the source archive and runtime files to their actual bytes.  The
# auditor checks all of these fields before it deserializes the RDS.
source_archive <- file.path(snapshot, "source.archive")
writeLines("fixture immutable source archive", source_archive, useBytes = TRUE)
source_archive <- normalizePath(source_archive, mustWork = TRUE)

pixi_toml <- file.path(tree, "pixi.toml")
pixi_lock <- file.path(tree, "pixi.lock")
writeLines("[workspace]\nname = 'fixture'", pixi_toml, useBytes = TRUE)
writeLines("fixture lock", pixi_lock, useBytes = TRUE)

sha256 <- function(path) {
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}
source_manifest <- file.path(identity_dir, "source.manifest")
writeLines(
  c(
    "FORMAT=1",
    paste0("SOURCE_ROOT=", source_root),
    "SOURCE_COMMIT=0123456789abcdef0123456789abcdef01234567",
    paste0("SOURCE_ARCHIVE_PATH=", source_archive),
    paste0("SOURCE_ARCHIVE_SHA256=", sha256(source_archive)),
    paste0("CONFIG_HELPER_SHA256=", sha256(contract_snapshot)),
    paste0("DATASETS_SHA256=", sha256(config_path)),
    paste0("PIXI_TOML_SHA256=", sha256(pixi_toml)),
    paste0("PIXI_LOCK_SHA256=", sha256(pixi_lock)),
    paste0("AUX_ROOT=", file.path(snapshot, "aux")),
    "SCGATE_DB_BRANCH=fixture"
  ),
  source_manifest,
  useBytes = TRUE
)
source_manifest <- normalizePath(source_manifest, mustWork = TRUE)

runtime_image <- file.path(runtime_dir, "fixture-runtime-image")
runtime_manifest <- file.path(runtime_dir, "fixture-runtime.manifest")
writeLines("fixture runtime image", runtime_image, useBytes = TRUE)
writeLines(c("FORMAT=1", "PROFILE=fixture"), runtime_manifest, useBytes = TRUE)
runtime_image <- normalizePath(runtime_image, mustWork = TRUE)
runtime_manifest <- normalizePath(runtime_manifest, mustWork = TRUE)

runtime_identity <- file.path(run_root, "manifests", "runtime.identity")
writeLines(
  c(
    paste0("RUNTIME_IMAGE=", runtime_image),
    paste0("RUNTIME_MANIFEST=", runtime_manifest),
    paste0("RUNTIME_IMAGE_SHA256=", sha256(runtime_image)),
    paste0("RUNTIME_MANIFEST_SHA256=", sha256(runtime_manifest)),
    paste0("RUNTIME_IMAGE_SIZE=", as.character(file.info(runtime_image)$size)),
    paste0("RUNTIME_MANIFEST_SIZE=", as.character(file.info(runtime_manifest)$size)),
    "RUNTIME_PROFILE=fixture"
  ),
  runtime_identity,
  useBytes = TRUE
)
runtime_identity <- normalizePath(runtime_identity, mustWork = TRUE)

run_auditor <- function(input_path, output_path) {
  output <- suppressWarnings(system2(
    "pixi",
    c(
      "run", "-e", "default", "Rscript", "--vanilla",
      file.path(root, "src", "utils", "r", "audit_corrected_source_metadata.R"),
      "--config", config_path,
      "--input-file", input_path,
      "--output", output_path,
      "--dataset", "Fixture",
      "--view", "batch_effect_corrected",
      "--source-root", source_root,
      "--source-manifest", source_manifest,
      "--runtime-identity", runtime_identity,
      "--run-root", run_root
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

report_path <- file.path(run_root, "reports", "Fixture_corrected_source.json")
positive <- run_auditor(source_rds, report_path)
stopifnot(positive$status == 0L, file.exists(report_path))

report <- jsonlite::fromJSON(report_path, simplifyVector = FALSE)
expected_source_md5 <- tolower(unname(as.character(tools::md5sum(source_rds))))
expected_source_sha256 <- sha256(source_rds)
expected_source_size <- as.numeric(file.info(source_rds)$size)

# Status and provenance establish that this is the RDS path, not an H5AD or a
# silently substituted source.  Compare the complete file identity fields.
stopifnot(
  identical(report$status, "SOURCE_METADATA_VALIDATED_RDS"),
  identical(report$source_type, "rds"),
  identical(report$source_path, source_rds),
  identical(report$source_identity$path, source_rds),
  as.numeric(report$source_identity$size) == expected_source_size,
  identical(tolower(report$source_identity$md5), expected_source_md5),
  identical(tolower(report$source_identity$sha256), expected_source_sha256),
  identical(report$provenance$source_manifest$SOURCE_ROOT, source_root),
  identical(report$provenance$runtime_identity$path, runtime_identity)
)

# The membership subset retains S1/S2/S3/S4 (eight cells) and drops S5 (two
# cells).  No repeated sample is split across the row-level mask.
stopifnot(
  as.integer(report$subset_audit$total_cells) == 10L,
  as.integer(report$subset_audit$retained_cells) == 8L,
  as.integer(report$subset_audit$dropped_cells) == 2L,
  as.integer(report$subset_audit$total_samples) == 5L,
  as.integer(report$subset_audit$retained_samples) == 4L,
  as.integer(report$subset_audit$dropped_samples) == 1L,
  as.integer(report$subset_audit$split_sample_count) == 0L,
  identical(as.character(unlist(report$subset_audit$retained_sample_ids, use.names = FALSE)), c("S1", "S2", "S3", "S4")),
  identical(as.character(unlist(report$subset_audit$dropped_sample_ids, use.names = FALSE)), "S5")
)

# Configured columns and the compact contract summary are part of the trust
# boundary: the selected metadata has two batch levels and an estimable design.
stopifnot(
  identical(report$config$sample_column, "sample_id"),
  identical(report$config$label_column, "label"),
  identical(as.character(unlist(report$config$batch_keys, use.names = FALSE)), "batch"),
  identical(report$config$subset_vars$subset_group$op, "in"),
  identical(as.character(unlist(report$config$subset_vars$subset_group$values, use.names = FALSE)), "keep"),
  isTRUE(report$validation_summary$valid),
  isTRUE(report$validation_summary$estimable),
  identical(report$validation_summary$sample_column, "sample_id"),
  identical(report$validation_summary$biological_column, "label"),
  identical(as.character(unlist(report$validation_summary$batch_keys, use.names = FALSE)), "batch"),
  as.integer(report$validation_summary$n_cells) == 8L,
  as.integer(report$validation_summary$n_samples) == 4L,
  as.integer(report$validation_summary$key_level_counts$batch) == 2L,
  as.integer(report$validation_summary$design_rank) == 2L,
  as.integer(report$validation_summary$design_columns) == 2L,
  is.character(report$validation_summary$fingerprint),
  nzchar(report$validation_summary$fingerprint)
)

# A source missing its configured batch column must fail closed before a
# report can be installed.  This is a distinct RDS with the same configured
# basename, so the negative path still reaches metadata extraction.
missing_batch_metadata <- metadata[, setdiff(names(metadata), "batch"), drop = FALSE]
missing_batch_rds <- file.path(tmp, "negative-input", "Fixture.rds")
saveRDS(missing_batch_metadata, missing_batch_rds)
missing_batch_rds <- normalizePath(missing_batch_rds, mustWork = TRUE)
negative_report <- file.path(run_root, "reports", "Fixture_missing_batch.json")
negative <- run_auditor(missing_batch_rds, negative_report)
stopifnot(negative$status != 0L, !file.exists(negative_report))

cat("corrected source metadata regression: OK\n")
}

main()
