#!/usr/bin/env Rscript
# Focused regression tests for benchmark_pipeline.R RDS loading boundaries.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
pipeline <- file.path(root, "src", "5_run_benchmark_methods", "benchmark_pipeline.R")

new_pipeline_env <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(pipeline, envir = environment)
  environment
}

assert_error <- function(expr, pattern, label) {
  condition <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (is.null(condition)) stop("expected error: ", label)
  if (!grepl(pattern, conditionMessage(condition), fixed = TRUE)) {
    stop("unexpected error for ", label, ": ", conditionMessage(condition))
  }
  invisible(condition)
}

capture_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    force(expr),
    warning = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = messages)
}

read_guard <- function(counter, message) {
  function(path) {
    counter$count <- counter$count + 1L
    stop(message)
  }
}

write_listed_sidecar <- function(fixture_root, files) {
  writeLines(
    vapply(
      files,
      function(path) paste0(
        unname(tools::md5sum(path)),
        "  ",
        file.path("results", basename(path))
      ),
      character(1)
    ),
    file.path(fixture_root, "checksums.md5")
  )
}

withTemporary <- function(code) {
  directory <- tempfile("ecoda-benchmark-loader-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  eval(substitute(code), envir = environment())
}

withTemporary({
  new_fixture <- function(label) {
    fixture_root <- file.path(directory, label)
    dir.create(file.path(fixture_root, "results"), recursive = TRUE)
    fixture_root
  }

  # A missing checksums.md5 is an explicit legacy skip, not permission to
  # deserialize an otherwise plausible-looking RDS file.
  missing_root <- new_fixture("missing-sidecar")
  missing_results <- file.path(missing_root, "results")
  missing_file <- file.path(missing_results, "Toy_gloscope.rds")
  saveRDS(list(marker = "unverified"), missing_file, compress = FALSE)
  missing_reads <- new.env(parent = emptyenv())
  missing_reads$count <- 0L
  missing_pipeline <- new_pipeline_env()
  missing_pipeline$readRDS <- read_guard(
    missing_reads,
    "readRDS was invoked for an unverified benchmark result"
  )
  missing <- capture_warnings(
    missing_pipeline$load_hpc_benchmark_results(
      list(), "Toy", missing_results, methods = "gloscope"
    )
  )
  if (!any(grepl("legacy_unverified", missing$messages, fixed = TRUE)) ||
      !any(grepl("readRDS was not attempted", missing$messages, fixed = TRUE))) {
    stop("missing checksums.md5 was not reported as legacy_unverified: ",
         paste(missing$messages, collapse = " | "))
  }
  if (!identical(missing_reads$count, 0L)) {
    stop("readRDS was invoked for a missing checksums.md5 sidecar")
  }
  if (length(missing$value$bmark$Toy) != 0L) {
    stop("legacy-unverified benchmark result was loaded")
  }

  # A present sidecar that omits the result is also unverified and must not be
  # treated as a legacy-compatible read.
  unlisted_root <- new_fixture("unlisted-sidecar")
  unlisted_results <- file.path(unlisted_root, "results")
  unlisted_file <- file.path(unlisted_results, "Toy_gloscope.rds")
  unrelated_file <- file.path(unlisted_results, "Toy_other.rds")
  saveRDS(list(marker = "unlisted"), unlisted_file, compress = FALSE)
  saveRDS(list(marker = "unrelated"), unrelated_file, compress = FALSE)
  write_listed_sidecar(unlisted_root, unrelated_file)
  unlisted_reads <- new.env(parent = emptyenv())
  unlisted_reads$count <- 0L
  unlisted_pipeline <- new_pipeline_env()
  unlisted_pipeline$readRDS <- read_guard(
    unlisted_reads,
    "readRDS was invoked for an unlisted benchmark result"
  )
  unlisted <- capture_warnings(
    unlisted_pipeline$load_hpc_benchmark_results(
      list(), "Toy", unlisted_results, methods = "gloscope"
    )
  )
  if (!any(grepl("legacy_unverified", unlisted$messages, fixed = TRUE)) ||
      !any(grepl("not listed in checksums.md5", unlisted$messages, fixed = TRUE))) {
    stop("unlisted benchmark result was not reported as legacy_unverified: ",
         paste(unlisted$messages, collapse = " | "))
  }
  if (!identical(unlisted_reads$count, 0L)) {
    stop("readRDS was invoked for an unlisted checksums.md5 entry")
  }
  if (length(unlisted$value$bmark$Toy) != 0L) {
    stop("unlisted benchmark result was loaded")
  }

  # A listed sidecar is checked immediately before deserialization and then
  # the result bundle is made available to the caller.
  valid_root <- new_fixture("valid-sidecar")
  valid_results <- file.path(valid_root, "results")
  valid_file <- file.path(valid_results, "Toy_gloscope.rds")
  saveRDS(list(valid_marker = "loaded"), valid_file, compress = FALSE)
  write_listed_sidecar(valid_root, valid_file)
  valid_reads <- new.env(parent = emptyenv())
  valid_reads$count <- 0L
  valid_pipeline <- new_pipeline_env()
  valid_pipeline$readRDS <- function(path) {
    valid_reads$count <- valid_reads$count + 1L
    base::readRDS(path)
  }
  valid <- capture_warnings(
    valid_pipeline$load_hpc_benchmark_results(
      list(), "Toy", valid_results, methods = "gloscope"
    )
  )
  if (!identical(valid_reads$count, 1L)) {
    stop("listed benchmark result was not read exactly once")
  }
  if (!identical(valid$value$bmark$Toy$valid_marker, "loaded")) {
    stop("listed benchmark result was not returned by the loader")
  }
  if (any(grepl("legacy_unverified", valid$messages, fixed = TRUE))) {
    stop("valid listed benchmark result was reported as legacy_unverified")
  }

  # Malformed checksums.md5 content is rejected while the RDS boundary remains
  # untouched, even when the payload itself would otherwise be readable.
  malformed_root <- new_fixture("malformed-sidecar")
  malformed_results <- file.path(malformed_root, "results")
  malformed_file <- file.path(malformed_results, "Toy_gloscope.rds")
  saveRDS(list(marker = "malformed-sidecar"), malformed_file, compress = FALSE)
  writeLines("not a checksums.md5 record", file.path(malformed_root, "checksums.md5"))
  malformed_reads <- new.env(parent = emptyenv())
  malformed_reads$count <- 0L
  malformed_pipeline <- new_pipeline_env()
  malformed_pipeline$readRDS <- read_guard(
    malformed_reads,
    "readRDS was invoked for a malformed checksums.md5 sidecar"
  )
  assert_error(
    malformed_pipeline$load_hpc_benchmark_results(
      list(), "Toy", malformed_results, methods = "gloscope"
    ),
    "Malformed checksums.md5 sidecar",
    "malformed checksums.md5"
  )
  if (!identical(malformed_reads$count, 0L)) {
    stop("readRDS was invoked before malformed sidecar rejection")
  }
})

cat("benchmark pipeline loader: OK\n")
