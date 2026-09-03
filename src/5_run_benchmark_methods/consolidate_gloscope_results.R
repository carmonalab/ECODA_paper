#!/usr/bin/env Rscript
# Consolidate run-owned GloScope parameter shards into the canonical method RDS.
# The shard workers write one validated bundle per dataset/parameter; this
# serialized tail runs only after the matrix watchdog has gated every shard.

project_root <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT is required")
}
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

args <- parse_flags(commandArgs(trailingOnly = TRUE))
for (required in c("manifest", "results_dir")) {
  if (is.null(args[[required]]) || identical(args[[required]], TRUE) || !nzchar(args[[required]])) {
    stop("Missing required --", required, " argument")
  }
}

manifest <- normalizePath(args$manifest, mustWork = TRUE)
results_dir <- normalizePath(args$results_dir, mustWork = TRUE)
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
  if (!all(vapply(shard_paths, artifact_checksum_ok, logical(1)))) {
    bad <- shard_paths[!vapply(shard_paths, artifact_checksum_ok, logical(1))]
    stop("GloScope shard checksum validation failed: ", paste(bad, collapse = ", "))
  }
  bundles <- lapply(shard_paths, readRDS)
  names(bundles) <- combo_names
  method_path <- file.path(results_dir, paste0(dataset, "_gloscope.rds"))
  save_rds_atomic(bundles, method_path)
  message("Consolidated ", length(bundles), " GloScope shards for ", dataset,
          " -> ", method_path)
}
