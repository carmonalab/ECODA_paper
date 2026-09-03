#!/usr/bin/env Rscript
set.seed(1)
root <- normalizePath(tempfile("ecoda-gloscope-consolidation-"), mustWork = FALSE)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
physical_results_dir <- file.path(root, "results")
dir.create(physical_results_dir)
results_dir <- file.path(root, "results_alias")
stopifnot(file.symlink(physical_results_dir, results_dir))
manifest <- normalizePath(file.path(root, "matrix.tsv"), mustWork = FALSE)

repo_root <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(PROJECT_ROOT = repo_root)
source(file.path(repo_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

combos <- c(
  "hvg2000_pcadims10",
  "hvg2000_pcadims30",
  "hvg2000_pcadims50",
  "hvg1000_pcadims30",
  "hvg3000_pcadims30"
)
writeLines(
  paste("Adams", "benchmark_analysis", "gloscope", combos, sep = "\t"),
  manifest
)
ids <- c("s1", "s2")
mat <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(ids, ids))
for (combo in combos) {
  bundle <- list(
    scores = c(silhouette = 0.1),
    feat_mat = mat,
    dist_mat = as.dist(mat),
    labels = structure(c("A", "B"), names = ids)
  )
  save_rds_atomic(
    bundle,
    file.path(results_dir, paste0("Adams_GloScope_", combo, ".rds"))
  )
}

script <- file.path(repo_root, "src/5_run_benchmark_methods/consolidate_gloscope_results.R")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c("--vanilla", script, "--manifest", manifest, "--results_dir", results_dir)
)
stopifnot(identical(status, 0L))
method_path <- file.path(results_dir, "Adams_gloscope.rds")
stopifnot(artifact_checksum_ok(method_path))
method <- readRDS(method_path)
stopifnot(identical(names(method), paste0("GloScope_", combos)))
cat("GloScope shard consolidation: OK\n")
