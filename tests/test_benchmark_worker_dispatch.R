#!/usr/bin/env Rscript
# Focused regression for the counts-free GloScope worker dispatch boundary.

script_arg <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(
  sub("^--file=", "", script_arg), mustWork = TRUE
)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
worker_path <- file.path(
  project_root,
  "src",
  "5_run_benchmark_methods",
  "run_r_sample_embedding_methods",
  "1.1.1_run_benchmark_methods_r.R"
)

with_temporary <- function(code) {
  fixture_root <- tempfile("ecoda-r-worker-dispatch-")
  dir.create(fixture_root, recursive = TRUE)
  on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)
  eval(substitute(code), envir = environment())
}

with_temporary({
input_path <- file.path(fixture_root, "Synthetic.h5ad")
writeBin(charToRaw("counts-backed fixture sentinel"), input_path)
config_path <- file.path(fixture_root, "datasets.json")
writeLines("{}", config_path)
results_dir <- file.path(fixture_root, "results")
log_file <- file.path(fixture_root, "execution.feather")

fake_obs <- data.frame(
  Sample = c("s1", "s1", "s2"),
  label = c("case", "case", "control"),
  row.names = c("cell1", "cell2", "cell3"),
  stringsAsFactors = FALSE
)
fake_embeddings <- list(
  X_pca_benchmark_analysis_hvg1000 = matrix(
    c(1, 2, 3, 4, 5, 6), nrow = 3L, byrow = TRUE,
    dimnames = list(rownames(fake_obs), c("PC1", "PC2"))
  ),
  X_pca_benchmark_analysis_hvg2000 = matrix(
    c(11, 12, 13, 14, 15, 16), nrow = 3L, byrow = TRUE,
    dimnames = list(rownames(fake_obs), c("PC1", "PC2"))
  ),
  X_pca_benchmark_analysis_hvg3000 = matrix(
    c(21, 22, 23, 24, 25, 26), nrow = 3L, byrow = TRUE,
    dimnames = list(rownames(fake_obs), c("PC1", "PC2"))
  )
)
fake_adata <- list(obs = fake_obs, obsm = fake_embeddings)
events <- new.env(parent = emptyenv())
events$counts_free_calls <- 0L
events$seurat_calls <- 0L

environment <- new.env(parent = globalenv())
environment$source <- function(file, ...) invisible(NULL)
environment$commandArgs <- function(trailingOnly = FALSE) {
  if (isTRUE(trailingOnly)) {
    c(
      "--config_path", config_path,
      "--ds_name", "Synthetic",
      "--view", "benchmark_analysis",
      "--method", "gloscope",
      "--input_dir", dirname(input_path),
      "--results_dir", results_dir,
      "--log_file", log_file,
      "--force"
    )
  } else {
    paste0("--file=", worker_path)
  }
}
environment$parse_flags <- function(raw_args) {
  parsed <- list()
  index <- 1L
  while (index <= length(raw_args)) {
    token <- raw_args[[index]]
    if (!startsWith(token, "--")) stop("unexpected worker argument: ", token)
    key <- sub("^--", "", token)
    if (identical(key, "force")) {
      parsed[[key]] <- TRUE
    } else {
      index <- index + 1L
      if (index > length(raw_args) || startsWith(raw_args[[index]], "--")) {
        stop("missing value for worker argument: ", token)
      }
      parsed[[key]] <- raw_args[[index]]
    }
    index <- index + 1L
  }
  parsed
}
environment$read_datasets_json <- function(path, view = NULL) {
  list(Synthetic = list(
    label_col = "label",
    cell_type_low_res = "ct_low",
    cell_type_high_res = "ct_high",
    batch_col = NULL,
    not_suitable_for_auto_annotation = character()
  ))
}
environment$get_h5ad_path <- function(config, ds, view, input_dir) input_path
environment$load_h5ad_counts_free <- function(
  h5ad_path,
  obs_columns,
  embedding_keys,
  ...
) {
  events$counts_free_calls <- events$counts_free_calls + 1L
  events$loader_path <- h5ad_path
  events$loader_obs_columns <- as.character(obs_columns)
  events$loader_embedding_keys <- as.character(embedding_keys)
  events$loader_method <- list(...)$method
  fake_adata
}
environment$load_benchmark_seurat <- function(...) {
  events$seurat_calls <- events$seurat_calls + 1L
  stop("count-backed Seurat branch was reached")
}
environment$py_to_r <- function(value) value
environment$get_hvg_rank_genes <- function(adata) character()
environment$collapse_sample_metadata <- function(obs, sample_col = "Sample") {
  selected <- !duplicated(obs[[sample_col]])
  reduced <- obs[selected, c(sample_col, "label"), drop = FALSE]
  rownames(reduced) <- reduced[[sample_col]]
  reduced
}
environment$run_gloscope_hpc <- function(
  seurat,
  metadata,
  ...,
  embedding_matrices = NULL,
  embedding_sample_ids = NULL
) {
  events$gloscope <- list(
    seurat = seurat,
    metadata = metadata,
    embedding_matrices = embedding_matrices,
    embedding_sample_ids = embedding_sample_ids
  )
  list(dispatch = "gloscope")
}
environment$save_rds_atomic <- function(object, path) {
  saveRDS(object, path)
  invisible(path)
}

Sys.setenv(PROJECT_ROOT = fixture_root)
Sys.unsetenv(c(
  "ECODA_RUN_ID",
  "ECODA_RUN_ROOT",
  "ECODA_ARTIFACT_PRODUCER_RUN_ID",
  "ANALYSIS_VARIANT",
  "ANALYSIS_PASS"
))
sys.source(worker_path, envir = environment, keep.source = FALSE)

stopifnot(identical(events$counts_free_calls, 1L))
stopifnot(identical(events$seurat_calls, 0L))
stopifnot(identical(events$loader_path, input_path))
stopifnot(identical(events$loader_method, "gloscope"))
stopifnot(identical(events$loader_obs_columns, c("Sample", "label")))
stopifnot(identical(
  events$loader_embedding_keys,
  c(
    "X_pca_benchmark_analysis_hvg1000",
    "X_pca_benchmark_analysis_hvg2000",
    "X_pca_benchmark_analysis_hvg3000"
  )
))
stopifnot(is.list(events$gloscope))
stopifnot(is.null(events$gloscope$seurat))
stopifnot(identical(
  names(events$gloscope$embedding_matrices),
  c("hvg1000", "hvg2000", "hvg3000")
))
stopifnot(identical(
  events$gloscope$embedding_matrices[["hvg2000"]],
  fake_embeddings[["X_pca_benchmark_analysis_hvg2000"]]
))
stopifnot(identical(
  events$gloscope$embedding_sample_ids,
  c("s1", "s1", "s2")
))
})

cat("GloScope counts-free worker dispatch: OK\n")
