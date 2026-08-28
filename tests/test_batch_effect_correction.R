#!/usr/bin/env Rscript

raw_args <- commandArgs(trailingOnly = FALSE)
script_arg <- raw_args[grepl("^--file=", raw_args)][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), ".."))
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_methods_r.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_pipeline.R"))

# Corrected CLR removes the fitted technical effect and preserves the CLR
# invariant exactly after row recentering.
set.seed(11)
n <- 12
feat <- matrix(
  rnorm(n * 4, sd = 0.2),
  nrow = n,
  dimnames = list(paste0("s", seq_len(n)), paste0("ct", seq_len(4)))
)
batch <- rep(c("A", "B"), each = n / 2)
feat[batch == "B", 1:2] <- feat[batch == "B", 1:2] + 4
feat <- feat - rowMeans(feat)
meta <- data.frame(Sample = rownames(feat), tech = factor(batch))
corrected <- correct_clr_batch_lmm(feat, meta, "tech")
stopifnot(max(abs(rowSums(corrected))) < 1e-8)
stopifnot(identical(dimnames(corrected), dimnames(feat)))

bad_order <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
order_error <- tryCatch(
  correct_clr_batch_lmm(feat, bad_order, "tech"),
  error = identity
)
stopifnot(inherits(order_error, "error"))

missing_batch <- meta
missing_batch$tech[1] <- NA
missing_error <- tryCatch(
  correct_clr_batch_lmm(feat, missing_batch, "tech"),
  error = identity
)
stopifnot(inherits(missing_error, "error"))

# The pseudobulk driver forwards the two modes without exposing labels to the
# DESeq2 normalizer.
captured <- list()
get_pb_deseq2 <- function(...) {
  captured <<- list(...)
  matrix(1, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), paste0("s", 1:3)))
}
invisible(prepare_pseudobulks_hpc(
  list(),
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000",
  batch_col = "tech",
  blind = FALSE,
  correct_batch = TRUE
))
stopifnot(identical(captured$batch_col, "tech"))
stopifnot(identical(captured$blind, FALSE))
stopifnot(identical(captured$correct_batch, TRUE))
invisible(prepare_pseudobulks_hpc(
  list(),
  hvg_rank_genes = c("g1", "g2"),
  variants = "hvg2000"
))
stopifnot(identical(captured$blind, TRUE))
stopifnot(identical(captured$correct_batch, FALSE))

single_batch_error <- tryCatch(
  correct_clr_batch_lmm(feat, transform(meta, tech = factor("A")), "tech"),
  error = identity
)
stopifnot(inherits(single_batch_error, "error"))

# The correction model is batch-only by construction.
correction_body <- paste(deparse(body(correct_clr_batch_lmm)), collapse = " ")
stopifnot(grepl("y ~ 1 \\+ \\(1 \\| batch\\)", correction_body))
stopifnot(!grepl("label|Status|disease", correction_body, ignore.case = TRUE))

# Batch pseudobulk result bundles use a pass-qualified stem and only the
# hvg2000 high-resolution result.
process_pseudobulk_fig <- function(feat_mat, labels, ...) {
  list(feat_mat = feat_mat, labels = labels)
}
tmp_results <- tempfile()
dir.create(tmp_results)
pb <- list(
  hvg2000 = list(
    pb = matrix(1, nrow = 3, ncol = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:3))),
    time_secs = 0
  )
)
batch_pb <- run_pseudobulk_hpc(
  list(),
  labels = structure(factor(c("A", "B", "A")), names = paste0("s", 1:3)),
  pb_variants = pb,
  results_dir = tmp_results,
  ds = "DS",
  batch_mode = TRUE,
  result_stem = "DS_batch_effect_uncorrected"
)

stopifnot(identical(names(batch_pb), "Pseudobulk_hvg2000"))
stopifnot(file.exists(file.path(
  tmp_results, "DS_batch_effect_uncorrected_Pseudobulk_hvg2000.rds"
)))

cat("batch-effect CLR and pseudobulk modes OK\n")
# CT pseudobulk must fail closed when every per-CT normalization fails, and
# successful runs must expose contribution counts in the result bundle.
ct_counts <- matrix(
  1,
  nrow = 4,
  ncol = 6,
  dimnames = list(paste0("g", 1:4), paste0("cell", 1:6))
)
ct_meta <- data.frame(
  Sample = rep(paste0("s", 1:3), each = 2),
  ct = rep(c("A", "B"), 3),
  row.names = colnames(ct_counts)
)
ct_seurat <- Seurat::CreateSeuratObject(
  counts = ct_counts,
  meta.data = ct_meta
)
ct_labels <- structure(
  factor(c("A", "B", "A")),
  names = paste0("s", 1:3)
)
get_pb_deseq2 <- function(...) {
  stop("synthetic per-CT failure")
}
all_failed <- tryCatch(
  process_pseudobulk_ct_fig(
    ct_seurat,
    ct_labels,
    ct_col = "ct",
    sample_col = "Sample",
    hvg = 2,
    min_cells = 1
  ),
  error = identity
)
stopifnot(
  inherits(all_failed, "error"),
  grepl("no successful cell-type pseudobulks", conditionMessage(all_failed))
)

create_result_bundle <- function(feat_mat, labels, dist_mat = NULL, extra = list()) {
  c(list(feat_mat = feat_mat, labels = labels), extra)
}
get_pb_deseq2 <- function(...) {
  matrix(
    seq_len(6),
    nrow = 3,
    ncol = 2,
    dimnames = list(paste0("s", 1:3), c("g1", "g2"))
  )
}
ct_success <- process_pseudobulk_ct_fig(
  ct_seurat,
  ct_labels,
  ct_col = "ct",
  sample_col = "Sample",
  hvg = 2,
  min_cells = 1
)
stopifnot(
  identical(ct_success$n_ct_success, 2L),
  identical(ct_success$n_sample_pairs_contributed, 3L),
  identical(ct_success$n_ct_pair_contributions, 6L),
  identical(ct_success$successful_cell_types, c("A", "B"))
)

cat("batch-effect CLR, pseudobulk modes, and CT contribution guard OK\n")
