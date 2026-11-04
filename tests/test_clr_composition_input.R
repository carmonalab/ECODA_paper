#!/usr/bin/env Rscript
# Exercise the actual CLR producer's numeric data.frame at the correction API.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
Sys.setenv(
  PROJECT_ROOT = root,
  RETICULATE_PYTHON = file.path(root, ".pixi/envs/default/bin/python")
)
suppressPackageStartupMessages(library(dplyr))
source(file.path(root, "src/utils/math_utils.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_methods_r.R"))
set.seed(19)
counts <- matrix(
  rpois(48L, 20) + 1,
  nrow = 12L,
  dimnames = list(paste0("s", 1:12), paste0("ct", 1:4))
)
counts[7:12, 1:2] <- counts[7:12, 1:2] * 5
features <- clr(as.data.frame(counts))
metadata <- data.frame(Sample = rownames(features), tech = rep(c("A", "B"), each = 6L))
validation <- ecoda_batch_validate_metadata(metadata, batch_keys = "tech", sample_col = "Sample")
actual <- correct_clr_batch_limma(
  features, metadata, batch_keys = "tech", metadata_validation = validation
)
design <- ecoda_batch_fixed_effect_design(
  metadata, batch_keys = "tech", validation = validation, sample_col = "Sample"
)
expected <- t(limma::removeBatchEffect(
  t(as.matrix(features)),
  covariates = design$design[, design$technical_columns, drop = FALSE],
  design = matrix(1, nrow = nrow(features), ncol = 1L)
))
expected <- sweep(expected, 1L, rowMeans(expected), "-")
stopifnot(
  isTRUE(all.equal(actual, expected, tolerance = 1e-8)),
  identical(dimnames(actual), dimnames(features)),
  all(rowSums(actual) == 0),
  max(abs(actual - as.matrix(features))) > 0.1
)
# Converting the container must not turn nonnumeric features into valid input.
bad <- features
bad[[1L]] <- as.character(bad[[1L]])
rejected <- tryCatch(
  correct_clr_batch_limma(bad, metadata, batch_keys = "tech", metadata_validation = validation),
  error = identity
)
stopifnot(inherits(rejected, "error"))
cat("PASS actual CLR data.frame correction matches limma, preserves sample identity and rejects nonnumeric features\n")
