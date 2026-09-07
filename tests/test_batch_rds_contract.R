#!/usr/bin/env Rscript

raw_args <- commandArgs(trailingOnly = FALSE)
script_arg <- raw_args[grep("^--file=", raw_args)][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), ".."))
validator <- file.path(
  root,
  "src/5_run_benchmark_methods/validate_benchmark_rds_contract.R"
)

combo <- function() {
  feat <- matrix(
    c(1, 2),
    nrow = 2,
    dimnames = list(c("s1", "s2"), "feature")
  )
  labels <- structure(factor(c("A", "B")), names = rownames(feat))
  list(
    scores = 1,
    feat_mat = feat,
    dist_mat = dist(feat),
    labels = labels
  )
}

write_bundle <- function(keys, path) {
  bundle <- setNames(lapply(keys, function(ignored) combo()), keys)
  saveRDS(bundle, path)
  writeLines(
    c(
      paste0("MD5=", unname(tools::md5sum(path))),
      paste0("SIZE=", file.info(path)$size),
      paste0("PATH=", path)
    ),
    paste0(path, ".md5")
  )
  path
}

run_validator <- function(path) {
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "--vanilla", validator,
      "--artifact", path,
      "--method", "composition",
      "--dataset", "Synthetic",
      "--view", "batch_effect_uncorrected",
      "--batch-pass", "uncorrected"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  list(status = attr(output, "status") %||% 0, output = output)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
tmp <- tempfile("batch-rds-contract-")
base <- c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2")

base_result <- run_validator(write_bundle(base, tmp))
stopifnot(base_result$status == 0)

legacy_result <- run_validator(write_bundle(
  c(base, "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR"),
  paste0(tmp, ".legacy")
))
stopifnot(legacy_result$status == 0)

unknown_result <- run_validator(write_bundle(
  c(base, "unexpected_legacy_combo"),
  paste0(tmp, ".unknown")
))
stopifnot(unknown_result$status != 0)

missing_result <- run_validator(write_bundle(
  c("ECODA_authors_HR", "ECODA_authors_HR_NULL"),
  paste0(tmp, ".missing")
))
stopifnot(missing_result$status != 0)

write_metadata <- function(path) {
  labels <- structure(factor(c("A", "B")), names = c("s1", "s2"))
  value <- list(
    labels = labels,
    n_cells = 4,
    n_samples = 2,
    cells_per_sample = structure(c(2, 2), names = c("s1", "s2"))
  )
  saveRDS(value, path)
  writeLines(
    c(
      paste0("MD5=", unname(tools::md5sum(path))),
      paste0("SIZE=", file.info(path)$size),
      paste0("PATH=", path)
    ),
    paste0(path, ".md5")
  )
}

artifact_root <- tempfile("batch-rds-root-")
dir.create(file.path(artifact_root, "results"), recursive = TRUE)
selection <- file.path(artifact_root, "selection.tsv")
writeLines(
  "Synthetic\tbatch_effect_uncorrected\tbatch_effect_uncorrected",
  selection
)
writeLines(
  c(
    paste0("MD5=", unname(tools::md5sum(selection))),
    paste0("SIZE=", file.info(selection)$size),
    paste0("PATH=", selection)
  ),
  paste0(selection, ".md5")
)
root_bundle <- file.path(
  artifact_root,
  "results",
  "Synthetic_batch_effect_uncorrected_composition.rds"
)
write_bundle(c(base, "ECODA_HiTME_HR_layer2"), root_bundle)
write_metadata(file.path(
  artifact_root,
  "results",
  "Synthetic_batch_effect_uncorrected_metadata.rds"
))
root_output <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "--vanilla", validator,
    "--root", artifact_root,
    "--selection", selection,
    "--labels", "composition",
    "--batch-pass", "uncorrected"
  ),
  stdout = TRUE,
  stderr = TRUE
)
stopifnot((attr(root_output, "status") %||% 0) == 0)

cat("batch RDS contract: OK\n")
