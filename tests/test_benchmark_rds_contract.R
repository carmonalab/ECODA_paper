#!/usr/bin/env Rscript
# Focused regression tests for method-specific benchmark RDS contracts.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
validator <- file.path(root, "src", "5_run_benchmark_methods", "validate_benchmark_rds_contract.R")
RDS_TEST_ENV <- c(
  paste0("RETICULATE_PYTHON=", file.path(root, ".pixi", "envs", "default", "bin", "python"))
)

write_checked <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(value, path)
  digest <- unname(tools::md5sum(path))
  writeLines(c(
    paste0("MD5=", digest),
    paste0("SIZE=", file.info(path)$size),
    paste0("PATH=", path)
  ), paste0(path, ".md5"))
}

write_checked_text <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  digest <- unname(tools::md5sum(path))
  writeLines(c(
    paste0("MD5=", digest),
    paste0("SIZE=", file.info(path)$size),
    paste0("PATH=", path)
  ), paste0(path, ".md5"))
}

run_validator <- function(arguments) {
  system2(
    "pixi", c("run", "Rscript", "--vanilla", validator, arguments),
    stdout = FALSE, stderr = FALSE, env = RDS_TEST_ENV
  )
}

run_validator_capture <- function(arguments) {
  output <- system2(
    "pixi", c("run", "Rscript", "--vanilla", validator, arguments),
    stdout = TRUE, stderr = TRUE, env = RDS_TEST_ENV
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = status, output = paste(output, collapse = "\n"))
}

expect_ok <- function(arguments, label) {
  status <- run_validator(arguments)
  if (!identical(status, 0L)) stop("expected validator success: ", label)
}

expect_fail <- function(arguments, label) {
  status <- run_validator(arguments)
  if (identical(status, 0L)) stop("expected validator failure: ", label)
}

combo <- function(ids = c("s1", "s2")) {
  feature <- matrix(
    c(1, 0, 0, 1), nrow = 2L,
    dimnames = list(ids, c("f1", "f2"))
  )
  distance <- as.matrix(dist(feature))
  dimnames(distance) <- list(ids, ids)
  labels <- factor(c("A", "B"))
  names(labels) <- ids
  list(
    scores = list(sil_score = 0.5),
    feat_mat = feature,
    dist_mat = distance,
    labels = labels
  )
}

withTemporary <- function(code) {
  directory <- tempfile("ecoda-rds-contract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  eval(substitute(code), envir = environment())
}

withTemporary({
  batch_args <- function(path, method) c(
    "--artifact", path,
    "--method", method,
    "--dataset", "Synthetic",
    "--view", "batch_effect_uncorrected",
    "--batch-pass", "uncorrected"
  )

  gloscope <- file.path(directory, "gloscope.rds")
  write_checked(gloscope, list(GloScope_hvg2000_pcadims30 = combo()))
  expect_ok(batch_args(gloscope, "gloscope"), "batch GloScope")

  pseudobulk <- file.path(directory, "pseudobulk_hvg2000.rds")
  pb <- matrix(c(1, 2, 3, 4), nrow = 2L, dimnames = list(c("s1", "s2"), c("g1", "g2")))
  write_checked(pseudobulk, pb)
  expect_ok(batch_args(pseudobulk, "pseudobulk"), "batch pseudobulk")
  artifact_list <- file.path(directory, "rds-artifacts.tsv")
  write_checked_text(artifact_list, c(
    paste(gloscope, "gloscope", "Synthetic", "benchmark_analysis", "0", sep = "\t"),
    paste(pseudobulk, "prepare_pseudobulk", "Synthetic", "benchmark_analysis", "0", sep = "\t")
  ))
  expect_ok(c("--artifact-list", artifact_list), "grouped RDS artifact list")

  composition <- file.path(directory, "composition.rds")
  composition_keys <- c(
    "ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2",
    "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR"
  )
  write_checked(composition, setNames(lapply(composition_keys, function(x) combo()), composition_keys))
  expect_ok(batch_args(composition, "composition"), "batch composition")

  missing_key <- file.path(directory, "gloscope-missing.rds")
  write_checked(missing_key, list(GloScope_hvg2000_pcadims30 = combo(), extra = combo()))
  expect_fail(batch_args(missing_key, "gloscope"), "extra GloScope key")
  write_checked(missing_key, list())
  expect_fail(batch_args(missing_key, "gloscope"), "missing GloScope key")

  trans <- file.path(directory, "trans.rds")
  write_checked(trans, data.frame(
    trans_method = "counts",
    ANOSIM_score = 0.1,
    Modularity_score = 0.2,
    Adjusted_Rand_Index = 0.3,
    stringsAsFactors = FALSE
  ))
  expect_ok(c("--artifact", trans, "--method", "trans"), "transformation")
  trans_bad <- file.path(directory, "trans-bad.rds")
  write_checked(trans_bad, data.frame(
    trans_method = "counts",
    ANOSIM_score = "not-numeric",
    Modularity_score = 0.2,
    Adjusted_Rand_Index = 0.3,
    stringsAsFactors = FALSE
  ))
  expect_fail(c("--artifact", trans_bad, "--method", "trans"), "nonnumeric transformation score")

  zeroimp <- file.path(directory, "zeroimp.rds")
  write_checked(zeroimp, list(counts_all_1 = list(score = 0.25)))
  expect_ok(c("--artifact", zeroimp, "--method", "zeroimp"), "zero-imputation")
  zeroimp_bad <- file.path(directory, "zeroimp-bad.rds")
  write_checked(zeroimp_bad, list(counts_all_1 = list(score = Inf)))
  expect_fail(c("--artifact", zeroimp_bad, "--method", "zeroimp"), "nonfinite zero-imputation score")

  nonfinite <- file.path(directory, "nonfinite.rds")
  bad_combo <- combo()
  bad_combo$feat_mat[[1L, 1L]] <- Inf
  write_checked(nonfinite, bad_combo)
  expect_fail(c("--artifact", nonfinite, "--method", "mofa"), "nonfinite feature matrix")

  mismatch <- file.path(directory, "mismatch.rds")
  mismatched <- combo()
  names(mismatched$labels) <- c("s2", "s1")
  write_checked(mismatch, mismatched)
  expect_fail(c("--artifact", mismatch, "--method", "mofa"), "label identifier mismatch")

  scratch <- file.path(directory, "scratch")
  h5ad <- file.path(scratch, "Synthetic", "output", "synthetic.h5ad")
  dir.create(dirname(h5ad), recursive = TRUE, showWarnings = FALSE)
  python <- paste(
    "import anndata as ad, numpy as np, pandas as pd, sys;",
    "a=ad.AnnData(X=np.ones((2,1),dtype='float32'),",
    "obs=pd.DataFrame({'Sample':['s1','s2']},index=['c1','c2']));",
    "a.write_h5ad(sys.argv[1])"
  )
  status <- system2("pixi", c("run", "python", "-c", shQuote(python), shQuote(h5ad)))
  if (!identical(status, 0L)) stop("could not create tiny source h5ad")
  writeLines(c(
    paste0("MD5=", unname(tools::md5sum(h5ad))),
    paste0("SIZE=", file.info(h5ad)$size),
    paste0("PATH=", h5ad)
  ), paste0(h5ad, ".md5"))
  config <- file.path(directory, "datasets.json")
  writeLines('{"Synthetic":{"views":{"benchmark_analysis":{"output_file_name":"synthetic.h5ad"}}}}', config)
  identity_selection <- file.path(directory, "identity-selection.tsv")
  writeLines("Synthetic\tbenchmark_analysis\tscitd", identity_selection)
  identity <- file.path(directory, "source-identity.json")
  status <- system2("pixi", c(
    "run", "python", file.path(root, "src", "utils", "py", "h5ad_source_identity.py"),
    "--output", identity, "--selection", identity_selection,
    "--input-root", scratch, "--config", config
  ))
  if (!identical(status, 0L)) stop("could not create source identity")
  writeLines(c(
    paste0("MD5=", unname(tools::md5sum(identity))),
    paste0("SIZE=", file.info(identity)$size),
    paste0("PATH=", identity)
  ), paste0(identity, ".md5"))

  scitd_subset <- list(
    scITD_hvg2000_factors5 = list(
      scores = list(sil_score = 0.5),
      feat_mat = matrix(1, nrow = 1L, ncol = 1L,
                        dimnames = list("s1", "f1")),
      dist_mat = matrix(0, nrow = 1L, ncol = 1L,
                        dimnames = list("s1", "s1")),
      labels = structure(factor("A"), names = "s1")
    )
  )
  scitd_artifact <- file.path(directory, "scitd-subset.rds")
  write_checked(scitd_artifact, scitd_subset)
  scitd_arguments <- c(
    "--artifact", scitd_artifact,
    "--method", "scitd",
    "--dataset", "Synthetic",
    "--view", "benchmark_analysis",
    "--input-root", scratch,
    "--config", config,
    "--source-identity", identity,
    "--source-identity-verified"
  )
  scitd_result <- run_validator_capture(scitd_arguments)
  if (!identical(scitd_result$status, 0L)) stop(scitd_result$output)
  if (!grepl("dropped sample IDs: s2", scitd_result$output, fixed = TRUE)) {
    stop("scITD dropped sample IDs were not reported: ", scitd_result$output)
  }
  mofa_arguments <- scitd_arguments
  mofa_arguments[match("--method", mofa_arguments) + 1L] <- "mofa"
  expect_fail(mofa_arguments, "non-scITD sample drop")

  scitd_root <- file.path(directory, "scitd-root")
  write_checked(file.path(scitd_root, "results", "Synthetic_scitd.rds"), scitd_subset)
  scitd_selection <- file.path(directory, "scitd-selection.tsv")
  writeLines("Synthetic\tbenchmark_analysis\tscitd", scitd_selection)
  writeLines(c(
    paste0("MD5=", unname(tools::md5sum(scitd_selection))),
    paste0("SIZE=", file.info(scitd_selection)$size),
    paste0("PATH=", scitd_selection)
  ), paste0(scitd_selection, ".md5"))
  scitd_root_result <- run_validator_capture(c(
    "--root", scitd_root,
    "--selection", scitd_selection,
    "--labels", "scitd",
    "--input-root", scratch,
    "--config", config
  ))
  if (!identical(scitd_root_result$status, 0L)) stop(scitd_root_result$output)
  if (!grepl("dropped sample IDs: s2", scitd_root_result$output, fixed = TRUE)) {
    stop("root scITD dropped sample IDs were not reported: ", scitd_root_result$output)
  }

  ordinary_root <- file.path(directory, "ordinary")
  ordinary_selection <- file.path(directory, "ordinary-selection.tsv")
  writeLines("Synthetic\tbenchmark_analysis\tgloscope", ordinary_selection)
  writeLines(c(
    paste0("MD5=", unname(tools::md5sum(ordinary_selection))),
    paste0("SIZE=", file.info(ordinary_selection)$size),
    paste0("PATH=", ordinary_selection)
  ), paste0(ordinary_selection, ".md5"))
  reordered <- file.path(ordinary_root, "results", "Synthetic_gloscope.rds")
  write_checked(reordered, list(combo = combo(c("s2", "s1"))))
  expect_fail(c(
    "--root", ordinary_root, "--selection", ordinary_selection, "--labels", "gloscope",
    "--input-root", scratch, "--config", config
  ), "reordered samples")

  writeLines("Synthetic\tbenchmark_analysis\tgloscope-modified", ordinary_selection)
  expect_fail(c(
    "--root", ordinary_root, "--selection", ordinary_selection, "--labels", "gloscope",
    "--input-root", scratch, "--config", config
  ), "selection checksum mismatch")

  writeBin(charToRaw("changed"), h5ad)
  expect_fail(c(
    "--root", ordinary_root, "--selection", ordinary_selection, "--labels", "gloscope",
    "--input-root", scratch, "--config", config
  ), "source h5ad checksum mismatch")
})

cat("benchmark RDS contract: OK\n")
