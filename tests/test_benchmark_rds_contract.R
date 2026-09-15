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

expect_fail_capture <- function(arguments, label, pattern) {
  result <- run_validator_capture(arguments)
  if (identical(result$status, 0L)) stop("expected validator failure: ", label)
  if (!grepl(pattern, result$output, fixed = TRUE)) {
    stop("validator failure did not report ", label, ": ", result$output)
  }
  invisible(result)
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

pseudobulk_result <- function(ids = c("s1", "s2")) {
  values <- matrix(
    seq_len(length(ids) * 2L),
    nrow = length(ids),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(ids, c("g1", "g2"))
  )
  list(pb = values, time_secs = 0, mem_GB = 0)
}

withTemporary <- function(code) {
  directory <- tempfile("ecoda-rds-contract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  eval(substitute(code), envir = environment())
}
sys.source(file.path(root, "src", "utils", "batch_contract.R"), envir = .GlobalEnv)
sys.source(
  file.path(root, "src", "5_run_benchmark_methods", "benchmark_hpc_utils.R"),
  envir = .GlobalEnv
)

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

  missing_sidecar <- file.path(directory, "gloscope-missing-sidecar.rds")
  writeBin(charToRaw("not an RDS stream"), missing_sidecar)
  expect_fail_capture(
    batch_args(missing_sidecar, "gloscope"),
    "missing result sidecar",
    "Missing or invalid result checksum"
  )

  malformed_sidecar <- file.path(directory, "gloscope-malformed-sidecar.rds")
  writeBin(charToRaw("not an RDS stream"), malformed_sidecar)
  writeLines("not a checksums.md5 record", paste0(malformed_sidecar, ".md5"))
  expect_fail_capture(
    batch_args(malformed_sidecar, "gloscope"),
    "malformed result sidecar",
    "Missing or invalid result checksum"
  )

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

  final_batch_args <- function(path, method) c(
    "--artifact", path,
    "--method", method,
    "--dataset", "Synthetic",
    "--view", "batch_effect_uncorrected",
    "--batch-pass", "uncorrected",
    "--analysis-variant", "final"
  )
  final_root <- file.path(directory, "uncorrected_final")
  final_results <- file.path(final_root, "results")
  final_cache <- file.path(final_root, "pseudobulks")
  dir.create(final_results, recursive = TRUE)
  dir.create(final_cache, recursive = TRUE)

  final_composition <- file.path(
    final_results,
    "Synthetic_batch_effect_uncorrected_final_composition.rds"
  )
  final_composition_keys <- c(
    "ECODA_authors_HR",
    "ECODA_seuratres_2",
    "ECODA_authors_HR_NULL"
  )
  write_checked(
    final_composition,
    setNames(lapply(final_composition_keys, function(x) combo()), final_composition_keys)
  )
  expect_ok(
    final_batch_args(final_composition, "composition"),
    "final composition shared bundle with explicit null key"
  )

  final_metadata <- file.path(
    final_results,
    "Synthetic_batch_effect_uncorrected_final_metadata.rds"
  )
  write_checked(
    final_metadata,
    list(
      labels = structure(factor(c("A", "B")), names = c("s1", "s2")),
      n_cells = 200,
      n_samples = 2,
      cells_per_sample = structure(c(100, 100), names = c("s1", "s2"))
    )
  )

  final_pseudobulk <- file.path(
    final_results,
    "Synthetic_batch_effect_uncorrected_final_pseudobulk.rds"
  )
  write_checked(
    final_pseudobulk,
    list(Pseudobulk_hvg2000 = pseudobulk_result())
  )
  final_pseudobulk_cache <- file.path(
    final_cache,
    "Synthetic_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds"
  )
  write_checked(final_pseudobulk_cache, pb)
  expect_ok(
    final_batch_args(final_pseudobulk, "pseudobulk"),
    "final pseudobulk result bundle"
  )

  final_selection <- file.path(directory, "final-selection.tsv")
  write_checked_text(
    final_selection,
    c(
      "Synthetic\tbatch_effect_uncorrected\tbatch_effect_uncorrected"
    )
  )
  expect_ok(
    c(
      "--root", final_root,
      "--selection", final_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "uncorrected",
      "--analysis-variant", "final"
    ),
    "final pseudobulk result path under results"
  )
  unlink(final_pseudobulk)
  expect_fail(
    c(
      "--root", final_root,
      "--selection", final_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "uncorrected",
      "--analysis-variant", "final"
    ),
    "final pseudobulk cache is not a result bundle"
  )
  write_checked(final_pseudobulk, list(Pseudobulk_hvg2000 = pseudobulk_result()))
  expect_ok(
    c(
      "--root", final_root,
      "--selection", final_selection,
      "--labels", "composition",
      "--batch-pass", "uncorrected",
      "--analysis-variant", "final"
    ),
    "final composition result and metadata paths"
  )
  corrected_sample_metadata <- data.frame(
    Sample = c("s1", "s2", "s3", "s4"),
    batch = factor(c("A", "B", "A", "B")),
    stringsAsFactors = FALSE
  )
  corrected_validation <- ecoda_batch_validate_metadata(
    corrected_sample_metadata,
    batch_keys = "batch",
    sample_col = "Sample"
  )
  corrected_identity <- function(method_id) {
    model_id <- if (identical(method_id, "Pseudobulk")) {
      "pseudobulk_limma_fixed_effects_v1"
    } else if (identical(method_id, "GloScope")) {
      "embedding_consumer_harmony_v1"
    } else {
      "limma_fixed_effects_v1"
    }
    ecoda_hpc_augment_batch_contract(
      identity = ecoda_hpc_batch_contract_identity(
        batch_keys = "batch",
        sample_col = "Sample",
        method_id = method_id,
        model_id = model_id
      ),
      validation = corrected_validation,
      method_id = method_id,
      batch_keys = "batch"
    )
  }
  corrected_combo <- function(method_id) {
    ids <- corrected_sample_metadata$Sample
    feature <- matrix(
      c(1, 0, 0, 1, 1, 2, 2, 3),
      nrow = length(ids),
      byrow = TRUE,
      dimnames = list(ids, c("f1", "f2"))
    )
    labels <- structure(factor(c("A", "B", "A", "B")), names = ids)
    result <- list(
      scores = list(sil_score = 0.5),
      feat_mat = feature,
      dist_mat = as.matrix(dist(feature)),
      labels = labels
    )
    result$batch_contract <- corrected_identity(method_id)
    result
  }
  corrected_root <- file.path(directory, "batch_effect", "corrected_final")
  corrected_results <- file.path(corrected_root, "results")
  corrected_cache <- file.path(corrected_root, "pseudobulks")
  dir.create(corrected_results, recursive = TRUE)
  dir.create(corrected_cache, recursive = TRUE)
  corrected_config <- file.path(directory, "corrected-final-datasets.json")
  writeLines(
    paste0(
      '{"Synthetic":{"columns":{"batch":["batch"]},"views":',
      '{"benchmark_analysis":{"output_file_name":"synthetic.h5ad"},',
      '"batch_effect_corrected":{"output_file_name":"synthetic-corrected.h5ad"}}}}'
    ),
    corrected_config
  )
  corrected_composition <- file.path(
    corrected_results,
    "Synthetic_batch_effect_corrected_final_composition.rds"
  )
  corrected_composition_bundle <- list(
    batch_contract = corrected_identity("ECODA_authors_HR"),
    ECODA_authors_HR = corrected_combo("ECODA_authors_HR"),
    ECODA_seuratres_2 = corrected_combo("ECODA_seuratres_2"),
    ECODA_authors_HR_NULL = corrected_combo("ECODA_authors_HR_NULL")
  )
  write_checked(corrected_composition, corrected_composition_bundle)
  corrected_batch_args <- function(path, method) c(
    "--artifact", path,
    "--method", method,
    "--dataset", "Synthetic",
    "--view", "batch_effect_corrected",
    "--batch-pass", "corrected",
    "--analysis-variant", "corrected_final",
    "--config", corrected_config
  )
  expect_ok(
    corrected_batch_args(corrected_composition, "composition"),
    "corrected-final composition shared bundle with explicit null key"
  )
  # Corrected-final GloScope uses its semantic result key and final-qualified
  # path, just like the worker and synchronization contracts.
  corrected_gloscope_identity <- ecoda_hpc_augment_batch_contract(
    identity = ecoda_hpc_batch_contract_identity(
      batch_keys = "batch",
      sample_col = "Sample",
      method_id = "GloScope",
      model_id = "embedding_consumer_harmony_v1"
    ),
    validation = corrected_validation,
    method_id = "GloScope",
    batch_keys = "batch"
  )
  corrected_gloscope_value <- corrected_combo("GloScope")
  corrected_gloscope_value$batch_contract <- corrected_gloscope_identity
  corrected_gloscope <- file.path(
    corrected_results,
    "Synthetic_batch_effect_corrected_final_gloscope.rds"
  )
  corrected_gloscope_bundle <- list(
    batch_contract = corrected_gloscope_identity,
    GloScope_hvg2000_pcadims30 = corrected_gloscope_value
  )
  write_checked(corrected_gloscope, corrected_gloscope_bundle)
  expect_ok(
    corrected_batch_args(corrected_gloscope, "gloscope"),
    "corrected-final GloScope result path and bundle key"
  )
  corrected_matrix_root <- file.path(
    directory, "batch_effect", "corrected_final", "recovery_35row"
  )
  corrected_matrix_results <- file.path(corrected_matrix_root, "results")
  corrected_matrix_gloscope <- file.path(
    corrected_matrix_results,
    "Synthetic_batch_effect_corrected_final_gloscope.rds"
  )
  write_checked(corrected_matrix_gloscope, corrected_gloscope_bundle)
  expect_ok(
    corrected_batch_args(corrected_matrix_gloscope, "gloscope"),
    "corrected-final matrix recovery root and direct GloScope path"
  )
  # Recovery roots must retain the same strict corrected identity checks as
  # the canonical corrected-final root; root binding is not an identity
  # validation bypass.
  bad_recovery_bundle <- corrected_gloscope_bundle
  bad_recovery_bundle$batch_contract$required_source_obs_columns <-
    c("batch", "Sample")
  write_checked(corrected_matrix_gloscope, bad_recovery_bundle)
  expect_fail_capture(
    corrected_batch_args(corrected_matrix_gloscope, "gloscope"),
    "corrected-final recovery root identity",
    "source/config identity"
  )
  write_checked(corrected_matrix_gloscope, corrected_gloscope_bundle)
  stopifnot(
    identical(
      corrected_composition_bundle$batch_contract$model_id,
      "limma_fixed_effects_v1"
    ),
    identical(
      corrected_composition_bundle$batch_contract$correction_mode,
      "limma_fixed_effects"
    ),
    identical(
      corrected_composition_bundle$batch_contract$fixed_effect_aliases,
      c(batch = "batch_key_1")
    ),
    identical(
      corrected_composition_bundle$batch_contract$correction_design_formula,
      "~1 + batch_key_1"
    ),
    identical(
      corrected_composition_bundle$batch_contract$effective_batch_keys,
      "batch"
    ),
    identical(
      corrected_composition_bundle$batch_contract$non_estimable_batch_keys,
      character()
    ),
    identical(
      corrected_composition_bundle$batch_contract$correction_state,
      "BATCH_CORRECTION"
    ),
    identical(
      corrected_composition_bundle$batch_contract$design_rank,
      2L
    ),
    identical(
      corrected_composition_bundle$batch_contract$design_columns,
      2L
    ),
    corrected_composition_bundle$batch_contract$design_residual_df > 0L,
    identical(
      corrected_composition_bundle$batch_contract$validation_summary$schema_version,
      1L
    ),
    all(c("batch", "Sample") %in%
      corrected_composition_bundle$batch_contract$required_source_obs_columns),
    isTRUE(corrected_composition_bundle$batch_contract$reserved_obs_absent),
    !grepl(
      "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
      paste(
        vapply(
          corrected_composition_bundle[
            c("ECODA_authors_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
          ],
          function(value) value$batch_contract$correction_formula,
          character(1)
        ),
        collapse = "\n"
      ),
      perl = TRUE
    )
  )

  # Active limma identities are rejected when their advertised design loses
  # rank or has no residual degrees of freedom.
  bad_rank_bundle <- corrected_composition_bundle
  bad_rank_bundle$ECODA_authors_HR$batch_contract$design_rank <- 1L
  write_checked(corrected_composition, bad_rank_bundle)
  expect_fail_capture(
    corrected_batch_args(corrected_composition, "composition"),
    "active limma design rank",
    "design"
  )
  bad_residual_bundle <- corrected_composition_bundle
  bad_residual_bundle$ECODA_authors_HR$batch_contract$design_residual_df <- 0L
  write_checked(corrected_composition, bad_residual_bundle)
  expect_fail_capture(
    corrected_batch_args(corrected_composition, "composition"),
    "active limma residual degrees of freedom",
    "design"
  )
  write_checked(corrected_composition, corrected_composition_bundle)

  # An active limma identity with no estimable technical key is a valid exact
  # no-op, provided its empty effective set and intercept-only dimensions are
  # explicit.
  no_correction_sample_metadata <- corrected_sample_metadata
  no_correction_sample_metadata$batch <- factor(rep("A", nrow(
    no_correction_sample_metadata
  )))
  no_correction_validation <- ecoda_batch_validate_metadata(
    no_correction_sample_metadata,
    batch_keys = "batch",
    sample_col = "Sample"
  )
  no_correction_identity <- function(method_id) {
    model_id <- if (identical(method_id, "Pseudobulk")) {
      "pseudobulk_limma_fixed_effects_v1"
    } else {
      "limma_fixed_effects_v1"
    }
    ecoda_hpc_augment_batch_contract(
      identity = ecoda_hpc_batch_contract_identity(
        batch_keys = "batch",
        sample_col = "Sample",
        method_id = method_id,
        model_id = model_id
      ),
      validation = no_correction_validation,
      method_id = method_id,
      batch_keys = "batch"
    )
  }
  no_correction_bundle <- corrected_composition_bundle
  no_correction_bundle$batch_contract <- no_correction_identity("ECODA_authors_HR")
  for (method_id in c(
    "ECODA_authors_HR",
    "ECODA_seuratres_2",
    "ECODA_authors_HR_NULL"
  )) {
    no_correction_bundle[[method_id]]$batch_contract <-
      no_correction_identity(method_id)
  }
  write_checked(corrected_composition, no_correction_bundle)
  expect_ok(
    corrected_batch_args(corrected_composition, "composition"),
    "active limma all-constant no-op"
  )
  stopifnot(
    identical(
      no_correction_bundle$batch_contract$effective_batch_keys,
      character()
    ),
    identical(
      no_correction_bundle$batch_contract$correction_state,
      "NO_CORRECTION"
    ),
    identical(
      no_correction_bundle$batch_contract$correction_formula,
      "NO_CORRECTION: no estimable technical batch key"
    )
  )
  write_checked(corrected_composition, corrected_composition_bundle)
  rejected_final_variant_args <- corrected_batch_args(
    corrected_composition,
    "composition"
  )
  rejected_final_variant_args[
    match("--analysis-variant", rejected_final_variant_args) + 1L
  ] <- "final"
  expect_fail(
    rejected_final_variant_args,
    "final variant rejects corrected pass"
  )
  corrected_metadata <- file.path(
    corrected_results,
    "Synthetic_batch_effect_corrected_final_metadata.rds"
  )
  write_checked(
    corrected_metadata,
    list(
      batch_contract = corrected_identity("ECODA_authors_HR"),
      labels = structure(
        factor(c("A", "B", "A", "B")),
        names = c("s1", "s2", "s3", "s4")
      ),
      n_cells = 400,
      n_samples = 4,
      cells_per_sample = structure(
        c(100, 100, 100, 100),
        names = c("s1", "s2", "s3", "s4")
      )
    )
  )
  corrected_pseudobulk <- file.path(
    corrected_results,
    "Synthetic_batch_effect_corrected_final_pseudobulk.rds"
  )
  corrected_pseudobulk_result <- list(
    pb = matrix(
      c(1, 2, 3, 4, 5, 6, 7, 8),
      nrow = 4L,
      byrow = TRUE,
      dimnames = list(c("s1", "s2", "s3", "s4"), c("g1", "g2"))
    ),
    time_secs = 0,
    mem_GB = 0,
    batch_contract = corrected_identity("Pseudobulk")
  )
  corrected_pseudobulk_bundle <- list(
    batch_contract = corrected_identity("Pseudobulk"),
    Pseudobulk_hvg2000 = corrected_pseudobulk_result
  )
  write_checked(corrected_pseudobulk, corrected_pseudobulk_bundle)
  corrected_pseudobulk_cache <- file.path(
    corrected_cache,
    "Synthetic_batch_effect_corrected_final_pseudobulk_hvg2000.rds"
  )
  write_checked(corrected_pseudobulk_cache, pb)
  corrected_selection <- file.path(directory, "corrected-final-selection.tsv")
  write_checked_text(
    corrected_selection,
    "Synthetic\tbatch_effect_corrected\tbatch_effect_corrected"
  )
  expect_ok(
    c(
      "--root", corrected_root,
      "--selection", corrected_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "corrected",
      "--analysis-variant", "corrected_final",
      "--config", corrected_config
    ),
    "corrected-final pseudobulk result bundle"
  )
  stopifnot(
    identical(
      corrected_pseudobulk_bundle$batch_contract$model_id,
      "pseudobulk_limma_fixed_effects_v1"
    ),
    identical(
      corrected_pseudobulk_result$batch_contract$correction_mode,
      "limma_fixed_effects_pseudobulk"
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$correction_mode,
      "limma_fixed_effects_pseudobulk"
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$fixed_effect_aliases,
      c(batch = "batch_key_1")
    ),
    identical(
      corrected_pseudobulk_result$batch_contract$correction_design_formula,
      "~1 + batch_key_1"
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$effective_batch_keys,
      "batch"
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$non_estimable_batch_keys,
      character()
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$correction_state,
      "BATCH_CORRECTION"
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$design_rank,
      2L
    ),
    identical(
      corrected_pseudobulk_bundle$batch_contract$design_columns,
      2L
    ),
    corrected_pseudobulk_bundle$batch_contract$design_residual_df > 0L,
    identical(
      corrected_pseudobulk_bundle$batch_contract$validation_summary$schema_version,
      1L
    ),
    all(is.finite(corrected_pseudobulk_result$pb)),
    identical(
      rownames(corrected_pseudobulk_result$pb),
      corrected_sample_metadata$Sample
    ),
    grepl(
      "DESeq2 design=~ 1; model.matrix(~ 1 + batch_key_1); limma::removeBatchEffect(covariates=technical_covariates, design=intercept)",
      corrected_pseudobulk_result$batch_contract$correction_formula,
      fixed = TRUE
    ),
    !grepl(
      "__ecoda_batch_combined_v1|lme4|lmer|\\(1 \\|",
      corrected_pseudobulk_result$batch_contract$correction_formula,
      perl = TRUE
    )
  )
  bad_pseudobulk_rank <- corrected_pseudobulk_bundle
  bad_pseudobulk_rank$Pseudobulk_hvg2000$batch_contract$design_rank <- 1L
  write_checked(corrected_pseudobulk, bad_pseudobulk_rank)
  expect_fail_capture(
    c(
      "--root", corrected_root,
      "--selection", corrected_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "corrected",
      "--analysis-variant", "corrected_final",
      "--config", corrected_config
    ),
    "active pseudobulk limma design rank",
    "design"
  )
  bad_pseudobulk_residual <- corrected_pseudobulk_bundle
  bad_pseudobulk_residual$Pseudobulk_hvg2000$batch_contract$design_residual_df <- 0L
  write_checked(corrected_pseudobulk, bad_pseudobulk_residual)
  expect_fail_capture(
    c(
      "--root", corrected_root,
      "--selection", corrected_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "corrected",
      "--analysis-variant", "corrected_final",
      "--config", corrected_config
    ),
    "active pseudobulk limma residual degrees of freedom",
    "design"
  )
  write_checked(corrected_pseudobulk, corrected_pseudobulk_bundle)
  unlink(corrected_pseudobulk)
  expect_fail(
    c(
      "--root", corrected_root,
      "--selection", corrected_selection,
      "--labels", "pseudobulk",
      "--batch-pass", "corrected",
      "--analysis-variant", "corrected_final",
      "--config", corrected_config
    ),
    "corrected-final pseudobulk cache is not a result bundle"
  )
  write_checked(corrected_pseudobulk, corrected_pseudobulk_bundle)
  wrong_corrected_stem <- file.path(
    corrected_results,
    "Synthetic_batch_effect_corrected_composition.rds"
  )
  write_checked(wrong_corrected_stem, corrected_composition_bundle)
  expect_fail(
    corrected_batch_args(wrong_corrected_stem, "composition"),
    "corrected-final result stem"
  )
  wrong_corrected_root <- file.path(directory, "batch_effect", "uncorrected_final")
  wrong_root_results <- file.path(wrong_corrected_root, "results")
  dir.create(wrong_root_results, recursive = TRUE)
  wrong_root_composition <- file.path(
    wrong_root_results,
    "Synthetic_batch_effect_corrected_final_composition.rds"
  )
  wrong_root_metadata <- file.path(
    wrong_root_results,
    "Synthetic_batch_effect_corrected_final_metadata.rds"
  )
  write_checked(wrong_root_composition, corrected_composition_bundle)
  write_checked(
    wrong_root_metadata,
    list(
      batch_contract = corrected_identity("ECODA_authors_HR"),
      labels = structure(
        factor(c("A", "B", "A", "B")),
        names = c("s1", "s2", "s3", "s4")
      ),
      n_cells = 400,
      n_samples = 4,
      cells_per_sample = structure(
        c(100, 100, 100, 100),
        names = c("s1", "s2", "s3", "s4")
      )
    )
  )
  expect_fail(
    c(
      "--root", wrong_corrected_root,
      "--selection", corrected_selection,
      "--labels", "composition",
      "--batch-pass", "corrected",
      "--analysis-variant", "corrected_final",
      "--config", corrected_config
    ),
    "corrected-final root"
  )
  bad_corrected_bundle <- corrected_composition_bundle
  bad_corrected_bundle$ECODA_authors_HR_NULL <- NULL
  write_checked(corrected_composition, bad_corrected_bundle)
  expect_fail(
    corrected_batch_args(corrected_composition, "composition"),
    "corrected-final composition requires explicit null bundle key"
  )
  write_checked(corrected_composition, corrected_composition_bundle)
  wrong_key_bundle <- corrected_composition_bundle
  wrong_key_bundle$ECODA_authors_HR_NULL$batch_contract <- corrected_identity(
    "ECODA_authors_HR"
  )
  write_checked(corrected_composition, wrong_key_bundle)
  expect_fail(
    corrected_batch_args(corrected_composition, "composition"),
    "corrected-final composition null key identity"
  )
  write_checked(corrected_composition, corrected_composition_bundle)


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
  unrelated_partial <- file.path(scitd_root, "unrelated", "nested", "stale.tmp.123")
  dir.create(dirname(unrelated_partial), recursive = TRUE)
  writeLines("stale", unrelated_partial)
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
  selected_partial <- file.path(
    scitd_root, "results", "Synthetic_scitd.rds.tmp.123"
  )
  writeLines("stale", selected_partial)
  selected_partial_result <- run_validator_capture(c(
    "--root", scitd_root,
    "--selection", scitd_selection,
    "--labels", "scitd",
    "--input-root", scratch,
    "--config", config
  ))
  if (identical(selected_partial_result$status, 0L) ||
      !grepl("partial benchmark artifacts remain",
             selected_partial_result$output, fixed = TRUE)) {
    stop("selected adjacent partial was accepted: ",
         selected_partial_result$output)
  }
  unlink(selected_partial)

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
