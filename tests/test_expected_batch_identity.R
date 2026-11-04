#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
Sys.setenv(PROJECT_ROOT = root,
           RETICULATE_PYTHON = file.path(root, ".pixi/envs/default/bin/python"))
source(file.path(root, "src/utils/batch_contract.R"))
source(file.path(root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
validator <- new.env(parent = globalenv())
expressions <- parse(file.path(root, "src/5_run_benchmark_methods/validate_benchmark_rds_contract.R"))
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      ((is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function"))) ||
       (is.symbol(expression[[2L]]) && startsWith(as.character(expression[[2L]]), ".batch_")))) {
    eval(expression, envir = validator)
  }
}
config <- tempfile("expected-identity-", fileext = ".json")
writeLines('{"Synthetic":{"columns":{"sample":"Sample","batch":"tech"},"views":{"batch_effect_corrected":{"output_file_name":"synthetic.h5ad"}}}}', config)
validator$config_path <- config
validator$config <- jsonlite::fromJSON(config, simplifyVector = FALSE)
validator$batch_pass <- "corrected"
validator$corrected_final_mode <- TRUE
validator$corrected_summary_required <- FALSE
metadata <- data.frame(Sample = paste0("s", 1:6), tech = rep(c("a", "b"), each = 3L))
validation <- ecoda_batch_validate_metadata(metadata, batch_keys = "tech", sample_col = "Sample")
for (method in c("prepare_pseudobulk", "pseudobulk", "composition")) {
  expected <- validator$.corrected_batch_row_identity("Synthetic", "batch_effect_corrected", method)
  supplied <- validator$.expected_row_batch_contract("Synthetic", "batch_effect_corrected", method, expected)
  recorded <- ecoda_hpc_augment_batch_contract(expected, validation,
      method_id = expected$method_id, batch_keys = "tech")
  checked <- validator$validate_batch_contract_identity(supplied, recorded,
      require_recorded = TRUE, require_summary = FALSE)
  stopifnot(identical(checked$fingerprint, recorded$fingerprint))
  missing_fit <- tryCatch(validator$validate_batch_contract_identity(supplied, expected,
      require_recorded = TRUE, require_summary = FALSE), error = identity)
  stopifnot(inherits(missing_fit, "error"))
  other <- ecoda_batch_contract_identity(batch_keys = "other_tech", sample_col = "Sample",
      method_id = expected$method_id, model_id = expected$model_id)
  mismatch <- tryCatch(validator$.expected_row_batch_contract(
      "Synthetic", "batch_effect_corrected", method, other), error = identity)
  stopifnot(inherits(mismatch, "error"))
}
unlink(config)
cat("PASS configuration-only expected identities, strict materialized fit metadata and source mismatch rejection\n")
