#!/usr/bin/env Rscript
# Validator-only corrected-final consumer contract barrier.
# It reads H5AD metadata/identity and exported sample metadata only; it never
# submits workers or materializes expression/count matrices.

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT is required")
}
project_root <- normalizePath(project_root, mustWork = TRUE)
setwd(project_root)

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))

args <- parse_flags(commandArgs(trailingOnly = TRUE))
required_args <- c("config", "selection", "analysis-root", "input-root", "output")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0L) {
  stop("missing required arguments: ", paste(missing_args, collapse = ", "))
}

config_path <- normalizePath(args[["config"]], mustWork = TRUE)
selection_path <- normalizePath(args[["selection"]], mustWork = TRUE)
analysis_root <- normalizePath(args[["analysis-root"]], mustWork = FALSE)
input_root <- normalizePath(args[["input-root"]], mustWork = FALSE)
output_path <- normalizePath(args[["output"]], mustWork = FALSE)
if (!grepl("/batch_effect/corrected_final$", analysis_root, perl = TRUE)) {
  stop("corrected-final consumer barrier requires batch_effect/corrected_final")
}
if (!grepl("^/", analysis_root, perl = TRUE) ||
    !grepl("^/", input_root, perl = TRUE) ||
    !grepl("^/", output_path, perl = TRUE)) {
  stop("consumer barrier paths must be absolute")
}

config <- read_datasets_json(config_path, view = "batch_effect_corrected")
selection_lines <- readLines(selection_path, warn = FALSE)
if (!length(selection_lines) || any(!nzchar(selection_lines))) {
  stop("corrected-final consumer selection is empty or malformed")
}
selection_parts <- strsplit(selection_lines, "\t", fixed = TRUE)
if (any(lengths(selection_parts) != 3L)) {
  stop("corrected-final consumer selection must have exactly three columns")
}
selection <- data.frame(
  dataset = vapply(selection_parts, `[[`, character(1), 1L),
  view = vapply(selection_parts, `[[`, character(1), 2L),
  label = vapply(selection_parts, `[[`, character(1), 3L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (any(selection$view != "batch_effect_corrected") ||
    any(selection$label != "batch_effect_corrected") ||
    any(!nzchar(selection$dataset)) ||
    anyDuplicated(selection$dataset)) {
  stop("corrected-final consumer selection has invalid rows")
}

method_specs <- list(
  prepare_pseudobulk = list(
    method_id = "Pseudobulk", model_id = "pseudobulk_composite_v1"
  ),
  pseudobulk = list(
    method_id = "Pseudobulk", model_id = "pseudobulk_composite_v1"
  ),
  gloscope = list(
    method_id = "GloScope", model_id = "embedding_consumer_harmony_v1"
  ),
  composition = list(
    method_id = "ECODA_authors_HR",
    model_id = "ecoda_additive_random_intercepts_v1"
  ),
  mrvi = list(
    method_id = "MrVI", model_id = "mrvi_composite_v1"
  ),
  pilot = list(
    method_id = "PILOT", model_id = "embedding_consumer_harmony_v1"
  ),
  qot = list(
    method_id = "QOT", model_id = "embedding_consumer_harmony_v1"
  )
)

rows <- vector("list", nrow(selection))
failures <- character()
for (index in seq_len(nrow(selection))) {
  ds <- selection$dataset[[index]]
  view <- selection$view[[index]]
  row_result <- list(
    dataset = ds,
    view = view,
    status = "FAIL"
  )
  tryCatch({
    entry <- config[[ds]]
    if (is.null(entry)) stop("dataset is absent from datasets.json")
    batch_keys <- ecoda_hpc_normalize_batch_keys(
      entry$batch_col,
      sample_col = "Sample",
      biological_label = entry$label_col
    )
    h5ad_path <- get_h5ad_path(
      config, ds, view, file.path(input_root, ds, "output")
    )
    h5ad_identity <- ecoda_hpc_batch_contract_identity(
      batch_keys,
      sample_col = "Sample",
      method_id = "preprocess",
      model_id = "hvg_composite_v1"
    )
    h5ad_metadata <- load_h5ad_pseudobulk_metadata(
      h5ad_path = h5ad_path,
      sample_col = "Sample",
      metadata_columns = unique(c(entry$label_col, batch_keys)),
      n_hvg = 1L,
      required_nonmissing_columns = unique(c(
        "Sample", entry$label_col, batch_keys
      )),
      expected_batch_contract = h5ad_identity,
      view = view,
      method = "preprocessing",
      allow_missing_summary = TRUE
    )
    h5ad_obs <- h5ad_metadata$obs
    if (!"Sample" %in% colnames(h5ad_obs)) {
      stop("H5AD metadata lacks standardized Sample")
    }
    expected_sample_ids <- unique(as.character(h5ad_obs[["Sample"]]))
    metadata_path <- file.path(
      analysis_root, "metadata", paste0(ds, "_sample_metadata.feather")
    )
    allow_unknown_keys <- if (identical(ds, "Breast_cancer")) {
      intersect(batch_keys, "suspension_dissociation_time")
    } else {
      character()
    }
    metadata_context <- ecoda_hpc_load_sample_metadata_contract(
      path = metadata_path,
      expected_sample_ids = expected_sample_ids,
      batch_keys = batch_keys,
      sample_col = "Sample",
      biological_label = entry$label_col,
      required_columns = entry$label_col,
      allow_unknown_keys = allow_unknown_keys
    )
    validation <- metadata_context$validation
    level_counts <- vapply(
      batch_keys,
      function(key) length(validation$per_key_levels[[key]]),
      integer(1)
    )
    effective_batch_keys <- validation[["effective_batch_keys"]]
    if (is.null(effective_batch_keys)) {
      effective_batch_keys <- unname(batch_keys[level_counts >= 2L])
    }
    effective_batch_keys <- unname(as.character(effective_batch_keys))
    non_estimable_batch_keys <- unname(
      setdiff(batch_keys, effective_batch_keys)
    )
    correction_state <- if (length(effective_batch_keys)) {
      "BATCH_CORRECTION"
    } else {
      "NO_CORRECTION"
    }
    estimable_keys <- effective_batch_keys
    checked_r_consumers <- character()
    consumer_contracts <- list()
    for (method in c("prepare_pseudobulk", "pseudobulk", "gloscope", "composition")) {
      spec <- method_specs[[method]]
      identity <- ecoda_hpc_batch_contract_identity(
        batch_keys,
        sample_col = "Sample",
        method_id = spec$method_id,
        model_id = spec$model_id
      )
      consumer_contracts[[method]] <- ecoda_hpc_augment_batch_contract(
        identity = identity,
        validation = validation,
        method_id = spec$method_id,
        batch_keys = batch_keys,
        scalar_batch_col = metadata_context$scalar_batch_col
      )
      checked_r_consumers <- c(checked_r_consumers, method)
    }
    method_identities <- lapply(method_specs, function(spec) {
      ecoda_hpc_batch_contract_identity(
        batch_keys,
        sample_col = "Sample",
        method_id = spec$method_id,
        model_id = spec$model_id
      )
    })
    row_result <- list(
      dataset = ds,
      view = view,
      status = "OK",
      h5ad_path = h5ad_path,
      metadata_path = metadata_path,
      batch_keys = unname(batch_keys),
      key_level_counts = as.list(level_counts),
      estimable_batch_keys = unname(estimable_keys),
      non_estimable_batch_keys = unname(non_estimable_batch_keys),
      correction_state = correction_state,
      composite_level_count = length(validation$composite_levels),
      checked_r_consumers = checked_r_consumers,
      correction_formulas = lapply(
        consumer_contracts,
        function(identity) identity$correction_formula
      ),
      method_identity_fingerprints = lapply(
        method_identities,
        function(identity) identity$fingerprint
      )
    )
  }, error = function(error) {
    message("corrected-final consumer contract failed for ", ds, ": ",
            conditionMessage(error))
    failures <<- c(failures, paste0(ds, ": ", conditionMessage(error)))
    row_result$error <- conditionMessage(error)
  })
  rows[[index]] <- row_result
}

report <- list(
  schema_version = 1L,
  status = if (length(failures)) "FAIL" else "CORRECTED_FINAL_CONSUMERS_VALIDATED",
  analysis_variant = "corrected_final",
  analysis_pass = "corrected",
  run_id = Sys.getenv("ECODA_RUN_ID", unset = ""),
  config_path = config_path,
  selection_path = selection_path,
  analysis_root = analysis_root,
  input_root = input_root,
  rows = rows,
  failures = unname(failures)
)

output_dir <- dirname(output_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(output_path, ".tmp.", Sys.getpid())
jsonlite::write_json(report, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null")
if (!file.rename(tmp, output_path)) {
  unlink(tmp, force = TRUE)
  stop("could not atomically install consumer contract report")
}
digest <- tolower(unname(tools::md5sum(output_path)))
sidecar <- paste0(output_path, ".md5")
sidecar_tmp <- paste0(sidecar, ".tmp.", Sys.getpid())
writeLines(c(
  paste0("MD5=", digest),
  paste0("SIZE=", file.info(output_path)$size),
  paste0("PATH=", output_path)
), sidecar_tmp, useBytes = TRUE)
if (!file.rename(sidecar_tmp, sidecar)) {
  unlink(sidecar_tmp, force = TRUE)
  stop("could not atomically install consumer contract checksum")
}

if (length(failures)) {
  quit(status = 1L, save = "no")
}
cat("corrected-final consumer contracts: OK\n")
