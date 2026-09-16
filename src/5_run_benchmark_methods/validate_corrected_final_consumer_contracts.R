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
if (!grepl("/batch_effect/corrected_final(/recovery_35row)?$", analysis_root, perl = TRUE)) {
  stop("corrected-final consumer barrier requires batch_effect/corrected_final or its recovery_35row replacement root")
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
method_matrix_mode <- "method-matrix" %in% names(args)
method_matrix_path <- NULL
method_matrix_rows <- NULL
method_matrix_md5 <- NULL
method_matrix_size <- NULL
method_matrix_sha256 <- NULL
method_matrix_identity <- NULL
method_matrix_declared_count <- NULL
if (method_matrix_mode) {
  raw_matrix_path <- args[["method-matrix"]]
  if (
    isTRUE(raw_matrix_path) ||
    length(raw_matrix_path) != 1L ||
    is.na(raw_matrix_path) ||
    !nzchar(as.character(raw_matrix_path))
  ) {
    stop("--method-matrix requires a readable path")
  }
  method_matrix_path <- normalizePath(
    as.character(raw_matrix_path), mustWork = TRUE
  )
  if (
    !grepl("^/", method_matrix_path, perl = TRUE) ||
    !file.exists(method_matrix_path) ||
    isTRUE(file.info(method_matrix_path)$isdir) ||
    nzchar(Sys.readlink(method_matrix_path))
  ) {
    stop("--method-matrix must be an absolute regular file, not a symlink")
  }
  matrix_sidecar <- paste0(method_matrix_path, ".md5")
  if (!file.exists(matrix_sidecar) || nzchar(Sys.readlink(matrix_sidecar))) {
    stop("run-owned method matrix checksum sidecar is missing or unsafe")
  }
  sidecar_lines <- readLines(matrix_sidecar, warn = FALSE)
  sidecar_md5_values <- sub("^MD5=", "", sidecar_lines[grepl("^MD5=", sidecar_lines)])
  sidecar_size_values <- sub("^SIZE=", "", sidecar_lines[grepl("^SIZE=", sidecar_lines)])
  sidecar_path_values <- sub("^PATH=", "", sidecar_lines[grepl("^PATH=", sidecar_lines)])
  if (
    length(sidecar_md5_values) != 1L ||
    length(sidecar_size_values) != 1L ||
    length(sidecar_path_values) != 1L
  ) {
    stop("run-owned method matrix checksum sidecar is malformed")
  }
  sidecar_md5 <- sidecar_md5_values[[1L]]
  sidecar_size <- sidecar_size_values[[1L]]
  sidecar_path <- sidecar_path_values[[1L]]
  method_matrix_md5 <- tolower(unname(tools::md5sum(method_matrix_path)))
  method_matrix_size <- as.numeric(file.info(method_matrix_path)$size)
  sidecar_path_normalized <- normalizePath(sidecar_path, mustWork = FALSE)
  if (
    !identical(sidecar_md5, method_matrix_md5) ||
    !identical(suppressWarnings(as.numeric(sidecar_size)), method_matrix_size) ||
    !identical(sidecar_path_normalized, method_matrix_path)
  ) {
    stop("run-owned method matrix checksum sidecar does not match the matrix")
  }
  if (!grepl("/batch_effect/corrected_final/recovery_35row$", analysis_root, perl = TRUE)) {
    stop("method-matrix corrected-final validation requires the recovery_35row root")
  }
  matrix_lines <- readLines(method_matrix_path, warn = FALSE)
  if (!length(matrix_lines) || any(!nzchar(matrix_lines))) {
    stop("method matrix is empty or contains blank rows")
  }
  matrix_parts <- strsplit(matrix_lines, "\t", fixed = TRUE)
  if (any(lengths(matrix_parts) != 3L)) {
    stop("method matrix must have exactly three columns")
  }
  matrix_values <- data.frame(
    dataset = vapply(matrix_parts, `[[`, character(1), 1L),
    view = vapply(matrix_parts, `[[`, character(1), 2L),
    method = vapply(matrix_parts, `[[`, character(1), 3L),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expected_matrix_methods <- c(
    "prepare_pseudobulk", "pseudobulk", "gloscope", "composition",
    "mrvi", "pilot", "qot"
  )
  if (
    any(!nzchar(matrix_values$dataset)) ||
    any(matrix_values$view != "batch_effect_corrected") ||
    any(!matrix_values$method %in% expected_matrix_methods) ||
    anyDuplicated(paste(matrix_values$dataset, matrix_values$view,
                        matrix_values$method, sep = "|"))
  ) {
    stop("method matrix has malformed, unsupported, or duplicate rows")
  }
  if (
    !setequal(unique(matrix_values$dataset), unique(selection$dataset)) ||
    any(!selection$dataset %in% matrix_values$dataset)
  ) {
    stop("method matrix does not cover exactly the selected datasets")
  }
  method_matrix_rows <- matrix_values
  method_matrix_declared_count <- nrow(matrix_values)
  method_matrix_md5 <- tolower(unname(tools::md5sum(method_matrix_path)))
  method_matrix_size <- as.numeric(file.info(method_matrix_path)$size)
  method_matrix_sha256_arg <- args[["method-matrix-sha256"]]
  method_matrix_identity_arg <- args[["method-matrix-identity"]]
  if (!is.null(args[["method-matrix-md5"]])) {
    supplied_md5 <- as.character(args[["method-matrix-md5"]])
    if (!identical(supplied_md5, method_matrix_md5)) {
      stop("method matrix MD5 does not match the bound file")
    }
  }
  if (!is.null(args[["method-matrix-size"]])) {
    supplied_size <- suppressWarnings(as.numeric(args[["method-matrix-size"]]))
    if (length(supplied_size) != 1L || is.na(supplied_size) ||
        supplied_size != method_matrix_size) {
      stop("method matrix size does not match the bound file")
    }
  }
  sha256_bin <- Sys.which("sha256sum")
  sha256_args <- character()
  if (!nzchar(sha256_bin)) {
    sha256_bin <- Sys.which("shasum")
    sha256_args <- "-a 256"
  }
  if (!nzchar(sha256_bin)) {
    stop("a SHA-256 utility is required for method matrix identity")
  }
  sha256_output <- system2(
    sha256_bin, c(sha256_args, method_matrix_path),
    stdout = TRUE, stderr = TRUE
  )
  computed_sha256 <- strsplit(trimws(sha256_output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[[:xdigit:]]{64}$", computed_sha256)) {
    stop("could not compute method matrix SHA-256 identity")
  }
  method_matrix_sha256 <- tolower(computed_sha256)
  if (!is.null(method_matrix_sha256_arg)) {
    supplied_sha256 <- as.character(method_matrix_sha256_arg)
    if (!identical(supplied_sha256, method_matrix_sha256)) {
      stop("method matrix SHA-256 does not match the bound file")
    }
  }
  if (!is.null(method_matrix_identity_arg)) {
    method_matrix_identity <- as.character(method_matrix_identity_arg)
    if (length(method_matrix_identity) != 1L ||
        !grepl("^[[:xdigit:]]{64}$", method_matrix_identity)) {
      stop("method matrix identity is malformed")
    }
    if (!is.null(method_matrix_sha256) &&
        !identical(method_matrix_identity, method_matrix_sha256)) {
      stop("method matrix identity disagrees with its SHA-256 identity")
    }
  } else {
    method_matrix_identity <- method_matrix_sha256
  }
  if (is.null(method_matrix_identity)) {
    # The wrapper records the SHA-256 identity.  An isolated validator can still
    # enforce the declared rows and file checksum without requiring a
    # platform-specific SHA utility.
    method_matrix_identity <- method_matrix_md5
  }
}
method_matrix_methods_for <- function(dataset) {
  if (!method_matrix_mode) {
    return(c("prepare_pseudobulk", "pseudobulk", "gloscope", "composition"))
  }
  unname(method_matrix_rows$method[method_matrix_rows$dataset == dataset])
}

method_specs <- list(
  prepare_pseudobulk = list(
    method_id = "Pseudobulk", model_id = "pseudobulk_limma_fixed_effects_v1"
  ),
  pseudobulk = list(
    method_id = "Pseudobulk", model_id = "pseudobulk_limma_fixed_effects_v1"
  ),
  gloscope = list(
    method_id = "GloScope", model_id = "embedding_consumer_harmony_v1"
  ),
  composition = list(
    method_id = "ECODA_authors_HR", model_id = "limma_fixed_effects_v1"
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
validate_fixed_effect_consumer <- function(
  contract,
  spec,
  batch_keys,
  biological_label,
  label
) {
  if (!is.list(contract)) stop(label, " contract is not a list")
  if (!identical(contract[["model_id"]], spec$model_id)) {
    stop(label, " has the wrong active corrected model identity")
  }
  required <- c(
    "effective_batch_keys", "non_estimable_batch_keys", "correction_state",
    "correction_mode", "correction_formula", "fixed_effect_aliases",
    "correction_design_formula", "design_rank", "design_columns",
    "design_residual_df"
  )
  missing <- setdiff(required, names(contract))
  if (length(missing)) {
    stop(label, " is missing fixed-effect metadata: ",
         paste(missing, collapse = ", "))
  }
  effective <- unname(as.character(contract[["effective_batch_keys"]]))
  non_estimable <- unname(as.character(contract[["non_estimable_batch_keys"]]))
  keys <- unname(as.character(batch_keys))
  if (
    anyNA(effective) || anyNA(non_estimable) ||
    any(!effective %in% keys) || any(!non_estimable %in% keys) ||
    anyDuplicated(effective) || anyDuplicated(non_estimable) ||
    !identical(effective, keys[keys %in% effective]) ||
    !identical(non_estimable, keys[keys %in% non_estimable]) ||
    !setequal(c(effective, non_estimable), keys) ||
    length(intersect(effective, non_estimable))
  ) {
    stop(label, " has invalid effective/non-estimable technical keys")
  }
  expected_aliases <- if (length(effective)) {
    setNames(paste0("batch_key_", match(effective, keys)), effective)
  } else {
    character()
  }
  aliases <- contract[["fixed_effect_aliases"]]
  if (
    !is.character(aliases) ||
    !identical(names(aliases), names(expected_aliases)) ||
    !identical(unname(aliases), unname(expected_aliases))
  ) {
    stop(label, " has the wrong fixed-effect aliases")
  }
  expected_spec <- ecoda_batch_correction_spec(
    method_id = spec$method_id,
    batch_keys = keys,
    effective_batch_keys = effective,
    non_estimable_batch_keys = non_estimable
  )
  expected_state <- if (length(effective)) "BATCH_CORRECTION" else "NO_CORRECTION"
  if (
    !identical(contract[["correction_state"]], expected_state) ||
    !identical(contract[["correction_mode"]], expected_spec$correction_mode) ||
    !identical(contract[["correction_formula"]], expected_spec$correction_formula)
  ) {
    stop(label, " has the wrong fixed-effect correction policy")
  }
  expected_design <- if (length(effective)) {
    paste0("~1 + ", paste(unname(expected_aliases), collapse = " + "))
  } else {
    "~1"
  }
  if (
    !is.character(contract[["correction_design_formula"]]) ||
    length(contract[["correction_design_formula"]]) != 1L ||
    gsub("[[:space:]]", "", contract[["correction_design_formula"]]) !=
      gsub("[[:space:]]", "", expected_design)
  ) {
    stop(label, " has the wrong separate-covariate design formula")
  }
  for (field in c("design_rank", "design_columns", "design_residual_df")) {
    value <- contract[[field]]
    if (
      !is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value != floor(value)
    ) {
      stop(label, " has invalid ", field)
    }
  }
  if (
    contract[["design_rank"]] != contract[["design_columns"]] ||
    contract[["design_residual_df"]] <= 0
  ) {
    stop(label, " design is not full rank with positive residual degrees of freedom")
  }
  if (contract[["design_rank"]] > contract[["design_columns"]]) {
    stop(label, " design rank exceeds design columns")
  }
  formula <- contract[["correction_formula"]]
  if (
    grepl(
      "__ecoda_batch_combined_v1|additive_random_intercepts|pseudobulk_composite|lmFit|remove technical contribution|\\(1 \\|",
      formula,
      perl = TRUE
    ) ||
    (nzchar(biological_label) && (
      grepl(biological_label, formula, fixed = TRUE) ||
      grepl(
        biological_label,
        contract[["correction_design_formula"]],
        fixed = TRUE
      )
    ))
  ) {
    stop(label, " advertises a historical/scalarized or biological design")
  }
  invisible(TRUE)
}

failures <- character()
for (index in seq_len(nrow(selection))) {
  ds <- selection$dataset[[index]]
  view <- selection$view[[index]]
  row_result <- list(
    dataset = ds,
    view = view,
    status = "FAIL"
  )
  consumer_contracts <- list()
  checked_r_consumers <- character()
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
    declared_methods <- method_matrix_methods_for(ds)
    if (
      !length(declared_methods) ||
      any(!declared_methods %in% names(method_specs)) ||
      anyDuplicated(declared_methods)
    ) {
      stop("method scope is missing, duplicated, or unsupported for ", ds)
    }
    for (method in declared_methods) {
      spec <- method_specs[[method]]
      identity <- ecoda_hpc_batch_contract_identity(
        batch_keys,
        sample_col = "Sample",
        method_id = spec$method_id,
        model_id = spec$model_id
      )
      consumer_contract <- ecoda_hpc_augment_batch_contract(
        identity = identity,
        validation = validation,
        method_id = spec$method_id,
        batch_keys = batch_keys,
        scalar_batch_col = metadata_context$scalar_batch_col
      )
      if (method %in% c("prepare_pseudobulk", "pseudobulk", "composition")) {
        validate_fixed_effect_consumer(
          consumer_contract,
          spec,
          batch_keys,
          entry$label_col,
          paste0(ds, "/", view, "/", method)
        )
      }
      consumer_contracts[[method]] <- consumer_contract
      checked_r_consumers <- c(checked_r_consumers, method)
    }
    method_identities <- lapply(
      if (method_matrix_mode) declared_methods else names(method_specs),
      function(method) {
        spec <- method_specs[[method]]
        ecoda_hpc_batch_contract_identity(
          batch_keys,
          sample_col = "Sample",
          method_id = spec$method_id,
          model_id = spec$model_id
        )
      }
    )
    row_result <- list(
      dataset = ds,
      view = view,
      status = "OK",
      h5ad_path = h5ad_path,
      metadata_path = metadata_path,
      declared_methods = declared_methods,
      estimable_batch_keys = unname(effective_batch_keys),
      key_level_counts = as.list(level_counts),
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
if (method_matrix_mode) {
  report$method_matrix <- list(
    path = method_matrix_path,
    md5 = method_matrix_md5,
    size = method_matrix_size,
    sha256 = method_matrix_sha256,
    identity = method_matrix_identity,
    declared_method_rows = method_matrix_declared_count,
    dataset_order = unname(unique(selection$dataset))
  )
}

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
