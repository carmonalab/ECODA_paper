# R adapter for the authoritative corrected-mode Python batch contract.
#
# Python owns normalization, canonicalization, validation, fingerprints, and
# serialized identities. R keeps only the fixed-effect model-matrix boundary
# needed by limma consumers.

.ecoda_batch_source_file <- tryCatch(
  sys.frame(1L)$ofile,
  error = function(error) ""
)

.ecoda_batch_python_module <- local({
  module <- NULL
  source_file <- .ecoda_batch_source_file
  function() {
    if (!is.null(module)) return(module)
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("reticulate is required for the authoritative Python batch contract")
    }
    project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
    source_dir <- if (
      is.character(source_file) &&
      length(source_file) == 1L &&
      nzchar(source_file)
    ) {
      dirname(normalizePath(source_file, mustWork = FALSE))
    } else {
      ""
    }
    if (!nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
      python_candidates <- c(
        if (nzchar(project_root)) file.path(
          project_root, ".pixi", "envs", "default", "bin", "python"
        ),
        if (nzchar(source_dir)) file.path(
          source_dir, "..", "..", ".pixi", "envs", "default", "bin", "python"
        )
      )
      python_candidates <- python_candidates[file.exists(python_candidates)]
      if (length(python_candidates)) {
        Sys.setenv(
          RETICULATE_PYTHON = normalizePath(
            python_candidates[[1L]], mustWork = TRUE
          )
        )
      }
    }
    candidates <- c(
      if (nzchar(source_dir)) file.path(source_dir, "py", "batch_contract.py"),
      if (nzchar(project_root)) {
        file.path(project_root, "src", "utils", "py", "batch_contract.py")
      },
      file.path(getwd(), "src", "utils", "py", "batch_contract.py")
    )
    candidates <- unique(candidates[file.exists(candidates)])
    if (!length(candidates)) {
      stop("authoritative Python batch contract not found")
    }
    module_dir <- dirname(normalizePath(candidates[[1L]], mustWork = TRUE))
    python_sys <- reticulate::import("sys", convert = FALSE)
    python_sys$path$insert(0L, module_dir)
    module <<- reticulate::import_from_path(
      "batch_contract",
      path = module_dir,
      convert = FALSE
    )
    module
  }
})

.ecoda_batch_python_keys <- function(batch_keys) {
  if (is.character(batch_keys) && length(batch_keys) > 1L) {
    return(as.list(unname(batch_keys)))
  }
  if (is.list(batch_keys)) return(unname(batch_keys))
  batch_keys
}

.ecoda_batch_python_value <- function(value) {
  if (is.factor(value)) as.character(value) else value
}


.ecoda_batch_identity_python <- function(identity) {
  if (!is.list(identity) || is.null(names(identity))) {
    stop("corrected batch identity must be a named list")
  }
  payload <- identity
  for (field in c("ordered_source_keys", "required_source_obs_columns")) {
    if (field %in% names(payload)) {
      payload[[field]] <- as.list(unname(as.character(payload[[field]])))
    }
  }
  reticulate::r_to_py(payload)
}
.ecoda_batch_summary_python <- function(summary) {
  if (!is.list(summary) || is.null(names(summary))) {
    stop("compact batch summary must be a named list")
  }
  payload <- summary
  for (field in c("sample_constancy", "per_key_levels", "key_level_counts")) {
    values <- payload[[field]]
    if (!is.list(values) || is.null(names(values))) {
      stop("compact batch summary has invalid ", field)
    }
    payload[[field]] <- setNames(lapply(values, function(value) {
      if (field == "sample_constancy") {
        isTRUE(as.logical(value[[1L]] %||% value))
      } else if (field == "key_level_counts") {
        as.integer(value[[1L]] %||% value)
      } else {
        as.list(unname(as.character(value)))
      }
    }), names(values))
  }
  payload[["composite_levels"]] <- as.list(unname(as.character(
    payload[["composite_levels"]]
  )))
  reticulate::r_to_py(payload)
}
.ecoda_batch_accepted_sentinels_python <- function(values) {
  if (is.null(values)) return(NULL)
  if (!is.list(values) || is.null(names(values))) {
    stop("accepted sentinel values must be a named list")
  }
  reticulate::r_to_py(setNames(lapply(values, function(value) {
    as.list(unname(as.character(value)))
  }), names(values)))
}

.ecoda_batch_to_r <- function(value) reticulate::py_to_r(value)

.ecoda_batch_validation_keys <- function(validation) {
  keys <- validation[["ordered_keys"]]
  if (is.null(keys)) keys <- validation[["ordered_source_keys"]]
  if (is.null(keys)) keys <- validation[["keys"]]
  if (is.null(keys)) stop("validated batch metadata has no ordered keys")
  unname(as.character(keys))
}

# Python owns the type and byte-level rules. These small adapters only preserve
# R's scalar/list calling convention at the language seam.
ecoda_batch_normalize_keys <- function(
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL
) {
  module <- .ecoda_batch_python_module()
  result <- module$normalize_batch_keys(
    reticulate::r_to_py(.ecoda_batch_python_keys(batch_keys)),
    sample_column = sample_col,
    biological_column = biological_label
  )
  unname(as.character(.ecoda_batch_to_r(result)))
}

ecoda_batch_canonical_value <- function(value, label = "batch value") {
  module <- .ecoda_batch_python_module()
  factor_value <- is.factor(value)
  value <- .ecoda_batch_python_value(value)
  unname(as.character(.ecoda_batch_to_r(module$canonicalize_batch_value(
    reticulate::r_to_py(value),
    factor = factor_value
  ))))
}

ecoda_batch_canonical_values <- function(
  values,
  label = "batch column",
  accepted_sentinel_values = NULL
) {
  module <- .ecoda_batch_python_module()
  values_for_python <- if (is.factor(values)) {
    as.list(as.character(values))
  } else {
    as.list(values)
  }
  result <- module$canonicalize_batch_values(
    reticulate::r_to_py(values_for_python),
    factor = is.factor(values),
    accepted_sentinel_values = if (is.null(accepted_sentinel_values)) {
      NULL
    } else {
      reticulate::r_to_py(as.list(unname(as.character(
        accepted_sentinel_values
      ))))
    }
  )
  unname(as.character(.ecoda_batch_to_r(result)))
}

ecoda_batch_composite_token <- function(
  batch_keys,
  values,
  accepted_sentinel_values = NULL
) {
  module <- .ecoda_batch_python_module()
  keys <- .ecoda_batch_python_keys(batch_keys)
  values_for_python <- if (is.list(values)) {
    lapply(values, .ecoda_batch_python_value)
  } else if (is.factor(values)) {
    as.list(as.character(values))
  } else {
    values
  }
  factors <- if (is.list(values)) {
    vapply(values, is.factor, logical(1L))
  } else if (is.factor(values)) {
    rep(TRUE, length(values))
  } else {
    rep(FALSE, length(values))
  }
  result <- module$composite_token(
    reticulate::r_to_py(keys),
    reticulate::r_to_py(values_for_python),
    factors = reticulate::r_to_py(as.list(factors)),
    accepted_sentinel_values = .ecoda_batch_accepted_sentinels_python(
      accepted_sentinel_values
    )
  )
  unname(as.character(.ecoda_batch_to_r(result)))
}

.ecoda_batch_design_from_validation <- function(
  metadata,
  validation,
  keys,
  sample_col = "Sample",
  label = "limma batch design"
) {
  if (!is.data.frame(metadata)) {
    stop(label, " requires a data.frame")
  }
  if (!is.character(sample_col) || length(sample_col) != 1L ||
      is.na(sample_col) || !nzchar(sample_col) ||
      !sample_col %in% colnames(metadata)) {
    stop(label, " is missing ", sample_col)
  }
  sample_ids <- unname(as.character(metadata[[sample_col]]))
  if (!length(sample_ids) || anyNA(sample_ids) ||
      any(!nzchar(trimws(sample_ids))) || anyDuplicated(sample_ids)) {
    stop(label, " requires unique sample identifiers")
  }
  if (!is.list(validation)) stop(label, " requires validation metadata")
  effective <- validation[["effective_batch_keys"]]
  if (is.null(effective)) {
    summary <- validation[["validation_summary"]]
    levels <- validation[["per_key_levels"]]
    if (is.null(levels)) levels <- validation[["levels"]]
    if (is.null(levels) && is.list(summary)) {
      levels <- summary[["per_key_levels"]]
    }
    if (!is.list(levels) || is.null(names(levels))) {
      stop(label, " is missing per-key levels")
    }
    effective <- keys[vapply(
      levels[keys],
      function(values) length(values) >= 2L,
      logical(1L)
    )]
  }
  effective <- unname(as.character(effective))
  if (length(effective) && (
    anyNA(effective) || any(!effective %in% keys) || anyDuplicated(effective) ||
    !identical(effective, keys[keys %in% effective])
  )) {
    stop(label, " has invalid effective batch keys")
  }
  non_estimable <- validation[["non_estimable_batch_keys"]]
  if (is.null(non_estimable)) non_estimable <- setdiff(keys, effective)
  non_estimable <- unname(as.character(non_estimable))
  if (!identical(non_estimable, unname(setdiff(keys, effective)))) {
    stop(label, " has invalid non-estimable batch keys")
  }

  model_data <- data.frame(
    row.names = sample_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  aliases <- if (length(effective)) {
    setNames(paste0("batch_key_", match(effective, keys)), effective)
  } else {
    character()
  }
  for (key in effective) {
    if (!key %in% colnames(metadata)) {
      stop(label, " metadata is missing key ", key)
    }
    values <- metadata[[key]]
    if (is.list(values) && !is.factor(values)) {
      stop(label, " key ", key, " is list-valued")
    }
    model_data[[unname(aliases[[key]])]] <- factor(as.character(values))
  }
  model_formula <- if (length(effective)) {
    stats::as.formula(paste0(
      "~ 1 + ", paste(unname(aliases), collapse = " + ")
    ))
  } else {
    stats::formula("~ 1")
  }
  design <- tryCatch(
    stats::model.matrix(model_formula, data = model_data),
    error = function(error) {
      stop(label, " could not be constructed: ", conditionMessage(error))
    }
  )
  if (!is.numeric(design) || any(!is.finite(design))) {
    stop(label, " contains non-finite design values")
  }
  rank <- qr(design, tol = 1e-10)$rank
  if (rank < ncol(design)) stop(label, " is rank deficient or confounded")
  residual_df <- nrow(design) - rank
  if (residual_df <= 0L) {
    stop(label, " is non-estimable: no residual degrees of freedom")
  }
  technical_columns <- which(colnames(design) != "(Intercept)")
  list(
    data = model_data,
    formula = model_formula,
    design = design,
    aliases = aliases,
    configured_batch_keys = unname(keys),
    effective_batch_keys = effective,
    non_estimable_batch_keys = non_estimable,
    technical_columns = as.integer(technical_columns),
    correction_state = if (length(effective)) {
      "BATCH_CORRECTION"
    } else {
      "NO_CORRECTION"
    },
    rank = as.integer(rank),
    columns = as.integer(ncol(design)),
    residual_df = as.integer(residual_df)
  )
}

# This is the only contract logic retained in R: fitting the separate
# technical fixed-effect design consumed by limma.
ecoda_batch_fixed_effect_design <- function(
  metadata,
  batch_keys,
  validation = NULL,
  sample_col = "Sample"
) {
  if (!is.data.frame(metadata)) {
    stop("limma batch metadata must be a data.frame")
  }
  keys <- ecoda_batch_normalize_keys(batch_keys, sample_col = sample_col)
  if (is.null(validation)) {
    validation <- ecoda_batch_validate_metadata(
      metadata = metadata,
      batch_keys = keys,
      sample_col = sample_col
    )
  }
  if (is.list(validation) && !is.null(validation[["validation_summary"]])) {
    summary <- validation[["validation_summary"]]
    validation[["validation_summary"]] <- ecoda_batch_build_validation_summary(
      validation,
      summary[["correction_mode"]],
      summary[["correction_formula"]]
    )
  }
  .ecoda_batch_design_from_validation(
    metadata = metadata,
    validation = validation,
    keys = keys,
    sample_col = sample_col
  )
}

.ecoda_batch_validation_from_python <- function(
  metadata,
  serialized,
  compact_identity,
  keys,
  sample_col,
  biological_label,
  accepted_sentinel_values = NULL
) {
  raw_sample_ids <- unname(as.character(metadata[[sample_col]]))
  sample_ids <- unname(as.character(serialized[["sample_ids"]]))
  first_indices <- match(sample_ids, raw_sample_ids)
  if (anyNA(first_indices)) {
    stop("Python batch validation returned unknown Sample identifiers")
  }
  sample_metadata <- data.frame(
    Sample = sample_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (key in keys) sample_metadata[[key]] <- metadata[[key]][first_indices]
  rownames(sample_metadata) <- sample_ids

  sample_group_ids <- unname(as.character(serialized[["sample_group_ids"]]))
  canonical_sample_metadata <- data.frame(
    Sample = sample_group_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  canonical_values <- setNames(vector("list", length(keys)), keys)
  for (key in keys) {
    values <- sample_metadata[[key]]
    canonical <- ecoda_batch_canonical_values(
      values,
      paste0(key, " samples"),
      if (is.null(accepted_sentinel_values)) {
        NULL
      } else {
        accepted_sentinel_values[[key]]
      }
    )
    canonical_values[[key]] <- canonical
    canonical_sample_metadata[[key]] <- canonical
  }
  rownames(canonical_sample_metadata) <- sample_group_ids

  levels <- serialized[["per_key_levels"]]
  if (!is.list(levels) || is.null(names(levels))) {
    stop("Python batch validation returned no per-key levels")
  }
  levels <- levels[keys]
  if (any(vapply(levels, is.null, logical(1L)))) {
    stop("Python batch validation omitted a configured key")
  }
  levels <- setNames(lapply(levels, function(values) {
    unname(as.character(unlist(values, use.names = FALSE)))
  }), keys)
  key_level_counts <- setNames(
    lapply(levels, function(values) as.integer(length(values))),
    keys
  )
  effective <- keys[vapply(levels, function(values) length(values) >= 2L, logical(1L))]
  non_estimable <- unname(setdiff(keys, effective))
  provisional <- list(
    valid = TRUE,
    ordered_keys = keys,
    sample_ids = sample_ids,
    canonical_sample_metadata = canonical_sample_metadata,
    per_key_levels = levels,
    effective_batch_keys = effective,
    non_estimable_batch_keys = non_estimable
  )
  design_info <- .ecoda_batch_design_from_validation(
    metadata = sample_metadata,
    validation = provisional,
    keys = keys,
    sample_col = "Sample",
    label = "Python-validated limma batch design"
  )
  summary <- compact_identity[["validation_summary"]]
  composite_levels <- unname(as.character(serialized[["composite_levels"]]))
  sample_composite_values <- unname(as.character(
    serialized[["sample_composite_values"]]
  ))
  list(
    valid = TRUE,
    token_version = as.character(serialized[["token_version"]]),
    encoding = as.character(serialized[["token_version"]]),
    encoding_version = as.character(serialized[["token_version"]]),
    sample_col = sample_col,
    sample_column = sample_col,
    biological_label = biological_label,
    biological_column = biological_label,
    ordered_keys = unname(keys),
    keys = unname(keys),
    key_count = as.integer(length(keys)),
    scalarization = as.character(serialized[["scalarization"]]),
    sample_ids = sample_ids,
    sample_group_ids = sample_group_ids,
    n_cells = as.integer(serialized[["n_obs"]]),
    n_obs = as.integer(serialized[["n_obs"]]),
    n_samples = as.integer(serialized[["n_samples"]]),
    sample_constancy = setNames(lapply(keys, function(key) TRUE), keys),
    sample_metadata = sample_metadata,
    canonical_sample_metadata = canonical_sample_metadata,
    canonical_values = canonical_values,
    canonical_cell_values = NULL,
    per_key_levels = levels,
    levels = levels,
    key_level_counts = key_level_counts,
    key_near_unique_fraction = setNames(
      lapply(levels, function(values) {
        as.numeric(length(values) / length(sample_ids))
      }),
      keys
    ),
    near_unique_fraction = as.numeric(serialized[["near_unique_fraction"]]),
    effective_batch_keys = effective,
    non_estimable_batch_keys = non_estimable,
    correction_state = design_info$correction_state,
    fixed_effect_aliases = design_info$aliases,
    correction_design_formula = paste(
      deparse(design_info$formula), collapse = ""
    ),
    effective_design = list(
      rank = design_info$rank,
      columns = design_info$columns,
      residual_df = design_info$residual_df
    ),
    effective_design_rank = design_info$rank,
    effective_design_columns = design_info$columns,
    effective_design_residual_df = design_info$residual_df,
    composite_values = sample_composite_values,
    composite_levels = composite_levels,
    composite_level_count = as.integer(length(composite_levels)),
    composite_near_unique_fraction = if (length(keys) >= 2L) {
      as.numeric(length(composite_levels) / length(sample_ids))
    } else {
      NA_real_
    },
    estimable = TRUE,
    additive_design = list(
      rank = design_info$rank,
      columns = design_info$columns,
      residual_df = design_info$residual_df
    ),
    design_rank = design_info$rank,
    design_columns = design_info$columns,
    design_residual_df = design_info$residual_df,
    composite_design = list(
      rank = as.integer(serialized[["composite_design_rank"]]),
      columns = as.integer(serialized[["composite_design_columns"]])
    ),
    composite_design_rank = as.integer(serialized[["composite_design_rank"]]),
    composite_design_columns = as.integer(serialized[["composite_design_columns"]]),
    validation_summary = summary,
    python_metadata = serialized,
    compact_identity = compact_identity
  )
}

# Validate full metadata in Python, then expose an R-only model view. The
# persisted identity is always the compact Python serialization, never this
# in-memory compatibility view.
ecoda_batch_validate_metadata <- function(
  metadata,
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50,
  accepted_sentinel_values = NULL,
  enforce_near_unique = TRUE,
  enforce_composite_near_unique = FALSE
) {
  if (!is.data.frame(metadata)) {
    stop("cell metadata must be a data.frame")
  }
  keys <- ecoda_batch_normalize_keys(
    batch_keys,
    sample_col = sample_col,
    biological_label = biological_label
  )
  module <- .ecoda_batch_python_module()
  py_validation <- module$validate_batch_metadata(
    reticulate::r_to_py(metadata),
    reticulate::r_to_py(as.list(keys)),
    sample_column = sample_col,
    biological_column = biological_label,
    near_unique_fraction = as.numeric(near_unique_fraction),
    accepted_sentinel_values = .ecoda_batch_accepted_sentinels_python(
      accepted_sentinel_values
    ),
    enforce_near_unique = enforce_near_unique,
    enforce_composite_near_unique = enforce_composite_near_unique
  )
  compact_identity <- .ecoda_batch_to_r(module$serialize_batch_contract_identity(
    py_validation,
    method_id = "Pseudobulk",
    model_id = "pseudobulk_limma_fixed_effects_v1"
  ))
  serialized <- list(
    token_version = as.character(.ecoda_batch_to_r(py_validation$token_version)),
    scalarization = as.character(.ecoda_batch_to_r(py_validation$scalarization)),
    sample_ids = as.character(.ecoda_batch_to_r(py_validation$sample_ids)),
    sample_group_ids = as.character(
      .ecoda_batch_to_r(py_validation$sample_group_ids)
    ),
    n_obs = as.integer(.ecoda_batch_to_r(py_validation$n_obs)),
    n_samples = as.integer(.ecoda_batch_to_r(py_validation$n_samples)),
    per_key_levels = .ecoda_batch_to_r(py_validation$levels),
    composite_levels = as.character(
      .ecoda_batch_to_r(py_validation$composite_levels)
    ),
    sample_composite_values = as.character(
      .ecoda_batch_to_r(py_validation$sample_composite_values)
    ),
    near_unique_fraction = as.numeric(
      .ecoda_batch_to_r(py_validation$near_unique_fraction)
    ),
    composite_design_rank = as.integer(
      .ecoda_batch_to_r(py_validation$composite_design_rank)
    ),
    composite_design_columns = as.integer(
      .ecoda_batch_to_r(py_validation$composite_design_columns)
    )
  )
  .ecoda_batch_validation_from_python(
    metadata = metadata,
    serialized = serialized,
    compact_identity = compact_identity,
    keys = keys,
    sample_col = sample_col,
    biological_label = biological_label,
    accepted_sentinel_values = accepted_sentinel_values
  )
}

# These wrappers keep existing callers on the Python contract while exposing
# the exact identity payload for cross-language parity tests.
ecoda_batch_fingerprint <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  module <- .ecoda_batch_python_module()
  value <- module$batch_contract_fingerprint(
    reticulate::r_to_py(.ecoda_batch_python_keys(batch_keys)),
    scalarization = scalarization,
    method_id = method_id,
    model_id = model_id
  )
  unname(as.character(.ecoda_batch_to_r(value)))
}

ecoda_batch_fingerprint_payload_hex <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  module <- .ecoda_batch_python_module()
  value <- module$batch_contract_fingerprint_payload_hex(
    reticulate::r_to_py(.ecoda_batch_python_keys(batch_keys)),
    scalarization = scalarization,
    method_id = method_id,
    model_id = model_id
  )
  unname(as.character(.ecoda_batch_to_r(value)))
}

ecoda_batch_fingerprint_payload <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  hex <- ecoda_batch_fingerprint_payload_hex(
    batch_keys, method_id, model_id, scalarization
  )
  if (!length(hex) || is.na(hex) || nchar(hex) %% 2L != 0L) {
    stop("Python batch fingerprint payload is malformed")
  }
  bytes <- vapply(seq(1L, nchar(hex), by = 2L), function(index) {
    value <- suppressWarnings(strtoi(substr(hex, index, index + 1L), base = 16L))
    if (is.na(value)) stop("Python batch fingerprint payload is not hexadecimal")
    value
  }, integer(1L))
  rawToChar(as.raw(bytes[bytes != 0L]))
}

ecoda_batch_contract_identity <- function(
  batch_keys,
  sample_col = "Sample",
  method_id = NULL,
  model_id = NULL
) {
  module <- .ecoda_batch_python_module()
  result <- module$build_batch_contract_identity(
    reticulate::r_to_py(.ecoda_batch_python_keys(batch_keys)),
    sample_column = sample_col,
    method_id = method_id,
    model_id = model_id
  )
  .ecoda_batch_to_r(result)
}

# Normalize and validate the compact Python summary before it crosses into R.
ecoda_batch_build_validation_summary <- function(
  validation,
  correction_mode,
  correction_formula
) {
  if (!is.list(validation)) stop("validated batch metadata must be a named list")
  summary <- validation[["validation_summary"]]
  keys <- .ecoda_batch_validation_keys(validation)
  if (is.null(summary)) {
    levels <- validation[["per_key_levels"]]
    if (is.null(levels)) levels <- validation[["levels"]]
    if (!is.list(levels) || is.null(names(levels))) {
      stop("validated batch metadata has no compact level summary")
    }
    levels <- levels[keys]
    summary <- list(
      schema_version = 1L,
      validated_before_reduction = TRUE,
      sample_constancy = setNames(lapply(keys, function(key) TRUE), keys),
      per_key_levels = levels,
      key_level_counts = setNames(
        lapply(levels, function(values) as.integer(length(values))), keys
      ),
      composite_levels = unname(as.character(
        validation[["composite_levels"]] %||% character()
      )),
      composite_level_count = as.integer(length(
        validation[["composite_levels"]] %||% character()
      )),
      n_cells = as.integer(validation[["n_cells"]] %||% validation[["n_obs"]]),
      n_samples = as.integer(validation[["n_samples"]]),
      correction_mode = correction_mode,
      correction_formula = correction_formula
    )
  } else {
    summary <- summary
    summary[["correction_mode"]] <- correction_mode
    summary[["correction_formula"]] <- correction_formula
  }
  module <- .ecoda_batch_python_module()
  normalized <- module$validate_batch_validation_summary(
    .ecoda_batch_summary_python(summary),
    reticulate::r_to_py(as.list(keys))
  )
  .ecoda_batch_to_r(normalized)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

ecoda_batch_augment_contract <- function(
  identity,
  validation,
  correction_mode,
  correction_formula
) {
  summary <- if (
    is.list(validation) && !is.null(validation[["validation_summary"]])
  ) {
    validation[["validation_summary"]]
  } else if (
    is.list(validation) && all(c(
      "schema_version", "validated_before_reduction", "sample_constancy",
      "per_key_levels", "key_level_counts", "composite_levels",
      "composite_level_count", "n_cells", "n_samples", "correction_mode",
      "correction_formula"
    ) %in% names(validation))
  ) {
    validation
  } else {
    ecoda_batch_build_validation_summary(
      validation, correction_mode, correction_formula
    )
  }
  summary[["correction_mode"]] <- correction_mode
  summary[["correction_formula"]] <- correction_formula
  module <- .ecoda_batch_python_module()
  result <- module$augment_batch_contract(
    .ecoda_batch_identity_python(identity),
    .ecoda_batch_summary_python(summary),
    correction_mode,
    correction_formula
  )
  .ecoda_batch_to_r(result)
}

ecoda_batch_correction_spec <- function(
  method_id,
  batch_keys,
  scalar_batch_col = NULL,
  effective_batch_keys = NULL,
  non_estimable_batch_keys = NULL
) {
  keys <- ecoda_batch_normalize_keys(batch_keys)
  if (!is.null(scalar_batch_col) && length(keys) == 1L &&
      !identical(as.character(scalar_batch_col), keys[[1L]])) {
    stop("scalar batch column differs from the configured key")
  }
  effective <- if (is.null(effective_batch_keys)) {
    unname(keys)
  } else {
    unname(as.character(effective_batch_keys))
  }
  non_estimable <- if (is.null(non_estimable_batch_keys)) {
    unname(setdiff(keys, effective))
  } else {
    unname(as.character(non_estimable_batch_keys))
  }
  module <- .ecoda_batch_python_module()
  result <- .ecoda_batch_to_r(module$batch_correction_spec_for_keys(
    method_id,
    reticulate::r_to_py(as.list(keys)),
    effective_batch_keys = reticulate::r_to_py(as.list(effective)),
    non_estimable_batch_keys = reticulate::r_to_py(as.list(non_estimable))
  ))
  aliases <- if (length(effective)) {
    setNames(paste0("batch_key_", match(effective, keys)), effective)
  } else {
    character()
  }
  list(
    correction_mode = as.character(result[[1L]]),
    correction_formula = as.character(result[[2L]]),
    correction_state = if (length(effective)) {
      "BATCH_CORRECTION"
    } else {
      "NO_CORRECTION"
    },
    effective_batch_keys = effective,
    non_estimable_batch_keys = non_estimable,
    aliases = aliases
  )
}

# Build a compact identity from Python-validated metadata. No vector-bearing
# field from the in-memory compatibility view is persisted.
ecoda_batch_build_composite <- function(
  metadata,
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50
) {
  validation <- ecoda_batch_validate_metadata(
    metadata, batch_keys, sample_col, biological_label, near_unique_fraction
  )
  keys <- validation$ordered_keys
  list(
    token_version = validation$token_version,
    encoding = validation$encoding,
    ordered_keys = keys,
    key_count = validation$key_count,
    scalarization = validation$scalarization,
    composite_name = if (length(keys) >= 2L) {
      as.character(validation$compact_identity$reserved_obs_name)
    } else {
      keys[[1L]]
    },
    reserved_name = as.character(validation$compact_identity$reserved_obs_name),
    sample_ids = validation$sample_ids,
    effective_batch_keys = validation$effective_batch_keys,
    non_estimable_batch_keys = validation$non_estimable_batch_keys,
    correction_state = validation$correction_state,
    fixed_effect_aliases = validation$fixed_effect_aliases,
    correction_design_formula = validation$correction_design_formula,
    composite_values = validation$composite_values,
    composite_levels = validation$composite_levels,
    composite_level_count = validation$composite_level_count,
    per_key_levels = validation$per_key_levels,
    levels = validation$per_key_levels,
    sample_metadata = validation$sample_metadata,
    canonical_sample_metadata = validation$canonical_sample_metadata,
    validation = validation
  )
}

ecoda_batch_serialize_metadata <- function(
  metadata,
  batch_keys,
  method_id,
  model_id,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50
) {
  if (!is.data.frame(metadata)) stop("cell metadata must be a data.frame")
  keys <- ecoda_batch_normalize_keys(
    batch_keys, sample_col = sample_col, biological_label = biological_label
  )
  module <- .ecoda_batch_python_module()
  validation <- module$validate_batch_metadata(
    reticulate::r_to_py(metadata),
    reticulate::r_to_py(as.list(keys)),
    sample_column = sample_col,
    biological_column = biological_label,
    near_unique_fraction = as.numeric(near_unique_fraction)
  )
  serialized_py <- module$serialize_batch_metadata(
    validation,
    method_id = method_id,
    model_id = model_id
  )
  serialized_py$pop("fingerprint_payload")
  serialized <- .ecoda_batch_to_r(serialized_py)
  serialized[["fingerprint_payload"]] <- ecoda_batch_fingerprint_payload(
    keys,
    method_id = method_id,
    model_id = model_id,
    scalarization = serialized[["scalarization"]]
  )
  serialized
}

# Explicit aliases retained for the small number of legacy in-memory callers.
ecoda_batch_validation_summary <- ecoda_batch_build_validation_summary
ecoda_batch_augment_contract_identity <- ecoda_batch_augment_contract

if (exists(".ecoda_batch_source_file", inherits = FALSE)) {
  rm(.ecoda_batch_source_file)
}
