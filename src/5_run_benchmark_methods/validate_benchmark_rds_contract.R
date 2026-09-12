#!/usr/bin/env Rscript
# Validate selected result bundles, matrix identifiers, and atomic artifacts.
args <- commandArgs(trailingOnly = TRUE)
value_for <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing ", flag)
  args[[i + 1L]]
}
has_flag <- function(flag) flag %in% args
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
validator_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
} else {
  getwd()
}
artifact_path <- value_for("--artifact", "")
artifact_list <- value_for("--artifact-list", "")
root <- value_for("--root")
selection <- value_for("--selection")
labels_arg <- value_for("--labels", "")
labels <- if (nzchar(labels_arg)) strsplit(labels_arg, ",", fixed = TRUE)[[1L]] else character()
batch_pass <- value_for("--batch-pass", "")
config_path <- value_for("--config", Sys.getenv("DATASETS_JSON_FILE", unset = ""))
input_root <- value_for("--input-root", "")
dataset_arg <- value_for("--dataset", "")
view_arg <- value_for("--view", "")
method_arg <- value_for("--method", "")
metadata_kind <- has_flag("--metadata")
source_identity_path <- value_for("--source-identity", "")
expected_batch_contract_arg <- value_for("--expected-batch-contract", "")
batch_contract_arg <- value_for("--batch-contract", "")
source_identity_verified <- has_flag("--source-identity-verified")
exact <- has_flag("--exact")
batch <- nzchar(batch_pass)
supported_labels <- c("gloscope", "mofa", "pseudobulk", "composition", "scitd",
                      "prepare_pseudobulk", "trans", "zeroimp")
if (nzchar(artifact_path) && nzchar(artifact_list)) {
  stop("--artifact and --artifact-list are mutually exclusive")
}
if (!nzchar(artifact_path) && !nzchar(artifact_list) &&
    any(!labels %in% supported_labels)) {
  stop("unsupported selected RDS label: ",
       paste(labels[!labels %in% supported_labels], collapse = ", "))
}
if (!nzchar(artifact_path) && !nzchar(artifact_list) &&
    (is.null(root) || is.null(selection) || !length(labels) || any(!nzchar(labels)))) {
  stop("--root, --selection, and --labels are required")
}

checksum_ok <- function(file) {
  sidecar <- paste0(file, ".md5")
  info <- if (file.exists(file)) file.info(file) else NULL
  if (is.null(info) || !isTRUE(info$isdir == FALSE) ||
      is.na(info$size) || info$size <= 0 || !file.exists(sidecar)) {
    return(FALSE)
  }
  lines <- tryCatch(readLines(sidecar, warn = FALSE), error = function(e) NULL)
  keys <- c("MD5", "SIZE", "PATH")
  if (is.null(lines) || length(lines) != length(keys) ||
      any(!startsWith(lines, paste0(keys, "=")))) {
    return(FALSE)
  }
  values <- substring(lines, nchar(keys) + 2L)
  digest <- values[[1L]]
  size <- values[[2L]]
  recorded_path <- values[[3L]]
  if (!grepl("^[0-9a-f]{32}$", digest, perl = TRUE) ||
      !grepl("^[1-9][0-9]*$", size, perl = TRUE) ||
      !identical(recorded_path, file) ||
      !identical(size, as.character(info$size))) {
    return(FALSE)
  }
  actual <- tryCatch(unname(tools::md5sum(file)), error = function(e) NA_character_)
  isTRUE(!is.na(actual) && identical(actual, digest))
}

partial_name_patterns <- function(path) {
  bases <- c(path, paste0(path, ".md5"))
  hidden <- file.path(dirname(bases), paste0(".", basename(bases)))
  bases <- unique(c(bases, hidden))
  suffixes <- c(
    ".tmp", ".tmp.*",
    ".build", ".build.*",
    ".partial", ".partial.*"
  )
  unlist(lapply(bases, paste0, suffixes), use.names = FALSE)
}

selected_partial_paths <- function(paths) {
  patterns <- unlist(
    lapply(unique(as.character(paths)), partial_name_patterns),
    use.names = FALSE
  )
  matches <- unique(unlist(lapply(patterns, Sys.glob), use.names = FALSE))
  matches[file.exists(matches)]
}

reject_selected_partials <- function(paths) {
  partials <- selected_partial_paths(paths)
  if (length(partials)) {
    stop("partial benchmark artifacts remain: ",
         paste(partials, collapse = ", "))
  }
}



artifact_list_parts <- list()
if (nzchar(artifact_path)) {
  selection_rows <- character()
  parts <- list()
} else if (nzchar(artifact_list)) {
  if (!checksum_ok(artifact_list)) {
    stop("artifact-list checksum is missing or invalid: ", artifact_list)
  }
  artifact_list_rows <- readLines(artifact_list, warn = FALSE)
  if (!length(artifact_list_rows) || any(!nzchar(artifact_list_rows))) {
    stop("artifact-list is empty or contains a blank row: ", artifact_list)
  }
  artifact_list_parts <- strsplit(artifact_list_rows, "\t", fixed = TRUE)
  if (any(vapply(artifact_list_parts, length, integer(1L)) != 5L) ||
      any(vapply(artifact_list_parts, function(x) any(!nzchar(x)), logical(1L))) ||
      any(vapply(artifact_list_parts, function(x) !x[[5L]] %in% c("0", "1"), logical(1L)))) {
    stop("artifact-list rows must be PATH<TAB>METHOD<TAB>DATASET<TAB>VIEW<TAB>METADATA")
  }
  parts <- list()
} else {
  if (!checksum_ok(selection)) stop("selection checksum is missing or invalid: ", selection)
  selection_rows <- readLines(selection, warn = FALSE)
  if (!length(selection_rows) || any(!nzchar(selection_rows))) {
    stop("selection is empty or contains a blank row: ", selection)
  }
  parts <- strsplit(selection_rows, "\t", fixed = TRUE)
  if (any(vapply(parts, length, integer(1L)) != 3L) ||
      any(vapply(parts, function(x) any(!nzchar(x)), logical(1L)))) {
    stop("selection rows must have three non-empty columns")
  }
  row_keys <- vapply(parts, paste, character(1L), collapse = "\t")
  if (anyDuplicated(row_keys)) stop("selection contains duplicate rows")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
.batch_identity_alias <- function(identity, fields, label) {
  if (!is.list(identity) || is.null(names(identity))) {
    stop("corrected batch contract identity for ", label, " must be a named list")
  }
  present <- fields[fields %in% names(identity)]
  if (!length(present)) {
    stop("corrected batch contract identity is missing ", label)
  }
  value <- identity[[present[[1L]]]]
  if (length(present) > 1L &&
      any(!vapply(present[-1L], function(field) {
        identical(identity[[field]], value)
      }, logical(1L)))) {
    stop("corrected batch contract identity has mismatched ", label, " aliases")
  }
  value
}

.batch_identity_string_list <- function(value, label, nonempty = TRUE) {
  if (is.character(value)) {
    values <- unname(value)
  } else if (is.list(value)) {
    if (length(value)) {
      values <- vapply(seq_along(value), function(index) {
        item <- value[[index]]
        if (!is.character(item) || length(item) != 1L || is.na(item)) {
          stop(
            "corrected batch contract identity ", label,
            " must be an ordered list of strings"
          )
        }
        unname(item)
      }, character(1L))
    } else {
      values <- character()
    }
  } else {
    stop(
      "corrected batch contract identity ", label,
      " must be an ordered list of strings"
    )
  }
  if (nonempty && !length(values)) {
    stop("corrected batch contract identity ", label, " must be nonempty")
  }
  if (anyNA(values) || any(!nzchar(values)) ||
      any(values != trimws(values)) || anyDuplicated(values)) {
    stop(
      "corrected batch contract identity ", label,
      " contains a blank, padded, or duplicate value"
    )
  }
  unname(values)
}

.batch_identity_scalar <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop("corrected batch contract identity ", label, " must be one nonblank string")
  }
  unname(value)
}

.batch_identity_vector_fields <- c(
  "composite_values", "sample_composite_values", "sample_ids",
  "sample_group_ids", "canonical_values", "canonical_cell_values",
  "canonical_sample_metadata", "sample_metadata", "row_tokens", "tokens",
  "scalarized_values", "sample_scalarized_values"
)
.batch_validation_summary_fields <- c(
  "schema_version", "validated_before_reduction", "sample_constancy",
  "per_key_levels", "key_level_counts", "composite_levels",
  "composite_level_count", "n_cells", "n_samples", "correction_mode",
  "correction_formula"
)

.batch_summary_named_map <- function(value, label) {
  if (!is.list(value) || is.data.frame(value) || is.null(names(value)) ||
      anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    stop(label, " must be an ordered named map")
  }
  value
}

.batch_summary_integer <- function(value, label, minimum = NULL) {
  if ((!is.integer(value) && !is.numeric(value)) ||
      length(value) != 1L || is.na(value) ||
      !is.finite(value) || floor(value) != value ||
      (!is.null(minimum) && value < minimum)) {
    stop(label, " must be an integer scalar")
  }
  unname(as.integer(value))
}

.batch_summary_string_list <- function(value, label) {
  if (is.character(value)) {
    values <- unname(value)
  } else if (is.list(value)) {
    values <- if (length(value)) {
      vapply(seq_along(value), function(index) {
        item <- value[[index]]
        if (!is.character(item) || length(item) != 1L || is.na(item)) {
          stop(label, " must contain scalar strings")
        }
        unname(item)
      }, character(1L))
    } else {
      character()
    }
  } else {
    stop(label, " must be an ordered list of strings")
  }
  if (anyNA(values) || any(!vapply(values, validUTF8, logical(1L)))) {
    stop(label, " must contain valid strings")
  }
  unname(values)
}

.batch_summary_raw_sort <- function(values, label) {
  if (!length(values)) return(character())
  if (!is.character(values) || anyNA(values) ||
      any(!vapply(values, validUTF8, logical(1L)))) {
    stop(label, " must contain valid strings")
  }
  raw_hex <- vapply(unname(values), function(value) {
    paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))), collapse = "")
  }, character(1L))
  unname(values)[order(raw_hex, method = "radix")]
}

validate_batch_validation_summary <- function(
  summary,
  batch_keys = NULL,
  label = "validation_summary"
) {
  if (!is.list(summary) || is.null(names(summary)) ||
      length(names(summary)) != length(.batch_validation_summary_fields) ||
      anyNA(names(summary)) || anyDuplicated(names(summary)) ||
      !setequal(names(summary), .batch_validation_summary_fields)) {
    stop(label, " has an invalid field set")
  }
  schema_version <- .batch_summary_integer(
    summary[["schema_version"]], paste0(label, "$schema_version")
  )
  if (!identical(schema_version, 1L)) {
    stop(label, "$schema_version must be 1")
  }
  if (!is.logical(summary[["validated_before_reduction"]]) ||
      length(summary[["validated_before_reduction"]]) != 1L ||
      is.na(summary[["validated_before_reduction"]]) ||
      !isTRUE(summary[["validated_before_reduction"]])) {
    stop(label, "$validated_before_reduction must be TRUE")
  }

  constancy <- .batch_summary_named_map(
    summary[["sample_constancy"]], paste0(label, "$sample_constancy")
  )
  levels <- .batch_summary_named_map(
    summary[["per_key_levels"]], paste0(label, "$per_key_levels")
  )
  counts <- .batch_summary_named_map(
    summary[["key_level_counts"]], paste0(label, "$key_level_counts")
  )
  constancy_keys <- names(constancy)
  level_keys <- names(levels)
  count_keys <- names(counts)
  if (!setequal(constancy_keys, level_keys) ||
      !setequal(level_keys, count_keys)) {
    stop(label, " key maps must contain the same configured keys")
  }

  if (!is.null(batch_keys)) {
    .batch_identity_load_contract()
    expected_keys <- tryCatch(
      ecoda_batch_normalize_keys(
        if (is.character(batch_keys) && length(batch_keys) > 1L) {
          as.list(batch_keys)
        } else {
          batch_keys
        }
      ),
      error = function(error) {
        stop(label, " keys do not match configured batch keys: ",
             conditionMessage(error))
      }
    )
    if (!setequal(level_keys, expected_keys)) {
      stop(label, " keys do not match configured batch keys")
    }
    keys <- unname(expected_keys)
  } else {
    if (!identical(level_keys, constancy_keys) ||
        !identical(level_keys, count_keys)) {
      stop(label, " key maps must use one configured order")
    }
    keys <- level_keys
  }
  if (!length(keys) || anyNA(keys) || any(!nzchar(keys)) ||
      any(keys != trimws(keys))) {
    stop(label, " contains invalid batch keys")
  }

  normalized_levels <- setNames(vector("list", length(keys)), keys)
  for (key in keys) {
    constant <- constancy[[key]]
    if (!is.logical(constant) || length(constant) != 1L ||
        is.na(constant) || !isTRUE(constant)) {
      stop(label, "$sample_constancy is false for ", key)
    }
    key_levels <- .batch_summary_string_list(
      levels[[key]], paste0(label, "$per_key_levels$", key)
    )
    if (anyDuplicated(key_levels) || length(key_levels) < 2L ||
        !identical(
          unname(key_levels),
          .batch_summary_raw_sort(key_levels, paste0(key, " levels"))
        )) {
      stop(label, " has invalid levels for ", key)
    }
    normalized_levels[[key]] <- key_levels
    level_count <- .batch_summary_integer(
      counts[[key]], paste0(label, "$key_level_counts$", key)
    )
    if (!identical(level_count, as.integer(length(key_levels)))) {
      stop(label, " has an invalid level count for ", key)
    }
  }

  composite_levels <- .batch_summary_string_list(
    summary[["composite_levels"]], paste0(label, "$composite_levels")
  )
  if (anyDuplicated(composite_levels) ||
      !identical(
        unname(composite_levels),
        .batch_summary_raw_sort(composite_levels, paste0(label, "$composite_levels"))
      )) {
    stop(label, " has invalid composite_levels")
  }

  composite_count <- .batch_summary_integer(
    summary[["composite_level_count"]],
    paste0(label, "$composite_level_count")
  )
  if (!identical(composite_count, as.integer(length(composite_levels)))) {
    stop(label, " has an invalid composite_level_count")
  }
  if (length(keys) == 1L && length(composite_levels)) {
    stop(label, " direct validation cannot carry composite levels")
  }
  if (length(keys) >= 2L && !length(composite_levels)) {
    stop(label, " composite validation requires composite levels")
  }

  n_cells <- .batch_summary_integer(
    summary[["n_cells"]], paste0(label, "$n_cells"), minimum = 1L
  )
  n_samples <- .batch_summary_integer(
    summary[["n_samples"]], paste0(label, "$n_samples"), minimum = 2L
  )
  for (field in c("correction_mode", "correction_formula")) {
    value <- summary[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value) || !identical(value, trimws(value)) ||
        !isTRUE(validUTF8(value))) {
      stop(label, " has invalid ", field)
    }
  }

  list(
    schema_version = schema_version,
    validated_before_reduction = TRUE,
    sample_constancy = setNames(
      lapply(keys, function(key) TRUE), keys
    ),
    per_key_levels = setNames(
      lapply(keys, function(key) normalized_levels[[key]]), keys
    ),
    key_level_counts = setNames(
      lapply(keys, function(key) .batch_summary_integer(
        counts[[key]], paste0(label, "$key_level_counts$", key)
      )), keys
    ),
    composite_levels = unname(composite_levels),
    composite_level_count = composite_count,
    n_cells = n_cells,
    n_samples = n_samples,
    correction_mode = unname(summary[["correction_mode"]]),
    correction_formula = unname(summary[["correction_formula"]])
  )
}

.batch_identity_summary <- function(
  identity,
  keys,
  label,
  method_id = NULL,
  required = FALSE,
  validate_optional = TRUE
) {
  if (!is.list(identity) || is.null(names(identity))) {
    stop("corrected batch contract identity for ", label, " must be a named list")
  }
  has_summary <- "validation_summary" %in% names(identity)
  vector_fields <- sort(intersect(.batch_identity_vector_fields, names(identity)))
  if ((isTRUE(required) || (isTRUE(validate_optional) && has_summary)) &&
      length(vector_fields)) {
    stop(
      label, " contains forbidden per-cell/sample vector fields: ",
      paste(vector_fields, collapse = ", ")
    )
  }
  if (!has_summary) {
    if (isTRUE(required)) stop(label, " is missing validation_summary")
    return(NULL)
  }
  if (!isTRUE(required) && !isTRUE(validate_optional)) return(NULL)
  normalized <- tryCatch(
    validate_batch_validation_summary(
      identity[["validation_summary"]],
      as.list(unname(keys)),
      paste0(label, "$validation_summary")
    ),
    error = function(error) {
      stop(label, " has an invalid validation_summary: ",
           conditionMessage(error))
    }
  )
  if (!is.null(method_id)) {
    .batch_identity_load_contract()
    spec <- tryCatch(
      ecoda_batch_correction_spec(
        method_id = method_id,
        batch_keys = as.list(unname(keys))
      ),
      error = function(error) {
        stop(label, " has an invalid validation_summary: ",
             conditionMessage(error))
      }
    )
    if (!identical(normalized[["correction_mode"]], spec[["correction_mode"]]) ||
        !identical(
          normalized[["correction_formula"]],
          spec[["correction_formula"]]
        )) {
      stop(
        label, " has the wrong correction policy for method ", method_id
      )
    }
  }
  normalized
}
.batch_identity_load_contract <- function() {
  if (exists("ecoda_batch_fingerprint", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }
  candidates <- unique(c(
    file.path(validator_dir, "..", "utils", "batch_contract.R"),
    file.path("src", "utils", "batch_contract.R")
  ))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("corrected batch contract fingerprint implementation is unavailable")
  }
  sys.source(candidates[[1L]], envir = .GlobalEnv)
  if (!exists("ecoda_batch_fingerprint", mode = "function", inherits = TRUE)) {
    stop("corrected batch contract fingerprint implementation is unavailable")
  }
  invisible(TRUE)
}

.batch_identity_normalize <- function(identity, label = "batch contract") {
  if (!is.list(identity) || is.null(names(identity))) {
    stop(
      "corrected batch contract identity for ", label,
      " must be a named list"
    )
  }
  keys <- .batch_identity_string_list(
    .batch_identity_alias(
      identity,
      c(
        "ordered_source_keys", "ordered_keys", "source_keys",
        "batch_keys", "keys"
      ),
      "ordered source keys"
    ),
    "ordered source keys"
  )
  reserved_name <- "__ecoda_batch_combined_v1"
  if (reserved_name %in% keys) {
    stop(
      "corrected batch contract identity ordered source keys contain reserved column ",
      reserved_name
    )
  }

  token_values <- .batch_identity_alias(
    identity,
    c("token_version", "encoding_version", "encoding"),
    "token/encoding version"
  )
  token_values <- .batch_identity_scalar(token_values, "token/encoding version")
  if (!identical(token_values, "ecoda_batch_composite_v1")) {
    stop(
      "corrected batch contract identity token/encoding version must be ",
      "ecoda_batch_composite_v1"
    )
  }

  scalarization <- .batch_identity_scalar(
    identity[["scalarization"]], "scalarization"
  )
  expected_scalarization <- if (length(keys) >= 2L) "composite_v1" else "direct_v1"
  if (!identical(scalarization, expected_scalarization)) {
    stop(
      "corrected batch contract identity scalarization does not match ",
      "the ordered source key count"
    )
  }

  method_id <- .batch_identity_scalar(
    .batch_identity_alias(
      identity, c("method_id", "method", "method_policy"), "method policy"
    ),
    "method policy"
  )
  model_id <- .batch_identity_scalar(
    .batch_identity_alias(
      identity, c("model_id", "model", "model_policy"), "model policy"
    ),
    "model policy"
  )
  method_ids <- c(
    "preprocess", "ECODA_authors_HR", "ECODA_seuratres_2",
    "ECODA_authors_HR_NULL", "Pseudobulk", "GloScope", "PILOT",
    "MrVI", "QOT"
  )
  model_ids <- c(
    "hvg_composite_v1", "harmony_native_list_v1",
    "ecoda_additive_random_intercepts_v1", "pseudobulk_composite_v1",
    "mrvi_composite_v1", "embedding_consumer_harmony_v1"
  )
  if (!method_id %in% method_ids) {
    stop("corrected batch contract identity has unsupported method policy: ", method_id)
  }
  if (!model_id %in% model_ids) {
    stop("corrected batch contract identity has unsupported model policy: ", model_id)
  }
  if ("contract_version" %in% names(identity) &&
      !identical(identity[["contract_version"]], "ecoda_batch_contract_v1")) {
    stop("corrected batch contract identity has an unsupported contract version")
  }

  required_obs_columns <- .batch_identity_string_list(
    .batch_identity_alias(
      identity,
      c(
        "required_source_obs_columns", "required_obs_columns",
        "source_obs_columns", "obs_columns"
      ),
      "required source obs columns"
    ),
    "required source obs columns"
  )
  if (reserved_name %in% required_obs_columns) {
    stop(
      "corrected batch contract identity required source obs columns contain ",
      "reserved column ", reserved_name
    )
  }
  missing_keys <- keys[!keys %in% required_obs_columns]
  if (length(missing_keys)) {
    stop(
      "corrected batch contract identity required source obs columns omit ",
      "ordered source keys: ", paste(missing_keys, collapse = ", ")
    )
  }

  reserved_absent <- .batch_identity_alias(
    identity,
    c("reserved_obs_absent", "reserved_column_absent", "reserved_absent"),
    "reserved-column absence"
  )
  if (!is.logical(reserved_absent) || length(reserved_absent) != 1L ||
      is.na(reserved_absent) || !isTRUE(reserved_absent)) {
    stop(
      "corrected batch contract identity must assert absence of ",
      reserved_name
    )
  }
  reserved_name_fields <- c("reserved_obs_name", "reserved_name")
  present_reserved_name_fields <- reserved_name_fields[
    reserved_name_fields %in% names(identity)
  ]
  if (length(present_reserved_name_fields)) {
    recorded_name <- .batch_identity_alias(
      identity, reserved_name_fields, "reserved-column name"
    )
    if (!identical(recorded_name, reserved_name)) {
      stop(
        "corrected batch contract identity has an unsupported ",
        "reserved-column name"
      )
    }
  }

  fingerprint <- .batch_identity_scalar(
    .batch_identity_alias(
      identity,
      c("fingerprint", "batch_contract_fingerprint", "key_set_fingerprint"),
      "fingerprint"
    ),
    "fingerprint"
  )
  if (!grepl("^[0-9a-f]{64}$", fingerprint, perl = TRUE)) {
    stop("corrected batch contract identity fingerprint is malformed")
  }
  .batch_identity_load_contract()
  expected_fingerprint <- ecoda_batch_fingerprint(
    batch_keys = as.list(keys),
    method_id = method_id,
    model_id = model_id,
    scalarization = scalarization
  )
  if (!identical(fingerprint, expected_fingerprint)) {
    stop(
      "corrected batch contract identity fingerprint does not match its ",
      "source/configuration fields"
    )
  }
  normalized <- list(
    ordered_source_keys = keys,
    token_version = token_values,
    scalarization = scalarization,
    method_id = method_id,
    model_id = model_id,
    fingerprint = fingerprint,
    required_source_obs_columns = required_obs_columns,
    reserved_obs_absent = TRUE,
    contract_version = "ecoda_batch_contract_v1"
  )
  if ("validation_summary" %in% names(identity)) {
    .batch_identity_summary(
      identity,
      keys,
      label,
      method_id = method_id
    )
  }
  normalized
}

validate_batch_contract_identity <- function(
  expected_batch_contract = NULL,
  batch_contract = NULL,
  source_obs_columns = NULL,
  reserved_absent = NULL,
  require_recorded = FALSE,
  label = "batch contract",
  require_summary = NULL
) {
  if (is.null(expected_batch_contract) && is.null(batch_contract)) {
    return(NULL)
  }
  expected_source <- if (!is.null(expected_batch_contract)) {
    expected_batch_contract
  } else {
    batch_contract
  }
  expected <- .batch_identity_normalize(
    expected_source, paste0(label, " (expected)")
  )
  summary_required <- if (is.null(require_summary)) {
    isTRUE(require_recorded)
  } else {
    isTRUE(require_summary)
  }
  expected_summary <- .batch_identity_summary(
    expected_source,
    expected$ordered_source_keys,
    paste0(label, " expected identity"),
    method_id = expected$method_id,
    validate_optional = TRUE
  )
  if (is.null(batch_contract)) {
    if (isTRUE(require_recorded)) {
      stop(label, " is missing recorded corrected batch contract identity")
    }
    recorded <- NULL
    recorded_summary <- NULL
  } else {
    recorded <- .batch_identity_normalize(
      batch_contract, paste0(label, " (recorded)")
    )
    recorded_summary <- .batch_identity_summary(
      batch_contract,
      recorded$ordered_source_keys,
      paste0(label, " recorded identity"),
      method_id = recorded$method_id,
      required = isTRUE(summary_required),
      validate_optional = TRUE
    )
    if (!identical(recorded, expected)) {
      stop(
        label,
        " source/config identity does not match the expected corrected batch contract"
      )
    }
    if (!is.null(expected_summary) && !is.null(recorded_summary) &&
        !identical(recorded_summary, expected_summary)) {
      stop(label, " validation_summary does not match the expected corrected metadata")
    }
  }
  if (!is.null(source_obs_columns)) {
    if (!is.character(source_obs_columns)) {
      stop(label, " source obs columns are not strings")
    }
    missing <- expected$required_source_obs_columns[
      !expected$required_source_obs_columns %in% source_obs_columns
    ]
    if (length(missing)) {
      stop(label, " source obs columns are missing: ", paste(missing, collapse = ", "))
    }
    if ("__ecoda_batch_combined_v1" %in% source_obs_columns) {
      stop(
        label, " source obs contains reserved temporary column ",
        "__ecoda_batch_combined_v1"
      )
    }
  }
  if (!is.null(reserved_absent) &&
      (!is.logical(reserved_absent) || length(reserved_absent) != 1L ||
       is.na(reserved_absent) || !isTRUE(reserved_absent))) {
    stop(
      label, " source obs contains reserved temporary column ",
      "__ecoda_batch_combined_v1"
    )
  }
  if (is.null(recorded)) expected else recorded
}

load_batch_contract_argument <- function(value) {
  if (is.null(value) || !nzchar(value)) return(NULL)
  text <- if (file.exists(value)) {
    paste(readLines(value, warn = FALSE), collapse = "\n")
  } else {
    value
  }
  identity <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = FALSE),
    error = function(error) {
      stop(
        "batch contract identity must be a JSON object or a readable JSON path: ",
        conditionMessage(error)
      )
    }
  )
  if (!is.list(identity) || is.null(names(identity))) {
    stop("batch contract identity JSON must be an object")
  }
  identity
}

expected_batch_contract <- load_batch_contract_argument(expected_batch_contract_arg)
batch_contract <- load_batch_contract_argument(batch_contract_arg)
if (
  batch_pass == "corrected" &&
  !nzchar(artifact_path) &&
  (!nzchar(config_path) || !file.exists(config_path))
) {
  stop("corrected batch artifact validation requires an existing --config")
}
source_identity_records <- list()
source_identity_verified_keys <- new.env(parent = emptyenv())
if (nzchar(source_identity_path)) {
  if (!checksum_ok(source_identity_path)) {
    stop("source identity checksum is missing or invalid: ", source_identity_path)
  }
  identity <- jsonlite::fromJSON(source_identity_path, simplifyVector = FALSE)
  if (is.null(identity$schema) || as.integer(identity$schema) != 1L ||
      !is.list(identity$entries) || !length(identity$entries)) {
    stop("source identity manifest has an invalid schema: ", source_identity_path)
  }
  for (entry in identity$entries) {
    required <- c("dataset", "view", "path", "size", "md5", "sample_ids")
    if (!is.list(entry) || !identical(sort(names(entry)), sort(required))) {
      stop("source identity entry has invalid fields: ", source_identity_path)
    }
    ids <- as.character(unlist(entry$sample_ids, use.names = FALSE))
    if (!nzchar(entry$dataset) || !nzchar(entry$view) || !nzchar(entry$path) ||
        !grepl("^[0-9]+$", as.character(entry$size)) ||
        !grepl("^[[:xdigit:]]{32}$", entry$md5) ||
        !length(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
      stop("source identity entry is malformed: ", source_identity_path)
    }
    key <- paste(entry$dataset, entry$view, sep = "\t")
    if (!is.null(source_identity_records[[key]])) {
      stop("source identity contains duplicate rows: ", key)
    }
    source_identity_records[[key]] <- list(
      path = entry$path,
      size = as.character(entry$size),
      md5 = tolower(entry$md5),
      sample_ids = ids
    )
  }
}

source_reader <- NULL
read_h5ad_sample_ids_h5py <- function(path) {
  if (is.null(source_reader)) {
    module_dir <- normalizePath(file.path(validator_dir, "..", "utils", "py"), mustWork = TRUE)
    source_reader <<- reticulate::import_from_path(
      "h5ad_source_identity", path = module_dir, convert = TRUE
    )
  }
  as.character(source_reader$read_h5ad_sample_ids(path))
}

expected_samples <- function(ds, view) {
  if (!nzchar(input_root)) return(NULL)
  if (!nzchar(config_path) || !file.exists(config_path)) {
    stop("--input-root requires an existing --config")
  }
  config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  entry <- config[[ds]]
  view_spec <- if (is.null(entry)) NULL else entry$views[[view]]
  output <- if (is.null(view_spec)) NULL else {
    view_spec$output_file_name %||% view_spec$output_file
  }
  if (is.null(output) || length(output) != 1L || !nzchar(output)) {
    stop("missing h5ad output contract for ", ds, "/", view)
  }
  path <- file.path(input_root, ds, "output", output)
  canonical_path <- gsub("/{2,}", "/", path)
  key <- paste(ds, view, sep = "\t")
  if (nzchar(source_identity_path)) {
    record <- source_identity_records[[key]]
    if (is.null(record) || !identical(record$path, canonical_path)) {
      stop("source identity is missing or mismatched for ", ds, "/", view)
    }
    if (!source_identity_verified && !exists(key, envir = source_identity_verified_keys, inherits = FALSE)) {
      if (!file.exists(path) || as.character(file.info(path)$size) != record$size ||
          tolower(unname(tools::md5sum(path))) != record$md5) {
        stop("source identity digest mismatch for ", ds, "/", view)
      }
      current_ids <- read_h5ad_sample_ids_h5py(path)
      if (!identical(current_ids, record$sample_ids)) {
        stop("source identity Sample order mismatch for ", ds, "/", view)
      }
      assign(key, TRUE, envir = source_identity_verified_keys)
    }
    return(record$sample_ids)
  }
  if (!checksum_ok(path)) stop("missing or checksum-invalid input h5ad: ", path)
  read_h5ad_sample_ids_h5py(path)
}

finite_numeric <- function(value) {
  if (is.numeric(value)) return(all(is.finite(value)))
  if (is.data.frame(value)) {
    return(all(vapply(value, function(column) {
      is.numeric(column) && all(is.finite(column))
    }, logical(1L))))
  }
  if (is.list(value)) return(all(vapply(value, finite_numeric, logical(1L))))
  FALSE
}

config <- if (nzchar(config_path) && file.exists(config_path)) {
  jsonlite::fromJSON(config_path, simplifyVector = FALSE)
} else {
  list()
}

.corrected_batch_method_model <- function(label) {
  if (identical(label, "composition")) {
    return(list(
      method_id = "ECODA_authors_HR",
      model_id = "ecoda_additive_random_intercepts_v1"
    ))
  }
  if (identical(label, "gloscope")) {
    return(list(
      method_id = "GloScope",
      model_id = "embedding_consumer_harmony_v1"
    ))
  }
  if (label %in% c("prepare_pseudobulk", "pseudobulk")) {
    return(list(
      method_id = "Pseudobulk",
      model_id = "pseudobulk_composite_v1"
    ))
  }
  stop("unsupported corrected batch method for identity: ", label)
}

.corrected_batch_row_identity <- function(ds, view, label) {
  if (!identical(batch_pass, "corrected")) return(NULL)
  if (!nzchar(config_path) || !file.exists(config_path)) {
    stop("corrected batch artifact validation requires an existing --config")
  }
  entry <- config[[ds]]
  if (is.null(entry) || !is.list(entry)) {
    stop("dataset ", ds, " is missing from the selected config")
  }
  views <- entry[["views"]]
  view_spec <- if (is.list(views)) views[[view]] else NULL
  if (is.null(view_spec) || !is.list(view_spec)) {
    stop("dataset ", ds, " is missing the selected view ", view)
  }
  columns <- entry[["columns"]] %||% list()
  view_columns <- view_spec[["columns"]] %||% list()
  if (!is.list(columns) || !is.list(view_columns)) {
    stop("batch column configuration is malformed for ", ds, "/", view)
  }
  columns <- modifyList(columns, view_columns)
  batch_keys <- columns[["batch"]]
  if (is.null(batch_keys)) {
    stop(
      "corrected batch config is missing columns.batch for ",
      ds, "/", view
    )
  }
  .batch_identity_load_contract()
  method_model <- .corrected_batch_method_model(label)
  tryCatch(
    ecoda_batch_contract_identity(
      batch_keys = batch_keys,
      sample_col = "Sample",
      method_id = method_model[["method_id"]],
      model_id = method_model[["model_id"]]
    ),
    error = function(error) {
      stop(
        "invalid corrected batch config for ", ds, "/", view, "/", label,
        ": ", conditionMessage(error)
      )
    }
  )
}

.expected_row_batch_contract <- function(ds, view, label, supplied = NULL) {
  derived <- .corrected_batch_row_identity(ds, view, label)
  if (is.null(derived)) return(supplied)
  if (!is.null(supplied)) {
    validate_batch_contract_identity(
      derived,
      supplied,
      require_recorded = TRUE,
      require_summary = FALSE,
      label = paste0(ds, "/", view, "/", label, " supplied identity")
    )
    return(supplied)
  }
  derived
}
batch_required_keys <- function(ds, label) {
  if (label == "gloscope") return("GloScope_hvg2000_pcadims30")
  if (label == "pseudobulk") return("Pseudobulk_hvg2000")
  if (label == "composition") {
    # Batch composition has a deliberately smaller contract than ordinary
    # benchmark composition. Optional annotation bundles are legacy extras.
    return(c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2"))
  }
  NULL
}
batch_allowed_extra_keys <- function(ds, label) {
  if (label == "composition" && !identical(batch_pass, "corrected")) {
    return(c("ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR"))
  }
  character()
}

.validate_corrected_composition_nested <- function(
  bundle,
  expected_batch_contract,
  recorded_batch_contract,
  file
) {
  identity_source <- expected_batch_contract %||% recorded_batch_contract
  if (is.null(identity_source)) {
    stop("corrected composition bundle is missing its source identity: ", file)
  }

  # Configuration-only expected identities are used by the corrected CLI, so
  # retain the recorded top-level summary when the expected identity has none.
  source_summary_identity <- expected_batch_contract
  source_summary <- if (is.list(source_summary_identity)) {
    source_summary_identity[["validation_summary"]]
  } else {
    NULL
  }
  if (is.null(source_summary)) {
    source_summary_identity <- recorded_batch_contract
    source_summary <- if (is.list(source_summary_identity)) {
      source_summary_identity[["validation_summary"]]
    } else {
      NULL
    }
  }
  if (is.null(source_summary)) {
    stop(
      "corrected composition bundle is missing its source validation_summary: ",
      file
    )
  }

  normalized_source <- .batch_identity_normalize(
    identity_source,
    paste0("RDS composition (", file, ")")
  )
  source_summary <- .batch_identity_summary(
    source_summary_identity,
    normalized_source$ordered_source_keys,
    paste0("RDS composition (", file, ") source identity"),
    method_id = normalized_source$method_id,
    required = TRUE
  )
  methods <- c(
    "ECODA_authors_HR",
    "ECODA_authors_HR_NULL",
    "ECODA_seuratres_2"
  )
  .batch_identity_load_contract()
  for (method_id in methods) {
    combo <- bundle[[method_id]]
    if (!is.list(combo)) {
      stop(
        "corrected composition bundle ", method_id,
        " is missing or not a list: ", file
      )
    }
    batch_keys <- as.list(unname(normalized_source$ordered_source_keys))
    expected_nested <- tryCatch(
      ecoda_batch_contract_identity(
        batch_keys = batch_keys,
        sample_col = "Sample",
        method_id = method_id,
        model_id = "ecoda_additive_random_intercepts_v1"
      ),
      error = function(error) {
        stop(
          "invalid corrected composition identity for ", method_id,
          " in ", file, ": ", conditionMessage(error)
        )
      }
    )
    correction_spec <- tryCatch(
      ecoda_batch_correction_spec(
        method_id = method_id,
        batch_keys = batch_keys
      ),
      error = function(error) {
        stop(
          "invalid corrected composition correction policy for ", method_id,
          " in ", file, ": ", conditionMessage(error)
        )
      }
    )
    nested_summary <- source_summary
    nested_summary[["correction_mode"]] <- correction_spec[["correction_mode"]]
    nested_summary[["correction_formula"]] <- correction_spec[["correction_formula"]]
    expected_nested[["validation_summary"]] <- nested_summary
    .rds_batch_contract_values(
      combo,
      expected_nested,
      NULL,
      paste0("RDS composition ", method_id, " (", file, ")")
    )
  }
  invisible(NULL)
}

# scITD may legitimately emit an ordered subset of source samples; every
# other method must retain the complete source sample universe.
report_sample_universe <- function(ids, expected, method, file) {
  if (is.null(expected)) return(invisible(NULL))
  ids <- as.character(ids)
  expected <- as.character(expected)
  if (identical(ids, expected)) return(invisible(NULL))
  if (!identical(method, "scitd")) {
    stop("sample identifiers do not match the ordered selected h5ad: ", file)
  }
  expected_subset <- expected[expected %in% ids]
  if (!identical(ids, expected_subset)) {
    stop("scITD sample identifiers must be an ordered subset of the selected h5ad: ", file)
  }
  dropped <- expected[!expected %in% ids]
  if (length(dropped)) {
    message(
      "scITD sample-universe exception: dropped sample IDs: ",
      paste(dropped, collapse = ", "), " (", file, ")"
    )
  }
  invisible(dropped)
}

validate_matrix <- function(matrix_value, labels_value, file, expected = NULL, method = "") {
  dimensions <- dim(matrix_value)
  if (is.null(dimensions) || length(dimensions) < 2L ||
      dimensions[[1L]] <= 0L || dimensions[[2L]] <= 0L) {
    stop("feat_mat is empty or not matrix-like: ", file)
  }
  if (is.data.frame(matrix_value)) {
    if (!all(vapply(matrix_value, is.numeric, logical(1L)))) {
      stop("feat_mat contains nonnumeric features: ", file)
    }
  } else if (!is.numeric(matrix_value)) {
    stop("feat_mat contains nonnumeric features: ", file)
  }
  if (!finite_numeric(matrix_value)) stop("feat_mat contains nonfinite values: ", file)
  ids <- rownames(matrix_value)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("feat_mat has missing or duplicate sample identifiers: ", file)
  }
  if (length(labels_value) != nrow(matrix_value) || anyNA(labels_value)) {
    stop("labels are missing or misaligned with feat_mat: ", file)
  }
  label_ids <- names(labels_value)
  if (is.null(label_ids) || !identical(as.character(label_ids), as.character(ids))) {
    stop("labels names do not exactly match feat_mat row names: ", file)
  }
  report_sample_universe(ids, expected, method, file)
}

validate_dist <- function(dist_value, n, ids, file) {
  if (inherits(dist_value, "dist")) {
    if (attr(dist_value, "Size") != n || !finite_numeric(as.numeric(dist_value))) {
      stop("dist_mat is invalid or nonfinite: ", file)
    }
    dist_ids <- attr(dist_value, "Labels")
    if (!is.null(dist_ids) && !identical(as.character(dist_ids), as.character(ids))) {
      stop("dist_mat labels are misaligned: ", file)
    }
    return(invisible(NULL))
  }
  dimensions <- dim(dist_value)
  if (is.null(dimensions) || length(dimensions) < 2L ||
      dimensions[[1L]] != n || dimensions[[2L]] != n ||
      !is.numeric(dist_value) || !finite_numeric(dist_value)) {
    stop("dist_mat dimensions or values are invalid: ", file)
  }
  if (!is.null(rownames(dist_value)) &&
      !identical(as.character(rownames(dist_value)), as.character(ids))) {
    stop("dist_mat row identifiers are misaligned: ", file)
  }
  if (!is.null(colnames(dist_value)) &&
      !identical(as.character(colnames(dist_value)), as.character(ids))) {
    stop("dist_mat column identifiers are misaligned: ", file)
  }
}

validate_combo <- function(combo, file, expected = NULL, method = "") {
  required <- c("scores", "feat_mat", "dist_mat", "labels")
  missing <- setdiff(required, names(combo))
  if (length(missing)) stop("result combo missing fields in ", file, ": ", paste(missing, collapse = ", "))
  if (is.null(combo$scores) || length(combo$scores) == 0L ||
      !finite_numeric(combo$scores)) stop("scores are empty or nonfinite: ", file)
  validate_matrix(combo$feat_mat, combo$labels, file, expected, method)
  validate_dist(combo$dist_mat, nrow(combo$feat_mat), rownames(combo$feat_mat), file)
}

.rds_embedded_batch_contract <- function(value, label) {
  if (!is.list(value)) return(NULL)
  fields <- c(
    "batch_contract", "batch_contract_identity",
    "ecoda_batch_contract", "_ecoda_batch_contract"
  )
  present <- fields[
    fields %in% names(value) &
      vapply(fields, function(field) !is.null(value[[field]]), logical(1L))
  ]
  if (!length(present)) return(NULL)
  embedded <- value[[present[[1L]]]]
  if (length(present) > 1L) {
    for (field in present[-1L]) {
      validate_batch_contract_identity(
        embedded,
        value[[field]],
        require_recorded = TRUE,
        require_summary = identical(batch_pass, "corrected"),
        label = paste0(label, " embedded identity aliases")
      )
    }
  }
  embedded
}

.rds_batch_contract_values <- function(
  value,
  expected_batch_contract,
  batch_contract,
  label,
  require_summary = NULL
) {
  embedded <- .rds_embedded_batch_contract(value, label)
  summary_required <- if (is.null(require_summary)) {
    identical(batch_pass, "corrected")
  } else {
    isTRUE(require_summary)
  }
  if (summary_required && is.null(embedded)) {
    stop(label, " is missing embedded corrected batch contract identity")
  }
  if (!is.null(batch_contract) && !is.null(embedded)) {
    validate_batch_contract_identity(
      batch_contract,
      embedded,
      require_recorded = TRUE,
      require_summary = summary_required,
      label = paste0(label, " embedded identity")
    )
  }
  recorded <- embedded %||% batch_contract
  validate_batch_contract_identity(
    expected_batch_contract,
    recorded,
    require_recorded = !is.null(expected_batch_contract),
    require_summary = summary_required,
    label = label
  )
  list(embedded = embedded, recorded = recorded)
}

validate_result_file <- function(
  file,
  expected = NULL,
  required_keys = NULL,
  allowed_extra_keys = character(),
  method = "",
  expected_batch_contract = NULL,
  batch_contract = NULL
) {
  if (identical(batch_pass, "corrected") && identical(method, "composition")) {
    allowed_extra_keys <- character()
  }
  if (!checksum_ok(file)) stop("Missing or invalid result checksum: ", file)
  bundle <- readRDS(file)
  identity_values <- .rds_batch_contract_values(
    bundle,
    expected_batch_contract,
    batch_contract,
    paste0("RDS ", method, " (", file, ")")
  )
  recorded_batch_contract <- identity_values[["recorded"]]
  identity_validation_active <- !is.null(expected_batch_contract) ||
    !is.null(batch_contract) || !is.null(recorded_batch_contract)
  if (identity_validation_active && is.list(bundle)) {
    bundle <- bundle[
      setdiff(
        names(bundle),
        c(
          "batch_contract", "batch_contract_identity",
          "ecoda_batch_contract", "_ecoda_batch_contract"
        )
      )
    ]
  }
  required <- c("scores", "feat_mat", "dist_mat", "labels")
  is_combo <- is.list(bundle) && all(required %in% names(bundle))
  if (is_combo) {
    if (!is.null(required_keys)) {
      stop("result bundle keys do not match the method contract: ", file)
    }
    combos <- list(bundle)
  } else {
    if (!is.list(bundle) || !length(bundle) || is.null(names(bundle)) ||
        any(!nzchar(names(bundle)))) {
      stop("Result bundle is empty or unnamed: ", file)
    }
    if (!is.null(required_keys)) {
      actual_keys <- if (
        !is.null(expected_batch_contract) || !is.null(batch_contract)
      ) {
        setdiff(
          names(bundle),
          c(
            "batch_contract", "batch_contract_identity",
            "ecoda_batch_contract", "_ecoda_batch_contract"
          )
        )
      } else {
        names(bundle)
      }
      allowed_keys <- unique(c(required_keys, allowed_extra_keys))
      missing_keys <- setdiff(required_keys, actual_keys)
      unexpected_keys <- setdiff(actual_keys, allowed_keys)
      if (length(missing_keys) > 0L || length(unexpected_keys) > 0L ||
          anyDuplicated(actual_keys)) {
        stop("result combo keys do not match the method contract: ", file)
      }
    }
    combos <- bundle
  }
  if (identical(batch_pass, "corrected") && identical(method, "pseudobulk")) {
    nested <- if (is.list(bundle)) bundle[["Pseudobulk_hvg2000"]] else NULL
    if (!is.list(nested)) {
      stop(
        "corrected pseudobulk bundle Pseudobulk_hvg2000 is missing or not a list: ",
        file
      )
    }
    identity_source <- expected_batch_contract %||% recorded_batch_contract
    if (is.null(identity_source)) {
      stop("corrected pseudobulk bundle is missing its source identity: ", file)
    }
    summary_source_identity <- expected_batch_contract
    source_summary <- if (is.list(summary_source_identity)) {
      summary_source_identity[["validation_summary"]]
    } else {
      NULL
    }
    if (is.null(source_summary)) {
      summary_source_identity <- recorded_batch_contract
      source_summary <- if (is.list(summary_source_identity)) {
        summary_source_identity[["validation_summary"]]
      } else {
        NULL
      }
    }
    if (is.null(source_summary)) {
      stop(
        "corrected pseudobulk bundle is missing its source validation_summary: ",
        file
      )
    }
    normalized_source <- .batch_identity_normalize(
      identity_source,
      paste0("RDS pseudobulk (", file, ")")
    )
    source_summary <- .batch_identity_summary(
      summary_source_identity,
      normalized_source$ordered_source_keys,
      paste0("RDS pseudobulk (", file, ") source identity"),
      method_id = "Pseudobulk",
      required = TRUE
    )
    .batch_identity_load_contract()
    expected_nested <- tryCatch(
      ecoda_batch_contract_identity(
        batch_keys = as.list(unname(normalized_source$ordered_source_keys)),
        sample_col = "Sample",
        method_id = "Pseudobulk",
        model_id = "pseudobulk_composite_v1"
      ),
      error = function(error) {
        stop(
          "invalid corrected pseudobulk nested identity in ", file, ": ",
          conditionMessage(error)
        )
      }
    )
    expected_nested[["validation_summary"]] <- source_summary
    .rds_batch_contract_values(
      nested,
      expected_nested,
      NULL,
      paste0("RDS pseudobulk Pseudobulk_hvg2000 (", file, ")"),
      require_summary = TRUE
    )
  }
  if (identical(batch_pass, "corrected") && identical(method, "composition")) {
    .validate_corrected_composition_nested(
      bundle,
      expected_batch_contract,
      recorded_batch_contract,
      file
    )
  }
  for (combo in combos) {
    if (!is.list(combo)) stop("Result combo is not a list: ", file)
    validate_combo(combo, file, expected, method)
  }
}
validate_timing_scalar <- function(value, field, file) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop("Pseudobulk ", field, " is invalid: ", file)
  }
  invisible(TRUE)
}

validate_memory_scalar <- function(value, field, file) {
  if (is.null(value)) return(invisible(TRUE))
  if (!is.numeric(value) || length(value) != 1L) {
    stop("Pseudobulk ", field, " is invalid: ", file)
  }
  if (is.na(value)) {
    if (is.nan(value)) stop("Pseudobulk ", field, " is invalid: ", file)
    return(invisible(TRUE))
  }
  if (!is.finite(value) || value < 0) {
    stop("Pseudobulk ", field, " is invalid: ", file)
  }
  invisible(TRUE)
}

validate_pseudobulk <- function(
  file,
  expected = NULL,
  expected_batch_contract = NULL,
  batch_contract = NULL
){
  if (!checksum_ok(file)) stop("Missing or invalid pseudobulk checksum: ", file)
  value <- readRDS(file)
  identity_values <- .rds_batch_contract_values(
    value,
    expected_batch_contract,
    batch_contract,
    paste0("RDS pseudobulk (", file, ")")
  )
  recorded_batch_contract <- identity_values[["recorded"]]
  if (
    (!is.null(expected_batch_contract) ||
     !is.null(batch_contract) ||
     !is.null(recorded_batch_contract)) &&
    is.list(value)
  ) {
    value <- value[
      setdiff(
        names(value),
        c(
          "batch_contract", "batch_contract_identity",
          "ecoda_batch_contract", "_ecoda_batch_contract"
        )
      )
    ]
  }
  timing <- NULL
  memory <- NULL
  if (is.list(value) && !is.data.frame(value) && !is.null(value$pb)) {
    timing <- value$time_secs
    memory <- value$mem_GB
    schema_fields <- c(
      "timing_schema", "aggregate_time_secs", "shared_fit_time_secs",
      "shared_time_secs", "variant_time_secs", "shared_mem_GB", "timing_id"
    )
    present_schema_fields <- intersect(names(value), schema_fields)
    if (length(present_schema_fields)) {
      required_schema2 <- c(
        "pb", "time_secs", "mem_GB", "aggregate_time_secs",
        "shared_fit_time_secs", "shared_time_secs", "variant_time_secs",
        "shared_mem_GB", "timing_id", "timing_schema"
      )
      actual_fields <- if (
        !is.null(expected_batch_contract) || !is.null(batch_contract)
      ) {
        setdiff(
          names(value),
          c(
            "batch_contract", "batch_contract_identity",
            "ecoda_batch_contract", "_ecoda_batch_contract"
          )
        )
      } else {
        names(value)
      }
      if (is.null(actual_fields) ||
          length(actual_fields) != length(required_schema2) ||
          anyDuplicated(actual_fields) ||
          !setequal(actual_fields, required_schema2) ||
          !is.numeric(value$timing_schema) ||
          length(value$timing_schema) != 1L ||
          is.na(value$timing_schema) ||
          !is.finite(value$timing_schema) ||
          value$timing_schema != 2 ||
          value$timing_schema != floor(value$timing_schema)) {
        stop("Pseudobulk schema-2 timing fields are invalid: ", file)
      }
      for (field in c(
        "time_secs", "aggregate_time_secs", "shared_fit_time_secs",
        "shared_time_secs", "variant_time_secs"
      )) {
        validate_timing_scalar(value[[field]], field, file)
      }
      for (field in c("mem_GB", "shared_mem_GB")) {
        validate_memory_scalar(value[[field]], field, file)
      }
      timing_id <- value$timing_id
      timing_id_parts <- if (is.character(timing_id) && length(timing_id) == 1L &&
                             !is.na(timing_id)) {
        strsplit(timing_id, ":", fixed = TRUE)[[1L]]
      } else {
        character()
      }
      if (!is.character(timing_id) || length(timing_id) != 1L ||
          is.na(timing_id) || !nzchar(trimws(timing_id)) ||
          length(timing_id_parts) != 4L ||
          any(!nzchar(trimws(timing_id_parts))) ||
          grepl("[[:cntrl:]]", timing_id, perl = TRUE)) {
        stop("Pseudobulk timing_id is invalid: ", file)
      }
      shared_time <- value$aggregate_time_secs +
        value$shared_fit_time_secs
      if (!is.finite(shared_time) ||
          !isTRUE(all.equal(
            as.numeric(value$shared_time_secs),
            as.numeric(shared_time),
            tolerance = 0
          )) ||
          !isTRUE(all.equal(
            as.numeric(value$time_secs),
            as.numeric(value$variant_time_secs),
            tolerance = 0
          ))) {
        stop("Pseudobulk schema-2 timing totals are inconsistent: ", file)
      }
    } else if (length(intersect(
      names(value),
      c(
        "aggregate_time_secs", "shared_fit_time_secs", "shared_time_secs",
        "variant_time_secs", "shared_mem_GB", "timing_id"
      )
    ))) {
      stop("Pseudobulk schema-2 timing fields are incomplete: ", file)
    }
    value <- value$pb
  }
  dimensions <- dim(value)
  if (is.null(dimensions) || length(dimensions) < 2L ||
      dimensions[[1L]] <= 0L || dimensions[[2L]] <= 0L) {
    stop("Pseudobulk artifact is empty: ", file)
  }
  if (is.data.frame(value)) {
    if (!all(vapply(value, is.numeric, logical(1L)))) {
      stop("Pseudobulk values are nonnumeric: ", file)
    }
  } else if (!is.numeric(value)) {
    stop("Pseudobulk values are nonnumeric: ", file)
  }
  ids <- rownames(value)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Pseudobulk identifiers invalid: ", file)
  }
  if (!is.null(expected) && !identical(as.character(ids), as.character(expected))) {
    stop("Pseudobulk sample identifiers do not match ordered selected h5ad: ", file)
  }
  if (!finite_numeric(value)) stop("Pseudobulk values are nonfinite: ", file)
  if (!is.null(timing)) validate_timing_scalar(timing, "timing", file)
  if (!is.null(memory)) validate_memory_scalar(memory, "memory", file)
}

validate_trans <- function(file) {
  if (!checksum_ok(file)) stop("Missing or invalid transformation checksum: ", file)
  value <- readRDS(file)
  required <- c("trans_method", "ANOSIM_score", "Modularity_score", "Adjusted_Rand_Index")
  if (!is.data.frame(value) || nrow(value) <= 0L || !all(required %in% names(value))) {
    stop("Transformation result has the wrong summarized schema: ", file)
  }
  scores <- value[required[2:4]]
  if (!all(vapply(scores, is.numeric, logical(1L))) || !finite_numeric(scores)) {
    stop("Transformation scores are nonnumeric or nonfinite: ", file)
  }
}
validate_zeroimp <- function(file) {
  if (!checksum_ok(file)) stop("Missing or invalid zero-imputation checksum: ", file)
  value <- readRDS(file)
  if (!is.list(value) || !length(value) || is.null(names(value)) ||
      any(!nzchar(names(value))) || anyDuplicated(names(value))) {
    stop("Zero-imputation result must be a nonempty named list: ", file)
  }
  if (!all(vapply(value, function(score) is.list(score) && length(score) > 0L &&
                  finite_numeric(score), logical(1L)))) {
    stop("Zero-imputation scores are malformed or nonfinite: ", file)
  }
}
validate_metadata <- function(
  file,
  expected = NULL,
  expected_batch_contract = NULL,
  batch_contract = NULL
) {
  if (!checksum_ok(file)) stop("Missing or invalid metadata checksum: ", file)
  value <- readRDS(file)
  identity_values <- .rds_batch_contract_values(
    value,
    expected_batch_contract,
    batch_contract,
    paste0("RDS metadata (", file, ")")
  )
  recorded_batch_contract <- identity_values[["recorded"]]
  required <- c("labels", "n_cells", "n_samples", "cells_per_sample")
  if (!is.list(value) || !all(required %in% names(value))) {
    stop("composition metadata bundle is malformed: ", file)
  }
  numeric_fields <- value[c("n_cells", "n_samples", "cells_per_sample")]
  if (!is.numeric(value$n_cells) || length(value$n_cells) != 1L ||
      !is.numeric(value$n_samples) || length(value$n_samples) != 1L ||
      !is.numeric(value$cells_per_sample) ||
      !finite_numeric(numeric_fields) ||
      value$n_cells <= 0 || value$n_samples <= 0) {
    stop("composition metadata numeric fields are malformed: ", file)
  }
  ids <- names(value$labels)
  if (length(value$labels) == 0L || is.null(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids) || anyNA(value$labels) ||
      any(!nzchar(trimws(as.character(value$labels))))) {
    stop("composition metadata labels are unnamed, blank, or duplicated: ", file)
  }
  if (value$n_samples != length(ids) ||
      length(value$cells_per_sample) != length(ids)) {
    stop("composition metadata sample counts do not match labels: ", file)
  }
  if (!is.null(expected) && !identical(as.character(ids), as.character(expected))) {
    stop("composition metadata sample identifiers do not match ordered selected h5ad: ", file)
  }
  if (!identical(as.character(names(value$cells_per_sample)), as.character(ids))) {
    stop("composition metadata sample order is invalid: ", file)
  }
}

validate_artifact_contract <- function(
  file,
  method,
  ds = "",
  view = "",
  metadata = FALSE,
  expected_batch_contract = NULL,
  batch_contract = NULL
) {
  row_expected_batch_contract <- expected_batch_contract
  if (
    identical(view, "batch_effect_corrected") &&
    nzchar(ds) && nzchar(view)
  ) {
    row_expected_batch_contract <- .expected_row_batch_contract(
      ds, view, method, expected_batch_contract
    )
  }
  if (
    identical(view, "batch_effect_corrected") &&
    is.null(row_expected_batch_contract) && is.null(batch_contract)
  ) {
    stop(
      "corrected batch artifact validation requires explicit contract identity"
    )
  }
  expected <- if (nzchar(ds) && nzchar(view)) {
    expected_samples(ds, view)
  } else {
    NULL
  }
  if (metadata) {
    validate_metadata(
      file,
      expected,
      row_expected_batch_contract,
      batch_contract
    )
  } else if (method == "trans") {
    validate_trans(file)
  } else if (method == "zeroimp") {
    validate_zeroimp(file)
  } else if (method == "prepare_pseudobulk" ||
             grepl("pseudobulk_", basename(file), fixed = TRUE)) {
    validate_pseudobulk(
      file,
      expected,
      row_expected_batch_contract,
      batch_contract
    )
  } else {
    required_keys <- if (batch && method %in% c("gloscope", "pseudobulk", "composition")) {
      batch_required_keys(ds, method)
    } else {
      NULL
    }
    allowed_extra_keys <- if (batch) batch_allowed_extra_keys(ds, method) else character()
    validate_result_file(
      file,
      expected,
      required_keys,
      allowed_extra_keys,
      method,
      row_expected_batch_contract,
      batch_contract
    )
  }
}

if (nzchar(artifact_path)) {
  validate_artifact_contract(
    artifact_path,
    method_arg,
    dataset_arg,
    view_arg,
    metadata_kind,
    expected_batch_contract,
    batch_contract
  )
  reject_selected_partials(artifact_path)
  cat("benchmark RDS artifact contract OK\n")
  quit(save = "no", status = 0)
}

if (nzchar(artifact_list)) {
  selected_artifacts <- artifact_list
  for (part in artifact_list_parts) {
    selected_artifacts <- c(selected_artifacts, part[[1L]])
    validate_artifact_contract(
      part[[1L]],
      part[[2L]],
      part[[3L]],
      part[[4L]],
      part[[5L]] == "1",
      expected_batch_contract,
      batch_contract
    )
  }
  reject_selected_partials(selected_artifacts)
  cat("benchmark RDS artifact-list contract OK\n")
  quit(save = "no", status = 0)
}

if (!dir.exists(root)) stop("benchmark result root is missing: ", root)
if (batch && !batch_pass %in% c("uncorrected", "corrected")) {
  stop("batch validation requires uncorrected or corrected pass")
}
if (
  batch &&
  batch_pass == "corrected" &&
  (!nzchar(config_path) || !file.exists(config_path))
) {
  stop("corrected batch artifact validation requires an existing --config")
}
if (batch && exact) {
  expected_rows <- paste(
    c("Alzheimer", "Breast_cancer", "Covid19_PBMC", "Kidney_KPMP_full",
      "Myocardial_infarction", "Diabetes", "Lupus_PBMC", "Lung",
      "Parkinson", "Joanito", "Stephenson", "CombinedPBMC"),
    "batch_effect_uncorrected", "batch_effect_uncorrected", sep = "\t"
  )
  if (!identical(selection_rows, expected_rows) || batch_pass != "uncorrected") {
    stop("batch exact selection is not the literal twelve-row uncorrected matrix")
  }
}

selected_artifacts <- selection
for (part in parts) {
  ds <- part[[1L]]
  view <- part[[2L]]
  scope <- part[[3L]]
  if (batch) {
    expected_view <- paste0("batch_effect_", batch_pass)
    if (view != expected_view) stop("batch selection view mismatch: ", ds, "/", view)
    if (scope != expected_view) stop("batch selection scope mismatch: ", ds, "/", scope)
    selected_labels <- labels
  } else {
    selected_labels <- if (exact) scope else labels
    if (exact && !scope %in% labels) stop("selection scope is not selected: ", scope)
  }
  expected <- expected_samples(ds, view)
  for (label in selected_labels) {
    row_expected_batch_contract <- if (
      batch && batch_pass == "corrected" &&
      !label %in% c("trans", "zeroimp")
    ) {
      .expected_row_batch_contract(ds, view, label, expected_batch_contract)
    } else {
      expected_batch_contract
    }
    if (label == "prepare_pseudobulk") {
      variants <- if (batch) "hvg2000" else c("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")
      stem <- if (batch) paste0(ds, "_batch_effect_", batch_pass) else ds
      for (variant in variants) {
        file <- file.path(root, "pseudobulks", paste0(stem, "_pseudobulk_", variant, ".rds"))
        selected_artifacts <- c(selected_artifacts, file)
        validate_pseudobulk(
          file,
          expected,
          row_expected_batch_contract,
          batch_contract
        )
      }
    } else if (label == "trans") {
      file <- file.path(root, "results", paste0(ds, "_trans.rds"))
      selected_artifacts <- c(selected_artifacts, file)
      validate_trans(file)
    } else if (label == "zeroimp") {
      file <- file.path(root, "results", paste0(ds, "_zeroimp.rds"))
      selected_artifacts <- c(selected_artifacts, file)
      validate_zeroimp(file)
    } else {
      stem <- if (batch) paste0(ds, "_batch_effect_", batch_pass) else ds
      file <- file.path(root, "results", paste0(stem, "_", label, ".rds"))
      selected_artifacts <- c(selected_artifacts, file)
      required_keys <- if (batch) batch_required_keys(ds, label) else NULL
      allowed_extra_keys <- if (batch) batch_allowed_extra_keys(ds, label) else character()
      validate_result_file(
        file,
        expected,
        required_keys,
        allowed_extra_keys,
        label,
        row_expected_batch_contract,
        batch_contract
      )
      if (label == "composition") {
        metadata_file <- file.path(root, "results", paste0(stem, "_metadata.rds"))
        selected_artifacts <- c(selected_artifacts, metadata_file)
        validate_metadata(
          metadata_file,
          expected,
          row_expected_batch_contract,
          batch_contract
        )
      }
    }
  }
}
reject_selected_partials(selected_artifacts)
cat("benchmark RDS bundle contract OK\n")
