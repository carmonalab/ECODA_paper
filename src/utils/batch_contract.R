# Corrected-mode batch contract.
#
# This module is deliberately independent of the raw configuration readers and
# of any corrected-mode consumer.  It validates the full cell-level metadata
# before reducing it to one record per Sample and provides byte-stable values,
# composite tokens, configuration fingerprints, and run-owned metadata.

.ecoda_batch_token_version <- "ecoda_batch_composite_v1"
.ecoda_batch_contract_version <- "ecoda_batch_contract_v1"
.ecoda_batch_reserved_name <- "__ecoda_batch_combined_v1"
.ecoda_batch_sentinels <- c("na", "nan", "none", "<na>", "n/a", "null", "unknown")

# Corrected method IDs are part of cache/run identity.  ``preprocess`` is the
# preprocessing boundary; the remaining entries are the exact Stage 5 labels.
.ecoda_batch_method_ids <- c(
  "preprocess",
  "ECODA_authors_HR",
  "ECODA_seuratres_2",
  "ECODA_authors_HR_NULL",
  "Pseudobulk",
  "GloScope",
  "PILOT",
  "MrVI",
  "QOT"
)
.ecoda_batch_model_ids <- c(
  "hvg_composite_v1",
  "harmony_native_list_v1",
  "ecoda_additive_random_intercepts_v1",
  "pseudobulk_composite_v1",
  "mrvi_composite_v1",
  "embedding_consumer_harmony_v1"
)

.ecoda_batch_stop <- function(...) {
  stop(..., call. = FALSE)
}

.ecoda_batch_utf8_text <- function(value, label = "value") {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    .ecoda_batch_stop(label, " must be one valid UTF-8 string")
  }
  value <- enc2utf8(unname(value))
  if (length(value) != 1L || is.na(value) || !isTRUE(validUTF8(value))) {
    .ecoda_batch_stop(label, " must be one valid UTF-8 string")
  }
  value
}

.ecoda_batch_utf8_raw <- function(value, label = "value") {
  charToRaw(.ecoda_batch_utf8_text(value, label))
}

.ecoda_batch_hex <- function(value) {
  if (length(value) == 0L) return("")
  tolower(paste(sprintf("%02x", as.integer(value)), collapse = ""))
}

.ecoda_batch_raw_sort <- function(values, label = "levels") {
  if (length(values) == 0L) return(character())
  values <- vapply(seq_along(values), function(index) {
    .ecoda_batch_utf8_text(values[[index]], paste0(label, "[", index, "]"))
  }, character(1))
  # Hexadecimal UTF-8 is ASCII and preserves unsigned raw-byte lexicographic
  # order.  Radix ordering avoids locale-dependent collation.
  raw_hex <- vapply(values, function(value) {
    .ecoda_batch_hex(.ecoda_batch_utf8_raw(value, label))
  }, character(1))
  values[order(raw_hex, method = "radix")]
}

.ecoda_batch_unique_raw_sorted <- function(values, label = "levels") {
  values <- unique(as.character(values))
  .ecoda_batch_raw_sort(values, label)
}

.ecoda_batch_config_name <- function(value, label, allow_null = FALSE) {
  if (is.null(value) && isTRUE(allow_null)) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    .ecoda_batch_stop(label, " must be one nonblank string")
  }
  value <- .ecoda_batch_utf8_text(value, label)
  if (!nzchar(value) || !identical(trimws(value), value)) {
    .ecoda_batch_stop(label, " must be one nonblank string without surrounding whitespace")
  }
  value
}

.ecoda_batch_key_vector <- function(batch_keys) {
  if (is.character(batch_keys) && length(batch_keys) == 1L) {
    keys <- unname(batch_keys)
  } else if (is.list(batch_keys) && length(batch_keys) > 0L) {
    keys <- vapply(seq_along(batch_keys), function(index) {
      value <- batch_keys[[index]]
      if (!is.character(value) || length(value) != 1L || is.na(value)) {
        .ecoda_batch_stop(
          "batch keys must be a scalar string or an ordered nonempty list of strings"
        )
      }
      unname(value)
    }, character(1))
  } else {
    .ecoda_batch_stop(
      "batch keys must be a scalar string or an ordered nonempty list of strings"
    )
  }
  if (length(keys) == 0L || anyNA(keys)) {
    .ecoda_batch_stop("batch keys must be a scalar string or an ordered nonempty list of strings")
  }
  keys <- vapply(seq_along(keys), function(index) {
    .ecoda_batch_config_name(keys[[index]], paste0("batch key ", index))
  }, character(1))
  if (anyDuplicated(keys)) {
    .ecoda_batch_stop("batch keys must not contain duplicates")
  }
  if (any(keys == .ecoda_batch_reserved_name)) {
    .ecoda_batch_stop(
      "batch key uses the reserved temporary name: ", .ecoda_batch_reserved_name
    )
  }
  unname(keys)
}

# Normalize the scalar/list configuration without trimming, reordering, or
# silently converting an atomic vector of several keys into a list.
ecoda_batch_normalize_keys <- function(
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL
) {
  sample_col <- .ecoda_batch_config_name(sample_col, "sample column")
  biological_label <- .ecoda_batch_config_name(
    biological_label, "biological label", allow_null = TRUE
  )
  if (identical(sample_col, .ecoda_batch_reserved_name)) {
    .ecoda_batch_stop("sample column uses the reserved temporary name")
  }
  keys <- .ecoda_batch_key_vector(batch_keys)
  if (any(keys == sample_col)) {
    .ecoda_batch_stop("batch keys must not include the standardized Sample column")
  }
  if (!is.null(biological_label) && any(keys == biological_label)) {
    .ecoda_batch_stop("batch keys must not include the biological label column")
  }
  keys
}

.ecoda_batch_is_missing_or_sentinel <- function(value) {
  if (is.null(value)) return(TRUE)
  if (length(value) != 1L) return(TRUE)
  if (is.factor(value)) {
    text <- as.character(value)
    if (length(text) != 1L || is.na(text)) return(TRUE)
  } else if (is.character(value)) {
    text <- value
    if (is.na(text)) return(TRUE)
  } else if (is.logical(value) || is.integer(value) || is.double(value)) {
    if (is.na(value) || (is.double(value) && is.nan(value))) return(TRUE)
    return(FALSE)
  } else {
    return(FALSE)
  }
  text <- enc2utf8(unname(text))
  if (is.na(text) || !isTRUE(validUTF8(text))) return(TRUE)
  trimmed <- trimws(text)
  !nzchar(trimmed) || tolower(trimmed) %in% .ecoda_batch_sentinels
}

.ecoda_batch_assert_scalar_supported <- function(value, label) {
  if (is.list(value) || is.data.frame(value)) {
    .ecoda_batch_stop(label, " has an unsupported list/object value")
  }
  if (is.object(value) && !is.factor(value)) {
    .ecoda_batch_stop(label, " has an unsupported date/object value")
  }
  if (!(is.factor(value) || is.character(value) || is.logical(value) ||
        is.integer(value) || is.double(value))) {
    .ecoda_batch_stop(label, " has an unsupported value type")
  }
  invisible(NULL)
}

.ecoda_batch_double_bits <- function(value, label) {
  connection <- rawConnection(raw(0L), open = "wb")
  on.exit(if (isOpen(connection)) close(connection), add = TRUE)
  writeBin(as.numeric(value), connection, size = 8L, endian = "big")
  bytes <- rawConnectionValue(connection)
  close(connection)
  if (length(bytes) != 8L) {
    .ecoda_batch_stop(label, " could not be encoded as binary64")
  }
  .ecoda_batch_hex(bytes)
}

# Canonicalize one nonmissing category value.  Character/factor values retain
# their exact UTF-8 label; numeric values retain their R scalar type.
ecoda_batch_canonical_value <- function(value, label = "batch value") {
  if (length(value) != 1L) {
    .ecoda_batch_stop(label, " must be one scalar value")
  }
  .ecoda_batch_assert_scalar_supported(value, label)
  if (.ecoda_batch_is_missing_or_sentinel(value)) {
    .ecoda_batch_stop(label, " is missing, blank, or a forbidden sentinel")
  }
  if (is.factor(value) || is.character(value)) {
    text <- if (is.factor(value)) as.character(value) else value
    text <- .ecoda_batch_utf8_text(text, label)
    return(paste0("s:", text))
  }
  if (is.logical(value)) {
    return(if (isTRUE(value)) "b:true" else "b:false")
  }
  if (is.integer(value)) {
    # R integer scalars are signed 32-bit values; as.character is exact
    # decimal and never introduces a floating-point spelling.
    return(paste0("i:", as.character(value)))
  }
  if (is.double(value)) {
    if (!isTRUE(is.finite(value))) {
      .ecoda_batch_stop(label, " must be finite when numeric")
    }
    return(paste0("f64:", .ecoda_batch_double_bits(value, label)))
  }
  .ecoda_batch_stop(label, " has an unsupported value type")
}

# Canonicalize every scalar in one metadata column without mutating the source.
ecoda_batch_canonical_values <- function(values, label = "batch column") {
  if (is.list(values) && !is.factor(values)) {
    .ecoda_batch_stop(label, " is list-valued and cannot be canonicalized")
  }
  if (length(values) == 0L) return(character())
  vapply(seq_along(values), function(index) {
    ecoda_batch_canonical_value(
      values[[index]], paste0(label, "[", index, "]")
    )
  }, character(1))
}

.ecoda_batch_canonical_token_pair <- function(key, canonical_value, key_index) {
  key_raw <- .ecoda_batch_utf8_raw(key, paste0("batch key ", key_index))
  value_raw <- .ecoda_batch_utf8_raw(
    canonical_value, paste0("canonical value for batch key ", key_index)
  )
  paste0(
    length(key_raw), ":", .ecoda_batch_hex(key_raw), ",",
    length(value_raw), ":", .ecoda_batch_hex(value_raw)
  )
}

.ecoda_batch_composite_token_canonical <- function(batch_keys, canonical_values) {
  keys <- if (is.character(batch_keys) && length(batch_keys) >= 2L) {
    # Internal callers already hold the normalized ordered character vector;
    # convert to a list solely to reuse strict key validation.
    .ecoda_batch_key_vector(as.list(unname(batch_keys)))
  } else {
    .ecoda_batch_key_vector(batch_keys)
  }
  if (length(keys) < 2L) {
    .ecoda_batch_stop("composite tokens require at least two batch keys")
  }
  if (!is.list(canonical_values) || length(canonical_values) != length(keys)) {
    .ecoda_batch_stop("canonical composite values must be one value per batch key")
  }
  pairs <- vapply(seq_along(keys), function(index) {
    value <- canonical_values[[index]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      .ecoda_batch_stop("canonical composite values must be nonmissing strings")
    }
    # Values supplied here have already been canonicalized.  Requiring one of
    # the fixed tags prevents accidental double-canonicalization.
    if (!grepl("^(s:|b:(true|false)$|i:-?[0-9]+$|f64:[0-9a-f]{16}$)", value, perl = TRUE)) {
      .ecoda_batch_stop("canonical composite value has an invalid type tag")
    }
    .ecoda_batch_canonical_token_pair(keys[[index]], value, index)
  }, character(1))
  paste0(.ecoda_batch_token_version, "|", length(keys), "|", paste(pairs, collapse = ";"))
}

# Encode one raw row in the versioned token format for two or more keys.
# Scalar/one-key corrected configurations remain direct and never call this
# encoder.
ecoda_batch_composite_token <- function(batch_keys, values) {
  keys <- ecoda_batch_normalize_keys(batch_keys)
  if (length(keys) < 2L) {
    .ecoda_batch_stop("composite tokens require at least two batch keys")
  }
  if (is.data.frame(values) || is.null(values) || length(values) != length(keys)) {
    .ecoda_batch_stop("composite values must contain one scalar per batch key")
  }
  if (is.list(values)) {
    value_names <- names(values)
    if (!is.null(value_names)) {
      if (length(value_names) != length(keys) || anyNA(value_names) ||
          any(!nzchar(value_names)) || anyDuplicated(value_names) ||
          !setequal(value_names, keys)) {
        .ecoda_batch_stop(
          "named composite values must contain exactly one value per batch key"
        )
      }
      cells <- unname(values[keys])
    } else {
      cells <- unname(values)
    }
  } else {
    value_names <- names(values)
    if (!is.null(value_names)) {
      if (length(value_names) != length(keys) || anyNA(value_names) ||
          any(!nzchar(value_names)) || anyDuplicated(value_names) ||
          !setequal(value_names, keys)) {
        .ecoda_batch_stop(
          "named composite values must contain exactly one value per batch key"
        )
      }
      cells <- lapply(keys, function(key) values[[key]])
    } else {
      cells <- lapply(seq_along(values), function(index) values[[index]])
    }
  }
  canonical <- lapply(seq_along(keys), function(index) {
    ecoda_batch_canonical_value(
      cells[[index]], paste0("batch key ", keys[[index]], " value")
    )
  })
  .ecoda_batch_composite_token_canonical(keys, canonical)
}

.ecoda_batch_metadata_names <- function(metadata) {
  if (!is.data.frame(metadata)) {
    .ecoda_batch_stop("cell metadata must be a data.frame")
  }
  columns <- colnames(metadata)
  if (is.null(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns)) {
    .ecoda_batch_stop("cell metadata must have unique nonblank column names")
  }
  columns <- vapply(seq_along(columns), function(index) {
    .ecoda_batch_utf8_text(columns[[index]], paste0("metadata column ", index))
  }, character(1))
  if (anyDuplicated(columns)) {
    .ecoda_batch_stop("cell metadata column names must be unique UTF-8 strings")
  }
  columns
}

.ecoda_batch_validate_fraction <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0 || value >= 1) {
    .ecoda_batch_stop("near_unique_fraction must be one finite value in (0, 1)")
  }
  as.numeric(value)
}

.ecoda_batch_sample_ids <- function(values) {
  if (is.list(values) && !is.factor(values)) {
    .ecoda_batch_stop("Sample metadata has an unsupported list/object value")
  }
  if (!(is.character(values) || is.factor(values))) {
    .ecoda_batch_stop("Sample metadata must be a character or factor column")
  }
  if (length(values) == 0L) {
    .ecoda_batch_stop("cell metadata must contain at least one cell")
  }
  vapply(seq_along(values), function(index) {
    ecoda_batch_canonical_value(
      values[[index]], paste0("Sample[", index, "]")
    )
  }, character(1))
}

.ecoda_batch_sample_labels <- function(values, first_indices) {
  labels <- values[first_indices]
  if (is.factor(labels)) labels <- as.character(labels)
  if (is.character(labels)) {
    labels <- vapply(seq_along(labels), function(index) {
      .ecoda_batch_utf8_text(labels[[index]], paste0("Sample label ", index))
    }, character(1))
  } else {
    labels <- as.character(labels)
  }
  unname(labels)
}

.ecoda_batch_sample_design <- function(sample_frame, label) {
  if (!is.data.frame(sample_frame) || nrow(sample_frame) == 0L) {
    .ecoda_batch_stop(label, " has no sample rows")
  }
  design <- tryCatch(
    stats::model.matrix(~ ., data = sample_frame),
    error = function(error) {
      .ecoda_batch_stop(label, " could not be constructed: ", conditionMessage(error))
    }
  )
  if (!is.numeric(design) || any(!is.finite(design))) {
    .ecoda_batch_stop(label, " contains non-finite design values")
  }
  rank <- qr(design, tol = 1e-10)$rank
  if (rank < ncol(design)) {
    .ecoda_batch_stop(label, " is rank deficient or perfectly confounded")
  }
  residual_df <- nrow(design) - rank
  if (residual_df < 1L) {
    .ecoda_batch_stop(label, " is non-estimable: no residual degrees of freedom")
  }
  list(
    rank = as.integer(rank),
    columns = as.integer(ncol(design)),
    residual_df = as.integer(residual_df)
  )
}

# Validate every cell before any sample-level collapse.  The returned
# sample_metadata uses the first row only after exact canonical equality within
# each Sample has been established; canonical_sample_metadata is the stable
# representation used by composite/fingerprint consumers.
ecoda_batch_validate_metadata <- function(
  metadata,
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50
) {
  columns <- .ecoda_batch_metadata_names(metadata)
  if (nrow(metadata) == 0L) {
    .ecoda_batch_stop("cell metadata must contain at least one cell")
  }
  sample_col <- .ecoda_batch_config_name(sample_col, "sample column")
  biological_label <- .ecoda_batch_config_name(
    biological_label, "biological label", allow_null = TRUE
  )
  keys <- ecoda_batch_normalize_keys(
    batch_keys,
    sample_col = sample_col,
    biological_label = biological_label
  )
  if (!sample_col %in% columns) {
    .ecoda_batch_stop("cell metadata is missing the standardized Sample column: ", sample_col)
  }
  if (!is.null(biological_label) && !biological_label %in% columns) {
    .ecoda_batch_stop("cell metadata is missing the biological label column: ", biological_label)
  }
  missing_keys <- keys[!keys %in% columns]
  if (length(missing_keys) > 0L) {
    .ecoda_batch_stop(
      "cell metadata is missing configured batch columns: ",
      paste(missing_keys, collapse = ", ")
    )
  }
  if (.ecoda_batch_reserved_name %in% columns) {
    .ecoda_batch_stop(
      "cell metadata already contains the reserved temporary column: ",
      .ecoda_batch_reserved_name
    )
  }
  near_unique_fraction <- .ecoda_batch_validate_fraction(near_unique_fraction)
  sample_values <- metadata[[sample_col]]
  sample_group_tokens <- .ecoda_batch_sample_ids(sample_values)
  sample_group_order <- unique(sample_group_tokens)
  sample_first_indices <- match(sample_group_order, sample_group_tokens)
  sample_ids <- .ecoda_batch_sample_labels(
    sample_values, sample_first_indices
  )
  sample_index <- match(sample_group_tokens, sample_group_order)
  n_samples <- length(sample_group_order)
  if (n_samples < 2L) {
    .ecoda_batch_stop("corrected batch metadata requires at least two Samples")
  }

  canonical_cells <- vector("list", length(keys))
  canonical_samples <- vector("list", length(keys))
  original_samples <- vector("list", length(keys))
  per_key_levels <- vector("list", length(keys))
  names(canonical_cells) <- keys
  names(canonical_samples) <- keys
  names(original_samples) <- keys
  names(per_key_levels) <- keys
  sample_constancy <- setNames(rep(TRUE, length(keys)), keys)
  key_level_counts <- setNames(integer(length(keys)), keys)
  key_near_unique_fraction <- setNames(numeric(length(keys)), keys)

  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    column <- metadata[[key]]
    if (is.list(column) && !is.factor(column)) {
      .ecoda_batch_stop("configured batch column is list-valued: ", key)
    }
    canonical <- vapply(seq_len(nrow(metadata)), function(cell_index) {
      cell <- column[[cell_index]]
      if (.ecoda_batch_is_missing_or_sentinel(cell)) {
        .ecoda_batch_stop(
          "configured batch column '", key,
          "' has missing, blank, or forbidden sentinel value at cell ", cell_index
        )
      }
      ecoda_batch_canonical_value(cell, paste0(key, "[", cell_index, "]"))
    }, character(1))
    canonical_cells[[key_index]] <- canonical

    collapsed <- character(n_samples)
    first_indices <- integer(n_samples)
    for (sample_index_value in seq_len(n_samples)) {
      cell_indices <- which(sample_index == sample_index_value)
      first_indices[[sample_index_value]] <- cell_indices[[1L]]
      values <- canonical[cell_indices]
      if (length(unique(values)) != 1L) {
        sample_constancy[[key]] <- FALSE
        .ecoda_batch_stop(
          "configured batch column '", key,
          "' disagrees within Sample '", sample_ids[[sample_index_value]], "'"
        )
      }
      collapsed[[sample_index_value]] <- values[[1L]]
    }
    canonical_samples[[key_index]] <- collapsed
    original_samples[[key_index]] <- column[first_indices]

    levels <- .ecoda_batch_unique_raw_sorted(collapsed, paste0(key, " levels"))
    if (length(levels) < 2L) {
      .ecoda_batch_stop(
        "configured batch column '", key, "' has fewer than two observed levels"
      )
    }
    key_level_counts[[key]] <- length(levels)
    key_near_unique_fraction[[key]] <- length(levels) / n_samples
    if (key_near_unique_fraction[[key]] > near_unique_fraction) {
      .ecoda_batch_stop(
        "configured batch column '", key, "' is near-unique: ",
        length(levels), "/", n_samples,
        " levels (threshold ", format(near_unique_fraction, trim = TRUE), ")"
      )
    }
    per_key_levels[[key_index]] <- levels
  }

  canonical_sample_metadata <- data.frame(
    Sample = sample_group_order,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  sample_metadata <- data.frame(
    Sample = sample_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    canonical_sample_metadata[[key]] <- canonical_samples[[key_index]]
    sample_metadata[[key]] <- original_samples[[key_index]]
  }
  rownames(canonical_sample_metadata) <- sample_group_order
  rownames(sample_metadata) <- sample_ids

  scalarization <- if (length(keys) >= 2L) "composite_v1" else "direct_v1"

  additive_frame <- canonical_sample_metadata[keys]
  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    additive_frame[[key]] <- factor(
      canonical_samples[[key_index]],
      levels = per_key_levels[[key_index]]
    )
  }
  key_design <- .ecoda_batch_sample_design(
    additive_frame, "additive batch design"
  )

  if (length(keys) >= 2L) {
    canonical_composite_values <- vapply(seq_len(n_samples), function(sample_index_value) {
      .ecoda_batch_composite_token_canonical(
        keys,
        lapply(canonical_samples, function(values) values[[sample_index_value]])
      )
    }, character(1))
    composite_values <- canonical_composite_values
    composite_levels <- .ecoda_batch_unique_raw_sorted(
      canonical_composite_values, "composite levels"
    )
    if (length(composite_levels) < 2L) {
      .ecoda_batch_stop("combined batch design has fewer than two observed levels")
    }
    composite_near_unique_fraction <- length(composite_levels) / n_samples
    if (composite_near_unique_fraction > near_unique_fraction) {
      .ecoda_batch_stop(
        "combined batch design is near-unique: ", length(composite_levels), "/",
        n_samples, " levels (threshold ", format(near_unique_fraction, trim = TRUE), ")"
      )
    }
    composite_frame <- data.frame(
      composite = factor(canonical_composite_values, levels = composite_levels),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    composite_design <- .ecoda_batch_sample_design(
      composite_frame, "combined batch design"
    )
  } else {
    canonical_composite_values <- canonical_samples[[1L]]
    composite_values <- canonical_composite_values
    # Direct scalarization does not create a composite factor.  Keep the
    # direct key's levels in per_key_levels and expose no composite levels.
    composite_levels <- character()
    composite_near_unique_fraction <- NA_real_
    composite_design <- key_design
  }

  list(
    valid = TRUE,
    token_version = .ecoda_batch_token_version,
    encoding = .ecoda_batch_token_version,
    encoding_version = .ecoda_batch_token_version,
    sample_col = sample_col,
    sample_column = sample_col,
    biological_label = biological_label,
    biological_column = biological_label,
    ordered_keys = keys,
    keys = keys,
    key_count = as.integer(length(keys)),
    scalarization = scalarization,
    sample_ids = sample_ids,
    sample_group_ids = sample_group_order,
    n_cells = as.integer(nrow(metadata)),
    n_obs = as.integer(nrow(metadata)),
    n_samples = as.integer(n_samples),
    sample_constancy = sample_constancy,
    sample_metadata = sample_metadata,
    canonical_sample_metadata = canonical_sample_metadata,
    canonical_cell_values = canonical_cells,
    canonical_values = canonical_cells,
    per_key_levels = per_key_levels,
    levels = per_key_levels,
    key_level_counts = key_level_counts,
    key_near_unique_fraction = key_near_unique_fraction,
    near_unique_fraction = as.numeric(near_unique_fraction),
    composite_values = composite_values,
    composite_levels = composite_levels,
    composite_level_count = as.integer(length(composite_levels)),
    composite_near_unique_fraction = as.numeric(composite_near_unique_fraction),
    estimable = TRUE,
    additive_design = key_design,
    design_rank = key_design$rank,
    design_columns = key_design$columns,
    composite_design = composite_design,
    composite_design_rank = composite_design$rank,
    composite_design_columns = composite_design$columns
  )
}

.ecoda_batch_ordered_key_encoding <- function(batch_keys) {
  keys <- if (is.character(batch_keys) && length(batch_keys) >= 2L) {
    .ecoda_batch_key_vector(as.list(unname(batch_keys)))
  } else {
    .ecoda_batch_key_vector(batch_keys)
  }
  encoded <- vapply(seq_along(keys), function(index) {
    raw <- .ecoda_batch_utf8_raw(keys[[index]], paste0("batch key ", index))
    paste0(length(raw), ":", .ecoda_batch_hex(raw), ";")
  }, character(1))
  paste0(length(keys), "|", paste(encoded, collapse = ""))
}

.ecoda_batch_field_encoding <- function(name, value) {
  name_raw <- .ecoda_batch_utf8_raw(name, "fingerprint field name")
  value_raw <- .ecoda_batch_utf8_raw(value, "fingerprint field value")
  paste0(
    length(name_raw), ":", .ecoda_batch_hex(name_raw), ",",
    length(value_raw), ":", .ecoda_batch_hex(value_raw), ";"
  )
}

.ecoda_batch_id <- function(value, label, allowed = NULL) {
  value <- .ecoda_batch_config_name(value, label)
  if (is.null(value)) .ecoda_batch_stop(label, " is required")
  if (!is.null(allowed) && !value %in% allowed) {
    .ecoda_batch_stop(label, " is not a recognized corrected-mode identifier")
  }
  value
}

.ecoda_batch_fingerprint_parts <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  keys <- ecoda_batch_normalize_keys(batch_keys)
  method_id <- .ecoda_batch_id(method_id, "method_id", .ecoda_batch_method_ids)
  model_id <- .ecoda_batch_id(model_id, "model_id", .ecoda_batch_model_ids)
  if (is.null(scalarization)) {
    scalarization <- if (length(keys) >= 2L) "composite_v1" else "direct_v1"
  }
  if (!is.character(scalarization) || length(scalarization) != 1L ||
      is.na(scalarization) ||
      !scalarization %in% c("direct_v1", "composite_v1")) {
    .ecoda_batch_stop("scalarization must be exactly direct_v1 or composite_v1")
  }
  expected <- if (length(keys) >= 2L) "composite_v1" else "direct_v1"
  if (!identical(scalarization, expected)) {
    .ecoda_batch_stop(
      "scalarization does not match the number of configured batch keys"
    )
  }
  list(
    keys = keys,
    method_id = method_id,
    model_id = model_id,
    scalarization = scalarization,
    ordered_key_encoding = .ecoda_batch_ordered_key_encoding(keys)
  )
}

# Return the exact text payload (including the embedded NUL after the prefix)
# used by ecoda_batch_fingerprint.  All fields after that NUL are ASCII hex.
ecoda_batch_fingerprint_payload <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  parts <- .ecoda_batch_fingerprint_parts(
    batch_keys, method_id, model_id, scalarization
  )
  fields <- paste0(
    .ecoda_batch_field_encoding("encoding", .ecoda_batch_token_version),
    .ecoda_batch_field_encoding("keys", parts$ordered_key_encoding),
    .ecoda_batch_field_encoding("scalarization", parts$scalarization),
    .ecoda_batch_field_encoding("method", parts$method_id),
    .ecoda_batch_field_encoding("model", parts$model_id)
  )
  paste0(
    "ecoda_batch_contract_v1",
    rawToChar(as.raw(0L)),
    fields
  )
}

.ecoda_batch_fingerprint_raw <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  parts <- .ecoda_batch_fingerprint_parts(
    batch_keys, method_id, model_id, scalarization
  )
  fields <- paste0(
    .ecoda_batch_field_encoding("encoding", .ecoda_batch_token_version),
    .ecoda_batch_field_encoding("keys", parts$ordered_key_encoding),
    .ecoda_batch_field_encoding("scalarization", parts$scalarization),
    .ecoda_batch_field_encoding("method", parts$method_id),
    .ecoda_batch_field_encoding("model", parts$model_id)
  )
  c(
    charToRaw(.ecoda_batch_contract_version),
    as.raw(0L),
    charToRaw(fields)
  )
}

# SHA-256 is taken over the exact bytes, not over R's serialized object form.
ecoda_batch_fingerprint <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  payload <- .ecoda_batch_fingerprint_raw(
    batch_keys, method_id, model_id, scalarization
  )
  if (!requireNamespace("digest", quietly = TRUE)) {
    .ecoda_batch_stop("the digest package is required for batch fingerprints")
  }

  fingerprint <- tolower(digest::digest(payload, algo = "sha256", serialize = FALSE))
  if (!is.character(fingerprint) || length(fingerprint) != 1L ||
      !grepl("^[0-9a-f]{64}$", fingerprint, perl = TRUE)) {
    .ecoda_batch_stop("batch fingerprint did not produce a lowercase SHA-256 digest")
  }
  fingerprint
}
# Lowercase hexadecimal transport form for the exact fingerprint bytes.
ecoda_batch_fingerprint_payload_hex <- function(
  batch_keys,
  method_id,
  model_id,
  scalarization = NULL
) {
  .ecoda_batch_hex(.ecoda_batch_fingerprint_raw(
    batch_keys, method_id, model_id, scalarization
  ))
}


# Build lightweight source/configuration identity without reading cell metadata.
# The returned named list is the strict, JSON-safe identity consumed by
# corrected artifact validators; it intentionally contains no cell/sample values.
ecoda_batch_contract_identity <- function(
  batch_keys,
  sample_col = "Sample",
  method_id = NULL,
  model_id = NULL
) {
  sample_col <- .ecoda_batch_config_name(sample_col, "sample column")
  if (identical(sample_col, .ecoda_batch_reserved_name)) {
    .ecoda_batch_stop("sample column uses the reserved temporary name")
  }
  keys <- ecoda_batch_normalize_keys(
    batch_keys,
    sample_col = sample_col
  )
  scalarization <- if (length(keys) >= 2L) {
    "composite_v1"
  } else {
    "direct_v1"
  }
  method_id <- .ecoda_batch_id(
    method_id, "method_id", .ecoda_batch_method_ids
  )
  model_id <- .ecoda_batch_id(
    model_id, "model_id", .ecoda_batch_model_ids
  )

  # Reproduce the existing fingerprint byte stream from normalized keys
  # directly.  This avoids materializing metadata and also permits a custom
  # sample column when it is not one of the source batch keys.
  fields <- paste0(
    .ecoda_batch_field_encoding("encoding", .ecoda_batch_token_version),
    .ecoda_batch_field_encoding(
      "keys", .ecoda_batch_ordered_key_encoding(keys)
    ),
    .ecoda_batch_field_encoding("scalarization", scalarization),
    .ecoda_batch_field_encoding("method", method_id),
    .ecoda_batch_field_encoding("model", model_id)
  )
  payload <- c(
    charToRaw(.ecoda_batch_contract_version),
    as.raw(0L),
    charToRaw(fields)
  )
  if (!requireNamespace("digest", quietly = TRUE)) {
    .ecoda_batch_stop("the digest package is required for batch fingerprints")
  }
  fingerprint <- tolower(
    digest::digest(payload, algo = "sha256", serialize = FALSE)
  )
  if (!is.character(fingerprint) || length(fingerprint) != 1L ||
      !grepl("^[0-9a-f]{64}$", fingerprint, perl = TRUE)) {
    .ecoda_batch_stop("batch fingerprint did not produce a lowercase SHA-256 digest")
  }

  list(
    contract_version = .ecoda_batch_contract_version,
    token_version = .ecoda_batch_token_version,
    ordered_source_keys = unname(keys),
    scalarization = scalarization,
    method_id = method_id,
    model_id = model_id,
    required_source_obs_columns = unname(c(sample_col, keys)),
    reserved_obs_name = .ecoda_batch_reserved_name,
    reserved_obs_absent = TRUE,
    fingerprint = fingerprint
  )
}
# Build the compact, validated summary carried by corrected artifact
# identities.  This boundary deliberately accepts only the full-cell
# validation result; it never derives values from a reduced sample table and
# never carries per-cell or per-sample vectors into the persisted identity.
ecoda_batch_build_validation_summary <- function(
  validation,
  correction_mode,
  correction_formula
) {
  if (!is.list(validation) || !isTRUE(validation[["valid"]])) {
    .ecoda_batch_stop(
      "compact corrected batch summary requires a validated full-cell result"
    )
  }
  keys <- validation[["ordered_keys"]]
  if (is.null(keys)) keys <- validation[["keys"]]
  if (!is.character(keys) || length(keys) == 0L ||
      anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    .ecoda_batch_stop(
      "compact corrected batch summary requires ordered validated keys"
    )
  }
  correction_mode <- .ecoda_batch_config_name(
    correction_mode, "correction mode"
  )
  correction_formula <- .ecoda_batch_config_name(
    correction_formula, "correction formula"
  )

  constancy <- validation[["sample_constancy"]]
  levels <- validation[["per_key_levels"]]
  if (is.null(levels)) levels <- validation[["levels"]]
  if (!is.list(levels) || is.null(names(levels)) ||
      !identical(names(levels), unname(keys))) {
    .ecoda_batch_stop(
      "compact corrected batch summary requires per-key levels in key order"
    )
  }
  if (!is.null(constancy)) {
    if (is.null(names(constancy)) ||
        !all(unname(keys) %in% names(constancy))) {
      .ecoda_batch_stop(
        "compact corrected batch summary requires per-key sample constancy"
      )
    }
    constancy_values <- vapply(unname(keys), function(key) {
      value <- constancy[[key]]
      is.logical(value) && length(value) == 1L && !is.na(value) &&
        isTRUE(value)
    }, logical(1))
    if (any(!constancy_values)) {
      .ecoda_batch_stop(
        "compact corrected batch summary requires constant configured keys"
      )
    }
  } else {
    .ecoda_batch_stop(
      "compact corrected batch summary requires per-key sample constancy"
    )
  }

  per_key_levels <- setNames(lapply(unname(keys), function(key) {
    values <- levels[[key]]
    if (!is.character(values) || anyNA(values) || any(!nzchar(values)) ||
        anyDuplicated(values)) {
      .ecoda_batch_stop(
        "compact corrected batch summary has invalid levels for key ", key
      )
    }
    # Validation already performs byte-order sorting.  Re-checking the order
    # here prevents a caller from attaching a hand-edited compact summary.
    expected <- .ecoda_batch_raw_sort(values, paste0(key, " levels"))
    if (!identical(unname(values), unname(expected))) {
      .ecoda_batch_stop(
        "compact corrected batch summary levels are not canonically sorted for key ",
        key
      )
    }
    unname(values)
  }), unname(keys))
  key_level_counts <- setNames(lapply(per_key_levels, function(values) {
    as.integer(length(values))
  }), unname(keys))

  composite_levels <- validation[["composite_levels"]]
  if (is.null(composite_levels)) composite_levels <- character()
  if (!is.character(composite_levels) || anyNA(composite_levels) ||
      any(!nzchar(composite_levels)) || anyDuplicated(composite_levels)) {
    .ecoda_batch_stop(
      "compact corrected batch summary has invalid composite levels"
    )
  }
  expected_composite_levels <- .ecoda_batch_raw_sort(
    composite_levels, "composite levels"
  )
  if (!identical(unname(composite_levels), unname(expected_composite_levels))) {
    .ecoda_batch_stop(
      "compact corrected batch summary composite levels are not canonically sorted"
    )
  }

  n_cells <- validation[["n_cells"]]
  if (is.null(n_cells)) n_cells <- validation[["n_obs"]]
  n_samples <- validation[["n_samples"]]
  if (!is.numeric(n_cells) || length(n_cells) != 1L || is.na(n_cells) ||
      !is.finite(n_cells) || n_cells < 1 || n_cells != floor(n_cells) ||
      !is.numeric(n_samples) || length(n_samples) != 1L || is.na(n_samples) ||
      !is.finite(n_samples) || n_samples < 1 || n_samples != floor(n_samples)) {
    .ecoda_batch_stop(
      "compact corrected batch summary requires positive integer cell/sample counts"
    )
  }
  if (length(composite_levels) > 0L && length(keys) < 2L) {
    .ecoda_batch_stop(
      "direct scalar corrected batch summary must have empty composite levels"
    )
  }
  if (length(composite_levels) == 0L && length(keys) >= 2L) {
    .ecoda_batch_stop(
      "composite corrected batch summary must include composite levels"
    )
  }

  list(
    schema_version = 1L,
    validated_before_reduction = TRUE,
    sample_constancy = setNames(
      lapply(unname(keys), function(key) TRUE),
      unname(keys)
    ),
    per_key_levels = per_key_levels,
    key_level_counts = key_level_counts,
    composite_levels = unname(composite_levels),
    composite_level_count = as.integer(length(composite_levels)),
    n_cells = as.integer(n_cells),
    n_samples = as.integer(n_samples),
    correction_mode = correction_mode,
    correction_formula = correction_formula
  )
}

# Drop fields that can encode cell/sample vectors before attaching the compact
# summary.  Full serializers remain available for synthetic/in-memory use;
# this helper is only for run-owned corrected artifact identities.
.ecoda_batch_compact_identity <- function(identity) {
  forbidden <- c(
    "composite_values", "sample_composite_values", "sample_ids",
    "sample_group_ids", "canonical_cell_values", "canonical_values",
    "sample_metadata", "canonical_sample_metadata", "tokens", "row_tokens",
    "scalarized_values"
  )
  identity[setdiff(names(identity), forbidden)]
}

# Attach a compact validated summary to a configuration/fingerprint identity.
# The source identity fields and fingerprint are retained byte-for-byte; only
# vector-bearing serializer fields are removed and validation_summary added.
ecoda_batch_augment_contract <- function(
  identity,
  validation,
  correction_mode,
  correction_formula
) {
  if (!is.list(identity) || is.null(names(identity)) ||
      anyNA(names(identity)) || any(!nzchar(names(identity)))) {
    .ecoda_batch_stop(
      "compact corrected batch identity must be a named list"
    )
  }
  summary <- ecoda_batch_build_validation_summary(
    validation = validation,
    correction_mode = correction_mode,
    correction_formula = correction_formula
  )
  identity <- .ecoda_batch_compact_identity(identity)
  identity[["validation_summary"]] <- summary
  identity
}

# Explicit aliases keep the independent R boundary discoverable alongside the
# Python ``build_batch_validation_summary``/``augment_batch_contract`` API.
ecoda_batch_validation_summary <- ecoda_batch_build_validation_summary
ecoda_batch_augment_contract_identity <- ecoda_batch_augment_contract

# Return the fixed correction policy for one corrected consumer.  The formula
# text is metadata, not an executable model formula; aliases are deliberately
# fixed so configured biological names cannot leak into the persisted identity.
ecoda_batch_correction_spec <- function(
  method_id,
  batch_keys,
  scalar_batch_col = NULL
) {
  batch_keys_input <- if (
    is.character(batch_keys) && length(batch_keys) > 1L
  ) {
    as.list(batch_keys)
  } else {
    batch_keys
  }
  keys <- ecoda_batch_normalize_keys(batch_keys_input)
  scalar <- if (is.null(scalar_batch_col)) {
    if (length(keys) >= 2L) .ecoda_batch_reserved_name else keys[[1L]]
  } else {
    .ecoda_batch_config_name(scalar_batch_col, "scalar batch column")
  }
  key_text <- paste(keys, collapse = ",")
  method_id <- .ecoda_batch_config_name(method_id, "method_id")
  if (identical(method_id, "preprocess")) {
    return(list(
      correction_mode = "preprocess_hvg_harmony",
      correction_formula = paste0(
        "HVG batch_key=", scalar, "; Harmony vars_use=[", key_text, "]"
      )
    ))
  }
  if (identical(method_id, "Pseudobulk")) {
    return(list(
      correction_mode = "batch_only_pseudobulk",
      correction_formula = paste0(
        "DESeq2 design=~ 1; limma removeBatchEffect(batch=", scalar, ")"
      )
    ))
  }
  if (method_id %in% c(
    "ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2"
  )) {
    aliases <- if (length(keys) == 1L) {
      "batch"
    } else {
      paste0("batch_key_", seq_along(keys))
    }
    return(list(
      correction_mode = "additive_random_intercepts",
      correction_formula = paste0(
        "y ~ 1 + ", paste0("(1 | ", aliases, ")", collapse = " + ")
      )
    ))
  }
  if (method_id %in% c("GloScope", "PILOT", "QOT")) {
    return(list(
      correction_mode = "native_harmony_embedding",
      correction_formula = "embedding=X_pca_harmony_batch_effect_corrected_hvg2000"
    ))
  }
  if (identical(method_id, "MrVI")) {
    return(list(
      correction_mode = "mrvi_composite_batch",
      correction_formula = paste0(
        "MRVI.setup_anndata(batch_key=", scalar, ")"
      )
    ))
  }
  .ecoda_batch_stop(
    "unsupported corrected batch correction spec method_id: ", method_id
  )
}



# Build a validated sample-level composite without modifying the input.  For a
# scalar/one-key configuration, composite_values is the canonical value of the
# direct key (sample_metadata retains its original direct column); for two or
# more keys it is the exact ecoda_batch_composite_v1 token.
ecoda_batch_build_composite <- function(
  metadata,
  batch_keys,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50
) {
  validation <- ecoda_batch_validate_metadata(
    metadata = metadata,
    batch_keys = batch_keys,
    sample_col = sample_col,
    biological_label = biological_label,
    near_unique_fraction = near_unique_fraction
  )
  keys <- validation$ordered_keys
  composite_name <- if (length(keys) >= 2L) {
    .ecoda_batch_reserved_name
  } else {
    keys[[1L]]
  }
  list(
    token_version = validation$token_version,
    encoding = validation$encoding,
    ordered_keys = keys,
    key_count = validation$key_count,
    scalarization = validation$scalarization,
    composite_name = composite_name,
    reserved_name = .ecoda_batch_reserved_name,
    sample_ids = validation$sample_ids,
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

# Serialize all semantic identity needed by a run-owned corrected artifact.
# This is an ordinary ordered R list (not JSON or R's binary serialization); its
# fingerprint payload is independently available for cross-language fixtures.
ecoda_batch_serialize_metadata <- function(
  metadata,
  batch_keys,
  method_id,
  model_id,
  sample_col = "Sample",
  biological_label = NULL,
  near_unique_fraction = 0.50
) {
  composite <- ecoda_batch_build_composite(
    metadata = metadata,
    batch_keys = batch_keys,
    sample_col = sample_col,
    biological_label = biological_label,
    near_unique_fraction = near_unique_fraction
  )
  validation <- composite$validation
  fingerprint <- ecoda_batch_fingerprint(
    batch_keys = composite$ordered_keys,
    method_id = method_id,
    model_id = model_id,
    scalarization = composite$scalarization
  )
  fingerprint_payload <- ecoda_batch_fingerprint_payload(
    batch_keys = composite$ordered_keys,
    method_id = method_id,
    model_id = model_id,
    scalarization = composite$scalarization
  )
  fingerprint_payload_hex <- ecoda_batch_fingerprint_payload_hex(
    batch_keys = composite$ordered_keys,
    method_id = method_id,
    model_id = model_id,
    scalarization = composite$scalarization
  )
  list(
    contract_version = .ecoda_batch_contract_version,
    token_version = composite$token_version,
    encoding = composite$encoding,
    encoding_version = composite$token_version,
    ordered_keys = composite$ordered_keys,
    keys = composite$ordered_keys,
    key_count = composite$key_count,
    scalarization = composite$scalarization,
    per_key_levels = composite$per_key_levels,
    levels = composite$levels,
    composite_levels = composite$composite_levels,
    composite_level_count = composite$composite_level_count,
    composite_values = composite$composite_values,
    method = .ecoda_batch_id(method_id, "method_id", .ecoda_batch_method_ids),
    method_id = .ecoda_batch_id(method_id, "method_id", .ecoda_batch_method_ids),
    model = .ecoda_batch_id(model_id, "model_id", .ecoda_batch_model_ids),
    model_id = .ecoda_batch_id(model_id, "model_id", .ecoda_batch_model_ids),
    fingerprint_payload = fingerprint_payload,
    fingerprint_payload_hex = fingerprint_payload_hex,
    fingerprint = fingerprint,
    sample_column = validation$sample_col,
    sample_col = validation$sample_col,
    biological_column = validation$biological_label,
    biological_label = validation$biological_label,
    required_source_obs_columns = c(validation$sample_col, composite$ordered_keys),
    reserved_obs_name = composite$reserved_name,
    reserved_name = composite$reserved_name,
    reserved_obs_absent = TRUE,
    n_obs = validation$n_cells,
    n_cells = validation$n_cells,
    n_samples = validation$n_samples,
    sample_ids = validation$sample_ids,
    sample_group_ids = validation$sample_group_ids,
    sample_constancy = validation$sample_constancy,
    estimable = TRUE,
    near_unique_fraction = as.numeric(near_unique_fraction),
    key_level_counts = validation$key_level_counts,
    key_near_unique_fraction = validation$key_near_unique_fraction,
    composite_near_unique_fraction = validation$composite_near_unique_fraction,
    composite_name = composite$composite_name,
    design_rank = validation$additive_design$rank,
    design_columns = validation$additive_design$columns,
    composite_design_rank = validation$composite_design$rank,
    composite_design_columns = validation$composite_design$columns,
    additive_design = validation$additive_design,
    composite_design = validation$composite_design
  )
}
