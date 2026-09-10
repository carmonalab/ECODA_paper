# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Apply method label recoding using centralized lookup
apply_method_labels <- function(data, label_map = method_label_map_main) {
  data %>% mutate(method = recode(method, !!!label_map))
}

# Execution timing contract.  A schema-2 pseudobulk method row carries only
# variant-local work in time_secs; separate global shared rows carry aggregate
# plus shared-fit work once per timing identity.
# Rows without timing_schema == 2 retain their historical inclusive timing.
.execution_shared_timing_method <- "prepare_pseudobulk_shared"
.execution_shared_timing_methods <- c(
  "prepare_pseudobulk_shared",
  "prepare_pseudobulk_ct_shared"
)
.execution_is_shared_timing_method <- function(method) {
  method <- as.character(method)
  !is.na(method) & (
    method %in% .execution_shared_timing_methods |
      grepl(
        "^prepare_pseudobulk_ct_shared_[^[:space:]]+$",
        method,
        perl = TRUE
      )
  )
}
.execution_shared_method_scalar <- function(
  value, label = "shared_timing_method"
) {
  value <- .execution_id_scalar(value, label)
  if (!.execution_is_shared_timing_method(value)) {
    stop(label, " is not a recognized shared timing method.")
  }
  value
}
.execution_timing_columns <- c(
  "aggregate_time_secs",
  "shared_fit_time_secs",
  "shared_time_secs",
  "variant_time_secs",
  "shared_mem_GB",
  "timing_id",
  "timing_schema"
)

.execution_scalar <- function(value, label = "timing value", allow_na = TRUE) {
  if (is.null(value) || length(value) == 0L) {
    if (allow_na) return(NA_real_)
    stop(label, " is missing.")
  }
  if (length(value) != 1L) {
    stop(label, " must be one finite nonnegative number.")
  }
  if (is.list(value)) {
    value <- value[[1L]]
    if (is.null(value) || length(value) != 1L) {
      stop(label, " must be one finite nonnegative number.")
    }
  }
  if (is.factor(value)) value <- as.character(value)
  missing_input <- is.na(value)
  nan_input <- is.nan(value)
  if (inherits(value, "difftime")) {
    value <- suppressWarnings(as.numeric(value, units = "secs"))
  } else {
    value <- suppressWarnings(as.numeric(value))
  }
  if (length(value) != 1L) {
    stop(label, " must be one finite nonnegative number.")
  }
  if (is.na(value)) {
    if (allow_na && isTRUE(missing_input) && !isTRUE(nan_input)) {
      return(NA_real_)
    }
    stop(label, " must be one finite nonnegative number.")
  }
  if (!is.finite(value) || value < 0) {
    stop(label, " must be one finite nonnegative number.")
  }
  as.numeric(value)
}

.execution_cell <- function(column, index) {
  if (is.list(column)) column[[index]] else column[index]
}

.execution_id_scalar <- function(value, label = "timing_id") {
  if (is.null(value) || length(value) == 0L) {
    stop(label, " is missing.")
  }
  if (length(value) != 1L) stop(label, " must be one nonblank string.")
  if (is.list(value)) {
    value <- value[[1L]]
    if (is.null(value) || length(value) != 1L) {
      stop(label, " must be one nonblank string.")
    }
  }
  if (is.factor(value)) value <- as.character(value)
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop(label, " must be one nonblank string.")
  }
  value
}

.execution_schema2_mask <- function(data) {
  if (!"timing_schema" %in% names(data)) {
    return(rep(FALSE, nrow(data)))
  }
  raw_schema <- data[["timing_schema"]]
  if (length(raw_schema) != nrow(data)) {
    stop("Execution timing schema must have one value per row.")
  }
  schema <- rep(NA_real_, nrow(data))
  if (nrow(data) > 0L) {
    for (index in seq_len(nrow(data))) {
      cell <- .execution_cell(raw_schema, index)
      if (is.null(cell) || length(cell) == 0L) next
      schema[[index]] <- .execution_scalar(
        cell, "timing_schema", allow_na = TRUE
      )
    }
  }
  invalid <- !is.na(schema) & (
    !is.finite(schema) | schema != floor(schema) | schema != 2
  )
  if (any(invalid)) {
    stop("Unsupported execution timing schema; only schema 2 is recognized.")
  }
  !is.na(schema) & schema == 2
}

.execution_nonblank <- function(value) {
  if (length(value) == 0L) return(logical())
  vapply(seq_along(value), function(index) {
    cell <- .execution_cell(value, index)
    if (is.null(cell) || length(cell) != 1L) return(FALSE)
    if (is.factor(cell)) cell <- as.character(cell)
    cell <- as.character(cell)
    length(cell) == 1L && !is.na(cell) && nzchar(trimws(cell))
  }, logical(1))
}

.execution_numeric_column <- function(
  data, column, rows, allow_na = FALSE, label = column
) {
  if (!column %in% names(data)) {
    stop("Execution timing rows are missing required column: ", column)
  }
  values <- rep(NA_real_, nrow(data))
  for (index in which(rows)) {
    values[[index]] <- .execution_scalar(
      .execution_cell(data[[column]], index),
      label,
      allow_na = allow_na
    )
  }
  values
}

.execution_validate_schema2 <- function(data, schema2) {
  if (!any(schema2)) return(invisible(NULL))
  method_values <- as.character(data[["method"]])
  shared_rows <- .execution_is_shared_timing_method(method_values)
  local_rows <- schema2 & !shared_rows
  required <- c(
    "timing_id",
    "shared_time_secs",
    "shared_mem_GB"
  )
  if (any(local_rows)) required <- c(required, "variant_time_secs")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      "Schema-2 execution rows are missing required timing fields: ",
      paste(missing, collapse = ", ")
    )
  }
  time <- .execution_numeric_column(
    data, "time_secs", schema2, allow_na = FALSE, label = "time_secs"
  )
  mem <- .execution_numeric_column(
    data, "mem_GB", schema2, allow_na = TRUE, label = "mem_GB"
  )
  shared_time <- .execution_numeric_column(
    data, "shared_time_secs", schema2, allow_na = FALSE,
    label = "shared_time_secs"
  )
  shared_mem <- .execution_numeric_column(
    data, "shared_mem_GB", schema2, allow_na = TRUE,
    label = "shared_mem_GB"
  )
  variant <- if ("variant_time_secs" %in% names(data)) {
    .execution_numeric_column(
      data, "variant_time_secs", schema2, allow_na = TRUE,
      label = "variant_time_secs"
    )
  } else {
    rep(NA_real_, nrow(data))
  }
  ids <- rep(NA_character_, nrow(data))
  for (index in which(schema2)) {
    ids[[index]] <- .execution_id_scalar(
      .execution_cell(data[["timing_id"]], index)
    )
  }
  if ("shared_timing_method" %in% names(data)) {
    for (index in which(schema2)) {
      cell <- .execution_cell(data[["shared_timing_method"]], index)
      if (is.null(cell) || length(cell) == 0L ||
          (length(cell) == 1L && is.na(cell))) {
        next
      }
      .execution_shared_method_scalar(
        cell, "shared_timing_method"
      )
    }
  }
  if (any(local_rows & is.na(variant))) {
    stop("Schema-2 execution rows have invalid variant-local timing.")
  }
  if (any(schema2 & shared_rows & !is.na(variant))) {
    stop(
      "Schema-2 shared execution rows must not carry variant-local timing."
    )
  }
  if (any(schema2 & shared_rows &
          abs(time - shared_time) > 1e-9)) {
    stop(
      "Schema-2 shared execution rows must use shared_time_secs in time_secs."
    )
  }
  if (any(local_rows & abs(time - variant) > 1e-9)) {
    stop(
      "Schema-2 execution rows must use variant_time_secs in time_secs."
    )
  }

  decomposition <- c("aggregate_time_secs", "shared_fit_time_secs")
  present_decomposition <- intersect(decomposition, names(data))
  if (length(present_decomposition) > 0L &&
      length(present_decomposition) != length(decomposition)) {
    stop("Schema-2 execution rows have an incomplete shared decomposition.")
  }
  if (length(present_decomposition) == length(decomposition)) {
    aggregate <- .execution_numeric_column(
      data, "aggregate_time_secs", schema2, allow_na = TRUE,
      label = "aggregate_time_secs"
    )
    shared_fit <- .execution_numeric_column(
      data, "shared_fit_time_secs", schema2, allow_na = TRUE,
      label = "shared_fit_time_secs"
    )
    incomplete <- schema2 & xor(is.na(aggregate), is.na(shared_fit))
    if (any(incomplete)) {
      stop("Schema-2 execution rows have an incomplete shared decomposition.")
    }
    complete <- schema2 & !is.na(aggregate) & !is.na(shared_fit)
    if (any(complete & abs(shared_time - aggregate - shared_fit) > 1e-9)) {
      stop(
        "Schema-2 shared_time_secs must equal aggregate_time_secs + ",
        "shared_fit_time_secs."
      )
    }
  }

  identity_key <- paste(
    as.character(data[["dataset"]]), ids, sep = "\r"
  )
  for (key in unique(identity_key[schema2])) {
    rows <- schema2 & identity_key == key
    values <- shared_time[rows]
    if (length(values) > 1L &&
        (max(values) - min(values)) > 1e-9) {
      stop(
        "Schema-2 rows with one timing_id have inconsistent shared timing."
      )
    }
  }
  invisible(list(
    time = time,
    mem = mem,
    shared_time = shared_time,
    shared_mem = shared_mem,
    variant = variant,
    ids = ids
  ))
}


# Normalize execution rows without inventing a decomposition for old records.
# Timing schema 2 is the sole marker for the extended timing contract.  A
# missing marker keeps a row legacy-inclusive, even when bookkeeping fields
# such as timing_id happen to be present.
normalize_exec_times <- function(exec_times) {
  if (!is.data.frame(exec_times)) {
    stop("exec_times must be a data.frame.")
  }
  required <- c("dataset", "method", "time_secs", "mem_GB")
  missing <- setdiff(required, names(exec_times))
  if (length(missing) > 0L) {
    stop("Execution times are missing required columns: ",
         paste(missing, collapse = ", "))
  }

  data <- exec_times
  data[["time_secs"]] <- .execution_numeric_column(
    data,
    "time_secs",
    rep(TRUE, nrow(data)),
    allow_na = FALSE,
    label = "time_secs"
  )
  data[["mem_GB"]] <- .execution_numeric_column(
    data,
    "mem_GB",
    rep(TRUE, nrow(data)),
    allow_na = TRUE,
    label = "mem_GB"
  )
  schema2 <- .execution_schema2_mask(data)
  validation <- .execution_validate_schema2(data, schema2)
  if (any(schema2)) {
    local_rows <- schema2 &
      !.execution_is_shared_timing_method(as.character(data[["method"]]))
    if (any(local_rows)) {
      local <- validation$variant
      data[["time_secs"]][local_rows] <- local[local_rows]
    }
  }
  data
}

# Keep ordinary rows unique by dataset/method.  If an implementation publishes
# timing_id in an extended log, shared rows are unique by dataset/timing_id so
# distinct timing contexts cannot be charged twice or silently collapsed.
deduplicate_exec_times <- function(exec_times, keep = c("last", "first")) {
  keep <- match.arg(keep)
  data <- normalize_exec_times(exec_times)
  if (nrow(data) == 0L) return(data)
  key <- paste(as.character(data[["dataset"]]),
               as.character(data[["method"]]), sep = "\r")
  shared <- .execution_is_shared_timing_method(
    as.character(data[["method"]])
  )
  if ("timing_id" %in% names(data)) {
    identified <- shared & .execution_nonblank(data[["timing_id"]])
    if (any(identified)) {
      key[identified] <- paste(
        key[identified],
        as.character(data[["timing_id"]][identified]),
        sep = "\r"
      )
    }
  }
  data[!duplicated(key, fromLast = identical(keep, "last")), , drop = FALSE]
}

# Convert one in-memory result bundle to canonical report rows.  This keeps
# schema-2 variant timing local and emits the explicit shared row once per
# timing identity; legacy bundles continue to use inclusive exec_time.
execution_time_rows_from_bundle <- function(dataset, method, bundle) {
  if (!is.list(bundle)) stop("Benchmark timing bundle must be a list.")

  # Only an explicit, valid timing_schema marker enables schema 2.  In
  # particular, timing_id by itself is not evidence of a decomposition.
  schema2 <- FALSE
  if ("timing_schema" %in% names(bundle)) {
    schema_value <- .execution_scalar(
      bundle[["timing_schema"]], "timing_schema", allow_na = TRUE
    )
    schema2 <- !is.na(schema_value) && schema_value == 2
    if (!is.na(schema_value) && !schema2) {
      stop("Unsupported execution timing schema in result bundle.")
    }
  }

  method_value <- as.character(method)
  shared_bundle <- .execution_is_shared_timing_method(method_value)
  for (field in intersect(
    c("exec_time", "time_secs", "variant_time_secs"), names(bundle)
  )) {
    .execution_scalar(
      bundle[[field]],
      field,
      allow_na = identical(field, "variant_time_secs")
    )
  }

  if (schema2) {
    required_timing <- c(
      "shared_time_secs",
      "shared_mem_GB",
      "timing_id"
    )
    if (!shared_bundle) {
      required_timing <- c(required_timing, "variant_time_secs")
    }
    missing_timing <- setdiff(required_timing, names(bundle))
    if (length(missing_timing) > 0L) {
      stop("Schema-2 timing bundle is missing: ",
           paste(missing_timing, collapse = ", "))
    }
  }

  timing_fields <- c("exec_time", "time_secs", "variant_time_secs")
  if (schema2 && shared_bundle) {
    timing_fields <- c(timing_fields, "shared_time_secs")
  }
  has_time <- any(timing_fields %in% names(bundle))
  if (!has_time) return(NULL)
  raw_time <- if ("exec_time" %in% names(bundle)) {
    .execution_scalar(bundle[["exec_time"]], "exec_time", allow_na = FALSE)
  } else if ("time_secs" %in% names(bundle)) {
    .execution_scalar(bundle[["time_secs"]], "time_secs", allow_na = FALSE)
  } else if (schema2 && shared_bundle) {
    .execution_scalar(
      bundle[["shared_time_secs"]], "shared_time_secs", allow_na = FALSE
    )
  } else {
    .execution_scalar(
      bundle[["variant_time_secs"]], "variant_time_secs", allow_na = FALSE
    )
  }
  local_time <- if (schema2 && !shared_bundle) {
    .execution_scalar(
      bundle[["variant_time_secs"]], "variant_time_secs", allow_na = FALSE
    )
  } else {
    raw_time
  }
  if (schema2 && !shared_bundle) {
    for (field in intersect(c("exec_time", "time_secs"), names(bundle))) {
      recorded <- .execution_scalar(
        bundle[[field]], field, allow_na = FALSE
      )
      if (abs(recorded - local_time) > 1e-9) {
        stop("Schema-2 timing fields must agree with variant_time_secs.")
      }
    }
  }
  mem <- if ("mem_GB" %in% names(bundle)) {
    .execution_scalar(bundle[["mem_GB"]], "mem_GB", allow_na = TRUE)
  } else {
    NA_real_
  }
  row <- data.frame(
    dataset = as.character(dataset),
    method = as.character(method),
    time_secs = local_time,
    mem_GB = mem,
    stringsAsFactors = FALSE
  )
  if (!schema2) return(row)

  shared_method <- if ("shared_timing_method" %in% names(bundle)) {
    .execution_shared_method_scalar(
      bundle[["shared_timing_method"]], "shared_timing_method"
    )
  } else if (shared_bundle) {
    method_value
  } else if (grepl("^Pseudobulk_CT_", method_value)) {
    "prepare_pseudobulk_ct_shared"
  } else {
    .execution_shared_timing_method
  }

  timing_id <- .execution_id_scalar(
    bundle[["timing_id"]], "timing_id"
  )
  shared_time <- .execution_scalar(
    bundle[["shared_time_secs"]], "shared_time_secs", allow_na = FALSE
  )
  shared_mem <- .execution_scalar(
    bundle[["shared_mem_GB"]], "shared_mem_GB", allow_na = TRUE
  )
  if (schema2 && shared_bundle) {
    for (field in intersect(c("exec_time", "time_secs"), names(bundle))) {
      recorded <- .execution_scalar(
        bundle[[field]], field, allow_na = FALSE
      )
      if (abs(recorded - shared_time) > 1e-9) {
        stop(
          "Schema-2 shared timing fields must agree with shared_time_secs."
        )
      }
    }
  }
  if (shared_bundle && "variant_time_secs" %in% names(bundle) &&
      !is.na(.execution_scalar(
        bundle[["variant_time_secs"]], "variant_time_secs", allow_na = TRUE
      ))) {
    stop(
      "Schema-2 shared timing bundles must not carry variant-local timing."
    )
  }
  if ("aggregate_time_secs" %in% names(bundle) ||
      "shared_fit_time_secs" %in% names(bundle)) {
    if (!all(c("aggregate_time_secs", "shared_fit_time_secs") %in%
             names(bundle))) {
      stop("Schema-2 timing bundle has an incomplete shared decomposition.")
    }
    aggregate_time <- .execution_scalar(
      bundle[["aggregate_time_secs"]], "aggregate_time_secs",
      allow_na = FALSE
    )
    shared_fit_time <- .execution_scalar(
      bundle[["shared_fit_time_secs"]], "shared_fit_time_secs",
      allow_na = FALSE
    )
    if (abs(shared_time - aggregate_time - shared_fit_time) > 1e-9) {
      stop("Schema-2 shared_time_secs must equal aggregate_time_secs + ",
           "shared_fit_time_secs.")
    }
  }

  for (column in .execution_timing_columns) {
    if (column %in% names(bundle)) {
      row[[column]] <- if (column == "timing_id") {
        timing_id
      } else if (column == "timing_schema") {
        2
      } else {
        .execution_scalar(bundle[[column]], column, allow_na = TRUE)
      }
    }
  }
  row[["timing_id"]] <- timing_id
  row[["timing_schema"]] <- 2
  row[["shared_time_secs"]] <- shared_time
  row[["shared_mem_GB"]] <- shared_mem
  row[["variant_time_secs"]] <- local_time

  shared_row <- row
  shared_row[["method"]] <- shared_method
  shared_row[["time_secs"]] <- shared_time
  shared_row[["mem_GB"]] <- shared_mem
  shared_row[["variant_time_secs"]] <- NA_real_
  if (.execution_is_shared_timing_method(method_value)) {
    return(shared_row)
  }
  dplyr::bind_rows(row, shared_row)
}

# Summarize rows with one shared charge per timing identity.  Rows carrying
# timing_schema == 2 are classified independently from legacy rows in the
# same dataset.  A missing schema marker keeps a row's historical inclusive
# timing, while the dedicated shared method remains one global charge.
summarize_exec_times <- function(exec_times) {
  data <- deduplicate_exec_times(exec_times)
  if (nrow(data) == 0L) {
    return(data.frame(
      dataset = character(),
      variant_time_secs = numeric(),
      shared_time_secs = numeric(),
      legacy_inclusive_time_secs = numeric(),
      total_time_secs = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  schema2 <- .execution_schema2_mask(data)
  method_values <- as.character(data[["method"]])
  shared <- .execution_is_shared_timing_method(method_values)

  # Timing-extended logs may carry the shared fields on every local row but
  # omit the dedicated row.  Reconstitute one shared row for each
  # (dataset, timing_id), not merely one row per dataset.
  if (any(schema2 & !shared)) {
    if (!all(c("timing_id", "shared_time_secs") %in% names(data))) {
      stop(
        "Schema-2 execution rows are missing timing_id/shared_time_secs."
      )
    }
    timing_ids <- .execution_nonblank(data[["timing_id"]])
    shared_values <- .execution_numeric_column(
      data,
      "shared_time_secs",
      schema2,
      allow_na = FALSE,
      label = "shared_time_secs"
    )
    derive <- schema2 & !shared & timing_ids
    if (any(derive)) {
      derived <- data[derive, , drop = FALSE]
      derived_key <- paste(
        as.character(derived[["dataset"]]),
        as.character(derived[["timing_id"]]),
        sep = "\r"
      )
      derived <- derived[!duplicated(derived_key), , drop = FALSE]
      explicit_key <- character()
      explicit <- shared & .execution_nonblank(data[["timing_id"]])
      if (any(explicit)) {
        explicit_key <- paste(
          as.character(data[["dataset"]][explicit]),
          as.character(data[["timing_id"]][explicit]),
          sep = "\r"
        )
      }
      derived <- derived[!derived_key %in% explicit_key, ,
                         drop = FALSE]
      if (nrow(derived) > 0L) {
        derived_method <- if ("shared_timing_method" %in% names(derived)) {
          vapply(seq_len(nrow(derived)), function(index) {
            cell <- .execution_cell(
              derived[["shared_timing_method"]], index
            )
            if (is.null(cell) || length(cell) == 0L ||
                (length(cell) == 1L && is.na(cell))) {
              return(NA_character_)
            }
            .execution_shared_method_scalar(
              cell, "shared_timing_method"
            )
          }, character(1))
        } else {
          rep(NA_character_, nrow(derived))
        }
        fallback_method <- ifelse(
          grepl("^Pseudobulk_CT_", as.character(derived[["method"]])),
          "prepare_pseudobulk_ct_shared",
          .execution_shared_timing_method
        )
        derived_method[is.na(derived_method)] <-
          fallback_method[is.na(derived_method)]
        derived[["method"]] <- derived_method
        derived[["time_secs"]] <- .execution_numeric_column(
          derived,
          "shared_time_secs",
          rep(TRUE, nrow(derived)),
          allow_na = FALSE,
          label = "shared_time_secs"
        )
        derived[["mem_GB"]] <- .execution_numeric_column(
          derived,
          "shared_mem_GB",
          rep(TRUE, nrow(derived)),
          allow_na = TRUE,
          label = "shared_mem_GB"
        )
        if ("variant_time_secs" %in% names(derived)) {
          derived[["variant_time_secs"]] <- NA_real_
        }
        data <- deduplicate_exec_times(
          dplyr::bind_rows(data, derived)
        )
      }
    }
  }

  schema2 <- .execution_schema2_mask(data)
  method_values <- as.character(data[["method"]])
  shared <- .execution_is_shared_timing_method(method_values)
  dataset_values <- as.character(data[["dataset"]])
  time_values <- as.numeric(data[["time_secs"]])
  datasets <- unique(dataset_values)
  # A missing schema marker is legacy-inclusive, including legacy
  # prepare_pseudobulk_* rows.  Only explicitly marked schema-2 local rows
  # enter the variant decomposition; the method name alone identifies the
  # global shared row.
  variant_rows <- schema2 & !shared
  legacy_rows <- !schema2 & !shared

  sum_by_dataset <- function(values) {
    result <- setNames(numeric(length(datasets)), datasets)
    totals <- tapply(values, dataset_values, sum, na.rm = TRUE)
    if (length(totals) > 0L) {
      result[names(totals)] <- as.numeric(totals)
    }
    result
  }

  variant_by_ds <- sum_by_dataset(ifelse(variant_rows, time_values, 0))
  legacy_by_ds <- sum_by_dataset(ifelse(legacy_rows, time_values, 0))

  shared_values <- time_values
  schema_shared <- schema2 & shared
  if (any(schema_shared)) {
    shared_values[schema_shared] <- .execution_numeric_column(
      data,
      "shared_time_secs",
      schema_shared,
      allow_na = FALSE,
      label = "shared_time_secs"
    )[schema_shared]
  }
  shared_by_ds <- sum_by_dataset(ifelse(shared, shared_values, 0))

  out <- data.frame(
    dataset = datasets,
    variant_time_secs = as.numeric(variant_by_ds[datasets]),
    shared_time_secs = as.numeric(shared_by_ds[datasets]),
    legacy_inclusive_time_secs = as.numeric(legacy_by_ds[datasets]),
    stringsAsFactors = FALSE
  )
  out[is.na(out)] <- 0
  out$total_time_secs <- out$variant_time_secs +
    out$shared_time_secs + out$legacy_inclusive_time_secs
  out
}

# Merge benchmark results with execution times cleanly.  Shared preparation
# rows are global accounting rows and must never join to every variant.
# `time_secs` remains variant-local for schema 2 and inclusive for legacy data.
merge_exec_times <- function(df_results, exec_times) {
  method_times <- deduplicate_exec_times(exec_times)
  method_times <- method_times[
    !.execution_is_shared_timing_method(
      as.character(method_times[["method"]])
    ),
    ,
    drop = FALSE
  ]
  df_results %>%
    mutate(temp_key = paste0(dataset, method)) %>%
    left_join(
      method_times %>%
        mutate(temp_key = paste0(dataset, method)) %>%
        dplyr::select(temp_key, time_secs),
      by = "temp_key"
    ) %>%
    dplyr::select(-temp_key)
}

# For ggplot x-axis label bolding (and optionally color)
highlight <- function(x, pat, color = "black", family = "") {
  ifelse(
    grepl(pat, x),
    glue::glue("<b style='font-family:{family}; color:{color}'>{x}</b>"),
    x
  )
}
