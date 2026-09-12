#!/usr/bin/env Rscript
# Read-only corrected Stage 3 source metadata audit.
#
# This boundary deliberately reads an RDS only far enough to obtain Seurat
# @meta.data (or a data.frame).  It never touches assays, expression matrices,
# count layers, or H5AD files.  The resulting report is written atomically
# below the caller's run root.

options(stringsAsFactors = FALSE)

.stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

`%||%` <- function(value, fallback) {
  if (is.null(value)) fallback else value
}

.safe_text <- function(value, label, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    .stop(label, " must be one string")
  }
  value <- enc2utf8(unname(value))
  if (!isTRUE(validUTF8(value)) || grepl("[\\r\\n\\t]", value, perl = TRUE)) {
    .stop(label, " contains invalid record-delimiter characters")
  }
  if (!allow_empty && (!nzchar(value) || !identical(trimws(value), value))) {
    .stop(label, " must be nonblank and free of surrounding whitespace")
  }
  value
}

.safe_absolute <- function(value, label) {
  value <- .safe_text(value, label)
  if (!grepl("^/", value)) .stop(label, " must be absolute")
  value
}

.no_symlink <- function(path, label) {
  link <- tryCatch(Sys.readlink(path), error = function(error) "")
  if (length(link) > 0L && nzchar(link[[1L]])) {
    .stop(label, " must not be a symlink: ", path)
  }
  invisible(path)
}

.existing_file <- function(path, label) {
  path <- .safe_absolute(path, label)
  info <- file.info(path)
  if (!file.exists(path) || is.na(info$isdir) || isTRUE(info$isdir)) {
    .stop(label, " is not an existing regular file: ", path)
  }
  if (is.na(info$size) || info$size <= 0) {
    .stop(label, " is empty: ", path)
  }
  .no_symlink(path, label)
  normalizePath(path, mustWork = TRUE)
}

.existing_directory <- function(path, label) {
  path <- .safe_absolute(path, label)
  if (!dir.exists(path)) .stop(label, " is not an existing directory: ", path)
  .no_symlink(path, label)
  normalizePath(path, mustWork = TRUE)
}

.path_below <- function(path, root) {
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

.run_owned_output <- function(path, run_root) {
  path <- .safe_absolute(path, "output path")
  run_root <- .existing_directory(run_root, "run root")
  parent <- dirname(path)
  if (!dir.exists(parent)) .stop("output parent directory does not exist: ", parent)
  parent <- normalizePath(parent, mustWork = TRUE)
  if (!.path_below(parent, run_root)) {
    .stop("output path escapes the run root: ", path)
  }
  if (file.exists(path)) .no_symlink(path, "output path")
  path
}

.read_kv_manifest <- function(path, label) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) == 0L) .stop(label, " is empty")
  result <- list()
  for (index in seq_along(lines)) {
    line <- lines[[index]]
    equals <- regexpr("=", line, fixed = TRUE)[[1L]]
    if (equals <= 1L) .stop(label, " has malformed line ", index)
    key <- substr(line, 1L, equals - 1L)
    value <- substr(line, equals + 1L, nchar(line))
    .safe_text(key, paste0(label, " key ", index))
    if (key %in% names(result)) .stop(label, " duplicates key ", key)
    result[[key]] <- .safe_text(value, paste0(label, " value ", key), allow_empty = TRUE)
  }
  result
}

.require_manifest_field <- function(manifest, key, label, allow_empty = FALSE) {
  if (is.null(manifest[[key]])) .stop(label, " is missing ", key)
  .safe_text(manifest[[key]], paste0(label, " ", key), allow_empty = allow_empty)
}

.hex_digest <- function(path, algorithm) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    .stop("the digest package is required for source identity")
  }
  tolower(digest::digest(path, algo = algorithm, file = TRUE, serialize = FALSE))
}

.file_identity <- function(path, label) {
  path <- .existing_file(path, label)
  info <- file.info(path)
  list(
    path = path,
    size = as.numeric(info$size),
    md5 = .hex_digest(path, "md5"),
    sha256 = .hex_digest(path, "sha256")
  )
}
.validate_identity_manifests <- function(
  source_root,
  source_manifest_path,
  runtime_identity_path,
  run_root
) {
  source_manifest_path <- .existing_file(source_manifest_path, "source manifest")
  runtime_identity_path <- .existing_file(runtime_identity_path, "runtime identity")
  run_root <- .existing_directory(run_root, "run root")
  if (!identical(dirname(dirname(runtime_identity_path)), run_root) ||
      !identical(basename(dirname(runtime_identity_path)), "manifests")) {
    .stop("runtime identity is not the run-bound manifests/runtime.identity")
  }
  expected_source_manifest <- normalizePath(
    file.path(dirname(source_root), "identity", "source.manifest"),
    mustWork = FALSE
  )
  if (!identical(source_manifest_path, expected_source_manifest)) {
    .stop("source manifest is not the immutable snapshot manifest")
  }
  expected_runtime_identity <- normalizePath(
    file.path(dirname(runtime_identity_path), "runtime.identity"),
    mustWork = FALSE
  )
  if (!identical(runtime_identity_path, expected_runtime_identity)) {
    .stop("runtime identity path is malformed")
  }

  source <- .read_kv_manifest(source_manifest_path, "source manifest")
  runtime <- .read_kv_manifest(runtime_identity_path, "runtime identity")
  if (!identical(.require_manifest_field(source, "FORMAT", "source manifest"), "1")) {
    .stop("source manifest FORMAT must be 1")
  }
  if (!identical(.require_manifest_field(source, "SOURCE_ROOT", "source manifest"), source_root)) {
    .stop("source manifest SOURCE_ROOT does not match source root")
  }
  for (key in c(
    "SOURCE_COMMIT", "SOURCE_ARCHIVE_PATH", "SOURCE_ARCHIVE_SHA256",
    "CONFIG_HELPER_SHA256", "DATASETS_SHA256", "PIXI_TOML_SHA256",
    "PIXI_LOCK_SHA256", "AUX_ROOT", "SCGATE_DB_BRANCH"
  )) {
    .require_manifest_field(source, key, "source manifest")
  }
  for (key in c("SOURCE_ARCHIVE_SHA256", "CONFIG_HELPER_SHA256", "DATASETS_SHA256",
                "PIXI_TOML_SHA256", "PIXI_LOCK_SHA256")) {
    if (!grepl("^[0-9A-Fa-f]{64}$", source[[key]], perl = TRUE)) {
      .stop("source manifest ", key, " is not a SHA-256 digest")
    }
  }
  source_archive <- .safe_absolute(source[["SOURCE_ARCHIVE_PATH"]], "source archive path")
  if (!file.exists(source_archive)) .stop("source archive is missing: ", source_archive)
  .no_symlink(source_archive, "source archive")
  if (!identical(.hex_digest(source_archive, "sha256"), tolower(source[["SOURCE_ARCHIVE_SHA256"]]))) {
    .stop("source archive SHA-256 does not match source manifest")
  }
  datasets_path <- file.path(source_root, "datasets.json")
  if (!identical(.hex_digest(datasets_path, "sha256"), tolower(source[["DATASETS_SHA256"]]))) {
    .stop("immutable datasets.json SHA-256 does not match source manifest")
  }

  runtime_image <- .safe_absolute(
    .require_manifest_field(runtime, "RUNTIME_IMAGE", "runtime identity"),
    "runtime image path"
  )
  runtime_manifest <- .safe_absolute(
    .require_manifest_field(runtime, "RUNTIME_MANIFEST", "runtime identity"),
    "runtime manifest path"
  )
  runtime_image_identity <- .file_identity(runtime_image, "runtime image")
  runtime_manifest_identity <- .file_identity(runtime_manifest, "runtime manifest")
  runtime_image_sha <- .require_manifest_field(runtime, "RUNTIME_IMAGE_SHA256", "runtime identity")
  runtime_manifest_sha <- .require_manifest_field(runtime, "RUNTIME_MANIFEST_SHA256", "runtime identity")
  runtime_image_size <- .require_manifest_field(runtime, "RUNTIME_IMAGE_SIZE", "runtime identity")
  runtime_manifest_size <- .require_manifest_field(runtime, "RUNTIME_MANIFEST_SIZE", "runtime identity")
  if (!grepl("^[0-9A-Fa-f]{64}$", runtime_image_sha, perl = TRUE) ||
      !identical(tolower(runtime_image_sha), runtime_image_identity$sha256)) {
    .stop("runtime image SHA-256 does not match runtime identity")
  }
  if (!grepl("^[0-9A-Fa-f]{64}$", runtime_manifest_sha, perl = TRUE) ||
      !identical(tolower(runtime_manifest_sha), runtime_manifest_identity$sha256)) {
    .stop("runtime manifest SHA-256 does not match runtime identity")
  }
  if (!grepl("^[1-9][0-9]*$", runtime_image_size, perl = TRUE) ||
      as.numeric(runtime_image_size) != runtime_image_identity$size) {
    .stop("runtime image size does not match runtime identity")
  }
  if (!grepl("^[1-9][0-9]*$", runtime_manifest_size, perl = TRUE) ||
      as.numeric(runtime_manifest_size) != runtime_manifest_identity$size) {
    .stop("runtime manifest size does not match runtime identity")
  }
  runtime_manifest_fields <- .read_kv_manifest(runtime_manifest, "runtime manifest")
  runtime_format <- .require_manifest_field(runtime_manifest_fields, "FORMAT", "runtime manifest")
  if (!runtime_format %in% c("1", "2")) .stop("runtime manifest FORMAT is unsupported")

  list(
    source_manifest = c(list(path = source_manifest_path), source),
    runtime_identity = c(list(path = runtime_identity_path), runtime),
    runtime_manifest = c(
      list(
        path = runtime_manifest,
        sha256 = runtime_manifest_identity$sha256,
        size = runtime_manifest_identity$size
      ),
      runtime_manifest_fields
    ),
    runtime_image = runtime_image_identity
  )
}

.scalar_string <- function(value, label) {
  .safe_text(value, label)
}

.merge_columns <- function(entry, view, dataset) {
  base <- entry[["columns"]] %||% list()
  override <- view[["columns"]] %||% list()
  if (!is.list(base) || !is.list(override)) .stop(dataset, " columns must be objects")
  columns <- base
  if (length(override) > 0L) {
    if (is.null(names(override)) || anyNA(names(override)) || any(!nzchar(names(override)))) {
      .stop(dataset, " view columns have malformed names")
    }
    for (key in names(override)) columns[[key]] <- override[[key]]
  }
  columns
}

.batch_keys <- function(value, dataset) {
  if (is.character(value) && length(value) == 1L) {
    keys <- value
  } else if (is.character(value) && length(value) > 1L) {
    keys <- value
  } else if (is.list(value) && length(value) > 0L) {
    keys <- vapply(seq_along(value), function(index) {
      .scalar_string(value[[index]], paste0(dataset, " batch key ", index))
    }, character(1))
  } else {
    .stop(dataset, " columns.batch must be a nonempty string or array of strings")
  }
  keys <- vapply(seq_along(keys), function(index) {
    .scalar_string(keys[[index]], paste0(dataset, " batch key ", index))
  }, character(1))
  if (anyDuplicated(keys) || any(keys %in% c("Sample", "__ecoda_batch_combined_v1"))) {
    .stop(dataset, " columns.batch contains duplicate or reserved keys")
  }
  unname(keys)
}

.resolve_config <- function(config_path, dataset, view_name, input_path) {
  config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  if (!is.list(config) || is.null(names(config))) .stop("datasets.json must be an object")
  if (!nzchar(dataset) || startsWith(dataset, "_")) .stop("dataset is not a production dataset: ", dataset)
  entry <- config[[dataset]]
  if (!is.list(entry)) .stop("dataset is not configured: ", dataset)
  if (!identical(entry[["use_for_batch_effect"]], TRUE)) {
    .stop("dataset is not enabled for batch-effect processing: ", dataset)
  }
  views <- entry[["views"]]
  if (!is.list(views) || !is.list(views[[view_name]])) {
    .stop("corrected view is not configured: ", dataset, "/", view_name)
  }
  view <- views[[view_name]]
  input_name <- view[["input_file_name"]] %||% view[["input_file"]]
  output_name <- view[["output_file_name"]] %||% view[["output_file"]]
  input_name <- .scalar_string(input_name, "configured input_file_name")
  output_name <- .scalar_string(output_name, "configured output_file_name")
  if (!identical(basename(input_path), basename(input_name))) {
    .stop("source path basename does not match configured input_file_name")
  }
  if (grepl("/", output_name, fixed = TRUE) || startsWith(output_name, ".")) {
    .stop("configured output_file_name is unsafe: ", output_name)
  }
  columns <- .merge_columns(entry, view, dataset)
  sample_col <- .scalar_string(columns[["sample"]], "configured sample column")
  label_col <- .scalar_string(columns[["label"]], "configured label column")
  if (identical(sample_col, label_col)) .stop("configured sample and label columns must differ")
  batch_keys <- .batch_keys(columns[["batch"]], dataset)
  if (any(batch_keys %in% c(sample_col, label_col))) {
    .stop("configured batch keys overlap sample or label column")
  }
  subset_vars <- view[["subset_vars"]] %||% list()
  if (!is.list(subset_vars)) .stop("subset_vars must be an object")
  if (length(subset_vars) > 0L && (is.null(names(subset_vars)) || anyNA(names(subset_vars)) ||
                                  any(!nzchar(names(subset_vars))))) {
    .stop("subset_vars has malformed column names")
  }
  list(
    config = config,
    entry = entry,
    view = view,
    input_file_name = input_name,
    output_file_name = output_name,
    columns = columns,
    sample_column = sample_col,
    label_column = label_col,
    batch_keys = batch_keys,
    subset_vars = subset_vars
  )
}

.normalize_rule_values <- function(value, label) {
  if (is.null(value) || is.data.frame(value) || is.environment(value)) {
    .stop(label, " must be a scalar or nonempty sequence")
  }
  if (is.list(value)) {
    if (length(value) == 0L) .stop(label, " must not be empty")
    return(unname(value))
  }
  if (length(value) == 0L) .stop(label, " must not be empty")
  unname(as.list(value))
}

.is_missing_value <- function(value) {
  if (length(value) != 1L || is.null(value)) return(TRUE)
  if (is.factor(value)) value <- as.character(value)
  if (is.character(value)) {
    return(is.na(value) || !nzchar(trimws(value)))
  }
  if (is.logical(value)) return(is.na(value))
  if (is.numeric(value)) return(is.na(value) || !is.finite(value))
  is.na(value)[[1L]]
}

.scalar_exact_equal <- function(left, right) {
  if (length(left) != 1L || length(right) != 1L) return(FALSE)
  if (is.factor(left)) left <- as.character(left)
  if (is.factor(right)) right <- as.character(right)
  if (anyNA(left) || anyNA(right)) return(FALSE)
  if (is.character(left) || is.character(right)) {
    return(is.character(left) && is.character(right) && identical(left, right))
  }
  if (is.logical(left) || is.logical(right)) {
    return(is.logical(left) && is.logical(right) && identical(left, right))
  }
  if (is.numeric(left) || is.numeric(right)) {
    return(is.numeric(left) && is.numeric(right) && isTRUE(left == right))
  }
  identical(left, right)
}

.exact_membership <- function(column, values) {
  vapply(seq_along(column), function(index) {
    any(vapply(values, function(value) .scalar_exact_equal(column[[index]], value), logical(1)))
  }, logical(1))
}

.numeric_threshold <- function(value, label) {
  if (is.logical(value) || length(value) != 1L || .is_missing_value(value)) {
    .stop(label, " must be one finite numeric value")
  }
  text <- trimws(as.character(value))
  number <- suppressWarnings(as.numeric(text))
  if (length(number) != 1L || is.na(number) || !is.finite(number)) {
    .stop(label, " must be one finite numeric value")
  }
  number
}

.evaluate_subset <- function(metadata, subset_vars) {
  n <- nrow(metadata)
  mask <- rep(TRUE, n)
  if (length(subset_vars) == 0L) return(mask)
  operators <- c("in", "notin", "<=", "<", ">=", ">")
  for (column in names(subset_vars)) {
    if (!column %in% names(metadata)) {
      .stop("subset_vars references missing metadata column: ", column)
    }
    rule <- subset_vars[[column]]
    if (!is.list(rule) || is.null(rule[["op"]])) {
      .stop("subset rule for ", column, " is malformed")
    }
    operator <- .scalar_string(
      rule[["op"]], paste0("subset rule operator for ", column)
    )
    if (!operator %in% operators) {
      .stop("unknown subset operator for ", column, ": ", operator)
    }
    if (is.null(rule[["values"]])) {
      .stop("subset rule for ", column, " is missing values")
    }
    values <- .normalize_rule_values(rule[["values"]], paste0(column, ".values"))
    has_include <- !is.null(rule[["include_values"]])
    if (has_include && operator %in% c("in", "notin")) {
      .stop("include_values is only valid for comparison rules: ", column)
    }
    column_values <- metadata[[column]]
    if (is.list(column_values) && !is.factor(column_values)) {
      .stop("subset metadata column is list-valued: ", column)
    }
    if (operator %in% c("in", "notin")) {
      column_mask <- .exact_membership(column_values, values)
      if (operator == "notin") column_mask <- !column_mask
      invalid_rows <- vapply(column_values, .is_missing_value, logical(1))
      # A fail-closed subset never turns an absent/blank/non-finite value into
      # a retained row merely because the operator is ``notin``.
      column_mask[invalid_rows] <- FALSE
    } else {
      if (length(values) != 1L) {
        .stop("comparison rule for ", column, " requires one threshold")
      }
      threshold <- .numeric_threshold(
        values[[1L]], paste0(column, ".values")
      )
      text <- as.character(column_values)
      text <- trimws(text)
      numeric <- suppressWarnings(as.numeric(text))
      finite <- !is.na(numeric) & is.finite(numeric)
      column_mask <- rep(FALSE, n)
      if (operator == "<=") column_mask <- finite & numeric <= threshold
      if (operator == "<") column_mask <- finite & numeric < threshold
      if (operator == ">=") column_mask <- finite & numeric >= threshold
      if (operator == ">") column_mask <- finite & numeric > threshold
      if (has_include) {
        include_values <- .normalize_rule_values(
          rule[["include_values"]], paste0(column, ".include_values")
        )
        include_text <- vapply(seq_along(include_values), function(index) {
          value <- include_values[[index]]
          if (.is_missing_value(value)) {
            .stop(column, ".include_values contains missing/blank value")
          }
          value <- trimws(as.character(value))
          if (!nzchar(value)) {
            .stop(column, ".include_values contains blank value")
          }
          tolower(value)
        }, character(1))
        trimmed_text <- tolower(text)
        column_mask <- column_mask |
          (!is.na(trimmed_text) & trimmed_text %in% include_text)
      }
    }
    mask <- mask & as.logical(column_mask)
  }
  as.logical(mask)
}

.subset_audit <- function(metadata, mask, sample_col, context) {
  if (!is.character(metadata[[sample_col]]) && !is.factor(metadata[[sample_col]])) {
    .stop(context, ": sample column must be character or factor")
  }
  sample_values <- as.character(metadata[[sample_col]])
  invalid <- is.na(sample_values) | !nzchar(trimws(sample_values))
  if (any(invalid)) .stop(context, ": sample metadata contains missing/blank IDs")
  ordered_samples <- unique(sample_values)
  retained_samples <- ordered_samples[vapply(ordered_samples, function(sample) {
    any(sample_values == sample & mask)
  }, logical(1))]
  dropped_samples <- ordered_samples[vapply(ordered_samples, function(sample) {
    any(sample_values == sample & !mask)
  }, logical(1))]
  split_samples <- intersect(retained_samples, dropped_samples)
  if (length(split_samples) > 0L) {
    .stop(context, ": subset splits sample(s): ", paste(head(split_samples, 5L), collapse = ", "))
  }
  list(
    context = context,
    sample_column = sample_col,
    total_cells = as.integer(length(mask)),
    retained_cells = as.integer(sum(mask)),
    dropped_cells = as.integer(sum(!mask)),
    total_samples = as.integer(length(ordered_samples)),
    retained_samples = as.integer(length(retained_samples)),
    dropped_samples = as.integer(length(dropped_samples)),
    split_sample_count = 0L,
    retained_sample_ids = unname(retained_samples),
    dropped_sample_ids = unname(dropped_samples),
    split_sample_ids = character()
  )
}

.extract_metadata <- function(source_path) {
  # This is intentionally the first operation involving the deserialized
  # object: obtain @meta.data and immediately release the parent object.
  parent <- readRDS(source_path)
  metadata <- if (is.data.frame(parent)) {
    parent
  } else if (isS4(parent) && "meta.data" %in% methods::slotNames(parent)) {
    methods::slot(parent, "meta.data")
  } else {
    .stop("RDS source is neither a data.frame nor a Seurat object with @meta.data")
  }
  rm(parent)
  invisible(gc())
  if (!is.data.frame(metadata)) .stop("Seurat @meta.data is not a data.frame")
  metadata
}

.as_json_object <- function(value) {
  if (is.null(value) || is.null(names(value))) return(value)
  as.list(value)
}

.compact_validation <- function(validation, contract_identity) {
  list(
    valid = isTRUE(validation$valid),
    sample_column = validation$sample_col %||% validation$sample_column,
    biological_column = validation$biological_label %||% validation$biological_column,
    batch_keys = validation$ordered_keys %||% validation$keys,
    n_cells = as.integer(validation$n_cells %||% validation$n_obs),
    n_samples = as.integer(validation$n_samples),
    key_level_counts = .as_json_object(validation$key_level_counts),
    key_near_unique_fraction = .as_json_object(validation$key_near_unique_fraction),
    near_unique_fraction = validation$near_unique_fraction,
    estimable = isTRUE(validation$estimable),
    design_rank = as.integer(validation$design_rank),
    design_columns = as.integer(validation$design_columns),
    composite_design_rank = as.integer(validation$composite_design_rank),
    composite_design_columns = as.integer(validation$composite_design_columns),
    fingerprint = contract_identity$fingerprint,
    method_id = contract_identity$method_id,
    model_id = contract_identity$model_id
  )
}

.parse_args <- function(args) {
  if (length(args) == 0L) .stop("arguments are required")
  result <- list()
  index <- 1L
  while (index <= length(args)) {
    flag <- args[[index]]
    if (!startsWith(flag, "--") || index == length(args)) .stop("malformed argument: ", flag)
    key <- switch(flag,
      "--config" = "config",
      "--config-path" = "config",
      "--input" = "input",
      "--input-file" = "input",
      "--output" = "output",
      "--output-file" = "output",
      "--dataset" = "dataset",
      "--ds-name" = "dataset",
      "--view" = "view",
      "--source-root" = "source_root",
      "--source-manifest" = "source_manifest",
      "--runtime-identity" = "runtime_identity",
      "--run-root" = "run_root",
      NULL
    )
    if (is.null(key)) .stop("unknown argument: ", flag)
    if (!is.null(result[[key]])) .stop("duplicate argument: ", flag)
    result[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  required <- c("config", "input", "output", "dataset", "view", "source_root",
                "source_manifest", "runtime_identity", "run_root")
  if (any(vapply(required, function(key) is.null(result[[key]]), logical(1)))) {
    .stop("missing required corrected source metadata audit argument")
  }
  result
}

.write_json_atomic <- function(payload, path) {
  tmp <- paste0(path, ".build.", Sys.getpid())
  if (file.exists(tmp)) unlink(tmp)
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  jsonlite::write_json(
    payload,
    path = tmp,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  if (!file.rename(tmp, path)) .stop("could not atomically install report: ", path)
  invisible(path)
}

main <- function() {
  args <- .parse_args(commandArgs(trailingOnly = TRUE))
  source_root <- .existing_directory(args$source_root, "source root")
  if (!identical(basename(source_root), "tree")) .stop("source root must end in /tree")
  run_root <- .existing_directory(args$run_root, "run root")
  input_path <- .existing_file(args$input, "RDS source")
  if (!grepl("\\.rds$", input_path, ignore.case = TRUE, perl = TRUE)) {
    .stop("RDS source path must end in .rds")
  }
  output_path <- .run_owned_output(args$output, run_root)
  config_path <- .existing_file(args$config, "configuration")
  expected_config <- normalizePath(file.path(source_root, "datasets.json"), mustWork = TRUE)
  if (!identical(config_path, expected_config)) .stop("configuration is not the immutable snapshot datasets.json")
  if (!identical(args$view, "batch_effect_corrected")) .stop("audit requires batch_effect_corrected view")
  dataset <- .safe_text(args$dataset, "dataset")
  view_name <- "batch_effect_corrected"
  provenance <- .validate_identity_manifests(
    source_root, args$source_manifest, args$runtime_identity, run_root
  )
  resolved <- .resolve_config(config_path, dataset, view_name, input_path)

  metadata <- .extract_metadata(input_path)
  if (is.null(names(metadata)) || anyNA(names(metadata)) || anyDuplicated(names(metadata))) {
    .stop("source metadata must have unique named columns")
  }
  required_columns <- unique(c(
    resolved$sample_column, resolved$label_column, resolved$batch_keys,
    names(resolved$subset_vars)
  ))
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns) > 0L) {
    .stop("source metadata is missing configured columns: ", paste(missing_columns, collapse = ", "))
  }
  subset_mask <- .evaluate_subset(metadata, resolved$subset_vars)
  subset_audit <- .subset_audit(
    metadata, subset_mask, resolved$sample_column,
    paste0(dataset, "/", view_name, " corrected source subset")
  )
  if (subset_audit$retained_cells < 1L) .stop("corrected source subset retained no cells")
  selected_metadata <- metadata[subset_mask, , drop = FALSE]

  contract_path <- .existing_file(
    file.path(source_root, "src", "utils", "batch_contract.R"),
    "snapshot batch contract"
  )
  contract_env <- new.env(parent = globalenv())
  sys.source(contract_path, envir = contract_env)
  validation <- contract_env$ecoda_batch_validate_metadata(
    metadata = selected_metadata,
    batch_keys = as.list(resolved$batch_keys),
    sample_col = resolved$sample_column,
    biological_label = resolved$label_column
  )
  contract_identity <- contract_env$ecoda_batch_contract_identity(
    batch_keys = as.list(resolved$batch_keys),
    sample_col = resolved$sample_column,
    method_id = "preprocess",
    model_id = "hvg_composite_v1"
  )
  source_identity <- .file_identity(input_path, "RDS source")
  report <- list(
    schema_version = 1L,
    status = "SOURCE_METADATA_VALIDATED_RDS",
    dataset = dataset,
    view = view_name,
    source_type = "rds",
    source_path = source_identity$path,
    source_identity = source_identity,
    config = list(
      dataset = dataset,
      view = view_name,
      input_file_name = resolved$input_file_name,
      output_file_name = resolved$output_file_name,
      sample_column = resolved$sample_column,
      label_column = resolved$label_column,
      batch_keys = unname(resolved$batch_keys),
      subset_vars = resolved$subset_vars
    ),
    subset_audit = subset_audit,
    validation_summary = .compact_validation(validation, contract_identity),
    provenance = provenance
  )
  if (!isTRUE(report$validation_summary$valid) || !isTRUE(report$validation_summary$estimable)) {
    .stop("corrected source metadata validation did not produce a valid estimable contract")
  }
  .write_json_atomic(report, output_path)
  invisible(TRUE)
}

tryCatch(
  main(),
  error = function(error) {
    cat("ERROR: ", conditionMessage(error), "\n", file = stderr())
    quit(save = "no", status = 1L)
  }
)
