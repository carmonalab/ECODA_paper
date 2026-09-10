#!/usr/bin/env Rscript
# ==============================================================================
# extract_derived_composition.R
#
# Build one compact, run-owned ECODA composition snapshot on the HPC.  The
# source H5ADs stay on scratch: this script reads only the requested obs
# columns through h5ad_obs_free.py, once per dataset, and never materializes
# X, layers["counts"], or an embedding.  The resulting RDS contains only
# sample labels and sample-by-cell-type counts for the three derived analyses.
#
# This file is deliberately sourceable.  The CLI is guarded at the bottom and
# does not dispatch any of the ordinary Pipeline 1--5 scripts.
# ==============================================================================

.ecoda_snapshot_missing_sentinels <- c(
  "NA", "nan", "None", "<NA>", "Unknown", "n/a", "null"
)
.ecoda_snapshot_state <- new.env(parent = emptyenv())

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}

.ecoda_snapshot_scalar <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(label, " must be one non-empty string")
  }
  value
}

.ecoda_snapshot_validate_run_id <- function(value) {
  value <- .ecoda_snapshot_scalar(as.character(value), "run_id")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", value, perl = TRUE)) {
    stop("run_id is invalid")
  }
  value
}

.ecoda_snapshot_path <- function(path, must_work = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("path must be one non-empty string")
  }
  path <- path.expand(path)
  normalizePath(path, mustWork = must_work)
}

.ecoda_snapshot_is_symlink <- function(path) {
  link <- tryCatch(Sys.readlink(path), error = function(error) "")
  is.character(link) && length(link) == 1L && !is.na(link) && nzchar(link)
}

.ecoda_snapshot_within <- function(path, root) {
  path <- .ecoda_snapshot_path(path)
  root <- .ecoda_snapshot_path(root)
  if (identical(path, root)) return(TRUE)
  boundary <- if (identical(root, .Platform$file.sep)) {
    root
  } else {
    paste0(root, .Platform$file.sep)
  }
  startsWith(path, boundary)
}

.ecoda_snapshot_md5 <- function(path) {
  path <- .ecoda_snapshot_path(path, must_work = TRUE)
  digest <- unname(tools::md5sum(path))
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[0-9a-fA-F]{32}$", digest, perl = TRUE)) {
    stop("Could not compute MD5 for: ", path)
  }
  tolower(as.character(digest))
}

.ecoda_snapshot_is_missing <- function(values) {
  values <- trimws(as.character(values))
  is.na(values) | !nzchar(values) |
    tolower(values) %in% tolower(.ecoda_snapshot_missing_sentinels)
}
.ecoda_snapshot_is_missing_high_res <- function(values) {
  values <- trimws(as.character(values))
  .ecoda_snapshot_is_missing(values) |
    tolower(values) == "unassigned"
}

.ecoda_snapshot_atomic_text <- function(path, lines) {
  path <- .ecoda_snapshot_path(path)
  if (.ecoda_snapshot_is_symlink(path)) {
    stop("refusing to replace symbolic link: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid(), ".", as.integer(runif(1L, 1L, 2147483647L)))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(as.character(lines), temporary, useBytes = TRUE)
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
      file.info(temporary)$size <= 0) {
    stop("empty temporary publication: ", temporary)
  }
  if (file.exists(path) || .ecoda_snapshot_is_symlink(path) ||
      !file.rename(temporary, path)) {
    stop("could not atomically publish: ", path)
  }
  invisible(path)
}

.ecoda_snapshot_to_json <- function(value) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to publish the composition snapshot")
  }
  jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null",
    dataframe = "rows"
  )
}

.ecoda_snapshot_atomic_json <- function(value, path) {
  .ecoda_snapshot_atomic_text(path, paste0(.ecoda_snapshot_to_json(value), "\n"))
}

.ecoda_snapshot_atomic_rds <- function(value, path) {
  path <- .ecoda_snapshot_path(path)
  if (.ecoda_snapshot_is_symlink(path)) {
    stop("refusing to replace symbolic link: ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid(), ".", as.integer(runif(1L, 1L, 2147483647L)))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
      file.info(temporary)$size <= 0) {
    stop("empty temporary RDS publication: ", temporary)
  }
  if (file.exists(path) || .ecoda_snapshot_is_symlink(path) ||
      !file.rename(temporary, path)) {
    stop("could not atomically publish RDS: ", path)
  }
  invisible(path)
}

.ecoda_snapshot_write_sidecar <- function(path) {
  path <- .ecoda_snapshot_path(path, must_work = TRUE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) stop("artifact is empty: ", path)
  sidecar <- paste0(path, ".md5")
  .ecoda_snapshot_atomic_text(
    sidecar,
    c(
      paste0("MD5=", .ecoda_snapshot_md5(path)),
      paste0("SIZE=", as.character(info$size)),
      paste0("PATH=", path)
    )
  )
  invisible(sidecar)
}

.ecoda_snapshot_verify_h5ad_sidecar <- function(path) {
  path <- .ecoda_snapshot_path(path, must_work = TRUE)
  info <- file.info(path)
  if (is.na(info$isdir) || isTRUE(info$isdir) || is.na(info$size) ||
      info$size <= 0 || is.na(info$mtime)) {
    stop("H5AD is missing, empty, or not a regular file: ", path)
  }
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || .ecoda_snapshot_is_symlink(sidecar)) {
    stop("H5AD checksum sidecar is missing: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("H5AD checksum sidecar has the wrong schema: ", sidecar)
  }
  expected_md5 <- sub("^MD5=", "", lines[[1L]])
  expected_size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[0-9a-fA-F]{32}$", expected_md5, perl = TRUE) ||
      !grepl("^[1-9][0-9]*$", expected_size, perl = TRUE) ||
      as.character(info$size) != expected_size || !nzchar(recorded_path)) {
    stop("H5AD checksum sidecar metadata mismatch: ", sidecar)
  }
  recorded_path <- tryCatch(
    .ecoda_snapshot_path(recorded_path),
    error = function(error) NA_character_
  )
  if (is.na(recorded_path) || !identical(recorded_path, path)) {
    stop("H5AD checksum sidecar PATH mismatch: ", sidecar)
  }
  actual_md5 <- .ecoda_snapshot_md5(path)
  if (!identical(tolower(expected_md5), actual_md5)) {
    stop("H5AD checksum mismatch: ", path)
  }
  list(
    MD5 = actual_md5,
    SIZE = as.numeric(info$size),
    MTIME = as.numeric(info$mtime),
    PATH = path
  )
}

.ecoda_snapshot_validate_output_dir <- function(output_dir) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty string")
  }
  candidate <- path.expand(output_dir)
  if (.ecoda_snapshot_is_symlink(candidate)) {
    stop("output_dir must not be a symbolic link: ", candidate)
  }
  if (file.exists(candidate) && !dir.exists(candidate)) {
    stop("output_dir must name a directory: ", candidate)
  }
  if (dir.exists(candidate) &&
      length(list.files(candidate, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir must be a new empty run-owned directory: ", candidate)
  }
  if (!dir.exists(candidate)) {
    dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(candidate) || .ecoda_snapshot_is_symlink(candidate)) {
    stop("could not create output_dir: ", candidate)
  }
  .ecoda_snapshot_path(candidate, must_work = TRUE)
}

.ecoda_snapshot_validate_path_component <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || value %in% c(".", "..") ||
      grepl("[/\\\\]", value) || grepl("[[:cntrl:]]", value)) {
    stop(label, " is invalid")
  }
  value
}

.ecoda_snapshot_resolve_h5ad <- function(input_dir, dataset, view_spec) {
  input_dir <- .ecoda_snapshot_path(input_dir, must_work = TRUE)
  dataset <- .ecoda_snapshot_validate_path_component(
    dataset, "benchmark dataset key"
  )
  output_file <- view_spec$output_file_name %||% view_spec$output_file
  if (!is.character(output_file) || length(output_file) != 1L ||
      is.na(output_file) || !nzchar(output_file) ||
      grepl("^/", output_file, perl = TRUE) ||
      grepl("^[A-Za-z]:[/\\\\]", output_file, perl = TRUE) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", output_file, perl = TRUE) ||
      grepl("[[:cntrl:]]", output_file)) {
    stop("benchmark_analysis output file is invalid for ", dataset)
  }
  candidate <- .ecoda_snapshot_path(
    file.path(input_dir, dataset, "output", output_file)
  )
  if (!.ecoda_snapshot_within(candidate, input_dir)) {
    stop("benchmark_analysis H5AD escapes input_dir: ", candidate)
  }
  if (!file.exists(candidate) || dir.exists(candidate)) {
    stop("declared benchmark_analysis H5AD is missing: ", candidate)
  }
  candidate
}

# Return the raw-config union without silently filtering malformed entries.
# A dataset is eligible when use_for_benchmark is true OR it declares a
# benchmark_analysis view.  Underscore-prefixed keys are reserved for debug.
ecoda_snapshot_raw_config_union <- function(raw_config, scope = "benchmark_union") {
  if (!is.list(raw_config) || is.null(names(raw_config))) {
    stop("datasets config must be a named object")
  }
  scope <- match.arg(as.character(scope), c("benchmark_union", "debug"))
  config_names <- names(raw_config)
  if (identical(scope, "debug")) {
    return(if ("_debug" %in% config_names) "_debug" else character())
  }
  eligible <- vapply(raw_config, function(entry) {
    if (!is.list(entry)) return(FALSE)
    has_view <- is.list(entry$views) &&
      "benchmark_analysis" %in% names(entry$views)
    isTRUE(entry$use_for_benchmark) || has_view
  }, logical(1L))
  config_names[eligible & !startsWith(config_names, "_")]
}

ecoda_snapshot_targets <- function() {
  c("all cells", "2000", "1000", "500", "400", "300", "200", "150", "100", "50")
}

ecoda_snapshot_parse_seeds <- function(value = NULL) {
  if (is.null(value)) return(101L:120L)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop("--seeds must be a non-empty integer sequence")
  }
  pieces <- strsplit(value, ",", fixed = TRUE)[[1L]]
  parsed <- integer()
  for (piece in pieces) {
    piece <- trimws(piece)
    if (!nzchar(piece)) stop("--seeds contains an empty element")
    if (grepl("^[+-]?[0-9]+:[+-]?[0-9]+$", piece, perl = TRUE)) {
      bounds <- suppressWarnings(as.integer(strsplit(piece, ":", fixed = TRUE)[[1L]]))
      if (anyNA(bounds) || bounds[[1L]] > bounds[[2L]]) {
        stop("--seeds range must be finite and ascending: ", piece)
      }
      parsed <- c(parsed, seq.int(bounds[[1L]], bounds[[2L]]))
    } else if (grepl("^[+-]?[0-9]+$", piece, perl = TRUE)) {
      parsed <- c(parsed, suppressWarnings(as.integer(piece)))
    } else {
      stop("--seeds must contain integers or ascending integer ranges: ", piece)
    }
  }
  if (!length(parsed) || anyNA(parsed) || anyDuplicated(parsed)) {
    stop("--seeds must contain at least one unique finite integer")
  }
  as.integer(parsed)
}

.ecoda_snapshot_parse_flags <- function(raw_args) {
  allowed <- c("config_path", "input_dir", "output_dir", "scope", "run_id", "seeds")
  result <- list()
  i <- 1L
  while (i <= length(raw_args)) {
    token <- raw_args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token)
    name <- sub("^--", "", token)
    if (identical(name, "help")) {
      if (length(raw_args) != 1L) stop("--help cannot be combined with other arguments")
      result$help <- TRUE
      i <- i + 1L
      next
    }
    if (grepl("=", name, fixed = TRUE)) {
      pieces <- strsplit(name, "=", fixed = TRUE)[[1L]]
      if (length(pieces) != 2L || !nzchar(pieces[[1L]]) ||
          pieces[[1L]] %in% names(result)) {
        stop("Malformed or repeated argument: ", token)
      }
      flag_name <- pieces[[1L]]
      value <- pieces[[2L]]
      if (!flag_name %in% allowed || !nzchar(value)) {
        stop("Unknown or empty argument: ", token)
      }
      result[[flag_name]] <- value
      i <- i + 1L
    } else {
      if (!name %in% allowed) stop("Unknown argument: ", token)
      if (name %in% names(result)) stop("Repeated argument: --", name)
      if (i >= length(raw_args) || startsWith(raw_args[[i + 1L]], "--")) {
        stop("Missing value for --", name)
      }
      result[[name]] <- raw_args[[i + 1L]]
      i <- i + 2L
    }
  }
  result
}

.ecoda_snapshot_script_path <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  token <- full[startsWith(full, "--file=")]
  if (length(token)) {
    return(.ecoda_snapshot_path(sub("^--file=", "", token[[1L]]), must_work = TRUE))
  }
  candidate <- file.path(getwd(), "src", "5_run_benchmark_methods", "extract_derived_composition.R")
  .ecoda_snapshot_path(candidate, must_work = TRUE)
}

.ecoda_snapshot_load_obs <- function(path, obs_columns, project_root) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required for the metadata-only H5AD reader")
  }
  module_dir <- .ecoda_snapshot_path(
    file.path(project_root, "src", "utils", "py"), must_work = TRUE
  )
  if (is.null(.ecoda_snapshot_state$loader) ||
      !identical(.ecoda_snapshot_state$module_dir, module_dir)) {
    python_sys <- reticulate::import("sys", convert = FALSE)
    # Prevent an import from creating __pycache__ under the source tree.  The
    # extractor's only intended writes are the fresh output directory files.
    python_sys$dont_write_bytecode <- TRUE
    .ecoda_snapshot_state$loader <- reticulate::import_from_path(
      "h5ad_obs_free", path = module_dir, convert = FALSE
    )
    .ecoda_snapshot_state$module_dir <- module_dir
  }
  value <- .ecoda_snapshot_state$loader$load_h5ad_obs_free(
    path, as.list(as.character(obs_columns))
  )
  value <- reticulate::py_to_r(value)
  if (!is.data.frame(value)) value <- as.data.frame(value, stringsAsFactors = FALSE)
  for (column in colnames(value)) {
    if (is.factor(value[[column]])) value[[column]] <- as.character(value[[column]])
  }
  value
}

.ecoda_snapshot_as_char_column <- function(obs, column) {
  if (!column %in% colnames(obs)) {
    stop("H5AD obs is missing requested column: ", column)
  }
  value <- obs[[column]]
  if (is.list(value)) {
    stop("H5AD obs column is not scalar-valued: ", column)
  }
  as.character(value)
}

.ecoda_snapshot_count_rows <- function(
  dataset,
  analysis,
  target,
  replicate,
  seed,
  sample_ids,
  category_values,
  categories,
  groups,
  total_cells,
  annotated_cells,
  original_cells,
  effective_cells
) {
  if (!length(categories) || length(groups) != length(sample_ids) ||
      length(total_cells) != length(sample_ids) ||
      length(annotated_cells) != length(sample_ids) ||
      length(original_cells) != length(sample_ids) ||
      length(effective_cells) != length(sample_ids)) {
    stop("invalid compact composition dimensions for ", dataset, "/", analysis)
  }
  rows <- vector("list", length(sample_ids) * length(categories))
  row_index <- 0L
  for (sample_index in seq_along(sample_ids)) {
    keep <- as.integer(groups[[sample_index]])
    if (length(keep) != effective_cells[[sample_index]]) {
      stop("effective cell count mismatch for ", dataset, "/", sample_ids[[sample_index]])
    }
    category_index <- match(category_values[keep], categories)
    if (anyNA(category_index)) stop("invalid category while composing ", dataset)
    counts <- tabulate(category_index, nbins = length(categories))
    if (sum(counts) != effective_cells[[sample_index]]) {
      stop("composition counts do not sum to effective cells for ", dataset)
    }
    for (category_index_value in seq_along(categories)) {
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        dataset = dataset,
        analysis = analysis,
        target = target,
        replicate = as.integer(replicate),
        seed = as.integer(seed),
        Sample = sample_ids[[sample_index]],
        category = categories[[category_index_value]],
        count = as.integer(counts[[category_index_value]]),
        total_cells = as.integer(total_cells[[sample_index]]),
        annotated_cells = as.integer(annotated_cells[[sample_index]]),
        original_cells = as.integer(original_cells[[sample_index]]),
        effective_cells = as.integer(effective_cells[[sample_index]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.ecoda_snapshot_compose_dataset <- function(obs, dataset, label_col, high_res_col, seeds) {
  required <- unique(c(
    "Sample",
    label_col,
    high_res_col,
    "leiden_res_50_benchmark_analysis_hvg2000",
    "leiden_res_2_benchmark_analysis_hvg2000_harmony"
  ))
  missing <- setdiff(required, colnames(obs))
  if (length(missing)) {
    stop("H5AD obs is missing required columns for ", dataset, ": ",
         paste(missing, collapse = ", "))
  }
  sample_values <- .ecoda_snapshot_as_char_column(obs, "Sample")
  if (!length(sample_values) || any(.ecoda_snapshot_is_missing(sample_values))) {
    stop("H5AD Sample contains missing or blank values for ", dataset)
  }
  sample_ids <- unique(sample_values)
  total_groups <- lapply(
    sample_ids,
    function(sample_id) which(sample_values == sample_id)
  )
  total_cells <- as.integer(lengths(total_groups))
  names(total_cells) <- sample_ids

  high_values <- .ecoda_snapshot_as_char_column(obs, high_res_col)
  high_missing <- .ecoda_snapshot_is_missing_high_res(high_values)
  usable_high_res <- vapply(
    total_groups,
    function(indices) any(!high_missing[indices]),
    logical(1L)
  )
  if (any(!usable_high_res)) {
    stop("H5AD high-resolution column '", high_res_col,
         "' has no usable cell type for sample(s): ",
         paste(sample_ids[!usable_high_res], collapse = ", "))
  }
  annotated_keep <- !high_missing
  annotated_sample_values <- sample_values[annotated_keep]
  annotated_high_values <- high_values[annotated_keep]
  groups <- lapply(
    sample_ids,
    function(sample_id) which(annotated_sample_values == sample_id)
  )
  annotated_cells <- as.integer(lengths(groups))
  names(annotated_cells) <- sample_ids
  if (any(annotated_cells <= 0L)) {
    stop("H5AD high-resolution filtering removed every annotated cell for sample(s): ",
         paste(sample_ids[annotated_cells <= 0L], collapse = ", "))
  }

  label_values <- .ecoda_snapshot_as_char_column(obs, label_col)
  label_rows <- vector("list", length(sample_ids))
  for (sample_index in seq_along(sample_ids)) {
    values <- label_values[total_groups[[sample_index]]]
    if (any(.ecoda_snapshot_is_missing(values))) {
      stop("biological labels for sample '", sample_ids[[sample_index]],
           "' are missing or blank in ", dataset)
    }
    if (length(unique(values)) != 1L) {
      stop("biological labels for sample '", sample_ids[[sample_index]],
           "' are conflicting in ", dataset)
    }
    label_rows[[sample_index]] <- data.frame(
      dataset = dataset,
      Sample = sample_ids[[sample_index]],
      label = values[[1L]],
      total_cells = total_cells[[sample_index]],
      annotated_cells = annotated_cells[[sample_index]],
      original_cells = annotated_cells[[sample_index]],
      stringsAsFactors = FALSE
    )
  }
  labels <- do.call(rbind, label_rows)

  res50_values <- .ecoda_snapshot_as_char_column(
    obs, "leiden_res_50_benchmark_analysis_hvg2000"
  )
  harmony_values <- .ecoda_snapshot_as_char_column(
    obs, "leiden_res_2_benchmark_analysis_hvg2000_harmony"
  )
  for (pair in list(
    list(name = "res50", values = res50_values),
    list(name = "harmony", values = harmony_values)
  )) {
    if (any(.ecoda_snapshot_is_missing(pair$values))) {
      stop("H5AD technical cluster column contains missing values for ",
           dataset, "/", pair$name)
    }
  }

  count_parts <- list()
  count_index <- 0L
  count_index <- count_index + 1L
  count_parts[[count_index]] <- .ecoda_snapshot_count_rows(
    dataset = dataset,
    analysis = "res50",
    target = "all cells",
    replicate = 0L,
    seed = NA_integer_,
    sample_ids = sample_ids,
    category_values = res50_values,
    categories = unique(res50_values),
    groups = total_groups,
    total_cells = total_cells,
    annotated_cells = annotated_cells,
    original_cells = total_cells,
    effective_cells = total_cells
  )
  count_index <- count_index + 1L
  count_parts[[count_index]] <- .ecoda_snapshot_count_rows(
    dataset = dataset,
    analysis = "harmony",
    target = "all cells",
    replicate = 0L,
    seed = NA_integer_,
    sample_ids = sample_ids,
    category_values = harmony_values,
    categories = unique(harmony_values),
    groups = total_groups,
    total_cells = total_cells,
    annotated_cells = annotated_cells,
    original_cells = total_cells,
    effective_cells = total_cells
  )

  targets <- ecoda_snapshot_targets()
  if (!identical(as.integer(seeds), 101L:120L)) {
    stop("cell-subsetting seeds must use the exact publication schedule 101:120")
  }
  high_categories <- unique(annotated_high_values)
  if (any(.ecoda_snapshot_is_missing_high_res(high_categories))) {
    stop("annotated high-resolution categories contain missing values for ", dataset)
  }
  schedule_index <- 0L
  for (target in targets) {
    if (identical(target, "all cells")) {
      schedule_index <- schedule_index + 1L
      count_index <- count_index + 1L
      count_parts[[count_index]] <- .ecoda_snapshot_count_rows(
        dataset = dataset,
        analysis = "cell_subsetting",
        target = target,
        replicate = 0L,
        seed = NA_integer_,
        sample_ids = sample_ids,
        category_values = annotated_high_values,
        categories = high_categories,
        groups = groups,
        total_cells = total_cells,
        annotated_cells = annotated_cells,
        original_cells = annotated_cells,
        effective_cells = annotated_cells
      )
      next
    }
    target_cells <- as.integer(target)
    for (replicate in seq_along(seeds)) {
      seed <- as.integer(seeds[[replicate]])
      set.seed(seed)
      effective_cells <- as.integer(pmin(lengths(groups), target_cells))
      sampled_groups <- lapply(seq_along(groups), function(sample_index) {
        indices <- groups[[sample_index]]
        n_keep <- effective_cells[[sample_index]]
        if (n_keep == length(indices)) {
          indices
        } else {
          # This is intentionally the base-R per-sample sampling operation;
          # do not replace it with a multinomial or hypergeometric draw.
          sample(indices, size = n_keep, replace = FALSE)
        }
      })
      schedule_index <- schedule_index + 1L
      count_index <- count_index + 1L
      count_parts[[count_index]] <- .ecoda_snapshot_count_rows(
        dataset = dataset,
        analysis = "cell_subsetting",
        target = target,
        replicate = replicate,
        seed = seed,
        sample_ids = sample_ids,
        category_values = annotated_high_values,
        categories = high_categories,
        groups = sampled_groups,
        total_cells = total_cells,
        annotated_cells = annotated_cells,
        original_cells = annotated_cells,
        effective_cells = effective_cells
      )
    }
  }
  if (!identical(schedule_index, 181L)) {
    stop("cell-subsetting schedule is not exactly 181 rows for ", dataset)
  }
  list(labels = labels, counts = do.call(rbind, count_parts))
}

.ecoda_snapshot_validate_tables <- function(snapshot, datasets) {
  if (!is.list(snapshot) || !all(c("labels", "counts", "metadata") %in% names(snapshot))) {
    stop("snapshot object has the wrong top-level fields")
  }
  labels <- snapshot$labels
  counts <- snapshot$counts
  labels_required <- c(
    "dataset", "Sample", "label", "total_cells", "annotated_cells",
    "original_cells"
  )
  counts_required <- c(
    "dataset", "analysis", "target", "replicate", "seed", "Sample",
    "category", "count", "total_cells", "annotated_cells",
    "original_cells", "effective_cells"
  )
  if (!is.data.frame(labels) || !all(labels_required %in% colnames(labels)) ||
      !is.data.frame(counts) || !all(counts_required %in% colnames(counts))) {
    stop("snapshot tables are missing required columns")
  }
  if (!nrow(labels) || !nrow(counts) || !length(datasets)) {
    stop("snapshot tables must be nonempty")
  }
  if (anyNA(labels$dataset) || anyNA(labels$Sample) ||
      any(!nzchar(as.character(labels$dataset))) ||
      any(!nzchar(as.character(labels$Sample))) ||
      anyNA(labels$label) || any(.ecoda_snapshot_is_missing(labels$label)) ||
      anyDuplicated(paste(labels$dataset, labels$Sample, sep = "\r"))) {
    stop("snapshot labels are invalid")
  }
  label_total <- suppressWarnings(as.numeric(labels$total_cells))
  label_annotated <- suppressWarnings(as.numeric(labels$annotated_cells))
  label_original <- suppressWarnings(as.numeric(labels$original_cells))
  if (any(!is.finite(label_total)) || any(label_total <= 0) ||
      any(label_total != as.integer(label_total)) ||
      any(!is.finite(label_annotated)) || any(label_annotated <= 0) ||
      any(label_annotated != as.integer(label_annotated)) ||
      any(!is.finite(label_original)) || any(label_original <= 0) ||
      any(label_original != as.integer(label_original)) ||
      any(label_annotated > label_total) ||
      any(label_original != label_annotated) ||
      any(!as.character(labels$dataset) %in% as.character(datasets))) {
    stop("snapshot labels contain invalid cell diagnostics")
  }
  if (!identical(
    sort(unique(as.character(labels$dataset))),
    sort(as.character(datasets))
  )) {
    stop("snapshot labels are not dataset-complete")
  }

  if (anyNA(counts$dataset) || anyNA(counts$analysis) ||
      anyNA(counts$target) || anyNA(counts$Sample) ||
      anyNA(counts$category) ||
      any(!nzchar(as.character(counts$dataset))) ||
      any(!nzchar(as.character(counts$Sample))) ||
      any(!nzchar(as.character(counts$category))) ||
      any(!as.character(counts$analysis) %in%
          c("res50", "harmony", "cell_subsetting")) ||
      anyNA(counts$count) || anyNA(counts$total_cells) ||
      anyNA(counts$annotated_cells) || anyNA(counts$original_cells) ||
      anyNA(counts$effective_cells) ||
      any(
        tolower(trimws(as.character(counts$category))) == "unassigned" &
          as.character(counts$analysis) == "cell_subsetting"
      )) {
    stop("snapshot counts are invalid")
  }
  count_values <- lapply(
    counts[c("count", "total_cells", "annotated_cells",
             "original_cells", "effective_cells")],
    function(value) suppressWarnings(as.numeric(value))
  )
  if (any(vapply(count_values, function(value) {
    any(!is.finite(value)) || any(value < 0) ||
      any(value != as.integer(value))
  }, logical(1L))) ||
      any(count_values$total_cells <= 0) ||
      any(count_values$annotated_cells <= 0) ||
      any(count_values$original_cells <= 0) ||
      any(count_values$effective_cells <= 0) ||
      any(count_values$annotated_cells > count_values$total_cells) ||
      any(count_values$count > count_values$effective_cells) ||
      any(!as.character(counts$dataset) %in% as.character(datasets))) {
    stop("snapshot counts contain invalid cell diagnostics")
  }
  for (dataset in datasets) {
    ds_labels <- labels[labels$dataset == dataset, , drop = FALSE]
    subset_counts <- counts[
      counts$dataset == dataset & counts$analysis == "cell_subsetting",
      , drop = FALSE
    ]
    schedule <- unique(subset_counts[, c("target", "replicate", "seed"), drop = FALSE])
    if (nrow(schedule) != 181L) {
      stop("cell-subsetting schedule is not exactly 181 rows for ", dataset)
    }
    for (analysis_name in c("res50", "harmony", "cell_subsetting")) {
      analysis_rows <- counts[
        counts$dataset == dataset & counts$analysis == analysis_name,
        , drop = FALSE
      ]
      if (!nrow(analysis_rows) ||
          !identical(
            sort(unique(as.character(analysis_rows$Sample))),
            sort(as.character(ds_labels$Sample))
          )) {
        stop("snapshot counts are not sample-complete for ", dataset, "/", analysis_name)
      }
      for (sample_id in as.character(ds_labels$Sample)) {
        sample_rows <- analysis_rows[
          analysis_rows$Sample == sample_id, , drop = FALSE
        ]
        label_row <- ds_labels[ds_labels$Sample == sample_id, , drop = FALSE]
        total <- unique(as.numeric(sample_rows$total_cells))
        annotated <- unique(as.numeric(sample_rows$annotated_cells))
        original <- unique(as.numeric(sample_rows$original_cells))
        expected_total <- as.numeric(label_row$total_cells)
        expected_annotated <- as.numeric(label_row$annotated_cells)
        if (length(total) != 1L || length(annotated) != 1L ||
            length(original) != 1L ||
            !identical(total, expected_total) ||
            !identical(annotated, expected_annotated) ||
            (analysis_name %in% c("res50", "harmony") &&
             !identical(original, total)) ||
            (identical(analysis_name, "cell_subsetting") &&
             !identical(original, annotated))) {
          stop("snapshot counts have inconsistent cell diagnostics for ",
               dataset, "/", analysis_name, "/", sample_id)
        }
        if (identical(analysis_name, "cell_subsetting") &&
            any(as.numeric(sample_rows$effective_cells) > original)) {
          stop("snapshot cell-subsetting effective counts exceed annotated cells for ",
               dataset, "/", sample_id)
        }
      }
    }
  }
  invisible(TRUE)
}

# Run the full source verification/read/composition/publish operation.  The
# input_dir is the full HPC scratch root and is never written by this function.
ecoda_snapshot_extract <- function(
  config_path,
  input_dir,
  output_dir,
  scope = "benchmark_union",
  run_id = NULL,
  seeds = 101L:120L,
  script_path = NULL
) {
  config_path <- .ecoda_snapshot_path(config_path, must_work = TRUE)
  input_dir <- .ecoda_snapshot_path(input_dir, must_work = TRUE)
  if (!dir.exists(input_dir)) stop("input_dir must name a directory: ", input_dir)
  output_dir <- .ecoda_snapshot_validate_output_dir(output_dir)
  scope <- match.arg(as.character(scope), c("benchmark_union", "debug"))
  seeds <- as.integer(seeds)
  if (!identical(seeds, 101L:120L)) {
    stop("cell-subsetting seeds must use the exact publication schedule 101:120")
  }
  if (is.null(run_id)) run_id <- basename(output_dir)
  run_id <- .ecoda_snapshot_validate_run_id(as.character(run_id))
  if (is.null(script_path)) script_path <- .ecoda_snapshot_script_path()
  script_path <- .ecoda_snapshot_path(script_path, must_work = TRUE)
  project_root <- dirname(dirname(dirname(script_path)))

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read datasets.json and publish the snapshot")
  }
  raw_config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  selected <- ecoda_snapshot_raw_config_union(raw_config, scope = scope)
  if (!length(selected)) stop("raw-config union is empty for scope ", scope)

  config_md5 <- .ecoda_snapshot_md5(config_path)
  extractor_md5 <- .ecoda_snapshot_md5(script_path)
  source_rows <- list()
  label_parts <- list()
  count_parts <- list()

  for (dataset in selected) {
    entry <- raw_config[[dataset]]
    if (!is.list(entry)) stop("raw config entry is not an object: ", dataset)
    views <- entry$views
    if (!is.list(views) || !is.list(views$benchmark_analysis)) {
      stop("benchmark_analysis view is missing or invalid: ", dataset)
    }
    view_spec <- views$benchmark_analysis
    merged_columns <- entry$columns %||% list()
    if (!is.list(merged_columns)) merged_columns <- list()
    if (is.list(view_spec$columns)) {
      merged_columns <- modifyList(merged_columns, view_spec$columns)
    }
    label_col <- merged_columns$label %||% entry$label_col
    high_res_col <- merged_columns$cell_type_high_res %||% entry$cell_type_high_res
    if (!is.character(label_col) || length(label_col) != 1L || is.na(label_col) ||
        !nzchar(label_col) || !is.character(high_res_col) ||
        length(high_res_col) != 1L || is.na(high_res_col) || !nzchar(high_res_col)) {
      stop("configured label/high-resolution columns are missing: ", dataset)
    }
    h5ad_path <- .ecoda_snapshot_resolve_h5ad(input_dir, dataset, view_spec)
    source_identity <- .ecoda_snapshot_verify_h5ad_sidecar(h5ad_path)
    obs_columns <- unique(c(
      "Sample",
      label_col,
      high_res_col,
      "leiden_res_50_benchmark_analysis_hvg2000",
      "leiden_res_2_benchmark_analysis_hvg2000_harmony"
    ))
    # Exactly one metadata-only H5AD read for this dataset.  All 181 subset
    # schedules and both ordinary analyses consume this returned obs table.
    obs <- .ecoda_snapshot_load_obs(h5ad_path, obs_columns, project_root)
    post_info <- file.info(h5ad_path)
    if (is.na(post_info$size) || as.numeric(post_info$size) != source_identity$SIZE ||
        is.na(post_info$mtime) || as.numeric(post_info$mtime) != source_identity$MTIME) {
      stop("H5AD changed during metadata-only read: ", h5ad_path)
    }
    composed <- .ecoda_snapshot_compose_dataset(
      obs,
      dataset = dataset,
      label_col = label_col,
      high_res_col = high_res_col,
      seeds = seeds
    )
    label_parts[[length(label_parts) + 1L]] <- composed$labels
    count_parts[[length(count_parts) + 1L]] <- composed$counts
    source_rows[[length(source_rows) + 1L]] <- list(
      dataset = dataset,
      view = "benchmark_analysis",
      h5ad_path = h5ad_path,
      h5ad_md5 = source_identity$MD5,
      h5ad_size = source_identity$SIZE,
      h5ad_mtime = source_identity$MTIME,
      label_col = label_col,
      high_res_col = high_res_col
    )
    rm(obs, composed, source_identity, post_info)
    gc(verbose = FALSE)
  }

  snapshot <- list(
    labels = do.call(rbind, label_parts),
    counts = do.call(rbind, count_parts),
    metadata = list(
      schema_version = 1L,
      source_mode = "composition_snapshot",
      run_id = run_id,
      scope = scope,
      config_path = config_path,
      config_md5 = config_md5,
      extractor_path = script_path,
      extractor_md5 = extractor_md5,
      analyses = c("res50", "harmony", "cell_subsetting"),
      targets = ecoda_snapshot_targets(),
      seeds = as.integer(seeds),
      cell_subsetting_schedule_rows_per_dataset = 181L,
      datasets = selected
    )
  )
  .ecoda_snapshot_validate_tables(snapshot, selected)

  snapshot_path <- .ecoda_snapshot_path(
    file.path(output_dir, "derived_composition_snapshot.rds")
  )
  snapshot_json_path <- .ecoda_snapshot_path(
    file.path(output_dir, "derived_composition_snapshot.json")
  )
  manifest_path <- .ecoda_snapshot_path(
    file.path(output_dir, "composition_snapshot_manifest.json")
  )
  # Serialize JSON before installing anything so a serialization failure cannot
  # leave an otherwise plausible but incomplete snapshot publication.
  snapshot_json_text <- paste0(.ecoda_snapshot_to_json(snapshot), "\n")
  .ecoda_snapshot_atomic_rds(snapshot, snapshot_path)
  .ecoda_snapshot_write_sidecar(snapshot_path)
  .ecoda_snapshot_atomic_text(snapshot_json_path, snapshot_json_text)
  .ecoda_snapshot_write_sidecar(snapshot_json_path)

  snapshot_info <- file.info(snapshot_path)
  snapshot_json_info <- file.info(snapshot_json_path)
  manifest <- list(
    schema_version = 1L,
    run_id = run_id,
    scope = scope,
    status = "COMPLETED",
    source_mode = "composition_snapshot",
    config_path = config_path,
    config_md5 = config_md5,
    extractor_path = script_path,
    extractor_md5 = extractor_md5,
    snapshot_path = snapshot_path,
    snapshot_md5 = .ecoda_snapshot_md5(snapshot_path),
    snapshot_size = as.numeric(snapshot_info$size),
    snapshot_json_path = snapshot_json_path,
    snapshot_json_md5 = .ecoda_snapshot_md5(snapshot_json_path),
    snapshot_json_size = as.numeric(snapshot_json_info$size),
    sources = source_rows
  )
  .ecoda_snapshot_atomic_json(manifest, manifest_path)
  message("Published composition snapshot: ", snapshot_path)
  invisible(manifest)
}

.ecoda_snapshot_usage <- function() {
  paste(
    "Usage:",
    "Rscript extract_derived_composition.R",
    "--config_path <datasets.json>",
    "--input_dir <HPC scratch root>",
    "--output_dir <new empty run directory>",
    "[--scope benchmark_union|debug] [--run_id <id>] [--seeds 101:120]"
  )
}

.ecoda_snapshot_cli <- function() {
  args <- .ecoda_snapshot_parse_flags(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    cat(.ecoda_snapshot_usage(), "\n")
    return(invisible(NULL))
  }
  for (required in c("config_path", "input_dir", "output_dir")) {
    if (is.null(args[[required]])) stop("Missing required --", required, " argument")
  }
  seeds <- ecoda_snapshot_parse_seeds(args$seeds)
  run_id <- args$run_id %||% NULL
  ecoda_snapshot_extract(
    config_path = args$config_path,
    input_dir = args$input_dir,
    output_dir = args$output_dir,
    scope = args$scope %||% "benchmark_union",
    run_id = run_id,
    seeds = seeds
  )
  invisible(NULL)
}

# Do not execute when sourced for focused tests.  Rscript supplies --file=;
# matching this basename keeps Rscript -e evaluation inert.
.ecoda_snapshot_runner_invocation <- any(vapply(
  commandArgs(trailingOnly = FALSE)[startsWith(
    commandArgs(trailingOnly = FALSE), "--file="
  )],
  function(token) {
    identical(
      basename(sub("^--file=", "", token)),
      "extract_derived_composition.R"
    )
  },
  logical(1L)
))
if (.ecoda_snapshot_runner_invocation) {
  tryCatch(
    .ecoda_snapshot_cli(),
    error = function(error) stop(conditionMessage(error), call. = FALSE)
  )
}
