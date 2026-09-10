# ==============================================================================
# run_local_ecoda_derived.R
#
# Local, run-owned ECODA derived analyses.  This file deliberately does not
# register methods with the ordinary Stage 5 matrix.  Each selector emits one
# standalone composition bundle: `res50` emits ordinary unsupervised Leiden
# resolution 50, `harmony` emits Harmony Leiden resolution 2, and
# `cell_subsetting` emits the fixed cell-depth experiment.
# ==============================================================================

# The definitions above the Rscript guard are intentionally sourceable.  They
# contain the selection/target/manifest contracts used by focused tests without
# loading reticulate or any benchmark package.

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}
.ecoda_missing_label_sentinels <- c("NA", "nan", "None", "<NA>", "Unknown", "n/a", "null")
.ecoda_is_missing_sentinel <- function(values) {
  values <- trimws(as.character(values))
  is.na(values) | !nzchar(values) |
    tolower(values) %in% tolower(.ecoda_missing_label_sentinels)
}
.ecoda_is_missing_high_res <- function(values) {
  values <- trimws(as.character(values))
  .ecoda_is_missing_sentinel(values) | tolower(values) == "unassigned"
}
.ecoda_in_process_owner <- new.env(parent = emptyenv())

.ecoda_clear_in_process_owner <- function() {
  .ecoda_in_process_owner$token <- NULL
  invisible(NULL)
}

.ecoda_record_in_process_owner <- function(output_dir, run_id) {
  .ecoda_in_process_owner$token <- list(
    output_dir = .ecoda_path(output_dir, must_work = TRUE),
    run_id = .ecoda_validate_run_id(run_id)
  )
  invisible(NULL)
}

.ecoda_in_process_owner_matches <- function(output_dir, run_id) {
  token <- .ecoda_in_process_owner$token
  if (!is.list(token) || length(token) != 2L ||
      !identical(names(token), c("output_dir", "run_id"))) {
    return(FALSE)
  }
  identical(token$output_dir, output_dir) &&
    identical(token$run_id, run_id)
}



.ecoda_scalar <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(label, " must be one non-empty string")
  }
  value
}

# Return the raw-config union.  This intentionally does not call
# read_datasets_json(): that helper filters entries which do not have a valid
# view, while this preflight must report those entries as blocked.
ecoda_raw_config_union <- function(
  raw_config,
  scope = c("benchmark_union", "debug")
) {
  scope <- match.arg(scope)
  if (!is.list(raw_config) || is.null(names(raw_config))) {
    stop("raw datasets config must be a named list")
  }
  names_config <- names(raw_config)
  if (identical(scope, "debug")) {
    return(if ("_debug" %in% names_config) "_debug" else character())
  }
  eligible <- vapply(raw_config, function(entry) {
    if (!is.list(entry)) return(FALSE)
    has_view <- is.list(entry$views) &&
      "benchmark_analysis" %in% names(entry$views)
    isTRUE(entry$use_for_benchmark) || has_view
  }, logical(1))
  names_config[eligible & !startsWith(names_config, "_")]
}

# Exact publication order for the cell-depth analysis.
ecoda_derived_targets <- function() {
  c("all cells", "2000", "1000", "500", "400", "300", "200", "150", "100", "50")
}

ecoda_parse_seeds <- function(value = NULL) {
  if (is.null(value)) return(101L:120L)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("--seeds must be a non-empty integer sequence")
  }
  text <- as.character(value)
  pieces <- strsplit(text, ",", fixed = TRUE)[[1L]]
  parsed <- integer()
  for (piece in pieces) {
    piece <- trimws(piece)
    if (!nzchar(piece)) stop("--seeds contains an empty element")
    if (grepl("^[+-]?[0-9]+:[+-]?[0-9]+$", piece, perl = TRUE)) {
      bounds <- as.integer(strsplit(piece, ":", fixed = TRUE)[[1L]])
      if (anyNA(bounds) || bounds[1L] > bounds[2L]) {
        stop("--seeds range must be finite and ascending: ", piece)
      }
      parsed <- c(parsed, seq.int(bounds[1L], bounds[2L]))
    } else if (grepl("^[+-]?[0-9]+$", piece, perl = TRUE)) {
      parsed <- c(parsed, as.integer(piece))
    } else {
      stop("--seeds must contain integers or ascending integer ranges: ", piece)
    }
  }
  if (!length(parsed) || anyDuplicated(parsed) || any(!is.finite(parsed))) {
    stop("--seeds must contain at least one unique finite integer")
  }
  as.integer(parsed)
}

# Build only the deterministic target/replicate schedule.  The execution path
# adds sample/effective-count diagnostics after each composition call.
ecoda_subsetting_plan <- function(
  sample_ids,
  counts,
  seeds = 101L:120L,
  targets = ecoda_derived_targets()
) {
  sample_ids <- as.character(sample_ids)
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("subsetting sample_ids must be nonempty and unique")
  }
  if (is.null(names(counts)) ||
      !identical(as.character(names(counts)), sample_ids) ||
      length(counts) != length(sample_ids) || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != as.integer(counts))) {
    stop("subsetting counts must be named, finite, nonnegative integers aligned to sample_ids")
  }
  seeds <- as.integer(seeds)
  if (!length(seeds) || anyNA(seeds) || anyDuplicated(seeds)) {
    stop("subsetting seeds must be unique finite integers")
  }
  if (!identical(seeds, 101L:120L)) {
    stop("subsetting seeds must use the exact publication schedule 101:120")
  }

  if (!identical(as.character(targets), ecoda_derived_targets())) {
    stop("subsetting targets must use the exact ordered target contract")
  }
  rows <- list()
  row_index <- 0L
  for (target in targets) {
    if (identical(target, "all cells")) {
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        target = target,
        target_cells = NA_integer_,
        replicate = 0L,
        seed = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      target_cells <- as.integer(target)
      for (replicate in seq_along(seeds)) {
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          target = target,
          target_cells = target_cells,
          replicate = as.integer(replicate),
          seed = as.integer(seeds[[replicate]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# Validate the stable portions of either derived manifest.  Source manifests
# have no self-referential source_manifest_md5; run manifests do.  This helper
# accepts both forms while checking every field that is present.
ecoda_validate_derived_manifest <- function(
  manifest,
  required_sources = TRUE
) {
  if (!is.list(manifest)) stop("derived manifest must be an object")
  required <- c(
    "run_id", "analysis", "scope", "status", "config_path", "config_md5",
    "runner_path", "runner_md5", "sources"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("derived manifest is missing: ", paste(missing, collapse = ", "))
  if (!is.character(manifest$run_id) || length(manifest$run_id) != 1L ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", manifest$run_id, perl = TRUE)) {
    stop("derived manifest run_id is invalid")
  }
  if (!identical(manifest$analysis, "res50") &&
      !identical(manifest$analysis, "harmony") &&
      !identical(manifest$analysis, "cell_subsetting")) {
    stop("derived manifest analysis is invalid")
  }
  if (!identical(manifest$scope, "benchmark_union") &&
      !identical(manifest$scope, "debug")) {
    stop("derived manifest scope is invalid")
  }
  if (!is.character(manifest$status) || length(manifest$status) != 1L ||
      !nzchar(manifest$status)) stop("derived manifest status is invalid")
  for (field in c("config_path", "runner_path")) {
    if (!is.character(manifest[[field]]) || length(manifest[[field]]) != 1L ||
        !nzchar(manifest[[field]])) stop("derived manifest ", field, " is invalid")
  }
  for (field in c("config_md5", "runner_md5")) {
    if (!is.character(manifest[[field]]) || length(manifest[[field]]) != 1L ||
        !grepl("^[0-9a-f]{32}$", tolower(manifest[[field]]), perl = TRUE)) {
      stop("derived manifest ", field, " is invalid")
    }
  }
  snapshot_config_fields <- c("snapshot_config_path", "snapshot_config_md5")
  snapshot_config_present <- snapshot_config_fields %in% names(manifest)
  if (any(snapshot_config_present) && !all(snapshot_config_present)) {
    stop("derived manifest snapshot config identity is incomplete")
  }
  if (all(snapshot_config_present)) {
    if (!is.character(manifest$snapshot_config_path) ||
        length(manifest$snapshot_config_path) != 1L ||
        !nzchar(manifest$snapshot_config_path) ||
        !is.character(manifest$snapshot_config_md5) ||
        length(manifest$snapshot_config_md5) != 1L ||
        !grepl("^[0-9a-f]{32}$",
               tolower(manifest$snapshot_config_md5), perl = TRUE) ||
        !identical(
          tolower(manifest$snapshot_config_md5),
          tolower(manifest$config_md5)
        )) {
      stop("derived manifest snapshot config identity is invalid")
    }
  }
  if (required_sources && !is.list(manifest$sources)) {
    stop("derived manifest sources must be an array")
  }
  if (is.list(manifest$sources) && length(manifest$sources)) {
    for (source in manifest$sources) {
      needed <- c("dataset", "view", "h5ad_path", "h5ad_md5", "analysis")
      if (!is.list(source) || !all(needed %in% names(source))) {
        stop("derived source entry is malformed")
      }
      if (!identical(source$view, "benchmark_analysis") ||
          !identical(source$analysis, manifest$analysis) ||
          !is.character(source$dataset) || length(source$dataset) != 1L ||
          !nzchar(source$dataset) || !is.character(source$h5ad_path) ||
          length(source$h5ad_path) != 1L || !nzchar(source$h5ad_path) ||
          !is.character(source$h5ad_md5) || length(source$h5ad_md5) != 1L ||
          !grepl("^[0-9a-f]{32}$", tolower(source$h5ad_md5), perl = TRUE)) {
        stop("derived source entry is invalid")
      }
    }
  }
  if ("artifacts" %in% names(manifest)) {
    if (!is.list(manifest$artifacts)) stop("derived manifest artifacts must be an array")
    if (length(manifest$artifacts)) {
      for (artifact in manifest$artifacts) {
        needed <- c("path", "size", "md5", "dataset", "method")
        if (!is.list(artifact) || !all(needed %in% names(artifact))) {
          stop("derived artifact entry is malformed")
        }
        if (!is.character(artifact$path) || length(artifact$path) != 1L ||
            !nzchar(artifact$path) || !is.numeric(artifact$size) ||
            length(artifact$size) != 1L || !is.finite(artifact$size) ||
            artifact$size <= 0 || !is.character(artifact$md5) ||
            length(artifact$md5) != 1L ||
            !grepl("^[0-9a-f]{32}$", tolower(artifact$md5), perl = TRUE) ||
            !is.character(artifact$dataset) || length(artifact$dataset) != 1L ||
            !nzchar(artifact$dataset) || !is.character(artifact$method) ||
            length(artifact$method) != 1L || !nzchar(artifact$method)) {
          stop("derived artifact entry is invalid")
        }
      }
    }
  }
  if ("artifacts" %in% names(manifest) &&
      identical(manifest$status, "COMPLETED")) {
    source_datasets <- vapply(manifest$sources, `[[`, character(1), "dataset")
    artifact_paths <- vapply(manifest$artifacts, `[[`, character(1), "path")
    artifact_datasets <- vapply(
      manifest$artifacts, `[[`, character(1), "dataset"
    )
    artifact_methods <- vapply(
      manifest$artifacts, `[[`, character(1), "method"
    )
    if (manifest$analysis %in% c("res50", "harmony")) {
      expected_method <- if (identical(manifest$analysis, "res50")) {
        "ECODA_seuratres_50"
      } else {
        "ECODA_seuratres_2_harmony"
      }
      expected_paths <- paste0(
        source_datasets, "_", expected_method, ".rds"
      )
      if (!identical(artifact_datasets, source_datasets) ||
          !identical(artifact_methods, rep(expected_method, length(source_datasets))) ||
          !identical(basename(artifact_paths), expected_paths)) {
        stop(
          "completed ", manifest$analysis,
          " manifest does not record its exact output set"
        )
      }
    } else {
      expected_paths <- c(
        "ECODA_authors_HR_cell_subsetting.rds",
        "ECODA_authors_HR_cell_subsetting.csv"
      )
      if (!identical(artifact_datasets, c("ALL", "ALL")) ||
          !identical(
            artifact_methods,
            c(
              "ECODA_authors_HR_cell_subsetting",
              "ECODA_authors_HR_cell_subsetting"
            )
          ) ||
          !identical(basename(artifact_paths), expected_paths)) {
        stop("completed cell_subsetting manifest does not record its exact output set")
      }
    }
  }
  if ("source_manifest_path" %in% names(manifest) &&
      !is.null(manifest$source_manifest_path) &&
      (!is.character(manifest$source_manifest_path) ||
       length(manifest$source_manifest_path) != 1L ||
       !nzchar(manifest$source_manifest_path))) {
    stop("derived source_manifest_path is invalid")
  }
  if ("source_manifest_md5" %in% names(manifest) &&
      !is.null(manifest$source_manifest_md5) &&
      (!is.character(manifest$source_manifest_md5) ||
       length(manifest$source_manifest_md5) != 1L ||
       !grepl("^[0-9a-f]{32}$", tolower(manifest$source_manifest_md5), perl = TRUE))) {
    stop("derived source_manifest_md5 is invalid")
  }
  TRUE
}

# Explicit internal mapping: this is the only place the new standalone method
# keys are registered.  It does not touch the ordinary Stage 5 method map.
ECODA_DERIVED_METHODS <- list(
  res50 = list(
    ECODA_seuratres_50 = list(
      display = "ECODA_Leiden_res_50",
      obs_col = "leiden_res_50_benchmark_analysis_hvg2000",
      embedding = "X_pca_benchmark_analysis_hvg2000"
    )
  ),
  harmony = list(
    ECODA_seuratres_2_harmony = list(
      display = "ECODA_Leiden_res_2_harmony",
      obs_col = "leiden_res_2_benchmark_analysis_hvg2000_harmony",
      embedding = "X_pca_harmony_benchmark_analysis_hvg2000"
    )
  ),
  cell_subsetting = list(
    ECODA_authors_HR_cell_subsetting = list(
      display = "ECODA_authors_HR_cell_subsetting"
    )
  )
)

.ecoda_runner_path <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  token <- full[startsWith(full, "--file=")]
  if (length(token)) {
    return(normalizePath(sub("^--file=", "", token[[1L]]), mustWork = TRUE))
  }
  candidate <- file.path(getwd(), "src", "5_run_benchmark_methods", "run_local_ecoda_derived.R")
  if (file.exists(candidate)) normalizePath(candidate, mustWork = TRUE) else NA_character_
}

.ecoda_md5 <- function(path) {
  if (!file.exists(path)) stop("Cannot checksum missing path: ", path)
  digest <- unname(tools::md5sum(path))
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[0-9a-fA-F]{32}$", digest, perl = TRUE)) {
    stop("Could not compute MD5 for: ", path)
  }
  tolower(as.character(digest))
}

.ecoda_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (must_work) return(normalizePath(path, mustWork = TRUE))
  if (length(path) != 1L) {
    return(vapply(
      path,
      .ecoda_path,
      character(1),
      must_work = FALSE
    ))
  }

  # normalizePath(mustWork = FALSE) cannot resolve symlinks above a
  # nonexistent descendant.  Canonicalize the deepest existing ancestor
  # first, then attach the unresolved components to that canonical parent.
  candidate <- path
  suffix <- character()
  repeat {
    if (file.exists(candidate) || dir.exists(candidate)) {
      canonical_parent <- normalizePath(candidate, mustWork = TRUE)
      if (!length(suffix)) return(canonical_parent)

      # Walk the suffix from the canonical parent so existing symlinks are
      # resolved before any following .. component is applied.  Once a
      # component is missing, descendants remain unresolved until a .. pops
      # back to the existing prefix.
      current <- canonical_parent
      unresolved <- character()
      for (component in suffix) {
        if (!nzchar(component) || identical(component, ".")) next
        if (identical(component, "..")) {
          if (length(unresolved)) {
            unresolved <- unresolved[-length(unresolved)]
          } else {
            current <- dirname(current)
          }
          next
        }
        if (length(unresolved)) {
          unresolved <- c(unresolved, component)
          next
        }
        candidate <- file.path(current, component)
        if (file.exists(candidate) || dir.exists(candidate)) {
          current <- normalizePath(candidate, mustWork = TRUE)
        } else {
          unresolved <- component
        }
      }
      if (!length(unresolved)) return(current)
      rebuilt <- do.call(file.path, c(list(current), as.list(unresolved)))
      return(normalizePath(rebuilt, mustWork = FALSE))
    }

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    suffix <- c(basename(candidate), suffix)
    candidate <- parent
  }
  normalizePath(path, mustWork = FALSE)
}

.ecoda_path_within <- function(path, root) {
  path <- .ecoda_path(path)
  root <- .ecoda_path(root)
  boundary <- if (identical(root, .Platform$file.sep)) {
    root
  } else {
    paste0(root, .Platform$file.sep)
  }
  identical(path, root) || startsWith(path, boundary)
}

.ecoda_derived_root <- function(project_root) {
  project_root <- .ecoda_path(project_root, must_work = TRUE)
  configured <- Sys.getenv("ECODA_DERIVED_ROOT", unset = "")
  if (!nzchar(configured)) {
    configured <- file.path(
      project_root, "data", "benchmark", "results", "derived"
    )
  }
  derived_root <- .ecoda_path(configured)
  if (file.exists(derived_root) && !dir.exists(derived_root)) {
    stop("ECODA_DERIVED_ROOT must name a directory: ", derived_root)
  }
  derived_root
}

.ecoda_validate_derived_output_dir <- function(output_dir, derived_root) {
  output_dir <- .ecoda_path(output_dir)
  derived_root <- .ecoda_path(derived_root)
  if (identical(output_dir, derived_root) ||
      !.ecoda_path_within(output_dir, derived_root)) {
    stop(
      "output_dir must be a strict descendant of the derived output root: ",
      derived_root
    )
  }
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stop("output_dir must name a directory: ", output_dir)
  }
  list(output_dir = output_dir, derived_root = derived_root)
}

.ecoda_validate_derived_output_for_args <- function(output_arg) {
  runner_path <- .ecoda_runner_path()
  if (is.na(runner_path) || !file.exists(runner_path)) {
    stop("runner path is unavailable")
  }
  runner_path <- .ecoda_path(runner_path, must_work = TRUE)
  project_root <- dirname(dirname(dirname(runner_path)))
  derived_root <- .ecoda_derived_root(project_root)
  validated <- .ecoda_validate_derived_output_dir(output_arg, derived_root)
  c(validated, list(project_root = project_root, runner_path = runner_path))
}


.ecoda_validate_path_component <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || value %in% c(".", "..") ||
      grepl("^[A-Za-z]:", value, perl = TRUE) ||
      grepl("[/\\\\]", value, perl = TRUE) ||
      grepl("[[:cntrl:]]", value, perl = TRUE)) {
    stop(label, " is invalid")
  }
  value
}

.ecoda_validate_run_id <- function(run_id) {
  if (!is.character(run_id) || length(run_id) != 1L ||
      is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id, perl = TRUE)) {
    stop("run_id is invalid")
  }
  run_id
}

.ecoda_run_owner_path <- function(output_dir) {
  .ecoda_output_path(
    file.path(output_dir, ".ecoda_run_owner"),
    output_dir
  )
}

.ecoda_is_symlink <- function(path) {
  isTRUE(nzchar(Sys.readlink(path), keepNA = TRUE))
}

.ecoda_write_run_owner <- function(output_dir, run_id) {
  run_id <- .ecoda_validate_run_id(run_id)
  path <- .ecoda_run_owner_path(output_dir)
  if (file.exists(path) || .ecoda_is_symlink(path)) {
    stop("run ownership marker already exists: ", path)
  }
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(paste0("RUN_ID=", run_id), temporary, useBytes = TRUE)
  if (file.exists(path) || .ecoda_is_symlink(path) ||
      !file.rename(temporary, path)) {
    stop("Could not atomically publish run ownership marker: ", path)
  }
  invisible(path)
}

.ecoda_run_owner_matches <- function(output_dir, run_id) {
  run_id <- tryCatch(
    .ecoda_validate_run_id(run_id),
    error = function(error) NULL
  )
  if (is.null(run_id)) return(FALSE)
  path <- tryCatch(
    .ecoda_run_owner_path(output_dir),
    error = function(error) NULL
  )
  if (is.null(path) || !file.exists(path) ||
      isTRUE(file.info(path)$isdir) || .ecoda_is_symlink(path)) {
    return(FALSE)
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(error) NULL)
  if (is.null(lines) || length(lines) != 1L ||
      !grepl("^RUN_ID=[A-Za-z0-9][A-Za-z0-9_-]*$", lines[[1L]],
             perl = TRUE)) {
    return(FALSE)
  }
  identical(sub("^RUN_ID=", "", lines[[1L]]), run_id)
}

.ecoda_is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path, perl = TRUE)
}

.ecoda_output_path <- function(path, output_dir) {
  path <- .ecoda_path(path)
  output_dir <- .ecoda_path(output_dir)
  prefix <- paste0(output_dir, .Platform$file.sep)
  if (!identical(path, output_dir) && !startsWith(path, prefix)) {
    stop("derived output escapes output_dir: ", path)
  }
  path
}

.ecoda_write_text_atomic <- function(path, text) {
  path <- .ecoda_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(as.character(text), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Could not atomically publish: ", path)
  invisible(path)
}

.ecoda_write_json_atomic <- function(value, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to publish derived manifests")
  }
  text <- jsonlite::toJSON(
    value, auto_unbox = TRUE, pretty = TRUE, null = "null",
    na = "null", dataframe = "rows"
  )
  .ecoda_write_text_atomic(path, paste0(text, "\n"))
}

.ecoda_validate_source_sidecar <- function(path) {
  path <- .ecoda_path(path, must_work = TRUE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) stop("H5AD is empty: ", path)
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || file.info(sidecar)$size <= 0) {
    stop("H5AD checksum sidecar is missing: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("H5AD checksum sidecar has the wrong schema: ", sidecar)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[0-9a-fA-F]{32}$", md5, perl = TRUE) ||
      !grepl("^[0-9]+$", size, perl = TRUE) ||
      as.character(info$size) != size ||
      !identical(.ecoda_path(recorded_path), path)) {
    stop("H5AD checksum sidecar metadata mismatch: ", sidecar)
  }
  actual <- .ecoda_md5(path)
  if (!identical(tolower(md5), actual)) stop("H5AD checksum mismatch: ", path)
  list(MD5 = actual, SIZE = as.character(info$size), PATH = path)
}
.ecoda_validate_output_sidecar <- function(path, require_path = TRUE) {
  path <- .ecoda_path(path, must_work = TRUE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) {
    stop("Derived output is missing or empty: ", path)
  }
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || is.na(file.info(sidecar)$size) ||
      file.info(sidecar)$size <= 0) {
    stop("Derived output checksum sidecar is missing: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("Derived output checksum sidecar has the wrong schema: ", sidecar)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[0-9a-fA-F]{32}$", md5, perl = TRUE) ||
      !grepl("^[0-9]+$", size, perl = TRUE) || !nzchar(recorded_path) ||
      !identical(size, as.character(info$size))) {
    stop("Derived output checksum sidecar metadata mismatch: ", sidecar)
  }
  path_matches <- if (isTRUE(require_path)) {
    canonical_recorded_path <- tryCatch(
      .ecoda_path(recorded_path),
      error = function(error) NA_character_
    )
    !is.na(canonical_recorded_path) &&
      identical(canonical_recorded_path, path)
  } else {
    identical(basename(recorded_path), basename(path))
  }
  if (!path_matches) {
    stop("Derived output checksum sidecar metadata mismatch: ", sidecar)
  }
  actual <- .ecoda_md5(path)
  if (!identical(tolower(md5), actual)) {
    stop("Derived output checksum mismatch: ", path)
  }
  list(MD5 = actual, SIZE = as.character(info$size), PATH = path)
}
# Validate and load the compact composition snapshot.  This path intentionally
# verifies the manifest and RDS sidecar before readRDS(); source H5AD paths in
# the manifest are provenance only and are never touched by this consumer.
.ecoda_validate_composition_snapshot <- function(
  snapshot_arg,
  config_path,
  config_md5,
  scope,
  selected
) {
  if (!is.character(snapshot_arg) || length(snapshot_arg) != 1L ||
      is.na(snapshot_arg) || !nzchar(snapshot_arg) ||
      identical(snapshot_arg, TRUE)) {
    stop("--composition_snapshot must be one non-empty RDS path")
  }
  snapshot_path <- .ecoda_path(snapshot_arg, must_work = TRUE)
  if (isTRUE(file.info(snapshot_path)$isdir)) {
    stop("--composition_snapshot must name an RDS file: ", snapshot_path)
  }
  snapshot_identity <- .ecoda_validate_output_sidecar(
    snapshot_path, require_path = FALSE
  )
  manifest_path <- .ecoda_path(
    file.path(dirname(snapshot_path), "composition_snapshot_manifest.json"),
    must_work = TRUE
  )
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read composition snapshot manifests")
  }
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  if (!is.list(manifest)) stop("composition snapshot manifest must be an object")

  required_manifest <- c(
    "schema_version", "run_id", "scope", "status",
    "config_path", "config_md5", "extractor_path", "extractor_md5",
    "snapshot_path", "snapshot_md5", "snapshot_size", "sources"
  )
  missing_manifest <- setdiff(required_manifest, names(manifest))
  if (length(missing_manifest)) {
    stop(
      "composition snapshot manifest is missing: ",
      paste(missing_manifest, collapse = ", ")
    )
  }
  schema_version <- manifest$schema_version
  if ((!is.character(schema_version) && !is.numeric(schema_version)) ||
      length(schema_version) != 1L || is.na(schema_version) ||
      !nzchar(as.character(schema_version))) {
    stop("composition snapshot schema_version is invalid")
  }
  .ecoda_validate_run_id(as.character(manifest$run_id))
  if (!identical(manifest$scope, scope)) {
    stop("composition snapshot scope does not match --scope")
  }
  if (!identical(manifest$status, "COMPLETED")) {
    stop("composition snapshot manifest is not COMPLETED")
  }
  if (!is.character(manifest$config_path) || length(manifest$config_path) != 1L ||
      is.na(manifest$config_path) || !nzchar(manifest$config_path)) {
    stop("composition snapshot config_path provenance is invalid")
  }
  if (!is.character(manifest$config_md5) || length(manifest$config_md5) != 1L ||
      !grepl("^[0-9a-f]{32}$", tolower(manifest$config_md5), perl = TRUE) ||
      !identical(tolower(manifest$config_md5), tolower(config_md5))) {
    stop("composition snapshot config_md5 does not match --config_path")
  }
  if (!is.character(manifest$extractor_path) ||
      length(manifest$extractor_path) != 1L ||
      !nzchar(manifest$extractor_path) ||
      !is.character(manifest$extractor_md5) ||
      length(manifest$extractor_md5) != 1L ||
      !grepl("^[0-9a-f]{32}$", tolower(manifest$extractor_md5), perl = TRUE)) {
    stop("composition snapshot extractor provenance is invalid")
  }
  snapshot_size <- suppressWarnings(as.numeric(manifest$snapshot_size))
  if (!is.character(manifest$snapshot_path) ||
      length(manifest$snapshot_path) != 1L ||
      is.na(manifest$snapshot_path) ||
      !nzchar(manifest$snapshot_path) ||
      !identical(
        basename(manifest$snapshot_path),
        basename(snapshot_path)
      ) ||
      !is.character(manifest$snapshot_md5) ||
      length(manifest$snapshot_md5) != 1L ||
      !identical(tolower(manifest$snapshot_md5), snapshot_identity$MD5) ||
      length(snapshot_size) != 1L ||
      !is.finite(snapshot_size) ||
      snapshot_size <= 0 ||
      !identical(as.character(as.integer(snapshot_size)),
                 snapshot_identity$SIZE)) {
    stop("composition snapshot manifest does not match its RDS sidecar")
  }

  if (!is.list(manifest$sources) || !length(manifest$sources)) {
    stop("composition snapshot manifest sources must be a nonempty array")
  }
  source_datasets <- vapply(manifest$sources, function(source) {
    if (!is.list(source)) stop("composition snapshot source entry is malformed")
    needed <- c(
      "dataset", "view", "h5ad_path", "h5ad_md5", "h5ad_size", "h5ad_mtime",
      "label_col", "high_res_col"
    )
    missing <- setdiff(needed, names(source))
    if (length(missing)) {
      stop(
        "composition snapshot source entry is missing: ",
        paste(missing, collapse = ", ")
      )
    }
    if (!is.character(source$dataset) || length(source$dataset) != 1L ||
        !nzchar(source$dataset) ||
        !identical(source$view, "benchmark_analysis") ||
        !is.character(source$h5ad_path) || length(source$h5ad_path) != 1L ||
        !nzchar(source$h5ad_path) || !.ecoda_is_absolute_path(source$h5ad_path) ||
        !is.character(source$h5ad_md5) || length(source$h5ad_md5) != 1L ||
        !grepl("^[0-9a-f]{32}$", tolower(source$h5ad_md5), perl = TRUE) ||
        !is.numeric(source$h5ad_size) || length(source$h5ad_size) != 1L ||
        !is.finite(source$h5ad_size) || source$h5ad_size <= 0 ||
        ((!is.numeric(source$h5ad_mtime) && !is.character(source$h5ad_mtime)) ||
         length(source$h5ad_mtime) != 1L || is.na(source$h5ad_mtime) ||
         !nzchar(as.character(source$h5ad_mtime))) ||
        !is.character(source$label_col) ||
        length(source$label_col) != 1L || !nzchar(source$label_col) ||
        !is.character(source$high_res_col) ||
        length(source$high_res_col) != 1L || !nzchar(source$high_res_col)) {
      stop("composition snapshot source entry is invalid")
    }
    source$dataset
  }, character(1))
  if (anyDuplicated(source_datasets) ||
      !identical(sort(source_datasets), sort(as.character(selected)))) {
    stop("composition snapshot sources do not match the configured dataset union")
  }

  snapshot <- readRDS(snapshot_path)
  if (!is.list(snapshot) || !all(c("labels", "counts") %in% names(snapshot))) {
    stop("composition snapshot RDS must contain labels and counts tables")
  }
  labels <- as.data.frame(snapshot$labels, stringsAsFactors = FALSE)
  counts <- as.data.frame(snapshot$counts, stringsAsFactors = FALSE)
  label_columns <- c(
    "dataset", "Sample", "label", "total_cells", "annotated_cells",
    "original_cells"
  )
  count_columns <- c(
    "dataset", "analysis", "target", "replicate", "seed", "Sample",
    "category", "count", "total_cells", "annotated_cells",
    "original_cells", "effective_cells"
  )
  if (!all(label_columns %in% colnames(labels))) {
    stop("composition snapshot labels table has the wrong schema")
  }
  if (!all(count_columns %in% colnames(counts))) {
    stop("composition snapshot counts table has the wrong schema")
  }
  if (!nrow(labels) || !nrow(counts)) {
    stop("composition snapshot labels/counts tables must be nonempty")
  }
  if (anyNA(labels$dataset) || anyNA(labels$Sample) ||
      any(!nzchar(as.character(labels$dataset))) ||
      any(!nzchar(as.character(labels$Sample))) ||
      anyNA(labels$label) || any(.ecoda_is_missing_sentinel(labels$label)) ||
      anyNA(labels$total_cells) || anyNA(labels$annotated_cells) ||
      anyNA(labels$original_cells)) {
    stop("composition snapshot labels table contains invalid values")
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
      any(label_original != label_annotated)) {
    stop("composition snapshot labels table contains invalid cell diagnostics")
  }
  if (anyDuplicated(paste(labels$dataset, labels$Sample, sep = "\r")) ||
      !identical(sort(unique(as.character(labels$dataset))),
                 sort(as.character(selected)))) {
    stop("composition snapshot labels are not dataset-complete and unique")
  }
  if (anyNA(counts$dataset) || anyNA(counts$analysis) ||
      anyNA(counts$target) || anyNA(counts$Sample) ||
      anyNA(counts$category) ||
      any(!nzchar(as.character(counts$dataset))) ||
      any(!nzchar(as.character(counts$Sample))) ||
      any(!nzchar(as.character(counts$category))) ||
      any(!as.character(counts$analysis) %in%
          c("res50", "harmony", "cell_subsetting")) ||
      anyNA(counts$replicate) ||
      any(!is.finite(as.numeric(counts$replicate))) ||
      any(as.numeric(counts$replicate) != as.integer(counts$replicate)) ||
      anyNA(counts$count) || anyNA(counts$total_cells) ||
      anyNA(counts$annotated_cells) || anyNA(counts$original_cells) ||
      anyNA(counts$effective_cells) ||
      any(
        tolower(trimws(as.character(counts$category))) == "unassigned" &
          as.character(counts$analysis) == "cell_subsetting"
      )) {
    stop("composition snapshot counts table contains invalid values")
  }
  count_values <- lapply(
    counts[c(
      "count", "total_cells", "annotated_cells", "original_cells",
      "effective_cells"
    )],
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
      any(!as.character(counts$dataset) %in% as.character(selected))) {
    stop("composition snapshot counts table contains invalid cell diagnostics")
  }
  duplicate_key <- paste(
    counts$dataset, counts$analysis, counts$target, counts$replicate,
    ifelse(is.na(counts$seed), "<NA>", counts$seed),
    counts$Sample, counts$category, sep = "\r"
  )
  if (anyDuplicated(duplicate_key)) {
    stop("composition snapshot counts contain duplicate rows")
  }

  for (dataset in selected) {
    ds_labels <- labels[labels$dataset == dataset, , drop = FALSE]
    ds_samples <- as.character(ds_labels$Sample)
    for (analysis_name in c("res50", "harmony", "cell_subsetting")) {
      ds_counts <- counts[
        counts$dataset == dataset & counts$analysis == analysis_name,
        , drop = FALSE
      ]
      if (!nrow(ds_counts) ||
          !identical(sort(unique(as.character(ds_counts$Sample))),
                     sort(ds_samples))) {
        stop(
          "composition snapshot counts are not sample-complete for ",
          dataset, "/", analysis_name
        )
      }
      for (sample_id in ds_samples) {
        sample_rows <- ds_counts[ds_counts$Sample == sample_id, , drop = FALSE]
        total <- unique(as.numeric(sample_rows$total_cells))
        annotated <- unique(as.numeric(sample_rows$annotated_cells))
        original <- unique(as.numeric(sample_rows$original_cells))
        label_row <- ds_labels[ds_labels$Sample == sample_id, , drop = FALSE]
        label_total <- as.numeric(label_row$total_cells)
        label_annotated <- as.numeric(label_row$annotated_cells)
        if (length(total) != 1L || length(annotated) != 1L ||
            length(original) != 1L ||
            !identical(total, label_total) ||
            !identical(annotated, label_annotated) ||
            (analysis_name %in% c("res50", "harmony") &&
             !identical(original, total)) ||
            (identical(analysis_name, "cell_subsetting") &&
             !identical(original, annotated))) {
          stop(
            "composition snapshot counts have inconsistent cell diagnostics for ",
            dataset, "/", analysis_name, "/", sample_id
          )
        }
        if (identical(analysis_name, "cell_subsetting")) {
          keys <- unique(sample_rows[c("target", "replicate", "seed")])
          for (key_row in seq_len(nrow(keys))) {
            key_mask <- as.character(sample_rows$target) ==
                as.character(keys$target[[key_row]]) &
              as.integer(sample_rows$replicate) ==
                as.integer(keys$replicate[[key_row]])
            key_seed <- keys$seed[[key_row]]
            if (is.na(key_seed)) {
              key_mask <- key_mask & is.na(sample_rows$seed)
            } else {
              key_mask <- key_mask & !is.na(sample_rows$seed) &
                as.integer(sample_rows$seed) == as.integer(key_seed)
            }
            key_values <- sample_rows[key_mask, , drop = FALSE]
            effective <- unique(as.numeric(key_values$effective_cells))
            target_value <- as.character(keys$target[[key_row]])
            expected_effective <- if (identical(target_value, "all cells")) {
              original
            } else {
              suppressWarnings(pmin(original, as.numeric(target_value)))
            }
            if (length(effective) != 1L || length(expected_effective) != 1L ||
                is.na(expected_effective) ||
                !identical(effective, expected_effective) ||
                sum(as.numeric(key_values$count)) != effective) {
              stop(
                "composition snapshot counts have inconsistent cell diagnostics for ",
                dataset, "/", analysis_name, "/", sample_id
              )
            }
          }
        } else {
          effective <- unique(as.numeric(sample_rows$effective_cells))
          if (length(effective) != 1L ||
              !identical(effective, total) ||
              sum(as.numeric(sample_rows$count)) != effective) {
            stop(
              "composition snapshot counts have inconsistent cell diagnostics for ",
              dataset, "/", analysis_name, "/", sample_id
            )
          }
        }
      }
    }
  }

  # Validate the exact all-cell and cell-depth key schedule without touching
  # the compact table's cell-level representation.
  for (dataset in selected) {
    ds_labels <- labels[labels$dataset == dataset, , drop = FALSE]
    ds_counts <- counts[counts$dataset == dataset, , drop = FALSE]
    for (analysis_name in c("res50", "harmony")) {
      rows <- ds_counts[ds_counts$analysis == analysis_name, , drop = FALSE]
      if (any(as.character(rows$target) != "all cells") ||
          any(as.integer(rows$replicate) != 0L) ||
          any(!is.na(rows$seed))) {
        stop("composition snapshot has invalid all-cell keys for ", dataset)
      }
    }
    rows <- ds_counts[ds_counts$analysis == "cell_subsetting", , drop = FALSE]
    expected <- ecoda_subsetting_plan(
      unique(as.character(ds_labels$Sample)),
      setNames(as.integer(ds_labels$original_cells), ds_labels$Sample)
    )
    observed <- unique(rows[c("target", "replicate", "seed")])
    if (nrow(observed) != nrow(expected) ||
        !all(vapply(seq_len(nrow(expected)), function(i) {
          same_seed <- if (is.na(expected$seed[[i]])) {
            is.na(observed$seed)
          } else {
            observed$seed == expected$seed[[i]]
          }
          any(
            observed$target == expected$target[[i]] &
              as.integer(observed$replicate) == expected$replicate[[i]] &
              same_seed
          )
        }, logical(1))) ||
        any(as.integer(rows$seed[!is.na(rows$seed)]) < 101L) ||
        any(as.integer(rows$seed[!is.na(rows$seed)]) > 120L)) {
      stop("composition snapshot cell-subsetting key schedule is incomplete")
    }
  }
  list(
    path = snapshot_path,
    md5 = snapshot_identity$MD5,
    manifest_path = manifest_path,
    manifest = manifest,
    snapshot_config_path = manifest$config_path,
    snapshot_config_md5 = tolower(manifest$config_md5),
    labels = labels,
    counts = counts,
    metadata = snapshot$metadata %||% list()
  )
}

.ecoda_snapshot_labels <- function(labels, dataset) {
  rows <- labels[labels$dataset == dataset, , drop = FALSE]
  sample_ids <- as.character(rows$Sample)
  result <- as.factor(as.character(rows$label))
  names(result) <- sample_ids
  result
}

.ecoda_snapshot_counts_matrix <- function(rows, sample_ids) {
  rows <- rows[order(match(as.character(rows$Sample), sample_ids)), , drop = FALSE]
  categories <- unique(as.character(rows$category))
  has_subset_rows <- !"analysis" %in% colnames(rows) ||
    any(as.character(rows$analysis) == "cell_subsetting")
  if (!length(categories) ||
      (has_subset_rows && any(tolower(trimws(categories)) == "unassigned"))) {
    stop("composition snapshot has invalid cell-type categories")
  }
  matrix_counts <- matrix(
    0,
    nrow = length(sample_ids),
    ncol = length(categories),
    dimnames = list(sample_ids, categories)
  )
  for (row_index in seq_len(nrow(rows))) {
    sample_id <- as.character(rows$Sample[[row_index]])
    category <- as.character(rows$category[[row_index]])
    matrix_counts[sample_id, category] <-
      matrix_counts[sample_id, category] +
      as.numeric(rows$count[[row_index]])
  }
  keep_categories <- colSums(matrix_counts) > 0
  if (!any(keep_categories)) {
    stop("composition snapshot count matrix has no observed categories")
  }
  matrix_counts <- matrix_counts[, keep_categories, drop = FALSE]
  if (any(!is.finite(matrix_counts)) || any(rowSums(matrix_counts) <= 0)) {
    stop("composition snapshot count matrix is nonpositive or nonfinite")
  }
  as.data.frame(matrix_counts, check.names = FALSE, stringsAsFactors = FALSE)
}

.ecoda_snapshot_bundle <- function(count_rows, labels, expected_samples, method) {
  df_counts <- .ecoda_snapshot_counts_matrix(count_rows, expected_samples)
  df_imp <- impute_zeros(
    df_counts,
    clr_zero_impute_method = "counts_all",
    clr_zero_impute_num = 0.5
  )
  feat_mat <- clr(df_imp)
  bundle <- create_result_bundle(
    feat_mat,
    labels,
    dist_mat = dist(feat_mat),
    extra = list(counts = df_imp)
  )
  .ecoda_validate_bundle(bundle, expected_samples, method)
  bundle
}



.ecoda_parse_flags <- function(raw_args) {
  args <- list()
  i <- 1L
  while (i <= length(raw_args)) {
    flag <- raw_args[[i]]
    if (!startsWith(flag, "--")) stop("Unexpected positional argument: ", flag)
    name <- sub("^--", "", flag)
    if (grepl("=", name, fixed = TRUE)) {
      kv <- strsplit(name, "=", fixed = TRUE)[[1L]]
      if (length(kv) != 2L || !nzchar(kv[[1L]]) || kv[[1L]] %in% names(args)) {
        stop("Malformed or repeated argument: ", flag)
      }
      args[[kv[[1L]]]] <- kv[[2L]]
      i <- i + 1L
    } else if (i < length(raw_args) && !startsWith(raw_args[[i + 1L]], "--")) {
      if (name %in% names(args)) stop("Repeated argument: --", name)
      args[[name]] <- raw_args[[i + 1L]]
      i <- i + 2L
    } else {
      if (name %in% names(args)) stop("Repeated argument: --", name)
      args[[name]] <- TRUE
      i <- i + 1L
    }
  }
  args
}

.ecoda_resolve_source <- function(input_dir, dataset, view_spec) {
  input_dir <- .ecoda_path(input_dir)
  dataset <- .ecoda_validate_path_component(
    dataset,
    "benchmark_analysis dataset key"
  )
  output_file <- view_spec$output_file_name %||% view_spec$output_file
  if (!is.character(output_file) || length(output_file) != 1L ||
      is.na(output_file) || !nzchar(output_file) ||
      output_file %in% c(".", "..") ||
      .ecoda_is_absolute_path(output_file) ||
      grepl("^[A-Za-z]:", output_file, perl = TRUE) ||
      grepl("[\\\\]", output_file, perl = TRUE) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", output_file, perl = TRUE) ||
      grepl("[[:cntrl:]]", output_file, perl = TRUE)) {
    stop("benchmark_analysis output file is invalid")
  }
  candidates <- unique(c(
    file.path(input_dir, dataset, "output", output_file),
    file.path(input_dir, output_file),
    file.path(input_dir, dataset, output_file)
  ))
  canonical_candidates <- unique(vapply(
    candidates,
    function(path) .ecoda_path(path),
    character(1)
  ))
  if (any(!vapply(
    canonical_candidates,
    function(path) .ecoda_path_within(path, input_dir),
    logical(1)
  ))) {
    stop("declared benchmark_analysis source escapes input_dir")
  }
  existing <- canonical_candidates[file.exists(canonical_candidates)]
  if (!length(existing)) stop("declared benchmark_analysis H5AD is missing")
  if (length(existing) > 1L) {
    # The canonical dataset/output layout wins.  Never silently choose between
    # two different existing source bytes in the same run.
    canonical <- existing[[1L]]
    hashes <- vapply(existing, function(p) .ecoda_validate_source_sidecar(p)$MD5, character(1))
    if (length(unique(hashes)) > 1L) stop("ambiguous benchmark_analysis H5AD paths")
    return(.ecoda_path(canonical, must_work = TRUE))
  }
  .ecoda_path(existing[[1L]], must_work = TRUE)
}

.ecoda_as_obs <- function(value) {
  if (!is.data.frame(value)) value <- as.data.frame(value, stringsAsFactors = FALSE)
  for (name in colnames(value)) {
    if (is.factor(value[[name]])) value[[name]] <- as.character(value[[name]])
  }
  value
}

.ecoda_check_obs <- function(obs, entry, analysis, required_columns) {
  obs <- .ecoda_as_obs(obs)
  missing <- setdiff(required_columns, colnames(obs))
  if (length(missing)) stop("H5AD obs is missing: ", paste(missing, collapse = ", "))
  sample_values <- as.character(obs[["Sample"]])
  if (!length(sample_values) || any(.ecoda_is_missing_sentinel(sample_values))) {
    stop("H5AD Sample contains missing or blank values")
  }
  sample_ids <- unique(sample_values)
  total_cells <- tabulate(
    match(sample_values, sample_ids),
    nbins = length(sample_ids)
  )
  names(total_cells) <- sample_ids
  label_col <- entry$columns$label %||% entry$label_col
  label_values <- as.character(obs[[label_col]])
  if (any(.ecoda_is_missing_sentinel(label_values))) {
    stop("H5AD biological label column '", label_col,
         "' has missing or blank values")
  }
  for (sample_id in sample_ids) {
    values <- label_values[sample_values == sample_id]
    if (length(unique(values)) != 1L) {
      stop("H5AD biological labels for sample '", sample_id, "' are conflicting")
    }
  }
  if (analysis == "cell_subsetting") {
    ct_cols <- entry$columns$cell_type_high_res %||% entry$cell_type_high_res
  } else {
    method <- names(ECODA_DERIVED_METHODS[[analysis]])[[1L]]
    ct_cols <- ECODA_DERIVED_METHODS[[analysis]][[method]]$obs_col
  }
  for (column in unique(ct_cols)) {
    values <- as.character(obs[[column]])
    missing_mask <- if (analysis == "cell_subsetting") {
      .ecoda_is_missing_high_res(values)
    } else {
      .ecoda_is_missing_sentinel(values)
    }
    if (analysis == "cell_subsetting") {
      has_nonmissing <- vapply(sample_ids, function(sample_id) {
        any(!missing_mask[sample_values == sample_id])
      }, logical(1L))
      if (any(!has_nonmissing)) {
        stop(
          "H5AD obs column '", column,
          "' has no annotated cells for sample(s): ",
          paste(sample_ids[!has_nonmissing], collapse = ", ")
        )
      }
      keep <- !missing_mask
      obs <- obs[keep, , drop = FALSE]
      sample_values <- sample_values[keep]
      ordering <- order(match(sample_values, sample_ids))
      obs <- obs[ordering, , drop = FALSE]
      sample_values <- sample_values[ordering]
      annotated_cells <- tabulate(
        match(sample_values, sample_ids),
        nbins = length(sample_ids)
      )
      names(annotated_cells) <- sample_ids
      if (any(annotated_cells <= 0L)) {
        stop(
          "H5AD high-resolution filtering removed every annotated cell for sample(s): ",
          paste(sample_ids[annotated_cells <= 0L], collapse = ", ")
        )
      }
      attr(obs, "ecoda_total_cells") <-
        setNames(as.integer(total_cells), sample_ids)
      attr(obs, "ecoda_annotated_cells") <-
        setNames(as.integer(annotated_cells), sample_ids)
    } else if (any(missing_mask)) {
      stop("H5AD obs column contains invalid cell types: ", column)
    }
  }
  label_col <- entry$columns$label %||% entry$label_col
  label_values <- as.character(obs[[label_col]])
  if (any(.ecoda_is_missing_sentinel(label_values))) {
    stop("H5AD biological label column '", label_col,
         "' has missing or blank values")
  }
  invisible(obs)
}

.ecoda_labels <- function(obs, label_col) {
  if (!is.character(label_col) || length(label_col) != 1L ||
      is.na(label_col) || !nzchar(label_col) ||
      !"Sample" %in% colnames(obs) || !label_col %in% colnames(obs)) {
    stop("biological label column or Sample is missing")
  }
  sample_values <- as.character(obs[["Sample"]])
  label_values <- as.character(obs[[label_col]])
  if (length(sample_values) != length(label_values)) {
    stop("biological labels are not aligned with Sample")
  }
  sample_ids <- unique(sample_values)
  labels <- vapply(sample_ids, function(sample_id) {
    values <- label_values[sample_values == sample_id]
    if (any(.ecoda_is_missing_sentinel(values))) {
      stop("biological labels for sample '", sample_id,
           "' are missing or blank")
    }

    if (any(values != values[[1L]])) {
      stop("biological labels for sample '", sample_id, "' are conflicting")
    }
    values[[1L]]
  }, character(1))
  labels <- as.factor(labels)
  names(labels) <- sample_ids
  labels
}
# Load and validate one derived-analysis source.  The caller owns the returned
# obs table and must release it before moving to another dataset.
# Counts-free embedding methods retain their exact embedding contract.
.ecoda_load_derived_obs <- function(
  source_config,
  analysis,
  source_identity = NULL,
  expected_md5 = NULL
) {
  h5ad_path <- .ecoda_path(source_config$h5ad_path, must_work = TRUE)
  if (is.null(source_identity)) {
    source_identity <- .ecoda_validate_source_sidecar(h5ad_path)
  }
  if (!is.null(expected_md5) &&
      !identical(source_identity$MD5, expected_md5)) {
    stop("H5AD checksum changed since source preflight: ", h5ad_path)
  }

  if (analysis == "cell_subsetting") {
    project_root <- Sys.getenv("PROJECT_ROOT")
    if (project_root == "") {
      stop("PROJECT_ROOT not set; cannot load metadata-only H5AD.")
    }
    module_dir <- normalizePath(
      file.path(project_root, "src", "utils", "py"),
      mustWork = TRUE
    )
    python_sys <- reticulate::import("sys", convert = FALSE)
    python_sys$path$insert(0L, module_dir)
    loader <- reticulate::import_from_path(
      "h5ad_obs_free",
      path = module_dir,
      convert = FALSE
    )
    obs <- reticulate::py_to_r(
      loader$load_h5ad_obs_free(
        h5ad_path,
        as.list(as.character(source_config$obs_columns))
      )
    )
  } else {
    adata <- load_h5ad_counts_free(
      h5ad_path,
      obs_columns = source_config$obs_columns,
      embedding_keys = source_config$embedding_keys,
      view = "benchmark_analysis",
      method = "ecoda_derived"
    )
    obs <- reticulate::py_to_r(adata$obs)
    rm(adata)
  }
  post_load_identity <- .ecoda_validate_source_sidecar(h5ad_path)
  if (!identical(post_load_identity$MD5, source_identity$MD5)) {
    stop("H5AD checksum changed during load: ", h5ad_path)
  }
  if (!is.null(expected_md5) &&
      !identical(post_load_identity$MD5, expected_md5)) {
    stop("H5AD checksum changed since source preflight: ", h5ad_path)
  }

  checked_obs <- .ecoda_check_obs(
    obs,
    source_config$entry,
    analysis,
    source_config$obs_columns
  )
  .ecoda_labels(checked_obs, source_config$label_col)
  checked_obs
}


.ecoda_validate_bundle <- function(bundle, expected_samples, method) {
  required <- c("feat_mat", "labels", "dist_mat", "scores", "counts")
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    stop(method, " bundle is missing required fields")
  }
  expected_samples <- as.character(expected_samples)
  feat <- as.matrix(bundle$feat_mat)
  feat_ids <- rownames(feat)
  if (!nrow(feat) || is.null(feat_ids) || anyNA(feat_ids) ||
      anyDuplicated(feat_ids) || !identical(as.character(feat_ids), expected_samples) ||
      any(!is.finite(feat))) {
    stop(method, " feature matrix is not finite or is not sample-complete/aligned")
  }
  labels <- bundle$labels
  label_ids <- names(labels)
  if (is.null(label_ids) || anyNA(label_ids) || anyDuplicated(label_ids) ||
      !identical(as.character(label_ids), expected_samples) || anyNA(labels)) {
    stop(method, " labels are not sample-complete or aligned")
  }
  dist_matrix <- as.matrix(bundle$dist_mat)
  dist_ids <- rownames(dist_matrix)
  if (nrow(dist_matrix) != length(expected_samples) ||
      ncol(dist_matrix) != length(expected_samples) ||
      is.null(dist_ids) || is.null(colnames(dist_matrix)) ||
      !identical(as.character(dist_ids), expected_samples) ||
      !identical(as.character(colnames(dist_matrix)), expected_samples) ||
      any(!is.finite(dist_matrix))) {
    stop(method, " distance matrix is not finite or sample-complete/aligned")
  }
  if (!is.list(bundle$scores) || !length(bundle$scores)) stop(method, " scores are empty")
  score_values <- unlist(bundle$scores, use.names = FALSE)
  if (!is.numeric(score_values) || any(!is.finite(score_values))) {
    stop(method, " scores contain nonfinite values")
  }
  counts <- as.data.frame(bundle$counts, stringsAsFactors = FALSE)
  if (!nrow(counts) || is.null(rownames(counts)) ||
      !identical(as.character(rownames(counts)), expected_samples) ||
      any(!vapply(counts, function(column) {
        is.numeric(column) && all(is.finite(column))
      }, logical(1)))) {
    stop(method, " counts are not sample-complete/aligned and finite")
  }
  invisible(TRUE)
}
.ecoda_validate_cell_subsetting_results <- function(result_table) {
  required <- c(
    "dataset", "target", "target_cells", "replicate", "seed", "sample_ids",
    "total_cells_per_sample", "original_cells_per_sample",
    "annotated_cells_per_sample", "effective_cells_per_sample",
    "original_total_cells", "annotated_cells", "total_cells",
    "effective_cells", "ANOSIM"
  )
  if (!is.data.frame(result_table) || !nrow(result_table) ||
      !all(required %in% colnames(result_table))) {
    stop("cell-subsetting result table is incomplete or nonfinite")
  }
  integer_fields <- c(
    "original_total_cells", "annotated_cells", "total_cells",
    "effective_cells"
  )
  for (field in integer_fields) {
    values <- suppressWarnings(as.numeric(result_table[[field]]))
    if (any(!is.finite(values)) || any(values <= 0) ||
        any(values != as.integer(values))) {
      stop("cell-subsetting result table has invalid ", field)
    }
  }
  if (any(result_table$annotated_cells > result_table$original_total_cells) ||
      any(result_table$total_cells > result_table$annotated_cells) ||
      any(result_table$effective_cells != result_table$total_cells) ||
      any(!is.finite(as.numeric(result_table$ANOSIM)))) {
    stop("cell-subsetting result table has inconsistent cell diagnostics")
  }
  for (dataset in unique(as.character(result_table$dataset))) {
    rows <- result_table[as.character(result_table$dataset) == dataset, , drop = FALSE]
    if (nrow(rows) != 181L ||
        !identical(unique(as.character(rows$target)), ecoda_derived_targets())) {
      stop("cell-subsetting result table is not exactly 181 rows for ", dataset)
    }
  }
  invisible(TRUE)
}

.ecoda_save_rds <- function(object, path) {
  # Explicit empty producer/run_id prevent the shared helper from creating an
  # artifact record under ECODA_RUNS_ROOT.  The only write is the requested
  # output RDS and its adjacent checksum sidecar.
  save_rds_atomic(object, path, producer = "", run_id = "")
  if (!isTRUE(artifact_checksum_ok(path, producer = "", run_id = ""))) {
    stop("published RDS failed checksum validation: ", path)
  }
  invisible(path)
}
.ecoda_save_csv <- function(object, path) {
  path <- .ecoda_path(path)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
      file.info(temporary)$size <= 0) {
    stop("Empty CSV temporary file: ", temporary)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not atomically publish CSV: ", path)
  }
  digest <- .ecoda_md5(path)
  size <- as.character(file.info(path)$size)
  .ecoda_write_text_atomic(
    paste0(path, ".md5"),
    c(
      paste0("MD5=", digest),
      paste0("SIZE=", size),
      paste0("PATH=", path)
    )
  )
  invisible(path)
}


.ecoda_artifact_record <- function(path, dataset, method) {
  path <- .ecoda_path(path, must_work = TRUE)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) stop("published artifact is empty: ", path)
  list(
    path = path,
    size = as.numeric(info$size),
    md5 = .ecoda_md5(path),
    dataset = dataset,
    method = method
  )
}

.ecoda_status <- function(run_id, analysis, scope, status, output_dir,
                          error = NULL, blocked = list(), derived_root) {
  list(
    run_id = run_id,
    analysis = analysis,
    scope = scope,
    status = status,
    output_dir = output_dir,
    derived_root = derived_root,
    error = error,
    blocked_datasets = blocked
  )
}


.ecoda_json_artifact <- function(path, dataset, method) {
  .ecoda_artifact_record(path, dataset = dataset, method = method)
}

.ecoda_mark_failure <- function(args, error) {
  output_arg <- args$output_dir
  if (is.null(output_arg) || identical(output_arg, TRUE) ||
      !is.character(output_arg) || length(output_arg) != 1L || is.na(output_arg)) {
    return(invisible(NULL))
  }
  validated <- tryCatch(
    .ecoda_validate_derived_output_for_args(output_arg),
    error = function(e) NULL
  )
  if (is.null(validated)) return(invisible(NULL))
  output_dir <- validated$output_dir
  derived_root <- validated$derived_root
  if (!dir.exists(output_dir)) return(invisible(NULL))
  attempted_run_id <- args$run_id %||% Sys.getenv("ECODA_RUN_ID", unset = "")
  if (!is.character(attempted_run_id) || length(attempted_run_id) != 1L ||
      is.na(attempted_run_id) || !nzchar(attempted_run_id)) {
    attempted_run_id <- basename(output_dir)
  }
  attempted_run_id <- tryCatch(
    .ecoda_validate_run_id(as.character(attempted_run_id)),
    error = function(e) NULL
  )
  if (is.null(attempted_run_id) ||
      !.ecoda_in_process_owner_matches(output_dir, attempted_run_id) ||
      !.ecoda_run_owner_matches(output_dir, attempted_run_id)) {
    return(invisible(NULL))
  }
  paths <- tryCatch({
    names <- c(
      source = "derived_source_manifest.json",
      run = "derived_run_manifest.json",
      status = "status.json"
    )
    paths <- vapply(names, function(name) {
      candidate <- file.path(output_dir, name)
      if (.ecoda_is_symlink(candidate)) {
        stop("failure status path is a symbolic link")
      }
      .ecoda_output_path(candidate, output_dir)
    }, character(1))
    paths
  }, error = function(e) NULL)
  if (is.null(paths)) return(invisible(NULL))
  source_path <- paths[["source"]]
  run_path <- paths[["run"]]
  status_path <- paths[["status"]]
  terminal <- FALSE
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    for (manifest_path in c(status_path, run_path)) {
      if (!file.exists(manifest_path)) next
      terminal <- terminal || isTRUE(tryCatch({
        status_value <- jsonlite::fromJSON(
          manifest_path, simplifyVector = FALSE
        )$status
        status_value %in% c("BLOCKED", "COMPLETED", "FAILED")
      }, error = function(e) FALSE))
    }
  }
  if (terminal) return(invisible(NULL))
  analysis <- if (is.character(args$analysis) && length(args$analysis) == 1L) {
    as.character(args$analysis)
  } else {
    "unknown"
  }
  scope <- if (is.character(args$scope) && length(args$scope) == 1L) {
    as.character(args$scope)
  } else {
    "unknown"
  }
  if (file.exists(run_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    tryCatch({
      manifest <- jsonlite::fromJSON(run_path, simplifyVector = FALSE)
      if (!identical(manifest$status, "BLOCKED") &&
          !identical(manifest$status, "COMPLETED")) {
        manifest$status <- "FAILED"
        manifest$error <- conditionMessage(error)
        .ecoda_write_json_atomic(manifest, run_path)
      }
    }, error = function(e) invisible(NULL))
  }
  tryCatch(
    .ecoda_write_json_atomic(
      .ecoda_status(
        attempted_run_id, analysis, scope, "FAILED", output_dir,
        error = conditionMessage(error), derived_root = derived_root
      ),
      status_path
    ),
    error = function(e) invisible(NULL)
  )
}
.ecoda_run <- function(args) {
  .ecoda_clear_in_process_owner()
  snapshot_mode <- !is.null(args$composition_snapshot)
  if (snapshot_mode &&
      (identical(args$composition_snapshot, TRUE) ||
       !is.character(args$composition_snapshot) ||
       length(args$composition_snapshot) != 1L ||
       is.na(args$composition_snapshot) ||
       !nzchar(args$composition_snapshot))) {
    stop("--composition_snapshot must be one non-empty RDS path")
  }
  required <- c("config_path", "output_dir", "analysis", "scope")
  if (!snapshot_mode) required <- c(required, "input_dir")
  for (name in required) {
    if (is.null(args[[name]]) || identical(args[[name]], TRUE)) {
      stop("Missing required --", name, " argument")
    }
  }
  analysis <- match.arg(
    as.character(args$analysis),
    c("res50", "harmony", "cell_subsetting")
  )
  scope <- match.arg(as.character(args$scope), c("benchmark_union", "debug"))
  if (analysis != "cell_subsetting" && !is.null(args$seeds)) {
    stop("--seeds is supported only for --analysis cell_subsetting")
  }
  seeds <- if (analysis == "cell_subsetting") ecoda_parse_seeds(args$seeds) else NULL
  if (analysis == "cell_subsetting" && !identical(seeds, 101L:120L)) {
    stop("subsetting seeds must use the exact publication schedule 101:120")
  }

  runner_path <- .ecoda_runner_path()
  if (is.na(runner_path) || !file.exists(runner_path)) {
    stop("runner path is unavailable")
  }
  runner_path <- .ecoda_path(runner_path, must_work = TRUE)
  project_root <- dirname(dirname(dirname(runner_path)))
  derived_root <- .ecoda_derived_root(project_root)

  config_path <- .ecoda_path(as.character(args$config_path), must_work = TRUE)
  input_dir <- if (snapshot_mode) {
    NULL
  } else {
    .ecoda_path(as.character(args$input_dir))
  }
  validated_output <- .ecoda_validate_derived_output_dir(
    args$output_dir, derived_root
  )
  output_dir <- validated_output$output_dir
  derived_root <- validated_output$derived_root
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) stop("Could not create output_dir: ", output_dir)
  existing_outputs <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing_outputs)) {
    stop("output_dir must be a new empty run-owned directory: ", output_dir)
  }
  run_id <- args$run_id %||% Sys.getenv("ECODA_RUN_ID", unset = "")
  if (!nzchar(run_id)) run_id <- basename(output_dir)
  run_id <- .ecoda_scalar(as.character(run_id), "run_id")
  # The output directory is the sole write root.  Any later exception is
  # converted into a terminal FAILED status by the Rscript guard below.
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id, perl = TRUE)) {
    stop("run_id is invalid")
  }
  .ecoda_validate_run_id(run_id)
  .ecoda_write_run_owner(output_dir, run_id)
  .ecoda_record_in_process_owner(output_dir, run_id)
  if (!nzchar(Sys.getenv("PROJECT_ROOT", unset = ""))) {
    Sys.setenv(PROJECT_ROOT = project_root)
  }
  oldwd <- getwd()
  tryCatch({
    setwd(project_root)
    if (snapshot_mode) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("dplyr is required for compact composition scoring")
      }
      suppressPackageStartupMessages(
        library("dplyr", character.only = TRUE)
      )
    } else {
      source(file.path(project_root, "src/utils/imports_worker_core.R"))
    }
    source(file.path(project_root, "src/utils/load_worker_functions.R"))
    source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
  }, finally = setwd(oldwd))
  config_md5 <- .ecoda_md5(config_path)
  runner_md5 <- .ecoda_md5(runner_path)
  source_manifest_path <- .ecoda_output_path(
    file.path(output_dir, "derived_source_manifest.json"), output_dir
  )
  run_manifest_path <- .ecoda_output_path(
    file.path(output_dir, "derived_run_manifest.json"), output_dir
  )
  status_path <- .ecoda_output_path(file.path(output_dir, "status.json"), output_dir)

  .ecoda_write_json_atomic(
    .ecoda_status(
      run_id, analysis, scope, "RUNNING", output_dir,
      derived_root = derived_root
    ),
    status_path
  )

  raw_config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  if (!is.list(raw_config) || is.null(names(raw_config))) stop("datasets config is not an object")
  selected <- ecoda_raw_config_union(raw_config, scope = scope)
  if (!length(selected)) {
    selected <- if (scope == "debug") "_debug" else character()
  }
  blocked <- list()
  source_rows <- list()
  source_configs <- list()
  snapshot <- NULL
  if (snapshot_mode) {
    snapshot <- .ecoda_validate_composition_snapshot(
      args$composition_snapshot,
      config_path = config_path,
      config_md5 = config_md5,
      scope = scope,
      selected = selected
    )
    snapshot_sources <- snapshot$manifest$sources
    snapshot_source_names <- vapply(
      snapshot_sources, `[[`, character(1), "dataset"
    )
    source_rows <- lapply(selected, function(dataset) {
      source <- snapshot_sources[[match(dataset, snapshot_source_names)]]
      list(
        dataset = source$dataset,
        view = source$view,
        h5ad_path = source$h5ad_path,
        h5ad_md5 = tolower(source$h5ad_md5),
        analysis = analysis,
        h5ad_size = as.numeric(source$h5ad_size),
        h5ad_mtime = source$h5ad_mtime,
        label_col = source$label_col,
        high_res_col = source$high_res_col
      )
    })
  }
  if (!length(selected)) {
    blocked[[1L]] <- list(dataset = "<selection>", view = "benchmark_analysis",
                           reason = "raw-config union is empty")
  }

  if (!snapshot_mode) {
  for (dataset in selected) {
    entry <- raw_config[[dataset]]
    block <- function(reason) {
      blocked[[length(blocked) + 1L]] <<- list(
        dataset = dataset, view = "benchmark_analysis", reason = as.character(reason)
      )
    }
    if (!is.list(entry)) {
      block("raw config entry is not an object")
      next
    }
    views <- entry$views
    if (!is.list(views) || !"benchmark_analysis" %in% names(views) ||
        !is.list(views$benchmark_analysis)) {
      block("benchmark_analysis view is missing or invalid; batch-effect views are not substitutes")
      next
    }
    view_spec <- views$benchmark_analysis
    entry$columns <- modifyList(
      entry$columns %||% list(),
      view_spec$columns %||% list()
    )
    h5ad_path <- tryCatch(
      .ecoda_resolve_source(input_dir, dataset, view_spec),
      error = function(error) {
        block(conditionMessage(error))
        NULL
      }
    )
    if (is.null(h5ad_path)) next
    source_identity <- tryCatch(
      .ecoda_validate_source_sidecar(h5ad_path),
      error = function(error) {
        block(conditionMessage(error))
        NULL
      }
    )
    if (is.null(source_identity)) next
    label_col <- entry$columns$label %||% entry$label_col
    high_col <- entry$columns$cell_type_high_res %||% entry$cell_type_high_res
    if (!is.character(label_col) || length(label_col) != 1L || !nzchar(label_col) ||
        (analysis == "cell_subsetting" &&
         (!is.character(high_col) || length(high_col) != 1L || !nzchar(high_col)))) {
      block("configured label/high-resolution cell-type columns are missing")
      next
    }
    required_obs <- if (analysis == "cell_subsetting") {
      c("Sample", label_col, high_col)
    } else {
      method <- names(ECODA_DERIVED_METHODS[[analysis]])[[1L]]
      c("Sample", label_col,
        ECODA_DERIVED_METHODS[[analysis]][[method]]$obs_col)
    }
    embedding_keys <- if (analysis == "cell_subsetting") {
      character(0)
    } else {
      method <- names(ECODA_DERIVED_METHODS[[analysis]])[[1L]]
      ECODA_DERIVED_METHODS[[analysis]][[method]]$embedding
    }
    source_config <- list(
      entry = list(
        columns = list(
          label = label_col,
          cell_type_high_res = high_col
        )
      ),
      h5ad_path = h5ad_path,
      h5ad_md5 = source_identity$MD5,
      obs_columns = unique(required_obs),
      embedding_keys = as.character(embedding_keys),
      label_col = label_col,
      high_col = high_col
    )
    preflight_obs <- tryCatch(
      .ecoda_load_derived_obs(
        source_config,
        analysis,
        source_identity = source_identity
      ),
      error = function(error) {
        block(paste0("H5AD content contract failed: ", conditionMessage(error)))
        NULL
      }
    )
    if (is.null(preflight_obs)) {
      rm(preflight_obs)
      gc(verbose = FALSE)
      next
    }
    source_rows[[length(source_rows) + 1L]] <- list(
      dataset = dataset,
      view = "benchmark_analysis",
      h5ad_path = h5ad_path,
      h5ad_md5 = source_identity$MD5,
      analysis = analysis
    )
    source_configs[[dataset]] <- source_config
    rm(preflight_obs, source_config, source_identity)
    gc(verbose = FALSE)
  }
  }

  source_manifest <- list(
    run_id = run_id,
    analysis = analysis,
    scope = scope,
    derived_root = derived_root,
    status = if (length(blocked)) "BLOCKED" else "SOURCES_VERIFIED",
    config_path = config_path,
    config_md5 = config_md5,
    runner_path = runner_path,
    runner_md5 = runner_md5,
    source_manifest_path = source_manifest_path,
    source_manifest_md5 = NULL,
    artifacts = list(),
    sources = source_rows
  )
  ecoda_validate_derived_manifest(source_manifest, required_sources = TRUE)
  if (snapshot_mode) {
    source_manifest$source_mode <- "composition_snapshot"
    source_manifest$snapshot_path <- snapshot$path
    source_manifest$snapshot_md5 <- snapshot$md5
    source_manifest$snapshot_config_path <- snapshot$snapshot_config_path
    source_manifest$snapshot_config_md5 <- snapshot$snapshot_config_md5
  }
  .ecoda_write_json_atomic(source_manifest, source_manifest_path)
  source_manifest_md5 <- .ecoda_md5(source_manifest_path)

  run_manifest <- list(
    run_id = run_id,
    analysis = analysis,
    scope = scope,
    status = if (length(blocked)) "BLOCKED" else "RUNNING",
    output_dir = output_dir,
    derived_root = derived_root,
    config_path = config_path,
    config_md5 = config_md5,
    runner_path = runner_path,
    runner_md5 = runner_md5,
    source_manifest_path = source_manifest_path,
    source_manifest_md5 = source_manifest_md5,
    artifacts = list(),
    sources = source_rows
  )
  if (snapshot_mode) {
    run_manifest$source_mode <- "composition_snapshot"
    run_manifest$snapshot_path <- snapshot$path
    run_manifest$snapshot_md5 <- snapshot$md5
    run_manifest$snapshot_config_path <- snapshot$snapshot_config_path
    run_manifest$snapshot_config_md5 <- snapshot$snapshot_config_md5
  }
  if (length(blocked)) {
    blocked_path <- .ecoda_output_path(
      file.path(output_dir, "blocked_datasets.json"), output_dir
    )
    .ecoda_write_json_atomic(
      list(run_id = run_id, analysis = analysis, scope = scope,
           union = selected, blocked = blocked), blocked_path
    )
    run_manifest$artifacts <- list(
      .ecoda_json_artifact(blocked_path, "ALL", "blocked_datasets")
    )
    ecoda_validate_derived_manifest(run_manifest, required_sources = TRUE)
    .ecoda_write_json_atomic(run_manifest, run_manifest_path)
    .ecoda_write_json_atomic(
      .ecoda_status(
        run_id, analysis, scope, "BLOCKED", output_dir,
        error = "one or more config/source rows are blocked",
        blocked = blocked, derived_root = derived_root
      ),
      status_path
    )
    stop("derived analysis blocked; see ", blocked_path)
  }

  artifacts <- list()
  if (analysis %in% c("res50", "harmony")) {
    method <- names(ECODA_DERIVED_METHODS[[analysis]])[[1L]]
    method_spec <- ECODA_DERIVED_METHODS[[analysis]][[method]]
    if (snapshot_mode) {
      snapshot_counts <- snapshot$counts[
        snapshot$counts$analysis == analysis, , drop = FALSE
      ]
      method <- names(ECODA_DERIVED_METHODS[[analysis]])[[1L]]
      for (dataset in selected) {
        ds_labels <- snapshot$labels[
          snapshot$labels$dataset == dataset, , drop = FALSE
        ]
        expected_samples <- as.character(ds_labels$Sample)
        labels <- .ecoda_snapshot_labels(snapshot$labels, dataset)
        ds_counts <- snapshot_counts[
          snapshot_counts$dataset == dataset, , drop = FALSE
        ]
        bundle <- .ecoda_snapshot_bundle(
          ds_counts, labels, expected_samples, method
        )
        artifact_path <- .ecoda_output_path(
          file.path(output_dir, paste0(dataset, "_", method, ".rds")),
          output_dir
        )
        .ecoda_save_rds(bundle, artifact_path)
        artifacts[[length(artifacts) + 1L]] <-
          .ecoda_artifact_record(artifact_path, dataset, method)
        rm(ds_labels, expected_samples, labels, ds_counts, bundle)
        gc(verbose = FALSE)
      }
    } else {
    for (dataset in selected) {
      source_config <- source_configs[[dataset]]
      obs <- .ecoda_load_derived_obs(
        source_config,
        analysis,
        expected_md5 = source_config$h5ad_md5
      )
      labels <- .ecoda_labels(obs, source_config$label_col)
      expected_samples <- unique(as.character(obs$Sample))
      bundle <- process_coda_fig(
        seurat = NULL,
        labels = labels,
        ct_col = method_spec$obs_col,
        obs = obs
      )
      .ecoda_validate_bundle(bundle, expected_samples, method)
      artifact_path <- .ecoda_output_path(
        file.path(output_dir, paste0(dataset, "_", method, ".rds")), output_dir
      )
      .ecoda_save_rds(bundle, artifact_path)
      artifacts[[length(artifacts) + 1L]] <-
        .ecoda_artifact_record(artifact_path, dataset, method)
      rm(obs, source_config, labels, expected_samples, bundle)
      gc(verbose = FALSE)
    }
    }
  } else {
    if (snapshot_mode) {
      result_rows <- list()
      row_index <- 0L
      targets <- ecoda_derived_targets()
      method <- "ECODA_authors_HR_cell_subsetting"
      for (dataset in selected) {
        ds_labels <- snapshot$labels[
          snapshot$labels$dataset == dataset, , drop = FALSE
        ]
        sample_ids <- as.character(ds_labels$Sample)
        total_counts <- setNames(
          as.integer(ds_labels$total_cells), sample_ids
        )
        annotated_counts <- setNames(
          as.integer(ds_labels$annotated_cells), sample_ids
        )
        original_counts <- setNames(
          as.integer(ds_labels$original_cells), sample_ids
        )
        if (!identical(original_counts, annotated_counts)) {
          stop("composition snapshot original cells must equal annotated cells")
        }
        ds_counts <- snapshot$counts[
          snapshot$counts$dataset == dataset &
            snapshot$counts$analysis == "cell_subsetting",
          , drop = FALSE
        ]
        labels <- .ecoda_snapshot_labels(snapshot$labels, dataset)
        plan <- ecoda_subsetting_plan(
          sample_ids, original_counts, seeds = seeds
        )
        for (plan_row in seq_len(nrow(plan))) {
          target <- plan$target[[plan_row]]
          target_cells <- plan$target_cells[[plan_row]]
          seed <- plan$seed[[plan_row]]
          key_mask <- as.character(ds_counts$target) == target &
            as.integer(ds_counts$replicate) ==
            as.integer(plan$replicate[[plan_row]])
          if (is.na(seed)) {
            key_mask <- key_mask & is.na(ds_counts$seed)
          } else {
            key_mask <- key_mask & !is.na(ds_counts$seed) &
              as.integer(ds_counts$seed) == as.integer(seed)
          }
          key_rows <- ds_counts[key_mask, , drop = FALSE]
          if (!nrow(key_rows)) {
            stop(
              "composition snapshot is missing cell-subsetting rows for ",
              dataset, "/", target, "/", plan$replicate[[plan_row]]
            )
          }
          effective_counts <- vapply(sample_ids, function(sample_id) {
            values <- unique(as.numeric(
              key_rows$effective_cells[key_rows$Sample == sample_id]
            ))
            if (length(values) != 1L) {
              stop("cell-subsetting effective count is not unique")
            }
            as.integer(values)
          }, integer(1))
          names(effective_counts) <- sample_ids
          bundle <- .ecoda_snapshot_bundle(
            key_rows, labels, sample_ids, method
          )
          score <- bundle$scores$anosim_score
          if (!is.numeric(score) || length(score) != 1L || !is.finite(score)) {
            stop("cell-subsetting ANOSIM is nonfinite for ", dataset, "/", target)
          }
          row_index <- row_index + 1L
          result_rows[[row_index]] <- data.frame(
            dataset = dataset,
            target = target,
            target_cells = as.integer(target_cells),
            replicate = as.integer(plan$replicate[[plan_row]]),
            seed = as.integer(seed),
            sample_ids = paste(sample_ids, collapse = "|"),
            total_cells_per_sample = paste(
              paste0(sample_ids, ":", as.integer(total_counts)),
              collapse = "|"
            ),
            original_cells_per_sample = paste(
              paste0(sample_ids, ":", as.integer(original_counts)),
              collapse = "|"
            ),
            annotated_cells_per_sample = paste(
              paste0(sample_ids, ":", as.integer(annotated_counts)),
              collapse = "|"
            ),
            effective_cells_per_sample = paste(
              paste0(sample_ids, ":", as.integer(effective_counts)),
              collapse = "|"
            ),
            original_total_cells = as.integer(sum(total_counts)),
            annotated_cells = as.integer(sum(annotated_counts)),
            total_cells = as.integer(sum(effective_counts)),
            effective_cells = as.integer(sum(effective_counts)),
            ANOSIM = as.numeric(score),
            stringsAsFactors = FALSE
          )
          rm(key_rows, effective_counts, bundle, score)
        }
        rm(
          ds_labels, sample_ids, total_counts, annotated_counts,
          original_counts, ds_counts, labels, plan
        )
        gc(verbose = FALSE)
      }
      result_table <- do.call(rbind, result_rows)
      .ecoda_validate_cell_subsetting_results(result_table)
      rds_path <- .ecoda_output_path(
        file.path(output_dir, "ECODA_authors_HR_cell_subsetting.rds"),
        output_dir
      )
      csv_path <- .ecoda_output_path(
        file.path(output_dir, "ECODA_authors_HR_cell_subsetting.csv"),
        output_dir
      )
      .ecoda_save_rds(result_table, rds_path)
      .ecoda_save_csv(result_table, csv_path)
      .ecoda_validate_output_sidecar(csv_path)
      artifacts[[1L]] <- .ecoda_artifact_record(
        rds_path, "ALL", "ECODA_authors_HR_cell_subsetting"
      )
      artifacts[[2L]] <- .ecoda_artifact_record(
        csv_path, "ALL", "ECODA_authors_HR_cell_subsetting"
      )
    } else {
    result_rows <- list()
    row_index <- 0L
    targets <- ecoda_derived_targets()
    for (dataset in selected) {
      source_config <- source_configs[[dataset]]
      obs <- .ecoda_load_derived_obs(
        source_config,
        analysis,
        expected_md5 = source_config$h5ad_md5
      )
      sample_values <- as.character(obs$Sample)
      sample_ids <- unique(sample_values)
      total_counts <- attr(obs, "ecoda_total_cells", exact = TRUE)
      annotated_counts <- attr(obs, "ecoda_annotated_cells", exact = TRUE)
      if (is.null(total_counts) || is.null(annotated_counts) ||
          is.null(names(total_counts)) || is.null(names(annotated_counts)) ||
          !identical(as.character(names(total_counts)), sample_ids) ||
          !identical(as.character(names(annotated_counts)), sample_ids)) {
        stop("cell-subsetting source is missing total/annotated cell diagnostics")
      }
      total_counts <- as.integer(total_counts)
      annotated_counts <- as.integer(annotated_counts)
      original_counts <- annotated_counts
      if (any(total_counts <= 0L) || any(annotated_counts <= 0L) ||
          any(annotated_counts > total_counts)) {
        stop("cell-subsetting source has invalid total/annotated cell diagnostics")
      }
      plan <- ecoda_subsetting_plan(sample_ids, original_counts, seeds = seeds)
      high_col <- source_config$high_col
      label_col <- source_config$label_col
      groups <- lapply(sample_ids, function(sample_id) which(sample_values == sample_id))
      names(groups) <- sample_ids
      for (plan_row in seq_len(nrow(plan))) {
        target <- plan$target[[plan_row]]
        target_cells <- plan$target_cells[[plan_row]]
        seed <- plan$seed[[plan_row]]
        if (identical(target, "all cells")) {
          keep <- seq_len(nrow(obs))
        } else {
          set.seed(seed)
          keep <- unlist(lapply(groups, function(indices) {
            n_keep <- min(length(indices), target_cells)
            if (n_keep == length(indices)) indices else sample(indices, n_keep, replace = FALSE)
          }), use.names = FALSE)
          keep <- sort(as.integer(keep))
        }
        obs_sub <- obs[keep, , drop = FALSE]
        effective_counts <- vapply(groups, function(indices) sum(keep %in% indices), integer(1))
        names(effective_counts) <- sample_ids
        labels <- .ecoda_labels(obs_sub, label_col)
        bundle <- process_coda_fig(
          seurat = NULL,
          labels = labels,
          ct_col = high_col,
          obs = obs_sub
        )
        .ecoda_validate_bundle(bundle, sample_ids, "ECODA_authors_HR_cell_subsetting")
        score <- bundle$scores$anosim_score
        if (!is.numeric(score) || length(score) != 1L || !is.finite(score)) {
          stop("cell-subsetting ANOSIM is nonfinite for ", dataset, "/", target)
        }
        row_index <- row_index + 1L
        result_rows[[row_index]] <- data.frame(
          dataset = dataset,
          target = target,
          target_cells = as.integer(target_cells),
          replicate = as.integer(plan$replicate[[plan_row]]),
          seed = as.integer(seed),
          sample_ids = paste(sample_ids, collapse = "|"),
          total_cells_per_sample = paste(
            paste0(sample_ids, ":", as.integer(total_counts)), collapse = "|"
          ),
          original_cells_per_sample = paste(
            paste0(sample_ids, ":", as.integer(original_counts)), collapse = "|"
          ),
          annotated_cells_per_sample = paste(
            paste0(sample_ids, ":", as.integer(annotated_counts)), collapse = "|"
          ),
          effective_cells_per_sample = paste(
            paste0(sample_ids, ":", as.integer(effective_counts)), collapse = "|"
          ),
          original_total_cells = as.integer(sum(total_counts)),
          annotated_cells = as.integer(sum(annotated_counts)),
          total_cells = as.integer(sum(effective_counts)),
          effective_cells = as.integer(sum(effective_counts)),
          ANOSIM = as.numeric(score),
          stringsAsFactors = FALSE
        )
        rm(obs_sub, labels, bundle, keep, score)
      }
      rm(
        obs, source_config, sample_values, sample_ids, total_counts,
        annotated_counts, original_counts, plan, groups, high_col, label_col,
        effective_counts, target, target_cells, seed
      )
      gc(verbose = FALSE)
    }
    result_table <- do.call(rbind, result_rows)
    .ecoda_validate_cell_subsetting_results(result_table)
    rds_path <- .ecoda_output_path(
      file.path(output_dir, "ECODA_authors_HR_cell_subsetting.rds"), output_dir
    )
    csv_path <- .ecoda_output_path(
      file.path(output_dir, "ECODA_authors_HR_cell_subsetting.csv"), output_dir
    )
    .ecoda_save_rds(result_table, rds_path)
    .ecoda_save_csv(result_table, csv_path)
    .ecoda_validate_output_sidecar(csv_path)
    artifacts[[1L]] <- .ecoda_artifact_record(
      rds_path, "ALL", "ECODA_authors_HR_cell_subsetting"
    )
    artifacts[[2L]] <- .ecoda_artifact_record(
      csv_path, "ALL", "ECODA_authors_HR_cell_subsetting"
    )
    }
  }

  run_manifest$status <- "COMPLETED"
  run_manifest$artifacts <- artifacts
  ecoda_validate_derived_manifest(run_manifest, required_sources = TRUE)
  .ecoda_write_json_atomic(run_manifest, run_manifest_path)
  .ecoda_write_json_atomic(
    .ecoda_status(
      run_id, analysis, scope, "COMPLETED", output_dir,
      derived_root = derived_root
    ),
    status_path
  )
  invisible(run_manifest)
}

# Do not execute when source()d for focused contract tests.  A test launched
# by Rscript also carries --file=..., so match the runner basename exactly.
runner_invocation <- any(vapply(
  commandArgs(trailingOnly = FALSE)[startsWith(
    commandArgs(trailingOnly = FALSE), "--file="
  )],
  function(token) {
    identical(
      basename(sub("^--file=", "", token)),
      "run_local_ecoda_derived.R"
    )
  },
  logical(1)
))
if (runner_invocation) {
  .ecoda_clear_in_process_owner()
  args <- .ecoda_parse_flags(commandArgs(trailingOnly = TRUE))
  tryCatch(
    .ecoda_run(args),
    error = function(error) {
      .ecoda_mark_failure(args, error)
      stop(conditionMessage(error), call. = FALSE)
    }
  )
}
