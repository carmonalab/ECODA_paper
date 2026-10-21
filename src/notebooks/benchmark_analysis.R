# Helpers for notebooks/benchmark_analysis.rmd.
# Source into the notebook evaluation environment after load_all_functions.R;
# scope helpers resolve the notebook's selected methods, datasets, and metrics
# when called, not when this file is sourced.

# Analysis scope
# Scope helpers (closed over analysis_methods, resolved after df_benchmark is
# built; only called from figure chunks, which run later).
methods_in_scope <- function(x) intersect(x, analysis_methods)
methods_in_scope_named <- function(m) m[names(m) %in% analysis_methods]

# Derived-run loading
# Strict loader for local, run-owned derived analyses.  A derived RDS is never
# trusted from its filename alone: the completed run record, source manifest,
# output record, output checksum, and (for ordinary H5AD mode) source H5AD
# digest are all checked before readRDS() is called.  Snapshot mode verifies
# the compact local composition snapshot and preserves remote H5AD paths only
# as provenance.
.derived_require_scalar <- function(x, label, allow_na = FALSE) {
  if (!is.character(x) || length(x) != 1L ||
      (!allow_na && (is.na(x) || !nzchar(x)))) {
    stop("Invalid ", label, " in derived-run manifest.")
  }
  x
}

.derived_resolve_path <- function(path, base_dir = NULL) {
  .derived_require_scalar(as.character(path), "path")
  if (!is.null(base_dir) &&
      !grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(base_dir, path)
  }
  normalizePath(path, mustWork = FALSE)
}

.derived_path_within <- function(path, root) {
  path <- .derived_resolve_path(path)
  root <- .derived_resolve_path(root)
  boundary <- if (identical(root, .Platform$file.sep)) {
    root
  } else {
    paste0(root, .Platform$file.sep)
  }
  identical(path, root) || startsWith(path, boundary)
}

.derived_run_artifact_path <- function(run_root, relative_path,
                                       label = "derived artifact") {
  root <- .derived_resolve_path(run_root)
  relative_path <- .derived_require_scalar(
    as.character(relative_path),
    paste0(label, " relative path")
  )
  if (!dir.exists(root) ||
      grepl("^(/|[A-Za-z]:[/\\\\])", relative_path, perl = TRUE) ||
      grepl("^[A-Za-z]:", relative_path, perl = TRUE) ||
      relative_path %in% c(".", "..") ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", relative_path, perl = TRUE) ||
      grepl("[\\\\]", relative_path, perl = TRUE) ||
      grepl("[[:cntrl:]]", relative_path, perl = TRUE)) {
    stop(label, " is not a safe path within the requested run root.")
  }
  path <- .derived_resolve_path(file.path(root, relative_path))
  if (!.derived_path_within(path, root) || identical(path, root)) {
    stop(label, " escapes the requested run root: ", path)
  }
  path
}

.derived_dataset_artifact_path <- function(run_root, dataset, suffix,
                                           label = "derived artifact") {
  dataset <- .derived_require_scalar(
    as.character(dataset),
    paste0(label, " dataset")
  )
  suffix <- .derived_require_scalar(
    as.character(suffix),
    paste0(label, " suffix")
  )
  if (dataset %in% c(".", "..") ||
      grepl("^[A-Za-z]:", dataset, perl = TRUE) ||
      grepl("[/\\\\]", dataset, perl = TRUE) ||
      grepl("[[:cntrl:]]", dataset, perl = TRUE) ||
      grepl("^(/|[A-Za-z]:[/\\\\])", suffix, perl = TRUE) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", suffix, perl = TRUE) ||
      grepl("[\\\\]", suffix, perl = TRUE) ||
      grepl("[[:cntrl:]]", suffix, perl = TRUE)) {
    stop(label, " has an unsafe dataset/path component.")
  }
  .derived_run_artifact_path(
    run_root,
    paste0(dataset, suffix),
    label
  )
}

.derived_md5 <- function(path, label = "artifact") {
  if (!file.exists(path) || !isTRUE(file.info(path)$isdir == FALSE) ||
      !isTRUE(file.info(path)$size > 0)) {
    stop("Missing or empty ", label, ": ", path)
  }
  digest <- tolower(unname(tools::md5sum(path)))
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[[:xdigit:]]{32}$", digest)) {
    stop("Could not compute a valid MD5 for ", label, ": ", path)
  }
  digest
}

.derived_read_json <- function(path, label) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    stop("Missing or empty ", label, ": ", path)
  }
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) {
      stop("Invalid ", label, ": ", path, ": ", conditionMessage(e))
    }
  )
}

.derived_require_fields <- function(x, fields, label) {
  if (!is.list(x) || !all(fields %in% names(x))) {
    stop(label, " is missing required fields: ",
         paste(setdiff(fields, names(x)), collapse = ", "))
  }
  invisible(x)
}

.derived_sidecar <- function(path, require_path = TRUE) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || !isTRUE(file.info(sidecar)$size > 0)) {
    stop("Missing derived-artifact checksum sidecar: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L || any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("Malformed derived-artifact checksum sidecar: ", sidecar)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[[:xdigit:]]{32}$", md5) ||
      !grepl("^[0-9]+$", size) || !nzchar(recorded_path)) {
    stop("Malformed derived-artifact checksum sidecar: ", sidecar)
  }
  actual_size <- as.character(file.info(path)$size)
  actual_md5 <- .derived_md5(path)
  if (!identical(tolower(md5), actual_md5) ||
      !identical(size, actual_size)) {
    stop("Derived-artifact checksum sidecar mismatch: ", path)
  }
  if (isTRUE(require_path)) {
    path_matches <- identical(
      .derived_resolve_path(recorded_path),
      .derived_resolve_path(path)
    )
  } else {
    path_matches <- identical(basename(recorded_path), basename(path))
  }
  if (!path_matches) {
    stop("Derived-artifact checksum sidecar PATH mismatch: ", path)
  }
  list(md5 = actual_md5, size = as.numeric(actual_size))
}
.derived_snapshot_context <- function(
  snapshot_path_ref,
  snapshot_md5_ref,
  snapshot_size_ref = NULL,
  run_dir,
  expected_scope = NULL,
  expected_config_md5 = NULL
) {
  snapshot_path <- .derived_resolve_snapshot_file(
    snapshot_path_ref,
    run_dir,
    "composition snapshot"
  )
  snapshot_md5_ref <- as.character(snapshot_md5_ref)
  if (length(snapshot_md5_ref) != 1L ||
      !grepl("^[[:xdigit:]]{32}$", snapshot_md5_ref)) {
    stop("Invalid composition snapshot MD5 in derived manifest.")
  }
  snapshot_sidecar <- .derived_sidecar(snapshot_path, require_path = FALSE)
  if (!identical(snapshot_sidecar$md5, tolower(snapshot_md5_ref))) {
    stop("Composition snapshot checksum mismatch: ", snapshot_path)
  }
  if (!is.null(snapshot_size_ref)) {
    snapshot_size_ref <- suppressWarnings(as.numeric(snapshot_size_ref))
    if (length(snapshot_size_ref) != 1L ||
        !is.finite(snapshot_size_ref) ||
        snapshot_size_ref <= 0 ||
        snapshot_size_ref != snapshot_sidecar$size) {
      stop("Composition snapshot size mismatch: ", snapshot_path)
    }
  }

  snapshot_manifest_path <- .derived_resolve_path(
    file.path(dirname(snapshot_path), "composition_snapshot_manifest.json")
  )
  snapshot <- .derived_read_json(
    snapshot_manifest_path,
    "composition_snapshot_manifest.json"
  )
  .derived_require_fields(
    snapshot,
    c(
      "schema_version", "run_id", "scope", "status", "source_mode",
      "config_path", "config_md5", "extractor_path", "extractor_md5",
      "snapshot_path", "snapshot_md5", "snapshot_size", "sources"
    ),
    "composition_snapshot_manifest.json"
  )
  schema_version <- suppressWarnings(as.numeric(snapshot$schema_version))
  if (length(schema_version) != 1L || !is.finite(schema_version) ||
      schema_version != 1) {
    stop("Composition snapshot schema_version is invalid.")
  }
  .derived_require_scalar(
    as.character(snapshot$run_id),
    "composition snapshot run_id"
  )
  .derived_require_scalar(
    as.character(snapshot$scope),
    "composition snapshot scope"
  )
  if (!identical(as.character(snapshot$status), "COMPLETED")) {
    stop("Composition snapshot manifest is not a completed snapshot.")
  }
  if (!identical(as.character(snapshot$source_mode), "composition_snapshot")) {
    stop("Composition snapshot manifest source_mode is invalid.")
  }
  if (!is.null(expected_scope) &&
      !identical(as.character(snapshot$scope), as.character(expected_scope))) {
    stop("Composition snapshot scope does not match the derived run.")
  }
  for (field in c("config_path", "extractor_path")) {
    .derived_require_scalar(
      as.character(snapshot[[field]]),
      paste0("composition snapshot ", field)
    )
  }
  if (!grepl("^/", as.character(snapshot$extractor_path))) {
    stop("Composition snapshot extractor_path must be absolute.")
  }

  for (field in c("config_md5", "extractor_md5", "snapshot_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(snapshot[[field]]))) {
      stop("Invalid composition snapshot ", field, ".")
    }
  }
  if (!is.null(expected_config_md5) &&
      !identical(
        tolower(as.character(snapshot$config_md5)),
        tolower(as.character(expected_config_md5))
      )) {
    stop("Composition snapshot config checksum does not match the derived run.")
  }
  manifest_snapshot_md5 <- tolower(as.character(snapshot$snapshot_md5))
  if (!identical(manifest_snapshot_md5, snapshot_sidecar$md5) ||
      !identical(manifest_snapshot_md5, tolower(snapshot_md5_ref))) {
    stop("Composition snapshot manifest checksum does not match the RDS.")
  }
  manifest_snapshot_size <- suppressWarnings(as.numeric(snapshot$snapshot_size))
  if (length(manifest_snapshot_size) != 1L ||
      !is.finite(manifest_snapshot_size) ||
      manifest_snapshot_size <= 0 ||
      manifest_snapshot_size != snapshot_sidecar$size) {
    stop("Composition snapshot manifest size does not match the RDS.")
  }
  manifest_snapshot_path <- .derived_require_scalar(
    as.character(snapshot$snapshot_path),
    "composition snapshot snapshot_path"
  )
  json_fields <- c(
    "snapshot_json_path", "snapshot_json_md5", "snapshot_json_size"
  )
  json_present <- json_fields %in% names(snapshot)
  snapshot_json_path <- NULL
  if (any(json_present) && !all(json_present)) {
    stop("Composition snapshot JSON identity is incomplete.")
  }
  if (all(json_present)) {
    snapshot_json_path <- .derived_resolve_snapshot_file(
      snapshot$snapshot_json_path,
      dirname(snapshot_path),
      "composition snapshot JSON"
    )
    json_sidecar <- .derived_sidecar(
      snapshot_json_path, require_path = FALSE
    )
    json_md5 <- tolower(as.character(snapshot$snapshot_json_md5))
    json_size <- suppressWarnings(as.numeric(snapshot$snapshot_json_size))
    if (!grepl("^[[:xdigit:]]{32}$", json_md5) ||
        !identical(json_md5, json_sidecar$md5) ||
        length(json_size) != 1L || !is.finite(json_size) ||
        json_size != json_sidecar$size) {
      stop("Composition snapshot JSON identity does not match its sidecar.")
    }
  }


  if (!is.list(snapshot$sources) || length(snapshot$sources) == 0L) {
    stop("Composition snapshot manifest has no source provenance rows.")
  }
  snapshot_sources <- lapply(seq_along(snapshot$sources), function(i) {
    row <- snapshot$sources[[i]]
    .derived_require_fields(
      row,
      c("dataset", "view", "h5ad_path", "h5ad_md5"),
      paste0("composition snapshot source entry ", i)
    )
    dataset <- .derived_require_scalar(
      as.character(row$dataset),
      paste0("composition snapshot source ", i, " dataset")
    )
    view <- .derived_require_scalar(
      as.character(row$view),
      paste0("composition snapshot source ", i, " view")
    )
    if (!identical(view, "benchmark_analysis")) {
      stop("Composition snapshot source entry ", i,
           " is not the benchmark_analysis view.")
    }
    h5ad_path <- as.character(row$h5ad_path)
    if (length(h5ad_path) != 1L || is.na(h5ad_path) ||
        !grepl("^/", h5ad_path) || !nzchar(h5ad_path)) {
      stop("Composition snapshot source entry ", i,
           " must retain an absolute HPC H5AD path.")
    }
    h5ad_path <- .derived_resolve_path(h5ad_path)
    h5ad_md5 <- tolower(as.character(row$h5ad_md5))
    if (length(h5ad_md5) != 1L ||
        !grepl("^[[:xdigit:]]{32}$", h5ad_md5)) {
      stop("Invalid H5AD MD5 in composition snapshot source entry ", i)
    }
    get_source_field <- function(candidates, label) {
      present <- candidates[candidates %in% names(row)]
      if (!length(present)) {
        stop("Composition snapshot source entry ", i,
             " is missing ", label, ".")
      }
      value <- row[[present[[1L]]]]
      if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
        stop("Invalid ", label, " in composition snapshot source entry ", i)
      }
      value
    }
    source_size <- suppressWarnings(as.numeric(
      get_source_field(c("h5ad_size", "size"), "H5AD size")
    ))
    if (length(source_size) != 1L || !is.finite(source_size) ||
        source_size <= 0 || source_size != floor(source_size)) {
      stop("Invalid H5AD size in composition snapshot source entry ", i)
    }
    source_mtime <- get_source_field(c("h5ad_mtime", "mtime"), "H5AD mtime")
    label_col <- as.character(get_source_field(
      c("label_col", "label_column", "label"),
      "biological-label column"
    ))
    high_col <- as.character(get_source_field(
      c(
        "high_res_col", "high_res_column", "cell_type_high_res",
        "cell_type_high_res_column"
      ),
      "high-resolution cell-type column"
    ))
    list(
      dataset = dataset,
      view = view,
      h5ad_path = h5ad_path,
      h5ad_md5 = h5ad_md5,
      size = source_size,
      mtime = as.character(source_mtime),
      label_column = label_col,
      high_res_column = high_col
    )
  })
  source_datasets <- vapply(snapshot_sources, `[[`, character(1), "dataset")
  if (anyDuplicated(source_datasets)) {
    stop("Composition snapshot source provenance contains duplicate datasets.")
  }
  list(
    path = snapshot_path,
    md5 = snapshot_sidecar$md5,
    size = snapshot_sidecar$size,
    manifest_path = snapshot_manifest_path,
    snapshot_config_path = as.character(snapshot$config_path),
    snapshot_config_md5 = tolower(as.character(snapshot$config_md5)),
    snapshot_remote_path = manifest_snapshot_path,
    snapshot_json_path = if (exists("snapshot_json_path")) {
      snapshot_json_path
    } else {
      NULL
    },
    manifest = snapshot,
    sources = snapshot_sources
  )
}
.derived_snapshot_source_field <- function(row, candidates, label) {
  present <- candidates[candidates %in% names(row)]
  if (!length(present)) {
    stop("Composition snapshot source is missing ", label, ".")
  }
  values <- vapply(present, function(name) {
    value <- row[[name]]
    if (length(value) != 1L || is.na(value) ||
        !nzchar(as.character(value))) {
      stop("Invalid ", label, " in composition snapshot source.")
    }
    as.character(value)
  }, character(1))
  if (length(unique(values)) != 1L) {
    stop("Conflicting ", label, " fields in composition snapshot source.")
  }
  values[[1L]]
}

.derived_compare_snapshot_sources <- function(snapshot_sources, source_rows) {
  if (length(snapshot_sources) != length(source_rows)) {
    stop("Derived source rows do not match the composition snapshot.")
  }
  for (i in seq_along(snapshot_sources)) {
    for (field in c(
      "dataset", "view", "h5ad_path", "h5ad_md5", "size", "mtime",
      "label_column", "high_res_column"
    )) {
      if (!field %in% names(source_rows[[i]]) ||
          !identical(
            as.character(snapshot_sources[[i]][[field]]),
            as.character(source_rows[[i]][[field]])
          )) {
        stop("Derived source row ", i,
             " does not match the composition snapshot provenance.")
      }
    }
  }
  invisible(TRUE)
}


.derived_resolve_bound_file <- function(path, run_dir, label) {
  path <- .derived_require_scalar(as.character(path), label)
  candidates <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path
  } else {
    c(file.path(run_dir, path), file.path(getwd(), path))
  }
  candidates <- unique(vapply(
    candidates,
    function(x) .derived_resolve_path(x),
    character(1)
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) != 1L) {
    stop("Could not resolve exactly one ", label, ": ", path)
  }
  existing[[1L]]
}

.derived_resolve_snapshot_file <- function(path_ref, run_dir, label) {
  path_ref <- .derived_require_scalar(as.character(path_ref), label)
  candidates <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path_ref)) {
    path_ref
  } else {
    c(file.path(run_dir, path_ref), file.path(getwd(), path_ref))
  }
  candidates <- c(candidates, file.path(run_dir, basename(path_ref)))
  candidates <- unique(vapply(
    candidates,
    function(x) .derived_resolve_path(x),
    character(1)
  ))
  info <- file.info(candidates)
  existing <- candidates[
    file.exists(candidates) & !is.na(info$isdir) & !info$isdir
  ]
  if (length(existing) != 1L) {
    stop("Could not resolve exactly one ", label, ": ", path_ref)
  }
  existing[[1L]]
}

.derived_validate_record <- function(record, run_dir, label) {
  .derived_require_fields(
    record,
    c("path", "size", "md5", "dataset", "method"),
    label
  )
  artifact_path <- .derived_resolve_path(record$path, run_dir)
  if (!.derived_path_within(artifact_path, run_dir) ||
      identical(artifact_path, .derived_resolve_path(run_dir))) {
    stop(label, " escapes the run-owned directory: ", artifact_path)
  }
  if (!grepl("^[[:xdigit:]]{32}$", as.character(record$md5)) ||
      !grepl("^[1-9][0-9]*$", as.character(record$size))) {
    stop("Invalid checksum/size in ", label, ": ", artifact_path)
  }
  actual_size <- as.character(file.info(artifact_path)$size)
  actual_md5 <- .derived_md5(artifact_path, label)
  if (!identical(as.character(record$size), actual_size) ||
      !identical(tolower(as.character(record$md5)), actual_md5)) {
    stop("Derived-run artifact manifest mismatch: ", artifact_path)
  }
  if (!is.character(record$dataset) || length(record$dataset) != 1L ||
      is.na(record$dataset) ||
      !is.character(record$method) || length(record$method) != 1L ||
      is.na(record$method)) {
    stop("Invalid dataset/method provenance in ", label)
  }
  list(
    path = artifact_path,
    size = as.numeric(actual_size),
    md5 = actual_md5,
    dataset = record$dataset,
    method = record$method
  )
}

.derived_run_context <- function(path, expected_analysis = NULL) {
  path <- .derived_resolve_path(path)
  run_dir <- dirname(path)
  run_manifest_path <- file.path(run_dir, "derived_run_manifest.json")
  source_manifest_default <- file.path(run_dir, "derived_source_manifest.json")
  status_path <- file.path(run_dir, "status.json")
  run <- .derived_read_json(run_manifest_path, "derived_run_manifest.json")
  .derived_require_fields(
    run,
    c(
      "run_id", "analysis", "scope", "status", "output_dir", "derived_root",
      "config_path",
      "config_md5", "runner_path", "runner_md5", "source_manifest_path",
      "source_manifest_md5", "artifacts", "sources"
    ),
    "derived_run_manifest.json"
  )
  source_mode <- if (is.null(run$source_mode)) {
    "h5ad"
  } else {
    .derived_require_scalar(
      as.character(run$source_mode),
      "derived_run_manifest.json source_mode"
    )
  }
  if (!source_mode %in%
      c("h5ad", "composition_snapshot", "local_composition")) {
    stop("Unsupported derived source_mode: ", source_mode)
  }
  snapshot_context <- NULL

  run_id <- .derived_require_scalar(
    as.character(run$run_id),
    "derived_run_manifest.json run_id"
  )
  owner_path <- file.path(run_dir, ".ecoda_run_owner")
  if (!file.exists(owner_path) ||
      !isTRUE(file_test("-f", owner_path)) ||
      isTRUE(nzchar(Sys.readlink(owner_path), keepNA = TRUE))) {
    stop("Missing or invalid .ecoda_run_owner marker: ", owner_path)
  }
  owner_lines <- tryCatch(
    readLines(owner_path, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(owner_lines) || length(owner_lines) != 1L ||
      !identical(owner_lines[[1L]], paste0("RUN_ID=", run_id))) {
    stop(".ecoda_run_owner marker does not match derived run: ", owner_path)
  }

  run_derived_root <- .derived_require_scalar(
    run$derived_root,
    "derived_run_manifest.json derived_root"
  )
  if (!identical(run_derived_root, .derived_resolve_path(run_derived_root))) {
    stop("derived_run_manifest.json derived_root must be a canonical path")
  }
  if (!.derived_path_within(run_dir, run_derived_root) ||
      identical(run_dir, run_derived_root)) {
    stop(
      "Derived run directory must be a strict descendant of derived_root: ",
      run_derived_root
    )
  }
  if (!identical(as.character(run$status), "COMPLETED") ||
      (!is.null(expected_analysis) &&
       !identical(as.character(run$analysis), expected_analysis))) {
    stop(
      "Derived run is not a completed ",
      if (is.null(expected_analysis)) "" else expected_analysis,
      " run: ", run_manifest_path
    )
  }
  status <- .derived_read_json(status_path, "status.json")
  .derived_require_fields(
    status,
    c("status", "run_id", "analysis", "scope", "output_dir", "derived_root"),
    "status.json"
  )
  if (!identical(as.character(status$status), "COMPLETED")) {
    stop("Derived run terminal status is not COMPLETED: ", status_path)
  }
  status_derived_root <- .derived_require_scalar(
    status$derived_root,
    "status.json derived_root"
  )
  if (!identical(status_derived_root,
                 .derived_resolve_path(status_derived_root))) {
    stop("status.json derived_root must be a canonical path")
  }
  if (!identical(status_derived_root, run_derived_root)) {
    stop("status.json derived_root does not match derived_run_manifest.json")
  }
  for (field in c("run_id", "analysis", "scope", "output_dir")) {
    .derived_require_scalar(
      as.character(status[[field]]),
      paste0("status.json ", field)
    )
    .derived_require_scalar(
      as.character(run[[field]]),
      paste0("derived_run_manifest.json ", field)
    )
    if (!identical(as.character(status[[field]]),
                   as.character(run[[field]]))) {
      stop("status.json ", field, " does not match derived_run_manifest.json")
    }
  }
  if (!identical(.derived_resolve_path(run$output_dir), run_dir)) {
    stop(
      "derived_run_manifest.json output_dir does not match the run-owned directory"
    )
  }
  for (field in c("config_md5", "runner_md5", "source_manifest_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(run[[field]]))) {
      stop("Invalid ", field, " in derived_run_manifest.json")
    }
  }

  config_path <- .derived_resolve_bound_file(
    run$config_path, run_dir, "config_path"
  )
  runner_path <- .derived_resolve_bound_file(
    run$runner_path, run_dir, "runner_path"
  )
  if (!identical(.derived_md5(config_path, "config_path"),
                 tolower(as.character(run$config_md5))) ||
      !identical(.derived_md5(runner_path, "runner_path"),
                 tolower(as.character(run$runner_md5)))) {
    stop("Derived run config/runner checksum mismatch: ", run_manifest_path)
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_require_fields(
      run,
      c(
        "snapshot_path", "snapshot_md5",
        "snapshot_config_path", "snapshot_config_md5"
      ),
      "composition_snapshot derived_run_manifest.json"
    )
    run_snapshot_config_path <- .derived_require_scalar(
      as.character(run$snapshot_config_path),
      "derived run snapshot_config_path"
    )
    run_snapshot_config_md5 <- .derived_require_scalar(
      as.character(run$snapshot_config_md5),
      "derived run snapshot_config_md5"
    )
    if (!grepl("^[[:xdigit:]]{32}$", run_snapshot_config_md5)) {
      stop("Invalid snapshot_config_md5 in derived_run_manifest.json")
    }
    snapshot_context <- .derived_snapshot_context(
      snapshot_path_ref = run$snapshot_path,
      snapshot_md5_ref = run$snapshot_md5,
      snapshot_size_ref = if ("snapshot_size" %in% names(run)) {
        run$snapshot_size
      } else {
        NULL
      },
      run_dir = run_dir,
      expected_scope = run$scope,
      expected_config_md5 = run$config_md5
    )
    if (!identical(
      run_snapshot_config_path,
      snapshot_context$snapshot_config_path
    ) ||
        !identical(
          tolower(run_snapshot_config_md5),
          snapshot_context$snapshot_config_md5
        )) {
      stop("Derived run snapshot config provenance does not match the snapshot.")
    }
  } else if (any(c(
    "snapshot_path", "snapshot_md5", "snapshot_size",
    "snapshot_config_path", "snapshot_config_md5"
  ) %in% names(run))) {
    stop("Snapshot identity fields require source_mode=composition_snapshot.")
  }


  source_manifest_ref <- .derived_require_scalar(
    as.character(run$source_manifest_path),
    "derived_run_manifest.json source_manifest_path"
  )
  source_manifest_path <- .derived_resolve_bound_file(
    source_manifest_ref,
    run_dir,
    "source_manifest_path"
  )
  canonical_source_manifest <- .derived_resolve_path(source_manifest_default)

  if (!identical(source_manifest_path, canonical_source_manifest)) {
    stop(
      "derived_run_manifest.json source_manifest_path must resolve to ",
      "the canonical run-root derived_source_manifest.json"
    )
  }

  if (!.derived_path_within(source_manifest_path, run_dir)) {
    stop("source_manifest_path escapes the run-owned directory.")
  }
  if (!identical(.derived_md5(source_manifest_path, "source manifest"),
                 tolower(as.character(run$source_manifest_md5)))) {
    stop("Derived source manifest checksum mismatch: ", source_manifest_path)
  }
  source <- .derived_read_json(source_manifest_path, "derived_source_manifest.json")
  .derived_require_fields(
    source,
    c(
      "run_id", "analysis", "scope", "derived_root", "config_path",
      "config_md5", "runner_path", "runner_md5", "sources"
    ),
    "derived_source_manifest.json"
  )
  source_mode_source <- if (is.null(source$source_mode)) {
    "h5ad"
  } else {
    .derived_require_scalar(
      as.character(source$source_mode),
      "derived_source_manifest.json source_mode"
    )
  }
  if (!source_mode_source %in%
      c("h5ad", "composition_snapshot", "local_composition") ||
      !identical(source_mode_source, source_mode)) {
    stop("Derived source manifest source_mode does not match the run.")
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_require_fields(
      source,
      c(
        "snapshot_path", "snapshot_md5",
        "snapshot_config_path", "snapshot_config_md5"
      ),
      "composition_snapshot derived_source_manifest.json"
    )
    source_snapshot_path <- .derived_require_scalar(
      as.character(source$snapshot_path),
      "derived source snapshot_path"
    )
    run_snapshot_path <- .derived_require_scalar(
      as.character(run$snapshot_path),
      "derived run snapshot_path"
    )
    source_snapshot_md5 <- .derived_require_scalar(
      as.character(source$snapshot_md5),
      "derived source snapshot_md5"
    )
    source_snapshot_config_path <- .derived_require_scalar(
      as.character(source$snapshot_config_path),
      "derived source snapshot_config_path"
    )
    source_snapshot_config_md5 <- .derived_require_scalar(
      as.character(source$snapshot_config_md5),
      "derived source snapshot_config_md5"
    )
    if (!grepl("^[[:xdigit:]]{32}$", source_snapshot_md5) ||
        !grepl("^[[:xdigit:]]{32}$", source_snapshot_config_md5)) {
      stop("Invalid snapshot identity in derived_source_manifest.json.")
    }
    if (!identical(
      basename(source_snapshot_path),
      basename(run_snapshot_path)
    ) ||
        !identical(
          tolower(source_snapshot_md5),
          tolower(as.character(run$snapshot_md5))
        )) {
      stop("Derived source/run snapshot identity does not match.")
    }
    if (!identical(
      source_snapshot_config_path,
      as.character(run$snapshot_config_path)
    ) ||
        !identical(
          tolower(source_snapshot_config_md5),
          tolower(as.character(run$snapshot_config_md5))
        ) ||
        !identical(
          source_snapshot_config_path,
          snapshot_context$snapshot_config_path
        ) ||
        !identical(
          tolower(source_snapshot_config_md5),
          snapshot_context$snapshot_config_md5
        )) {
      stop("Derived source/run snapshot config provenance does not match.")
    }
    snapshot_size_matches <- function(value, label) {
      value <- suppressWarnings(as.numeric(value))
      if (length(value) != 1L || !is.finite(value) ||
          value != snapshot_context$size) {
        stop(label, " does not match the local snapshot.")
      }
      invisible(TRUE)
    }
    if ("snapshot_size" %in% names(run)) {
      snapshot_size_matches(run$snapshot_size, "Derived run snapshot size")
    }
    if ("snapshot_size" %in% names(source)) {
      snapshot_size_matches(
        source$snapshot_size,
        "Derived source snapshot size"
      )
    }
  } else if (any(c(
    "snapshot_path", "snapshot_md5", "snapshot_size",
    "snapshot_config_path", "snapshot_config_md5"
  ) %in% names(source))) {
    stop("Snapshot identity fields require source_mode=composition_snapshot.")
  }

  source_derived_root <- .derived_require_scalar(
    source$derived_root,
    "derived_source_manifest.json derived_root"
  )
  if (!identical(source_derived_root,
                 .derived_resolve_path(source_derived_root))) {
    stop("derived_source_manifest.json derived_root must be a canonical path")
  }
  for (field in c("run_id", "analysis", "scope", "config_path", "runner_path")) {
    .derived_require_scalar(
      as.character(source[[field]]),
      paste0("derived source manifest ", field)
    )
  }
  for (field in c("config_md5", "runner_md5")) {
    if (!grepl("^[[:xdigit:]]{32}$", as.character(source[[field]]))) {
      stop("Invalid ", field, " in derived_source_manifest.json")
    }
  }
  if (!identical(as.character(source$run_id), as.character(run$run_id)) ||
      !identical(as.character(source$analysis), as.character(run$analysis)) ||
      !identical(as.character(source$scope), as.character(run$scope)) ||
      !identical(source_derived_root, run_derived_root) ||
      !identical(as.character(source$config_path), as.character(run$config_path)) ||
      !identical(tolower(as.character(source$config_md5)),
                 tolower(as.character(run$config_md5))) ||
      !identical(as.character(source$runner_path), as.character(run$runner_path)) ||
      !identical(tolower(as.character(source$runner_md5)),
                 tolower(as.character(run$runner_md5)))) {
    stop("Derived source manifest run provenance does not match the run.")
  }
  if (!identical(as.character(source$analysis), as.character(run$analysis)) ||
      !is.list(source$sources) || length(source$sources) == 0L) {
    stop("Derived source manifest analysis/sources do not match the run.")
  }
  source_rows <- lapply(seq_along(source$sources), function(i) {
    row <- source$sources[[i]]
    dataset <- .derived_require_scalar(
      as.character(row$dataset),
      paste0("source entry ", i, " dataset")
    )
    if (identical(source_mode, "local_composition")) {
      local_fields <- c(
        "dataset", "view", "composition_path", "composition_md5",
        "composition_size", "composition_mtime", "metadata_path",
        "metadata_md5", "metadata_size", "metadata_mtime", "analysis"
      )
      .derived_require_fields(
        row, local_fields, paste0("local composition source entry ", i)
      )
      if (!identical(as.character(row$view), "local_composition") ||
          !identical(as.character(row$analysis), as.character(run$analysis)) ||
          !grepl("^[[:xdigit:]]{32}$", as.character(row$composition_md5)) ||
          !grepl("^[[:xdigit:]]{32}$", as.character(row$metadata_md5))) {
        stop("Invalid local composition source entry ", i)
      }
      validate_local_input <- function(
        path_ref,
        expected_md5,
        expected_size,
        expected_mtime,
        label
      ) {
        path <- .derived_resolve_path(path_ref)
        info <- file.info(path)
        if (!file.exists(path) || isTRUE(info$isdir) ||
            is.na(info$size) || info$size <= 0 || is.na(info$mtime)) {
          stop("Invalid ", label, ": ", path)
        }
        actual_md5 <- .derived_md5(path, label)
        if (!identical(actual_md5, tolower(as.character(expected_md5))) ||
            !identical(
              as.character(info$size),
              as.character(expected_size)
            ) ||
            !identical(
              as.character(as.numeric(info$mtime)),
              as.character(as.numeric(expected_mtime))
            )) {
          stop(label, " identity mismatch: ", path)
        }
        list(
          path = path,
          md5 = actual_md5,
          size = as.numeric(info$size),
          mtime = as.numeric(info$mtime)
        )
      }
      composition_identity <- validate_local_input(
        row$composition_path,
        row$composition_md5,
        row$composition_size,
        row$composition_mtime,
        "local HR composition bundle"
      )
      metadata_identity <- validate_local_input(
        row$metadata_path,
        row$metadata_md5,
        row$metadata_size,
        row$metadata_mtime,
        "local metadata bundle"
      )
      return(list(
        dataset = dataset,
        view = as.character(row$view),
        composition_path = composition_identity$path,
        composition_md5 = composition_identity$md5,
        composition_size = composition_identity$size,
        composition_mtime = composition_identity$mtime,
        metadata_path = metadata_identity$path,
        metadata_md5 = metadata_identity$md5,
        metadata_size = metadata_identity$size,
        metadata_mtime = metadata_identity$mtime,
        analysis = as.character(row$analysis)
      ))
    }
    .derived_require_fields(
      row, c("dataset", "view", "h5ad_path", "h5ad_md5", "analysis"),
      paste0("derived source entry ", i)
    )
    if (!identical(as.character(row$view), "benchmark_analysis") ||
        !identical(as.character(row$analysis), as.character(run$analysis)) ||
        !grepl("^[[:xdigit:]]{32}$", as.character(row$h5ad_md5))) {
      stop("Invalid derived source entry ", i)
    }
    h5ad_path <- as.character(row$h5ad_path)
    if (length(h5ad_path) != 1L || is.na(h5ad_path) || !nzchar(h5ad_path)) {
      stop("Invalid H5AD path in derived source entry ", i)
    }
    h5ad_path <- .derived_resolve_path(h5ad_path)
    h5ad_md5 <- tolower(as.character(row$h5ad_md5))
    if (identical(source_mode, "composition_snapshot")) {
      if (!grepl("^/", h5ad_path)) {
        stop("Snapshot source H5AD path is not absolute: ", h5ad_path)
      }
    } else {
      .derived_sidecar(h5ad_path)
      actual_h5ad_md5 <- .derived_md5(h5ad_path, "source H5AD")
      if (!identical(actual_h5ad_md5, h5ad_md5)) {
        stop("Source H5AD checksum mismatch: ", h5ad_path)
      }
      h5ad_md5 <- actual_h5ad_md5
    }
    source_row <- list(
      dataset = dataset,
      view = as.character(row$view),
      h5ad_path = h5ad_path,
      h5ad_md5 = h5ad_md5,
      analysis = as.character(row$analysis)
    )
    if (identical(source_mode, "composition_snapshot")) {
      source_size <- suppressWarnings(as.numeric(
        .derived_snapshot_source_field(
          row, c("size", "h5ad_size"), "H5AD size"
        )
      ))
      if (length(source_size) != 1L || !is.finite(source_size) ||
          source_size <= 0 || source_size != floor(source_size)) {
        stop("Invalid H5AD size in derived source entry ", i)
      }
      source_row$size <- source_size
      source_row$mtime <- as.character(.derived_snapshot_source_field(
        row, c("mtime", "h5ad_mtime"), "H5AD mtime"
      ))
      source_row$label_column <- as.character(
        .derived_snapshot_source_field(
          row, c("label_column", "label_col"), "biological-label column"
        )
      )
      source_row$high_res_column <- as.character(
        .derived_snapshot_source_field(
          row,
          c("high_res_column", "high_res_col"),
          "high-resolution cell-type column"
        )
      )
    }
    source_row
  })
  source_datasets <- vapply(source_rows, `[[`, character(1), "dataset")
  if (anyDuplicated(source_datasets)) {
    stop("Derived source manifest contains duplicate datasets.")
  }
  if (!is.list(run$sources) || length(run$sources) != length(source_rows)) {
    stop("derived_run_manifest.json sources do not match source manifest.")
  }
  if (identical(source_mode, "composition_snapshot")) {
    .derived_compare_snapshot_sources(
      snapshot_context$sources,
      source_rows
    )
  }
  for (i in seq_along(source_rows)) {
    run_source <- run$sources[[i]]
    if (identical(source_mode, "local_composition")) {
      local_fields <- c(
        "dataset", "view", "composition_path", "composition_md5",
        "composition_size", "composition_mtime", "metadata_path",
        "metadata_md5", "metadata_size", "metadata_mtime", "analysis"
      )
      .derived_require_fields(
        run_source, local_fields,
        paste0("derived run local source entry ", i)
      )
      for (field in local_fields) {
        if (!identical(
          as.character(run_source[[field]]),
          as.character(source_rows[[i]][[field]])
        )) {
          stop("Derived run/source local provenance mismatch at entry ", i)
        }
      }
    } else {
      .derived_require_fields(
        run_source,
        c("dataset", "view", "h5ad_path", "h5ad_md5", "analysis"),
        paste0("derived run source entry ", i)
      )
      for (field in c(
        "dataset", "view", "h5ad_path", "h5ad_md5", "analysis"
      )) {
        if (!identical(
          as.character(run_source[[field]]),
          as.character(source_rows[[i]][[field]])
        )) {
          stop("derived run/source manifest source mismatch at entry ", i)
        }
      }
      if (identical(source_mode, "composition_snapshot")) {
        provenance_fields <- list(
          size = c("size", "h5ad_size"),
          mtime = c("mtime", "h5ad_mtime"),
          label_column = c("label_column", "label_col"),
          high_res_column = c("high_res_column", "high_res_col")
        )
        for (field in names(provenance_fields)) {
          run_value <- .derived_snapshot_source_field(
            run_source,
            provenance_fields[[field]],
            field
          )
          if (!identical(
            as.character(run_value),
            as.character(source_rows[[i]][[field]])
          )) {
            stop("Derived run/source snapshot provenance mismatch at entry ", i)
          }
        }
      }
    }
  }

  if (!is.list(run$artifacts) || length(run$artifacts) == 0L) {
    stop("derived_run_manifest.json has no output artifact records.")
  }
  artifact_rows <- lapply(seq_along(run$artifacts), function(i) {
    .derived_validate_record(
      run$artifacts[[i]], run_dir,
      paste0("derived output artifact ", i)
    )
  })
  artifact_paths <- vapply(artifact_rows, `[[`, character(1), "path")
  if (anyDuplicated(artifact_paths)) {
    stop("Derived run manifest contains duplicate output paths.")
  }
  artifact_datasets <- vapply(artifact_rows, `[[`, character(1), "dataset")
  artifact_methods <- vapply(artifact_rows, `[[`, character(1), "method")
  if (run$analysis %in% c("res50", "harmony")) {
    expected_method <- if (identical(run$analysis, "res50")) {
      "ECODA_seuratres_50"
    } else {
      "ECODA_seuratres_2_harmony"
    }
    expected_paths <- paste0(
      source_datasets, "_", expected_method, ".rds"
    )
    if (!identical(artifact_datasets, source_datasets) ||
        !identical(
          artifact_methods,
          rep(expected_method, length(source_datasets))
        ) ||
        !identical(basename(artifact_paths), expected_paths)) {
      stop(
        "Derived run manifest does not record the exact ",
        run$analysis, " output set."
      )
    }
  } else if (identical(run$analysis, "cell_subsetting")) {
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
      stop("Derived run manifest does not record the exact cell-subsetting output set.")
    }
  }
  list(
    run_dir = run_dir,
    run = run,
    status = status,
    source = source,
    source_mode = source_mode,
    snapshot_path = if (is.null(snapshot_context)) NULL else snapshot_context$path,
    snapshot_md5 = if (is.null(snapshot_context)) NULL else snapshot_context$md5,
    snapshot_size = if (is.null(snapshot_context)) NULL else snapshot_context$size,
    snapshot_config_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_config_path
    },
    snapshot_config_md5 = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_config_md5
    },
    snapshot_remote_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_remote_path
    },
    snapshot_json_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$snapshot_json_path
    },
    snapshot_manifest_path = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$manifest_path
    },
    snapshot_manifest = if (is.null(snapshot_context)) {
      NULL
    } else {
      snapshot_context$manifest
    },
    sources = source_rows,
    artifacts = artifact_rows,
    config_path = config_path,
    runner_path = runner_path,
    source_manifest_path = source_manifest_path,
    derived_root = run_derived_root
  )
}

load_derived_bundle <- function(
  path,
  dataset,
  method,
  expected_samples = NULL,
  allow_dropped = FALSE
) {
  if (isTRUE(allow_dropped) &&
      length(method) == 1L &&
      !is.na(method) &&
      method %in% c(
        "ECODA_seuratres_2_harmony",
        "ECODA_seuratres_50"
      )) {
    stop(
      "allow_dropped=TRUE is not permitted for complete-sample method: ",
      method
    )
  }
  path <- .derived_resolve_path(path)
  expected_analysis <- switch(
    method,
    ECODA_seuratres_2_harmony = "harmony",
    ECODA_seuratres_50 = "res50",
    NULL
  )
  context <- .derived_run_context(path, expected_analysis = expected_analysis)
  records <- Filter(function(x) identical(x$path, path), context$artifacts)
  if (length(records) != 1L) {
    stop("Derived bundle is not uniquely recorded in the run manifest: ", path)
  }
  record <- records[[1L]]
  if (!identical(record$dataset, dataset) || !identical(record$method, method)) {
    stop("Derived bundle provenance mismatch: ", path)
  }
  source_rows <- Filter(function(x) identical(x$dataset, dataset), context$sources)
  if (length(source_rows) != 1L) {
    stop("Derived bundle has no unique source H5AD entry: ", dataset)
  }
  .derived_sidecar(path)
  bundle <- tryCatch(
    readRDS(path),
    error = function(e) stop("Could not read derived bundle ", path, ": ",
                             conditionMessage(e))
  )
  required <- c("feat_mat", "labels", "dist_mat", "scores", "counts")
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    stop("Derived bundle missing required fields: ", path)
  }
  feat_mat <- bundle$feat_mat
  feat_values <- suppressWarnings(as.numeric(as.matrix(feat_mat)))
  if ((!is.matrix(feat_mat) && !is.data.frame(feat_mat)) ||
      nrow(feat_mat) == 0L || ncol(feat_mat) == 0L ||
      is.null(rownames(feat_mat)) || anyNA(rownames(feat_mat)) ||
      any(!nzchar(rownames(feat_mat))) || anyDuplicated(rownames(feat_mat)) ||
      any(!is.finite(feat_values))) {
    stop("Invalid feat_mat in derived bundle: ", path)
  }
  actual_samples <- rownames(feat_mat)
  labels <- bundle$labels
  if (length(labels) != length(actual_samples) ||
      is.null(names(labels)) || anyDuplicated(names(labels)) ||
      !identical(as.character(names(labels)), as.character(actual_samples))) {
    stop("Derived bundle labels do not match feat_mat sample IDs: ", path)
  }
  distances <- bundle$dist_mat
  distance_matrix <- as.matrix(distances)
  if (!is.matrix(distance_matrix) ||
      nrow(distance_matrix) != length(actual_samples) ||
      ncol(distance_matrix) != length(actual_samples) ||
      is.null(rownames(distance_matrix)) ||
      is.null(colnames(distance_matrix)) ||
      !identical(as.character(rownames(distance_matrix)),
                 as.character(actual_samples)) ||
      !identical(as.character(colnames(distance_matrix)),
                 as.character(actual_samples)) ||
      any(!is.finite(distance_matrix))) {
    stop("Invalid dist_mat in derived bundle: ", path)
  }
  if (!is.list(bundle$scores) && !is.numeric(bundle$scores)) {
    stop("Invalid scores in derived bundle: ", path)
  }
  score_values <- unlist(bundle$scores, use.names = FALSE)
  if (!is.numeric(score_values) ||
      length(score_values) == 0L || any(!is.finite(score_values))) {
    stop("Derived bundle scores are empty or nonfinite: ", path)
  }
  counts <- bundle$counts
  counts_values <- suppressWarnings(as.numeric(as.matrix(counts)))
  if ((!is.matrix(counts) && !is.data.frame(counts)) ||
      nrow(counts) == 0L || ncol(counts) == 0L ||
      is.null(rownames(counts)) || anyDuplicated(rownames(counts)) ||
      !identical(as.character(rownames(counts)),
                 as.character(actual_samples)) ||
      any(!is.finite(counts_values))) {
    stop("Invalid counts/sample IDs in derived bundle: ", path)
  }
  if (!is.null(expected_samples)) {
    expected_samples <- as.character(expected_samples)
    if (length(expected_samples) == 0L || anyNA(expected_samples) ||
        anyDuplicated(expected_samples)) {
      stop("Invalid expected sample IDs for derived bundle: ", path)
    }
    extras <- setdiff(actual_samples, expected_samples)
    missing <- setdiff(expected_samples, actual_samples)
    if (length(extras) > 0L || (!allow_dropped && length(missing) > 0L)) {
      stop(
        "Derived bundle sample policy failed for ", dataset, "/", method,
        "; missing=", paste(missing, collapse = ","),
        "; extras=", paste(extras, collapse = ",")
      )
    }
    if (!allow_dropped && !identical(actual_samples, expected_samples)) {
      stop("Derived bundle sample order differs from expected samples: ", path)
    }
    if (allow_dropped && length(missing) > 0L) {
      attr(bundle, "derived_dropped_samples") <- missing
    }
  }
  bundle
}

load_derived_cell_subsetting <- function(
  rds_path,
  csv_path,
  expected_datasets,
  expected_samples_by_dataset
) {
  rds_path <- .derived_resolve_path(rds_path)
  csv_path <- .derived_resolve_path(csv_path)
  context <- .derived_run_context(rds_path, expected_analysis = "cell_subsetting")
  expected_datasets <- as.character(expected_datasets)
  if (length(expected_datasets) == 0L || anyDuplicated(expected_datasets)) {
    stop("Invalid expected datasets for cell-subsetting results.")
  }
  source_datasets <- vapply(context$sources, `[[`, character(1), "dataset")
  if (length(source_datasets) == 0L || anyDuplicated(source_datasets) ||
      !setequal(source_datasets, expected_datasets)) {
    stop("Cell-subsetting source manifest does not cover the expected datasets.")
  }
  rds_records <- Filter(function(x) identical(x$path, rds_path), context$artifacts)
  csv_records <- Filter(function(x) identical(x$path, csv_path), context$artifacts)
  if (length(rds_records) != 1L || length(csv_records) != 1L) {
    stop("Cell-subsetting RDS/CSV are not uniquely recorded in the run manifest.")
  }
  if (!identical(rds_records[[1L]]$method,
                 "ECODA_authors_HR_cell_subsetting") ||
      !identical(csv_records[[1L]]$method,
                 "ECODA_authors_HR_cell_subsetting")) {
    stop("Cell-subsetting artifact method provenance is invalid.")
  }
  .derived_sidecar(rds_path)
  .derived_sidecar(csv_path)
  raw <- tryCatch(
    readRDS(rds_path),
    error = function(e) stop("Could not read cell-subsetting RDS: ",
                             conditionMessage(e))
  )
  results <- if (is.data.frame(raw)) {
    raw
  } else if (is.list(raw) && is.data.frame(raw$results)) {
    raw$results
  } else {
    stop("Cell-subsetting RDS must contain a results data.frame.")
  }
  rds_table <- results
  required <- c(
    "dataset", "target", "target_cells", "replicate", "seed", "sample_ids",
    "original_cells_per_sample", "effective_cells_per_sample",
    "total_cells", "ANOSIM"
  )
  if (!all(required %in% names(results))) {
    stop("Cell-subsetting results missing fields: ",
         paste(setdiff(required, names(results)), collapse = ", "))
  }
  if (identical(context$source_mode, "local_composition")) {
    methodology_fields <- c(
      "source_mode", "sampling_method", "source_limitation"
    )
    if (!all(methodology_fields %in% names(results)) ||
        any(!nzchar(as.character(results$source_mode))) ||
        any(as.character(results$source_mode) != "local_composition") ||
        any(!nzchar(as.character(results$sampling_method))) ||
        any(!nzchar(as.character(results$source_limitation)))) {
      stop(
        "Local composition subsetting methodology is missing or invalid."
      )
    }
  }
  if (!setequal(unique(as.character(results$dataset)), expected_datasets)) {
    stop("Cell-subsetting results do not cover the expected datasets.")
  }
  target_levels <- c("all cells", "2000", "1000", "500", "400",
                     "300", "200", "150", "100", "50")
  results$dataset <- as.character(results$dataset)
  results$target <- as.character(results$target)
  if (any(!results$target %in% target_levels)) {
    stop("Cell-subsetting results contain an unexpected target.")
  }
  target_cells_raw <- results$target_cells
  results$target_cells <- suppressWarnings(as.numeric(target_cells_raw))
  target_cells_expected <- ifelse(
    results$target == "all cells",
    NA_real_,
    suppressWarnings(as.numeric(results$target))
  )
  if (length(results$target_cells) != nrow(results) ||
      any(is.na(results$target_cells) & !is.na(target_cells_raw)) ||
      any(is.nan(results$target_cells)) ||
      any(is.na(results$target_cells) != is.na(target_cells_expected)) ||
      any(!is.na(results$target_cells) &
          (!is.finite(results$target_cells) |
           results$target_cells < 1 |
           results$target_cells != floor(results$target_cells))) ||
      any(!is.na(results$target_cells) &
          results$target_cells != target_cells_expected)) {
    stop("Cell-subsetting target_cells values are invalid.")
  }
  if (any(!is.finite(as.numeric(results$replicate))) ||
      any(as.integer(results$replicate) != as.numeric(results$replicate)) ||
      any(!is.finite(as.numeric(results$ANOSIM)))) {
    stop("Cell-subsetting replicate/ANOSIM values are invalid.")
  }
  results$replicate <- as.integer(results$replicate)
  seed_raw <- results$seed
  results$seed <- suppressWarnings(as.numeric(seed_raw))
  if (length(results$seed) != nrow(results) ||
      any(is.na(results$seed) & !is.na(seed_raw)) ||
      any(is.nan(results$seed)) ||
      any(!is.na(results$seed) & !is.finite(results$seed))) {
    stop("Cell-subsetting seed values are invalid.")
  }
  results$anosim <- as.numeric(results$ANOSIM)
  as_ids <- function(x) {
    if (is.list(x)) {
      if (length(x) == 1L && is.list(x[[1L]])) x <- x[[1L]]
      else x <- unlist(x, use.names = FALSE)
    }
    if (length(x) == 1L && is.character(x)) {
      x <- strsplit(x, "|", fixed = TRUE)[[1L]]
    }
    trimws(as.character(x))
  }
  as_named_counts <- function(x) {
    if (is.list(x)) {
      if (length(x) == 1L && is.list(x[[1L]])) x <- x[[1L]]
      else x <- unlist(x, use.names = FALSE)
    }
    if (length(x) == 1L && is.character(x)) {
      tokens <- strsplit(x, "|", fixed = TRUE)[[1L]]
      fields <- strsplit(tokens, ":", fixed = TRUE)
      if (any(lengths(fields) != 2L)) return(numeric())
      values <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1), 2L)))
      names(values) <- trimws(vapply(fields, `[[`, character(1), 1L))
      return(values)
    }
    as.numeric(x)
  }
  sample_lists <- lapply(results$sample_ids, as_ids)
  original_count_lists <- lapply(
    results$original_cells_per_sample,
    as_named_counts
  )
  count_lists <- lapply(
    results$effective_cells_per_sample,
    as_named_counts
  )
  for (i in seq_len(nrow(results))) {
    ds <- results$dataset[[i]]
    expected_ids <- as.character(expected_samples_by_dataset[[ds]])
    ids <- sample_lists[[i]]
    original_counts <- original_count_lists[[i]]
    counts <- count_lists[[i]]
    if (length(expected_ids) == 0L || anyDuplicated(expected_ids) ||
        !identical(ids, expected_ids) || anyDuplicated(ids) ||
        length(original_counts) != length(ids) ||
        length(counts) != length(ids) ||
        !identical(names(original_counts), ids) ||
        !identical(names(counts), ids) ||
        any(!is.finite(original_counts)) || any(!is.finite(counts)) ||
        any(original_counts < 1) || any(counts < 1) ||
        any(original_counts != as.integer(original_counts)) ||
        any(counts != as.integer(counts)) ||
        !isTRUE(as.numeric(results$total_cells[[i]]) == sum(counts))) {
      stop("Cell-subsetting sample/count diagnostics are invalid at row ", i)
    }
  }
  for (ds in expected_datasets) {
    ds_rows <- results[results$dataset == ds, , drop = FALSE]
    if (!identical(sort(unique(ds_rows$target)), sort(target_levels))) {
      stop("Cell-subsetting targets are incomplete for dataset ", ds)
    }
    baseline <- ds_rows[ds_rows$target == "all cells", , drop = FALSE]
    if (nrow(baseline) != 1L || baseline$replicate[[1L]] != 0L ||
        !is.na(baseline$seed[[1L]])) {
      stop("Cell-subsetting baseline policy failed for dataset ", ds)
    }
    baseline_i <- which(rownames(results) == rownames(baseline)[[1L]])
    baseline_original <- original_count_lists[[baseline_i]]
    baseline_counts <- count_lists[[baseline_i]]
    if (!identical(baseline_counts, baseline_original)) {
      stop("Cell-subsetting all-cells baseline changed effective counts.")
    }
    for (target in target_levels[-1L]) {
      rows <- ds_rows[ds_rows$target == target, , drop = FALSE]
      if (nrow(rows) != 20L ||
          !setequal(rows$replicate, seq_len(20L)) ||
          !setequal(rows$seed, seq(101, 120))) {
        stop("Cell-subsetting replicate/seed policy failed for ", ds, "/", target)
      }
      target_counts <- as.numeric(target)
      for (row_name in rownames(rows)) {
        row_i <- which(rownames(results) == row_name)
        if (!identical(
          original_count_lists[[row_i]],
          baseline_original
        ) || !identical(
          count_lists[[row_i]],
          pmin(baseline_counts, target_counts)
        )) {
          stop("Cell-subsetting effective counts violate min(original,target).")
        }
      }
    }
  }
  csv <- tryCatch(
    utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Could not read cell-subsetting CSV: ",
                             conditionMessage(e))
  )
  csv_required <- c(
    "dataset", "target", "target_cells", "replicate", "seed", "sample_ids",
    "original_cells_per_sample", "effective_cells_per_sample",
    "total_cells", "ANOSIM"
  )
  if (!identical(names(csv), names(rds_table)) ||
      !all(csv_required %in% names(csv)) ||
      nrow(csv) != nrow(rds_table)) {
    stop("Cell-subsetting CSV does not match the RDS result table.")
  }
  exact_columns <- setdiff(names(rds_table), c("seed", "ANOSIM"))
  for (column in exact_columns) {
    if (!identical(as.character(csv[[column]]),
                   as.character(rds_table[[column]]))) {
      stop("Cell-subsetting CSV column does not match the RDS result table: ",
           column)
    }
  }
  csv_seed_text <- as.character(csv$seed)
  csv_seed <- suppressWarnings(as.numeric(csv_seed_text))
  rds_seed <- results$seed
  allowed_blank_seed <- is.na(csv_seed_text) |
    !nzchar(trimws(csv_seed_text)) |
    trimws(csv_seed_text) == "NA"
  if (any(is.na(csv_seed) & !allowed_blank_seed) ||
      any(!is.na(csv_seed) & !is.finite(csv_seed)) ||
      any(is.na(csv_seed) != is.na(rds_seed)) ||
      any(!is.na(csv_seed) & csv_seed != rds_seed)) {
    stop("Cell-subsetting CSV seed diagnostics do not match the RDS result table.")
  }
  csv_anosim <- suppressWarnings(as.numeric(csv$ANOSIM))
  if (any(!is.finite(csv_anosim))) {
    stop("Cell-subsetting CSV ANOSIM diagnostics are invalid.")
  }
  anosim_tolerance <- 1e-8
  if (any(abs(csv_anosim - results$anosim) >
          anosim_tolerance * pmax(1, abs(results$anosim)))) {
    stop("Cell-subsetting CSV ANOSIM diagnostics do not match the RDS result table.")
  }
  list(results = results, raw = raw, context = context)
}

# Optional derived-run roots
.derived_requested_root <- function(env_name, label) {
  value <- Sys.getenv(env_name, unset = NA_character_)
  if (is.na(value)) return(NULL)
  if (!nzchar(value)) {
    stop(env_name, " must be a non-empty directory when set.")
  }
  path <- normalizePath(path.expand(value), mustWork = FALSE)
  if (!dir.exists(path)) {
    stop("Completed ", label, " derived-run directory not found: ", path)
  }
  path
}

# Annotation eligibility
# datasets.json `not_suitable_for_auto_annotation` flag (see AGENTS.md):
# datasets flagged for the annotation method families are dropped from the
# HiTME/scATOMIC annotation-method figures (Figure 3 A + Supp fig 18) with
# this documented message. The dataset stays in Figure 2a, where the author
# CT columns drive the methods. Absent/empty flag = suitable for all.
flag_status_for <- function(ds_key) {
  e <- datasets[[ds_key]]
  if (is.null(e) || is.null(e$not_suitable_for_auto_annotation)) {
    character(0)
  } else {
    e$not_suitable_for_auto_annotation
  }
}

flag_drop_annotation_methods <- function(df) {
  df %>%
    filter(!(dataset %in% hitme_flagged_ds &
               method %in% c("ECODA_HiTME_HR_layer2", "ECODA_HiTME_HR_layer3"))) %>%
    filter(!(dataset %in% scatomic_flagged_ds & method == "ECODA_scATOMIC_HR"))
}

# Annotation cell-type counts
annotation_cell_type_count <- function(bundle, method, dataset) {
  if (is.null(bundle) || !"feat_mat" %in% names(bundle)) {
    stop(
      "Missing feature matrix for annotation method '", method,
      "' in dataset '", dataset, "'."
    )
  }

  if (identical(method, "ECODA_deconv")) {
    # ECODA_deconv has inferred EPIC components but no counts field.
    cell_type_names <- colnames(bundle$feat_mat)
    count_source <- "feat_mat component columns"
  } else {
    if (!"counts" %in% names(bundle) || is.null(bundle$counts)) {
      stop(
        "Missing top-level counts field for annotation method '", method,
        "' in dataset '", dataset, "'."
      )
    }
    cell_type_names <- colnames(bundle$counts)
    count_source <- "counts columns"
  }

  if (is.null(cell_type_names) || length(cell_type_names) == 0L) {
    stop(
      "Empty ", count_source, " for annotation method '", method,
      "' in dataset '", dataset, "'."
    )
  }

  cell_type_names <- trimws(as.character(cell_type_names))
  invalid <- is.na(cell_type_names) |
    !nzchar(cell_type_names) |
    cell_type_names %in% c("NA", "nan", "None", "Unknown")
  valid_names <- cell_type_names[!invalid]
  unique_names <- unique(valid_names)
  if (length(valid_names) == 0L || length(unique_names) != length(valid_names)) {
    stop(
      "Invalid or duplicated cell-type identities in ", count_source,
      " for annotation method '", method, "' in dataset '", dataset, "'."
    )
  }
  length(unique_names)
}

# Execution-log keys
# Normalize keys so the merge with df_benchmark (paste0(dataset, method)
# join, helpers.R::merge_exec_times) hits:
# - GongSharma: internal key Gongsharma_cmv_young_males is displayed as
#   GongSharma (same rename as bmark_result_list above);
# - PILOT/QOT/PILOT-GM-VAE: the feather rows are <method>_hvg<n>_highres,
#   the bundles <method>_hvg<n>.
normalize_exec_keys <- function(df) {
  df %>%
    mutate(
      dataset = recode(dataset, "Gongsharma_cmv_young_males" = "GongSharma"),
      method = gsub("^(PILOT_hvg[0-9]+)_highres$", "\\1", method),
      method = gsub("^(QOT_hvg[0-9]+)_highres$", "\\1", method),
      method = gsub("^(PILOT-GM-VAE_hvg[0-9]+)_highres$", "\\1", method)
    )
}

# Feather rows win for ordinary dataset/method records.  Shared preparation
# rows are different when their explicit timing identities differ; an
# anti_join on only (dataset, method) would silently discard those contexts.
execution_union_key <- function(df) {
  key <- paste(
    as.character(df[["dataset"]]),
    as.character(df[["method"]]),
    sep = "\r"
  )
  shared <- .execution_is_shared_timing_method(
    as.character(df[["method"]])
  )
  if ("timing_id" %in% names(df)) {
    identified <- shared & .execution_nonblank(df[["timing_id"]])
    if (any(identified)) {
      key[identified] <- paste(
        key[identified],
        as.character(df[["timing_id"]][identified]),
        sep = "\r"
      )
    }
  }
  key
}

# Execution-time and memory plots
scientific_A_10n <- function(x) {
  s <- scales::scientific_format(digits = 2)(x) 
  s_converted <- gsub("e\\+", " %*% 10^", s) 
  return(parse(text = s_converted)) 
}
# Fit and draw the regression in log10(number of cells) / log10(response)
# coordinates. The raw axes are log-scaled below, so the displayed line is
# straight in the same log-log space used for the model.
fit_log_log_lines <- function(data, y_col) {
  groups <- split(data, data$method)
  dplyr::bind_rows(lapply(groups, function(df) {
    x <- as.numeric(df$n_cells)
    y <- as.numeric(df[[y_col]])
    keep <- is.finite(x) & x > 0 & is.finite(y) & y > 0
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 2L) return(NULL)
    fit <- stats::lm(log10(y) ~ log10(x))
    x_grid <- 10^seq(
      log10(min(x)),
      log10(max(x)),
      length.out = 100L
    )
    data.frame(
      method = as.character(df$method[[1L]]),
      n_cells = x_grid,
      fitted = 10^stats::predict(
        fit,
        newdata = data.frame(x = x_grid)
      ),
      stringsAsFactors = FALSE
    )
  }))
}
# Place labels in a deterministic right-hand lane. Repulsion in raw
# coordinates is unstable after log transformation; spacing labels in the
# already-log10 y coordinate guarantees no text overlap.
spread_log_label_y <- function(data, min_gap = 0.2, x_offset = 0.18) {
  data <- data[order(data$log_pred_y, decreasing = TRUE), , drop = FALSE]
  data$log_label_y <- data$log_pred_y
  if (nrow(data) > 1L) {
    for (i in 2:nrow(data)) {
      data$log_label_y[i] <- min(
        data$log_label_y[i],
        data$log_label_y[i - 1L] - min_gap
      )
    }
  }
  data$log_label_x <- max(data$log_max_x) + x_offset
  data
}
spread_log_label_y_centered <- function(
  data,
  min_gap = 0.15,
  x_offset = 0.18
) {
  data <- data[order(data$log_pred_y, decreasing = TRUE), , drop = FALSE]
  n_labels <- nrow(data)
  if (n_labels == 0L) return(data)
  center_y <- mean(range(data$log_pred_y))
  span_y <- min_gap * max(n_labels - 1L, 0L)
  data$log_label_y <- seq(
    center_y + span_y / 2,
    center_y - span_y / 2,
    length.out = n_labels
  )
  data$log_label_x <- max(data$log_max_x) + x_offset
  data
}

# Rank scaling
scale_ranks_global <- function(x, n_meth) {
  if(all(is.na(x))) return(x)
  # BEST rank (1) should map to 1.0
  # WORST rank (n_meth) should map to 0.0
  # Formula: (Worst - Current) / (Worst - Best)
  (n_meth - x) / (n_meth - 1)
}

# Global benchmark averages
# Define a function to get the global average over the canonical benchmark
# metrics only. Alternate modularity neighborhoods are diagnostics-only.
get_global_avg <- function(score_list) {
  # 1. Calculate the mean for each dataset
  ds_means <- sapply(score_list, function(ds) {
    relevant_vals <- unlist(ds[names(ds) %in% benchmark_metrics])
    # Return the average of ANOSIM, Modularity [kNN3], and ARI.
    mean(relevant_vals, na.rm = TRUE)
  })
  
  # 2. Return the average of all dataset means
  return(mean(ds_means, na.rm = TRUE))
}

# Metric axes
# --- 2. REUSABLE AXIS AND PLOTTING HELPERS ---
# Use each variable's observed range rather than one broad range for every
# metric. This keeps narrow metrics such as Silhouette from looking flat while
# retaining a 5% margin around every observed point.
metric_axis_range <- function(x, metric_name) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    stop("Supp fig 1: no finite values available for ", metric_name)
  }

  observed <- range(x)
  span <- diff(observed)
  padding <- if (span > 0) 0.05 * span else max(0.05, 0.05 * abs(observed[1]))
  observed + c(-padding, padding)
}

metric_axis_breaks <- function(x, limits) {
  observed <- x[is.finite(x)]
  if (length(observed) == 0) {
    stop("Supp fig 1: no finite values available for tick breaks")
  }

  breaks <- scales::breaks_extended(n = 4)(limits)

  # Keep a negative boundary tick visible when the observed data cross zero
  # but the extended-break algorithm starts at zero.
  if (min(observed) < 0 && min(breaks) >= 0) {
    negative_break <- floor(min(observed) / 0.1) * 0.1
    breaks <- sort(unique(c(negative_break, breaks)))
  }

  if (length(breaks) < 3L || length(breaks) > 5L) {
    stop("Supp fig 1: expected 3-5 meaningful ticks, got ", length(breaks))
  }
  breaks
}

build_metric_axis_spec <- function(data, variables) {
  metric_values <- setNames(
    lapply(variables, function(metric) data[[metric]]),
    variables
  )
  axis_ranges <- setNames(
    lapply(variables, function(metric) {
      metric_axis_range(metric_values[[metric]], metric)
    }),
    variables
  )
  axis_breaks <- Map(metric_axis_breaks, metric_values, axis_ranges)
  axis_ranges <- Map(
    function(lims, breaks) range(c(lims, breaks)),
    axis_ranges,
    axis_breaks
  )
  stopifnot(
    all(vapply(
      axis_breaks,
      function(breaks) length(breaks) >= 3L && length(breaks) <= 5L,
      logical(1)
    )),
    all(vapply(names(axis_ranges), function(metric) {
      values <- metric_values[[metric]]
      values <- values[is.finite(values)]
      all(
        values >= axis_ranges[[metric]][1] &
          values <= axis_ranges[[metric]][2]
      )
    }, logical(1)))
  )
  list(
    values = metric_values,
    ranges = axis_ranges,
    breaks = axis_breaks
  )
}

format_metric_ticks <- function(x) {
  x[abs(x) < 0.005] <- 0
  sprintf("%.2f", x)
}

# Pairwise metric plots
# Add variable-specific meaningful 3-5 tick scales to each GGally panel. The
# coordinate limits are applied with coord_cartesian() so smoothing/statistical
# layers still see all rows; limits include every finite observation.
add_pairwise_axes <- function(plot, mapping, axis_spec, include_y = TRUE) {
  x_name <- rlang::as_name(mapping$x)
  plot <- plot +
    scale_x_continuous(
      breaks = axis_spec$breaks[[x_name]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    )

  if (include_y) {
    y_name <- rlang::as_name(mapping$y)
    plot <- plot +
      scale_y_continuous(
        breaks = axis_spec$breaks[[y_name]],
        labels = format_metric_ticks,
        expand = c(0, 0)
      ) +
      coord_cartesian(
        xlim = axis_spec$ranges[[x_name]],
        ylim = axis_spec$ranges[[y_name]],
        expand = FALSE
      )
  } else {
    plot <- plot +
      coord_cartesian(
        xlim = axis_spec$ranges[[x_name]],
        expand = FALSE
      )
  }

  plot
}

make_pairwise_upper <- function(digits = 2, text_size = 6.5) {
  function(data, mapping, ...) {
    ggally_cor(
      data,
      mapping,
      ...,
      stars = FALSE,
      digits = digits,
      title_args = list(size = text_size),
      group_args = list(size = text_size)
    )
  }
}

make_pairwise_lower <- function(axis_spec) {
  function(data, mapping, ...) {
    plot <- ggplot(data, mapping) +
      geom_point(alpha = 0.4, color = "#000000") +
      geom_smooth(method = "lm", color = "black", se = TRUE)
    add_pairwise_axes(plot, mapping, axis_spec)
  }
}

make_pairwise_diag <- function(axis_spec) {
  function(data, mapping, ...) {
    add_pairwise_axes(
      ggally_densityDiag(data, mapping, ...),
      mapping,
      axis_spec,
      include_y = FALSE
    )
  }
}


# Metric-versus-mean plots
# --- 3. DYNAMICALLY GENERATE "METRIC VS MEAN" PLOTS ---
# purrr::map iterates over your metrics and creates a list of plots
plot_mean_corr <- function(data, metric_var, axis_spec) {
  ggplot(data, aes(x = Mean_Score, y = .data[[metric_var]])) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.8
    ) +
    geom_point(alpha = 0.4, color = "#000000", size = 1.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 1, se = FALSE) +
    stat_cor(
      aes(label = after_stat(rr.label)),
      label.x.npc = "right",
      label.y.npc = "bottom",
      hjust = 1,
      vjust = 0,
      size = 5
    ) +
    labs(x = NULL, y = metric_var) +
    scale_x_continuous(
      breaks = axis_spec$breaks[["Mean_Score"]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = axis_spec$breaks[[metric_var]],
      labels = format_metric_ticks,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = axis_spec$ranges[["Mean_Score"]],
      ylim = axis_spec$ranges[[metric_var]],
      expand = FALSE
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 13),
      aspect.ratio = 1
    )
}

