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
  if (!file.exists(file) || file.info(file)$size <= 0L || !file.exists(sidecar)) return(FALSE)
  lines <- readLines(sidecar, warn = FALSE)
  get <- function(prefix) {
    value <- lines[startsWith(lines, prefix)]
    if (!length(value)) return("")
    sub(paste0("^", prefix), "", value[[1L]])
  }
  identical(get("PATH="), file) &&
    identical(get("SIZE="), as.character(file.info(file)$size)) &&
    identical(get("MD5="), unname(tools::md5sum(file)))
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
batch_required_keys <- function(ds, label) {
  if (label == "gloscope") return("GloScope_hvg2000_pcadims30")
  if (label == "pseudobulk") return("Pseudobulk_hvg2000")
  if (label == "composition") {
    excluded <- unlist(config[[ds]]$not_suitable_for_auto_annotation %||% character())
    keys <- c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_seuratres_2")
    if (!all(c("hitme", "scatomic") %in% excluded)) {
      keys <- c(keys, "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR")
    }
    return(keys)
  }
  NULL
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

validate_result_file <- function(file, expected = NULL, required_keys = NULL, method = "") {
  if (!checksum_ok(file)) stop("Missing or invalid result checksum: ", file)
  bundle <- readRDS(file)
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
    if (!is.null(required_keys) &&
        (!identical(sort(names(bundle)), sort(required_keys)))) {
      stop("result combo keys do not match the method contract: ", file)
    }
    combos <- bundle
  }
  for (combo in combos) {
    if (!is.list(combo)) stop("Result combo is not a list: ", file)
    validate_combo(combo, file, expected, method)
  }
}
validate_pseudobulk <- function(file, expected = NULL) {
  if (!checksum_ok(file)) stop("Missing or invalid pseudobulk checksum: ", file)
  value <- readRDS(file)
  timing <- NULL
  memory <- NULL
  if (is.list(value) && !is.data.frame(value) && !is.null(value$pb)) {
    timing <- value$time_secs
    memory <- value$mem_GB
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
  if (!is.null(timing) && (!is.numeric(timing) || length(timing) != 1L || !finite_numeric(timing))) {
    stop("Pseudobulk timing is invalid: ", file)
  }
  if (!is.null(memory) && (!is.numeric(memory) || length(memory) != 1L || !finite_numeric(memory))) {
    stop("Pseudobulk memory is invalid: ", file)
  }
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
validate_metadata <- function(file, expected = NULL) {
  if (!checksum_ok(file)) stop("Missing or invalid metadata checksum: ", file)
  value <- readRDS(file)
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

validate_artifact_contract <- function(file, method, ds = "", view = "", metadata = FALSE) {
  expected <- if (nzchar(ds) && nzchar(view)) {
    expected_samples(ds, view)
  } else {
    NULL
  }
  if (metadata) {
    validate_metadata(file, expected)
  } else if (method == "trans") {
    validate_trans(file)
  } else if (method == "zeroimp") {
    validate_zeroimp(file)
  } else if (method == "prepare_pseudobulk" ||
             grepl("pseudobulk_", basename(file), fixed = TRUE)) {
    validate_pseudobulk(file, expected)
  } else {
    required_keys <- if (batch && method %in% c("gloscope", "pseudobulk", "composition")) {
      batch_required_keys(ds, method)
    } else {
      NULL
    }
    validate_result_file(file, expected, required_keys, method)
  }
}

if (nzchar(artifact_path)) {
  validate_artifact_contract(
    artifact_path, method_arg, dataset_arg, view_arg, metadata_kind
  )
  cat("benchmark RDS artifact contract OK\n")
  quit(save = "no", status = 0)
}

if (nzchar(artifact_list)) {
  for (part in artifact_list_parts) {
    validate_artifact_contract(
      part[[1L]], part[[2L]], part[[3L]], part[[4L]], part[[5L]] == "1"
    )
  }
  cat("benchmark RDS artifact-list contract OK\n")
  quit(save = "no", status = 0)
}

if (!dir.exists(root)) stop("benchmark result root is missing: ", root)
if (batch && !batch_pass %in% c("uncorrected", "corrected")) {
  stop("batch validation requires uncorrected or corrected pass")
}
if (batch && exact) {
  expected_rows <- paste(
    c("Alzheimer", "Breast_cancer", "Covid19_PBMC", "Kidney_KPMP",
      "Myocardial_infarction", "Diabetes", "Lupus_PBMC", "Lung",
      "Parkinson", "Joanito", "Stephenson", "CombinedPBMC"),
    "batch_effect_uncorrected", "batch_effect_uncorrected", sep = "\t"
  )
  if (!identical(selection_rows, expected_rows) || batch_pass != "uncorrected") {
    stop("batch exact selection is not the literal twelve-row uncorrected matrix")
  }
}

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
    if (label == "prepare_pseudobulk") {
      variants <- if (batch) "hvg2000" else c("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")
      stem <- if (batch) paste0(ds, "_batch_effect_", batch_pass) else ds
      for (variant in variants) {
        validate_pseudobulk(file.path(root, "pseudobulks", paste0(stem, "_pseudobulk_", variant, ".rds")), expected)
      }
    } else if (label == "trans") {
      validate_trans(file.path(root, "results", paste0(ds, "_trans.rds")))
    } else if (label == "zeroimp") {
      validate_zeroimp(file.path(root, "results", paste0(ds, "_zeroimp.rds")))
    } else {
      stem <- if (batch) paste0(ds, "_batch_effect_", batch_pass) else ds
      file <- file.path(root, "results", paste0(stem, "_", label, ".rds"))
      validate_result_file(file, expected, if (batch) batch_required_keys(ds, label) else NULL, label)
      if (label == "composition") {
        validate_metadata(file.path(root, "results", paste0(stem, "_metadata.rds")), expected)
      }
    }
  }
}

partials <- list.files(root, recursive = TRUE, full.names = TRUE)
partials <- partials[grepl("\\.(tmp|partial|build)(\\.|$)", basename(partials))]
if (length(partials)) stop("partial benchmark artifacts remain: ", paste(partials, collapse = ", "))
cat("benchmark RDS bundle contract OK\n")
