# Shared strict contract for the uncorrected batch-effect analysis pass.
# Biological labels are evaluation metadata only; these helpers consume persisted
# distance artifacts and never rerun preprocessing or embedding methods.

.batch_stop <- function(...) {
  stop(..., call. = FALSE)
}

.batch_or <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

.batch_nonblank <- function(x) {
  !is.na(x) && nzchar(trimws(as.character(x)))
}

.batch_values <- function(x) {
  value <- as.character(x)
  value[is.na(value) | trimws(value) %in% c("", "<NA>", "NA", "nan", "None")] <- NA_character_
  value
}

.batch_method_names <- function() {
  c(
    "ECODA_authors_HR",
    "ECODA_seuratres_2",
    "Pseudobulk_hvg2000",
    "GloScope_hvg2000_pcadims30",
    "MrVI_hvg2000",
    "PILOT_hvg2000",
    "QOT_hvg2000"
  )
}

.batch_join_warnings <- function(...) {
  values <- as.character(unlist(list(...), use.names = FALSE))
  values <- values[!is.na(values) & nzchar(values)]
  paste(unique(values), collapse = "; ")
}

batch_uncorrected_method_specs <- function() {
  specs <- list(
    ECODA_authors_HR = list(
      method = "ECODA_authors_HR",
      kind = "rds",
      relative_path = "results/{dataset}_batch_effect_uncorrected_composition.rds",
      bundle_key = "ECODA_authors_HR"
    ),
    ECODA_seuratres_2 = list(
      method = "ECODA_seuratres_2",
      kind = "rds",
      relative_path = "results/{dataset}_batch_effect_uncorrected_composition.rds",
      bundle_key = "ECODA_seuratres_2"
    ),
    Pseudobulk_hvg2000 = list(
      method = "Pseudobulk_hvg2000",
      kind = "rds",
      relative_path = "results/{dataset}_batch_effect_uncorrected_pseudobulk.rds",
      bundle_key = "Pseudobulk_hvg2000"
    ),
    GloScope_hvg2000_pcadims30 = list(
      method = "GloScope_hvg2000_pcadims30",
      kind = "rds",
      relative_path = "results/{dataset}_batch_effect_uncorrected_gloscope.rds",
      bundle_key = "GloScope_hvg2000_pcadims30"
    ),
    MrVI_hvg2000 = list(
      method = "MrVI_hvg2000",
      kind = "feather",
      relative_path = "embeddings/{dataset}_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"
    ),
    PILOT_hvg2000 = list(
      method = "PILOT_hvg2000",
      kind = "feather",
      relative_path = "embeddings/{dataset}_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
    ),
    QOT_hvg2000 = list(
      method = "QOT_hvg2000",
      kind = "feather",
      relative_path = "embeddings/{dataset}_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather"
    )
  )
  if (!identical(names(specs), .batch_method_names())) {
    .batch_stop("uncorrected method specification order changed")
  }
  specs
}

.batch_sidecar_fields <- function(path) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || file.info(sidecar)$size <= 0) {
    .batch_stop("missing checksum sidecar for ", path, ": ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  fields <- list()
  for (line in lines) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    key <- parts[[1L]]
    if (key %in% c("MD5", "SIZE", "PATH")) {
      if (!is.null(fields[[key]])) .batch_stop("duplicate ", key, " in checksum sidecar: ", sidecar)
      fields[[key]] <- paste(parts[-1L], collapse = "=")
    }
  }
  missing <- setdiff(c("MD5", "SIZE", "PATH"), names(fields))
  if (length(missing) > 0L) {
    .batch_stop("checksum sidecar missing fields ", paste(missing, collapse = ", "), ": ", sidecar)
  }
  fields
}

validate_batch_artifact <- function(path, kind = "artifact") {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    .batch_stop(kind, " path must be one non-empty string")
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) .batch_stop(kind, " is missing: ", path)
  info <- file.info(path)
  if (isTRUE(info$isdir) || is.na(info$size) || info$size <= 0) {
    .batch_stop(kind, " is a directory or empty: ", path)
  }
  fields <- .batch_sidecar_fields(path)
  actual_md5 <- unname(tools::md5sum(path))
  actual_size <- as.character(file.info(path)$size)
  if (!identical(as.character(fields$PATH), path)) {
    .batch_stop(kind, " PATH mismatch: recorded '", fields$PATH, "', expected '", path, "'")
  }
  if (!identical(tolower(as.character(fields$MD5)), tolower(actual_md5))) {
    .batch_stop(kind, " MD5 mismatch: ", path)
  }
  if (!identical(as.character(fields$SIZE), actual_size)) {
    .batch_stop(kind, " SIZE mismatch: ", path)
  }
  invisible(list(path = path, md5 = actual_md5, size = as.numeric(actual_size), kind = kind))
}

.batch_write_checksum <- function(path) {
  sidecar <- paste0(path, ".md5")
  tmp <- paste0(sidecar, ".tmp.", Sys.getpid(), ".", as.integer(stats::runif(1, 1, 1e9)))
  writeLines(c(
    paste0("MD5=", unname(tools::md5sum(path))),
    paste0("SIZE=", file.info(path)$size),
    paste0("PATH=", normalizePath(path, mustWork = FALSE))
  ), tmp)
  if (!file.rename(tmp, sidecar)) {
    if (file.exists(sidecar)) unlink(sidecar)
    if (!file.rename(tmp, sidecar)) {
      if (file.exists(tmp)) unlink(tmp)
      .batch_stop("could not atomically install checksum sidecar: ", sidecar)
    }
  }
  invisible(sidecar)
}

write_batch_metadata_sidecar <- function(
  metadata,
  path,
  expected_sample_ids = NULL
) {
  if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"Sample" %in% colnames(metadata)) .batch_stop("metadata sidecar requires Sample column")
  sample_ids <- as.character(metadata[["Sample"]])
  if (length(sample_ids) == 0L || anyNA(sample_ids) || any(!nzchar(trimws(sample_ids))) || anyDuplicated(sample_ids)) {
    .batch_stop("metadata sidecar Sample values must be unique and nonblank")
  }
  if (!is.null(expected_sample_ids) && !identical(sample_ids, as.character(expected_sample_ids))) {
    .batch_stop("metadata sidecar Sample order does not match expected sample IDs")
  }
  for (column in names(metadata)) {
    if (is.factor(metadata[[column]])) metadata[[column]] <- as.character(metadata[[column]])
    if (is.list(metadata[[column]])) .batch_stop("metadata sidecar column is list-valued: ", column)
  }
  path <- normalizePath(path, mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", as.integer(stats::runif(1, 1, 1e9)))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  arrow::write_feather(metadata, tmp)
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) .batch_stop("empty metadata sidecar temporary file: ", tmp)
  if (!file.rename(tmp, path)) {
    if (file.exists(path)) unlink(path)
    if (!file.rename(tmp, path)) .batch_stop("could not atomically install metadata sidecar: ", path)
  }
  .batch_write_checksum(path)
  validate_batch_artifact(path, "metadata sidecar")
  invisible(path)
}

read_batch_metadata_sidecar <- function(path, expected_sample_ids, required_label) {
  validate_batch_artifact(path, "metadata sidecar")
  metadata <- tryCatch(
    arrow::read_feather(path),
    error = function(error) .batch_stop("malformed metadata sidecar ", path, ": ", conditionMessage(error))
  )
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!"Sample" %in% colnames(metadata)) .batch_stop("metadata sidecar lacks Sample: ", path)
  if (length(required_label) != 1L || is.na(required_label) || !nzchar(required_label)) {
    .batch_stop("required primary label is missing")
  }
  if (!required_label %in% colnames(metadata)) {
    .batch_stop("metadata sidecar lacks primary biological label '", required_label, "': ", path)
  }
  sample_ids <- as.character(metadata[["Sample"]])
  expected_sample_ids <- as.character(expected_sample_ids)
  if (length(sample_ids) == 0L || anyNA(sample_ids) || any(!nzchar(trimws(sample_ids))) || anyDuplicated(sample_ids)) {
    .batch_stop("metadata sidecar has invalid Sample values: ", path)
  }
  if (!identical(sample_ids, expected_sample_ids)) {
    .batch_stop("metadata sidecar Sample membership/order mismatch: ", path)
  }
  labels <- .batch_values(metadata[[required_label]])
  if (anyNA(labels)) .batch_stop("metadata sidecar primary biological label is incomplete: ", path)
  rownames(metadata) <- sample_ids
  metadata
}

.batch_py_list <- function(value) {
  if (is.null(value)) return(character())
  as.character(unlist(value, use.names = FALSE))
}

batch_candidate_registry <- function(config, dataset_order) {
  if (!requireNamespace("reticulate", quietly = TRUE)) .batch_stop("reticulate is required for the onboarding registry")
  dataset_order <- as.character(dataset_order)
  if (length(dataset_order) != 12L || anyDuplicated(dataset_order)) {
    .batch_stop("batch candidate registry requires exactly twelve unique datasets")
  }
  spec_path <- normalizePath(file.path("notebooks", "dataset_onboarding", "dataset_specs.py"), mustWork = TRUE)
  spec_module <- reticulate::import_from_path(
    "dataset_specs",
    path = dirname(spec_path),
    convert = TRUE
  )
  authoritative_order <- as.character(spec_module$BATCH_EFFECT_DATASET_ORDER)
  if (!identical(authoritative_order, dataset_order)) {
    .batch_stop("dataset order does not match BATCH_EFFECT_DATASET_ORDER")
  }
  technical_specs <- spec_module$BATCH_EFFECT_SPECS
  if (!setequal(names(technical_specs), dataset_order)) {
    .batch_stop("BATCH_EFFECT_SPECS does not cover exactly the twelve datasets")
  }
  rows <- list()
  row_index <- 0L
  for (dataset in dataset_order) {
    entry <- config[[dataset]]
    if (is.null(entry)) .batch_stop("datasets.json lacks batch-effect dataset: ", dataset)
    primary <- as.character(.batch_or(entry$label_col, ""))
    if (length(primary) != 1L || is.na(primary) || !nzchar(primary)) {
      .batch_stop(dataset, ": missing primary biological label in datasets.json")
    }
    full_spec <- spec_module$DATASET_SPECS[[dataset]]
    secondary <- character()
    if (!is.null(full_spec)) {
      bio_col <- as.character(.batch_or(full_spec$bio_col, ""))
      stable <- .batch_py_list(full_spec$sample_stable_cols)
      if (nzchar(bio_col)) secondary <- stable[stable != bio_col]
    }
    secondary <- unique(secondary[nzchar(secondary) & secondary != primary])
    technical <- unique(.batch_py_list(technical_specs[[dataset]]))
    technical <- technical[nzchar(technical) & technical != primary]
    secondary_kept <- secondary[!secondary %in% technical]
    candidates <- c(primary, secondary_kept, technical)
    candidate_class <- c(
      "biological_primary",
      rep("biological_secondary", length(secondary_kept)),
      rep("technical", length(technical))
    )
    ordered <- data.frame(
      dataset = dataset,
      candidate = candidates,
      candidate_class = candidate_class,
      is_primary = c(TRUE, rep(FALSE, length(candidates) - 1L)),
      is_secondary_biology = candidates %in% secondary,
      is_technical = candidates %in% technical,
      stringsAsFactors = FALSE
    )
    ordered$formula_alias <- paste0("batch_term_", seq_len(nrow(ordered)))
    row_index <- row_index + 1L
    rows[[row_index]] <- ordered
  }
  registry <- do.call(rbind, rows)
  rownames(registry) <- NULL
  registry
}

.batch_registry_for_dataset <- function(candidate_registry, dataset = NULL) {
  if (is.null(candidate_registry) || !is.data.frame(candidate_registry)) {
    .batch_stop("candidate registry must be a data.frame")
  }
  required <- c("candidate", "candidate_class", "is_primary", "formula_alias")
  missing <- setdiff(required, colnames(candidate_registry))
  if (length(missing) > 0L) .batch_stop("candidate registry lacks: ", paste(missing, collapse = ", "))
  if (!is.null(dataset)) {
    if (!"dataset" %in% colnames(candidate_registry)) .batch_stop("candidate registry lacks dataset")
    candidate_registry <- candidate_registry[candidate_registry$dataset == dataset, , drop = FALSE]
  }
  if (nrow(candidate_registry) == 0L) .batch_stop("candidate registry is empty")
  candidate_registry
}

.batch_method_path <- function(input_root, dataset, method_spec) {
  relative <- gsub("\\{dataset\\}", dataset, method_spec$relative_path, fixed = FALSE)
  normalizePath(file.path(input_root, relative), mustWork = FALSE)
}

.batch_validate_matrix <- function(matrix, expected_sample_ids, kind) {
  if (!is.matrix(matrix)) matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  expected_sample_ids <- as.character(expected_sample_ids)
  if (nrow(matrix) != ncol(matrix) || nrow(matrix) != length(expected_sample_ids)) {
    .batch_stop(kind, " distance matrix is not square or has unexpected size")
  }
  row_ids <- rownames(matrix)
  col_ids <- colnames(matrix)
  if (is.null(row_ids) || is.null(col_ids) || anyNA(row_ids) || anyNA(col_ids) ||
      any(!nzchar(row_ids)) || any(!nzchar(col_ids)) || anyDuplicated(row_ids) || anyDuplicated(col_ids) ||
      !identical(as.character(row_ids), as.character(col_ids)) ||
      !identical(as.character(row_ids), expected_sample_ids)) {
    .batch_stop(kind, " distance matrix sample IDs/order mismatch")
  }
  if (any(!is.finite(matrix))) .batch_stop(kind, " distance matrix contains non-finite values")
  matrix
}

read_batch_method_dist <- function(input_root, dataset, method_spec, expected_sample_ids) {
  if (is.null(method_spec$kind) || is.null(method_spec$relative_path)) {
    .batch_stop("malformed method specification for ", method_spec$method)
  }
  path <- .batch_method_path(input_root, dataset, method_spec)
  validate_batch_artifact(path, paste0(method_spec$method, " artifact"))
  scores <- NULL
  if (identical(method_spec$kind, "rds")) {
    bundle <- tryCatch(readRDS(path), error = function(error) {
      .batch_stop("malformed RDS distance artifact ", path, ": ", conditionMessage(error))
    })
    result <- bundle[[method_spec$bundle_key]]
    if (is.null(result) || is.null(result$dist_mat)) {
      .batch_stop("RDS distance bundle lacks key ", method_spec$bundle_key, ": ", path)
    }
    matrix <- .batch_validate_matrix(as.matrix(result$dist_mat), expected_sample_ids, "RDS")
    scores <- result$scores
  } else if (identical(method_spec$kind, "feather")) {
    frame <- tryCatch(arrow::read_feather(path), error = function(error) {
      .batch_stop("malformed Feather distance artifact ", path, ": ", conditionMessage(error))
    })
    frame <- as.data.frame(frame, stringsAsFactors = FALSE)
    if (ncol(frame) < 2L) .batch_stop("Feather distance artifact lacks an index column: ", path)
    index <- as.character(frame[[ncol(frame)]])
    column_ids <- as.character(colnames(frame)[seq_len(ncol(frame) - 1L)])
    if (anyNA(index) || any(!nzchar(index)) || anyDuplicated(index) || !identical(column_ids, index)) {
      .batch_stop("Feather distance artifact columns/index are not in identical order: ", path)
    }
    matrix <- as.matrix(frame[, seq_len(ncol(frame) - 1L), drop = FALSE])
    storage.mode(matrix) <- "double"
    rownames(matrix) <- index
    colnames(matrix) <- index
    matrix <- .batch_validate_matrix(matrix, expected_sample_ids, "Feather")
  } else {
    .batch_stop("unsupported batch distance artifact kind: ", method_spec$kind)
  }
  list(
    method = as.character(method_spec$method),
    path = path,
    kind = as.character(method_spec$kind),
    matrix = matrix,
    dist_mat = stats::as.dist(matrix),
    sample_ids = rownames(matrix),
    scores = scores
  )
}

load_batch_uncorrected_dataset <- function(input_root, metadata_root, dataset, registry) {
  registry <- .batch_registry_for_dataset(registry, dataset)
  metadata_summary_path <- normalizePath(
    file.path(input_root, "results", paste0(dataset, "_batch_effect_uncorrected_metadata.rds")),
    mustWork = FALSE
  )
  validate_batch_artifact(metadata_summary_path, "metadata summary")
  summary <- tryCatch(readRDS(metadata_summary_path), error = function(error) {
    .batch_stop("malformed metadata summary ", metadata_summary_path, ": ", conditionMessage(error))
  })
  labels <- summary$labels
  if (is.null(labels) || is.null(names(labels)) || anyNA(names(labels)) || any(!nzchar(names(labels))) ||
      anyDuplicated(names(labels)) || anyNA(labels)) {
    .batch_stop(dataset, ": metadata summary labels are malformed")
  }
  sample_ids <- as.character(names(labels))
  if (!identical(as.integer(summary$n_samples), length(sample_ids))) {
    .batch_stop(dataset, ": metadata summary sample count mismatch")
  }
  metadata_path <- normalizePath(
    file.path(metadata_root, paste0(dataset, "_sample_metadata.feather")),
    mustWork = FALSE
  )
  metadata <- read_batch_metadata_sidecar(metadata_path, sample_ids, registry$candidate[registry$is_primary][[1L]])
  label_col <- registry$candidate[registry$is_primary][[1L]]
  if (!identical(.batch_values(metadata[[label_col]]), .batch_values(labels))) {
    .batch_stop(dataset, ": metadata sidecar labels do not match metadata summary")
  }
  specs <- batch_uncorrected_method_specs()
  methods <- lapply(specs, function(spec) read_batch_method_dist(input_root, dataset, spec, sample_ids))
  names(methods) <- names(specs)
  list(
    dataset = dataset,
    metadata = metadata,
    sample_ids = sample_ids,
    label_col = label_col,
    registry = registry,
    methods = methods,
    metadata_summary_path = metadata_summary_path,
    metadata_path = metadata_path,
    artifact_paths = vapply(methods, function(x) x$path, character(1))
  )
}

.batch_validity <- function(metadata, candidate, label_col = NULL) {
  if (!candidate %in% colnames(metadata)) {
    return(list(valid = FALSE, status = "MISSING", warning = paste0("missing candidate column ", candidate), values = NULL))
  }
  values <- .batch_values(metadata[[candidate]])
  valid <- !is.na(values)
  warnings <- character()
  if (!all(valid)) warnings <- c(warnings, "incomplete candidate")
  levels <- unique(values[valid])
  if (length(levels) < 2L) warnings <- c(warnings, "constant candidate")
  if (length(levels) == sum(valid) && sum(valid) == nrow(metadata)) {
    warnings <- c(warnings, "sample-unique candidate")
  }
  if (!is.null(label_col) && label_col %in% colnames(metadata)) {
    labels <- .batch_values(metadata[[label_col]])
    pair <- valid & !is.na(labels)
    if (sum(pair) > 0L) {
      tab <- table(labels[pair], values[pair])
      perfect <- nrow(tab) == ncol(tab) && all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
      if (perfect) warnings <- c(warnings, "perfectly confounded with biology")
    }
  }
  if (length(warnings) > 0L) {
    status <- if (any(grepl("perfectly confounded", warnings, fixed = TRUE))) "NON_ESTIMABLE" else "INVALID"
    return(list(valid = FALSE, status = status, warning = paste(unique(warnings), collapse = "; "), values = values))
  }
  list(valid = TRUE, status = "VALID", warning = "", values = values)
}

.batch_prepare_distance <- function(dist_mat, sample_ids = NULL) {
  matrix <- if (inherits(dist_mat, "dist")) as.matrix(dist_mat) else as.matrix(dist_mat)
  storage.mode(matrix) <- "double"
  ids <- rownames(matrix)
  if (is.null(ids)) ids <- colnames(matrix)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(matrix)))
  if (is.null(sample_ids)) sample_ids <- ids
  matrix <- .batch_validate_matrix(matrix, as.character(sample_ids), "analysis")
  stats::as.dist(matrix)
}

.batch_complete_metadata <- function(metadata, columns) {
  keep <- rep(TRUE, nrow(metadata))
  for (column in columns) {
    if (!column %in% colnames(metadata)) return(integer())
    value <- .batch_values(metadata[[column]])
    keep <- keep & !is.na(value)
  }
  which(keep)
}

.batch_empty_metric <- function(status, warning, n_samples = 0L, candidate = NA_character_) {
  data.frame(
    candidate = candidate,
    r2 = NA_real_,
    pseudo_f = NA_real_,
    p_value = NA_real_,
    status = status,
    warning = warning,
    n_samples = as.integer(n_samples),
    stringsAsFactors = FALSE
  )
}

.batch_adonis_term <- function(
  dist_mat,
  data,
  aliases,
  term_alias,
  permutations,
  parallel = getOption("mc.cores")
) {
  distance <- dist_mat
  formula <- stats::reformulate(aliases, response = "distance")
  environment(formula) <- environment()
  warnings <- character()
  result <- tryCatch(
    withCallingHandlers(
      vegan::adonis2(
        formula,
        data = data,
        permutations = permutations,
        by = "margin",
        parallel = parallel
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    return(list(result = NULL, warning = paste(unique(c(warnings, conditionMessage(result))), collapse = "; ")))
  }
  row <- which(rownames(result) == term_alias)
  if (length(row) != 1L) {
    return(list(result = NULL, warning = paste(unique(c(warnings, "term not returned by adonis2")), collapse = "; ")))
  }
  list(result = result[row, , drop = FALSE], warning = paste(unique(warnings), collapse = "; "))
}

compute_batch_anosim <- function(
  dist_mat,
  grouping,
  permutations = 999L,
  parallel = getOption("mc.cores")
) {
  grouping <- .batch_values(grouping)
  valid <- !is.na(grouping)
  if (sum(valid) < 3L || length(unique(grouping[valid])) < 2L) {
    return(list(statistic = NA_real_, p_value = NA_real_, status = "INVALID", warning = "ANOSIM requires at least three samples and two levels", n_samples = sum(valid)))
  }
  distance_matrix <- as.matrix(dist_mat)
  if (nrow(distance_matrix) != length(grouping) || ncol(distance_matrix) != length(grouping)) {
    .batch_stop("ANOSIM grouping length does not match distance matrix")
  }
  distance <- stats::as.dist(distance_matrix[valid, valid, drop = FALSE])
  result <- tryCatch(
    suppressWarnings(vegan::anosim(
      distance,
      grouping[valid],
      permutations = as.integer(permutations),
      parallel = parallel
    )),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    return(list(statistic = NA_real_, p_value = NA_real_, status = "NON_ESTIMABLE", warning = conditionMessage(result), n_samples = sum(valid)))
  }
  list(
    statistic = as.numeric(result$statistic),
    p_value = as.numeric(result$signif),
    status = "ESTIMABLE",
    warning = "",
    n_samples = sum(valid)
  )
}

compute_batch_permanova <- function(
  dist_mat,
  metadata,
  candidate,
  biology_col = NULL,
  permutations = 999L,
  parallel = getOption("mc.cores")
) {
  if (is.null(biology_col)) {
    validity <- .batch_validity(metadata, candidate)
    model_columns <- candidate
  } else {
    if (!biology_col %in% colnames(metadata)) {
      return(.batch_empty_metric("MISSING", paste0("missing biology column ", biology_col), candidate = candidate))
    }
    validity <- .batch_validity(metadata, candidate, biology_col)
    model_columns <- c(biology_col, candidate)
  }
  if (!isTRUE(validity$valid)) {
    return(.batch_empty_metric(validity$status, validity$warning, candidate = candidate))
  }
  complete <- .batch_complete_metadata(metadata, model_columns)
  if (length(complete) < 3L) {
    return(.batch_empty_metric("INVALID", "fewer than three complete-case samples", candidate = candidate))
  }
  data <- data.frame(row.names = rownames(metadata)[complete], stringsAsFactors = FALSE)
  aliases <- character()
  for (index in seq_along(model_columns)) {
    alias <- paste0("batch_model_", index)
    data[[alias]] <- factor(.batch_values(metadata[[model_columns[[index]]]])[complete])
    aliases <- c(aliases, alias)
  }
  if (any(vapply(data, function(x) nlevels(x) < 2L, logical(1)))) {
    return(.batch_empty_metric("INVALID", "complete-case model contains a constant term", length(complete), candidate))
  }
  distance_matrix <- as.matrix(if (inherits(dist_mat, "dist")) dist_mat else dist_mat)
  distance_matrix <- distance_matrix[complete, complete, drop = FALSE]
  distance <- stats::as.dist(distance_matrix)
  term_alias <- aliases[[length(aliases)]]
  adonis <- .batch_adonis_term(distance, data, aliases, term_alias, permutations, parallel)
  if (is.null(adonis$result)) {
    return(.batch_empty_metric("NON_ESTIMABLE", adonis$warning, length(complete), candidate))
  }
  row <- adonis$result
  data.frame(
    candidate = candidate,
    r2 = as.numeric(row$R2[[1L]]),
    pseudo_f = as.numeric(row$F[[1L]]),
    p_value = as.numeric(row$`Pr(>F)`[[1L]]),
    status = "ESTIMABLE",
    warning = adonis$warning,
    n_samples = as.integer(length(complete)),
    stringsAsFactors = FALSE
  )
}

batch_joint_nested_rules <- function(dataset) {
  rules <- list()
  if (identical(dataset, "Lung")) {
    rules[["origin_fine"]] <- list(
      keep = "origin",
      reason = "origin_fine is a documented finer-grained origin field; retain origin"
    )
  }
  if (identical(dataset, "Kidney_KPMP")) {
    rules[["condition.l2"]] <- list(
      keep = "condition.l1",
      reason = "condition.l2 is finer-grained than condition.l1; retain condition.l1"
    )
    rules[["condition.long"]] <- list(
      keep = "condition.l1",
      reason = "condition.long is finer-grained than condition.l1; retain condition.l1"
    )
  }
  rules
}

.batch_partition_equivalent <- function(left, right) {
  left <- .batch_values(left)
  right <- .batch_values(right)
  if (length(left) != length(right)) return(FALSE)
  keep <- !is.na(left) & !is.na(right)
  if (!any(keep) || any(is.na(left)) || any(is.na(right))) return(FALSE)
  identical(
    outer(left[keep], left[keep], "=="),
    outer(right[keep], right[keep], "==")
  )
}

.batch_joint_model_class <- function(registry, source_indices) {
  classes <- unique(registry$candidate_class[source_indices])
  if (length(classes) == 1L) classes else "shared_biological_technical"
}

.batch_joint_composite_values <- function(metadata, candidates) {
  values <- lapply(candidates, function(candidate) .batch_values(metadata[[candidate]]))
  Reduce(function(left, right) paste(left, right, sep = "__"), values)
}

.batch_joint_model_info <- function(metadata, entries, active) {
  if (length(active) == 0L) return(NULL)
  complete <- rep(TRUE, nrow(metadata))
  for (entry_index in active) {
    values <- entries[[entry_index]]$values
    complete <- complete & !is.na(values)
  }
  complete_rows <- which(complete)
  if (length(complete_rows) < 3L) return(NULL)
  aliases <- paste0("joint_term_", seq_along(active))
  data <- data.frame(
    row.names = rownames(metadata)[complete_rows],
    stringsAsFactors = FALSE
  )
  for (index in seq_along(active)) {
    data[[aliases[[index]]]] <- factor(entries[[active[[index]]] ]$values[complete_rows])
  }
  formula <- stats::reformulate(aliases)
  model_matrix <- tryCatch(
    stats::model.matrix(formula, data = data),
    error = function(error) error
  )
  if (inherits(model_matrix, "error")) {
    return(list(error = model_matrix, complete = complete_rows, data = data, aliases = aliases))
  }
  assignment <- attr(model_matrix, "assign")
  rank_full <- qr(model_matrix)$rank
  term_estimable <- vapply(seq_along(aliases), function(index) {
    columns <- which(assignment == index)
    if (length(columns) == 0L) return(FALSE)
    rank_without <- qr(
      model_matrix[, setdiff(seq_len(ncol(model_matrix)), columns), drop = FALSE]
    )$rank
    rank_full - rank_without == length(columns)
  }, logical(1))
  list(
    error = NULL,
    complete = complete_rows,
    data = data,
    aliases = aliases,
    formula = formula,
    model_matrix = model_matrix,
    rank = rank_full,
    term_estimable = term_estimable
  )
}

prepare_batch_joint_design <- function(
  metadata,
  candidate_registry,
  near_unique_fraction = 0.50
) {
  registry <- .batch_registry_for_dataset(candidate_registry)
  if (length(near_unique_fraction) != 1L ||
      !is.finite(near_unique_fraction) ||
      near_unique_fraction <= 0 ||
      near_unique_fraction >= 1) {
    .batch_stop("near_unique_fraction must be one finite value in (0, 1)")
  }
  dataset <- if ("dataset" %in% colnames(registry)) registry$dataset[[1L]] else ""
  primary_row <- which(registry$is_primary)
  if (length(primary_row) != 1L) {
    .batch_stop("registry must contain exactly one primary biology column")
  }
  primary <- registry$candidate[[primary_row]]
  rows <- data.frame(
    candidate = registry$candidate,
    candidate_class = registry$candidate_class,
    is_primary = as.logical(registry$is_primary),
    formula_alias = registry$formula_alias,
    status = NA_character_,
    warning = NA_character_,
    joint_status = NA_character_,
    joint_warning = NA_character_,
    joint_model_term = NA_character_,
    joint_model_class = NA_character_,
    joint_reduction_group = NA_character_,
    n_samples = as.integer(nrow(metadata)),
    n_levels = integer(nrow(registry)),
    near_unique_fraction = NA_real_,
    stringsAsFactors = FALSE
  )
  validities <- vector("list", nrow(registry))
  for (index in seq_len(nrow(registry))) {
    candidate <- registry$candidate[[index]]
    validities[[index]] <- if (isTRUE(registry$is_primary[[index]])) {
      .batch_validity(metadata, candidate)
    } else {
      .batch_validity(metadata, candidate, primary)
    }
    values <- if (candidate %in% colnames(metadata)) .batch_values(metadata[[candidate]]) else character()
    rows$n_levels[[index]] <- length(unique(values[!is.na(values)]))
    rows$near_unique_fraction[[index]] <- if (nrow(metadata) > 0L) {
      rows$n_levels[[index]] / nrow(metadata)
    } else {
      NA_real_
    }
    rows$status[[index]] <- validities[[index]]$status
    rows$warning[[index]] <- validities[[index]]$warning
  }

  processing_order <- c(
    primary_row,
    which(!registry$is_primary & registry$candidate_class == "technical"),
    which(!registry$is_primary & registry$candidate_class != "technical")
  )
  entries <- list()
  nested_rules <- batch_joint_nested_rules(dataset)
  for (index in processing_order) {
    candidate <- registry$candidate[[index]]
    validity <- validities[[index]]
    values <- if (candidate %in% colnames(metadata)) .batch_values(metadata[[candidate]]) else character()
    if (!isTRUE(validity$valid) && !isTRUE(registry$is_primary[[index]])) {
      if (identical(validity$status, "NON_ESTIMABLE") &&
          grepl("perfectly confounded", validity$warning, fixed = TRUE)) {
        rows$joint_status[[index]] <- "DROPPED_CONFOUNDED_PRIMARY"
        rows$joint_warning[[index]] <- "candidate is perfectly confounded with primary biology"
      } else {
        rows$joint_status[[index]] <- "DROPPED_INVALID"
        rows$joint_warning[[index]] <- validity$warning
      }
      rows$status[[index]] <- rows$joint_status[[index]]
      next
    }
    if (!isTRUE(validity$valid) && isTRUE(registry$is_primary[[index]])) {
      rows$joint_status[[index]] <- "INVALID_PRIMARY"
      rows$joint_warning[[index]] <- validity$warning
      rows$status[[index]] <- rows$joint_status[[index]]
      next
    }
    if (!isTRUE(registry$is_primary[[index]]) &&
        rows$n_levels[[index]] / nrow(metadata) > near_unique_fraction) {
      rows$joint_status[[index]] <- "DROPPED_NEAR_UNIQUE"
      rows$joint_warning[[index]] <- paste0(
        "near-unique candidate: ", rows$n_levels[[index]], "/", nrow(metadata),
        " levels (threshold ", format(near_unique_fraction, trim = TRUE), ")"
      )
      rows$status[[index]] <- rows$joint_status[[index]]
      next
    }
    nested <- nested_rules[[candidate]]
    if (!is.null(nested) && any(vapply(
      entries,
      function(entry) nested$keep %in% entry$source_candidates,
      logical(1)
    ))) {
      rows$joint_status[[index]] <- "DROPPED_NESTED"
      rows$joint_warning[[index]] <- nested$reason
      rows$status[[index]] <- rows$joint_status[[index]]
      next
    }
    matching <- which(vapply(
      entries,
      function(entry) .batch_partition_equivalent(values, entry$partition_reference),
      logical(1)
    ))
    if (length(matching) > 0L) {
      entry_index <- matching[[1L]]
      entry <- entries[[entry_index]]
      if (isTRUE(entry$is_primary)) {
        rows$joint_status[[index]] <- "DROPPED_DUPLICATE_PRIMARY"
        rows$joint_warning[[index]] <- paste0(
          "exact metadata partition already represented by primary biology '",
          entry$term, "'"
        )
        rows$status[[index]] <- rows$joint_status[[index]]
        next
      }
      source_indices <- c(entry$source_indices, index)
      source_candidates <- c(entry$source_candidates, candidate)
      composite_term <- paste(source_candidates, collapse = "__")
      entry$source_indices <- source_indices
      entry$source_candidates <- source_candidates
      entry$term <- composite_term
      entry$values <- .batch_joint_composite_values(metadata, source_candidates)
      entry$class <- .batch_joint_model_class(registry, source_indices)
      entries[[entry_index]] <- entry
      rows$joint_status[source_indices] <- c(
        "RETAINED_COMPOSITE",
        rep("MERGED_EXACT_PARTITION", length(source_indices) - 1L)
      )
      rows$joint_warning[source_indices] <- paste0(
        "exact metadata partition merged as ", composite_term,
        "; original fields are not modeled separately"
      )
      rows$joint_model_term[source_indices] <- composite_term
      rows$joint_model_class[source_indices] <- entry$class
      rows$joint_reduction_group[source_indices] <- composite_term
      rows$status[source_indices] <- rows$joint_status[source_indices]
      next
    }
    entry_index <- length(entries) + 1L
    entries[[entry_index]] <- list(
      source_indices = index,
      source_candidates = candidate,
      term = candidate,
      values = values,
      partition_reference = values,
      class = registry$candidate_class[[index]],
      is_primary = isTRUE(registry$is_primary[[index]])
    )
    rows$joint_status[[index]] <- "RETAINED"
    rows$joint_model_term[[index]] <- candidate
    rows$joint_model_class[[index]] <- registry$candidate_class[[index]]
    rows$joint_reduction_group[[index]] <- candidate
    rows$status[[index]] <- rows$joint_status[[index]]
  }

  primary_entry <- which(vapply(entries, function(entry) isTRUE(entry$is_primary), logical(1)))
  if (length(primary_entry) != 1L) {
    rows$joint_design_status <- "UNAVAILABLE"
    rows$joint_design_warning <- "primary biology is not estimable"
    return(list(
      rows = rows,
      entries = entries,
      active = integer(),
      model_info = NULL,
      status = "UNAVAILABLE",
      warning = rows$joint_design_warning[[1L]],
      n_samples = 0L,
      reduced = TRUE
    ))
  }

  active <- seq_along(entries)
  repeat {
    info <- .batch_joint_model_info(metadata, entries, active)
    if (is.null(info)) break
    if (is.null(info$error) && all(info$term_estimable)) break
    if (is.null(info$error)) {
      bad <- which(!info$term_estimable)
    } else {
      bad <- seq_along(active)
    }
    bad <- bad[!vapply(active[bad], function(entry_index) isTRUE(entries[[entry_index]]$is_primary), logical(1))]
    if (length(bad) == 0L) break
    # Preserve the authoritative registry order (primary, technical, then
    # secondary biology). When an unrelated term is algebraically aliased,
    # drop the latest retained term rather than claiming semantic nesting.
    drop_position <- bad[[length(bad)]]
    drop_entry_index <- active[[drop_position]]
    entry <- entries[[drop_entry_index]]
    retained_terms <- vapply(active[active != drop_entry_index], function(entry_index) entries[[entry_index]]$term, character(1))
    warning_text <- paste0(
      "dropped from reduced joint model because it was aliased/rank-deficient; ",
      "retained terms: ", paste(retained_terms, collapse = ", ")
    )
    rows$joint_status[entry$source_indices] <- "DROPPED_ALIASED"
    rows$status[entry$source_indices] <- rows$joint_status[entry$source_indices]
    rows$joint_warning[entry$source_indices] <- warning_text
    active <- active[active != drop_entry_index]
  }
  info <- .batch_joint_model_info(metadata, entries, active)
  if (is.null(info) || is.null(info$error) == FALSE || any(!info$term_estimable)) {
    rows$joint_design_status <- "NON_ESTIMABLE"
    rows$joint_design_warning <- "reduced joint design remains rank-deficient after deterministic reductions"
    return(list(
      rows = rows,
      entries = entries,
      active = active,
      model_info = info,
      status = "NON_ESTIMABLE",
      warning = rows$joint_design_warning[[1L]],
      n_samples = if (is.null(info)) 0L else length(info$complete),
      reduced = TRUE
    ))
  }
  for (entry_index in active) {
    entry <- entries[[entry_index]]
    rows$joint_model_term[entry$source_indices] <- entry$term
    rows$joint_model_class[entry$source_indices] <- entry$class
    rows$joint_reduction_group[entry$source_indices] <- entry$term
  }
  reduced <- any(!rows$joint_status %in% c("RETAINED", "RETAINED_COMPOSITE"))
  design_status <- if (reduced) "ESTIMABLE_REDUCED" else "ESTIMABLE"
  design_warning <- .batch_join_warnings(rows$joint_warning)
  rows$joint_design_status <- design_status
  rows$joint_design_warning <- design_warning
  list(
    rows = rows,
    entries = entries,
    active = active,
    model_info = info,
    status = design_status,
    warning = design_warning,
    n_samples = length(info$complete),
    reduced = reduced
  )
}
.batch_joint_design <- function(
  dist_mat,
  metadata,
  registry,
  permutations = 999L,
  parallel = getOption("mc.cores")
) {
  reduction <- prepare_batch_joint_design(metadata, registry)
  status_table <- reduction$rows
  if (identical(reduction$status, "UNAVAILABLE")) {
    return(list(
      rows = status_table,
      decomposition = NULL,
      status = reduction$status,
      warning = reduction$warning,
      n_samples = reduction$n_samples,
      full_model_r2 = NA_real_
    ))
  }
  if (identical(reduction$status, "NON_ESTIMABLE")) {
    return(list(
      rows = status_table,
      decomposition = NULL,
      status = reduction$status,
      warning = reduction$warning,
      n_samples = reduction$n_samples,
      full_model_r2 = NA_real_
    ))
  }
  info <- reduction$model_info
  active <- reduction$active
  distance_matrix <- as.matrix(if (inherits(dist_mat, "dist")) dist_mat else dist_mat)
  distance <- stats::as.dist(distance_matrix[info$complete, info$complete, drop = FALSE])
  term_formula <- stats::reformulate(info$aliases, response = "distance")
  environment(term_formula) <- environment()
  warnings <- character()
  full_model <- tryCatch(
    withCallingHandlers(
      vegan::adonis2(
        term_formula,
        data = info$data,
        permutations = as.integer(permutations),
        parallel = parallel
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  marginal_model <- tryCatch(
    withCallingHandlers(
      vegan::adonis2(
        term_formula,
        data = info$data,
        permutations = as.integer(permutations),
        by = "margin",
        parallel = parallel
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  if (inherits(full_model, "error") || inherits(marginal_model, "error")) {
    error <- if (inherits(full_model, "error")) full_model else marginal_model
    status_table$status <- "NON_ESTIMABLE"
    status_table$joint_status <- "NON_ESTIMABLE"
    status_table$warning <- .batch_join_warnings(status_table$warning, conditionMessage(error))
    return(list(
      rows = status_table,
      decomposition = NULL,
      status = "NON_ESTIMABLE",
      warning = conditionMessage(error),
      n_samples = length(info$complete),
      full_model_r2 = NA_real_
    ))
  }

  full_model_r2 <- as.numeric(full_model$R2[[1L]])
  status_table$r2 <- NA_real_
  status_table$pseudo_f <- NA_real_
  status_table$p_value <- NA_real_
  status_table$p_adjusted_holm <- NA_real_
  for (index in seq_along(active)) {
    entry <- reduction$entries[[active[[index]]]]
    source_row <- entry$source_indices[[1L]]
    result_row <- which(rownames(marginal_model) == info$aliases[[index]])
    if (length(result_row) != 1L || !is.finite(marginal_model$R2[[result_row]])) {
      status_table$status[entry$source_indices] <- "NON_ESTIMABLE"
      status_table$joint_status[entry$source_indices] <- "NON_ESTIMABLE"
      status_table$joint_warning[entry$source_indices] <- "term was not returned by adonis2"
      next
    }
    status_table$status[[source_row]] <- "ESTIMABLE"
    status_table$r2[[source_row]] <- as.numeric(marginal_model$R2[[result_row]])
    status_table$pseudo_f[[source_row]] <- as.numeric(marginal_model$F[[result_row]])
    status_table$p_value[[source_row]] <- as.numeric(marginal_model$`Pr(>F)`[[result_row]])
  }
  status_table$n_samples <- as.integer(length(info$complete))
  status_table$full_model_r2 <- full_model_r2
  adjusted <- status_table$status == "ESTIMABLE" & is.finite(status_table$p_value)
  if (any(adjusted)) {
    status_table$p_adjusted_holm[adjusted] <- stats::p.adjust(
      status_table$p_value[adjusted],
      method = "holm"
    )
  }
  status_table$warning <- mapply(
    .batch_join_warnings,
    status_table$warning,
    status_table$joint_warning,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  status_table$warning[is.na(status_table$warning)] <- ""
  term_r2 <- vapply(active, function(entry_index) {
    source_row <- reduction$entries[[entry_index]]$source_indices[[1L]]
    status_table$r2[[source_row]]
  }, numeric(1))
  term_classes <- vapply(active, function(entry_index) reduction$entries[[entry_index]]$class, character(1))
  unique_biological <- sum(term_r2[term_classes %in% c("biological_primary", "biological_secondary")], na.rm = TRUE)
  unique_technical <- sum(term_r2[term_classes == "technical"], na.rm = TRUE)
  shared <- max(0, full_model_r2 - unique_biological - unique_technical)
  residual <- max(0, 1 - full_model_r2)
  component_warning <- .batch_join_warnings(reduction$warning, paste(warnings, collapse = "; "))
  components <- data.frame(
    component = c("Unique Biological", "Unique Technical", "Shared / Confounded", "Residual / Unexplained"),
    r2 = c(unique_biological, unique_technical, shared, residual),
    status = reduction$status,
    warning = component_warning,
    n_samples = as.integer(length(info$complete)),
    full_model_r2 = full_model_r2,
    stringsAsFactors = FALSE
  )
  status_table$joint_design_warning <- component_warning
  list(
    rows = status_table,
    decomposition = components,
    status = reduction$status,
    warning = component_warning,
    n_samples = length(info$complete),
    full_model_r2 = full_model_r2
  )
}

compute_batch_joint_permanova <- function(
  dist_mat,
  metadata,
  candidate_registry,
  permutations = 999L,
  parallel = getOption("mc.cores")
) {
  .batch_joint_design(dist_mat, metadata, candidate_registry, permutations, parallel)
}

compute_batch_nmi <- function(metadata, candidate_registry) {
  registry <- .batch_registry_for_dataset(candidate_registry)
  candidates <- registry$candidate
  matrix <- matrix(NA_real_, nrow = length(candidates), ncol = length(candidates), dimnames = list(candidates, candidates))
  .nmi <- function(x, y) {
    x <- .batch_values(x); y <- .batch_values(y)
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]; y <- y[keep]
    if (length(x) == 0L || length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
    table_xy <- table(x, y)
    pxy <- table_xy / sum(table_xy)
    px <- rowSums(pxy); py <- colSums(pxy)
    nonzero <- pxy > 0
    mutual_information <- sum(pxy[nonzero] * log(pxy[nonzero] / outer(px, py)[nonzero]))
    entropy_x <- -sum(px[px > 0] * log(px[px > 0]))
    entropy_y <- -sum(py[py > 0] * log(py[py > 0]))
    if (entropy_x <= 0 || entropy_y <= 0) return(NA_real_)
    max(0, min(1, mutual_information / sqrt(entropy_x * entropy_y)))
  }
  for (i in seq_along(candidates)) {
    if (candidates[[i]] %in% colnames(metadata)) matrix[i, i] <- 1
    if (i < length(candidates)) {
      for (j in seq.int(i + 1L, length(candidates))) {
        score <- if (candidates[[i]] %in% colnames(metadata) && candidates[[j]] %in% colnames(metadata)) {
          .nmi(metadata[[candidates[[i]]]], metadata[[candidates[[j]]]])
        } else {
          NA_real_
        }
        matrix[i, j] <- score
        matrix[j, i] <- score
      }
    }
  }
  matrix
}

plot_batch_nmi_heatmap <- function(nmi_matrix, dataset = "", threshold = 0.70) {
  if (!is.matrix(nmi_matrix) || nrow(nmi_matrix) != ncol(nmi_matrix)) .batch_stop("NMI heatmap requires a square matrix")
  ordered <- rownames(nmi_matrix)
  if (is.null(ordered) || !identical(ordered, colnames(nmi_matrix))) .batch_stop("NMI matrix row/column labels differ")
  frame <- expand.grid(row = ordered, column = ordered, stringsAsFactors = FALSE)
  frame$nmi <- as.vector(nmi_matrix[cbind(match(frame$row, ordered), match(frame$column, ordered))])
  frame$label <- ifelse(is.na(frame$nmi), "NA", sprintf("%.2f", frame$nmi))
  frame$severe <- !is.na(frame$nmi) & frame$nmi > threshold
  ggplot2::ggplot(frame, ggplot2::aes(x = column, y = row, fill = nmi)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2, na.rm = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5, na.rm = FALSE) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#B2182B", limits = c(0, 1), na.value = "grey80",
      name = "NMI"
    ) +
    ggplot2::scale_x_discrete(limits = ordered, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = rev(ordered), drop = FALSE) +
    ggplot2::labs(
      title = paste0(dataset, " metadata NMI"),
      subtitle = paste0("Metadata-only; NMI > ", format(threshold, nsmall = 2), " indicates severe collinearity"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
}

make_batch_metric_table <- function(
  dataset,
  method,
  metadata,
  candidate_registry,
  dist_mat,
  artifact_path,
  permutations = 999L,
  parallel = getOption("mc.cores"),
  permanova_permutations = permutations,
  precomputed_anosim = NULL
) {
  registry <- .batch_registry_for_dataset(candidate_registry, dataset)
  primary <- registry$candidate[registry$is_primary][[1L]]
  rows <- lapply(seq_len(nrow(registry)), function(index) {
    candidate <- registry$candidate[[index]]
    cached_primary_anosim <- identical(candidate, primary) &&
      length(precomputed_anosim) == 1L &&
      is.finite(precomputed_anosim)
    anosim <- if (cached_primary_anosim) {
      list(
        statistic = as.numeric(precomputed_anosim),
        p_value = NA_real_,
        status = "ESTIMABLE",
        warning = "primary ANOSIM statistic reused from persisted/cached heatmap metric; permutation p-value not persisted",
        n_samples = nrow(metadata)
      )
    } else if (candidate %in% colnames(metadata)) {
      compute_batch_anosim(dist_mat, metadata[[candidate]], permutations, parallel)
    } else {
      list(statistic = NA_real_, p_value = NA_real_, status = "MISSING", warning = paste0("missing candidate column ", candidate), n_samples = 0L)
    }
    univariate <- compute_batch_permanova(dist_mat, metadata, candidate, NULL, permanova_permutations, parallel)
    informed <- if (identical(candidate, primary)) {
      .batch_empty_metric("NOT_A_CANDIDATE", "primary biology omitted from biology-informed candidate panel", candidate = candidate)
    } else {
      compute_batch_permanova(dist_mat, metadata, candidate, primary, permanova_permutations, parallel)
    }
    data.frame(
      dataset = dataset,
      method = method,
      method_label = if (exists("method_label_map_main", inherits = TRUE) && method %in% names(method_label_map_main)) unname(method_label_map_main[[method]]) else method,
      candidate = candidate,
      candidate_class = registry$candidate_class[[index]],
      is_primary = isTRUE(registry$is_primary[[index]]),
      sample_count = nrow(metadata),
      anosim_statistic = anosim$statistic,
      anosim_p_value = anosim$p_value,
      anosim_status = anosim$status,
      anosim_warning = anosim$warning,
      univariate_r2 = univariate$r2,
      univariate_pseudo_f = univariate$pseudo_f,
      univariate_p_value = univariate$p_value,
      univariate_status = univariate$status,
      univariate_warning = univariate$warning,
      biology_informed_r2 = informed$r2,
      biology_informed_pseudo_f = informed$pseudo_f,
      biology_informed_p_value = informed$p_value,
      biology_informed_status = informed$status,
      biology_informed_warning = informed$warning,
      artifact_path = artifact_path,
      warning = .batch_join_warnings(anosim$warning, univariate$warning, informed$warning),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$biology_informed_p_adjusted_holm <- NA_real_
  valid_univariate <- result$univariate_status == "ESTIMABLE" & is.finite(result$univariate_p_value)
  valid_informed <- result$biology_informed_status == "ESTIMABLE" & is.finite(result$biology_informed_p_value)
  if (any(valid_univariate)) result$univariate_p_adjusted_holm[valid_univariate] <- stats::p.adjust(result$univariate_p_value[valid_univariate], "holm")
  if (any(valid_informed)) result$biology_informed_p_adjusted_holm[valid_informed] <- stats::p.adjust(result$biology_informed_p_value[valid_informed], "holm")
  result
}

make_batch_joint_table <- function(dataset, method, joint_result, artifact_path) {
  rows <- joint_result$rows
  numeric_columns <- c(
    "r2", "pseudo_f", "p_value", "p_adjusted_holm", "n_samples",
    "full_model_r2"
  )
  for (column in numeric_columns) {
    if (!column %in% colnames(rows)) rows[[column]] <- NA_real_
  }
  character_columns <- c(
    "joint_status", "joint_warning", "joint_model_term",
    "joint_model_class", "joint_reduction_group",
    "joint_design_status", "joint_design_warning"
  )
  for (column in character_columns) {
    if (!column %in% colnames(rows)) rows[[column]] <- NA_character_
  }
  rows$dataset <- dataset
  rows$method <- method
  rows$method_label <- if (exists("method_label_map_main", inherits = TRUE) && method %in% names(method_label_map_main)) unname(method_label_map_main[[method]]) else method
  rows$sample_count <- as.integer(joint_result$n_samples)
  rows$full_model_r2 <- ifelse(
    is.finite(rows$full_model_r2),
    rows$full_model_r2,
    as.numeric(joint_result$full_model_r2)
  )
  rows$reduction_status <- rows$joint_status
  rows$reduction_warning <- rows$joint_warning
  rows$joint_design_status <- ifelse(
    nzchar(rows$joint_design_status),
    rows$joint_design_status,
    as.character(joint_result$status)
  )
  rows$joint_design_warning <- ifelse(
    nzchar(rows$joint_design_warning),
    rows$joint_design_warning,
    as.character(joint_result$warning)
  )
  rows$artifact_path <- artifact_path
  rows
}

make_batch_nmi_table <- function(dataset, nmi_matrix, candidate_registry) {
  registry <- .batch_registry_for_dataset(candidate_registry)
  ordered <- registry$candidate
  pairs <- list()
  row_index <- 0L
  for (i in seq_along(ordered)) {
    for (j in seq_along(ordered)) {
      row_index <- row_index + 1L
      score <- nmi_matrix[ordered[[i]], ordered[[j]]]
      pairs[[row_index]] <- data.frame(
        dataset = dataset,
        candidate_a = ordered[[i]],
        candidate_b = ordered[[j]],
        candidate_a_class = registry$candidate_class[[i]],
        candidate_b_class = registry$candidate_class[[j]],
        nmi = as.numeric(score),
        severe_collinearity = !identical(ordered[[i]], ordered[[j]]) &&
          !is.na(score) && score > 0.70,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(pairs) == 0L) return(data.frame())
  do.call(rbind, pairs)
}
write_batch_analysis_table <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", as.integer(stats::runif(1, 1, 1e9)))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) .batch_stop("empty analysis table temporary file: ", tmp)
  if (!file.rename(tmp, path)) {
    if (file.exists(path)) unlink(path)
    if (!file.rename(tmp, path)) .batch_stop("could not atomically install analysis table: ", path)
  }
  .batch_write_checksum(path)
  validate_batch_artifact(path, "analysis table")
  invisible(path)
}
