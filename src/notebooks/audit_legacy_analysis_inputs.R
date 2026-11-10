# Evidence only: no fits, score calculation, publication receipts or source edits.
# Source this file for the focused fixture; the CLI is invoked only on HPC.
.legacy_need <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
.legacy_fields <- function(frame, fields) {
  .legacy_need(is.data.frame(frame) && nrow(frame) > 0L && all(fields %in% names(frame)),
    paste("Missing finite manifest fields:", paste(fields, collapse = ", ")))
}
.legacy_ids <- function(ids) {
  !is.null(ids) && length(ids) > 0L && !anyNA(ids) &&
    all(nzchar(trimws(ids))) && !anyDuplicated(ids)
}
.legacy_checked_read <- function(path, md5, size, reader) {
  .analysis_checked_file(path, md5)
  .legacy_need(identical(as.numeric(file.info(path)$size), as.numeric(size)), "Manifest size mismatch")
  value <- reader(path)
  .analysis_checked_file(path, md5)
  .legacy_need(identical(as.numeric(file.info(path)$size), as.numeric(size)), "Input changed during read")
  value
}
.legacy_input_witness <- function(row) {
  path <- as.character(row$current_input_path); proof <- as.character(row$current_input_proof)
  for (p in c(path, proof)) .legacy_need(file.exists(p) && !dir.exists(p) &&
    identical(normalizePath(p, mustWork = TRUE), p), "Noncanonical current-input path/proof")
  .legacy_need(file.info(path)$size == as.numeric(row$current_input_size) &&
    file.info(proof)$size > 0 && file.info(proof)$size < 16384, "Invalid current-input stat/proof size")
  lines <- readLines(proof, warn = FALSE)
  .legacy_need(length(lines) == 3L && identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH")) &&
    identical(sub("^MD5=", "", lines[1]), as.character(row$current_input_md5)) &&
    identical(sub("^SIZE=", "", lines[2]), as.character(row$current_input_size)) &&
    identical(sub("^PATH=", "", lines[3]), path), "Current-input checksum proof disagrees with manifest")
  list(path = path, md5 = as.character(row$current_input_md5), size = as.numeric(row$current_input_size),
    proof_path = proof, proof_md5 = unname(tools::md5sum(proof)),
    mtime = as.numeric(file.info(path)$mtime), ctime = as.numeric(file.info(path)$ctime), mode = as.character(file.info(path)$mode),
    identity_guarantee = "strict proof and size checked; no full H5AD rehash here; prior content-validation witness required")
}
.legacy_read_obs <- function(path, columns, module_dir) {
  # Reuse the h5py column decoder directly: no X/layers/raw/var/obsm values.
  python <- Sys.getenv("RETICULATE_PYTHON")
  .legacy_need(grepl("^/", python) && file.exists(python), "Explicit configured RETICULATE_PYTHON required")
  h5py <- reticulate::import("h5py", convert = FALSE)
  decoder <- reticulate::import_from_path("h5ad_source_identity", path = module_dir, convert = FALSE)
  handle <- h5py$File(path, "r"); on.exit(handle$close())
  group <- handle$get("obs")
  .legacy_need(!reticulate::py_is_null_xptr(group), "Missing H5AD obs")
  # Stage5 uses the existing canonical Sample derived from the configured
  # source column. This additional witness is metadata, never expression.
  if (!"Sample" %in% columns && isTRUE(reticulate::py_to_r(group$`__contains__`("Sample"))))
    columns <- c(columns, "Sample")
  values <- setNames(lapply(columns, function(column)
    as.character(reticulate::py_to_r(decoder$read_obs_column(group, column)))), columns)
  .legacy_need(length(unique(lengths(values))) == 1L && lengths(values)[1] > 0L, "Invalid obs column lengths")
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}
.legacy_canonical_obs <- function(obs, entry) {
  raw <- as.character(obs[[entry$sample_col]])
  .legacy_need(length(raw) == nrow(obs) && !anyNA(raw) && all(nzchar(trimws(raw))), "Invalid configured sample IDs")
  original <- unique(raw)
  # Configured Sample is already authoritative, including legitimate hyphens.
  # Only a distinct configured raw alias needs the producer's mapping witness.
  if (identical(entry$sample_col, "Sample"))
    return(list(obs = obs, aliases = data.frame(raw = original, canonical = original, stringsAsFactors = FALSE)))
  expected <- unname(standardize_sample_names(original))
  .legacy_need(.legacy_ids(expected), "Configured sample IDs collide after production canonicalization")
  if ("Sample" %in% names(obs)) {
    aliases <- unique(data.frame(raw = raw, canonical = as.character(obs$Sample), stringsAsFactors = FALSE))
    .legacy_need(.legacy_ids(aliases$raw) && .legacy_ids(aliases$canonical) &&
      identical(unname(aliases$canonical), expected), "Configured-to-canonical Sample witness is not the production bijection")
    obs[[entry$sample_col]] <- as.character(obs$Sample)
  } else {
    .legacy_need(identical(original, expected), "Changed sample IDs require an existing canonical Sample witness")
    aliases <- data.frame(raw = original, canonical = original, stringsAsFactors = FALSE)
  }
  list(obs = obs, aliases = aliases)
}
.legacy_labels <- function(obs, entry) {
  .legacy_need(all(c(entry$sample_col, entry$label_col) %in% names(obs)), "Missing configured sample/label column")
  sample <- as.character(obs[[entry$sample_col]]); label <- as.character(obs[[entry$label_col]])
  .legacy_need(!anyNA(sample) && all(nzchar(trimws(sample))) &&
    !anyNA(label) && all(nzchar(trimws(label))) &&
    !any(sample %in% c("NA", "nan", "None")) && !any(label %in% c("NA", "nan", "None")), "Missing sample/biological label")
  ids <- unique(sample)
  labels <- vapply(ids, function(id) {
    x <- unique(label[sample == id]); .legacy_need(length(x) == 1L, "Within-sample biological labels disagree"); x
  }, character(1))
  setNames(labels, ids)
}
.legacy_bundle <- function(value, row, labels) {
  .legacy_need(is.list(value), "Selected result is not a bundle")
  distance <- evaluation_distance(value$dist_mat); ids <- unname(attr(distance, "Labels"))
  .legacy_need(all(ids %in% names(labels)), "Distance has extra sample IDs")
  dropped <- setdiff(names(labels), ids); scitd <- grepl("^scITD(_|$)", row$method_key)
  .legacy_need(scitd || !length(dropped), "Only scITD may drop samples")
  .legacy_need(length(unique(labels[ids])) >= 2L, "Retained samples have fewer than two biological groups")
  .legacy_need(.legacy_ids(names(value$labels)) && identical(unname(names(value$labels)), ids) &&
    !anyNA(value$labels) && all(nzchar(trimws(as.character(value$labels)))), "Stored named labels must match distance order")
  null <- grepl("(^|_)NULL($|_)", row$method_key)
  if (null) .legacy_need(identical(sort(as.character(value$labels)), sort(unname(labels[ids]))),
    "Stored NULL shuffle does not preserve biological category counts")
  else .legacy_need(identical(as.character(value$labels), unname(labels[ids])), "Stored biological labels disagree")
  scores <- value$scores; present <- intersect(evaluation_metrics(), names(scores)); policy <- value$scoring_policy
  valid <- present[vapply(scores[present], function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))]
  if (identical(policy, evaluation_scoring_policy()) && "lisi_score" %in% valid &&
      (scores$lisi_score < 0 || scores$lisi_score > 1)) valid <- setdiff(valid, "lisi_score")
  if (anyDuplicated(names(scores))) valid <- character()
  timing_names <- intersect(c("exec_time", "time_secs", "mem_GB", "timing_schema", "timing_id", "shared_timing_id",
    "shared_time_secs", "shared_mem_GB", "variant_time_secs", "variant_mem_GB", "aggregate_time_secs", "shared_fit_time_secs"), names(value))
  timing_valid <- timing_names[vapply(value[timing_names], function(x) length(x) == 1L &&
    ((is.numeric(x) && is.finite(x) && x >= 0) || (is.character(x) && !is.na(x) && nzchar(x))), logical(1))]
  list(sample_ids = ids, sample_order_matches_current = identical(ids, unname(names(labels))),
    stored_labels = setNames(as.character(value$labels), ids), labels_md5 = .analysis_fingerprint(value$labels),
    distance_md5 = .analysis_fingerprint(value$dist_mat), n_samples = length(ids),
    n_features = if (is.null(value$feat_mat)) NULL else ncol(value$feat_mat),
    dropped_sample_ids = dropped, dropped_ids_require_explicit_review = length(dropped) > 0L,
    score_names = names(scores), missing_core_scores = setdiff(evaluation_metrics(), present),
    invalid_core_scores = setdiff(present, valid), stored_core_scores = scores[valid], scoring_policy = policy,
    original_timing = value[timing_valid], invalid_timing_fields = setdiff(timing_names, timing_valid),
    score_disposition = if (length(valid) == 5L && identical(policy, evaluation_scoring_policy()))
      "stored_policy_claim_requires_root_review" else "needs_validated_distance_rescore",
    historical_producer = "unproven", historical_runtime = "unproven", origin = "legacy_reuse")
}
.legacy_pair <- function(author, metadata, obs, entry, labels) {
  counts <- as.matrix(get_ct_comp_df(obs, entry$sample_col, entry$cell_type_high_res))
  .legacy_need(setequal(rownames(counts), names(labels)) && ncol(counts) > 0L, "Author annotations omit whole samples")
  stored <- as.matrix(author$counts); features <- as.matrix(author$feat_mat)
  .legacy_need(.legacy_ids(rownames(stored)) && .legacy_ids(colnames(stored)) &&
    setequal(rownames(stored), rownames(counts)) && setequal(colnames(stored), colnames(counts)), "Author count axes differ")
  counts <- counts[rownames(stored), colnames(stored), drop = FALSE]
  .legacy_need(is.numeric(stored) && all(is.finite(stored)) &&
    isTRUE(all.equal(unname(stored), unname(counts + .5), tolerance = 1e-12)), "Author counts are not current raw counts +0.5")
  .legacy_need(.legacy_ids(rownames(features)) && .legacy_ids(colnames(features)) &&
    setequal(rownames(features), rownames(stored)) && setequal(colnames(features), colnames(stored)), "Author CLR feature axes differ")
  clr <- log(stored) - rowMeans(log(stored))
  # Align the expected view only; original axes, values and fingerprints remain.
  clr <- clr[rownames(features), colnames(features), drop = FALSE]
  .legacy_need(isTRUE(all.equal(unname(features), unname(clr), tolerance = 1e-10)), "Author CLR features differ")
  distance <- as.matrix(evaluation_distance(author$dist_mat))
  .legacy_need(identical(unname(rownames(distance)), unname(rownames(features))) &&
    isTRUE(all.equal(unname(distance), unname(as.matrix(stats::dist(clr))), tolerance = 1e-10)), "Author distances differ from current CLR")
  sample <- as.character(obs[[entry$sample_col]])
  cells <- setNames(as.numeric(table(factor(sample, levels = names(labels)))), names(labels))
  .legacy_need(setequal(names(metadata$labels), names(labels)) &&
    identical(as.character(metadata$labels[names(labels)]), unname(labels)) &&
    identical(as.integer(metadata$n_samples), as.integer(length(labels))) &&
    identical(as.numeric(metadata$n_cells), as.numeric(nrow(obs))) &&
    setequal(names(metadata$cells_per_sample), names(labels)) &&
    identical(as.numeric(metadata$cells_per_sample[names(labels)]), unname(cells)), "Dataset metadata disagrees with current obs")
  categories <- length(unique(obs[[entry$cell_type_high_res]]))
  .legacy_need(identical(as.integer(metadata$n_cell_types_high_res), as.integer(categories)), "Author category metadata differs")
  list(status = "COMPATIBLE_PENDING_REVIEW", n_cells = nrow(obs), n_samples = length(labels),
    sample_ids = unname(names(metadata$cells_per_sample)),
    valid_annotated_cells = sum(counts), n_valid_categories = ncol(counts), n_recorded_categories = categories,
    counts_policy = "observed current get_ct_comp_df +0.5; CLR and Euclidean distance checked; no scores calculated",
    raw_counts_md5 = .analysis_fingerprint(counts), cell_subsetting_source_eligible_pending_review = TRUE)
}

.legacy_header_scope <- function(rows, observation) {
  selected <- rows[rows$kind %in% c("distance_feather", "embedding_feather"), , drop = FALSE]
  .legacy_need(identical(observation$schema, "ecoda_p6_feather_header_observation_v1") &&
    length(observation$rows) == nrow(selected), "Feather header observation schema/scope differs")
  observed_paths <- vapply(observation$rows, `[[`, character(1), "artifact_path")
  .legacy_need(!anyDuplicated(observed_paths) && setequal(observed_paths, selected$artifact_path), "Feather header file union differs")
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]; found <- observation$rows[[match(row$artifact_path, observed_paths)]]
    for (field in c("dataset", "method_key", "kind", "md5", "size"))
      .legacy_need(identical(as.character(found[[field]]), as.character(row[[field]])), "Feather header row identity differs")
    fields <- vapply(found$fields, `[[`, character(1), "name")
    .legacy_need(identical(unname(unlist(found$pandas$index_columns)), as.character(row$sample_id_col)) &&
      sum(fields == row$sample_id_col) == 1L && !anyDuplicated(fields) && identical(found$record_batches_read, 0L),
      "Declared Feather sample ID lacks exact pandas-index schema evidence")
  }
  invisible(TRUE)
}

audit_legacy_analysis_inputs <- function(rows, files, pairs, datasets, obs_reader) {
  .legacy_fields(rows, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "sample_id_col", "md5", "size",
    "current_input_path", "current_input_md5", "current_input_size", "current_input_proof", "origin", "status"))
  .legacy_fields(files, c("dataset", "artifact_path", "md5", "size", "needed_bundle_keys"))
  .legacy_fields(pairs, c("dataset", "author_path", "author_md5", "author_bundle_key", "metadata_path", "metadata_md5", "current_input_path", "current_input_md5"))
  .legacy_need(!anyDuplicated(.analysis_key(rows)) && !anyDuplicated(files$artifact_path) &&
    !anyDuplicated(pairs$dataset) && all(rows$origin == "legacy_reuse") &&
    all(rows$status == "BLOCKED_UNTIL_SCIENTIFIC_REVIEW") && all(rows$view == "benchmark_analysis"), "Invalid/approved/duplicate legacy scope")
  .legacy_need(all(rows$kind %in% c("rds_bundle", "standalone_result", "distance_feather", "embedding_feather")) &&
    all((rows$kind == "rds_bundle") == nzchar(rows$bundle_key)) && all(pairs$dataset %in% rows$dataset), "Unknown artifact kind/pair cohort")
  feather_rows <- rows$kind %in% c("distance_feather", "embedding_feather")
  .legacy_need(!anyNA(rows$sample_id_col[feather_rows]) && all(nzchar(rows$sample_id_col[feather_rows])),
    "Feather scope requires explicit nonempty sample_id_col from reviewed schema evidence")
  .legacy_need(all(rows$dataset %in% names(datasets)) &&
    !any(grepl("zeroimp|_trans\\.rds|\\.h5ad$|/pseudobulks/", rows$artifact_path)) &&
    !any(grepl("_bl$|^ECODA_PB_combo_", rows$method_key)), "Forbidden artifact in legacy scope")
  .legacy_need(setequal(files$artifact_path[files$dataset != "__all__"], rows$artifact_path), "Physical/logical scope mismatch")
  for (i in seq_len(nrow(files))) {
    f <- files[i, ]; selected <- rows[rows$artifact_path == f$artifact_path, , drop = FALSE]
    .legacy_need(all(selected$md5 == f$md5 & selected$size == f$size & selected$dataset == f$dataset), "Physical/logical identity mismatch")
    expected <- unique(selected$bundle_key[selected$kind == "rds_bundle"])
    hints <- if (nzchar(f$needed_bundle_keys)) strsplit(f$needed_bundle_keys, ";", fixed = TRUE)[[1]] else character()
    .legacy_need(!anyDuplicated(hints) && setequal(hints, expected), "Physical needed_bundle_keys differs from exact selected union")
  }
  report <- list(schema = "ecoda_legacy_analysis_compatibility_v1", scientifically_approved_rows = 0L,
    model_fits = 0L, score_calculations = 0L, h5ad_expression_values_read = FALSE,
    h5ad_content_rehashed_here = FALSE, prior_validated_content_witness_required = TRUE,
    physical_files = nrow(files), logical_rows = nrow(rows), rows = list(), author_hr_pairs = list(), cohorts = list())
  for (dataset in unique(rows$dataset)) {
    selected <- rows[rows$dataset == dataset, , drop = FALSE]; entry <- datasets[[dataset]]
    .legacy_need(isTRUE(entry$use_for_benchmark), "Dataset not enabled for benchmark")
    input_fields <- c("current_input_path", "current_input_md5", "current_input_size", "current_input_proof")
    .legacy_need(nrow(unique(selected[input_fields])) == 1L, "Cohort rows disagree on current input")
    before <- .legacy_input_witness(selected[1, ])
    has_pair <- dataset %in% pairs$dataset
    columns <- unique(c(entry$sample_col, entry$label_col, if (has_pair) entry$cell_type_high_res))
    obs <- obs_reader(before$path, columns)
    canonical <- .legacy_canonical_obs(obs, entry); obs <- canonical$obs
    labels <- .legacy_labels(obs, entry)
    report$cohorts[[dataset]] <- c(before, list(sample_col = entry$sample_col, label_col = entry$label_col,
      author_hr_col = if (has_pair) entry$cell_type_high_res else NULL,
      configured_sample_ids = canonical$aliases$raw, canonical_sample_mapping = canonical$aliases,
      canonical_sample_column = if ("Sample" %in% names(obs)) "Sample" else entry$sample_col,
      sample_ids = names(labels), n_cells = nrow(obs)))
    objects <- list(); author <- metadata <- NULL
    for (path in unique(selected$artifact_path)) {
      f <- files[match(path, files$artifact_path), ]; subset <- selected[selected$artifact_path == path, , drop = FALSE]
      reader <- if (grepl("\\.rds$", path)) readRDS else function(p) as.data.frame(arrow::read_feather(p))
      object <- .legacy_checked_read(path, f$md5, f$size, reader)
      value <- NULL
      if (any(subset$kind == "rds_bundle")) .legacy_need(is.list(object) && !is.null(names(object)) &&
        !anyDuplicated(names(object)), "Malformed aggregate bundle keys")
      for (i in seq_len(nrow(subset))) {
        row <- subset[i, ]; key <- paste(dataset, row$method_key, sep = "/")
        result <- tryCatch({
          if (row$kind == "rds_bundle") {
            .legacy_need(row$bundle_key %in% names(object), "Missing selected aggregate key")
            value <- object[[row$bundle_key]]
          } else value <- object
          if (row$method_key == "__dataset_metadata__") {
            metadata <- value; list(status = "PAIR_CHECK_PENDING")
          } else if (row$kind %in% c("distance_feather", "embedding_feather")) {
            .legacy_need(is.data.frame(value) && ncol(value) >= 2L && !anyDuplicated(names(value)) &&
              row$sample_id_col %in% names(value), "Feather lacks its exact declared sample_id_col")
            id_col <- row$sample_id_col; ids <- as.character(value[[id_col]])
            .legacy_need(.legacy_ids(ids) && setequal(ids, names(labels)), "Feather sample IDs differ")
            x <- as.matrix(value[, names(value) != id_col, drop = FALSE]); rownames(x) <- ids
            .legacy_need(is.numeric(x) && all(is.finite(x)), "Nonfinite Feather values")
            if (row$kind == "distance_feather") evaluation_distance(x)
            list(status = "COMPATIBLE_PENDING_REVIEW", sample_id_col = id_col, sample_ids = ids,
              sample_order_matches_current = identical(ids, names(labels)), shape = dim(x),
              value_md5 = .analysis_fingerprint(x), stored_labels_claim = FALSE,
              score_disposition = "needs_validated_distance_rescore", origin = "legacy_reuse")
          } else {
            details <- .legacy_bundle(value, row, labels)
            if (row$method_key == "ECODA_authors_HR") author <- value
            if (grepl("(^|_)NULL($|_)", row$method_key)) objects[[row$method_key]] <- value
            c(list(status = "COMPATIBLE_PENDING_REVIEW"), details)
          }
        }, error = function(e) list(status = "INCOMPATIBLE", reason = conditionMessage(e)))
        report$rows[[key]] <- c(list(dataset = dataset, view = row$view, method_key = row$method_key,
          artifact_path = path, artifact_md5 = row$md5, actual_aggregate_keys = if (row$kind == "rds_bundle") names(object) else NULL), result)
      }
      .analysis_checked_file(path, f$md5); rm(object, value); invisible(gc(FALSE))
    }
    for (key in names(objects)) {
      same <- !is.null(author) && identical(objects[[key]]$feat_mat, author$feat_mat) && identical(objects[[key]]$dist_mat, author$dist_mat)
      name <- paste(dataset, key, sep = "/"); report$rows[[name]]$null_features_and_distances_match_author <- same
      if (!same) report$rows[[name]]$status <- "INCOMPATIBLE"
    }
    if (has_pair) {
      pair <- pairs[pairs$dataset == dataset, ]
      a <- selected[selected$method_key == pair$author_bundle_key, , drop = FALSE]
      m <- selected[selected$method_key == "__dataset_metadata__", , drop = FALSE]
      .legacy_need(nrow(a) == 1L && nrow(m) == 1L && identical(a$artifact_path, pair$author_path) &&
        identical(a$md5, pair$author_md5) && identical(m$artifact_path, pair$metadata_path) &&
        identical(m$md5, pair$metadata_md5) && identical(before$path, pair$current_input_path) &&
        identical(before$md5, pair$current_input_md5), "Author/metadata pair differs from selected scope")
      report$author_hr_pairs[[dataset]] <- tryCatch(.legacy_pair(author, metadata, obs, entry, labels),
        error = function(e) list(status = "INCOMPATIBLE", reason = conditionMessage(e)))
      report$rows[[paste(dataset, "__dataset_metadata__", sep = "/")]]$status <- report$author_hr_pairs[[dataset]]$status
    }
    .legacy_need(identical(before, .legacy_input_witness(selected[1, ])), "Current input/proof changed during cohort read")
    rm(obs, labels, objects, author, metadata); invisible(gc(FALSE))
  }
  timing <- files[files$dataset == "__all__", , drop = FALSE]
  .legacy_need(nrow(timing) == 1L, "Exactly one selected merged timing file required")
  table <- .legacy_checked_read(timing$artifact_path, timing$md5, timing$size, function(p) as.data.frame(arrow::read_feather(p)))
  .legacy_fields(table, c("dataset", "method", "time_secs", "mem_GB"))
  timing_keys <- paste(table$dataset, table$method, sep = "/")
  wanted <- paste(rows$dataset[rows$method_key != "__dataset_metadata__"], rows$method_key[rows$method_key != "__dataset_metadata__"], sep = "/")
  selected_timing <- timing_keys %in% wanted
  valid_timing <- if (is.numeric(table$time_secs) && is.numeric(table$mem_GB))
    is.finite(table$time_secs) & table$time_secs >= 0 & (is.na(table$mem_GB) | (is.finite(table$mem_GB) & table$mem_GB >= 0)) else rep(FALSE, nrow(table))
  report$timing <- list(path = timing$artifact_path, md5 = timing$md5, columns = names(table),
    missing_exact_method_keys = setdiff(wanted, timing_keys), duplicate_exact_keys = unique(timing_keys[duplicated(timing_keys)]),
    selected_rows = table[selected_timing, , drop = FALSE], selected_costs_valid = all(valid_timing[selected_timing]),
    invalid_selected_rows = which(selected_timing & !valid_timing), invalid_unrelated_rows = which(!selected_timing & !valid_timing),
    unrelated_method_keys = unique(timing_keys[!selected_timing]),
    alias_policy = "no guessed aliases; missing coverage is a row-specific limitation, not a refit trigger")
  report$compatible <- report$timing$selected_costs_valid && !any(vapply(c(report$rows, report$author_hr_pairs), function(x)
    !identical(x$status, "COMPATIBLE_PENDING_REVIEW"), logical(1)))
  report$status <- if (report$compatible) "COMPATIBILITY_EVIDENCE_COMPLETE_PENDING_REVIEW" else "COMPATIBILITY_MISMATCHES"
  report
}

.legacy_sha_file <- function(path, sha) {
  .legacy_need(length(path) == 1L && grepl("^[a-f0-9]{64}$", sha) &&
    file.exists(path) && !dir.exists(path) && identical(normalizePath(path), path) &&
    identical(digest::digest(file = path, algo = "sha256"), sha), "Invalid external SHA256 file identity")
  .analysis_checked_file(path, unname(tools::md5sum(path)))
}
.legacy_json_metadata <- function(path, sha) {
  .legacy_sha_file(path, sha)
  .legacy_need(file.info(path)$size <= 16 * 1024^2, "Evidence/control must be bounded metadata")
  value <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  .legacy_sha_file(path, sha)
  value
}
legacy_analysis_audit_main <- function(args = commandArgs(TRUE)) {
  .legacy_need(length(args) %% 2L == 0L && length(args) > 0L, "Expected explicit flag/value arguments")
  opts <- setNames(as.list(args[seq.int(2L, length(args), 2L)]), args[seq.int(1L, length(args), 2L)])
  required <- c("--control", "--control-sha256", "--config", "--config-md5", "--source-identity", "--source-identity-md5",
    "--runtime-identity", "--runtime-identity-md5", "--observation", "--observation-sha256", "--zeroimp-schema",
    "--zeroimp-schema-sha256", "--zeroimp-policy", "--zeroimp-policy-sha256", "--feather-headers", "--feather-headers-sha256", "--report")
  .legacy_need(!anyDuplicated(names(opts)) && setequal(names(opts), required), "Missing or unknown explicit audit argument")
  output <- opts[["--report"]]
  .legacy_need(grepl("^/", output) && !.analysis_exists(output) && dir.exists(dirname(output)) &&
    identical(normalizePath(dirname(output)), dirname(output)), "Report must be create-only under an existing canonical directory")
  for (name in c("config", "source-identity", "runtime-identity"))
    .analysis_checked_file(opts[[paste0("--", name)]], opts[[paste0("--", name, "-md5")]])
  control <- .legacy_json_metadata(opts[["--control"]], opts[["--control-sha256"]])
  .legacy_need(identical(control$schema, "ecoda_p6_legacy_science_audit_proposal_v4") &&
    control$total_logical_rows == 209L && control$physical_files == 106L && control$author_hr_metadata_pairs == 9L &&
    control$scientifically_approved_rows == 0L, "Unexpected finite audit control")
  root <- dirname(opts[["--control"]]); frames <- list()
  for (name in c("rows.209.tsv", "files.106.tsv", "author_hr_metadata_pairs.9.tsv")) {
    path <- file.path(root, name); sha <- control$controls[[name]]
    .legacy_sha_file(path, sha)
    frames[[name]] <- utils::read.delim(path, colClasses = "character", check.names = FALSE, quote = "", comment.char = "", na.strings = NULL)
    .legacy_sha_file(path, sha)
  }
  .legacy_need(nrow(frames[[1]]) == 209L && nrow(frames[[2]]) == 106L && nrow(frames[[3]]) == 9L,
    "Finite manifest row counts differ")
  .legacy_need(identical(opts[["--observation-sha256"]], control$source_observation$sha256), "Original observation identity differs")
  observation <- .legacy_json_metadata(opts[["--observation"]], opts[["--observation-sha256"]])
  .legacy_need(identical(opts[["--feather-headers-sha256"]], control$feather_header_observation$sha256), "Feather header witness identity differs")
  headers <- .legacy_json_metadata(opts[["--feather-headers"]], opts[["--feather-headers-sha256"]])
  .legacy_header_scope(frames[[1]], headers)
  zero_schema <- .legacy_json_metadata(opts[["--zeroimp-schema"]], opts[["--zeroimp-schema-sha256"]])
  zero_policy <- .legacy_json_metadata(opts[["--zeroimp-policy"]], opts[["--zeroimp-policy-sha256"]])
  .legacy_need(identical(zero_policy$source_schema_sha256, opts[["--zeroimp-schema-sha256"]]) &&
    length(zero_schema$rows) == 9L && length(zero_policy$datasets) == 9L &&
    zero_policy$scientifically_approved_rows == 0L &&
    setequal(vapply(zero_schema$rows, `[[`, character(1), "dataset"), frames[[3]]$dataset) &&
    setequal(vapply(zero_policy$datasets, `[[`, character(1), "dataset"), frames[[3]]$dataset),
    "Zero-imputation evidence identities/scope differ")
  datasets <- read_datasets_json(opts[["--config"]], view = "benchmark_analysis")
  module_dir <- normalizePath(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "src/utils/py"), mustWork = TRUE)
  result <- audit_legacy_analysis_inputs(frames[[1]], frames[[2]], frames[[3]], datasets,
    function(path, columns) .legacy_read_obs(path, columns, module_dir))
  result$validation_identities <- opts[setdiff(names(opts), "--report")]
  result$original_observation_status <- observation$status
  result$zero_imputation <- list(schema_evidence = opts[["--zeroimp-schema"]], schema_sha256 = opts[["--zeroimp-schema-sha256"]],
    policy_evidence = opts[["--zeroimp-policy"]], policy_sha256 = opts[["--zeroimp-policy-sha256"]],
    datasets = vapply(zero_policy$datasets, `[[`, character(1), "dataset"),
    status = "HISTORY_AND_SCHEMA_ONLY; CURRENT_COUNTS_PAIR_AND_ROOT_REVIEW_REQUIRED; NO_PAYLOAD_LABEL_DISTANCE_OR_POLICY_WITNESS",
    science_approved = FALSE, payloads_read_by_this_audit = 0L, lisi_inverted = FALSE)
  for (name in c("config", "source-identity", "runtime-identity"))
    .analysis_checked_file(opts[[paste0("--", name)]], opts[[paste0("--", name, "-md5")]])
  .legacy_sha_file(opts[["--control"]], opts[["--control-sha256"]])
  tmp <- tempfile(".legacy-audit-", tmpdir = dirname(output)); on.exit(unlink(tmp))
  jsonlite::write_json(result, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  .legacy_need(!.analysis_exists(output) && file.link(tmp, output), "Cannot atomically publish create-only audit report")
  message("LEGACY_ANALYSIS_COMPATIBILITY=", result$status)
  invisible(result)
}

if (sys.nframe() == 0L) {
  project <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
  for (path in c("analysis_release.R", "evaluation_scores.R", "datasets_io.R", "seurat_utils.R"))
    source(file.path(project, "src/utils", path))
  suppressPackageStartupMessages(library(dplyr))
  result <- legacy_analysis_audit_main()
  if (!result$compatible) quit(status = 1L)
}
