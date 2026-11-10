# An explicit, finite selection of reviewed artifacts; no directory discovery or
# fallback to another generation. Receipts are evidence, not producer records.
.analysis_key <- function(x) {
  do.call(paste, c(x[c("dataset", "view", "method_key")], sep = "\r"))
}

.analysis_fields <- function(x, fields, label) {
  if (!is.data.frame(x) || !all(fields %in% names(x)) || !nrow(x))
    stop(label, " requires nonempty rows and fields: ", paste(fields, collapse = ", "))
  for (field in fields) {
    if (anyNA(x[[field]]) || any(!nzchar(trimws(as.character(x[[field]])))))
      stop(label, " has missing ", field)
  }
}

.analysis_checked_file <- function(path, md5) {
  if (length(path) != 1L || !grepl("^/", path) || !file.exists(path) ||
      dir.exists(path) || !is.finite(file.info(path)$size) || file.info(path)$size <= 0 ||
      !grepl("^[a-f0-9]{32}$", md5) ||
      !identical(unname(tools::md5sum(path)), md5))
    stop("Missing, empty or checksum-invalid analysis artifact: ", path)
  canonical <- normalizePath(path, mustWork = TRUE)
  if (!identical(canonical, path)) stop("Analysis paths must be canonical, without aliases or symlinks: ", path)
  canonical
}

validate_analysis_release <- function(ledger, expected_scope, require_projection = FALSE) {
  if (!is.logical(require_projection) || length(require_projection) != 1L || is.na(require_projection))
    stop("require_projection must be an explicit logical scalar")
  keys <- c("dataset", "view", "method_key")
  .analysis_fields(ledger, c(keys, "kind", "artifact_path", "run_root", "md5",
                            "receipt_path", "receipt_md5"), "release ledger")
  .analysis_fields(expected_scope, c(keys, "current_input_md5", "source_policy", "status"), "expected scope")
  if (!"require_new" %in% names(expected_scope) ||
      !is.logical(expected_scope$require_new) || anyNA(expected_scope$require_new))
    stop("expected scope requires explicit logical require_new")
  if (!"bundle_key" %in% names(ledger) || anyNA(ledger$bundle_key))
    stop("ledger requires explicit bundle_key (empty for standalone/Feather)")
  if (anyDuplicated(.analysis_key(ledger)) || anyDuplicated(.analysis_key(expected_scope)) ||
      !setequal(.analysis_key(ledger), .analysis_key(expected_scope)))
    stop("Release keys must exactly equal the unique expected scope")
  ledger <- ledger[match(.analysis_key(expected_scope), .analysis_key(ledger)), , drop = FALSE]
  allowed <- c("rds_bundle", "standalone_result", "distance_feather", "embedding_feather", "metadata_feather")
  if (any(!ledger$kind %in% allowed)) stop("Unknown release artifact kind")
  if (any((ledger$kind == "rds_bundle") != nzchar(ledger$bundle_key)))
    stop("Only aggregate RDS rows require a bundle_key")
  if (any(!expected_scope$status %in% c("approved", "blocked"))) stop("Unknown expected-scope review status")
  if (any(expected_scope$status == "blocked")) {
    if (!"blocked_reason" %in% names(expected_scope)) stop("Blocked scope requires blocked_reason")
    stop("Blocked release scope: ", paste(expected_scope$blocked_reason[expected_scope$status == "blocked"], collapse = "; "))
  }
  if (any((ledger$kind == "metadata_feather") != (ledger$method_key == "__sample_metadata__")))
    stop("Metadata Feather requires the reserved __sample_metadata__ key")
  receipts <- vector("list", nrow(ledger))
  for (i in seq_len(nrow(ledger))) {
    row <- ledger[i, , drop = FALSE]
    if (!grepl("^/", row$run_root)) stop("run_root must be absolute")
    root <- normalizePath(row$run_root, mustWork = TRUE)
    if (!identical(root, row$run_root)) stop("run_root must be canonical, without aliases or symlinks")
    if (!dir.exists(root) || !startsWith(row$artifact_path, paste0(root, "/")))
      stop("Artifact is outside its explicit run_root: ", row$artifact_path)
    receipt_path <- .analysis_checked_file(row$receipt_path, row$receipt_md5)
    if (file.info(receipt_path)$size > 16 * 1024^2) stop("Release receipt must be small metadata")
    receipt <- readRDS(receipt_path)
    .analysis_checked_file(row$receipt_path, row$receipt_md5)
    if (!is.list(receipt) || !identical(receipt$schema, "ecoda_analysis_release_receipt_v1") ||
        !identical(receipt$status, "reviewed") ||
        !receipt$origin %in% c("new", "legacy_reuse") ||
        !is.character(receipt$evidence) || length(receipt$evidence) != 1L ||
        !nzchar(trimws(receipt$evidence))) stop("Invalid reviewed release receipt: ", receipt_path)
    .analysis_checked_file(receipt$evidence_path, receipt$evidence_md5)
    if (!identical(receipt$source_policy, as.character(expected_scope$source_policy[i])))
      stop("Receipt does not establish the approved source policy")
    for (field in c(keys, "kind", "artifact_path", "bundle_key", "run_root", "md5")) {
      if (!identical(receipt[[field]], as.character(row[[field]])))
        stop("Receipt disagrees with ledger field ", field, ": ", receipt_path)
    }
    if (!identical(receipt$input_md5, as.character(expected_scope$current_input_md5[i])) ||
        !grepl("^[a-f0-9]{32}$", receipt$input_md5))
      stop("Receipt does not establish the approved current input identity")
    if (expected_scope$require_new[i] && !identical(receipt$origin, "new"))
      stop("Required new result cannot reuse a legacy artifact: ", row$method_key)
    # This check precedes even opening the selected RDS. A scientific reuse
    # receipt alone does not authorize transferring an RNA aggregate locally.
    compact_screen <- identical(receipt$metadata_schema, "ecoda_cell_depth_five_v1") &&
      identical(as.character(row$kind), "standalone_result") &&
      identical(as.character(row$method_key), "ECODA_authors_HR_cell_subsetting")
    if ((!is.null(receipt$projection_schema)) ||
        (require_projection && !compact_screen && !row$kind %in% c("distance_feather", "metadata_feather")))
      .analysis_projection_proof(row, receipt)
    .analysis_checked_file(row$artifact_path, row$md5)
    receipts[[i]] <- receipt
  }
  attr(ledger, "receipts") <- receipts
  ledger
}

load_analysis_release <- function(ledger, expected_scope, require_projection = FALSE) {
  ledger <- validate_analysis_release(ledger, expected_scope, require_projection)
  receipts <- attr(ledger, "receipts")
  results <- vector("list", nrow(ledger))
  names(results) <- .analysis_key(ledger)
  for (i in seq_len(nrow(ledger))) {
    row <- ledger[i, , drop = FALSE]
    # Recheck at the read boundary, including when one aggregate supplies many keys.
    .analysis_checked_file(row$artifact_path, row$md5)
    if (row$kind %in% c("rds_bundle", "standalone_result")) {
      result <- readRDS(row$artifact_path)
      if (row$kind == "rds_bundle") {
        if (!is.list(result) || is.null(names(result)) || anyDuplicated(names(result)) ||
            !row$bundle_key %in% names(result)) stop("Missing explicit aggregate bundle key")
        result <- result[[row$bundle_key]]
      }
      if (!is.list(result)) stop("Selected result must be a bundle list")
    } else {
      if (!requireNamespace("arrow", quietly = TRUE)) stop("Feather release needs installed arrow")
      if (!"sample_id_col" %in% names(row) || is.na(row$sample_id_col) || !nzchar(row$sample_id_col))
        stop("Feather rows require explicit sample_id_col")
      frame <- as.data.frame(arrow::read_feather(row$artifact_path))
      if (!row$sample_id_col %in% names(frame)) stop("Feather lacks declared sample identity")
      ids <- as.character(frame[[row$sample_id_col]])
      if (row$kind == "metadata_feather") {
        if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) || !length(ids))
          stop("Metadata Feather requires unique nonmissing sample IDs")
        results[[i]] <- list(result = list(metadata = frame), source = row, receipt = receipts[[i]],
          identity = list(dataset = row$dataset, view = row$view, method_key = row$method_key,
            artifact_md5 = row$md5, receipt_md5 = row$receipt_md5,
            input_md5 = receipts[[i]]$input_md5, source_policy = receipts[[i]]$source_policy))
        .analysis_checked_file(row$artifact_path, row$md5)
        next
      }
      frame[[row$sample_id_col]] <- NULL
      if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
          !nrow(frame) || !ncol(frame) || !all(vapply(frame, is.numeric, logical(1))))
        stop("Malformed Feather samples or numeric data")
      values <- as.matrix(frame)
      if (any(!is.finite(values))) stop("Nonfinite Feather data")
      rownames(values) <- ids
      if (row$kind == "distance_feather") {
        if (!identical(colnames(values), ids)) stop("Feather distances require identical named axes")
        result <- list(dist_mat = values)
      } else result <- list(feat_mat = values, dist_mat = stats::dist(values))
    }
    .analysis_checked_file(row$artifact_path, row$md5)
    if (!is.null(receipts[[i]]$projection_schema)) {
      proof <- .analysis_projection_proof(row, receipts[[i]])
      if (!identical(.analysis_fingerprint(result$dist_mat), proof$distance_md5) ||
          !identical(.analysis_fingerprint(result$labels), proof$labels_md5))
        stop("Projected distances or stored labels differ from the parent witness")
    }
    results[[i]] <- list(result = result, source = row, receipt = receipts[[i]],
      identity = list(dataset = row$dataset, view = row$view, method_key = row$method_key,
                      artifact_md5 = row$md5, receipt_md5 = row$receipt_md5,
                      input_md5 = receipts[[i]]$input_md5, source_policy = receipts[[i]]$source_policy))
  }
  results
}

.analysis_fingerprint <- function(value) {
  path <- tempfile("analysis-fingerprint-")
  on.exit(unlink(path))
  # Version 2 materializes ALTREP vectors, so a Feather->RDS round trip does
  # not change the witness solely through an internal vector representation.
  saveRDS(value, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

# The producer's counts/source checks occur before projection. This compact
# contract retains every recipe outcome and only sample-level PASS geometry.
.analysis_zeroimp_screen <- function(value) {
  metrics <- evaluation_metrics()
  counts <- c("0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")
  percent <- c("0.001", "0.01", "0.1", "1", "2", "5")
  menu <- c(unlist(lapply(counts, function(x) paste0(c("counts_zeros_", "counts_all_"), x)), use.names = FALSE),
    unlist(lapply(percent, function(x) paste0(c("percentage_all_", "percentage_zeros_", "multLN_", "multRepl_"), x, "%")), use.names = FALSE), "asinsqrt")
  need <- function(ok, message) if (!isTRUE(ok)) stop(message)
  need(identical(value$schema, "ecoda_zeroimp_counts_five_v1") &&
    identical(value$policy$schema, value$schema) && identical(value$policy$version, 1L) &&
    identical(value$policy$recipe_menu, menu) && identical(value$policy$metrics, metrics) &&
    identical(value$policy$scoring_policy, evaluation_scoring_policy()) &&
    identical(value$policy$parent_counts_offset, 0.5) && identical(value$policy$seed, 123L),
    "Unknown counts-only zero-imputation policy")
  # These are the original producer's fingerprints, not hashes of a new RDS
  # serialization: the checked A receipt binds the entire unchanged package.
  hash <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{32}$", x)
  need(all(vapply(value[c("policy_md5", "labels_md5", "counts_md5", "recipe_body_md5",
    "config_md5", "runner_md5", "release_md5")], hash, logical(1))), "Zero-imputation producer identities missing")
  ids <- unname(names(value$labels))
  need(length(ids) > 1L && !anyNA(ids) && !anyDuplicated(ids) && all(nzchar(ids)) &&
    !anyNA(value$labels) && all(nzchar(as.character(value$labels))), "Invalid zero-imputation sample labels")
  if (!is.null(value[["counts"]])) need(is.matrix(value[["counts"]]) && is.numeric(value[["counts"]]) &&
    all(is.finite(value[["counts"]])) && all(value[["counts"]] >= 0) && all(value[["counts"]] == round(value[["counts"]])) &&
    identical(unname(rownames(value[["counts"]])), ids), "Zero-imputation parent count identity differs")
  else need(identical(value$projection_schema, "ecoda_zeroimp_counts_projection_v1"),
    "Unprojected zero-imputation package lacks parent counts")
  table <- value$results
  need(is.data.frame(table) && identical(table$setting, menu) && identical(names(value$recipes), menu) &&
    all(c("status", "reason", "scoring_sample_ids", "dropped_sample_ids", metrics) %in% names(table)) &&
    !anyNA(table$status) && all(table$status %in% c("PASS", "FAILED_OR_DROPPED_BY_ORDINARY_PRODUCER")),
    "Zero-imputation screen requires all 43 explicit recipe outcomes")
  for (i in seq_along(menu)) {
    record <- value$recipes[[i]]; row <- table[i, , drop = FALSE]
    need(identical(record$setting, menu[i]) && identical(record$status, row$status) &&
      identical(record$reason, row$reason) &&
      identical(row$scoring_sample_ids, as.character(jsonlite::toJSON(record$sample_ids, auto_unbox = FALSE, null = "null"))) &&
      identical(row$dropped_sample_ids, as.character(jsonlite::toJSON(record$dropped_ids, auto_unbox = FALSE, null = "null"))),
      "Zero-imputation recipe outcome evidence differs")
    if (row$status == "PASS") {
      d <- as.matrix(record$dist_mat)
      need(inherits(record$dist_mat, "dist") && identical(record$sample_ids, ids) &&
        !length(record$dropped_ids) && !length(record$extra_ids) && identical(record$labels, value$labels) &&
        identical(unname(rownames(d)), ids) && identical(unname(colnames(d)), ids) &&
        all(is.finite(d)) && all(d >= -1e-10) && max(abs(d - t(d))) <= 1e-10 && max(abs(diag(d))) <= 1e-10 &&
        hash(record$distance_md5), "PASS zero-imputation geometry/sample identity differs")
      need(is.list(record$scores) && !anyDuplicated(names(record$scores)) &&
        all(vapply(metrics, function(m) is.numeric(record$scores[[m]]) && length(record$scores[[m]]) == 1L &&
          is.finite(record$scores[[m]]) && identical(record$scores[[m]], row[[m]]), logical(1))) &&
        row$lisi_score >= 0 && row$lisi_score <= 1, "PASS zero-imputation canonical five differ")
    } else need(is.null(record$scores) && is.null(record$dist_mat) && is.null(record$labels) &&
      all(is.na(row[metrics])) && length(record$reason) == 1L && !is.na(record$reason) && nzchar(record$reason),
      "Failed zero-imputation recipe cannot carry invented scores/geometry")
  }
  need(identical(value$recipes$counts_all_0.5$status, "PASS"), "Zero-imputation reference must pass")
  invisible(value)
}

.analysis_projection_proof <- function(row, receipt) {
  if (!identical(receipt$projection_schema, "ecoda_analysis_projection_v1"))
    stop("Local result intake requires a reviewed metadata projection before reading RDS: ", row$method_key)
  path <- .analysis_checked_file(receipt$evidence_path, receipt$evidence_md5)
  if (file.info(path)$size > 16 * 1024^2) stop("Projection witness must be small metadata")
  proof <- readRDS(path)
  .analysis_checked_file(path, receipt$evidence_md5)
  if (!is.list(proof) || !identical(proof$schema, receipt$projection_schema) ||
      !identical(proof$dataset, as.character(row$dataset)) ||
      !identical(proof$view, as.character(row$view)) ||
      !identical(proof$method_key, as.character(row$method_key)) ||
      !identical(proof$output_md5, as.character(row$md5)) ||
      !identical(proof$bundle_key, as.character(row$bundle_key)) ||
      !identical(proof$origin, receipt$origin)) stop("Invalid projection identity or scientific origin")
  for (field in c("parent_receipt_copy", "parent_evidence_copy")) {
    relative <- proof[[field]]
    if (!is.character(relative) || length(relative) != 1L ||
        !grepl("^parents/[0-9]+\\.(receipt|evidence)$", relative)) stop("Invalid copied parent proof path")
  }
  parent <- readRDS(.analysis_checked_file(file.path(row$run_root, proof$parent_receipt_copy), proof$parent$receipt_md5))
  .analysis_checked_file(file.path(row$run_root, proof$parent_evidence_copy), parent$evidence_md5)
  fields <- c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")
  if (!is.list(parent) || !identical(parent$schema, "ecoda_analysis_release_receipt_v1") ||
      !identical(parent$status, "reviewed") || !identical(parent$origin, receipt$origin) ||
      !identical(parent$input_md5, receipt$input_md5) || !identical(parent$source_policy, receipt$source_policy) ||
      any(!vapply(fields, function(field) identical(parent[[field]], proof$parent[[field]]), logical(1))))
    stop("Projection cannot change the parent scientific review, origin or input policy")
  for (field in c("source", "runtime", "control"))
    .analysis_checked_file(file.path(row$run_root, "parents", paste0("validation.", field)), proof$validation_md5[[field]])
  control_path <- file.path(row$run_root, "parents", "validation.control")
  if (file.info(control_path)$size > 16 * 1024^2) stop("Frozen parent control must be small metadata")
  control <- readRDS(control_path)
  key <- .analysis_key(row)
  original <- control$ledger[.analysis_key(control$ledger) == key, , drop = FALSE]
  expected <- control$expected_scope[.analysis_key(control$expected_scope) == key, , drop = FALSE]
  if (nrow(original) != 1L || nrow(expected) != 1L || expected$status != "approved" ||
      !identical(expected$current_input_md5, receipt$input_md5) ||
      !identical(expected$source_policy, receipt$source_policy) ||
      any(!vapply(c(fields, "receipt_path", "receipt_md5"), function(field)
        identical(as.character(original[[field]]), proof$parent[[field]]), logical(1))))
    stop("Projection parent must match the original approved control row, including its scientific receipt")
  proof
}

# Only these sample-by-cell-type matrices feed maintained numbered PCA/ratio
# panels. Other composition methods retain scalar/category metadata; RNA
# features and pseudobulk/cache objects are never included in a projection.
analysis_projection_ct_methods <- function() c("ECODA_authors_HR", "ECODA_authors_LR",
  "ECODA_authors_HR_NULL", "Freq_highres", "ECODA_HiTME_HR_layer2",
  "ECODA_scATOMIC_HR", "ECODA_seuratres_0.4", "ECODA_authors_HR_top_varexp0")

.analysis_thin_result <- function(item) {
  value <- item$result; key <- as.character(item$source$method_key)
  if (key == "__dataset_metadata__") {
    fields <- c("labels", "n_samples", "n_cells", "cells_per_sample", "n_cell_types_high_res")
    if (!all(fields %in% names(value))) stop("Incomplete dataset metadata")
    return(value[fields])
  }
  if (key == "__trans__") {
    if (!exists("validate_transformation_result", mode = "function")) stop("Projection requires the pure transformation validator")
    validate_transformation_result(value)
    fields <- c("trans_method", unname(transformation_analysis_policy()$score_columns))
    thin <- as.data.frame(lapply(fields, function(field) value[[field]]), stringsAsFactors = FALSE)
    names(thin) <- fields
    for (field in c("transformation_policy", "transformation_replicates", "composition_counts", "sample_labels"))
      attr(thin, field) <- attr(value, field)
    validate_transformation_result(thin)
    return(thin)
  }
  if (key == "__zeroimp__") {
    if (!is.null(value$schema)) {
      .analysis_zeroimp_screen(value)
      if (!identical(value$dataset, as.character(item$source$dataset))) stop("Counts-only screen cohort differs from its ledger row")
      if (!identical(item$receipt$origin, "new")) stop("Counts-only screens require their new reviewed producer")
      fields <- c("schema", "dataset", "run_id", "results", "recipes", "labels", "labels_md5", "counts_md5",
        "policy", "policy_md5", "source", "config_md5", "runner_md5", "release_md5", "recipe_body_md5",
        "rng_before", "rng_after", "warnings", "csv_md5")
      thin <- value[intersect(fields, names(value))]
      thin$projection_schema <- "ecoda_zeroimp_counts_projection_v1"
      .analysis_zeroimp_screen(thin)
      return(thin)
    }
    clean <- function(x, depth = 0L) {
      if (is.list(x) && length(x) && !is.null(names(x)) &&
          all(names(x) %in% c(evaluation_metrics(), "mod_knnsqrtn_score", "mod_knn6_score", "mod_knn9_score"))) {
        if (anyDuplicated(names(x)) || !all(evaluation_metrics() %in% names(x)) ||
            any(!vapply(x, function(v) is.numeric(v) && length(v) == 1L, logical(1))) ||
            any(!vapply(x[evaluation_metrics()], is.finite, logical(1))))
          stop("Invalid zero-imputation scalar scores")
        return(x)
      }
      if (!is.list(x) || !length(x) || depth >= 2L || is.null(names(x)) || anyDuplicated(names(x)))
        stop("Zero-imputation projection accepts only named scalar-score tables")
      lapply(x, clean, depth = depth + 1L)
    }
    return(clean(value))
  }
  distance <- as.matrix(value$dist_mat)
  # Incidental names on an ID vector are not sample identity; preserve value objects.
  ids <- unname(rownames(distance))
  if (!is.numeric(distance) || !length(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
      !identical(unname(colnames(distance)), ids) || any(!is.finite(distance)) ||
      any(distance < 0) || any(abs(distance - t(distance)) > 1e-10) || any(abs(diag(distance)) > 1e-10))
    stop("Projection requires finite named symmetric distances")
  if (grepl("(^|_)NULL($|_)", key) && is.null(value$labels)) stop("NULL projections require actual stored labels")
  if (!is.null(value$labels) && (is.null(names(value$labels)) || !identical(unname(names(value$labels)), ids) ||
      anyNA(value$labels) || any(!nzchar(trimws(as.character(value$labels))))))
    stop("Stored labels must retain exact distance sample order")
  fields <- c("dist_mat", "labels", "scores", "exec_time", "time_secs", "mem_GB", "timing_schema", "timing_id",
    "variant_time_secs", "variant_mem_GB", "shared_time_secs", "shared_mem_GB", "shared_timing_id",
    "aggregate_time_secs", "shared_fit_time_secs", "allowed_dropped_sample_ids")
  fields <- intersect(fields, names(value))
  thin <- setNames(lapply(fields, function(field) value[[field]]), fields)
  thin$n_features <- if (is.null(value$feat_mat)) {
    if (is.null(value$n_features)) NA_integer_ else value$n_features
  } else ncol(value$feat_mat)
  if (!is.null(thin$n_features) && (!is.numeric(thin$n_features) || length(thin$n_features) != 1L ||
      (!is.na(thin$n_features) && (!is.finite(thin$n_features) || thin$n_features < 1 || thin$n_features != as.integer(thin$n_features)))))
    stop("Original feature count must be a positive scalar")
  if (!is.null(value$scores) && (!is.list(value$scores) || anyDuplicated(names(value$scores)) ||
      any(!vapply(value$scores, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1)))))
    stop("Projection scores must be finite scalars")
  # Count annotation categories without retaining every annotation matrix.
  if (grepl("^ECODA_(authors|HiTME|scATOMIC|seuratres|deconv)", key) || key == "Freq_highres") {
    categories <- if (key == "ECODA_deconv") colnames(value$feat_mat) else colnames(value$counts)
    if (!is.null(categories)) thin$cell_type_names <- categories
  }
  if (item$source$view == "benchmark_analysis" && key %in% analysis_projection_ct_methods()) {
    if (is.null(value$feat_mat) || !identical(unname(rownames(value$feat_mat)), ids) ||
        any(!is.finite(as.matrix(value$feat_mat)))) stop("Invalid whitelisted CT feature matrix")
    thin$feat_mat <- value$feat_mat
    if (!is.null(value$counts)) {
      if (!setequal(rownames(value$counts), ids) || any(!is.finite(as.matrix(value$counts)))) stop("Invalid CT counts")
      thin$counts <- value$counts
    }
  }
  thin
}

.analysis_save_new <- function(value, path) {
  if (.analysis_exists(path)) stop("Projection output already exists: ", path)
  tmp <- tempfile(".projection-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(value, tmp, version = 3)
  if (!file.link(tmp, path)) stop("Cannot atomically publish new projection: ", path)
  invisible(path)
}

.analysis_exists <- function(path) file.exists(path) ||
  (!is.na(Sys.readlink(path)) && nzchar(Sys.readlink(path)))

.analysis_projection_paths <- function(control, root) {
  ledger <- control$ledger
  c(ledger$artifact_path, ledger$receipt_path,
    file.path(root, "receipts", paste0(seq_len(nrow(ledger)), ".projection.rds")),
    file.path(root, "parents", paste0(seq_len(nrow(ledger)), ".receipt")),
    file.path(root, "parents", paste0(seq_len(nrow(ledger)), ".evidence")),
    file.path(root, "parents", paste0("validation.", c("source", "runtime", "control"))),
    file.path(root, "analysis_release.rds"),
    if (!is.null(control$execution_times)) control$execution_times$artifact_path)
}

# The producer's raw copied root remains immutable. Only fresh local access
# receipts change paths; no scientific parent identity or thin payload changes.
mirror_analysis_projection <- function(copied_root, output_root, control_md5) {
  if (!identical(normalizePath(copied_root, mustWork = TRUE), copied_root) ||
      !grepl("^/", output_root) || .analysis_exists(output_root)) stop("Explicit canonical copied root and fresh local root required")
  path <- .analysis_checked_file(file.path(copied_root, "analysis_release.rds"), control_md5)
  if (file.info(path)$size > 16 * 1024^2) stop("Projection control must be small metadata")
  control <- readRDS(path)
  roots <- unique(control$ledger$run_root)
  if (length(roots) != 1L || !grepl("^/", roots)) stop("Expected one explicit producer projection root")
  original_root <- roots[[1L]]
  rebase <- function(paths, target) {
    if (anyNA(paths) || any(!startsWith(paths, paste0(original_root, "/")))) stop("Projection path escapes its declared producer root")
    file.path(target, substring(paths, nchar(original_root) + 2L))
  }
  inventory <- utils::read.delim(file.path(copied_root, "files.tsv"), header = FALSE,
    col.names = c("path", "md5"), colClasses = "character", quote = "", comment.char = "")
  expected_paths <- .analysis_projection_paths(control, original_root)
  if (!nrow(inventory) || anyDuplicated(inventory$path) ||
      !setequal(inventory$path, expected_paths)) stop("Copied projection inventory differs from its finite control")
  source_paths <- rebase(inventory$path, copied_root)
  for (i in seq_along(source_paths)) .analysis_checked_file(source_paths[i], inventory$md5[i])
  # Bind each copied witness to the immutable parent control before copying
  # anything, even when the original HPC paths are unavailable locally.
  for (i in seq_len(nrow(control$ledger))) {
    row <- control$ledger[i, , drop = FALSE]
    receipt <- readRDS(.analysis_checked_file(rebase(row$receipt_path, copied_root), row$receipt_md5))
    row$run_root <- copied_root; row$artifact_path <- rebase(row$artifact_path, copied_root)
    receipt$evidence_path <- rebase(receipt$evidence_path, copied_root)
    .analysis_projection_proof(row, receipt)
  }
  if (!dir.create(output_root) || !identical(normalizePath(output_root), output_root)) stop("Cannot create canonical local projection root")
  dir.create(file.path(output_root, "parents")); dir.create(file.path(output_root, "receipts"))
  destination_paths <- rebase(inventory$path, output_root)
  # Receipts and the control are written afresh below, never overwritten.
  access_paths <- c(control$ledger$receipt_path, file.path(original_root, "analysis_release.rds"))
  for (i in which(!inventory$path %in% access_paths))
    if (!file.copy(source_paths[i], destination_paths[i], overwrite = FALSE)) stop("Cannot copy finite projection metadata")
  for (i in seq_len(nrow(control$ledger))) {
    row <- control$ledger[i, , drop = FALSE]
    receipt <- readRDS(rebase(row$receipt_path, copied_root))
    control$ledger$artifact_path[i] <- rebase(row$artifact_path, output_root)
    control$ledger$receipt_path[i] <- rebase(row$receipt_path, output_root)
    control$ledger$run_root[i] <- output_root
    receipt$artifact_path <- control$ledger$artifact_path[i]; receipt$run_root <- output_root
    receipt$evidence_path <- rebase(receipt$evidence_path, output_root)
    .analysis_save_new(receipt, control$ledger$receipt_path[i])
    control$ledger$receipt_md5[i] <- unname(tools::md5sum(control$ledger$receipt_path[i]))
  }
  if (!is.null(control$execution_times)) control$execution_times$artifact_path <- rebase(control$execution_times$artifact_path, output_root)
  validate_analysis_release(control$ledger, control$expected_scope, require_projection = TRUE)
  .analysis_save_new(control, file.path(output_root, "analysis_release.rds"))
  local_inventory <- data.frame(path = destination_paths, md5 = unname(tools::md5sum(destination_paths)))
  tmp <- tempfile(".projection-inventory-", tmpdir = output_root); on.exit(unlink(tmp))
  utils::write.table(local_inventory, tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (!file.link(tmp, file.path(output_root, "files.tsv"))) stop("Cannot publish local access inventory")
  for (i in seq_along(source_paths)) .analysis_checked_file(source_paths[i], inventory$md5[i])
  invisible(control)
}

analysis_release_input_inventory <- function(control) {
  ledger <- validate_analysis_release(control$ledger, control$expected_scope)
  receipts <- attr(ledger, "receipts")
  inventory <- data.frame(path = c(ledger$artifact_path, ledger$receipt_path,
    vapply(receipts, `[[`, character(1), "evidence_path")),
    md5 = c(ledger$md5, ledger$receipt_md5, vapply(receipts, `[[`, character(1), "evidence_md5")))
  if (!is.null(control$execution_times)) {
    timing <- control$execution_times
    .analysis_checked_file(timing$artifact_path, timing$md5)
    inventory <- rbind(inventory, data.frame(path = timing$artifact_path, md5 = timing$md5))
  }
  inventory <- unique(inventory)
  if (anyDuplicated(inventory$path) || any(grepl("[\t\r\n]", inventory$path))) stop("Conflicting or malformed explicit input paths")
  inventory
}

project_analysis_release <- function(control, output_root, validation_files) {
  if (!is.list(control) || !all(c("ledger", "expected_scope") %in% names(control))) stop("Explicit approved A-control required")
  allowed_control <- c("ledger", "expected_scope", "reference_method_key", "execution_times",
    "schema", "status", "release_id", "datasets", "view", "expected_grid", "control_md5", "evaluation_axes")
  if (any(!names(control) %in% allowed_control) || anyDuplicated(names(control)))
    stop("Projection accepts a finite metadata control, never a loaded release/result object")
  for (frame in list(control$ledger, control$expected_scope))
    if (any(!vapply(frame, function(column) is.atomic(column) && is.null(dim(column)), logical(1))))
      stop("Control ledger/scope columns must be metadata vectors")
  parent <- validate_analysis_release(control$ledger, control$expected_scope)
  if (!identical(names(validation_files), c("source", "runtime", "control"))) stop("Explicit validation source/runtime/control files required")
  validation_md5 <- vapply(validation_files, function(path) {
    hash <- unname(tools::md5sum(path)); .analysis_checked_file(path, hash); hash
  }, character(1))
  if (file.info(validation_files[["control"]])$size > 16 * 1024^2) stop("Parent control must be small metadata")
  if (!identical(readRDS(validation_files[["control"]]), control)) stop("Approved parent control differs from frozen control")
  if (.analysis_exists(output_root) || !grepl("^/", output_root)) stop("Projection root must be new and absolute")
  if (!dir.create(output_root) || !identical(normalizePath(output_root), output_root)) stop("Projection root is not canonical")
  dir.create(file.path(output_root, "parents")); dir.create(file.path(output_root, "receipts"))
  for (field in names(validation_files))
    if (!file.copy(validation_files[[field]], file.path(output_root, "parents", paste0("validation.", field)), overwrite = FALSE))
      stop("Cannot freeze validation identity")
  ledger <- parent; attr(ledger, "receipts") <- NULL
  receipts <- attr(parent, "receipts")
  for (dataset in unique(parent$dataset)) {
    selected <- which(parent$dataset == dataset)
    # Whole RNA-bearing aggregates are deserialized only on HPC, one cohort
    # at a time. No original aggregate or cache is copied into this root.
    loaded <- load_analysis_release(parent[selected, , drop = FALSE],
      control$expected_scope[match(.analysis_key(parent[selected, , drop = FALSE]), .analysis_key(control$expected_scope)), , drop = FALSE])
    for (j in seq_along(selected)) {
      i <- selected[j]; item <- loaded[[j]]; original <- receipts[[i]]
      parent_row <- parent[i, , drop = FALSE]
      suffix <- if (parent_row$kind == "metadata_feather") ".feather" else ".rds"
      path <- file.path(output_root, paste0(i, suffix))
      thin <- if (parent_row$kind == "metadata_feather") item$result else .analysis_thin_result(item)
      if (parent_row$kind == "metadata_feather") {
        if (!file.copy(parent_row$artifact_path, path, overwrite = FALSE)) stop("Cannot copy selected sample metadata")
        ledger$bundle_key[i] <- ""
      } else {
        .analysis_save_new(setNames(list(thin), parent_row$method_key), path)
        ledger$kind[i] <- "rds_bundle"; ledger$bundle_key[i] <- parent_row$method_key
      }
      ledger$artifact_path[i] <- path; ledger$run_root[i] <- output_root
      ledger$md5[i] <- unname(tools::md5sum(path))
      parent_receipt_copy <- paste0("parents/", i, ".receipt")
      parent_evidence_copy <- paste0("parents/", i, ".evidence")
      for (pair in list(c(parent_row$receipt_path, parent_receipt_copy), c(original$evidence_path, parent_evidence_copy))) {
        if (file.info(pair[1])$size > 16 * 1024^2 || !file.copy(pair[1], file.path(output_root, pair[2]), overwrite = FALSE))
          stop("Parent review evidence must be bounded metadata")
      }
      proof <- list(schema = "ecoda_analysis_projection_v1", dataset = parent_row$dataset, view = parent_row$view,
        method_key = parent_row$method_key, origin = original$origin, parent = as.list(parent_row),
        parent_receipt_copy = parent_receipt_copy, parent_evidence_copy = parent_evidence_copy,
        output_md5 = ledger$md5[i], bundle_key = ledger$bundle_key[i],
        distance_md5 = .analysis_fingerprint(item$result$dist_mat), labels_md5 = .analysis_fingerprint(item$result$labels),
        original_feature_count = thin$n_features, retained_fields = names(thin),
        validation_md5 = validation_md5)
      proof_path <- file.path(output_root, "receipts", paste0(i, ".projection.rds"))
      .analysis_save_new(proof, proof_path)
      receipt <- original
      for (field in c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5"))
        receipt[[field]] <- as.character(ledger[[field]][i])
      receipt$projection_schema <- proof$schema
      receipt$evidence <- "Metadata-only projection; parent scientific origin unchanged; not a new fit"
      receipt$evidence_path <- proof_path; receipt$evidence_md5 <- unname(tools::md5sum(proof_path))
      ledger$receipt_path[i] <- file.path(output_root, "receipts", paste0(i, ".receipt.rds"))
      .analysis_save_new(receipt, ledger$receipt_path[i]); ledger$receipt_md5[i] <- unname(tools::md5sum(ledger$receipt_path[i]))
      .analysis_checked_file(parent_row$artifact_path, parent_row$md5)
      .analysis_checked_file(parent_row$receipt_path, parent_row$receipt_md5)
      .analysis_checked_file(original$evidence_path, original$evidence_md5)
    }
    rm(loaded, item, thin)
    invisible(gc(FALSE))
  }
  control$ledger <- ledger
  if (!is.null(control$execution_times)) {
    timing <- control$execution_times
    .analysis_checked_file(timing$artifact_path, timing$md5)
    timing_frame <- as.data.frame(arrow::read_feather(timing$artifact_path))
    if (!all(c("dataset", "method", "time_secs", "mem_GB") %in% names(timing_frame)) ||
        any(!vapply(timing_frame, function(column) is.atomic(column) && is.null(dim(column)), logical(1))))
      stop("Execution-time input must be a compact timing table")
    path <- file.path(output_root, "execution_times.feather")
    if (!file.copy(timing$artifact_path, path, overwrite = FALSE)) stop("Cannot copy explicit timing metadata")
    .analysis_checked_file(path, timing$md5)
    control$execution_times$artifact_path <- path
    control$execution_times$parent_artifact_path <- timing$artifact_path
  }
  validate_analysis_release(ledger, control$expected_scope, require_projection = TRUE)
  for (dataset in unique(ledger$dataset)) {
    selected <- ledger$dataset == dataset
    checked <- load_analysis_release(ledger[selected, , drop = FALSE],
      control$expected_scope[match(.analysis_key(ledger[selected, , drop = FALSE]), .analysis_key(control$expected_scope)), , drop = FALSE],
      require_projection = TRUE)
    rm(checked)
  }
  .analysis_save_new(control, file.path(output_root, "analysis_release.rds"))
  paths <- .analysis_projection_paths(control, output_root)
  inventory <- data.frame(path = paths, md5 = unname(tools::md5sum(paths)))
  tmp <- tempfile(".projection-inventory-", tmpdir = output_root)
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(inventory, tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (!file.link(tmp, file.path(output_root, "files.tsv"))) stop("Cannot publish projection inventory")
  invisible(control)
}
