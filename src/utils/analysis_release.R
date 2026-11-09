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

validate_analysis_release <- function(ledger, expected_scope) {
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
    path <- .analysis_checked_file(row$artifact_path, row$md5)
    if (!grepl("^/", row$run_root)) stop("run_root must be absolute")
    root <- normalizePath(row$run_root, mustWork = TRUE)
    if (!identical(root, row$run_root)) stop("run_root must be canonical, without aliases or symlinks")
    if (!dir.exists(root) || !startsWith(path, paste0(root, "/")))
      stop("Artifact is outside its explicit run_root: ", path)
    receipt_path <- .analysis_checked_file(row$receipt_path, row$receipt_md5)
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
    receipts[[i]] <- receipt
  }
  attr(ledger, "receipts") <- receipts
  ledger
}

load_analysis_release <- function(ledger, expected_scope) {
  ledger <- validate_analysis_release(ledger, expected_scope)
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
    results[[i]] <- list(result = result, source = row, receipt = receipts[[i]],
      identity = list(dataset = row$dataset, view = row$view, method_key = row$method_key,
                      artifact_md5 = row$md5, receipt_md5 = row$receipt_md5,
                      input_md5 = receipts[[i]]$input_md5, source_policy = receipts[[i]]$source_policy))
  }
  results
}
