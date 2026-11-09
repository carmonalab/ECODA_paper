# These are the core separation metrics, all higher-is-better. In particular,
# calc_sep_score's LISI is already (N - LISI)/(N - 1); never invert it again.
evaluation_metrics <- function() {
  c("anosim_score", "mod_knn3_score", "cluster_score", "sil_score", "lisi_score")
}

evaluation_scoring_policy <- function() "calc_sep_score_five_v1_lisi_separation"

.evaluation_named_labels <- function(labels, ids, allow_drops = FALSE) {
  if (is.null(labels) || is.null(names(labels)) || anyDuplicated(names(labels)) ||
      anyNA(names(labels)) || any(!nzchar(trimws(names(labels)))) ||
      anyNA(labels) || any(!nzchar(trimws(as.character(labels)))))
    stop("Evaluation labels require unique, nonmissing sample names and values")
  if (!all(ids %in% names(labels)) || (!allow_drops && !setequal(names(labels), ids)))
    stop("Evaluation labels do not match the distance sample universe")
  labels[ids]
}

evaluation_distance <- function(distance) {
  matrix <- as.matrix(distance)
  ids <- rownames(matrix)
  if (!is.numeric(matrix) || !nrow(matrix) || nrow(matrix) != ncol(matrix) ||
      is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
      !identical(ids, colnames(matrix)) || any(!is.finite(matrix)) ||
      any(matrix < 0) || any(abs(diag(matrix)) > 1e-10) ||
      !isTRUE(all.equal(unname(matrix), unname(t(matrix)), tolerance = 1e-10)))
    stop("Evaluation requires finite symmetric distances with identical unique named axes")
  stats::as.dist(matrix)
}

.evaluation_fingerprint <- function(value) {
  path <- tempfile("ecoda-evaluation-")
  on.exit(unlink(path))
  saveRDS(value, path, compress = FALSE, version = 3)
  unname(tools::md5sum(path))
}

# No method fitting: score immutable distances against explicitly supplied axes.
# Primary NULL labels always come from the original bundle. Extra NULL axes
# require their own recorded seed (Kidney tissue uses 124).
evaluate_release_result <- function(
    result, dataset, view, method_key, evaluations, primary_evaluation,
    is_null = FALSE, null_seeds = integer(), stored_policy = NULL,
    score_fn = calc_sep_score, source_identity, allowed_dropped_sample_ids = character()) {
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) saved_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) assign(".Random.seed", saved_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.list(result) || is.null(result$dist_mat)) stop("Validated distance bundle is required")
  if (!is.list(evaluations) || is.null(names(evaluations)) || anyDuplicated(names(evaluations)) ||
      any(!nzchar(names(evaluations))) || !primary_evaluation %in% names(evaluations))
    stop("Explicit named evaluation axes and primary_evaluation are required")
  if (!is.list(source_identity) || !all(c("artifact_md5", "receipt_md5", "input_md5", "source_policy") %in% names(source_identity)) ||
      any(!vapply(source_identity[c("artifact_md5", "receipt_md5", "input_md5")],
        function(x) is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{32}$", x), logical(1))) ||
      !is.character(source_identity$source_policy) || length(source_identity$source_policy) != 1L ||
      is.na(source_identity$source_policy) || !nzchar(source_identity$source_policy))
    stop("Explicit validated artifact/receipt/input digests and source policy are required")
  if (!identical(source_identity$dataset, dataset) || !identical(source_identity$view, view) ||
      !identical(source_identity$method_key, method_key)) stop("Release identity disagrees with dataset/view/method")
  distance <- evaluation_distance(result$dist_mat)
  ids <- attr(distance, "Labels")
  distance_md5 <- .evaluation_fingerprint(distance)
  features_md5 <- if (is.null(result$feat_mat)) NULL else .evaluation_fingerprint(result$feat_mat)
  scitd <- grepl("^scITD(_|$)", method_key)
  dropped <- list()
  labels_used <- list()
  scores <- list()
  provenance <- list()
  for (evaluation in names(evaluations)) {
    true_labels <- .evaluation_named_labels(evaluations[[evaluation]], ids, scitd)
    dropped[[evaluation]] <- setdiff(names(evaluations[[evaluation]]), ids)
    if (!setequal(dropped[[evaluation]], allowed_dropped_sample_ids) ||
        (length(allowed_dropped_sample_ids) && !scitd))
      stop("Dropped samples require exact documented scITD-only permission")
    seed <- NA_integer_
    if (identical(evaluation, primary_evaluation)) {
      # Even complete historical NULL scores cannot justify fabricating labels.
      labels <- if (is.null(result$labels) && !is_null) true_labels else .evaluation_named_labels(result$labels, ids)
      if (is_null && !identical(sort(as.character(labels)), sort(as.character(true_labels))))
        stop("Stored primary NULL labels do not preserve the true axis category counts")
      if (!is_null && !identical(as.character(labels), as.character(true_labels)))
        stop("Stored primary labels disagree with the approved evaluation labels")
    } else if (is_null) {
      if (is.null(names(null_seeds)) || !evaluation %in% names(null_seeds))
        stop("Additional NULL axis requires its own explicit shuffle seed")
      seed <- null_seeds[[evaluation]]
      if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed))
        stop("NULL shuffle seed must be one finite integer")
      # Restore the caller's RNG state; this is evaluation-only shuffling.
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
      shuffled <- tryCatch({ set.seed(seed); sample(seq_along(true_labels)) },
                           finally = {
                             if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
                             else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                               rm(".Random.seed", envir = .GlobalEnv)
                           })
      labels <- true_labels[shuffled]
      names(labels) <- ids
    } else labels <- true_labels
    if (length(unique(labels)) < 2L || length(unique(labels)) >= length(labels))
      stop("Separation evaluation requires 2..n-1 biological groups")
    metrics <- evaluation_metrics()
    values <- setNames(rep(NA_real_, length(metrics)), metrics)
    reused <- rep(FALSE, length(metrics))
    if (identical(evaluation, primary_evaluation) &&
        identical(stored_policy, evaluation_scoring_policy()) && !is.null(result$labels) && !is.null(result$scores)) {
      for (j in seq_along(metrics)) {
        value <- result$scores[[metrics[j]]]
        if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
            (metrics[j] != "lisi_score" || (value >= 0 && value <= 1))) {
          values[j] <- value
          reused[j] <- TRUE
        }
      }
    }
    if (any(!reused)) {
      calculated <- score_fn(distance, labels)
      for (j in which(!reused)) {
        value <- calculated[[metrics[j]]]
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            (metrics[j] == "lisi_score" && (value < 0 || value > 1)))
          stop("Core scorer did not produce finite metric ", metrics[j])
        values[j] <- value
      }
    }
    labels_used[[evaluation]] <- labels
    scores[[evaluation]] <- data.frame(
      dataset = dataset, view = view, evaluation_id = evaluation,
      method_key = method_key, metric = metrics, value = unname(values),
      score_source = ifelse(reused, "stored_compatible", "rescored_distance"),
      stringsAsFactors = FALSE)
    provenance[[evaluation]] <- list(
      scoring_policy = evaluation_scoring_policy(), source_identity = source_identity,
      lisi_definition = "mean(clamp((N - standard_LISI)/(N - 1), 0, 1)); higher is better",
      primary = identical(evaluation, primary_evaluation), is_null = is_null,
      shuffle_seed = seed, labels = labels, labels_md5 = .evaluation_fingerprint(labels),
      distance_md5 = distance_md5, features_md5 = features_md5,
      dropped_sample_ids = dropped[[evaluation]])
  }
  features <- result$feat_mat
  if (!is.null(features) && (is.null(rownames(features)) || anyDuplicated(rownames(features)) ||
      !setequal(rownames(features), ids))) stop("Feature coordinates do not match distance samples")
  list(scores = do.call(rbind, scores), labels = labels_used, provenance = provenance,
       dist_mat = result$dist_mat, feat_mat = features)
}

# expected_grid explicitly enumerates every required dataset/view/axis/method/
# metric. Validate the whole finite grid before scaling or averaging anything.
aggregate_evaluation_scores <- function(scores, expected_grid) {
  keys <- c("dataset", "view", "evaluation_id", "method_key", "metric")
  key <- function(frame) do.call(paste, c(frame[keys], sep = "\r"))
  for (frame in list(scores, expected_grid)) {
    if (!is.data.frame(frame) || !nrow(frame) || !all(keys %in% names(frame)) ||
        anyNA(frame[keys]) || any(vapply(frame[keys], function(x) any(!nzchar(trimws(x))), logical(1))) ||
        anyDuplicated(key(frame))) stop("Evaluation grid has missing or duplicate keys")
    groups <- split(frame$metric, interaction(frame[keys[1:4]], drop = TRUE))
    if (any(!vapply(groups, function(x) setequal(x, evaluation_metrics()) && length(x) == 5L, logical(1))))
      stop("Every evaluation entry requires exactly the five approved metrics")
  }
  if (!setequal(key(scores), key(expected_grid)) || !is.numeric(scores$value) ||
      any(!is.finite(scores$value))) stop("Incomplete or nonfinite required evaluation grid")
  scores <- scores[match(key(expected_grid), key(scores)), , drop = FALSE]
  group <- interaction(scores[c("dataset", "view", "evaluation_id", "metric")], drop = TRUE)
  scores$scaled_value <- ave(scores$value, group, FUN = function(x) {
    limits <- range(x)
    if (limits[1] == limits[2]) rep(0.5, length(x)) else (x - limits[1]) / diff(limits)
  })
  entries <- stats::aggregate(scores$scaled_value,
    scores[c("dataset", "view", "evaluation_id", "method_key")], mean)
  names(entries)[ncol(entries)] <- "mean_score"
  # All methods must cover the same evaluation entries before a fair mean.
  for (view in unique(entries$view)) {
    selected <- entries[entries$view == view, , drop = FALSE]
    scopes <- split(paste(selected$dataset, selected$evaluation_id, sep = "\r"), selected$method_key)
    if (any(!vapply(scopes, function(x) setequal(x, scopes[[1]]), logical(1))))
      stop("Ranking methods must cover the same evaluation entries")
  }
  # Each evaluation entry has equal weight, including Kidney's two axes.
  summary <- stats::aggregate(entries$mean_score, entries[c("view", "method_key")], mean)
  names(summary)[ncol(summary)] <- "mean_score"
  summary$rank <- ave(-summary$mean_score, summary$view, FUN = function(x) rank(x, ties.method = "min"))
  summary <- summary[order(summary$view, summary$rank, summary$method_key), , drop = FALSE]
  list(scores = scores, entries = entries, ranking = summary,
       policy = "per-dataset/view/evaluation/metric minmax; constant=0.5; five equal weights; evaluation entries equal; descending arithmetic mean")
}
