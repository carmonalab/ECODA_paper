library(dplyr)
source("src/utils/scoring_metrics.R")
source("src/utils/evaluation_scores.R")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
release_identity <- list(dataset = "Kidney_KPMP_full", view = "batch_effect_uncorrected", method_key = "ECODA_authors_HR",
                 artifact_md5 = strrep("a", 32), receipt_md5 = strrep("b", 32),
                 input_md5 = strrep("c", 32), source_policy = "input-v1")
ids <- paste0("s", 1:24)
features <- cbind(x = rep(c(0, 4, 8), each = 8) + sin(1:24)/5, y = cos(1:24))
rownames(features) <- ids
labels <- setNames(rep(c("A", "B", "C"), each = 8), ids)
tissue <- setNames(rep(c("kidney", "cortex", "medulla"), 8), ids)
distance <- dist(features)
# Actual core parity, including transformed LISI and mean Ward/PAM ARI.
core <- calc_sep_score(distance, labels, parallel = 1)
stopifnot(all(is.finite(unlist(core[evaluation_metrics()]))))
result <- list(dist_mat = distance, feat_mat = features, labels = labels, scores = core)
count <- 0L
scorer <- function(d, l) { count <<- count + 1L; calc_sep_score(d, l, parallel = 1) }
evaluate <- function(result = result, axes = list(condition = labels), ...) {
  evaluate_release_result(result, "Kidney_KPMP_full", "batch_effect_uncorrected", "ECODA_authors_HR",
    axes, "condition", source_identity = release_identity, score_fn = scorer, ...)
}
# Compatible values preserved exactly, without rescoring/inverting LISI.
observed <- evaluate(result, stored_policy = evaluation_scoring_policy())
stopifnot(count == 0L, identical(observed$scores$value, unname(unlist(core[evaluation_metrics()]))))
missing <- result; missing$scores$sil_score <- NULL; missing$scores$lisi_score <- NULL
observed <- evaluate(missing, stored_policy = evaluation_scoring_policy())
stopifnot(count == 1L, sum(observed$scores$score_source == "rescored_distance") == 2L,
          isTRUE(all.equal(observed$scores$value, unname(unlist(core[evaluation_metrics()])))))
feather_result <- list(dist_mat = distance)
observed <- evaluate(feather_result)
stopifnot(identical(observed$labels$condition, labels), all(observed$scores$score_source == "rescored_distance"))
count <- 0L
observed <- evaluate(result, stored_policy = "legacy_lisi_unknown")
stopifnot(count == 1L, all(observed$scores$score_source == "rescored_distance"))
# Axis relabeling reuses exactly the original features/distances.
observed <- evaluate(result, list(condition = labels[rev(ids)], tissue = tissue), stored_policy = evaluation_scoring_policy())
stopifnot(identical(observed$dist_mat, distance), identical(observed$feat_mat, features),
          identical(observed$labels$condition, labels), identical(observed$labels$tissue, tissue),
          identical(observed$provenance$condition$distance_md5, observed$provenance$tissue$distance_md5),
          identical(observed$provenance$condition$features_md5, observed$provenance$tissue$features_md5))
# Preserve the actual historical primary NULL shuffle, never regenerate it.
set.seed(92); null_labels <- labels[sample(seq_along(labels))]; names(null_labels) <- ids
null <- result; null$labels <- null_labels; null$scores <- calc_sep_score(distance, null_labels, parallel = 1)
set.seed(777); rng_before <- .Random.seed
null_eval <- evaluate(null, list(condition = labels, tissue = tissue), is_null = TRUE,
                      null_seeds = c(tissue = 124), stored_policy = evaluation_scoring_policy())
stopifnot(identical(.Random.seed, rng_before), identical(null_eval$labels$condition, null_labels),
          null_eval$provenance$tissue$shuffle_seed == 124,
          identical(null_eval$dist_mat, distance), identical(null_eval$feat_mat, features))
set.seed(124); expected_tissue <- tissue[sample(seq_along(tissue))]; names(expected_tissue) <- ids
stopifnot(identical(null_eval$labels$tissue, expected_tissue),
          identical(null_eval$scores$value[1:5], unname(unlist(null$scores[evaluation_metrics()]))))
fail(evaluate(null, list(condition = labels, tissue = tissue), is_null = TRUE))
no_labels <- null; no_labels$labels <- NULL
fail(evaluate(no_labels, is_null = TRUE, stored_policy = evaluation_scoring_policy()))
wrong <- result; wrong$labels[1] <- "other"; fail(evaluate(wrong))
bad_null <- null; bad_null$labels[1] <- "other"
fail(evaluate(bad_null, is_null = TRUE, stored_policy = evaluation_scoring_policy()))
bad_null <- null; bad_null$labels[1] <- bad_null$labels[which(bad_null$labels != bad_null$labels[1])[1]]
fail(evaluate(bad_null, is_null = TRUE, stored_policy = evaluation_scoring_policy()))
wrong_lisi <- result; wrong_lisi$scores$lisi_score <- 2
fixed_lisi <- evaluate(wrong_lisi, stored_policy = evaluation_scoring_policy())
stopifnot(fixed_lisi$scores$value[5] == core$lisi_score,
          fixed_lisi$scores$score_source[5] == "rescored_distance")
bad_axes <- labels; names(bad_axes)[1] <- names(bad_axes)[2]; fail(evaluate(result, list(condition = bad_axes)))
bad <- result; bad$dist_mat <- as.matrix(distance); bad$dist_mat[1,2] <- NA_real_; fail(evaluate(bad))
# Only scITD may omit exactly the explicitly documented IDs.
dropped <- result; dropped$dist_mat <- as.dist(as.matrix(distance)[-24, -24]); dropped$labels <- labels[-24]; dropped$feat_mat <- features[-24, ]
fail(evaluate(dropped))
scitd_identity <- release_identity; scitd_identity$dataset <- "D"; scitd_identity$view <- "V"; scitd_identity$method_key <- "scITD_hvg2000_factors5"
fail(evaluate_release_result(dropped, "D", "V", "scITD_hvg2000_factors5", list(condition = labels), "condition",
  source_identity = scitd_identity, score_fn = scorer))
scitd <- evaluate_release_result(dropped, "D", "V", "scITD_hvg2000_factors5", list(condition = labels), "condition",
  source_identity = scitd_identity, score_fn = scorer, allowed_dropped_sample_ids = "s24")
stopifnot(identical(scitd$provenance$condition$dropped_sample_ids, "s24"))
# A complete five-metric grid: Kidney's two axes each receive the same weight
# as another dataset evaluation, rather than averaging Kidney first.
grid <- expand.grid(dataset = c("Other", "Kidney"), evaluation_id = c("condition", "tissue"),
  method_key = c("a", "b"), metric = evaluation_metrics(), stringsAsFactors = FALSE)
grid <- grid[grid$dataset == "Kidney" | grid$evaluation_id == "condition", ]
grid$view <- "uncorrected"
scores <- grid
scores$value <- as.numeric((scores$method_key == "a") == (scores$dataset == "Other"))
ranked <- aggregate_evaluation_scores(scores, grid)
stopifnot(ranked$ranking$method_key[1] == "b", ranked$ranking$mean_score[1] == 2/3)
constant <- scores; constant$value <- 4
neutral <- aggregate_evaluation_scores(constant, grid)
stopifnot(all(neutral$scores$scaled_value == .5), all(neutral$ranking$mean_score == .5), all(neutral$ranking$rank == 1))
fail(aggregate_evaluation_scores(scores[-1, ], grid))
fail(aggregate_evaluation_scores(rbind(scores, scores[1, ]), grid))
bad <- scores; bad$value[1] <- Inf; fail(aggregate_evaluation_scores(bad, grid))
bad <- grid; bad$metric[1] <- "unknown"; fail(aggregate_evaluation_scores(scores, bad))
cat("evaluation five-metric policies: OK\n")
