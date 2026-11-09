source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/constants.R")
source("src/notebooks/benchmark_analysis.R")
source("src/utils/plotting.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(ggplot2) })
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = base::identity), "error"))
tmp <- tempfile("benchmark-release-"); dir.create(tmp); tmp <- normalizePath(tmp)
md5 <- function(path) unname(tools::md5sum(path))
ids <- paste0("s", 1:6)
labels <- setNames(rep(c("a", "b"), each = 3), ids)
features <- matrix(seq_len(18), 6, dimnames = list(ids, letters[1:3]))
bundle <- list(feat_mat = features, dist_mat = dist(features), labels = labels)
null <- bundle; null$labels <- labels[c(1, 4, 2, 5, 3, 6)]; names(null$labels) <- ids
metadata <- list(labels = labels, n_samples = 6L, n_cells = 60L,
  cells_per_sample = setNames(rep(10L, 6), ids), n_cell_types_high_res = 3L)
methods <- c("ECODA_authors_HR", "ECODA_authors_HR_NULL", "Pseudobulk_hvg2000", "Pseudobulk_hvg2000_bl")
aggregate <- file.path(tmp, "old.rds")
saveRDS(setNames(list(bundle, null, bundle), methods[1:3]), aggregate)
standalone <- file.path(tmp, "bl.rds"); saveRDS(bundle, standalone)
meta <- file.path(tmp, "metadata.rds"); saveRDS(metadata, meta)
evidence <- file.path(tmp, "review.txt"); writeLines("Bounded current-input/artifact review", evidence)
ledger <- data.frame(dataset = "D", view = "benchmark_analysis",
  method_key = c(methods, "__dataset_metadata__"),
  kind = c(rep("rds_bundle", 3), rep("standalone_result", 2)),
  artifact_path = c(rep(aggregate, 3), standalone, meta), bundle_key = c(methods[1:3], "", ""),
  run_root = tmp, md5 = "", receipt_path = file.path(tmp, paste0(seq_len(5), ".receipt.rds")), receipt_md5 = "")
ledger$md5 <- vapply(ledger$artifact_path, md5, character(1))
scope <- ledger[c("dataset", "view", "method_key")]
scope$require_new <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
scope$current_input_md5 <- strrep("a", 32); scope$source_policy <- "current-input-v1"
scope$status <- "approved"; scope$blocked_reason <- ""
receipt <- function(i, origin = if (i == 4) "new" else "legacy_reuse") {
  value <- c(as.list(ledger[i, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
    list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = origin,
      input_md5 = scope$current_input_md5[i], source_policy = scope$source_policy[i],
      evidence = "bounded review", evidence_path = evidence, evidence_md5 = md5(evidence)))
  saveRDS(value, ledger$receipt_path[i]); ledger$receipt_md5[i] <<- md5(ledger$receipt_path[i])
}
for (i in seq_len(nrow(ledger))) receipt(i)
control <- function() list(ledger = ledger, expected_scope = scope, reference_method_key = c(D = "ECODA_authors_HR"))
load_raw_benchmark_release <- load_benchmark_analysis_release
load_benchmark_analysis_release <- function(control, ...) {
  path <- tempfile("parent-control-", tmpdir = tmp); saveRDS(control, path)
  source <- file.path(tmp, "source.manifest"); runtime <- file.path(tmp, "runtime.identity")
  writeLines("test source", source); writeLines("test runtime", runtime)
  projected <- project_analysis_release(control, tempfile("projection-", tmpdir = tmp),
    c(source = source, runtime = runtime, control = path))
  load_raw_benchmark_release(projected, ...)
}
control_path <- file.path(tmp, "control.rds"); saveRDS(control(), control_path)
stopifnot(identical(read_benchmark_analysis_control(control_path, md5(control_path)), control()))
fail(read_benchmark_analysis_control(control_path, strrep("0", 32)))
fail(read_benchmark_analysis_control(file.path(tmp, ".", "control.rds"), md5(control_path)))
datasets <- list(D = list(label_col = "biology", evaluation_labels = c(primary = "biology")))
calls <- 0L
scorer <- function(distance, labels) {
  calls <<- calls + 1L
  c(setNames(as.list(c(0.6, 0.4, 0.2, 0.3, 0.7)), evaluation_metrics()), list(mod_knn6_score = 0.4))
}
before <- tools::md5sum(c(aggregate, standalone, meta)); set.seed(19); rng <- .Random.seed
release <- load_benchmark_analysis_release(control(), datasets, scorer)
stopifnot(calls == 4L, identical(.Random.seed, rng), nrow(release$core_scores) == 20L,
  identical(release$results$bmark$D$ECODA_authors_HR_NULL$labels, null$labels),
  identical(tools::md5sum(c(aggregate, standalone, meta)), before), nrow(release$statistics) == 1L,
  all(evaluation_metrics() %in% names(release$results$bmark$D$ECODA_authors_HR$scores)))
bad <- control(); bad$ledger <- bad$ledger[-4, ]; fail(load_benchmark_analysis_release(bad, datasets, scorer))
bad <- control(); bad$expected_scope$require_new[4] <- FALSE; fail(load_benchmark_analysis_release(bad, datasets, scorer))
receipt(4, "legacy_reuse"); fail(load_benchmark_analysis_release(control(), datasets, scorer)); receipt(4)
bad <- control(); bad$ledger <- bad$ledger[-1, ]; fail(load_benchmark_analysis_release(bad, datasets, scorer))
bad <- control(); bad$expected_scope$status[1] <- "blocked"; bad$expected_scope$blocked_reason[1] <- "Explicit pending incident"
fail(load_benchmark_analysis_release(bad, datasets, scorer))
fail(load_benchmark_analysis_release(control(), datasets, function(...) setNames(as.list(rep(NA_real_, 5)), evaluation_metrics())))
# Rehashed source changes still need scientifically consistent labels.
bad_bundle <- bundle; bad_bundle$labels[1] <- "b"; saveRDS(bad_bundle, standalone)
ledger$md5[4] <- md5(standalone); receipt(4)
fail(load_benchmark_analysis_release(control(), datasets, scorer))
saveRDS(bundle, standalone); ledger$md5[4] <- md5(standalone); receipt(4)
# scITD is the sole permitted sample-dropping method, and only with exact IDs.
reduced <- bundle; reduced$feat_mat <- features[-6, ]; reduced$dist_mat <- dist(reduced$feat_mat); reduced$labels <- labels[-6]
stopifnot(is.null(release$loaded), is.list(release$control))
identity <- release$provenance[[paste("D", methods[1], sep = "/")]]$primary$source_identity
identity$method_key <- "scITD_hvg2000_factors5"
fail(evaluate_release_result(reduced, "D", "benchmark_analysis", identity$method_key,
  list(primary = labels), "primary", score_fn = scorer, source_identity = identity))
stopifnot(nrow(evaluate_release_result(reduced, "D", "benchmark_analysis", identity$method_key,
  list(primary = labels), "primary", score_fn = scorer, source_identity = identity,
  allowed_dropped_sample_ids = "s6")$scores) == 5L)
identity$method_key <- "ordinary"
fail(evaluate_release_result(reduced, "D", "benchmark_analysis", "ordinary",
  list(primary = labels), "primary", score_fn = scorer, source_identity = identity, allowed_dropped_sample_ids = "s6"))
# Five metrics stay equally weighted. Constant comparisons receive 0.5.
frame <- release$core_scores; frame$method <- frame$method_key; frame$score <- frame$metric
clean <- prepare_evaluation_plot_data(frame, methods, expected_grid = release$expected_grid)
stopifnot(all(clean$scaled_value == 0.5), get_global_avg(frame, methods[1], methods) == 0.5)
fail(prepare_evaluation_plot_data(frame[-1, ], methods, expected_grid = release$expected_grid))
fail(prepare_evaluation_plot_data(frame, c(methods, "missing"), expected_grid = release$expected_grid))
bad <- frame; bad$value[1] <- Inf; fail(prepare_evaluation_plot_data(bad, methods))
fail(prepare_evaluation_plot_data(frame, c(methods, "ECODA_PB_combo_old")))
stopifnot(length(unique(method_label_map_main[c("ecoda_unsupervised_harmony", "ecoda_unsupervised_harmony_limma")])) == 2L)
analysis_methods <- methods; analysis_requested_methods <- NULL; analysis_scope <- "all"
stopifnot(identical(methods_in_scope(methods), methods))
fail(methods_in_scope(c(methods, "missing")))
analysis_requested_methods <- methods[1]
stopifnot(identical(methods_in_scope(c(methods, "missing")), methods[1]))
# Unequal raw ranges cannot cause one metric to dominate the normalized mean.
frame$value <- rep(c(-10, 3, -100, -0.2, 0.1), 4) + rep(0:3, each = 5)
clean <- prepare_evaluation_plot_data(frame, methods)
stopifnot(isTRUE(all.equal(as.numeric(tapply(clean$scaled_value, clean$method, mean)[methods]), c(0, 1/3, 2/3, 1))))
# The NULL comparator changes the range even if it is absent from drawn MDS
# points. The subtitle uses the same full pool as its heatmap entry.
frame$value[frame$method == "ECODA_authors_HR_NULL"] <- -1000
normalized <- prepare_evaluation_plot_data(frame, methods)
heatmap_entry <- mean(normalized$scaled_value[normalized$method == methods[1]])
mds_mean <- benchmark_normalized_mean(frame, "D", methods[1], normalization_methods = methods)
plotted_only <- benchmark_normalized_mean(frame, "D", methods[1], normalization_methods = setdiff(methods, "ECODA_authors_HR_NULL"))
stopifnot(identical(mds_mean, heatmap_entry), mds_mean != plotted_only,
  "ECODA_authors_HR_NULL" %in% benchmark_main_method_keys,
  "ECODA_authors_HR_NULL" %in% benchmark_extended_method_keys)
p <- plot_mds(dist(features), labels, cluster_score = FALSE, mod_score = FALSE, anosim_score = FALSE,
  separation_score = 0.5)
stopifnot(identical(p$labels$subtitle, "Five-metric min–max mean: 0.50"))
# Parse every maintained R chunk without executing its analyses or producing figures.
for (path in c("src/notebooks/benchmark_analysis.R", "src/utils/plotting.R", "src/utils/constants.R")) invisible(parse(path))
lines <- readLines("notebooks/benchmark_analysis.rmd"); active <- FALSE; chunks <- 0L; code <- character()
for (line in lines) {
  if (startsWith(line, "```{r")) { active <- TRUE; code <- character(); next }
  if (active && startsWith(line, "```")) { invisible(parse(text = code)); chunks <- chunks + 1L; active <- FALSE; next }
  if (active) code <- c(code, line)
}
stopifnot(!active, chunks > 50L)
# Cohort streaming is an observable boundary: neither the source projector nor
# maintained intake may pass more than one cohort to the payload reader.
two <- control(); extra_ledger <- ledger; extra_ledger$dataset <- "E"
extra_scope <- scope; extra_scope$dataset <- "E"
for (i in seq_len(nrow(extra_ledger))) {
  old_receipt <- readRDS(ledger$receipt_path[i]);old_receipt$dataset <- "E"
  extra_ledger$receipt_path[i] <- file.path(tmp,paste0("E",i,".receipt.rds"))
  saveRDS(old_receipt,extra_ledger$receipt_path[i]);extra_ledger$receipt_md5[i] <- md5(extra_ledger$receipt_path[i])
}
two$ledger <- rbind(ledger,extra_ledger);two$expected_scope <- rbind(scope,extra_scope)
two$reference_method_key <- c(D="ECODA_authors_HR",E="ECODA_authors_HR")
reader <- load_analysis_release; cohorts_read <- list()
load_analysis_release <- function(ledger,...) {
  cohorts_read[[length(cohorts_read)+1L]] <<- unique(ledger$dataset)
  reader(ledger,...)
}
streamed <- load_benchmark_analysis_release(two,c(datasets,list(E=datasets$D)),scorer)
load_analysis_release <- reader
stopifnot(length(cohorts_read)>=4L,all(lengths(cohorts_read)==1L),is.null(streamed$loaded),
  nrow(streamed$core_scores)==40L,is.null(streamed$results$bmark$E$Pseudobulk_hvg2000$feat_mat),
  streamed$results$bmark$E$Pseudobulk_hvg2000$n_features==3L)
unlink(tmp, recursive = TRUE)
cat("benchmark release reports: PASS; parsed", chunks, "R chunks\n")
