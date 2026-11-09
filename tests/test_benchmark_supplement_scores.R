# Compact reviewed summaries only; no model, RNA or production figure execution.
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/constants.R")
source("src/notebooks/benchmark_analysis.R")
calc_sep_score <- datrans <- function(...) stop("Reporting must not fit or recalculate scores")
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = base::identity), "error"))
tmp <- tempfile("ecoda-supplements-"); dir.create(tmp); tmp <- normalizePath(tmp)
md5 <- function(path) unname(tools::md5sum(path))
metrics <- evaluation_metrics()
values <- setNames(as.list(c(0.7, 0.4, 0.6, -0.2, 0.8)), metrics)
summary <- list(counts_all_0.5 = c(values, list(mod_knn6_score = 123, mod_knn9_score = Inf, mod_knnsqrtn_score = 456)),
                counts_all_1 = list(multLN = values))
evidence <- file.path(tmp, "review.txt"); writeLines("Bounded current-input summary review", evidence)
ledger <- data.frame(dataset = c("A", "B"), view = "benchmark_analysis", method_key = "__zeroimp__",
  kind = "standalone_result", artifact_path = file.path(tmp, c("A.rds", "B.rds")), bundle_key = "", run_root = tmp,
  md5 = "", receipt_path = file.path(tmp, c("A.receipt.rds", "B.receipt.rds")), receipt_md5 = "")
scope <- ledger[c("dataset", "view", "method_key")]
scope$require_new <- FALSE; scope$current_input_md5 <- strrep("a", 32)
scope$source_policy <- "current-input-v1"; scope$status <- "approved"
for (i in seq_len(nrow(ledger))) {
  saveRDS(summary, ledger$artifact_path[i]); ledger$md5[i] <- md5(ledger$artifact_path[i])
  receipt <- c(as.list(ledger[i, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
    list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = "legacy_reuse",
      input_md5 = scope$current_input_md5[i], source_policy = scope$source_policy[i], evidence = "bounded review",
      evidence_path = evidence, evidence_md5 = md5(evidence)))
  saveRDS(receipt, ledger$receipt_path[i]); ledger$receipt_md5[i] <- md5(ledger$receipt_path[i])
}
project_screen_control <- function(ledger, scope) {
  control <- list(ledger=ledger, expected_scope=scope)
  path <- tempfile("parent-control-", tmpdir=tmp); saveRDS(control,path)
  source <- file.path(tmp,"source.manifest"); runtime <- file.path(tmp,"runtime.identity")
  writeLines("fixture source",source);writeLines("fixture runtime",runtime)
  list(control=project_analysis_release(control,tempfile("projection-",tmpdir=tmp),
    c(source=source,runtime=runtime,control=path)))
}
reject_corrupt_screen <- function(release, analysis, mutate) {
  path <- release$control$ledger$artifact_path[1]
  old <- readRDS(path)
  on.exit(saveRDS(old,path))
  altered <- old; altered[[1L]] <- mutate(altered[[1L]])
  saveRDS(altered,path)
  if(analysis=="trans") reject(benchmark_transformation_scores(release,c("A","B")))
  else reject(benchmark_zeroimp_scores(release,c("A","B")))
}
release <- project_screen_control(ledger,scope)
before <- tools::md5sum(ledger$artifact_path)
frame <- benchmark_zeroimp_scores(release, c("A", "B"))
stopifnot(nrow(frame) == 20L, setequal(unique(frame$score), metrics),
  all(frame$value[frame$score == "lisi_score"] == 0.8),
  all(frame$value[frame$score == "mod_knn3_score"] == 0.4),
  identical(tools::md5sum(ledger$artifact_path), before))
bad <- release; bad$control$ledger <- bad$control$ledger[-1,]; reject(benchmark_zeroimp_scores(bad, c("A", "B")))
reject_corrupt_screen(release,"zeroimp",function(x){x$counts_all_0.5$lisi_score<-NULL;x})
reject_corrupt_screen(release,"zeroimp",function(x){x$counts_all_0.5$sil_score<-NA_real_;x})
reject_corrupt_screen(release,"zeroimp",function(x){x$counts_all_1<-NULL;x})
reject(benchmark_validate_screen_scores(frame[-1, ]))
reject(benchmark_validate_screen_scores(rbind(frame, frame[1, ])))
panels <- lapply(metrics, function(metric) {
  panel <- plot_benchmark_raw_metric(frame, metric)
  stopifnot(nrow(panel$data) == 4L, all(panel$data$score == metric),
    identical(panel$labels$y, paste("Raw", score_facet_labels[[metric]])))
  ggplot2::ggplot_build(panel)
  panel
})
# Actual small temporary plot exercises the shared screen rendering seam.
ggplot2::ggsave(file.path(tmp, "silhouette.pdf"), panels[[4]], width = 5, height = 4)
stopifnot(file.info(file.path(tmp, "silhouette.pdf"))$size > 0)
reject(plot_benchmark_raw_metric(frame, "mod_knn9_score"))
# The reporter delegates the versioned scientific contract to the producer's
# pure validator. Its independent tests cover reconstruction and ALR semantics.
aliases <- setNames(c("ANOSIM_score", "Modularity_score", "Adjusted_Rand_Index", "Silhouette_score", "LISI_score"), metrics)
contract <- new.env(parent = globalenv())
sys.source("src/5_run_benchmark_methods/benchmark_pipeline.R", envir = contract)
transformation_analysis_policy <- contract$transformation_analysis_policy
stopifnot(identical(transformation_analysis_policy()$score_columns, aliases))
validated <- 0L
validate_transformation_result <- function(table) {
  validated <<- validated + 1L
  contract$validate_transformation_result(table)
}
trans <- data.frame(trans_method = c("counts", "freq", "arcsine_sqrt", "alr_randref", "alr_mincvref", "clr"))
for (metric in metrics) trans[[aliases[[metric]]]] <- values[[metric]]
attr(trans, "transformation_policy") <- transformation_analysis_policy()
counts <- matrix(seq_len(18), 6, dimnames = list(paste0("s", 1:6), c("A", "B", "C")))
labels <- setNames(rep(c("a", "b"), each = 3), rownames(counts))
draws <- data.frame(trans_method = rep(trans$trans_method, each = 20), zero_imp_method = "counts_all__1",
  n_celltypes = 0, Amount_of_perturbation = 0, bootstrap_id = rep(1:20, 6), alr_reference = NA_character_)
features <- lapply(draws$trans_method, function(method) {
  matrix <- counts + 1
  frequency <- matrix / rowSums(matrix)
  switch(method, counts = matrix, freq = 100 * frequency, arcsine_sqrt = asin(sqrt(frequency)),
    clr = log(matrix) - rowMeans(log(matrix)),
    alr_randref = log(matrix[, -1, drop = FALSE] / matrix[, 1]),
    alr_mincvref = log(matrix[, -1, drop = FALSE] / matrix[, 1]))
})
draws$alr_reference[grepl("^alr_", draws$trans_method)] <- "A"
draws$feat_mat <- I(features); draws$sample_labels <- I(rep(list(labels), nrow(draws)))
set.seed(123, kind = "L'Ecuyer-CMRG")
draws$rng_feature_state <- I(rep(list(.Random.seed), nrow(draws)))
draws$rng_score_state <- I(rep(list(.Random.seed), nrow(draws)))
for (metric in metrics) draws[[aliases[[metric]]]] <- values[[metric]]
attr(trans, "transformation_replicates") <- draws
attr(trans, "composition_counts") <- counts
attr(trans, "sample_labels") <- labels
trans_ledger <- ledger; trans_ledger$method_key <- "__trans__"
trans_ledger$artifact_path <- file.path(tmp, c("A_trans.rds", "B_trans.rds"))
trans_ledger$receipt_path <- paste0(trans_ledger$artifact_path, ".receipt.rds")
trans_scope <- scope; trans_scope$method_key <- "__trans__"; trans_scope$require_new <- TRUE
for (i in seq_len(nrow(trans_ledger))) {
  saveRDS(trans, trans_ledger$artifact_path[i]); trans_ledger$md5[i] <- md5(trans_ledger$artifact_path[i])
  receipt <- c(as.list(trans_ledger[i, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
    list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = "new",
      input_md5 = trans_scope$current_input_md5[i], source_policy = trans_scope$source_policy[i], evidence = "new producer review",
      evidence_path = evidence, evidence_md5 = md5(evidence)))
  saveRDS(receipt, trans_ledger$receipt_path[i]); trans_ledger$receipt_md5[i] <- md5(trans_ledger$receipt_path[i])
}
trans_release <- project_screen_control(trans_ledger,trans_scope)
validated <- 0L
trans_frame <- benchmark_transformation_scores(trans_release, c("A", "B"))
stopifnot(nrow(trans_frame) == 60L, validated == 2L, setequal(unique(trans_frame$score), metrics),
  all(trans_frame$value[trans_frame$score == "mod_knn3_score"] == 0.4),
  all(trans_frame$value[trans_frame$score == "lisi_score"] == 0.8))
bad <- trans_release
receipt <- readRDS(bad$control$ledger$receipt_path[1]);receipt$origin<-"legacy_reuse"
receipt_path<-tempfile("bad-origin-",tmpdir=tmp);saveRDS(receipt,receipt_path)
bad$control$ledger$receipt_path[1]<-receipt_path;bad$control$ledger$receipt_md5[1]<-md5(receipt_path)
reject(benchmark_transformation_scores(bad, c("A", "B")))
reject_corrupt_screen(trans_release,"trans",function(x){attr(x,"transformation_policy")<-NULL;x})
reject_corrupt_screen(trans_release,"trans",function(x){x$LISI_score<-NULL;x})
reject_corrupt_screen(trans_release,"trans",function(x){x$Silhouette_score[1]<-NA_real_;x})
reject_corrupt_screen(release,"zeroimp",function(x){x$counts_all_0.5$lisi_score<-2;x})
# Parse maintained screens without executing notebook or historical diagnostics.
lines <- readLines("notebooks/benchmark_analysis.rmd"); active <- FALSE; chunks <- 0L; code <- character()
for (line in lines) {
  if (startsWith(line, "```{r")) { active <- TRUE; code <- character(); next }
  if (active && startsWith(line, "```")) { invisible(parse(text = code)); chunks <- chunks + 1L; active <- FALSE; next }
  if (active) code <- c(code, line)
}
stopifnot(!active, chunks >= 62L)
unlink(tmp, recursive = TRUE)
cat("PASS: reviewed five-metric supplemental grids, raw panels, no diagnostic substitution or fitting;", chunks, "chunks parsed\n")
