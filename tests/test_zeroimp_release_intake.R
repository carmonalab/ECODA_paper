# Actual counts-screen package -> checked A projection -> supplemental reader.
# No H5AD, RNA, fits, score calculation or figure writes during intake.
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/constants.R")
source("src/5_run_benchmark_methods/run_local_ecoda_derived.R")
source("src/notebooks/run_zeroimp_count_screen.R")
source("src/notebooks/benchmark_analysis.R")
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
root <- tempfile("zeroimp-intake-"); dir.create(root); root <- normalizePath(root)
md5 <- function(p) unname(tools::md5sum(p))
ids <- paste0("s", 1:6)
counts <- matrix(as.double(1:18), 6, dimnames = list(ids, paste0("CT", 1:3)))
labels <- setNames(factor(rep(c("A", "B"), each = 3)), ids)
recipe <- function(df, labels) list(counts_all_0.5 = calc_sep_score(stats::dist(as.matrix(df) + .5), labels),
  counts_all_1 = calc_sep_score(stats::dist(as.matrix(df) + 1), labels))
score <- function(...) setNames(as.list(c(.4, .3, .2, -.1, .8, 5, 6, 9)),
  c(evaluation_metrics(), "mod_knnsqrtn_score", "mod_knn6_score", "mod_knn9_score"))
a <- zeroimp_run_cohort(list(counts = counts, labels = labels), score_fun = score, recipe_fun = recipe)
a$schema <- zeroimp_policy()$schema; a$dataset <- "A"; a$source <- list(dataset = "A", parent_receipts = list(origin = "legacy_reuse"))
a$run_id <- "new-screen"; a$config_md5 <- strrep("a", 32); a$runner_md5 <- strrep("b", 32); a$release_md5 <- strrep("c", 32)
b <- a; b$dataset <- "B"; b$source$dataset <- "B"
# A different expected-invalid recipe does not change the comparison cohort denominator.
b$recipes$counts_all_1 <- b$recipes$asinsqrt; b$recipes$counts_all_1$setting <- "counts_all_1"
i <- match("counts_all_1", b$results$setting); j <- match("asinsqrt", b$results$setting)
b$results[i, ] <- b$results[j, ]; b$results$setting[i] <- "counts_all_1"
.analysis_zeroimp_screen(a); .analysis_zeroimp_screen(b)
evidence <- file.path(root, "review.txt"); writeLines("Reviewed current all43 counts screen; parent origin stays legacy", evidence)
ledger <- data.frame(dataset = c("A", "B"), view = "benchmark_analysis", method_key = "__zeroimp__",
  kind = "standalone_result", artifact_path = file.path(root, c("A.rds", "B.rds")), bundle_key = "", run_root = root,
  md5 = "", receipt_path = file.path(root, c("A.receipt.rds", "B.receipt.rds")), receipt_md5 = "")
scope <- ledger[c("dataset", "view", "method_key")]; scope$require_new <- TRUE
scope$current_input_md5 <- strrep("d", 32); scope$source_policy <- "reviewed-all11-count-screen-v1"; scope$status <- "approved"
for (i in 1:2) {
  saveRDS(list(a, b)[[i]], ledger$artifact_path[i]); ledger$md5[i] <- md5(ledger$artifact_path[i])
  receipt <- c(as.list(ledger[i, c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
    list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = "new", input_md5 = scope$current_input_md5[i],
      source_policy = scope$source_policy[i], evidence = "complete current recipe review", evidence_path = evidence, evidence_md5 = md5(evidence)))
  saveRDS(receipt, ledger$receipt_path[i]); ledger$receipt_md5[i] <- md5(ledger$receipt_path[i])
}
control <- list(ledger = ledger, expected_scope = scope)
parent <- file.path(root, "control.rds"); saveRDS(control, parent)
source <- file.path(root, "source.manifest"); runtime <- file.path(root, "runtime.identity")
writeLines("source", source); writeLines("runtime", runtime)
before <- tools::md5sum(c(ledger$artifact_path, ledger$receipt_path))
# Intake and projection cannot calculate any scores or refit recipes.
calc_sep_score <- run_zeroimp_analysis <- zeroimp_run_cohort <- function(...) stop("Intake must not compute")
projected <- project_analysis_release(control, file.path(root, "projection"), c(source = source, runtime = runtime, control = parent))
items <- load_analysis_release(projected$ledger, projected$expected_scope, require_projection = TRUE)
for (item in items) stopifnot(is.null(item$result[["counts"]]), identical(item$result$source$parent_receipts$origin, "legacy_reuse"),
  identical(item$result$policy$recipe_menu, zeroimp_menu()), length(item$result$recipes) == 43L,
  is.null(item$result$recipes$asinsqrt$scores), is.null(item$result$recipes$asinsqrt$dist_mat),
  length(item$result$recipes$counts_all_0.5$scores) == 8L)
warnings <- character()
frame <- withCallingHandlers(benchmark_zeroimp_scores(list(control = projected), c("A", "B")),
  warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
stopifnot(nrow(frame) == 10L, identical(unique(frame$method), "counts_all_0.5"), length(warnings) == 1L,
  nrow(attr(frame, "recipe_status")) == 86L, nrow(attr(frame, "all_recipe_scores")) == 430L,
  "counts_all_1" %in% attr(frame, "excluded_methods"), length(attr(frame, "excluded_methods")) == 42L,
  all(frame$value[frame$score == "lisi_score"] == .8))
full <- attr(frame, "all_recipe_scores")
stopifnot(all(is.na(full$value[full$status != "PASS"])),
  sum(full$dataset == "A" & full$method == "counts_all_1" & full$status == "PASS") == 5L,
  sum(full$dataset == "B" & full$method == "counts_all_1" & is.na(full$value)) == 5L)
# All original package/receipt bytes stay unchanged, and failures cannot be fabricated away.
stopifnot(identical(before, tools::md5sum(names(before))))
x <- a; x$results <- x$results[-1, ]; reject(.analysis_zeroimp_screen(x))
x <- a; x$recipes$asinsqrt$scores <- score(); reject(.analysis_zeroimp_screen(x))
x <- a; x$recipes$counts_all_0.5$labels <- rev(x$labels); reject(.analysis_zeroimp_screen(x))
x <- a; x$recipes$counts_all_0.5$sample_ids <- rev(ids); reject(.analysis_zeroimp_screen(x))
x <- a; x$policy$seed <- 124L; reject(.analysis_zeroimp_screen(x))
x <- a; x$results$lisi_score[1] <- .2; reject(.analysis_zeroimp_screen(x))
x <- a; x$schema <- "unknown"; reject(.analysis_thin_result(list(result = x, source = list(method_key = "__zeroimp__"), receipt = list(origin = "new"))))
reject(.analysis_thin_result(list(result = a, source = list(method_key = "__zeroimp__"), receipt = list(origin = "legacy_reuse"))))
reject(.analysis_thin_result(list(result = a, source = list(method_key = "__zeroimp__", dataset = "B"), receipt = list(origin = "new"))))
unlink(root, recursive = TRUE)
cat("PASS typed zeroimp projection/intake, all43 statuses, common full-cohort menu, no scores/refits, unchanged origins and source bytes\n")
