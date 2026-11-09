suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(arrow); library(funkyheatmap)})
source("src/utils/datasets_io.R")
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/batch_effect_analysis.R")
source("src/utils/constants.R")
source("src/utils/plotting.R")
fail <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
root <- tempfile("batch-release-"); dir.create(root); root <- normalizePath(root)
md5 <- function(path) unname(tools::md5sum(path))
evidence <- file.path(root, "review.txt"); writeLines("Synthetic current-input/release identity review", evidence)
ids <- paste0("s", 1:12)
truth <- setNames(rep(c("a", "b", "c"), each = 4), ids)
features <- cbind(x = 1:12, y = cos(1:12)); rownames(features) <- ids
shuffled <- truth[c(12:1)]; names(shuffled) <- ids
ledger <- list(); expected <- list(); index <- 0L
for (view in c("batch_effect_uncorrected", "batch_effect_corrected")) {
  datasets <- read_datasets_json(view = view)
  datasets <- datasets[setdiff(names(datasets)[vapply(datasets, function(ds) isTRUE(ds$use_for_batch_effect), logical(1))], "_debug")]
  stopifnot(length(datasets) == 8L, !"Diabetes" %in% names(datasets))
  for (dataset in names(datasets)) {
    ds <- datasets[[dataset]]
    metadata <- data.frame(Sample = ids)
    for (column in unlist(ds$evaluation_labels)) metadata[[column]] <- unname(truth)
    if (dataset == "Kidney_KPMP_full") metadata$tissue <- rep(c("kidney", "cortex", "medulla"), 4)
    for (column in unlist(ds$batch_col)) metadata[[column]] <- rep(c("batch1", "batch2"), 6)
    metadata_path <- file.path(root, paste0(dataset, "_", view, "_metadata.feather"))
    arrow::write_feather(metadata, metadata_path)
    method_keys <- batch_release_method_keys(view)
    bundles <- setNames(lapply(method_keys, function(method) list(dist_mat = dist(features), feat_mat = features,
      labels = if (method == "ECODA_authors_HR_NULL") shuffled else truth)), method_keys)
    aggregate_path <- file.path(root, paste0(dataset, "_", view, ".rds")); saveRDS(bundles, aggregate_path)
    for (method in c(method_keys, "__sample_metadata__")) {
      index <- index + 1L
      row <- data.frame(dataset = dataset, view = view, method_key = method,
        kind = if (method == "__sample_metadata__") "metadata_feather" else "rds_bundle",
        artifact_path = if (method == "__sample_metadata__") metadata_path else aggregate_path,
        bundle_key = if (method == "__sample_metadata__") "" else method,
        sample_id_col = "Sample", run_root = root, md5 = "",
        receipt_path = file.path(root, paste0("receipt_", index, ".rds")), receipt_md5 = "")
      if (method == "ECODA_authors_HR") {
        row$kind <- "standalone_result"; row$bundle_key <- ""
        row$artifact_path <- file.path(root, paste0(dataset, "_", view, "_newHR.rds"))
        saveRDS(bundles[[method]], row$artifact_path)
      }
      row$md5 <- md5(row$artifact_path)
      receipt <- c(as.list(row[c("dataset", "view", "method_key", "kind", "artifact_path", "bundle_key", "run_root", "md5")]),
        list(schema = "ecoda_analysis_release_receipt_v1", status = "reviewed", origin = "new",
             input_md5 = strrep("a", 32), source_policy = "cohort-v1", evidence = "synthetic reviewed identities",
             evidence_path = evidence, evidence_md5 = md5(evidence)))
      saveRDS(receipt, row$receipt_path); row$receipt_md5 <- md5(row$receipt_path)
      ledger[[index]] <- row
      scope <- row[c("dataset", "view", "method_key")]
      scope$require_new <- TRUE; scope$current_input_md5 <- strrep("a", 32)
      scope$source_policy <- "cohort-v1"; scope$status <- "approved"; scope$blocked_reason <- ""
      expected[[index]] <- scope
    }
  }
}
control <- list(schema = "ecoda_batch_analysis_release_v1", status = "approved", release_id = "tiny-v1",
                ledger = do.call(rbind, ledger), expected_scope = do.call(rbind, expected))
parent_control_path <- file.path(root, "parent-release.rds"); saveRDS(control, parent_control_path)
# Maintained local intake rejects original bundles before result deserialization.
raw_error <- tryCatch(read_batch_release_control(parent_control_path, md5(parent_control_path),
  "batch_effect_uncorrected"), error = identity)
stopifnot(inherits(raw_error, "error"),
          grepl("reviewed metadata projection", conditionMessage(raw_error), fixed = TRUE))
validation_source <- file.path(root, "validation-source.txt")
validation_runtime <- file.path(root, "validation-runtime.txt")
writeLines("Synthetic immutable source identity", validation_source)
writeLines("Synthetic immutable runtime identity", validation_runtime)
control <- project_analysis_release(control, file.path(root, "projection"),
  c(source = validation_source, runtime = validation_runtime, control = parent_control_path))
control_path <- file.path(root, "projection", "analysis_release.rds")
seen <- list()
scorer <- function(distance, labels) {
  seen[[length(seen) + 1L]] <<- labels
  as.list(setNames(c(.2, .3, .4, .5, .6), evaluation_metrics()))
}
left_control <- read_batch_release_control(control_path, md5(control_path), "batch_effect_uncorrected")
right_control <- read_batch_release_control(control_path, md5(control_path), "batch_effect_corrected")
left <- score_batch_release(left_control, scorer); right <- score_batch_release(right_control, scorer)
stopifnot(nrow(left$scores) == 360L, nrow(right$scores) == 405L,
          nrow(left$entries) == 72L, nrow(right$entries) == 81L,
          all(left$ranking$mean_score == .5), all(right$ranking$mean_score == .5),
          !"ECODA_seuratres_2" %in% right$scores$method_key,
          all(c("ecoda_unsupervised_harmony", "ecoda_unsupervised_harmony_limma") %in% right$scores$method_key),
          identical(left$provenance$Kidney_KPMP_full$ECODA_authors_HR_NULL$condition$labels, shuffled),
          left$provenance$Kidney_KPMP_full$ECODA_authors_HR_NULL$tissue$shuffle_seed == 124L)
unc <- load_batch_release_dataset(left_control, "Kidney_KPMP_full", scorer)
corr <- load_batch_release_dataset(right_control, "Kidney_KPMP_full", scorer)
labels <- setNames(batch_release_method_keys("batch_effect_uncorrected"), batch_release_method_keys("batch_effect_uncorrected"))
panels <- plot_batch_release_mds(unc, left_control$datasets$Kidney_KPMP_full$batch_col, labels, scorer)
stopifnot(all(c("biological", "evaluation_tissue") %in% names(panels)))
for (i in seq_len(length(labels) - 1L)) stopifnot(
  identical(panels$biological[[i]]$data[c("Dim1", "Dim2")], panels$evaluation_tissue[[i]]$data[c("Dim1", "Dim2")]),
  identical(as.character(panels$evaluation_tissue[[i]]$data$Cluster), unc$metadata$tissue))
# Biological MDS summaries use the exact heatmap entry mean, including NULL
# in each min-max range despite excluding it from the plotted method panels.
varying_scorer <- function(distance, labels) {
  base <- if (identical(labels, truth)) .4 else .2
  as.list(setNames(base + seq(0, .4, .1), evaluation_metrics()))
}
scaled_unc <- load_batch_release_dataset(left_control, "Kidney_KPMP_full", varying_scorer)
scaled_panels <- plot_batch_release_mds(scaled_unc, left_control$datasets$Kidney_KPMP_full$batch_col, labels, varying_scorer)
for (axis in c("condition", "tissue")) {
  panel <- if (axis == "condition") "biological" else "evaluation_tissue"
  entry <- scaled_unc$entry_means[scaled_unc$entry_means$evaluation_id == axis &
                                 scaled_unc$entry_means$method_key == "ECODA_authors_HR", ]
  stopifnot(grepl(sprintf("Scaled mean (five metrics): %.3f", entry$mean_score),
                 scaled_panels[[panel]][[1]]$labels$subtitle, fixed = TRUE))
}
stopifnot(scaled_unc$entry_means$mean_score[scaled_unc$entry_means$evaluation_id == "condition" &
          scaled_unc$entry_means$method_key == "ECODA_authors_HR"] == 1)
scaled_unc$metadata$constant <- "one level"
constant_panels <- plot_batch_release_mds(scaled_unc, "constant", labels, varying_scorer)
stopifnot(identical(constant_panels$batch_constant[[1]]$labels$subtitle,
                    "Technical labels: unscorable (one level)"))
corrected_labels <- setNames(batch_release_method_keys("batch_effect_corrected"), batch_release_method_keys("batch_effect_corrected"))
pairs <- plot_batch_release_comparison(unc, corr, right_control$datasets$Kidney_KPMP_full$batch_col, corrected_labels, scorer)
stopifnot(length(pairs$biological) == 16L,
          identical(pairs$biological[[3]]$data[c("Dim1", "Dim2")], pairs$biological[[5]]$data[c("Dim1", "Dim2")]))
# New method and metadata absence never invokes a historical manifest fallback.
bad <- control; bad$ledger <- bad$ledger[-1, ]; saveRDS(bad, control_path)
fail(read_batch_release_control(control_path, md5(control_path), "batch_effect_uncorrected"))
bad <- control; bad$ledger <- bad$ledger[!(bad$ledger$dataset == "Kidney_KPMP_full" & bad$ledger$method_key == "__sample_metadata__"), ]; saveRDS(bad, control_path)
fail(read_batch_release_control(control_path, md5(control_path), "batch_effect_uncorrected"))
bad <- control; bad$expected_scope$require_new[1] <- FALSE; saveRDS(bad, control_path)
fail(read_batch_release_control(control_path, md5(control_path), "batch_effect_uncorrected"))
bad <- right_control; bad$ledger$method_key[bad$ledger$method_key == "ecoda_unsupervised_harmony"] <- "ECODA_seuratres_2"
fail(load_batch_release_dataset(bad, "Kidney_KPMP_full", scorer))
# Create-only derived snapshots preserve validated old bytes on disagreement.
output <- file.path(root, "derived", "scores.rds")
write_batch_release_artifact(left, output); original <- md5(output)
write_batch_release_artifact(left, output)
fail(write_batch_release_artifact(right, output)); stopifnot(md5(output) == original)
orphan <- file.path(root, "derived", "orphan.rds"); writeLines("preserved receipt", paste0(orphan, ".md5"))
orphan_md5 <- md5(paste0(orphan, ".md5")); fail(write_batch_release_artifact(left, orphan))
stopifnot(!file.exists(orphan), md5(paste0(orphan, ".md5")) == orphan_md5)
# Parse all maintained notebook R chunks; never render or generate production figures.
for (path in c("notebooks/batch_effect_analysis_uncorrected.rmd", "notebooks/batch_effect_analysis_corrected.rmd")) {
  text <- readLines(path); inside <- FALSE; code <- character()
  for (line in text) {
    if (grepl("^```\\{r", line)) { inside <- TRUE; next }
    if (grepl("^```", line)) { inside <- FALSE; next }
    if (inside) code <- c(code, line)
  }
  parse(text = code)
}
# Exercise the actual shared five-metric/no-runtime builder only on this tiny
# synthetic release. Never render production notebook figures in a test.
for (report in list(left, right)) {
  frame <- report$scores; frame$method <- frame$method_key; frame$score <- frame$metric
  methods <- unique(frame$method_key); datasets <- unique(frame$dataset)
  grid <- frame[c("dataset", "view", "evaluation_id", "method_key", "metric")]
  prepared <- prepare_benchmark_data(frame, NULL, methods, setNames(methods, methods),
    setNames(datasets, datasets), evaluation_metrics(), expected_grid = grid)
  stopifnot(length(unique(prepared$score)) == 5L, length(unique(prepared$dataset)) == 9L,
            all(prepared$scaled_value == .5))
  filename <- paste0(unique(frame$view), ".pdf")
  plot <- build_funky_heatmap(frame, NULL, methods, setNames(methods, methods),
    setNames(datasets, datasets), evaluation_metrics(), output_file = filename,
    output_dir = root, expected_grid = grid)
  stopifnot(file.exists(file.path(root, filename)), file.info(file.path(root, filename))$size > 0)
}
stopifnot(all(vapply(control$ledger$artifact_path, md5, character(1)) == control$ledger$md5))
unlink(root, recursive = TRUE)
cat("current batch release report contracts: OK\n")
