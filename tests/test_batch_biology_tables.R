#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
setwd(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
source("src/utils/batch_effect_analysis.R")
metadata <- data.frame(
  Sample = paste0("sample_", 1:6),
  biology = c("normal", "normal", "cancer", "normal", "cancer", NA),
  batch = c("A", "A", "B", "B", "B", NA),
  assay = c("X", "Y", "X", "X", "X", "Y")
)
tables <- batch_biology_tables(metadata, "biology", c("batch", "assay"))
counts <- tables$batch_batch
stopifnot(
  counts["A", "normal"] == 2L, counts["A", "cancer"] == 0L,
  counts["B", "normal"] == 1L, counts["B", "cancer"] == 2L,
  counts["<NA>", "<NA>"] == 1L,
  all(vapply(tables, sum, integer(1)) == nrow(metadata)),
  tables$batch_combined["A | X", "normal"] == 1L,
  tables$batch_combined["A | Y", "normal"] == 1L,
  tables$batch_combined["<NA> | Y", "<NA>"] == 1L
)
root <- tempfile("bio-batch-contract-")
write_batch_biology_tables(metadata, "Synthetic", "biology", c("batch", "assay"), root, root, "test")
composition <- read.csv(file.path(root, "Synthetic_bio_vs_batch_batch_composition.csv"))
stopifnot(
  composition$pct_within_batch[composition$batch == "A" & composition$biology == "normal"] == 100,
  isTRUE(all.equal(composition$pct_within_batch[composition$batch == "B" & composition$biology == "cancer"], 200 / 3)),
  isTRUE(all.equal(composition$pct_within_biology[composition$batch == "A" & composition$biology == "normal"], 200 / 3))
)
unlink(root, recursive = TRUE)
error <- tryCatch(batch_biology_tables(metadata[c(1L, 1L), ], "biology", "batch"), error = conditionMessage)
stopifnot(is.character(error), grepl("unique", error, fixed = TRUE))
error <- tryCatch(batch_biology_tables(metadata, "biology", "absent"), error = conditionMessage)
stopifnot(is.character(error), grepl("configured", error, fixed = TRUE))

# Kidney anatomy is diagnostic only: retain unknown regions and do not add it
# to the correction keys or their combined-batch table.
root <- tempfile("kidney-region-contract-")
dir.create(file.path(root, "metadata"), recursive = TRUE)
dir.create(file.path(root, "results"))
metadata$region.l1 <- c("Cortex", "Cortex", "Medulla", "Papilla", "NA", NA)
config <- list(sample_col = "specimen", label_col = "biology", batch_col = "batch")
for (dataset in c("Kidney_KPMP_full", "Other")) {
  metadata_path <- file.path(root, "metadata", paste0(dataset, "_sample_metadata.feather"))
  arrow::write_feather(metadata, metadata_path)
  .batch_write_checksum(metadata_path)
  summary_path <- file.path(root, "results", paste0(dataset, "_batch_effect_corrected_final_metadata.rds"))
  saveRDS(list(labels = stats::setNames(metadata$biology, metadata$Sample)), summary_path)
  .batch_write_checksum(summary_path)
}
provenance <- write_batch_biology_report(
  root, normalizePath(root),
  list(Kidney_KPMP_full = config, Other = config),
  c("Kidney_KPMP_full", "Other"), file.path(root, "plots"), file.path(root, "tables")
)
region_counts <- read.csv(file.path(root, "tables", "Kidney_KPMP_full_bio_vs_batch_region.l1_counts.csv"),
                          check.names = FALSE, na.strings = "")
stopifnot(
  sum(as.matrix(region_counts[, -1L])) == 6L,
  region_counts$normal[region_counts$batch_level == "Cortex"] == 2L,
  region_counts$cancer[region_counts$batch_level == "Papilla"] == 0L,
  region_counts$cancer[region_counts$batch_level == "<NA>"] == 1L,
  region_counts[region_counts$batch_level == "<NA>", "<NA>"] == 1L,
  all(provenance$batch_columns == "batch"),
  !file.exists(file.path(root, "tables", "Kidney_KPMP_full_bio_vs_batch_combined_counts.csv")),
  !file.exists(file.path(root, "tables", "Other_bio_vs_batch_region.l1_counts.csv"))
)
metadata$region.l1 <- NULL
metadata_path <- file.path(root, "metadata", "Kidney_KPMP_full_sample_metadata.feather")
arrow::write_feather(metadata, metadata_path)
.batch_write_checksum(metadata_path)
error <- tryCatch(write_batch_biology_report(
  root, normalizePath(root), list(Kidney_KPMP_full = config), "Kidney_KPMP_full",
  file.path(root, "plots"), file.path(root, "tables")
), error = conditionMessage)
stopifnot(is.character(error), grepl("configured label/batch columns", error, fixed = TRUE))
unlink(root, recursive = TRUE)
cat("PASS: Kidney region counts retain missing/zero combinations, preserve correction keys, and reject absent anatomy\n")
cat("PASS: sample counts, absent combinations, missing labels, combined batches, and duplicate rejection\n")
