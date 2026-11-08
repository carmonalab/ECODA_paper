#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
setwd(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
source("src/utils/datasets_io.R")
datasets <- read_datasets_json()
stopifnot(
  identical(datasets$Kidney_KPMP_full$label_col, "condition.l1"),
  identical(datasets$Kidney_KPMP_full$batch_col, "suspension_type"),
  identical(datasets$Kidney_KPMP_full$evaluation_labels,
            list(condition = "condition.l1", tissue = "tissue")),
  identical(datasets$Adams$evaluation_labels, list(primary = "Disease_Identity")),
  identical(datasets$Breast_cancer$batch_metadata_policy$sample_aggregation, "majority_v1"),
  identical(datasets$Breast_cancer$batch_metadata_policy$majority_keys,
            list("suspension_dissociation_time")),
  identical(datasets$Kidney_KPMP_full$views$batch_effect_corrected$columns,
            datasets$Kidney_KPMP_full$views$batch_effect_uncorrected$columns)
)
fixture <- jsonlite::fromJSON("datasets.json", simplifyVector = FALSE)
fixture$Kidney_KPMP_full$columns$evaluation_labels <- list(tissue = list("invalid"))
path <- tempfile(fileext = ".json")
jsonlite::write_json(fixture, path, auto_unbox = TRUE)
error <- tryCatch(read_datasets_json(path), error = conditionMessage)
unlink(path)
stopifnot(is.character(error), grepl("evaluation_labels", error, fixed = TRUE))
cat("PASS: dataset-level evaluation labels, primary compatibility, and technical policy exposure\n")
