#!/usr/bin/env Rscript

# Exercise the shared notebook panel builder without production artifacts.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
setwd(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
suppressPackageStartupMessages(library(ggplot2))
source("src/utils/plotting.R")
source("src/utils/batch_effect_analysis.R")
source("src/utils/scoring_metrics.R")

selected_methods <- c("first", "second")
method_label_map_main <- c(first = "First", second = "Second")
metadata <- data.frame(
  Sample = paste0("sample_", 1:6),
  assay = c("A", "A", "B", "B", "A", NA),
  site = c("X", "Y", "X", "Y", "X", "Y"),
  biology = c("case", "case", "control", "control", "case", "control"),
  row.names = paste0("sample_", 1:6)
)
features <- cbind(c(0, 1, 3, 2, 5, 4), c(0, 2, 1, 4, 3, 6))
rownames(features) <- rownames(metadata)
distance <- dist(features)
result <- list(
  dataset = "Synthetic", metadata = metadata, label_col = "biology",
  methods = list(
    first = list(dist_mat = distance, scores = list(anosim_score = 0.42, cluster_score = 0.12, mod_knn3_score = 0.72)),
    second = list(dist_mat = distance * 2)
  )
)
panels <- plot_batch_mds_panels(result, list("assay", "site"), method_label_map_main[selected_methods])
displayed_scores <- function(plot) {
  values <- sub(".*: ", "", strsplit(plot$labels$subtitle, "\n", fixed = TRUE)[[1L]])
  values[values == "NA"] <- NA_character_
  as.numeric(values)
}
expected_scores <- function(distance, labels) {
  keep <- !is.na(labels)
  distance <- as.dist(as.matrix(distance)[keep, keep])
  labels <- factor(labels[keep])
  components <- c(
    vegan::anosim(distance, labels, permutations = 0)$statistic,
    clust_eval(distance, labels),
    calc_modularity(distance, labels, 3)
  )
  round(c(components[[1L]], mean(components)), 3)
}
for (method_index in seq_along(selected_methods)) {
  assay <- panels$batch_assay[[method_index]]
  site <- panels$batch_site[[method_index]]
  combined <- panels$batch_combined[[method_index]]
  stopifnot(
    identical(as.character(assay$data$Cluster), metadata$assay),
    identical(as.character(site$data$Cluster), metadata$site),
    identical(as.character(combined$data$Cluster), c("A | X", "A | Y", "B | X", "B | Y", "A | X", NA_character_)),
    identical(assay$data[, c("Dim1", "Dim2")], site$data[, c("Dim1", "Dim2")]),
    identical(assay$data[, c("Dim1", "Dim2")], combined$data[, c("Dim1", "Dim2")])
  )
  # Incomplete metadata must not silently drop a sample or reuse biological scores.
  built <- ggplot_build(assay)
  stopifnot(all(is.finite(built$data[[1L]]$x)), all(is.finite(built$data[[1L]]$y)),
            all(!is.na(built$data[[1L]]$colour)), all(!is.na(built$data[[1L]]$shape)))
  distance <- result$methods[[method_index]]$dist_mat
  stopifnot(
    isTRUE(all.equal(displayed_scores(assay), expected_scores(distance, metadata$assay))),
    isTRUE(all.equal(displayed_scores(site), expected_scores(distance, metadata$site))),
    isTRUE(all.equal(displayed_scores(combined), expected_scores(distance, as.character(combined$data$Cluster))))
  )
}
# Biological panels retain persisted components; other panels use their own labels.
stopifnot(identical(displayed_scores(panels$biological[[1L]]), c(0.42, 0.42)),
          isTRUE(all.equal(displayed_scores(panels$biological[[2L]]),
                          expected_scores(result$methods$second$dist_mat, metadata$biology))))
single <- plot_batch_mds_panels(result, "site", method_label_map_main[selected_methods])
stopifnot(identical(as.character(single$batch_site[[1L]]$data$Cluster), metadata$site),
          is.null(single$batch_combined))
result$metadata$constant <- "one batch"
constant <- plot_batch_mds_panels(result, "constant", method_label_map_main[selected_methods])
stopifnot(all(is.na(displayed_scores(constant$batch_constant[[1L]]))))
missing_error <- tryCatch(plot_batch_mds_panels(result, "missing", method_label_map_main[selected_methods]), error = conditionMessage)
stopifnot(identical(missing_error, "Synthetic: metadata lacks configured batch columns: missing"))

# Reject incompatible identities instead of intersecting or padding the lanes.
expect_identity_rejection <- function(left, right) {
  error <- tryCatch({
    plot_batch_mds_comparison(left, right, c("assay", "site"), method_label_map_main[selected_methods])
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"))
}
shorter <- result
shorter$metadata <- metadata[1:5, , drop = FALSE]
shorter$methods <- lapply(result$methods, function(method) {
  list(dist_mat = as.dist(as.matrix(method$dist_mat)[1:5, 1:5]))
})
expect_identity_rejection(result, shorter)
reordered <- result
reordered$metadata <- metadata[6:1, , drop = FALSE]
reordered$methods <- lapply(result$methods, function(method) {
  list(dist_mat = as.dist(as.matrix(method$dist_mat)[6:1, 6:1]))
})
expect_identity_rejection(result, reordered)
different <- result
different$metadata$Sample[[6L]] <- "other_sample"
expect_identity_rejection(result, different)
different$dataset <- "Other"
different$metadata <- metadata
expect_identity_rejection(result, different)
for (invalid_dataset in list(NULL, character(), NA_character_, "", c("Synthetic", "Other"))) {
  invalid <- result
  invalid$dataset <- invalid_dataset
  expect_identity_rejection(invalid, invalid)
}
for (invalid_samples in list(NULL, rep("", 6), rep(" ", 6),
                             c(NA_character_, metadata$Sample[-1L]),
                             rep("duplicate", 6))) {
  invalid <- result
  invalid$metadata$Sample <- invalid_samples
  expect_identity_rejection(result, invalid)
  expect_identity_rejection(invalid, result)
}
empty <- result
empty$metadata <- metadata[FALSE, , drop = FALSE]
expect_identity_rejection(empty, empty)

# Valid pairing retains every sample and each lane's distances/scores, with shared colours.
corrected <- result
corrected$metadata$assay[[2L]] <- "C"
corrected_features <- features
corrected_features[, 2L] <- corrected_features[, 2L] * 3
corrected$methods <- list(
  first = list(dist_mat = dist(corrected_features),
               scores = list(anosim_score = 0.21, cluster_score = 0.31, mod_knn3_score = 0.41)),
  second = list(dist_mat = dist(corrected_features) * 2)
)
right_panels <- plot_batch_mds_panels(corrected, c("assay", "site"), method_label_map_main[selected_methods])
comparison <- plot_batch_mds_comparison(result, corrected, c("assay", "site"), method_label_map_main[selected_methods])
for (grouping in names(comparison)) {
  for (index in seq_along(selected_methods)) {
    left <- comparison[[grouping]][[2L * index - 1L]]
    right <- comparison[[grouping]][[2L * index]]
    stopifnot(
      nrow(left$data) == nrow(metadata),
      nrow(right$data) == nrow(metadata),
      identical(left$data, panels[[grouping]][[index]]$data),
      identical(right$data, right_panels[[grouping]][[index]]$data),
      identical(displayed_scores(left), displayed_scores(panels[[grouping]][[index]])),
      identical(displayed_scores(right), displayed_scores(right_panels[[grouping]][[index]]))
    )
  }
}
stopifnot(
  !isTRUE(all.equal(panels$biological[[1L]]$data[, c("Dim1", "Dim2")],
                    right_panels$biological[[1L]]$data[, c("Dim1", "Dim2")])),
  identical(displayed_scores(comparison$biological[[1L]]), c(0.42, 0.42)),
  identical(displayed_scores(comparison$biological[[2L]]), c(0.21, 0.31))
)
for (index in seq_along(selected_methods)) {
  left <- comparison$batch_assay[[2L * index - 1L]]
  right <- comparison$batch_assay[[2L * index]]
  stopifnot(
    identical(left$data, panels$batch_assay[[index]]$data),
    identical(right$data, right_panels$batch_assay[[index]]$data),
    identical(displayed_scores(left), displayed_scores(panels$batch_assay[[index]])),
    identical(displayed_scores(right), displayed_scores(right_panels$batch_assay[[index]]))
  )
  left_built <- ggplot_build(left)
  right_built <- ggplot_build(right)
  stopifnot(
    identical(left_built$data[[1L]]$colour[c(1L, 3L)], right_built$data[[1L]]$colour[c(1L, 3L)]),
    all(c("A", "B", "C") %in% left_built$plot$scales$get_scales("colour")$get_limits()),
    identical(left_built$plot$scales$get_scales("colour")$get_limits(),
              right_built$plot$scales$get_scales("colour")$get_limits())
  )
}
cat("PASS: separate/combined labels, complete-case scores, rejected identity mismatches, paired lane identity, and shared comparison colours\n")
