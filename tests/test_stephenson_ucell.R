# Run from repository root: Rscript tests/test_stephenson_ucell.R
source("src/notebooks/stephenson_ucell.R")
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
ids <- paste0("S", 1:4)
mat <- cbind(typeA = 1:4, typeB = c(2, 4, 1, 3),
             rare = c(NA, 1, 2, 3), constant = rep(1, 4))
rownames(mat) <- ids
selected <- prepare_stephenson_ucell_features(mat)
stopifnot(identical(colnames(selected$features), c("typeA", "typeB")),
          identical(rownames(selected$features), ids),
          identical(selected$dropped_incomplete, "rare"),
          identical(selected$dropped_constant, "constant"))
labels <- factor(c("Covid", "Healthy", "Covid", "Healthy"))
names(labels) <- ids
bundle <- list(feat_mat = mat[, 1:2, drop = FALSE][c(4, 2, 1, 3), ],
               labels = labels[c(3, 4, 2, 1)])
aligned <- align_stephenson_bundle(bundle, ids, labels)
stopifnot(identical(aligned$features, mat[, 1:2, drop = FALSE]))
bad <- bundle
rownames(bad$feat_mat)[1] <- "unexpected"
expect_error(align_stephenson_bundle(bad, ids, labels), "sample IDs")
bad <- bundle
bad$feat_mat <- bad$feat_mat[-1, , drop = FALSE]
expect_error(align_stephenson_bundle(bad, ids, labels), "sample IDs")
bad <- bundle
bad$labels[1] <- if (bad$labels[1] == "Covid") "Healthy" else "Covid"
expect_error(align_stephenson_bundle(bad, ids, labels), "labels disagree")

# Exercise the cached-feature comparison with independently reordered rows and labels.
library(ggplot2)
source("src/utils/plotting.R")
source("src/utils/scoring_metrics.R")
pseudobulk_mat <- cbind(IFI27 = 1:4, RPS4Y1 = c(2, 4, 1, 3),
                       ISG15 = c(4, 1, 3, 2), constant = rep(1, 4),
                       IGHG1 = c(20, 1, 15, 2), IGLC2 = c(1, 30, 2, 20))
rownames(pseudobulk_mat) <- ids
pseudobulk <- list(feat_mat = pseudobulk_mat, labels = labels)
comparison_artifacts <- list(samples = data.frame(Sample = ids), labels = labels)
reference <- compare_stephenson_pseudobulk(comparison_artifacts, pseudobulk, n_loadings = 10L)
pseudobulk$feat_mat <- pseudobulk$feat_mat[c(4, 2, 1, 3), , drop = FALSE]
pseudobulk$labels <- pseudobulk$labels[c(3, 4, 2, 1)]
reordered <- compare_stephenson_pseudobulk(comparison_artifacts, pseudobulk, n_loadings = 20L)
stopifnot(
  identical(reordered$features, pseudobulk_mat[, 1:3, drop = FALSE]),
  identical(reordered$labels, labels),
  identical(reordered$table, reference$table),
  reordered$table$n_samples == length(ids), reordered$table$n_features == 3L,
  reordered$table$n_excluded_ig == 2L,
  identical(reordered$excluded_genes, c("IGHG1", "IGLC2")),
  identical(pseudobulk$feat_mat, pseudobulk_mat[c(4, 2, 1, 3), , drop = FALSE])
)
all_ig <- pseudobulk
all_ig$feat_mat <- all_ig$feat_mat[, c("IGHG1", "IGLC2"), drop = FALSE]
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, all_ig),
             "fewer than two")
bad <- pseudobulk
rownames(bad$feat_mat)[1] <- "unexpected"
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, bad), "sample IDs")
bad <- pseudobulk
bad$feat_mat <- bad$feat_mat[-1, , drop = FALSE]
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, bad), "sample IDs")
bad <- pseudobulk
bad$labels <- bad$labels[-1]
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, bad), "sample IDs")
bad <- pseudobulk
names(bad$labels)[1] <- "unexpected"
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, bad), "sample IDs")
bad <- pseudobulk
bad$labels[1] <- if (bad$labels[1] == "Covid") "Healthy" else "Covid"
expect_error(compare_stephenson_pseudobulk(comparison_artifacts, bad), "labels disagree")
for (limit in list(0L, -1L, 1.5, NA_real_, Inf, "15", TRUE, numeric(), c(8L, 15L))) {
  expect_error(compare_stephenson_ucell(NULL, NULL, n_loadings = limit), "positive integer")
  expect_error(compare_stephenson_pseudobulk(NULL, NULL, n_loadings = limit), "positive integer")
}

# Build a complete compact artifact so checksum rejection exercises the loader.
root <- tempfile("stephenson_ucell_test_")
dir.create(root)
files <- list()
write_table <- function(table, name) {
  path <- file.path(root, name)
  arrow::write_feather(table, path)
  files[[name]] <<- digest::digest(path, algo = "sha256", file = TRUE)
}
for (annotation in stephenson_ucell_annotations) {
  for (score in stephenson_ucell_scores) {
    write_table(data.frame(Sample = ids, typeA = 1:4, typeB = c(2, 4, 1, 3)),
                paste0(annotation, "__", score, ".feather"))
  }
  write_table(data.frame(Sample = ids, typeA = rep(1L, 4), typeB = rep(1L, 4)),
              paste0(annotation, "__n_cells.feather"))
}
write_table(data.frame(Sample = ids, Status = as.character(labels)), "samples.feather")
manifest <- list(
  schema_version = 1L, dataset = "Stephenson", view = "benchmark_analysis",
  annotations = as.list(c(authors_lowres = "initial_clustering", authors_highres = "full_clustering",
                          hitme_layer1 = "layer1", hitme_layer2 = "layer2")),
  scores = stephenson_ucell_scores, sample_ids = ids, label_column = "Status", n_cells = 8L,
  excluded_annotation_cells = as.list(setNames(rep(0L, 4), stephenson_ucell_annotations)),
  cell_type_columns = setNames(rep(list(c("typeA", "typeB")), 4), stephenson_ucell_annotations),
  files = files
)
jsonlite::write_json(manifest, file.path(root, "manifest.json"), auto_unbox = TRUE)
artifacts <- load_stephenson_ucell(root)
stopifnot(length(artifacts$matrices) == 16L, nrow(artifacts$coverage) == 16L,
          identical(names(artifacts$labels), ids))
corrupted <- file.path(root, "authors_lowres__IFN_UCell.feather")
connection <- file(corrupted, "ab")
writeBin(as.raw(1L), connection)
close(connection)
expect_error(load_stephenson_ucell(root), "checksum mismatch")
unlink(root, recursive = TRUE)
cat("Stephenson UCell artifact/feature/alignment checks passed\n")

# Real compact CD14 export: fixed 55-sample/2000-gene contract, including IG genes.
cd14_root <- tempfile("stephenson_cd14_test_")
dir.create(cd14_root)
cd14_ids <- sprintf("S%02d", seq_len(55L))
cd14_labels <- c(rep("Covid", 43L), rep("Healthy", 12L))
cd14_samples <- data.frame(Sample = cd14_ids, Status = cd14_labels,
                           n_cells = seq_len(55L) + 5L)
cd14_matrix <- outer(seq_len(55L), seq_len(2000L), function(i, j) sin(i * j / 29))
colnames(cd14_matrix) <- c("IGHG1", paste0("gene", seq_len(1999L)))
rownames(cd14_matrix) <- cd14_ids
cd14_expression <- data.frame(Sample = cd14_ids, cd14_matrix, check.names = FALSE)
cd14_manifest <- list(
  schema_version = 1L, dataset = "Stephenson", view = "benchmark_analysis",
  annotation_column = "initial_clustering", cell_type = "CD14",
  sample_ids = cd14_ids, n_selected_cells = sum(cd14_samples$n_cells), n_features = 2000L,
  normalization = list(blind = TRUE, batch_col = NULL, correct_batch = FALSE,
                       n_hvg = 2000L, black_list = "none"),
  files = list()
)
write_cd14_table <- function(table, name) {
  path <- file.path(cd14_root, name)
  arrow::write_feather(table, path)
  cd14_manifest$files[[name]] <<- digest::digest(path, algo = "sha256", file = TRUE)
}
write_cd14_manifest <- function(value = cd14_manifest) {
  jsonlite::write_json(value, file.path(cd14_root, "manifest.json"),
                       auto_unbox = TRUE, null = "null")
}
write_cd14_table(cd14_expression, "expression.feather")
write_cd14_table(cd14_samples, "samples.feather")
write_cd14_manifest()
cd14 <- load_stephenson_cd14_pseudobulk(cd14_root)
stopifnot(identical(cd14$feat_mat, cd14_matrix),
          identical(as.character(cd14$labels), cd14_labels),
          identical(names(cd14$labels), cd14_ids),
          identical(cd14$metadata, cd14_samples),
          "IGHG1" %in% colnames(cd14$feat_mat))
cd14_reference <- list(feat_mat = cd14_matrix[55:1, , drop = FALSE],
                       labels = cd14$labels[55:1])
stopifnot(identical(align_stephenson_bundle(
  cd14_reference, cd14_ids, cd14$labels
)$features, cd14_matrix))
bad_cd14 <- cd14_reference
bad_cd14$labels[1L] <- "Covid"
expect_error(align_stephenson_bundle(bad_cd14, cd14_ids, cd14$labels), "labels disagree")
bad_cd14 <- cd14_reference
rownames(bad_cd14$feat_mat)[1L] <- "unexpected"
expect_error(align_stephenson_bundle(bad_cd14, cd14_ids, cd14$labels), "sample IDs")
for (field in c("cell_type", "annotation_column", "normalization")) {
  invalid <- cd14_manifest
  if (field == "normalization") {
    invalid$normalization$black_list <- "IG"
  } else {
    invalid[[field]] <- if (field == "cell_type") "CD14_mono" else "full_clustering"
  }
  write_cd14_manifest(invalid)
  expect_error(load_stephenson_cd14_pseudobulk(cd14_root), "identity/normalization")
}
write_cd14_manifest()
connection <- file(file.path(cd14_root, "expression.feather"), "ab")
writeBin(as.raw(1L), connection)
close(connection)
expect_error(load_stephenson_cd14_pseudobulk(cd14_root), "checksum mismatch")
write_cd14_table(cd14_expression[55:1, ], "expression.feather")
write_cd14_manifest()
expect_error(load_stephenson_cd14_pseudobulk(cd14_root), "axes")
write_cd14_table(cd14_expression, "expression.feather")
invalid_samples <- cd14_samples
invalid_samples$Status[1L] <- "unknown"
write_cd14_table(invalid_samples, "samples.feather")
write_cd14_manifest()
expect_error(load_stephenson_cd14_pseudobulk(cd14_root), "sample labels")
invalid_samples <- cd14_samples
invalid_samples$n_cells[1L] <- 4L
write_cd14_table(invalid_samples, "samples.feather")
write_cd14_manifest()
expect_error(load_stephenson_cd14_pseudobulk(cd14_root), "cell counts")
write_cd14_table(cd14_samples, "samples.feather")
write_cd14_manifest()
cd14_swapped <- cd14_samples
cd14_swapped$Status <- rev(cd14_swapped$Status)
write_cd14_table(cd14_swapped, "samples.feather")
write_cd14_manifest()
cd14_relabeled <- load_stephenson_cd14_pseudobulk(cd14_root)
stopifnot(identical(cd14_relabeled$feat_mat, cd14$feat_mat))
expect_error(align_stephenson_bundle(
  cd14_reference, cd14_ids, cd14_relabeled$labels
), "labels disagree")
unlink(cd14_root, recursive = TRUE)
cat("Stephenson CD14 compact-export checks passed\n")
