#!/usr/bin/env Rscript

# Deterministic synthetic contract test for the uncorrected batch-effect
# analysis helpers. No production cohort, scheduler state, or datasets.json
# mutation is allowed here.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Cannot determine repository root")
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
setwd(repo_root)

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
  library(vegan)
})
source("src/utils/batch_effect_analysis.R")
source("src/utils/plotting.R")

expect_error <- function(expr, pattern = NULL) {
  error_text <- NULL
  tryCatch(
    eval(substitute(expr), envir = parent.frame()),
    error = function(error) error_text <<- conditionMessage(error)
  )
  if (is.null(error_text)) stop("expected an error")
  if (!is.null(pattern) && !grepl(pattern, error_text, fixed = TRUE)) {
    stop("error did not contain expected text: ", pattern, "\n", error_text)
  }
  invisible(error_text)
}

set.seed(42)
root <- tempfile("batch-analysis-contract-")
input_root <- file.path(root, "uncorrected")
metadata_root <- file.path(input_root, "metadata")
dir.create(file.path(input_root, "results"), recursive = TRUE)
dir.create(file.path(input_root, "embeddings"), recursive = TRUE)
dir.create(metadata_root, recursive = TRUE)

sample_ids <- paste0("sample_", seq_len(8L))
labels <- rep(c("mild", "severe"), each = 4L)
metadata <- data.frame(
  Sample = sample_ids,
  `CoVID-19 severity` = labels,
  `secondary biology` = rep(c("A", "B"), 4L),
  technical = rep(c("batch_1", "batch_2"), each = 2L, times = 2L),
  technical_duplicate = rep(c("dup_1", "dup_2"), each = 2L, times = 2L),
  near_unique = c("n1", "n2", "n3", "n4", "n5", "n1", "n2", "n3"),
  confounded = labels,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(metadata) <- sample_ids
registry <- data.frame(
  dataset = "Synthetic",
  candidate = c("CoVID-19 severity", "secondary biology", "technical"),
  candidate_class = c("biological_primary", "biological_secondary", "technical"),
  is_primary = c(TRUE, FALSE, FALSE),
  is_secondary_biology = c(FALSE, TRUE, FALSE),
  is_technical = c(FALSE, FALSE, TRUE),
  formula_alias = paste0("batch_term_", 1:3),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

summary_path <- file.path(input_root, "results", "Synthetic_batch_effect_uncorrected_metadata.rds")
summary <- list(
  labels = structure(factor(labels, levels = c("mild", "severe")), names = sample_ids),
  n_cells = 800L,
  n_samples = length(sample_ids),
  cells_per_sample = structure(rep(100L, length(sample_ids)), names = sample_ids),
  n_cell_types_high_res = 4L
)
saveRDS(summary, summary_path)
.batch_write_checksum(summary_path)

make_r_bundle <- function(path, key, matrix) {
  saveRDS(setNames(list(list(dist_mat = stats::as.dist(matrix))), key), path)
  .batch_write_checksum(path)
}

make_feather <- function(path, matrix) {
  frame <- as.data.frame(matrix, check.names = FALSE)
  frame$Sample <- rownames(matrix)
  arrow::write_feather(frame, path)
  .batch_write_checksum(path)
}

base_matrix <- as.matrix(dist(matrix(c(
  0, 1, 0, 1, 0, 1, 0, 1,
  1, 0, 1, 0, 1, 0, 1, 0,
  0, 1, 0, 1, 0, 1, 0, 1
), nrow = 8L, byrow = TRUE)))
# Use a full-rank deterministic sample distance matrix rather than a feature
# matrix with repeated rows.
base_matrix <- as.matrix(dist(matrix(seq_len(24L), nrow = 8L, ncol = 3L)))
rownames(base_matrix) <- sample_ids
colnames(base_matrix) <- sample_ids
composition_path <- file.path(input_root, "results", "Synthetic_batch_effect_uncorrected_composition.rds")
saveRDS(
  list(
    ECODA_authors_HR = list(dist_mat = stats::as.dist(base_matrix)),
    ECODA_seuratres_2 = list(dist_mat = stats::as.dist(base_matrix * 1.1))
  ),
  composition_path
)
.batch_write_checksum(composition_path)
make_r_bundle(
  file.path(input_root, "results", "Synthetic_batch_effect_uncorrected_pseudobulk.rds"),
  "Pseudobulk_hvg2000",
  base_matrix * 1.2
)
make_r_bundle(
  file.path(input_root, "results", "Synthetic_batch_effect_uncorrected_gloscope.rds"),
  "GloScope_hvg2000_pcadims30",
  base_matrix * 1.3
)
for (suffix in c("mrvi", "pilot", "qot")) {
  make_feather(
    file.path(input_root, "embeddings", paste0("Synthetic_batch_effect_uncorrected_hvg2000_highres_", suffix, "_dists.feather")),
    base_matrix * (1 + match(suffix, c("mrvi", "pilot", "qot")) / 10)
  )
}
write_batch_metadata_sidecar(
  metadata,
  file.path(metadata_root, "Synthetic_sample_metadata.feather"),
  expected_sample_ids = sample_ids
)

specs <- batch_uncorrected_method_specs()
stopifnot(identical(names(specs), c(
  "ECODA_authors_HR", "ECODA_seuratres_2", "Pseudobulk_hvg2000",
  "GloScope_hvg2000_pcadims30", "MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000"
)))
loaded <- load_batch_uncorrected_dataset(input_root, metadata_root, "Synthetic", registry)
stopifnot(
  identical(names(loaded$methods), names(specs)),
  identical(loaded$sample_ids, sample_ids),
  identical(rownames(loaded$methods[["MrVI_hvg2000"]]$matrix), sample_ids)
)

expect_error(
  read_batch_method_dist(input_root, "Synthetic", specs[["MrVI_hvg2000"]], rev(sample_ids)),
  "sample IDs/order mismatch"
)

# Strict sidecars: missing, stale, and corrected/legacy fallback names all fail.
metadata_sidecar <- file.path(metadata_root, "Synthetic_sample_metadata.feather")
file.remove(paste0(metadata_sidecar, ".md5"))
expect_error(read_batch_metadata_sidecar(metadata_sidecar, sample_ids, "CoVID-19 severity"), "missing checksum sidecar")
write_batch_metadata_sidecar(metadata, metadata_sidecar, expected_sample_ids = sample_ids)
writeLines(c(
  "MD5=00000000000000000000000000000000",
  paste0("SIZE=", file.info(metadata_sidecar)$size),
  paste0("PATH=", normalizePath(metadata_sidecar))
), paste0(metadata_sidecar, ".md5"))
expect_error(read_batch_metadata_sidecar(metadata_sidecar, sample_ids, "CoVID-19 severity"), "MD5 mismatch")
write_batch_metadata_sidecar(metadata, metadata_sidecar, expected_sample_ids = sample_ids)
unlink(file.path(input_root, "embeddings", "Synthetic_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"))
file.copy(
  file.path(input_root, "embeddings", "Synthetic_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"),
  file.path(input_root, "embeddings", "Synthetic_batch_effect_corrected_hvg2000_highres_mrvi_dists.feather")
)
expect_error(load_batch_uncorrected_dataset(input_root, metadata_root, "Synthetic", registry), "missing")

# Recreate the selected artifact for the remaining behavioral checks.
make_feather(
  file.path(input_root, "embeddings", "Synthetic_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather"),
  base_matrix * 1.1
)

mds <- plot_mds(
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  factor(labels),
  title = "Synthetic",
  cluster_score = FALSE,
  mod_score = FALSE,
  sil_score = FALSE,
  anosim_score = TRUE
)
stopifnot(grepl("ANOSIM score:", mds$labels$title, fixed = TRUE))

anosim <- compute_batch_anosim(loaded$methods[["ECODA_authors_HR"]]$dist_mat, metadata[["CoVID-19 severity"]], permutations = 9L)
univariate <- compute_batch_permanova(loaded$methods[["ECODA_authors_HR"]]$dist_mat, metadata, "technical", permutations = 9L)
informed <- compute_batch_permanova(loaded$methods[["ECODA_authors_HR"]]$dist_mat, metadata, "technical", "CoVID-19 severity", permutations = 9L)
stopifnot(
  is.finite(anosim$statistic),
  is.finite(univariate$r2),
  is.finite(informed$r2)
)

metric_table <- make_batch_metric_table(
  "Synthetic", "ECODA_authors_HR", metadata, registry,
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  loaded$methods[["ECODA_authors_HR"]]$path,
  permutations = 9L
)
stopifnot(
  nrow(metric_table) == nrow(registry),
  all(metric_table$biology_informed_status[metric_table$is_primary] == "NOT_A_CANDIDATE"),
  any(is.finite(metric_table$univariate_r2))
)

joint <- compute_batch_joint_permanova(
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  metadata,
  registry,
  permutations = 9L
)
joint_table <- make_batch_joint_table("Synthetic", "ECODA_authors_HR", joint, loaded$methods[["ECODA_authors_HR"]]$path)
stopifnot(
  all(registry$candidate %in% joint$rows$candidate),
  any(joint$rows$status == "ESTIMABLE"),
  all(is.finite(joint$rows$p_adjusted_holm[joint$rows$status == "ESTIMABLE"])),
  identical(joint$decomposition$component, c("Unique Biological", "Unique Technical", "Shared / Confounded", "Residual / Unexplained")),
  abs(sum(joint$decomposition$r2) - 1) < 1e-8,
  nrow(joint_table) == nrow(registry)
)

aliased_registry <- registry
aliased_registry$candidate[3L] <- "confounded"
aliased_registry$candidate_class[3L] <- "technical"
aliased_registry$is_technical[3L] <- TRUE
aliased <- compute_batch_joint_permanova(
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  metadata,
  aliased_registry,
  permutations = 9L
)
stopifnot(
  aliased$rows$status[aliased$rows$candidate == "confounded"] == "DROPPED_CONFOUNDED_PRIMARY",
  identical(aliased$status, "ESTIMABLE_REDUCED"),
  !is.null(aliased$decomposition),
  abs(sum(aliased$decomposition$r2) - 1) < 1e-8
)

duplicate_registry <- rbind(registry, data.frame(
  dataset = "Synthetic",
  candidate = "technical_duplicate",
  candidate_class = "technical",
  is_primary = FALSE,
  is_secondary_biology = FALSE,
  is_technical = TRUE,
  formula_alias = "batch_term_4",
  check.names = FALSE,
  stringsAsFactors = FALSE
))
duplicate <- compute_batch_joint_permanova(
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  metadata,
  duplicate_registry,
  permutations = 9L
)
stopifnot(
  any(duplicate$rows$joint_model_term == "technical__technical_duplicate"),
  duplicate$rows$status[duplicate$rows$candidate == "technical_duplicate"] == "MERGED_EXACT_PARTITION"
)

near_unique_registry <- rbind(registry, data.frame(
  dataset = "Synthetic",
  candidate = "near_unique",
  candidate_class = "technical",
  is_primary = FALSE,
  is_secondary_biology = FALSE,
  is_technical = TRUE,
  formula_alias = "batch_term_4",
  check.names = FALSE,
  stringsAsFactors = FALSE
))
near_unique <- compute_batch_joint_permanova(
  loaded$methods[["ECODA_authors_HR"]]$dist_mat,
  metadata,
  near_unique_registry,
  permutations = 9L
)
stopifnot(
  near_unique$rows$status[near_unique$rows$candidate == "near_unique"] == "DROPPED_NEAR_UNIQUE"
)

nested_metadata <- data.frame(
  Sample = sample_ids,
  disease = labels,
  origin = rep(c("A", "A", "B", "B"), 2L),
  origin_fine = rep(c("A1", "A2", "B1", "B2"), 2L),
  site = rep(c("site1", "site2", "site1", "site2"), 2L),
  dissociation_protocol = rep(c("p1", "p2", "p2", "p1"), 2L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(nested_metadata) <- sample_ids
nested_registry <- data.frame(
  dataset = "Lung",
  candidate = c("disease", "origin", "origin_fine", "site", "dissociation_protocol"),
  candidate_class = c("biological_primary", "biological_secondary", "biological_secondary", "technical", "technical"),
  is_primary = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  is_secondary_biology = c(FALSE, TRUE, TRUE, FALSE, FALSE),
  is_technical = c(FALSE, FALSE, FALSE, TRUE, TRUE),
  formula_alias = paste0("batch_term_", 1:5),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
nested <- prepare_batch_joint_design(nested_metadata, nested_registry)
stopifnot(
  nested$rows$joint_status[nested$rows$candidate == "origin_fine"] == "DROPPED_NESTED",
  nested$rows$joint_status[nested$rows$candidate == "site"] == "RETAINED",
  nested$rows$joint_status[nested$rows$candidate == "dissociation_protocol"] == "RETAINED"
)

nmi_registry <- rbind(registry, data.frame(
  dataset = "Synthetic",
  candidate = "confounded",
  candidate_class = "technical",
  is_primary = FALSE,
  is_secondary_biology = FALSE,
  is_technical = TRUE,
  formula_alias = "batch_term_4",
  check.names = FALSE,
  stringsAsFactors = FALSE
))
nmi <- compute_batch_nmi(metadata, nmi_registry)
nmi_table <- make_batch_nmi_table("Synthetic", nmi, nmi_registry)
stopifnot(
  isTRUE(all.equal(nmi, t(nmi), check.attributes = TRUE)),
  all(diag(nmi)[seq_len(nrow(nmi) - 1L)] == 1),
  nmi["CoVID-19 severity", "confounded"] > 0.70,
  nrow(nmi_table) == nrow(nmi) * ncol(nmi),
  any(nmi_table$severe_collinearity)
)

cat("test_batch_effect_analysis.R: all checks passed\n")
