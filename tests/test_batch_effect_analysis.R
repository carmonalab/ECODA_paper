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

final_datasets <- c(
  "Alzheimer", "Breast_cancer", "Covid19_PBMC", "Kidney_KPMP_full",
  "Diabetes", "Lupus_PBMC", "Lung", "Joanito", "Stephenson"
)
final_methods <- c(names(specs), "ECODA_authors_HR_NULL")
final_analysis_root <- file.path(root, "uncorrected_final")
final_results_root <- file.path(final_analysis_root, "results")
final_embeddings_root <- file.path(final_analysis_root, "embeddings")
final_metadata_root <- file.path(final_analysis_root, "metadata")
dir.create(final_results_root, recursive = TRUE)
dir.create(final_embeddings_root, recursive = TRUE)
dir.create(final_metadata_root, recursive = TRUE)

final_summary_path <- file.path(
  final_results_root,
  "Synthetic_batch_effect_uncorrected_final_metadata.rds"
)
saveRDS(summary, final_summary_path)
.batch_write_checksum(final_summary_path)
final_metadata_path <- file.path(final_metadata_root, "Synthetic_sample_metadata.feather")
write_batch_metadata_sidecar(
  metadata,
  final_metadata_path,
  expected_sample_ids = sample_ids
)

final_combo <- function(scale) {
  list(
    scores = list(sil_score = 0.5),
    dist_mat = stats::as.dist(base_matrix * scale)
  )
}
final_composition_path <- file.path(
  final_results_root,
  "Synthetic_batch_effect_uncorrected_final_composition.rds"
)
saveRDS(
  list(
    ECODA_authors_HR = final_combo(1.0),
    ECODA_seuratres_2 = final_combo(1.1),
    ECODA_authors_HR_NULL = final_combo(1.2)
  ),
  final_composition_path
)
.batch_write_checksum(final_composition_path)
final_pseudobulk_path <- file.path(
  final_results_root,
  "Synthetic_batch_effect_uncorrected_final_pseudobulk.rds"
)
saveRDS(list(Pseudobulk_hvg2000 = final_combo(1.3)), final_pseudobulk_path)
.batch_write_checksum(final_pseudobulk_path)
final_gloscope_path <- file.path(
  final_results_root,
  "Synthetic_batch_effect_uncorrected_final_gloscope.rds"
)
saveRDS(list(GloScope_hvg2000_pcadims30 = final_combo(1.4)), final_gloscope_path)
.batch_write_checksum(final_gloscope_path)
for (suffix in c("mrvi", "pilot", "qot")) {
  make_feather(
    file.path(
      final_embeddings_root,
      paste0(
        "Synthetic_batch_effect_uncorrected_final_hvg2000_highres_",
        suffix, "_dists.feather"
      )
    ),
    base_matrix * (1 + match(suffix, c("mrvi", "pilot", "qot")) / 20)
  )
}

manifest_write <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, mustWork = TRUE)
}
manifest_row <- function(values) paste(values, collapse = "\t")
final_metadata_manifest <- manifest_write(
  file.path(root, "final-analysis-metadata.tsv"),
  c(
    "dataset\tsummary_lane\tfeather_lane\tmetadata_summary_path\tmetadata_feather_path",
    manifest_row(c(
      "Synthetic", "final", "final",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_metadata.rds",
      "uncorrected_final/metadata/Synthetic_sample_metadata.feather"
    ))
  )
)
final_artifact_manifest <- manifest_write(
  file.path(root, "final-analysis-artifacts.tsv"),
  c(
    "dataset\tmethod\tlane\tartifact_kind\tartifact_path\tbundle_key",
    manifest_row(c(
      "Synthetic", "ECODA_authors_HR", "final", "rds_bundle",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_composition.rds",
      "ECODA_authors_HR"
    )),
    manifest_row(c(
      "Synthetic", "ECODA_seuratres_2", "final", "rds_bundle",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_composition.rds",
      "ECODA_seuratres_2"
    )),
    manifest_row(c(
      "Synthetic", "Pseudobulk_hvg2000", "final", "rds_bundle",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_pseudobulk.rds",
      "Pseudobulk_hvg2000"
    )),
    manifest_row(c(
      "Synthetic", "GloScope_hvg2000_pcadims30", "final", "rds_bundle",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_gloscope.rds",
      "GloScope_hvg2000_pcadims30"
    )),
    manifest_row(c(
      "Synthetic", "MrVI_hvg2000", "final", "distance_feather",
      "uncorrected_final/embeddings/Synthetic_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather",
      ""
    )),
    manifest_row(c(
      "Synthetic", "PILOT_hvg2000", "final", "distance_feather",
      "uncorrected_final/embeddings/Synthetic_batch_effect_uncorrected_final_hvg2000_highres_pilot_dists.feather",
      ""
    )),
    manifest_row(c(
      "Synthetic", "QOT_hvg2000", "final", "distance_feather",
      "uncorrected_final/embeddings/Synthetic_batch_effect_uncorrected_final_hvg2000_highres_qot_dists.feather",
      ""
    )),
    manifest_row(c(
      "Synthetic", "ECODA_authors_HR_NULL", "final", "rds_bundle",
      "uncorrected_final/results/Synthetic_batch_effect_uncorrected_final_composition.rds",
      "ECODA_authors_HR_NULL"
    ))
  )
)

final_metadata_frame <- read_batch_final_manifest(
  final_metadata_manifest,
  expected_datasets = "Synthetic",
  expected_methods = final_methods,
  repository_root = root
)
final_artifact_frame <- read_batch_final_manifest(
  final_artifact_manifest,
  expected_datasets = "Synthetic",
  expected_methods = final_methods,
  repository_root = root
)
stopifnot(
  identical(attr(final_metadata_frame, "manifest_kind"), "metadata"),
  identical(attr(final_artifact_frame, "manifest_kind"), "artifacts"),
  all(vapply(final_metadata_frame, is.character, logical(1L))),
  all(vapply(final_artifact_frame, is.character, logical(1L))),
  identical(
    final_artifact_frame$bundle_key[
      final_artifact_frame$artifact_kind == "distance_feather"
    ],
    rep("", 3L)
  )
)
composition_rows <- final_artifact_frame[
  final_artifact_frame$method %in% c(
    "ECODA_authors_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL"
  ),
  ,
  drop = FALSE
]
stopifnot(
  length(unique(composition_rows$artifact_path)) == 1L,
  identical(
    composition_rows$bundle_key,
    c("ECODA_authors_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
  )
)
loaded_final <- load_batch_uncorrected_dataset_from_manifest(
  final_metadata_frame,
  final_artifact_frame,
  "Synthetic",
  registry,
  root
)
stopifnot(
  identical(names(loaded_final$methods), final_methods),
  identical(loaded_final$sample_ids, sample_ids),
  identical(
    loaded_final$methods[["Pseudobulk_hvg2000"]]$path,
    normalizePath(final_pseudobulk_path, mustWork = TRUE)
  ),
  identical(
    loaded_final$methods[["ECODA_authors_HR_NULL"]]$bundle_key,
    "ECODA_authors_HR_NULL"
  ),
  identical(
    loaded_final$methods[["ECODA_authors_HR_NULL"]]$sample_ids,
    sample_ids
  )
)

# A legacy null score artifact is aggregate-only and intentionally has no
# sample-ID axis.  The explicit standalone_scores mapping permits it while
# retaining strict final RDS/Feather checks above.
legacy_null_path <- file.path(
  input_root,
  "results",
  "Synthetic_batch_effect_uncorrected_ECODA_authors_HR_NULL.rds"
)
saveRDS(
  list(scores = list(
    anosim_score = 0.1,
    mod_knn3_score = 0.2,
    cluster_score = 0.3
  )),
  legacy_null_path
)
legacy_metadata_manifest <- manifest_write(
  file.path(root, "legacy-analysis-metadata.tsv"),
  c(
    "dataset\tsummary_lane\tfeather_lane\tmetadata_summary_path\tmetadata_feather_path",
    manifest_row(c(
      "LegacySynthetic", "legacy", "legacy",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_metadata.rds",
      "uncorrected/metadata/Synthetic_sample_metadata.feather"
    ))
  )
)
legacy_artifact_manifest <- manifest_write(
  file.path(root, "legacy-analysis-artifacts.tsv"),
  c(
    "dataset\tmethod\tlane\tartifact_kind\tartifact_path\tbundle_key",
    manifest_row(c(
      "LegacySynthetic", "ECODA_authors_HR", "legacy", "rds_bundle",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_composition.rds",
      "ECODA_authors_HR"
    )),
    manifest_row(c(
      "LegacySynthetic", "ECODA_seuratres_2", "legacy", "rds_bundle",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_composition.rds",
      "ECODA_seuratres_2"
    )),
    manifest_row(c(
      "LegacySynthetic", "Pseudobulk_hvg2000", "legacy", "rds_bundle",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_pseudobulk.rds",
      "Pseudobulk_hvg2000"
    )),
    manifest_row(c(
      "LegacySynthetic", "GloScope_hvg2000_pcadims30", "legacy", "rds_bundle",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_gloscope.rds",
      "GloScope_hvg2000_pcadims30"
    )),
    manifest_row(c(
      "LegacySynthetic", "MrVI_hvg2000", "legacy", "distance_feather",
      "uncorrected/embeddings/Synthetic_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather",
      ""
    )),
    manifest_row(c(
      "LegacySynthetic", "PILOT_hvg2000", "legacy", "distance_feather",
      "uncorrected/embeddings/Synthetic_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather",
      ""
    )),
    manifest_row(c(
      "LegacySynthetic", "QOT_hvg2000", "legacy", "distance_feather",
      "uncorrected/embeddings/Synthetic_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather",
      ""
    )),
    manifest_row(c(
      "LegacySynthetic", "ECODA_authors_HR_NULL", "legacy", "standalone_scores",
      "uncorrected/results/Synthetic_batch_effect_uncorrected_ECODA_authors_HR_NULL.rds",
      "scores"
    ))
  )
)
legacy_metadata_frame <- read_batch_final_manifest(
  legacy_metadata_manifest,
  expected_datasets = "LegacySynthetic",
  expected_methods = final_methods,
  repository_root = root
)
legacy_artifact_frame <- read_batch_final_manifest(
  legacy_artifact_manifest,
  expected_datasets = "LegacySynthetic",
  expected_methods = final_methods,
  repository_root = root
)
legacy_registry <- registry
legacy_registry$dataset <- "LegacySynthetic"
loaded_legacy <- load_batch_uncorrected_dataset_from_manifest(
  legacy_metadata_frame,
  legacy_artifact_frame,
  "LegacySynthetic",
  legacy_registry,
  root
)
stopifnot(
  identical(
    loaded_legacy$methods[["ECODA_authors_HR_NULL"]]$artifact_kind,
    "standalone_scores"
  ),
  is.null(loaded_legacy$methods[["ECODA_authors_HR_NULL"]]$sample_ids),
  identical(
    names(loaded_legacy$methods[["ECODA_authors_HR_NULL"]]$scores),
    c("anosim_score", "mod_knn3_score", "cluster_score")
  )
)

# Final RDS and distance-Feather bundles must carry the ordered sample IDs.
bad_rds_path <- file.path(final_results_root, "bad-composition.rds")
bad_rds_combo <- final_combo(1.0)
bad_rds_combo$dist_mat <- stats::as.dist(unname(as.matrix(bad_rds_combo$dist_mat)))
saveRDS(
  setNames(
    list(bad_rds_combo, bad_rds_combo, bad_rds_combo),
    c("ECODA_authors_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
  ),
  bad_rds_path
)
.batch_write_checksum(bad_rds_path)
bad_rds_lines <- readLines(final_artifact_manifest, warn = FALSE)
for (index in c(2L, 3L, 9L)) {
  bad_rds_parts <- strsplit(bad_rds_lines[[index]], "\t", fixed = TRUE)[[1L]]
  bad_rds_parts[[5L]] <- "uncorrected_final/results/bad-composition.rds"
  bad_rds_lines[[index]] <- manifest_row(bad_rds_parts)
}
bad_rds_manifest <- manifest_write(
  file.path(root, "bad-rds-artifacts.tsv"),
  bad_rds_lines
)
bad_rds_frame <- read_batch_final_manifest(
  bad_rds_manifest,
  expected_datasets = "Synthetic",
  expected_methods = final_methods,
  repository_root = root
)
expect_error(
  load_batch_uncorrected_dataset_from_manifest(
    final_metadata_frame,
    bad_rds_frame,
    "Synthetic",
    registry,
    root
  ),
  "sample IDs/order mismatch"
)

bad_feather_path <- file.path(final_embeddings_root, "bad-mrvi.feather")
arrow::write_feather(as.data.frame(base_matrix), bad_feather_path)
.batch_write_checksum(bad_feather_path)
bad_feather_lines <- readLines(final_artifact_manifest, warn = FALSE)
bad_feather_parts <- strsplit(bad_feather_lines[[6L]], "\t", fixed = TRUE)[[1L]]
bad_feather_parts[[5L]] <- "uncorrected_final/embeddings/bad-mrvi.feather"
bad_feather_lines[[6L]] <- paste0(manifest_row(bad_feather_parts), "\t")
bad_feather_manifest <- manifest_write(
  file.path(root, "bad-feather-artifacts.tsv"),
  bad_feather_lines
)
bad_feather_frame <- read_batch_final_manifest(
  bad_feather_manifest,
  expected_datasets = "Synthetic",
  expected_methods = final_methods,
  repository_root = root
)
expect_error(
  load_batch_uncorrected_dataset_from_manifest(
    final_metadata_frame,
    bad_feather_frame,
    "Synthetic",
    registry,
    root
  ),
  "distance Feather schema is not square"
)

# The physical sixth TSV field is mandatory even when its value is empty.
bad_key_lines <- readLines(final_artifact_manifest, warn = FALSE)
bad_key_parts <- strsplit(bad_key_lines[[8L]], "\t", fixed = TRUE)[[1L]]
bad_key_parts[[6L]] <- "not-empty"
bad_key_lines[[8L]] <- manifest_row(bad_key_parts)
expect_error(
  read_batch_final_manifest(
    manifest_write(file.path(root, "bad-key-artifacts.tsv"), bad_key_lines),
    expected_datasets = "Synthetic",
    expected_methods = final_methods,
    repository_root = root
  ),
  "empty bundle_key"
)
bad_physical_lines <- readLines(final_artifact_manifest, warn = FALSE)
bad_physical_fields <- strsplit(bad_physical_lines[[7L]], "\t", fixed = TRUE)[[1L]]
bad_physical_lines[[7L]] <- paste(bad_physical_fields[-length(bad_physical_fields)], collapse = "\t")
expect_error(
  read_batch_final_manifest(
    manifest_write(file.path(root, "bad-physical-artifacts.tsv"), bad_physical_lines),
    expected_datasets = "Synthetic",
    expected_methods = final_methods,
    repository_root = root
  ),
  "wrong physical field count"
)
bad_path_lines <- readLines(final_metadata_manifest, warn = FALSE)
bad_path_parts <- strsplit(bad_path_lines[[2L]], "\t", fixed = TRUE)[[1L]]
bad_path_parts[[5L]] <- "../outside.feather"
bad_path_lines[[2L]] <- manifest_row(bad_path_parts)
expect_error(
  read_batch_final_manifest(
    manifest_write(file.path(root, "bad-path-metadata.tsv"), bad_path_lines),
    expected_datasets = "Synthetic",
    expected_methods = final_methods,
    repository_root = root
  ),
  "escapes repository root"
)

# The production final manifest has an exact mixed-source nine-dataset order.
scope_metadata_lines <- c(
  "dataset\tsummary_lane\tfeather_lane\tmetadata_summary_path\tmetadata_feather_path",
  vapply(
    final_datasets,
    function(dataset) {
      frozen <- dataset %in% c(
        "Alzheimer", "Breast_cancer", "Lupus_PBMC", "Stephenson"
      )
      summary_lane <- if (frozen || dataset == "Kidney_KPMP_full") {
        "legacy"
      } else {
        "final"
      }
      feather_lane <- if (frozen) "legacy" else "final"
      manifest_row(c(
        dataset,
        summary_lane,
        feather_lane,
        paste0("scope/", dataset, "/metadata_summary.rds"),
        paste0("scope/", dataset, "/metadata.feather")
      ))
    },
    character(1L)
  )
)
scope_artifact_lines <- c(
  "dataset\tmethod\tlane\tartifact_kind\tartifact_path\tbundle_key",
  unlist(lapply(final_datasets, function(dataset) {
    lane <- if (dataset %in% c("Alzheimer", "Breast_cancer", "Lupus_PBMC", "Stephenson")) "legacy" else "final"
    vapply(final_methods, function(method) {
      kind <- if (method %in% c("MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000")) {
        "distance_feather"
      } else if (method == "ECODA_authors_HR_NULL" && lane == "legacy") {
        "standalone_scores"
      } else {
        "rds_bundle"
      }
      key <- if (kind == "distance_feather") "" else if (kind == "standalone_scores") "scores" else method
      path <- if (
        lane == "final" &&
        method %in% c("ECODA_authors_HR", "ECODA_seuratres_2", "ECODA_authors_HR_NULL")
      ) {
        paste0("scope/", dataset, "/composition.rds")
      } else {
        paste0(
          "scope/", dataset, "/", method,
          if (kind == "distance_feather") ".feather" else ".rds"
        )
      }
      manifest_row(c(dataset, method, lane, kind, path, key))
    }, character(1L))
  }), use.names = FALSE)
)
scope_metadata_manifest <- manifest_write(
  file.path(root, "scope-metadata.tsv"),
  scope_metadata_lines
)
scope_artifact_manifest <- manifest_write(
  file.path(root, "scope-artifacts.tsv"),
  scope_artifact_lines
)
scope_metadata_frame <- read_batch_final_manifest(
  scope_metadata_manifest,
  expected_datasets = final_datasets,
  expected_methods = final_methods,
  repository_root = root
)
scope_artifact_frame <- read_batch_final_manifest(
  scope_artifact_manifest,
  expected_datasets = final_datasets,
  expected_methods = final_methods,
  repository_root = root
)
stopifnot(
  identical(scope_metadata_frame$dataset, final_datasets),
  identical(
    scope_metadata_frame$summary_lane,
    c("legacy", "legacy", "final", "legacy", "final", "legacy", "final", "final", "legacy")
  ),
  identical(
    scope_metadata_frame$feather_lane,
    c("legacy", "legacy", "final", "final", "final", "legacy", "final", "final", "legacy")
  ),
  identical(
    scope_artifact_frame$dataset,
    rep(final_datasets, each = length(final_methods))
  ),
  identical(
    scope_artifact_frame$method,
    rep(final_methods, times = length(final_datasets))
  )
)

cat("test_batch_effect_analysis.R: all checks passed\n")

