#!/usr/bin/env Rscript
# Build all 6-cohort publication batch figures (excluding Alzheimer) and derived
# batch tables under data/derived/batch_6cohort_release_20261126/:
#   - Figure 5A: plots/Figure_5_A_batch_6cohorts_funkyheatmap.pdf (combined 16-row
#     uncorrected + corrected + Shuffled_baseline funkyheatmap with joint min-max
#     scaling per (dataset, metric), bold corrected / muted-gray uncorrected row
#     labels, "Corrected" X column, "Correction" text column, and sLISI label)
#   - Figure 5B: plots/Figure_5_B_batch_correction_arrows_ANOSIM_6cohorts.pdf
#     (2x3 ANOSIM arrows for the 6 cohorts) + 5-metric companions
#   - Composite Figure 5: plots/Figure_5.pdf (+ plots/Figure_5.png)
#   - SF26: plots/Supp_fig_TMP_ECODA_becm_comp_ANOSIM.{pdf,png}
#     (2x3 |ΔCLR| >= 0..4 trajectories + limma for 6 cohorts; droplow50pct hidden)
#   - BatchMDS (SF27-SF32):
#     1) Internal (with 5-metric score subtitles):
#        plots/batch_mds_internal_with_scores_20261126/BatchMDS_<Dataset>_with_scores.pdf
#     2) Publication (no score subtitles, larger legends, short <Cohort> header):
#        plots/Supp_fig_TMP_BatchMDS_<Dataset>.pdf (6 retained cohorts only)

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".", mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(cowplot)
  library(ggpubr)
  library(funkyheatmap)
  library(arrow)
  library(grid)
  library(tibble)
})

source("src/utils/load_all_functions.R")
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/batch_effect_analysis.R")

set.seed(42)
options(mc.cores = 1L)

retained_cohorts <- c(
  "Breast_cancer",
  "Covid19_PBMC",
  "Joanito",
  "Kidney_KPMP_full",
  "Lupus_PBMC",
  "Stephenson"
)

cohort_display_map <- c(
  Breast_cancer = "Breast cancer",
  Covid19_PBMC = "COVID-19 PBMC",
  Joanito = "Joanito",
  Kidney_KPMP_full = "Kidney (KPMP)",
  Lupus_PBMC = "Lupus PBMC",
  Stephenson = "Stephenson"
)

cohort_tech_col <- c(
  Breast_cancer = "suspension_dissociation_time",
  Covid19_PBMC = "datasets",
  Joanito = "seqtec",
  Kidney_KPMP_full = "suspension_type",
  Lupus_PBMC = "batch_cov",
  Stephenson = "Site"
)

cohort_bio_col <- c(
  Breast_cancer = "disease",
  Covid19_PBMC = "CoVID-19 severity",
  Joanito = "sample.origin",
  Kidney_KPMP_full = "condition.l1",
  Lupus_PBMC = "Status",
  Stephenson = "Status"
)

# Local metric labels using sLISI (never modifying src/utils/constants.R)
batch_score_label_map <- c(
  anosim_score = "ANOSIM",
  mod_knn3_score = "Modularity",
  cluster_score = "ARI",
  sil_score = "Silhouette",
  lisi_score = "sLISI"
)

clean_short_batch <- function(dataset, x) {
  x <- as.character(x)
  if (dataset == "Breast_cancer") {
    return(ifelse(x == "12 hour", "12 h", ifelse(x == "15 minute", "15 min", x)))
  }
  if (dataset == "Lupus_PBMC") {
    x <- sub("^dmx_", "", x)
    x <- sub("^YS-JY-", "YS-", x)
    x <- sub("_pool", " p", x)
    return(x)
  }
  x
}

check_local_sidecar <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  sc <- paste0(path, ".md5")
  stopifnot(file.exists(sc))
  fields <- list()
  for (ln in readLines(sc, warn = FALSE)) {
    p <- strsplit(ln, "=", fixed = TRUE)[[1L]]
    if (length(p) >= 2L) fields[[p[[1L]]]] <- paste(p[-1L], collapse = "=")
  }
  stopifnot(
    identical(tolower(fields$MD5), tolower(unname(tools::md5sum(path)))),
    identical(as.character(fields$SIZE), as.character(file.info(path)$size)),
    identical(basename(fields$PATH), basename(path))
  )
  path
}

read_dist_feather_strict <- function(path, sample_ids) {
  check_local_sidecar(path)
  tbl <- as.data.frame(arrow::read_feather(path), check.names = FALSE)
  stopifnot("__index_level_0__" %in% names(tbl))
  rids <- as.character(tbl[["__index_level_0__"]])
  tbl[["__index_level_0__"]] <- NULL
  cids <- as.character(colnames(tbl))
  mat <- as.matrix(tbl)
  storage.mode(mat) <- "double"
  rownames(mat) <- rids
  colnames(mat) <- cids
  stopifnot(
    identical(rids, cids),
    setequal(rids, sample_ids),
    length(rids) == length(sample_ids),
    all(is.finite(mat)),
    all(mat >= 0),
    max(abs(diag(mat))) <= 1e-10,
    max(abs(mat - t(mat))) <= 1e-10
  )
  stats::as.dist(mat[sample_ids, sample_ids, drop = FALSE])
}

read_rds_bundle_strict <- function(path, key, sample_ids) {
  check_local_sidecar(path)
  obj <- readRDS(path)
  bundle <- if (!is.null(obj[[key]]) && is.list(obj[[key]])) obj[[key]] else obj
  stopifnot(is.list(bundle), !is.null(bundle$dist_mat))
  dmat <- as.matrix(bundle$dist_mat)
  ids <- rownames(dmat)
  stopifnot(
    !is.null(ids),
    identical(ids, colnames(dmat)),
    setequal(ids, sample_ids),
    all(is.finite(dmat)),
    all(dmat >= 0),
    max(abs(diag(dmat))) <= 1e-10,
    max(abs(dmat - t(dmat))) <= 1e-10
  )
  bundle$dist_mat <- stats::as.dist(dmat[sample_ids, sample_ids, drop = FALSE])
  bundle
}

load_strict_lupus_view <- function(view_name, lupus_dir) {
  subdir <- if (view_name == "batch_effect_uncorrected") "uncorrected_final" else "corrected_final"
  base <- file.path(lupus_dir, subdir)
  meta_path <- check_local_sidecar(file.path(base, "metadata", "Lupus_PBMC_sample_metadata.feather"))
  meta <- as.data.frame(arrow::read_feather(meta_path), stringsAsFactors = FALSE)
  ids <- as.character(meta$Sample)
  stopifnot(
    length(ids) == 94L,
    !anyDuplicated(ids),
    all(c("Sample", "Status", "batch_cov") %in% names(meta))
  )
  rownames(meta) <- ids
  true_bio <- setNames(as.character(meta$Status), ids)

  methods_vec <- batch_release_method_keys(view_name)
  method_objs <- list()
  for (mk in methods_vec) {
    if (mk %in% c("MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000")) {
      short_py <- tolower(sub("_hvg2000$", "", mk))
      fpath <- file.path(base, "embeddings", paste0("Lupus_PBMC_", view_name, "_final_hvg2000_highres_", short_py, "_dists.feather"))
      d <- read_dist_feather_strict(fpath, ids)
      bundle <- list(dist_mat = d, labels = true_bio, scores = NULL)
    } else {
      rpath <- file.path(base, "results", paste0("Lupus_PBMC_", view_name, "_final_", mk, ".rds"))
      bundle <- read_rds_bundle_strict(rpath, mk, ids)
    }
    method_objs[[mk]] <- list(
      dist_mat = bundle$dist_mat,
      labels = true_bio,
      path = if (mk %in% c("MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000")) fpath else rpath
    )
  }

  list(
    dataset = "Lupus_PBMC",
    view = view_name,
    metadata = meta,
    sample_ids = ids,
    label_col = "Status",
    batch_col = "batch_cov",
    evaluations = list(primary = true_bio),
    primary_evaluation = "primary",
    methods = method_objs
  )
}

message("=== Step 1: Load 6-cohort inputs & compute CLR >= 2 drop 5-metric scores ===")
becm_cache_path <- check_local_sidecar("data/derived/ECODA_becm_comp_20261121_v2/comparison_cache.rds")
becm_cache <- readRDS(becm_cache_path)

metrics <- evaluation_metrics()
clr2_variants <- c(
  "ECODA_authors_sigdiff_ctdrop_2",
  "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_2"
)

clr2_bio_rows <- list()
clr2_detail_rows <- list()

for (ds in retained_cohorts) {
  entry <- becm_cache$datasets[[ds]]
  stopifnot(!is.null(entry), is.data.frame(entry$metadata))
  meta <- entry$metadata
  ids <- as.character(meta$Sample)
  stopifnot(
    identical(names(meta), c("Sample", "bio", "batch")),
    !anyNA(meta),
    !anyDuplicated(ids),
    length(unique(meta$bio)) >= 2L,
    length(unique(meta$batch)) >= 2L
  )

  for (vk in clr2_variants) {
    v <- entry$variants[[vk]]
    stopifnot(!is.null(v), identical(v$status, "PASS"), !is.null(v$dist_mat))
    dmat <- as.matrix(v$dist_mat)
    stopifnot(
      identical(rownames(dmat), ids),
      identical(colnames(dmat), ids),
      all(is.finite(dmat)),
      all(dmat >= 0),
      max(abs(diag(dmat))) <= 1e-10,
      max(abs(dmat - t(dmat))) <= 1e-10,
      min(dmat[upper.tri(dmat)]) > 0
    )
    dist_ordered <- stats::as.dist(dmat[ids, ids, drop = FALSE])

    sc_bio <- calc_sep_score(dist_ordered, meta$bio)
    sc_bat <- calc_sep_score(dist_ordered, meta$batch)
    vals_bio <- vapply(metrics, function(m) as.numeric(sc_bio[[m]]), numeric(1))
    vals_bat <- vapply(metrics, function(m) as.numeric(sc_bat[[m]]), numeric(1))
    stopifnot(all(is.finite(vals_bio)), all(is.finite(vals_bat)))

    # Cross-check ANOSIM against cached BECM score table
    cached_row <- becm_cache$scores[
      becm_cache$scores$dataset == ds & becm_cache$scores$method_id == vk,
      ,
      drop = FALSE
    ]
    stopifnot(
      nrow(cached_row) == 1L,
      abs(vals_bio[1] - cached_row$biological_anosim) < 1e-10,
      abs(vals_bat[1] - cached_row$batch_anosim) < 1e-10
    )

    clr2_bio_rows[[paste(ds, vk, sep = "::")]] <- data.frame(
      dataset = ds,
      view = "batch_effect_corrected",
      evaluation_id = "primary",
      method_key = vk,
      metric = metrics,
      value = unname(vals_bio),
      score_source = "becm_cache_clr2_calc_sep_score",
      stringsAsFactors = FALSE
    )

    clr2_detail_rows[[paste(ds, vk, sep = "::")]] <- data.frame(
      dataset = ds,
      method_id = vk,
      metric = metrics,
      biological_value = unname(vals_bio),
      batch_value = unname(vals_bat),
      n_initial = v$n_initial,
      n_significant = v$n_significant,
      n_dropped = v$n_dropped,
      n_remaining = v$n_remaining,
      stringsAsFactors = FALSE
    )
  }
}

clr2_bio_df <- do.call(rbind, clr2_bio_rows)
clr2_detail_df <- do.call(rbind, clr2_detail_rows)
rownames(clr2_bio_df) <- rownames(clr2_detail_df) <- NULL

# Load 7-cohort CSVs, filter out Alzheimer and ECODA limma methods
scores_7_unc <- read.csv(
  "data/derived/batch_7cohort_release_20261121/batch_7cohort_scores_uncorrected.csv",
  stringsAsFactors = FALSE
)
scores_7_cor <- read.csv(
  "data/derived/batch_7cohort_release_20261121/batch_7cohort_scores_corrected.csv",
  stringsAsFactors = FALSE
)
arrow_7 <- read.csv(
  "data/derived/batch_7cohort_release_20261121/batch_7cohort_arrow_scores.csv",
  stringsAsFactors = FALSE
)

scores_6_unc <- scores_7_unc[scores_7_unc$dataset %in% retained_cohorts, , drop = FALSE]
scores_6_cor_nonlimma <- scores_7_cor[
  scores_7_cor$dataset %in% retained_cohorts &
    !scores_7_cor$method_key %in% c("ECODA_authors_HR", "ecoda_unsupervised_harmony_limma", "ECODA_authors_HR_NULL"),
  ,
  drop = FALSE
]

# Define unified 16-row method table (7 uncorrected + 8 corrected + 1 Shuffled_baseline)
ctdrop_2_label <- "CT-drop (|\u0394CLR| \u2265 2)"
harmony_ctdrop_2_label <- "Harmony + CT-drop (|\u0394CLR| \u2265 2)"

row_spec <- tibble::tribble(
  ~row_id,                                     ~method_label,       ~state,        ~corrected_flag, ~correction_method,    ~source_view,
  "MrVI_hvg2000__corrected",                   "MrVI",              "corrected",   "X",             "MrVI",                "corrected",
  "GloScope_hvg2000_pcadims30__corrected",     "GloScope",          "corrected",   "X",             "Harmony",             "corrected",
  "Pseudobulk_hvg2000__corrected",             "Pseudobulk",        "corrected",   "X",             "limma",               "corrected",
  "QOT_hvg2000__corrected",                    "QOT",               "corrected",   "X",             "Harmony",             "corrected",
  "PILOT_hvg2000__corrected",                  "PILOT",             "corrected",   "X",             "Harmony",             "corrected",
  "ECODA_authors_sigdiff_ctdrop_2",            "ECODA_authors",     "corrected",   "X",             ctdrop_2_label,        "clr2",
  "ecoda_unsupervised_harmony__corrected",     "ECODA_Leiden",      "corrected",   "X",             "Harmony",             "corrected",
  "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_2", "ECODA_Leiden",      "corrected",   "X",             harmony_ctdrop_2_label,"clr2",
  "MrVI_hvg2000__uncorrected",                 "MrVI",              "uncorrected", "",              "",                    "uncorrected",
  "GloScope_hvg2000_pcadims30__uncorrected",   "GloScope",          "uncorrected", "",              "",                    "uncorrected",
  "Pseudobulk_hvg2000__uncorrected",           "Pseudobulk",        "uncorrected", "",              "",                    "uncorrected",
  "QOT_hvg2000__uncorrected",                  "QOT",               "uncorrected", "",              "",                    "uncorrected",
  "PILOT_hvg2000__uncorrected",                "PILOT",             "uncorrected", "",              "",                    "uncorrected",
  "ECODA_authors_HR__uncorrected",             "ECODA_authors",     "uncorrected", "",              "",                    "uncorrected",
  "ECODA_seuratres_2__uncorrected",            "ECODA_Leiden",      "uncorrected", "",              "",                    "uncorrected",
  "Shuffled_baseline",                         "Shuffled_baseline", "baseline",    "",              "",                    "uncorrected_null"
)

unc_nonnull <- scores_6_unc[scores_6_unc$method_key != "ECODA_authors_HR_NULL", , drop = FALSE]
unc_nonnull$row_id <- paste0(unc_nonnull$method_key, "__uncorrected")

unc_null <- scores_6_unc[scores_6_unc$method_key == "ECODA_authors_HR_NULL", , drop = FALSE]
unc_null$row_id <- "Shuffled_baseline"

cor_nonlimma <- scores_6_cor_nonlimma
cor_nonlimma$row_id <- paste0(cor_nonlimma$method_key, "__corrected")

clr2_part <- clr2_bio_df
clr2_part$row_id <- clr2_part$method_key

cols_keep <- c("dataset", "view", "evaluation_id", "method_key", "row_id", "metric", "value", "score_source")
combined_scores <- rbind(
  unc_nonnull[, cols_keep],
  cor_nonlimma[, cols_keep],
  clr2_part[, cols_keep],
  unc_null[, cols_keep]
) %>%
  left_join(row_spec, by = "row_id")

stopifnot(
  nrow(combined_scores) == 16L * 6L * 5L,
  !anyNA(combined_scores$method_label),
  !anyNA(combined_scores$value),
  all(is.finite(combined_scores$value)),
  !anyDuplicated(paste(combined_scores$dataset, combined_scores$row_id, combined_scores$metric, sep = "::"))
)

# Joint min-max scaling per (dataset, metric) across ALL 16 rows
combined_scores$scaled_value <- ave(
  combined_scores$value,
  interaction(combined_scores$dataset, combined_scores$metric, drop = TRUE),
  FUN = function(x) {
    r <- range(x)
    if (r[1] == r[2]) rep(0.5, length(x)) else (x - r[1]) / diff(r)
  }
)

# Build per-dataset, per-metric, and overall 6-cohort means & rank descending
ds_means <- combined_scores %>%
  group_by(row_id, method_label, state, corrected_flag, correction_method, dataset) %>%
  summarise(mean_score = mean(scaled_value), .groups = "drop")

metric_means <- combined_scores %>%
  group_by(row_id, metric) %>%
  summarise(mean_score = mean(scaled_value), .groups = "drop") %>%
  mutate(metric_label = unname(batch_score_label_map[metric]))

overall_means <- ds_means %>%
  group_by(row_id, method_label, state, corrected_flag, correction_method) %>%
  summarise(overall_mean = mean(mean_score), .groups = "drop") %>%
  mutate(rank = rank(-overall_mean, ties.method = "min")) %>%
  arrange(rank, method_label, correction_method)

ds_wide <- ds_means %>%
  mutate(col = paste0("Dataset__", dataset)) %>%
  select(row_id, col, mean_score) %>%
  pivot_wider(names_from = col, values_from = mean_score)

metric_wide <- metric_means %>%
  mutate(col = paste0("Metric__", metric_label)) %>%
  select(row_id, col, mean_score) %>%
  pivot_wider(names_from = col, values_from = mean_score)

ranking_combined <- overall_means %>%
  left_join(ds_wide, by = "row_id") %>%
  left_join(metric_wide, by = "row_id") %>%
  arrange(rank, desc(overall_mean))

# Build 6-cohort arrow table (8 non-NULL method pairs x 6 cohorts x 5 metrics = 240 rows)
arrow_nonlimma <- arrow_7[
  arrow_7$dataset %in% retained_cohorts &
    !arrow_7$method_key %in% c("ECODA_authors_HR", "ecoda_unsupervised_harmony_limma"),
  ,
  drop = FALSE
]
arrow_label_map <- c(
  ecoda_unsupervised_harmony = "ECODA_Leiden: Harmony",
  Pseudobulk_hvg2000 = "Pseudobulk",
  GloScope_hvg2000_pcadims30 = "GloScope",
  MrVI_hvg2000 = "MrVI",
  PILOT_hvg2000 = "PILOT",
  QOT_hvg2000 = "QOT"
)
arrow_nonlimma$method_label <- unname(arrow_label_map[arrow_nonlimma$method_key])
stopifnot(!anyNA(arrow_nonlimma$method_label))

# Build CT-drop (|ΔCLR| >= 2) arrow rows paired with respective uncorrected baselines
clr2_arrow_pairs <- tibble::tribble(
  ~method_key,                                 ~uncorrected_method_key, ~method_label,
  "ECODA_authors_sigdiff_ctdrop_2",            "ECODA_authors_HR",      "ECODA_authors",
  "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_2", "ECODA_seuratres_2",     "ECODA_Leiden: Harmony + CT-drop"
)

clr2_arrow_rows <- list()
for (i in seq_len(nrow(clr2_arrow_pairs))) {
  mk <- clr2_arrow_pairs$method_key[i]
  umk <- clr2_arrow_pairs$uncorrected_method_key[i]
  mlab <- clr2_arrow_pairs$method_label[i]
  ref_sub <- arrow_7[
    arrow_7$dataset %in% retained_cohorts &
      arrow_7$uncorrected_method_key == umk,
    c("dataset", "evaluation_id", "batch_column", "metric", "biological_uncorrected", "batch_uncorrected", "uncorrected_path"),
    drop = FALSE
  ]
  # Keep one unique uncorrected baseline row per (dataset, metric)
  ref_sub <- ref_sub[!duplicated(paste(ref_sub$dataset, ref_sub$metric, sep = "::")), , drop = FALSE]
  cor_sub <- clr2_detail_df[clr2_detail_df$method_id == mk, , drop = FALSE]
  merged <- merge(ref_sub, cor_sub, by = c("dataset", "metric"))
  stopifnot(nrow(merged) == 6L * 5L)
  merged$method_key <- mk
  merged$uncorrected_method_key <- umk
  merged$method_label <- mlab
  merged$biological_corrected <- merged$biological_value
  merged$batch_corrected <- merged$batch_value
  merged$corrected_path <- becm_cache_path
  clr2_arrow_rows[[mk]] <- merged[, names(arrow_nonlimma), drop = FALSE]
}

arrow_6 <- rbind(arrow_nonlimma, do.call(rbind, clr2_arrow_rows))
rownames(arrow_6) <- NULL
stopifnot(
  nrow(arrow_6) == 6L * 8L * 5L,
  setequal(unique(arrow_6$dataset), retained_cohorts),
  all(is.finite(as.matrix(arrow_6[, c("biological_uncorrected", "biological_corrected", "batch_uncorrected", "batch_corrected")])))
)

out_release_dir <- file.path("data", "derived", "batch_6cohort_release_20261126")
dir.create(out_release_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(ranking_combined, file.path(out_release_dir, "batch_6cohort_ranking_combined.csv"), row.names = FALSE)
write.csv(combined_scores, file.path(out_release_dir, "batch_6cohort_scores_combined.csv"), row.names = FALSE)
write.csv(arrow_6, file.path(out_release_dir, "batch_6cohort_arrow_scores.csv"), row.names = FALSE)
write.csv(clr2_detail_df, file.path(out_release_dir, "batch_6cohort_clr2_drop_scores.csv"), row.names = FALSE)
saveRDS(
  list(
    schema = "ecoda_batch_6cohort_release_20261126_v1",
    retained_cohorts = retained_cohorts,
    excluded = c("Alzheimer", "Lung", "Kidney_KPMP_full:tissue", "ECODA_authors_limma", "ECODA_Harmony_Leiden_limma", "droplow50pct"),
    ranking_combined = ranking_combined,
    scores_combined = combined_scores,
    arrow_scores = arrow_6,
    clr2_drop_scores = clr2_detail_df
  ),
  file.path(out_release_dir, "batch_6cohort_summary.rds")
)

cat("\n=== 6-Cohort Combined Ranking (Joint Min-Max Scaled Across All 16 Rows) ===\n")
print(
  ranking_combined[, c("rank", "method_label", "corrected_flag", "correction_method", "overall_mean",
                       paste0("Dataset__", retained_cohorts),
                       paste0("Metric__", unname(batch_score_label_map)))],
  row.names = FALSE,
  digits = 4
)

message("\n=== Step 2: Figure 5A (16-row Combined Ranked Funkyheatmap) ===")
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

build_batch_6cohort_funkyheatmap <- function(ranking_df, output_pdf, plot_width = 14.2, plot_height = 7.8) {
  ncolors <- 11L
  fh_data <- ranking_df %>%
    transmute(
      id = row_id,
      method = method_label,
      Corrected = corrected_flag,
      Correction = correction_method,
      Overall_Score = overall_mean,
      Dataset__Breast_cancer = Dataset__Breast_cancer,
      Dataset__Covid19_PBMC = Dataset__Covid19_PBMC,
      Dataset__Joanito = Dataset__Joanito,
      Dataset__Kidney_KPMP_full = Dataset__Kidney_KPMP_full,
      Dataset__Lupus_PBMC = Dataset__Lupus_PBMC,
      Dataset__Stephenson = Dataset__Stephenson,
      Metric__ANOSIM = Metric__ANOSIM,
      Metric__Modularity = Metric__Modularity,
      Metric__ARI = Metric__ARI,
      Metric__Silhouette = Metric__Silhouette,
      Metric__sLISI = Metric__sLISI
    )

  # Rename dataset columns to clean display names
  ds_display_order <- unname(cohort_display_map[retained_cohorts])
  for (ds in retained_cohorts) {
    colnames(fh_data)[colnames(fh_data) == paste0("Dataset__", ds)] <- paste0("Dataset__", cohort_display_map[[ds]])
  }

  info_method <- tibble::tribble(
    ~id,          ~group, ~name,        ~geom,  ~palette,           ~options,
    "method",     "",     "",           "text", NA_character_,      list(hjust = 0, width = 5.6),
    "Corrected",  "",     "Corrected",  "text", NA_character_,      list(hjust = 0.5, width = 2.0),
    "Correction", "",     "Correction", "text", NA_character_,      list(hjust = 0, width = 9.8)
  )
  info_overall <- tibble::tribble(
    ~id,             ~group,    ~name,   ~geom, ~palette,           ~options,
    "Overall_Score", "Overall", "Mean",  "bar", "palette_avg_mean", list(width = 4.0)
  )
  info_datasets <- tibble::tibble(
    id = paste0("Dataset__", sort(ds_display_order)),
    group = "Datasets",
    name = sort(ds_display_order),
    geom = "circle",
    palette = "Datasets",
    options = rep(list(list()), length(ds_display_order))
  )
  metric_order <- unname(batch_score_label_map[evaluation_metrics()])
  info_metrics <- tibble::tibble(
    id = paste0("Metric__", metric_order),
    group = "Metrics",
    name = metric_order,
    geom = "bar",
    palette = "Metrics",
    options = rep(list(list(width = 2.2)), length(metric_order))
  )

  column_info <- dplyr::bind_rows(info_method, info_overall, info_datasets, info_metrics)
  column_groups <- tibble::tibble(group = c("Overall", "Datasets", "Metrics"))

  palettes <- list(
    palette_avg_mean = RColorBrewer::brewer.pal(n = ncolors, name = "RdBu"),
    Datasets = RColorBrewer::brewer.pal(n = ncolors, name = "RdBu"),
    Metrics = RColorBrewer::brewer.pal(n = ncolors, name = "RdBu")
  )
  legends <- list(
    list(
      palette = "palette_avg_mean",
      enabled = TRUE,
      title = "Overall Mean",
      labels = c(". 0", rep("", 5), "0.5", rep("", 4), "1")
    ),
    list(
      palette = "Datasets",
      enabled = TRUE,
      title = "Dataset separation",
      labels = c(".     Min", rep("", ncolors - 2), "Max   .")
    ),
    list(
      palette = "Metrics",
      enabled = TRUE,
      title = "Metric separation",
      labels = c(". 0", rep("", round((ncolors - 2) / 2) + 1), "0.5",
                 rep("", round((ncolors - 2) / 2)), "1")
    )
  )

  ns_fh <- asNamespace("funkyheatmap")
  position_args <- funkyheatmap::position_arguments(
    col_width = 1.5,
    col_space = 0.25,
    col_annot_offset = 5.4,
    expand_ymax = 0.15
  )

  data_v <- ns_fh$verify_data(fh_data)
  col_v <- ns_fh$verify_column_info(column_info, data_v)
  row_v <- ns_fh$verify_row_info(NULL, data_v)
  cg_v <- ns_fh$verify_column_groups(column_groups, col_v)
  rg_v <- ns_fh$verify_row_groups(NULL, row_v)
  pal_v <- ns_fh$verify_palettes(palettes, col_v, data_v)
  leg_v <- ns_fh$verify_legends(legends, pal_v, col_v, data_v)

  geom_positions <- ns_fh$calculate_geom_positions(
    data_v, col_v, row_v, cg_v, rg_v, pal_v, position_args,
    scale_column = FALSE, add_abc = TRUE
  )

  # Style row labels: bold #000000 for corrected methods, plain muted gray #555555 for uncorrected/baseline
  is_cor <- setNames(ranking_df$state == "corrected", ranking_df$row_id)
  if ("text_data" %in% names(geom_positions) && nrow(geom_positions$text_data) > 0) {
    td <- geom_positions$text_data
    if (!"fontface" %in% names(td)) td$fontface <- "plain"
    if (!"colour" %in% names(td)) td$colour <- "black"
    td$fontface[is.na(td$fontface)] <- "plain"
    td$colour[is.na(td$colour)] <- "black"

    idx_method <- which(td$column_id == "method")
    if (length(idx_method)) {
      rids <- td$row_id[idx_method]
      cor_vec <- unname(is_cor[rids])
      td$fontface[idx_method] <- ifelse(cor_vec, "bold", "plain")
      td$colour[idx_method] <- ifelse(cor_vec, "#000000", "#555555")
    }

    idx_cor <- which(td$column_id == "Corrected")
    if (length(idx_cor)) {
      td$fontface[idx_cor] <- ifelse(td$label_value[idx_cor] == "X", "bold", "plain")
      td$colour[idx_cor] <- ifelse(td$label_value[idx_cor] == "X", "#000000", "#555555")
    }

    idx_corr_col <- which(td$column_id == "Correction")
    if (length(idx_corr_col)) {
      rids <- td$row_id[idx_corr_col]
      cor_vec <- unname(is_cor[rids])
      td$fontface[idx_corr_col] <- "plain"
      td$colour[idx_corr_col] <- ifelse(cor_vec, "#000000", "#555555")
    }
    geom_positions$text_data <- td
  }

  main_plot <- ns_fh$compose_ggplot(geom_positions, position_args)

  geom_legend_funs <- list(
    funkyrect = ns_fh$create_funkyrect_legend,
    circle = ns_fh$create_circle_legend,
    rect = ns_fh$create_rect_legend,
    pie = ns_fh$create_pie_legend,
    text = ns_fh$create_text_legend,
    image = ns_fh$create_image_legend,
    bar = ns_fh$create_bar_legend
  )
  legend_plots <- list()
  for (legend in leg_v) {
    if (isTRUE(legend$enabled)) {
      legend_fun <- geom_legend_funs[[legend$geom]]
      legend_args <- legend
      legend_args$geom <- NULL
      legend_args$enabled <- NULL
      legend_args$palette <- NULL
      legend_args$position_args <- position_args
      legend_plots <- c(legend_plots, list(do.call(legend_fun, legend_args)))
    }
  }
  legend_widths <- vapply(legend_plots, function(x) x$width, numeric(1))
  legend_heights <- vapply(legend_plots, function(x) x$height, numeric(1))
  out <- patchwork::wrap_plots(
    main_plot,
    patchwork::plot_spacer(),
    patchwork::wrap_plots(legend_plots, nrow = 1, widths = legend_widths),
    ncol = 1,
    heights = c(main_plot$height, 0.1, max(legend_heights))
  )
  out$width <- max(main_plot$width, sum(legend_widths))
  out$height <- main_plot$height + 0.1 + max(legend_heights)

  ggsave(output_pdf, plot = out, width = plot_width, height = plot_height, device = cairo_pdf)
  invisible(out)
}

p5A <- build_batch_6cohort_funkyheatmap(
  ranking_combined,
  file.path("plots", "Figure_5_A_batch_6cohorts_funkyheatmap.pdf")
)
message("Wrote plots/Figure_5_A_batch_6cohorts_funkyheatmap.pdf")

message("\n=== Step 3: Figure 5B (2x3 ANOSIM Arrows) & Composite Figure 5 ===")
arrow_companion_dir <- file.path("plots", "batch_effect_corrected_6cohorts_20261126")
dir.create(arrow_companion_dir, recursive = TRUE, showWarnings = FALSE)

arrow_method_order <- c(
  "ECODA_authors",
  "ECODA_Leiden: Harmony",
  "ECODA_Leiden: Harmony + CT-drop",
  "MrVI",
  "GloScope",
  "Pseudobulk",
  "QOT",
  "PILOT"
)

palette_arrows <- setNames(
  c("#1B9E77", "#0072B2", "#984EA3", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7", "#4D4D4D"),
  arrow_method_order
)

build_arrow_panel_figure_6cohort <- function(arrow_df, metric_name) {
  sub <- arrow_df[arrow_df$metric == metric_name, , drop = FALSE]
  sub$dataset_label <- factor(
    cohort_display_map[sub$dataset],
    levels = unname(cohort_display_map[retained_cohorts])
  )
  sub$method <- factor(sub$method_label, levels = arrow_method_order)
  stopifnot(!anyNA(sub$dataset_label), !anyNA(sub$method))

  metric_nice <- c(
    anosim_score = "ANOSIM",
    mod_knn3_score = "Modularity (kNN=3)",
    cluster_score = "ARI",
    sil_score = "Silhouette",
    lisi_score = "sLISI"
  )[[metric_name]]

  x_rng <- range(c(sub$biological_uncorrected, sub$biological_corrected))
  y_rng <- range(c(sub$batch_uncorrected, sub$batch_corrected))
  x_lim <- c(min(-0.08, x_rng[1] - 0.04), max(0.68, x_rng[2] + 0.05))
  y_lim <- c(min(-0.08, y_rng[1] - 0.04), max(1.02, y_rng[2] + 0.04))

  facet_panels <- lapply(seq_along(retained_cohorts), function(i) {
    ds <- retained_cohorts[i]
    df_ds <- sub[sub$dataset == ds, , drop = FALSE]
    pts_ds <- rbind(
      transform(df_ds, biological = biological_uncorrected, batch = batch_uncorrected, state = "Uncorrected"),
      transform(df_ds, biological = biological_corrected, batch = batch_corrected, state = "Corrected")
    )
    pts_ds$state <- factor(pts_ds$state, levels = c("Uncorrected", "Corrected"))
    row_idx <- ceiling(i / 3)
    col_idx <- ((i - 1) %% 3) + 1

    p <- ggplot(df_ds, aes(colour = method)) +
      geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_segment(
        aes(x = biological_uncorrected, y = batch_uncorrected,
            xend = biological_corrected, yend = batch_corrected),
        linewidth = 0.6, alpha = 0.85,
        arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "closed"),
        show.legend = FALSE
      ) +
      geom_point(data = pts_ds, aes(x = biological, y = batch, shape = state), size = 2.2, stroke = 0.75) +
      facet_wrap(~dataset_label) +
      scale_colour_manual(values = palette_arrows, drop = FALSE) +
      scale_shape_manual(values = c(Uncorrected = 1, Corrected = 16), drop = FALSE) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_x_continuous(
        breaks = scales::breaks_pretty(n = 4),
        labels = function(z) format(signif(z, 2), trim = TRUE, scientific = FALSE)
      ) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF3", colour = "grey70"),
        strip.text = element_text(size = 10, face = "bold", margin = margin(t = 3.5, b = 3.5)),
        axis.title = element_blank(),
        axis.text = element_text(size = 8.8),
        plot.margin = margin(4, 5, 4, 5)
      ) +
      labs(colour = "Method", shape = "State") +
      guides(
        colour = guide_legend(nrow = 2, byrow = TRUE, order = 1, override.aes = list(size = 2.6)),
        shape = guide_legend(nrow = 1, order = 2, override.aes = list(size = 2.6))
      )
    if (col_idx != 1) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    if (row_idx == 1) {
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    p
  })

  core_grid <- patchwork::wrap_plots(facet_panels, ncol = 3, nrow = 2, guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(size = 9.5, face = "bold"),
      legend.text = element_text(size = 8.8),
      legend.key.height = grid::unit(4.0, "mm"),
      legend.spacing.y = grid::unit(1.5, "mm"),
      legend.margin = margin(t = 2, b = 0)
    )

  # Place outer axis labels cleanly around the 2x3 grid (bottom legend sits below x-axis title or grid)
  # Use patchwork plot_annotation/wrap or cowplot with space for bottom legend:
  leg <- cowplot::get_legend(
    facet_panels[[1]] + theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 9.5, face = "bold"),
      legend.text = element_text(size = 8.8),
      legend.key.height = grid::unit(4.0, "mm"),
      legend.spacing.x = grid::unit(4.0, "mm")
    )
  )
  panels_noleg <- patchwork::wrap_plots(
    lapply(facet_panels, function(p) p + theme(legend.position = "none")),
    ncol = 3, nrow = 2
  )
  plot_with_axes <- cowplot::ggdraw() +
    cowplot::draw_plot(panels_noleg, x = 0.038, y = 0.06, width = 0.955, height = 0.93) +
    cowplot::draw_label(
      paste0("Biological separation (", metric_nice, ")"),
      x = 0.515, y = 0.025, size = 11, fontface = "bold"
    ) +
    cowplot::draw_label(
      paste0("Batch separation (", metric_nice, ")"),
      x = 0.015, y = 0.525, angle = 90, size = 11, fontface = "bold"
    )

  cowplot::plot_grid(plot_with_axes, leg, ncol = 1, rel_heights = c(1, 0.14))
}

p5B <- build_arrow_panel_figure_6cohort(arrow_6, "anosim_score")
ggsave(file.path("plots", "Figure_5_B_batch_correction_arrows_ANOSIM_6cohorts.pdf"), p5B, width = 13.5, height = 7.2, device = cairo_pdf)
for (m in evaluation_metrics()) {
  pm <- build_arrow_panel_figure_6cohort(arrow_6, m)
  ggsave(file.path(arrow_companion_dir, paste0("batch_effect_corrected_6cohorts_", m, "_arrows.pdf")),
         pm, width = 13.5, height = 7.2, device = cairo_pdf)
}

fig5_combined <- cowplot::plot_grid(
  p5A,
  p5B,
  ncol = 1,
  rel_heights = c(1.02, 0.95),
  labels = c("A", "B"),
  label_size = 15,
  label_fontface = "bold"
)
ggsave(file.path("plots", "Figure_5.pdf"), fig5_combined, width = 14.2, height = 14.6, device = cairo_pdf)
ggsave(file.path("plots", "Figure_5.png"), fig5_combined, width = 14.2, height = 14.6, dpi = 220)
message("Wrote plots/Figure_5.{pdf,png} and Figure_5_A/B standalone PDFs.")

message("\n=== Step 4: SF26 (Supp_fig_TMP_ECODA_becm_comp_ANOSIM, 2x3, with limma, droplow50pct hidden) ===")
build_becm_anosim_6cohort <- function(scores_df) {
  # Filter to the 6 retained cohorts and hide droplow50pct only (limma restored)
  data <- scores_df[
    scores_df$dataset %in% retained_cohorts &
      scores_df$method_id != "ECODA_authors_sigdiff_ctdrop_2_droplow50pct",
    ,
    drop = FALSE
  ]
  author_drop <- grepl("^ECODA_authors_sigdiff_ctdrop_", data$method_id) & !grepl("droplow", data$method_id)
  harmony_drop <- grepl("^ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_", data$method_id)

  auth_traj_label <- "ECODA_authors: |\u0394CLR| \u2265 0..4"
  harm_traj_label <- "ECODA_Leiden (Harmony): |\u0394CLR| \u2265 0..4"

  labels <- c(
    ECODA_authors_uncorrected = "ECODA_authors (uncorrected)",
    ECODA_Leiden_CLR_uncorrected = "ECODA_Leiden (uncorrected)",
    ECODA_authors_limma = "ECODA_authors (limma)",
    ECODA_Harmony_Leiden_CLR = "ECODA_Leiden (Harmony)",
    ECODA_Harmony_Leiden_CLR_limma = "ECODA_Leiden (Harmony + limma)"
  )
  series_levels <- c(
    "ECODA_authors (uncorrected)",
    "ECODA_Leiden (uncorrected)",
    "ECODA_authors (limma)",
    "ECODA_Leiden (Harmony)",
    "ECODA_Leiden (Harmony + limma)",
    auth_traj_label,
    harm_traj_label
  )
  data$series <- unname(labels[data$method_id])
  data$series[author_drop] <- auth_traj_label
  data$series[harmony_drop] <- harm_traj_label
  stopifnot(!anyNA(data$series))
  data$series <- factor(data$series, levels = series_levels)
  data$valid <- is.finite(data$biological_anosim) & is.finite(data$batch_anosim)

  palette <- setNames(
    c("#333333", "#888888", "#1B9E77", "#377EB8", "#984EA3", "#D73027", "#E68A00"),
    series_levels
  )
  shapes <- setNames(
    c(1, 2, 17, 15, 18, 16, 16),
    series_levels
  )
  arrow_pairs <- list(
    c("ECODA_authors_uncorrected", "ECODA_authors_limma"),
    c("ECODA_Leiden_CLR_uncorrected", "ECODA_Harmony_Leiden_CLR"),
    c("ECODA_Leiden_CLR_uncorrected", "ECODA_Harmony_Leiden_CLR_limma")
  )

  x_lim <- c(-0.04, 0.72)
  y_lim <- c(-0.09, 1.05)

  panels <- lapply(seq_along(retained_cohorts), function(i) {
    ds <- retained_cohorts[i]
    d_ds <- data[data$dataset == ds, , drop = FALSE]
    n_auth <- d_ds$n_initial[d_ds$method_id == "ECODA_authors_uncorrected"]
    n_harm <- d_ds$n_initial[d_ds$method_id == "ECODA_Harmony_Leiden_CLR"]
    stopifnot(length(n_auth) == 1L, length(n_harm) == 1L)
    strip_txt <- paste0(cohort_display_map[[ds]], "\n(", n_auth, " author / ", n_harm, " Harmony cell types)")
    d_ds$strip <- factor(strip_txt, levels = strip_txt)

    traj <- d_ds[grepl("^ECODA_authors_sigdiff_ctdrop_", d_ds$method_id) |
                   grepl("^ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_", d_ds$method_id), , drop = FALSE]
    traj <- traj[order(traj$series, traj$threshold), ]
    groups <- interaction(traj$series, drop = TRUE)
    traj$segment <- ave(as.integer(!traj$valid), groups, FUN = cumsum)
    traj$line <- interaction(traj$series, traj$segment, drop = TRUE)
    line_ds <- traj[traj$valid, , drop = FALSE]
    pts_ds <- d_ds[d_ds$valid, , drop = FALSE]

    arr_ds <- do.call(rbind, lapply(arrow_pairs, function(pair) {
      from <- pts_ds[pts_ds$method_id == pair[1], c("biological_anosim", "batch_anosim")]
      to <- pts_ds[pts_ds$method_id == pair[2], c("biological_anosim", "batch_anosim", "series")]
      if (!nrow(from) || !nrow(to)) return(NULL)
      cbind(
        stats::setNames(from, c("x_from", "y_from")),
        stats::setNames(to, c("x_to", "y_to", "series"))
      )
    }))

    row_idx <- ceiling(i / 3)
    col_idx <- ((i - 1) %% 3) + 1
    p <- ggplot(pts_ds, aes(x = biological_anosim, y = batch_anosim, colour = series, shape = series)) +
      geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_segment(
        data = arr_ds,
        aes(x = x_from, y = y_from, xend = x_to, yend = y_to, colour = series),
        inherit.aes = FALSE, linewidth = 0.55, alpha = 0.85,
        arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "closed"),
        show.legend = FALSE
      ) +
      geom_path(data = line_ds, aes(group = line), linewidth = 0.65, show.legend = FALSE) +
      geom_point(size = 2.2) +
      ggrepel::geom_text_repel(
        data = line_ds, aes(label = threshold),
        size = 2.6, max.overlaps = Inf, box.padding = 0.36, point.padding = 0.24,
        force = 5.5, force_pull = 0.3, max.iter = 8000,
        min.segment.length = 0.1, segment.size = 0.25, segment.alpha = 0.65,
        seed = 42, show.legend = FALSE
      ) +
      facet_wrap(~strip) +
      scale_colour_manual(values = palette, limits = series_levels, drop = FALSE) +
      scale_shape_manual(values = shapes, limits = series_levels, drop = FALSE) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), labels = c("0", "0.2", "0.4", "0.6")) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF3", colour = "grey70"),
        strip.text = element_text(size = 9.5, face = "bold", margin = margin(t = 3.5, b = 3.5)),
        axis.title = element_blank(),
        axis.text = element_text(size = 8.8),
        plot.margin = margin(4, 5, 4, 5)
      ) +
      labs(colour = NULL, shape = NULL) +
      guides(
        colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.7)),
        shape = guide_legend(nrow = 2, byrow = TRUE)
      )
    if (col_idx != 1) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (row_idx == 1) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p
  })

  leg <- cowplot::get_legend(
    panels[[1]] + theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9.5),
      legend.key.height = grid::unit(4.5, "mm"),
      legend.spacing.x = grid::unit(5, "mm")
    )
  )
  panels_noleg <- patchwork::wrap_plots(
    lapply(panels, function(p) p + theme(legend.position = "none")),
    ncol = 3, nrow = 2
  )
  plot_with_axes <- cowplot::ggdraw() +
    cowplot::draw_plot(panels_noleg, x = 0.038, y = 0.06, width = 0.955, height = 0.935) +
    cowplot::draw_label(
      "Biological ANOSIM (higher = stronger biological separation)",
      x = 0.515, y = 0.025, size = 11, fontface = "bold"
    ) +
    cowplot::draw_label(
      "Batch ANOSIM (higher = stronger batch effect)",
      x = 0.015, y = 0.525, angle = 90, size = 11, fontface = "bold"
    )

  cowplot::plot_grid(plot_with_axes, leg, ncol = 1, rel_heights = c(1, 0.12))
}

p_becm <- build_becm_anosim_6cohort(becm_cache$scores)
ggsave(file.path("plots", "Supp_fig_TMP_ECODA_becm_comp_ANOSIM.pdf"), p_becm, width = 13.5, height = 8.2, device = cairo_pdf)
ggsave(file.path("plots", "Supp_fig_TMP_ECODA_becm_comp_ANOSIM.png"), p_becm, width = 13.5, height = 8.2, dpi = 220)
message("Wrote plots/Supp_fig_TMP_ECODA_becm_comp_ANOSIM.{pdf,png}.")

message("\n=== Step 5: BatchMDS (SF27-SF32, internal with scores + publication clean) ===")
release_path <- normalizePath("data/derived/batch_access_v2/batch_analysis_release.rds", mustWork = TRUE)
release_md5 <- unname(tools::md5sum(release_path))
ctrl_unc <- read_batch_release_control(release_path, release_md5, "batch_effect_uncorrected")
ctrl_cor <- read_batch_release_control(release_path, release_md5, "batch_effect_corrected")
lupus_local_dir <- normalizePath("data/derived/lupus_strict_stage5_20261121", mustWork = TRUE)

unc_datasets <- list()
cor_datasets <- list()
for (ds in setdiff(retained_cohorts, "Lupus_PBMC")) {
  left <- load_batch_release_dataset(ctrl_unc, ds)
  right <- load_batch_release_dataset(ctrl_cor, ds)
  prim <- right$primary_evaluation
  stopifnot(identical(left$primary_evaluation, prim))
  if (length(right$evaluations) > 1L || !identical(prim, "primary")) {
    left$evaluations <- list(primary = left$evaluations[[prim]])
    left$primary_evaluation <- "primary"
    right$evaluations <- list(primary = right$evaluations[[prim]])
    right$primary_evaluation <- "primary"
  }
  unc_datasets[[ds]] <- left
  cor_datasets[[ds]] <- right
}
unc_datasets[["Lupus_PBMC"]] <- load_strict_lupus_view("batch_effect_uncorrected", lupus_local_dir)
cor_datasets[["Lupus_PBMC"]] <- load_strict_lupus_view("batch_effect_corrected", lupus_local_dir)

# Attach CLR >= 2 drop dist_mat from comparison_cache.rds (verified sample order)
for (ds in retained_cohorts) {
  ids <- cor_datasets[[ds]]$sample_ids
  cm <- becm_cache$datasets[[ds]]$metadata
  stopifnot(setequal(cm$Sample, ids), length(cm$Sample) == length(ids))
  for (vk in clr2_variants) {
    dmat <- as.matrix(becm_cache$datasets[[ds]]$variants[[vk]]$dist_mat)
    cor_datasets[[ds]]$methods[[vk]] <- list(
      dist_mat = stats::as.dist(dmat[ids, ids, drop = FALSE])
    )
  }
}

mds_rows <- tibble::tribble(
  ~row_label,                         ~unc_key,            ~cor_key,
  "ECODA_authors",                    "ECODA_authors_HR",  "ECODA_authors_sigdiff_ctdrop_2",
  "ECODA_Leiden: Harmony",            "ECODA_seuratres_2", "ecoda_unsupervised_harmony",
  "ECODA_Leiden: Harmony + CT-drop",  "ECODA_seuratres_2", "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_2",
  "Pseudobulk",                       "Pseudobulk_hvg2000","Pseudobulk_hvg2000",
  "GloScope",                         "GloScope_hvg2000_pcadims30", "GloScope_hvg2000_pcadims30",
  "MrVI",                             "MrVI_hvg2000",      "MrVI_hvg2000",
  "PILOT",                            "PILOT_hvg2000",     "PILOT_hvg2000",
  "QOT",                              "QOT_hvg2000",       "QOT_hvg2000"
)

compute_mds_coords <- function(dist_obj, ids) {
  d <- as.matrix(dist_obj)[ids, ids, drop = FALSE]
  res <- stats::cmdscale(stats::as.dist(d), k = 2, eig = TRUE)
  eig <- res$eig
  pos <- sum(eig[eig > 0])
  pct <- if (pos > 0) signif(100 * eig[1:2] / pos, 2L) else c(0, 0)
  pts <- res$points
  if (ncol(pts) < 2L) pts <- cbind(pts, 0)
  list(
    coords = data.frame(Sample = ids, MDS1 = pts[, 1], MDS2 = pts[, 2], stringsAsFactors = FALSE),
    pct = pct
  )
}

build_cohort_4col_mds <- function(ds, include_scores = TRUE) {
  left <- unc_datasets[[ds]]
  right <- cor_datasets[[ds]]
  ids <- right$sample_ids
  stopifnot(identical(left$sample_ids, ids))
  tech_col <- cohort_tech_col[[ds]]
  bio_col <- cohort_bio_col[[ds]]

  bio_vals <- as.character(right$evaluations$primary[ids])
  tech_raw <- as.character(right$metadata[ids, tech_col])
  tech_short <- clean_short_batch(ds, tech_raw)

  bio_levels <- sort(unique(bio_vals))
  tech_levels <- sort(unique(tech_short))
  bio_fac <- factor(bio_vals, levels = bio_levels)
  tech_fac <- factor(tech_short, levels = tech_levels)

  bio_pal <- setNames(scales::hue_pal()(length(bio_levels)), bio_levels)
  tech_pal <- setNames(scales::brewer_pal(palette = "Dark2")(max(3, length(tech_levels))), tech_levels)
  if (length(tech_levels) > 8L) {
    tech_pal <- setNames(scales::hue_pal(h = c(30, 330), l = 50)(length(tech_levels)), tech_levels)
  }

  fmt_metrics <- function(dist_obj, grp_vec) {
    sc <- calc_sep_score(dist_obj, grp_vec)
    sprintf("ANOSIM %.2f | Mod %.2f | ARI %.2f\nSil %.2f | sLISI %.2f",
            sc$anosim_score, sc$mod_knn3_score, sc$cluster_score, sc$sil_score, sc$lisi_score)
  }

  make_cell_plot <- function(mds_res, grp_fac, pal, sub_txt = NULL) {
    df <- cbind(mds_res$coords, Group = grp_fac)
    xr <- range(df$MDS1); yr <- range(df$MDS2)
    padx <- if (diff(xr) > 0) diff(xr) * 0.08 else 0.5
    pady <- if (diff(yr) > 0) diff(yr) * 0.08 else 0.5
    ggplot(df, aes(x = MDS1, y = MDS2, colour = Group)) +
      geom_point(size = if (include_scores) 1.9 else 2.2, alpha = 0.88) +
      scale_colour_manual(values = pal, drop = FALSE) +
      coord_cartesian(xlim = c(xr[1] - padx, xr[2] + padx), ylim = c(yr[1] - pady, yr[2] + pady)) +
      labs(
        subtitle = sub_txt,
        x = paste0("Dim1 (", mds_res$pct[1], "%)"),
        y = paste0("Dim2 (", mds_res$pct[2], "%)")
      ) +
      theme_classic(base_size = if (include_scores) 10.5 else 11.5) +
      theme(
        legend.position = "none",
        plot.subtitle = if (include_scores) {
          element_text(size = 9.8, lineheight = 1.08, hjust = 0.5, margin = margin(b = 3))
        } else {
          element_blank()
        },
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = if (include_scores) 9.5 else 10.5, colour = "grey25"),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
        plot.margin = margin(3, 4, 3, 4)
      )
  }

  leg_title_size <- if (include_scores) 12 else 14
  leg_text_size <- if (include_scores) 11 else 13
  leg_pt_size <- if (include_scores) 3.4 else 4.2

  legend_tech_plot <- ggplot(data.frame(x = 1, y = 1, g = factor(tech_levels, levels = tech_levels)),
                             aes(x, y, colour = g)) +
    geom_point(size = leg_pt_size) +
    scale_colour_manual(values = tech_pal, name = paste0("Batch (", tech_col, "): "), drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = leg_title_size, face = "bold"),
      legend.text = element_text(size = leg_text_size),
      legend.spacing.x = grid::unit(3.5, "mm")
    ) +
    guides(colour = guide_legend(nrow = if (length(tech_levels) > 4L || sum(nchar(c(tech_col, tech_levels))) > 36L) 2 else 1, byrow = TRUE))

  legend_bio_plot <- ggplot(data.frame(x = 1, y = 1, g = factor(bio_levels, levels = bio_levels)),
                            aes(x, y, colour = g)) +
    geom_point(size = leg_pt_size) +
    scale_colour_manual(values = bio_pal, name = paste0("Biology (", bio_col, "): "), drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = leg_title_size, face = "bold"),
      legend.text = element_text(size = leg_text_size),
      legend.spacing.x = grid::unit(3.5, "mm")
    ) +
    guides(colour = guide_legend(nrow = if (length(bio_levels) > 4L || sum(nchar(c(bio_col, bio_levels))) > 36L) 2 else 1, byrow = TRUE))

  leg_tech <- cowplot::get_legend(legend_tech_plot)
  leg_bio <- cowplot::get_legend(legend_bio_plot)

  banner_box <- function(txt, fill_hex = "#ECEFF3") {
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = fill_hex, colour = "grey70", linewidth = 0.4) +
      annotate("text", x = 0.5, y = 0.5, label = txt, fontface = "bold", size = if (include_scores) 4.8 else 5.2) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(2, 3, 2, 3))
  }

  row_label_box <- function(txt) {
    pretty <- sub(": ", ":\n", txt, fixed = TRUE)
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = pretty, fontface = "bold",
               size = if (include_scores) 4.1 else 4.3, lineheight = 0.95, angle = 90) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(2, 1, 2, 1))
  }

  divider_col <- ggplot() +
    geom_vline(xintercept = 0.5, colour = "grey35", linewidth = 0.7) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  header_row <- patchwork::wrap_plots(
    patchwork::plot_spacer(),
    banner_box("Batch | Uncorrected", "#EAECEF"),
    banner_box("Batch | Corrected", "#EAECEF"),
    divider_col,
    banner_box("Biology | Uncorrected", "#E6EFF6"),
    banner_box("Biology | Corrected", "#E6EFF6"),
    ncol = 6,
    widths = c(0.64, 1, 1, 0.07, 1, 1)
  )

  row_grobs <- vector("list", nrow(mds_rows))
  for (r in seq_len(nrow(mds_rows))) {
    mk_unc <- mds_rows$unc_key[r]
    mk_cor <- mds_rows$cor_key[r]
    mds_unc <- compute_mds_coords(left$methods[[mk_unc]]$dist_mat, ids)
    mds_cor <- compute_mds_coords(right$methods[[mk_cor]]$dist_mat, ids)

    sub_c1 <- if (include_scores) fmt_metrics(left$methods[[mk_unc]]$dist_mat, as.character(tech_fac)) else NULL
    sub_c2 <- if (include_scores) fmt_metrics(right$methods[[mk_cor]]$dist_mat, as.character(tech_fac)) else NULL
    sub_c3 <- if (include_scores) fmt_metrics(left$methods[[mk_unc]]$dist_mat, as.character(bio_fac)) else NULL
    sub_c4 <- if (include_scores) fmt_metrics(right$methods[[mk_cor]]$dist_mat, as.character(bio_fac)) else NULL

    c1 <- make_cell_plot(mds_unc, tech_fac, tech_pal, sub_c1)
    c2 <- make_cell_plot(mds_cor, tech_fac, tech_pal, sub_c2)
    c3 <- make_cell_plot(mds_unc, bio_fac, bio_pal, sub_c3)
    c4 <- make_cell_plot(mds_cor, bio_fac, bio_pal, sub_c4)

    row_grobs[[r]] <- patchwork::wrap_plots(
      row_label_box(mds_rows$row_label[r]),
      c1, c2, divider_col, c3, c4,
      ncol = 6,
      widths = c(0.64, 1, 1, 0.07, 1, 1)
    )
  }

  body_grid <- patchwork::wrap_plots(c(list(header_row), row_grobs), ncol = 1, heights = c(0.30, rep(1, nrow(mds_rows))))
  legends_row <- cowplot::plot_grid(NULL, leg_tech, NULL, leg_bio, ncol = 4, rel_widths = c(0.64, 2.0, 0.07, 2.0))
  title_txt <- if (include_scores) {
    paste0(cohort_display_map[[ds]], " (n = ", length(ids), " samples) | Batch vs Biological MDS")
  } else {
    cohort_display_map[[ds]]
  }
  title_grob <- cowplot::ggdraw() + cowplot::draw_label(
    title_txt, fontface = "bold", size = if (include_scores) 18 else 19, x = 0.5, y = 0.5
  )

  cowplot::plot_grid(
    title_grob,
    body_grid,
    legends_row,
    ncol = 1,
    rel_heights = if (include_scores) c(0.036, 0.896, 0.068) else c(0.034, 0.884, 0.082)
  )
}

internal_mds_dir <- file.path("plots", "batch_mds_internal_with_scores_20261126")
dir.create(internal_mds_dir, recursive = TRUE, showWarnings = FALSE)

# Move Alzheimer BatchMDS out of publication plots/ if present
alz_pub <- file.path("plots", "Supp_fig_TMP_BatchMDS_Alzheimer.pdf")
if (file.exists(alz_pub)) {
  file.rename(alz_pub, file.path(internal_mds_dir, "Supp_fig_TMP_BatchMDS_Alzheimer.pdf"))
  message("Moved plots/Supp_fig_TMP_BatchMDS_Alzheimer.pdf to ", internal_mds_dir)
}

for (ds in retained_cohorts) {
  fig_internal <- build_cohort_4col_mds(ds, include_scores = TRUE)
  out_int <- file.path(internal_mds_dir, paste0("BatchMDS_", ds, "_with_scores.pdf"))
  ggsave(out_int, fig_internal, width = 13.5, height = 21.5, device = cairo_pdf)
  message("Wrote ", out_int)

  fig_pub <- build_cohort_4col_mds(ds, include_scores = FALSE)
  out_pub <- file.path("plots", paste0("Supp_fig_TMP_BatchMDS_", ds, ".pdf"))
  ggsave(out_pub, fig_pub, width = 13.5, height = 19.5, device = cairo_pdf)
  message("Wrote ", out_pub)
}

message("ALL_6COHORT_BATCH_FIGURES_COMPLETE")
