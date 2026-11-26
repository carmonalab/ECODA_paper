#!/usr/bin/env Rscript
# Refresh 7-cohort batch-effect release scores (integrating strict majority_only
# Lupus_PBMC Stage 5 results while retaining the 6 unchanged cohorts from
# batch_access_v2, excluding Lung and Kidney tissue), refresh the
# ECODA_becm_comp cache, and generate all Phase 2B publication-grade figures:
#   - Figure 5 (A: 7-cohort funkyheatmap, B: 7-cohort ANOSIM arrows) + companions
#   - Supp_fig_TMP_ECODA_becm_comp_ANOSIM.{pdf,png}
#   - Supp_fig_TMP_BatchMDS_<Dataset>.pdf (7 files, 4 columns x 8 rows + divider)

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
})

source("src/utils/load_all_functions.R")
source("src/utils/analysis_release.R")
source("src/utils/evaluation_scores.R")
source("src/utils/batch_effect_analysis.R")

set.seed(42)
options(mc.cores = 1L)

retained_cohorts <- c(
  "Alzheimer",
  "Breast_cancer",
  "Covid19_PBMC",
  "Joanito",
  "Kidney_KPMP_full",
  "Lupus_PBMC",
  "Stephenson"
)

cohort_display_map <- c(
  Alzheimer = "Alzheimer",
  Breast_cancer = "Breast cancer",
  Covid19_PBMC = "COVID-19 PBMC",
  Joanito = "Joanito",
  Kidney_KPMP_full = "Kidney (KPMP)",
  Lupus_PBMC = "Lupus PBMC",
  Stephenson = "Stephenson"
)

cohort_tech_col <- c(
  Alzheimer = "assay",
  Breast_cancer = "suspension_dissociation_time",
  Covid19_PBMC = "datasets",
  Joanito = "seqtec",
  Kidney_KPMP_full = "suspension_type",
  Lupus_PBMC = "batch_cov",
  Stephenson = "Site"
)

cohort_bio_col <- c(
  Alzheimer = "Cognitive status",
  Breast_cancer = "disease",
  Covid19_PBMC = "CoVID-19 severity",
  Joanito = "sample.origin",
  Kidney_KPMP_full = "condition.l1",
  Lupus_PBMC = "Status",
  Stephenson = "Status"
)

batch_method_display_map <- c(
  "ECODA_authors_HR" = "ECODA_authors (limma)",
  "ecoda_unsupervised_harmony" = "ECODA_Leiden (Harmony)",
  "ecoda_unsupervised_harmony_limma" = "ECODA_Leiden (Harmony + limma)",
  "MrVI_hvg2000" = "MrVI (MrVI)",
  "GloScope_hvg2000_pcadims30" = "GloScope (Harmony)",
  "Pseudobulk_hvg2000" = "Pseudobulk (limma)",
  "QOT_hvg2000" = "QOT (Harmony)",
  "PILOT_hvg2000" = "PILOT (Harmony)",
  "ECODA_authors_HR_NULL" = "Shuffled_baseline"
)

clean_short_batch <- function(dataset, x) {
  x <- as.character(x)
  if (dataset == "Alzheimer") {
    return(ifelse(x == "10x 3' v3", "3' v3", ifelse(x == "10x Multiome", "Multiome", x)))
  }
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
  true_batch <- setNames(as.character(meta$batch_cov), ids)

  methods_vec <- batch_release_method_keys(view_name)
  method_objs <- list()
  score_list <- list()
  metrics <- evaluation_metrics()
  prefix <- if (view_name == "batch_effect_uncorrected") {
    "Lupus_PBMC_batch_effect_uncorrected_final"
  } else {
    "Lupus_PBMC_batch_effect_corrected_final"
  }

  for (mk in methods_vec) {
    if (mk %in% c("MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000")) {
      short_py <- tolower(sub("_hvg2000$", "", mk))
      fpath <- file.path(base, "embeddings", paste0(prefix, "_hvg2000_highres_", short_py, "_dists.feather"))
      d <- read_dist_feather_strict(fpath, ids)
      bundle <- list(dist_mat = d, labels = true_bio, scores = NULL)
    } else {
      rpath <- file.path(base, "results", paste0(prefix, "_", mk, ".rds"))
      bundle <- read_rds_bundle_strict(rpath, mk, ids)
    }

    eval_labels <- true_bio
    if (identical(mk, "ECODA_authors_HR_NULL")) {
      stopifnot(!is.null(bundle$labels), setequal(names(bundle$labels), ids))
      eval_labels <- setNames(as.character(bundle$labels[ids]), ids)
      stopifnot(identical(sort(unname(eval_labels)), sort(unname(true_bio))))
    }

    sc <- calc_sep_score(bundle$dist_mat, eval_labels)
    vals <- vapply(metrics, function(m) as.numeric(sc[[m]]), numeric(1))
    stopifnot(all(is.finite(vals)))
    if (!is.null(bundle$scores)) {
      for (m in metrics) {
        if (!is.null(bundle$scores[[m]]) && is.finite(bundle$scores[[m]])) {
          stopifnot(abs(as.numeric(bundle$scores[[m]]) - sc[[m]]) < 1e-6)
        }
      }
    }

    method_objs[[mk]] <- list(
      dist_mat = bundle$dist_mat,
      labels = eval_labels,
      scores = sc,
      path = if (mk %in% c("MrVI_hvg2000", "PILOT_hvg2000", "QOT_hvg2000")) fpath else rpath
    )
    score_list[[mk]] <- data.frame(
      dataset = "Lupus_PBMC",
      view = view_name,
      evaluation_id = "primary",
      method_key = mk,
      metric = metrics,
      value = unname(vals),
      score_source = "lupus_strict_20261121",
      stringsAsFactors = FALSE
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
    methods = method_objs,
    evaluated = lapply(methods_vec, function(mk) list(scores = score_list[[mk]])) |> setNames(methods_vec),
    scores = do.call(rbind, score_list)
  )
}

message("=== Step 1: Load 6 unchanged cohorts from batch_access_v2 & strict Lupus_PBMC ===")
release_path <- normalizePath("data/derived/batch_access_v2/batch_analysis_release.rds", mustWork = TRUE)
release_md5 <- unname(tools::md5sum(release_path))
ctrl_unc <- read_batch_release_control(release_path, release_md5, "batch_effect_uncorrected")
ctrl_cor <- read_batch_release_control(release_path, release_md5, "batch_effect_corrected")

lupus_local_dir <- normalizePath("data/derived/lupus_strict_stage5_20261121", mustWork = TRUE)
unchanged_6 <- setdiff(retained_cohorts, "Lupus_PBMC")

unc_datasets <- list()
cor_datasets <- list()
old_lupus_unc <- load_batch_release_dataset(ctrl_unc, "Lupus_PBMC")
old_lupus_cor <- load_batch_release_dataset(ctrl_cor, "Lupus_PBMC")

for (ds in unchanged_6) {
  left <- load_batch_release_dataset(ctrl_unc, ds)
  right <- load_batch_release_dataset(ctrl_cor, ds)
  prim <- right$primary_evaluation
  stopifnot(identical(left$primary_evaluation, prim))
  if (length(right$evaluations) > 1L || !identical(prim, "primary")) {
    left$evaluations <- list(primary = left$evaluations[[prim]])
    left$primary_evaluation <- "primary"
    left$scores <- left$scores[left$scores$evaluation_id == prim, , drop = FALSE]
    left$scores$evaluation_id <- "primary"
    for (mk in names(left$evaluated)) {
      df_m <- left$evaluated[[mk]]$scores
      df_m <- df_m[df_m$evaluation_id == prim, , drop = FALSE]
      df_m$evaluation_id <- "primary"
      left$evaluated[[mk]]$scores <- df_m
    }

    right$evaluations <- list(primary = right$evaluations[[prim]])
    right$primary_evaluation <- "primary"
    right$scores <- right$scores[right$scores$evaluation_id == prim, , drop = FALSE]
    right$scores$evaluation_id <- "primary"
    for (mk in names(right$evaluated)) {
      df_m <- right$evaluated[[mk]]$scores
      df_m <- df_m[df_m$evaluation_id == prim, , drop = FALSE]
      df_m$evaluation_id <- "primary"
      right$evaluated[[mk]]$scores <- df_m
    }
  }
  unc_datasets[[ds]] <- left
  cor_datasets[[ds]] <- right
}
unc_datasets[["Lupus_PBMC"]] <- load_strict_lupus_view("batch_effect_uncorrected", lupus_local_dir)
cor_datasets[["Lupus_PBMC"]] <- load_strict_lupus_view("batch_effect_corrected", lupus_local_dir)

scores_unc_raw <- do.call(rbind, lapply(retained_cohorts, function(ds) unc_datasets[[ds]]$scores))
scores_cor_raw <- do.call(rbind, lapply(retained_cohorts, function(ds) cor_datasets[[ds]]$scores))
rownames(scores_unc_raw) <- rownames(scores_cor_raw) <- NULL

stopifnot(
  nrow(scores_unc_raw) == 7L * 8L * 5L,
  nrow(scores_cor_raw) == 7L * 9L * 5L,
  setequal(unique(scores_unc_raw$dataset), retained_cohorts),
  setequal(unique(scores_cor_raw$dataset), retained_cohorts),
  all(scores_unc_raw$evaluation_id == "primary"),
  all(scores_cor_raw$evaluation_id == "primary")
)

expected_grid_cor <- expand.grid(
  dataset = retained_cohorts,
  view = "batch_effect_corrected",
  evaluation_id = "primary",
  method_key = batch_release_method_keys("batch_effect_corrected"),
  metric = evaluation_metrics(),
  stringsAsFactors = FALSE
)
expected_grid_unc <- expand.grid(
  dataset = retained_cohorts,
  view = "batch_effect_uncorrected",
  evaluation_id = "primary",
  method_key = batch_release_method_keys("batch_effect_uncorrected"),
  metric = evaluation_metrics(),
  stringsAsFactors = FALSE
)

agg_cor <- aggregate_evaluation_scores(scores_cor_raw, expected_grid_cor)
agg_unc <- aggregate_evaluation_scores(scores_unc_raw, expected_grid_unc)

cat("\n=== 7-Cohort Corrected Overall Ranking (mean scaled 5 metrics) ===\n")
print(agg_cor$ranking, row.names = FALSE)

# Build Lupus old (majority_v1) vs new strict (majority_only) comparison table
old_all <- rbind(
  transform(old_lupus_unc$scores, state = "old_majority_v1"),
  transform(old_lupus_cor$scores, state = "old_majority_v1")
)
new_all <- rbind(
  transform(unc_datasets$Lupus_PBMC$scores, state = "strict_majority_only"),
  transform(cor_datasets$Lupus_PBMC$scores, state = "strict_majority_only")
)
cmp_lupus <- merge(
  old_all[, c("view", "method_key", "metric", "value")],
  new_all[, c("view", "method_key", "metric", "value")],
  by = c("view", "method_key", "metric"),
  suffixes = c("_old_v1", "_strict")
)
cmp_lupus$delta <- cmp_lupus$value_strict - cmp_lupus$value_old_v1
cat("\n=== Lupus_PBMC Old (majority_v1) vs Strict (majority_only) ANOSIM comparison ===\n")
print(cmp_lupus[cmp_lupus$metric == "anosim_score", ], row.names = FALSE)

# Compute paired biological vs batch arrow scores for all 7 retained cohorts
cor_labels_non_null <- batch_method_display_map[setdiff(batch_release_method_keys("batch_effect_corrected"), "ECODA_authors_HR_NULL")]
arrow_tables <- lapply(retained_cohorts, function(ds) {
  batch_release_arrow_scores(unc_datasets[[ds]], cor_datasets[[ds]], ctrl_cor$datasets[[ds]]$batch_col, cor_labels_non_null)
})
arrow_scores <- do.call(rbind, arrow_tables)
stopifnot(
  nrow(arrow_scores) == 7L * 8L * 5L,
  setequal(unique(arrow_scores$dataset), retained_cohorts),
  all(arrow_scores$evaluation_id == "primary")
)

out_release_dir <- file.path("data", "derived", "batch_7cohort_release_20261121")
dir.create(out_release_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(agg_unc$scores, file.path(out_release_dir, "batch_7cohort_scores_uncorrected.csv"), row.names = FALSE)
write.csv(agg_cor$scores, file.path(out_release_dir, "batch_7cohort_scores_corrected.csv"), row.names = FALSE)
write.csv(agg_cor$ranking, file.path(out_release_dir, "batch_7cohort_ranking_corrected.csv"), row.names = FALSE)
write.csv(arrow_scores, file.path(out_release_dir, "batch_7cohort_arrow_scores.csv"), row.names = FALSE)
write.csv(cmp_lupus, file.path(out_release_dir, "lupus_old_vs_strict_comparison.csv"), row.names = FALSE)
saveRDS(
  list(
    schema = "ecoda_batch_7cohort_release_20261121_v1",
    retained_cohorts = retained_cohorts,
    excluded = c("Lung", "Kidney_KPMP_full:tissue"),
    ranking_corrected = agg_cor$ranking,
    scores_uncorrected = agg_unc$scores,
    scores_corrected = agg_cor$scores,
    arrow_scores = arrow_scores,
    lupus_comparison = cmp_lupus
  ),
  file.path(out_release_dir, "batch_7cohort_summary.rds")
)

message("\n=== Step 2: Refresh ECODA_becm_comp cache with strict Lupus_PBMC ===")
worker_lines <- readLines("data/derived/ECODA_becm_comp_20261116_v1/worker.R")
becm_env <- new.env(parent = globalenv())
eval(parse(text = worker_lines[1:223]), envir = becm_env)

old_becm_path <- normalizePath("data/derived/ECODA_becm_comp_20261116_v1/comparison_cache.rds", mustWork = TRUE)
stopifnot(identical(unname(tools::md5sum(old_becm_path)), "d22ed3205de4b439f4063f1b88f37ebf"))
becm_cache <- readRDS(old_becm_path)

comp_u <- readRDS(check_local_sidecar(file.path(
  lupus_local_dir, "uncorrected_final/results/Lupus_PBMC_batch_effect_uncorrected_final_composition.rds"
)))
comp_c <- readRDS(check_local_sidecar(file.path(
  lupus_local_dir, "corrected_final/results/Lupus_PBMC_batch_effect_corrected_final_composition.rds"
)))

extract_becm_pair <- function(bundle, ids, verify_clr = FALSE) {
  clr <- as.matrix(bundle$feat_mat)
  padded <- as.matrix(bundle$counts)
  stopifnot(
    identical(rownames(clr), ids),
    identical(rownames(padded), ids),
    identical(dim(clr), dim(padded)),
    all(is.finite(padded)),
    all(padded >= 0.5),
    max(abs(padded - 0.5 - round(padded - 0.5))) < 1e-8,
    identical(colnames(clr), colnames(padded)) || identical(colnames(clr), make.names(colnames(padded), unique = TRUE))
  )
  colnames(clr) <- colnames(padded)
  counts <- padded - 0.5
  if (verify_clr) {
    stopifnot(max(abs(unname(clr) - unname(log(padded) - rowMeans(log(padded))))) < 1e-8)
  }
  d <- as.matrix(bundle$dist_mat)
  stopifnot(max(abs(as.matrix(stats::dist(clr)) - d[ids, ids, drop = FALSE])) < 1e-8)
  list(clr = clr, counts = counts)
}

lupus_meta_df <- unc_datasets$Lupus_PBMC$metadata
lupus_ids <- unc_datasets$Lupus_PBMC$sample_ids
becm_lupus_inputs <- list(
  authors = extract_becm_pair(comp_u$ECODA_authors_HR, lupus_ids, verify_clr = TRUE),
  leiden = extract_becm_pair(comp_u$ECODA_seuratres_2, lupus_ids, verify_clr = TRUE),
  authors_limma = extract_becm_pair(comp_c$ECODA_authors_HR, lupus_ids, verify_clr = FALSE),
  harmony = extract_becm_pair(comp_c$ecoda_unsupervised_harmony, lupus_ids, verify_clr = TRUE),
  harmony_limma = extract_becm_pair(comp_c$ecoda_unsupervised_harmony_limma, lupus_ids, verify_clr = FALSE)
)
stopifnot(
  identical(becm_lupus_inputs$authors$counts, becm_lupus_inputs$authors_limma$counts),
  identical(becm_lupus_inputs$harmony$counts, becm_lupus_inputs$harmony_limma$counts)
)

becm_meta <- data.frame(
  Sample = lupus_ids,
  bio = as.character(lupus_meta_df$Status),
  batch = as.character(lupus_meta_df$batch_cov),
  stringsAsFactors = FALSE
)
lupus_built <- becm_env$becm_build_variants(becm_lupus_inputs, becm_meta$batch)
lupus_becm_scores <- becm_env$becm_score_variants(lupus_built$variants, becm_meta, "Lupus_PBMC")
lupus_becm_boxes <- becm_env$becm_box_summary(
  becm_lupus_inputs$authors$clr,
  becm_lupus_inputs$authors$counts,
  becm_meta$batch,
  lupus_built$screens$authors,
  "Lupus_PBMC"
)

becm_cache$datasets$Lupus_PBMC$metadata <- becm_meta
becm_cache$datasets$Lupus_PBMC$variants <- lupus_built$variants
becm_cache$datasets$Lupus_PBMC$screens <- lupus_built$screens
becm_cache$datasets$Lupus_PBMC$policy$lupus <- "strict cell-level majority_only, 374,642 cells / 94 samples"
becm_cache$boxplots$Lupus_PBMC <- lupus_becm_boxes

non_lupus_idx <- becm_cache$scores$dataset != "Lupus_PBMC"
becm_cache$scores <- rbind(
  becm_cache$scores[non_lupus_idx, , drop = FALSE],
  lupus_becm_scores
)
becm_cache$scores <- becm_cache$scores[order(match(becm_cache$scores$dataset, retained_cohorts)), , drop = FALSE]
rownames(becm_cache$scores) <- NULL

becm_out_dir <- file.path("data", "derived", "ECODA_becm_comp_20261121_v2")
dir.create(becm_out_dir, recursive = TRUE, showWarnings = FALSE)
becm_cache_out <- file.path(becm_out_dir, "comparison_cache.rds")
saveRDS(becm_cache, becm_cache_out, version = 3)
writeLines(
  paste0("MD5=", unname(tools::md5sum(becm_cache_out)), "\nSIZE=", file.info(becm_cache_out)$size, "\nPATH=", normalizePath(becm_cache_out)),
  paste0(becm_cache_out, ".md5")
)
write.csv(becm_cache$scores, file.path(becm_out_dir, "ANOSIM_scores_and_cell_type_counts.csv"), row.names = FALSE)

cat("\n=== Lupus_PBMC refreshed ECODA_becm_comp key variants ===\n")
print(lupus_becm_scores[lupus_becm_scores$method_id %in% c(
  "ECODA_authors_uncorrected", "ECODA_Leiden_CLR_uncorrected", "ECODA_authors_limma",
  "ECODA_Harmony_Leiden_CLR", "ECODA_Harmony_Leiden_CLR_limma",
  "ECODA_authors_sigdiff_ctdrop_1", "ECODA_authors_sigdiff_ctdrop_2",
  "ECODA_authors_sigdiff_ctdrop_2_droplow50pct", "ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_2"
), c("method_id", "biological_anosim", "batch_anosim", "n_initial", "n_significant", "n_dropped", "n_remaining")], row.names = FALSE)

message("\n=== Step 3A: Main Figure 5 (Funkyheatmap + 7-cohort ANOSIM Arrows) ===")
dir.create("plots", recursive = TRUE, showWarnings = FALSE)
arrow_companion_dir <- file.path("plots", "batch_effect_corrected_7cohorts_20261121")
dir.create(arrow_companion_dir, recursive = TRUE, showWarnings = FALSE)

selected_cor_methods <- batch_release_method_keys("batch_effect_corrected")
plot_data_7 <- scores_cor_raw
plot_data_7$method <- plot_data_7$method_key
plot_data_7$score <- plot_data_7$metric

# Stand-alone Panel A + ggplot object for composite Figure 5
p5A <- build_funky_heatmap(
  plot_data_7,
  exec_times = NULL,
  filter_methods = selected_cor_methods,
  method_recode = batch_method_display_map,
  dataset_recode = cohort_display_map,
  score_set = evaluation_metrics(),
  output_file = "Figure_5_A_batch_corrected_funkyheatmap_7cohorts.pdf",
  output_dir = "plots",
  plot_width = 13,
  plot_height = 6.2,
  method_width = 9.8,
  top_padding = 0.6,
  expected_grid = expected_grid_cor
)

# Panel B (and companion metrics): 7 facets in 2x4 grid with legend in empty slot 8
build_arrow_panel_figure <- function(arrow_df, metric_name) {
  sub <- arrow_df[arrow_df$metric == metric_name, , drop = FALSE]
  sub$dataset_label <- factor(cohort_display_map[sub$dataset], levels = unname(cohort_display_map[retained_cohorts]))
  methods_order <- names(cor_labels_non_null)
  method_display <- unname(cor_labels_non_null)
  sub$method <- factor(sub$method_key, levels = methods_order, labels = method_display)

  palette_methods <- setNames(
    c("#1B9E77", "#0072B2", "#984EA3", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7", "#4D4D4D"),
    method_display
  )
  metric_nice <- c(
    anosim_score = "ANOSIM",
    mod_knn3_score = "Modularity (kNN=3)",
    cluster_score = "ARI",
    sil_score = "Silhouette",
    lisi_score = "LISI"
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
    row_idx <- ceiling(i / 4)
    col_idx <- ((i - 1) %% 4) + 1
    p <- ggplot(df_ds, aes(colour = method)) +
      geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_segment(
        aes(x = biological_uncorrected, y = batch_uncorrected,
            xend = biological_corrected, yend = batch_corrected),
        linewidth = 0.55, alpha = 0.85,
        arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed"),
        show.legend = FALSE
      ) +
      geom_point(data = pts_ds, aes(x = biological, y = batch, shape = state), size = 2.0, stroke = 0.7) +
      facet_wrap(~dataset_label) +
      scale_colour_manual(values = palette_methods, drop = FALSE) +
      scale_shape_manual(values = c(Uncorrected = 1, Corrected = 16), drop = FALSE) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 4), labels = clean_becm_axis <- function(z) format(signif(z, 2), trim = TRUE, scientific = FALSE)) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      theme_bw(base_size = 9.5) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF3", colour = "grey70"),
        strip.text = element_text(size = 9.2, face = "bold", margin = margin(t = 3, b = 3)),
        axis.title = element_blank(),
        plot.margin = margin(3, 4, 3, 4)
      ) +
      labs(colour = "Method", shape = "State") +
      guides(
        colour = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2.4)),
        shape = guide_legend(ncol = 2, order = 2, override.aes = list(size = 2.4))
      )
    if (col_idx != 1) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    if (row_idx == 1 && i <= 3) {
      # Panels 1..3 have panel below in row 2; Panel 4 has empty guide slot below so keep x-ticks on Panel 4!
      p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    p
  })

  core_grid <- patchwork::wrap_plots(
    c(facet_panels, list(patchwork::guide_area())),
    ncol = 4, nrow = 2, guides = "collect"
  ) & theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.4),
    legend.key.height = grid::unit(3.6, "mm"),
    legend.spacing.y = grid::unit(2, "mm")
  )

  # Wrap grid inside outer axis labels
  cowplot::ggdraw() +
    cowplot::draw_plot(core_grid, x = 0.032, y = 0.045, width = 0.965, height = 0.95) +
    cowplot::draw_label(paste0("Biological separation (", metric_nice, ")"),
      x = 0.515, y = 0.018, size = 10.5, fontface = "bold") +
    cowplot::draw_label(paste0("Batch separation (", metric_nice, ")"),
      x = 0.015, y = 0.52, angle = 90, size = 10.5, fontface = "bold")
}

p5B <- build_arrow_panel_figure(arrow_scores, "anosim_score")
ggsave(file.path("plots", "Figure_5_B_batch_correction_arrows_ANOSIM_7cohorts.pdf"), p5B, width = 13, height = 6.2, device = "pdf")
for (m in evaluation_metrics()) {
  pm <- build_arrow_panel_figure(arrow_scores, m)
  ggsave(file.path(arrow_companion_dir, paste0("batch_effect_corrected_7cohorts_", m, "_arrows.pdf")),
         pm, width = 13, height = 6.2, device = "pdf")
}

fig5_combined <- cowplot::plot_grid(
  p5A,
  p5B,
  ncol = 1,
  rel_heights = c(0.88, 1.05),
  labels = c("A", "B"),
  label_size = 15,
  label_fontface = "bold"
)
ggsave(file.path("plots", "Figure_5.pdf"), fig5_combined, width = 13.5, height = 12.2, device = "pdf")
message("Wrote plots/Figure_5.pdf and Figure_5_A/B standalone PDFs.")

message("\n=== Step 3B: Supp_fig_TMP_ECODA_becm_comp_ANOSIM ===")
build_becm_anosim_7cohort <- function(scores_df) {
  data <- scores_df[scores_df$dataset %in% retained_cohorts, , drop = FALSE]
  author_drop <- grepl("^ECODA_authors_sigdiff_ctdrop_", data$method_id) & !grepl("droplow", data$method_id)
  harmony_drop <- grepl("^ECODA_Harmony_Leiden_CLR_sigdiff_ctdrop_", data$method_id)
  labels <- c(
    ECODA_authors_uncorrected = "ECODA_authors (uncorrected)",
    ECODA_Leiden_CLR_uncorrected = "ECODA_Leiden (uncorrected)",
    ECODA_authors_limma = "ECODA_authors (limma)",
    ECODA_Harmony_Leiden_CLR = "ECODA_Leiden (Harmony)",
    ECODA_Harmony_Leiden_CLR_limma = "ECODA_Leiden (Harmony + limma)",
    ECODA_authors_sigdiff_ctdrop_2_droplow50pct = "ECODA_authors: cutoff 2 + zero >=50%"
  )
  series_levels <- c(
    "ECODA_authors (uncorrected)", "ECODA_Leiden (uncorrected)",
    "ECODA_authors (limma)", "ECODA_Leiden (Harmony)", "ECODA_Leiden (Harmony + limma)",
    "ECODA_authors: drop trajectory", "ECODA_Leiden (Harmony): drop trajectory",
    "ECODA_authors: cutoff 2 + zero >=50%"
  )
  data$series <- unname(labels[data$method_id])
  data$series[author_drop] <- "ECODA_authors: drop trajectory"
  data$series[harmony_drop] <- "ECODA_Leiden (Harmony): drop trajectory"
  stopifnot(!anyNA(data$series))
  data$series <- factor(data$series, levels = series_levels)
  data$valid <- is.finite(data$biological_anosim) & is.finite(data$batch_anosim)

  palette <- c(
    "ECODA_authors (uncorrected)" = "#333333",
    "ECODA_Leiden (uncorrected)" = "#888888",
    "ECODA_authors (limma)" = "#1B9E77",
    "ECODA_Leiden (Harmony)" = "#377EB8",
    "ECODA_Leiden (Harmony + limma)" = "#984EA3",
    "ECODA_authors: drop trajectory" = "#D73027",
    "ECODA_Leiden (Harmony): drop trajectory" = "#E68A00",
    "ECODA_authors: cutoff 2 + zero >=50%" = "#A65628"
  )
  shapes <- c(
    "ECODA_authors (uncorrected)" = 1,
    "ECODA_Leiden (uncorrected)" = 2,
    "ECODA_authors (limma)" = 17,
    "ECODA_Leiden (Harmony)" = 15,
    "ECODA_Leiden (Harmony + limma)" = 18,
    "ECODA_authors: drop trajectory" = 16,
    "ECODA_Leiden (Harmony): drop trajectory" = 16,
    "ECODA_authors: cutoff 2 + zero >=50%" = 8
  )
  arrow_pairs <- list(
    c("ECODA_Leiden_CLR_uncorrected", "ECODA_Harmony_Leiden_CLR"),
    c("ECODA_Leiden_CLR_uncorrected", "ECODA_Harmony_Leiden_CLR_limma"),
    c("ECODA_authors_uncorrected", "ECODA_authors_limma"),
    c("ECODA_authors_uncorrected", "ECODA_authors_sigdiff_ctdrop_2_droplow50pct")
  )

  x_lim <- c(-0.04, 0.72)
  y_lim <- c(-0.07, 1.05)

  panels <- lapply(seq_along(retained_cohorts), function(i) {
    ds <- retained_cohorts[i]
    d_ds <- data[data$dataset == ds, , drop = FALSE]
    n_auth <- d_ds$n_initial[d_ds$method_id == "ECODA_authors_uncorrected"]
    n_harm <- d_ds$n_initial[d_ds$method_id == "ECODA_Harmony_Leiden_CLR"]
    stopifnot(length(n_auth) == 1L, length(n_harm) == 1L)
    strip_txt <- paste0(cohort_display_map[[ds]], "\n(", n_auth, " author / ", n_harm, " Harmony)")
    d_ds$strip <- factor(strip_txt, levels = strip_txt)

    traj <- d_ds[(grepl("^ECODA_authors_sigdiff_ctdrop_", d_ds$method_id) & !grepl("droplow", d_ds$method_id)) |
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

    row_idx <- ceiling(i / 4)
    col_idx <- ((i - 1) %% 4) + 1
    p <- ggplot(pts_ds, aes(x = biological_anosim, y = batch_anosim, colour = series, shape = series)) +
      geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
      geom_segment(
        data = arr_ds,
        aes(x = x_from, y = y_from, xend = x_to, yend = y_to, colour = series),
        inherit.aes = FALSE, linewidth = 0.48, alpha = 0.85,
        arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed"),
        show.legend = FALSE
      ) +
      geom_path(data = line_ds, aes(group = line), linewidth = 0.6, show.legend = FALSE) +
      geom_point(size = 2.0) +
      ggrepel::geom_text_repel(
        data = line_ds, aes(label = threshold),
        size = 2.35, max.overlaps = Inf, box.padding = 0.22, point.padding = 0.14,
        force = 3.2, force_pull = 0.4, max.iter = 5000,
        min.segment.length = 0.12, segment.size = 0.25, segment.alpha = 0.65,
        seed = 42, show.legend = FALSE
      ) +
      facet_wrap(~strip) +
      scale_colour_manual(values = palette, limits = series_levels, drop = FALSE) +
      scale_shape_manual(values = shapes, limits = series_levels, drop = FALSE) +
      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), labels = c("0", "0.2", "0.4", "0.6")) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      theme_bw(base_size = 9.5) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF3", colour = "grey70"),
        strip.text = element_text(size = 8.8, face = "bold", margin = margin(t = 3, b = 3)),
        axis.title = element_blank(),
        plot.margin = margin(3, 4, 3, 4)
      ) +
      labs(colour = "Method / trajectory", shape = "Method / trajectory") +
      guides(
        colour = guide_legend(ncol = 1, override.aes = list(size = 2.5)),
        shape = guide_legend(ncol = 1)
      )
    if (col_idx != 1) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (row_idx == 1 && i <= 3) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p
  })

  core_grid <- patchwork::wrap_plots(
    c(panels, list(patchwork::guide_area())),
    ncol = 4, nrow = 2, guides = "collect"
  ) & theme(
    legend.position = "right",
    legend.title = element_text(size = 9.2, face = "bold"),
    legend.text = element_text(size = 8.1),
    legend.key.height = grid::unit(4.4, "mm")
  )

  cowplot::ggdraw() +
    cowplot::draw_plot(core_grid, x = 0.035, y = 0.05, width = 0.96, height = 0.945) +
    cowplot::draw_label("Biological ANOSIM (higher = stronger biological separation)",
      x = 0.515, y = 0.02, size = 10.5, fontface = "bold") +
    cowplot::draw_label("Batch ANOSIM (higher = stronger batch effect)",
      x = 0.015, y = 0.52, angle = 90, size = 10.5, fontface = "bold")
}

p_becm <- build_becm_anosim_7cohort(becm_cache$scores)
ggsave(file.path("plots", "Supp_fig_TMP_ECODA_becm_comp_ANOSIM.pdf"), p_becm, width = 13.5, height = 7.5, device = "pdf")
ggsave(file.path("plots", "Supp_fig_TMP_ECODA_becm_comp_ANOSIM.png"), p_becm, width = 13.5, height = 7.5, dpi = 220)
message("Wrote plots/Supp_fig_TMP_ECODA_becm_comp_ANOSIM.{pdf,png}.")

message("\n=== Step 3C: Generate 7 clean 4-column per-cohort Batch MDS figures ===")
mds_row_order <- c(
  "ECODA_authors_HR",
  "ecoda_unsupervised_harmony",
  "ecoda_unsupervised_harmony_limma",
  "Pseudobulk_hvg2000",
  "GloScope_hvg2000_pcadims30",
  "MrVI_hvg2000",
  "PILOT_hvg2000",
  "QOT_hvg2000"
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

build_cohort_4col_mds <- function(ds) {
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

  # Distinct color palettes so batch colors and bio colors cannot visually blur
  bio_pal <- setNames(scales::hue_pal()(length(bio_levels)), bio_levels)
  tech_pal <- setNames(scales::brewer_pal(palette = "Dark2")(max(3, length(tech_levels))), tech_levels)
  if (length(tech_levels) > 8L) {
    tech_pal <- setNames(scales::hue_pal(h = c(30, 330), l = 50)(length(tech_levels)), tech_levels)
  }

  fmt_metrics <- function(dist_obj, grp_vec) {
    sc <- calc_sep_score(dist_obj, grp_vec)
    sprintf("ANOSIM %.2f | Mod %.2f | ARI %.2f\nSil %.2f | LISI %.2f",
            sc$anosim_score, sc$mod_knn3_score, sc$cluster_score, sc$sil_score, sc$lisi_score)
  }

  make_cell_plot <- function(mds_res, grp_fac, pal, sub_txt) {
    df <- cbind(mds_res$coords, Group = grp_fac)
    xr <- range(df$MDS1); yr <- range(df$MDS2)
    padx <- if (diff(xr) > 0) diff(xr) * 0.08 else 0.5
    pady <- if (diff(yr) > 0) diff(yr) * 0.08 else 0.5
    ggplot(df, aes(x = MDS1, y = MDS2, colour = Group)) +
      geom_point(size = 1.9, alpha = 0.88) +
      scale_colour_manual(values = pal, drop = FALSE) +
      coord_cartesian(xlim = c(xr[1] - padx, xr[2] + padx), ylim = c(yr[1] - pady, yr[2] + pady)) +
      labs(
        subtitle = sub_txt,
        x = paste0("Dim1 (", mds_res$pct[1], "%)"),
        y = paste0("Dim2 (", mds_res$pct[2], "%)")
      ) +
      theme_classic(base_size = 10.5) +
      theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9.8, lineheight = 1.08, hjust = 0.5, margin = margin(b = 3)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 9.5, colour = "grey25"),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
        plot.margin = margin(3, 4, 3, 4)
      )
  }

  # Extract shared legends from dummy plots
  legend_tech_plot <- ggplot(data.frame(x = 1, y = 1, g = factor(tech_levels, levels = tech_levels)),
                             aes(x, y, colour = g)) +
    geom_point(size = 3.4) +
    scale_colour_manual(values = tech_pal, name = paste0("Batch (", tech_col, "): "), drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.spacing.x = grid::unit(3, "mm")
    ) +
    guides(colour = guide_legend(nrow = if (length(tech_levels) > 4L || sum(nchar(c(tech_col, tech_levels))) > 42L) 2 else 1, byrow = TRUE))

  legend_bio_plot <- ggplot(data.frame(x = 1, y = 1, g = factor(bio_levels, levels = bio_levels)),
                            aes(x, y, colour = g)) +
    geom_point(size = 3.4) +
    scale_colour_manual(values = bio_pal, name = paste0("Biology (", bio_col, "): "), drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.spacing.x = grid::unit(3, "mm")
    ) +
    guides(colour = guide_legend(nrow = if (length(bio_levels) > 4L || sum(nchar(c(bio_col, bio_levels))) > 42L) 2 else 1, byrow = TRUE))

  leg_tech <- cowplot::get_legend(legend_tech_plot)
  leg_bio <- cowplot::get_legend(legend_bio_plot)

  banner_box <- function(txt, fill_hex = "#ECEFF3") {
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = fill_hex, colour = "grey70", linewidth = 0.4) +
      annotate("text", x = 0.5, y = 0.5, label = txt, fontface = "bold", size = 4.8) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(2, 3, 2, 3))
  }

  row_label_box <- function(txt) {
    pretty <- sub(" \\(", "\n(", txt)
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = pretty, fontface = "bold", size = 4.3, lineheight = 1.0, angle = 90) +
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
    banner_box("Batch — Uncorrected", "#EAECEF"),
    banner_box("Batch — Corrected", "#EAECEF"),
    divider_col,
    banner_box("Biology — Uncorrected", "#E6EFF6"),
    banner_box("Biology — Corrected", "#E6EFF6"),
    ncol = 6,
    widths = c(0.56, 1, 1, 0.07, 1, 1)
  )

  row_grobs <- vector("list", length(mds_row_order))
  for (r in seq_along(mds_row_order)) {
    mk_cor <- mds_row_order[r]
    mk_unc <- if (mk_cor %in% c("ecoda_unsupervised_harmony", "ecoda_unsupervised_harmony_limma")) {
      "ECODA_seuratres_2"
    } else {
      mk_cor
    }
    mds_unc <- compute_mds_coords(left$methods[[mk_unc]]$dist_mat, ids)
    mds_cor <- compute_mds_coords(right$methods[[mk_cor]]$dist_mat, ids)

    c1 <- make_cell_plot(mds_unc, tech_fac, tech_pal, fmt_metrics(left$methods[[mk_unc]]$dist_mat, as.character(tech_fac)))
    c2 <- make_cell_plot(mds_cor, tech_fac, tech_pal, fmt_metrics(right$methods[[mk_cor]]$dist_mat, as.character(tech_fac)))
    c3 <- make_cell_plot(mds_unc, bio_fac, bio_pal, fmt_metrics(left$methods[[mk_unc]]$dist_mat, as.character(bio_fac)))
    c4 <- make_cell_plot(mds_cor, bio_fac, bio_pal, fmt_metrics(right$methods[[mk_cor]]$dist_mat, as.character(bio_fac)))

    row_grobs[[r]] <- patchwork::wrap_plots(
      row_label_box(unname(batch_method_display_map[[mk_cor]])),
      c1, c2, divider_col, c3, c4,
      ncol = 6,
      widths = c(0.56, 1, 1, 0.07, 1, 1)
    )
  }

  body_grid <- patchwork::wrap_plots(c(list(header_row), row_grobs), ncol = 1, heights = c(0.30, rep(1, length(mds_row_order))))
  legends_row <- cowplot::plot_grid(NULL, leg_tech, NULL, leg_bio, ncol = 4, rel_widths = c(0.56, 2.0, 0.07, 2.0))
  title_grob <- cowplot::ggdraw() + cowplot::draw_label(
    paste0(cohort_display_map[[ds]], " (n = ", length(ids), " samples) — Batch vs Biological MDS"),
    fontface = "bold", size = 18, x = 0.5, y = 0.5
  )

  cowplot::plot_grid(title_grob, body_grid, legends_row, ncol = 1, rel_heights = c(0.036, 0.896, 0.068))
}

for (ds in retained_cohorts) {
  fig_mds <- build_cohort_4col_mds(ds)
  out_pdf <- file.path("plots", paste0("Supp_fig_TMP_BatchMDS_", ds, ".pdf"))
  ggsave(out_pdf, fig_mds, width = 13.5, height = 21.5, device = "pdf")
  message("Wrote ", out_pdf)
}

message("ALL_7COHORT_BATCH_FIGURES_COMPLETE")
