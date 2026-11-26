#!/usr/bin/env Rscript
# Refresh Stephenson 53-sample transcriptional complementarity (CD14+ pseudobulk
# HVG2000 PCA + IFN_UCell authors_highres PCA) without label leakage and plot
# plots/Supp_fig_TMP_Stephenson_IFN.pdf.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = root)

suppressPackageStartupMessages({
  source("src/utils/load_all_functions.R")
  source("src/utils/analysis_release.R")
  source("src/utils/evaluation_scores.R")
  source("src/notebooks/stephenson_ucell.R")
})

cd14_dir <- file.path(root, "data/derived/stephenson_cd14_pseudobulk_20261121")
ucell_dir <- file.path(root, "data/derived/stephenson_ucell_20261121")
path_plots <- file.path(root, "plots")

# 1. Verify exact 53 samples & labels against access_v2
rel <- readRDS(file.path(root, "data/derived/access_v2/analysis_release.rds"))
row_hr <- rel$ledger[rel$ledger$dataset == "Stephenson" & rel$ledger$method_key == "ECODA_authors_HR", ][1, ]
bundle_hr <- readRDS(row_hr$artifact_path)[[row_hr$bundle_key]]

samples_df <- read.csv(file.path(cd14_dir, "samples_53.csv"), stringsAsFactors = FALSE)
samples_df <- samples_df[order(samples_df$Sample), , drop = FALSE]
stopifnot(nrow(samples_df) == 53L, setequal(samples_df$Sample, names(bundle_hr$labels)))
labels <- factor(samples_df$Status, levels = c("Covid", "Healthy"))
names(labels) <- samples_df$Sample
stopifnot(identical(as.character(labels[names(bundle_hr$labels)]), as.character(bundle_hr$labels)))

# 2. Label-blind DESeq2 VST + top-2000 variable genes on 53-sample CD14 raw counts
pb_raw <- read.csv(file.path(cd14_dir, "cd14_raw_counts_53.csv"), row.names = 1, check.names = FALSE)
pb_raw <- pb_raw[samples_df$Sample, , drop = FALSE]
counts_g_by_s <- t(as.matrix(pb_raw))
storage.mode(counts_g_by_s) <- "integer"

# Metadata passed to DESeq2 contains Sample ONLY (zero Status, zero Site -> no label leakage)
pb_meta <- data.frame(Sample = samples_df$Sample, row.names = samples_df$Sample, stringsAsFactors = FALSE)
cd14_norm <- get_pb_deseq2_from_counts(
  counts = counts_g_by_s,
  metadata = pb_meta,
  hvg = NULL,
  n_hvg = 2000L,
  black_list = "none",
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
)
stopifnot(identical(rownames(cd14_norm), samples_df$Sample), ncol(cd14_norm) == 2000L, all(is.finite(cd14_norm)))

expr_out <- cbind(data.frame(Sample = rownames(cd14_norm), check.names = FALSE), as.data.frame(cd14_norm, check.names = FALSE))
samp_out <- data.frame(Sample = samples_df$Sample, Status = samples_df$Status, n_cells = as.integer(samples_df$cd14_n_cells), stringsAsFactors = FALSE)
arrow::write_feather(expr_out, file.path(cd14_dir, "expression.feather"))
arrow::write_feather(samp_out, file.path(cd14_dir, "samples.feather"))

# 3. Prepare 53-sample UCell IFN authors_highres features (complete nonconstant CTs)
ifn_raw <- read.csv(file.path(ucell_dir, "authors_highres__IFN_UCell_53.csv"), row.names = 1, check.names = FALSE)
ifn_raw <- as.matrix(ifn_raw[samples_df$Sample, , drop = FALSE])
ifn_prep <- prepare_stephenson_ucell_features(ifn_raw, require_pca = TRUE)
cat(sprintf("UCell IFN authors_highres: %d observed CTs -> %d retained complete nonconstant CTs (dropped incomplete: %s)\n",
            ncol(ifn_raw), ncol(ifn_prep$features), paste(ifn_prep$dropped_incomplete, collapse = ", ")))

# Save compact Feather + manifest
ifn_feather_df <- cbind(data.frame(Sample = rownames(ifn_raw), check.names = FALSE), as.data.frame(ifn_raw, check.names = FALSE))
arrow::write_feather(ifn_feather_df, file.path(ucell_dir, "authors_highres__IFN_UCell.feather"))
arrow::write_feather(samp_out[, c("Sample", "Status")], file.path(ucell_dir, "samples.feather"))

# 4. Score both representations + ECODA_authors_HR across all 5 canonical metrics
score_all <- function(name, feat_mat) {
  d <- stats::dist(feat_mat)
  set.seed(20261121L)
  sc <- calc_sep_score(d, labels[rownames(feat_mat)], parallel = 1L)
  data.frame(
    representation = name,
    n_samples = nrow(feat_mat),
    n_features = ncol(feat_mat),
    ANOSIM = unname(sc$anosim_score),
    Modularity_k3 = unname(sc$mod_knn3_score),
    ARI = unname(sc$cluster_score),
    Silhouette = unname(sc$sil_score),
    LISI = unname(sc$lisi_score),
    stringsAsFactors = FALSE
  )
}

scores_df <- rbind(
  score_all("CD14_LR_pseudobulk_HVG2000", cd14_norm),
  score_all("Authors_HR_IFN_UCell", ifn_prep$features),
  score_all("ECODA_authors_HR", bundle_hr$feat_mat[samples_df$Sample, , drop = FALSE])
)
print(scores_df, row.names = FALSE)
write.csv(scores_df, file.path(cd14_dir, "stephenson_53_complementarity_scores.csv"), row.names = FALSE)

# 5. Plot single-panel CD14+ monocytes pseudobulk PCA
comp_cd14 <- stephenson_pca_comparison(
  cd14_norm, labels,
  title = "CD14+ monocytes pseudobulk",
  n_loadings = 8L, seed = 20261022L
)

# Inspect exact top loadings for text reporting
pca_cd14 <- prcomp(cd14_norm, scale. = FALSE)
cd14_top <- names(sort(rowSums(factoextra::get_pca_var(pca_cd14)$contrib[, 1:2]), decreasing = TRUE))[1:8]
cat("Top CD14 HVG2000 PC1-2 loadings:", paste(cd14_top, collapse = ", "), "\n")

p_cd14_1 <- comp_cd14$plot +
  ggplot2::labs(title = "CD14+ monocytes pseudobulk", subtitle = NULL) +
  theme(
    plot.title = element_text(size = 11.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(path_plots, "Supp_fig_TMP_Stephenson_IFN.pdf"), plot = p_cd14_1, width = 5.8, height = 5.2)
ggsave(file.path(path_plots, "Supp_fig_TMP_Stephenson_IFN.png"), plot = p_cd14_1, width = 5.8, height = 5.2, dpi = 300)
cat("Wrote plots/Supp_fig_TMP_Stephenson_IFN.{pdf,png}\n")
