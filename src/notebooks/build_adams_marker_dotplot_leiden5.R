#!/usr/bin/env Rscript
# Build Supplementary Adams canonical marker-gene dotplot comparing Authors_HR,
# Leiden_res_5, HiTME (layer2), and Authors_LR across 305,701 Adams IPF cells.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

df_other <- read.csv(
  file.path(root, "data/derived/adams_marker_dotplot_20261120/canonical_markers_by_annotation.csv"),
  stringsAsFactors = FALSE
) %>% filter(strategy %in% c("Authors_HR", "HiTME", "Authors_LR"))

df_l5 <- read.csv(
  file.path(root, "data/derived/adams_leiden5_20261121/leiden5_canonical_markers.csv"),
  stringsAsFactors = FALSE
) %>% mutate(cluster = as.character(cluster))
df_all <- bind_rows(df_other, df_l5)

# Data-verified correspondence table for Leiden_res_5 cluster labels
corr <- read.csv(
  file.path(root, "data/derived/adams_annotation_correspondence_20261023/counts.csv"),
  row.names = 1, check.names = FALSE
)
l5_selected <- c("69", "70", "76", "79", "57", "75", "56", "77", "55", "54", "86",
                 "47", "17", "44", "60", "40", "8", "7", "88", "27", "65", "80")

format_l5_label <- function(cl) {
  vals <- sort(stats::setNames(as.numeric(corr[, cl]), rownames(corr)), decreasing = TRUE)
  tot <- sum(vals)
  p1 <- round(100 * vals[[1]] / tot)
  top1 <- names(vals)[1]
  # For endothelial/mural clusters merging close subtypes, note top 2 subtypes
  if (p1 < 75 && length(vals) >= 2 && vals[[2]] / tot >= 0.20) {
    p2 <- round(100 * vals[[2]] / tot)
    return(sprintf("cl_%s (%s %d%% / %s %d%%)", cl, top1, p1, names(vals)[2], p2))
  }
  sprintf("cl_%s (%s %d%%)", cl, top1, p1)
}
l5_label_map <- setNames(vapply(l5_selected, format_l5_label, character(1)), l5_selected)

authors_hr_selected <- c(
  "VE_Peribronchial", "VE_Arterial", "VE_Venous", "VE_Capillary_B", "VE_Capillary_A",
  "ATI", "ATII", "Aberrant_Basaloid", "Club", "Basal", "Ciliated",
  "Myofibroblast", "Fibroblast", "SMC", "Pericyte",
  "Macrophage_Alveolar", "Macrophage", "cMonocyte", "pDC", "T_Cytotoxic", "B", "B_Plasma"
)
hitme_selected <- c(
  "Endothelial", "Epithelial", "Fibroblast",
  "Mac_MARCO", "Mac_SPP1", "Mono_THBS1", "Mono_CD14", "pDC",
  "CD8.CM", "NK", "Bcell", "PlasmaCell"
)
authors_lr_selected <- c("Endothelial", "Epithelial", "Stromal", "Myeloid", "Lymphoid")

selected_list <- list(
  Authors_HR = setNames(authors_hr_selected, authors_hr_selected),
  Leiden_res_5 = l5_label_map,
  HiTME = setNames(hitme_selected, hitme_selected),
  Authors_LR = setNames(authors_lr_selected, authors_lr_selected)
)

gene_groups <- list(
  Endothelial = c("PECAM1", "VWF", "ACKR1", "COL15A1", "PLVAP", "VWA1", "FCN3", "EDNRB", "CA4", "GJA5"),
  Epithelial  = c("EPCAM", "AGER", "SFTPC", "SFTPB", "ABCA3", "SCGB1A1", "SCGB3A2", "FOXJ1", "KRT5", "KRT17", "TP63", "MMP7"),
  Mesenchymal = c("COL1A1", "COL1A2", "DCN", "LUM", "CTHRC1", "POSTN", "ACTA2", "MYH11", "PDGFRB", "RGS5"),
  Myeloid     = c("PTPRC", "MARCO", "FABP4", "SPP1", "CD14", "FCGR3A", "LILRA4"),
  Lymphoid    = c("CD3D", "CD8A", "NKG7", "MS4A1", "JCHAIN")
)
genes_ordered <- unlist(gene_groups, use.names = FALSE)
gene_lineage <- rep(names(gene_groups), lengths(gene_groups))
names(gene_lineage) <- genes_ordered

df_sub <- bind_rows(lapply(names(selected_list), function(strat) {
  lmap <- selected_list[[strat]]
  sub <- df_all %>%
    filter(strategy == strat, cluster %in% names(lmap), gene %in% genes_ordered)
  stopifnot(setequal(unique(sub$cluster), names(lmap)))
  sub %>%
    mutate(
      cluster_label = unname(lmap[as.character(cluster)]),
      cluster_key = paste0(strat, "|", cluster_label),
      strategy_label = factor(
        strat,
        levels = c("Authors_HR", "Leiden_res_5", "HiTME", "Authors_LR"),
        labels = c("Authors_HR (22 of 38)", "Leiden_res_5 (22 of 92 clusters)", "HiTME layer2 (12 of 44)", "Authors_LR (5 of 5)")
      )
    )
}))

ordered_keys <- unlist(lapply(names(selected_list), function(strat) {
  paste0(strat, "|", unname(selected_list[[strat]]))
}), use.names = FALSE)

# Z-score mean_expr across selected clusters within each strategy & gene (capped [-2, 2.5])
df_sub <- df_sub %>%
  group_by(strategy, gene) %>%
  mutate(
    mu = mean(mean_expr),
    sdv = sd(mean_expr),
    z_expr = ifelse(sdv > 1e-6, (mean_expr - mu) / sdv, 0),
    z_expr_capped = pmax(-2.0, pmin(2.5, z_expr))
  ) %>%
  ungroup() %>%
  mutate(
    cluster_key = factor(cluster_key, levels = rev(ordered_keys)),
    gene = factor(gene, levels = genes_ordered),
    lineage = factor(gene_lineage[as.character(gene)], levels = names(gene_groups))
  )

p <- ggplot(df_sub, aes(x = gene, y = cluster_key)) +
  geom_point(aes(size = pct_expr, color = z_expr_capped)) +
  scale_y_discrete(labels = function(x) sub("^[^|]+\\|", "", x)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2.0, 2.5),
    name = "Scaled mean\nexpr (z-score)"
  ) +
  scale_size_continuous(range = c(0.3, 4.2), limits = c(0, 100), breaks = c(10, 25, 50, 75, 100),
                        name = "Expressing\ncells (%)") +
  facet_grid(strategy_label ~ lineage, scales = "free", space = "free") +
  labs(
    x = "Canonical Adams et al. (2020) lineage and disease marker genes",
    y = NULL,
    title = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major = element_line(color = "grey94", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EAEEF3", color = "grey65"),
    strip.text.x = element_text(face = "bold", size = 8.5),
    strip.text.y = element_text(face = "bold", size = 8.5, angle = 0),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.5),
    axis.text.y = element_text(size = 7.5),
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

ggsave(file.path(root, "plots/Supp_fig_TMP_Adams_marker_dotplot_leiden5.pdf"), plot = p, width = 14.5, height = 10.5)
ggsave(file.path(root, "plots/Supp_fig_TMP_Adams_marker_dotplot_leiden5.png"), plot = p, width = 14.5, height = 10.5, dpi = 250)
cat("Wrote plots/Supp_fig_TMP_Adams_marker_dotplot_leiden5.{pdf,png}\n")
