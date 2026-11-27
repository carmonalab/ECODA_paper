#!/usr/bin/env Rscript
# Fast standalone regeneration of Figure 4 (panels A-F) from core_scores.csv
# and the 4 required Adams/Kfoury bundles in data/derived/access_v2/, without
# calling load_benchmark_analysis_release() across all 1,100 bundles.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = root)

suppressPackageStartupMessages({
  source("src/utils/load_all_functions.R")
  source("src/utils/analysis_release.R")
  source("src/utils/evaluation_scores.R")
  source("src/notebooks/benchmark_analysis.R")
  library(ggpubr)
  library(gtools)
})

path_plots <- file.path(root, "plots")
dir.create(path_plots, recursive = TRUE, showWarnings = FALSE)

# 1. Load core scores for panels A & B
core_csv <- file.path(root, "data/derived/benchmark_reporting_20261110_v2/core_scores.csv")
stopifnot(file.exists(core_csv))
df_benchmark <- read.csv(core_csv, stringsAsFactors = FALSE)
df_benchmark$dataset[df_benchmark$dataset == "Gongsharma_cmv_young_males"] <- "GongSharma"
benchmark_metrics <- evaluation_metrics()

df_benchmark_sub <- df_benchmark %>%
  filter(score %in% benchmark_metrics) %>%
  filter(grepl("ECODA_authors_HR", method) & !grepl("HiTME|scATOMIC|PCA|NULL|most|least", method)) %>%
  filter(!grepl("cmv|sex", dataset)) %>%
  filter(!method == "ECODA_authors_HR_top_varexp0") %>%
  mutate(method = recode(method,
                         "ECODA_authors_HR" = "100",
                         "ECODA_authors_HR_top_varexp0" = "2most_varcts")) %>%
  mutate(method = gsub("ECODA_authors_HR", "", method)) %>%
  mutate(method = if_else(
    grepl("^_top_varexp0\\.\\d", method),
    paste0(as.integer(gsub("_top_varexp0\\.", "", method)) * 10),
    method
  )) %>%
  mutate(method = fct_relevel(method)) %>%
  mutate(method = factor(method, levels = rev(mixedsort(levels(method))))) %>%
  mutate(dataset = recode(dataset, "GongSharma_all_age" = "GongSharma"))

df_benchmark_anosim <- df_benchmark_sub %>% filter(score == "anosim_score")
stopifnot(nrow(df_benchmark_anosim) == 110L,
          length(unique(df_benchmark_anosim$dataset)) == 11L,
          !anyNA(df_benchmark_anosim$n_features))

pointsize <- 1.5

p1 <- df_benchmark_anosim %>%
  ggplot(aes(x = method, y = n_features)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "bar", fill = "white", color = "black") +
  geom_point(aes(shape = dataset, color = dataset), size = pointsize) +
  scale_shape_manual(values = 1:length(unique(df_benchmark_sub$dataset))) +
  geom_line(aes(group = dataset, linetype = dataset, color = dataset), alpha = 1) +
  labs(x = "Variance explained (%)", y = "Number of cell types") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1),
        legend.position = "none")

p2 <- df_benchmark_anosim %>%
  ggplot(aes(x = method, y = value)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "bar", fill = "white", color = "black") +
  geom_point(aes(shape = dataset, color = dataset), size = pointsize) +
  scale_shape_manual(values = 1:length(unique(df_benchmark_sub$dataset))) +
  geom_line(aes(group = dataset, linetype = dataset, color = dataset), alpha = 1) +
  labs(x = "Variance explained (%)", y = "Separation score\n(ANOSIM)") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1),
        legend.position = "none")

figure_4_ab <- list(p1, p2)

# 2. Load ONLY the 4 needed Adams/Kfoury bundles from checksum-verified access_v2
release_rds <- file.path(root, "data/derived/access_v2/analysis_release.rds")
release_md5 <- "fa146125e805593a70810229e47028ec"
control <- read_benchmark_analysis_control(release_rds, release_md5)

load_bundle <- function(ds, method_key) {
  row <- control$ledger[control$ledger$dataset == ds &
                        control$ledger$view == "benchmark_analysis" &
                        control$ledger$method_key == method_key, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  path <- .analysis_checked_file(row$artifact_path[[1L]], row$md5[[1L]])
  obj <- readRDS(path)
  bundle <- obj[[row$bundle_key[[1L]]]]
  stopifnot(!is.null(bundle), !is.null(bundle$feat_mat), !is.null(bundle$labels))
  bundle
}

build_pca_and_boxplot <- function(ds) {
  b_hr <- load_bundle(ds, "ECODA_authors_HR")
  b_top <- load_bundle(ds, "ECODA_authors_HR_top_varexp0")
  feat_mat_recalc <- b_top$feat_mat
  feat_mat <- b_hr$feat_mat[, colnames(feat_mat_recalc), drop = FALSE]
  labels <- b_hr$labels

  df_logr <- as.data.frame(feat_mat_recalc[, 1, drop = FALSE] * 2)
  df_logr$Groups <- labels

  p_pca <- plot_pca(
    feat_mat = feat_mat,
    labels = labels,
    cluster_score = FALSE, mod_score = FALSE, anosim_score = FALSE,
    coord_equal = TRUE,
    title = ds,
    n_ct_show = 4, invisible = "quali", repel = TRUE, pointsize = 2
  )

  p_box <- ggplot(df_logr, aes(x = Groups, y = !!sym(colnames(df_logr)[1]), fill = Groups)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = NULL,
      y = paste0("Log ratio\n", colnames(feat_mat_recalc)[1], " / ", colnames(feat_mat_recalc)[2])
    ) &
    NoLegend()

  list(pca = p_pca, box = p_box)
}

adams <- build_pca_and_boxplot("Adams")
kfoury <- build_pca_and_boxplot("Kfoury")

figure_4_c <- adams$pca
figure_4_e <- adams$box
figure_4_d <- kfoury$pca
figure_4_f <- kfoury$box

# Combine with smaller overall dimensions so labels appear larger naturally,
# enlarged middle PCA panels (C & D) with preserved coord_equal 1:1 PC aspect
# ratio, and narrower boxplot panels (E & F).
figure_4_left <- wrap_plots(list(figure_4_ab[[1]] + labs(tag = "A"),
  figure_4_ab[[2]] + labs(tag = "B")), ncol = 1, guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 6, 4, 6)) &
  labs(colour = NULL, shape = NULL, linetype = NULL) &
  guides(colour = guide_legend(ncol = 3), shape = guide_legend(ncol = 3),
    linetype = guide_legend(ncol = 3))

figure_4_middle <- wrap_plots(list(figure_4_c + labs(tag = "C"),
  figure_4_d + labs(tag = "D")), ncol = 1, guides = "keep") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 9),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0),
        plot.margin = margin(4, 6, 4, 6)) &
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

figure_4_right <- wrap_plots(list(figure_4_e + labs(tag = "E"),
  figure_4_f + labs(tag = "F")), ncol = 1) &
  theme(legend.position = "none",
        plot.margin = margin(4, 6, 4, 6))

figure_4 <- wrap_plots(list(figure_4_left, figure_4_middle, figure_4_right),
  nrow = 1, widths = c(1.05, 1.45, 0.65)) &
  theme(plot.tag = element_text(size = 13, face = "bold"))

out_pdf <- file.path(path_plots, "Figure_4.pdf")
ggsave(out_pdf, plot = figure_4, width = 10.0, height = 5.6)
cat("Wrote:", out_pdf, "\n")
