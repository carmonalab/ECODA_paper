#!/usr/bin/env Rscript
# Extend the zero-handling count-offset screen (counts_all_* and counts_zeros_*)
# downward to include 0.001, 0.01, and 0.1 across all 11 clean benchmark cohorts,
# verify exact numerical match on all 43 overlapping settings against data/derived/access_v2/,
# and regenerate Supp_fig_16B_* and Supp_fig_16AB_trans_zero_anosim.pdf.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = root)

suppressPackageStartupMessages({
  source("src/utils/load_all_functions.R")
  source("src/utils/imports_worker_transzeroimp.R")
  source("src/utils/analysis_release.R")
  source("src/utils/evaluation_scores.R")
  source("src/5_run_benchmark_methods/run_local_ecoda_derived.R")
  source("src/notebooks/benchmark_analysis.R")
  source("src/5_run_benchmark_methods/benchmark_pipeline.R")
  source("src/notebooks/run_zeroimp_count_screen.R")
  library(ggpubr)
})

release_rds <- file.path(root, "data/derived/access_v2/analysis_release.rds")
release_md5 <- "fa146125e805593a70810229e47028ec"
control <- read_benchmark_analysis_control(release_rds, release_md5)
release_stub <- list(control = control)

out_dir <- file.path(root, "data/derived/zeroimp_count_sweep_20261127")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
path_plots <- file.path(root, "plots")
dir.create(path_plots, recursive = TRUE, showWarnings = FALSE)

datasets_cfg <- read_datasets_json(view = "benchmark_analysis")
datasets_cfg <- datasets_cfg[vapply(datasets_cfg, function(x) isTRUE(x$use_for_benchmark), logical(1))]
datasets_cfg <- datasets_cfg[!grepl("^_", names(datasets_cfg))]
ds_names <- names(datasets_cfg)
stopifnot(length(ds_names) == 11L)

trans_rows <- benchmark_reviewed_screen_rows(release_stub, "trans", ds_names)
zi_rows <- benchmark_reviewed_screen_rows(release_stub, "zeroimp", ds_names)
metrics <- evaluation_metrics()
menu <- zeroimp_menu()
stopifnot(length(menu) == 49L)

all_frames <- list()
overlap_diffs <- numeric()

for (j in seq_along(ds_names)) {
  ds <- ds_names[j]
  t_row <- trans_rows[trans_rows$dataset == ds, , drop = FALSE]
  z_row <- zi_rows[zi_rows$dataset == ds, , drop = FALSE]
  t_item <- benchmark_load_screen_row(release_stub, t_row)
  z_item <- benchmark_load_screen_row(release_stub, z_row)

  ct_comps <- as.matrix(attr(t_item$result, "composition_counts"))
  t_labels <- attr(t_item$result, "sample_labels")
  z_labels <- z_item$result$labels

  stopifnot(
    !is.null(ct_comps), !is.null(t_labels), !is.null(z_labels),
    all(ct_comps >= 0), all(ct_comps == floor(ct_comps)),
    setequal(rownames(ct_comps), names(z_labels)),
    identical(as.character(t_labels[names(z_labels)]), as.character(z_labels))
  )

  # Match exact sample ordering, double storage mode, and labels of the reviewed __zeroimp__ input
  ct_comps <- ct_comps[names(z_labels), , drop = FALSE]
  storage.mode(ct_comps) <- "double"
  stopifnot(
    identical(.ecoda_screen_digest(ct_comps), z_item$result$counts_md5),
    identical(.ecoda_screen_digest(z_labels), z_item$result$labels_md5)
  )
  input <- list(counts = ct_comps, labels = z_labels)

  rds_path <- file.path(out_dir, paste0(ds, "_zeroimp_counts.rds"))
  csv_path <- file.path(out_dir, paste0(ds, "_zeroimp_counts.csv"))
  md5_path <- paste0(rds_path, ".md5")

  pkg <- NULL
  if (file.exists(rds_path) && file.exists(md5_path)) {
    recorded <- sub("^MD5=", "", grep("^MD5=", readLines(md5_path, warn = FALSE), value = TRUE)[1])
    if (identical(unname(tools::md5sum(rds_path)), recorded)) {
      cached <- readRDS(rds_path)
      if (identical(cached$results$setting, menu) &&
          identical(cached$counts_md5, .ecoda_screen_digest(ct_comps)) &&
          identical(cached$labels_md5, .ecoda_screen_digest(z_labels))) {
        pkg <- cached
      }
    }
  }

  if (is.null(pkg)) {
    pkg <- zeroimp_run_cohort(input)
    pkg$schema <- zeroimp_policy()$schema
    pkg$dataset <- ds
    pkg$run_id <- "zeroimp_count_sweep_20261127"
    write.csv(pkg$results, csv_path, row.names = FALSE)
    saveRDS(pkg, rds_path)
    writeLines(
      c(paste0("MD5=", unname(tools::md5sum(rds_path))),
        paste0("SIZE=", file.info(rds_path)$size),
        paste0("PATH=", rds_path)),
      md5_path
    )
  }

  # Mandatory audit: exact match on all 43 overlapping settings vs access_v2
  # (4 deterministic/RNG-exact metrics < 1e-12; lisi_score < 1e-6 for x86_64 vs ARM64 C++ perplexity rounding)
  old_res <- z_item$result$results
  common_settings <- old_res$setting
  new_sub <- pkg$results[match(common_settings, pkg$results$setting), , drop = FALSE]
  stopifnot(identical(new_sub$setting, old_res$setting),
            identical(new_sub$status, old_res$status))
  pass_idx <- which(old_res$status == "PASS")
  for (m in metrics) {
    diff_m <- abs(new_sub[[m]][pass_idx] - old_res[[m]][pass_idx])
    tol <- if (m == "lisi_score") 1e-6 else 1e-12
    stopifnot(all(is.finite(diff_m)), max(diff_m) < tol)
    overlap_diffs <- c(overlap_diffs, max(diff_m))
  }
  new_only <- paste0(rep(c("counts_all_", "counts_zeros_"), each = 3L), c("0.001", "0.01", "0.1"))
  new_sub2 <- pkg$results[match(new_only, pkg$results$setting), , drop = FALSE]
  stopifnot(all(new_sub2$status == "PASS"))

  for (m in metrics) {
    all_frames[[length(all_frames) + 1L]] <- data.frame(
      dataset = ds,
      method = pkg$results$setting,
      score = m,
      value = pkg$results[[m]],
      status = pkg$results$status,
      reason = pkg$results$reason,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("[%02d/11] %-26s OK (49 settings, max overlap diff = %.2e)\n",
              j, ds, max(overlap_diffs)))
}

cat(sprintf("\nOVERLAP AUDIT PASSED: max abs diff across all 11 cohorts x 43 settings x 5 metrics = %.3e\n",
            max(overlap_diffs)))

df_zeroimp_all <- do.call(rbind, all_frames)
write.csv(df_zeroimp_all, file.path(out_dir, "zeroimp_scores_all11.csv"), row.names = FALSE)

# Build Supp Fig 16B using exact notebook logic with extended zeroimp_counts_order
zeroimp_counts_order <- c("0.001", "0.01", "0.1", "0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")
zeroimp_order <- c(
  paste0("counts_all_", zeroimp_counts_order),
  paste0("counts_zeros_", zeroimp_counts_order),
  paste0("multLN_", c("0.001", "0.01", "0.1", "1", "2", "5"), "%")
)

# Verify every displayed method PASSed in all 11 cohorts
disp_rows <- df_zeroimp_all %>% filter(method %in% zeroimp_order)
stopifnot(all(disp_rows$status == "PASS"),
          all(is.finite(disp_rows$value)))

zeroimp_companion_scores <- disp_rows %>%
  mutate(
    method = factor(method, levels = zeroimp_order),
    family = factor(
      sub("_(?:[^_]+)$", "", as.character(method), perl = TRUE),
      levels = c("counts_all", "counts_zeros", "multLN")
    )
  )
stopifnot(setequal(unique(as.character(zeroimp_companion_scores$method)), zeroimp_order),
          !anyNA(zeroimp_companion_scores$family))
benchmark_validate_screen_scores(zeroimp_companion_scores)

# Print & save summary statistics (mean, SEM, two-sided paired Wilcoxon p vs counts_all_0.5)
anosim_rows <- zeroimp_companion_scores %>% filter(score == "anosim_score")
ref_anosim <- anosim_rows[anosim_rows$method == "counts_all_0.5", c("dataset", "value")]
summary_list <- lapply(zeroimp_order, function(m_key) {
  vals <- anosim_rows[as.character(anosim_rows$method) == m_key, ]
  base <- ref_anosim$value[match(vals$dataset, ref_anosim$dataset)]
  stopifnot(identical(vals$dataset, ref_anosim$dataset[match(vals$dataset, ref_anosim$dataset)]))
  p_val <- if (m_key == "counts_all_0.5" || all(vals$value == base)) 1.0 else
    suppressWarnings(stats::wilcox.test(vals$value, base, paired = TRUE)$p.value)
  data.frame(
    method = m_key,
    family = as.character(vals$family[[1L]]),
    mean_anosim = mean(vals$value),
    sem_anosim = sd(vals$value) / sqrt(length(vals$value)),
    p_vs_counts_all_0.5 = p_val,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "zeroimp_summary_stats.csv"), row.names = FALSE)
cat("\n=== Extended Zero-Handling ANOSIM Summary (11 cohorts) ===\n")
print(summary_df, digits = 4, row.names = FALSE)

zeroimp_metric_panel <- function(metric) {
  rows <- zeroimp_companion_scores %>% filter(score == metric)
  reference <- rows[rows$method == "counts_all_0.5", c("dataset", "value")]
  tests <- bind_rows(lapply(setdiff(zeroimp_order, "counts_all_0.5"), function(method_key) {
    values <- rows[as.character(rows$method) == method_key, ]
    baseline <- reference$value[match(values$dataset, reference$dataset)]
    stopifnot(length(baseline) == nrow(reference), !anyNA(baseline),
              !anyDuplicated(values$dataset), setequal(values$dataset, reference$dataset))
    p_value <- if (all(values$value == baseline)) 1 else
      suppressWarnings(stats::wilcox.test(values$value, baseline, paired = TRUE)$p.value)
    stopifnot(is.finite(p_value))
    data.frame(
      method = method_key, family = as.character(values$family[[1]]),
      label = as.character(cut(p_value, c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                               labels = c("****", "***", "**", "*", "")))
    )
  }))
  tests$y <- max(rows$value) + 0.08 * max(diff(range(rows$value)), 0.1)
  plot_benchmark_raw_metric(zeroimp_companion_scores, metric, method_order = zeroimp_order, monochrome = TRUE) +
    facet_grid(cols = vars(family), scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = function(x) sub("^(counts_all_|counts_zeros_|multLN_)", "", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    geom_text(data = tests, aes(x = method, y = y, label = label), inherit.aes = FALSE, size = 3) +
    labs(y = paste0("Separation score (", score_label_map[[metric]], ")")) +
    theme(legend.position = "none")
}

p2 <- zeroimp_metric_panel("anosim_score")
ggsave(file.path(path_plots, "Supp_fig_16B_zero_replacement_anosim.pdf"),
       plot = p2, width = 13, height = 5)

zeroimp_companions <- lapply(setdiff(evaluation_metrics(), "anosim_score"), function(metric) {
  panel <- zeroimp_metric_panel(metric)
  ggsave(file.path(path_plots, paste0("Supp_fig_16B_zero_replacement_", metric, ".pdf")),
         plot = panel, width = 13, height = 5)
  panel
})
ggsave(file.path(path_plots, "Supp_fig_16B_zero_replacement_metric_companions.pdf"),
       plot = wrap_plots(zeroimp_companions, ncol = 2, guides = "collect"),
       width = 26, height = 10)

# Load Panel A (transformation screen under unified counts_all__0.5) from verified cache
trans_dir <- file.path(root, "data/derived/transformation_counts_all_0.5_20261121")
columns <- transformation_analysis_policy()$score_columns
trans_frames <- list()
for (ds in ds_names) {
  rds_p <- file.path(trans_dir, paste0(ds, "_trans_counts_all_0.5.rds"))
  md5_p <- paste0(rds_p, ".md5")
  rec <- sub("^MD5=", "", grep("^MD5=", readLines(md5_p, warn = FALSE), value = TRUE)[1])
  stopifnot(identical(unname(tools::md5sum(rds_p)), rec))
  tbl <- readRDS(rds_p)[["__trans__"]]
  for (i in seq_len(nrow(tbl))) {
    m <- as.character(tbl$trans_method[i])
    vals <- vapply(unname(columns), function(col) tbl[[col]][i], numeric(1))
    trans_frames[[length(trans_frames) + 1L]] <- data.frame(
      dataset = ds, method = m, score = names(columns), value = vals, stringsAsFactors = FALSE
    )
  }
}
df_norm_raw <- benchmark_validate_screen_scores(do.call(rbind, trans_frames))
df_normalization <- df_norm_raw %>%
  mutate(trans_method = recode(method,
    "clr" = "CLR", "alr_mincvref" = "ALR_minCVref", "alr_randref" = "ALR_randref",
    "arcsine_sqrt" = "Arcsine_sqrt", "counts" = "Counts", "freq" = "Frequency(%)"))
df_norm_anosim <- df_normalization %>%
  filter(score == "anosim_score") %>%
  group_by(trans_method) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  arrange(desc(mean))
df_normalization$trans_method <- factor(df_normalization$trans_method, levels = df_norm_anosim$trans_method)

p1 <- df_normalization %>%
  filter(score == "anosim_score") %>%
  ggplot(aes(x = trans_method, y = value)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "bar", fill = "white", color = "black") +
  geom_point(position = "identity", alpha = 0.6, size = 0.8) +
  geom_line(aes(group = dataset), color = "grey20", alpha = 0.2) +
  stat_compare_means(
    ref.group = "CLR", method = "wilcox.test", method.args = list(exact = FALSE),
    paired = TRUE, hide.ns = TRUE, label.y = 0.96,
    aes(label = after_stat(p.signif)), size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  coord_cartesian(ylim = c(-0.06, 1.05)) +
  labs(x = "", y = "Separation score (ANOSIM)") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1))

combined_16ab <- wrap_plots(list(p1, p2)) +
  plot_layout(nrow = 1, widths = c(1, 3.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(file.path(path_plots, "Supp_fig_16AB_trans_zero_anosim.pdf"),
       plot = combined_16ab, width = 11.5, height = 4.5)

cat("DONE: Extended Supp_fig_16B_* and Supp_fig_16AB_trans_zero_anosim.pdf regenerated.\n")
