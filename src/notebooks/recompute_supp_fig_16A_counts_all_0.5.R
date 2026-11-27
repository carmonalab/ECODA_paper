#!/usr/bin/env Rscript
# Recompute Supp Fig 16A (transformation screen) locally under the unified
# counts_all__0.5 zero-imputation policy from existing integer composition_counts
# in data/derived/access_v2/, and regenerate Supp_fig_16A_* and 16AB PDFs.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = root)

suppressPackageStartupMessages({
  source("src/utils/load_all_functions.R")
  source("src/utils/imports_worker_transzeroimp.R")
  source("src/utils/analysis_release.R")
  source("src/utils/evaluation_scores.R")
  source("src/notebooks/benchmark_analysis.R")
  source("src/5_run_benchmark_methods/benchmark_pipeline.R")
  source("src/5_run_benchmark_methods/benchmark_hpc_utils.R")
})

release_rds <- file.path(root, "data/derived/access_v2/analysis_release.rds")
release_md5 <- "fa146125e805593a70810229e47028ec"
out_dir <- file.path(root, "data/derived/transformation_counts_all_0.5_20261121")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
path_plots <- file.path(root, "plots")

policy <- transformation_analysis_policy()
stopifnot(identical(policy$zero_imputation, "counts_all__0.5"))

benchmark_control <- read_benchmark_analysis_control(release_rds, release_md5)
release_stub <- list(control = benchmark_control)

datasets_cfg <- read_datasets_json(view = "benchmark_analysis")
datasets_cfg <- datasets_cfg[vapply(datasets_cfg, function(x) isTRUE(x$use_for_benchmark), logical(1))]
datasets_cfg <- datasets_cfg[!grepl("^_", names(datasets_cfg))]
ds_names <- names(datasets_cfg)
stopifnot(length(ds_names) == 11L)

trans_rows <- benchmark_reviewed_screen_rows(release_stub, "trans", ds_names)
columns <- policy$score_columns

old_frames <- list()
new_frames <- list()
cmp_rows <- list()

for (j in seq_len(nrow(trans_rows))) {
  row <- trans_rows[j, , drop = FALSE]
  ds <- row$dataset[[1L]]
  item <- benchmark_load_screen_row(release_stub, row)
  old_table <- item$result
  ct_comps <- as.data.frame(attr(old_table, "composition_counts"), check.names = FALSE)
  labels <- attr(old_table, "sample_labels")
  stopifnot(!is.null(ct_comps), !is.null(labels), isTRUE(all(as.matrix(ct_comps) == floor(as.matrix(ct_comps)))))

  rds_path <- file.path(out_dir, paste0(ds, "_trans_counts_all_0.5.rds"))
  md5_path <- paste0(rds_path, ".md5")
  new_table <- NULL
  if (file.exists(rds_path) && file.exists(md5_path)) {
    recorded <- sub("^MD5=", "", grep("^MD5=", readLines(md5_path, warn = FALSE), value = TRUE)[1])
    if (identical(unname(tools::md5sum(rds_path)), recorded)) {
      cached <- readRDS(rds_path)[["__trans__"]]
      if (isTRUE(tryCatch({ validate_transformation_result(cached, expected_samples = rownames(ct_comps)); TRUE }, error = function(e) FALSE))) {
        new_table <- cached
      }
    }
  }
  if (is.null(new_table)) {
    new_table <- run_transformation_analysis(ct_comps, labels)
    validate_transformation_result(new_table, expected_samples = rownames(ct_comps))
    bundle <- list(`__trans__` = new_table)
    saveRDS(bundle, rds_path)
    md5_val <- unname(tools::md5sum(rds_path))
    writeLines(
      c(paste0("MD5=", md5_val), paste0("SIZE=", file.info(rds_path)$size), paste0("PATH=", rds_path)),
      md5_path
    )
  }

  for (i in seq_len(nrow(new_table))) {
    m <- as.character(new_table$trans_method[i])
    new_vals <- vapply(unname(columns), function(col) new_table[[col]][i], numeric(1))
    old_idx <- match(m, as.character(old_table$trans_method))
    old_vals <- vapply(unname(columns), function(col) old_table[[col]][old_idx], numeric(1))

    new_frames[[length(new_frames) + 1L]] <- data.frame(
      dataset = ds, method = m, score = names(columns), value = new_vals, stringsAsFactors = FALSE
    )
    old_frames[[length(old_frames) + 1L]] <- data.frame(
      dataset = ds, method = m, score = names(columns), value = old_vals, stringsAsFactors = FALSE
    )
    cmp_rows[[length(cmp_rows) + 1L]] <- data.frame(
      dataset = ds, trans_method = m, score = names(columns),
      old_counts_all_1 = old_vals, new_counts_all_0.5 = new_vals,
      diff = new_vals - old_vals, stringsAsFactors = FALSE
    )
  }
  cat(sprintf("[%02d/11] %s: OK (%d samples x %d CTs)\n", j, ds, nrow(ct_comps), ncol(ct_comps)))
}

df_cmp <- do.call(rbind, cmp_rows)
write.csv(df_cmp, file.path(out_dir, "transformation_scores_summary.csv"), row.names = FALSE)

agg_cmp <- aggregate(cbind(old_counts_all_1, new_counts_all_0.5, diff) ~ trans_method + score, df_cmp, mean)
cat("\n=== Mean Score Comparison (counts_all__1 vs counts_all__0.5, 11 datasets) ===\n")
print(agg_cmp[order(agg_cmp$score, -agg_cmp$new_counts_all_0.5), ], row.names = FALSE)

df_norm_raw <- benchmark_validate_screen_scores(do.call(rbind, new_frames))
df_normalization <- df_norm_raw %>%
  mutate(
    trans_method = recode(
      method,
      "clr" = "CLR",
      "alr_mincvref" = "ALR_minCVref",
      "alr_randref" = "ALR_randref",
      "arcsine_sqrt" = "Arcsine_sqrt",
      "counts" = "Counts",
      "freq" = "Frequency(%)"
    )
  )

df_normalization_anosim <- df_normalization %>%
  filter(score == "anosim_score") %>%
  group_by(trans_method) %>%
  summarize(mean = mean(value), sem = sd(value) / sqrt(length(value)), .groups = "drop") %>%
  arrange(desc(mean))

df_normalization_anosim$trans_method <- factor(
  df_normalization_anosim$trans_method,
  levels = df_normalization_anosim$trans_method
)
df_normalization$trans_method <- factor(
  df_normalization$trans_method,
  levels = df_normalization_anosim$trans_method
)

p1 <- df_normalization %>%
  filter(score == "anosim_score") %>%
  ggplot(aes(x = trans_method, y = value)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "bar", fill = "white", color = "black") +
  geom_point(position = "identity", alpha = 0.6, size = 0.8) +
  geom_line(aes(group = dataset), color = "grey20", alpha = 0.2) +
  stat_compare_means(
    ref.group = "CLR",
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    paired = TRUE,
    hide.ns = TRUE,
    label.y = 0.96,
    aes(label = after_stat(p.signif)),
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  coord_cartesian(ylim = c(-0.06, 1.05)) +
  labs(x = "", y = "Separation score (ANOSIM)") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1))

ggsave(file.path(path_plots, "Supp_fig_16A_transformations_anosim.pdf"), plot = p1, width = 2.5, height = 3.5)

transformation_unscaled_means <- benchmark_unscaled_transformation_means(df_normalization)
p_transformation_unscaled_mean <- (p1 + transformation_unscaled_means) +
  coord_cartesian(ylim = c(NA, max(transformation_unscaled_means$value) * 1.15)) +
  labs(y = "Separation score\n(mean, unscaled)") +
  theme(axis.title.y = element_text(size = 10))
ggsave(
  file.path(path_plots, "Supp_fig_16A_transformations_mean_unscaled.pdf"),
  plot = p_transformation_unscaled_mean, width = 3.5, height = 3.5
)

transformation_companion_scores <- df_normalization %>% mutate(method = as.character(trans_method))
transformation_companions <- lapply(setdiff(evaluation_metrics(), "anosim_score"), function(metric) {
  panel <- plot_benchmark_raw_metric(
    transformation_companion_scores, metric,
    method_order = levels(df_normalization$trans_method), reference = "CLR", monochrome = TRUE
  )
  ggsave(file.path(path_plots, paste0("Supp_fig_16A_transformations_", metric, ".pdf")),
         plot = panel, width = 5, height = 4)
  panel
})
ggsave(
  file.path(path_plots, "Supp_fig_16A_transformations_metric_companions.pdf"),
  plot = wrap_plots(transformation_companions, ncol = 2, guides = "collect"),
  width = 10, height = 8
)

# Regenerate Supp_fig_16AB_trans_zero_anosim.pdf combining updated p1 with extended p2
sweep_csv <- file.path(root, "data/derived/zeroimp_count_sweep_20261127/zeroimp_scores_all11.csv")
df_zeroimp <- if (file.exists(sweep_csv)) {
  read.csv(sweep_csv, stringsAsFactors = FALSE)
} else {
  benchmark_zeroimp_scores(release_stub, ds_names)
}
zeroimp_counts_order <- c("0.001", "0.01", "0.1", "0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")
zeroimp_order <- c(
  paste0("counts_all_", zeroimp_counts_order),
  paste0("counts_zeros_", zeroimp_counts_order),
  paste0("multLN_", c("0.001", "0.01", "0.1", "1", "2", "5"), "%")
)
zeroimp_companion_scores <- df_zeroimp %>%
  filter(method %in% zeroimp_order) %>%
  mutate(
    method = factor(method, levels = zeroimp_order),
    family = factor(
      sub("_(?:[^_]+)$", "", as.character(method), perl = TRUE),
      levels = c("counts_all", "counts_zeros", "multLN")
    )
  )
stopifnot(setequal(unique(as.character(zeroimp_companion_scores$method)), zeroimp_order))
benchmark_validate_screen_scores(zeroimp_companion_scores)

zeroimp_metric_panel <- function(metric) {
  rows <- zeroimp_companion_scores %>% filter(score == metric)
  reference <- rows[rows$method == "counts_all_0.5", c("dataset", "value")]
  tests <- bind_rows(lapply(setdiff(zeroimp_order, "counts_all_0.5"), function(method_key) {
    values <- rows[as.character(rows$method) == method_key, ]
    baseline <- reference$value[match(values$dataset, reference$dataset)]
    p_value <- if (all(values$value == baseline)) 1 else
      suppressWarnings(stats::wilcox.test(values$value, baseline, paired = TRUE)$p.value)
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
combined_16ab <- wrap_plots(list(p1, p2)) +
  plot_layout(nrow = 1, widths = c(1, 3.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(file.path(path_plots, "Supp_fig_16AB_trans_zero_anosim.pdf"),
       plot = combined_16ab, width = 11.5, height = 4.5)

cat("DONE: Supp_fig_16A and Supp_fig_16AB regenerated with counts_all__0.5\n")
