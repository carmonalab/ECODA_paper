#!/usr/bin/env Rscript
# Regenerate Supp_fig_14A, Supp_fig_14B, and composite Supp_fig_14.{pdf,png}
# by overlaying the scoped 16-row Bamboo classic-full-H5AD profiling table onto
# data/derived/access_v2/execution_times.feather in memory (leaving access_v2
# immutable). Also emits 11-cohort summary statistics and log-log regression
# residual diagnostics.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")

Sys.setenv(
  PROJECT_ROOT = root,
  ECODA_BENCHMARK_RELEASE_RDS = file.path(root, "data/derived/access_v2/analysis_release.rds"),
  ECODA_BENCHMARK_RELEASE_MD5 = "fa146125e805593a70810229e47028ec"
)

source(file.path(root, "src/utils/load_all_functions.R"))
source(file.path(root, "src/utils/analysis_release.R"))
source(file.path(root, "src/utils/evaluation_scores.R"))
source(file.path(root, "src/notebooks/benchmark_analysis.R"))

profile_csv <- file.path(root, "data/derived/sf15_classic_profile_20261126/profile_16rows.csv")
if (!file.exists(profile_csv)) {
  stop("Missing 16-row classic profile CSV: ", profile_csv)
}
profile_16 <- utils::read.csv(profile_csv, stringsAsFactors = FALSE)
stopifnot(
  nrow(profile_16) == 16L,
  all(c("dataset", "method", "time_secs", "mem_GB") %in% names(profile_16)),
  all(is.finite(profile_16$time_secs) & profile_16$time_secs > 0),
  all(is.finite(profile_16$mem_GB) & profile_16$mem_GB > 0)
)
profile_16 <- normalize_exec_keys(profile_16)
stopifnot(!anyDuplicated(paste(profile_16$dataset, profile_16$method, sep = "\r")))

datasets <- read_datasets_json(file.path(root, "datasets.json"), view = "benchmark_analysis")
datasets <- datasets[vapply(datasets, function(x) isTRUE(x$use_for_benchmark), logical(1))]
datasets <- datasets[!grepl("^_", names(datasets))]
stopifnot(length(datasets) == 11L)

benchmark_control <- read_benchmark_analysis_control(
  Sys.getenv("ECODA_BENCHMARK_RELEASE_RDS"),
  Sys.getenv("ECODA_BENCHMARK_RELEASE_MD5")
)
benchmark_release <- load_benchmark_analysis_release(benchmark_control, datasets)
benchmark_release_times <- read_benchmark_execution_times(benchmark_control)

result_list <- benchmark_release$results
release_statistics <- benchmark_release$statistics
release_statistics$dataset[release_statistics$dataset == "Gongsharma_cmv_young_males"] <- "GongSharma"
n_cells_by_ds <- as.list(setNames(release_statistics$n_cells, release_statistics$dataset))
n_cells_lookup <- unlist(n_cells_by_ds)

feather_exec_times <- normalize_exec_keys(benchmark_release_times)
bundle_exec_times <- purrr::imap_dfr(result_list[["bmark"]], ~ {
  dataset_name <- .y
  purrr::imap_dfr(.x, function(bundle, method_name) {
    execution_time_rows_from_bundle(dataset_name, method_name, bundle)
  })
})
if (ncol(bundle_exec_times) == 0L) {
  bundle_exec_times <- data.frame(
    dataset = character(),
    method = character(),
    time_secs = numeric(),
    mem_GB = numeric(),
    stringsAsFactors = FALSE
  )
}
bundle_exec_times <- normalize_exec_keys(bundle_exec_times)
feather_exec_times$timing_source <- "reviewed Feather"
bundle_exec_times$timing_source <- "selected result bundle"

feather_exec_times[[".execution_union_key"]] <- execution_union_key(feather_exec_times)
bundle_exec_times[[".execution_union_key"]] <- execution_union_key(bundle_exec_times)

bundle_only <- bundle_exec_times %>%
  dplyr::anti_join(
    feather_exec_times %>% dplyr::select(.execution_union_key),
    by = ".execution_union_key"
  )
exec_times <- dplyr::bind_rows(
  feather_exec_times %>% dplyr::select(-.execution_union_key),
  bundle_only %>% dplyr::select(-.execution_union_key)
) %>%
  deduplicate_exec_times()
exec_times <- benchmark_complete_execution_memory(exec_times, bundle_exec_times)
exec_times$n_cells <- n_cells_lookup[exec_times$dataset]
exec_times <- exec_times[exec_times$dataset %in% names(n_cells_by_ds), ]

exec_times_method_cost <- benchmark_method_execution_costs(exec_times, bundle_exec_times)

# Overlay the 16 classic full-H5AD measurements in memory
overlay_keys <- paste(profile_16$dataset, profile_16$method, sep = "\r")
cost_keys <- paste(exec_times_method_cost$dataset, exec_times_method_cost$method, sep = "\r")
idx <- match(overlay_keys, cost_keys)
if (anyNA(idx)) {
  stop("Profiled rows missing from execution table: ", paste(overlay_keys[is.na(idx)], collapse = ", "))
}
exec_times_method_cost$time_secs[idx] <- as.numeric(profile_16$time_secs)
exec_times_method_cost$mem_GB[idx] <- as.numeric(profile_16$mem_GB)
exec_times_method_cost$timing_source[idx] <- "classic_full_h5ad_profile_20261126"
exec_times_method_cost$memory_source[idx] <- "classic_full_h5ad_profile_20261126"
exec_times_method_cost$cost_accounting[idx] <- "classic full-H5AD profile"

filter_methods <- c(
  "ECODA_authors_HR",
  "Pseudobulk_hvg2000",
  "MOFA_hvg2000_factors15",
  "scITD_hvg2000_factors5",
  "GloScope_hvg2000_pcadims30",
  "GloProp",
  "MrVI_hvg2000",
  "scPoli_hvg2000_dims15_highres",
  "PILOT_hvg2000",
  "QOT_hvg2000",
  "PILOT-GM-VAE_hvg2000"
)

method_recode <- c(
  "ECODA_authors_HR" = "ECODA",
  "Pseudobulk_hvg2000" = "Pseudobulk",
  "MOFA_hvg2000_factors15" = "MOFA",
  "scITD_hvg2000_factors5" = "scITD",
  "GloScope_hvg2000_pcadims30" = "GloScope",
  "GloProp" = "GloProp",
  "MrVI_hvg2000" = "MrVI",
  "scPoli_hvg2000_dims15_highres" = "scPoli",
  "PILOT_hvg2000" = "PILOT",
  "QOT_hvg2000" = "QOT",
  "PILOT-GM-VAE_hvg2000" = "PILOT-GM-VAE"
)

plot_data <- exec_times_method_cost %>%
  dplyr::filter(
    method %in% filter_methods,
    !.execution_is_shared_timing_method(method)
  ) %>%
  dplyr::mutate(
    raw_method = method,
    method = dplyr::recode(method, !!!method_recode)
  )

stopifnot(nrow(plot_data) == 121L)
ram_data <- plot_data %>% dplyr::filter(!is.na(mem_GB))
stopifnot(nrow(ram_data) == 121L)

# Log-log regression diagnostics across all 121 points
compute_residuals <- function(df, y_col) {
  split_df <- split(df, df$method)
  dplyr::bind_rows(lapply(split_df, function(sub) {
    fit <- stats::lm(log10(sub[[y_col]]) ~ log10(sub$n_cells))
    pred <- stats::predict(fit)
    sub$metric_col <- y_col
    sub$log10_pred <- pred
    sub$log10_resid <- log10(sub[[y_col]]) - pred
    sub$ratio_obs_to_pred <- sub[[y_col]] / (10^pred)
    sub$r_squared <- summary(fit)$r.squared
    sub$slope <- stats::coef(fit)[[2L]]
    sub
  }))
}

res_time <- compute_residuals(plot_data, "time_secs")
res_mem <- compute_residuals(ram_data, "mem_GB")
res_all <- dplyr::bind_rows(res_time, res_mem)
utils::write.csv(
  res_all %>%
    dplyr::select(
      metric_col, method, raw_method, dataset, n_cells, time_secs, mem_GB,
      timing_source, memory_source, log10_pred, log10_resid, ratio_obs_to_pred,
      slope, r_squared
    ),
  file.path(root, "data/derived/sf15_classic_profile_20261126/sf15_regression_residuals_121rows.csv"),
  row.names = FALSE
)

summary_11 <- plot_data %>%
  dplyr::group_by(method, raw_method) %>%
  dplyr::summarise(
    n_cohorts = dplyr::n(),
    mean_time_secs = mean(time_secs),
    median_time_secs = stats::median(time_secs),
    min_time_secs = min(time_secs),
    max_time_secs = max(time_secs),
    mean_mem_GB = mean(mem_GB),
    median_mem_GB = stats::median(mem_GB),
    min_mem_GB = min(mem_GB),
    max_mem_GB = max(mem_GB),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    res_time %>% dplyr::distinct(method, time_r2 = r_squared, time_slope = slope),
    by = "method"
  ) %>%
  dplyr::left_join(
    res_mem %>% dplyr::distinct(method, mem_r2 = r_squared, mem_slope = slope),
    by = "method"
  ) %>%
  dplyr::arrange(mean_time_secs)

utils::write.csv(
  summary_11,
  file.path(root, "data/derived/sf15_classic_profile_20261126/sf15_summary_11cohorts.csv"),
  row.names = FALSE
)

cat("=== 11-COHORT METHOD SUMMARY (CLASSIC-ALIGNED) ===\n")
print(as.data.frame(summary_11), digits = 4, row.names = FALSE)
cat("\n=== STEPHENSON & BASSEZ MEMORY RESIDUALS (log10 & ratio) ===\n")
print(
  as.data.frame(
    res_mem %>%
      dplyr::filter(dataset %in% c("Stephenson", "Bassez")) %>%
      dplyr::select(dataset, method, n_cells, mem_GB, log10_resid, ratio_obs_to_pred) %>%
      dplyr::arrange(dataset, method)
  ),
  digits = 4,
  row.names = FALSE
)
cat("\n=== TOP 10 ABS LOG10 MEM RESIDUALS ACROSS ALL 121 POINTS ===\n")
print(
  as.data.frame(
    res_mem %>%
      dplyr::arrange(dplyr::desc(abs(log10_resid))) %>%
      dplyr::slice_head(n = 10L) %>%
      dplyr::select(dataset, method, n_cells, mem_GB, log10_resid, ratio_obs_to_pred)
  ),
  digits = 4,
  row.names = FALSE
)

# --- Plotting (exact match to notebooks/benchmark_analysis.rmd) ---
regression_lines <- fit_log_log_lines(plot_data, "time_secs")
line_ends <- regression_lines %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    max_x = max(n_cells),
    pred_y = fitted[which.max(n_cells)],
    .groups = "drop"
  )
line_ends_log <- spread_log_label_y(
  line_ends %>%
    dplyr::mutate(
      log_max_x = log10(max_x),
      log_pred_y = log10(pred_y)
    ),
  min_gap = 0.22
)
line_ends_raw <- line_ends_log %>%
  dplyr::mutate(
    endpoint_x = 10^log_max_x,
    endpoint_y = 10^log_pred_y,
    label_x = 10^log_label_x,
    label_y = 10^log_label_y
  )
plot_x_breaks <- c(1e5, 3e5, 1e6)
plot_y_breaks <- c(1, 10, 100, 1000, 10000)

supp_fig_14a <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = n_cells, y = time_secs, color = method)
) +
  ggplot2::annotation_logticks(sides = "bl", color = "grey80") +
  ggplot2::geom_point(shape = 4, size = 1, stroke = 1.25, alpha = 0.5) +
  ggplot2::geom_line(
    data = regression_lines,
    ggplot2::aes(x = n_cells, y = fitted, color = method),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  ggplot2::geom_segment(
    data = line_ends_raw,
    ggplot2::aes(
      x = endpoint_x,
      y = endpoint_y,
      xend = label_x,
      yend = label_y
    ),
    inherit.aes = FALSE,
    color = "grey50",
    linewidth = 0.2
  ) +
  ggplot2::geom_text(
    data = line_ends_raw,
    ggplot2::aes(
      x = label_x,
      y = label_y,
      label = method,
      color = method
    ),
    hjust = 0,
    size = 4,
    show.legend = FALSE
  ) +
  ggplot2::labs(x = "Number of cells", y = "Time (s)") +
  ggplot2::theme_minimal() +
  ggplot2::scale_x_log10(
    breaks = plot_x_breaks,
    labels = function(x) sapply(x, scientific_A_10n),
    limits = c(
      min(plot_data$n_cells),
      max(plot_data$n_cells) * 10^0.65
    )
  ) +
  ggplot2::scale_y_log10(
    breaks = plot_y_breaks,
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme(legend.position = "none")

ram_regression_lines <- fit_log_log_lines(ram_data, "mem_GB")
ram_line_ends <- ram_regression_lines %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    max_x = max(n_cells),
    pred_y = fitted[which.max(n_cells)],
    .groups = "drop"
  )
ram_line_ends_log <- spread_log_label_y_centered(
  ram_line_ends %>%
    dplyr::mutate(
      log_max_x = log10(max_x),
      log_pred_y = log10(pred_y)
    ),
  min_gap = 0.15
)
ram_line_ends_raw <- ram_line_ends_log %>%
  dplyr::mutate(
    endpoint_x = 10^log_max_x,
    endpoint_y = 10^log_pred_y,
    label_x = 10^log_label_x,
    label_y = 10^log_label_y
  )
ram_x_breaks <- c(1e5, 3e5, 1e6)
ram_y_breaks <- c(1, 10, 100)

supp_fig_14b <- ggplot2::ggplot(
  ram_data,
  ggplot2::aes(x = n_cells, y = mem_GB, color = method)
) +
  ggplot2::annotation_logticks(sides = "bl", color = "grey80") +
  ggplot2::geom_point(shape = 4, size = 1, stroke = 1.25, alpha = 0.5) +
  ggplot2::geom_line(
    data = ram_regression_lines,
    ggplot2::aes(x = n_cells, y = fitted, color = method),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  ggplot2::geom_segment(
    data = ram_line_ends_raw,
    ggplot2::aes(
      x = endpoint_x,
      y = endpoint_y,
      xend = label_x,
      yend = label_y
    ),
    inherit.aes = FALSE,
    color = "grey50",
    linewidth = 0.2
  ) +
  ggplot2::geom_text(
    data = ram_line_ends_raw,
    ggplot2::aes(
      x = label_x,
      y = label_y,
      label = method,
      color = method
    ),
    hjust = 0,
    size = 4,
    show.legend = FALSE
  ) +
  ggplot2::labs(x = "Number of cells", y = "Peak RAM usage (GB)") +
  ggplot2::theme_minimal() +
  ggplot2::scale_x_log10(
    breaks = ram_x_breaks,
    labels = function(x) sapply(x, scientific_A_10n),
    limits = c(
      min(ram_data$n_cells),
      max(ram_data$n_cells) * 10^0.65
    )
  ) +
  ggplot2::scale_y_log10(
    breaks = ram_y_breaks,
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ggplot2::theme(legend.position = "none")

path_plots <- file.path(root, "plots")
dir.create(path_plots, showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  file.path(path_plots, "Supp_fig_14A_benchmark_exec_times.pdf"),
  plot = supp_fig_14a,
  width = 6.0,
  height = 3.5,
  limitsize = FALSE
)
ggplot2::ggsave(
  file.path(path_plots, "Supp_fig_14B_benchmark_mem_GB.pdf"),
  plot = supp_fig_14b,
  width = 6.0,
  height = 3.5,
  limitsize = FALSE
)

supp_fig_14 <- patchwork::wrap_plots(list(supp_fig_14a, supp_fig_14b), nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  file.path(path_plots, "Supp_fig_14.pdf"),
  plot = supp_fig_14,
  width = 12,
  height = 3.5,
  limitsize = FALSE
)
ggplot2::ggsave(
  file.path(path_plots, "Supp_fig_14.png"),
  plot = supp_fig_14,
  width = 12,
  height = 3.5,
  dpi = 300,
  limitsize = FALSE
)
cat("DONE: Regenerated Supp_fig_14A, Supp_fig_14B, Supp_fig_14.pdf, and Supp_fig_14.png.\n")
