# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

plot_pca <- function(
  feat_mat,
  labels,
  scale. = FALSE,
  pca_dims = NULL,
  knn_k = 3,
  title = NULL,
  cluster_score = TRUE,
  mod_score = TRUE,
  sil_score = FALSE,
  anosim_score = TRUE,
  digits = 3,
  pointsize = 3,
  labelsize = 4,
  coord_equal = TRUE,
  axes = c(1, 2),
  plotly_3d = FALSE,
  invisible = c("var", "quali"),
  n_ct_show = Inf,
  repel = FALSE
) {
  res.pca <- prcomp(feat_mat, scale. = scale., rank. = pca_dims)
  dist_mat <- dist(feat_mat)
  format_str <- paste0("%.", digits, "f")

  if (anosim_score) {
    anosim_score <- round(
      vegan::anosim(x = dist_mat, grouping = labels, distance = "euclidean")[[
        "statistic"
      ]],
      3
    )
    title <- paste0(
      title,
      "\nANOSIM score: ",
      sprintf(format_str, anosim_score)
    )
  }
  if (cluster_score) {
    cluster_score <- clust_eval(dist_mat, labels)
    title <- paste0(title, "\nARI: ", sprintf(format_str, cluster_score))
  }
  if (mod_score) {
    mod_score <- round(calc_modularity(dist_mat, labels, knn_k), 3)
    title <- paste0(
      title,
      "\nModularity score: ",
      sprintf(format_str, mod_score)
    )
  }
  if (sil_score) {
    sil_score <- round(calc_sil(dist_mat, labels), 3)
    title <- paste0(
      title,
      "\nSilhouette score: ",
      sprintf(format_str, sil_score)
    )
  }

  if (plotly_3d) {
    df <- as.data.frame(res.pca$x)
    df$id <- seq_len(nrow(df))
    df$vs <- factor(labels)
    ms <- replicate(2, df, simplify = F)
    ms[[2]]$PC3 <- min(df$PC3)
    m <- ms %>% bind_rows() %>% plotly::group2NA("id", "vs")
    p <- plotly::plot_ly(color = ~vs) %>%
      plotly::add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
      plotly::add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 0.2)
  } else {
    eig.val <- factoextra::get_eig(res.pca)
    pc1_var <- round(eig.val[1, 2], 1)
    pc2_var <- round(eig.val[2, 2], 1)
    p <- factoextra::fviz_pca(
      res.pca,
      axes = axes,
      habillage = labels,
      label = "var",
      pointsize = pointsize,
      labelsize = labelsize,
      invisible = invisible,
      select.var = list(contrib = n_ct_show),
      repel = repel,
      geom = "point",
      axes.linetype = NA
    ) +
      ggtitle(title) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      labs(
        x = paste0("PCA dim1 (", pc1_var, "%)"),
        y = paste0("PCA dim2 (", pc2_var, "%)")
      ) +
      scale_shape_manual(values = rep(19, length(unique(labels))))
    if (coord_equal) p <- p + coord_equal()
  }
  return(p)
}

plot_pca_contributions_horizontal <- function(
  res.pca,
  pcs = c("PC1", "PC2"),
  n = 3,
  absolute = TRUE
) {
  loadings <- as.data.frame(res.pca$rotation)
  missing_pcs <- setdiff(pcs, colnames(loadings))
  if (length(missing_pcs) > 0) {
    stop(paste("Missing PCs:", paste(missing_pcs, collapse = ", ")))
  }

  plot_data <- lapply(pcs, function(pc) {
    df <- loadings %>%
      rownames_to_column(var = "Feature") %>%
      dplyr::select(Feature, Loading = !!sym(pc))
    if (absolute) {
      top_bottom <- df %>%
        mutate(Loading = abs(.data$Loading)) %>%
        arrange(desc(.data$Loading)) %>%
        head(n) %>%
        mutate(PC = pc)
    } else {
      df <- df %>% arrange(desc(.data$Loading))
      top_bottom <- bind_rows(head(df, n), tail(df, n)) %>%
        distinct() %>%
        mutate(PC = pc)
    }
    return(top_bottom)
  }) %>%
    bind_rows()

  plot_data <- plot_data %>%
    mutate(Feature_PC = factor(paste(.data$Feature, .data$PC, sep = "___"))) %>%
    arrange(.data$PC, .data$Loading) %>%
    mutate(Feature_PC = factor(.data$Feature_PC, levels = .data$Feature_PC))

  x_label <- if (absolute) "Absolute PCA Loading" else "PCA Loading"
  ggplot(plot_data, aes(y = Feature_PC, x = Loading)) +
    geom_col(fill = "gray70", color = "black", linewidth = 0.4) +
    scale_y_discrete(labels = function(x) gsub("___.*", "", x)) +
    facet_wrap(~PC, nrow = length(pcs), scales = "free_y") +
    theme_classic() +
    labs(title = "Main Contributing Features", y = "Feature", x = x_label) +
    theme(
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      strip.background = element_blank(),
      panel.spacing = unit(1.5, "lines"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold", color = "black"),
      plot.title = element_text(face = "bold")
    )
}

plot_mds <- function(
  dist_mat,
  labels,
  knn_k = 3,
  title = NULL,
  cluster_score = TRUE,
  mod_score = TRUE,
  sil_score = FALSE,
  anosim_score = TRUE,
  digits = 3,
  pointsize = 3,
  labelsize = 4,
  coord_equal = TRUE,
  axes = c(1, 2)
) {
  mds_res <- cmdscale(dist_mat, k = max(axes), eig = TRUE)
  format_str <- paste0("%.", digits, "f")
  mds_df <- data.frame(
    Dim1 = mds_res$points[, axes[1]],
    Dim2 = mds_res$points[, axes[2]],
    Cluster = labels
  )

  if (anosim_score) {
    anosim_score_val <- round(
      vegan::anosim(x = dist_mat, grouping = labels, distance = "euclidean")[[
        "statistic"
      ]],
      3
    )
    title <- paste0(
      title,
      "\nANOSIM score: ",
      sprintf(format_str, anosim_score_val)
    )
  }
  if (cluster_score) {
    cluster_score_val <- clust_eval(dist_mat, labels)
    title <- paste0(title, "\nARI: ", sprintf(format_str, cluster_score_val))
  }
  if (mod_score) {
    mod_score_val <- round(calc_modularity(dist_mat, labels, knn_k), 3)
    title <- paste0(
      title,
      "\nModularity score: ",
      sprintf(format_str, mod_score_val)
    )
  }
  if (sil_score) {
    sil_score_val <- round(calc_sil(dist_mat, labels), 3)
    title <- paste0(
      title,
      "\nSilhouette score: ",
      sprintf(format_str, sil_score_val)
    )
  }

  eig <- mds_res$eig
  perc_var <- round(100 * eig / sum(eig[eig > 0]), 1)
  lab_x <- paste0("MDS dim", axes[1], " (", perc_var[axes[1]], "%)")
  lab_y <- paste0("MDS dim", axes[2], " (", perc_var[axes[2]], "%)")

  p <- ggplot2::ggplot(
    mds_df,
    aes(x = Dim1, y = Dim2, color = Cluster, shape = Cluster)
  ) +
    geom_point(size = pointsize) +
    scale_shape_manual(values = rep(19, length(unique(labels)))) +
    labs(
      x = lab_x,
      y = lab_y,
      title = title,
      color = "Groups",
      shape = "Groups"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  if (coord_equal) {
    p <- p + coord_equal()
  }
  return(p)
}

# Shared cleaning step for the funky-heatmap figures (Figure 2A, Supp fig 15,
# "For presentation"): merge exec times, filter methods, recode method/score/
# dataset display names. `score_set` optionally restricts to the configured
# metric names (benchmark_metrics); NULL keeps every score column (needed by
# Supp fig 1, whose correlation matrix includes Silhouette and LISI).
prepare_benchmark_data <- function(
  df_results,
  exec_times,
  filter_methods,
  method_recode,
  dataset_recode,
  score_set = NULL
) {
  clean_data <- merge_exec_times(df_results, exec_times) %>%
    filter(method %in% filter_methods) %>%
    mutate(method = recode(method, !!!method_recode)) %>%
    mutate(score = recode(score, !!!score_label_map)) %>%
    mutate(dataset = recode(dataset, !!!dataset_recode))
  if (!is.null(score_set)) {
    # dplyr::recode keeps unmatched values, so without this filter the raw
    # mod_knnsqrtn_score / mod_knn6_score / mod_knn9_score leftovers would
    # leak into the metric columns of the heatmap.
    clean_data <- clean_data %>%
      filter(score %in% unname(score_label_map[score_set]))
  }
  return(clean_data)
}

# Parameterized funky-heatmap builder (refactor of the previously duplicated
# Figure 2A / Supp fig 15 / "For presentation" blocks). Mirrors the original
# pipeline: min-max per (dataset, score) -> rank/overall/dataset/metric/runtime
# aggregates -> column_info -> funky_heatmap -> ggsave. Outputs are identical
# to the legacy blocks except that only the configured `score_set` metrics
# are ranked/displayed. Returns the funkyheatmap grob invisibly.
build_funky_heatmap <- function(
  df_results,
  exec_times,
  filter_methods,
  method_recode,
  dataset_recode,
  score_set = c("anosim_score", "mod_knn3_score", "cluster_score"),
  output_file,
  output_dir,
  plot_width = 10,
  plot_height = 6,
  method_width = 6,
  overall_col_name = "Mean",
  include_rank = FALSE,
  overall_legend_title = "Overall Mean",
  overall_legend_labels = c(". 0", rep("", 5), "0.5", rep("", 4), "1")
) {
  ncolors <- 11

  clean_data <- prepare_benchmark_data(
    df_results, exec_times, filter_methods, method_recode, dataset_recode,
    score_set = score_set
  )

  # --- 2. CALCULATE GEOMETRIC MEAN OF RANKS ---
  df_ranks <- clean_data %>%
    group_by(dataset, score) %>%
    mutate(rank = rank(-value, ties.method = "min")) %>%
    ungroup() %>%
    group_by(method) %>%
    dplyr::summarise(geo_mean_rank = exp(mean(log(rank), na.rm = TRUE))) %>%
    ungroup() %>%
    # Normalize so it plots nicely (0 to 1). Invert: best rank (lowest
    # number) = 1.0 (biggest bar).
    mutate(Rank_Score = (max(geo_mean_rank) - geo_mean_rank) /
             (max(geo_mean_rank) - min(geo_mean_rank))) %>%
    dplyr::select(method, Rank_Score)

  # --- 3. SCALE DATA (0-1) FOR HEATMAP VISUALS ---
  scaled_data <- clean_data %>%
    group_by(dataset, score) %>%
    mutate(value = min_max(value)) %>%
    ungroup()

  # --- 4. AGGREGATE SCORES ---
  df_overall <- scaled_data %>%
    group_by(method) %>%
    dplyr::summarise(Overall_Score = mean(value, na.rm = TRUE))

  df_datasets_wide <- scaled_data %>%
    group_by(method, dataset) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = method, names_from = dataset, values_from = value,
      names_prefix = "Dataset__"
    )

  df_metrics_wide <- scaled_data %>%
    group_by(method, score) %>%
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      id_cols = method, names_from = score, values_from = value,
      names_prefix = "Metric__"
    )

  df_runtime <- clean_data %>%
    distinct(method, dataset, time_secs) %>%
    group_by(method) %>%
    dplyr::summarise(raw_time = mean(time_secs, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      log_time = log10(raw_time + 1),
      # Inverting so SHORTER time (faster) = BIGGER bar/BETTER score
      Runtime = min_max(log_time)
    ) %>%
    dplyr::select(method, Runtime)

  # --- 5. MERGE EVERYTHING ---
  final_data <- df_overall %>%
    left_join(df_ranks, by = "method") %>%
    left_join(df_datasets_wide, by = "method") %>%
    left_join(df_metrics_wide, by = "method") %>%
    left_join(df_runtime, by = "method") %>%
    arrange(desc(Overall_Score))

  # --- 6. BUILD COLUMN INFO ---
  info_method <- tribble(
    ~id,      ~group, ~name, ~geom,  ~palette, ~options,
    "method", "",     "",    "text", NA,       list(hjust = 0, width = method_width)
  )

  info_score_avg <- tribble(
    ~id,             ~group,    ~name,            ~geom, ~palette,           ~options,
    "Overall_Score", "Overall", overall_col_name, "bar", "palette_avg_mean", list(width = 4)
  )

  palettes <- list(palette_avg_mean = brewer.pal(n = ncolors, name = "RdBu"))
  legends <- list()
  column_info_parts <- list(info_method, info_score_avg)

  if (include_rank) {
    info_rank_avg <- tribble(
      ~id,          ~group,    ~name,      ~geom, ~palette,          ~options,
      "Rank_Score", "Overall", "Rank Mean", "bar", "palette_avg_rank", list(width = 4, legend = FALSE)
    )
    column_info_parts <- c(column_info_parts, list(info_rank_avg))
    palettes[["palette_avg_rank"]] <- brewer.pal(n = ncolors, name = "RdBu")
    legends <- c(legends, list(
      list(palette = "palette_avg_rank", enabled = FALSE, title = "REMOVE PLACEHOLDER")
    ))
  }

  dataset_cols <- colnames(final_data)[startsWith(colnames(final_data), "Dataset__")]
  info_datasets <- tibble(id = dataset_cols) %>%
    mutate(
      group = "Datasets",
      name = str_remove(id, "Dataset__"),
      geom = "circle",
      palette = "Datasets",
      options = list(list())
    ) %>%
    arrange(name)

  metric_cols <- colnames(final_data)[startsWith(colnames(final_data), "Metric__")]
  target_metric_order <- unname(score_label_map[score_set])
  info_metrics_main <- tibble(id = metric_cols) %>%
    mutate(
      group = "Metrics",
      name = str_remove(id, "Metric__"),
      geom = "bar",
      palette = "Metrics",
      options = list(list(width = 2)),
      name_f = factor(name, levels = target_metric_order)
    ) %>%
    arrange(name_f) %>%
    dplyr::select(-name_f)

  info_runtime <- tribble(
    ~id,       ~group,    ~name,     ~geom, ~palette,          ~options,
    "Runtime", "Runtime", "Runtime", "bar", "palette_runtime", list(width = 3)
  )

  column_info <- bind_rows(c(
    column_info_parts,
    list(info_datasets),
    list(info_metrics_main),
    list(info_runtime)
  ))

  palettes[["Datasets"]] <- brewer.pal(n = ncolors, name = "RdBu")
  palettes[["Metrics"]] <- brewer.pal(n = ncolors, name = "RdBu")
  palettes[["palette_runtime"]] <- rev(brewer.pal(n = ncolors, name = "RdYlGn"))

  legends <- c(legends, list(
    list(
      palette = "palette_avg_mean",
      enabled = TRUE,
      title = overall_legend_title,
      labels = overall_legend_labels
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
      labels = c(". 0", rep("", round((ncolors - 2) / 2) + 1), "0.5",
                 rep("", round((ncolors - 2) / 2)), "1")
    ),
    list(
      palette = "palette_runtime",
      enabled = TRUE,
      title = "Runtime log10(s)",
      labels = c("secs", "mins", "hours")
    )
  ))

  # --- 7. PLOT ---
  p <- funky_heatmap(
    final_data,
    column_info = column_info,
    scale_column = FALSE,
    column_groups = tibble(
      group = c("Overall", "Datasets", "Metrics", "Runtime")
    ),
    palettes = palettes,
    legends = legends,
    position_args = position_arguments(
      col_width = 1.5,
      col_space = 0.2,
      col_annot_offset = 5
    )
  )
  ggsave(file.path(output_dir, output_file),
         width = plot_width,
         height = plot_height)
  return(invisible(p))
}
