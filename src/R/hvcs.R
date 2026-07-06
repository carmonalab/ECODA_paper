# ============================================================
# HVC (HIGHLY VARIABLE CELL TYPE) FUNCTIONS
# ============================================================

# get_ct_var and helpers ----

get_ct_var <- function(
  df,
  show_plot = TRUE,
  plot_title = "",
  smooth_method = "lm",
  descending = TRUE
) {
  df_var <- df %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "celltype",
      values_to = "values"
    ) %>%
    dplyr::group_by(.data$celltype) %>%
    dplyr::summarize(
      Relative_abundance = mean(values, na.rm = TRUE),
      Variance = var(values, na.rm = TRUE)
    )

  if (descending) {
    df_var <- df_var %>%
      dplyr::arrange(dplyr::desc(.data$Variance))
  } else {
    df_var <- df_var %>%
      dplyr::arrange(.data$Variance)
  }

  if (show_plot) {
    p <- varmeanplot(
      data = df_var,
      title = plot_title,
      smooth_method = smooth_method
    )
    print(p)
  }

  return(df_var)
}


get_hvcs <- function(df_var, top_n_hvcs = NULL, variance_threshold = 0.8) {
  if (!is.null(top_n_hvcs)) {
    if (top_n_hvcs < 1) {
      top_hvcs <- df_var %>%
        slice_head(prop = top_n_hvcs) %>%
        pull(celltype)
    } else {
      top_hvcs <- df_var %>%
        slice_head(n = top_n_hvcs) %>%
        pull(celltype)
    }
  } else {
    top_hvcs <- select_by_variance_explained(
      df_var,
      variance_threshold = variance_threshold
    )
  }

  # Select at least two cell types
  if (length(top_hvcs) < 2) {
    top_hvcs <- df_var %>%
      slice_head(n = 2) %>%
      pull(celltype)
  }

  return(top_hvcs)
}


select_by_variance_explained <- function(df_var, variance_threshold) {
  # Ensure data is sorted by variance in descending order
  df_var_sorted <- df_var %>%
    dplyr::arrange(dplyr::desc(.data$Variance))

  # Calculate the total variance
  total_variance <- sum(df_var_sorted$Variance)

  # Calculate the cumulative variance and the percentage of variance explained
  df_with_cumulative_var <- df_var_sorted %>%
    dplyr::mutate(
      cumulative_variance = cumsum(Variance),
      variance_explained = (cumulative_variance / total_variance)
    )

  # Select the cell types that meet the variance threshold
  selected_celltypes <- df_with_cumulative_var %>%
    dplyr::filter(.data$variance_explained <= variance_threshold) %>%
    dplyr::pull(celltype)

  return(selected_celltypes)
}


varmeanplot <- function(
  data,
  title = "",
  smooth_method = "lm",
  label_points = FALSE
) {
  p <- ggplot(data, aes(x = Relative_abundance, y = Variance)) +
    geom_point() +
    geom_smooth(
      method = smooth_method,
      color = "red",
      fill = "#69b3a2",
      se = TRUE
    ) +
    labs(title = paste(title)) +
    theme_classic() +
    xlab("Mean") +
    ylab("Variance")

  if (label_points) {
    p <- p +
      ggrepel::geom_text_repel(data = data, aes(label = celltype), vjust = -0.5)
  }

  return(p)
}


plot_varmean <- function(
  df_var,
  highlight_celltypes = NULL,
  plot_title = "",
  highlight_hvcs = TRUE,
  labels = c("only_hvc", "all", "none"),
  plot_fit_line = FALSE,
  smooth_method = "lm"
) {
  labels <- match.arg(labels)

  # --- 1. Create a highlighting factor column ---
  if (highlight_hvcs) {
    df_var <- df_var %>%
      mutate(
        is_highlighted = if_else(
          .data$celltype %in% highlight_celltypes,
          "HVC selected",
          "Not selected"
        )
      )
    color_map <- c("HVC selected" = "red", "Not selected" = "black")

    p <- ggplot(
      df_var,
      aes(
        x = Relative_abundance,
        y = Variance,
        color = is_highlighted
      )
    ) +
      scale_color_manual(values = color_map, name = "Cell Type Group")
  } else {
    p <- ggplot(df_var, aes(x = avg_clr_abundance, y = Variance))
  }

  p <- p +
    geom_point() +
    labs(title = paste(plot_title)) +
    theme_classic() +
    xlab("Mean (CLR)") +
    ylab("Variance (CLR)")

  if (plot_fit_line) {
    p <- p +
      geom_smooth(
        method = smooth_method,
        color = "red",
        fill = "#69b3a2",
        se = TRUE
      )
  }

  # Ensure the legend is not shown if we are not highlighting anything
  if (!highlight_hvcs) {
    p <- p + theme(legend.position = "none")
  }

  if (labels != "none") {
    # Filter data for labeling based on user choice
    label_df <- df_var
    if (labels == "only_hvc") {
      label_df <- df_var %>%
        dplyr::filter(.data$celltype %in% highlight_celltypes)
    }

    # Use the color aesthetic inside ggrepel so
    # highlighted labels match the point color
    if (highlight_hvcs) {
      p <- p +
        geom_text_repel(
          data = label_df,
          aes(label = celltype, color = is_highlighted),
          size = 3,
          show.legend = FALSE
        )
    } else {
      p <- p +
        geom_text_repel(
          data = label_df,
          aes(label = celltype),
          size = 3,
          color = "black"
        )
    }
  }

  return(p)
}