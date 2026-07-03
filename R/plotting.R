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
      select(Feature, Loading = !!sym(pc))
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
