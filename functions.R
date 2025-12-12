suppressPackageStartupMessages({
  # install.packages("devtools")
  library(compositions)
  library(doParallel)
  library(factoextra)
  library(foreach)
  library(ggplot2)
  library(ggpubr)
  library(ggtext)
  library(mclust)
  library(plotly)
  # library("robCompositions")
  library(Seurat)
  library(limma)
  # devtools::install_github('satijalab/seurat-data')
  library(SeuratData)
  # install.packages("hdf5r")
  # remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
  library(SignatuR)
  library(stringr)
  library(tidyr)
  library(tidyverse)
  # install.packages("zCompositions")
  # library(zCompositions)
  library(rstatix)

  library(parallelly)

  # remotes::install_github("carmonalab/scooter", ref="f31eab3")
  library(scooter)

  # Import at the end, otherwise the "select" function gets re-assigned by another package
  library(dplyr)

  # library(arrow)
})


clr <- function(df) {
  geometric_mean <- apply(df, 1, function(row) exp(mean(log(row))))
  clr_df <- apply(df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()

  return(clr_df)
}


# Calculate coefficient of variation
cv <- function(x) {
  sd_x <- sd(x) # Standard deviation
  mean_x <- mean(x) # Mean
  cv_value <- abs(sd_x / mean_x) * 100 # Coefficient of variation (%)
  return(cv_value)
}


calc_perc_df <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>% as.data.frame()
  return(df)
}


calc_sep_score <- function(feat_mat,
                           labels,
                           knn_k = NULL) {
  # sil_score <- round(calc_sil(df, labels), 3)
  mod_score <- unlist(round(calc_modularity(feat_mat, labels, knn_k), 3))
  mod_knn3_score <- unlist(round(calc_modularity(feat_mat, labels, 3), 3))
  cluster_score <- clust_eval(matrix = feat_mat, labels)
  anosim_score <- vegan::anosim(x = feat_mat, grouping = labels, distance = "euclidean")[["statistic"]]

  res <- list(
    # sil_score = sil_score,
    mod_score = mod_score,
    mod_knn3_score = mod_knn3_score,
    anosim_score = anosim_score,
    cluster_score = cluster_score
  )

  return(res)
}


# Calculate average silhouette width
calc_sil <- function(feat_mat,
                     labels) {
  sils <- cluster::silhouette(
    x = as.numeric(factor(labels)),
    dist = dist(feat_mat)
  ) %>%
    as.data.frame()
  score <- mean(sils[["sil_width"]])
  return(score)
}


#----------------------------------------------------------->
calc_modularity <- function(feat_mat,
                            labels,
                            knn_k = NULL) {
  ngroups <- length(unique(labels))

  if (is.null(knn_k)) {
    # min_group_size <- min(table(labels))
    # half_group_size <- round(nrow(feat_mat) / ngroups / 2)
    # # knn_k is equal to half average group size or at least 3
    # knn_k <- max(half_group_size, 3)

    knn_k <- max(3, round(sqrt(nrow(feat_mat))))
  }

  # Create a graph object
  g <- compute_snn_graph(feat_mat = feat_mat, knn_k = knn_k)
  # Compute modularity
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))

  # NOTE:
  # Maximum modularity depends on the number of groups: max(mod) = 1 - 1 / (number of groups)
  # see Brandes, Ulrik, et al. On finding graph clusterings with maximum modularity.

  # Adjust modularity score for number of groups
  maximum_modularity_score <- 1 - (1 / ngroups)
  adjusted_modularity_score <- modularity_score / maximum_modularity_score

  return(adjusted_modularity_score)
}


calc_avg_sep_score <- function(feat_mat, labels, digits = 2) {
  sep_scores <- calc_sep_score(feat_mat = feat_mat, labels = labels)
  avg_sep_score <- sep_scores[c("mod_knn3_score", "anosim_score", "cluster_score")] %>%
    unlist() %>%
    mean() %>%
    round(digits)
}


compute_KNN <- function(feat_mat, knn_k) {
  # Compute KNN
  knn <- RANN::nn2(as.matrix(feat_mat), k = knn_k + 1)$nn.idx
  knn <- knn[, -1] # Remove self-neighbor
  return(knn)
}


compute_snn_graph <- function(feat_mat, knn_k) {
  knn <- compute_KNN(feat_mat = feat_mat, knn_k = knn_k)
  # Initialize adjacency matrix
  n <- nrow(as.matrix(feat_mat))
  adj_matrix <- matrix(0, n, n)
  # Count shared neighbors
  for (i in seq_len(n)) {
    for (j in knn[i, ]) {
      shared_neighbors <- length(intersect(knn[i, ], knn[j, ]))
      adj_matrix[i, j] <- shared_neighbors
      adj_matrix[j, i] <- shared_neighbors # Ensure symmetry
    }
  }
  # adj_matrix <- spdep::knearneigh(feat_mat, k=5, longlat = NULL)
  # Create graph object
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  return(g)
}


#-----------------------------------------------------------<


# Cluster samples and compare to original annotation
clust_eval <- function(matrix,
                       labels,
                       nclusts = NULL,
                       digits = 3,
                       return_mean = TRUE) {
  results <- list()
  dist_mat <- dist(matrix)

  if (is.null(nclusts)) {
    nclusts <- length(unique(labels))
  }

  # Perform hierarchical clustering
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  clust_labels <- stats::cutree(hc, k = nclusts)
  results[["hclust_accuracy"]] <- mclust::adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  # Perform PAM clustering
  clust_labels <- cluster::pam(matrix, k = nclusts)$cluster
  results[["pamclust_accuracy"]] <- mclust::adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  if (return_mean) {
    return(round(mean(unlist(results)), digits))
  } else {
    results[["hclust_accuracy"]] <- round(results[["hclust_accuracy"]], digits)
    results[["pamclust_accuracy"]] <- round(results[["pamclust_accuracy"]], digits)
    return(results)
  }
}


convert_seurat_to_h5ad <- function(seurat,
                                   ds,
                                   shared_storage_path = "/Users/christianhalter/Library/CloudStorage/OneDrive-Personal/Temp") {
  temp_file_path <- file.path(getwd(), "data")
  temp_file_path_h5ad <- file.path(temp_file_path, paste0(ds, ".h5ad"))
  file_path_h5seurat <- file.path(temp_file_path, paste0(ds, ".h5Seurat"))
  file_path_h5ad <- file.path(shared_storage_path, paste0(ds, ".h5ad"))
  if (!file.exists(file_path_h5ad)) {
    SaveH5Seurat(seurat, filename = file_path_h5seurat)
    Convert(file_path_h5seurat, dest = "h5ad")
    unlink(file_path_h5seurat)
    file.rename(temp_file_path_h5ad, file_path_h5ad)
    print(paste0("File saved to: ", file_path_h5ad))
  } else {
    print("File already processed")
  }
}


create_clean_seuratv5_object <- function(seurat) {
  # Create new v5 seurat object
  # Remove any possible processing, only saving counts and metadata
  seurat <- CreateSeuratObject(counts = seurat@assays$RNA$counts, meta.data = seurat@meta.data)

  return(seurat)
}


# Test data transformation methods for ECODA
datrans <- function(count_mat,
                    labels = NULL,
                    Amount_of_perturbation, # Percent cell abundance difference (e.g. 100 equals one cell type being twice as abundant)
                    n_ct_to_select, # Number of randomly selected cell types to be differentially abundant
                    cts = NULL, # Cell types to be differentially abundant. If NULL, randomly select a specified number of cell types (n_ct_to_select)
                    reps = 20, # Number of random shuffling to calculate separation using different cell types and samples for DA
                    trans_method = c(
                      "counts",
                      # "counts_imputed",
                      # "counts_pca",
                      "freq",
                      # "freq_imputed",
                      # "freq_pca",
                      "arcsine_sqrt",
                      # "arcsine_sqrt_pca",
                      "alr_randref",
                      # "alr_randref_pca",
                      "alr_mincvref",
                      # "alr_mincvref_pca",
                      # "ilr", "ilr_pca",
                      "clr" # ,
                      # "clr_pca"
                    ),
                    zero_imp_method = "counts_all__1",
                    n_cores = 8) {
  colnames(count_mat) <- make.names(colnames(count_mat), unique = TRUE)

  n_half_samples <- round(dim(count_mat)[1] / 2)

  if (!is.null(cts)) {
    n_ct_to_select[n_ct_to_select >= length(cts)] <- length(cts)
  }

  # Register cluster
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)

  rets <- foreach(
    da = Amount_of_perturbation,
    .export = c(
      "calc_perc_df",
      "impute_zeros",
      "calc_sil",
      "calc_modularity", "compute_snn_graph", "compute_KNN",
      "clust_eval", "adjustedRandIndex",
      "cv",
      "clr"
    ),
    .packages = c("dplyr", "robCompositions", "zCompositions"),
    .errorhandling = "pass",
    .combine = rbind
  ) %dopar% {
    res <- data.frame(
      trans_method = character(),
      zero_imp_method = character(),
      n_celltypes = numeric(),
      Amount_of_perturbation = numeric(),
      Silhouette_score = numeric(),
      ANOSIM_score = numeric(),
      Modularity_score = numeric(),
      Adjusted_Rand_Index = numeric(),
      bootstrap_id = numeric(),
      diff_abu_cts = list()
    )

    print(paste0("da: ", da))
    for (nct in n_ct_to_select) {
      print(paste0("nct: ", nct))
      for (rep in 1:reps) {
        # Prepare data
        df_counts_temp <- count_mat

        if (!is.null(cts)) {
          ct_da <- sample(cts, size = nct)
        } else {
          ct_da <- sample(colnames(df_counts_temp), size = nct)
        }

        if (is.null(labels)) {
          half_samples_da <- sample(row.names(df_counts_temp), size = n_half_samples)
          labels_random <- as.numeric(row.names(count_mat) %in% half_samples_da)

          # Simulate differential abundance
          rsums_before <- rowSums(df_counts_temp)
          df_counts_temp[half_samples_da, ct_da] <- round(df_counts_temp[half_samples_da, ct_da] * da)
          rsums_after <- rowSums(df_counts_temp)
          df_counts_temp <- round(df_counts_temp / (rsums_after / rsums_before))
        }

        # Remove columns if they contain all zeros
        df_counts_temp <- df_counts_temp %>%
          select_if(colSums(.) != 0) %>%
          mutate_all(as.numeric)

        for (zmet in zero_imp_method) {
          df_freq <- df_counts_temp %>% calc_perc_df()
          df_arcsine_sqrt <- asin(sqrt(df_freq / 100))

          if (grepl("percentage|counts|multRepl", zmet)) {
            zero_imp_method_split <- strsplit(zmet, "__")
            if (grepl("percentage|counts", zmet)) {
              df_counts_temp_imputed <- df_counts_temp %>% impute_zeros(clr_zero_impute_method = zero_imp_method_split[[1]][1], clr_zero_impute_num = as.numeric(zero_imp_method_split[[1]][2]))
              df_freq_imputed <- df_counts_temp_imputed %>% calc_perc_df()
            } else if (grepl("multRepl", zmet)) {
              df_freq_imputed <- df_freq %>% zCompositions::multRepl(label = 0, dl = rep(as.numeric(zero_imp_method_split[[1]][2]), ncol(df_freq)), z.warning = 1, frac = 1)
            }
          } else if (zmet == "multLN") {
            df_freq_imputed <- df_freq %>% zCompositions::multLN(label = 0, dl = rep(0.1, ncol(df_freq)), z.warning = 0.9)
          }

          for (met in trans_method) {
            if (grepl("counts", met)) {
              df <- df_counts_temp
              # } else if (grepl("counts_imputed", met)) {df <- df_counts_temp_imputed
            } else if (grepl("freq", met)) {
              df <- df_freq
              # } else if (grepl("freq_imputed", met)) {df <- df_freq_imputed
            } else if (grepl("arcsine_sqrt", met)) {
              df <- df_arcsine_sqrt
            } else if (grepl("alr_mincvref", met)) {
              ct_ref <- sample(colnames(df_freq_imputed)[!colnames(df_freq_imputed) %in% ct_da], size = 1)
              df <- Hotelling::alr(as.formula(paste0(ct_ref, "~.")), df_freq_imputed)
            } else if (grepl("alr_randref", met)) {
              cvs <- apply(df_freq_imputed, 2, cv)
              min_cv <- min(cvs[!colnames(df_freq_imputed) %in% ct_da])
              ct_ref_mincv <- colnames(df_freq_imputed)[which(cvs == min_cv)][1]
              df <- Hotelling::alr(as.formula(paste0(ct_ref_mincv, "~.")), df_freq_imputed)
            } else if (grepl("ilr", met)) {
              df <- compositions::ilr(df_freq)
            } else if (grepl("clr", met)) {
              df <- clr(df_freq_imputed)
              # } else if (met == "clr_centered") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = FALSE)
              # } else if (met == "clr_centered_scaled") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = TRUE)
            }

            if (grepl("pca", met)) {
              df <- prcomp(df)$x[, 1:2]
            }


            # Calculate scores

            if (is.null(labels)) {
              avg_sil <- calc_sil(feat_mat = df, labels_random)
              # Reduced number of permutations to speed up calculation time
              anosim_score <- vegan::anosim(x = df, grouping = labels_random, distance = "euclidean", permutations = 49)[["statistic"]]
              mod <- calc_modularity(feat_mat = df, labels_random)
              cluster_score <- clust_eval(matrix = df, labels_random)
            } else {
              avg_sil <- calc_sil(feat_mat = df, labels)
              # Reduced number of permutations to speed up calculation time
              anosim_score <- vegan::anosim(x = df, grouping = labels, distance = "euclidean", permutations = 49)[["statistic"]]
              mod <- calc_modularity(feat_mat = df, labels)
              cluster_score <- clust_eval(matrix = df, labels)
            }


            # Append results
            new_row_df <- data.frame(
              trans_method = met,
              zero_imp_method = zmet,
              n_celltypes = nct,
              Amount_of_perturbation = da,
              Silhouette_score = avg_sil,
              ANOSIM_score = anosim_score,
              Modularity_score = mod,
              Adjusted_Rand_Index = cluster_score,
              bootstrap_id = rep
            )

            new_row_df$diff_abu_cts <- list(ct_da)
            new_row_df$feat_mat <- list(df)

            res <- rbind(res, new_row_df)
          }
        }
      }
    }

    return(res)
  }
  # stop cluster
  stopCluster(cluster)

  return(rets)
}


DESeq2.normalize <- function(matrix,
                             metadata,
                             n_hvg = 2000) {
  suppressMessages({
    suppressWarnings({
      # Normalize pseudobulk data using DESeq2
      # do formula for design with the cluster_by elements in order
      matrix <- DESeq2::DESeqDataSetFromMatrix(
        countData = matrix,
        colData = metadata,
        design = stats::formula(paste("~ 1"))
      )

      matrix <- DESeq2::estimateSizeFactors(matrix)

      # Set minimum number of counts per gene
      nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(matrix, normalized = TRUE)) > 10))

      # transform counts using vst
      matrix <- DESeq2::vst(matrix, blind = T, nsub = nsub)
      matrix <- SummarizedExperiment::assay(matrix)

      # get top variable genes
      rv <- MatrixGenerics::rowVars(matrix)
      select <- order(rv, decreasing = TRUE)[seq_len(min(n_hvg, length(rv)))]
      select <- row.names(matrix)[select]

      matrix <- matrix[select[select %in% row.names(matrix)], ]
    })
  })

  return(matrix)
}




exec_time <- function(fun) {
  start_time <- Sys.time()
  fun
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  time_taken
  return(time_taken)
}




FindClusters_multi <- function(seurat,
                               res_broad = c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100),
                               res_broad_2sub = c(), # c(0.1, 0.5),
                               res_sub = c() # c(0.1)
) {
  for (rb in res_broad) {
    print(paste0("Running broad clustering at resolution: ", rb))
    Idents(seurat) <- "Sample" # Resetting Idents for initial broad clustering

    # Run broad clustering
    seurat <- FindClusters(seurat, resolution = rb)

    broad_clus_col_name <- paste0("RNA_snn_res.", rb)

    # Sub-clustering based on res_broad_2sub
    if (rb %in% res_broad_2sub) {
      print(paste0("Proceeding with sub-clustering for broad resolution: ", rb))

      # Get the names of your main clusters
      main_clusters <- levels(seurat@meta.data[, broad_clus_col_name])

      # Inner loop for different sub-clustering resolutions
      for (rs in res_sub) {
        print(paste0("  Sub-clustering at resolution: ", rs))

        # Create a new metadata column to store the final combined subclusters for this resolution
        subcluster_col_name <- paste0("subcluster_broad_res", rb, "_sub_res", rs)
        seurat@meta.data[, subcluster_col_name] <- as.character(seurat@meta.data[, broad_clus_col_name]) # Initialize

        for (cluster_id in main_clusters) {
          print(paste("    Subclustering broad cluster:", cluster_id, " / ", length(main_clusters)))

          # Subset the main Seurat object for the current broad cluster
          Idents(seurat) <- broad_clus_col_name # Set broad clusters as active identity
          sub_seurat <- subset(seurat, idents = cluster_id)

          # Skip if the subset is empty (e.g., a cluster might not exist after subsetting)
          if (ncol(sub_seurat) < 50) {
            print(paste("    Skipping subclustering for broad cluster", cluster_id, "as it contains too few cells."))
            next
          }

          # Re-processing for sub-clustering: FindVariableFeatures, ScaleData, RunPCA
          # These steps should be applied to the subsetted object for optimal sub-clustering.
          sub_seurat <- sub_seurat |>
            FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
            ScaleData(verbose = FALSE) |>
            RunPCA(dims = 1:50, verbose = FALSE)

          # Determine dimensions for sub-clustering. Ensure it doesn't exceed available PCs.
          # It's important to check if PCA actually returned components.
          if (is.null(sub_seurat@reductions$pca)) {
            print(paste("    Skipping subclustering for broad cluster", cluster_id, "due to PCA failure."))
            next
          }

          n_pcs_available <- ncol(sub_seurat@reductions$pca@cell.embeddings)
          if (n_pcs_available == 0) {
            print(paste("    Skipping subclustering for broad cluster", cluster_id, "due to 0 available PCs."))
            next
          }

          # Run FindNeighbors and FindClusters on the subset
          sub_seurat <- FindNeighbors(sub_seurat, dims = dims, verbose = FALSE)
          sub_seurat <- FindClusters(sub_seurat, resolution = rs, verbose = FALSE)

          sub_cluster_idents_col <- paste0("RNA_snn_res.", rs)
          sub_cluster_ids <- sub_seurat@meta.data[[sub_cluster_idents_col]]

          if (!is.factor(sub_cluster_ids)) {
            sub_cluster_ids <- factor(sub_cluster_ids)
          }

          # Rename subclusters to be unique (e.g., "BroadCluster_SubclusterID")
          # Ensure levels are handled correctly to preserve order and avoid issues with factor conversions.
          new_subcluster_names <- paste0(cluster_id, "_", levels(sub_cluster_ids))
          names(new_subcluster_names) <- levels(sub_cluster_ids)

          # Map the old levels to the new names
          mapped_sub_cluster_ids <- factor(sub_cluster_ids, levels = levels(sub_cluster_ids), labels = new_subcluster_names)

          # Store the subcluster information back into the main Seurat object's new column
          seurat@meta.data[Cells(sub_seurat), subcluster_col_name] <- as.character(mapped_sub_cluster_ids)
        }

        # After iterating through all broad clusters for a given 'rs',
        # set the new subcluster annotation as the active identity for plotting this specific result
        Idents(seurat) <- subcluster_col_name

        # Convert the new column to a factor if it's not already, for consistent plotting
        seurat@meta.data[, subcluster_col_name] <- as.factor(seurat@meta.data[, subcluster_col_name])
      }
    }
  }

  return(seurat)
}


get_ct_comp_df_seurat <- function(seurat, sample_col, ct_col) {
  ct_comp_df <- table(seurat@meta.data[[sample_col]], seurat@meta.data[[ct_col]]) %>%
    t() %>%
    as.data.frame.matrix() %>%
    t() %>%
    as.data.frame()
  ct_comp_df <- ct_comp_df[rowSums(ct_comp_df) != 0, ]

  return(ct_comp_df)
}


# get_ct_var and helpers ----


get_ct_var <- function(df,
                       show_plot = TRUE,
                       plot_title = "",
                       smooth_method = "lm",
                       descending = TRUE) {
  df_var <- df %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "celltype",
      values_to = "values"
    ) %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarize(
      Relative_abundance = mean(values, na.rm = TRUE),
      Variance = var(values, na.rm = TRUE)
    )

  if (descending) {
    df_var <- df_var %>%
      dplyr::arrange(dplyr::desc(Variance))
  } else {
    df_var <- df_var %>%
      dplyr::arrange(Variance)
  }

  if (show_plot) {
    p <- varmeanplot(data = df_var, title = plot_title, smooth_method = smooth_method)
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
    top_hvcs <- select_by_variance_explained(df_var, variance_threshold = variance_threshold)
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
    dplyr::arrange(dplyr::desc(Variance))

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
    dplyr::filter(variance_explained <= variance_threshold) %>%
    dplyr::pull(celltype)

  return(selected_celltypes)
}



varmeanplot <- function(data, plot_title = "", smooth_method = "lm", label_points = FALSE) {
  p <- ggplot(data, aes(x = Relative_abundance, y = Variance)) +
    geom_point() +
    geom_smooth(method = smooth_method, color = "red", fill = "#69b3a2", se = TRUE) +
    labs(title = paste(plot_title)) +
    theme_classic() +
    xlab("Mean") +
    ylab("Variance")

  if (label_points) {
    p <- p + ggrepel::geom_text_repel(data = data, aes(label = celltype), vjust = -0.5)
  }

  return(p)
}


# ----


get_pb <- function(seurat, sample_col = "Sample", hvg = NULL) {
  pb <- as.matrix(AggregateExpression(seurat, group.by = sample_col, assays = "RNA")[["RNA"]])
  colnames(pb) <- gsub("-", "_", colnames(pb))
  if (!is.null(hvg)) {
    pb <- pb[hvg, ]
  }
  return(pb)
}


get_pb_deseq2 <- function(seurat, sample_col = "Sample", hvg = NULL, n_hvg = 2000, black_list = "none") {
  pb <- get_pb(seurat, sample_col = sample_col, hvg = hvg)

  if (is.null(hvg) & black_list == "default") {
    black_list <- default_black_list # From SignatuR package
    black_list <- black_list[!names(black_list) %in% c("Xgenes", "Ygenes")]
    black_list <- unlist(black_list)

    pb <- pb[!rownames(pb) %in% black_list, ]
  }

  metadata <- get_metadata(seurat, sample_col = sample_col)
  metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- t(DESeq2.normalize(pb, metadata = metadata, n_hvg = n_hvg))
  return(pb_norm)
}



get_metadata <- function(seurat, sample_col = "Sample") {
  metadata <- seurat@meta.data %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::slice(1)

  return(metadata)
}

get_labels <- function(seurat, label_col) {
  metadata <- get_metadata(seurat)
  labels <- as.factor(metadata[[label_col]])
  names(labels) <- metadata[["Sample"]]

  return(labels)
}


# For ggplot x-axis label bolding (and optionally color)
library(glue)
highlight <- function(x, pat, color = "black", family = "") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}



impute_zeros <- function(df,
                         clr_zero_impute_method = c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all"),
                         clr_zero_impute_num = 1) {
  if (!clr_zero_impute_method %in% c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all")) {
    stop("clr_zero_impute_method not found")
  }

  # Apply specified zero imputation method
  if (clr_zero_impute_method == "percentage_zeros") {
    for (row in 1:nrow(df)) {
      df[row, ][df[row, ] == 0] <- sum(df[row, ]) / 100 * clr_zero_impute_num
    }
  } else if (clr_zero_impute_method == "percentage_all") {
    for (row in 1:nrow(df)) {
      df[row, ] <- df[row, ] + sum(df[row, ]) / 100 * clr_zero_impute_num
    }
  } else if (clr_zero_impute_method == "counts_zeros") {
    # Impute zeros by replacing them with a small non-zero value (1 in this case)
    df[df == 0] <- clr_zero_impute_num
  } else if (clr_zero_impute_method == "counts_all") {
    # Impute by adding a fixed count to all values
    df <- df + clr_zero_impute_num
  }

  return(df)
}


# Need to load h5ad without metadata and add it manually, otherwise R errors
load_h5ad_to_seurat <- function(file_name) {
  file_name_short <- strsplit(file_name, "\\.")[[1]][1]
  file_name_seurat <- paste0(file_name_short, ".h5seurat")

  if (file.exists(file_name) & !file.exists(file_name_seurat)) {
    Convert(file_name, dest = "h5seurat", overwrite = TRUE)
  }

  seurat <- SeuratDisk::LoadH5Seurat(file_name_seurat, meta.data = F)
  meta <- arrow::read_feather(paste0(file_name_short, ".feather"))
  seurat <- Seurat::AddMetaData(seurat, meta)

  return(seurat)
}


# Min-max-scaling score metrics
min_max <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  if (max(x, na.rm = T) == min(x, na.rm = T)) {
    return(x)
  }
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(feat_mat,
                     labels,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = 3,
                     title = NULL,
                     cluster_score = TRUE,
                     mod_score = TRUE,
                     sil_score = FALSE,
                     anosim_score = TRUE,
                     pointsize = 3,
                     labelsize = 4,
                     coord_equal = TRUE,
                     axes = c(1, 2),
                     plotly_3d = FALSE,
                     invisible = c("var", "quali"),
                     n_ct_show = Inf,
                     repel = FALSE) {
  res.pca <- prcomp(feat_mat, scale. = scale., rank. = pca_dims)


  if (anosim_score) {
    anosim_score <- round(vegan::anosim(x = feat_mat, grouping = labels, distance = "euclidean")[["statistic"]], 3)
    title <- paste0(title, "\nANOSIM score: ", anosim_score)
  }
  if (cluster_score) {
    cluster_score <- clust_eval(feat_mat, labels)
    title <- paste0(title, "\nARI: ", cluster_score)
  }
  if (mod_score) {
    mod_score <- round(calc_modularity(feat_mat, labels, knn_k), 3)
    title <- paste0(title, "\nModularity score: ", mod_score)
  }
  if (sil_score) {
    sil_score <- round(calc_sil(feat_mat, labels), 3)
    title <- paste0(title, "\nSilhouette score: ", sil_score)
  }

  if (plotly_3d) {
    df <- as.data.frame(res.pca$x)
    df$id <- seq_len(nrow(df))
    df$vs <- factor(labels)
    ms <- replicate(2, df, simplify = F)
    ms[[2]]$PC3 <- min(df$PC3)
    m <- ms %>%
      bind_rows() %>%
      plotly::group2NA("id", "vs")
    # Plotting with plotly
    p <- plotly::plot_ly(color = ~vs) %>%
      plotly::add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
      plotly::add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 0.2)
  } else {
    p <- factoextra::fviz_pca(res.pca,
      axes = axes,
      habillage = labels,
      label = "var",
      pointsize = pointsize,
      labelsize = labelsize,
      invisible = invisible,
      select.var = list(contrib = n_ct_show),
      repel = repel,
      geom = "point"
    ) +
      ggtitle(title) +
      scale_shape_manual(values = rep(19, length(unique(labels))))
    if (coord_equal) {
      p <- p + coord_equal()
    }
  }

  return(p)
}


# Plot rocauc curve separating samples by PCA loading 1
plot_roc_pcad1 <- function(feat_mat, labels) {
  res.pca <- prcomp(feat_mat)

  pca_scores_pc1 <- data.frame(
    PC1 = data.frame(res.pca$x)$PC1,
    label = factor(labels)
  )
  p <- ggplot(pca_scores_pc1, aes(x = PC1, fill = label)) +
    geom_density(alpha = 0.4) +
    geom_histogram(aes(y = ..density..),
      alpha = 0.3, colour = "black",
      position = "identity"
    ) +
    labs(
      title = "Density Plot of PC1 Scores",
      x = "PC1 Scores",
      y = "Density"
    ) +
    theme_classic()
  print(p)

  # Calculate the ROC curve
  roc_curve <- roc(pca_scores_pc1$label, ifelse(pca_scores_pc1$PC1 > 0, 1, 0))

  # Plot the ROC curve
  plot(roc_curve, col = "blue", main = "ROC Curve", print.auc = TRUE)
}


remove_low_cellcount_samples <- function(seurat,
                                         sample_col = "Sample",
                                         min_cells_per_sample = 500,
                                         show_plot = FALSE) {
  Idents(seurat) <- sample_col
  cells_per_sample <- table(seurat@meta.data[[sample_col]])
  seurat <- subset(seurat, idents = names(cells_per_sample[cells_per_sample > min_cells_per_sample]))

  if (show_plot) {
    barplot(sort(cells_per_sample))
    abline(v = min_cells_per_sample)
  }

  return(seurat)
}




# To replace HiTME layer3 annotation based on other UCell score threshold
replace_HiTMElayer3_annot <- function(seurat, thresh = 0.1) {
  seurat$layer3 <- gsub(pattern = "_resting|_IFN|_cellCycle.G1S|_cellCycle.G2M|_HeatShock|_Heatshock|_Prolif", replacement = "", seurat$layer3)

  # Define the simple conditions to loop through
  simple_conditions <- list(
    IFN_UCell = list(threshold = thresh, suffix = "_IFN") # ,
    # HeatShock_UCell = list(threshold = thresh, suffix = "_Heatshock")
  )

  # Loop through the simple conditions
  for (cond_col in names(simple_conditions)) {
    threshold <- simple_conditions[[cond_col]]$threshold
    suffix <- simple_conditions[[cond_col]]$suffix

    seurat@meta.data <- seurat@meta.data %>%
      mutate(
        layer3 = if_else(
          !is.na(layer3) & .data[[cond_col]] > threshold,
          paste0(layer3, suffix),
          layer3
        )
      )
  }

  # Now, handle the combined "Proliferation" condition separately
  seurat@meta.data <- seurat@meta.data %>%
    mutate(
      layer3 = if_else(
        !is.na(layer3) & (cellCycle.G1S_UCell > thresh | cellCycle.G2M_UCell > thresh),
        paste0(layer3, "_Prolif"),
        layer3
      )
    )

  return(seurat)
}


# run_analyses and helpers ----

run_analyses <- function(result_list,
                         ds,
                         seurat,
                         path_data,
                         path_plots) {
  print(paste("Running benchmark analysis for dataset: ", ds))
  result_list[["bmark"]][[ds]] <- run_benchmark_analysis(
    res_list = result_list[["bmark"]][[ds]],
    ds = ds,
    seurat = seurat,
    path_data = path_data,
    path_plots = path_plots
  )

  labels <- get_labels(seurat, seurat@misc$label_col)
  ct_comps <- get_ct_comp_df_seurat(seurat, sample_col = "Sample", ct_col = seurat@misc$hi_res_ct_col)

  if (is.null(result_list[["trans"]][[ds]])) {
    print(paste("Running transformation analysis for dataset: ", ds))
    result_list[["trans"]][[ds]] <- run_transformation_analysis(ct_comps, labels)
  }

  if (is.null(result_list[["zeroimp"]][[ds]])) {
    print(paste("Running zero imputation analysis for dataset: ", ds))
    result_list[["zeroimp"]][[ds]] <- run_zeroimp_analysis(ct_comps, labels)
  }

  saveRDS(result_list, file = "result_list.rds")

  return(result_list)
}


#----------------------------------------------------------->
run_benchmark_analysis <- function(res_list,
                                   ds,
                                   seurat,
                                   sample_col = "Sample",
                                   factors_test = c(2, 3, 5, 10, 15),
                                   path_data,
                                   path_plots,
                                   seurat_res = c(0.4, 2, 5, 20),
                                   HVGs = c(1000, 2000, 3000),
                                   # ECODA_select_top_hvct = TRUE,
                                   # ECODA_top_n_hvct = 0.25,
                                   ECODA_top_varexp_hvct = seq(0, 0.9, 0.1),
                                   gloscope_n_pca_dims = c(10, 30, 50)) {
  if (grepl("GongSharma", ds)) {
    ds_filename <- "GongSharma_all"
  } else {
    ds_filename <- ds
  }

  # Files preprocessed with python
  for (i in HVGs) {
    if (i == 2000) {
      scpoli_dims <- factors_test
    } else {
      scpoli_dims <- 5
    }

    file_mrvi <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_mrvi_dists.feather"))
    file_pilot <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_pilot_dists.feather"))
    files_scpoli <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_scpoli_dims", scpoli_dims, "_embs.feather"))

    files_to_check <- c(file_mrvi, file_pilot, files_scpoli)
    missing_files <- files_to_check[!file.exists(files_to_check)]

    if (length(missing_files) > 0) {
      stop("The following file(s) are missing:\n", paste(missing_files, collapse = "\n"))
    }
  }


  # Sample names starting with digits are not allowed in seurat
  seurat@meta.data[[sample_col]] <- standardize_sample_names(seurat@meta.data[[sample_col]])

  metadata <- get_metadata(seurat)

  label_col <- seurat@misc$label_col
  labels <- get_labels(seurat, label_col)

  hvg <- get_current_hvgs(seurat)


  if (!"Pseudobulk_hvg2000" %in% names(res_list)) {
    exec_time_pb_norm <- exec_time(
      pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
    )
    res_list[["Pseudobulk_hvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_hvg2000"]] <- process_pseudobulk_fig(pb_norm, labels)
    ) + exec_time_pb_norm
  }


  # These methods need pb_norm (with unsupervised HVGs)
  test_items <- c(
    "ECODA_deconv",
    paste0("Pseudobulk_", factors_test, "_PCA_dims"),
    paste0("MOFA_hvg2000_factors", factors_test)
  )
  exec_time_pb_norm <- exec_time(
    pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = NULL, n_hvg = 2000)
  )
  # So need to check if any of those methods need to be run
  # (mainly whether to calculate pb_norm or not)
  if (any(!"Pseudobulk_unsup_hvg2000" %in% names(res_list))) {
    res_list[["Pseudobulk_unsup_hvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_unsup_hvg2000"]] <- process_pseudobulk_fig(pb_norm, labels)
    ) + exec_time_pb_norm
  }

  exec_time_pb_norm_bl <- exec_time(
    pb_norm_bl <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = NULL, n_hvg = 2000)
  )
  # So need to check if any of those methods need to be run
  # (mainly whether to calculate pb_norm or not)
  if (any(!"Pseudobulk_unsup_hvg2000_bl" %in% names(res_list))) {
    res_list[["Pseudobulk_unsup_hvg2000_bl"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_unsup_hvg2000_bl"]] <- process_pseudobulk_fig(pb_norm_bl, labels)
    ) + exec_time_pb_norm_bl
  }


  if (!"Avg_PCA_embedding" %in% names(res_list)) {
    res_list[["Avg_PCA_embedding"]][["exec_time"]] <- exec_time(
      res_list[["Avg_PCA_embedding"]] <- process_avg_pca_embedding_fig(seurat, labels)
    ) + exec_time_pb_norm
  }

  # Deconvolute using EPIC
  if (!"ECODA_deconv" %in% names(res_list)) {
    res_list[["ECODA_deconv"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_deconv"]] <- process_deconv_fig(t(pb_norm), labels)
    )
  }


  # CoDA

  ## layer1: low res. cell types
  if (!is.null(seurat@misc$low_res_ct_col)) {
    res_list[["ECODA_authors_LR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_LR"]] <- process_coda_fig(
        seurat, labels,
        ct_col = seurat@misc$low_res_ct_col,
        title = "ECODA\nlow res."
      )
    )
  }

  ## layer2: high res. cell types
  if (!is.null(seurat@misc$hi_res_ct_col)) {
    res_list[["ECODA_authors_HR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR"]] <- process_coda_fig(
        seurat, labels,
        ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\nhigh res."
      )
    )
    res_list[["ECODA_authors_HR_NULL"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR_NULL"]] <- process_coda_fig(
        seurat, labels,
        ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\nhigh res.", shuffle_labels = TRUE
      )
    )

    for (varexp_hvc in ECODA_top_varexp_hvct) {
      ECODA_authors_HR_top_varexp_hvc <- paste0("ECODA_authors_HR_top_varexp", varexp_hvc)
      res_list[[ECODA_authors_HR_top_varexp_hvc]][["exec_time"]] <- exec_time(
        res_list[[ECODA_authors_HR_top_varexp_hvc]] <-
          process_coda_fig(
            seurat, labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = seurat@misc$hi_res_ct_col,
            title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%")
          )
      )

      ECODA_HiTME_HR_layer2_top_varexp_hvc <- paste0("ECODA_HiTME_HR_layer2_top_varexp", varexp_hvc)
      res_list[[ECODA_HiTME_HR_layer2_top_varexp_hvc]][["exec_time"]] <- exec_time(
        res_list[[ECODA_HiTME_HR_layer2_top_varexp_hvc]] <-
          process_coda_fig(
            seurat, labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = "layer2",
            title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%")
          )
      )

      ECODA_HiTME_HR_layer3_top_varexp_hvc <- paste0("ECODA_HiTME_HR_layer3_top_varexp", varexp_hvc)
      res_list[[ECODA_HiTME_HR_layer3_top_varexp_hvc]][["exec_time"]] <- exec_time(
        res_list[[ECODA_HiTME_HR_layer3_top_varexp_hvc]] <-
          process_coda_fig(
            seurat, labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = "layer3",
            title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%")
          )
      )
    }

    res_list[["Freq_highres"]][["exec_time"]] <- exec_time(
      res_list[["Freq_highres"]] <- process_coda_fig(
        seurat, labels,
        calc_clr = FALSE,
        ct_col = seurat@misc$hi_res_ct_col,
        title = "Cell type composition (%)\nhigh res."
      )
    )
  }

  res_list[["ECODA_authors_HR_3most_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3most_varcts"]] <-
      process_coda_fig(
        seurat, labels,
        ECODA_top_n_hvct = 3,
        var_ct_desc = TRUE, ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\n2 least var. cell types"
      )
  )

  res_list[["ECODA_authors_HR_2least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_2least_varcts"]] <-
      process_coda_fig(
        seurat, labels,
        ECODA_top_n_hvct = 2,
        var_ct_desc = FALSE, ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\n2 least var. cell types"
      )
  )

  res_list[["ECODA_authors_HR_3least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3least_varcts"]] <-
      process_coda_fig(
        seurat, labels,
        ECODA_top_n_hvct = 3,
        var_ct_desc = FALSE, ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\n2 least var. cell types"
      )
  )

  res_list[["ECODA_HiTME_HR_layer2"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer2"]] <-
      process_coda_fig(
        seurat, labels,
        ct_col = "layer2",
        title = "ECODA\nHiTME layer2"
      )
  )
  res_list[["ECODA_HiTME_HR_layer3"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer3"]] <-
      process_coda_fig(
        seurat, labels,
        ct_col = "layer3",
        title = "ECODA\nHiTME layer3"
      )
  )
  res_list[["ECODA_scATOMIC_HR"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_scATOMIC_HR"]] <-
      process_coda_fig(
        seurat, labels,
        ct_col = "scATOMIC_pred",
        title = "ECODA\nscATOMIC"
      )
  )


  # Analyze for all resolutions

  ## Ultra high res. cell type clusters based on Leiden clustering to artificially increase the number of cell types (clusters), e.g. to 250 cell types (clusters)
  for (r in seurat_res) {
    res_col_name <- paste0("RNA_snn_res.", r)
    ECODA_seuratres_r <- paste0("ECODA_seuratres_", r)
    res_list[[ECODA_seuratres_r]][["exec_time"]] <- exec_time(
      res_list[[ECODA_seuratres_r]] <-
        process_coda_fig(
          seurat, labels,
          ct_col = res_col_name,
          title = paste0("ECODA\nLeiden clustering ", r)
        )
    )
  }


  # Methods that use different number of factors (e.g. PCA or dims)

  ### Required for MOFA to run
  seurat@version <- package_version("3.1.5")

  for (i in factors_test) {
    # Pseudobulk with PCA
    pb_pca_i <- paste0("Pseudobulk_", i, "_PCA_dims")
    if (!pb_pca_i %in% names(res_list)) {
      res_list[[pb_pca_i]][["exec_time"]] <- exec_time(
        res_list[[pb_pca_i]] <-
          process_pseudobulk_fig(
            pb_norm, labels,
            pca_dims = i,
            title = paste0("Pseudobulk gene expression\n+ PCA (", i, " dims)")
          )
      ) + exec_time_pb_norm
    }

    # Hires CODA with PCA
    ecoda_pca_i <- paste0("ECODA_authors_HR_", i, "_PCA_dims")
    if (!ecoda_pca_i %in% names(res_list)) {
      res_list[[ecoda_pca_i]][["exec_time"]] <- exec_time(
        res_list[[ecoda_pca_i]] <-
          process_coda_fig(
            seurat, labels,
            pca_dims = i,
            ct_col = seurat@misc$hi_res_ct_col,
            title = paste0("ECODA\nhigh res. + PCA (", i, " dims)")
          )
      )
    }

    # MOFA
    mofa_factor_i <- paste0("MOFA_hvg2000_factors", i)
    if (!mofa_factor_i %in% names(res_list)) {
      res_list[[mofa_factor_i]][["exec_time"]] <- exec_time(
        res_list[[mofa_factor_i]] <-
          process_mofa_bulk_fig(
            pb_norm,
            metadata = metadata, labels,
            num_factors = i
          )
      ) + exec_time_pb_norm
    }

    # scITD
    scitd_factor_i <- paste0("scITD_hvg2000_factors", i)
    if (!scitd_factor_i %in% names(res_list)) {
      res_list[[scitd_factor_i]][["exec_time"]] <- exec_time(
        res_list[[scitd_factor_i]] <-
          process_scitd_fig(
            seurat,
            ct_col = seurat@misc$low_res_ct_col,
            label_col = label_col, hvg,
            num_factors = i
          )
      )
    }
  }

  # GloScope
  ## With different numbers of PCA dims

  for (i in gloscope_n_pca_dims) {
    gloscope_dist_file <- file.path(path_data, paste0(ds_filename, "_gloscope_hvg2000_pcadims", i, "_dists.rds"))
    gloscope_pca_i <- paste0("GloScope_hvg2000_pcadims", i)
    if (!gloscope_pca_i %in% names(res_list)) {
      res_list[[gloscope_pca_i]][["exec_time"]] <- exec_time(
        res_list[[gloscope_pca_i]] <-
          process_gloscope_fig(
            seurat, metadata, label_col,
            gloscope_dist_file = gloscope_dist_file,
            n_pca_dims = i
          )
      )
      res_list[[gloscope_pca_i]][["exec_time"]] <- exec_time(
        res_list[[paste0(gloscope_pca_i, "_sqrtmat")]] <-
          process_gloscope_sqrtmat_fig(
            metadata, label_col,
            gloscope_dist_file = gloscope_dist_file
          )
      ) + res_list[[gloscope_pca_i]][["exec_time"]]
    }
  }


  # Methods that use different number of HVGs

  ### Do only for non-default HVGs (default = 2000)
  for (i in HVGs[!HVGs %in% 2000]) {
    Pseudobulk_hvg_i <- paste0("Pseudobulk_hvg", i)
    Pseudobulk_unsup_hvg_i <- paste0("Pseudobulk_unsup_hvg", i)
    MOFA_hvg_i_15_factors <- paste0("MOFA_hvg", i, "_15_factors")
    GloScope_hvg_i_pcadims30 <- paste0("GloScope_hvg", i, "_pcadims30")
    GloScope_hvg_i_pcadims30_sqrtmat <- paste0(GloScope_hvg_i_pcadims30, "_sqrtmat")
    test_items <- c(
      Pseudobulk_hvg_i,
      Pseudobulk_unsup_hvg_i,
      MOFA_hvg_i_15_factors,
      GloScope_hvg_i_pcadims30,
      GloScope_hvg_i_pcadims30_sqrtmat
    )

    if (any(!test_items %in% names(res_list))) {
      # Memory critical steps
      gc()
      seurat <- create_clean_seuratv5_object(seurat)
      gc()
      seurat <- NormalizeData(seurat)
      gc()
      seurat <- FindVariableFeatures(seurat, nfeatures = i)
      gc()
      seurat <- ScaleData(seurat)
      gc()
      # Needs a lot of memory for 3000 HVGs
      seurat <- RunPCA(seurat, dims = 1:50, verbose = FALSE)
      gc()
      hvg <- get_current_hvgs(seurat)

      if (!Pseudobulk_hvg_i %in% names(res_list)) {
        exec_time_pb_norm <- exec_time(
          pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
        )
        res_list[[Pseudobulk_hvg_i]][["exec_time"]] <- exec_time(
          res_list[[Pseudobulk_hvg_i]] <- process_pseudobulk_fig(pb_norm, labels)
        ) + exec_time_pb_norm
      }

      test_items <- c(
        Pseudobulk_unsup_hvg_i,
        MOFA_hvg_i_15_factors
      )
      if (any(!test_items %in% names(res_list))) {
        exec_time_pb_norm <- exec_time(
          pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = NULL, n_hvg = i)
        )
        res_list[[Pseudobulk_unsup_hvg_i]][["exec_time"]] <- exec_time(
          res_list[[Pseudobulk_unsup_hvg_i]] <- process_pseudobulk_fig(pb_norm, labels)
        ) + exec_time_pb_norm

        res_list[[MOFA_hvg_i_15_factors]][["exec_time"]] <- exec_time(
          res_list[[MOFA_hvg_i_15_factors]] <-
            process_mofa_bulk_fig(
              pb_norm,
              metadata = metadata, labels,
              num_factors = 15
            )
        ) + exec_time_pb_norm
      }

      test_items <- c(
        GloScope_hvg_i_pcadims30,
        GloScope_hvg_i_pcadims30_sqrtmat
      )
      if (any(!test_items %in% names(res_list))) {
        gloscope_dist_file <- file.path(path_data, paste0(ds_filename, "_gloscope_hvg", i, "_pcadims50_dists.rds"))
        res_list[[GloScope_hvg_i_pcadims30]][["exec_time"]] <- exec_time(
          res_list[[GloScope_hvg_i_pcadims30]] <-
            process_gloscope_fig(
              seurat, metadata, label_col,
              gloscope_dist_file = gloscope_dist_file
            )
        )
        res_list[[GloScope_hvg_i_pcadims30_sqrtmat]][["exec_time"]] <- exec_time(
          res_list[[GloScope_hvg_i_pcadims30_sqrtmat]] <-
            process_gloscope_sqrtmat_fig(
              metadata, label_col,
              gloscope_dist_file = gloscope_dist_file
            )
        ) + res_list[[GloScope_hvg_i_pcadims30]][["exec_time"]]
      }
    }
  }


  for (i in HVGs) {
    # --- MrVI (Runs once per HVG) ---
    mrvi_dist_file <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_mrvi_dists.feather"))
    res_list[[paste0("MrVI_hvg", i)]] <- process_mrvi_fig(mrvi_dist_file = mrvi_dist_file, labels)

    # --- PILOT (Runs once per HVG) ---
    pilot_dist_file <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_pilot_dists.feather"))
    res_list[[paste0("PILOT_hvg", i)]] <- process_pilot_fig(pilot_dist_file = pilot_dist_file, labels)

    # --- scPoli (Runs once OR multiple times depending on HVG) ---
    if (i == 2000) {
      target_dims <- factors_test
    } else {
      target_dims <- 5
    }

    for (f in target_dims) {
      scpoli_emb_file <- file.path(path_data, paste0(ds_filename, "_hvg", i, "_scpoli_dims", f, "_embs.feather"))
      res_list[[paste0("scPoli_hvg", i, "_dims", f)]] <- process_scpoli_fig(scpoli_emb_file = scpoli_emb_file, labels)
    }
  }


  # ECODA + PB
  for (i in c(0, 0.25, 0.5, 0.75, 1)) {
    for (norm in c("max", "median", "zscore", "quantile")) {
      res_list[[paste0("ECODA_PB_combo_norm", norm, "_ecodaweight", i)]] <- process_ecodapb_fig(
        dist_mat_ecoda = res_list[["ECODA_authors_HR"]][["dist_mat"]],
        dist_mat_pb = res_list[["Pseudobulk_unsup_hvg2000"]][["dist_mat"]],
        norm_method = norm,
        ecoda_weight = i,
        labels = res_list[["ECODA_authors_HR"]][["labels"]],
      )
      res_list[[paste0("ECODA_PB_combo_norm", norm, "_ecodaweight", i)]][["exec_time"]] <-
        res_list[["ECODA_authors_HR"]][["exec_time"]] + res_list[["Pseudobulk_unsup_hvg2000"]][["exec_time"]]
    }
  }


  return(res_list)
}


get_current_hvgs <- function(seurat) {
  if ("var.features" %in% slotNames(seurat@assays[["RNA"]])) {
    return(seurat@assays[["RNA"]]@var.features)
  } else if ("var.features" %in% colnames(seurat@assays[["RNA"]]@meta.data)) {
    vf <- seurat@assays[["RNA"]]@meta.data[["var.features"]]
    return(vf[!is.na(vf)])
  } else {
    stop("Could not find variable features.")
  }
}


standardize_sample_names <- function(sample_names) {
  # Sample names starting with digits are not allowed in seurat
  sample_names_starting_with_digit <- grepl("^\\d", unique(sample_names))
  if (any(sample_names_starting_with_digit)) {
    sample_names[sample_names_starting_with_digit] <- paste0("g", sample_names[sample_names_starting_with_digit])
  }
  sample_names <- gsub("-", "_", sample_names)
}


process_deconv_fig <- function(pseudobulk,
                               labels,
                               title = "Deconvoluted cell type composition") {
  out <- EPIC(pseudobulk, BRef)
  deconv_ct_comps <- as.data.frame(out[["mRNAProportions"]])
  row.names(deconv_ct_comps) <- colnames(pseudobulk)

  feat_mat <- clr(deconv_ct_comps + 0.001)
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}


process_coda_fig <- function(seurat,
                             labels,
                             ECODA_top_n_hvct = NULL,
                             ECODA_top_varexp_hvct = NULL,
                             hvct_recalc_clr = TRUE,
                             calc_clr = TRUE,
                             pca_dims = NULL,
                             sample_col = "Sample",
                             ct_col,
                             title,
                             clr_zero_impute_method = "counts_all",
                             clr_zero_impute_num = 1,
                             feat_mat = NULL,
                             var_ct_desc = TRUE,
                             shuffle_labels = FALSE) {
  if (is.null(feat_mat)) {
    df_counts <- get_ct_comp_df_seurat(seurat, sample_col = sample_col, ct_col)

    df_imp <- df_counts %>%
      impute_zeros(
        clr_zero_impute_method = clr_zero_impute_method,
        clr_zero_impute_num = clr_zero_impute_num
      ) %>%
      calc_perc_df()

    if (calc_clr) {
      feat_mat <- df_imp %>%
        clr()
    } else {
      feat_mat <- df_imp
    }

    top_hvct <- NULL
    if (!is.null(ECODA_top_n_hvct)) {
      top_hvct <- get_ct_var(feat_mat, show_plot = FALSE, descending = var_ct_desc) %>%
        get_hvcs(top_n_hvcs = ECODA_top_n_hvct, variance_threshold = NULL)
      feat_mat <- feat_mat[, top_hvct]
    }
    if (!is.null(ECODA_top_varexp_hvct)) {
      top_hvct <- get_ct_var(feat_mat, show_plot = FALSE, descending = var_ct_desc) %>%
        get_hvcs(top_n_hvcs = NULL, variance_threshold = ECODA_top_varexp_hvct)
      feat_mat <- feat_mat[, top_hvct]
    }

    if (hvct_recalc_clr & !is.null(top_hvct)) {
      feat_mat <- df_counts[, top_hvct] %>%
        impute_zeros(clr_zero_impute_method = clr_zero_impute_method, clr_zero_impute_num = clr_zero_impute_num) %>%
        calc_perc_df() %>%
        clr()
    }
  }

  if (!is.null(pca_dims)) {
    feat_mat <- prcomp(feat_mat, rank. = pca_dims)[["x"]]
  }

  labels <- labels[rownames(feat_mat)]

  if (shuffle_labels) {
    set.seed(123)
    labels <- sample(labels)
  }

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}


process_pseudobulk_fig <- function(feat_mat,
                                   labels,
                                   pca_dims = NULL,
                                   title = "Pseudobulk gene expression",
                                   knn_k = NULL) {
  if (!is.null(pca_dims)) {
    feat_mat <- prcomp(feat_mat, rank. = pca_dims)[["x"]]
  }
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}


process_avg_pca_embedding_fig <- function(seurat,
                                          labels,
                                          sample_col = "Sample",
                                          title = "Avg_PCA_embedding") {
  feat_mat <- as.data.frame(seurat@reductions$pca@cell.embeddings)
  feat_mat$Sample <- seurat@meta.data[[sample_col]]

  feat_mat <- feat_mat %>%
    group_by(Sample) %>%
    summarise(across(starts_with("PC_"), mean)) %>%
    ungroup() %>%
    column_to_rownames(var = "Sample")
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}



process_mofa_bulk_fig <- function(pb_norm,
                                  metadata,
                                  labels,
                                  title = "MOFA",
                                  num_factors = 5,
                                  maxiter = 1000) {
  title <- paste0(title, " (", num_factors, " factors)")

  pb_list <- list(pb = t(pb_norm))

  MOFAobject <- create_mofa(pb_list)

  metadata$sample <- metadata$Sample

  samples_metadata(MOFAobject) <- metadata

  # Default data options
  data_opts <- get_default_data_options(MOFAobject)

  # Default model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- num_factors

  # Training options
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "fast"
  train_opts$seed <- 42
  train_opts$maxiter <- maxiter

  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  MOFAobject <- run_mofa(MOFAobject, use_basilisk = FALSE, save_data = FALSE)

  # feat_mat <- get_factors(MOFAobject, as.data.frame = T)
  feat_mat <- as.data.frame(MOFAobject@expectations[["Z"]])
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title, coord_equal = FALSE)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}



process_scitd_fig <- function(seurat,
                              ct_col,
                              label_col,
                              hvg,
                              num_factors = 5,
                              title = "scITD") {
  title <- paste0(title, " (", num_factors, " factors)")

  seurat$donors <- seurat$Sample
  seurat$ctypes <- as.character(seurat@meta.data[[ct_col]])

  # if (!file.exists(file.path(path_output, "pbmc_container_cmv.rds"))) {
  ctypes <- unique(seurat$ctypes)
  ctypes <- ctypes[!is.na(ctypes)]

  # Step 1: Calculate table of donors by ctypes
  cell_counts <- table(seurat$donors, seurat$ctypes)

  # Step 2: Determine low abundant cell types (less than 5 cells in more than 20% of donors)
  ctypes_to_drop <- names(which(colMeans(cell_counts < 5) > 0.2))

  ctypes <- ctypes[!ctypes %in% ctypes_to_drop]

  # set up project parameters
  param_list <- initialize_params(
    ctypes_use = as.character(ctypes),
    ncores = 6, rand_seed = 10
  )

  # create project container
  pbmc_container <- make_new_container(
    count_data = seurat@assays$RNA$counts,
    meta_data = seurat@meta.data[, c("donors", "ctypes", label_col, "Sample")],
    gn_convert = NULL,
    params = param_list
  )

  # form the tensor from the data
  pbmc_container <- form_tensor(pbmc_container,
    custom_genes = hvg
  )

  # run the tensor decomposition
  pbmc_container <- run_tucker_ica(pbmc_container,
    ranks = c(num_factors, num_factors + 5)
  )

  feat_mat <- pbmc_container[["tucker_results"]][[1]]
  labels_scITD <- seurat@meta.data %>%
    filter(donors %in% row.names(feat_mat)) %>%
    distinct(donors, .keep_all = TRUE) %>%
    dplyr::select(donors, !!sym(label_col))
  labels <- as.factor(labels_scITD[[label_col]])
  names(labels) <- labels_scITD[["donors"]]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}


process_mrvi_fig <- function(mrvi_dist_file, labels, title = "MrVI") {
  feat_mat <- arrow::read_feather(mrvi_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()

  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)
  res[["labels"]] <- labels

  return(res)
}


process_gloscope_fig <- function(seurat,
                                 metadata,
                                 label_col,
                                 gloscope_dist_file,
                                 n_pca_dims = 30,
                                 dens = "KNN",
                                 dist_mat = c("KL"),
                                 k = 25,
                                 BPPARAM = BiocParallel::MulticoreParam(workers = parallelly::availableCores() - 2, progressbar = TRUE),
                                 title = "GloScope") {
  if (!file.exists(gloscope_dist_file)) {
    feat_mat <- GloScope::gloscope(
      embedding_matrix = seurat@reductions$pca@cell.embeddings[, 1:n_pca_dims],
      cell_sample_ids = seurat$Sample,
      dens = dens,
      dist_mat = dist_mat,
      k = k,
      BPPARAM = BPPARAM
    )
    saveRDS(feat_mat, file = gloscope_dist_file)
  } else {
    feat_mat <- readRDS(gloscope_dist_file)
  }

  samples <- row.names(feat_mat)
  samples <- standardize_sample_names(samples)

  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  names(labels) <- metadata[match(samples, metadata$Sample), ][["Sample"]]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- feat_mat
  res[["labels"]] <- labels

  return(res)
}


process_gloscope_sqrtmat_fig <- function(metadata,
                                         label_col,
                                         gloscope_dist_file,
                                         title = "GloScope") {
  if (!file.exists(gloscope_dist_file)) {
    stop(paste(gloscope_dist_file, "not found!"))
  } else {
    feat_mat <- readRDS(gloscope_dist_file)
  }

  feat_mat <- sqrt(feat_mat)

  # As suggested here:
  # https://github.com/epurdom/GloScope/issues/3
  feat_mat[is.na(feat_mat)] <- 0

  samples <- row.names(feat_mat)
  samples <- standardize_sample_names(samples)

  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  names(labels) <- metadata[match(samples, metadata$Sample), ][["Sample"]]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- feat_mat
  res[["labels"]] <- labels

  return(res)
}


process_scpoli_fig <- function(scpoli_emb_file, labels, title = "scPoli") {
  feat_mat <- arrow::read_feather(scpoli_emb_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()

  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  res[["labels"]] <- labels

  return(res)
}


process_pilot_fig <- function(pilot_dist_file, labels, title = "PILOT") {
  feat_mat <- arrow::read_feather(pilot_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()

  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  labels <- labels[rownames(feat_mat)]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)
  res[["labels"]] <- labels

  return(res)
}


process_ecodapb_fig <- function(
  dist_mat_ecoda,
  dist_mat_pb,
  feat_mat_ecoda = NULL,
  feat_mat_pb = NULL,
  ecoda_weight = 0.5,
  norm_method = c("max", "median", "zscore", "quantile"),
  labels
) {
  norm_method <- match.arg(norm_method)

  if (!is.null(feat_mat_ecoda) & !is.null(feat_mat_pb)) {
    # Cosine distance matrix
    dist_mat_ecoda_normed <- proxy::dist(feat_mat_ecoda, method = "cosine")
    dist_mat_pb_normed <- proxy::dist(feat_mat_pb, method = "cosine")
  } else {
    # Euclidian distance matrix
    if (norm_method == "max") {
      dist_mat_ecoda_normed <- dist_mat_ecoda / max(dist_mat_ecoda)
      dist_mat_pb_normed <- dist_mat_pb / max(dist_mat_pb)
    } else if (norm_method == "median") {
      dist_mat_ecoda_normed <- dist_mat_ecoda / median(dist_mat_ecoda)
      dist_mat_pb_normed <- dist_mat_pb / median(dist_mat_pb)
    } else if (norm_method == "zscore") {
      dist_mat_ecoda_normed <- zscore_transform(dist_mat_ecoda)
      dist_mat_pb_normed <- zscore_transform(dist_mat_pb)
    } else if (norm_method == "quantile") {
      dist_mat_ecoda_normed <- global_quantile_norm_gaussian(dist_mat_ecoda)
      dist_mat_pb_normed <- global_quantile_norm_gaussian(dist_mat_pb)
    }
  }


  feat_mat <- dist_mat_ecoda_normed * ecoda_weight + dist_mat_pb_normed * (1 - ecoda_weight)

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)
  res[["labels"]] <- labels

  return(res)
}


zscore_transform <- function(dist_mat) {
  mu <- mean(dist_mat)
  sigma <- sd(dist_mat)
  z_score_matrix <- (dist_mat - mu) / sigma

  return(z_score_matrix)
}


global_quantile_norm_gaussian <- function(dist_mat) {
  # 1. Flatten the matrix to a single vector
  v <- as.vector(dist_mat)

  # 2. Calculate ranks (handling ties by averaging)
  r <- rank(v, ties.method = "average")

  # 3. Convert ranks to probabilities (Blom's method to avoid Inf)
  # (r - 0.5) / n is a standard formula to center probabilities
  probs <- (r - 0.5) / length(v)

  # 4. Apply Inverse Normal CDF (qnorm)
  v_norm <- qnorm(probs)

  # 5. Reshape back to original matrix dimensions
  mat_norm <- matrix(v_norm, nrow = nrow(dist_mat), ncol = ncol(dist_mat))

  # Restore names if they existed
  dimnames(mat_norm) <- dimnames(dist_mat)

  return(mat_norm)
}

#-----------------------------------------------------------<


run_transformation_analysis <- function(ct_comps, labels) {
  res_list <- datrans(ct_comps, labels,
    Amount_of_perturbation = 0,
    n_ct_to_select = 0,
    reps = 20,
    n_cores = 1,
    trans_method = c(
      "counts",
      "freq",
      "arcsine_sqrt",
      "alr_randref",
      "alr_mincvref",
      "clr"
    )
  )
  res_list <- res_list %>%
    dplyr::group_by(trans_method) %>%
    summarize(
      ANOSIM_score = mean(ANOSIM_score),
      # Silhouette_score = mean(Silhouette_score),
      Modularity_score = mean(Modularity_score),
      Adjusted_Rand_Index = mean(Adjusted_Rand_Index)
    ) %>%
    ungroup()
  res_list$dataset <- ds

  return(res_list)
}


run_zeroimp_analysis <- function(ct_comps, labels) {
  df <- ct_comps %>%
    select_if(colSums(.) != 0) %>%
    mutate_all(as.numeric)

  res_list <- list()
  res_list[["counts_zeros_1"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "counts_zeros", clr_zero_impute_num = 1) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["counts_all_1"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "counts_all", clr_zero_impute_num = 1) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.001%"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.001) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.01%"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.01) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.1%"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.1) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["perc_all_1%"]] <- df %>%
    impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 1) %>%
    calc_perc_df() %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["asinsqrt"]] <- df %>%
    calc_perc_df() %>%
    mutate(across(everything(), ~ . / 100)) %>%
    sqrt() %>%
    asin() %>%
    calc_sep_score(labels)
  # res_list[["impRZilr"]] <- df %>% {impRZilr(.)[["x"]]} %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  # Need to check if multLN drops samples
  df_multLN <- df %>%
    calc_perc_df() %>%
    zCompositions::multLN(label = 0, dl = rep(0.1, ncol(df)), z.warning = 0.9)
  labels_multLN <- labels[row.names(df) %in% row.names(df_multLN)]
  res_list[["multLN"]] <- df_multLN %>%
    clr() %>%
    calc_sep_score(labels_multLN)
  res_list[["multRepl_0.01%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(label = 0, dl = rep(0.01, ncol(df)), z.warning = 1, frac = 1) %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["multRepl_0.1%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(label = 0, dl = rep(0.1, ncol(df)), z.warning = 1, frac = 1) %>%
    clr() %>%
    calc_sep_score(labels)
  res_list[["multRepl_1%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(label = 0, dl = rep(1, ncol(df)), z.warning = 1, frac = 1) %>%
    clr() %>%
    calc_sep_score(labels)

  return(res_list)
}
