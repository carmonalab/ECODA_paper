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
  # devtools::install_github('satijalab/seurat-data')
  library(SeuratData)
  # install.packages("hdf5r")
  # remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
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


calc_sep_score <- function(df,
                           labels,
                           knn_k = NULL) {
  # sil_score <- round(calc_sil(df, labels), 3)
  mod_score <- unlist(round(calc_modularity(df, labels, knn_k), 3))
  cluster_score <- clust_eval(matrix = df, labels)
  anosim_score <- vegan::anosim(x = df, grouping = labels, distance = "euclidean")[["statistic"]]

  res <- list(
    # sil_score = sil_score,
    mod_score = mod_score,
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


get_pb <- function(seurat, sample_col = "Sample", hvg = NULL) {
  pb <- as.matrix(AggregateExpression(seurat, group.by = sample_col, assays = "RNA")[["RNA"]])
  colnames(pb) <- gsub("-", "_", colnames(pb))
  if (!is.null(hvg)) {
    pb <- pb[hvg, ]
  }
  return(pb)
}


get_pb_deseq2 <- function(seurat, sample_col = "Sample", hvg = NULL, n_hvg = 2000) {
  pb <- get_pb(seurat, sample_col = sample_col, hvg = hvg)
  metadata <- get_metadata(seurat)
  metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- t(DESeq2.normalize(pb, metadata = metadata, n_hvg = n_hvg))
  return(pb_norm)
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


get_metadata <- function(seurat, sample_col = "Sample") {
  metadata <- seurat@meta.data %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::slice(1)

  return(metadata)
}

get_labels <- function(seurat, label_col) {
  metadata <- get_metadata(seurat)
  labels <- as.factor(metadata[[label_col]])

  return(labels)
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


# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(feat_mat,
                     labels,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = NULL,
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


run_analyses <- function(result_list,
                         ds,
                         seurat,
                         path_data,
                         path_plots,
                         factors_test) {
  if (is.null(result_list[["bmark"]][[ds]])) {
    print(paste("Running benchmark analysis for dataset: ", ds))
    result_list[["bmark"]][[ds]] <- run_benchmark_analysis(
      seurat = seurat,
      ds = ds,
      path_data = path_data,
      path_plots = path_plots,
      factors_test = factors_test
    )
  }

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
run_benchmark_analysis <- function(seurat,
                                   ds,
                                   sample_col = "Sample",
                                   factors_test,
                                   path_data,
                                   path_plots,
                                   Pseudobulk = TRUE,
                                   ECODA_deconv = TRUE,
                                   ECODA_low_res = TRUE,
                                   ECODA_high_res = TRUE,
                                   ECODA_seuratcluster_res = TRUE,
                                   # ECODA_select_top_hvct = TRUE,
                                   # ECODA_top_n_hvct = 0.25,
                                   ECODA_top_varexp_hvct = seq(0, 0.9, 0.1),
                                   Pseudobulk_PCA = TRUE,
                                   ECODA_high_res_PCA = TRUE,
                                   MOFA = TRUE,
                                   scITD = TRUE,
                                   MrVI = TRUE,
                                   GloScope = TRUE,
                                   gloscope_n_pca_dims = c(10, 30, 50),
                                   scPoli = TRUE,
                                   PILOT = TRUE,
                                   show_pca_plots = FALSE,
                                   save_pca_plots = TRUE) {
  res_list <- list()

  # Files preprocessed with python
  for (i in c(1000, 2000, 3000)) {
    if (grepl("GongSharma", ds)) {
      ds_unif <- "GongSharma_all"
      files <- list(
        mrvi_dist_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_mrvi_dists.feather")),
        scpoli_emb_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_scpoli_embs.feather")),
        pilot_dist_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_pilot_dists.feather"))
      )
    } else {
      files <- list(
        mrvi_dist_file = file.path(path_data, paste0(ds, "_hvg", i, "_mrvi_dists.feather")),
        scpoli_emb_file = file.path(path_data, paste0(ds, "_hvg", i, "_scpoli_embs.feather")),
        pilot_dist_file = file.path(path_data, paste0(ds, "_hvg", i, "_pilot_dists.feather"))
      )
    }

    for (file in files) {
      if (!file.exists(file)) {
        stop("File not found: ", file)
      }
    }
  }

  # Sample names starting with digits are not allowed in seurat
  sample_names_starting_with_digit <- grepl("^\\d", unique(seurat@meta.data[[sample_col]]))
  if (any(sample_names_starting_with_digit)) {
    seurat@meta.data[[sample_col]][grepl("^\\d", seurat@meta.data[[sample_col]])] <- paste0("g", seurat@meta.data[[sample_col]][grepl("^\\d", seurat@meta.data[[sample_col]])])
  }

  metadata <- get_metadata(seurat)
  metadata[[sample_col]] <- gsub("-", "_", metadata[[sample_col]])

  label_col <- seurat@misc$label_col
  labels <- get_labels(seurat, label_col)
  names(labels) <- metadata[[sample_col]]

  if ("var.features" %in% slotNames(seurat@assays[["RNA"]])) {
    # This path is for older Seurat objects (primarily v2/v3)
    hvg <- seurat@assays[["RNA"]]@var.features
  } else if ("var.features" %in% colnames(seurat@assays[["RNA"]]@meta.data)) {
    varf <- seurat@assays[["RNA"]]@meta.data[["var.features"]]
    hvg <- varf[!is.na(varf)]
  } else {
    stop("Could not find variable features in expected locations.")
  }

  if (Pseudobulk) {
    exec_time_1 <- exec_time(
      pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
    )
    exec_time_2 <- exec_time(
      res_list[["Pseudobulk_hvg2000"]] <- process_pseudobulk_fig(pb_norm, labels)
    )
    res_list[["Pseudobulk_hvg2000"]][["exec_time"]] <- exec_time_1 + exec_time_2

    pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = NULL, n_hvg = 2000)
    res_list[["Pseudobulk_unsup_hvg2000"]] <- process_pseudobulk_fig(pb_norm, labels)
  }

  res_list[["Avg_PCA_embedding"]] <- process_avg_pca_embedding_fig(seurat, labels)

  if (ECODA_deconv) {
    # Deconvolute using EPIC
    res_list[["ECODA_deconv"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_deconv"]] <- process_deconv_fig(t(pb_norm), labels)
    )
  }

  # CoDA
  if (ECODA_low_res) {
    ## layer1: low res. cell types
    if (!is.null(seurat@misc$low_res_ct_col)) {
      res_list[["ECODA_authors_LR"]] <- process_coda_fig(seurat, labels, ct_col = seurat@misc$low_res_ct_col, title = "ECODA\nlow res.")
    }
  }

  if (ECODA_high_res) {
    ## layer2: high res. cell types
    if (!is.null(seurat@misc$hi_res_ct_col)) {
      res_list[["ECODA_authors_HR"]][["exec_time"]] <- exec_time(
        res_list[["ECODA_authors_HR"]] <- process_coda_fig(seurat, labels, ct_col = seurat@misc$hi_res_ct_col, title = "ECODA\nhigh res.")
      )
      res_list[["ECODA_authors_HR_NULL"]] <- process_coda_fig(seurat, labels, ct_col = seurat@misc$hi_res_ct_col, title = "ECODA\nhigh res.", shuffle_labels = TRUE)

      for (varexp_hvc in ECODA_top_varexp_hvct) {
        res_list[[paste0("ECODA_authors_HR_top_varexp", varexp_hvc)]] <-
          process_coda_fig(seurat, labels, ECODA_top_varexp_hvct = varexp_hvc, ct_col = seurat@misc$hi_res_ct_col, title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%"))

        res_list[[paste0("ECODA_HiTME_HR_layer2_top_varexp", varexp_hvc)]] <-
          process_coda_fig(seurat, labels, ECODA_top_varexp_hvct = varexp_hvc, ct_col = "layer2", title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%"))
        res_list[[paste0("ECODA_HiTME_HR_layer3_top_varexp", varexp_hvc)]] <-
          process_coda_fig(seurat, labels, ECODA_top_varexp_hvct = varexp_hvc, ct_col = "layer3", title = paste0("ECODA\nhigh res. var. exp. ", varexp_hvc * 100, "%"))
      }

      res_list[["Freq_highres"]] <- process_coda_fig(seurat, labels, clr = FALSE, ct_col = seurat@misc$hi_res_ct_col, title = "Cell type composition (%)\nhigh res.")
    }

    res_list[["ECODA_authors_HR_3most_varcts"]] <-
      process_coda_fig(seurat, labels, ECODA_top_n_hvct = 3, var_ct_desc = TRUE, ct_col = seurat@misc$hi_res_ct_col, title = "ECODA\n2 least var. cell types")

    res_list[["ECODA_authors_HR_2least_varcts"]] <-
      process_coda_fig(seurat, labels, ECODA_top_n_hvct = 2, var_ct_desc = FALSE, ct_col = seurat@misc$hi_res_ct_col, title = "ECODA\n2 least var. cell types")

    res_list[["ECODA_authors_HR_3least_varcts"]] <-
      process_coda_fig(seurat, labels, ECODA_top_n_hvct = 3, var_ct_desc = FALSE, ct_col = seurat@misc$hi_res_ct_col, title = "ECODA\n2 least var. cell types")

    res_list[["ECODA_HiTME_HR"]] <- process_coda_fig(seurat, labels, ct_col = "layer2", title = "ECODA\nHiTME")
    res_list[["ECODA_HiTME_HR_layer2"]] <- process_coda_fig(seurat, labels, ct_col = "layer2", title = "ECODA\nHiTME layer2")
    res_list[["ECODA_HiTME_HR_layer3"]] <- process_coda_fig(seurat, labels, ct_col = "layer3", title = "ECODA\nHiTME layer3")
    res_list[["ECODA_scATOMIC_HR"]] <- process_coda_fig(seurat, labels, ct_col = "scATOMIC_pred", title = "ECODA\nscATOMIC")
  }

  # Analyze for all resolutions

  if (ECODA_seuratcluster_res) {
    ## Ultra high res. cell type clusters based on Leiden clustering to artificially increase the number of cell types (clusters), e.g. to 250 cell types (clusters)
    # Resolutions for seurat clustering
    res <- c(0.4, 2, 5, 20)

    for (r in res) {
      res_col_name <- paste0("RNA_snn_res.", r)
      res_list[[paste0("ECODA_seuratres_", r)]] <- process_coda_fig(seurat, labels, ct_col = res_col_name, title = paste0("ECODA\nLeiden clustering ", r))
    }
  }

  if (any(Pseudobulk_PCA, ECODA_high_res_PCA, MOFA, scITD)) {
    for (i in factors_test) {
      if (Pseudobulk_PCA) {
        # Pseudobulk with PCA
        res_list[[paste0("Pseudobulk_", i, "_PCA_dims")]] <- process_pseudobulk_fig(pb_norm, labels, pca_dims = i, title = paste0("Pseudobulk gene expression\n+ PCA (", i, " dims)"))
      }

      if (ECODA_high_res_PCA) {
        # Hires CODA with PCA
        res_list[[paste0("ECODA_authors_HR_", i, "_PCA_dims")]] <- process_coda_fig(seurat, labels, pca_dims = i, ct_col = seurat@misc$hi_res_ct_col, title = paste0("ECODA\nhigh res. + PCA (", i, " dims)"))
      }

      if (MOFA) {
        # MOFA
        # Required for MOFA to run
        seurat@version <- package_version("3.1.5")
        res_list[[paste0("MOFA_hvg2000_factors", i)]][["exec_time"]] <- exec_time(
          res_list[[paste0("MOFA_hvg2000_", "factors", i)]] <- process_mofa_bulk_fig(pb_norm, metadata = metadata, labels, num_factors = i)
        )
      }

      if (scITD) {
        # scITD
        res_list[[paste0("scITD_hvg2000_factors", i)]][["exec_time"]] <- exec_time(
          res_list[[paste0("scITD_hvg2000_factors", i)]] <- process_scitd_fig(seurat, ct_col = seurat@misc$low_res_ct_col, label_col = label_col, hvg, num_factors = i)
        )
      }
    }
  }

  if (GloScope) {
    for (n_pca_dims in gloscope_n_pca_dims) {
      if (grepl("GongSharma", ds)) {
        ds_unif <- "GongSharma_all"
        gloscope_dist_file <- file.path(path_data, paste0(ds_unif, "_gloscope_hvg2000_pcadims", n_pca_dims, "_dists.rds"))
      } else {
        gloscope_dist_file <- file.path(path_data, paste0(ds, "_gloscope_hvg2000_pcadims", n_pca_dims, "_dists.rds"))
      }

      res_list[[paste0("GloScope_hvg2000_pcadims", n_pca_dims)]][["exec_time"]] <- exec_time(
        res_list[[paste0("GloScope_hvg2000_pcadims", n_pca_dims)]] <- process_gloscope_fig(seurat, metadata, label_col, gloscope_dist_file = gloscope_dist_file, n_pca_dims = n_pca_dims)
      )
      res_list[[paste0("GloScope_hvg2000_pcadims", n_pca_dims, "_sqrtmat")]] <- process_gloscope_sqrtmat_fig(metadata, label_col, gloscope_dist_file = gloscope_dist_file)
    }
  }


  for (i in c(1000, 3000)) {
    if (i > 2000) {
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
    }

    if ("var.features" %in% slotNames(seurat@assays[["RNA"]])) {
      # This path is for older Seurat objects (primarily v2/v3)
      hvg <- seurat@assays[["RNA"]]@var.features
    } else if ("var.features" %in% colnames(seurat@assays[["RNA"]]@meta.data)) {
      varf <- seurat@assays[["RNA"]]@meta.data[["var.features"]]
      hvg <- varf[!is.na(varf)]
    } else {
      stop("Could not find variable features in expected locations.")
    }

    pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
    res_list[[paste0("Pseudobulk_hvg", i)]] <- process_pseudobulk_fig(pb_norm, labels)

    pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = NULL, n_hvg = i)
    res_list[[paste0("Pseudobulk_unsup_hvg", i)]] <- process_pseudobulk_fig(pb_norm, labels)

    res_list[[paste0("MOFA_hvg", i, "_15_factors")]] <- process_mofa_bulk_fig(pb_norm, metadata = metadata, labels, num_factors = 15)

    if (GloScope) {
      if (grepl("GongSharma", ds)) {
        ds_unif <- "GongSharma_all"
        gloscope_dist_file <- file.path(path_data, paste0(ds_unif, "_gloscope_hvg", i, "_pcadims50_dists.rds"))
      } else {
        gloscope_dist_file <- file.path(path_data, paste0(ds, "_gloscope_hvg", i, "_pcadims50_dists.rds"))
      }

      res_list[[paste0("GloScope_hvg", i, "_pcadims50")]][["exec_time"]] <- exec_time(
        res_list[[paste0("GloScope_hvg", i, "_pcadims50")]] <- process_gloscope_fig(seurat, metadata, label_col, gloscope_dist_file = gloscope_dist_file)
      )
      res_list[[paste0("GloScope_hvg", i, "_pcadims50", "_sqrtmat")]] <- process_gloscope_sqrtmat_fig(metadata, label_col, gloscope_dist_file = gloscope_dist_file)
    }
  }


  for (i in c(1000, 2000, 3000)) {
    if (grepl("GongSharma", ds)) {
      ds_unif <- "GongSharma_all"
      files <- list(
        mrvi_dist_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_mrvi_dists.feather")),
        scpoli_emb_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_scpoli_embs.feather")),
        pilot_dist_file = file.path(path_data, paste0(ds_unif, "_hvg", i, "_pilot_dists.feather"))
      )
    } else {
      files <- list(
        mrvi_dist_file = file.path(path_data, paste0(ds, "_hvg", i, "_mrvi_dists.feather")),
        scpoli_emb_file = file.path(path_data, paste0(ds, "_hvg", i, "_scpoli_embs.feather")),
        pilot_dist_file = file.path(path_data, paste0(ds, "_hvg", i, "_pilot_dists.feather"))
      )
    }

    if (MrVI) {
      # MrVI
      res_list[[paste0("MrVI_hvg", i)]] <- process_mrvi_fig(files[["mrvi_dist_file"]], labels)
    }

    if (scPoli) {
      res_list[[paste0("scPoli_hvg", i)]] <- process_scpoli_fig(files[["scpoli_emb_file"]], labels)
    }

    if (PILOT) {
      res_list[[paste0("PILOT_hvg", i)]] <- process_pilot_fig(files[["pilot_dist_file"]], labels)
    }
  }


  return(res_list)
}


process_deconv_fig <- function(pseudobulk,
                               labels,
                               title = "Deconvoluted cell type composition") {
  out <- EPIC(pseudobulk, BRef)
  deconv_ct_comps <- as.data.frame(out[["mRNAProportions"]])
  row.names(deconv_ct_comps) <- colnames(pseudobulk)

  feat_mat <- clr(deconv_ct_comps + 0.001)

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

  return(res)
}


process_coda_fig <- function(seurat,
                             labels,
                             ECODA_top_n_hvct = NULL,
                             ECODA_top_varexp_hvct = NULL,
                             hvct_recalc_clr = TRUE,
                             clr = TRUE,
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

    df_freq <- df_counts %>%
      impute_zeros(clr_zero_impute_method = clr_zero_impute_method, clr_zero_impute_num = clr_zero_impute_num) %>%
      calc_perc_df()

    if (clr) {
      feat_mat <- df_freq %>%
        clr()
    } else {
      feat_mat <- df_freq
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

  if (shuffle_labels) {
    set.seed(123)
    labels <- sample(labels)
  }

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

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

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

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

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

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

  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title, coord_equal = FALSE)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

  return(res)
}

# process_mofa_fig <- function(seurat,
#                              labels,
#                              hvg,
#                              pseudobulk_hvg,
#                              save_file_path,
#                              title = "MOFA",
#                              num_factors = 5,
#                              log_abs_transform = FALSE,
#                              maxiter = 22) {
#   title <- paste0(title, " (", num_factors, " factors)")
#
#   if (!file.exists(save_file_path)) {
#     MOFAobject <- create_mofa(seurat, assays = "RNA", features = hvg)
#
#     # Default data options
#     data_opts <- get_default_data_options(MOFAobject)
#
#     # Default model options
#     model_opts <- get_default_model_options(MOFAobject)
#     model_opts$num_factors <- num_factors
#
#     # Training options
#     train_opts <- get_default_training_options(MOFAobject)
#     train_opts$convergence_mode <- "fast"
#     train_opts$seed <- 42
#     train_opts$maxiter <- maxiter
#
#     MOFAobject <- prepare_mofa(
#       object = MOFAobject,
#       data_options = data_opts,
#       model_options = model_opts,
#       training_options = train_opts
#     )
#
#     MOFAobject <- run_mofa(MOFAobject, use_basilisk = FALSE, save_data = FALSE)
#
#     saveRDS(MOFAobject, save_file_path)
#   } else {
#     MOFAobject <- readRDS(save_file_path)
#   }
#
#   feat_mat <- t(t(MOFAobject@expectations[["W"]]$RNA) %*% pseudobulk_hvg)
#
#   if (log_abs_transform) {
#     feat_mat <- log(abs(feat_mat))
#   }
#
#   res <- list()
#   res[["plot"]] <- plot_pca(feat_mat, labels, plotly_3d = FALSE, title = title, coord_equal = FALSE)
#   res[["scores"]] <- calc_sep_score(feat_mat, labels)
#   res[["feat_mat"]] <- feat_mat
#   res[["dist_mat"]] <- as.matrix(dist(feat_mat))
#
#   return(res)
# }


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

  # get donor scores-metadata associations
  pbmc_container <- get_meta_associations(pbmc_container, vars_test = c("Sample"))

  # get significant genes
  # pbmc_container <- get_lm_pvals(pbmc_container)

  feat_mat <- pbmc_container[["tucker_results"]][[1]]
  labels_scITD <- seurat@meta.data %>%
    filter(donors %in% row.names(feat_mat)) %>%
    distinct(donors, .keep_all = TRUE) %>%
    dplyr::select(donors, !!sym(label_col))
  labels <- as.factor(labels_scITD[[label_col]])

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

  return(res)
}


process_mrvi_fig <- function(mrvi_dist_file, labels, title = "MrVI") {
  feat_mat <- arrow::read_feather(mrvi_dist_file) %>%
    dplyr::select(-last_col())

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)

  return(res)
}


process_gloscope_fig <- function(seurat,
                                 metadata,
                                 label_col,
                                 gloscope_dist_file,
                                 n_pca_dims = ncol(seurat@reductions$pca@cell.embeddings),
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
  sample_names_starting_with_digit <- grepl("^\\d", samples)
  if (any(sample_names_starting_with_digit)) {
    samples[grepl("^\\d", samples)] <- paste0("g", samples[grepl("^\\d", samples)])
  }
  samples <- gsub("-", "_", samples)

  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- feat_mat

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
  sample_names_starting_with_digit <- grepl("^\\d", samples)
  if (any(sample_names_starting_with_digit)) {
    samples[grepl("^\\d", samples)] <- paste0("g", samples[grepl("^\\d", samples)])
  }
  samples <- gsub("-", "_", samples)

  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- feat_mat

  return(res)
}


process_scpoli_fig <- function(scpoli_emb_file, labels, title = "scPoli") {
  feat_mat <- arrow::read_feather(scpoli_emb_file) %>%
    dplyr::select(-last_col())

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))

  return(res)
}


process_pilot_fig <- function(pilot_dist_file, labels, title = "PILOT") {
  feat_mat <- arrow::read_feather(pilot_dist_file) %>%
    dplyr::select(-last_col())

  res <- list()
  # res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)

  return(res)
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


#' Create an example Seurat object for single-cell data simulation.
#'
#' This function generates a Seurat object with random count data and
#' simulated metadata for samples, cell types, and sample groups.
#'
#' @param n_samples Integer, the number of unique sample identifiers to create.
#' @param n_cells Integer, the total number of cells in the dataset.
#' @param n_genes Integer, the total number of genes in the dataset.
#' @param n_cell_types Integer, the number of distinct cell types to simulate.
#' @param n_groups Integer, the number of distinct groups to assign samples to.
#' @param min_cells Integer, minimum number of cells a gene must be detected in
#'   to be included in the Seurat object (passed to `CreateSeuratObject`).
#' @param min_features Integer, minimum number of features (genes) a cell must
#'   express to be included in the Seurat object (passed to `CreateSeuratObject`).
#' @param project_name Character string, the name for the Seurat project.
#' @return A Seurat object.
#' @examples
#' # Create a small Seurat object with groups
#' seurat_obj_small_groups <- create_example_seurat_object(
#'   n_samples = 4, n_cells = 200, n_genes = 50, n_cell_types = 3, n_groups = 2
#' )
#' print(seurat_obj_small_groups)
#' head(seurat_obj_small_groups@meta.data)
#' table(seurat_obj_small_groups$Sample, seurat_obj_small_groups$Group) # Check sample-group assignment
#'
#' # Create a larger Seurat object with groups
#' seurat_obj_large_groups <- create_example_seurat_object(
#'   n_samples = 10, n_cells = 1000, n_genes = 100, n_cell_types = 5, n_groups = 3,
#'   min_cells = 5, min_features = 20
#' )
#' print(seurat_obj_large_groups)
#' head(seurat_obj_large_groups@meta.data)
#' table(seurat_obj_large_groups$Sample, seurat_obj_large_groups$Group) # Check sample-group assignment
create_example_seurat_object <- function(
  n_samples = 10,
  n_cells = 1000,
  n_genes = 100,
  n_cell_types = 5,
  n_groups = 2, # <--- NEW: Number of distinct groups for samples
  min_cells = 3,
  min_features = 200,
  project_name = "ExampleProject"
) {
  # Ensure n_groups is not more than n_samples, as each sample needs a group
  if (n_groups > n_samples) {
    warning("Number of groups (n_groups) is greater than number of samples (n_samples). Setting n_groups = n_samples.")
    n_groups <- n_samples
  }

  # 1. Create a random count matrix with integers between 0 and 100
  count_matrix <- matrix(
    sample(0:100, size = n_cells * n_genes, replace = TRUE),
    nrow = n_genes, # Seurat expects genes as rows, cells as columns
    ncol = n_cells,
    dimnames = list(
      paste0("Gene_", 1:n_genes), # Gene names
      paste0("Cell_", 1:n_cells) # Cell names
    )
  )

  # Convert the count matrix to a sparse matrix (recommended for scRNA-seq)
  count_matrix_sparse <- as(count_matrix, "sparseMatrix")

  # 2. Create the 'Sample' metadata column
  sample_names <- paste0("Sample_", 1:n_samples)
  sample_assignments <- sample(sample_names, size = n_cells, replace = TRUE)

  # 3. Create cell type annotations based on n_cell_types
  cell_type_names <- paste0("Type_", LETTERS[1:n_cell_types])
  cell_type_assignments <- sample(cell_type_names, size = n_cells, replace = TRUE)

  # 4. Create Group assignments for each unique Sample
  group_names <- paste0("Group_", 1:n_groups)
  # Create a mapping from sample to group
  sample_to_group_map <- data.frame(
    Sample = sample_names,
    Group = sample(group_names, size = n_samples, replace = TRUE) # Assign each unique sample to a group
  )

  # 5. Create a data frame for the cell metadata (equivalent to AnnData's .obs)
  # Row names must be cell names to match the count matrix columns
  cell_metadata <- data.frame(
    Sample = sample_assignments,
    CellType = cell_type_assignments,
    row.names = paste0("Cell_", 1:n_cells) # Ensure row names match cell IDs in count_matrix
  )

  # 6. Merge the Group information into the cell_metadata based on Sample
  cell_metadata <- dplyr::left_join(cell_metadata, sample_to_group_map, by = "Sample")
  # Restore row names after join, as left_join often drops them
  rownames(cell_metadata) <- paste0("Cell_", 1:n_cells)


  # 7. Create the Seurat object
  seurat_object <- CreateSeuratObject(
    counts = count_matrix_sparse,
    project = project_name,
    min.cells = min_cells,
    min.features = min_features,
    meta.data = cell_metadata
  )

  # Optional: Print some summary information
  message(paste0(
    "Created Seurat object with ", ncol(seurat_object), " cells and ",
    nrow(seurat_object), " genes."
  ))
  message(paste0("Number of unique samples: ", length(unique(seurat_object$Sample))))
  message(paste0("Number of unique cell types: ", length(unique(seurat_object$CellType))))
  message(paste0("Number of unique groups: ", length(unique(seurat_object$Group))))

  # Return the Seurat object
  return(seurat_object)
}


exec_time <- function(fun) {
  start_time <- Sys.time()
  fun
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  time_taken
  return(time_taken)
}
