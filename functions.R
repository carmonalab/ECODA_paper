library(foreach)
library(doParallel)


# Calculate coefficient of variation
cv <- function(x) {
  sd_x <- sd(x)  # Standard deviation
  mean_x <- mean(x)  # Mean
  cv_value <- abs(sd_x / mean_x) * 100  # Coefficient of variation (%)
  return(cv_value)
}


# Calculate average silhouette width
calc_sil <- function(feat_mat,
                     labels) {
  sils <- cluster::silhouette(x = as.numeric(factor(labels)),
                              dist = dist(feat_mat)) %>%
    as.data.frame()
  score <- mean(sils[["sil_width"]])
  return(score)
}



compute_KNN <- function(dist, knn_k){
  # Compute KNN
  knn <- RANN::nn2(as.matrix(dist), k = knn_k + 1)$nn.idx
  knn <- knn[, -1]  # Remove self-neighbor
  return(knn)
}
# moudularity
compute_snn_graph <- function(dist, knn_k) {
  
  knn <- compute_KNN(dist = dist, knn_k = knn_k)
  # Initialize adjacency matrix
  n <- nrow(as.matrix(dist))
  adj_matrix <- matrix(0, n, n)
  # Count shared neighbors
  for (i in seq_len(n)) {
    for (j in knn[i, ]) {
      shared_neighbors <- length(intersect(knn[i, ], knn[j, ]))
      adj_matrix[i, j] <- shared_neighbors
      adj_matrix[j, i] <- shared_neighbors  # Ensure symmetry
    }
  }
  # Create graph object
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)
  return(g)
}
calc_modularity <- function(feat_mat,
                            labels,
                            knn_k = NULL) {
  if (is.null(knn_k)) {
    knn_k <- round(nrow(feat_mat) / 4)
  }
  
  # Create a graph object
  g <- compute_snn_graph(dist = feat_mat, knn_k = knn_k)
  # Compute modularity
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))
  return(modularity_score)
}


# # Calculate modularity score
# calc_modularity <- function(feat_mat,
#                             labels,
#                             d = 5) {
#   graph <- scran::buildKNNGraph(as.matrix(feat_mat),
#                                 d = d,
#                                 transposed = TRUE,
#                                 k = 5)
#   score <- igraph::modularity(graph, membership = as.numeric(factor(labels)))
#   return(score)
# }




# Test data transformation methods for ECODA
datrans <- function(count_mat,
                    labels = NULL,
                    Percent_difference, # Percent cell abundance difference (e.g. 100 equals one cell type being twice as abundant)
                    n_ct_to_select, # Number of randomly selected cell types to be differentially abundant
                    cts = NULL, # Cell types to be differentially abundant. If NULL, randomly select a specified number of cell types (n_ct_to_select)
                    reps = 20, # Number of random shuffling to calculate separation using different cell types and samples for DA
                    methods = c(
                      "counts",
                      # "counts_imputed",
                      # "counts_pca",
                      "freq",
                      # "freq_imputed",
                      # "freq_pca",
                      "asin_sqrt",
                      "alr_randref",
                      # "alr_randref_pca",
                      "alr_mincvref",
                      # "alr_mincvref_pca",
                      # "ilr", "ilr_pca",
                      "clr",
                      "clr_pca"
                    ),
                    n_cores = 8
) {
  
  labs <- labels
  
  colnames(count_mat) <- make.names(colnames(count_mat), unique = TRUE)
  
  n_half_samples <- round(dim(count_mat)[1] / 2)
  
  if (!is.null(cts)) {
    n_ct_to_select[n_ct_to_select >= length(cts)] <- length(cts)
  }
  
  # Register cluster
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  rets <- foreach (pd = Percent_difference,
                   .export = c(
                     "calc_perc_df",
                     "impute_zeros",
                     "calc_sil",
                     "calc_modularity", "compute_snn_graph", "compute_KNN",
                     "clust_eval", "adjustedRandIndex",
                     "cv",
                     "clr"),
                   .packages = c("dplyr"),
                   .combine = rbind) %dopar% {
                     
                     res <- data.frame(
                       method = character(),
                       n_celltypes = numeric(),
                       Percent_difference = numeric(),
                       Silhouette_score = numeric(),
                       Modularity_score = numeric(),
                       Adjusted_Rand_Index = numeric()
                     )
                     
                     print(paste0("pd: ", pd))
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
                           labs <- as.numeric(row.names(count_mat) %in% half_samples_da)
                           
                           # Simulate differential abundance
                           rsums_before <- rowSums(df_counts_temp)
                           df_counts_temp[half_samples_da, ct_da] <- round(df_counts_temp[half_samples_da, ct_da] * (1 + pd/100))
                           rsums_after <- rowSums(df_counts_temp)
                           df_counts_temp <- round(df_counts_temp / (rsums_after / rsums_before))
                         }
                         
                         
                         df_freq <- df_counts_temp %>% calc_perc_df()
                         df_asin_sqrt <- asin(sqrt(df_freq/100))
                         df_counts_temp_imputed <- df_counts_temp %>% impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.1)
                         df_freq_imputed <- df_counts_temp_imputed %>% calc_perc_df()
                         
                         for (met in methods) {
                           if (met == "counts") {df <- df_counts_temp}
                           else if (met == "counts_imputed") {df <- df_counts_temp_imputed}
                           else if (met == "freq") {df <- df_freq}
                           else if (met == "freq_imputed") {df <- df_freq_imputed}
                           else if (met == "asin_sqrt") {df <- df_asin_sqrt}
                           else if (met == "asin_sqrt_imputed") {df <- df_asin_sqrt}
                           else if (met == "alr_mincvref") {
                             ct_ref <- sample(colnames(df_freq_imputed)[!colnames(df_freq_imputed) %in% ct_da], size = 1)
                             df <- Hotelling::alr(as.formula(paste0(ct_ref, "~.")), df_freq_imputed)
                           }
                           else if (met == "alr_randref") {
                             cvs <- apply(df_freq_imputed, 2, cv)
                             min_cv <- min(cvs[!colnames(df_freq_imputed) %in% ct_da])
                             ct_ref_mincv <- colnames(df_freq_imputed)[which(cvs == min_cv)][1]
                             df <- Hotelling::alr(as.formula(paste0(ct_ref_mincv, "~.")), df_freq_imputed)
                           }
                           else if (met == "ilr") {df <- compositions::ilr(df_freq)}
                           else if (met == "clr") {df <- clr(df_freq_imputed)}
                           else if (met == "clr_centered") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = FALSE)}
                           else if (met == "clr_centered_scaled") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = TRUE)}
                           
                           if (grepl("pca", met)) {df <- prcomp(df)$x[, 1:3]}
                           
                           
                           # Calculate scores
                           avg_sil <- calc_sil(feat_mat = df, labels = labs)
                           mod <- calc_modularity(feat_mat = df, labels = labs)
                           cluster_score <- clust_eval(matrix = df, labels = labs)
                           
                           
                           # Append results
                           new_row <- list(
                             method = met,
                             n_celltypes = nct,
                             Percent_difference = pd,
                             Silhouette_score = avg_sil,
                             Modularity_score = mod,
                             Adjusted_Rand_Index = cluster_score
                           )
                           res <- rbind(res, new_row)
                         }
                       }
                     }
                     
                     return(res)
                   }
  #stop cluster
  stopCluster(cluster)
  
  return(rets)
}


# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(feat_mat,
                     labels,
                     scale. = FALSE,
                     knn_k = NULL,
                     title = NULL,
                     sil_score = TRUE,
                     mod_score = TRUE,
                     cluster_score = TRUE,
                     pointsize = 3,
                     coord_equal = TRUE,
                     plotly_3d = FALSE) {
  
  res.pca <- prcomp(feat_mat, scale. = scale.)
  
  if (cluster_score) {
    cluster_score <- clust_eval(feat_mat, labels)
    title <- paste0(title, "\nCluster score: ", cluster_score)
  }
  if (mod_score) {
    mod_score <- round(calc_modularity(feat_mat, labels, knn_k), 3)
    title <- paste0(title, "\nModularity score: ", mod_score)
  }
  if (sil_score) {
    sil_score <- round(calc_sil(feat_mat, labels), 3)
    title <- paste0(title, "\nSilhouette score: ", sil_score)
  }
  
  p <- factoextra::fviz_pca(res.pca,
                            habillage = labels,
                            label = "var",
                            pointsize = pointsize,
                            invisible = c("var", "quali"),
                            geom = "point") +
    ggtitle(title) +
    scale_shape_manual(values = rep(19, length(unique(labels))))
  if (coord_equal) {
    p <- p + coord_equal()
  }
  print(p)
  
  if (plotly_3d) {
    df <- as.data.frame(res.pca$x)
    df$id <- seq_len(nrow(df))
    df$vs <- factor(labels)
    ms <- replicate(2, df, simplify = F)
    ms[[2]]$PC3 <- min(df$PC3)
    m <- ms %>%
      bind_rows() %>%
      group2NA("id", "vs")
    # Plotting with plotly
    plotly <- plot_ly(color = ~vs) %>%
      add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
      add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity=0.2)
    print(plotly)
  }
  
  return(p)
}


# Plot rocauc curve separating samples by PCA loading 1
plot_roc_pcad1 <- function(feat_mat, labels) {
  
  res.pca <- prcomp(feat_mat)
  
  pca_scores_pc1 <- data.frame(PC1 = data.frame(res.pca$x)$PC1,
                               label = factor(labels))
  p <- ggplot(pca_scores_pc1, aes(x = PC1, fill = label)) +
    geom_density(alpha=0.4) + 
    geom_histogram(aes(y=..density..), alpha=0.3, colour="black", 
                   position="identity") +
    labs(title = "Density Plot of PC1 Scores",
         x = "PC1 Scores",
         y = "Density") +
    theme_classic()
  print(p)
  
  # Calculate the ROC curve
  roc_curve <- roc(pca_scores_pc1$label, ifelse(pca_scores_pc1$PC1 > 0, 1, 0))
  
  # Plot the ROC curve
  plot(roc_curve, col = "blue", main = "ROC Curve", print.auc = TRUE)
}


# Cluster samples and compare to original annotation
clust_eval <- function(matrix,
                       labels,
                       nclusts = 2,
                       digits = 3,
                       return_mean = TRUE) {
  results <- list()
  dist_mat <- dist(matrix)
  
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


DESeq2.normalize <- function(matrix,
                             metadata,
                             nvar_genes = 2000) {
  suppressMessages({
    suppressWarnings({
      
      # Normalize pseudobulk data using DESeq2
      # do formula for design with the cluster_by elements in order
      matrix <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                               colData = metadata,
                                               design = stats::formula(paste("~ 1")))
      
      matrix <- DESeq2::estimateSizeFactors(matrix)
      
      # Set minimum number of counts per gene
      nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(matrix, normalized=TRUE)) > 10 ))
      
      # transform counts using vst
      matrix <- DESeq2::vst(matrix, blind = T, nsub = nsub)
      matrix <- SummarizedExperiment::assay(matrix)
      
      # get top variable genes
      rv <- MatrixGenerics::rowVars(matrix)
      select <- order(rv, decreasing=TRUE)[seq_len(min(nvar_genes, length(rv)))]
      select <- row.names(matrix)[select]
      
      matrix <- matrix[select[select %in% row.names(matrix)],]
    })
  })
  
  return(matrix)
}


load_seurat_object <- function(file_name) {
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

get_metadata <- function(seurat) {
  metadata <- seurat@meta.data %>%
    group_by(Sample) %>%
    slice(1)

  return(metadata)
}

get_labels <- function(seurat, label_col) {
  metadata <- get_metadata(seurat)
  labels <- metadata[[label_col]]
  
  return(labels)
}


process_pseudobulk_fig <- function(pb_norm,
                                   labels,
                                   pca_dims = NULL,
                                   title = "Pseudobulk gene expression",
                                   knn_k = NULL) {
  metadata <- get_metadata(seurat)
  
  if (!is.null(pca_dims)) {
    res.pca <- prcomp(pb_norm, center = TRUE, scale. = FALSE)
    feat_mat <- res.pca[["x"]][, 1:pca_dims]
  } else {
    feat_mat <- pb_norm
  }
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, plotly_3d = FALSE, knn_k = knn_k, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_deconv_fig <- function(pseudobulk,
                               labels,
                               title = "Deconvoluted cell type composition") {
  out <- EPIC(pseudobulk, BRef)
  deconv_ct_comps <- as.data.frame(out[["mRNAProportions"]])
  row.names(deconv_ct_comps) <- colnames(pseudobulk)
  
  feat_mat <- clr(deconv_ct_comps + 0.001)
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_coda_fig <- function(seurat,
                             labels,
                             pca_dims = NULL,
                             ct_col,
                             title,
                             knn_k = NULL) {
  
  df_counts <- get_ct_comp_df_seurat(seurat, sample_col = "Sample", ct_col = ct_col)
  
  feat_mat <- df_counts %>%
    impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.1) %>%
    calc_perc_df() %>%
    clr()
  
  if (!is.null(pca_dims)) {
    res.pca <- prcomp(feat_mat, center = TRUE, scale. = FALSE)
    feat_mat <- res.pca[["x"]][, 1:pca_dims]
  }
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, plotly_3d = FALSE, knn_k = knn_k, title = title)
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
  
  feat_mat <- as.data.frame(MOFAobject@expectations[["Z"]])

  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, plotly_3d = FALSE, title = title, coord_equal = FALSE)
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
  seurat$ctypes <- as.character(seurat[[ct_col]])
  
  # if (!file.exists(file.path(path_output, "pbmc_container_cmv.rds"))) {
  ctypes <- unique(seurat$ctypes)
  ctypes <- ctypes[!is.na(ctypes)]
  
  # Step 1: Calculate table of donors by ctypes
  cell_counts <- table(seurat$donors, seurat$ctypes)
  
  # Step 2: Determine low abundant cell types (less than 5 cells in more than 20% of donors)
  ctypes_to_drop <- names(which(colMeans(cell_counts < 5) > 0.2))
  
  ctypes <- ctypes[!ctypes %in% ctypes_to_drop]
  
  # set up project parameters
  param_list <- initialize_params(ctypes_use = as.character(ctypes),
                                  ncores = 6, rand_seed = 10)
  
  # create project container
  pbmc_container <- make_new_container(count_data = seurat@assays$RNA$counts,
                                       meta_data = seurat@meta.data[, c("donors", "ctypes", label_col, "Sample")],
                                       gn_convert = NULL,
                                       params = param_list)
  
  # form the tensor from the data
  pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5,
                                norm_method='trim', scale_factor=10000,
                                custom_genes=hvg,
                                scale_var = TRUE, var_scale_power = 2)
  
  # run the tensor decomposition
  pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(num_factors, num_factors + 5),
                                   tucker_type = 'regular', rotation_type = 'hybrid')
  
  # get donor scores-metadata associations
  pbmc_container <- get_meta_associations(pbmc_container, vars_test=c("Sample"), stat_use='pval')
  
  # get significant genes
  # pbmc_container <- get_lm_pvals(pbmc_container)
  
  feat_mat <- pbmc_container[["tucker_results"]][[1]]
  labels_scITD <- seurat@meta.data %>%
    filter(donors %in% row.names(feat_mat)) %>%
    distinct(donors, .keep_all = TRUE) %>%
    select(donors, !!sym(label_col))
  labels <- as.numeric(as.factor(labels_scITD[[label_col]]))
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_mrvi_fig <- function(seurat,
                             dist_file,
                             ct_col,
                             label_col,
                             title = "MrVI") {
  dist_list <- list()
  
  # Open the NetCDF file
  dists <- nc_open(dist_file)
  
  # List the available layers (cell types)
  layers <- ncvar_get(dists, paste0(ct_col, "_name"))
  
  samples <- ncvar_get(dists, "sample_x")
  metadata <- seurat@meta.data %>%
    group_by(Sample) %>%
    slice(1)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  
  for (l in layers) {
    # Find the index corresponding to "CD4T"
    index <- which(layers == l)
    
    # If you want to get the distance matrix for "CD4T":
    dist_mat <- ncvar_get(dists, ct_col, start = c(1, 1, index), count = c(-1, -1, 1))
    
    dist_list[[l]] <- dist_mat
  }
  
  # Close the NetCDF file handle when done
  nc_close(dists)
  
  feat_mat <- Reduce(`+`, dist_list) / length(dist_list)

  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)
  
  return(res)
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


calc_perc_df <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>% as.data.frame()
  return(df)
}


impute_zeros <- function(df,
                         clr_zero_impute_method = c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all"),
                         clr_zero_impute_num = 1) {
  if (!clr_zero_impute_method %in% c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all")) {
    stop("clr_zero_impute_method not found")
  }
  
  # Apply specified zero imputation method
  if (clr_zero_impute_method == "percentage_zeros") {
    for(row in 1:nrow(df)){
      df[row,][df[row,] == 0] <- sum(df[row,])/100 * clr_zero_impute_num
    }
  }
  else if (clr_zero_impute_method == "percentage_all") {
    for(row in 1:nrow(df)){
      df[row,] <- df[row,] + sum(df[row,])/100 * clr_zero_impute_num
    }
  }
  else if (clr_zero_impute_method == "counts_zeros") {
    # Impute zeros by replacing them with a small non-zero value (1 in this case)
    df[df == 0] <- clr_zero_impute_num
  }
  else if (clr_zero_impute_method == "counts_all") {
    # Impute by adding a fixed count to all values
    df <- df + clr_zero_impute_num
  }
  
  return(df)
}


clr <- function(df) {
  percentage_df <- calc_perc_df(df)
  
  geometric_mean <- apply(percentage_df, 1, function(row) exp(mean(log(row))))
  clr_df <- apply(percentage_df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()
  
  return(clr_df)
}



calc_sep_score <- function(df,
                           labels,
                           knn_k = NULL) {
  sil_score <- round(calc_sil(df, labels), 3)
  mod_score <- unlist(round(calc_modularity(df, labels, knn_k), 3))
  cluster_score <- clust_eval(matrix = df, labels = labels)
  anosim_score <- vegan::anosim(x = df, grouping = labels, distance = "euclidean")[["statistic"]]
  
  res <- list(sil_score = sil_score,
              mod_score = mod_score,
              cluster_score = cluster_score)
  
  return(res)
}