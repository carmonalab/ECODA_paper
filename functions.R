# Calculate coefficient of variation
cv <- function(x) {
  sd_x <- sd(x)  # Standard deviation
  mean_x <- mean(x)  # Mean
  cv_value <- sd_x / mean_x * 100  # Coefficient of variation (%)
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



compute_KNN <- function(dist, KNNGraph_k){
  # Compute KNN
  knn <- RANN::nn2(as.matrix(dist), k = KNNGraph_k + 1)$nn.idx
  knn <- knn[, -1]  # Remove self-neighbor
  return(knn)
}
# moudularity
compute_snn_graph <- function(dist, KNNGraph_k = 5, knn = NULL) {
  if(is.null(knn)){
    knn <- compute_KNN(dist = dist, KNNGraph_k = KNNGraph_k)
  }
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
                            KNNGraph_k = 5,
                            knn = NULL) {
  # Create a graph object
  g <- compute_snn_graph(dist = feat_mat,
                         KNNGraph_k = KNNGraph_k,
                         knn = knn)
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
                    perc_diffs, # Percent cell abundance difference (e.g. 100 equals one cell type being twice as abundant)
                    n_ct_to_select, # Number of randomly selected cell types to be differentially abundant
                    cts = NULL, # Cell types to be differentially abundant. If NULL, randomly select a specified number of cell types (n_ct_to_select)
                    reps = 20, # Number of random shuffling to calculate separation using different cell types and samples for DA
                    methods = c(
                      "counts",
                      # "counts_pca",
                      "freq",
                      # "freq_pca",
                      "asin_sqrt",
                      "alr_randref",
                      # "alr_randref_pca",
                      "alr_mincvref",
                      # "alr_mincvref_pca",
                      # "ilr", "ilr_pca",
                      "clr",
                      "clr_pca"
                    )
) {
  
  n_half_samples <- round(dim(df_counts)[1] / 2)
  
  results <- data.frame(
    method = character(),
    n_celltypes = numeric(),
    perc_diff = numeric(),
    silhouette = numeric(),
    modularity = numeric(),
    clustering = numeric()
  )
  
  if (!is.null(cts)) {
    n_ct_to_select[n_ct_to_select >= length(cts)] <- length(cts)
  }
  
  for (pd in perc_diffs) {
    print(paste0("pd: ", pd))
    for (nct in n_ct_to_select) {
      
      print(paste0("nct: ", nct))
      for (rep in 1:reps) {
        
        # Prepare data
        df_counts_temp <- df_counts
        
        if (!is.null(cts)) {
          ct_da <- sample(cts, size = nct)
        } else {
          ct_da <- sample(colnames(df_counts_temp), size = nct)
        }
        
        half_samples_da <- sample(row.names(df_counts_temp), size = n_half_samples)
        labels <- as.numeric(row.names(df_counts) %in% half_samples_da)
        
        # Simulate differential abundance
        rsums_before <- rowSums(df_counts_temp)
        df_counts_temp[half_samples_da, ct_da] <- round(df_counts_temp[half_samples_da, ct_da] * (1 + pd/100))
        rsums_after <- rowSums(df_counts_temp)
        df_counts_temp <- round(df_counts_temp / (rsums_after / rsums_before)) + 1
        df_freq <- na.omit(df_counts_temp/rowSums(df_counts_temp))
        df_asin_sqrt <- asin(sqrt(df_freq))
        
        for (met in methods) {
          if (grepl("counts", met)) {
            df <- df_counts_temp
          } else {
             if (grepl("freq", met)) {
              df <- df_freq
            } else {
              if (grepl("asin_sqrt", met)) {
                df <- df_asin_sqrt
              } else if (grepl("alr_randref", met)) {
                ct_ref <- sample(colnames(df_freq)[!colnames(df_freq) %in% ct_da], size = 1)
                df <- Hotelling::alr(as.formula(paste0(ct_ref, "~.")), df_freq)
              } else if (grepl("alr_mincvref", met)) {
                cvs <- apply(df_freq, 2, cv)
                min_cv <- min(cvs[!colnames(df_freq) %in% ct_da])
                ct_ref_mincv <- colnames(df_freq)[which(cvs == min_cv)]
                df <- Hotelling::alr(as.formula(paste0(ct_ref_mincv, "~.")), df_freq)
                # } else if (grepl("ilr", met)) {
                #   df <- compositions::ilr(df_freq)
              } else if (grepl("clr", met)) {
                df <- Hotelling::clr(df_freq)
              }
            }
            
          }
          
          if (grepl("pca", met)) {
            df <- prcomp(df)$x[, 1:3]
          }
          
          
          # Calculate scores
          avg_sil <- calc_sil(feat_mat = df, labels = labels)
          mod <- calc_modularity(feat_mat = df, labels = labels)
          cluster_score <- clust_eval(matrix = df, labels = labels)
          
          
          # Append results
          new_row <- list(
            method = met,
            n_celltypes = nct,
            perc_diff = pd,
            silhouette = avg_sil,
            modularity = mod,
            clustering = cluster_score
          )
          results <- rbind(results, new_row)
        }
      }
    }
  }
  
  return(results)
}


# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(feat_mat,
                     labels,
                     knn_k = 5,
                     title = NULL,
                     coord_equal = TRUE,
                     plotly_3d = TRUE) {
  
  res.pca <- prcomp(feat_mat)
  
  sil_score <- round(calc_sil(feat_mat, labels), 3)
  mod_score <- round(calc_modularity(feat_mat, labels, knn_k), 3)
  
  p <- factoextra::fviz_pca(res.pca,
                            habillage = labels,
                            label = "var",
                            pointsize = 3,
                            invisible = c("var", "quali"),
                            geom = "point") +
    ggtitle(paste(title, "\nSilhouette score: ", sil_score, "\nModularity score: ", mod_score)) +
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
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  clust_labels <- stats::cutree(hc, k = nclusts)
  acc <- sum(as.numeric(as.factor(labels)) == clust_labels) / length(labels)
  results[["hclust_accuracy"]] <- max(acc, 1-acc)
  # Perform PAM clustering
  pam_result <- cluster::pam(matrix, k = nclusts)
  # View cluster memberships
  clust_labels <- pam_result$cluster
  acc <- sum(as.numeric(as.factor(labels)) == clust_labels) / length(labels)
  results[["pamclust_accuracy"]] <- max(acc, 1-acc)
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
                             nvar_genes = 1500) {
  suppressMessages({
    suppressWarnings({
      
      # Normalize pseudobulk data using DESeq2
      # do formula for design with the cluster_by elements in order
      matrix <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                               colData = metadata,
                                               design = ~1)
      
      matrix <- DESeq2::estimateSizeFactors(matrix)
      
      # transform counts using vst
      matrix <- DESeq2::vst(matrix, blind = T)
      matrix <- SummarizedExperiment::assay(matrix)
      
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
  
  seurat <- LoadH5Seurat(file_name_seurat, meta.data = F)
  meta <- arrow::read_feather(paste0(file_name_short, ".feather"))
  seurat <- AddMetaData(seurat, meta)
  
  return(seurat)
}


get_labels <- function(seurat, label_col) {
  metadata <- seurat@meta.data %>%
    group_by(Sample) %>%
    slice(1)
  labels <- metadata[[label_col]]
  
  return(labels)
}


process_pseudobulk_fig <- function(pseudobulk_hvg,
                                   labels,
                                   hvg,
                                   title = "Pseudobulk gene expression",
                                   knn_k = 5) {
  metadata <- seurat@meta.data %>%
    group_by(Sample) %>%
    slice(1)
  
  pseudobulk_normalized <- t(DESeq2.normalize(pseudobulk_hvg, metadata))
  
  fig <- plot_pca(pseudobulk_normalized, labels, plotly_3d = FALSE, knn_k = knn_k, title = title)
  
  return(fig)
}


process_deconv_fig <- function(pseudobulk,
                               labels,
                               title = "Pseudobulk - deconvoluted cell type composition") {
  out <- EPIC(pseudobulk, BRef)
  deconv_ct_comps <- as.data.frame(out[["mRNAProportions"]])
  row.names(deconv_ct_comps) <- colnames(pseudobulk)
  fig <- plot_pca(Hotelling::clr(deconv_ct_comps + 0.01), labels, title = title)
  
  return(fig)
}


process_coda_fig <- function(seurat,
                             labels,
                             ct_col,
                             title,
                             knn_k = 5) {
  feat_mat <- table(seurat$Sample, seurat@meta.data[[ct_col]]) + 1
  feat_mat <- as.data.frame.matrix(t(prop.table(feat_mat, margin = 1)))
  feat_mat <- feat_mat %>%
    t() %>%
    Hotelling::clr() %>%
    t() %>%
    scale(center = TRUE, scale = FALSE) %>%
    t()
  fig <- plot_pca(feat_mat, labels, plotly_3d = FALSE, knn_k = knn_k, title = title)
}


process_mofa_fig <- function(seurat,
                             labels,
                             hvg,
                             pseudobulk_hvg,
                             save_file_path,
                             title = "MOFA2",
                             num_factors = 5,
                             maxiter = 22) {
  if (!file.exists(save_file_path)) {
    MOFAobject <- create_mofa(seurat, assays = "RNA", features = hvg)
    
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
    
    MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE, save_data = FALSE)
    
    saveRDS(MOFAobject, save_file_path)
  } else {
    MOFAobject <- readRDS(save_file_path)
  }
  
  feat_mat <- t(t(MOFAobject@expectations[["W"]]$RNA) %*% pseudobulk_hvg)
  fig <- plot_pca(feat_mat, labels, plotly_3d = FALSE, title = title, coord_equal = FALSE)
  
  return(fig)
}


process_scitd_fig <- function(seurat,
                              ctype_col,
                              label_col,
                              hvg,
                              title = "scITD") {
  
  seurat$donors <- seurat$Sample
  seurat$ctypes <- as.character(seurat[[ctype_col]])
  
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
                                       meta_data = seurat@meta.data[, c("donors", "ctypes", "subject.cmv", "Sample")],
                                       gn_convert = NULL,
                                       params = param_list)
  
  # form the tensor from the data
  pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5,
                                norm_method='trim', scale_factor=10000,
                                custom_genes=hvg,
                                scale_var = TRUE, var_scale_power = 2)
  
  # run the tensor decomposition
  pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10),
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
  labels <- as.numeric(labels_scITD[[label_col]])
  
  fig <- plot_pca(feat_mat, labels, title = title)
  
  return(fig)
}


process_mrvi_fig <- function(seurat,
                             dist_file,
                             ctype_col,
                             label_col,
                             title = "MrVI") {
  dist_list <- list()
  
  # Open the NetCDF file
  dists <- nc_open(dist_file)
  
  # List the available layers (cell types)
  layers <- ncvar_get(dists, paste0(ctype_col, "_name"))
  
  samples <- ncvar_get(dists, "sample_x")
  metadata <- seurat@meta.data %>%
    group_by(Sample) %>%
    slice(1)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  
  for (l in layers) {
    # Find the index corresponding to "CD4T"
    index <- which(layers == l)
    
    # If you want to get the distance matrix for "CD4T":
    dist_mat <- ncvar_get(dists, ctype_col, start = c(1, 1, index), count = c(-1, -1, 1))
    
    dist_list[[l]] <- dist_mat
  }
  
  # Close the NetCDF file handle when done
  nc_close(dists)
  
  dist_sum <- Reduce(`+`, dist_list) / length(dist_list)
  fig <- plot_pca(dist_sum, labels, title = "MrVI")
  
  return(fig)
}


get_ct_comp_df_seurat <- function(seurat, sample_col, ct_col) {
  ct_comp_df <- table(seurat@meta.data[[sample_col]], seurat@meta.data[[ct_col]]) %>%
    t() %>%
    as.data.frame.matrix() %>%
    t() %>%
    as.data.frame()

  return(ct_comp_df)
}


calc_perc_df <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>% as.data.frame()
  return(df)
}


clr <- function(df,
                clr_zero_impute_method = c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all"),
                clr_zero_impute_num = 1) {
  if (!clr_zero_impute_method %in% c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all")) {
    stop("clr_zero_impute_method not found")
  }
  
  # Check if clr_zero_impute_method is valid
  if (!clr_zero_impute_method %in% c("percentage_zeros", "percentage_all", "counts_zeros", "counts_all")) {
    stop("clr_zero_impute_method not found")
  }
  
  # Apply specified zero imputation method
  if (clr_zero_impute_method == "percentage_zeros") {
    for(row in 1:nrow(df)){
      df[row,][df[row,] == 0] <- sum(df[row,])/100
    }
  }
  if (clr_zero_impute_method == "percentage_all") {
    for(row in 1:nrow(df)){
      df[row,] <- sum(df[row,])/100
    }
  }
  if (clr_zero_impute_method == "counts_zeros") {
    # Impute zeros by replacing them with a small non-zero value (1 in this case)
    df[df == 0] <- clr_zero_impute_num
  }
  if (clr_zero_impute_method == "counts_all") {
    # Impute by adding a fixed count to all values
    df <- df + clr_zero_impute_num
  }
  
  percentage_df <- calc_perc_df(df)
  
  geometric_mean <- apply(percentage_df, 1, function(row) exp(mean(log(row))))
  clr_transformed <- apply(percentage_df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()

  return(clr_df)
}



calc_sep_score <- function(df,
                           labels,
                           ct_col,
                           knn_k = 5) {
  sil_score <- round(calc_sil(df, labels), 3)
  mod_score <- unlist(round(calc_modularity(df, labels, knn_k), 3))
  
  res <- list(sil_score = sil_score,
              mod_score = mod_score)
  
  return(res)
}
