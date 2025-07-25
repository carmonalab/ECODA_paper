
suppressPackageStartupMessages({
  # install.packages("devtools")
  library(compositions)
  library(doParallel)
  library(factoextra)
  library(foreach)
  library(ggplot2)
  library(ggpubr)
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
  
  # devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
  library(EPIC)
  # BiocManager::install("MOFA2")
  library(reticulate)
  py_install("mofapy2", envname = "my_env", method="auto")
  import("mofapy2")
  library(basilisk)
  library("MOFA2")
  # BiocManager::install("fgsea")
  # BiocManager::install("sva")
  # devtools::install_github("kharchenkolab/scITD")
  library(scITD)
  # install.packages("ncdf4")
  # install.packages("GloScope")
  # library(GloScope)
  library(parallelly)
  
  # remotes::install_github("carmonalab/scooter", ref="f31eab3")
  library(scooter)
  
  # Import at the end, otherwise the "select" function gets re-assigned by another package
  library(dplyr)
  
  # library(arrow)
})


clr <- function(df) {
  percentage_df <- calc_perc_df(df)
  
  geometric_mean <- apply(percentage_df, 1, function(row) exp(mean(log(row))))
  clr_df <- apply(percentage_df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()
  
  return(clr_df)
}


# Calculate coefficient of variation
cv <- function(x) {
  sd_x <- sd(x)  # Standard deviation
  mean_x <- mean(x)  # Mean
  cv_value <- abs(sd_x / mean_x) * 100  # Coefficient of variation (%)
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
  sils <- cluster::silhouette(x = as.numeric(factor(labels)),
                              dist = dist(feat_mat)) %>%
    as.data.frame()
  score <- mean(sils[["sil_width"]])
  return(score)
}


#----------------------------------------------------------->
calc_modularity <- function(feat_mat,
                            labels,
                            knn_k = NULL) {
  if (is.null(knn_k)) {
    # min_group_size <- min(table(labels))
    ngroups <- length(unique(labels))
    half_group_size <- round(nrow(feat_mat) / ngroups / 2)
    # # knn_k is equal to half average group size or at least 3
    knn_k <- max(half_group_size, 3)
  }
  
  # Create a graph object
  g <- compute_snn_graph(feat_mat = feat_mat, knn_k = knn_k)
  # Compute modularity
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))
  
  # NOTE:
  # Maximum modularity depends on the number of groups: max(mod) = 1 - 1 / (number of groups)
  # see Brandes, Ulrik, et al. On finding graph clusterings with maximum modularity.
  
  return(modularity_score)
}


compute_KNN <- function(feat_mat, knn_k){
  # Compute KNN
  knn <- RANN::nn2(as.matrix(feat_mat), k = knn_k + 1)$nn.idx
  knn <- knn[, -1]  # Remove self-neighbor
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
      adj_matrix[j, i] <- shared_neighbors  # Ensure symmetry
    }
  }
  # adj_matrix <- spdep::knearneigh(feat_mat, k=5, longlat = NULL)
  # Create graph object
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)
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
  file_path_h5seurat <- file.path(temp_file_path, paste0(ds, ".h5Seurat"))
  file_path_h5ad <- file.path(shared_storage_path, paste0(ds, ".h5ad"))
  if (!file.exists(file_path_h5ad)) {
    SaveH5Seurat(seurat, filename = file_path_h5seurat)
    Convert(file_path_h5seurat, dest = "h5ad")
    unlink(file_path_h5seurat)
    print(paste0("File saved to: ", file_path_h5ad))
  } else {
    print("File already processed")
  }
}


# Test data transformation methods for ECODA
datrans <- function(count_mat,
                    labels = NULL,
                    Amount_of_perturbation, # Percent cell abundance difference (e.g. 100 equals one cell type being twice as abundant)
                    n_ct_to_select, # Number of randomly selected cell types to be differentially abundant
                    cts = NULL, # Cell types to be differentially abundant. If NULL, randomly select a specified number of cell types (n_ct_to_select)
                    reps = 20, # Number of random shuffling to calculate separation using different cell types and samples for DA
                    trans_method = c(
                      # "counts",
                      # "counts_imputed",
                      "counts_pca",
                      # "freq",
                      # "freq_imputed",
                      "freq_pca",
                      # "arcsine_sqrt",
                      "arcsine_sqrt_pca",
                      # "alr_randref",
                      "alr_randref_pca",
                      # "alr_mincvref",
                      "alr_mincvref_pca",
                      # "ilr", "ilr_pca",
                      # "clr",
                      "clr_pca"
                    ),
                    zero_imp_method = "percentage_all__0.1",
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
  
  rets <- foreach (da = Amount_of_perturbation,
                   .export = c(
                     "calc_perc_df",
                     "impute_zeros",
                     "calc_sil",
                     "calc_modularity", "compute_snn_graph", "compute_KNN",
                     "clust_eval", "adjustedRandIndex",
                     "cv",
                     "clr"),
                   .packages = c("dplyr", "robCompositions", "zCompositions"),
                   .errorhandling="pass",
                   .combine = rbind) %dopar% {
                     
                     res <- data.frame(
                       trans_method = character(),
                       zero_imp_method = character(),
                       n_celltypes = numeric(),
                       Amount_of_perturbation = numeric(),
                       Silhouette_score = numeric(),
                       Modularity_score = numeric(),
                       Adjusted_Rand_Index = numeric()
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
                           labs <- as.numeric(row.names(count_mat) %in% half_samples_da)
                           
                           # Simulate differential abundance
                           rsums_before <- rowSums(df_counts_temp)
                           df_counts_temp[half_samples_da, ct_da] <- round(df_counts_temp[half_samples_da, ct_da] * da)
                           rsums_after <- rowSums(df_counts_temp)
                           df_counts_temp <- round(df_counts_temp / (rsums_after / rsums_before))
                         }
                         
                         # Remove columns if they contain all zeros
                         df_counts_temp <- df_counts_temp %>% select_if(colSums(.) != 0) %>% mutate_all(as.numeric)
                         
                         for (zmet in zero_imp_method) {
                           df_freq <- df_counts_temp %>% calc_perc_df()
                           df_arcsine_sqrt <- asin(sqrt(df_freq/100))
                           
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
                             if (grepl("counts", met)) {df <- df_counts_temp
                             # } else if (grepl("counts_imputed", met)) {df <- df_counts_temp_imputed
                             } else if (grepl("freq", met)) {df <- df_freq
                             # } else if (grepl("freq_imputed", met)) {df <- df_freq_imputed
                             } else if (grepl("arcsine_sqrt", met)) {df <- df_arcsine_sqrt
                             } else if (grepl("alr_mincvref", met)) {
                               ct_ref <- sample(colnames(df_freq_imputed)[!colnames(df_freq_imputed) %in% ct_da], size = 1)
                               df <- Hotelling::alr(as.formula(paste0(ct_ref, "~.")), df_freq_imputed)
                             } else if (grepl("alr_randref", met)) {
                               cvs <- apply(df_freq_imputed, 2, cv)
                               min_cv <- min(cvs[!colnames(df_freq_imputed) %in% ct_da])
                               ct_ref_mincv <- colnames(df_freq_imputed)[which(cvs == min_cv)][1]
                               df <- Hotelling::alr(as.formula(paste0(ct_ref_mincv, "~.")), df_freq_imputed)
                             } else if (grepl("ilr", met)) {df <- compositions::ilr(df_freq)
                             } else if (grepl("clr", met)) {df <- clr(df_freq_imputed)
                             # } else if (met == "clr_centered") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = FALSE)
                             # } else if (met == "clr_centered_scaled") {df <- scale(clr(df_freq_imputed), center = TRUE, scale = TRUE)
                             }
                             
                             if (grepl("pca", met)) {df <- prcomp(df)$x[, 1:2]}
                             
                             
                             # Calculate scores
                             avg_sil <- calc_sil(feat_mat = df, labs)
                             mod <- calc_modularity(feat_mat = df, labs)
                             cluster_score <- clust_eval(matrix = df, labs)
                             
                             
                             # Append results
                             new_row <- list(
                               trans_method = met,
                               zero_imp_method = zmet,
                               n_celltypes = nct,
                               Amount_of_perturbation = da,
                               Silhouette_score = avg_sil,
                               Modularity_score = mod,
                               Adjusted_Rand_Index = cluster_score,
                               bootstrap_id = rep
                             )
                             res <- rbind(res, new_row)
                           }
                         }
                         
                       }
                     }
                     
                     return(res)
                   }
  #stop cluster
  stopCluster(cluster)
  
  return(rets)
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


FindClusters_multi <- function(seurat,
                               res_broad = c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100),
                               res_broad_2sub = c(), # c(0.1, 0.5),
                               res_sub = c() #c(0.1)
                               ) {
  for (rb in res_broad){
    print(paste0("Running broad clustering at resolution: ", rb))
    Idents(seurat) <- "Sample" # Resetting Idents for initial broad clustering
    
    # Run broad clustering
    seurat <- FindClusters(seurat, resolution = rb)
    
    broad_clus_col_name <- paste0("RNA_snn_res.", rb)
    
    # Sub-clustering based on res_broad_2sub
    if (rb %in% res_broad_2sub){
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
                       plot_title = "") {
  df_var <- df %>%
    pivot_longer(
      cols = everything(),
      names_to = "celltype",
      values_to = "values"
    ) %>%
    group_by(celltype) %>%
    dplyr::summarize(
      Relative_abundance = mean(values, na.rm=TRUE),
      Variance = var(values, na.rm=TRUE)
    ) %>%
    arrange(desc(Variance))
  
  if (show_plot == TRUE) {
    p <- varmeanplot(data = df_var, title = plot_title)
    print(p)
  }
  
  return(df_var)
}


get_pb <- function(seurat, sample_col = "Sample", hvg = NULL) {
  pb <- as.matrix(AggregateExpression(seurat, group.by = sample_col, assays = "RNA")[["RNA"]])
  colnames(pb) <- gsub("-", "_", colnames(pb))
  if(!is.null(hvg)) {
    pb <- pb[hvg,]
  }
  return(pb)
}


get_pb_deseq2 <- function(seurat, sample_col = "Sample", hvg) {
  pb <- get_pb(seurat, sample_col = sample_col, hvg = hvg)
  metadata <- get_metadata(seurat)
  metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- t(DESeq2.normalize(pb, metadata = metadata))
  return(pb_norm)
}


varmeanplot <- function(data, title) {
  
  lmod <- lm(Variance ~ Relative_abundance, data = data)
  
  ggplot(data, aes(x = Relative_abundance, y = Variance)) + 
    geom_point() +
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
    labs(title = paste(title)) +
    theme_classic() +
    xlab("CLR-transformed relative abundance") + ylab("Variance")
}



get_metadata <- function(seurat, sample_col = "Sample") {
  metadata <- seurat@meta.data %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    slice(1)
  
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





# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(feat_mat,
                     labels,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = NULL,
                     title = NULL,
                     sil_score = TRUE,
                     mod_score = TRUE,
                     cluster_score = TRUE,
                     pointsize = 3,
                     labelsize = 1,
                     coord_equal = TRUE,
                     axes = c(1, 2),
                     plotly_3d = FALSE,
                     invisible = c("var", "quali"),
                     repel = FALSE) {
  
  res.pca <- prcomp(feat_mat, scale. = scale., rank. = pca_dims)
  
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
  
  p <- factoextra::fviz_pca(res.pca,
                            axes = axes,
                            habillage = labels,
                            label = "var",
                            pointsize = pointsize,
                            labelsize = labelsize,
                            invisible = invisible,
                            repel = repel,
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


run_analyses <- function(result_list,
                         ds,
                         seurat,
                         path_data,
                         factors_test) {
  
  if (is.null(result_list[["bmark"]][[ds]])) {
    print(paste("Running benchmark analysis for dataset: ", ds))
    result_list[["bmark"]][[ds]] <- run_benchmark_analysis(seurat = seurat,
                                                           path_data = path_data,
                                                           factors_test = factors_test)
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
                                   sample_col = "Sample",
                                   n_hvg = 2000,
                                   factors_test,
                                   path_data,
                                   Pseudobulk = TRUE,
                                   ECODA_deconv = TRUE,
                                   ECODA_low_res = TRUE,
                                   ECODA_high_res = TRUE,
                                   ECODA_ultrahigh_res = TRUE,
                                   ECODA_ultrahigh_res_resparam = 20,
                                   ECODA_selectg_top_hvct = TRUE,
                                   ECODA_top_hvct = 0.25,
                                   Pseudobulk_PCA = TRUE,
                                   ECODA_high_res_PCA = TRUE,
                                   MOFA = TRUE,
                                   scITD = TRUE,
                                   MrVI = TRUE,
                                   GloScope = TRUE,
                                   scPoli = TRUE,
                                   PILOT = TRUE,
                                   show_pca_plots = FALSE,
                                   save_pca_plots = TRUE) {
  
  res_list <- list()
  
  files <- list(
    mrvi_dist_file = file.path(path_data, paste0(ds, "_mrvi_dists.feather")),
    scpoli_emb_file = file.path(path_data, paste0(ds, "_scpoli_embs.feather")),
    pilot_dist_file = file.path(path_data, paste0(ds, "_pilot_dists.feather"))
  )
  
  for (file in files) {
    if (!file.exists(file)) {
      stop("File not found: ", file)
    }
  }
  
  
  metadata <- get_metadata(seurat)
  metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  
  labels <- get_labels(seurat, seurat@misc$label_col)
  
  
  if (Pseudobulk | Pseudobulk_PCA | scITD | GloScope | ECODA_ultrahigh_res) {
    seurat <- seurat |>
      NormalizeData() |>
      FindVariableFeatures(nfeatures = n_hvg)

    if ("var.features" %in% slotNames(seurat@assays[["RNA"]])) {
      # This path is for older Seurat objects (primarily v2/v3)
      hvg <- seurat@assays[["RNA"]]@var.features
    } else if ("var.features" %in% colnames(seurat@assays[["RNA"]]@meta.data)) {
      varf <- seurat@assays[["RNA"]]@meta.data[["var.features"]]
      hvg <- varf[!is.na(varf)]
    } else {
      stop("Could not find variable features in expected locations.")
    }
    
    if (GloScope | ECODA_ultrahigh_res) {
      seurat <- seurat |>
        ScaleData() |>
        RunPCA(dims = 1:50, verbose = FALSE)
      
      if (ECODA_ultrahigh_res) {
        seurat <- seurat |>
          FindNeighbors() |>
          FindClusters(resolution = ECODA_ultrahigh_res_resparam)
        ECODA_ultrahigh_res_col_name <- paste0("RNA_snn_res.", ECODA_ultrahigh_res_resparam)
      }
    }
  }
  
  # As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
  # First group.by variable `ident` starts with a number, appending `g` to ensure valid variable names
  # This message is displayed once per session.
  sample_names_starting_with_digit <- grepl("^\\d", unique(seurat[[sample_col]]))
  if (any(sample_names_starting_with_digit)) {
    seurat[[sample_col]][grepl("^\\d", seurat[[sample_col]])] <- paste0("g", seurat[[sample_col]][grepl("^\\d", seurat[[sample_col]])])
  }
  
  if (Pseudobulk | Pseudobulk_PCA) {
    pb_norm <- get_pb_deseq2(seurat, sample_col = "Sample", hvg = hvg)
  }
  
  if (Pseudobulk) {
    res_list[["Pseudobulk"]] <- process_pseudobulk_fig(pb_norm, labels)
  }
  
  if (ECODA_deconv) {
    # Deconvolute using EPIC
    res_list[["ECODA_deconv"]] <- process_deconv_fig(pseudobulk, labels)
  }
  
  # CoDA
  if (ECODA_low_res) {
    ## layer1: low res. cell types
    res_list[["CtCompECODA_lowres"]] <- process_coda_fig(seurat, labels, ct_col = seurat@misc$low_res_ct_col, title = "Cell type composition (CLR)\nlow res.")
  }
  
  if (ECODA_high_res) {
    ## layer2: high res. cell types
    res_list[["CtCompECODA_highres"]] <- process_coda_fig(seurat, labels, ct_col = seurat@misc$hi_res_ct_col, title = "Cell type composition (CLR)\nhigh res.")
    if (ECODA_selectg_top_hvct) {
      res_list[[paste0("CtCompECODA_highres_top", ECODA_top_hvct)]] <-
        process_coda_fig(seurat, labels, ECODA_top_hvct, ct_col = seurat@misc$hi_res_ct_col, title = "Cell type composition (CLR)\nhigh res. top", ECODA_top_hvct*100, "%")
    }
    res_list[["CtCompFreq_highres"]] <- process_coda_fig(seurat, labels, clr = FALSE, ct_col = seurat@misc$hi_res_ct_col, title = "Cell type composition (%)\nhigh res.")
  }
  
  if (ECODA_ultrahigh_res) {
    ## Ultra high res. cell type clusters based on Leiden clustering to artificially increase the number of cell types (clusters), e.g. to 250 cell types (clusters)
    res_list[["CtCompECODA_ultrahighres"]] <- process_coda_fig(seurat, labels, ct_col = ECODA_ultrahigh_res_col_name, title = "Cell type composition (CLR)\nultra-high res.")
    if (ECODA_selectg_top_hvct) {
      res_list[[paste0("CtCompECODA_ultrahighres_top", ECODA_top_hvct)]] <-
        process_coda_fig(seurat, labels, ECODA_top_hvct, ct_col = seurat@misc$hi_res_ct_col, title = "Cell type composition (CLR)\nultra-high res. top", ECODA_top_hvct*100, "%")
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
        res_list[[paste0("ECODA_high_res_", i, "_PCA_dims")]] <- process_coda_fig(seurat, labels, pca_dims = i, ct_col = seurat@misc$hi_res_ct_col, title = paste0("Cell type composition\nhigh res. + PCA (", i, " dims)"))
      }

      if (MOFA) {
        # MOFA
        # Required for MOFA to run
        seurat@version <- package_version('3.1.5')
        res_list[[paste0("MOFA_", i, "_factors")]] <- process_mofa_bulk_fig(pb_norm, metadata = metadata, labels, num_factors = i)
      }

      if (scITD) {
        # scITD
        res_list[[paste0("scITD_", i, "_factors")]] <- process_scitd_fig(seurat, ct_col = seurat@misc$low_res_ct_col, label_col = seurat@misc$label_col, hvg, num_factors = i)
      }
    }
  }
  
  if (MrVI) {
    # MrVI
    res_list[["MrVI"]] <- process_mrvi_fig(seurat, files[["mrvi_dist_file"]], metadata, seurat@misc$label_col)
  }
  
  if (scPoli) {
    res_list[["scPoli"]] <- process_scpoli_fig(files[["scpoli_emb_file"]], metadata, seurat@misc$label_col)
  }
  
  if (PILOT) {
    res_list[["PILOT"]] <- process_pilot_fig(files[["pilot_dist_file"]], metadata, seurat@misc$label_col)
  }
  
  if (GloScope) {
    res_list[["GloScope"]] <- process_gloscope_fig(seurat, metadata, seurat@misc$label_col)
  }
  
  
  if (show_pca_plots | save_pca_plots) {
    plot_list <- lapply(res_list, function(method) method[["plot"]])
    p <- cowplot::plot_grid(plotlist = plot_list)
    res_list[["PCA_plots"]] <- p
  }
  if (show_pca_plots) {print(p)}
  if (save_pca_plots) {
    path <- str_split(mrvi_dist_file, "/")[[1]]
    ds <- path[length(path)]
    ds <- gsub(".nc", "", ds)
    path <- path[1:(length(path)-3)]
    path_plots <- file.path(do.call(file.path, as.list(path)), "plots")
    ggsave(file.path(path_plots, paste0("benchmark_pca_", ds, ".pdf")),
           plot = p,
           width = 15,
           height = 10)
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
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_coda_fig <- function(seurat,
                             labels,
                             ECODA_top_hvct = NULL,
                             clr = TRUE,
                             pca_dims = NULL,
                             ct_col,
                             title,
                             clr_zero_impute_method = "percentage_all",
                             clr_zero_impute_num = 0.1,
                             feat_mat = NULL) {
  
  if (is.null(feat_mat)) {
    df_counts <- get_ct_comp_df_seurat(seurat, sample_col = "Sample", ct_col)
    
    feat_mat <- df_counts %>%
      impute_zeros(clr_zero_impute_method = clr_zero_impute_method, clr_zero_impute_num = clr_zero_impute_num) %>%
      calc_perc_df()
    
    if (clr) {
      feat_mat <- feat_mat %>%
        clr()
    }
    
    if (!is.null(ECODA_top_hvct)) {
      top_hvct <- get_ct_var(feat_mat, show_plot = FALSE) %>%
        slice_head(prop = ECODA_top_hvct) %>%
        pull(celltype)
      feat_mat <- feat_mat[, top_hvct]
    }
  }
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, pca_dims = pca_dims, title = title)
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
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, pca_dims = pca_dims, title = title)
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
    dplyr::select(donors, !!sym(label_col))
  labels <- as.factor(labels_scITD[[label_col]])
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_mrvi_fig <- function(mrvi_dist_file, metadata, label_col, title = "MrVI") {
  
  feat_mat <- arrow::read_feather(mrvi_dist_file)
  samples <- row.names(feat_mat)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(feat_mat)
  
  return(res)
}


process_gloscope_fig <- function(seurat,
                                 metadata,
                                 label_col,
                                 dens = "KNN",
                                 dist_mat = c("KL"),
                                 k = 9,
                                 BPPARAM = BiocParallel::MulticoreParam(workers = parallelly::availableCores() - 2, progressbar = TRUE),
                                 title = "GloScope") {
  
  feat_mat <- GloScope::gloscope(
    embedding_matrix = seurat@reductions$pca@cell.embeddings,
    cell_sample_ids = seurat$Sample,
    dens = dens,
    dist_mat = dist_mat,
    k = k,
    BPPARAM = BPPARAM)
  samples <- row.names(feat_mat)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- feat_mat
  
  return(res)
}


process_scpoli_fig <- function(scpoli_emb_file, metadata, label_col, title = "scPoli") {
  
  feat_mat <- arrow::read_feather(scpoli_emb_file)
  samples <- row.names(feat_mat)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]
  
  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.matrix(dist(feat_mat))
  
  return(res)
}


process_pilot_fig <- function(pilot_dist_file, metadata, label_col, title = "PILOT") {
  
  feat_mat <- arrow::read_feather(pilot_dist_file)
  samples <- row.names(feat_mat)
  labels <- metadata[match(samples, metadata$Sample), ][[label_col]]

  res <- list()
  res[["plot"]] <- plot_pca(feat_mat, labels, title = title)
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
                      methods = c(
                        "counts",
                        "freq",
                        "arcsine_sqrt",
                        "arcsine_sqrt_pca",
                        "alr_randref",
                        "alr_mincvref",
                        "clr",
                        "clr_pca"
                      ))
  res_list <- res_list %>%
    group_by(method) %>%
    summarize(Silhouette_score = mean(Silhouette_score),
              Modularity_score = mean(Modularity_score),
              Adjusted_Rand_Index = mean(Adjusted_Rand_Index)) %>%
    ungroup()
  res_list$dataset <- ds
  
  return(res_list)
}


run_zeroimp_analysis <- function(ct_comps, labels) {
  
  test <- colSums(df) == 0
  
  df <- ct_comps %>% select_if(colSums(.) != 0) %>% mutate_all(as.numeric)
  
  res_list <- list()
  res_list[["counts_zeros_1"]] <- df %>% impute_zeros(clr_zero_impute_method = "counts_zeros", clr_zero_impute_num = 1) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["counts_all_1"]] <- df %>% impute_zeros(clr_zero_impute_method = "counts_all", clr_zero_impute_num = 1) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["perc_all_0.001%"]] <- df %>% impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.001) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["perc_all_0.01%"]] <- df %>% impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.01) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["perc_all_0.1%"]] <- df %>% impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 0.1) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["perc_all_1%"]] <- df %>% impute_zeros(clr_zero_impute_method = "percentage_all", clr_zero_impute_num = 1) %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  res_list[["asinsqrt"]] <- df %>% calc_perc_df() %>% mutate(across(everything(), ~ . / 100)) %>% sqrt() %>% asin() %>% calc_sep_score(labels)
  # res_list[["impRZilr"]] <- df %>% {impRZilr(.)[["x"]]} %>% calc_perc_df() %>% clr() %>% calc_sep_score(labels)
  # Need to check if multLN drops samples
  df_multLN <- df %>% calc_perc_df() %>% zCompositions::multLN(label = 0, dl = rep(0.1, ncol(df)), z.warning = 0.9)
  labels_multLN <- labels[row.names(df) %in% row.names(df_multLN)]
  res_list[["multLN"]] <- df_multLN %>% clr() %>% calc_sep_score(labels_multLN)
  res_list[["multRepl_0.01%"]] <- df %>% calc_perc_df() %>% zCompositions::multRepl(label = 0, dl = rep(0.01, ncol(df)), z.warning = 1, frac = 1) %>% clr() %>% calc_sep_score(labels)
  res_list[["multRepl_0.1%"]] <- df %>% calc_perc_df() %>% zCompositions::multRepl(label = 0, dl = rep(0.1, ncol(df)), z.warning = 1, frac = 1) %>% clr() %>% calc_sep_score(labels)
  res_list[["multRepl_1%"]] <- df %>% calc_perc_df() %>% zCompositions::multRepl(label = 0, dl = rep(1, ncol(df)), z.warning = 1, frac = 1) %>% clr() %>% calc_sep_score(labels)
  
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
      paste0("Cell_", 1:n_cells)  # Cell names
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
  message(paste0("Created Seurat object with ", ncol(seurat_object), " cells and ",
                 nrow(seurat_object), " genes."))
  message(paste0("Number of unique samples: ", length(unique(seurat_object$Sample))))
  message(paste0("Number of unique cell types: ", length(unique(seurat_object$CellType))))
  message(paste0("Number of unique groups: ", length(unique(seurat_object$Group))))
  
  # Return the Seurat object
  return(seurat_object)
}

