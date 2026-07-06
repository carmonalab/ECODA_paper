# ============================================================
# SEURAT UTILITY FUNCTIONS
# ============================================================

# Create clean v5 seurat object (only counts + metadata)
create_clean_seuratv5_object <- function(seurat) {
  if (!is.null(seurat@assays$RNA$counts)) {
    seurat <- CreateSeuratObject(
      counts = seurat@assays$RNA$counts,
      meta.data = seurat@meta.data
    )
  } else if (!is.null(seurat@assays[["RNA"]]@layers[["X"]])) {
    seurat <- CreateSeuratObject(
      counts = seurat@assays$RNA$counts,
      meta.data = seurat@meta.data
    )
  }
  return(seurat)
}

# Multi-resolution clustering with optional sub-clustering
FindClusters_multi <- function(
  seurat,
  res_broad = c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100),
  res_broad_2sub = c(),
  res_sub = c()
) {
  for (rb in res_broad) {
    print(paste0("Running broad clustering at resolution: ", rb))
    Idents(seurat) <- "Sample"
    seurat <- FindClusters(seurat, resolution = rb)
    broad_clus_col_name <- paste0("RNA_snn_res.", rb)

    if (rb %in% res_broad_2sub) {
      print(paste0("Proceeding with sub-clustering for broad resolution: ", rb))
      main_clusters <- levels(seurat@meta.data[, broad_clus_col_name])

      for (rs in res_sub) {
        print(paste0("  Sub-clustering at resolution: ", rs))
        subcluster_col_name <- paste0(
          "subcluster_broad_res",
          rb,
          "_sub_res",
          rs
        )
        seurat@meta.data[, subcluster_col_name] <-
          as.character(seurat@meta.data[, broad_clus_col_name])

        for (cluster_id in main_clusters) {
          print(paste(
            "    Subclustering broad cluster:",
            cluster_id,
            " / ",
            length(main_clusters)
          ))
          Idents(seurat) <- broad_clus_col_name
          sub_seurat <- subset(seurat, idents = cluster_id)

          if (ncol(sub_seurat) < 50) {
            print(paste(
              "    Skipping subclustering for broad cluster",
              cluster_id,
              "as it contains too few cells."
            ))
            next
          }

          sub_seurat <- sub_seurat |>
            FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
            ScaleData(verbose = FALSE) |>
            RunPCA(dims = 1:50, verbose = FALSE)

          if (is.null(sub_seurat@reductions$pca)) {
            print(paste(
              "    Skipping subclustering for broad cluster",
              cluster_id,
              "due to PCA failure."
            ))
            next
          }

          n_pcs_available <- ncol(sub_seurat@reductions$pca@cell.embeddings)
          if (n_pcs_available == 0) {
            print(paste(
              "    Skipping subclustering for broad cluster",
              cluster_id,
              "due to 0 available PCs."
            ))
            next
          }

          sub_seurat <- FindNeighbors(sub_seurat, dims = dims, verbose = FALSE)
          sub_seurat <- FindClusters(
            sub_seurat,
            resolution = rs,
            verbose = FALSE
          )

          sub_cluster_idents_col <- paste0("RNA_snn_res.", rs)
          sub_cluster_ids <- sub_seurat@meta.data[[sub_cluster_idents_col]]
          if (!is.factor(sub_cluster_ids)) {
            sub_cluster_ids <- factor(sub_cluster_ids)
          }

          new_subcluster_names <- paste0(
            cluster_id,
            "_",
            levels(sub_cluster_ids)
          )
          names(new_subcluster_names) <- levels(sub_cluster_ids)
          mapped_sub_cluster_ids <- factor(
            sub_cluster_ids,
            levels = levels(sub_cluster_ids),
            labels = new_subcluster_names
          )

          seurat@meta.data[Cells(sub_seurat), subcluster_col_name] <-
            as.character(mapped_sub_cluster_ids)
        }

        Idents(seurat) <- subcluster_col_name
        seurat@meta.data[, subcluster_col_name] <-
          as.factor(seurat@meta.data[, subcluster_col_name])
      }
    }
  }
  return(seurat)
}

# Get cell type composition dataframe from seurat
get_ct_comp_df_seurat <- function(seurat, sample_col, ct_col) {
  ct_comp_df <- table(
    seurat@meta.data[[sample_col]],
    seurat@meta.data[[ct_col]]
  ) %>%
    t() %>%
    as.data.frame.matrix() %>%
    t() %>%
    as.data.frame()
  ct_comp_df <- ct_comp_df[rowSums(ct_comp_df) != 0, ]
  return(ct_comp_df)
}

# Remove samples with low cell counts
remove_low_cellcount_samples <- function(
  seurat,
  sample_col = "Sample",
  min_cells_per_sample = 500,
  show_plot = FALSE
) {
  Idents(seurat) <- sample_col
  cells_per_sample <- table(seurat@meta.data[[sample_col]])
  seurat <- subset(
    seurat,
    idents = names(cells_per_sample[cells_per_sample > min_cells_per_sample])
  )
  if (show_plot) {
    barplot(sort(cells_per_sample))
    abline(v = min_cells_per_sample)
  }
  return(seurat)
}

# Replace HiTME layer3 annotation based on UCell score threshold
replace_HiTMElayer3_annot <- function(seurat, thresh = 0.1) {
  seurat$layer3 <- gsub(
    pattern = "_resting|_IFN|_cellCycle.G1S|_cellCycle.G2M|_HeatShock|_Heatshock|_Prolif",
    replacement = "",
    seurat$layer3
  )

  simple_conditions <- list(
    IFN_UCell = list(threshold = thresh, suffix = "_IFN")
  )

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

  seurat@meta.data <- seurat@meta.data %>%
    mutate(
      layer3 = if_else(
        !is.na(layer3) &
          (cellCycle.G1S_UCell > thresh | cellCycle.G2M_UCell > thresh),
        paste0(layer3, "_Prolif"),
        layer3
      )
    )
  return(seurat)
}

# Load h5ad file to seurat object
load_h5ad_to_seurat <- function(file_name) {
  seurat <- read_h5ad(file_name, as = "Seurat")
  if (!is.null(seurat@assays$RNA$X)) {
    seurat[["RNA"]]$counts <- seurat[["RNA"]]$X
    seurat[["RNA"]]$X <- NULL
    gc()
  }

  file_name_short <- tools::file_path_sans_ext(file_name)
  feather_path <- paste0(file_name_short, ".feather")

  if (file.exists(feather_path)) {
    arrow::set_cpu_count(1)
    meta <- arrow::read_feather(feather_path)
    arrow::set_cpu_count(parallelly::availableCores() - 2)
    seurat <- Seurat::AddMetaData(seurat, metadata = as.data.frame(meta))
  }
  return(seurat)
}

# Get current variable features from seurat object
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

# Standardize sample names (prevent leading digits, replace hyphens)
standardize_sample_names <- function(sample_names) {
  sample_names_starting_with_digit <- grepl("^\\d", unique(sample_names))
  if (any(sample_names_starting_with_digit)) {
    sample_names[sample_names_starting_with_digit] <- paste0(
      "g",
      sample_names[sample_names_starting_with_digit]
    )
  }
  sample_names <- gsub("-", "_", sample_names)
}

# Get unique metadata per sample
get_metadata <- function(seurat, sample_col = "Sample") {
  metadata <- seurat@meta.data %>%
    dplyr::group_by(!!sym(sample_col)) %>%
    dplyr::slice(1)

  return(metadata)
}

# Get labels from seurat metadata
get_labels <- function(seurat, label_col) {
  metadata <- get_metadata(seurat)
  labels <- as.factor(metadata[[label_col]])
  names(labels) <- metadata[["Sample"]]

  return(labels)
}
