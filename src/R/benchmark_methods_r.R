# ============================================================
# BENCHMARK METHOD PROCESSING FUNCTIONS
# Each function computes a feature matrix and returns a result bundle
# ============================================================

# EPIC deconvolution
process_deconv_fig <- function(
  pseudobulk,
  labels
) {
  out <- EPIC(pseudobulk, BRef)
  deconv_ct_comps <- as.data.frame(out[["mRNAProportions"]])
  row.names(deconv_ct_comps) <- colnames(pseudobulk)
  deconv_ct_comps[deconv_ct_comps == 0] <- deconv_ct_comps[
    deconv_ct_comps == 0
  ] +
    (2 / 3) * min(deconv_ct_comps)
  feat_mat <- clr(deconv_ct_comps)
  return(create_result_bundle(feat_mat, labels))
}

# CoDA (compositional data analysis)
process_coda_fig <- function(
  seurat,
  labels,
  ECODA_top_n_hvct = NULL,
  ECODA_top_varexp_hvct = NULL,
  hvct_recalc_clr = TRUE,
  calc_clr = TRUE,
  pca_dims = NULL,
  sample_col = "Sample",
  ct_col,
  clr_zero_impute_method = "counts_all",
  clr_zero_impute_num = 0.5,
  feat_mat = NULL,
  var_ct_desc = TRUE,
  shuffle_labels = FALSE
) {
  if (is.null(feat_mat)) {
    df_counts <- get_ct_comp_df_seurat(seurat, sample_col = sample_col, ct_col)
    df_imp <- df_counts %>%
      impute_zeros(
        clr_zero_impute_method = clr_zero_impute_method,
        clr_zero_impute_num = clr_zero_impute_num
      )
    if (calc_clr) {
      feat_mat <- df_imp %>% clr()
    } else {
      feat_mat <- df_imp %>% calc_perc_df()
    }
    top_hvct <- NULL
    if (!is.null(ECODA_top_n_hvct)) {
      top_hvct <- get_ct_var(
        feat_mat,
        show_plot = FALSE,
        descending = var_ct_desc
      ) %>%
        get_hvcs(top_n_hvcs = ECODA_top_n_hvct, variance_threshold = NULL)
      feat_mat <- feat_mat[, top_hvct]
    }
    if (!is.null(ECODA_top_varexp_hvct)) {
      top_hvct <- get_ct_var(
        feat_mat,
        show_plot = FALSE,
        descending = var_ct_desc
      ) %>%
        get_hvcs(top_n_hvcs = NULL, variance_threshold = ECODA_top_varexp_hvct)
      feat_mat <- feat_mat[, top_hvct]
    }
    if (hvct_recalc_clr & !is.null(top_hvct)) {
      feat_mat <- df_counts[, top_hvct] %>%
        impute_zeros(
          clr_zero_impute_method = clr_zero_impute_method,
          clr_zero_impute_num = clr_zero_impute_num
        ) %>%
        clr()
    }
  }
  if (!is.null(pca_dims)) {
    feat_mat <- prcomp(feat_mat, rank. = pca_dims)[["x"]]
  }
  dist_mat <- dist(feat_mat)
  labels <- labels[rownames(feat_mat)]
  if (shuffle_labels) {
    set.seed(123)
    labels <- sample(labels)
  }
  res <- list()
  res[["scores"]] <- calc_sep_score(dist_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- dist_mat
  res[["labels"]] <- labels
  res[["counts"]] <- df_imp
  return(res)
}

# Pseudobulk processing
process_pseudobulk_fig <- function(
  feat_mat,
  labels,
  pca_dims = NULL,
  knn_k = NULL
) {
  if (!is.null(pca_dims)) {
    feat_mat <- prcomp(feat_mat, rank. = pca_dims)[["x"]]
  }
  return(create_result_bundle(feat_mat, labels))
}

# Cell type pseudobulk processing
process_pseudobulk_ct_fig <- function(
  seurat,
  labels,
  hvg = 500,
  sample_col = "Sample",
  ct_col,
  min_cells = 5
) {
  all_samples <- sort(unique(seurat[[sample_col, drop = TRUE]]))
  n_samples <- length(all_samples)
  total_dist <- matrix(
    0,
    nrow = n_samples,
    ncol = n_samples,
    dimnames = list(all_samples, all_samples)
  )
  count_mat <- matrix(
    0,
    nrow = n_samples,
    ncol = n_samples,
    dimnames = list(all_samples, all_samples)
  )
  cell_types <- unique(seurat[[ct_col, drop = TRUE]])
  cell_types <- cell_types[!is.na(cell_types)]

  for (ct in cell_types) {
    sub <- subset(x = seurat, subset = !!sym(ct_col) == ct)
    counts_per_sample <- table(sub@meta.data[, sample_col])
    keep_samples <- names(counts_per_sample)[counts_per_sample >= min_cells]
    if (length(keep_samples) < 2) {
      next
    }
    sub <- subset(sub, subset = !!sym(sample_col) %in% keep_samples)
    pb_norm <- get_pb_deseq2(sub, sample_col = sample_col, n_hvg = hvg)
    dist_mat_ct <- as.matrix(dist(pb_norm))
    total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] + dist_mat_ct
    count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] + 1
  }
  final_dist_mat <- total_dist / count_mat
  final_dist_mat[is.nan(final_dist_mat)] <- 0
  return(create_result_bundle(
    feat_mat = final_dist_mat,
    labels,
    dist_mat = as.dist(final_dist_mat)
  ))
}

# Average PCA embedding
process_avg_pca_embedding_fig <- function(
  seurat,
  labels,
  sample_col = "Sample"
) {
  feat_mat <- as.data.frame(seurat@reductions$pca@cell.embeddings)
  feat_mat$Sample <- seurat@meta.data[[sample_col]]
  feat_mat <- feat_mat %>%
    group_by(.data$Sample) %>%
    summarise(across(starts_with("PC_"), mean)) %>%
    ungroup() %>%
    column_to_rownames(var = "Sample")
  return(create_result_bundle(feat_mat, labels))
}

# MOFA processing
process_mofa_bulk_fig <- function(
  pb_norm,
  metadata,
  labels,
  num_factors = 5,
  maxiter = 1000
) {
  pb_list <- list(pb = t(pb_norm))
  MOFAobject <- create_mofa(pb_list)
  metadata$sample <- metadata$Sample
  samples_metadata(MOFAobject) <- metadata
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- num_factors
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
  return(create_result_bundle(feat_mat, labels))
}

# scITD processing
process_scitd_fig <- function(
  seurat,
  ct_col,
  label_col,
  hvg,
  num_factors = 5
) {
  seurat$donors <- seurat$Sample
  seurat$ctypes <- as.character(seurat@meta.data[[ct_col]])
  ctypes <- unique(seurat$ctypes)
  ctypes <- ctypes[!is.na(ctypes)]
  cell_counts <- table(seurat$donors, seurat$ctypes)
  ctypes_to_drop <- names(which(colMeans(cell_counts < 5) > 0.2))
  ctypes <- ctypes[!ctypes %in% ctypes_to_drop]
  param_list <- initialize_params(
    ctypes_use = as.character(ctypes),
    ncores = parallelly::availableCores() - 2,
    rand_seed = 10
  )
  pbmc_container <- make_new_container(
    count_data = seurat@assays$RNA$counts,
    meta_data = seurat@meta.data[, c("donors", "ctypes", label_col, "Sample")],
    gn_convert = NULL,
    params = param_list
  )
  pbmc_container <- form_tensor(pbmc_container, custom_genes = hvg)
  pbmc_container <- run_tucker_ica(
    pbmc_container,
    ranks = c(num_factors, num_factors + 5)
  )
  feat_mat <- pbmc_container[["tucker_results"]][[1]]
  labels_scITD <- seurat@meta.data %>%
    filter(.data$donors %in% row.names(feat_mat)) %>%
    distinct(.data$donors, .keep_all = TRUE) %>%
    dplyr::select(.data$donors, !!sym(label_col))
  labels <- as.factor(labels_scITD[[label_col]])
  names(labels) <- labels_scITD[["donors"]]
  return(create_result_bundle(feat_mat, labels))
}

# MrVI processing
process_mrvi_fig <- function(mrvi_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(mrvi_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# GloScope processing
process_gloscope_fig <- function(
  seurat,
  metadata,
  label_col,
  gloscope_dist_file,
  n_pca_dims = 30,
  dens = "KNN",
  dist_metric = c("KL"),
  k = 25
) {
  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(
      workers = parallelly::availableCores() - 2,
      progressbar = TRUE
    )
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(
      workers = parallelly::availableCores() - 2,
      progressbar = TRUE
    )
  }
  if (!file.exists(gloscope_dist_file)) {
    feat_mat <- GloScope::gloscope(
      embedding_matrix = seurat@reductions$pca@cell.embeddings[, 1:n_pca_dims],
      cell_sample_ids = seurat$Sample,
      dens = dens,
      dist_metric = dist_metric,
      k = k,
      BPPARAM = BPPARAM
    )
    saveRDS(feat_mat, file = gloscope_dist_file)
  } else {
    feat_mat <- readRDS(gloscope_dist_file)
  }
  row.names(feat_mat) <- standardize_sample_names(row.names(feat_mat))
  labels <- metadata[match(row.names(feat_mat), metadata$Sample), ][[label_col]]
  names(labels) <- metadata[match(row.names(feat_mat), metadata$Sample), ][[
    "Sample"
  ]]
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# GloScope sqrt matrix processing
process_gloscope_sqrtmat_fig <- function(
  metadata,
  label_col,
  gloscope_dist_file
) {
  if (!file.exists(gloscope_dist_file)) {
    stop(paste(gloscope_dist_file, "not found!"))
  }
  feat_mat <- readRDS(gloscope_dist_file)
  feat_mat <- sqrt(feat_mat)
  feat_mat[is.na(feat_mat)] <- 0
  row.names(feat_mat) <- standardize_sample_names(row.names(feat_mat))
  labels <- metadata[match(row.names(feat_mat), metadata$Sample), ][[label_col]]
  names(labels) <- metadata[match(row.names(feat_mat), metadata$Sample), ][[
    "Sample"
  ]]
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# GloProp processing
process_gloprop_fig <- function(
  seurat,
  metadata,
  ct_col,
  label_col,
  sample_col = "Sample",
  dist_metric = c("KL")
) {
  sample_id <- seurat@meta.data[[sample_col]]
  cluster_id <- seurat@meta.data[[ct_col]]
  dist_result <- gloscopeProp(
    sample_id,
    cluster_id,
    ep = 0.5,
    dist_metric = dist_metric
  )
  feat_mat <- sqrt(dist_result)
  row.names(feat_mat) <- standardize_sample_names(row.names(feat_mat))
  labels <- metadata[match(row.names(feat_mat), metadata$Sample), ][[label_col]]
  names(labels) <- metadata[match(row.names(feat_mat), metadata$Sample), ][[
    "Sample"
  ]]
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# scPoli processing
process_scpoli_fig <- function(scpoli_emb_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(scpoli_emb_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  return(create_result_bundle(feat_mat, labels))
}

# PILOT processing
process_pilot_fig <- function(pilot_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(pilot_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  rownames(feat_mat) <- standardize_sample_names(rownames(feat_mat))
  return(create_result_bundle(feat_mat, labels))
}

# ECODA-PB combo processing
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
    dist_mat_ecoda_normed <- proxy::dist(feat_mat_ecoda, method = "cosine")
    dist_mat_pb_normed <- proxy::dist(feat_mat_pb, method = "cosine")
  } else {
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
  feat_mat <- dist_mat_ecoda_normed *
    ecoda_weight +
    dist_mat_pb_normed * (1 - ecoda_weight)
  res <- list()
  res[["scores"]] <- calc_sep_score(feat_mat, labels)
  res[["feat_mat"]] <- feat_mat
  res[["dist_mat"]] <- as.dist(feat_mat)
  res[["labels"]] <- labels

  return(res)
}

# Create result bundle with scores, feature matrix, distance matrix, labels
create_result_bundle <- function(
  feat_mat,
  labels,
  dist_mat = NULL,
  extra = list()
) {
  if (is.null(dist_mat)) {
    dist_mat <- dist(feat_mat)
  } else if (is.matrix(dist_mat)) {
    dist_mat <- as.dist(dist_mat)
  }
  labels <- labels[rownames(feat_mat)]
  result <- list(
    scores = calc_sep_score(dist_mat, labels),
    feat_mat = feat_mat,
    dist_mat = dist_mat,
    labels = labels
  )
  result <- c(result, extra)

  return(result)
}
