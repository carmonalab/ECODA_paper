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
  shuffle_labels = FALSE,
  obs = NULL
) {
  if (is.null(feat_mat)) {
    if (!is.null(obs)) {
      df_counts <- get_ct_comp_df(obs, sample_col = sample_col, ct_col)
    } else {
      df_counts <- get_ct_comp_df_seurat(seurat, sample_col = sample_col, ct_col)
    }
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
  if (shuffle_labels) {
    label_ids <- names(labels)
    set.seed(123)
    labels <- labels[sample(seq_along(labels))]
    names(labels) <- label_ids
  }
  res <- create_result_bundle(
    feat_mat,
    labels,
    dist_mat = dist_mat,
    extra = list(counts = df_imp)
  )
  return(res)
}

# Correct CLR composition by subtracting only the fitted technical batch
# random effect. Biological labels are deliberately absent from this model.
correct_clr_batch_lmm <- function(feat_mat, sample_meta, batch_col) {
  if (is.null(batch_col) || !nzchar(batch_col)) {
    stop("correct_clr_batch_lmm: batch_col is required")
  }
  feat_mat <- as.matrix(feat_mat)
  if (is.null(rownames(feat_mat)) || anyDuplicated(rownames(feat_mat))) {
    stop("correct_clr_batch_lmm: feature matrix needs unique sample rownames")
  }
  if (is.null(sample_meta) || !batch_col %in% colnames(sample_meta)) {
    stop("correct_clr_batch_lmm: batch column missing from sample metadata")
  }
  if (!"Sample" %in% colnames(sample_meta)) {
    stop("correct_clr_batch_lmm: sample metadata needs a Sample column")
  }
  sample_ids <- as.character(sample_meta[["Sample"]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("correct_clr_batch_lmm: missing sample or batch IDs")
  }
  if (anyDuplicated(sample_ids)) {
    stop("correct_clr_batch_lmm: duplicate sample IDs in metadata")
  }
  if (!identical(sample_ids, rownames(feat_mat))) {
    stop("correct_clr_batch_lmm: sample-order mismatch")
  }
  batch <- sample_meta[[batch_col]]
  valid <- !is.na(batch) & nzchar(as.character(batch))
  if (!all(valid)) {
    stop("correct_clr_batch_lmm: missing sample or batch IDs")
  }
  batch <- factor(batch)
  if (nlevels(batch) < 2) {
    stop("correct_clr_batch_lmm: fewer than two batch levels")
  }

  corrected <- feat_mat
  model_data <- data.frame(batch = batch)
  for (feature in seq_len(ncol(feat_mat))) {
    model_data$y <- as.numeric(feat_mat[, feature])
    fit <- tryCatch(
      lme4::lmer(y ~ 1 + (1 | batch), data = model_data, REML = TRUE),
      error = function(e) {
        stop("correct_clr_batch_lmm: nonconvergence for feature ",
             colnames(feat_mat)[feature], ": ", conditionMessage(e))
      }
    )
    convergence <- fit@optinfo$conv$lme4$messages
    if (!is.null(convergence)) {
      stop("correct_clr_batch_lmm: nonconvergence for feature ",
           colnames(feat_mat)[feature], ": ",
           paste(convergence, collapse = "; "))
    }
    random_effects <- lme4::ranef(fit)$batch[["(Intercept)"]]
    names(random_effects) <- rownames(lme4::ranef(fit)$batch)
    corrected[, feature] <- model_data$y - random_effects[as.character(batch)]
  }
  corrected <- corrected - rowMeans(corrected)
  dimnames(corrected) <- dimnames(feat_mat)
  if (any(abs(rowSums(corrected)) > 1e-8)) {
    stop("correct_clr_batch_lmm: row recentering failed to restore zero sums")
  }
  corrected
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
  successful_cts <- character()
  n_ct_pair_contributions <- 0L

  for (ct in cell_types) {
    sub <- subset(x = seurat, subset = !!sym(ct_col) == ct)
    counts_per_sample <- table(sub@meta.data[, sample_col])
    keep_samples <- names(counts_per_sample)[counts_per_sample >= min_cells]
    if (length(keep_samples) < 2) {
      next
    }
    sub <- subset(sub, subset = !!sym(sample_col) %in% keep_samples)
    # One sparse CT must never kill the whole method: skip on error
    # (label alignment is centralized in create_result_bundle()).
    pb_norm <- tryCatch(
      get_pb_deseq2(sub, sample_col = sample_col, n_hvg = hvg),
      error = function(e) {
        warning(
          "process_pseudobulk_ct_fig: pseudobulk failed for cell type '", ct,
          "': ", conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(pb_norm)) {
      next
    }
    dist_mat_ct <- as.matrix(dist(pb_norm))
    total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] + dist_mat_ct
    count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] + 1
    successful_cts <- c(successful_cts, as.character(ct))
    n_ct_pair_contributions <- n_ct_pair_contributions +
      sum(upper.tri(dist_mat_ct))
  }

  n_sample_pairs_contributed <- sum(
    count_mat[upper.tri(count_mat)] > 0
  )
  if (length(successful_cts) == 0L || n_sample_pairs_contributed == 0L) {
    stop(
      "process_pseudobulk_ct_fig: no successful cell-type pseudobulks ",
      "contributed sample-pair distances for ct_col='", ct_col,
      "', hvg=", hvg, " (", length(cell_types), " cell types inspected)."
    )
  }

  final_dist_mat <- total_dist / count_mat
  final_dist_mat[is.nan(final_dist_mat)] <- 0
  return(create_result_bundle(
    feat_mat = final_dist_mat,
    labels,
    dist_mat = as.dist(final_dist_mat),
    extra = list(
      n_ct_success = length(successful_cts),
      successful_cell_types = successful_cts,
      n_sample_pairs_contributed = n_sample_pairs_contributed,
      n_ct_pair_contributions = n_ct_pair_contributions
    )
  ))
}

# Average PCA embedding
process_avg_pca_embedding_fig <- function(
  seurat,
  labels,
  sample_col = "Sample",
  pca_emb = NULL,
  obs = NULL
) {
  if (!is.null(pca_emb)) {
    feat_mat <- as.data.frame(pca_emb)
    # obsm numpy arrays carry no column names; the per-sample mean step
    # selects starts_with("PC_") (Seurat's reduction naming convention).
    if (is.null(colnames(feat_mat))) {
      colnames(feat_mat) <- paste0("PC_", seq_len(ncol(feat_mat)))
    }
    # Align the embedding cells to their sample ids via cell barcodes
    # (rownames). Both come from the same h5ad, so positions match; a missing
    # rownames fallback covers matrices materialized without the index.
    if (is.null(rownames(feat_mat)) && !is.null(obs) &&
        nrow(feat_mat) == nrow(obs)) {
      rownames(feat_mat) <- rownames(obs)
    }
    if (!is.null(obs) && !is.null(rownames(feat_mat))) {
      sample_ids <- obs[[sample_col]]
      names(sample_ids) <- rownames(obs)
      feat_mat$Sample <- sample_ids[rownames(feat_mat)]
    } else {
      stop("process_avg_pca_embedding_fig: cannot align pca_emb to obs ",
           "(mismatched rownames/lengths)")
    }
  } else {
    feat_mat <- as.data.frame(seurat@reductions$pca@cell.embeddings)
    feat_mat$Sample <- seurat@meta.data[[sample_col]]
  }
  feat_mat <- feat_mat %>%
    group_by(.data$Sample) %>%
    dplyr::summarise(across(starts_with("PC_"), mean)) %>%
    ungroup() %>%
    tibble::column_to_rownames(var = "Sample")
  return(create_result_bundle(feat_mat, labels))
}

# MOFA requires row names, not merely a `sample` metadata column, to match
# the sample names in each view. Benchmark metadata is cell-indexed upstream,
# so set the explicit per-sample row-name contract immediately before the
# samples_metadata<- assignment.
prepare_mofa_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata)
  if (!"Sample" %in% colnames(metadata)) {
    stop("MOFA metadata requires a 'Sample' column.")
  }
  sample_ids <- as.character(metadata$Sample)
  if (anyNA(sample_ids) || any(!nzchar(trimws(sample_ids))) ||
      anyDuplicated(sample_ids)) {
    stop("MOFA metadata sample IDs must be nonmissing, non-empty, and unique.")
  }
  rownames(metadata) <- sample_ids
  metadata$sample <- sample_ids
  metadata
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
  metadata <- prepare_mofa_metadata(metadata)
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
  # Feather rownames already carry the upstream-standardized sample names
  # (written from obs); do NOT re-standardize (see process_gloprop_fig).
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# GloScope processing (sqrtmat variant merged; the sqrt transform + NA->0 is
# now always applied, matching the legacy GloScope_*_sqrtmat results)
process_gloscope_fig <- function(
  embedding_matrix,
  sample_ids,
  metadata,
  label_col,
  gloscope_dist_file,
  n_pca_dims = 30,
  dens = "KNN",
  dist_metric = c("KL"),
  k = 25,
  force = FALSE
) {
  n_samples <- length(unique(sample_ids))
  if (k >= n_samples) {
    k <- min(k, n_samples - 1)
    warning(paste(
      "GloScope: k adjusted to", k,
      "(n_samples - 1 =", n_samples - 1, ") for the", n_samples,
      "sample dataset."
    ))
  }
  if (n_pca_dims > ncol(embedding_matrix)) {
    warning(paste(
      "GloScope: n_pca_dims adjusted to", ncol(embedding_matrix),
      "(available PCs) for the requested", n_pca_dims, "dims."
    ))
    n_pca_dims <- ncol(embedding_matrix)
  }
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
  if (force || !artifact_checksum_ok(gloscope_dist_file)) {
    feat_mat <- GloScope::gloscope(
      embedding_matrix = embedding_matrix[, 1:n_pca_dims],
      cell_sample_ids = sample_ids,
      dens = dens,
      dist_metric = dist_metric,
      k = k,
      BPPARAM = BPPARAM
    )
    save_rds_atomic(feat_mat, gloscope_dist_file)
  } else {
    feat_mat <- readRDS(gloscope_dist_file)
  }
  feat_mat <- sqrt(feat_mat)
  feat_mat[is.na(feat_mat)] <- 0
  # GloScope preserves the exact canonical sample IDs supplied in sample_ids.
  # Keep them unchanged so labels retain the canonical metadata identifiers.
  # Keep labels in canonical obs order; create_result_bundle reorders the
  # distance matrix/feature rows rather than adopting GloScope's order.
  labels <- as.factor(metadata[[label_col]])
  names(labels) <- as.character(metadata[["Sample"]])
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# GloProp processing
process_gloprop_fig <- function(
  seurat,
  metadata,
  ct_col,
  label_col,
  sample_col = "Sample",
  dist_metric = c("KL"),
  obs = NULL
) {
  if (!is.null(obs)) {
    sample_id <- obs[[sample_col]]
    cluster_id <- as.character(obs[[ct_col]])
  } else {
    sample_id <- seurat@meta.data[[sample_col]]
    cluster_id <- as.character(seurat@meta.data[[ct_col]])
  }
  valid_mask <- !is.na(cluster_id) & !cluster_id %in% c("NA", "nan", "None", "Unknown") & cluster_id != ""
  sample_id <- sample_id[valid_mask]
  cluster_id <- cluster_id[valid_mask]

  dist_result <- gloscopeProp(
    sample_id,
    cluster_id,
    ep = 0.5,
    dist_metric = dist_metric
  )
  feat_mat <- sqrt(dist_result)
  # Sample names are standardized upstream (1.1.1_preprocess.py): do NOT
  # re-standardize here (hyphen->underscore) or names diverge from the obs
  # labels for h5ads that predate the python change (e.g. Adams).
  # Keep labels in canonical obs order; create_result_bundle enforces the
  # complete source universe and aligns the returned matrix.
  labels <- as.factor(metadata[[label_col]])
  names(labels) <- as.character(metadata[["Sample"]])
  return(create_result_bundle(feat_mat, labels, dist_mat = as.dist(feat_mat)))
}

# scPoli processing
process_scpoli_fig <- function(scpoli_emb_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(scpoli_emb_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  # Feather rownames already carry the upstream-standardized sample names
  # (written from obs); do NOT re-standardize (see process_gloprop_fig).
  return(create_result_bundle(feat_mat, labels))
}

# PILOT processing
process_pilot_fig <- function(pilot_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(pilot_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  # Feather rownames already carry the upstream-standardized sample names
  # (written from obs); do NOT re-standardize (see process_gloprop_fig).
  return(create_result_bundle(feat_mat, labels))
}

# QOT processing (same layout as PILOT: sample x sample distance matrix,
# plain DataFrame.to_feather() with the pandas index = sample names)
process_qot_fig <- function(qot_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(qot_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  # Feather rownames already carry the upstream-standardized sample names
  # (written from obs); do NOT re-standardize (see process_gloprop_fig).
  return(create_result_bundle(feat_mat, labels))
}

# PILOT-GM-VAE processing (same layout as PILOT: sample x sample distance
# matrix, plain DataFrame.to_feather() with the pandas index = sample names)
process_pilotgm_fig <- function(pilotgm_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- arrow::read_feather(pilotgm_dist_file) %>%
    tibble::column_to_rownames(var = names(.)[ncol(.)]) %>%
    as.data.frame()
  arrow::set_cpu_count(parallelly::availableCores() - 2)
  # Feather rownames already carry the upstream-standardized sample names
  # (written from obs); do NOT re-standardize (see process_gloprop_fig).
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

# Align a result's feature rows and optional distance matrix to the named
# sample order supplied by the canonical preprocessed obs metadata.
align_result_samples <- function(feat_mat, labels, dist_mat = NULL) {
  sample_ids <- rownames(feat_mat)
  label_ids <- names(labels)
  if (is.null(sample_ids) || length(sample_ids) == 0L ||
      anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) ||
      is.null(label_ids) || length(label_ids) == 0L ||
      anyNA(label_ids) || any(!nzchar(label_ids)) ||
      anyDuplicated(label_ids) || anyNA(labels)) {
    stop(
      "Result sample identifiers must be nonmissing, non-empty, and unique."
    )
  }
  if (length(sample_ids) != length(label_ids) ||
      !setequal(as.character(sample_ids), as.character(label_ids))) {
    stop(
      "Result feature sample names do not match the labels vector. Check ",
      "sample-name standardization consistency between the input artifact ",
      "(h5ad/feather/pseudobulk) and the obs labels."
    )
  }

  reorder <- match(label_ids, sample_ids)
  if (anyNA(reorder)) {
    stop("Result sample identifiers cannot be aligned to labels.")
  }
  if (!identical(reorder, seq_along(reorder))) {
    square_feature_matrix <- !is.null(colnames(feat_mat)) &&
      ncol(feat_mat) == length(sample_ids) &&
      !anyNA(colnames(feat_mat)) &&
      !anyDuplicated(colnames(feat_mat)) &&
      setequal(as.character(colnames(feat_mat)), as.character(sample_ids))
    if (square_feature_matrix) {
      feat_mat <- feat_mat[reorder, reorder, drop = FALSE]
    } else {
      feat_mat <- feat_mat[reorder, , drop = FALSE]
    }
    if (!is.null(dist_mat)) {
      dist_matrix <- as.matrix(dist_mat)
      if (nrow(dist_matrix) != length(sample_ids) ||
          ncol(dist_matrix) != length(sample_ids)) {
        stop("Result distance matrix dimensions do not match feature rows.")
      }
      dist_matrix <- dist_matrix[reorder, reorder, drop = FALSE]
      dimnames(dist_matrix) <- list(label_ids, label_ids)
      dist_mat <- dist_matrix
    }
  }
  names(labels) <- label_ids
  list(feat_mat = feat_mat, labels = labels, dist_mat = dist_mat)
}

# Create result bundle with scores, feature matrix, distance matrix, labels
create_result_bundle <- function(
  feat_mat,
  labels,
  dist_mat = NULL,
  extra = list()
) {
  aligned <- align_result_samples(feat_mat, labels, dist_mat = dist_mat)
  feat_mat <- aligned$feat_mat
  labels <- aligned$labels
  dist_mat <- aligned$dist_mat
  if (is.null(dist_mat)) {
    dist_mat <- dist(feat_mat)
  } else if (is.matrix(dist_mat)) {
    dist_mat <- as.dist(dist_mat)
  }
  result <- list(
    scores = calc_sep_score(dist_mat, labels),
    feat_mat = feat_mat,
    dist_mat = dist_mat,
    labels = labels
  )
  result <- c(result, extra)

  return(result)
}
