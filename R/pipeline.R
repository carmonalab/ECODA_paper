# ============================================================
# PIPELINE FUNCTIONS
# ============================================================

datrans <- function(
  count_mat,
  labels = NULL,
  Amount_of_perturbation,
  n_ct_to_select,
  cts = NULL,
  reps = 20,
  trans_method = c(
    "counts",
    "freq",
    "arcsine_sqrt",
    "alr_randref",
    "alr_mincvref",
    "clr"
  ),
  zero_imp_method = "counts_all__1",
  n_cores = 8
) {
  colnames(count_mat) <- make.names(colnames(count_mat), unique = TRUE)
  n_half_samples <- round(dim(count_mat)[1] / 2)
  if (!is.null(cts)) {
    n_ct_to_select[n_ct_to_select >= length(cts)] <- length(cts)
  }
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  rets <- foreach(
    da = Amount_of_perturbation,
    .export = c(
      "calc_perc_df",
      "impute_zeros",
      "calc_sil",
      "calc_modularity",
      "compute_snn_graph",
      "compute_KNN_from_dist",
      "clust_eval",
      "adjustedRandIndex",
      "cv",
      "clr"
    ),
    .packages = c("dplyr", "mclust", "robCompositions", "zCompositions"),
    .errorhandling = "pass",
    .combine = rbind
  ) %dopar%
    {
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
      for (nct in n_ct_to_select) {
        for (rep in 1:reps) {
          df_counts_temp <- count_mat
          if (!is.null(cts)) {
            ct_da <- sample(cts, size = nct)
          } else {
            ct_da <- sample(colnames(df_counts_temp), size = nct)
          }
          if (is.null(labels)) {
            half_samples_da <- sample(
              row.names(df_counts_temp),
              size = n_half_samples
            )
            labels_random <- as.numeric(
              row.names(count_mat) %in% half_samples_da
            )
            rsums_before <- rowSums(df_counts_temp)
            df_counts_temp[half_samples_da, ct_da] <- round(
              df_counts_temp[half_samples_da, ct_da] * da
            )
            rsums_after <- rowSums(df_counts_temp)
            df_counts_temp <- round(
              df_counts_temp / (rsums_after / rsums_before)
            )
          }
          df_counts_temp <- df_counts_temp %>%
            select_if(colSums(.) != 0) %>%
            mutate_all(as.numeric)
          for (zmet in zero_imp_method) {
            df_freq <- df_counts_temp %>% calc_perc_df()
            df_arcsine_sqrt <- asin(sqrt(df_freq / 100))
            if (grepl("percentage|counts|multRepl", zmet)) {
              zero_imp_method_split <- strsplit(zmet, "__")
              if (grepl("percentage|counts", zmet)) {
                df_counts_temp_imputed <- df_counts_temp %>%
                  impute_zeros(
                    clr_zero_impute_method = zero_imp_method_split[[1]][1],
                    clr_zero_impute_num = as.numeric(eval(parse(
                      text = zero_imp_method_split[[1]][2]
                    )))
                  )
                df_freq_imputed <- df_counts_temp_imputed %>% calc_perc_df()
              } else if (grepl("multRepl", zmet)) {
                df_freq_imputed <- df_freq %>%
                  zCompositions::multRepl(
                    label = 0,
                    dl = rep(
                      as.numeric(zero_imp_method_split[[1]][2]),
                      ncol(df_freq)
                    ),
                    z.warning = 1,
                    frac = 1
                  )
              }
            } else if (zmet == "multLN") {
              df_freq_imputed <- df_freq %>%
                zCompositions::multLN(
                  label = 0,
                  dl = rep(0.1, ncol(df_freq)),
                  z.warning = 0.9
                )
            }
            for (met in trans_method) {
              if (grepl("counts", met)) {
                df <- df_counts_temp
              } else if (grepl("freq", met)) {
                df <- df_freq
              } else if (grepl("arcsine_sqrt", met)) {
                df <- df_arcsine_sqrt
              } else if (grepl("alr_mincvref", met)) {
                ct_ref <- sample(
                  colnames(df_freq_imputed)[
                    !colnames(df_freq_imputed) %in% ct_da
                  ],
                  size = 1
                )
                df <- Hotelling::alr(
                  as.formula(paste0(ct_ref, "~.")),
                  df_freq_imputed
                )
              } else if (grepl("alr_randref", met)) {
                cvs <- apply(df_freq_imputed, 2, cv)
                ct_ref_mincv <- colnames(df_freq_imputed)[which(
                  cvs == min(cvs[!colnames(df_freq_imputed) %in% ct_da])
                )][1]
                df <- Hotelling::alr(
                  as.formula(paste0(ct_ref_mincv, "~.")),
                  df_freq_imputed
                )
              } else if (grepl("clr", met)) {
                df <- clr(df_freq_imputed)
              }
              dist_mat <- dist(df)
              if (is.null(labels)) {
                avg_sil <- calc_sil(dist_mat, labels_random)
                anosim_score <- vegan::anosim(
                  x = dist_mat,
                  grouping = labels_random,
                  distance = "euclidean",
                  permutations = 99
                )[["statistic"]]
                mod <- calc_modularity(dist_mat, labels_random)
                cluster_score <- clust_eval(dist_mat, labels_random)
              } else {
                avg_sil <- calc_sil(dist_mat, labels)
                anosim_score <- vegan::anosim(
                  x = dist_mat,
                  grouping = labels,
                  distance = "euclidean",
                  permutations = 99
                )[["statistic"]]
                mod <- calc_modularity(dist_mat, labels)
                cluster_score <- clust_eval(dist_mat, labels)
              }
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
              new_row_df$dist_mat <- list(df)
              res <- rbind(res, new_row_df)
            }
          }
        }
      }
      return(res)
    }
  stopCluster(cluster)
  return(rets)
}

run_analyses <- function(result_list, ds, seurat, path_data, path_plots) {
  result_list[["bmark"]][[ds]] <- run_benchmark_analysis(
    res_list = result_list[["bmark"]][[ds]],
    ds = ds,
    seurat = seurat,
    path_data = path_data,
    path_plots = path_plots
  )
  labels <- get_labels(seurat, seurat@misc$label_col)
  ct_comps <- get_ct_comp_df_seurat(
    seurat,
    sample_col = "Sample",
    ct_col = seurat@misc$hi_res_ct_col
  )
  if (is.null(result_list[["trans"]][[ds]])) {
    result_list[["trans"]][[ds]] <- run_transformation_analysis(
      ct_comps,
      labels
    )
  }
  if (is.null(result_list[["zeroimp"]][[ds]])) {
    result_list[["zeroimp"]][[ds]] <- run_zeroimp_analysis(ct_comps, labels)
  }
  saveRDS(result_list, file = "result_list.rds")
  return(result_list)
}

run_benchmark_analysis <- function(
  res_list,
  ds,
  seurat,
  sample_col = "Sample",
  factors_test = c(2, 3, 5, 10, 15),
  path_data,
  path_plots,
  seurat_res = c(0.1, 0.4, 2, 5, 20),
  HVGs = c(1000, 2000, 3000),
  ECODA_top_varexp_hvct = seq(0, 0.9, 0.1),
  gloscope_n_pca_dims = c(10, 30, 50)
) {
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
      scpoli_dims <- 15
    }

    file_mrvi <- file.path(
      path_data,
      paste0(ds_filename, "_hvg", i, "_mrvi_dists.feather")
    )
    file_pilot <- file.path(
      path_data,
      paste0(ds_filename, "_hvg", i, "_highres_pilot_dists.feather")
    )
    files_scpoli <- file.path(
      path_data,
      paste0(
        ds_filename,
        "_hvg",
        i,
        "_highres_scpoli_dims",
        scpoli_dims,
        "_embs.feather"
      )
    )
    file_highres_pilot <- file.path(
      path_data,
      paste0(ds_filename, "_hvg2000_highres_pilot_dists.feather")
    )

    files_to_check <- c(file_mrvi, file_pilot, files_scpoli, file_highres_pilot)
    missing_files <- files_to_check[!file.exists(files_to_check)]

    if (length(missing_files) > 0) {
      stop(
        "The following file(s) are missing:\n",
        paste(missing_files, collapse = "\n")
      )
    }
  }

  # Sample names starting with digits are not allowed in seurat
  seurat@meta.data[[sample_col]] <- standardize_sample_names(seurat@meta.data[[
    sample_col
  ]])
  metadata <- get_metadata(seurat)
  label_col <- seurat@misc$label_col
  labels <- get_labels(seurat, label_col)
  hvg <- get_current_hvgs(seurat)

  if (!"Pseudobulk_schvg2000" %in% names(res_list)) {
    exec_time_pb_norm <- exec_time(
      pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
    )
    res_list[["Pseudobulk_schvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_schvg2000"]] <- process_pseudobulk_fig(
        pb_norm,
        labels
      )
    ) +
      exec_time_pb_norm
  }
  exec_time_pb_norm <- exec_time(
    pb_norm <- get_pb_deseq2(
      seurat,
      sample_col = sample_col,
      hvg = NULL,
      n_hvg = 2000
    )
  )
  if (!"Pseudobulk_hvg2000" %in% names(res_list)) {
    res_list[["Pseudobulk_hvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_hvg2000"]] <- process_pseudobulk_fig(
        pb_norm,
        labels
      )
    ) +
      exec_time_pb_norm
  }

  exec_time_pb_norm_hvg500 <- exec_time(
    pb_norm_hvg500 <- get_pb_deseq2(
      seurat,
      sample_col = sample_col,
      hvg = NULL,
      n_hvg = 500
    )
  )

  if (!"Pseudobulk_hvg500" %in% names(res_list)) {
    res_list[["Pseudobulk_hvg500"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_hvg500"]] <- process_pseudobulk_fig(
        pb_norm_hvg500,
        labels
      )
    ) +
      exec_time_pb_norm_hvg500
  }

  exec_time_pb_norm_bl <- exec_time(
    pb_norm_bl <- get_pb_deseq2(
      seurat,
      sample_col = sample_col,
      hvg = NULL,
      n_hvg = 2000,
      black_list = "default_without_sex_genes"
    )
  )

  if (!"Pseudobulk_hvg2000_bl" %in% names(res_list)) {
    res_list[["Pseudobulk_hvg2000_bl"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_hvg2000_bl"]] <- process_pseudobulk_fig(
        pb_norm_bl,
        labels
      )
    ) +
      exec_time_pb_norm_bl
  }

  if (!"Avg_PCA_embedding" %in% names(res_list)) {
    res_list[["Avg_PCA_embedding"]][["exec_time"]] <- exec_time(
      res_list[["Avg_PCA_embedding"]] <- process_avg_pca_embedding_fig(
        seurat,
        labels
      )
    )
  }

  # Calculate ct pseudobulks
  if (
    !"Pseudobulk_CT_LR_hvg2000" %in% names(res_list) &&
      !is.null(seurat@misc$low_res_ct_col)
  ) {
    res_list[["Pseudobulk_CT_LR_hvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_CT_LR_hvg2000"]] <- process_pseudobulk_ct_fig(
        seurat,
        labels,
        ct_col = seurat@misc$low_res_ct_col,
        sample_col = sample_col,
        hvg = 2000
      )
    )
  }
  if (
    !"Pseudobulk_CT_HR_hvg2000" %in% names(res_list) &&
      !is.null(seurat@misc$hi_res_ct_col)
  ) {
    res_list[["Pseudobulk_CT_HR_hvg2000"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_CT_HR_hvg2000"]] <- process_pseudobulk_ct_fig(
        seurat,
        labels,
        ct_col = seurat@misc$hi_res_ct_col,
        sample_col = sample_col,
        hvg = 2000
      )
    )
  }

  if (
    !"Pseudobulk_CT_LR_hvg500" %in% names(res_list) &&
      !is.null(seurat@misc$hi_res_ct_col)
  ) {
    res_list[["Pseudobulk_CT_LR_hvg500"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_CT_LR_hvg500"]] <- process_pseudobulk_ct_fig(
        seurat,
        labels,
        ct_col = seurat@misc$low_res_ct_col,
        sample_col = sample_col,
        hvg = 500
      )
    )
  }

  if (
    !"Pseudobulk_CT_HR_hvg500" %in% names(res_list) &&
      !is.null(seurat@misc$hi_res_ct_col)
  ) {
    res_list[["Pseudobulk_CT_HR_hvg500"]][["exec_time"]] <- exec_time(
      res_list[["Pseudobulk_CT_HR_hvg500"]] <- process_pseudobulk_ct_fig(
        seurat,
        labels,
        ct_col = seurat@misc$hi_res_ct_col,
        sample_col = sample_col,
        hvg = 500
      )
    )
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
        seurat,
        labels,
        ct_col = seurat@misc$low_res_ct_col
      )
    )
  }

  ## layer2: high res. cell types
  if (!is.null(seurat@misc$hi_res_ct_col)) {
    res_list[["ECODA_authors_HR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = seurat@misc$hi_res_ct_col
      )
    )
    res_list[["ECODA_authors_HR_NULL"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR_NULL"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = seurat@misc$hi_res_ct_col,
        shuffle_labels = TRUE
      )
    )
    res_list[["GloProp"]][["exec_time"]] <- exec_time(
      res_list[["GloProp"]] <- process_gloprop_fig(
        seurat,
        metadata,
        ct_col = seurat@misc$hi_res_ct_col
      )
    )

    for (varexp_hvc in ECODA_top_varexp_hvct) {
      ECODA_authors_HR_top_varexp_hvc <- paste0(
        "ECODA_authors_HR_top_varexp",
        varexp_hvc
      )
      res_list[[ECODA_authors_HR_top_varexp_hvc]][["exec_time"]] <- exec_time(
        res_list[[ECODA_authors_HR_top_varexp_hvc]] <-
          process_coda_fig(
            seurat,
            labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = seurat@misc$hi_res_ct_col
          )
      )

      ECODA_HiTME_HR_layer2_top_varexp_hvc <- paste0(
        "ECODA_HiTME_HR_layer2_top_varexp",
        varexp_hvc
      )
      res_list[[ECODA_HiTME_HR_layer2_top_varexp_hvc]][[
        "exec_time"
      ]] <- exec_time(
        res_list[[ECODA_HiTME_HR_layer2_top_varexp_hvc]] <-
          process_coda_fig(
            seurat,
            labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = "layer2"
          )
      )

      ECODA_HiTME_HR_layer3_top_varexp_hvc <- paste0(
        "ECODA_HiTME_HR_layer3_top_varexp",
        varexp_hvc
      )
      res_list[[ECODA_HiTME_HR_layer3_top_varexp_hvc]][[
        "exec_time"
      ]] <- exec_time(
        res_list[[ECODA_HiTME_HR_layer3_top_varexp_hvc]] <-
          process_coda_fig(
            seurat,
            labels,
            ECODA_top_varexp_hvct = varexp_hvc,
            ct_col = "layer3"
          )
      )
    }

    res_list[["Freq_highres"]][["exec_time"]] <- exec_time(
      res_list[["Freq_highres"]] <- process_coda_fig(
        seurat,
        labels,
        calc_clr = FALSE,
        ct_col = seurat@misc$hi_res_ct_col
      )
    )
  }
  res_list[["ECODA_authors_HR_3most_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3most_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 3,
      var_ct_desc = TRUE,
      ct_col = seurat@misc$hi_res_ct_col
    )
  )
  res_list[["ECODA_authors_HR_2least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_2least_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 2,
      var_ct_desc = FALSE,
      ct_col = seurat@misc$hi_res_ct_col,
    )
  )

  res_list[["ECODA_authors_HR_3least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3least_varcts"]] <-
      process_coda_fig(
        seurat,
        labels,
        ECODA_top_n_hvct = 3,
        var_ct_desc = FALSE,
        ct_col = seurat@misc$hi_res_ct_col,
        title = "ECODA\n2 least var. cell types"
      )
  )
  res_list[["ECODA_HiTME_HR_layer2"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer2"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "layer2"
    )
  )
  res_list[["ECODA_HiTME_HR_layer3"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer3"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "layer3"
    )
  )
  res_list[["ECODA_scATOMIC_HR"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_scATOMIC_HR"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "scATOMIC_pred"
    )
  )

  # Analyze for all resolutions

  ## Ultra high res. cell type clusters based on Leiden clustering to artificially increase the number of cell types (clusters), e.g. to 250 cell types (clusters)
  for (r in seurat_res) {
    res_col_name <- paste0("RNA_snn_res.", r)
    nm <- paste0("ECODA_seuratres_", r)
    res_list[[nm]][["exec_time"]] <- exec_time(
      res_list[[nm]] <- process_coda_fig(seurat, labels, ct_col = res_col_name)
    )
  }

  # Methods that use different number of factors (e.g. PCA or dims)

  ### Required for MOFA to run
  seurat@version <- package_version("3.1.5")
  for (i in factors_test) {
    # Pseudobulk with PCA
    nm1 <- paste0("Pseudobulk_", i, "_PCA_dims")
    if (!nm1 %in% names(res_list)) {
      res_list[[nm1]][["exec_time"]] <- exec_time(
        res_list[[nm1]] <- process_pseudobulk_fig(pb_norm, labels, pca_dims = i)
      ) +
        exec_time_pb_norm
    }
    # Hires CODA with PCA
    nm2 <- paste0("ECODA_authors_HR_", i, "_PCA_dims")
    if (!nm2 %in% names(res_list)) {
      res_list[[nm2]][["exec_time"]] <- exec_time(
        res_list[[nm2]] <- process_coda_fig(
          seurat,
          labels,
          pca_dims = i,
          ct_col = seurat@misc$hi_res_ct_col
        )
      )
    }
    # MOFA
    nm3 <- paste0("MOFA_hvg2000_factors", i)
    if (!nm3 %in% names(res_list)) {
      res_list[[nm3]][["exec_time"]] <- exec_time(
        res_list[[nm3]] <- process_mofa_bulk_fig(
          pb_norm,
          metadata,
          labels,
          num_factors = i
        )
      ) +
        exec_time_pb_norm
    }
    # scITD
    nm4 <- paste0("scITD_hvg2000_factors", i)
    if (!nm4 %in% names(res_list)) {
      res_list[[nm4]][["exec_time"]] <- exec_time(
        res_list[[nm4]] <- process_scitd_fig(
          seurat,
          ct_col = seurat@misc$low_res_ct_col,
          label_col = label_col,
          hvg,
          num_factors = i
        )
      )
    }
  }

  # GloScope
  ## With different numbers of PCA dims
  for (i in gloscope_n_pca_dims) {
    gf <- file.path(
      path_data,
      paste0(ds_filename, "_gloscope_hvg2000_pcadims", i, "_dists.rds")
    )
    nm <- paste0("GloScope_hvg2000_pcadims", i)
    if (!nm %in% names(res_list)) {
      res_list[[nm]][["exec_time"]] <- exec_time(
        res_list[[nm]] <- process_gloscope_fig(
          seurat,
          metadata,
          label_col,
          gf,
          n_pca_dims = i
        )
      )
      res_list[[paste0(nm, "_sqrtmat")]][["exec_time"]] <- exec_time(
        res_list[[paste0(nm, "_sqrtmat")]] <- process_gloscope_sqrtmat_fig(
          metadata,
          label_col,
          gf
        )
      ) +
        res_list[[nm]][["exec_time"]]
    }
  }

  # Methods that use different number of HVGs
  ### Do only for non-default HVGs (default = 2000)
  for (i in HVGs[!HVGs %in% 2000]) {
    # Pseudobulk_schvg_i <- paste0("Pseudobulk_schvg", i)
    Pseudobulk_hvg_i <- paste0("Pseudobulk_hvg", i)
    MOFA_hvg_i_factors15 <- paste0("MOFA_hvg", i, "_factors15")
    scITD_hvg_i_factors5 <- paste0("scITD_hvg", i, "_factors5")
    GloScope_hvg_i_pcadims30 <- paste0("GloScope_hvg", i, "_pcadims30")
    GloScope_hvg_i_pcadims30_sqrtmat <- paste0(
      GloScope_hvg_i_pcadims30,
      "_sqrtmat"
    )
    test_items <- c(
      # Pseudobulk_schvg_i,
      Pseudobulk_hvg_i,
      MOFA_hvg_i_factors15,
      scITD_hvg_i_factors5,
      GloScope_hvg_i_pcadims30,
      GloScope_hvg_i_pcadims30_sqrtmat
    )

    if (any(!test_items %in% names(res_list))) {
      # Memory critical steps
      gc()
      misc <- seurat@misc
      seurat <- create_clean_seuratv5_object(seurat)
      seurat@misc <- misc
      gc()
      seurat <- NormalizeData(seurat)
      gc()
      seurat <- FindVariableFeatures(seurat, nfeatures = i)
      hvg <- get_current_hvgs(seurat)
      gc()
      if (
        any(
          !c(GloScope_hvg_i_pcadims30, GloScope_hvg_i_pcadims30_sqrtmat) %in%
            names(res_list)
        )
      ) {
        seurat <- ScaleData(seurat)
        gc()
        # Needs a lot of memory for 3000 HVGs
        seurat <- RunPCA(seurat, dims = 1:50, verbose = FALSE)
        gc()
      }

      if (!Pseudobulk_hvg_i %in% names(res_list)) {
        exec_time_pb_norm <- exec_time(
          pb_norm <- get_pb_deseq2(seurat, sample_col = sample_col, hvg = hvg)
        )
        res_list[[Pseudobulk_hvg_i]][["exec_time"]] <- exec_time(
          res_list[[Pseudobulk_hvg_i]] <- process_pseudobulk_fig(
            pb_norm,
            labels
          )
        ) +
          exec_time_pb_norm
      }

      test_items <- c(
        Pseudobulk_hvg_i,
        MOFA_hvg_i_factors15
      )
      if (any(!test_items %in% names(res_list))) {
        exec_time_pb_norm <- exec_time(
          pb_norm <- get_pb_deseq2(
            seurat,
            sample_col = sample_col,
            hvg = NULL,
            n_hvg = i
          )
        )
        res_list[[Pseudobulk_hvg_i]][["exec_time"]] <- exec_time(
          res_list[[Pseudobulk_hvg_i]] <- process_pseudobulk_fig(
            pb_norm,
            labels
          )
        ) +
          exec_time_pb_norm

        res_list[[MOFA_hvg_i_factors15]][["exec_time"]] <- exec_time(
          res_list[[MOFA_hvg_i_factors15]] <-
            process_mofa_bulk_fig(
              pb_norm,
              metadata = metadata,
              labels,
              num_factors = 15
            )
        ) +
          exec_time_pb_norm
      }

      if (!scITD_hvg_i_factors5 %in% names(res_list)) {
        res_list[[scITD_hvg_i_factors5]][["exec_time"]] <- exec_time(
          res_list[[scITD_hvg_i_factors5]] <- process_scitd_fig(
            seurat,
            ct_col = seurat@misc$low_res_ct_col,
            label_col = label_col,
            hvg,
            num_factors = 5
          )
        )
      }

      test_items <- c(
        GloScope_hvg_i_pcadims30,
        GloScope_hvg_i_pcadims30_sqrtmat
      )
      if (any(!test_items %in% names(res_list))) {
        gloscope_dist_file <- file.path(
          path_data,
          paste0(ds_filename, "_gloscope_hvg", i, "_pcadims30_dists.rds")
        )
        res_list[[GloScope_hvg_i_pcadims30]][["exec_time"]] <- exec_time(
          res_list[[GloScope_hvg_i_pcadims30]] <-
            process_gloscope_fig(
              seurat,
              metadata,
              label_col,
              gloscope_dist_file = gloscope_dist_file
            )
        )
        res_list[[GloScope_hvg_i_pcadims30_sqrtmat]][[
          "exec_time"
        ]] <- exec_time(
          res_list[[GloScope_hvg_i_pcadims30_sqrtmat]] <-
            process_gloscope_sqrtmat_fig(
              metadata,
              label_col,
              gloscope_dist_file = gloscope_dist_file
            )
        ) +
          res_list[[GloScope_hvg_i_pcadims30]][["exec_time"]]
      }
    }
  }

  for (i in HVGs) {
    # --- MrVI (Runs once per HVG) ---
    mrvi_dist_file <- file.path(
      path_data,
      paste0(ds_filename, "_hvg", i, "_mrvi_dists.feather")
    )
    res_list[[paste0("MrVI_hvg", i)]] <- process_mrvi_fig(
      mrvi_dist_file = mrvi_dist_file,
      labels
    )

    # --- PILOT (Runs once per HVG) ---
    pilot_dist_file <- file.path(
      path_data,
      paste0(ds_filename, "_hvg", i, "_highres_pilot_dists.feather")
    )
    res_list[[paste0("PILOT_hvg", i)]] <- process_pilot_fig(
      pilot_dist_file = pilot_dist_file,
      labels
    )

    scpoli_emb_file <- file.path(
      path_data,
      paste0(ds_filename, "_hvg", i, "_highres_scpoli_dims15_embs.feather")
    )
    res_list[[paste0(
      "scPoli_hvg",
      i,
      "_dims15_highres"
    )]] <- process_scpoli_fig(scpoli_emb_file = scpoli_emb_file, labels)

    # --- scPoli (Runs once OR multiple times depending on HVG) ---
    if (i == 2000) {
      target_dims <- factors_test
      for (f in target_dims) {
        scpoli_emb_file <- file.path(
          path_data,
          paste0(
            ds_filename,
            "_hvg",
            i,
            "_highres_scpoli_dims",
            f,
            "_embs.feather"
          )
        )
        res_list[[paste0(
          "scPoli_hvg",
          i,
          "_dims",
          f,
          "_highres"
        )]] <- process_scpoli_fig(scpoli_emb_file = scpoli_emb_file, labels)
      }

      pilot_dist_file <- file.path(
        path_data,
        paste0(ds_filename, "_hvg", i, "_lowres_pilot_dists.feather")
      )
      res_list[[paste0("PILOT_hvg", i, "_lowres")]] <- process_pilot_fig(
        pilot_dist_file = pilot_dist_file,
        labels
      )

      scpoli_emb_file <- file.path(
        path_data,
        paste0(ds_filename, "_hvg", i, "_lowres_scpoli_dims15_embs.feather")
      )
      res_list[[paste0(
        "scPoli_hvg",
        i,
        "_dims15_lowres"
      )]] <- process_scpoli_fig(scpoli_emb_file = scpoli_emb_file, labels)
    }
  }

  # # ECODA + PB
  # for (i in c(0, 0.25, 0.5, 0.75, 1)) {
  #   for (norm in c("max", "median", "zscore", "quantile")) {
  #     res_list[[paste0("ECODA_PB_combo_norm", norm, "_ecodaweight", i)]] <- process_ecodapb_fig(
  #       dist_mat_ecoda = res_list[["ECODA_authors_HR"]][["dist_mat"]],
  #       dist_mat_pb = res_list[["Pseudobulk_hvg2000"]][["dist_mat"]],
  #       norm_method = norm,
  #       ecoda_weight = i,
  #       labels = res_list[["ECODA_authors_HR"]][["labels"]],
  #     )
  #     res_list[[paste0("ECODA_PB_combo_norm", norm, "_ecodaweight", i)]][["exec_time"]] <-
  #       res_list[["ECODA_authors_HR"]][["exec_time"]] + res_list[["Pseudobulk_hvg2000"]][["exec_time"]]
  #   }
  # }

  return(res_list)
}


run_transformation_analysis <- function(ct_comps, labels) {
  res_list <- datrans(
    ct_comps,
    labels,
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
    dplyr::group_by(.data$trans_method) %>%
    summarize(
      ANOSIM_score = mean(ANOSIM_score),
      Modularity_score = mean(Modularity_score),
      Adjusted_Rand_Index = mean(Adjusted_Rand_Index)
    ) %>%
    ungroup()
  return(res_list)
}

run_zeroimp_analysis <- function(ct_comps, labels) {
  df <- ct_comps %>% select_if(colSums(.) != 0) %>% mutate_all(as.numeric)
  res_list <- list()
  res_list[["counts_zeros_2/3min"]] <- df %>%
    impute_zeros("counts_zeros", 2 / 3) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["counts_zeros_1"]] <- df %>%
    impute_zeros("counts_zeros", 1) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["counts_all_1"]] <- df %>%
    impute_zeros("counts_all", 1) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["counts_zeros_0.5"]] <- df %>%
    impute_zeros("counts_zeros", 0.5) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["counts_all_0.5"]] <- df %>%
    impute_zeros("counts_all", 0.5) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.001%"]] <- df %>%
    impute_zeros("percentage_all", 0.001) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.01%"]] <- df %>%
    impute_zeros("percentage_all", 0.01) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["perc_all_0.1%"]] <- df %>%
    impute_zeros("percentage_all", 0.1) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["perc_all_1%"]] <- df %>%
    impute_zeros("percentage_all", 1) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["asinsqrt"]] <- df %>%
    calc_perc_df() %>%
    mutate(across(everything(), ~ . / 100)) %>%
    sqrt() %>%
    asin() %>%
    dist() %>%
    calc_sep_score(labels)
  df_multLN <- df %>%
    calc_perc_df() %>%
    zCompositions::multLN(label = 0, dl = rep(0.1, ncol(df)), z.warning = 0.9)
  res_list[["multLN"]] <- df_multLN %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels[row.names(df) %in% row.names(df_multLN)])
  res_list[["multRepl_0.01%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(
      label = 0,
      dl = rep(0.01, ncol(df)),
      z.warning = 1,
      frac = 1
    ) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["multRepl_0.1%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(
      label = 0,
      dl = rep(0.1, ncol(df)),
      z.warning = 1,
      frac = 1
    ) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  res_list[["multRepl_1%"]] <- df %>%
    calc_perc_df() %>%
    zCompositions::multRepl(
      label = 0,
      dl = rep(1, ncol(df)),
      z.warning = 1,
      frac = 1
    ) %>%
    clr() %>%
    dist() %>%
    calc_sep_score(labels)
  return(res_list)
}


# Timing wrapper
exec_time <- function(fun) {
  start_time <- Sys.time()
  fun
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  time_taken
  return(time_taken)
}
