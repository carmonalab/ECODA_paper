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

# Read the checksums.md5 sidecar written by the benchmark submitters (paths
# relative to the benchmark results root, "<md5>  <path>" lines, GNU md5sum
# format). Returns a named character vector (path -> hash), or NULL when no
# sidecar exists (e.g. results that predate the sidecar).
read_md5_sidecar <- function(checksum_file) {
  if (!file.exists(checksum_file)) return(NULL)
  lines <- readLines(checksum_file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  parts <- strsplit(lines, "  ", fixed = TRUE)
  hashes <- vapply(parts, `[`, character(1), 1)
  paths <- vapply(parts, `[`, character(1), 2)
  names(hashes) <- paths
  return(hashes)
}

# Best-effort integrity check before readRDS of NAS-produced bundles (RDS
# deserialization can execute code): files listed in the sidecar must match
# their md5; files not listed (legacy results) are read unverified.
verify_md5_sidecar <- function(file_path, checksums) {
  if (is.null(checksums) || length(checksums) == 0) {
    return(invisible(NULL))
  }
  key <- file.path(basename(dirname(file_path)), basename(file_path))
  expected <- unname(checksums[key])
  if (is.na(expected)) return(invisible(NULL))
  actual <- unname(tools::md5sum(file_path))
  if (is.na(actual) || actual != expected) {
    stop("Checksum mismatch for ", file_path, " (expected ", expected,
         ", got ", actual, "). The HPC result bundle was modified or ",
         "corrupted; re-run the benchmark pipeline (or remove the file) ",
         "before loading it.")
  }
  return(invisible(NULL))
}

# Load HPC-computed benchmark results (Pipeline A methods + Pipeline B
# trans/zeroimp) into the notebook's result_list, preserving the legacy rerun
# semantics: entries already present in result_list$bmark[[ds]] are kept; only
# missing methods are read from ${path_results_nas}/<ds>_<method>.rds (warning
# on missing files). Bundles are verified against the checksums.md5 sidecar
# (written by the submit scripts) before deserialization. Also reads
# <ds>_trans.rds / <ds>_zeroimp.rds into result_list$trans[[ds]] /
# result_list$zeroimp[[ds]] (only if not yet present). Saves result_list.rds
# and returns the list.
load_hpc_benchmark_results <- function(
  result_list,
  ds,
  path_results_nas,
  methods = c("gloscope", "mofa", "pseudobulk", "scitd")
) {
  if (is.null(result_list[["bmark"]])) result_list[["bmark"]] <- list()
  if (is.null(result_list[["trans"]])) result_list[["trans"]] <- list()
  if (is.null(result_list[["zeroimp"]])) result_list[["zeroimp"]] <- list()
  if (is.null(result_list[["bmark"]][[ds]])) {
    result_list[["bmark"]][[ds]] <- list()
  }

  checksums <- read_md5_sidecar(file.path(path_results_nas, "checksums.md5"))

  for (method in methods) {
    method_file <- file.path(
      path_results_nas,
      paste0(ds, "_", method, ".rds")
    )
    if (!file.exists(method_file)) {
      warning("HPC benchmark result file not found: ", method_file)
      next
    }
    verify_md5_sidecar(method_file, checksums)
    bundles <- readRDS(method_file)
    for (nm in names(bundles)) {
      if (!nm %in% names(result_list[["bmark"]][[ds]])) {
        result_list[["bmark"]][[ds]][[nm]] <- bundles[[nm]]
      }
    }
  }

  trans_file <- file.path(path_results_nas, paste0(ds, "_trans.rds"))
  if (is.null(result_list[["trans"]][[ds]])) {
    if (file.exists(trans_file)) {
      verify_md5_sidecar(trans_file, checksums)
      result_list[["trans"]][[ds]] <- readRDS(trans_file)
    } else {
      warning("HPC transformation result file not found: ", trans_file)
    }
  }

  zeroimp_file <- file.path(path_results_nas, paste0(ds, "_zeroimp.rds"))
  if (is.null(result_list[["zeroimp"]][[ds]])) {
    if (file.exists(zeroimp_file)) {
      verify_md5_sidecar(zeroimp_file, checksums)
      result_list[["zeroimp"]][[ds]] <- readRDS(zeroimp_file)
    } else {
      warning("HPC zero-imputation result file not found: ", zeroimp_file)
    }
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
  seurat_res = c(0.1, 0.4, 2, 5, 20),
  HVGs = c(1000, 2000, 3000),
  ECODA_top_varexp_hvct = seq(0, 0.9, 0.1)
) {
  # Files preprocessed with python. Feather names always use the datasets.json
  # key (both the HPC Python pipeline and the legacy qmd write them that way),
  # so the legacy GongSharma -> GongSharma_all remap is dropped.
  for (i in HVGs) {
    if (i == 2000) {
      scpoli_dims <- factors_test
    } else {
      scpoli_dims <- 15
    }

    file_mrvi <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_mrvi_dists.feather")
    )
    file_pilot <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_highres_pilot_dists.feather")
    )
    files_scpoli <- file.path(
      path_data,
      paste0(
        ds,
        "_hvg",
        i,
        "_highres_scpoli_dims",
        scpoli_dims,
        "_embs.feather"
      )
    )
    file_highres_pilot <- file.path(
      path_data,
      paste0(ds, "_hvg2000_highres_pilot_dists.feather")
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

  # Pseudobulk is computed by the HPC pipeline (run_pseudobulk_hpc). A local
  # DESeq2 pseudobulk is still needed by ECODA_deconv below; as in the legacy
  # code it is NOT wrapped in exec_time.
  pb_norm <- get_pb_deseq2(
    seurat,
    sample_col = sample_col,
    hvg = NULL,
    n_hvg = 2000
  )

  if (!"Avg_PCA_embedding" %in% names(res_list)) {
    res_list[["Avg_PCA_embedding"]][["exec_time"]] <- exec_time(
      res_list[["Avg_PCA_embedding"]] <- process_avg_pca_embedding_fig(
        seurat,
        labels
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
  if (!is.null(seurat@misc$cell_type_low_res)) {
    res_list[["ECODA_authors_LR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_LR"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = seurat@misc$cell_type_low_res
      )
    )
  }

  ## layer2: high res. cell types
  if (!is.null(seurat@misc$cell_type_high_res)) {
    res_list[["ECODA_authors_HR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = seurat@misc$cell_type_high_res
      )
    )
    res_list[["ECODA_authors_HR_NULL"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR_NULL"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = seurat@misc$cell_type_high_res,
        shuffle_labels = TRUE
      )
    )
    res_list[["GloProp"]][["exec_time"]] <- exec_time(
      res_list[["GloProp"]] <- process_gloprop_fig(
        seurat,
        metadata,
        ct_col = seurat@misc$cell_type_high_res
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
            ct_col = seurat@misc$cell_type_high_res
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
        ct_col = seurat@misc$cell_type_high_res
      )
    )
  }
  res_list[["ECODA_authors_HR_3most_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3most_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 3,
      var_ct_desc = TRUE,
      ct_col = seurat@misc$cell_type_high_res
    )
  )
  res_list[["ECODA_authors_HR_2least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_2least_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 2,
      var_ct_desc = FALSE,
      ct_col = seurat@misc$cell_type_high_res,
    )
  )

  res_list[["ECODA_authors_HR_3least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3least_varcts"]] <-
      process_coda_fig(
        seurat,
        labels,
        ECODA_top_n_hvct = 3,
        var_ct_desc = FALSE,
        ct_col = seurat@misc$cell_type_high_res,
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

  # Pseudobulk with PCA, MOFA and scITD moved to the HPC pipeline
  # (run_pseudobulk_hpc / run_mofa_hpc / run_scitd_hpc).
  for (i in factors_test) {
    # Hires CODA with PCA
    nm2 <- paste0("ECODA_authors_HR_", i, "_PCA_dims")
    if (!nm2 %in% names(res_list)) {
      res_list[[nm2]][["exec_time"]] <- exec_time(
        res_list[[nm2]] <- process_coda_fig(
          seurat,
          labels,
          pca_dims = i,
          ct_col = seurat@misc$cell_type_high_res
        )
      )
    }
  }

  # Methods that use different number of HVGs
  ### Pseudobulk/MOFA/scITD/GloScope for hvg1000/3000 moved to the HPC
  ### pipeline (run_pseudobulk_hpc / run_mofa_hpc / run_scitd_hpc /
  ### run_gloscope_hpc).
  for (i in HVGs) {
    # --- MrVI (Runs once per HVG) ---
    mrvi_dist_file <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_mrvi_dists.feather")
    )
    res_list[[paste0("MrVI_hvg", i)]] <- process_mrvi_fig(
      mrvi_dist_file = mrvi_dist_file,
      labels
    )

    # --- PILOT (Runs once per HVG) ---
    pilot_dist_file <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_highres_pilot_dists.feather")
    )
    res_list[[paste0("PILOT_hvg", i)]] <- process_pilot_fig(
      pilot_dist_file = pilot_dist_file,
      labels
    )

    scpoli_emb_file <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_highres_scpoli_dims15_embs.feather")
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
            ds,
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
        paste0(ds, "_hvg", i, "_lowres_pilot_dists.feather")
      )
      res_list[[paste0("PILOT_hvg", i, "_lowres")]] <- process_pilot_fig(
        pilot_dist_file = pilot_dist_file,
        labels
      )

      scpoli_emb_file <- file.path(
        path_data,
        paste0(ds, "_hvg", i, "_lowres_scpoli_dims15_embs.feather")
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

# ============================================================
# HPC DRIVERS (called by src/5_run_benchmark_methods/
# run_r_sample_embedding_methods/ workers; not used by the notebook)
#
# Each driver computes its combos, times every combo with exec_time(),
# appends the numeric-seconds exec_time to each result bundle, handles the
# per-combo cache files (<ds>_<combo>.rds, skip-if-exists unless
# --force; combo names are method-prefixed, so no method infix) and writes
# per-combo exec-log rows. Returns a named list of
# result bundles (legacy result names, minus the GloScope _sqrtmat suffix).
# ============================================================

# Precompute the shared DESeq2 pseudobulks used by MOFA and Pseudobulk.
# Returns per-variant list(pb, time_secs). Variants (legacy combo list):
#   schvg2000  get_pb_deseq2(hvg = top-2000 hvg_rank genes)
#   hvg2000    get_pb_deseq2(n_hvg = 2000)
#   hvg500     get_pb_deseq2(n_hvg = 500)
#   hvg2000_bl get_pb_deseq2(n_hvg = 2000, black_list = "default_without_sex_genes")
#   hvg1000    get_pb_deseq2(n_hvg = 1000)
#   hvg3000    get_pb_deseq2(n_hvg = 3000)
# `variants` restricts the computed subset (default: all PB_VARIANT_NAMES);
# failure-resume callers pass only the missing set so existing caches are
# neither recomputed nor overwritten.
# NOTE: the legacy get_pb_deseq2 black-list behavior (incl. the pre-existing
# `%in% black_list` no-op typo at pseudobulk.R:84) is preserved on purpose,
# to keep Pseudobulk_hvg2000_bl legacy-equivalent.
prepare_pseudobulks_hpc <- function(
  seurat,
  sample_col = "Sample",
  hvg_rank_genes = NULL,
  variants = PB_VARIANT_NAMES
) {
  n_ranked <- length(hvg_rank_genes)
  pb_specs <- list(
    schvg2000 = list(
      hvg = hvg_rank_genes[seq_len(min(2000, n_ranked))],
      n_hvg = NULL,
      black_list = "none"
    ),
    hvg2000 = list(hvg = NULL, n_hvg = 2000, black_list = "none"),
    hvg500 = list(hvg = NULL, n_hvg = 500, black_list = "none"),
    hvg2000_bl = list(
      hvg = NULL,
      n_hvg = 2000,
      black_list = "default_without_sex_genes"
    ),
    hvg1000 = list(hvg = NULL, n_hvg = 1000, black_list = "none"),
    hvg3000 = list(hvg = NULL, n_hvg = 3000, black_list = "none")
  )

  variants <- intersect(variants, names(pb_specs))
  if (length(variants) == 0) {
    stop("No valid pseudobulk variant names requested (valid: ",
         paste(names(pb_specs), collapse = ", "), ")")
  }

  results <- list()
  for (variant in variants) {
    spec <- pb_specs[[variant]]
    time_secs <- exec_time(
      pb <- get_pb_deseq2(
        seurat,
        sample_col = sample_col,
        hvg = spec$hvg,
        n_hvg = spec$n_hvg,
        black_list = spec$black_list
      )
    )
    results[[variant]] <- list(
      pb = pb,
      time_secs = as.numeric(time_secs, units = "secs")
    )
  }
  return(results)
}

# GloScope combos: hvg2000 x pcadims {10,30,50}; hvg1000, hvg3000 x pcadims 30.
# Raw GloScope distances are cached at <gloscope_cache_dir>/<ds>_gloscope_hvg<n>
# _pcadims<d>_dists.rds (sqrt + NA->0 applied by process_gloscope_fig); on a
# cache miss the combo time includes the distance computation, on a hit it is
# sqrt + read only — matching the legacy path_data dist-cache semantics.
run_gloscope_hpc <- function(
  seurat,
  metadata,
  label_col,
  sample_col = "Sample",
  gloscope_cache_dir,
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL
) {
  combos <- list(
    list(hvg = 2000, pcadims = 10),
    list(hvg = 2000, pcadims = 30),
    list(hvg = 2000, pcadims = 50),
    list(hvg = 1000, pcadims = 30),
    list(hvg = 3000, pcadims = 30)
  )

  results <- list()
  for (combo in combos) {
    n_hvg <- combo$hvg
    n_pca_dims <- combo$pcadims
    nm <- paste0("GloScope_hvg", n_hvg, "_pcadims", n_pca_dims)
    bundle_file <- file.path(
      results_dir,
      paste0(ds, "_", nm, ".rds")
    )
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }

    emb_key <- paste0("pca_benchmark_analysis_hvg", n_hvg)
    if (!emb_key %in% names(seurat@reductions)) {
      warning("Embedding '", emb_key, "' not found in seurat; skipping ", nm)
      next
    }

    dist_file <- file.path(
      gloscope_cache_dir,
      paste0(ds, "_gloscope_hvg", n_hvg, "_pcadims", n_pca_dims, "_dists.rds")
    )
    time_secs <- exec_time(
      res <- process_gloscope_fig(
        embedding_matrix = seurat@reductions[[emb_key]]@cell.embeddings,
        sample_ids = seurat@meta.data[[sample_col]],
        metadata = metadata,
        label_col = label_col,
        gloscope_dist_file = dist_file,
        n_pca_dims = n_pca_dims
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }
  return(results)
}

# MOFA combos: MOFA_hvg2000_factors{2,3,5,10,15}, MOFA_hvg{1000,3000}_factors15.
# Each combo time = pb creation time (pb_variants$time_secs) + MOFA runtime.
# Combos with num_factors >= n_samples are skipped with a warning (the model
# cannot fit that many factors; makes the 5-sample _debug dataset usable).
run_mofa_hpc <- function(
  metadata,
  labels,
  pb_variants,
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL
) {
  combos <- c(
    paste0("MOFA_hvg2000_factors", c(2, 3, 5, 10, 15)),
    "MOFA_hvg1000_factors15",
    "MOFA_hvg3000_factors15"
  )
  n_samples <- length(labels)

  results <- list()
  for (nm in combos) {
    n_hvg <- as.integer(sub("MOFA_hvg(\\d+)_factors.*", "\\1", nm))
    num_factors <- as.integer(sub(".*factors", "", nm))
    if (num_factors >= n_samples) {
      warning(nm, " skipped: num_factors (", num_factors,
              ") >= n_samples (", n_samples, ")")
      next
    }
    pb_variant <- pb_variants[[paste0("hvg", n_hvg)]]
    if (is.null(pb_variant)) {
      stop("Pseudobulk variant 'hvg", n_hvg, "' missing for ", nm)
    }

    bundle_file <- file.path(
      results_dir,
      paste0(ds, "_", nm, ".rds")
    )
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }

    time_secs <- exec_time(
      res <- process_mofa_bulk_fig(
        pb_variant$pb,
        metadata = metadata,
        labels,
        num_factors = num_factors
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs") +
      pb_variant$time_secs
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }
  return(results)
}

# Pseudobulk combos: Pseudobulk_schvg2000, Pseudobulk_hvg2000,
# Pseudobulk_hvg500, Pseudobulk_hvg2000_bl, Pseudobulk_CT_LR_hvg2000,
# Pseudobulk_CT_HR_hvg2000, Pseudobulk_CT_LR_hvg500, Pseudobulk_CT_HR_hvg500,
# Pseudobulk_{2,3,5,10,15}_PCA_dims, Pseudobulk_hvg1000, Pseudobulk_hvg3000.
# Plain variants reuse the precomputed pb_variants (time = pb creation +
# processing); CT variants compute their per-cell-type pseudobulk internally
# (each gated on its own ct col being non-null).
run_pseudobulk_hpc <- function(
  seurat,
  labels,
  pb_variants,
  sample_col = "Sample",
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL
) {
  plain_combos <- list(
    Pseudobulk_schvg2000 = "schvg2000",
    Pseudobulk_hvg2000 = "hvg2000",
    Pseudobulk_hvg500 = "hvg500",
    Pseudobulk_hvg2000_bl = "hvg2000_bl",
    Pseudobulk_hvg1000 = "hvg1000",
    Pseudobulk_hvg3000 = "hvg3000"
  )
  pca_combos <- paste0("Pseudobulk_", c(2, 3, 5, 10, 15), "_PCA_dims")
  ct_combos <- list(
    Pseudobulk_CT_LR_hvg2000 = list(
      ct_col = seurat@misc$cell_type_low_res,
      hvg = 2000
    ),
    Pseudobulk_CT_HR_hvg2000 = list(
      ct_col = seurat@misc$cell_type_high_res,
      hvg = 2000
    ),
    Pseudobulk_CT_LR_hvg500 = list(
      ct_col = seurat@misc$cell_type_low_res,
      hvg = 500
    ),
    Pseudobulk_CT_HR_hvg500 = list(
      ct_col = seurat@misc$cell_type_high_res,
      hvg = 500
    )
  )

  results <- list()

  for (nm in names(plain_combos)) {
    pb_variant <- pb_variants[[plain_combos[[nm]]]]
    if (is.null(pb_variant)) {
      stop("Pseudobulk variant '", plain_combos[[nm]], "' missing for ", nm)
    }
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }
    time_secs <- exec_time(
      res <- process_pseudobulk_fig(pb_variant$pb, labels)
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs") +
      pb_variant$time_secs
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }

  pb_hvg2000 <- pb_variants[["hvg2000"]]
  if (is.null(pb_hvg2000)) {
    stop("Pseudobulk variant 'hvg2000' missing for the PCA-dims combos")
  }
  for (nm in pca_combos) {
    n_pca_dims <- as.integer(sub("Pseudobulk_(\\d+)_PCA_dims", "\\1", nm))
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }
    time_secs <- exec_time(
      res <- process_pseudobulk_fig(
        pb_hvg2000$pb,
        labels,
        pca_dims = n_pca_dims
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs") +
      pb_hvg2000$time_secs
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }

  for (nm in names(ct_combos)) {
    spec <- ct_combos[[nm]]
    if (is.null(spec$ct_col)) {
      warning(nm, " skipped: ct column is null for this dataset")
      next
    }
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }
    time_secs <- exec_time(
      res <- process_pseudobulk_ct_fig(
        seurat,
        labels,
        ct_col = spec$ct_col,
        sample_col = sample_col,
        hvg = spec$hvg
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }

  return(results)
}

# scITD combos: scITD_hvg2000_factors{2,3,5,10,15}, scITD_hvg{1000,3000}_factors5.
# Uses the cell_type_low_res ct column (as the legacy code). Combos where
# num_factors + 5 >= n_samples are skipped with a warning: the tucker rank
# c(f, f+5) must be < n_samples, so the 5-sample _debug dataset cannot run
# scITD at all.
run_scitd_hpc <- function(
  seurat,
  label_col,
  hvg_sets,
  sample_col = "Sample",
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL
) {
  combos <- c(
    paste0("scITD_hvg2000_factors", c(2, 3, 5, 10, 15)),
    "scITD_hvg1000_factors5",
    "scITD_hvg3000_factors5"
  )
  n_samples <- length(unique(seurat@meta.data[[sample_col]]))

  results <- list()
  for (nm in combos) {
    n_hvg <- as.integer(sub("scITD_hvg(\\d+)_factors.*", "\\1", nm))
    num_factors <- as.integer(sub(".*factors", "", nm))
    if (num_factors + 5 >= n_samples) {
      warning(nm, " skipped: tucker rank c(", num_factors, ", ",
              num_factors + 5, ") must be < n_samples (", n_samples, ")")
      next
    }
    hvg <- hvg_sets[[paste0("hvg", n_hvg)]]
    if (is.null(hvg)) {
      stop("HVG set 'hvg", n_hvg, "' missing for ", nm)
    }

    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (file.exists(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file)
      }
      next
    }

    time_secs <- exec_time(
      res <- process_scitd_fig(
        seurat,
        ct_col = seurat@misc$cell_type_low_res,
        label_col = label_col,
        hvg = hvg,
        num_factors = num_factors
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file)
    results[[nm]] <- res
  }
  return(results)
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
  perc_df <- df %>% calc_perc_df()
  res_list <- list()

  for (i in c(0.5, 2/3, 1, 5, 10, 20, 50, 100, 200)) {
    res_list[[paste0("counts_zeros", i)]] <- df %>%
      impute_zeros("counts_zeros", i) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)

    res_list[[paste0("counts_all", i)]] <- df %>%
      impute_zeros("counts_all", i) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)
  }

  for (i in c(0.001, 0.01, 0.1, 1, 2, 5)) {
    res_list[[paste0("percentage_all", i, "%")]] <- df %>%
      impute_zeros("percentage_all", i) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)

    df_multLN <- perc_df %>%
      zCompositions::multLN(label = 0, dl = rep(i, ncol(df)), z.warning = 0.9)
    try(
      res_list[[paste0("multLN", i, "%")]][["multLN"]] <- df_multLN %>%
        clr() %>%
        dist() %>%
        calc_sep_score(labels[row.names(df) %in% row.names(df_multLN)])
    )

    try(
      res_list[[paste0("multRepl_", i, "%")]] <- perc_df %>%
        zCompositions::multRepl(
          label = 0,
          dl = rep(i, ncol(df)),
          z.warning = 1,
          frac = 1
        ) %>%
        clr() %>%
        dist() %>%
        calc_sep_score(labels)
    )
  }

  res_list[["asinsqrt"]] <- perc_df %>%
    mutate(across(everything(), ~ . / 100)) %>%
    sqrt() %>%
    asin() %>%
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
