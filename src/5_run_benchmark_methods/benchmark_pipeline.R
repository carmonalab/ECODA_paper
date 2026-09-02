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
      "cv",
      "clr"
    ),
    # Only dplyr is attached in the worker: the loop body uses its bare verbs
    # (select_if, mutate_all, %>%). Every other package is called
    # namespace-qualified (zCompositions::, Hotelling::, vegan::, mclust::,
    # cluster::, igraph::, Matrix::), so it loads on demand and must NOT be
    # listed here -- attaching packages that are never called bare only adds
    # worker startup time and failure points (e.g. the robCompositions worker
    # abort that previously killed the whole job).
    .packages = c("dplyr"),
    # "stop" (default): a worker error must surface in the .err log with the
    # real message. "pass" returned the error condition object as the result,
    # which crashed the downstream dplyr pipe ("no applicable method for
    # 'group_by'") and masked the actual failure.
    .errorhandling = "stop",
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
# trans/zeroimp) into the notebook's result_list. Every knit starts from a
# fresh list() and loads ALL bundles anew (no result_list.rds persistence,
# no "entries already present are kept" rerun semantics): the feather
# methods recompute in seconds and the stats come from <ds>_metadata.rds,
# so a stale session list can never silently clobber a good file. Bundles
# are verified against the checksums.md5 sidecar (written by the submit
# scripts) before deserialization. Missing files warn and are skipped.
# Returns the (unmodified, in-memory) result_list.
load_hpc_benchmark_results <- function(
  result_list,
  ds,
  path_results_nas,
  methods = c("gloscope", "mofa", "pseudobulk", "scitd", "composition")
) {
  if (is.null(result_list[["bmark"]])) result_list[["bmark"]] <- list()
  if (is.null(result_list[["trans"]])) result_list[["trans"]] <- list()
  if (is.null(result_list[["zeroimp"]])) result_list[["zeroimp"]] <- list()
  if (is.null(result_list[["bmark"]][[ds]])) {
    result_list[["bmark"]][[ds]] <- list()
  }

  checksums <- read_md5_sidecar(file.path(dirname(path_results_nas), "checksums.md5"))

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

  return(result_list)
}

# DEPRECATED (notebook-only caller, kept for reference; no deletion): the
# composition-based methods (ECODA_*, GloProp, Freq_highres, Avg_PCA_embedding,
# ECODA_deconv) moved to the HPC composition worker (run_composition_methods_hpc)
# on 2026-08-16; the python-feather methods below are now called directly by
# benchmark_analysis.rmd with labels from the <ds>_metadata.rds bundle.
run_benchmark_analysis <- function(
  res_list,
  ds,
  seurat,
  sample_col = "Sample",
  factors_test = c(2, 3, 5, 10, 15),
  path_data,
  seurat_res = c(0.1, 0.4, 2, 5, 20),
  HVGs = c(1000, 2000, 3000),
  ECODA_top_varexp_hvct = seq(0, 0.9, 0.1),
  obs = NULL,          # cell-level metadata data.frame (py_to_r(adata$obs));
                       # new-pipeline input mode: overrides the Seurat path
  adata = NULL,        # backed anndata handle; used only to fetch pca_emb
                       # from obsm when pca_emb is not given explicitly
  pca_emb = NULL,      # cells x PCs matrix (obsm X_pca_<view>_hvg2000 read);
                       # Avg_PCA_embedding input on the obs path
  pb_norm = NULL,      # DESeq2-normalized pseudobulk (samples x genes) on the
                       # obs path; expected from the NAS hvg2000 bundle
                       # (<ds>_pseudobulk_hvg2000.rds$pb), which the legacy
                       # get_pb_deseq2(seurat, n_hvg = 2000) call reproduces
  view = "benchmark_analysis",
  label_col = NULL,    # obs-path only (seurat path: seurat@misc$label_col)
  ct_col_low_res = NULL,  # obs-path only (seurat path: seurat@misc$cell_type_low_res)
  ct_col_high_res = NULL  # obs-path only (seurat path: seurat@misc$cell_type_high_res)
) {
  # Files preprocessed with python. Feather names always use the datasets.json
  # key (both the HPC Python pipeline and the legacy qmd write them that way),
  # so the legacy GongSharma -> GongSharma_all remap is dropped. On the
  # new-pipeline path `path_data` points at the NAS benchmark/embeddings dir.
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

  # Input mode: new-pipeline obs data.frame (backed h5ad obs read, no counts
  # access) vs legacy Seurat object. Sample names are already standardized in
  # the preprocessed obs (1.1.1_preprocess.py), so no standardize_sample_names
  # on the obs path.
  use_obs <- !is.null(obs)
  if (use_obs) {
    if (is.null(label_col)) {
      stop("run_benchmark_analysis: label_col is required when obs is used")
    }
    for (req_col in c(sample_col, label_col)) {
      if (!req_col %in% colnames(obs)) {
        stop("run_benchmark_analysis: column '", req_col,
             "' not found in obs")
      }
    }
    # Per-sample metadata + labels (get_metadata/get_labels equivalents on a
    # data.frame: slice(1) per sample, names = Sample).
    metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
    labels <- as.factor(metadata[[label_col]])
    names(labels) <- metadata[[sample_col]]
    # Map the new-pipeline Leiden columns
    # (leiden_res_<r>_<view>_hvg2000) to the legacy RNA_snn_res.* names used
    # by the ECODA_seuratres_* methods below.
    obs <- rename_leiden_cols(obs, view = view)
    # Pseudobulk: the HPC pipeline (prepare_pseudobulks_hpc) precomputed it;
    # the notebook passes the NAS hvg2000 variant. No counts access here.
    if (is.null(pb_norm)) {
      stop("run_benchmark_analysis: pb_norm (NAS pseudobulk hvg2000 bundle ",
           "$pb) is required when obs is used (ECODA_deconv)")
    }
    # Avg_PCA_embedding: per-sample means of the benchmark PCA embedding.
    # Prefer the explicitly passed matrix; fall back to an obsm read when an
    # adata handle is available.
    if (is.null(pca_emb) && !is.null(adata)) {
      obsm_keys <- py_to_r(import_builtins(convert = FALSE)$list(
        adata$obsm$keys()
      ))
      emb_key <- paste0("X_pca_", view, "_hvg2000")
      if (emb_key %in% obsm_keys) {
        pca_emb <- py_to_r(adata$obsm[[emb_key]])
        if (is.null(rownames(pca_emb))) rownames(pca_emb) <- rownames(obs)
      } else {
        warning("run_benchmark_analysis: embedding '", emb_key,
                "' not found in adata.obsm; skipping Avg_PCA_embedding")
      }
    }
  } else {
    if (is.null(seurat)) {
      stop("run_benchmark_analysis: either obs or seurat must be provided")
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
  }

  if (use_obs) {
    if (is.null(ct_col_high_res)) {
      stop("run_benchmark_analysis: ct_col_high_res is required when obs is used")
    }
    ct_low_res <- ct_col_low_res
    ct_high_res <- ct_col_high_res
  } else {
    ct_low_res <- seurat@misc$cell_type_low_res
    ct_high_res <- seurat@misc$cell_type_high_res
  }

  if (!"Avg_PCA_embedding" %in% names(res_list)) {
    res_list[["Avg_PCA_embedding"]][["exec_time"]] <- exec_time(
      res_list[["Avg_PCA_embedding"]] <- process_avg_pca_embedding_fig(
        seurat,
        labels,
        pca_emb = pca_emb,
        obs = obs
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
  if (!is.null(ct_low_res)) {
    res_list[["ECODA_authors_LR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_LR"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = ct_low_res,
        obs = obs
      )
    )
  }

  ## layer2: high res. cell types
  if (!is.null(ct_high_res)) {
    res_list[["ECODA_authors_HR"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = ct_high_res,
        obs = obs
      )
    )
    res_list[["ECODA_authors_HR_NULL"]][["exec_time"]] <- exec_time(
      res_list[["ECODA_authors_HR_NULL"]] <- process_coda_fig(
        seurat,
        labels,
        ct_col = ct_high_res,
        shuffle_labels = TRUE,
        obs = obs
      )
    )
    res_list[["GloProp"]][["exec_time"]] <- exec_time(
      res_list[["GloProp"]] <- process_gloprop_fig(
        seurat,
        metadata,
        ct_col = ct_high_res,
        label_col = label_col,
        obs = obs
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
            ct_col = ct_high_res,
            obs = obs
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
            ct_col = "layer2",
            obs = obs
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
            ct_col = "layer3",
            obs = obs
          )
      )
    }

    res_list[["Freq_highres"]][["exec_time"]] <- exec_time(
      res_list[["Freq_highres"]] <- process_coda_fig(
        seurat,
        labels,
        calc_clr = FALSE,
        ct_col = ct_high_res,
        obs = obs
      )
    )
  }
  res_list[["ECODA_authors_HR_3most_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3most_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 3,
      var_ct_desc = TRUE,
      ct_col = ct_high_res,
      obs = obs
    )
  )
  res_list[["ECODA_authors_HR_2least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_2least_varcts"]] <- process_coda_fig(
      seurat,
      labels,
      ECODA_top_n_hvct = 2,
      var_ct_desc = FALSE,
      ct_col = ct_high_res,
      obs = obs,
    )
  )

  res_list[["ECODA_authors_HR_3least_varcts"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_authors_HR_3least_varcts"]] <-
      process_coda_fig(
        seurat,
        labels,
        ECODA_top_n_hvct = 3,
        var_ct_desc = FALSE,
        ct_col = ct_high_res,
        obs = obs
      )
  )
  res_list[["ECODA_HiTME_HR_layer2"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer2"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "layer2",
      obs = obs
    )
  )
  res_list[["ECODA_HiTME_HR_layer3"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_HiTME_HR_layer3"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "layer3",
      obs = obs
    )
  )
  res_list[["ECODA_scATOMIC_HR"]][["exec_time"]] <- exec_time(
    res_list[["ECODA_scATOMIC_HR"]] <- process_coda_fig(
      seurat,
      labels,
      ct_col = "scATOMIC_pred",
      obs = obs
    )
  )

  # Analyze for all resolutions

  ## Ultra high res. cell type clusters based on Leiden clustering to artificially increase the number of cell types (clusters), e.g. to 250 cell types (clusters)
  for (r in seurat_res) {
    res_col_name <- paste0("RNA_snn_res.", r)
    nm <- paste0("ECODA_seuratres_", r)
    res_list[[nm]][["exec_time"]] <- exec_time(
      res_list[[nm]] <- process_coda_fig(
        seurat, labels, ct_col = res_col_name, obs = obs
      )
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
          ct_col = ct_high_res,
          obs = obs
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

    # --- QOT (Runs once per HVG) ---
    # Pending method (TODO.md Phase 3): feathers may be absent -> skip with a
    # message instead of failing the whole dataset.
    qot_dist_file <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_highres_qot_dists.feather")
    )
    if (file.exists(qot_dist_file)) {
      res_list[[paste0("QOT_hvg", i)]] <- process_qot_fig(
        qot_dist_file = qot_dist_file,
        labels
      )
    } else {
      message("QOT_hvg", i, " skipped for ", ds,
              ": feather not found (pending method, TODO.md Phase 3)")
    }

    # --- PILOT-GM-VAE (Runs once per HVG) ---
    # Pending method (TODO.md Phase 3): feathers may be absent -> skip with a
    # message instead of failing the whole dataset.
    pilotgm_dist_file <- file.path(
      path_data,
      paste0(ds, "_hvg", i, "_highres_pilotgm_dists.feather")
    )
    if (file.exists(pilotgm_dist_file)) {
      res_list[[paste0("PILOT-GM-VAE_hvg", i)]] <- process_pilotgm_fig(
        pilotgm_dist_file = pilotgm_dist_file,
        labels
      )
    } else {
      message("PILOT-GM-VAE_hvg", i, " skipped for ", ds,
              ": feather not found (pending method, TODO.md Phase 3)")
    }

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

      qot_dist_file <- file.path(
        path_data,
        paste0(ds, "_hvg", i, "_lowres_qot_dists.feather")
      )
      if (file.exists(qot_dist_file)) {
        res_list[[paste0("QOT_hvg", i, "_lowres")]] <- process_qot_fig(
          qot_dist_file = qot_dist_file,
          labels
        )
      } else {
        message("QOT_hvg", i, "_lowres skipped for ", ds,
                ": feather not found (pending method, TODO.md Phase 3)")
      }

      pilotgm_dist_file <- file.path(
        path_data,
        paste0(ds, "_hvg", i, "_lowres_pilotgm_dists.feather")
      )
      if (file.exists(pilotgm_dist_file)) {
        res_list[[paste0("PILOT-GM-VAE_hvg", i, "_lowres")]] <- process_pilotgm_fig(
          pilotgm_dist_file = pilotgm_dist_file,
          labels
        )
      } else {
        message("PILOT-GM-VAE_hvg", i, "_lowres skipped for ", ds,
                ": feather not found (pending method, TODO.md Phase 3)")
      }

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
# per-combo exec-log rows. Since 2026-08-16 each bundle also stores
# mem_GB = peak_rss_gb() (VmHWM at combo completion; NA_real_ off-Linux), so
# the re-emit paths replay the ORIGINAL peak on cache reuse instead of the
# live cumulative peak (which would overstate a resumed combo's RAM). Returns
# a named list of result bundles (legacy result names, minus the GloScope
# _sqrtmat suffix).
# ============================================================

# Precompute the shared DESeq2 pseudobulks used by MOFA and Pseudobulk.
# Returns per-variant list(pb, time_secs, mem_GB) (mem_GB = peak_rss_gb() at
# variant completion; old bundles without it re-emit NA). Variants (legacy
# combo list):
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
  variants = PB_VARIANT_NAMES,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
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
        black_list = spec$black_list,
        batch_col = batch_col,
        blind = blind,
        correct_batch = correct_batch
      )
    )
    results[[variant]] <- list(
      pb = pb,
      time_secs = as.numeric(time_secs, units = "secs"),
      mem_GB = peak_rss_gb()
    )
  }
  return(results)
}

# GloScope combos: hvg2000 x pcadims {10,30,50}; hvg1000, hvg3000 x pcadims 30.
# Raw GloScope distances are cached at <gloscope_cache_dir>/<ds>_gloscope_hvg<n>
# _pcadims<d>_dists.rds (sqrt + NA->0 applied by process_gloscope_fig); on a
# cache miss the combo time includes the distance computation, on a hit it is
# sqrt + read only — matching the legacy path_data dist-cache semantics.
# GloScope combos: hvg2000 x pcadims {10,30,50}; hvg1000, hvg3000 x pcadims 30.
# Batch-effect mode deliberately runs only the high-resolution hvg2000/30-PC
# result. Raw GloScope distances are cached under the pass-qualified stem.
run_gloscope_hpc <- function(
  seurat,
  metadata,
  label_col,
  sample_col = "Sample",
  gloscope_cache_dir,
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL,
  batch_mode = FALSE,
  result_stem = ds,
  embedding_name = NULL
) {
  combos <- if (batch_mode) {
    list(list(hvg = 2000, pcadims = 30))
  } else {
    list(
      list(hvg = 2000, pcadims = 10),
      list(hvg = 2000, pcadims = 30),
      list(hvg = 2000, pcadims = 50),
      list(hvg = 1000, pcadims = 30),
      list(hvg = 3000, pcadims = 30)
    )
  }
  artifact_stem <- if (batch_mode) result_stem else ds

  results <- list()
  for (combo in combos) {
    n_hvg <- combo$hvg
    n_pca_dims <- combo$pcadims
    nm <- paste0("GloScope_hvg", n_hvg, "_pcadims", n_pca_dims)
    emb_key <- NULL
    if (batch_mode) {
      if (is.null(embedding_name) || !nzchar(embedding_name)) {
        stop("Batch GloScope requires an exact embedding reduction name")
      }
      emb_key <- embedding_name
      if (!emb_key %in% names(seurat@reductions)) {
        stop("Embedding reduction '", emb_key, "' not found in seurat")
      }
    }
    bundle_file <- file.path(
      results_dir,
      paste0(artifact_stem, "_", nm, ".rds")
    )
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run. The stored mem_GB
      # is replayed too, rather than using the live cumulative peak.
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    if (!batch_mode) {
      emb_key <- paste0("pca_benchmark_analysis_hvg", n_hvg)
      if (!emb_key %in% names(seurat@reductions)) {
        warning("Embedding '", emb_key, "' not found in seurat; skipping ", nm)
        next
      }
    }

    dist_file <- file.path(
      gloscope_cache_dir,
      paste0(
        artifact_stem, "_gloscope_hvg", n_hvg,
        "_pcadims", n_pca_dims, "_dists.rds"
      )
    )
    time_secs <- exec_time(
      res <- process_gloscope_fig(
        embedding_matrix = seurat@reductions[[emb_key]]@cell.embeddings,
        sample_ids = seurat@meta.data[[sample_col]],
        metadata = metadata,
        label_col = label_col,
        gloscope_dist_file = dist_file,
        n_pca_dims = n_pca_dims,
        force = force
      )
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
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
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets). The stored mem_GB
      # is replayed too (the live cumulative peak would overstate the
      # combo's RAM on a resume).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
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
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
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
  log_file = NULL,
  batch_mode = FALSE,
  result_stem = ds
) {
  if (batch_mode) {
    pb_variant <- pb_variants[["hvg2000"]]
    if (is.null(pb_variant)) {
      stop("Pseudobulk variant 'hvg2000' missing for batch-effect run")
    }
    nm <- "Pseudobulk_hvg2000"
    bundle_file <- file.path(results_dir, paste0(result_stem, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file, mem_gb = cached$mem_GB)
      }
      return(setNames(list(cached), nm))
    }
    time_secs <- exec_time(res <- process_pseudobulk_fig(pb_variant$pb, labels))
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs") + pb_variant$time_secs
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file, mem_gb = res[["mem_GB"]])
    return(setNames(list(res), nm))
  }
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
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets). The stored mem_GB
      # is replayed too (the live cumulative peak would overstate the
      # combo's RAM on a resume).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    time_secs <- exec_time(
      res <- process_pseudobulk_fig(pb_variant$pb, labels)
    )
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs") +
      pb_variant$time_secs
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
    results[[nm]] <- res
  }

  pb_hvg2000 <- pb_variants[["hvg2000"]]
  if (is.null(pb_hvg2000)) {
    stop("Pseudobulk variant 'hvg2000' missing for the PCA-dims combos")
  }
  for (nm in pca_combos) {
    n_pca_dims <- as.integer(sub("Pseudobulk_(\\d+)_PCA_dims", "\\1", nm))
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets). The stored mem_GB
      # is replayed too (the live cumulative peak would overstate the
      # combo's RAM on a resume).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
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
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
    results[[nm]] <- res
  }

  for (nm in names(ct_combos)) {
    spec <- ct_combos[[nm]]
    if (is.null(spec$ct_col)) {
      warning(nm, " skipped: ct column is null for this dataset")
      next
    }
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets). The stored mem_GB
      # is replayed too (the live cumulative peak would overstate the
      # combo's RAM on a resume).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
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
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
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
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      # Re-emit the stored timing on cache reuse: failure-resume runs must
      # not lose exec-log rows computed in an aborted run (the merge is
      # scoped to the current run's labels x datasets). The stored mem_GB
      # is replayed too (the live cumulative peak would overstate the
      # combo's RAM on a resume).
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
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
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
    results[[nm]] <- res
  }
  return(results)
}

# Composition-based benchmark methods (the former notebook-local set) as one
# HPC method ("composition"): Avg_PCA_embedding, ECODA_deconv, the
# ECODA_authors_* CoDA family (LR/HR, HR_NULL, top-varexp, top-n/least-varct,
# PCA-dims, seurat-res, HiTME layer2/3, scATOMIC), GloProp and Freq_highres.
# Obs-only: consumes the backed h5ad obs (cell-level metadata), the hvg2000
# obsm PCA embedding (Avg_PCA_embedding) and the hvg2000 DESeq2 pseudobulk
# variant (ECODA_deconv) — no Seurat materialization, no counts access.
# Defaults mirror run_benchmark_analysis (factors_test, seurat_res,
# ECODA_top_varexp_hvct). One worker process per dataset, set.seed(123) at
# driver start -> deterministic across HPC re-runs. NOTE:
# ECODA_authors_HR_NULL (shuffle_labels) will NOT be bit-identical to the
# old notebook runs (the notebook consumed the RNG stream sequentially
# before this method); it is a null control, so the difference is
# inconsequential (documented in benchmark_analysis.rmd).
# HiTME/scATOMIC combos are guarded: skipped with a warning when the ct
# column is absent from obs (annotation produces them, availability varies;
# the old notebook would crash on missing columns).
# Also emits <ds>_metadata.rds = list(labels = named factor, n_cells,
# n_samples, cells_per_sample = named int, n_cell_types_high_res) — replaces
# the notebook's per-dataset obs reads (stats, Supp table 1, exec-times
# n_cells, feather-method labels).
run_composition_methods_hpc <- function(
  labels,
  metadata,
  pca_emb,
  pb_hvg2000,
  obs,
  label_col,
  ct_col_low_res = NULL,
  ct_col_high_res = NULL,
  sample_col = "Sample",
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL,
  factors_test = c(2, 3, 5, 10, 15),
  seurat_res = c(0.1, 0.4, 2, 5, 20),
  ECODA_top_varexp_hvct = seq(0, 0.9, 0.1),
  not_suitable_for_auto_annotation = character(0),
  batch_mode = FALSE,
  result_stem = ds,
  batch_col = NULL,
  corrected = FALSE
) {
  set.seed(123)

  if (corrected && !batch_mode) {
    stop("Corrected CLR composition requires batch-effect mode")
  }
  if (batch_mode && corrected &&
      (is.null(batch_col) || !nzchar(batch_col))) {
    stop("Corrected CLR composition requires a confirmed technical batch")
  }

  skip_hitme <- "hitme" %in% not_suitable_for_auto_annotation
  skip_scatomic <- "scatomic" %in% not_suitable_for_auto_annotation
  if (skip_hitme || skip_scatomic) {
    warning("Dataset '", ds, "' flagged not_suitable_for_auto_annotation = [",
            paste(not_suitable_for_auto_annotation, collapse = ", "),
            "]; skipping annotation-driven composition combos")
  }

  cells_per_sample <- table(factor(
    as.character(obs[[sample_col]]), levels = as.character(names(labels))
  ))
  metadata_bundle <- list(
    labels = labels,
    n_cells = nrow(obs),
    n_samples = length(labels),
    cells_per_sample = cells_per_sample,
    n_cell_types_high_res = if (!is.null(ct_col_high_res) &&
                                 ct_col_high_res %in% colnames(obs)) {
      length(unique(obs[[ct_col_high_res]]))
    } else {
      NA_integer_
    }
  )
  artifact_stem <- if (batch_mode) result_stem else ds
  save_rds_atomic(
    metadata_bundle,
    file.path(results_dir, paste0(artifact_stem, "_metadata.rds"))
  )

  if (batch_mode) {
    if (is.null(ct_col_high_res) ||
        !ct_col_high_res %in% colnames(obs)) {
      stop("Batch composition requires the configured high-resolution cell-type column")
    }
    if (!"RNA_snn_res.2" %in% colnames(obs)) {
      stop("Batch composition requires RNA_snn_res.2 from the selected pass view")
    }

    correction_metadata <- metadata
    if (!"Sample" %in% colnames(correction_metadata)) {
      if (!sample_col %in% colnames(correction_metadata)) {
        stop("Batch composition metadata is missing the sample identifier")
      }
      correction_metadata[["Sample"]] <- correction_metadata[[sample_col]]
    }
    if (corrected && !batch_col %in% colnames(correction_metadata)) {
      stop("Confirmed technical batch column '", batch_col,
           "' is missing from composition metadata")
    }

    correct_result <- function(res) {
      if (!corrected) return(res)
      corrected_feat <- correct_clr_batch_lmm(
        res[["feat_mat"]], correction_metadata, batch_col
      )
      corrected_dist <- dist(corrected_feat)
      corrected_labels <- res[["labels"]]
      if (length(corrected_labels) != nrow(corrected_feat)) {
        stop("Corrected CLR composition labels do not cover all samples")
      }
      names(corrected_labels) <- rownames(corrected_feat)
      res[["feat_mat"]] <- corrected_feat
      res[["dist_mat"]] <- corrected_dist
      res[["labels"]] <- corrected_labels
      res[["scores"]] <- calc_sep_score(corrected_dist, corrected_labels)
      res
    }

    combos <- list()
    add_coda <- function(name, ct_col, shuffle_labels = FALSE) {
      combos[[name]] <<- function() {
        process_coda_fig(
          NULL, labels, ct_col = ct_col, obs = obs,
          clr_zero_impute_method = "counts_all",
          clr_zero_impute_num = 0.5,
          shuffle_labels = shuffle_labels
        )
      }
    }

    add_coda("ECODA_authors_HR", ct_col_high_res)
    add_coda("ECODA_authors_HR_NULL", ct_col_high_res, shuffle_labels = TRUE)

    if (!skip_hitme && "layer2" %in% colnames(obs)) {
      add_coda("ECODA_HiTME_HR_layer2", "layer2")
    } else if (!skip_hitme) {
      warning("ECODA_HiTME_HR_layer2 skipped for ", ds,
              ": obs has no 'layer2' column")
    }

    if (!skip_scatomic && "scATOMIC_pred" %in% colnames(obs)) {
      add_coda("ECODA_scATOMIC_HR", "scATOMIC_pred")
    } else if (!skip_scatomic) {
      warning("ECODA_scATOMIC_HR skipped for ", ds,
              ": obs has no 'scATOMIC_pred' column")
    }

    add_coda("ECODA_seuratres_2", "RNA_snn_res.2")

    results <- list()
    for (nm in names(combos)) {
      bundle_file <- file.path(
        results_dir, paste0(artifact_stem, "_", nm, ".rds")
      )
      if (artifact_checksum_ok(bundle_file) && !force) {
        cached <- readRDS(bundle_file)
        results[[nm]] <- cached
        if (!is.null(cached$exec_time)) {
          log_exec_row(ds, nm, cached$exec_time, log_file,
                       mem_gb = cached$mem_GB)
        }
        next
      }
      time_secs <- exec_time(res <- combos[[nm]]())
      res <- correct_result(res)
      res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
      res[["mem_GB"]] <- peak_rss_gb()
      save_rds_atomic(res, bundle_file)
      log_exec_row(ds, nm, res[["exec_time"]], log_file,
                   mem_gb = res[["mem_GB"]])
      results[[nm]] <- res
    }
    return(results)
  }

  # Ordinary benchmark composition behavior remains unchanged below.
  results <- list()
  combos <- list()

  combos[["Avg_PCA_embedding"]] <- function() {
    process_avg_pca_embedding_fig(NULL, labels, pca_emb = pca_emb, obs = obs)
  }

  if (is.null(pb_hvg2000)) {
    stop("run_composition_methods_hpc: pb_hvg2000 (hvg2000 pseudobulk ",
         "variant) is required for ECODA_deconv")
  }
  combos[["ECODA_deconv"]] <- function() {
    process_deconv_fig(t(pb_hvg2000$pb), labels)
  }

  if (!is.null(ct_col_low_res)) {
    combos[["ECODA_authors_LR"]] <- function() {
      process_coda_fig(NULL, labels, ct_col = ct_col_low_res, obs = obs)
    }
  } else {
    warning("ECODA_authors_LR skipped: cell_type_low_res is null for ", ds)
  }

  if (!is.null(ct_col_high_res)) {
    combos[["ECODA_authors_HR"]] <- function() {
      process_coda_fig(NULL, labels, ct_col = ct_col_high_res, obs = obs)
    }
    combos[["ECODA_authors_HR_NULL"]] <- function() {
      process_coda_fig(NULL, labels, ct_col = ct_col_high_res,
                       shuffle_labels = TRUE, obs = obs)
    }
    combos[["GloProp"]] <- function() {
      process_gloprop_fig(NULL, metadata, ct_col = ct_col_high_res,
                          label_col = label_col, obs = obs)
    }
    combos[["Freq_highres"]] <- function() {
      process_coda_fig(NULL, labels, calc_clr = FALSE,
                       ct_col = ct_col_high_res, obs = obs)
    }

    varexp_combos <- lapply(ECODA_top_varexp_hvct, function(v) {
      function() {
        process_coda_fig(NULL, labels, ECODA_top_varexp_hvct = v,
                         ct_col = ct_col_high_res, obs = obs)
      }
    })
    names(varexp_combos) <- paste0(
      "ECODA_authors_HR_top_varexp", ECODA_top_varexp_hvct
    )
    combos <- c(combos, varexp_combos)

    for (hitme_ct in c("layer2", "layer3")) {
      if (skip_hitme) {
        warning("ECODA_HiTME_HR_", hitme_ct, "* combos skipped for ", ds,
                ": not_suitable_for_auto_annotation incl. 'hitme'")
        next
      }
      if (!hitme_ct %in% colnames(obs)) {
        warning("ECODA_HiTME_HR_", hitme_ct, "* combos skipped for ", ds,
                ": obs has no '", hitme_ct, "' column")
        next
      }
      hitme_combos <- lapply(ECODA_top_varexp_hvct, function(v) {
        function() {
          process_coda_fig(NULL, labels, ECODA_top_varexp_hvct = v,
                           ct_col = hitme_ct, obs = obs)
        }
      })
      names(hitme_combos) <- paste0(
        "ECODA_HiTME_HR_", hitme_ct, "_top_varexp", ECODA_top_varexp_hvct
      )
      combos <- c(combos, hitme_combos)
    }

    combos[["ECODA_authors_HR_3most_varcts"]] <- function() {
      process_coda_fig(NULL, labels, ECODA_top_n_hvct = 3,
                       var_ct_desc = TRUE, ct_col = ct_col_high_res,
                       obs = obs)
    }
    combos[["ECODA_authors_HR_2least_varcts"]] <- function() {
      process_coda_fig(NULL, labels, ECODA_top_n_hvct = 2,
                       var_ct_desc = FALSE, ct_col = ct_col_high_res,
                       obs = obs)
    }
    combos[["ECODA_authors_HR_3least_varcts"]] <- function() {
      process_coda_fig(NULL, labels, ECODA_top_n_hvct = 3,
                       var_ct_desc = FALSE, ct_col = ct_col_high_res,
                       obs = obs)
    }

    for (hitme_ct in c("layer2", "layer3")) {
      if (skip_hitme || !hitme_ct %in% colnames(obs)) next
      hitme_plain <- lapply(hitme_ct, function(ct) {
        function() {
          process_coda_fig(NULL, labels, ct_col = ct, obs = obs)
        }
      })
      names(hitme_plain) <- paste0("ECODA_HiTME_HR_", hitme_ct)
      combos <- c(combos, hitme_plain)
    }
    if (skip_scatomic) {
      warning("ECODA_scATOMIC_HR combos skipped for ", ds,
              ": not_suitable_for_auto_annotation incl. 'scatomic'")
    } else if ("scATOMIC_pred" %in% colnames(obs)) {
      combos[["ECODA_scATOMIC_HR"]] <- function() {
        process_coda_fig(NULL, labels, ct_col = "scATOMIC_pred", obs = obs)
      }
    } else {
      warning("ECODA_scATOMIC_HR skipped for ", ds,
              ": obs has no 'scATOMIC_pred' column")
    }

    pca_combos <- lapply(factors_test, function(i) {
      function() {
        process_coda_fig(NULL, labels, pca_dims = i,
                         ct_col = ct_col_high_res, obs = obs)
      }
    })
    names(pca_combos) <- paste0("ECODA_authors_HR_", factors_test, "_PCA_dims")
    combos <- c(combos, pca_combos)
  } else {
    warning("ECODA_authors_HR* / GloProp / Freq_highres combos skipped for ",
            ds, ": cell_type_high_res is null")
  }

  seuratres_combos <- lapply(seurat_res, function(r) {
    res_col_name <- paste0("RNA_snn_res.", r)
    function() {
      process_coda_fig(NULL, labels, ct_col = res_col_name, obs = obs)
    }
  })
  names(seuratres_combos) <- paste0("ECODA_seuratres_", seurat_res)
  combos <- c(combos, seuratres_combos)

  for (nm in names(combos)) {
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- readRDS(bundle_file)
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    time_secs <- exec_time(res <- combos[[nm]]())
    res[["exec_time"]] <- as.numeric(time_secs, units = "secs")
    res[["mem_GB"]] <- peak_rss_gb()
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
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

  # Method keys use underscore separators with string-formatted values
  # (counts_all_0.5, counts_all_2/3, percentage_all_0.1%, multLN_0.1%, ...):
  # the notebook's legacy references (counts_all_0.5, counts_all_1) match
  # these names as-is, and the paste0-without-separator artifacts
  # (counts_all0.666666666666667) are gone. NOTE: this is a breaking change
  # for already-computed zeroimp bundles — re-run the zeroimp array with
  # --force (see TODO.md).
  counts_vals <- c(
    "0.5" = 0.5, "2/3" = 2 / 3, "1" = 1, "5" = 5, "10" = 10,
    "20" = 20, "50" = 50, "100" = 100, "200" = 200
  )
  for (i in names(counts_vals)) {
    i_val <- unname(counts_vals[[i]])
    res_list[[paste0("counts_zeros_", i)]] <- df %>%
      impute_zeros("counts_zeros", i_val) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)

    res_list[[paste0("counts_all_", i)]] <- df %>%
      impute_zeros("counts_all", i_val) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)
  }

  for (i in c("0.001", "0.01", "0.1", "1", "2", "5")) {
    i_val <- as.numeric(i)
    res_list[[paste0("percentage_all_", i, "%")]] <- df %>%
      impute_zeros("percentage_all", i_val) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)

    res_list[[paste0("percentage_zeros_", i, "%")]] <- df %>%
      impute_zeros("percentage_zeros", i_val) %>%
      clr() %>%
      dist() %>%
      calc_sep_score(labels)

    df_multLN <- perc_df %>%
      zCompositions::multLN(label = 0, dl = rep(i_val, ncol(df)), z.warning = 0.9)
    try(
      res_list[[paste0("multLN_", i, "%")]][["multLN"]] <- df_multLN %>%
        clr() %>%
        dist() %>%
        calc_sep_score(labels[row.names(df) %in% row.names(df_multLN)])
    )

    try(
      res_list[[paste0("multRepl_", i, "%")]] <- perc_df %>%
        zCompositions::multRepl(
          label = 0,
          dl = rep(i_val, ncol(df)),
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
