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

# Read checksums.md5 written by benchmark submitters. The paths are relative
# to the benchmark results root and use GNU md5sum's "<md5>  <path>" format.
read_md5_sidecar <- function(checksum_file) {
  if (!file.exists(checksum_file)) return(NULL)
  lines <- readLines(checksum_file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  parts <- strsplit(lines, "  ", fixed = TRUE)
  if (any(lengths(parts) != 2L)) {
    stop("Malformed checksums.md5 sidecar: ", checksum_file)
  }
  hashes <- vapply(parts, `[`, character(1), 1)
  paths <- vapply(parts, `[`, character(1), 2)
  if (any(!grepl("^[0-9a-fA-F]{32}$", hashes)) ||
      any(!nzchar(paths)) || anyDuplicated(paths)) {
    stop("Malformed checksums.md5 sidecar: ", checksum_file)
  }
  names(hashes) <- paths
  return(hashes)
}

md5_sidecar_key <- function(file_path) {
  file.path(basename(dirname(file_path)), basename(file_path))
}

md5_sidecar_listed <- function(file_path, checksums) {
  if (is.null(checksums) || length(checksums) == 0L) return(FALSE)
  expected <- unname(checksums[md5_sidecar_key(file_path)])
  length(expected) == 1L && !is.na(expected) && nzchar(expected)
}

# Validate a listed file immediately before readRDS. FALSE means the file is
# legacy-unverified (missing/unlisted manifest entry), not trusted content.
verify_md5_sidecar <- function(file_path, checksums) {
  if (!md5_sidecar_listed(file_path, checksums)) return(invisible(FALSE))
  expected <- unname(checksums[md5_sidecar_key(file_path)])
  actual <- tolower(unname(tools::md5sum(file_path)))
  if (length(actual) != 1L || is.na(actual) ||
      tolower(expected) != actual) {
    stop("Checksum mismatch for ", file_path, " (expected ", expected,
         ", got ", actual, "). The HPC result bundle was modified or ",
         "corrupted; re-run the benchmark pipeline (or remove the file) ",
         "before loading it.")
  }
  invisible(TRUE)
}

report_legacy_unverified <- function(file_path, reason) {
  warning(
    "Skipping ", file_path, " (legacy_unverified: ", reason,
    "); readRDS was not attempted."
  )
  invisible(NULL)
}

read_rds_sidecar_checked <- function(file_path, checksums) {
  if (!md5_sidecar_listed(file_path, checksums)) {
    report_legacy_unverified(file_path, "file is not listed in checksums.md5")
    return(NULL)
  }
  if (!isTRUE(verify_md5_sidecar(file_path, checksums))) {
    report_legacy_unverified(file_path, "checksums.md5 entry is unusable")
    return(NULL)
  }
  readRDS(file_path)
}

# Validate timing metadata carried by a cached benchmark result bundle before
# the bundle is reused or replayed.  Pseudobulk cache records have a stricter
# validator in benchmark_hpc_utils.R; this companion accepts both those
# records and downstream method bundles, whose payload omits `pb` and the
# aggregate decomposition.  A missing timing_schema is the legacy contract:
# timing_id alone must not promote an old bundle to schema 2.
validate_hpc_timing_bundle <- function(value, label = "Benchmark result") {
  if (!is.list(value)) stop(label, " is not a list.")
  if (!"timing_schema" %in% names(value)) {
    return(invisible(FALSE))
  }

  schema <- value[["timing_schema"]]
  if (!is.numeric(schema) || length(schema) != 1L ||
      is.na(schema) || !is.finite(schema) ||
      schema != 2 || schema != floor(schema)) {
    stop(label, " has an invalid timing_schema.")
  }

  require_fields <- c(
    "shared_time_secs", "variant_time_secs", "shared_mem_GB", "timing_id"
  )
  missing <- setdiff(require_fields, names(value))
  if (length(missing) > 0L) {
    stop(label, " is missing timing fields: ",
         paste(missing, collapse = ", "))
  }

  validate_nonnegative <- function(raw, field, allow_na = FALSE) {
    if (!is.numeric(raw) || length(raw) != 1L) {
      stop(label, " has invalid ", field, ".")
    }
    if (is.na(raw)) {
      if (allow_na && !is.nan(raw)) return(invisible(TRUE))
      stop(label, " has invalid ", field, ".")
    }
    if (!is.finite(raw) || raw < 0) {
      stop(label, " has invalid ", field, ".")
    }
    invisible(TRUE)
  }
  for (field in c("shared_time_secs", "variant_time_secs")) {
    validate_nonnegative(value[[field]], field)
  }
  for (field in intersect(
    c("time_secs", "exec_time", "aggregate_time_secs",
      "shared_fit_time_secs"),
    names(value)
  )) {
    validate_nonnegative(value[[field]], field)
  }
  for (field in intersect(c("mem_GB", "shared_mem_GB"), names(value))) {
    validate_nonnegative(value[[field]], field, allow_na = TRUE)
  }

  timing_id <- value[["timing_id"]]
  timing_id_parts <- if (is.character(timing_id) && length(timing_id) == 1L &&
                         !is.na(timing_id)) {
    strsplit(timing_id, ":", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  if (!is.character(timing_id) || length(timing_id) != 1L ||
      is.na(timing_id) || !nzchar(trimws(timing_id)) ||
      length(timing_id_parts) != 4L ||
      any(!nzchar(trimws(timing_id_parts))) ||
      grepl("[[:cntrl:]]", timing_id, perl = TRUE)) {
    stop(label, " has an invalid timing_id.")
  }

  if ("pb" %in% names(value)) {
    pb <- value[["pb"]]
    if (!is.matrix(pb) || length(dim(pb)) != 2L ||
        any(dim(pb) <= 0L) || !is.numeric(pb) ||
        any(!is.finite(pb))) {
      stop(label, " has a non-matrix or invalid pb payload.")
    }
    if (!all(c("time_secs", "mem_GB") %in% names(value))) {
      stop(label, " pseudobulk timing record is incomplete.")
    }
  }
  if ("time_secs" %in% names(value) &&
      !isTRUE(all.equal(
        as.numeric(value[["time_secs"]]),
        as.numeric(value[["variant_time_secs"]]),
        tolerance = 0
      ))) {
    stop(label, " time_secs does not equal variant_time_secs.")
  }
  if ("exec_time" %in% names(value)) {
    expected_exec_field <- if ("variant_time_secs" %in% names(value)) {
      "variant_time_secs"
    } else {
      "shared_time_secs"
    }
    if (!isTRUE(all.equal(
      as.numeric(value[["exec_time"]]),
      as.numeric(value[[expected_exec_field]]),
      tolerance = 0
    ))) {
      stop(
        label, " exec_time does not equal ", expected_exec_field, "."
      )
    }
  }
  has_aggregate <- "aggregate_time_secs" %in% names(value)
  has_fit <- "shared_fit_time_secs" %in% names(value)
  if (has_aggregate != has_fit) {
    stop(label, " has an incomplete shared timing decomposition.")
  }
  if (has_aggregate) {
    if (!isTRUE(all.equal(
      as.numeric(value[["shared_time_secs"]]),
      as.numeric(value[["aggregate_time_secs"]]) +
        as.numeric(value[["shared_fit_time_secs"]]),
      tolerance = 0
    ))) {
      stop(label, " shared_time_secs is not aggregate + shared_fit.")
    }
  }
  if ("shared_timing_method" %in% names(value)) {
    shared_method <- value[["shared_timing_method"]]
    if (!is.character(shared_method) || length(shared_method) != 1L ||
        is.na(shared_method) ||
        !grepl("^prepare_pseudobulk_ct_shared_[A-Za-z0-9_-]+$",
               shared_method)) {
      stop(label, " has an invalid shared_timing_method.")
    }
  }
  invisible(TRUE)
}

# Load HPC-computed benchmark results (Pipeline A methods + Pipeline B
# trans/zeroimp) into the notebook's result_list. Every knit starts from a
# fresh list() and loads ALL bundles anew (no result_list.rds persistence,
# no "entries already present are kept" rerun semantics): the feather
# methods recompute in seconds and the stats come from <ds>_metadata.rds,
# so a stale session list can never silently clobber a good file. A listed
# bundle is checksum-verified immediately before deserialization. A missing
# or unlisted checksums.md5 entry is reported as legacy_unverified and skipped.
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

  checksums <- read_md5_sidecar(
    file.path(dirname(path_results_nas), "checksums.md5")
  )
  legacy_reason <- if (is.null(checksums)) {
    "checksums.md5 is missing"
  } else if (length(checksums) == 0L) {
    "checksums.md5 is empty"
  } else {
    NULL
  }

  for (method in methods) {
    method_file <- file.path(
      path_results_nas,
      paste0(ds, "_", method, ".rds")
    )
    if (!file.exists(method_file)) {
      warning("HPC benchmark result file not found: ", method_file)
      next
    }
    if (!is.null(legacy_reason)) {
      report_legacy_unverified(method_file, legacy_reason)
      next
    }
    if (!md5_sidecar_listed(method_file, checksums)) {
      report_legacy_unverified(method_file, "file is not listed in checksums.md5")
      next
    }
    if (!isTRUE(verify_md5_sidecar(method_file, checksums))) {
      report_legacy_unverified(method_file, "checksums.md5 entry is unusable")
      next
    }
    bundles <- read_rds_sidecar_checked(method_file, checksums)
    if (!is.list(bundles)) {
      stop("HPC benchmark result bundle is not a list: ", method_file)
    }
    for (nm in names(bundles)) {
      if (is.list(bundles[[nm]])) {
        validate_hpc_timing_bundle(
          bundles[[nm]],
          label = paste0("HPC benchmark ", ds, "/", method, "/", nm)
        )
      }
      if (!nm %in% names(result_list[["bmark"]][[ds]])) {
        result_list[["bmark"]][[ds]][[nm]] <- bundles[[nm]]
      }
    }
  }

  trans_file <- file.path(path_results_nas, paste0(ds, "_trans.rds"))
  if (is.null(result_list[["trans"]][[ds]])) {
    if (file.exists(trans_file)) {
      if (!is.null(legacy_reason)) {
        report_legacy_unverified(trans_file, legacy_reason)
      } else if (!md5_sidecar_listed(trans_file, checksums)) {
        report_legacy_unverified(trans_file, "file is not listed in checksums.md5")
      } else if (!isTRUE(verify_md5_sidecar(trans_file, checksums))) {
        report_legacy_unverified(trans_file, "checksums.md5 entry is unusable")
      } else {
        result_list[["trans"]][[ds]] <- read_rds_sidecar_checked(
          trans_file, checksums
        )
      }
    } else {
      warning("HPC transformation result file not found: ", trans_file)
    }
  }

  zeroimp_file <- file.path(path_results_nas, paste0(ds, "_zeroimp.rds"))
  if (is.null(result_list[["zeroimp"]][[ds]])) {
    if (file.exists(zeroimp_file)) {
      if (!is.null(legacy_reason)) {
        report_legacy_unverified(zeroimp_file, legacy_reason)
      } else if (!md5_sidecar_listed(zeroimp_file, checksums)) {
        report_legacy_unverified(zeroimp_file, "file is not listed in checksums.md5")
      } else if (!isTRUE(verify_md5_sidecar(zeroimp_file, checksums))) {
        report_legacy_unverified(zeroimp_file, "checksums.md5 entry is unusable")
      } else {
        result_list[["zeroimp"]][[ds]] <- read_rds_sidecar_checked(
          zeroimp_file, checksums
        )
      }
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

# Focused tests may source this pipeline file without the HPC utility file.
# Keep the canonical timing-ID construction available in that isolated
# context; the production load order resolves the shared helper instead.
if (!exists("pb_timing_id", mode = "function", inherits = TRUE)) {
  pb_timing_id <- function(
    cache_stem, view = "benchmark_analysis", analysis_pass = NULL,
    run_id = NULL
  ) {
    if (is.null(run_id) || !nzchar(as.character(run_id))) {
      run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
    }
    if (!nzchar(as.character(run_id))) run_id <- "unbound"
    pass <- if (is.null(analysis_pass)) "none" else as.character(analysis_pass)
    paste(as.character(run_id), as.character(cache_stem), as.character(view),
          pass, sep = ":")
  }
}

# Focused tests may source this pipeline file without the HPC utility file.
# Keep CT method/token construction available in that isolated context; the
# production load order resolves the shared helpers from benchmark_hpc_utils.R.
if (!exists("ct_timing_token", mode = "function", inherits = TRUE)) {
  .ct_timing_validate_token <- function(
    token, label = "CT timing token"
  ) {
    if (is.factor(token)) token <- as.character(token)
    if (!is.character(token) || length(token) != 1L ||
        is.na(token) || !nzchar(token) ||
        !grepl("^[A-Za-z0-9_-]+$", token, perl = TRUE)) {
      stop(label, " must be a non-empty safe token.")
    }
    token
  }
  ct_timing_token <- function(ct_col) {
    if (is.factor(ct_col)) ct_col <- as.character(ct_col)
    if (!is.character(ct_col) || length(ct_col) != 1L ||
        is.na(ct_col) || !nzchar(ct_col)) {
      stop("CT column must be one non-empty string.")
    }
    ct_col <- tryCatch(
      enc2utf8(ct_col),
      error = function(error) {
        stop("CT column is not valid UTF-8: ", conditionMessage(error))
      }
    )
    if (grepl("[[:cntrl:]]", ct_col, perl = TRUE)) {
      stop("CT column contains control characters.")
    }
    raw <- charToRaw(ct_col)
    token <- paste(sprintf("%02x", as.integer(raw)), collapse = "")
    .ct_timing_validate_token(token)
  }
  ct_shared_timing_method_from_token <- function(token) {
    token <- .ct_timing_validate_token(token)
    paste0("prepare_pseudobulk_ct_shared_", token)
  }
  ct_shared_timing_method <- function(ct_col) {
    ct_shared_timing_method_from_token(ct_timing_token(ct_col))
  }
  ct_shared_timing_method_from_timing_id <- function(
    timing_id, fallback_method = NULL
  ) {
    if (is.factor(timing_id)) timing_id <- as.character(timing_id)
    if (!is.character(timing_id) || length(timing_id) != 1L ||
        is.na(timing_id) || !nzchar(trimws(timing_id))) {
      stop("CT timing_id must be one non-empty string.")
    }
    parts <- strsplit(timing_id, ":", fixed = TRUE)[[1L]]
    token <- NULL
    if (length(parts) == 4L && grepl("_ct_", parts[[2L]], fixed = TRUE)) {
      candidate <- sub("^.*_ct_", "", parts[[2L]])
      if (length(candidate) == 1L &&
          !is.na(candidate) && nzchar(candidate) &&
          grepl("^[A-Za-z0-9_-]+$", candidate, perl = TRUE)) {
        token <- candidate
      }
    }
    if (is.null(token) && !is.null(fallback_method)) {
      if (is.factor(fallback_method)) {
        fallback_method <- as.character(fallback_method)
      }
      if (is.character(fallback_method) && length(fallback_method) == 1L &&
          !is.na(fallback_method) &&
          grepl("^Pseudobulk_CT_[^_]+_.*$", fallback_method)) {
        candidate <- sub(
          "^Pseudobulk_CT_([^_]+)_.*$", "\\1", fallback_method
        )
        if (nzchar(candidate) &&
            grepl("^[A-Za-z0-9_-]+$", candidate, perl = TRUE)) {
          token <- candidate
        }
      }
    }
    if (is.null(token)) {
      stop("CT timing identity does not contain a recoverable token.")
    }
    ct_shared_timing_method_from_token(token)
  }
}


# Precompute shared DESeq2 pseudobulks from one raw H5AD Sample aggregate.
# Full-gene HVG variants share one fit; schvg2000 intentionally uses a
# separate pre-filtered fit.  The returned cache payload retains the legacy
# $pb/$time_secs/$mem_GB fields and adds schema-2 shared timing fields.
prepare_pseudobulks_hpc <- function(
  h5ad_path,
  sample_col = "Sample",
  hvg_rank_genes = NULL,
  variants = PB_VARIANT_NAMES,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE,
  cache_stem = basename(h5ad_path),
  view = "benchmark_analysis",
  analysis_pass = NULL,
  run_id = NULL,
  source_identity = NULL,
  chunk_size = 4096L
) {
  if (!is.character(h5ad_path) || length(h5ad_path) != 1L ||
      is.na(h5ad_path) || !nzchar(h5ad_path) || !file.exists(h5ad_path)) {
    stop("prepare_pseudobulks_hpc: H5AD path is missing: ", h5ad_path)
  }
  variants <- unique(as.character(variants))
  unknown <- setdiff(variants, PB_VARIANT_NAMES)
  if (length(unknown) > 0L) {
    stop("Unknown pseudobulk variant requested: ",
         paste(unknown, collapse = ", "))
  }
  if (length(variants) == 0L) {
    stop("No pseudobulk variants requested.")
  }
  if (is.null(hvg_rank_genes) || length(hvg_rank_genes) == 0L) {
    stop("prepare_pseudobulks_hpc requires ranked HVG genes.")
  }
  hvg_rank_genes <- as.character(hvg_rank_genes)
  if (anyNA(hvg_rank_genes) || any(!nzchar(hvg_rank_genes)) ||
      anyDuplicated(hvg_rank_genes)) {
    stop("prepare_pseudobulks_hpc received invalid ranked HVG genes.")
  }

  specs <- list(
    schvg2000 = list(n_hvg = 2000L, black_list = "none"),
    hvg2000 = list(n_hvg = 2000L, black_list = "none"),
    hvg500 = list(n_hvg = 500L, black_list = "none"),
    hvg2000_bl = list(
      n_hvg = 2000L, black_list = "default_without_sex_genes"
    ),
    hvg1000 = list(n_hvg = 1000L, black_list = "none"),
    hvg3000 = list(n_hvg = 3000L, black_list = "none")
  )
  timing_id <- pb_timing_id(
    cache_stem = cache_stem,
    view = view,
    analysis_pass = analysis_pass,
    run_id = run_id
  )

  # This call is the only raw Sample aggregation for the entire requested
  # variant set.  Biological labels are intentionally not requested.
  aggregate_time <- exec_time(
    aggregated <- load_h5ad_sample_aggregate(
      h5ad_path,
      sample_col = sample_col,
      metadata_columns = unique(c(sample_col, batch_col)),
      chunk_size = chunk_size,
      max_value = .Machine$integer.max
    )
  )
  aggregate_time <- as.numeric(aggregate_time, units = "secs")

  full_names <- c("hvg500", "hvg1000", "hvg2000", "hvg2000_bl", "hvg3000")
  shared_names <- intersect(variants, full_names)
  shared_fit <- NULL
  shared_fit_time <- 0
  shared_mem <- NA_real_
  if (length(shared_names) > 0L) {
    shared_fit_time <- exec_time(
      shared_fit <- fit_pseudobulk_deseq2(
        aggregated$counts,
        metadata = aggregated$metadata,
        batch_col = batch_col,
        blind = blind,
        correct_batch = correct_batch
      )
    )
    shared_fit_time <- as.numeric(shared_fit_time, units = "secs")
    shared_mem <- peak_rss_gb()
  }
  shared_time <- aggregate_time + shared_fit_time

  results <- list()
  for (variant in variants) {
    spec <- specs[[variant]]
    variant_start <- Sys.time()
    if (identical(variant, "schvg2000")) {
      ranked <- hvg_rank_genes[seq_len(min(2000L, length(hvg_rank_genes)))]
      positions <- match(ranked, rownames(aggregated$counts))
      positions <- positions[!is.na(positions)]
      if (length(positions) == 0L) {
        stop("schvg2000 has no ranked genes in the H5AD gene universe.")
      }
      schvg_fit_time <- exec_time(
        selected_fit <- fit_pseudobulk_deseq2(
          aggregated$counts[positions, , drop = FALSE],
          metadata = aggregated$metadata,
          batch_col = batch_col,
          blind = blind,
          correct_batch = correct_batch
        )
      )
      selected_fit_time <- as.numeric(schvg_fit_time, units = "secs")
      selected_time <- exec_time(
        selected <- select_pseudobulk_deseq2(
          selected_fit,
          n_hvg = spec$n_hvg,
          black_list = spec$black_list
        )
      )
      selected_time <- as.numeric(selected_time, units = "secs")
      variant_time <- selected_fit_time + selected_time
    } else {
      if (is.null(shared_fit)) {
        stop("Shared full-gene fit is missing for ", variant)
      }
      selected_time <- exec_time(
        selected <- select_pseudobulk_deseq2(
          shared_fit,
          n_hvg = spec$n_hvg,
          black_list = spec$black_list
        )
      )
      selected_time <- as.numeric(selected_time, units = "secs")
      variant_time <- selected_time
    }
    if (!is.matrix(selected)) selected <- as.matrix(selected)
    pb <- t(selected)
    canonical_samples <- as.character(aggregated$sample_ids)
    if (is.null(rownames(pb)) || is.null(colnames(pb))) {
      stop("Direct pseudobulk selection returned unnamed dimensions for ", variant)
    }
    pb <- pb[canonical_samples, , drop = FALSE]
    variant_elapsed <- as.numeric(
      difftime(Sys.time(), variant_start, units = "secs")
    )
    # Include any unmeasured R-side selection bookkeeping, while keeping the
    # shared aggregate/fit out of every variant-local charge.
    variant_time <- max(variant_time, variant_elapsed)
    variant_mem <- peak_rss_gb()
    results[[variant]] <- list(
      pb = pb,
      time_secs = variant_time,
      mem_GB = variant_mem,
      aggregate_time_secs = aggregate_time,
      shared_fit_time_secs = shared_fit_time,
      shared_time_secs = shared_time,
      variant_time_secs = variant_time,
      shared_mem_GB = shared_mem,
      timing_id = timing_id,
      timing_schema = 2L
    )
  }
  results
}

# GloScope combos: hvg2000 x pcadims {10,30,50}; hvg1000, hvg3000 x pcadims 30.
# Raw GloScope distances are cached at <gloscope_cache_dir>/<ds>_gloscope_hvg<n>
# _pcadims<d>_dists.rds (sqrt + NA->0 applied by process_gloscope_fig); on a
# cache miss the combo time includes the distance computation, on a hit it is
# sqrt + read only — matching the legacy path_data dist-cache semantics.
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
  embedding_name = NULL,
  combo_token = NULL,
  embedding_matrices = NULL,
  embedding_sample_ids = NULL
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
  if (!is.null(combo_token)) {
    if (batch_mode ||
        !is.character(combo_token) || length(combo_token) != 1L ||
        !nzchar(combo_token)) {
      stop("--combo is only supported for ordinary GloScope runs")
    }
    combo_names <- vapply(
      combos,
      function(combo) paste0("hvg", combo$hvg, "_pcadims", combo$pcadims),
      character(1)
    )
    matches <- which(combo_names == combo_token)
    if (length(matches) != 1L) {
      stop(
        "Unknown or ambiguous GloScope combo '", combo_token,
        "'; expected one of ", paste(combo_names, collapse = ", ")
      )
    }
    combos <- combos[matches]
  }
  artifact_stem <- if (batch_mode) result_stem else ds

  results <- list()
  for (combo in combos) {
    n_hvg <- combo$hvg
    n_pca_dims <- combo$pcadims
    nm <- paste0("GloScope_hvg", n_hvg, "_pcadims", n_pca_dims)
    bundle_file <- file.path(
      results_dir,
      paste0(artifact_stem, "_", nm, ".rds")
    )
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("GloScope result ", ds, "/", nm)
      )
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }

    emb_key <- NULL
    embedding_matrix <- NULL
    sample_ids <- NULL
    if (!is.null(seurat)) {
      if (batch_mode) {
        if (is.null(embedding_name) || !nzchar(embedding_name)) {
          stop("Batch GloScope requires an exact embedding reduction name")
        }
        emb_key <- embedding_name
      } else {
        emb_key <- paste0("pca_benchmark_analysis_hvg", n_hvg)
      }
      if (!emb_key %in% names(seurat@reductions)) {
        warning("Embedding '", emb_key, "' not found in seurat; skipping ", nm)
        next
      }
      embedding_matrix <- seurat@reductions[[emb_key]]@cell.embeddings
      sample_ids <- seurat@meta.data[[sample_col]]
    } else {
      if (!is.list(embedding_matrices) ||
          is.null(embedding_matrices[[paste0("hvg", n_hvg)]]) ||
          is.null(embedding_sample_ids)) {
        stop(
          "GloScope requires embedding matrices and sample IDs when Seurat is absent"
        )
      }
      embedding_matrix <- embedding_matrices[[paste0("hvg", n_hvg)]]
      sample_ids <- embedding_sample_ids
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
        embedding_matrix = embedding_matrix,
        sample_ids = sample_ids,
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
  results
}

# MOFA combos: MOFA_hvg2000_factors{2,3,5,10,15}, MOFA_hvg{1000,3000}_factors15.
# Shared pseudobulk preparation is not charged to each combo.  A schema-2
# result's exec_time is MOFA runtime plus the selected variant-local time;
# preparation's shared row is emitted by load_pb_variants().
run_mofa_hpc <- function(
  metadata,
  labels,
  pb_variants,
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL
) {
  local_variant_time <- function(value) {
    if (!is.list(value)) {
      stop("Pseudobulk variant timing is invalid.")
    }
    if ("timing_schema" %in% names(value)) {
      validate_hpc_timing_bundle(value, "Pseudobulk variant")
      raw <- value[["variant_time_secs"]]
    } else {
      raw <- value[["time_secs"]]
    }
    if (!is.numeric(raw) || length(raw) != 1L ||
        is.na(raw) || !is.finite(raw) || raw < 0) {
      stop("Pseudobulk variant timing is invalid.")
    }
    as.numeric(raw)
  }
  local_shared_time <- function(value) {
    if (!is.list(value) || !"timing_schema" %in% names(value)) return(0)
    validate_hpc_timing_bundle(value, "Pseudobulk variant")
    as.numeric(value[["shared_time_secs"]])
  }
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
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("MOFA result ", ds, "/", nm)
      )
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }

    process_time <- exec_time(
      res <- process_mofa_bulk_fig(
        pb_variant$pb,
        metadata = metadata,
        labels,
        num_factors = num_factors
      )
    )
    process_time <- as.numeric(process_time, units = "secs")
    variant_time <- local_variant_time(pb_variant)
    res[["exec_time"]] <- process_time + variant_time
    if ("timing_schema" %in% names(pb_variant)) {
      res[["variant_time_secs"]] <- process_time + variant_time
      res[["shared_time_secs"]] <- local_shared_time(pb_variant)
      res[["timing_schema"]] <- 2L
      res[["timing_id"]] <- as.character(pb_variant[["timing_id"]])
      res[["shared_mem_GB"]] <- pb_variant[["shared_mem_GB"]]
    }
    res[["mem_GB"]] <- peak_rss_gb()
    validate_hpc_timing_bundle(
      res, label = paste0("MOFA result ", ds, "/", nm)
    )
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
    results[[nm]] <- res
  }
  results
}

# Pseudobulk combos: Pseudobulk_schvg2000, Pseudobulk_hvg2000,
# Pseudobulk_hvg500, Pseudobulk_hvg2000_bl, Pseudobulk_CT_LR_hvg2000,
# Pseudobulk_CT_HR_hvg2000, Pseudobulk_CT_LR_hvg500, Pseudobulk_CT_HR_hvg500,
# Pseudobulk_{2,3,5,10,15}_PCA_dims, Pseudobulk_hvg1000, Pseudobulk_hvg3000.
# Plain variants reuse precomputed direct-matrix pseudobulks.  CT variants
# invoke the store-backed H5AD processor (the multi-variant helper when one
# CT column has multiple pending HVGs); the Seurat CT routine remains only as
# an explicit legacy fallback.
run_pseudobulk_hpc <- function(
  seurat = NULL,
  labels,
  pb_variants,
  sample_col = "Sample",
  results_dir,
  ds,
  force = FALSE,
  log_file = NULL,
  batch_mode = FALSE,
  result_stem = ds,
  h5ad_path = NULL,
  ct_col_low_res = NULL,
  ct_col_high_res = NULL,
  view = "benchmark_analysis",
  analysis_pass = NULL,
  run_id = NULL,
  temp_root = NULL,
  source_identity = NULL,
  chunk_size = 4096L
) {
  local_variant_time <- function(value) {
    if (!is.list(value)) {
      stop("Pseudobulk variant timing is invalid.")
    }
    if ("timing_schema" %in% names(value)) {
      validate_hpc_timing_bundle(value, "Pseudobulk variant")
      raw <- value[["variant_time_secs"]]
    } else {
      raw <- value[["time_secs"]]
    }
    if (!is.numeric(raw) || length(raw) != 1L ||
        is.na(raw) || !is.finite(raw) || raw < 0) {
      stop("Pseudobulk variant timing is invalid.")
    }
    as.numeric(raw)
  }
  local_shared_time <- function(value) {
    if (!is.list(value) || !"timing_schema" %in% names(value)) return(0)
    validate_hpc_timing_bundle(value, "Pseudobulk variant")
    as.numeric(value[["shared_time_secs"]])
  }
  add_pb_timing <- function(res, process_time, pb_variant) {
    process_time <- as.numeric(process_time)
    if (!is.finite(process_time) || process_time < 0) {
      stop("Pseudobulk process timing is invalid.")
    }
    res[["exec_time"]] <- process_time + local_variant_time(pb_variant)
    if ("timing_schema" %in% names(pb_variant)) {
      res[["variant_time_secs"]] <- res[["exec_time"]]
      res[["shared_time_secs"]] <- local_shared_time(pb_variant)
      res[["timing_schema"]] <- 2L
      res[["timing_id"]] <- as.character(pb_variant[["timing_id"]])
      res[["shared_mem_GB"]] <- pb_variant[["shared_mem_GB"]]
    }
    res[["mem_GB"]] <- peak_rss_gb()
    validate_hpc_timing_bundle(res, "Pseudobulk result")
    res
  }

  if (batch_mode) {
    pb_variant <- pb_variants[["hvg2000"]]
    if (is.null(pb_variant)) {
      stop("Pseudobulk variant 'hvg2000' missing for batch-effect run")
    }
    nm <- "Pseudobulk_hvg2000"
    bundle_file <- file.path(results_dir, paste0(result_stem, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("Pseudobulk result ", ds, "/", nm)
      )
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file, mem_gb = cached$mem_GB)
      }
      return(setNames(list(cached), nm))
    }
    process_time <- exec_time(
      res <- process_pseudobulk_fig(pb_variant$pb, labels)
    )
    res <- add_pb_timing(
      res, as.numeric(process_time, units = "secs"), pb_variant
    )
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
  if (is.null(ct_col_low_res) && !is.null(seurat)) {
    ct_col_low_res <- seurat@misc$cell_type_low_res
  }
  if (is.null(ct_col_high_res) && !is.null(seurat)) {
    ct_col_high_res <- seurat@misc$cell_type_high_res
  }
  ct_combos <- list(
    Pseudobulk_CT_LR_hvg2000 = list(ct_col = ct_col_low_res, hvg = 2000),
    Pseudobulk_CT_HR_hvg2000 = list(ct_col = ct_col_high_res, hvg = 2000),
    Pseudobulk_CT_LR_hvg500 = list(ct_col = ct_col_low_res, hvg = 500),
    Pseudobulk_CT_HR_hvg500 = list(ct_col = ct_col_high_res, hvg = 500)
  )

  results <- list()
  for (nm in names(plain_combos)) {
    pb_variant <- pb_variants[[plain_combos[[nm]]]]
    if (is.null(pb_variant)) {
      stop("Pseudobulk variant '", plain_combos[[nm]], "' missing for ", nm)
    }
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    if (artifact_checksum_ok(bundle_file) && !force) {
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("Pseudobulk result ", ds, "/", nm)
      )
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    process_time <- exec_time(
      res <- process_pseudobulk_fig(pb_variant$pb, labels)
    )
    res <- add_pb_timing(
      res, as.numeric(process_time, units = "secs"), pb_variant
    )
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
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("Pseudobulk result ", ds, "/", nm)
      )
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    process_time <- exec_time(
      res <- process_pseudobulk_fig(
        pb_hvg2000$pb,
        labels,
        pca_dims = n_pca_dims
      )
    )
    res <- add_pb_timing(
      res, as.numeric(process_time, units = "secs"), pb_hvg2000
    )
    save_rds_atomic(res, bundle_file)
    log_exec_row(ds, nm, res[["exec_time"]], log_file,
                 mem_gb = res[["mem_GB"]])
    results[[nm]] <- res
  }

  # CT combinations sharing a column also share one composite H5AD store.
  # Resolve/cache-hit bundles first so a complete CT result set never opens
  # the counts layer.  A partial set passes only missing HVG variants to the
  # reusable store-backed processor, preserving missing-only publication.
  ct_names <- names(ct_combos)
  # Validate configured CT-column eligibility before touching any CT bundle.
  # A stale result for a dataset whose CT column is now null is not a valid
  # cache hit and must not even be deserialized.
  eligible_ct_names <- ct_names[vapply(ct_names, function(nm) {
    value <- ct_combos[[nm]]$ct_col
    !is.null(value) && length(value) == 1L &&
      !is.na(value) && nzchar(as.character(value))
  }, logical(1))]
  ct_cache_hit <- setNames(rep(FALSE, length(ct_names)), ct_names)
  if (length(eligible_ct_names) > 0L) {
    ct_cache_hit[eligible_ct_names] <- vapply(
      eligible_ct_names,
      function(nm) {
        bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
        artifact_checksum_ok(bundle_file) && !force
      },
      logical(1)
    )
  }
  ct_shared_replay <- list()
  record_ct_shared <- function(timing, context) {
    if (is.null(timing)) return(invisible(NULL))
    required <- c(
      "shared_time_secs", "shared_mem_GB", "timing_id",
      "timing_schema", "shared_timing_method"
    )
    if (!all(required %in% names(timing))) {
      stop("CT timing metadata is incomplete for ", context)
    }
    shared_time <- timing[["shared_time_secs"]]
    shared_mem <- timing[["shared_mem_GB"]]
    timing_id <- timing[["timing_id"]]
    shared_method <- timing[["shared_timing_method"]]
    if (!is.numeric(shared_time) || length(shared_time) != 1L ||
        is.na(shared_time) || !is.finite(shared_time) || shared_time < 0 ||
        !is.numeric(shared_mem) || length(shared_mem) != 1L ||
        (is.na(shared_mem) && is.nan(shared_mem)) ||
        (!is.na(shared_mem) &&
         (!is.finite(shared_mem) || shared_mem < 0)) ||
        !is.character(timing_id) || length(timing_id) != 1L ||
        is.na(timing_id) || !nzchar(trimws(timing_id)) ||
        !is.character(shared_method) || length(shared_method) != 1L ||
        is.na(shared_method) ||
        !grepl("^prepare_pseudobulk_ct_shared_[A-Za-z0-9_-]+$",
               shared_method)) {
      stop("CT timing metadata is invalid for ", context)
    }
    prior <- ct_shared_replay[[timing_id]]
    if (!is.null(prior)) {
      if (!isTRUE(all.equal(prior$time, as.numeric(shared_time),
                            tolerance = 0)) ||
          !identical(prior$mem, shared_mem) ||
          !identical(prior$method, shared_method)) {
        stop("CT shared timing differs for timing_id ", timing_id)
      }
      return(invisible(NULL))
    }
    ct_shared_replay[[timing_id]] <<- list(
      time = as.numeric(shared_time), mem = shared_mem, method = shared_method
    )
    invisible(NULL)
  }

  # Deserialize only eligible CT hits, validate schema-2 timing before use,
  # and collect one shared-store row per timing identity.  The shared rows are
  # emitted below, after all hit bundles have been read, so one grouped store
  # can never be charged once per HVG.
  for (nm in eligible_ct_names[ct_cache_hit[eligible_ct_names]]) {
    bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
    cached <- read_rds_checked(bundle_file)
    validate_hpc_timing_bundle(
      cached, label = paste0("Pseudobulk result ", ds, "/", nm)
    )
    results[[nm]] <- cached
    if ("timing_schema" %in% names(cached)) {
      record_ct_shared(
        list(
          shared_time_secs = cached[["shared_time_secs"]],
          shared_mem_GB = cached[["shared_mem_GB"]],
          timing_id = cached[["timing_id"]],
          timing_schema = cached[["timing_schema"]],
          shared_timing_method = if (
            "shared_timing_method" %in% names(cached)
          ) {
            cached[["shared_timing_method"]]
          } else {
            ct_shared_timing_method(ct_combos[[nm]][["ct_col"]])
          }
        ),
        paste0(ds, "/", nm)
      )
    }
    if (!is.null(cached$exec_time)) {
      log_exec_row(ds, nm, cached$exec_time, log_file,
                   mem_gb = cached$mem_GB)
    }
  }
  for (nm in ct_names[!ct_cache_hit]) {
    if (!nm %in% eligible_ct_names) {
      warning(nm, " skipped: ct column is null for this dataset")
    }
  }
  if (length(eligible_ct_names) > 0L) {
    ct_groups <- split(
      eligible_ct_names,
      vapply(
        eligible_ct_names,
        function(nm) as.character(ct_combos[[nm]]$ct_col),
        character(1)
      )
    )
    for (ct_col_name in names(ct_groups)) {
      group_names <- ct_groups[[ct_col_name]]
      pending_names <- group_names[!ct_cache_hit[group_names]]
      if (length(pending_names) == 0L) next
      if (!is.null(h5ad_path)) {
        hvgs <- unique(vapply(
          pending_names,
          function(nm) as.integer(ct_combos[[nm]]$hvg),
          integer(1)
        ))
        ct_token <- ct_timing_token(ct_col_name)
        ct_timing_id <- pb_timing_id(
          cache_stem = paste0(result_stem, "_ct_", ct_token),
          view = view,
          analysis_pass = analysis_pass,
          run_id = run_id
        )
        variant_formals <- tryCatch(
          names(formals(process_pseudobulk_ct_h5ad_variants_fig)),
          error = function(error) character()
        )
        wrapper_formals <- tryCatch(
          names(formals(process_pseudobulk_ct_h5ad_fig)),
          error = function(error) character()
        )
        supports_timing <- "timing_id" %in% variant_formals &&
          (length(hvgs) > 1L || "timing_id" %in% wrapper_formals)
        call_timed <- function(fun, call_args, formals_names) {
          if ("timing_id" %in% formals_names) {
            call_args[["timing_id"]] <- ct_timing_id
          }
          do.call(fun, call_args)
        }

        # The current reducer exposes one shared store timing and one local
        # timing per HVG.  Keep a compatibility path for an older sourced
        # reducer by measuring each single-HVG call independently; it avoids
        # assigning a grouped elapsed time to every legacy result.
        if (supports_timing) {
          if (length(hvgs) > 1L) {
            ct_processed <- call_timed(
              process_pseudobulk_ct_h5ad_variants_fig,
              list(
                h5ad_path = h5ad_path,
                labels = labels,
                sample_col = sample_col,
                ct_col = ct_col_name,
                hvgs = hvgs,
                min_cells = 5,
                chunk_size = chunk_size,
                run_id = run_id,
                temp_root = temp_root,
                source_identity = source_identity
              ),
              variant_formals
            )
          } else {
            ct_processed <- setNames(
              list(call_timed(
                process_pseudobulk_ct_h5ad_fig,
                list(
                  h5ad_path = h5ad_path,
                  labels = labels,
                  sample_col = sample_col,
                  ct_col = ct_col_name,
                  hvg = hvgs[[1L]],
                  min_cells = 5,
                  chunk_size = chunk_size,
                  run_id = run_id,
                  temp_root = temp_root,
                  source_identity = source_identity
                ),
                wrapper_formals
              )),
              paste0("hvg", hvgs[[1L]])
            )
          }
        } else {
          ct_processed <- setNames(lapply(hvgs, function(n_hvg) {
            result <- NULL
            elapsed <- exec_time(
              result <- process_pseudobulk_ct_h5ad_fig(
                h5ad_path,
                labels,
                sample_col = sample_col,
                ct_col = ct_col_name,
                hvg = n_hvg,
                min_cells = 5,
                chunk_size = chunk_size,
                run_id = run_id,
                temp_root = temp_root,
                source_identity = source_identity
              )
            )
            attr(result, "legacy_ct_elapsed") <- as.numeric(
              elapsed, units = "secs"
            )
            result
          }), paste0("hvg", hvgs))
        }

        extract_ct_timing <- function(result, context) {
          if (!is.list(result)) {
            stop("CT processor returned a non-list result for ", context)
          }
          if ("timing_schema" %in% names(result)) {
            validate_hpc_timing_bundle(result, context)
            return(list(
              local = as.numeric(result[["variant_time_secs"]]),
              shared = as.numeric(result[["shared_time_secs"]]),
              mem = result[["shared_mem_GB"]],
              id = as.character(result[["timing_id"]]),
              schema2 = TRUE
            ))
          }
          timing <- attr(result, "ct_timing", exact = TRUE)
          if (is.list(timing) && !is.null(timing[["timing_id"]])) {
            required <- c(
              "shared_time_secs", "variant_time_secs", "shared_mem_GB",
              "timing_id", "timing_schema"
            )
            if (!all(required %in% names(timing))) {
              stop("CT timing metadata is incomplete for ", context)
            }
            timing_bundle <- list(
              shared_time_secs = timing[["shared_time_secs"]],
              variant_time_secs = timing[["variant_time_secs"]],
              shared_mem_GB = timing[["shared_mem_GB"]],
              timing_id = timing[["timing_id"]],
              timing_schema = timing[["timing_schema"]]
            )
            timing_bundle[["exec_time"]] <- timing_bundle[[
              "variant_time_secs"
            ]]
            validate_hpc_timing_bundle(timing_bundle, context)
            return(list(
              local = as.numeric(timing[["variant_time_secs"]]),
              shared = as.numeric(timing[["shared_time_secs"]]),
              mem = timing[["shared_mem_GB"]],
              id = as.character(timing[["timing_id"]]),
              schema2 = TRUE
            ))
          }
          elapsed <- attr(result, "legacy_ct_elapsed", exact = TRUE)
          if (!is.numeric(elapsed) || length(elapsed) != 1L ||
              is.na(elapsed) || !is.finite(elapsed) || elapsed < 0) {
            stop("CT processor did not return timing metadata for ", context)
          }
          list(
            local = as.numeric(elapsed),
            shared = NA_real_,
            mem = NA_real_,
            id = NULL,
            schema2 = FALSE
          )
        }

        ct_timings <- lapply(names(ct_processed), function(key) {
          extract_ct_timing(
            ct_processed[[key]],
            paste0(ds, "/", key, " (ct_col=", ct_col_name, ")")
          )
        })
        names(ct_timings) <- names(ct_processed)
        if (any(vapply(ct_timings, function(value) is.null(value), logical(1)))) {
          stop(
            "CT processor did not return per-variant timing metadata for ",
            ct_col_name
          )
        }
        schema2_flags <- vapply(
          ct_timings, function(value) isTRUE(value[["schema2"]]), logical(1)
        )
        if (any(schema2_flags) && !all(schema2_flags)) {
          stop("CT timing metadata mixes legacy and schema-2 results for ",
               ct_col_name)
        }
        if (all(schema2_flags)) {
          shared_times <- unique(vapply(
            ct_timings,
            function(value) as.numeric(value[["shared"]]),
            numeric(1)
          ))
          if (length(shared_times) != 1L ||
              !is.finite(shared_times[[1L]]) || shared_times[[1L]] < 0) {
            stop("CT shared timing differs across variants for ", ct_col_name)
          }
          shared_mem <- ct_timings[[1L]][["mem"]]
          if (!all(vapply(
            ct_timings,
            function(value) identical(value[["mem"]], shared_mem),
            logical(1)
          ))) {
            stop("CT shared memory differs across variants for ", ct_col_name)
          }
          timing_ids <- unique(vapply(
            ct_timings, function(value) as.character(value[["id"]]),
            character(1)
          ))
          if (length(timing_ids) != 1L) {
            stop("CT variants do not share one timing_id for ", ct_col_name)
          }
          record_ct_shared(
            list(
              shared_time_secs = shared_times[[1L]],
              shared_mem_GB = shared_mem,
              timing_id = timing_ids[[1L]],
              timing_schema = 2L,
              shared_timing_method = ct_shared_timing_method(ct_col_name)
            ),
            paste0(ds, " (ct_col=", ct_col_name, ")")
          )
        }
        for (nm in pending_names) {
          key <- paste0("hvg", as.integer(ct_combos[[nm]]$hvg))
          res <- ct_processed[[key]]
          if (is.null(res)) {
            stop("CT processor returned no result for ", key,
                 " (ct_col=", ct_col_name, ")")
          }
          timing <- ct_timings[[key]]
          if (isTRUE(timing[["schema2"]]) &&
              !"timing_schema" %in% names(res)) {
            res[["shared_time_secs"]] <- timing[["shared"]]
            res[["variant_time_secs"]] <- timing[["local"]]
            res[["shared_mem_GB"]] <- timing[["mem"]]
            res[["timing_id"]] <- timing[["id"]]
            res[["timing_schema"]] <- 2L
          }
          if (isTRUE(timing[["schema2"]])) {
            res[["shared_timing_method"]] <-
              ct_shared_timing_method(ct_col_name)
          }
          res[["exec_time"]] <- as.numeric(timing[["local"]])
          res[["mem_GB"]] <- peak_rss_gb()
          if (isTRUE(timing[["schema2"]])) {
            validate_hpc_timing_bundle(
              res, label = paste0("Pseudobulk result ", ds, "/", nm)
            )
          }
          bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
          save_rds_atomic(res, bundle_file)
          log_exec_row(ds, nm, res[["exec_time"]], log_file,
                       mem_gb = res[["mem_GB"]])
          results[[nm]] <- res
        }
      } else {
        if (is.null(seurat)) {
          stop("Canonical CT pseudobulk requires an H5AD path.")
        }
        for (nm in pending_names) {
          spec <- ct_combos[[nm]]
          bundle_file <- file.path(results_dir, paste0(ds, "_", nm, ".rds"))
          ct_time <- exec_time(
            res <- process_pseudobulk_ct_fig(
              seurat,
              labels,
              ct_col = as.character(spec$ct_col),
              sample_col = sample_col,
              hvg = spec$hvg
            )
          )
          res[["exec_time"]] <- as.numeric(ct_time, units = "secs")
          res[["mem_GB"]] <- peak_rss_gb()
          validate_hpc_timing_bundle(
            res, label = paste0("Pseudobulk result ", ds, "/", nm)
          )
          save_rds_atomic(res, bundle_file)
          log_exec_row(ds, nm, res[["exec_time"]], log_file,
                       mem_gb = res[["mem_GB"]])
          results[[nm]] <- res
        }
      }
    }
  }
  for (timing in ct_shared_replay) {
    log_exec_row(
      ds, timing$method, timing$time, log_file,
      mem_gb = timing$mem
    )
  }
  # Restore the historical CT combination order after grouping by cell-type
  # column; this keeps the published result-list order stable.
  results <- c(
    results[setdiff(names(results), ct_names)],
    results[ct_names[ct_names %in% names(results)]]
  )
  results
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
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("scITD result ", ds, "/", nm)
      )
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

    # Batch composition requires only the configured author annotation and
    # selected Leiden resolution. Existing HiTME/scATOMIC bundles are legacy
    # extras and remain valid, but are not regenerated or required here.
    add_coda("ECODA_authors_HR", ct_col_high_res)
    add_coda("ECODA_authors_HR_NULL", ct_col_high_res, shuffle_labels = TRUE)
    add_coda("ECODA_seuratres_2", "RNA_snn_res.2")

    results <- list()
    for (nm in names(combos)) {
      bundle_file <- file.path(
        results_dir, paste0(artifact_stem, "_", nm, ".rds")
      )
      if (artifact_checksum_ok(bundle_file) && !force) {
        cached <- read_rds_checked(bundle_file)
        validate_hpc_timing_bundle(
          cached, label = paste0("Composition result ", ds, "/", nm)
        )
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
  deconv_pb_timing <- NULL
  if ("timing_schema" %in% names(pb_hvg2000)) {
    validate_hpc_timing_bundle(
      pb_hvg2000, label = paste0("Pseudobulk cache ", ds, "/hvg2000")
    )
    deconv_pb_timing <- list(
      variant = as.numeric(pb_hvg2000[["variant_time_secs"]]),
      shared = as.numeric(pb_hvg2000[["shared_time_secs"]]),
      shared_mem = pb_hvg2000[["shared_mem_GB"]],
      timing_id = as.character(pb_hvg2000[["timing_id"]])
    )
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
      cached <- read_rds_checked(bundle_file)
      validate_hpc_timing_bundle(
        cached, label = paste0("Composition result ", ds, "/", nm)
      )
      results[[nm]] <- cached
      if (!is.null(cached$exec_time)) {
        log_exec_row(ds, nm, cached$exec_time, log_file,
                     mem_gb = cached$mem_GB)
      }
      next
    }
    time_secs <- exec_time(res <- combos[[nm]]())
    process_time <- as.numeric(time_secs, units = "secs")
    if (identical(nm, "ECODA_deconv") && !is.null(deconv_pb_timing)) {
      res[["exec_time"]] <- process_time + deconv_pb_timing[["variant"]]
      res[["variant_time_secs"]] <- res[["exec_time"]]
      res[["shared_time_secs"]] <- deconv_pb_timing[["shared"]]
      res[["shared_mem_GB"]] <- deconv_pb_timing[["shared_mem"]]
      res[["timing_id"]] <- deconv_pb_timing[["timing_id"]]
      res[["timing_schema"]] <- 2L
    } else {
      res[["exec_time"]] <- process_time
    }
    res[["mem_GB"]] <- peak_rss_gb()
    validate_hpc_timing_bundle(
      res, label = paste0("Composition result ", ds, "/", nm)
    )
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
