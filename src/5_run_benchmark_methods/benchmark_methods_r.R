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
    labels <- shuffle_labels_deterministic(labels)
  }
  res <- create_result_bundle(
    feat_mat,
    labels,
    dist_mat = dist_mat,
    extra = list(counts = df_imp)
  )
  return(res)
}

# Keep null controls deterministic without changing the feature matrix.  The
# caller supplies the original sample order; only the label values are
# permuted and their sample names remain attached to the same samples.
shuffle_labels_deterministic <- function(labels, seed = 123) {
  label_ids <- names(labels)
  set.seed(seed)
  labels <- labels[sample(seq_along(labels))]
  names(labels) <- label_ids
  labels
}

# Reuse a computed result bundle for a deterministic null control.  In
# particular, this must not re-run composition or batch correction.
shuffle_result_labels_deterministic <- function(res, seed = 123) {
  if (!is.list(res) || is.null(res[["labels"]])) {
    stop("shuffle_result_labels_deterministic: result labels are required")
  }
  res[["labels"]] <- shuffle_labels_deterministic(res[["labels"]], seed = seed)
  if (!is.null(res[["dist_mat"]])) {
    res[["scores"]] <- calc_sep_score(res[["dist_mat"]], res[["labels"]])
  }
  res
}

# Load the corrected-mode contract lazily.  Ordinary benchmark and
# uncorrected batch paths do not need this module and retain their existing
# package/source behavior.
.ecoda_require_batch_contract <- function() {
  required <- c(
    "ecoda_batch_normalize_keys",
    "ecoda_batch_validate_metadata"
  )
  if (all(vapply(
    required,
    function(name) exists(name, mode = "function", inherits = TRUE),
    logical(1)
  ))) {
    return(invisible(TRUE))
  }

  project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
  candidates <- character()
  if (nzchar(project_root)) {
    candidates <- c(
      candidates,
      file.path(project_root, "src", "utils", "batch_contract.R")
    )
  }
  candidates <- c(
    candidates,
    file.path(getwd(), "src", "utils", "batch_contract.R")
  )
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) {
    stop(
      "Corrected CLR composition requires src/utils/batch_contract.R"
    )
  }
  source(candidates[[1L]], local = .GlobalEnv)
  if (!all(vapply(
    required,
    function(name) exists(name, mode = "function", inherits = TRUE),
    logical(1)
  ))) {
    stop("src/utils/batch_contract.R did not define the corrected-mode contract")
  }
  invisible(TRUE)
}


# Correct CLR composition by subtracting only the fitted technical batch
# random effect. Biological labels are deliberately absent from this model.
correct_clr_batch_lmm <- function(
  feat_mat,
  sample_meta,
  batch_col = NULL,
  sample_col = "Sample",
  metadata_validation = NULL
) {
  if (is.null(batch_col)) {
    stop("correct_clr_batch_lmm: batch_col is required")
  }
  .ecoda_require_batch_contract()
  if (!is.data.frame(sample_meta)) {
    stop("correct_clr_batch_lmm: sample metadata must be a data.frame")
  }

  batch_keys <- ecoda_batch_normalize_keys(
    batch_col,
    sample_col = sample_col
  )
  feat_mat <- as.matrix(feat_mat)
  if (length(dim(feat_mat)) != 2L ||
      nrow(feat_mat) == 0L || ncol(feat_mat) == 0L) {
    stop("correct_clr_batch_lmm: feature matrix must be non-empty")
  }
  if (is.null(rownames(feat_mat)) || anyDuplicated(rownames(feat_mat))) {
    stop("correct_clr_batch_lmm: feature matrix needs unique sample rownames")
  }
  if (!is.numeric(feat_mat) || any(!is.finite(feat_mat))) {
    stop("correct_clr_batch_lmm: feature matrix must contain finite numeric values")
  }

  validation <- metadata_validation
  if (is.null(validation)) {
    # This call deliberately receives the cell-level table.  The contract
    # checks every row before creating its sample-level representation.
    validation <- ecoda_batch_validate_metadata(
      metadata = sample_meta,
      batch_keys = if (length(batch_keys) >= 2L) {
        as.list(unname(batch_keys))
      } else {
        batch_keys
      },
      sample_col = sample_col
    )
  }
  if (!is.list(validation) || !isTRUE(validation[["valid"]]) ||
      !identical(as.character(validation[["ordered_keys"]]), batch_keys)) {
    stop("correct_clr_batch_lmm: invalid batch metadata validation")
  }
  sample_metadata <- validation[["sample_metadata"]]
  canonical_metadata <- validation[["canonical_sample_metadata"]]
  if (!is.data.frame(sample_metadata) ||
      !is.data.frame(canonical_metadata) ||
      !all(batch_keys %in% colnames(sample_metadata)) ||
      !all(batch_keys %in% colnames(canonical_metadata))) {
    stop("correct_clr_batch_lmm: validated batch metadata is incomplete")
  }

  sample_ids <- as.character(validation[["sample_ids"]])
  if (length(sample_ids) != nrow(feat_mat) ||
      anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      !identical(sample_ids, rownames(feat_mat))) {
    stop("correct_clr_batch_lmm: sample-order mismatch")
  }
  if (nrow(sample_metadata) != nrow(feat_mat) ||
      nrow(canonical_metadata) != nrow(feat_mat)) {
    stop("correct_clr_batch_lmm: sample metadata does not cover all features")
  }

  # Never interpolate user-configured names into a formula.  The aliases are
  # fixed syntactic names and the original key order is retained separately.
  aliases <- if (length(batch_keys) == 1L) {
    "batch"
  } else {
    paste0("batch_key_", seq_along(batch_keys))
  }
  model_data <- data.frame(
    row.names = seq_len(nrow(feat_mat)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (key_index in seq_along(batch_keys)) {
    key <- batch_keys[[key_index]]
    values <- canonical_metadata[[key]]
    levels <- validation[["per_key_levels"]][[key]]
    if (length(values) != nrow(feat_mat) ||
        length(levels) < 2L ||
        anyNA(values) || any(!values %in% levels)) {
      stop("correct_clr_batch_lmm: invalid levels for batch key ", key)
    }
    model_data[[aliases[[key_index]]]] <- factor(
      as.character(values),
      levels = as.character(levels)
    )
  }

  # Keep the scalar formula literal for the established one-key path.  The
  # multi-key branch is generated only from the fixed aliases above.
  model_formula <- if (length(aliases) == 1L) {
    y ~ 1 + (1 | batch)
  } else {
    stats::as.formula(paste0(
      "y ~ 1 + ",
      paste0("(1 | ", aliases, ")", collapse = " + ")
    ))
  }

  corrected <- feat_mat
  for (feature in seq_len(ncol(feat_mat))) {
    model_data$y <- as.numeric(feat_mat[, feature])
    feature_name <- if (is.null(colnames(feat_mat))) {
      as.character(feature)
    } else {
      colnames(feat_mat)[[feature]]
    }
    fit <- tryCatch(
      lme4::lmer(model_formula, data = model_data, REML = TRUE),
      error = function(e) {
        stop("correct_clr_batch_lmm: nonconvergence for feature ",
             feature_name, ": ", conditionMessage(e))
      }
    )
    convergence <- fit@optinfo$conv$lme4$messages
    if (!is.null(convergence)) {
      stop("correct_clr_batch_lmm: nonconvergence for feature ",
           feature_name, ": ",
           paste(convergence, collapse = "; "))
    }

    random_effects <- numeric(nrow(feat_mat))
    random_tables <- lme4::ranef(fit)
    for (alias in aliases) {
      random_table <- random_tables[[alias]]
      if (is.null(random_table) ||
          !"(Intercept)" %in% colnames(random_table)) {
        stop("correct_clr_batch_lmm: nonconvergence for feature ",
             feature_name, ": missing random effect for ", alias)
      }
      effect_levels <- rownames(random_table)
      effect <- random_table[["(Intercept)"]][match(
        as.character(model_data[[alias]]), effect_levels
      )]
      if (length(effect) != nrow(feat_mat) ||
          anyNA(effect) || any(!is.finite(effect))) {
        stop("correct_clr_batch_lmm: nonconvergence for feature ",
             feature_name, ": invalid random effect for ", alias)
      }
      random_effects <- random_effects + as.numeric(effect)
    }
    corrected[, feature] <- as.numeric(model_data$y) - random_effects
  }

  # Recenter rows (not columns) and make the final coordinate the exact
  # negative sum of its prefix, restoring the CLR zero-sum invariant.
  corrected <- sweep(corrected, 1L, rowMeans(corrected), FUN = "-")
  final_column <- ncol(corrected)
  if (final_column == 1L) {
    corrected[, final_column] <- 0
  } else {
    prefix <- seq_len(final_column - 1L)
    for (row_index in seq_len(nrow(corrected))) {
      corrected[row_index, final_column] <- -sum(
        corrected[row_index, prefix, drop = FALSE]
      )
    }
  }
  dimnames(corrected) <- dimnames(feat_mat)
  row_sums <- rowSums(corrected)
  if (any(!is.finite(row_sums)) || any(row_sums != 0)) {
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
if (!exists("DEFAULT_CHUNK_SIZE", inherits = TRUE)) {
  DEFAULT_CHUNK_SIZE <- 4096L
}
# Canonical cell-type pseudobulk processing from a raw H5AD CSR store.
#
# Unlike process_pseudobulk_ct_fig(), this boundary never accepts a Seurat
# object.  The Python helper performs the metadata-only eligibility pass and
# the single raw counts pass, retaining only sparse Sample x cell-type
# contribution rows in a run-owned store.  One cell type is materialized at a
# time at the direct DESeq2 boundary.
process_pseudobulk_ct_h5ad_fig <- function(
  h5ad_path,
  labels,
  sample_col = "Sample",
  ct_col,
  hvg = 500,
  min_cells = 5,
  chunk_size = DEFAULT_CHUNK_SIZE,
  run_id = NULL,
  temp_root = NULL,
  source_identity = NULL,
  timing_id = NULL
) {
  variants <- process_pseudobulk_ct_h5ad_variants_fig(
    h5ad_path = h5ad_path,
    labels = labels,
    sample_col = sample_col,
    ct_col = ct_col,
    hvgs = hvg,
    min_cells = min_cells,
    chunk_size = chunk_size,
    run_id = run_id,
    temp_root = temp_root,
    source_identity = source_identity,
    timing_id = timing_id
  )
  if (length(hvg) == 1L) variants[[1L]] else variants
}


# Internal helpers for the canonical H5AD cell-type store boundary.  These
# deliberately live outside the public wrapper so one prepared store can serve
# all requested HVG variants without a second raw H5AD pass.
.ecoda_ct_scalar_character <- function(value, name) {
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !nzchar(value)) {
    stop(name, " must be one non-empty string")
  }
  value
}

.ecoda_ct_scalar_positive_integer <- function(value, name) {
  if (length(value) != 1L || is.na(value) ||
      !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value <= 0 ||
      value > .Machine$integer.max) {
    stop(name, " must be a positive integer")
  }
  as.integer(value)
}

.ecoda_ct_canonical_path <- function(value, name) {
  value <- .ecoda_ct_scalar_character(as.character(value), name)
  path <- path.expand(value)
  initial <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (length(initial) != 1L || is.na(initial) ||
      !grepl("^(/|[A-Za-z]:[/\\\\])", initial, perl = TRUE)) {
    stop(name, " must resolve to an absolute path")
  }

  # normalizePath(mustWork = FALSE) cannot resolve symlinks above a
  # nonexistent descendant.  Canonicalize the deepest existing ancestor
  # first, then attach unresolved components to that canonical parent.
  candidate <- path
  suffix <- character()
  repeat {
    if (file.exists(candidate) || dir.exists(candidate)) {
      current <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      if (!length(suffix)) {
        path <- current
        break
      }

      # Walk the suffix from the canonical parent so existing symlinks are
      # resolved before any following .. component is applied.  Once a
      # component is missing, descendants remain unresolved until a .. pops
      # back to the existing prefix.
      unresolved <- character()
      for (component in suffix) {
        if (!nzchar(component) || identical(component, ".")) next
        if (identical(component, "..")) {
          if (length(unresolved)) {
            unresolved <- unresolved[-length(unresolved)]
          } else {
            current <- dirname(current)
          }
          next
        }
        if (length(unresolved)) {
          unresolved <- c(unresolved, component)
          next
        }
        next_path <- file.path(current, component)
        if (file.exists(next_path) || dir.exists(next_path)) {
          current <- normalizePath(next_path, winslash = "/", mustWork = TRUE)
        } else {
          unresolved <- component
        }
      }
      if (!length(unresolved)) {
        path <- current
      } else {
        rebuilt <- do.call(file.path, c(list(current), as.list(unresolved)))
        path <- normalizePath(rebuilt, winslash = "/", mustWork = FALSE)
      }
      break
    }

    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    suffix <- c(basename(candidate), suffix)
    candidate <- parent
  }
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path, perl = TRUE)) {
    stop(name, " must resolve to an absolute path")
  }
  path
}

.ecoda_ct_scalar_path <- function(value, name) {
  .ecoda_ct_canonical_path(value, name)
}

.ecoda_ct_paths_equal <- function(left, right,
                                  left_name = "left path",
                                  right_name = "right path") {
  identical(
    .ecoda_ct_canonical_path(left, left_name),
    .ecoda_ct_canonical_path(right, right_name)
  )
}

.ecoda_ct_py_raw_member <- function(object, name) {
  value <- tryCatch(object[[name]], error = function(e) NULL)
  if (is.null(value)) {
    value <- switch(
      name,
      shape = tryCatch(object$shape, error = function(e) NULL),
      data = tryCatch(object$data, error = function(e) NULL),
      indices = tryCatch(object$indices, error = function(e) NULL),
      indptr = tryCatch(object$indptr, error = function(e) NULL),
      NULL
    )
  }
  if (is.null(value) && exists(
    "py_get_attr",
    envir = asNamespace("reticulate"),
    inherits = FALSE
  )) {
    value <- tryCatch(
      reticulate::py_get_attr(object, name),
      error = function(e) NULL
    )
  }
  value
}

.ecoda_ct_py_member <- function(object, name) {
  value <- .ecoda_ct_py_raw_member(object, name)
  if (is.null(value)) return(NULL)
  if (!requireNamespace("reticulate", quietly = TRUE)) return(value)
  reticulate::py_to_r(value)
}

.ecoda_ct_csr_to_matrix <- function(csr, expected_nrow, expected_ncol) {
  if (inherits(csr, "Matrix") || is.matrix(csr)) {
    return(csr)
  }
  shape <- .ecoda_ct_py_member(csr, "shape")
  data <- .ecoda_ct_py_member(csr, "data")
  indices <- .ecoda_ct_py_member(csr, "indices")
  indptr <- .ecoda_ct_py_member(csr, "indptr")
  if (is.null(shape) || is.null(data) || is.null(indices) ||
      is.null(indptr)) {
    stop("Python CT store reader returned no usable CSR arrays")
  }
  shape <- as.numeric(shape)
  if (length(shape) != 2L || any(!is.finite(shape)) ||
      any(shape != floor(shape)) ||
      shape[1L] != expected_nrow || shape[2L] != expected_ncol) {
    stop("selected CT CSR shape is inconsistent")
  }
  data <- as.numeric(data)
  indices <- as.numeric(indices)
  indptr <- as.numeric(indptr)
  if (any(!is.finite(data)) || any(data < 0) ||
      any(data != floor(data)) ||
      any(data > .Machine$integer.max) ||
      any(!is.finite(indices)) || any(indices != floor(indices)) ||
      any(indices < 0) || any(indices >= expected_ncol) ||
      any(!is.finite(indptr)) || any(indptr != floor(indptr)) ||
      any(indptr < 0) || any(indptr > .Machine$integer.max) ||
      length(indptr) != expected_nrow + 1L ||
      indptr[1L] != 0 || any(diff(indptr) < 0) ||
      indptr[length(indptr)] != length(data) ||
      length(indices) != length(data)) {
    stop("selected CT CSR arrays are invalid")
  }
  if (length(data) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(), j = integer(), x = numeric(),
      dims = c(expected_nrow, expected_ncol)
    ))
  }
  row_ids <- rep.int(seq_len(expected_nrow), diff(indptr))
  if (anyDuplicated(paste(row_ids, indices, sep = ":"))) {
    stop("selected CT CSR rows contain duplicate coordinates")
  }
  Matrix::sparseMatrix(
    i = as.integer(row_ids),
    j = as.integer(indices + 1),
    x = data,
    dims = c(expected_nrow, expected_ncol),
    repr = "C"
  )
}

.ecoda_ct_audit_field <- function(audit, name) {
  if (is.null(audit) || !is.list(audit) ||
      is.null(names(audit)) || !name %in% names(audit)) {
    return(NULL)
  }
  audit[[name]]
}

prepare_pseudobulk_ct_h5ad_store_context <- function(
  h5ad_path,
  sample_col = "Sample",
  ct_col,
  chunk_size = DEFAULT_CHUNK_SIZE,
  run_id = NULL,
  temp_root = NULL,
  source_identity = NULL
) {
  h5ad_path <- .ecoda_ct_scalar_path(h5ad_path, "h5ad_path")
  if (!file.exists(h5ad_path) ||
      !isTRUE(file.info(h5ad_path)$size > 0)) {
    stop("H5AD path is missing or empty: ", h5ad_path)
  }
  sample_col <- .ecoda_ct_scalar_character(sample_col, "sample_col")
  ct_col <- .ecoda_ct_scalar_character(ct_col, "ct_col")
  chunk_size <- .ecoda_ct_scalar_positive_integer(chunk_size, "chunk_size")
  if (is.character(source_identity) &&
      (length(source_identity) != 1L ||
       is.na(source_identity) || !nzchar(source_identity))) {
    stop("source_identity character value must be one non-empty string")
  }

  resolved_run_id <- run_id
  if (is.null(resolved_run_id) || length(resolved_run_id) == 0L ||
      !nzchar(as.character(resolved_run_id))) {
    resolved_run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  }
  if (!nzchar(resolved_run_id)) {
    resolved_run_id <- paste0(
      "local-", Sys.getpid(), "-", as.integer(as.numeric(Sys.time()))
    )
  }
  resolved_run_id <- .ecoda_ct_scalar_character(
    as.character(resolved_run_id),
    "run_id"
  )
  if (!grepl(
    "^[A-Za-z0-9][A-Za-z0-9_-]*$",
    resolved_run_id,
    perl = TRUE
  )) {
    stop("run_id contains unsafe path characters")
  }

  max_age_seconds <- suppressWarnings(as.numeric(Sys.getenv(
    "ECODA_PSEUDOBULK_CT_MAX_AGE_SECONDS", unset = "86400"
  )))
  if (length(max_age_seconds) != 1L || is.na(max_age_seconds) ||
      !is.finite(max_age_seconds) || max_age_seconds < 0) {
    stop("ECODA_PSEUDOBULK_CT_MAX_AGE_SECONDS must be nonnegative and finite")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required for canonical H5AD cell-type pseudobulk")
  }
  if (!exists("get_pb_deseq2_from_counts", mode = "function")) {
    stop("get_pb_deseq2_from_counts is required for canonical H5AD cell-type pseudobulk")
  }
  project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (!nzchar(project_root)) {
    stop("PROJECT_ROOT is required for canonical H5AD cell-type pseudobulk")
  }
  module_dir <- normalizePath(
    file.path(project_root, "src", "utils", "py"),
    mustWork = TRUE
  )
  python_sys <- reticulate::import("sys", convert = FALSE)
  python_sys$path$insert(0L, module_dir)
  loader <- reticulate::import_from_path(
    "h5ad_pseudobulk",
    path = module_dir,
    convert = FALSE
  )
  if (is.null(loader$prepare_h5ad_ct_group_store) ||
      is.null(loader$read_h5ad_ct_group_store) ||
      is.null(loader$audit_h5ad_ct_group_store)) {
    stop("h5ad_pseudobulk is missing the canonical cell-type store API")
  }

  if (is.null(temp_root) || length(temp_root) == 0L ||
      !nzchar(as.character(temp_root))) {
    configured_run_root <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
    if (nzchar(configured_run_root)) {
      temp_root <- file.path(configured_run_root, "scratch")
    } else {
      temp_root <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
      if (!nzchar(temp_root)) temp_root <- tempdir()
    }
  }
  temp_root <- .ecoda_ct_scalar_path(temp_root, "temp_root")
  if (!dir.exists(temp_root) &&
      !dir.create(temp_root, recursive = TRUE, showWarnings = FALSE)) {
    stop("cannot create temporary root: ", temp_root)
  }
  ct_parent <- file.path(temp_root, "pseudobulk_ct", resolved_run_id)
  if (!dir.exists(ct_parent) &&
      !dir.create(ct_parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("cannot create run-owned CT root: ", ct_parent)
  }
  if (nzchar(Sys.readlink(ct_parent))) {
    stop("run-owned CT root must not be a symlink")
  }
  unique_root <- tempfile(
    pattern = ".store-",
    tmpdir = ct_parent,
    fileext = ""
  )
  if (file.exists(unique_root) || dir.exists(unique_root) ||
      !dir.create(unique_root, recursive = FALSE, showWarnings = FALSE)) {
    stop("cannot acquire a unique run-owned CT store root")
  }
  store_path <- file.path(unique_root, "groups.h5")
  manifest_path <- paste0(store_path, ".manifest.json")
  lock_path <- paste0(store_path, ".lock")

  # Keep the Python ownership context available for the failure cleanup below.
  # A writing-state store is removable only while this exact process is still
  # inside the Python preparation call.  Once preparation returns, normal
  # ready-state cleanup remains the only permitted path.
  preparation_state <- "not_started"
  expected_source_identity <- source_identity
  if (is.null(expected_source_identity)) {
    expected_source_identity <- tryCatch({
      source_identity_fn <- reticulate::py_get_attr(loader, "_source_identity")
      reticulate::py_to_r(source_identity_fn(h5ad_path, NULL))
    }, error = function(error) {
      NULL
    })
  }
  scheduler_keys <- c(
    "SLURM_JOB_ID",
    "SLURM_STEP_ID",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "PBS_JOBID",
    "JOB_ID",
    "ECODA_JOB_ID"
  )
  scheduler_values <- Sys.getenv(scheduler_keys, unset = "")
  scheduler_values <- scheduler_values[nzchar(scheduler_values)]
  expected_scheduler_identity <- structure(
    as.list(scheduler_values),
    names = names(scheduler_values)
  )
  expected_owner_host <- tryCatch({
    socket <- reticulate::import("socket", convert = TRUE)
    as.character(socket$gethostname())
  }, error = function(error) {
    as.character(Sys.info()[["nodename"]])
  })

  .ecoda_ct_json_equal <- function(left, right) {
    if (is.null(left) || is.null(right)) {
      return(is.null(left) && is.null(right))
    }
    if (is.factor(left)) left <- as.character(left)
    if (is.factor(right)) right <- as.character(right)
    if (is.atomic(left) && is.atomic(right)) {
      if (length(left) != length(right)) return(FALSE)
      if (is.numeric(left) && is.numeric(right)) {
        if (anyNA(left) || anyNA(right) ||
            any(!is.finite(left)) || any(!is.finite(right))) {
          return(FALSE)
        }
        return(all(left == right))
      }
      return(identical(left, right))
    }
    if (!is.list(left) || !is.list(right) ||
        length(left) != length(right)) {
      return(FALSE)
    }
    left_names <- names(left)
    right_names <- names(right)
    if (is.null(left_names) != is.null(right_names)) return(FALSE)
    if (!is.null(left_names)) {
      if (anyDuplicated(left_names) || anyDuplicated(right_names) ||
          !setequal(left_names, right_names)) {
        return(FALSE)
      }
      return(all(vapply(
        left_names,
        function(name) .ecoda_ct_json_equal(left[[name]], right[[name]]),
        logical(1)
      )))
    }
    all(vapply(
      seq_along(left),
      function(index) .ecoda_ct_json_equal(left[[index]], right[[index]]),
      logical(1)
    ))
  }

  .ecoda_ct_path_is_symlink <- function(path) {
    link <- tryCatch(Sys.readlink(path), error = function(error) NA_character_)
    length(link) == 1L && !is.na(link) && nzchar(link)
  }

  .ecoda_ct_existing_path <- function(path) {
    file.exists(path) || dir.exists(path)
  }

  .ecoda_ct_read_manifest <- function() {
    if (!isTRUE(file_test("-f", manifest_path)) ||
        .ecoda_ct_path_is_symlink(manifest_path) ||
        !requireNamespace("jsonlite", quietly = TRUE)) {
      return(NULL)
    }
    info <- tryCatch(file.info(manifest_path), error = function(error) NULL)
    if (is.null(info) || nrow(info) != 1L || isTRUE(info$isdir) ||
        is.na(info$size) || info$size <= 0) {
      return(NULL)
    }
    manifest <- tryCatch(
      jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
      error = function(error) NULL
    )
    if (!is.list(manifest) || is.null(names(manifest))) return(NULL)
    manifest
  }

  .ecoda_ct_manifest_matches <- function(manifest, expected_state,
                                          require_fresh = FALSE) {
    required <- c(
      "schema", "stage", "state", "run_id", "pid", "owner_host",
      "scheduler_identity", "source_identity", "source_checksum",
      "store_path", "lock_path", "created_at"
    )
    if (!is.list(manifest) || is.null(names(manifest)) ||
        !all(required %in% names(manifest))) {
      return(FALSE)
    }
    schema <- suppressWarnings(as.numeric(manifest[["schema"]]))
    pid <- suppressWarnings(as.numeric(manifest[["pid"]]))
    created_at <- suppressWarnings(as.numeric(manifest[["created_at"]]))
    state <- as.character(manifest[["state"]])
    run_value <- as.character(manifest[["run_id"]])
    host_value <- as.character(manifest[["owner_host"]])
    manifest_store <- as.character(manifest[["store_path"]])
    manifest_lock <- as.character(manifest[["lock_path"]])
    if (length(schema) != 1L || is.na(schema) || !is.finite(schema) ||
        schema != 1 || length(pid) != 1L || is.na(pid) ||
        !is.finite(pid) || pid != floor(pid) ||
        length(created_at) != 1L || is.na(created_at) ||
        !is.finite(created_at) || created_at < 0 ||
        length(state) != 1L || is.na(state) ||
        length(run_value) != 1L || is.na(run_value) ||
        length(host_value) != 1L || is.na(host_value) ||
        !nzchar(run_value) || !nzchar(host_value) ||
        length(manifest_store) != 1L || is.na(manifest_store) ||
        length(manifest_lock) != 1L || is.na(manifest_lock) ||
        !nzchar(manifest_store) || !nzchar(manifest_lock) ||
        !identical(as.character(manifest[["stage"]]), "pseudobulk_ct") ||
        !identical(state, expected_state) ||
        !identical(run_value, resolved_run_id) ||
        !identical(host_value, expected_owner_host) ||
        !identical(pid, as.numeric(Sys.getpid())) ||
        !identical(manifest_store, store_path) ||
        !.ecoda_ct_paths_equal(
          manifest_store,
          store_path,
          "manifest store path",
          "run-owned store path"
        ) ||
        !identical(manifest_lock, lock_path) ||
        !.ecoda_ct_paths_equal(
          manifest_lock,
          lock_path,
          "manifest lock path",
          "run-owned lock path"
        ) ||
        !is.list(manifest[["scheduler_identity"]]) ||
        !.ecoda_ct_json_equal(
          manifest[["scheduler_identity"]],
          expected_scheduler_identity
        ) ||
        is.null(expected_source_identity) ||
        !.ecoda_ct_json_equal(
          manifest[["source_identity"]],
          expected_source_identity
        ) ||
        .ecoda_ct_existing_path(lock_path)) {
      return(FALSE)
    }
    if (isTRUE(require_fresh)) {
      age <- as.numeric(Sys.time()) - created_at
      if (!is.finite(age) || age < 0 || age >= max_age_seconds) {
        return(FALSE)
      }
    }
    TRUE
  }

  cleanup_done <- FALSE
  cleanup_store <- function(strict = FALSE, allow_writing = FALSE) {
    reject <- function(message) {
      if (strict) stop("CT store cleanup failed closed: ", message)
      warning("CT store cleanup skipped: ", message)
      invisible(FALSE)
    }
    if (isTRUE(cleanup_done)) return(invisible(TRUE))
    if (isTRUE(allow_writing) &&
        !identical(preparation_state, "running")) {
      return(reject("writing-state cleanup is not tied to active preparation"))
    }
    if (!dir.exists(unique_root) ||
        .ecoda_ct_path_is_symlink(unique_root)) {
      return(reject("CT store cleanup refused a missing or symlinked store root"))
    }
    protected_parents <- c(temp_root, dirname(ct_parent), ct_parent)
    if (any(vapply(
      protected_parents,
      function(path) !dir.exists(path) || .ecoda_ct_path_is_symlink(path),
      logical(1)
    ))) {
      return(reject("CT store cleanup refused an unsafe parent path"))
    }

    entries <- tryCatch(
      list.files(unique_root, all.files = TRUE, no.. = TRUE),
      error = function(error) NULL
    )
    if (is.null(entries)) {
      return(reject("CT store cleanup could not inspect the run-owned root"))
    }
    allowed_entries <- c(basename(store_path), basename(manifest_path))
    if (length(entries) > 0L &&
        any(!entries %in% allowed_entries)) {
      return(reject("CT store cleanup found unexpected root entries"))
    }
    entry_paths <- file.path(unique_root, entries)
    if (length(entry_paths) > 0L &&
        any(vapply(entry_paths, .ecoda_ct_path_is_symlink, logical(1)))) {
      return(reject("CT store cleanup refused symlinked store entries"))
    }
    if (.ecoda_ct_existing_path(store_path) &&
        !isTRUE(file_test("-f", store_path))) {
      return(reject("CT store cleanup refused a non-regular store file"))
    }
    if (.ecoda_ct_existing_path(manifest_path) &&
        !isTRUE(file_test("-f", manifest_path))) {
      return(reject("CT store cleanup refused a non-regular manifest file"))
    }

    has_store <- .ecoda_ct_existing_path(store_path) ||
      .ecoda_ct_existing_path(manifest_path) ||
      .ecoda_ct_existing_path(lock_path)
    if (!has_store) {
      if (length(entries) == 0L) {
        unlink(unique_root, recursive = TRUE, force = FALSE)
        if (.ecoda_ct_existing_path(unique_root)) {
          return(reject("CT store cleanup did not remove the empty root"))
        }
        cleanup_done <<- TRUE
        return(invisible(TRUE))
      }
      return(reject("CT store cleanup found an unmanifested store"))
    }

    manifest <- .ecoda_ct_read_manifest()
    if (isTRUE(allow_writing)) {
      # Do not require a valid CSR layout here: an exception can leave a
      # torn HDF5 file after Python closes its handles.  The manifest is the
      # ownership proof, and every path/owner/source field is checked above.
      if (is.null(manifest) ||
          !.ecoda_ct_manifest_matches(
            manifest, expected_state = "writing", require_fresh = TRUE
          )) {
        return(reject(
          "writing-state manifest is unknown, stale, or not current-run-owned"
        ))
      }
      if (.ecoda_ct_existing_path(store_path) &&
          (dir.exists(store_path) || .ecoda_ct_path_is_symlink(store_path))) {
        return(reject("CT store cleanup refused an unsafe store file"))
      }
    } else {
      audit <- tryCatch(
        reticulate::py_to_r(loader$audit_h5ad_ct_group_store(
          store_path,
          expected_run_id = resolved_run_id,
          max_age_seconds = max_age_seconds,
          cleanup = FALSE
        )),
        error = function(error) error
      )
      if (inherits(audit, "error") || is.null(audit)) {
        message <- if (inherits(audit, "error")) {
          conditionMessage(audit)
        } else {
          "audit returned no result"
        }
        return(reject(paste("CT store cleanup audit failed:", message)))
      }
      valid <- .ecoda_ct_audit_field(audit, "valid")
      audited_path <- .ecoda_ct_audit_field(audit, "store_path")
      owner_status <- tolower(as.character(
        .ecoda_ct_audit_field(audit, "owner_status")
      ))
      expired <- .ecoda_ct_audit_field(audit, "expired")
      state <- tolower(as.character(.ecoda_ct_audit_field(audit, "state")))
      if (!isTRUE(as.logical(valid)) ||
          length(audited_path) != 1L ||
          !.ecoda_ct_paths_equal(
            as.character(audited_path),
            store_path,
            "audited store path",
            "run-owned store path"
          ) ||
          length(owner_status) != 1L || owner_status != "active" ||
          length(expired) != 1L || is.na(as.logical(expired)) ||
          length(state) != 1L || !identical(state, "ready") ||
          is.null(manifest) ||
          !.ecoda_ct_manifest_matches(
            manifest, expected_state = "ready", require_fresh = FALSE
          )) {
        message <- .ecoda_ct_audit_field(audit, "reason")
        if (is.null(message) || !nzchar(as.character(message))) {
          message <- "ownership manifest or store audit is not safe"
        }
        return(reject(as.character(message)))
      }
    }

    # Remove only the two expected files, then the exact empty run root.
    # Never recursively delete an unexpected or symlinked entry.
    unlink(c(store_path, manifest_path), recursive = FALSE, force = FALSE)
    if (.ecoda_ct_existing_path(store_path) ||
        .ecoda_ct_existing_path(manifest_path) ||
        .ecoda_ct_existing_path(lock_path)) {
      return(reject("CT store cleanup did not remove the owned files"))
    }
    unlink(unique_root, recursive = TRUE, force = FALSE)
    if (.ecoda_ct_existing_path(unique_root)) {
      return(reject("CT store cleanup did not remove the owned root"))
    }
    cleanup_done <<- TRUE
    invisible(TRUE)
  }
  ownership_transferred <- FALSE
  on.exit({
    if (!isTRUE(ownership_transferred) && !isTRUE(cleanup_done)) {
      tryCatch(
        cleanup_store(
          strict = FALSE,
          allow_writing = identical(preparation_state, "running")
        ),
        error = function(error) {
          warning("CT store cleanup failed closed: ", conditionMessage(error))
        }
      )
    }
  }, add = TRUE)

  pre_audit <- tryCatch(
    reticulate::py_to_r(loader$audit_h5ad_ct_group_store(
      store_path,
      expected_run_id = resolved_run_id,
      max_age_seconds = max_age_seconds,
      cleanup = TRUE
    )),
    error = function(error) error
  )
  if (inherits(pre_audit, "error") || is.null(pre_audit)) {
    message <- if (inherits(pre_audit, "error")) {
      conditionMessage(pre_audit)
    } else {
      "audit returned no result"
    }
    stop("CT store reuse audit failed closed: ", message)
  }
  pre_valid <- .ecoda_ct_audit_field(pre_audit, "valid")
  if (length(pre_valid) != 1L || is.na(pre_valid)) {
    stop("CT store reuse audit returned an invalid validity flag")
  }
  if (isTRUE(as.logical(pre_valid)) &&
      (.ecoda_ct_existing_path(store_path) ||
       .ecoda_ct_existing_path(manifest_path) ||
       .ecoda_ct_existing_path(lock_path))) {
    stop("CT store path remains active after stale-store audit")
  }
  if (!isTRUE(as.logical(pre_valid)) &&
      (.ecoda_ct_existing_path(store_path) ||
       .ecoda_ct_existing_path(manifest_path) ||
       .ecoda_ct_existing_path(lock_path))) {
    stop("CT store path remains ambiguous after stale-store audit")
  }
  pre_entries <- tryCatch(
    list.files(unique_root, all.files = TRUE, no.. = TRUE),
    error = function(error) NULL
  )
  if (is.null(pre_entries) || length(pre_entries) > 0L) {
    stop("CT store temporary root is not empty before preparation")
  }

  preparation_state <- "running"
  prepared <- loader$prepare_h5ad_ct_group_store(
    path = h5ad_path,
    sample_col = sample_col,
    cell_type_col = ct_col,
    metadata_columns = as.list(as.character(unique(c(sample_col, ct_col)))),
    chunk_size = as.integer(chunk_size),
    max_value = as.integer(.Machine$integer.max),
    store_path = store_path,
    run_id = resolved_run_id,
    source_identity = source_identity
  )
  preparation_state <- "succeeded"
  reported_store <- .ecoda_ct_py_member(prepared, "store_path")
  if (!is.null(reported_store) &&
      !.ecoda_ct_paths_equal(
        as.character(reported_store),
        store_path,
        "store_path",
        "run-owned store path"
      )) {
    stop("Python CT store returned a path different from the run-owned target")
  }
  group_ids <- as.character(.ecoda_ct_py_member(prepared, "group_ids"))
  sample_ids <- as.character(.ecoda_ct_py_member(prepared, "sample_ids"))
  cell_type_ids <- as.character(.ecoda_ct_py_member(prepared, "cell_type_ids"))
  group_cell_counts <- as.numeric(
    .ecoda_ct_py_member(prepared, "group_cell_counts")
  )
  gene_names <- as.character(.ecoda_ct_py_member(prepared, "gene_names"))
  n_vars <- as.integer(.ecoda_ct_py_member(prepared, "n_vars"))
  all_sample_ids <- as.character(
    .ecoda_ct_py_member(prepared, "all_sample_ids")
  )
  n_groups <- length(group_ids)
  if (n_groups == 0L || length(sample_ids) != n_groups ||
      length(cell_type_ids) != n_groups ||
      length(group_cell_counts) != n_groups ||
      length(all_sample_ids) == 0L ||
      length(n_vars) != 1L || is.na(n_vars) || n_vars <= 0L ||
      length(gene_names) != n_vars ||
      anyNA(group_ids) || any(!nzchar(group_ids)) ||
      anyDuplicated(group_ids) ||
      anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyNA(cell_type_ids) || any(!nzchar(cell_type_ids)) ||
      any(!is.finite(group_cell_counts)) ||
      any(group_cell_counts < 0) ||
      any(group_cell_counts != floor(group_cell_counts)) ||
      anyNA(all_sample_ids) || any(!nzchar(all_sample_ids)) ||
      anyDuplicated(all_sample_ids)) {
    stop("Python CT store metadata has invalid dimensions or identifiers")
  }
  if (anyNA(gene_names) || any(!nzchar(gene_names)) ||
      anyDuplicated(gene_names)) {
    stop("Python CT store gene names are invalid")
  }
  context <- list(
    loader = loader,
    store_path = store_path,
    manifest_path = manifest_path,
    unique_root = unique_root,
    run_id = resolved_run_id,
    max_age_seconds = max_age_seconds,
    group_ids = group_ids,
    sample_ids = sample_ids,
    all_sample_ids = all_sample_ids,
    cell_type_ids = cell_type_ids,
    group_cell_counts = group_cell_counts,
    gene_names = gene_names,
    n_vars = n_vars,
    cleanup = cleanup_store
  )
  ownership_transferred <- TRUE
  context
}

process_pseudobulk_ct_store_fig <- function(
  store_context,
  labels,
  sample_col = "Sample",
  ct_col,
  hvg = 500,
  min_cells = 5
) {
  if (!is.list(store_context) ||
      !is.function(store_context[["cleanup"]]) ||
      is.null(store_context[["loader"]]) ||
      is.null(store_context[["store_path"]])) {
    stop("process_pseudobulk_ct_store_fig received an invalid store context")
  }
  sample_col <- .ecoda_ct_scalar_character(sample_col, "sample_col")
  ct_col <- .ecoda_ct_scalar_character(ct_col, "ct_col")
  hvg <- .ecoda_ct_scalar_positive_integer(hvg, "hvg")
  min_cells <- .ecoda_ct_scalar_positive_integer(min_cells, "min_cells")
  group_ids <- as.character(store_context[["group_ids"]])
  sample_ids <- as.character(store_context[["sample_ids"]])
  all_sample_ids <- as.character(store_context[["all_sample_ids"]])
  cell_type_ids <- as.character(store_context[["cell_type_ids"]])
  group_cell_counts <- as.numeric(store_context[["group_cell_counts"]])
  gene_names <- as.character(store_context[["gene_names"]])
  n_vars <- as.integer(store_context[["n_vars"]])
  n_groups <- length(group_ids)
  if (n_groups == 0L || length(sample_ids) != n_groups ||
      length(cell_type_ids) != n_groups ||
      length(group_cell_counts) != n_groups ||
      length(all_sample_ids) == 0L ||
      length(n_vars) != 1L || is.na(n_vars) || n_vars <= 0L ||
      length(gene_names) != n_vars) {
    stop("process_pseudobulk_ct_store_fig received invalid store metadata")
  }
  all_samples <- sort(unique(all_sample_ids))
  canonical_samples <- all_sample_ids
  total_dist <- matrix(
    0,
    nrow = length(all_samples),
    ncol = length(all_samples),
    dimnames = list(all_samples, all_samples)
  )
  count_mat <- matrix(
    0L,
    nrow = length(all_samples),
    ncol = length(all_samples),
    dimnames = list(all_samples, all_samples)
  )
  cell_types <- unique(cell_type_ids)
  cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]
  successful_cts <- character()
  n_ct_pair_contributions <- 0L

  for (ct in cell_types) {
    eligible <- which(cell_type_ids == ct & group_cell_counts >= min_cells)
    if (length(eligible) < 2L) next
    ct_result <- tryCatch(
      (function() {
        eligible_group_ids <- group_ids[eligible]
        eligible_sample_ids <- sample_ids[eligible]
        sample_rank <- match(eligible_sample_ids, canonical_samples)
        if (anyNA(sample_rank) || anyDuplicated(eligible_sample_ids)) {
          stop("Python CT store has duplicate or unalignable eligible sample groups")
        }
        eligible_order <- order(sample_rank)
        eligible_group_ids <- eligible_group_ids[eligible_order]
        eligible_sample_ids <- eligible_sample_ids[eligible_order]
        selected <- store_context$loader$read_h5ad_ct_group_store(
          store_context$store_path,
          as.list(as.character(eligible_group_ids))
        )
        selected_group_ids <- as.character(
          .ecoda_ct_py_member(selected, "group_ids")
        )
        selected_sample_ids <- as.character(
          .ecoda_ct_py_member(selected, "sample_ids")
        )
        if (length(selected_group_ids) != length(eligible_group_ids) ||
            anyNA(selected_group_ids) || anyDuplicated(selected_group_ids)) {
          stop("selected CT contribution IDs are invalid")
        }
        selected_rows <- match(eligible_group_ids, selected_group_ids)
        if (anyNA(selected_rows)) {
          stop("selected CT contribution IDs do not match the requested groups")
        }
        selected_counts <- .ecoda_ct_py_raw_member(selected, "counts")
        if (is.null(selected_counts)) {
          selected_counts <- .ecoda_ct_py_raw_member(selected, "rows")
        }
        if (is.null(selected_counts)) {
          stop("Python CT store reader returned no sparse contribution rows")
        }
        selected_counts <- .ecoda_ct_csr_to_matrix(
          selected_counts,
          expected_nrow = length(selected_group_ids),
          expected_ncol = n_vars
        )
        if (length(dim(selected_counts)) != 2L ||
            nrow(selected_counts) != length(selected_group_ids) ||
            ncol(selected_counts) != n_vars) {
          stop("selected CT contribution rows have inconsistent dimensions")
        }
        selected_counts <- selected_counts[selected_rows, , drop = FALSE]
        if (!is.null(selected_sample_ids) &&
            (length(selected_sample_ids) != length(selected_group_ids) ||
             anyNA(selected_sample_ids) ||
             !identical(
               selected_sample_ids[selected_rows],
               eligible_sample_ids
             ))) {
          stop("selected CT sample IDs do not match canonical contribution order")
        }
        dimnames(selected_counts) <- list(eligible_group_ids, gene_names)
        counts_ct <- Matrix::t(selected_counts)
        dimnames(counts_ct) <- list(gene_names, eligible_sample_ids)
        metadata_ct <- data.frame(
          stringsAsFactors = FALSE,
          check.names = FALSE,
          sample_id = eligible_sample_ids
        )
        colnames(metadata_ct)[1L] <- sample_col
        if (!identical(sample_col, "Sample")) {
          metadata_ct[["Sample"]] <- eligible_sample_ids
        }
        rownames(metadata_ct) <- eligible_sample_ids
        norm <- get_pb_deseq2_from_counts(
          counts = counts_ct,
          metadata = metadata_ct,
          hvg = NULL,
          n_hvg = as.integer(hvg),
          black_list = "none",
          batch_col = NULL,
          blind = TRUE,
          correct_batch = FALSE
        )
        if (is.data.frame(norm)) norm <- as.matrix(norm)
        if (length(dim(norm)) != 2L ||
            nrow(norm) != length(eligible_sample_ids) ||
            is.null(rownames(norm)) ||
            anyNA(rownames(norm)) ||
            anyDuplicated(rownames(norm)) ||
            !setequal(rownames(norm), eligible_sample_ids)) {
          stop("normalized matrix has invalid sample identifiers")
        }
        norm <- norm[eligible_sample_ids, , drop = FALSE]
        dist_mat_ct <- as.matrix(stats::dist(norm))
        dimnames(dist_mat_ct) <- list(
          eligible_sample_ids,
          eligible_sample_ids
        )
        list(dist_mat = dist_mat_ct)
      })(),
      error = function(error) {
        warning(
          "process_pseudobulk_ct_h5ad_fig: pseudobulk failed for cell type '",
          ct, "': ", conditionMessage(error)
        )
        NULL
      }
    )
    if (is.null(ct_result)) {
      invisible(gc(verbose = FALSE))
      next
    }
    dist_mat_ct <- ct_result[["dist_mat"]]
    total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      total_dist[rownames(dist_mat_ct), colnames(dist_mat_ct)] + dist_mat_ct
    count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] <-
      count_mat[rownames(dist_mat_ct), colnames(dist_mat_ct)] + 1L
    successful_cts <- c(successful_cts, as.character(ct))
    n_ct_pair_contributions <- n_ct_pair_contributions +
      as.integer(sum(upper.tri(dist_mat_ct)))
    rm(ct_result, dist_mat_ct)
    invisible(gc(verbose = FALSE))
  }

  n_sample_pairs_contributed <- as.integer(sum(
    count_mat[upper.tri(count_mat)] > 0
  ))
  if (length(successful_cts) == 0L || n_sample_pairs_contributed == 0L) {
    stop(
      "process_pseudobulk_ct_h5ad_fig: no successful cell-type pseudobulks ",
      "contributed sample-pair distances for ct_col='", ct_col,
      "', hvg=", hvg, " (", length(cell_types), " cell types inspected)."
    )
  }
  final_dist_mat <- total_dist / count_mat
  final_dist_mat[is.nan(final_dist_mat)] <- 0
  create_result_bundle(
    feat_mat = final_dist_mat,
    labels,
    dist_mat = stats::as.dist(final_dist_mat),
    extra = list(
      n_ct_success = length(successful_cts),
      successful_cell_types = successful_cts,
      n_sample_pairs_contributed = n_sample_pairs_contributed,
      n_ct_pair_contributions = n_ct_pair_contributions
    )
  )
}

process_pseudobulk_ct_h5ad_variants_fig <- function(
  h5ad_path,
  labels,
  sample_col = "Sample",
  ct_col,
  hvgs = c(500L, 2000L),
  min_cells = 5,
  chunk_size = DEFAULT_CHUNK_SIZE,
  run_id = NULL,
  temp_root = NULL,
  source_identity = NULL,
  timing_id = NULL
) {
  if (!is.numeric(hvgs) || length(hvgs) == 0L ||
      anyNA(hvgs) || any(!is.finite(hvgs)) ||
      any(hvgs != floor(hvgs)) || any(hvgs <= 0) ||
      any(hvgs > .Machine$integer.max) ||
      anyDuplicated(as.integer(hvgs))) {
    stop("hvgs must be unique positive integers")
  }
  hvgs <- as.integer(hvgs)
  if (!is.null(timing_id)) {
    timing_id <- .ecoda_ct_scalar_character(timing_id, "timing_id")
  }
  shared_start <- proc.time()[["elapsed"]]
  context <- prepare_pseudobulk_ct_h5ad_store_context(
    h5ad_path = h5ad_path,
    sample_col = sample_col,
    ct_col = ct_col,
    chunk_size = chunk_size,
    run_id = run_id,
    temp_root = temp_root,
    source_identity = source_identity
  )
  shared_time_secs <- as.numeric(
    proc.time()[["elapsed"]] - shared_start
  )
  if (length(shared_time_secs) != 1L || is.na(shared_time_secs) ||
      !is.finite(shared_time_secs) || shared_time_secs < 0) {
    stop("CT shared preparation timing is invalid")
  }
  shared_mem_GB <- NA_real_
  if (exists("peak_rss_gb", mode = "function", inherits = TRUE)) {
    measured_mem <- tryCatch(
      as.numeric(peak_rss_gb()),
      error = function(error) NA_real_
    )
    if (length(measured_mem) == 1L && !is.nan(measured_mem) &&
        (is.na(measured_mem) ||
         (is.finite(measured_mem) && measured_mem >= 0))) {
      shared_mem_GB <- measured_mem
    }
  }
  cleanup_done <- FALSE
  on.exit({
    if (!isTRUE(cleanup_done)) {
      tryCatch(
        context$cleanup(strict = FALSE),
        error = function(error) {
          warning("CT store cleanup failed closed: ", conditionMessage(error))
        }
      )
    }
  }, add = TRUE)
  results <- lapply(hvgs, function(n_hvg) {
    variant_start <- proc.time()[["elapsed"]]
    result <- process_pseudobulk_ct_store_fig(
      store_context = context,
      labels = labels,
      sample_col = sample_col,
      ct_col = ct_col,
      hvg = n_hvg,
      min_cells = min_cells
    )
    variant_time_secs <- as.numeric(
      proc.time()[["elapsed"]] - variant_start
    )
    if (length(variant_time_secs) != 1L || is.na(variant_time_secs) ||
        !is.finite(variant_time_secs) || variant_time_secs < 0) {
      stop("CT variant processing timing is invalid")
    }
    timing_metadata <- list(
      shared_time_secs = shared_time_secs,
      variant_time_secs = variant_time_secs,
      shared_mem_GB = shared_mem_GB,
      timing_id = timing_id,
      timing_schema = if (is.null(timing_id)) NULL else 2L
    )
    if (!is.null(timing_id)) {
      result[["shared_time_secs"]] <- shared_time_secs
      result[["variant_time_secs"]] <- variant_time_secs
      result[["shared_mem_GB"]] <- shared_mem_GB
      result[["timing_id"]] <- timing_id
      result[["timing_schema"]] <- 2L
    }
    attr(result, "ct_timing") <- timing_metadata
    result
  })
  context$cleanup(strict = TRUE)
  cleanup_done <- TRUE
  names(results) <- paste0("hvg", hvgs)
  results
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

# Feather deserialization must be preceded by a strict sidecar/content check.
# benchmark_hpc_utils.R replaces this with the record-aware implementation
# when an HPC worker sources it; the fallback keeps notebook loading safe.
read_feather_checked <- function(path) {
  if (exists("artifact_record_for_load", mode = "function")) {
    artifact_record_for_load(path)
  } else {
    sidecar <- paste0(path, ".md5")
    if (!file.exists(path) || !isTRUE(file.info(path)$size > 0) ||
        !file.exists(sidecar)) {
      stop("Feather checksum validation failed: ", path)
    }
    lines <- readLines(sidecar, warn = FALSE)
    if (length(lines) != 3L ||
        !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
      stop("Feather checksum sidecar has the wrong schema: ", sidecar)
    }
    md5 <- sub("^MD5=", "", lines[[1L]])
    size <- sub("^SIZE=", "", lines[[2L]])
    recorded <- sub("^PATH=", "", lines[[3L]])
    actual <- tolower(unname(tools::md5sum(path)))
    if (!identical(recorded, path) ||
        !grepl("^[0-9a-f]{32}$", md5) ||
        !identical(size, as.character(file.info(path)$size)) ||
        length(actual) != 1L || is.na(actual) ||
        !identical(tolower(md5), actual)) {
      stop("Feather checksum validation failed: ", path)
    }
  }
  arrow::read_feather(path)
}

# MrVI processing
process_mrvi_fig <- function(mrvi_dist_file, labels) {
  arrow::set_cpu_count(1)
  feat_mat <- read_feather_checked(mrvi_dist_file) %>%
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
    feat_mat <- read_rds_checked(gloscope_dist_file)
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
  feat_mat <- read_feather_checked(scpoli_emb_file) %>%
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
  feat_mat <- read_feather_checked(pilot_dist_file) %>%
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
  feat_mat <- read_feather_checked(qot_dist_file) %>%
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
  feat_mat <- read_feather_checked(pilotgm_dist_file) %>%
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
