# ============================================================
# PSEUDOBULK FUNCTIONS
# ============================================================

# Validate the raw genes-by-samples count matrix at the R/DESeq2 boundary.
# DESeq2 ultimately consumes an ordinary R integer matrix; validating the
# complete double/int64-like input before coercion prevents NA-producing
# overflow casts from reaching DESeq2.
validate_pseudobulk_counts_matrix <- function(
  counts,
  context = "pseudobulk",
  max_value = .Machine$integer.max
) {
  if (!is.character(context) || length(context) != 1L ||
      is.na(context) || !nzchar(context)) {
    stop("pseudobulk validation context must be one non-empty string")
  }
  if (!is.numeric(max_value) || length(max_value) != 1L ||
      is.na(max_value) || !is.finite(max_value) ||
      max_value < 0 || max_value != floor(max_value) ||
      max_value > .Machine$integer.max) {
    stop(
      context,
      " max_value must be an integer in [0, .Machine$integer.max]"
    )
  }
  if (is.null(dim(counts)) || length(dim(counts)) != 2L) {
    stop(context, " counts must be a numeric two-dimensional matrix")
  }
  # Sparse Matrix objects are valid matrix inputs at this boundary.  DESeq2
  # consumes a dense R integer matrix, so materialize only here after the
  # dimensionality check; the H5AD reducer remains sparse/on-disk upstream.
  if (inherits(counts, "Matrix")) counts <- as.matrix(counts)
  if (!is.numeric(counts)) {
    stop(context, " counts must be a numeric two-dimensional matrix")
  }

  count_dim <- dim(counts)
  if (anyNA(count_dim) || any(count_dim < 1L)) {
    stop(context, " counts must have at least one gene and one sample")
  }
  count_names <- dimnames(counts)
  if (is.null(count_names) || length(count_names) != 2L ||
      is.null(count_names[[1L]]) || is.null(count_names[[2L]])) {
    stop(context, " counts must have gene and sample identifiers")
  }
  gene_ids <- as.character(count_names[[1L]])
  sample_ids <- as.character(count_names[[2L]])
  if (length(gene_ids) != count_dim[[1L]] ||
      length(sample_ids) != count_dim[[2L]]) {
    stop(context, " counts identifiers do not match matrix dimensions")
  }
  if (anyNA(gene_ids) || anyNA(sample_ids) ||
      any(!nzchar(trimws(gene_ids))) || any(!nzchar(trimws(sample_ids))) ||
      anyDuplicated(gene_ids) || anyDuplicated(sample_ids)) {
    stop(context, " counts gene and sample identifiers must be nonmissing, ",
         "non-empty, and unique")
  }

  # Check every value while it is still in its source representation.  In
  # particular, do not coerce first: as.integer() can turn an oversized
  # double into NA and would conceal the actual DESeq2 ceiling violation.
  if (anyNA(counts) || any(!is.finite(counts))) {
    stop(context, " counts must contain only finite values")
  }
  if (any(counts < 0)) {
    stop(context, " counts must be nonnegative")
  }
  if (any(counts != round(counts))) {
    stop(context, " counts must be integer-valued")
  }
  if (any(counts > max_value)) {
    stop(
      context,
      " counts exceed the supported integer ceiling (.Machine$integer.max)"
    )
  }

  counts_integer <- matrix(
    as.integer(counts),
    nrow = count_dim[[1L]],
    ncol = count_dim[[2L]],
    dimnames = list(gene_ids, sample_ids)
  )
  if (anyNA(counts_integer)) {
    stop(context, " counts could not be represented as R integers")
  }
  counts_integer
}


# Normalize the public scalar compatibility argument into an ordered vector of
# original metadata column names. Corrected callers may provide several keys.
.pseudobulk_batch_keys <- function(batch_col, reject_reserved = FALSE) {
  if (is.null(batch_col)) return(character())
  keys <- if (is.character(batch_col)) {
    unname(batch_col)
  } else if (is.list(batch_col) && length(batch_col) > 0L) {
    vapply(seq_along(batch_col), function(index) {
      value <- batch_col[[index]]
      if (!is.character(value) || length(value) != 1L || is.na(value)) {
        stop("pseudobulk batch_col must contain one or more column names")
      }
      unname(value)
    }, character(1))
  } else {
    stop("pseudobulk batch_col must be NULL, a column name, or column names")
  }
  if (length(keys) == 0L || anyNA(keys) ||
      any(!nzchar(trimws(keys)))) {
    stop("pseudobulk batch_col must contain non-empty column names")
  }
  if (anyDuplicated(keys)) {
    stop("pseudobulk batch_col must not contain duplicate column names")
  }
  if (isTRUE(reject_reserved) &&
      any(keys == "__ecoda_batch_combined_v1")) {
    stop(
      "pseudobulk batch_col must use original technical columns; ",
      "__ecoda_batch_combined_v1 is reserved"
    )
  }
  keys
}


# Reorder metadata to the count matrix columns. H5AD aggregation can retain
# first-seen sample order while its canonical metadata bundle uses another
# order; fitting always follows the count columns, and direct publication
# restores the metadata order below.
.pseudobulk_align_metadata <- function(
  metadata,
  sample_ids,
  batch_col = NULL,
  reject_reserved = FALSE
) {


  if (is.null(metadata)) {
    stop("pseudobulk metadata is required")
  }
  if (!is.data.frame(metadata)) {
    metadata <- tryCatch(
      as.data.frame(metadata),
      error = function(e) {
        stop("pseudobulk metadata must be coercible to a data.frame")
      }
    )
  }
  if (nrow(metadata) != length(sample_ids)) {
    stop("pseudobulk metadata row count does not match sample count")
  }

  if ("Sample" %in% colnames(metadata)) {
    metadata_ids <- as.character(metadata[["Sample"]])
  } else {
    metadata_ids <- rownames(metadata)
  }
  if (is.null(metadata_ids) || length(metadata_ids) != nrow(metadata)) {
    stop("pseudobulk metadata must contain canonical Sample identifiers")
  }
  metadata_ids <- as.character(metadata_ids)
  if (anyNA(metadata_ids) || any(!nzchar(trimws(metadata_ids))) ||
      anyDuplicated(metadata_ids)) {
    stop("pseudobulk metadata Sample identifiers must be nonmissing, ",
         "non-empty, and unique")
  }

  sample_ids <- as.character(sample_ids)
  metadata_index <- match(sample_ids, metadata_ids)
  if (anyNA(metadata_index) && exists(
    "standardize_sample_names",
    mode = "function",
    inherits = TRUE
  )) {
    sample_alias <- standardize_sample_names(sample_ids)
    metadata_alias <- standardize_sample_names(metadata_ids)
    if (!anyDuplicated(sample_alias) && !anyDuplicated(metadata_alias)) {
      metadata_index <- match(sample_alias, metadata_alias)
    }
  }
  if (anyNA(metadata_index)) {
    stop("pseudobulk metadata Sample identifiers do not match count columns")
  }

  metadata <- metadata[metadata_index, , drop = FALSE]
  rownames(metadata) <- sample_ids
  batch_keys <- .pseudobulk_batch_keys(
    batch_col,
    reject_reserved = reject_reserved
  )
  if (length(batch_keys) > 0L &&
      any(!batch_keys %in% colnames(metadata))) {
    stop(
      "pseudobulk batch_col is missing from metadata: ",
      paste(batch_keys[!batch_keys %in% colnames(metadata)], collapse = ", ")
    )
  }
  if (isTRUE(reject_reserved) &&
      "__ecoda_batch_combined_v1" %in% colnames(metadata)) {
    stop(
      "pseudobulk metadata contains the reserved combined batch column: ",
      "__ecoda_batch_combined_v1"
    )
  }
  metadata
}


# Validate and factorize original technical columns for corrected fitting.
# Numeric-looking and suspension-like categorical values are deliberately kept
# as categorical levels; corrected DESeq2 fitting itself remains intercept-only.
.pseudobulk_factor_column <- function(values, key) {
  if (is.list(values) && !is.factor(values)) {
    stop("pseudobulk batch column is list-valued: ", key)
  }
  if (length(values) == 0L) {
    stop("pseudobulk batch column is empty: ", key)
  }
  if (anyNA(values)) {
    stop("pseudobulk batch column has missing values: ", key)
  }
  if (is.numeric(values) && any(!is.finite(values))) {
    stop("pseudobulk batch column has non-finite values: ", key)
  }
  values <- as.character(values)
  trimmed <- trimws(values)
  lower <- tolower(trimmed)
  invalid_sentinels <- c(
    "na", "nan", "none", "null", "<na>", "n/a", "unknown",
    "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"
  )
  # Breast's configured suspension-duration field permits only the literal
  # ``unknown`` token as an ordinary categorical level. Case variants and
  # surrounding whitespace remain invalid sentinels, as do unknown values on
  # every other technical key.
  invalid_mask <- lower %in% invalid_sentinels
  if (identical(key, "suspension_dissociation_time")) {
    invalid_mask <- invalid_mask & values != "unknown"
  }
  if (anyNA(values) || any(!nzchar(trimmed)) || any(invalid_mask)) {
    stop("pseudobulk batch column has missing or blank values: ", key)
  }
  factor(values, levels = unique(values))
}


# Build the additive fixed-effect design used only by corrected pseudobulk.
# Constants remain valid metadata but are omitted from the technical
# covariates. A rank-deficient or saturated design is rejected before limma.
.pseudobulk_corrected_design <- function(metadata, batch_keys, sample_ids) {
  if (length(batch_keys) == 0L) {
    stop("corrected pseudobulk fitting requires one or more technical keys")
  }
  if (nrow(metadata) != length(sample_ids) || length(sample_ids) < 2L) {
    stop("corrected pseudobulk fitting requires at least two samples")
  }

  factor_metadata <- metadata
  key_levels <- vector("list", length(batch_keys))
  names(key_levels) <- batch_keys
  for (index in seq_along(batch_keys)) {
    key <- batch_keys[[index]]
    factor_metadata[[key]] <- .pseudobulk_factor_column(
      factor_metadata[[key]],
      key
    )
    key_levels[[index]] <- levels(factor_metadata[[key]])
  }

  effective_keys <- unname(batch_keys[
    vapply(key_levels, length, integer(1)) >= 2L
  ])
  non_estimable_keys <- unname(setdiff(batch_keys, effective_keys))
  model_data <- data.frame(
    row.names = sample_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  aliases <- setNames(
    paste0("batch_key_", match(effective_keys, batch_keys)),
    effective_keys
  )
  for (key in effective_keys) {
    model_data[[unname(aliases[[key]])]] <- factor_metadata[[key]]
  }

  if (length(effective_keys) == 0L) {
    design <- matrix(
      1,
      nrow = length(sample_ids),
      ncol = 1L,
      dimnames = list(sample_ids, "(Intercept)")
    )
  } else {
    design_formula <- stats::as.formula(paste0(
      "~ 1 + ", paste(unname(aliases), collapse = " + ")
    ))
    design <- tryCatch(
      stats::model.matrix(design_formula, data = model_data),
      error = function(error) {
        stop(
          "corrected pseudobulk batch design could not be constructed: ",
          conditionMessage(error)
        )
      }
    )
    rownames(design) <- sample_ids
  }
  if (is.null(dim(design)) || nrow(design) != length(sample_ids) ||
      any(!is.finite(design))) {
    stop("corrected pseudobulk batch design is invalid")
  }
  rank <- qr(design, tol = 1e-10)$rank
  if (rank < ncol(design)) {
    stop("corrected pseudobulk batch design is rank deficient")
  }
  residual_df <- nrow(design) - rank
  if (residual_df <= 0L) {
    stop(
      "corrected pseudobulk batch design has nonpositive residual degrees ",
      "of freedom"
    )
  }
  technical_columns <- which(colnames(design) != "(Intercept)")
  technical_covariates <- design[
    , technical_columns,
    drop = FALSE
  ]
  preserve_design <- matrix(
    1,
    nrow = nrow(design),
    ncol = 1L,
    dimnames = list(sample_ids, "(Intercept)")
  )
  list(
    metadata = factor_metadata,
    configured_batch_keys = unname(batch_keys),
    effective_batch_keys = effective_keys,
    non_estimable_batch_keys = non_estimable_keys,
    key_levels = key_levels,
    design = design,
    technical_covariates = technical_covariates,
    preserve_design = preserve_design,
    aliases = aliases,
    rank = as.integer(rank),
    columns = as.integer(ncol(design)),
    residual_df = as.integer(residual_df),
    no_op = length(effective_keys) == 0L
  )
}


# Fit and transform every gene once.  HVG selection is deliberately kept out
# of this function so hvg500/1000/2000/3000 can share one DESeq2/VST pass.
# This private path requires counts to have already passed
# validate_pseudobulk_counts_matrix().
.fit_pseudobulk_deseq2_validated <- function(
  counts,
  metadata,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
) {
  if (!is.logical(blind) || length(blind) != 1L || is.na(blind) ||
      !is.logical(correct_batch) || length(correct_batch) != 1L ||
      is.na(correct_batch)) {
    stop("pseudobulk blind and correct_batch must be scalar logical values")
  }
  batch_keys <- .pseudobulk_batch_keys(
    batch_col,
    reject_reserved = isTRUE(correct_batch)
  )
  metadata <- .pseudobulk_align_metadata(
    metadata,
    colnames(counts),
    batch_col = if (length(batch_keys)) batch_keys else NULL,
    reject_reserved = isTRUE(correct_batch)
  )

  corrected_design <- NULL
  if (isTRUE(correct_batch)) {
    corrected_design <- .pseudobulk_corrected_design(
      metadata = metadata,
      batch_keys = batch_keys,
      sample_ids = colnames(counts)
    )
    metadata <- corrected_design[["metadata"]]
  }
  # Corrected DESeq2 normalization sees only Sample and the requested original
  # technical keys. Biological labels and every other metadata field remain
  # outside the normalization boundary.
  metadata_for_fit <- metadata
  if (isTRUE(correct_batch)) {
    fit_columns <- unique(c(
      if ("Sample" %in% colnames(metadata)) "Sample" else character(),
      batch_keys
    ))
    metadata_for_fit <- metadata[, fit_columns, drop = FALSE]
    rownames(metadata_for_fit) <- colnames(counts)
  }

  # Corrected mode removes separate technical fixed effects after a
  # blind=FALSE transformation, so DESeq2 itself always fits ~1. Ordinary and
  # uncorrected callers retain their historical batch-aware design behavior.
  design_formula <- if (isTRUE(correct_batch) || length(batch_keys) == 0L) {
    stats::formula("~ 1")
  } else {
    stats::reformulate(batch_keys)
  }
  dds <- suppressMessages(
    suppressWarnings(
      DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata_for_fit,
        design = design_formula
      )
    )
  )
  dds <- suppressMessages(suppressWarnings(DESeq2::estimateSizeFactors(dds)))

  # vst() hard-stops when fewer than nsub genes have mean normalized count > 5
  # (an unconditional check in DESeq2::vst). Keep the historical fallback
  # chain for sparse/tiny cell-type pseudobulks.
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  n_gt5 <- sum(MatrixGenerics::rowMeans(normalized_counts) > 5)
  if (n_gt5 == 0) {
    message(
      "DESeq2.normalize: no gene with mean normalized count > 5; ",
      "using log2(counts+1)"
    )
    norm_matrix <- log2(normalized_counts + 1)
  } else {
    norm_matrix <- tryCatch(
      suppressMessages(
        suppressWarnings(
          SummarizedExperiment::assay(
            DESeq2::vst(dds, blind = blind, nsub = min(1000, n_gt5))
          )
        )
      ),
      error = function(e) tryCatch(
        suppressMessages(
          suppressWarnings(
            SummarizedExperiment::assay(
              DESeq2::varianceStabilizingTransformation(
                dds, blind = blind, fitType = "mean"
              )
            )
          )
        ),
        error = function(e2) tryCatch(
          suppressMessages(
            suppressWarnings(
              SummarizedExperiment::assay(
                DESeq2::varianceStabilizingTransformation(
                  dds, blind = TRUE, fitType = "mean"
                )
              )
            )
          )
        ),
        error = function(e3) {
          warning(
            "DESeq2.normalize: vst failed (", conditionMessage(e),
            "); using log2(counts+1)"
          )
          log2(normalized_counts + 1)
        }
      )
    )
  }
  if (is.null(dimnames(norm_matrix))) {
    dimnames(norm_matrix) <- dimnames(counts)
  } else {
    if (is.null(rownames(norm_matrix))) rownames(norm_matrix) <- rownames(counts)
    if (is.null(colnames(norm_matrix))) colnames(norm_matrix) <- colnames(counts)
  }

  # Batch correction is technical-only. Build one dummy-variable model from
  # every varying original factor and protect only an intercept in limma.
  if (isTRUE(correct_batch) && !isTRUE(corrected_design[["no_op"]])) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("limma is required for corrected pseudobulk batch correction")
    }
    corrected <- tryCatch(
      limma::removeBatchEffect(
        x = norm_matrix,
        covariates = corrected_design[["technical_covariates"]],
        design = corrected_design[["preserve_design"]]
      ),
      error = function(error) {
        stop(
          "DESeq2.normalize: limma batch correction failed: ",
          conditionMessage(error)
        )
      }
    )
    if (is.null(dim(corrected)) ||
        !identical(dim(corrected), dim(norm_matrix)) ||
        any(!is.finite(corrected))) {
      stop("DESeq2.normalize: limma batch correction produced invalid output")
    }
    dimnames(corrected) <- dimnames(norm_matrix)
    norm_matrix <- corrected
  }

  row_variances <- MatrixGenerics::rowVars(norm_matrix)
  variance_order <- rownames(norm_matrix)[
    order(row_variances, decreasing = TRUE)
  ]
  result <- list(
    norm_matrix = norm_matrix,
    normalized_matrix = norm_matrix,
    variance_order = variance_order,
    variance_ordering = variance_order,
    row_variances = row_variances,
    counts = counts,
    metadata = metadata_for_fit,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  )
  if (isTRUE(correct_batch)) {
    result[["batch_correction"]] <- list(
      correction_mode = "limma_fixed_effects_pseudobulk",
      correction_formula = if (isTRUE(corrected_design[["no_op"]])) {
        "NO_CORRECTION: no estimable technical batch key"
      } else {
        paste0(
          "DESeq2 design=~ 1; model.matrix(~ 1 + ",
          paste(unname(corrected_design[["aliases"]]), collapse = " + "),
          "); limma::removeBatchEffect(",
          "covariates=technical_covariates, design=intercept)"
        )
      },
      correction_state = if (isTRUE(corrected_design[["no_op"]])) {
        "NO_CORRECTION"
      } else {
        "BATCH_CORRECTION"
      },
      configured_batch_keys = corrected_design[["configured_batch_keys"]],
      effective_batch_keys = corrected_design[["effective_batch_keys"]],
      non_estimable_batch_keys = corrected_design[["non_estimable_batch_keys"]],
      aliases = corrected_design[["aliases"]],
      design_rank = corrected_design[["rank"]],
      design_columns = corrected_design[["columns"]],
      design_residual_df = corrected_design[["residual_df"]]
    )
  }
  result
}

# Public fitting boundary.  Independent callers must always receive the full
# count-matrix validation/coercion before entering the private fit path.
fit_pseudobulk_deseq2 <- function(
  counts,
  metadata,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
) {
  counts <- validate_pseudobulk_counts_matrix(
    counts,
    context = "pseudobulk",
    max_value = .Machine$integer.max
  )
  .fit_pseudobulk_deseq2_validated(
    counts,
    metadata = metadata,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  )
}


# Select an HVG prefix from a shared full-gene fit.  The historical
# hvg2000_bl path passed the literal "default_without_sex_genes" into a
# `%in%` expression and therefore intentionally did not remove genes; retain
# that no-op here.
select_pseudobulk_deseq2 <- function(
  fit,
  n_hvg,
  black_list = "none"
) {
  if (!is.list(fit) || is.null(fit[["norm_matrix"]]) ||
      is.null(fit[["variance_order"]])) {
    stop("pseudobulk DESeq2 fit is missing normalized matrix/order")
  }
  norm_matrix <- fit[["norm_matrix"]]
  if (is.null(dim(norm_matrix)) || length(dim(norm_matrix)) != 2L ||
      is.null(rownames(norm_matrix)) || is.null(colnames(norm_matrix))) {
    stop("pseudobulk DESeq2 fit normalized matrix has invalid dimensions")
  }
  if (is.null(n_hvg)) {
    n_hvg <- nrow(norm_matrix)
  }
  if (!is.numeric(n_hvg) || length(n_hvg) != 1L || is.na(n_hvg) ||
      !is.finite(n_hvg) || n_hvg < 0 || n_hvg != floor(n_hvg)) {
    stop("pseudobulk n_hvg must be a nonnegative integer or NULL")
  }
  genes <- as.character(fit[["variance_order"]])
  if (anyNA(genes) || anyDuplicated(genes) ||
      any(!genes %in% rownames(norm_matrix))) {
    stop("pseudobulk DESeq2 fit variance order is invalid")
  }
  if (length(black_list) == 1L &&
      is.character(black_list) &&
      black_list %in% c("none", "default_without_sex_genes")) {
    # Intentional no-op for both the unfiltered and legacy hvg2000_bl flags.
  } else if (!is.null(black_list)) {
    if (!is.character(black_list)) {
      stop("pseudobulk black_list must be a character vector")
    }
    genes <- genes[!genes %in% black_list]
  }
  n_keep <- min(n_hvg, length(genes))
  genes <- genes[seq_len(n_keep)]
  norm_matrix[genes, , drop = FALSE]
}


# Direct matrix/metadata boundary for the canonical H5AD path.  The input is
# genes-by-samples; publication is samples-by-genes in canonical metadata
# order.  No Seurat object or AggregateExpression() call is involved.
get_pb_deseq2_from_counts <- function(
  counts,
  metadata,
  hvg = NULL,
  n_hvg = 2000,
  black_list = "none",
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
) {
  counts <- validate_pseudobulk_counts_matrix(
    counts,
    context = "pseudobulk",
    max_value = .Machine$integer.max
  )
  if (!is.data.frame(metadata)) {
    metadata <- tryCatch(
      as.data.frame(metadata),
      error = function(e) {
        stop("pseudobulk metadata must be coercible to a data.frame")
      }
    )
  }
  if ("Sample" %in% colnames(metadata)) {
    canonical_sample_ids <- as.character(metadata[["Sample"]])
  } else {
    canonical_sample_ids <- rownames(metadata)
  }
  if (is.null(canonical_sample_ids) ||
      length(canonical_sample_ids) != nrow(metadata) ||
      anyNA(canonical_sample_ids) ||
      any(!nzchar(trimws(canonical_sample_ids))) ||
      anyDuplicated(canonical_sample_ids)) {
    stop("pseudobulk metadata Sample identifiers must be nonmissing, ",
         "non-empty, and unique")
  }
  canonical_sample_ids <- as.character(canonical_sample_ids)
  if (nrow(metadata) != ncol(counts)) {
    stop("pseudobulk metadata row count does not match sample count")
  }

  fit_counts <- counts
  if (!is.null(hvg)) {
    if (!is.character(hvg) || length(hvg) < 1L || anyNA(hvg) || any(!nzchar(hvg)) ||
        anyDuplicated(hvg) || any(!hvg %in% rownames(counts))) {
      stop("pseudobulk hvg identifiers must be present, nonmissing, ",
           "non-empty, and unique")
    }
    # schvg2000 intentionally fits its prefiltered raw gene universe
    # separately rather than sharing the full-gene fit.
    fit_counts <- counts[hvg, , drop = FALSE]
  }
  fit <- .fit_pseudobulk_deseq2_validated(
    fit_counts,
    metadata = metadata,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  )
  selected <- select_pseudobulk_deseq2(
    fit,
    n_hvg = n_hvg,
    # Explicit hvg input is already a prefilter; the legacy wrapper bypassed
    # blacklist handling in that case as well.
    black_list = if (is.null(hvg)) black_list else "none"
  )
  published <- t(selected)

  # Restore canonical metadata order after fitting in count-column order.
  count_sample_ids <- colnames(counts)
  output_index <- match(canonical_sample_ids, count_sample_ids)
  if (anyNA(output_index) && exists(
    "standardize_sample_names",
    mode = "function",
    inherits = TRUE
  )) {
    count_alias <- standardize_sample_names(count_sample_ids)
    canonical_alias <- standardize_sample_names(canonical_sample_ids)
    if (!anyDuplicated(count_alias) && !anyDuplicated(canonical_alias)) {
      output_index <- match(canonical_alias, count_alias)
    }
  }
  if (anyNA(output_index)) {
    stop("pseudobulk metadata Sample identifiers do not match count columns")
  }
  rownames(published) <- count_sample_ids
  published <- published[output_index, , drop = FALSE]
  rownames(published) <- canonical_sample_ids
  published
}


# Compatibility core used by maintained legacy Seurat callers.  It retains the
# historical genes-by-samples return orientation; direct callers should use
# get_pb_deseq2_from_counts(), which publishes samples-by-genes.
DESeq2.normalize <- function(
  matrix,
  metadata,
  n_hvg = 2000,
  batch_col = NULL,
  blind = TRUE,
  correct_batch = FALSE
) {
  fit <- fit_pseudobulk_deseq2(
    matrix,
    metadata = metadata,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  )
  select_pseudobulk_deseq2(fit, n_hvg = n_hvg)
}


# Get pseudobulk matrix from seurat object
get_pb <- function(seurat, sample_col = "Sample", hvg = NULL) {
  pb <- as.matrix(AggregateExpression(
    seurat,
    group.by = sample_col,
    assays = "RNA"
  )[["RNA"]])
  # Seurat may sanitize separators while aggregating (e.g.
  # "BIOKEY-2-Pre" vs "BIOKEY_2_Pre"). get_pb_deseq2() reconciles the
  # resulting sample IDs against the canonical preprocessed obs IDs before
  # returning a samples-by-genes matrix.
  if (!is.null(hvg)) {
    pb <- pb[hvg, ]
  }
  return(pb)
}
# Reconcile Seurat/AggregateExpression sample IDs with the canonical IDs from
# the preprocessed obs. Exact matches win; alias matching is accepted only when
# it is one-to-one and covers every sample.
align_pseudobulk_sample_names <- function(pb, sample_ids) {
  pb_ids <- rownames(pb)
  sample_ids <- as.character(sample_ids)
  if (is.null(pb_ids) ||
      length(pb_ids) != length(sample_ids) ||
      anyNA(pb_ids) ||
      anyNA(sample_ids) ||
      any(!nzchar(pb_ids)) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(pb_ids) ||
      anyDuplicated(sample_ids)) {
    stop("Pseudobulk sample IDs must be nonmissing and unique.")
  }

  exact_match <- match(pb_ids, sample_ids)
  if (all(!is.na(exact_match))) {
    rownames(pb) <- sample_ids[exact_match]
    # AggregateExpression commonly sorts its sample columns. Reindex the
    # normalized samples to the first-appearance order from the canonical
    # preprocessed obs before any downstream result bundle is written.
    return(pb[sample_ids, , drop = FALSE])
  }

  if (!exists("standardize_sample_names", mode = "function")) {
    stop("Cannot reconcile pseudobulk sample IDs: standardize_sample_names is unavailable.")
  }
  pb_alias <- standardize_sample_names(pb_ids)
  sample_alias <- standardize_sample_names(sample_ids)
  if (anyDuplicated(pb_alias) || anyDuplicated(sample_alias)) {
    stop("Pseudobulk sample-ID aliasing is ambiguous.")
  }
  alias_match <- match(pb_alias, sample_alias)
  if (anyNA(alias_match)) {
    stop(
      "Pseudobulk sample IDs do not match canonical metadata IDs after ",
      "standardization."
    )
  }
  rownames(pb) <- sample_ids[alias_match]
  # Alias reconciliation changes names but must not preserve AggregateExpression
  # order; all consumers use the canonical obs order.
  pb[sample_ids, , drop = FALSE]
}



# Get DESeq2-normalized pseudobulk
get_pb_deseq2 <- function(
  seurat,
  sample_col = "Sample",
  hvg = NULL,
  n_hvg = 2000,
  black_list = "none",
  batch_col = NULL,      # batch column (batch-effect analysis only)
  blind = TRUE,          # benchmark = TRUE (legacy-equivalent)
  correct_batch = FALSE  # batch-only limma::removeBatchEffect, no design protection
) {
  pb <- get_pb(seurat, sample_col = sample_col, hvg = hvg)

  # Get default black list from STACAS. Plain `data("default_black_list")`
  # only searches ATTACHED packages; on worker nodes STACAS is installed
  # (pixi.toml, install_github) but never library()-ed, so the data set was
  # "not found" -> `object 'black.list' not found` -> prepare_pseudobulk task
  # crashed (observed 2026-08-17: Adams, the first dataset to recompute after
  # the Aug-12 cache; all other datasets had cached variants and skipped).
  # Load explicitly from the package when available; fall back to an empty
  # list otherwise — safe because the hvg2000_bl filter is a documented no-op
  # anyway (`pb[!rownames(pb) %in% black_list, ]` tests the literal flag
  # string, not this object).
  black.list <- NULL
  if (requireNamespace("STACAS", quietly = TRUE)) {
    suppressWarnings(
      data("default_black_list", package = "STACAS", envir = environment())
    )
  } else {
    message("get_pb_deseq2: STACAS unavailable; default_black_list fallback = empty")
  }
  if (is.null(black.list)) black.list <- character(0)
  default_black_list <- black.list

  if (is.null(hvg) & black_list == "default") {
    default_black_list <- unlist(default_black_list)
    pb <- pb[!rownames(pb) %in% default_black_list, ]
  } else if (is.null(hvg) & black_list == "default_without_sex_genes") {
    default_black_list <- default_black_list[
      !names(default_black_list) %in% c("Xgenes", "Ygenes")
    ]
    default_black_list <- unlist(default_black_list)
    pb <- pb[!rownames(pb) %in% black_list, ]
  }

  metadata <- get_metadata(seurat, sample_col = sample_col)
  pb_norm <- t(DESeq2.normalize(
    pb,
    metadata = metadata,
    n_hvg = n_hvg,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  ))
  pb_norm <- align_pseudobulk_sample_names(
    pb_norm,
    metadata[[sample_col]]
  )
  return(pb_norm)
}