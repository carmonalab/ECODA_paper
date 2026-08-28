# ============================================================
# PSEUDOBULK FUNCTIONS
# ============================================================

# DESeq2 normalization for pseudobulk data WITH batch correction
# Defaults preserve legacy behavior for every existing caller (benchmark):
# batch_col=NULL, blind=TRUE, correct_batch=FALSE -> design ~ 1, vst(blind=TRUE).
# Batch-effect analysis: batch_col + blind=FALSE + correct_batch=TRUE (batch-only
# limma::removeBatchEffect(x, batch), NO design argument -> no-leakage by construction).
DESeq2.normalize <- function(
  matrix,
  metadata,
  n_hvg = 2000,
  batch_col = NULL,      # batch column (batch-effect analysis only)
  blind = TRUE,          # explicit; benchmark = TRUE (legacy-equivalent)
  correct_batch = FALSE  # apply limma::removeBatchEffect(x, batch) — batch-only, no design protection
) {
  # 1. Setup design formula based on whether batch is provided
  design_formula <- if (!is.null(batch_col)) paste("~", batch_col) else "~ 1"

  # 2. Create DESeq object
  dds <- suppressMessages(
    suppressWarnings(
      DESeq2::DESeqDataSetFromMatrix(
        countData = matrix,
        colData = metadata,
        design = stats::formula(design_formula)
      )
    )
  )

  dds <- suppressMessages(suppressWarnings(DESeq2::estimateSizeFactors(dds)))

  # 3. Transform counts using batch-aware vst, with sparsity fallback chain.
  # vst() hard-stops when < nsub genes have mean normalized count > 5
  # (unconditional check in DESeq2::vst); tiny per-cell-type pseudobulk subsets
  # (e.g. Kfoury) need nsub capped at the actual count and a log2(counts+1)
  # fallback when nothing passes the threshold.
  n_gt5 <- sum(MatrixGenerics::rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)
  if (n_gt5 == 0) {
    message("DESeq2.normalize: no gene with mean normalized count > 5; using log2(counts+1)")
    norm_matrix <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
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
          ),
          error = function(e3) {
            warning(
              "DESeq2.normalize: vst failed (", conditionMessage(e),
              "); using log2(counts+1)"
            )
            log2(DESeq2::counts(dds, normalized = TRUE) + 1)
          }
        )
      )
    )
  }

  # 4. Remove the batch effect from the matrix (batch-only, no design protection)
  if (correct_batch) {
    if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
      norm_matrix <- limma::removeBatchEffect(
        x = norm_matrix,
        batch = metadata[[batch_col]]
      )
    } else {
      warning(
        "DESeq2.normalize: correct_batch=TRUE but batch_col is NULL or not in ",
        "metadata; skipping batch correction"
      )
    }
  }

  # 5. Get top variable genes
  rv <- MatrixGenerics::rowVars(norm_matrix)
  select <- order(rv, decreasing = TRUE)[seq_len(min(n_hvg, length(rv)))]
  select <- row.names(norm_matrix)[select]

  final_matrix <- norm_matrix[select[select %in% row.names(norm_matrix)], ]

  return(final_matrix)
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