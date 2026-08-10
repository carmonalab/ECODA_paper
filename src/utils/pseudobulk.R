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
  colnames(pb) <- gsub("-", "_", colnames(pb))
  if (!is.null(hvg)) {
    pb <- pb[hvg, ]
  }
  return(pb)
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

  # Get default black list from STACAS
  data("default_black_list")
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
  metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- t(DESeq2.normalize(
    pb,
    metadata = metadata,
    n_hvg = n_hvg,
    batch_col = batch_col,
    blind = blind,
    correct_batch = correct_batch
  ))
  return(pb_norm)
}