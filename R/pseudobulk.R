# ============================================================
# PSEUDOBULK FUNCTIONS
# ============================================================

# DESeq2 normalization for pseudobulk data
DESeq2.normalize <- function(matrix, metadata, n_hvg = 2000) {
  suppressMessages({
    suppressWarnings({
      # Normalize pseudobulk data using DESeq2
      matrix <- DESeq2::DESeqDataSetFromMatrix(
        countData = matrix,
        colData = metadata,
        design = stats::formula(paste("~ 1"))
      )

      matrix <- DESeq2::estimateSizeFactors(matrix)

      # transform counts using vst
      matrix <- DESeq2::vst(matrix)
      matrix <- SummarizedExperiment::assay(matrix)

      # get top variable genes
      rv <- MatrixGenerics::rowVars(matrix)
      select <- order(rv, decreasing = TRUE)[seq_len(min(n_hvg, length(rv)))]
      select <- row.names(matrix)[select]

      matrix <- matrix[select[select %in% row.names(matrix)], ]
    })
  })

  return(matrix)
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
  black_list = "none"
) {
  pb <- get_pb(seurat, sample_col = sample_col, hvg = hvg)

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
  pb_norm <- t(DESeq2.normalize(pb, metadata = metadata, n_hvg = n_hvg))
  return(pb_norm)
}