# ============================================================
# PSEUDOBULK FUNCTIONS
# ============================================================

# DESeq2 normalization for pseudobulk data WITH batch correction
DESeq2.normalize <- function(matrix, metadata, n_hvg = 2000, batch_col = NULL) {
  suppressMessages({
    suppressWarnings({
      
      # 1. Setup design formula based on whether batch is provided
      design_formula <- if (!is.null(batch_col)) paste("~", batch_col) else "~ 1"
      
      # 2. Create DESeq object
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = matrix,
        colData = metadata,
        design = stats::formula(design_formula)
      )

      dds <- DESeq2::estimateSizeFactors(dds)

      # 3. Transform counts using vst
      vst_data <- DESeq2::vst(dds, blind = FALSE)
      norm_matrix <- SummarizedExperiment::assay(vst_data)

      # 4. Remove the batch effect from the matrix
      if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
        norm_matrix <- limma::removeBatchEffect(
          x = norm_matrix, 
          batch = metadata[[batch_col]]
        )
      }

      # 5. Get top variable genes
      rv <- MatrixGenerics::rowVars(norm_matrix)
      select <- order(rv, decreasing = TRUE)[seq_len(min(n_hvg, length(rv)))]
      select <- row.names(norm_matrix)[select]

      final_matrix <- norm_matrix[select[select %in% row.names(norm_matrix)], ]
    })
  })

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