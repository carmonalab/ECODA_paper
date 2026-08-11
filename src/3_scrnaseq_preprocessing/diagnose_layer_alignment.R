#!/usr/bin/env Rscript
# ============================================================
# DIAGNOSTIC: anndataR layer-alignment failures in preprocessing
# ============================================================
# Background: some source RDSes (e.g. Wu) carry Assay5 layers whose gene names
# are NOT identical to the assay features (mismatch in the middle of the gene
# list; head/tail identical), which makes anndataR's write_h5ad fail with:
#   Error in private$.validate_aligned_mapping(value, "layers", ...)
#   Expected column names: ... Provided column names: ...
#
# This script prints the layer/dimnames layout of an RDS (raw and after
# create_clean_seuratv5_object) so the exact quirk can be confirmed before
# relying on the generalized alignment in preprocess_utils.py
# (convert_rds_to_raw_h5ad).
#
# Usage (HPC, working dir = PROJECT_ROOT; use the pixi R, never bare Rscript).
# NOTE: slurm_config.sh vars (HPC_SCRATCH_DIR) are NOT set in interactive login
# shells — use the concrete scratch path below or `source src/slurm_config.sh`
# first:
#   pixi run -e py-cuda13 Rscript src/3_scrnaseq_preprocessing/diagnose_layer_alignment.R \
#     ${HOME}/scratch/ECODA_paper/Wu/data/WuS_2021_34493872.rds
#
# For large RDS files run it via srun on the debug-cpu partition instead of the
# login node:
#   srun --partition=debug-cpu --mem=16G --cpus-per-task=2 --time=00:30:00 \
#     pixi run -e py-cuda13 Rscript src/3_scrnaseq_preprocessing/diagnose_layer_alignment.R \
#     ${HOME}/scratch/ECODA_paper/Wu/data/WuS_2021_34493872.rds
#
# Expected output interpretation:
#   - "rownames vs features: identical=FALSE setequal=TRUE"  -> middle-of-list
#     misalignment; the generalized alignment reindexes it (message: "reindexing").
#   - "colnames vs features: setequal=TRUE" with NULL rownames -> cell-major
#     (transposed) layer; the generalized alignment transposes it back
#     (message: "transposing cell-major").
#   - "names attr: ... NON-NULL!" on rownames/colnames or on the assay features
#     dimnames -> the Wu quirk: an attribute-only difference that identical()
#     (anndataR layer validation) still fails on. Element-wise checks
#     (setequal / identical on unname()d dimnames) cannot catch it — this
#     script now reports it explicitly.
#   - "setequal=FALSE" with missing/extra genes -> genuinely different gene
#     set; alignment fails closed (stop) and the RDS needs manual inspection.
# ============================================================
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
source("src/utils/seurat_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript diagnose_layer_alignment.R <input.rds>")
}
input_path <- args[1]
if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path)
}

describe_layer <- function(m, features, assay_name, lyr, label) {
  cat(sprintf("  [%s] assay '%s' layer '%s': class=%s dims=%dx%d\n",
      label, assay_name, lyr, class(m)[1], nrow(m), ncol(m)))
  if (is.null(rownames(m))) {
    cat("    rownames: NULL\n")
  } else {
    cat(sprintf("    rownames: n=%d head=%s\n",
        length(rownames(m)),
        paste(head(rownames(m), 3), collapse = ", ")))
    cat(sprintf("    rownames vs features: identical=%s setequal=%s missing=%d extra=%d\n",
        identical(rownames(m), features),
        setequal(rownames(m), features),
        sum(!features %in% rownames(m)),
        sum(!rownames(m) %in% features)))
    # Wu quirk: a 'names' attribute on the dimnames makes identical() (anndataR
    # layer validation) fail although all elements are equal; identical() on
    # unname()d dimnames is the reliable check.
    cat(sprintf("    names attr: rownames=%s colnames=%s\n",
        if (is.null(names(rownames(m)))) "NULL" else "NON-NULL!",
        if (is.null(names(colnames(m)))) "NULL" else "NON-NULL!"))
    cat(sprintf("    unname(rownames) vs features: identical=%s\n",
        identical(unname(rownames(m)), features)))
  }
  if (is.null(colnames(m))) {
    cat("    colnames: NULL\n")
  } else {
    cat(sprintf("    colnames: n=%d head=%s\n",
        length(colnames(m)),
        paste(head(colnames(m), 3), collapse = ", ")))
    cat(sprintf("    colnames vs features: setequal=%s (TRUE => genes in colnames, i.e. cell-major layer?)\n",
        setequal(colnames(m), features)))
  }
}

inspect <- function(seurat, label) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("R %s | Seurat %s | Matrix %s\n",
      R.version.string,
      as.character(packageVersion("Seurat")),
      as.character(packageVersion("Matrix"))))
  cat(sprintf("n cells: %d | assays: %s\n",
      ncol(seurat),
      paste(names(seurat@assays), collapse = ", ")))
  for (assay_name in names(seurat@assays)) {
    a <- seurat[[assay_name]]
    cat(sprintf("assay '%s': class=%s n_features=%d\n",
        assay_name, class(a)[1], nrow(a)))
    if (inherits(a, "Assay5")) {
      # Wu quirk lives here: a 'names' attribute on the features dimnames
      # (attr(slot(a, "features"), "dimnames")[[1]]) propagates into layer
      # rownames but is invisible to element-wise comparisons.
      feat_dn <- attr(slot(a, "features"), "dimnames")[[1]]
      cat(sprintf("    features dimnames names attr: rownames(a)=%s features-slot=%s\n",
          if (is.null(names(rownames(a)))) "NULL" else "NON-NULL!",
          if (is.null(names(feat_dn))) "NULL" else "NON-NULL!"))
      for (lyr in names(a@layers)) {
        describe_layer(a@layers[[lyr]], rownames(a), assay_name, lyr, label)
      }
    } else {
      cat("  (not Assay5; no @layers slot)\n")
    }
  }
}

cat("Reading:", input_path, "\n")
seurat_raw <- readRDS(input_path)
inspect(seurat_raw, "RAW RDS")

cat("\n--- create_clean_seuratv5_object() ---\n")
seurat_clean <- create_clean_seuratv5_object(seurat_raw)
inspect(seurat_clean, "AFTER create_clean_seuratv5_object")
