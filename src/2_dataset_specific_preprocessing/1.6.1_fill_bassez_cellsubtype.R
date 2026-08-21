# ==============================================================================
# 1.6.1_fill_bassez_cellsubtype.R — Fill missing cellSubType with cellType
# ==============================================================================
# In the Bassez dataset, original author subclustering (cellSubType) was only
# performed for T cells and Myeloid cells. All other cells (Cancer cells,
# Fibroblasts, B cells, Endothelial, Mast, etc. — ~68% of the dataset) have
# cellSubType == NA.
#
# Following the original preprocessing pipeline, missing cellSubType annotations
# are filled with the corresponding broad cellType annotation to preserve
# distinct non-immune and B-cell populations in high-resolution compositions.
#
# Must run AFTER 1_stage_data.sh and BEFORE the preprocess array
# (src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh).
# ==============================================================================

library(Seurat)

scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
if (scratch_dir == "") {
  # Fallback to local data dir if running locally
  scratch_dir <- file.path(getwd(), "data")
}

input <- file.path(scratch_dir, "Bassez", "data", "BassezA_2021_33958794whole.rds")
if (!file.exists(input)) {
  # Check alternative local / staging layout
  alt_input <- file.path(scratch_dir, "BassezA_2021_33958794", "output", "BassezA_2021_33958794whole.rds")
  if (file.exists(alt_input)) {
    input <- alt_input
  } else {
    stop("Bassez input RDS not found at: ", input, " or ", alt_input)
  }
}

seurat <- readRDS(input)

is_na <- is.na(seurat$cellSubType) | seurat$cellSubType == "NA" | seurat$cellSubType == ""
n_na <- sum(is_na)
message("Filling ", n_na, "/", ncol(seurat), " missing cellSubType entries with cellType in Bassez...")

seurat$cellSubType <- as.character(seurat$cellSubType)
seurat$cellSubType[is_na] <- as.character(seurat$cellType[is_na])
seurat$cellSubType <- as.factor(seurat$cellSubType)

saveRDS(seurat, input)
message("Saved updated Bassez Seurat object with filled cellSubType: ", input)
