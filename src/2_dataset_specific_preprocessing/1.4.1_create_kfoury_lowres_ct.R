# ==============================================================================
# 1.4.1_create_kfoury_lowres_ct.R — Create the Kfoury `cells_lowres` column
# ==============================================================================
# Ported from Preprocess_datasets.Rmd ("Create low res cell types for Kfoury"):
# collapses the author `cells` annotations into the low-res groups declared in
# datasets.json (columns.cell_type_low_res = "cells_lowres"):
#   Tcells, NKcells, Bcells, MoMac, DCcells (other labels are kept as-is).
# In-place idempotent saveRDS (recomputes cells_lowres each run).
#
# Must run AFTER 1_stage_data.sh and BEFORE the preprocess array
# (src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh).
# ==============================================================================

library(Seurat)

scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
if (scratch_dir == "") {
  stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")
}

input <- file.path(scratch_dir, "Kfoury", "data", "Kfoury_2021_34719426.rds")
seurat <- readRDS(input)

cells_lowres <- as.character(seurat$cells)
cells_lowres[
  cells_lowres %in%
    c("CD4+ Naive", "CD8+ Naive", "CTL-1", "CTL-2",
      "Th1/17", "Treg Active", "Treg Resting")
] <- "Tcells"
cells_lowres[cells_lowres %in% c("NK", "NKT")] <- "NKcells"
cells_lowres[
  cells_lowres %in% c("Immature B cells", "Mature B", "Pro-B", "memBcell")
] <- "Bcells"
cells_lowres[
  cells_lowres %in% c("Mono1", "Mono2", "Mono3", "Monocyte prog", "TAM", "TIM")
] <- "MoMac"
cells_lowres[cells_lowres %in% c("mDC", "pDC")] <- "DCcells"
seurat$cells_lowres <- as.factor(cells_lowres)

saveRDS(seurat, input)
message("Recomputed 'cells_lowres' for Kfoury and saved in place: ", input)
