library(Seurat)

scratch_dir <- Sys.getenv("HPC_SCRATCH_DIR")
if (scratch_dir == "") {
  stop("CRITICAL: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")
}

input <- file.path(scratch_dir, "Joanito", "data", "JoaI_2022_35773407_Nofilt_whole.rds")
seurat <- readRDS(input)

seurat$seqtec <- ifelse(
  seurat$dataset %in% c("CRC-SG1", "KUL5"),
  "5' seq",
  "3' seq"
)

# In-place save is idempotent (recomputes seqtec each run).
# Must run AFTER 1_stage_data.sh and BEFORE the preprocess array.
saveRDS(seurat, input)
