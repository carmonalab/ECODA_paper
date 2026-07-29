#!/bin/bash
#
# 2.1_process_chunk.sh — Process a specific subset of samples sequentially
#

set -euo pipefail

TASK_ID="$1"
export PROJECT_ROOT="$2"
CHUNK_FILE="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"

log_msg() {
  local msg="$1"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[${timestamp}] [TASK ${TASK_ID}] ${msg}"
}

log_msg "============================================"
log_msg "Processing Chunk ID: ${TASK_ID}"
log_msg "============================================"

# --- MODULE PURGE REMOVED ---
log_msg "Activating local Pixi context on worker node..."

if [[ ! -f "${CHUNK_FILE}" ]]; then
  log_msg "ERROR: Chunk file ${CHUNK_FILE} not found."
  exit 1
fi

cd "${PROJECT_ROOT}"

export RENV_CONFIG_EXTERNAL_LIBRARIES="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"
export R_LIBS_SITE="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library:${R_LIBS_SITE:-}"

# Run the inline R script via Pixi instead of default Rscript
"${HOME}/.pixi/bin/pixi" run Rscript --vanilla - <<'RS' "chunk_file__${CHUNK_FILE}"
# ==============================================================================
# PHASE 1: BOOTSTRAP RENV PROFILE (MUST HAPPEN BEFORE ANYTHING ELSE)
# ==============================================================================
project_root <- Sys.getenv("PROJECT_ROOT")

# DISABLE RENV SANDBOX: Force renv to allow fallback access to Pixi packages
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "false")

# Dynamically construct package execution environments based on cluster R build
r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
local_renv_lib <- file.path(project_root, "renv", "library", paste0("R-", r_mm), R.version$platform)
pixi_lib <- file.path(project_root, ".pixi", "envs", "default", "lib", "R", "library")

.libPaths(unique(c(local_renv_lib, pixi_lib, .libPaths())))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("CRITICAL: JSON tracking capabilities missing. Check path maps: ", paste(.libPaths(), collapse = ", "))
}

# ==============================================================================
# PHASE 2: ENVIRONMENT READY -> SAFE TO LOAD CONFIG & LIBRARIES
# ==============================================================================
source(file.path(project_root, "config_helper.R"))
paths <- get_pipeline_config()

library(scGate)
library(ProjecTILs)
library(SignatuR)
library(HiTME)
library(Seurat)
library(arrow)
library(scATOMIC)
library(cutoff.scATOMIC)
library(R.utils)

# Import anndata WITHOUT automatic R conversion
library(reticulate)

# # Explicitly bind reticulate to the Pixi environment's Python
# pixi_python <- file.path(project_root, ".pixi", "envs", "default", "bin", "python")
# reticulate::use_python(pixi_python, required = TRUE)

ad <- import("anndata", convert = FALSE)
set.seed(123)

# ==============================================================================
# HELPER FUNCTIONS & PARAMETER PARSING
# ==============================================================================
get_sample_seurat_obj <- function(adata, r_obs, target_sample, sample_colname) {
  # Get indices for target sample
  sample_indices <- which(r_obs[[sample_colname]] == target_sample) - 1 # -1 for 0-based Python indexing

  # Extract the matrix:
  # slice, load to memory, cast to float (required in case counts were saved as integers), and convert to CSC format
  subset_py <- adata[as.integer(sample_indices)]
  raw_X_py <- subset_py$X$astype("float64")$tocsc()

  # Bring it into R
  # py_to_r will create a standard dgCMatrix
  counts_matrix <- py_to_r(raw_X_py)

  # Create Seurat object
  # Transpose because Python is (Cells x Genes) and Seurat is (Genes x Cells)
  counts_matrix <- t(counts_matrix)
  rownames(counts_matrix) <- as.character(py_to_r(subset_py$var_names$values))
  colnames(counts_matrix) <- as.character(py_to_r(subset_py$obs_names$values))

  seurat_obj <- CreateSeuratObject(
    counts = counts_matrix,
    meta.data = as.data.frame(py_to_r(subset_py$obs))
  )
  return(seurat_obj)
}

# Parse custom environment entries (falling back to original configurations if empty)
env_sample_col    <- Sys.getenv("SAMPLE_COLNAME")
env_tissue        <- Sys.getenv("TISSUE_TYPE")
env_auth_annots   <- Sys.getenv("AUTHOR_ANNOT_COLNAMES")
env_normal_tissue <- Sys.getenv("NORMAL_TISSUE")

defaults <- list(
  chunk_file            = NULL,
  sample_colname        = if (env_sample_col != "") env_sample_col else "sample",
  tissue_type           = if (env_tissue != "") env_tissue else "Tumor",
  author_annot_colnames = if (env_auth_annots != "") unlist(strsplit(env_auth_annots, ",")) else character(),
  normal_tissue         = if (env_normal_tissue != "") as.logical(env_normal_tissue) else TRUE
)

# Column whitelists for annotation output
hitme_cols_keep <- c("IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                     "cellCycle.G2M_UCell", "layer1", "layer2", "layer3")
scatomic_cols <- c("layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                   "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence")
annot_cols <- c(hitme_cols_keep, scatomic_cols)

# Merge with incoming workflow command-line overrides if present
raw_args <- commandArgs(trailingOnly = TRUE)
args <- defaults

if (length(raw_args) > 0) {
  parsed_args_list <- unlist(strsplit(raw_args[1], "__"))
  keys <- parsed_args_list[seq(1, length(parsed_args_list), by = 2)]
  vals <- parsed_args_list[seq(2, length(parsed_args_list), by = 2)]
  
  args_list <- as.list(vals)
  names(args_list) <- keys
  args <- modifyList(defaults, args_list)
}

if (is.null(args$chunk_file) || !file.exists(args$chunk_file)) {
  stop("Valid 'chunk_file' parameter not parsed from execution context!")
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load ref maps ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Load scGate models ####
scGate_models_DB <- get_scGateDB(branch = "dev", verbose = T, force_update = TRUE)
scGate_models_blood <- scGate_models_DB$human$PBMC
scGate_models_blood$MoMac <- scGate_models_blood$Monocyte
scGate_models_blood$Monocyte <- NULL
scGate_models_tumor <- scGate_models_DB$human$HiTME
# Need to have same model names in both lists for HiTME downstream processing)
scGate_models_blood <- c(scGate_models_blood, scGate_models_tumor[!names(scGate_models_tumor) %in% names(scGate_models_blood)])

### Load ProjecTILs ref maps ####
ref.maps_sketched <- list(
  CD8 = load.reference.map(file.path(paths$path_ref, "sketched_CD8T_human_ref_v1.rds")),
  CD4 = load.reference.map(file.path(paths$path_ref, "sketched_CD4T_human_ref_v2.rds")),
  DC = load.reference.map(file.path(paths$path_ref, "sketched_DC_human_ref_v2.rds")),
  MoMac = load.reference.map(file.path(paths$path_ref, "sketched_MoMac_human_v1.rds"))
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Process data ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunk_lines <- readLines(args$chunk_file)
h5ad_file <- chunk_lines[1]
samples_to_process <- chunk_lines[-1]


### Tap h5ad file ####
adata <- ad$read_h5ad(h5ad_file, backed = "r")
obs <- py_to_r(adata$obs)

if (!args$sample_colname %in% colnames(obs)) {
  stop(args$sample_colname, " not found in h5ad obs colnames!")
}

# ==============================================================================
# PROCESSING LOOP
# ==============================================================================
message(paste("--- Starting processing for chunk file:", args$chunk_file, "---"))

chunk_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
annot_file <- file.path(paths$path_output, paste0("annotations_chunk_", chunk_id, ".feather"))
if (file.exists(annot_file)) {
  message(paste("Chunk already processed. Annotations exist at:", annot_file))
} else {
  annotations_list <- list()
  for (target_sample in samples_to_process) {
    message(paste("--- Processing sample:", target_sample, "---"))

    seurat_obj <- get_sample_seurat_obj(
      adata, obs, target_sample, args$sample_colname
    )

    # About 10 minutes per 10,000 cells
    timeout <- max(60, ncol(seurat_obj) / 10000 * 60 * 10)

    ### scATOMIC annotation ####
    if (is.null(seurat_obj@meta.data[["layer_1"]])) {
      for (a in 1:5) {
        message(paste("  scATOMIC attempt", a, "with", round(timeout), "s timeout"))
        result <- tryCatch({
          withTimeout({
            sca_preds <- run_scATOMIC(seurat_obj@assays$RNA$counts)
            sca_results <- create_summary_matrix(
              prediction_list = sca_preds,
              raw_counts = seurat_obj@assays$RNA$counts,
              normal_tissue = args$normal_tissue
            )
            "Complete"
          }, timeout = timeout)
        }, TimeoutException = function(te) {
          message("  scATOMIC timeout, retrying...")
          NULL
        }, error = function(er) {
          message(paste("  scATOMIC error:", er$message, "- retrying..."))
          NULL
        })
        if (!is.null(result)) {
          sca_cols <- intersect(scatomic_cols, colnames(sca_results))
          seurat_obj <- AddMetaData(seurat_obj, sca_results[, sca_cols, drop = FALSE])
          break
        }
      }
    }

    ### HiTME annotation ####
    for (a in 1:5) {
      message(paste("  HiTME attempt", a, "with", timeout, "s timeout"))
      result <- tryCatch({
        withTimeout({
          if (args$tissue_type == "Blood") {
            seurat_obj <- Run.HiTME(
              object = seurat_obj,
              scGate.model = scGate_models_blood,
              ref.maps = ref.maps_sketched,
              verbose = FALSE,
              ncores = 1
            )
          } else {
            seurat_obj <- Run.HiTME(
              object = seurat_obj,
              scGate.model = scGate_models_tumor,
              ref.maps = ref.maps_sketched,
              verbose = FALSE,
              ncores = 1
            )
          }
          "Complete"
        }, timeout = timeout)
      }, TimeoutException = function(te) {
        message("  HiTME timeout, retrying...")
        NULL
      }, error = function(er) {
        message(paste("  HiTME error:", er$message, "- retrying..."))
        NULL
      })
      if (!is.null(result)) break
    }

    ### Extract annotations ####
    meta <- seurat_obj@meta.data
    keep_cols <- intersect(annot_cols, colnames(meta))
    annot <- meta[, keep_cols, drop = FALSE]
    annot$cell_barcode <- rownames(annot)
    annot$sample <- target_sample
    annotations_list[[target_sample]] <- annot

    rm(seurat_obj)
    gc()
  }

  annotations_df <- do.call(rbind, annotations_list)
  rownames(annotations_df) <- NULL
  write_feather(annotations_df, annot_file)
  message(paste("Wrote annotations to:", annot_file))
}

message(paste("---", args$chunk_file, "processing complete! ---"))
RS

log_msg "✓ Chunk ${TASK_ID} processing complete."