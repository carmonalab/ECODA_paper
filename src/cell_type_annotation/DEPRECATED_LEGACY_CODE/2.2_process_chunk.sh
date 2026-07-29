#!/bin/bash
#
# 2.1_process_chunk.sh — Process a specific subset of samples sequentially
#

set -euo pipefail

TASK_ID="$1"
export PROJECT_ROOT="$2"
CHUNK_FILE="$3"

source "${PROJECT_ROOT}/config.env"

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
library(STACAS)
library(ProjecTILs)
library(SignatuR)
library(HiTME)
library(Seurat)

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
env_sample_col  <- Sys.getenv("SAMPLE_COLNAME")
env_tissue      <- Sys.getenv("TISSUE_TYPE")
env_get_comp    <- Sys.getenv("GET_CELLTYPE_COMP")
env_auth_annots <- Sys.getenv("AUTHOR_ANNOT_COLNAMES")

defaults <- list(
  chunk_file            = NULL, 
  sample_colname        = if (env_sample_col != "") env_sample_col else "sample", 
  tissue_type           = if (env_tissue != "") env_tissue else "Tumor", 
  get_celltype_comp     = if (env_get_comp != "") as.logical(env_get_comp) else TRUE,
  author_annot_colnames = if (env_auth_annots != "") unlist(strsplit(env_auth_annots, ",")) else character()
)

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

### Load gene symbols ref ####
geneRef <- data.table::fread(paths$gene_ref)


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
for (target_sample in samples_to_process) {
  message(paste("--- Processing sample:", target_sample, "---"))

  processed_file_path <- file.path(paths$path_output_samples, paste0(target_sample, ".rds"))
  if (file.exists(processed_file_path)) {
    message("Sample already processed. Skipping.")
    next
  }

  seurat_obj <- get_sample_seurat_obj(
    adata, obs, target_sample, args$sample_colname
  )

  ### Standardize gene symbols ####
  seurat_obj <- STACAS::StandardizeGeneSymbols(
    seurat_obj, geneRef,
    assay = "RNA", slots = "counts"
  )

  ### Annotate cells ####
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

  saveRDS(seurat_obj, file = processed_file_path)

  if (args$get_celltype_comp) {
    meta <- seurat_obj@meta.data
    meta <- meta[, !grepl("_UCell|is.pure_|.pseudocounts", colnames(meta))]
    
    # --- Part A: Constant Metadata ---
    # Identify columns where there is only 1 unique value
    is_constant <- sapply(meta, function(x) length(unique(x)) == 1)
    constant_df <- meta[1, is_constant, drop = FALSE]
    constant_df$sample <- target_sample

    # --- Part B: Cell Type Composition ---
    # HiTME
    # Calculate counts for layer1, layer2 and layer3
    # Convert to a dataframe for easier merging later
    comp_1 <- as.data.frame(table(layer1 = meta$layer1))
    comp_2 <- as.data.frame(table(layer2 = meta$layer2))
    comp_3 <- as.data.frame(table(layer3 = meta$layer3))

    ecoda <- list(
      metadata = constant_df,
      cell_counts = list(
        HiTME = list(
          layer1 = comp_1,
          layer2 = comp_2,
          layer3 = comp_3
        )
      )
    )

    # Optional: author annotations
    if (length(args$author_annot_colnames) > 0) {
      for (annot_col in args$author_annot_colnames) {
        if (annot_col %in% colnames(meta)) {
          comp <- as.data.frame(table(temp_name = meta[, annot_col]))
          colnames(comp)[1] <- paste0("authors_", annot_col)
          ecoda$cell_counts[[annot_col]] <- comp
        }
      }
    }
    saveRDS(ecoda, file = file.path(paths$path_output_ecoda, paste0(target_sample, ".rds")))
  }

  rm(seurat_obj)
  gc()
}
message(paste("---", args$chunk_file, "processing complete! ---"))
RS

log_msg "✓ Chunk ${TASK_ID} processing complete."