#!/bin/bash
#
# Centralized SLURM and project environment variables
#
# Source this from any bash script that needs project paths or SLURM config:
#   source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
# (adjust relative path based on script depth)
#

# --- Project Paths ---
# NOTE: All path/env vars below are exported so that srun/sbatch children and
# Sys.getenv()/os.environ consumers (R config_helper.R, Python scripts) see them.
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
export LOGS_DIR="${PROJECT_ROOT}/logs"

# --- NAS Paths ---
export NAS_PREFIX="/srv/smednas515.unige.ch/carmona_smb"
NAS_BASE_DIR="${NAS_PREFIX}/DataCollections"
export NAS_SC_DIR="${NAS_BASE_DIR}/Standardized_SingleCell_Datasets"
export NAS_TARGET_DIR="${NAS_PREFIX}/Projects/ECODA_paper"

# --- HPC Scratch Paths ---
export HPC_SCRATCH_DIR="${HOME}/scratch/ECODA_paper"
export SCRATCH_OUTPUT_DIR="${HPC_SCRATCH_DIR}/output"

# --- Pixi R library path ---
export PIXI_R_LIB="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"

# --- Reference atlas paths (cell type annotation) ---
export NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
export HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"

# --- Sample column name (cell type annotation) ---
# 1.1.1_preprocess.py standardizes every dataset's sample column to "Sample".
export SAMPLE_COLNAME="Sample"

# --- SLURM Configuration ---
SLURM_ACCOUNT=""
SLURM_PARTITION="shared-cpu" # TODO: Adapt for specific pipelines

# --- User Info ---
USER_EMAIL="${USER}@unige.ch"

# --- Parallelism ---
MAX_NUM_CHUNKS_PARALLEL=500

# --- Per-dataset Chunk Directory ---
# HOME_CHUNKS_DIR is dataset-specific; set it before calling array submission scripts.
# Typically: HOME_CHUNKS_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks" (see 2_submit_hpc_array.sh)
HOME_CHUNKS_DIR=""
