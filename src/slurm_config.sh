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

# --- Pixi-managed interpreter commands ---
# PYTHON_BIN: pixi python (has scanpy/anndata etc.); plain `python` on worker
# nodes may resolve to a bare system python without these packages.
# PIXI_RSCRIPT: pixi Rscript entry point (word-splits into the pixi command).
export PYTHON_BIN="${PROJECT_ROOT}/.pixi/envs/default/bin/python"
export PIXI_RSCRIPT="${HOME}/.pixi/bin/pixi run Rscript --vanilla"

# --- Reference atlas paths (cell type annotation) ---
export NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
export HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"

# --- scGate model DB cache (cell type annotation) ---
# Created once by 2.0_create_scgate_db.R (via 2_submit_hpc_array.sh); loaded by
# 2.1.1_process_chunk.R workers so they do not download in parallel.
# SCGATE_DB_BRANCH is the single source of truth for the model DB version.
export SCGATE_DB_PATH="${PROJECT_ROOT}/aux/scGateDB.rds"
export SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"

# --- Sample column name (cell type annotation) ---
# 1.1.1_preprocess.py standardizes every dataset's sample column to "Sample".
export SAMPLE_COLNAME="Sample"

# --- SLURM Configuration ---
# Passed at submit time via `--partition="${SLURM_PARTITION}"` (sbatch
# directives do not expand variables). Override per pipeline if needed.
SLURM_PARTITION="shared-cpu" # TODO: Adapt for specific pipelines

# --- User Info ---
USER_EMAIL="${USER}@unige.ch"

# --- Parallelism ---
MAX_NUM_CHUNKS_PARALLEL=500
