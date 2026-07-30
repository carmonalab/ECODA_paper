#!/bin/bash
#
# Centralized SLURM and project environment variables
#
# Source this from any bash script that needs project paths or SLURM config:
#   source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
# (adjust relative path based on script depth)
#

# --- Project Paths ---
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"

# --- NAS Paths ---
export NAS_PREFIX="/srv/smednas515.unige.ch/carmona_smb"
NAS_BASE_DIR="${NAS_PREFIX}/DataCollections"
NAS_SC_DIR="${NAS_BASE_DIR}/Standardized_SingleCell_Datasets"
NAS_TARGET_DIR="${NAS_BASE_DIR}/AnalysisResults/ECODA"

# --- HPC Scratch Paths ---
HPC_SCRATCH_DIR="${HOME}/scratch/ECODA_paper"
SCRATCH_OUTPUT_DIR="${HPC_SCRATCH_DIR}/output"

# --- Pixi R library path ---
PIXI_R_LIB="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"

# --- Reference atlas paths (cell type annotation) ---
NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"

# --- Gene reference (cell type annotation) ---
# Used for gene standardization with STACAS; now implemented in src/preprocess/1.1.1_preprocess.py
# which runs before cell type annotation. Kept for potential fallback.
GENE_REF_FILE="${PROJECT_ROOT}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
GENE_REF_URL="https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"

# --- SLURM Configuration ---
SLURM_ACCOUNT=""
SLURM_PARTITION="shared-cpu" # TODO: Adapt for specific pipelines

# --- User Info ---
USER_EMAIL="${USER}@unige.ch"

# --- Parallelism ---
MAX_NUM_CHUNKS_PARALLEL=500

# --- Per-dataset Chunk Directory ---
# HOME_CHUNKS_DIR is dataset-specific; set it before calling array submission scripts.
# Typically: HOME_CHUNKS_DIR="${HOME}/${DS_NAME}/output/chunks"
HOME_CHUNKS_DIR=""
