#!/bin/bash
#
# 1_prepare_chunks.sh — Build dataset chunks (Supports production or test modes)
#
# Usage:
#   ./1_prepare_chunks.sh        # Runs in production mode (chunk_size = 5)
#   ./1_prepare_chunks.sh test   # Runs in test mode (chunk_size = 1)
#

set -euo pipefail

# 1. Load central config
source "$(dirname "${BASH_SOURCE[0]}")/config.env"
cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/logs"

# Check if the user requested test mode
MODE_ARG="${1:-"production"}"
R_PASS_ARG="test__false"

if [ "$MODE_ARG" = "test" ]; then
  echo ">>> CONFIGURING PIPELINE IN TEST MODE <<<"
  R_PASS_ARG="test__true"
else
  echo ">>> CONFIGURING PIPELINE IN PRODUCTION MODE <<<"
fi

# -------------------------------------------------------------------------
# STAGE DATA: Copy raw files from UNIGE NAS to Cluster Scratch
# -------------------------------------------------------------------------

# Check if directory exists and is not empty
if [ -d "${HOME_DATA_DIR}" ] && [ "$(ls -A "${HOME_DATA_DIR}" 2>/dev/null)" ]; then
  echo ">>> Data files already exist in ${HOME_DATA_DIR}. Skipping rsync file transfer. <<<"
else
  echo "Staging raw files from NAS to home directory..."
  mkdir -p "${HOME_DATA_DIR}"
  rsync -av --progress "${NAS_DATA_DIR}" "${HOME_DATA_DIR}"
fi

# -------------------------------------------------------------------------
# STAGE REFERENCE MAPS: Copy ref files from UNIGE NAS to Cluster Scratch
# -------------------------------------------------------------------------

if [ -d "${HOME_REF_DIR}" ] && [ "$(ls -A "${HOME_REF_DIR}" 2>/dev/null)" ]; then
  echo ">>> Reference maps already exist in ${HOME_REF_DIR}. Skipping rsync. <<<"
else
  echo "Staging reference maps from NAS to home directory..."
  mkdir -p "${HOME_REF_DIR}"
  rsync -av --progress "${NAS_REF_DIR}" "${HOME_REF_DIR}"
fi

# -------------------------------------------------------------------------
# STAGE GENE ANNOTATIONS: Download Ensembl gene reference if missing
# -------------------------------------------------------------------------

if [ -f "${GENE_REF_FILE}" ]; then
  echo ">>> Gene reference file already exists. Skipping download. <<<"
else
  echo "Downloading gene reference annotations from GitHub..."
  curl -sSL "${GENE_REF_URL}" -o "${GENE_REF_FILE}"
fi
# -------------------------------------------------------------------------

echo "Allocating short-lived compute session to build dataset chunks..."
LOG_FILE="${PROJECT_ROOT}/logs/prepare_chunks_${MODE_ARG}.log"
ENV_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/default/bin/Rscript"

# Dynamically construct the Pixi library path using your config variables
PIXI_R_LIB="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"

# Export variables so renv and R recognize the Pixi environment during srun
export RENV_CONFIG_EXTERNAL_LIBRARIES="${PIXI_R_LIB}"
export R_LIBS_SITE="${PIXI_R_LIB}:${R_LIBS_SITE:-}"

srun --partition=shared-cpu \
     --time=00:10:00 \
     --ntasks=1 \
     --cpus-per-task=1 \
     --mem=4G \
     --output="${LOG_FILE}" \
     --error="${LOG_FILE}" \
     "${ENV_RSCRIPT}" --vanilla 1_prepare_chunks.r "${R_PASS_ARG}"

echo "✓ Chunk generation script finished executing. Log saved to: logs/prepare_chunks_${MODE_ARG}.log"