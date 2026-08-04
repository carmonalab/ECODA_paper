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
source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
mkdir -p "${LOGS_DIR}"

# Check if the user requested test mode
MODE_ARG="${1:-"production"}"
R_PASS_ARG="test__false"

if [ "$MODE_ARG" = "test" ]; then
  echo ">>> CONFIGURING PIPELINE IN TEST MODE <<<"
  R_PASS_ARG="test__true"
else
  echo ">>> CONFIGURING PIPELINE IN PRODUCTION MODE <<<"
fi

# DS_NAME must be set (exported by 2_submit_hpc_array.sh or set manually)
if [[ -z "${DS_NAME:-}" ]]; then
  echo "ERROR: DS_NAME is not set. Export it before running this script."
  exit 1
fi

# -------------------------------------------------------------------------
# STAGE REFERENCE MAPS: Copy ref files from UNIGE NAS to Cluster Scratch
# -------------------------------------------------------------------------

if [ -d "${HOME_REF_DIR}" ] && [ "$(ls -A "${HOME_REF_DIR}" 2>/dev/null)" ]; then
  echo ">>> Reference maps already exist in ${HOME_REF_DIR}. Skipping rsync. <<<"
else
  echo "Staging reference maps from NAS to home directory..."
  mkdir -p "${HOME_REF_DIR}"
  rsync -av --progress "${NAS_REF_DIR}" "${HOME_REF_DIR}/"
fi

echo "Allocating short-lived compute session to build dataset chunks..."
LOG_FILE="${LOGS_DIR}/prepare_chunks_${MODE_ARG}.log"
ENV_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/default/bin/Rscript"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export R_LIBS_SITE="${PIXI_R_LIB}:${R_LIBS_SITE:-}"

srun --partition=shared-cpu \
     --time=00:10:00 \
     --ntasks=1 \
     --cpus-per-task=1 \
     --mem=4G \
     --output="${LOG_FILE}" \
     --error="${LOG_FILE}" \
     "${ENV_RSCRIPT}" --vanilla "${SCRIPT_DIR}/1.1_prepare_chunks.r" "${R_PASS_ARG}"

echo "✓ Chunk generation script finished executing. Log saved to: logs/prepare_chunks_${MODE_ARG}.log"
