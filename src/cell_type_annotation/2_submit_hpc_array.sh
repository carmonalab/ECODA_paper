#!/bin/bash
set -euo pipefail

# IMPORTANT: DS_NAME must be exported before calling this script, e.g.:
#   export DS_NAME="Stephenson"
#   ./2_submit_hpc_array.sh

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

if [[ -z "${DS_NAME:-}" ]]; then
  echo "ERROR: DS_NAME is not set. Export it before running this script."
  echo "Usage: export DS_NAME=\"DatasetName\" && $0"
  exit 1
fi

export HOME_CHUNKS_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks"
NUM_CHUNKS=$(ls -1 "${HOME_CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)

if [[ ${NUM_CHUNKS} -eq 0 ]]; then
  echo "ERROR: No chunk files found in ${HOME_CHUNKS_DIR}! Run 1_prepare_chunks.sh first."
  exit 1
fi

echo "Found ${NUM_CHUNKS} chunks. Submitting job array range 1-${NUM_CHUNKS} to SLURM..."
SUBMIT_MSG=$(sbatch \
    --array=1-${NUM_CHUNKS}%${MAX_NUM_CHUNKS_PARALLEL} \
    --output="${PROJECT_ROOT}/logs/chunk_%A_%a.log" \
    --error="${PROJECT_ROOT}/logs/chunk_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    2.1_run_worker.sh)

ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
echo "Array Job ID allocated: ${ARRAY_JOB_ID}"

# ==============================================================================
# Post-Pipeline Sync: Run locally on Login Node because compute nodes lack NAS access
# ==============================================================================
echo "Monitoring Array Job ${ARRAY_JOB_ID} for completion..."

while squeue -u "$USER" 2>/dev/null | grep -q "${ARRAY_JOB_ID}"; do
    sleep 60
done

echo "Array Job ${ARRAY_JOB_ID} finished. Starting local sync to NAS from login node..."

if ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    mkdir -p "${NAS_TARGET_DIR}/output"
    rsync -rlptDv "${SCRATCH_OUTPUT_DIR}/" "${NAS_TARGET_DIR}/output/"
    echo "Success: Data safely synchronized to the NAS."
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable even from this login node."
    exit 1
fi
