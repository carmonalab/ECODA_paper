#!/bin/bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.env"
cd "${PROJECT_ROOT}"

HOME_CHUNKS_DIR="${HOME}/${DS_NAME}/output/chunks"
NUM_CHUNKS=$(ls -1 "${HOME_CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)

if [[ ${NUM_CHUNKS} -eq 0 ]]; then
  echo "ERROR: No chunk files found in ${HOME_CHUNKS_DIR}! Run 1_prepare_chunks.r first."
  exit 1
fi

CONFIG_FILE="${PROJECT_ROOT}/pipeline_config.json"
TEST_MODE=$(grep -o '"test_mode": [^,]*' "${CONFIG_FILE}" 2>/dev/null | cut -d' ' -f2 || echo "false")

if [ "$TEST_MODE" = "true" ]; then
  echo "[TEST MODE DETECTED] Restricting array execution to first 3 chunks."
  ARRAY_LIMIT=3
else
  ARRAY_LIMIT=${NUM_CHUNKS}
fi

echo "Submitting job array range 1-${ARRAY_LIMIT} to SLURM using 2.1_run_worker.sh..."
SUBMIT_MSG=$(sbatch \
    --array=1-${ARRAY_LIMIT}%${MAX_NUM_CHUNKS_PARALLEL} \
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

# Wait until the SLURM array job is no longer active (Pending or Running)
# A more robust check that won't break on array suffix changes
while squeue -u "$USER" 2>/dev/null | grep -q "${ARRAY_JOB_ID}"; do
    sleep 60
done

echo "Array Job ${ARRAY_JOB_ID} finished. Starting local sync to NAS from login node..."

# Run the sync directly on the login node where the mount lives
if ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    mkdir -p "${NAS_TARGET_DIR}/output"
    rsync -rlptDv "${SCRATCH_OUTPUT_DIR}/" "${NAS_TARGET_DIR}/output/"
    echo "✓ Success: Data safely synchronized to the NAS."
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable even from this login node."
    exit 1
fi