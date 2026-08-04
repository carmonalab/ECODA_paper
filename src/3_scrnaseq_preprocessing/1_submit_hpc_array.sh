#!/bin/bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load GCCcore/12.2.0
module load jq/1.6

# ---------------------------------------------------------------------------
# Submit preprocessing array job
# (Raw-data staging now lives in src/1_stage_data/1_stage_data.sh)
# ---------------------------------------------------------------------------
echo "=== Submitting preprocessing array job ==="

DATASET_NAMES=()
while IFS= read -r name; do
  DATASET_NAMES+=("$name")
done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")

NUM_DATASETS=${#DATASET_NAMES[@]}
if [[ ${NUM_DATASETS} -eq 0 ]]; then
  echo "ERROR: No datasets found in ${DATASETS_JSON_FILE}."
  exit 1
fi

echo "Found ${NUM_DATASETS} datasets."
echo "Datasets: ${DATASET_NAMES[*]}"

SUBMIT_MSG=$(sbatch \
    --array=1-${NUM_DATASETS}%${MAX_NUM_CHUNKS_PARALLEL} \
    --output="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.log" \
    --error="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    "$(dirname "${BASH_SOURCE[0]}")/1.1_run_worker.sh")

ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
echo "Array Job ID allocated: ${ARRAY_JOB_ID}"

# ---------------------------------------------------------------------------
# Monitor & sync results back to NAS
# ---------------------------------------------------------------------------
echo "=== Monitoring job completion ==="
while squeue -u "$USER" 2>/dev/null | grep -q "${ARRAY_JOB_ID}"; do
    sleep 60
done

echo "Array Job ${ARRAY_JOB_ID} finished. Syncing results to NAS..."
if ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    mkdir -p "${NAS_TARGET_DIR}/output"
    rsync -rlptDv "${SCRATCH_OUTPUT_DIR}/" "${NAS_TARGET_DIR}/output/"
    echo "Results synchronized to ${NAS_TARGET_DIR}/output/"
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
    exit 1
fi
