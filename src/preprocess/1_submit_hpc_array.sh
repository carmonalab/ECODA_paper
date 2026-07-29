#!/bin/bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

# ---------------------------------------------------------------------------
# Phase 1: Stage data from NAS to HPC scratch (login node has NAS access)
# ---------------------------------------------------------------------------
echo "=== Phase 1: Staging data from NAS to HPC scratch ==="
mkdir -p "${HPC_SCRATCH_DIR}"

jq -r '
  to_entries[] |
  .key as $key |
  .value.folder_name as $folder |
  .value.views |
  to_entries[] |
  .value.input_file_name |
  if type == "array" then .[] else . end |
  "\($folder) \(.)"
' "${DATASETS_JSON_FILE}" | sort -u | \
while read -r FOLDER_NAME RAW_FILE_NAME; do
    if [ "$FOLDER_NAME" == "null" ] || [ -z "$FOLDER_NAME" ]; then
        continue
    fi

    NAS_FILE_PATH="${NAS_SC_DIR}/${FOLDER_NAME}/output/${RAW_FILE_NAME}"
    echo "Dataset folder: ${FOLDER_NAME}, file: ${RAW_FILE_NAME}"

    if [ -f "$NAS_FILE_PATH" ]; then
        rsync -ah --progress "$NAS_FILE_PATH" "$HPC_SCRATCH_DIR"
    else
        echo "  -> [WARNING] Source not found: ${NAS_FILE_PATH}"
    fi
done

echo "Data staging complete."
echo ""

# ---------------------------------------------------------------------------
# Phase 2: Submit preprocessing array job
# ---------------------------------------------------------------------------
echo "=== Phase 2: Submitting preprocessing array job ==="

module load GCCcore/12.2.0
module load jq/1.6

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
    --output="${PROJECT_ROOT}/logs/preprocess_%A_%a.log" \
    --error="${PROJECT_ROOT}/logs/preprocess_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    "$(dirname "${BASH_SOURCE[0]}")/1.1_run_worker.sh")

ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
echo "Array Job ID allocated: ${ARRAY_JOB_ID}"

# ---------------------------------------------------------------------------
# Phase 3: Monitor & sync results back to NAS
# ---------------------------------------------------------------------------
echo "=== Phase 3: Monitoring job completion ==="
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
