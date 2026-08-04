#!/bin/bash
set -euo pipefail

# Usage:
#   ./2_submit_hpc_array.sh                # all datasets with chunks
#   ./2_submit_hpc_array.sh <DS_NAME>      # single dataset (must have chunks)

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

DS_NAME_ARG="${1:-}"

module load jq/1.6
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

# Build the list of datasets (all keys of datasets.json, or a single validated key)
DATASET_NAMES=()
if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" >/dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  DATASET_NAMES+=("${DS_NAME_ARG}")
else
  while IFS= read -r name; do
    DATASET_NAMES+=("$name")
  done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")
fi

# ---------------------------------------------------------------------------
# Build the global chunk manifest: one tab-separated line per chunk across all
# datasets:  DS_NAME<TAB>absolute_chunk_path
# SLURM_ARRAY_TASK_ID maps 1:1 to manifest line numbers, so task IDs are
# globally unique and per-chunk feather names never collide across datasets.
# The manifest is regenerated on every run; chunk dirs are recreated fresh by
# 1.1_prepare_chunks.py, so a stale manifest cannot be misread.
# ---------------------------------------------------------------------------
export CHUNKS_MANIFEST="${SCRATCH_OUTPUT_DIR}/chunks_manifest.txt"
: > "${CHUNKS_MANIFEST}"

shopt -s nullglob

SKIPPED_DATASETS=()
TOTAL_CHUNKS=0

for DS_NAME in "${DATASET_NAMES[@]}"; do
  CHUNKS_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks"
  CHUNK_FILES=("${CHUNKS_DIR}"/chunk_*.txt)
  NUM_CHUNKS=${#CHUNK_FILES[@]}

  if [[ ${NUM_CHUNKS} -eq 0 ]]; then
    if [[ -n "${DS_NAME_ARG}" ]]; then
      echo "ERROR: No chunk files found in ${CHUNKS_DIR}! Run 1_prepare_chunks.sh first."
      exit 1
    fi
    echo "WARNING: No chunk files found in ${CHUNKS_DIR}; skipping ${DS_NAME}."
    SKIPPED_DATASETS+=("${DS_NAME}")
    continue
  fi

  for CHUNK_FILE in "${CHUNK_FILES[@]}"; do
    printf '%s\t%s\n' "${DS_NAME}" "${CHUNK_FILE}" >> "${CHUNKS_MANIFEST}"
  done
  TOTAL_CHUNKS=$((TOTAL_CHUNKS + NUM_CHUNKS))
done

if [[ ${TOTAL_CHUNKS} -eq 0 ]]; then
  echo "ERROR: No chunk files found in any dataset. Run 1_prepare_chunks.sh first."
  exit 1
fi

echo "Manifest written to ${CHUNKS_MANIFEST} with ${TOTAL_CHUNKS} chunks."
if [[ ${#SKIPPED_DATASETS[@]} -gt 0 ]]; then
  echo "Skipped datasets (no chunks): ${SKIPPED_DATASETS[*]}"
fi

echo "Found ${TOTAL_CHUNKS} chunks. Submitting job array range 1-${TOTAL_CHUNKS} to SLURM..."
SUBMIT_MSG=$(sbatch \
    --array=1-${TOTAL_CHUNKS}%${MAX_NUM_CHUNKS_PARALLEL} \
    --output="${LOGS_DIR}/4_cell_type_annotation_%A_%a.log" \
    --error="${LOGS_DIR}/4_cell_type_annotation_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    "$(dirname "${BASH_SOURCE[0]}")/2.1_run_worker.sh")

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
