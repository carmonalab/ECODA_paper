#!/bin/bash
#
# 1_prepare_chunks.sh — Build dataset chunks (Supports production or test modes)
#
# Usage:
#   ./1_prepare_chunks.sh                        # production, all datasets (chunk_size = 5)
#   ./1_prepare_chunks.sh test                   # test mode, all datasets (chunk_size = 1)
#   ./1_prepare_chunks.sh production <DS_NAME>   # production, single dataset
#   ./1_prepare_chunks.sh test <DS_NAME>         # test mode, single dataset
#

set -euo pipefail

# 1. Load central config
source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
mkdir -p "${LOGS_DIR}"

# 2. Parse arguments: MODE (production/test) + optional single dataset
MODE_ARG="${1:-production}"
DS_NAME_ARG="${2:-}"
PY_ARGS=""

if [ "$MODE_ARG" = "test" ]; then
  echo ">>> CONFIGURING PIPELINE IN TEST MODE <<<"
  PY_ARGS="--test"
elif [ "$MODE_ARG" = "production" ]; then
  echo ">>> CONFIGURING PIPELINE IN PRODUCTION MODE <<<"
else
  echo "ERROR: Unknown mode '${MODE_ARG}'. Use 'production' or 'test'."
  exit 1
fi

module load jq/1.6
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

# 3. Build the list of datasets (all keys of datasets.json, or a single
#    validated key when DS_NAME_ARG is given)
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

if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: No datasets found in ${DATASETS_JSON_FILE}."
  exit 1
fi
echo "Datasets to process: ${DATASET_NAMES[*]}"

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

ENV_PYTHON="${PROJECT_ROOT}/.pixi/envs/default/bin/python"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 4. Build chunks per dataset (sequential, one short-lived compute session each;
#    --time=00:30:00 covers the full dataset loop, per-dataset runs are usually
#    much shorter)
FAILED_DATASETS=()
SKIPPED_DATASETS=()

for DS_NAME in "${DATASET_NAMES[@]}"; do
  echo ""
  echo ">>> Building chunks for dataset: ${DS_NAME} <<<"

  # Skip datasets without preprocessed .h5ad input (e.g. Zhu has no views)
  if ! ls "${HPC_SCRATCH_DIR}/${DS_NAME}/output"/*.h5ad >/dev/null 2>&1; then
    echo "WARNING: No preprocessed .h5ad files in ${HPC_SCRATCH_DIR}/${DS_NAME}/output; skipping ${DS_NAME}."
    SKIPPED_DATASETS+=("${DS_NAME}")
    continue
  fi

  LOG_FILE="${LOGS_DIR}/prepare_chunks_${MODE_ARG}_${DS_NAME}.log"
  export DS_NAME

  echo "Allocating short-lived compute session to build ${DS_NAME} chunks..."
  if ! srun --partition=shared-cpu \
       --time=00:30:00 \
       --ntasks=1 \
       --cpus-per-task=1 \
       --mem=4G \
       --output="${LOG_FILE}" \
       --error="${LOG_FILE}" \
       "${ENV_PYTHON}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" ${PY_ARGS}; then
    echo "ERROR: Chunk generation failed for ${DS_NAME}. See ${LOG_FILE}."
    FAILED_DATASETS+=("${DS_NAME}")
    continue
  fi

  echo "✓ Chunks generated for ${DS_NAME}. Log saved to: ${LOG_FILE}"
done

echo ""
echo "=== Chunk preparation summary ==="
echo "Processed: ${#DATASET_NAMES[@]} datasets"
echo "Skipped (no preprocessed .h5ad): ${#SKIPPED_DATASETS[@]}"
if [[ ${#SKIPPED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${SKIPPED_DATASETS[*]}"
fi
echo "Failed: ${#FAILED_DATASETS[@]}"
if [[ ${#FAILED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${FAILED_DATASETS[*]}"
  echo "Resolve the failures, then rerun this script before 2_submit_hpc_array.sh."
  exit 1
fi
