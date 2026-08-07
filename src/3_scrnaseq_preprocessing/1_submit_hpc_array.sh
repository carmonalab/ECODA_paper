#!/bin/bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

# ---------------------------------------------------------------------------
# Submit preprocessing array job
# (Raw-data staging now lives in src/1_stage_data/1_stage_data.sh)
#
# Usage:
#   ./1_submit_hpc_array.sh                  # all datasets (skips keys starting with _)
#   ./1_submit_hpc_array.sh --ds_name _debug # single dataset (explicitly allowed)
#   ./1_submit_hpc_array.sh --ds_name _debug --force   # recompute existing outputs
#
# Array task IDs map 1:1 to jq 'keys[]' line numbers (see 1.1_run_worker.sh).
# Convention: default-all skips keys starting with "_" (e.g. _debug) unless
# explicitly requested via --ds_name. Invariant relied on by the task->dataset
# mapping: every non-"_" key sorts BEFORE any "_" key (jq sorts by codepoint;
# "_" = 0x5F > A-Z = 0x41-0x5A, and all real dataset keys start with a capital
# letter), so the non-underscore keys are exactly a PREFIX of the sorted keys:
#   - default mode: array 1..N over the prefix (worker sed -n resolves the
#     i-th non-underscore key automatically)
#   - single mode: array INDEX..INDEX where INDEX is the dataset's position in
#     the FULL sorted key list (worker sed -n resolves it without any change)
# ---------------------------------------------------------------------------
echo "=== Submitting preprocessing array job ==="

DS_NAME_ARG=""
FORCE_ARG=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ds_name)
      DS_NAME_ARG="${2:-}"
      shift 2
      ;;
    --ds_name=*)
      DS_NAME_ARG="${1#*=}"
      shift
      ;;
    --force)
      FORCE_ARG=1
      shift
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Passed to workers via the environment (sbatch propagates the submit script's
# environment); 1.1_run_worker.sh forwards it to 1.1.1_preprocess.py --force.
export FORCE_PREPROCESS="${FORCE_ARG}"

DATASET_NAMES=()
if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  DATASET_NAMES+=("${DS_NAME_ARG}")
else
  while IFS= read -r name; do
    DATASET_NAMES+=("$name")
  done < <(jq -r 'keys[] | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
fi

NUM_DATASETS=${#DATASET_NAMES[@]}
if [[ ${NUM_DATASETS} -eq 0 ]]; then
  echo "ERROR: No datasets found in ${DATASETS_JSON_FILE}."
  exit 1
fi

echo "Found ${NUM_DATASETS} datasets."
echo "Datasets: ${DATASET_NAMES[*]}"

# Map single-dataset mode to a 1-task array at the dataset's position in the
# FULL sorted key list (includes "_" keys, which sort last).
if [[ -n "${DS_NAME_ARG}" ]]; then
  DS_INDEX="$(jq -r --arg ds "${DS_NAME_ARG}" '[keys[]] | index($ds) + 1' "${DATASETS_JSON_FILE}")"
  ARRAY_SPEC="${DS_INDEX}-${DS_INDEX}"
else
  ARRAY_SPEC="1-${NUM_DATASETS}"
fi
echo "Array specification: ${ARRAY_SPEC}"

mkdir -p "${LOGS_DIR}"
SUBMIT_MSG=$(sbatch \
    --array="${ARRAY_SPEC}%${MAX_NUM_CHUNKS_PARALLEL}" \
    --partition="${SLURM_PARTITION}" \
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

# Give sacct a moment to record final states (mirrors 1_submit_hpc.sh)
sleep 30

echo "Array Job ${ARRAY_JOB_ID} finished. Checking task states..."
STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
if [[ -z "${STATES//[[:space:]]/}" ]]; then
    echo "ERROR: sacct returned no states for Array Job ${ARRAY_JOB_ID}; NOT syncing to NAS."
    exit 1
fi
# Fail-closed: every row (array master + tasks + batch steps) must be COMPLETED.
if grep -qvE '^ *COMPLETED *$' <<< "${STATES}"; then
    echo "ERROR: Array Job ${ARRAY_JOB_ID} had non-COMPLETED tasks; NOT syncing to NAS."
    sacct -j "${ARRAY_JOB_ID}" --format=JobID,JobName,State,ExitCode
    exit 1
fi

echo "All tasks completed successfully. Syncing results to NAS..."
SYNCED_COUNT=0
if ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    if [[ -n "${DS_NAME_ARG}" ]]; then
      # Single-dataset mode: sync only the requested dataset's output dir.
      SYNC_DIRS=("${HPC_SCRATCH_DIR}/${DS_NAME_ARG}/output")
    else
      SYNC_DIRS=("${HPC_SCRATCH_DIR}"/*/output)
    fi
    for DS_DIR in "${SYNC_DIRS[@]}"; do
      [[ -d "${DS_DIR}" ]] || continue
      DS_NAME="$(basename "$(dirname "${DS_DIR}")")"
      mkdir -p "${NAS_TARGET_DIR}/${DS_NAME}/output"
      rsync -rlptDv "${DS_DIR}/" "${NAS_TARGET_DIR}/${DS_NAME}/output/"
      SYNCED_COUNT=$((SYNCED_COUNT + 1))
    done
    if [[ ${SYNCED_COUNT} -eq 0 ]]; then
        echo "ERROR: No dataset output dirs found under ${HPC_SCRATCH_DIR}; nothing to sync."
        exit 1
    fi
    echo "Results synchronized to ${NAS_TARGET_DIR}/<DS_NAME>/output/ (${SYNCED_COUNT} datasets)"
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
    exit 1
fi
