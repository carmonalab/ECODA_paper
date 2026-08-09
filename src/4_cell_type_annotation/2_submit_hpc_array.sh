#!/bin/bash
set -euo pipefail

# Usage:
#   ./2_submit_hpc_array.sh                # all datasets with chunks
#   ./2_submit_hpc_array.sh <DS_NAME>      # single dataset (must have chunks)
#   ./2_submit_hpc_array.sh <DS_NAME> --sync-only <job-id>  # resume: skip submission, re-check + sync

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils/bash/sync_status_email.sh"
mkdir -p "${LOGS_DIR}"

DS_NAME_ARG="${1:-}"
shift 2>/dev/null || true

SYNC_ONLY_IDS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sync-only)
      SYNC_ONLY_IDS="${2:-}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires a job id."
        exit 1
      fi
      shift 2
      ;;
    --sync-only=*)
      SYNC_ONLY_IDS="${1#*=}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires a job id."
        exit 1
      fi
      shift
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

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

if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: job ${SYNC_ONLY_IDS} (no scGate cache, no manifest, no submission) ==="
  ARRAY_JOB_ID="${SYNC_ONLY_IDS}"
else
# -------------------------------------------------------------------------
# STAGE scGate MODEL DB CACHE: Download the scGate model DB once into aux/ so
# the annotation array workers load it from disk instead of downloading in
# parallel (up to MAX_NUM_CHUNKS_PARALLEL concurrent downloads). Failure is
# non-fatal: workers fall back to download + persist themselves (see
# 2.1.1_process_chunk.R).
# -------------------------------------------------------------------------
if [[ -f "${SCGATE_DB_PATH}" ]]; then
  echo ">>> scGate DB cache already exists at ${SCGATE_DB_PATH}. Skipping. <<<"
else
  echo "Creating scGate DB cache at ${SCGATE_DB_PATH} (one-time download)..."
  if ! srun --partition="${SLURM_PARTITION}" \
       --time=00:30:00 \
       --ntasks=1 \
       --cpus-per-task=1 \
       --mem=4G \
       --output="${LOGS_DIR}/prepare_scgatedb.log" \
       --error="${LOGS_DIR}/prepare_scgatedb.log" \
       ${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.0_create_scgate_db.R"; then
    echo "WARNING: scGate DB cache creation failed (see ${LOGS_DIR}/prepare_scgatedb.log)."
    echo "         Continuing: annotation workers will download the model DB themselves (slower)."
  else
    echo "✓ scGate DB cache created. Log saved to: ${LOGS_DIR}/prepare_scgatedb.log"
  fi
fi

# ---------------------------------------------------------------------------
# Build the global chunk manifest: one tab-separated line per chunk across all
# datasets:  DS_NAME<TAB>absolute_chunk_path
# SLURM_ARRAY_TASK_ID maps 1:1 to manifest line numbers, so task IDs are
# globally unique and per-chunk feather names never collide across datasets.
# The manifest is regenerated on every run; chunk dirs are recreated fresh by
# 1.1_prepare_chunks.py, so a stale manifest cannot be misread.
# ---------------------------------------------------------------------------
export CHUNKS_MANIFEST="${HPC_SCRATCH_DIR}/chunks_manifest.txt"
: > "${CHUNKS_MANIFEST}"

shopt -s nullglob

SKIPPED_DATASETS=()
TOTAL_CHUNKS=0

for DS_NAME in "${DATASET_NAMES[@]}"; do
  CHUNKS_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks"
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
    --partition="${SLURM_PARTITION}" \
    --output="${LOGS_DIR}/4_cell_type_annotation_%A_%a.log" \
    --error="${LOGS_DIR}/4_cell_type_annotation_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    "${SCRIPT_DIR}/2.1_run_worker.sh")

ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
echo "Array Job ID allocated: ${ARRAY_JOB_ID}"
fi

# ==============================================================================
# Post-Pipeline Sync: Run locally on Login Node because compute nodes lack NAS access
# ==============================================================================
echo "Monitoring Array Job ${ARRAY_JOB_ID} for completion..."

# Block until the job leaves the scheduler. scontrol has NO plain `wait`
# command (only `wait_job`, which waits for node-ready — not completion —
# and is documented as unusable with SLURM_ARRAY_JOB_ID), so poll squeue
# for the exact job id (`-o %A` prints the array master id for every task,
# or the job id for plain jobs). The fail-closed sacct gate below is the
# authoritative check (covers cancellation, failure, purged controller
# records).
while squeue -u "$USER" -h -o "%A" 2>/dev/null | grep -qx "${ARRAY_JOB_ID}"; do
    sleep 60
done
echo "Array Job ${ARRAY_JOB_ID} left the scheduler."

# sacct may lag a few seconds behind the job leaving the scheduler; poll
# (bounded) until every state row is terminal instead of a blind fixed sleep.
# The 180-iteration cap (15 min) plus a 60-iteration grace window (5 min)
# covers pathological SlurmDBD accounting lag (scheduler said done, sacct
# still reports RUNNING); the fail-closed gate below is unchanged.
echo "Waiting for sacct to record terminal states for job ${ARRAY_JOB_ID} (bounded, max 20 min)..."
TAIL_ITER=0
while (( TAIL_ITER < 180 )); do  # max 15 min at 5s
    STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
    if [[ -n "${STATES//[[:space:]]/}" ]] \
       && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
        break
    fi
    sleep 5
    TAIL_ITER=$((TAIL_ITER + 1))
done
if (( TAIL_ITER >= 180 )); then
    echo "WARNING: sacct still reporting non-terminal states after 15 min; extending wait by a 5 min grace window..." >&2
    TAIL_ITER=0
    while (( TAIL_ITER < 60 )); do
        STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
        if [[ -n "${STATES//[[:space:]]/}" ]] \
           && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
            break
        fi
        sleep 5
        TAIL_ITER=$((TAIL_ITER + 1))
    done
fi

echo "Array Job ${ARRAY_JOB_ID} finished. Checking task states..."
STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
if [[ -z "${STATES//[[:space:]]/}" ]]; then
    echo "ERROR: sacct returned no states for Array Job ${ARRAY_JOB_ID}; NOT syncing to NAS."
    notify_sync_status \
        "ECODA: annotation NOT synced (job ${ARRAY_JOB_ID})" \
        "Annotation sync to NAS failed for job ${ARRAY_JOB_ID} (datasets: ${DATASET_NAMES[*]}): sacct returned no states (job purged or unknown id)."
    exit 1
fi
# Fail-closed: every row (array master + tasks + batch steps) must be COMPLETED.
if grep -qvE '^ *COMPLETED *$' <<< "${STATES}"; then
    echo "ERROR: Array Job ${ARRAY_JOB_ID} had non-COMPLETED tasks; NOT syncing to NAS."
    sacct -j "${ARRAY_JOB_ID}" --format=JobID,JobName,State,ExitCode
    notify_sync_status \
        "ECODA: annotation NOT synced (job ${ARRAY_JOB_ID})" \
        "Annotation sync to NAS failed for job ${ARRAY_JOB_ID} (datasets: ${DATASET_NAMES[*]}): non-COMPLETED tasks.
$(sacct -j "${ARRAY_JOB_ID}" --format=JobID,JobName,State,ExitCode -n 2>/dev/null || true)"
    exit 1
fi

echo "All tasks completed successfully. Starting local sync to NAS from login node..."
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
        notify_sync_status \
            "ECODA: annotation NOT synced (job ${ARRAY_JOB_ID})" \
            "Annotation sync to NAS failed for job ${ARRAY_JOB_ID}: no dataset output dirs found under ${HPC_SCRATCH_DIR}."
        exit 1
    fi
    echo "Results synchronized to ${NAS_TARGET_DIR}/<DS_NAME>/output/ (${SYNCED_COUNT} datasets)"
    echo "Success: Data safely synchronized to the NAS."
    notify_sync_status \
        "ECODA: annotation synced to NAS (job ${ARRAY_JOB_ID})" \
        "Annotation results for job ${ARRAY_JOB_ID} synced to NAS (${SYNCED_COUNT} datasets: ${DATASET_NAMES[*]})."
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable even from this login node."
    notify_sync_status \
        "ECODA: annotation NOT synced (job ${ARRAY_JOB_ID})" \
        "Annotation sync to NAS failed for job ${ARRAY_JOB_ID}: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
    exit 1
fi
