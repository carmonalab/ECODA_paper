#!/bin/bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/sync_status_email.sh"
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
#   ./1_submit_hpc_array.sh --ds_name Alzheimer --view batch_effect_uncorrected --mem 256G
#   ./1_submit_hpc_array.sh --partition shared-cpu # override scheduler partition
#                                                    # (repeat the original --ds_name flag)
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
# The completion/sync emails include a per-task report (task -> dataset,
# state, elapsed, exit code) + array wall time, built from the FULL sorted key
# list (DS_BY_TASK below) so single-dataset and _debug runs map correctly.
# ---------------------------------------------------------------------------
echo "=== Submitting preprocessing array job ==="

DS_NAME_ARG=""
FORCE_ARG=0
VIEW_ARG=""
MEMORY="128G"
PARTITION_ARG=""
SYNC_ONLY_IDS=""
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
    --view)
      VIEW_ARG="${2:-}"
      if [[ -z "${VIEW_ARG}" ]]; then
        echo "ERROR: --view requires a view name."
        exit 1
      fi
      shift 2
      ;;
    --view=*)
      VIEW_ARG="${1#*=}"
      if [[ -z "${VIEW_ARG}" ]]; then
        echo "ERROR: --view requires a view name."
        exit 1
      fi
      shift
      ;;
    --mem)
      MEMORY="${2:-}"
      shift 2
      ;;
    --partition)
      PARTITION_ARG="${2:-}"
      shift 2
      ;;
    --force)
      FORCE_ARG=1
      shift
      ;;
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

if [[ -n "${SYNC_ONLY_IDS}" && ${FORCE_ARG} -eq 1 ]]; then
  echo "ERROR: --sync-only cannot be combined with --force."
  exit 1
fi

if [[ -n "${SYNC_ONLY_IDS}" && -z "${VIEW_ARG}" ]]; then
  echo "ERROR: --sync-only requires the original --view when view-scoped."
  exit 1
fi

# Passed to workers via the environment (sbatch propagates the submit script's
# environment); 1.1_run_worker.sh forwards them to 1.1.1_preprocess.py.
export FORCE_PREPROCESS="${FORCE_ARG}"
export PREPROCESS_VIEW="${VIEW_ARG}"
# ---------------------------------------------------------------------------
# Task -> dataset mapping for the email reports. Task ids map 1:1 to the FULL
# sorted jq key list (includes "_" keys, which sort last): single-dataset mode
# submits a 1-task array at the dataset's position in this full list, so
# building DS_BY_TASK from the non-underscore prefix only would mis-map _debug
# tasks. Built once; used by preprocess_report below.
# ---------------------------------------------------------------------------
DS_BY_TASK=()
while IFS= read -r ds_name; do
  DS_BY_TASK+=("${ds_name}")
done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")

task_dataset() {
  local TASK_ID="$1"
  if [[ ${TASK_ID} -ge 1 && ${TASK_ID} -le ${#DS_BY_TASK[@]} ]]; then
    printf '%s' "${DS_BY_TASK[$((TASK_ID - 1))]}"
  else
    printf '%s' "?"
  fi
}

# Renders the per-task report block for email bodies: task -> dataset, state,
# elapsed, exit code (exit code only shown for non-COMPLETED tasks) + array
# wall time. Prints an n/a line when sacct is unavailable or purged.
preprocess_report() {
  local JOB_ID="$1"
  local TASKS
  TASKS="$(array_task_report "${JOB_ID}")"
  if [[ -z "${TASKS}" ]]; then
    printf '%s' "Per-task report: n/a (sacct unavailable or job purged)."
    return 0
  fi
  printf '%s\n' "Per-task report (task -> dataset, state, elapsed, exit code):"
  local task state elapsed exitcode extra
  while IFS=$'\t' read -r task state elapsed exitcode; do
    extra=""
    [[ "${state}" == "COMPLETED" ]] || extra=" (${exitcode})"
    printf '  %s -> %-30s %-14s%s  %s\n' "${task}" "$(task_dataset "${task}")" "${state}" "${extra}" "${elapsed}"
  done <<< "${TASKS}"
  printf 'Array wall time: %s\n' "$(array_wall_time "${JOB_ID}")"
}

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
echo "Memory: ${MEMORY}"
PARTITION="${PARTITION_ARG:-${SLURM_PARTITION}}"
echo "Partition: ${PARTITION}"

if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  echo "=== Sync-only resume mode: job ${SYNC_ONLY_IDS} (no submission) ==="
  ARRAY_JOB_ID="${SYNC_ONLY_IDS}"
else
  mkdir -p "${LOGS_DIR}"
  SUBMIT_MSG=$(sbatch \
      --mem="${MEMORY}" \
      --partition="${PARTITION}" \
      --output="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.log" \
      --error="${LOGS_DIR}/3_scrnaseq_preprocessing_%A_%a.err" \
      --mail-user="${USER_EMAIL}" \
      "$(dirname "${BASH_SOURCE[0]}")/1.1_run_worker.sh")

  ARRAY_JOB_ID=$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')
  echo "Array Job ID allocated: ${ARRAY_JOB_ID}"
fi

# ---------------------------------------------------------------------------
# Monitor & sync results back to NAS
# ---------------------------------------------------------------------------
echo "=== Monitoring job completion ==="
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
        "ECODA: preprocess NOT synced (job ${ARRAY_JOB_ID})" \
        "Preprocess sync to NAS failed for job ${ARRAY_JOB_ID} (datasets: ${DATASET_NAMES[*]}): sacct returned no states (job purged or unknown id)."
    exit 1
fi
# Fail-closed: every row (array master + tasks + batch steps) must be COMPLETED.
if grep -qvE '^ *COMPLETED *$' <<< "${STATES}"; then
    echo "ERROR: Array Job ${ARRAY_JOB_ID} had non-COMPLETED tasks; NOT syncing to NAS."
    sacct -j "${ARRAY_JOB_ID}" --format=JobID,JobName,State,ExitCode
    notify_sync_status \
        "ECODA: preprocess NOT synced (job ${ARRAY_JOB_ID})" \
        "Preprocess sync to NAS failed for job ${ARRAY_JOB_ID} (datasets: ${DATASET_NAMES[*]}): non-COMPLETED tasks.
$(preprocess_report "${ARRAY_JOB_ID}")"
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
        notify_sync_status \
            "ECODA: preprocess NOT synced (job ${ARRAY_JOB_ID})" \
            "Preprocess sync to NAS failed for job ${ARRAY_JOB_ID}: no dataset output dirs found under ${HPC_SCRATCH_DIR}."
        exit 1
    fi
    echo "Results synchronized to ${NAS_TARGET_DIR}/<DS_NAME>/output/ (${SYNCED_COUNT} datasets)"
    notify_sync_status \
        "ECODA: preprocess synced to NAS (job ${ARRAY_JOB_ID})" \
        "Preprocess results for job ${ARRAY_JOB_ID} synced to NAS (${SYNCED_COUNT} datasets: ${DATASET_NAMES[*]}).
$(preprocess_report "${ARRAY_JOB_ID}")"
else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
    notify_sync_status \
        "ECODA: preprocess NOT synced (job ${ARRAY_JOB_ID})" \
        "Preprocess sync to NAS failed for job ${ARRAY_JOB_ID}: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
    exit 1
fi
