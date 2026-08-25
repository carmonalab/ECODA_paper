#!/bin/bash
# Compute-node watchdog for the stage-wise onboarding preprocessing array.
# It owns terminal accounting, OOM-only retries, artifact validation, and NAS
# synchronization so an SSH/login-session loss cannot interrupt escalation.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Slurm may execute this script from /var/spool/slurmd; recover the repository
# path before sourcing the canonical configuration.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

if [[ $# -ne 6 ]]; then
  echo "Usage: 1.2_preprocess_watchdog.sh ARRAY_ID MANIFEST INITIAL_MEM MAX_MEM PARTITION THROTTLE" >&2
  exit 2
fi

ARRAY_ID="$1"
MANIFEST="$2"
ROOT_MANIFEST="${MANIFEST}"
CURRENT_MANIFEST="${MANIFEST}"
CURRENT_MEMORY="$3"
MAX_MEMORY="$4"
PARTITION="$5"
THROTTLE="$6"
VIEW="${PREPROCESS_VIEW:-batch_effect_uncorrected}"
FORCE_PREPROCESS="${FORCE_PREPROCESS:-0}"
WORKER_SCRIPT="${SCRIPT_DIR}/1.1_run_worker.sh"
STATUS_DIR="${HPC_SCRATCH_DIR}/_preprocessing_watchdog"
STATUS_FILE="${STATUS_DIR}/${SLURM_JOB_ID}.status"
mkdir -p "${STATUS_DIR}"

ARRAY_JOB_IDS="${ARRAY_ID}"
OOM_TASKS=()
FAILED_ROWS=()
TASK_COUNT=0

write_status() {
  local STATE="$1"
  local FAIL_REASON="$2"
  local TMP_FILE="${STATUS_FILE}.tmp.$$"
  {
    printf 'STATE=%s\n' "${STATE}"
    printf 'WATCHDOG_JOB_ID=%s\n' "${SLURM_JOB_ID}"
    printf 'ARRAY_JOB_IDS=%s\n' "${ARRAY_JOB_IDS}"
    printf 'VIEW=%s\n' "${VIEW}"
    printf 'MANIFEST=%s\n' "${ROOT_MANIFEST}"
    printf 'DATASETS=%s\n' "$(tr '\n' ' ' < "${ROOT_MANIFEST}" | sed 's/[[:space:]]*$//')"
    printf 'FAIL_REASON=%s\n' "${FAIL_REASON}"
    if [[ ${#FAILED_ROWS[@]} -gt 0 ]]; then
      printf 'FAILED_ROWS=%s\n' "${FAILED_ROWS[*]}"
    fi
  } > "${TMP_FILE}"
  mv -f "${TMP_FILE}" "${STATUS_FILE}"
}

scheduler_rows() {
  sacct -n -P -X -j "$1" \
    --format=JobID,State,ExitCode,MaxRSS 2>/dev/null || true
}

wait_for_terminal() {
  local JOB_ID="$1"
  local ATTEMPT=0
  local ROWS=""
  # The watchdog itself is the durable wait owner. Bound the wait by the
  # watchdog's Slurm time limit rather than leaving a login-node poll alive.
  while (( ATTEMPT < 2880 )); do
    ROWS="$(scheduler_rows "${JOB_ID}")"
    if [[ -n "${ROWS//[[:space:]]/}" ]] \
        && ! [[ "${ROWS}" =~ (PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING) ]]; then
      printf '%s\n' "${ROWS}"
      return 0
    fi
    sleep 15
    ATTEMPT=$((ATTEMPT + 1))
  done
  echo "ERROR: scheduler accounting did not reach terminal state for ${JOB_ID}." >&2
  return 1
}

analyze_rows() {
  local ROWS="$1"
  OOM_TASKS=()
  FAILED_ROWS=()
  TASK_COUNT=0
  local JOB_ID STATE EXIT_CODE MAX_RSS TASK_ID
  while IFS='|' read -r JOB_ID STATE EXIT_CODE MAX_RSS; do
    [[ -n "${JOB_ID}" ]] || continue
    [[ "${JOB_ID}" == "${ARRAY_ID}_"* ]] || continue
    TASK_ID="${JOB_ID#${ARRAY_ID}_}"
    if ! [[ "${TASK_ID}" =~ ^[0-9]+$ ]]; then
      FAILED_ROWS+=("${JOB_ID}|${STATE}|${EXIT_CODE}|${MAX_RSS}")
      continue
    fi
    TASK_COUNT=$((TASK_COUNT + 1))
    if [[ "${STATE}" == "COMPLETED" && "${EXIT_CODE}" == "0:0" ]]; then
      continue
    fi
    if [[ "${STATE}" == OUT_OF_MEMORY* ]]; then
      OOM_TASKS+=("${TASK_ID}")
    else
      FAILED_ROWS+=("${JOB_ID}|${STATE}|${EXIT_CODE}|${MAX_RSS}")
    fi
  done <<< "${ROWS}"
}

memory_gb() {
  case "$1" in
    *G) printf '%s' "${1%G}" ;;
    *T) printf '%s' "$(( ${1%T} * 1024 ))" ;;
    *) return 1 ;;
  esac
}

next_memory() {
  local CURRENT_GB MAX_GB NEXT_GB
  CURRENT_GB="$(memory_gb "${CURRENT_MEMORY}")" || return 1
  MAX_GB="$(memory_gb "${MAX_MEMORY}")" || return 1
  if (( CURRENT_GB >= MAX_GB )); then
    return 1
  fi
  NEXT_GB=$((CURRENT_GB * 2))
  (( NEXT_GB > MAX_GB )) && NEXT_GB="${MAX_GB}"
  printf '%sG' "${NEXT_GB}"
}

submit_retry_array() {
  local RETRY_MANIFEST="$1"
  local RETRY_COUNT
  RETRY_COUNT="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  local EXPORTS="ALL,PREPROCESS_DATASETS_FILE=${RETRY_MANIFEST},PREPROCESS_VIEW=${VIEW},FORCE_PREPROCESS=${FORCE_PREPROCESS}"
  local RETRY_ID
  RETRY_ID="$(sbatch --parsable \
    --array="1-${RETRY_COUNT}%${THROTTLE}" \
    --mem="${CURRENT_MEMORY}" \
    --partition="${PARTITION}" \
    --output="${LOGS_DIR}/3_scrnaseq_batch_effect_retry_%A_%a.log" \
    --error="${LOGS_DIR}/3_scrnaseq_batch_effect_retry_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    --export="${EXPORTS}" \
    "${WORKER_SCRIPT}")"
  RETRY_ID="${RETRY_ID%%;*}"
  if ! [[ "${RETRY_ID}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid retry array id: ${RETRY_ID}" >&2
    return 1
  fi
  printf '%s' "${RETRY_ID}"
}

validate_and_sync() {
  local DS_NAME OUTPUT_NAME OUTPUT_FILE DEST_DIR
  while IFS= read -r DS_NAME; do
    [[ -n "${DS_NAME}" ]] || continue
    OUTPUT_NAME="$(jq -er --arg ds "${DS_NAME}" --arg view "${VIEW}" \
      '.[$ds].views[$view].output_file_name' "${DATASETS_JSON_FILE}")"
    OUTPUT_FILE="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${OUTPUT_NAME}"
    if [[ ! -s "${OUTPUT_FILE}" ]]; then
      echo "ERROR: expected preprocessing artifact is missing or empty: ${OUTPUT_FILE}" >&2
      return 1
    fi
  done < "${ROOT_MANIFEST}"

  if [[ ! -d "${NAS_TARGET_DIR}" ]]; then
    echo "ERROR: NAS path is unreachable: ${NAS_TARGET_DIR}" >&2
    return 1
  fi
  while IFS= read -r DS_NAME; do
    [[ -n "${DS_NAME}" ]] || continue
    DEST_DIR="${NAS_TARGET_DIR}/${DS_NAME}/output"
    mkdir -p "${DEST_DIR}"
    rsync -rlptD "${HPC_SCRATCH_DIR}/${DS_NAME}/output/" "${DEST_DIR}/"
  done < "${ROOT_MANIFEST}"
}

if ! command -v jq >/dev/null 2>&1; then
  write_status FAIL "jq is unavailable on the watchdog node"
  exit 1
fi
if [[ ! -r "${MANIFEST}" ]]; then
  write_status FAIL "dataset manifest is unreadable: ${MANIFEST}"
  exit 1
fi

while :; do
  echo "Waiting for preprocessing array ${ARRAY_ID} to reach terminal accounting state."
  if ! ROWS="$(wait_for_terminal "${ARRAY_ID}")"; then
    write_status FAIL "scheduler accounting did not reach terminal state for ${ARRAY_ID}"
    exit 1
  fi
  analyze_rows "${ROWS}"
  if (( TASK_COUNT == 0 )); then
    write_status FAIL "no array task rows found for ${ARRAY_ID}"
    exit 1
  fi
  if [[ ${#FAILED_ROWS[@]} -gt 0 ]]; then
    echo "ERROR: preprocessing array ${ARRAY_ID} has non-OOM failures: ${FAILED_ROWS[*]}" >&2
    write_status FAIL "non-OOM array task failure"
    exit 1
  fi
  if [[ ${#OOM_TASKS[@]} -eq 0 ]]; then
    break
  fi

  NEXT_MEMORY="$(next_memory || true)"
  if [[ -z "${NEXT_MEMORY}" ]]; then
    echo "ERROR: OOM retry ceiling reached for ${ARRAY_ID}." >&2
    write_status FAIL "OOM retry ceiling reached at ${CURRENT_MEMORY}"
    exit 1
  fi
  RETRY_MANIFEST="${CURRENT_MANIFEST}.retry_${CURRENT_MEMORY}_to_${NEXT_MEMORY}"
  : > "${RETRY_MANIFEST}"
  for TASK_ID in "${OOM_TASKS[@]}"; do
    DS_NAME="$(sed -n "${TASK_ID}p" "${CURRENT_MANIFEST}")"
    if [[ -z "${DS_NAME}" ]]; then
      write_status FAIL "OOM task ${TASK_ID} has no dataset manifest row"
      exit 1
    fi
    printf '%s\n' "${DS_NAME}" >> "${RETRY_MANIFEST}"
  done
  CURRENT_MEMORY="${NEXT_MEMORY}"
  RETRY_ARRAY_ID="$(submit_retry_array "${RETRY_MANIFEST}")"
  ARRAY_JOB_IDS="${ARRAY_JOB_IDS},${RETRY_ARRAY_ID}"
  printf 'PREPROCESS_ARRAY_JOB_ID=%s\n' "${RETRY_ARRAY_ID}" >&2
  echo "Submitted OOM-only preprocessing retry ${RETRY_ARRAY_ID} at ${CURRENT_MEMORY}."
  ARRAY_ID="${RETRY_ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"
done

if ! validate_and_sync; then
  write_status FAIL "preprocessing artifacts or NAS synchronization failed"
  exit 1
fi
write_status OK "all preprocessing tasks completed and artifacts synchronized"
echo "Preprocessing stage completed successfully: ${ARRAY_JOB_IDS}"
