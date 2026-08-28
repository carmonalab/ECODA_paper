#!/bin/bash
# Compute-node watchdog for one Stage 5 method/analysis matrix array.
set -euo pipefail
SCRIPT_DIR=""
if [[ -n "${PROJECT_ROOT:-}" &&
      -f "${PROJECT_ROOT}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${PROJECT_ROOT}/src/5_run_benchmark_methods"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        -f "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}/src/5_run_benchmark_methods"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  if [[ -n "${SLURM_JOB_ID:-}" ]] &&
     command -v scontrol >/dev/null 2>&1; then
    submitted_command="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null |
      sed -n 's/.* Command=\([^ ]*\).*/\1/p' | head -1 || true)"
    submitted_dir="$(dirname "${submitted_command}")"
    if [[ -n "${submitted_command}" &&
          -f "${submitted_dir}/../slurm_config.sh" ]]; then
      SCRIPT_DIR="$(cd "${submitted_dir}" && pwd)"
    fi
  fi
fi
if [[ -z "${SCRIPT_DIR}" || ! -f "${SCRIPT_DIR}/../slurm_config.sh" ]]; then
  echo "ERROR: could not recover the repository source directory." >&2
  exit 1
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"
[[ $# -ge 9 ]] || { echo "Usage: matrix_watchdog.sh RUN_ROOT LABEL MANIFEST ARRAY_ID MEM MAX_MEM PARTITION THROTTLE WORKER [flags...]" >&2; exit 2; }
RUN_ROOT="$1"; LABEL="$2"; ROOT_MANIFEST="$3"; ARRAY_ID="$4"; CURRENT_MEMORY="$5"; MAX_MEMORY="$6"; PARTITION="$7"; THROTTLE="$8"; WORKER_SCRIPT="$9"; shift 9
WORKER_FLAGS=("$@")
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: matrix run root is missing: ${RUN_ROOT}" >&2; exit 1; }
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  { echo "ERROR: matrix manifest is outside the run root." >&2; exit 1; }
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 ||
  { echo "ERROR: matrix manifest is invalid." >&2; exit 1; }
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  unset BENCHMARK_MANIFEST
  # ROOT_MANIFEST is method-scoped (column 3 is the method); the immutable
  # selection's view/scope pair is validated by the submitter and root gates.
  if awk -F '\t' -v expected="batch_effect_${ANALYSIS_PASS}" \
      '$2 != expected {exit 1}' "${ROOT_MANIFEST}"; then
    :
  else
    echo "ERROR: pass matrix manifest contains a wrong view." >&2
    exit 1
  fi
fi
STATUS_DIR="${RUN_ROOT}/status/watchdogs"
STATUS_FILE="${STATUS_DIR}/${LABEL}.status"
CURRENT_MANIFEST="${ROOT_MANIFEST}"
RETRY_INDEX=0
SCHEDULER_IDS=("${ARRAY_ID}")
safe_label="$(printf '%s' "${LABEL}" | tr '/:,\t ' '_____')"
STATUS_FILE="${STATUS_DIR}/${safe_label}.status"
LOGS_DIR="${RUN_ROOT}/logs"
mkdir -p "${STATUS_DIR}" "${LOGS_DIR}"
status_write() {
  local state="$1" reason="${2:-}" tmp="${STATUS_FILE}.tmp.$$"
  {
    printf 'STATE=%s\nLABEL=%s\nREASON=%s\nARRAY_JOB_ID=%s\nRETRY_INDEX=%s\n' "${state}" "${LABEL}" "${reason}" "${ARRAY_ID}" "${RETRY_INDEX}"
    local scheduler_id
    for scheduler_id in "${SCHEDULER_IDS[@]}"; do
      printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
    done
  } > "${tmp}"
  mv -f "${tmp}" "${STATUS_FILE}"
}
fail() { status_write FAIL "$1"; exit 1; }
bump_mem() { [[ "$1" =~ ^([0-9]+)([GT])$ ]] || return 1; printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"; }
mem_ge() { local a="$1" b="$2" an as bn bs; [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1; an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"; [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1; bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"; [[ "${as}" == T ]] && an=$((an * 1024)); [[ "${bs}" == T ]] && bn=$((bn * 1024)); (( an >= bn )); }
classify() {
  local job="$1" expected="$2" rows jid state task
  OOM_TASKS=()
  FAILED_TASKS=()
  ecoda_wait_array_accounting "${job}" "${expected}" "${MATRIX_WATCHDOG_POLL_SECONDS:-30}" || return 1
  rows="${ECODA_ACCOUNTING_ROWS}"
  while IFS='|' read -r jid state exitcode; do
    [[ "${jid}" =~ ^${job}_[0-9]+$ ]] || continue
    task="${jid#${job}_}"
    state="${state%%+*}"
    case "${state}" in
      COMPLETED) [[ -z "${exitcode}" || "${exitcode}" == 0:0* ]] || FAILED_TASKS+=("${task}:${state}:${exitcode}") ;;
      OUT_OF_MEMORY) OOM_TASKS+=("${task}") ;;
      *) FAILED_TASKS+=("${task}:${state}:${exitcode}") ;;
    esac
  done <<< "${rows}"
}
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "benchmark matrix manifest empty"
while :; do
  classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal matrix task rows"
  [[ ${#FAILED_TASKS[@]} -eq 0 ]] || fail "non-OOM matrix task failure: ${FAILED_TASKS[*]}"
  [[ ${#OOM_TASKS[@]} -eq 0 ]] && break
  mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}" && fail "matrix OOM at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  NEXT_MEMORY="$(bump_mem "${CURRENT_MEMORY}")" || fail "unparseable matrix memory"; mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}" && NEXT_MEMORY="${MAX_MEMORY}"
  RETRY_INDEX=$((RETRY_INDEX + 1)); [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded matrix OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/${safe_label}.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    sed -n "${task}p" "${CURRENT_MANIFEST}" >> "${RETRY_TMP}"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 3; then
    fail "failed to install matrix retry manifest atomically"
  fi
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" ||
    fail "matrix retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 3 || fail "matrix retry manifest is invalid"
  retry_count="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  retry_export="ALL,ANALYSIS_MANIFEST=${RETRY_MANIFEST},MATRIX_RETRY=1,JOB_LOG_PREFIX=${LOGS_DIR}/5_matrix_${safe_label}_retry${RETRY_INDEX}"
  if [[ -n "${ANALYSIS_PASS:-}" ]]; then
    unset BENCHMARK_MANIFEST
    retry_export="${retry_export},ANALYSIS_PASS=${ANALYSIS_PASS}"
  else
    retry_export="${retry_export},BENCHMARK_MANIFEST=${RETRY_MANIFEST}"
  fi
  set +e
  retry_msg="$(sbatch --parsable --array="1-${retry_count}%${THROTTLE}" --partition="${PARTITION}" "${WORKER_FLAGS[@]}" --mem="${NEXT_MEMORY}" \
    --output="${LOGS_DIR}/5_matrix_${safe_label}_retry${RETRY_INDEX}_%A_%a.log" --error="${LOGS_DIR}/5_matrix_${safe_label}_retry${RETRY_INDEX}_%A_%a.err" \
    --mail-user="${USER_EMAIL}" --export="${retry_export}" "${WORKER_SCRIPT}")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected matrix OOM retry"
  ARRAY_ID="${retry_msg%%;*}"
  [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || fail "invalid matrix retry array id"
  SCHEDULER_IDS+=("${ARRAY_ID}")
  if [[ -n "${ANALYSIS_PASS:-}" ]]; then
    echo "BATCH_EFFECT_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  else
    echo "MATRIX_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  fi
  CURRENT_MANIFEST="${RETRY_MANIFEST}"; expected="${retry_count}"; CURRENT_MEMORY="${NEXT_MEMORY}"
done
status_write OK "matrix array and OOM retries completed"
