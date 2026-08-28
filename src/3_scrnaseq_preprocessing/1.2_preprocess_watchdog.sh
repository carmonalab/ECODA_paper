#!/bin/bash
# Compute-node watchdog for the manifest-driven Stage 3 array.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

if [[ $# -eq 7 ]]; then
  RUN_ID="$1"; ROOT_MANIFEST="$2"; ARRAY_ID="$3"; CURRENT_MEMORY="$4"; MAX_MEMORY="$5"; PARTITION="$6"; THROTTLE="$7"
elif [[ $# -eq 6 ]]; then
  # Compatibility for the former batch-effect submitter. New callers always
  # pass the explicit run id and root selection manifest.
  ARRAY_ID="$1"; ROOT_MANIFEST="$2"; CURRENT_MEMORY="$3"; MAX_MEMORY="$4"; PARTITION="$5"; THROTTLE="$6"
  RUN_ID="${PREPROCESS_RUN_ID:-legacy_${ARRAY_ID}}"
else
  echo "Usage: 1.2_preprocess_watchdog.sh RUN_ID MANIFEST ARRAY_ID MEM MAX_MEM PARTITION THROTTLE" >&2
  exit 2
fi
ecoda_validate_run_id "${RUN_ID}" || exit 1
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
if [[ -n "${PREPROCESS_RUN_ROOT:-}" ]]; then
  expected_root="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
  supplied_root="$(ecoda_realpath_existing "${PREPROCESS_RUN_ROOT}" 2>/dev/null || true)"
  [[ -n "${expected_root}" && "${expected_root}" == "${supplied_root}" ]] ||
    { echo "ERROR: PREPROCESS_RUN_ROOT is not the canonical Stage 3 run root." >&2; exit 1; }
fi
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: Stage 3 run root is missing: ${RUN_ROOT}" >&2; exit 1; }
STATUS_FILE="${RUN_ROOT}/status/watchdog"
CURRENT_MANIFEST="${PREPROCESS_PENDING_MANIFEST:-${ROOT_MANIFEST}}"
RETRY_INDEX=0
SCHEDULER_IDS=("${ARRAY_ID}")
atomic_status() {
  local state="$1" reason="${2:-}" tmp="${STATUS_FILE}.tmp.$$"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\nRUN_ID=%s\nREASON=%s\nRETRY_INDEX=%s\n' "${state}" "${RUN_ID}" "${reason}" "${RETRY_INDEX}"
    printf 'ARRAY_JOB_ID=%s\n' "${ARRAY_ID}"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      printf 'SCHEDULER_ID=%s\n' "${SLURM_JOB_ID}"
    fi
    if [[ ${#SCHEDULER_IDS[@]} -gt 0 ]]; then
      local scheduler_id
      for scheduler_id in "${SCHEDULER_IDS[@]}"; do
        printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
      done
    fi
  } > "${tmp}" || return 1
  mv -f "${tmp}" "${STATUS_FILE}"
}

set_owner_state() {
  local state="$1" reason="$2" owners_file="${RUN_ROOT}/manifests/owners.tsv"
  local row owner expected_owner owner_key rc=0 count=0
  [[ -r "${owners_file}" ]] || return 1
  while IFS=$'\t' read -r row owner; do
    [[ -n "${row}" && -n "${owner}" ]] || { rc=1; continue; }
    owner_key="${row}"
    [[ "${owner_key}" == */* ]] || owner_key="${owner_key}/batch_effect_uncorrected"
    expected_owner="$(ecoda_owner_dir stage3 "${owner_key}")"
    [[ "${owner}" == "${expected_owner}" ]] || { rc=1; continue; }
    count=$((count + 1))
    if ! ecoda_owner_set_state "${owner}" "${state}" "${reason}"; then
      rc=1
    fi
  done < "${owners_file}"
  [[ ${count} -gt 0 ]] || rc=1
  return "${rc}"
}

fail() {
  local reason="$1"
  local owner_rc=0
  set_owner_state FAIL "${reason}" || owner_rc=1
  atomic_status FAIL "${reason}" || owner_rc=1
  exit 1
}

bump_memory() {
  local value="$1"
  [[ "${value}" =~ ^([0-9]+)([GT])$ ]] || return 1
  printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"
}

mem_ge() {
  local a="$1" b="$2" an as bn bs
  [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1
  an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"
  [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1
  bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"
  [[ "${as}" == "T" ]] && an=$((an * 1024))
  [[ "${bs}" == "T" ]] && bn=$((bn * 1024))
  (( an >= bn ))
}

wait_and_classify() {
  local job="$1" expected="$2" rows jid state exitcode task
  OOM_TASKS=()
  FAILED_TASKS=()
  ecoda_wait_array_accounting "${job}" "${expected}" "${STAGE3_WATCHDOG_POLL_SECONDS:-30}" || return 1
  rows="${ECODA_ACCOUNTING_ROWS}"
  while IFS='|' read -r jid state exitcode; do
    [[ "${jid}" =~ ^${job}_[0-9]+$ ]] || continue
    task="${jid#${job}_}"
    state="${state%%+*}"
    case "${state}" in
      COMPLETED) [[ -z "${exitcode}" || "${exitcode}" == "0:0"* ]] || FAILED_TASKS+=("${task}:${state}:${exitcode}") ;;
      OUT_OF_MEMORY) OOM_TASKS+=("${task}") ;;
      *) FAILED_TASKS+=("${task}:${state}:${exitcode}") ;;
    esac
  done <<< "${rows}"
}

validate_manifest_outputs() {
  local manifest="$1" ds view output path
  while IFS=$'\t' read -r ds view; do
    [[ -n "${ds}" ]] || continue
    output="$(ecoda_view_output_name "${ds}" "${view}")"
    [[ -n "${output}" ]] || return 1
    path="${HPC_SCRATCH_DIR}/${ds}/output/${output}"
    [[ -s "${path}" ]] || return 1
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
      --path "${path}" --view "${view}" --method "Stage 3 watchdog" >/dev/null 2>&1 || return 1
    ecoda_validate_checksum "${path}" || return 1
  done < "${manifest}"
}

ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  fail "Stage 3 selection manifest is outside run root"
ecoda_validate_run_owned_path "${CURRENT_MANIFEST}" "${RUN_ROOT}" ||
  fail "Stage 3 array manifest is outside run root"
ecoda_validate_manifest "${ROOT_MANIFEST}" 2 ||
  fail "Stage 3 selection manifest is invalid"
ecoda_validate_manifest "${CURRENT_MANIFEST}" 2 ||
  fail "Stage 3 array manifest is invalid"
ecoda_validate_run_owned_path "${RUN_ROOT}/manifests/owners.tsv" "${RUN_ROOT}" ||
  fail "Stage 3 owner manifest is missing or outside run root"
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "Stage 3 array manifest is empty"

while :; do
  wait_and_classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal Stage 3 task rows"
  if [[ ${#FAILED_TASKS[@]} -gt 0 ]]; then
    fail "non-OOM Stage 3 task failure: ${FAILED_TASKS[*]}"
  fi
  if [[ ${#OOM_TASKS[@]} -eq 0 ]]; then break; fi
  if mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}"; then
    fail "OUT_OF_MEMORY Stage 3 tasks at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  fi
  NEXT_MEMORY="$(bump_memory "${CURRENT_MEMORY}")" || fail "unparseable Stage 3 memory '${CURRENT_MEMORY}'"
  if mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}"; then NEXT_MEMORY="${MAX_MEMORY}"; fi
  RETRY_INDEX=$((RETRY_INDEX + 1))
  [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded Stage 3 OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/selection.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    line="$(sed -n "${task}p" "${CURRENT_MANIFEST}")"
    [[ -n "${line}" ]] || fail "OOM task ${task} is absent from Stage 3 manifest"
    printf '%s\n' "${line}" >> "${RETRY_TMP}"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 2; then
    fail "failed to install Stage 3 retry manifest atomically"
  fi
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" ||
    fail "Stage 3 retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 2 ||
    fail "Stage 3 retry manifest is invalid"
  RETRY_COUNT="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  set +e
  RETRY_MSG="$(sbatch --parsable --array="1-${RETRY_COUNT}%${THROTTLE}" \
    --mem="${NEXT_MEMORY}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX}_%A_%a.log" \
    --error="${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX}_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    --export="ALL,PREPROCESS_SELECTION_FILE=${RETRY_MANIFEST},PREPROCESS_RUN_ROOT=${RUN_ROOT},FORCE_PREPROCESS=1,PREPROCESS_ERROR_PREFIX=${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX}" \
    "${SCRIPT_DIR}/1.1_run_worker.sh")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected Stage 3 OOM retry"
  ARRAY_ID="${RETRY_MSG%%;*}"
  [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || fail "invalid Stage 3 retry array id"
  SCHEDULER_IDS+=("${ARRAY_ID}")
  echo "PREPROCESS_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"
  expected="${RETRY_COUNT}"
  CURRENT_MEMORY="${NEXT_MEMORY}"
done

validate_manifest_outputs "${ROOT_MANIFEST}" || fail "Stage 3 h5ad schema/checksum validation failed"
set_owner_state OK "Stage 3 preprocessing artifacts validated" ||
  fail "failed to finalize Stage 3 owners"
if ! atomic_status OK "all selected Stage 3 tasks completed and artifacts validated"; then
  fail "failed to write Stage 3 watchdog success status"
fi
printf 'Stage 3 watchdog completed for run %s\n' "${RUN_ID}"
