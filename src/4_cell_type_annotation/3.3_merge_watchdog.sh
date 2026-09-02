#!/bin/bash
# Compute-node watchdog for dataset-parallel annotation merges.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
[[ $# -eq 7 ]] || { echo "Usage: 3.3_merge_watchdog.sh RUN_ID MANIFEST ARRAY_ID MEM MAX_MEM PARTITION THROTTLE" >&2; exit 2; }
RUN_ID="$1"; ROOT_MANIFEST="$2"; ARRAY_ID="$3"; CURRENT_MEMORY="$4"; MAX_MEMORY="$5"; PARTITION="$6"; THROTTLE="$7"
ANNOTATION_WORKER_TIME_LIMIT="${ANNOTATION_WORKER_TIME_LIMIT:-12:00:00}"
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
STATUS_FILE="${RUN_ROOT}/status/merge_watchdog"
CURRENT_MANIFEST="${ROOT_MANIFEST}"
RETRY_INDEX=0
SCHEDULER_IDS=("${ARRAY_ID}")
RUNTIME_EXPORT=""

status_write() {
  local state="$1" reason="${2:-}" tmp="${STATUS_FILE}.tmp.$$"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\nRUN_ID=%s\nREASON=%s\nARRAY_JOB_ID=%s\nRETRY_INDEX=%s\n' "${state}" "${RUN_ID}" "${reason}" "${ARRAY_ID}" "${RETRY_INDEX}"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then printf 'SCHEDULER_ID=%s\n' "${SLURM_JOB_ID}"; fi
    local scheduler_id
    for scheduler_id in "${SCHEDULER_IDS[@]}"; do printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"; done
  } > "${tmp}"
  mv -f "${tmp}" "${STATUS_FILE}"
}
fail() { echo "ERROR: $1" >&2; status_write FAIL "$1"; exit 1; }
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE}" || \
  fail "Stage 4 immutable runtime validation failed before merge retry handling"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage4 0)" || \
  fail "Stage 4 runtime export construction failed"
bump_mem() { [[ "$1" =~ ^([0-9]+)([GT])$ ]] || return 1; printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"; }
mem_ge() {
  local a="$1" b="$2" an as bn bs
  [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1; an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"
  [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1; bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"
  [[ "${as}" == T ]] && an=$((an * 1024)); [[ "${bs}" == T ]] && bn=$((bn * 1024)); (( an >= bn ))
}
classify() {
  local job="$1" expected="$2" rows jid state task
  OOM_TASKS=()
  FAILED_TASKS=()
  ecoda_wait_array_accounting "${job}" "${expected}" "${ANNOTATION_WATCHDOG_POLL_SECONDS:-30}" || return 1
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
validate_markers() {
  local ds views run_root marker view name path expected_sources expected_records
  local union_path union_md5 union_size marker_field
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && "${run_root}" == "${RUN_ROOT}" ]] || return 1
    marker="${RUN_ROOT}/datasets/${ds}/merge.ok"
    [[ -s "${marker}" ]] || { echo "ERROR: merge marker missing: ${marker}" >&2; return 1; }
    marker_field() { sed -n "s/^${1}=//p" "${marker}" | head -1; }
    [[ "$(marker_field STATE)" == "OK" && "$(marker_field DATASET)" == "${ds}" &&
       "$(marker_field VIEWS)" == "${views}" ]] || return 1
    union_path="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    [[ "$(marker_field UNION_PATH)" == "${union_path}" ]] || return 1
    ecoda_validate_checksum "${union_path}" || return 1
    union_md5="$(ecoda_md5_file "${union_path}")"
    union_size="$(wc -c < "${union_path}" | tr -d '[:space:]')"
    [[ "$(marker_field UNION_MD5)" == "${union_md5}" &&
       "$(marker_field UNION_SIZE)" == "${union_size}" ]] || return 1
    expected_sources=""
    expected_records=""
    IFS=',' read -r -a view_list <<< "${views}"
    for view in "${view_list[@]}"; do
      name="$(jq -r --arg ds "${ds}" --arg view "${view}" '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' "${DATASETS_JSON_FILE}")"
      [[ -n "${name}" ]] || return 1
      path="${HPC_SCRATCH_DIR}/${ds}/output/${name}"
      [[ -s "${path}" ]] || { echo "ERROR: merged h5ad missing: ${path}" >&2; return 1; }
      ecoda_validate_checksum "${path}" || return 1
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
        --h5ad "${path}" --require-sidecar >/dev/null 2>&1 || return 1
      record="${path}|$(ecoda_md5_file "${path}")|$(wc -c < "${path}" | tr -d '[:space:]')"
      [[ -z "${expected_sources}" ]] && expected_sources="${path}" ||
        expected_sources="${expected_sources};${path}"
      [[ -z "${expected_records}" ]] && expected_records="${record}" ||
        expected_records="${expected_records};${record}"
    done
    [[ "$(marker_field SOURCE_H5ADS)" == "${expected_sources}" &&
       "$(marker_field SOURCE_RECORDS)" == "${expected_records}" ]] || return 1
  done < "${ROOT_MANIFEST}"
}
ecoda_validate_run_id "${RUN_ID}" || exit 1
[[ -d "${RUN_ROOT}" ]] || fail "Stage 4 run root is missing"
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  fail "merge manifest is outside run root"
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 || fail "merge manifest is invalid"
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "merge manifest empty"
while :; do
  classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal merge rows"
  [[ ${#FAILED_TASKS[@]} -eq 0 ]] || fail "non-OOM merge failure: ${FAILED_TASKS[*]}"
  [[ ${#OOM_TASKS[@]} -eq 0 ]] && break
  mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}" && fail "merge OOM at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  NEXT_MEMORY="$(bump_mem "${CURRENT_MEMORY}")" || fail "unparseable merge memory"
  mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}" && NEXT_MEMORY="${MAX_MEMORY}"
  RETRY_INDEX=$((RETRY_INDEX + 1)); [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded merge OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/merge.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    sed -n "${task}p" "${CURRENT_MANIFEST}" >> "${RETRY_TMP}"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 3; then
    fail "failed to install merge retry manifest atomically"
  fi
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" ||
    fail "merge retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 3 || fail "merge retry manifest is invalid"
  retry_count="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  set +e
  retry_msg="$(sbatch --parsable --array="1-${retry_count}%${THROTTLE}" --partition="${PARTITION}" --time="${ANNOTATION_WORKER_TIME_LIMIT}" --mem="${NEXT_MEMORY}" \
    --output="${LOGS_DIR}/4_annotation_merge_retry${RETRY_INDEX}_%A_%a.log" --error="${LOGS_DIR}/4_annotation_merge_retry${RETRY_INDEX}_%A_%a.err" \
    --mail-user="${USER_EMAIL}" --export="ALL,ANNOTATION_MERGE_MANIFEST=${RETRY_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},FORCE_ANNOTATION=1,${RUNTIME_EXPORT}" \
    "${SCRIPT_DIR}/3.2_merge_worker.sh")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected merge OOM retry"
  ARRAY_ID="${retry_msg%%;*}"
  [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || fail "invalid merge retry array id"
  SCHEDULER_IDS+=("${ARRAY_ID}")
  echo "ANNOTATION_MERGE_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"; expected="${retry_count}"; CURRENT_MEMORY="${NEXT_MEMORY}"
done
validate_markers "${ROOT_MANIFEST}" || fail "merged annotation artifacts failed validation"
status_write OK "all selected dataset merges completed and validated"
