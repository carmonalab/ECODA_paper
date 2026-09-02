#!/bin/bash
# Compute-node watchdog for the global annotation chunk array.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
cd "${PROJECT_ROOT}"
[[ $# -eq 7 ]] || { echo "Usage: 1.2_annotation_watchdog.sh RUN_ID MANIFEST ARRAY_ID MEM MAX_MEM PARTITION THROTTLE" >&2; exit 2; }
RUN_ID="$1"; ROOT_MANIFEST="$2"; ARRAY_ID="$3"; CURRENT_MEMORY="$4"; MAX_MEMORY="$5"; PARTITION="$6"; THROTTLE="$7"
ANNOTATION_WORKER_TIME_LIMIT="${ANNOTATION_WORKER_TIME_LIMIT:-02:00:00}"
ecoda_validate_run_id "${RUN_ID}" || exit 1
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: Stage 4 run root is missing: ${RUN_ROOT}" >&2; exit 1; }
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  { echo "ERROR: annotation chunk manifest is outside run root." >&2; exit 1; }
STATUS_FILE="${RUN_ROOT}/status/annotation_watchdog"
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
fail() { status_write FAIL "$1"; exit 1; }
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_submission "${ECODA_RUNTIME_MODE}" || \
  fail "Stage 4 immutable runtime validation failed before annotation retry handling"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage4 0)" || \
  fail "Stage 4 runtime export construction failed"
bump_mem() { [[ "$1" =~ ^([0-9]+)([GT])$ ]] || return 1; printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"; }
mem_ge() { local a="$1" b="$2" an as bn bs; [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1; an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"; [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1; bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"; [[ "${as}" == T ]] && an=$((an * 1024)); [[ "${bs}" == T ]] && bn=$((bn * 1024)); (( an >= bn )); }
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
validate_feathers() {
  local manifest="$1" ds chunk feather_dir chunk_num feather union_path sample expected_union expected_chunk_dir
  local expected_args
  while IFS=$'\t' read -r ds chunk feather_dir; do
    [[ -n "${ds}" && -n "${chunk}" && -n "${feather_dir}" ]] || return 1
    ecoda_validate_run_owned_path "${chunk}" "${RUN_ROOT}" || return 1
    expected_chunk_dir="${RUN_ROOT}/datasets/${ds}/chunks"
    expected_union="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    [[ "${chunk}" == "${expected_chunk_dir}/chunk_"*.txt ]] || return 1
    [[ "${feather_dir}" == "${RUN_ROOT}/datasets/${ds}/annotations" ]] || return 1
    chunk_num="${chunk##*/chunk_}"; chunk_num="${chunk_num%.txt}"
    [[ "${chunk_num}" =~ ^[1-9][0-9]*$ ]] || return 1
    feather="${feather_dir}/annotations_chunk_${chunk_num}.feather"
    [[ -s "${feather}" ]] || return 1
    union_path="$(sed -n '1p' "${chunk}")"
    [[ "${union_path}" == "${expected_union}" ]] || return 1
    ecoda_validate_checksum "${union_path}" || return 1
    expected_args=(--expected-union "${union_path}")
    while IFS= read -r sample; do
      [[ -n "${sample}" ]] || return 1
      expected_args+=(--expected-sample "${sample}")
    done < <(sed -n '2,$p' "${chunk}")
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${feather}" --require-sidecar "${expected_args[@]}" >/dev/null 2>&1 || return 1
  done < "${manifest}"
}
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  fail "annotation chunk manifest is outside run root"
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 || fail "annotation chunk manifest is invalid"
manifest="${ANNOTATION_CHUNK_MANIFEST:-${ROOT_MANIFEST}}"
ecoda_validate_run_owned_path "${manifest}" "${RUN_ROOT}" ||
  fail "annotation validation manifest is outside run root"
ecoda_validate_manifest "${manifest}" 3 || fail "annotation validation manifest is invalid"
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "annotation array manifest empty"
while :; do
  classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal annotation rows"
  [[ ${#FAILED_TASKS[@]} -eq 0 ]] || fail "non-OOM annotation failure: ${FAILED_TASKS[*]}"
  [[ ${#OOM_TASKS[@]} -eq 0 ]] && break
  mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}" && fail "annotation OOM at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  NEXT_MEMORY="$(bump_mem "${CURRENT_MEMORY}")" || fail "unparseable annotation memory"; mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}" && NEXT_MEMORY="${MAX_MEMORY}"
  RETRY_INDEX=$((RETRY_INDEX + 1)); [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded annotation OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/chunks.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    sed -n "${task}p" "${CURRENT_MANIFEST}" >> "${RETRY_TMP}" ||
      fail "failed to build annotation retry manifest"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 3; then
    fail "failed to install annotation retry manifest atomically"
  fi
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" ||
    fail "annotation retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 3 || fail "annotation retry manifest is invalid"
  retry_count="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  set +e
  retry_msg="$(sbatch --parsable --array="1-${retry_count}%${THROTTLE}" --mem="${NEXT_MEMORY}" --time="${ANNOTATION_WORKER_TIME_LIMIT}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX}_%A_%a.log" --error="${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX}_%A_%a.err" --mail-user="${USER_EMAIL}" \
    --export="ALL,CHUNKS_MANIFEST=${RETRY_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},ANNOTATION_ERROR_PREFIX=${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX},${RUNTIME_EXPORT}" "${SCRIPT_DIR}/2.1_run_worker.sh")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected annotation OOM retry"
  ARRAY_ID="${retry_msg%%;*}"
  [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || fail "invalid annotation retry array id"
  SCHEDULER_IDS+=("${ARRAY_ID}")
  echo "ANNOTATION_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"; expected="${retry_count}"; CURRENT_MEMORY="${NEXT_MEMORY}"
done
validate_feathers "${manifest}" || fail "annotation feather schema/key validation failed"
status_write OK "all annotation chunks completed and validated"
