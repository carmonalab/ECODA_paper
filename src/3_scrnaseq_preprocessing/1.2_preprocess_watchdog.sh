#!/bin/bash
# Compute-node watchdog for the manifest-driven Stage 3 array.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
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
ecoda_open_run "${RUN_ID}" || exit 1
RUN_ROOT="${ECODA_RUN_ROOT}"
export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT
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
RUNTIME_EXPORT=""
stage3_manifest_value() {
  local manifest="$1" key="$2" value
  [[ -f "${manifest}" && ! -L "${manifest}" && -r "${manifest}" ]] || return 1
  value="$(awk -v wanted="${key}" '
    index($0, wanted "=") == 1 {
      count++
      result=substr($0, length(wanted) + 2)
    }
    END {
      if (count != 1 || result == "") exit 1
      print result
    }
  ' "${manifest}")" || return 1
  printf '%s\n' "${value}"
}

stage3_load_bound_run() {
  local source_copy="${RUN_ROOT}/manifests/source.manifest"
  local runtime_identity="${RUN_ROOT}/manifests/runtime.identity"
  local source_root source_manifest_original runtime_image runtime_manifest
  local snapshot_root runtime_format identity_count identity_image_sha identity_manifest_sha
  local identity_image_size identity_manifest_size
  [[ -s "${source_copy}" && ! -L "${source_copy}" && -r "${source_copy}" ]] || return 2
  [[ -s "${runtime_identity}" && ! -L "${runtime_identity}" && -r "${runtime_identity}" ]] || return 2
  ecoda_validate_run_owned_path "${source_copy}" "${RUN_ROOT}" || return 1
  ecoda_validate_run_owned_path "${runtime_identity}" "${RUN_ROOT}" || return 1
  [[ "$(stage3_manifest_value "${source_copy}" FORMAT)" == "1" ]] || return 1
  source_root="$(stage3_manifest_value "${source_copy}" SOURCE_ROOT)" || return 1
  [[ "${source_root}" = /* && "${source_root}" == */tree ]] || return 1
  snapshot_root="${source_root%/tree}"
  source_manifest_original="${snapshot_root}/identity/source.manifest"
  [[ -f "${source_manifest_original}" && ! -L "${source_manifest_original}" &&
     -r "${source_manifest_original}" ]] || return 1
  cmp -s "${source_copy}" "${source_manifest_original}" || {
    echo "ERROR: run-owned Stage 3 source manifest differs from immutable snapshot." >&2
    return 1
  }
  if [[ -n "${ECODA_SOURCE_MANIFEST_RUN:-}" ]]; then
    [[ "${ECODA_SOURCE_MANIFEST_RUN}" == "${source_copy}" ]] || return 1
  fi
  runtime_image="$(stage3_manifest_value "${runtime_identity}" RUNTIME_IMAGE)" || return 1
  identity_image_sha="$(stage3_manifest_value "${runtime_identity}" RUNTIME_IMAGE_SHA256)" || return 1
  identity_manifest_sha="$(stage3_manifest_value "${runtime_identity}" RUNTIME_MANIFEST_SHA256)" || return 1
  identity_image_size="$(stage3_manifest_value "${runtime_identity}" RUNTIME_IMAGE_SIZE)" || return 1
  identity_manifest_size="$(stage3_manifest_value "${runtime_identity}" RUNTIME_MANIFEST_SIZE)" || return 1
  runtime_manifest="$(stage3_manifest_value "${runtime_identity}" RUNTIME_MANIFEST)" || return 1
  [[ "${runtime_image}" = /* && "${runtime_manifest}" = /* &&
     -f "${runtime_manifest}" && ! -L "${runtime_manifest}" ]] || return 1
  if [[ -n "${ECODA_RUNTIME_IDENTITY:-}" ]]; then
    [[ "${ECODA_RUNTIME_IDENTITY}" == "${runtime_identity}" ]] || return 1
  fi
  identity_count="$(wc -l < "${runtime_identity}" | tr -d '[:space:]')" || return 1
  runtime_format="$(_ecoda_runtime_manifest_value "${runtime_manifest}" FORMAT 2>/dev/null || true)"
  case "${runtime_format}" in
    1) [[ "${identity_count}" == "6" ]] || return 1 ;;
    2) [[ "${identity_count}" == "8" ]] || return 1 ;;
    *) return 1 ;;
  esac
  if [[ -n "${ECODA_SOURCE_MANIFEST:-}" ]]; then
    cmp -s "${source_manifest_original}" "${ECODA_SOURCE_MANIFEST}" || return 1
  fi
  if [[ -n "${ECODA_SOURCE_ROOT:-}" ]]; then
    [[ "${ECODA_SOURCE_ROOT}" == "${source_root}" ]] || return 1
  fi
  if [[ -n "${ECODA_RUNTIME_IMAGE:-}" ]]; then
    [[ "${ECODA_RUNTIME_IMAGE}" == "${runtime_image}" ]] || return 1
  fi
  if [[ -n "${ECODA_RUNTIME_MANIFEST:-}" ]]; then
    [[ "${ECODA_RUNTIME_MANIFEST}" == "${runtime_manifest}" ]] || return 1
  fi
  if [[ "${runtime_format}" == "2" ]]; then
    [[ "${ECODA_RUNTIME_IMAGE_SHA256:-}" == "${identity_image_sha}" &&
       "${ECODA_RUNTIME_MANIFEST_SHA256:-}" == "${identity_manifest_sha}" &&
       "${ECODA_RUNTIME_IMAGE_SIZE:-}" == "${identity_image_size}" &&
       "${ECODA_RUNTIME_MANIFEST_SIZE:-}" == "${identity_manifest_size}" ]] || return 1
  fi
  SOURCE_ROOT="${source_root}"
  SOURCE_MANIFEST_ORIGINAL="${source_manifest_original}"
  SOURCE_MANIFEST_RUN="${source_copy}"
  RUNTIME_IDENTITY="${runtime_identity}"
  export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
  export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST_ORIGINAL}"
  export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
  export ECODA_AUX_ROOT="${source_root%/}/aux"
  export ECODA_RUNTIME_IMAGE="${runtime_image}"
  export ECODA_RUNTIME_MANIFEST="${runtime_manifest}"
  export ECODA_RUNTIME_IDENTITY="${runtime_identity}"
  export ECODA_RUN_ID="${RUN_ID}"
  export ECODA_RUN_ROOT="${RUN_ROOT}"
  PROJECT_ROOT="${source_root}"
  DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
  SCRIPT_DIR="${PROJECT_ROOT}/src/3_scrnaseq_preprocessing"
  export PROJECT_ROOT DATASETS_JSON_FILE
  LOGS_DIR="${ECODA_LOGS_DIR:-${LOGS_DIR:-${RUN_ROOT}/logs}}"
  export LOGS_DIR ECODA_LOGS_DIR="${LOGS_DIR}"
}

stage3_artifact_record_valid() {
  local path="$1" producer="$2" record
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  [[ -n "${record}" && -f "${record}" && ! -L "${record}" ]] || return 1
  ecoda_validate_artifact_record "${path}" "${producer}" "${RUN_ID}"
}

stage3_artifact_record_any() {
  local path="$1"
  stage3_artifact_record_valid "${path}" stage3 ||
    stage3_artifact_record_valid "${path}" stage3_preflight
}
STAGE3_OUTPUT_OWNER_DIRS_ALL=()

stage3_remember_output_owner() {
  local owner="$1" existing
  [[ -n "${owner}" ]] || return 1
  for existing in "${STAGE3_OUTPUT_OWNER_DIRS_ALL[@]:-}"; do
    [[ "${existing}" == "${owner}" ]] && return 0
  done
  STAGE3_OUTPUT_OWNER_DIRS_ALL+=("${owner}")
}

stage3_remember_output_owners() {
  local path owner
  for path in "${ECODA_OUTPUT_PATHS[@]:-}"; do
    [[ -n "${path}" ]] || continue
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    stage3_remember_output_owner "${owner}" || return 1
  done
  for owner in "${ECODA_OUTPUT_OWNER_DIRS[@]:-}"; do
    [[ -n "${owner}" ]] || continue
    stage3_remember_output_owner "${owner}" || return 1
  done
}
stage3_validate_scratch_output_ownership() {
  local selection="$1" saved_nas="${NAS_TARGET_DIR-}" had_nas=0 rc
  [[ -n "${NAS_TARGET_DIR+x}" ]] && had_nas=1
  unset NAS_TARGET_DIR
  ecoda_validate_output_ownership stage3 "${selection}" "${RUN_ID}"
  rc=$?
  if [[ ${had_nas} -eq 1 ]]; then
    export NAS_TARGET_DIR="${saved_nas}"
  else
    unset NAS_TARGET_DIR
  fi
  return "${rc}"
}

stage3_validate_output_owner() {
  local owner="$1" expected_owner
  _ecoda_artifact_owner_validate_dir "${owner}" || return 1
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${RUN_ID}" &&
     "${ECODA_ARTIFACT_OWNER_STAGE}" == stage3 ]] || return 2
  expected_owner="$(ecoda_artifact_owner_dir \
    "${ECODA_ARTIFACT_OWNER_CANONICAL_PATH}")" || return 1
  [[ "${owner}" == "${expected_owner}" ]]
}

stage3_finalize_output_owners() {
  local state="$1" reason="$2" owner rc=0
  if ! stage3_remember_output_owners; then
    rc=1
  fi
  for owner in "${STAGE3_OUTPUT_OWNER_DIRS_ALL[@]:-}"; do
    [[ -n "${owner}" ]] || continue
    if stage3_validate_output_owner "${owner}"; then
      ecoda_owner_set_state "${owner}" "${state}" "${reason}" || rc=1
    else
      # Never rewrite an owner belonging to another run or stage.
      rc=1
    fi
  done
  return "${rc}"
}


stage3_require_source_script() {
  local candidate="$1"
  [[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" &&
     "${SOURCE_ROOT:-${ECODA_SOURCE_ROOT:-}}" = /* ]] || return 1
  ecoda_require_source_script_path "${candidate}" "${SOURCE_ROOT:-${ECODA_SOURCE_ROOT}}"
}
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
  stage3_finalize_output_owners FAIL "${reason}" || owner_rc=1
  set_owner_state FAIL "${reason}" || owner_rc=1
  atomic_status FAIL "${reason}" || owner_rc=1
  exit 1
}
export ECODA_RUNTIME_PROFILE=stage3
set +e
stage3_load_bound_run
bound_rc=$?
set -e
if [[ ${bound_rc} -eq 2 ]]; then
  fail "legacy_source_unpinned"
fi
[[ ${bound_rc} -eq 0 ]] || fail "Stage 3 run-bound source/runtime identity is invalid"
ecoda_runtime_validate_bound_run ||
  fail "Stage 3 run-bound runtime validation failed before retry handling"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST_ORIGINAL}"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage3 0)" ||
  fail "Stage 3 runtime export construction failed"
RUNTIME_EXPORT="${RUNTIME_EXPORT},ECODA_RUNTIME_IDENTITY=${RUNTIME_IDENTITY},ECODA_SOURCE_MANIFEST_RUN=${SOURCE_MANIFEST_RUN},ECODA_RUNTIME_RUN_ID=${RUN_ID}"

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
    stage3_artifact_record_any "${path}" || return 1
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
stage3_validate_scratch_output_ownership "${CURRENT_MANIFEST}" ||
  fail "Stage 3 output ownership validation failed before watchdog handling"
stage3_remember_output_owners ||
  fail "failed to track Stage 3 global output owners"

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
stage3_validate_scratch_output_ownership "${RETRY_MANIFEST}" ||
    fail "Stage 3 output ownership validation failed before OOM retry"
  stage3_remember_output_owners ||
    fail "failed to track Stage 3 global output owners"
  retry_worker_script="$(stage3_require_source_script \
    "${SCRIPT_DIR}/1.1_run_worker.sh")" ||
    fail "Stage 3 retry worker script escaped immutable source root"
  set +e
  RETRY_MSG="$(sbatch --parsable --array="1-${RETRY_COUNT}%${THROTTLE}" \
    --mem="${NEXT_MEMORY}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX}_%A_%a.log" \
    --error="${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX}_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    --export="ALL,PREPROCESS_SELECTION_FILE=${RETRY_MANIFEST},PREPROCESS_RUN_ROOT=${RUN_ROOT},FORCE_PREPROCESS=1,PREPROCESS_ERROR_PREFIX=${LOGS_DIR}/3_scrnaseq_preprocessing_retry${RETRY_INDEX},${RUNTIME_EXPORT}" \
    "${retry_worker_script}")"
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
stage3_finalize_output_owners OK "Stage 3 preprocessing artifacts validated" ||
  fail "failed to finalize Stage 3 global artifact owners"
set_owner_state OK "Stage 3 preprocessing artifacts validated" ||
  fail "failed to finalize Stage 3 owners"
if ! atomic_status OK "all selected Stage 3 tasks completed and artifacts validated"; then
  fail "failed to write Stage 3 watchdog success status"
fi
printf 'Stage 3 watchdog completed for run %s\n' "${RUN_ID}"
