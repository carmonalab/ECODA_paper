#!/bin/bash
# Parameterized compute-node watchdog for all Stage 4 worker arrays.
# Scheduler accounting, OOM-only retries, and output ownership live here;
# phase-specific artifact checks are delegated to the validator command.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/../slurm_config.sh" ]]; then
  if [[ -n "${ECODA_SOURCE_ROOT:-}" &&
        -f "${ECODA_SOURCE_ROOT}/src/slurm_config.sh" ]]; then
    SCRIPT_DIR="${ECODA_SOURCE_ROOT}/src/4_cell_type_annotation"
  elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
          -f "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR}/src/4_cell_type_annotation"
  fi
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
cd "${PROJECT_ROOT}"

usage() {
  cat >&2 <<'EOF'
Usage: stage4_watchdog.sh PHASE RUN_ID ARRAY_MANIFEST OWNER_SELECTION ARRAY_ID MEM MAX_MEM PARTITION THROTTLE \
       RETRY_WORKER RETRY_MANIFEST_ENV RETRY_MANIFEST_STEM RETRY_LOG_PREFIX RETRY_TIME_LIMIT \
       RETRY_EXPORTS VALIDATOR_SCRIPT [VALIDATOR_ARGS...]

PHASE is preparation, annotation, or merge.  VALIDATOR_SCRIPT receives:
  PHASE RUN_ID ARRAY_MANIFEST OWNER_SELECTION [VALIDATOR_ARGS...]
EOF
}
[[ $# -ge 16 ]] || { usage; exit 2; }

PHASE="$1"
RUN_ID="$2"
ROOT_MANIFEST="$3"
OWNER_SELECTION="$4"
ARRAY_ID="$5"
CURRENT_MEMORY="$6"
MAX_MEMORY="$7"
PARTITION="$8"
THROTTLE="$9"
RETRY_WORKER_SCRIPT="${10}"
RETRY_MANIFEST_ENV="${11}"
RETRY_MANIFEST_STEM="${12}"
RETRY_LOG_PREFIX="${13}"
RETRY_TIME_LIMIT="${14}"
RETRY_EXPORTS="${15}"
VALIDATOR_SCRIPT="${16}"
shift 16
VALIDATOR_ARGS=("$@")
PHASE_UPPER="$(printf '%s' "${PHASE}" | tr '[:lower:]' '[:upper:]')"

case "${PHASE}" in
  preparation) STATUS_REASON="all selected dataset preparation artifacts validated" ;;
  annotation) STATUS_REASON="all annotation chunks completed and validated" ;;
  merge) STATUS_REASON="all selected dataset merges completed and validated" ;;
  *) echo "ERROR: unsupported Stage 4 watchdog phase: ${PHASE}" >&2; exit 2 ;;
esac
[[ "${RETRY_MANIFEST_ENV}" =~ ^[A-Z][A-Z0-9_]*$ ]] || {
  echo "ERROR: invalid Stage 4 retry manifest environment variable: ${RETRY_MANIFEST_ENV}" >&2
  exit 2
}
[[ "${RETRY_MANIFEST_STEM}" =~ ^[A-Za-z0-9_-]+$ ]] || {
  echo "ERROR: invalid Stage 4 retry manifest stem: ${RETRY_MANIFEST_STEM}" >&2
  exit 2
}
[[ "${ARRAY_ID}" =~ ^[0-9]+$ && "${THROTTLE}" =~ ^[1-9][0-9]*$ ]] || {
  echo "ERROR: invalid Stage 4 scheduler array configuration" >&2
  exit 2
}
[[ "${RETRY_TIME_LIMIT}" == "-" ]] && RETRY_TIME_LIMIT=""
[[ "${RETRY_EXPORTS}" == "-" ]] && RETRY_EXPORTS=""

RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
STATUS_FILE="${RUN_ROOT}/status/${PHASE}_watchdog"
CURRENT_MANIFEST="${ROOT_MANIFEST}"
RETRY_INDEX=0
SCHEDULER_IDS=("${ARRAY_ID}")
RUNTIME_EXPORT=""
STAGE4_OUTPUT_OWNER_DIRS=()
STAGE4_OUTPUT_OWNER_PATHS=()

stage4_validate_scratch_output_ownership() {
  local selection="$1" saved_nas="${NAS_TARGET_DIR-}" had_nas=0 rc
  [[ -n "${NAS_TARGET_DIR+x}" ]] && had_nas=1
  unset NAS_TARGET_DIR
  ecoda_validate_output_ownership stage4 "${selection}" "${RUN_ID}"
  rc=$?
  if [[ ${had_nas} -eq 1 ]]; then
    export NAS_TARGET_DIR="${saved_nas}"
  else
    unset NAS_TARGET_DIR
  fi
  return "${rc}"
}

stage4_validate_owner_scope() {
  local ds views run_root row_ds row_view extra view existing row matched
  local expected_root_real run_root_real
  local root_datasets=() expected_rows=() actual_rows=()
  if [[ "${PHASE}" != "annotation" ]]; then
    expected_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
    [[ -n "${expected_root_real}" ]] || return 1
  fi
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && -n "${views}" && -n "${run_root}" ]] || return 1
    matched=0
    for existing in "${root_datasets[@]-}"; do
      [[ "${existing}" == "${ds}" ]] && { matched=1; break; }
    done
    if [[ "${PHASE}" == "annotation" ]]; then
      ecoda_validate_run_owned_path "${views}" "${RUN_ROOT}" || return 1
      [[ "${views}" == "${RUN_ROOT}/datasets/${ds}/chunks/chunk_"*.txt &&
         "${run_root}" == "${RUN_ROOT}/datasets/${ds}/annotations" ]] || return 1
    else
      [[ ${matched} -eq 0 ]] || return 1
      run_root_real="$(ecoda_realpath_existing "${run_root}" 2>/dev/null || true)"
      [[ -n "${run_root_real}" && "${run_root_real}" == "${expected_root_real}" ]] || return 1
      ecoda_split_csv "${views}" || return 1
      for view in "${ECODA_ARRAY[@]}"; do
        expected_rows+=("${ds}"$'\t'"${view}")
      done
    fi
    [[ ${matched} -eq 1 ]] || root_datasets+=("${ds}")
  done < "${ROOT_MANIFEST}"
  while IFS=$'\t' read -r row_ds row_view extra; do
    [[ -n "${row_ds}" && -n "${row_view}" && -z "${extra}" ]] || return 1
    if [[ "${PHASE}" == "annotation" ]]; then
      matched=0
      for existing in "${root_datasets[@]-}"; do
        [[ "${existing}" == "${row_ds}" ]] && { matched=1; break; }
      done
      [[ ${matched} -eq 1 ]] || return 1
      matched=0
      for existing in "${actual_rows[@]-}"; do
        [[ "${existing}" == "${row_ds}" ]] && { matched=1; break; }
      done
      [[ ${matched} -eq 1 ]] || actual_rows+=("${row_ds}")
    else
      matched=0
      for existing in "${actual_rows[@]-}"; do
        [[ "${existing}" == "${row_ds}"$'\t'"${row_view}" ]] && { matched=1; break; }
      done
      [[ ${matched} -eq 0 ]] || return 1
      actual_rows+=("${row_ds}"$'\t'"${row_view}")
    fi
  done < "${OWNER_SELECTION}"
  if [[ "${PHASE}" == "annotation" ]]; then
    [[ ${#root_datasets[@]} -eq ${#actual_rows[@]} ]] || return 1
    for existing in "${root_datasets[@]}"; do
      matched=0
      for row in "${actual_rows[@]-}"; do
        [[ "${existing}" == "${row}" ]] && { matched=1; break; }
      done
      [[ ${matched} -eq 1 ]] || return 1
    done
    return 0
  fi
  [[ ${#expected_rows[@]} -eq ${#actual_rows[@]} ]] || return 1
  for existing in "${expected_rows[@]}"; do
    matched=0
    for row in "${actual_rows[@]-}"; do
      [[ "${existing}" == "${row}" ]] && { matched=1; break; }
    done
    [[ ${matched} -eq 1 ]] || return 1
  done
  return 0
}

stage4_collect_output_owners() {
  local path owner_dir owner_path write_flag existing seen
  local output_index=0
  ecoda_validate_run_owned_path "${OWNER_SELECTION}" "${RUN_ROOT}" || return 1
  ecoda_validate_manifest "${OWNER_SELECTION}" 2 || return 1
  stage4_validate_owner_scope || return 1
  stage4_validate_scratch_output_ownership "${OWNER_SELECTION}" || return 1
  for path in "${ECODA_OUTPUT_PATHS[@]}"; do
    write_flag="${ECODA_OUTPUT_WRITE_FLAGS[${output_index}]:-1}"
    output_index=$((output_index + 1))
    [[ "${write_flag}" == "1" ]] || continue
    owner_dir="$(ecoda_artifact_owner_dir "${path}")" || return 1
    [[ -d "${owner_dir}" ]] || return 1
    _ecoda_artifact_owner_validate_dir "${owner_dir}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${RUN_ID}" &&
       "${ECODA_ARTIFACT_OWNER_STAGE}" == "stage4" &&
       "${ECODA_ARTIFACT_OWNER_STATE}" == "ACTIVE" ]] || return 1
    seen=0
    if [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -gt 0 ]]; then
      for existing in "${STAGE4_OUTPUT_OWNER_DIRS[@]}"; do
        if [[ "${existing}" == "${owner_dir}" ]]; then
          seen=1
          break
        fi
      done
    fi
    if [[ ${seen} -eq 0 ]]; then
      STAGE4_OUTPUT_OWNER_DIRS+=("${owner_dir}")
      STAGE4_OUTPUT_OWNER_PATHS+=("${path}")
    fi
  done
}

stage4_finalize_output_owners() {
  local state="$1" reason="${2:-}" owner owner_path owner_state owner_run owner_stage
  local rc=0 output_index=0
  [[ "${state}" == "OK" || "${state}" == "FAIL" ]] || return 1
  [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -gt 0 ]] || return 0
  [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -eq ${#STAGE4_OUTPUT_OWNER_PATHS[@]} ]] || return 1
  for owner in "${STAGE4_OUTPUT_OWNER_DIRS[@]}"; do
    owner_path="${STAGE4_OUTPUT_OWNER_PATHS[${output_index}]:-}"
    output_index=$((output_index + 1))
    [[ -n "${owner_path}" ]] || { rc=1; continue; }
    _ecoda_artifact_owner_validate_dir "${owner}" "${owner_path}" || { rc=1; continue; }
    owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
    owner_stage="${ECODA_ARTIFACT_OWNER_STAGE}"
    owner_state="${ECODA_ARTIFACT_OWNER_STATE}"
    [[ "${owner_run}" == "${RUN_ID}" && "${owner_stage}" == "stage4" ]] || {
      rc=1
      continue
    }
    if [[ "${state}" == "OK" ]]; then
      [[ "${owner_state}" == "ACTIVE" ]] || { rc=1; continue; }
    else
      [[ "${owner_state}" == "ACTIVE" || "${owner_state}" == "FAIL" ]] || {
        rc=1
        continue
      }
    fi
  done
  [[ ${rc} -eq 0 ]] || return 1
  output_index=0
  for owner in "${STAGE4_OUTPUT_OWNER_DIRS[@]}"; do
    owner_path="${STAGE4_OUTPUT_OWNER_PATHS[${output_index}]:-}"
    output_index=$((output_index + 1))
    _ecoda_artifact_owner_validate_dir "${owner}" "${owner_path}" || { rc=1; continue; }
    owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
    owner_stage="${ECODA_ARTIFACT_OWNER_STAGE}"
    owner_state="${ECODA_ARTIFACT_OWNER_STATE}"
    [[ "${owner_run}" == "${RUN_ID}" && "${owner_stage}" == "stage4" &&
       "${owner_state}" == "ACTIVE" ]] || { rc=1; continue; }
    ecoda_owner_set_state "${owner}" "${state}" "${reason}" || rc=1
  done
  return "${rc}"
}

stage4_assert_output_owners_active() {
  local owner owner_path owner_run owner_stage owner_state
  local output_index=0
  [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -gt 0 &&
     ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -eq ${#STAGE4_OUTPUT_OWNER_PATHS[@]} ]] || return 1
  for owner in "${STAGE4_OUTPUT_OWNER_DIRS[@]}"; do
    owner_path="${STAGE4_OUTPUT_OWNER_PATHS[${output_index}]:-}"
    output_index=$((output_index + 1))
    [[ -n "${owner_path}" ]] || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${owner_path}" || return 1
    owner_run="${ECODA_ARTIFACT_OWNER_RUN}"
    owner_stage="${ECODA_ARTIFACT_OWNER_STAGE}"
    owner_state="${ECODA_ARTIFACT_OWNER_STATE}"
    [[ "${owner_run}" == "${RUN_ID}" &&
       "${owner_stage}" == "stage4" &&
       "${owner_state}" == "ACTIVE" ]] || return 1
  done
}

status_write() {
  local state="$1" reason="${2:-}" tmp="${STATUS_FILE}.tmp.$$"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\nPHASE=%s\nRUN_ID=%s\nSOURCE_ROOT=%s\nSOURCE_MANIFEST=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nREASON=%s\nARRAY_JOB_ID=%s\nRETRY_INDEX=%s\n' \
      "${state}" "${PHASE}" "${RUN_ID}" "${ECODA_SOURCE_ROOT:-}" "${ECODA_SOURCE_MANIFEST:-}" \
      "${ECODA_RUNTIME_IMAGE:-}" "${ECODA_RUNTIME_MANIFEST:-}" "${reason}" "${ARRAY_ID}" "${RETRY_INDEX}"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      printf 'SCHEDULER_ID=%s\n' "${SLURM_JOB_ID}"
    fi
    local scheduler_id
    for scheduler_id in "${SCHEDULER_IDS[@]}"; do
      printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"
    done
  } > "${tmp}"
  mv -f "${tmp}" "${STATUS_FILE}"
}

fail() {
  local reason="$1"
  ecoda_owner_finalize_tracked FAIL "${reason}" >/dev/null 2>&1 || true
  stage4_finalize_output_owners FAIL "${reason}" >/dev/null 2>&1 || true
  echo "ERROR: ${reason}" >&2
  status_write FAIL "${reason}"
  exit 1
}

load_bound_identity() {
  local source_copy="${RUN_ROOT}/manifests/source.manifest"
  local runtime_identity="${RUN_ROOT}/manifests/runtime.identity"
  local source_root expected_source_manifest runtime_image runtime_manifest
  local metadata_stage metadata_run
  [[ -s "${source_copy}" && ! -L "${source_copy}" &&
     -s "${runtime_identity}" && ! -L "${runtime_identity}" ]] || {
    echo "ERROR: legacy_source_unpinned: Stage 4 run lacks source/runtime identity manifests." >&2
    return 1
  }
  metadata_stage="$(sed -n 's/^STAGE=//p' "${RUN_ROOT}/metadata" | head -1 || true)"
  metadata_run="$(sed -n 's/^RUN_ID=//p' "${RUN_ROOT}/metadata" | head -1 || true)"
  [[ "${metadata_stage}" == "stage4" && "${metadata_run}" == "${RUN_ID}" ]] || return 1
  source_root="$(sed -n 's/^SOURCE_ROOT=//p' "${source_copy}" | head -1)"
  [[ "${source_root}" = /* && "${source_root##*/}" == tree ]] || return 1
  expected_source_manifest="${source_root%/tree}/identity/source.manifest"
  [[ -f "${expected_source_manifest}" && -r "${expected_source_manifest}" ]] || return 1
  cmp -s "${source_copy}" "${expected_source_manifest}" || return 1
  runtime_image="$(sed -n 's/^RUNTIME_IMAGE=//p' "${runtime_identity}" | head -1)"
  runtime_manifest="$(sed -n 's/^RUNTIME_MANIFEST=//p' "${runtime_identity}" | head -1)"
  [[ "${runtime_image}" = /* && "${runtime_manifest}" = /* ]] || return 1
  export ECODA_SOURCE_ROOT="${source_root}"
  export ECODA_SOURCE_MANIFEST="${expected_source_manifest}"
  export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
  export ECODA_AUX_ROOT="${source_root}/aux"
  export ECODA_RUNTIME_IMAGE="${runtime_image}"
  export ECODA_RUNTIME_MANIFEST="${runtime_manifest}"
  export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT="${RUN_ROOT}"
  export PROJECT_ROOT="${source_root}"
  export DATASETS_JSON_FILE="${source_root}/datasets.json"
  export SCGATE_DB_PATH="${ECODA_AUX_ROOT}/scGateDB.rds"
}

bump_mem() {
  [[ "$1" =~ ^([0-9]+)([GT])$ ]] || return 1
  printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"
}

mem_ge() {
  local a="$1" b="$2" an as bn bs
  [[ "${a}" =~ ^([0-9]+)([GT])$ ]] || return 1
  an="${BASH_REMATCH[1]}"; as="${BASH_REMATCH[2]}"
  [[ "${b}" =~ ^([0-9]+)([GT])$ ]] || return 1
  bn="${BASH_REMATCH[1]}"; bs="${BASH_REMATCH[2]}"
  [[ "${as}" == T ]] && an=$((an * 1024))
  [[ "${bs}" == T ]] && bn=$((bn * 1024))
  (( an >= bn ))
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

submit_retry() {
  local retry_manifest="$1" retry_count export_spec retry_exports retry_msg retry_rc
  local sbatch_args=()
  retry_count="$(wc -l < "${retry_manifest}" | tr -d '[:space:]')"
  [[ "${retry_count}" =~ ^[1-9][0-9]*$ ]] || return 1
  export_spec="ALL,${RETRY_MANIFEST_ENV}=${retry_manifest}"
  retry_exports="${RETRY_EXPORTS//@RETRY_INDEX@/${RETRY_INDEX}}"
  [[ -n "${retry_exports}" ]] && export_spec="${export_spec},${retry_exports}"
  export_spec="${export_spec},ANNOTATION_RUN_ID=${RUN_ID}"
  [[ -n "${RUNTIME_EXPORT}" ]] && export_spec="${export_spec},${RUNTIME_EXPORT}"
  sbatch_args=(--parsable --array="1-${retry_count}%${THROTTLE}"
    --mem="${NEXT_MEMORY}" --partition="${PARTITION}")
  [[ -n "${RETRY_TIME_LIMIT}" ]] &&
    sbatch_args+=(--time="${RETRY_TIME_LIMIT}")
  sbatch_args+=(--output="${RETRY_LOG_PREFIX}${RETRY_INDEX}_%A_%a.log"
    --error="${RETRY_LOG_PREFIX}${RETRY_INDEX}_%A_%a.err"
    --mail-user="${USER_EMAIL}" --export="${export_spec}" "${RETRY_WORKER_SCRIPT}")
  set +e
  retry_msg="$(sbatch "${sbatch_args[@]}")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || return 1
  retry_msg="${retry_msg%%;*}"
  [[ "${retry_msg}" =~ ^[0-9]+$ ]] || return 1
  printf '%s' "${retry_msg}"
}

export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT="${RUN_ROOT}"
ecoda_validate_run_id "${RUN_ID}" || exit 1
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: Stage 4 run root is missing: ${RUN_ROOT}" >&2; exit 1; }
load_bound_identity || fail "legacy_source_unpinned: Stage 4 ${PHASE} run is not source/runtime pinned"
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_bound_run || \
  fail "Stage 4 bound runtime validation failed before ${PHASE} retry handling"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage4 0)" || \
  fail "Stage 4 runtime export construction failed"
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" || \
  fail "Stage 4 ${PHASE} array manifest is outside the run root"
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 || \
  fail "Stage 4 ${PHASE} array manifest is invalid"
ecoda_validate_run_owned_path "${OWNER_SELECTION}" "${RUN_ROOT}" || \
  fail "Stage 4 ${PHASE} ownership selection is outside the run root"
ecoda_validate_manifest "${OWNER_SELECTION}" 2 || \
  fail "Stage 4 ${PHASE} ownership selection is invalid"
RETRY_WORKER_SCRIPT="$(ecoda_require_source_script_path \
  "${RETRY_WORKER_SCRIPT}" "${ECODA_SOURCE_ROOT}")" || \
  fail "Stage 4 ${PHASE} retry worker is outside the immutable source root"
VALIDATOR_SCRIPT="$(ecoda_require_source_script_path \
  "${VALIDATOR_SCRIPT}" "${ECODA_SOURCE_ROOT}")" || \
  fail "Stage 4 ${PHASE} validator is outside the immutable source root"
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "Stage 4 ${PHASE} array manifest is empty"
stage4_collect_output_owners || \
  fail "Stage 4 ${PHASE} output ownership validation failed before watchdog accounting"

while :; do
  classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal ${PHASE} rows"
  [[ ${#FAILED_TASKS[@]} -eq 0 ]] || fail "non-OOM ${PHASE} failure: ${FAILED_TASKS[*]}"
  [[ ${#OOM_TASKS[@]} -eq 0 ]] && break
  mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}" && \
    fail "${PHASE} OOM at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  NEXT_MEMORY="$(bump_mem "${CURRENT_MEMORY}")" || fail "unparseable ${PHASE} memory"
  mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}" && NEXT_MEMORY="${MAX_MEMORY}"
  RETRY_INDEX=$((RETRY_INDEX + 1))
  [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded ${PHASE} OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/${RETRY_MANIFEST_STEM}.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    line="$(sed -n "${task}p" "${CURRENT_MANIFEST}")" || \
      fail "failed to read ${PHASE} retry row ${task}"
    [[ -n "${line}" ]] || fail "${PHASE} OOM task ${task} is absent from the manifest"
    printf '%s\n' "${line}" >> "${RETRY_TMP}" || \
      fail "failed to build ${PHASE} retry manifest"
  done
  ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 3 || \
    fail "failed to install ${PHASE} retry manifest atomically"
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" || \
    fail "${PHASE} retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 3 || \
    fail "${PHASE} retry manifest is invalid"
  stage4_collect_output_owners || \
    fail "Stage 4 ${PHASE} retry output ownership validation failed"
  RETRY_ARRAY_ID="$(submit_retry "${RETRY_MANIFEST}")" || \
    fail "sbatch rejected ${PHASE} OOM retry"
  SCHEDULER_IDS+=("${RETRY_ARRAY_ID}")
  echo "STAGE4_${PHASE_UPPER}_RETRY_ARRAY_JOB_ID=${RETRY_ARRAY_ID}"
  ARRAY_ID="${RETRY_ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"
  expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
  CURRENT_MEMORY="${NEXT_MEMORY}"
done

if [[ ${#VALIDATOR_ARGS[@]} -gt 0 ]]; then
  bash "${VALIDATOR_SCRIPT}" "${PHASE}" "${RUN_ID}" "${ROOT_MANIFEST}" \
    "${OWNER_SELECTION}" "${VALIDATOR_ARGS[@]}"
else
  bash "${VALIDATOR_SCRIPT}" "${PHASE}" "${RUN_ID}" "${ROOT_MANIFEST}" \
    "${OWNER_SELECTION}"
fi || fail "Stage 4 ${PHASE} output validation failed"
case "${PHASE}" in
  merge)
    stage4_finalize_output_owners OK "validated by Stage 4 merge watchdog" || \
      fail "failed to finalize Stage 4 merge output owners"
    status_write OK "${STATUS_REASON}"
    ;;
  preparation|annotation)
    stage4_assert_output_owners_active || \
      fail "Stage 4 ${PHASE} output owners are not active after validation"
    status_write OK "${STATUS_REASON}"
    ;;
esac
