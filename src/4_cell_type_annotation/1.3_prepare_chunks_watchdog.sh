#!/bin/bash
# Compute-node watchdog for dataset-parallel annotation preparation.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"

[[ $# -eq 7 ]] || { echo "Usage: 1.3_prepare_chunks_watchdog.sh RUN_ID MANIFEST ARRAY_ID MEM MAX_MEM PARTITION THROTTLE" >&2; exit 2; }
RUN_ID="$1"; ROOT_MANIFEST="$2"; ARRAY_ID="$3"; CURRENT_MEMORY="$4"; MAX_MEMORY="$5"; PARTITION="$6"; THROTTLE="$7"
ecoda_validate_run_id "${RUN_ID}" || exit 1
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: Stage 4 run root is missing: ${RUN_ROOT}" >&2; exit 1; }
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  { echo "ERROR: annotation prep manifest is outside run root." >&2; exit 1; }
STATUS_FILE="${RUN_ROOT}/status/preparation_watchdog"
CURRENT_MANIFEST="${ROOT_MANIFEST}"
RETRY_INDEX=0
SCHEDULER_IDS=("${ARRAY_ID}")
RUNTIME_EXPORT=""
STAGE4_OUTPUT_OWNER_DIRS=()
STAGE4_OUTPUT_OWNER_PATHS=()

stage4_build_owner_selection() {
  local manifest="$1" output="$2"
  local ds views run_root view expected_root_real run_root_real existing
  local seen_datasets=""
  : > "${output}" || return 1
  expected_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
  [[ -n "${expected_root_real}" ]] || return 1
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && -n "${views}" && -n "${run_root}" ]] || return 1
    run_root_real="$(ecoda_realpath_existing "${run_root}" 2>/dev/null || true)"
    [[ "${run_root_real}" == "${expected_root_real}" ]] || return 1
    case " ${seen_datasets} " in
      *" ${ds} "*) return 1 ;;
    esac
    seen_datasets="${seen_datasets} ${ds}"
    ecoda_split_csv "${views}" || return 1
    for view in "${ECODA_ARRAY[@]}"; do
      printf '%s\t%s\n' "${ds}" "${view}" >> "${output}" || return 1
    done
  done < "${manifest}"
  [[ -s "${output}" ]]
}

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

stage4_collect_output_owners() {
  local manifest="$1" selection path owner_dir existing owner_path write_flag
  local output_index=0 seen
  selection="${RUN_ROOT}/manifests/.stage4_preparation_selection_$$-${RANDOM}.tsv"
  stage4_build_owner_selection "${manifest}" "${selection}" || {
    rm -f "${selection}"
    return 1
  }
  ecoda_validate_run_owned_path "${selection}" "${RUN_ROOT}" || {
    rm -f "${selection}"
    return 1
  }
  ecoda_validate_manifest "${selection}" 2 || {
    rm -f "${selection}"
    return 1
  }
  if ! stage4_validate_scratch_output_ownership "${selection}"; then
    rm -f "${selection}"
    return 1
  fi
  rm -f "${selection}" || return 1
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
    [[ "${owner_run}" == "${RUN_ID}" && "${owner_stage}" == "stage4" ]] ||
      { rc=1; continue; }
    if [[ "${state}" == "OK" ]]; then
      [[ "${owner_state}" == "ACTIVE" ]] || { rc=1; continue; }
    else
      [[ "${owner_state}" == "ACTIVE" || "${owner_state}" == "FAIL" ]] ||
        { rc=1; continue; }
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
  [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -gt 0 ]] || return 1
  [[ ${#STAGE4_OUTPUT_OWNER_DIRS[@]} -eq ${#STAGE4_OUTPUT_OWNER_PATHS[@]} ]] || return 1
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
    printf 'STATE=%s\nRUN_ID=%s\nSOURCE_ROOT=%s\nSOURCE_MANIFEST=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nREASON=%s\nARRAY_JOB_ID=%s\nRETRY_INDEX=%s\n' \
      "${state}" "${RUN_ID}" "${ECODA_SOURCE_ROOT:-}" "${ECODA_SOURCE_MANIFEST:-}" \
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

export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT="${RUN_ROOT}"
load_bound_identity || fail "legacy_source_unpinned: Stage 4 preparation run is not source/runtime pinned"
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_bound_run || \
  fail "Stage 4 bound runtime validation failed before preparation retry handling"
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage4 0)" || \
  fail "Stage 4 runtime export construction failed"
STAGE4_SELECTION_MANIFEST="${RUN_ROOT}/manifests/runnable_selection.tsv"
if [[ ! -r "${STAGE4_SELECTION_MANIFEST}" ]]; then
  STAGE4_SELECTION_MANIFEST="${RUN_ROOT}/manifests/selection.tsv"
fi
ecoda_validate_run_owned_path "${STAGE4_SELECTION_MANIFEST}" "${RUN_ROOT}" ||
  fail "Stage 4 ownership selection manifest is outside the run root"
ecoda_validate_manifest "${STAGE4_SELECTION_MANIFEST}" 2 ||
  fail "Stage 4 ownership selection manifest is invalid"
bump_mem() { [[ "$1" =~ ^([0-9]+)([GT])$ ]] || return 1; printf '%s%s' "$((BASH_REMATCH[1] * 2))" "${BASH_REMATCH[2]}"; }
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
validate_artifact() {
  local path="$1" record=""
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  if [[ -n "${record}" && -e "${record}" ]]; then
    ecoda_validate_artifact_record "${path}" "${RUN_ID}" "${RUN_ID}"
  else
    ecoda_validate_checksum "${path}"
  fi
}

validate_rows() {
  local manifest="$1" ds views run_root expected_root_real run_root_real
  local prepare_script="${ECODA_SOURCE_ROOT}/src/4_cell_type_annotation/1.1_prepare_chunks.py"
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && -n "${views}" && -n "${run_root}" ]] || return 1
    expected_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
    run_root_real="$(ecoda_realpath_existing "${run_root}" 2>/dev/null || true)"
    [[ -n "${expected_root_real}" && "${expected_root_real}" == "${run_root_real}" ]] || return 1
    export DS_NAME="${ds}" ANNOTATION_VIEWS="${views}" ANNOTATION_RUN_ROOT="${run_root}" ANNOTATION_RUN_ID="${RUN_ID}"
    "${PYTHON_BIN}" "${prepare_script}" --views "${views}" --run-root "${run_root}" --validate-only || return 1
    union_path="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    validate_artifact "${union_path}" || return 1
  done < "${manifest}"
}

ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" ||
  fail "annotation prep manifest is outside run root"
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 ||
  fail "annotation prep manifest is invalid"
expected="$(wc -l < "${CURRENT_MANIFEST}" | tr -d '[:space:]')"
[[ "${expected}" =~ ^[1-9][0-9]*$ ]] || fail "annotation prep manifest is empty"
stage4_collect_output_owners "${ROOT_MANIFEST}" ||
  fail "Stage 4 preparation output ownership validation failed before watchdog accounting"

while :; do
  classify "${ARRAY_ID}" "${expected}" || fail "sacct did not provide terminal preparation rows"
  [[ ${#FAILED_TASKS[@]} -eq 0 ]] || fail "non-OOM preparation failure: ${FAILED_TASKS[*]}"
  [[ ${#OOM_TASKS[@]} -eq 0 ]] && break
  mem_ge "${CURRENT_MEMORY}" "${MAX_MEMORY}" && fail "preparation OOM at ${MAX_MEMORY} ceiling: ${OOM_TASKS[*]}"
  NEXT_MEMORY="$(bump_mem "${CURRENT_MEMORY}")" || fail "unparseable preparation memory"
  mem_ge "${NEXT_MEMORY}" "${MAX_MEMORY}" && NEXT_MEMORY="${MAX_MEMORY}"
  RETRY_INDEX=$((RETRY_INDEX + 1)); [[ ${RETRY_INDEX} -le 4 ]] || fail "exceeded preparation OOM retry attempts"
  RETRY_MANIFEST="${RUN_ROOT}/manifests/preparation.retry_${RETRY_INDEX}.tsv"
  RETRY_TMP="${RETRY_MANIFEST}.build.$$"
  : > "${RETRY_TMP}"
  for task in "${OOM_TASKS[@]}"; do
    sed -n "${task}p" "${CURRENT_MANIFEST}" >> "${RETRY_TMP}"
  done
  if ! ecoda_atomic_install_manifest "${RETRY_TMP}" "${RETRY_MANIFEST}" 3; then
    fail "failed to install preparation retry manifest atomically"
  fi
  rm -f "${RETRY_TMP}"
  ecoda_validate_run_owned_path "${RETRY_MANIFEST}" "${RUN_ROOT}" ||
    fail "preparation retry manifest escaped the run root"
  ecoda_validate_manifest "${RETRY_MANIFEST}" 3 ||
    fail "preparation retry manifest is invalid"
  retry_count="$(wc -l < "${RETRY_MANIFEST}" | tr -d '[:space:]')"
  stage4_collect_output_owners "${ROOT_MANIFEST}" ||
    fail "Stage 4 preparation retry output ownership validation failed"
  retry_script="$(ecoda_require_source_script_path \
    "${SCRIPT_DIR}/1.2_prepare_chunks_worker.sh" "${ECODA_SOURCE_ROOT}")" ||
    fail "Stage 4 preparation retry worker is outside the immutable source root"
  set +e
  retry_msg="$(sbatch --parsable --array="1-${retry_count}%${THROTTLE}" --mem="${NEXT_MEMORY}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/4_annotation_prepare_retry${RETRY_INDEX}_%A_%a.log" \
    --error="${LOGS_DIR}/4_annotation_prepare_retry${RETRY_INDEX}_%A_%a.err" --mail-user="${USER_EMAIL}" \
    --export="ALL,ANNOTATION_PREP_MANIFEST=${RETRY_MANIFEST},ANNOTATION_TEST_MODE=0,FORCE_ANNOTATION=1,ANNOTATION_RUN_ID=${RUN_ID},${RUNTIME_EXPORT}" \
    "${retry_script}")"
  retry_rc=$?
  set -e
  [[ ${retry_rc} -eq 0 ]] || fail "sbatch rejected preparation OOM retry"
  ARRAY_ID="${retry_msg%%;*}"
  [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || fail "invalid preparation retry array id"
  SCHEDULER_IDS+=("${ARRAY_ID}")
  echo "ANNOTATION_PREP_RETRY_ARRAY_JOB_ID=${ARRAY_ID}"
  CURRENT_MANIFEST="${RETRY_MANIFEST}"; expected="${retry_count}"; CURRENT_MEMORY="${NEXT_MEMORY}"
done
validate_rows "${ROOT_MANIFEST}" || fail "annotation preparation artifact validation failed"
# Preparation is non-terminal: the submitter's final H5AD owners stay ACTIVE
# until the merge watchdog validates the merged outputs.
stage4_assert_output_owners_active ||
  fail "Stage 4 preparation output owners are not active after validation"
status_write OK "all selected dataset preparation artifacts validated"
