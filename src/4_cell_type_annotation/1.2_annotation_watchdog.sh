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
STAGE4_OUTPUT_OWNER_DIRS=()
STAGE4_OUTPUT_OWNER_PATHS=()

stage4_build_owner_selection() {
  local manifest="$1" output="$2"
  local ds chunk feather_dir sel_ds sel_view extra existing
  local selected_datasets=() matched_datasets=() matched
  : > "${output}" || return 1
  while IFS=$'\t' read -r ds chunk feather_dir; do
    [[ -n "${ds}" && -n "${chunk}" && -n "${feather_dir}" ]] || return 1
    ecoda_validate_run_owned_path "${chunk}" "${RUN_ROOT}" || return 1
    [[ "${chunk}" == "${RUN_ROOT}/datasets/${ds}/chunks/chunk_"*.txt &&
       "${feather_dir}" == "${RUN_ROOT}/datasets/${ds}/annotations" ]] || return 1
    matched=0
    if [[ ${#selected_datasets[@]} -gt 0 ]]; then
      for existing in "${selected_datasets[@]}"; do
        [[ "${existing}" == "${ds}" ]] && { matched=1; break; }
      done
    fi
    [[ ${matched} -eq 1 ]] || selected_datasets+=("${ds}")
  done < "${manifest}"
  while IFS=$'\t' read -r sel_ds sel_view extra; do
    [[ -n "${sel_ds}" && -n "${sel_view}" && -z "${extra}" ]] || return 1
    matched=0
    if [[ ${#selected_datasets[@]} -gt 0 ]]; then
      for ds in "${selected_datasets[@]}"; do
        if [[ "${ds}" == "${sel_ds}" ]]; then
          printf '%s\t%s\n' "${sel_ds}" "${sel_view}" >> "${output}" || return 1
          matched=1
          break
        fi
      done
    fi
    [[ ${matched} -eq 0 ]] && continue
    matched=0
    if [[ ${#matched_datasets[@]} -gt 0 ]]; then
      for existing in "${matched_datasets[@]}"; do
        [[ "${existing}" == "${sel_ds}" ]] && { matched=1; break; }
      done
    fi
    [[ ${matched} -eq 1 ]] || matched_datasets+=("${sel_ds}")
  done < "${STAGE4_SELECTION_MANIFEST}"
  for ds in "${selected_datasets[@]}"; do
    matched=0
    if [[ ${#matched_datasets[@]} -gt 0 ]]; then
      for existing in "${matched_datasets[@]}"; do
        [[ "${existing}" == "${ds}" ]] && { matched=1; break; }
      done
    fi
    [[ ${matched} -eq 1 ]] || return 1
  done
  [[ -s "${output}" ]]
}

stage4_collect_output_owners() {
  local manifest="$1" selection path owner_dir existing owner_path write_flag
  local output_index=0 seen
  selection="${RUN_ROOT}/manifests/.stage4_annotation_selection_$$-${RANDOM}.tsv"
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
  if ! ecoda_validate_output_ownership stage4 "${selection}" "${RUN_ID}"; then
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



STAGE4_ARTIFACT_RECORD_USED=0
status_write() {
  local state="$1" reason="${2:-}" tmp="${STATUS_FILE}.tmp.$$"
  mkdir -p "$(dirname "${STATUS_FILE}")"
  {
    printf 'STATE=%s\nRUN_ID=%s\nSOURCE_ROOT=%s\nSOURCE_MANIFEST=%s\nRUNTIME_IMAGE=%s\nRUNTIME_MANIFEST=%s\nREASON=%s\nARRAY_JOB_ID=%s\nRETRY_INDEX=%s\n' \
      "${state}" "${RUN_ID}" "${ECODA_SOURCE_ROOT:-}" "${ECODA_SOURCE_MANIFEST:-}" \
      "${ECODA_RUNTIME_IMAGE:-}" "${ECODA_RUNTIME_MANIFEST:-}" "${reason}" "${ARRAY_ID}" "${RETRY_INDEX}"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then printf 'SCHEDULER_ID=%s\n' "${SLURM_JOB_ID}"; fi
    local scheduler_id
    for scheduler_id in "${SCHEDULER_IDS[@]}"; do printf 'SCHEDULER_ID=%s\n' "${scheduler_id}"; done
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
load_bound_identity || fail "legacy_source_unpinned: Stage 4 annotation run is not source/runtime pinned"
export ECODA_RUNTIME_PROFILE=stage4
ecoda_runtime_validate_bound_run || \
  fail "Stage 4 bound runtime validation failed before annotation retry handling"
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
validate_artifact() {
  local path="$1" producer="${2:-${RUN_ID}}" record=""
  STAGE4_ARTIFACT_RECORD_USED=0
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  if [[ -n "${record}" && -e "${record}" ]]; then
    ecoda_validate_artifact_record "${path}" "${producer}" "${RUN_ID}" || return 1
    STAGE4_ARTIFACT_RECORD_USED=1
  else
    ecoda_validate_checksum "${path}"
  fi
}
validate_feathers() {
  local manifest="$1" ds chunk feather_dir chunk_num feather union_path sample expected_union expected_chunk_dir
  local feather_sidecar_flag
  local expected_args validated_unions=""
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
    case " ${validated_unions} " in
      *" ${union_path} "*) ;;
      *)
        validate_artifact "${union_path}" || return 1
        validated_unions="${validated_unions} ${union_path}"
        ;;
    esac
    expected_args=(--expected-union "${union_path}")
    while IFS= read -r sample; do
      [[ -n "${sample}" ]] || return 1
      expected_args+=(--expected-sample "${sample}")
    done < <(sed -n '2,$p' "${chunk}")
    validate_artifact "${feather}" stage4_annotation || return 1
    feather_sidecar_flag="--require-sidecar"
    [[ ${STAGE4_ARTIFACT_RECORD_USED} -eq 1 ]] && feather_sidecar_flag="--sidecar-validated"
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${feather}" "${feather_sidecar_flag}" "${expected_args[@]}" >/dev/null 2>&1 || return 1
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
stage4_collect_output_owners "${ROOT_MANIFEST}" ||
  fail "Stage 4 annotation output ownership validation failed before watchdog accounting"

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
  stage4_collect_output_owners "${ROOT_MANIFEST}" ||
    fail "Stage 4 annotation retry output ownership validation failed"
  retry_script="$(ecoda_require_source_script_path \
    "${SCRIPT_DIR}/2.1_run_worker.sh" "${ECODA_SOURCE_ROOT}")" ||
    fail "Stage 4 annotation retry worker is outside the immutable source root"
  set +e
  retry_msg="$(sbatch --parsable --array="1-${retry_count}%${THROTTLE}" --mem="${NEXT_MEMORY}" --time="${ANNOTATION_WORKER_TIME_LIMIT}" --partition="${PARTITION}" \
    --output="${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX}_%A_%a.log" --error="${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX}_%A_%a.err" --mail-user="${USER_EMAIL}" \
    --export="ALL,CHUNKS_MANIFEST=${RETRY_MANIFEST},ANNOTATION_RUN_ID=${RUN_ID},ANNOTATION_ERROR_PREFIX=${LOGS_DIR}/4_cell_type_annotation_retry${RETRY_INDEX},${RUNTIME_EXPORT}" "${retry_script}")"
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
# Annotation is non-terminal: the submitter's final H5AD owners stay ACTIVE
# until the merge watchdog validates the merged outputs.
stage4_assert_output_owners_active ||
  fail "Stage 4 annotation output owners are not active after validation"
status_write OK "all annotation chunks completed and validated"
