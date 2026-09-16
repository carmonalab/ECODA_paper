#!/bin/bash
# Phase-specific Stage 4 artifact validation used by stage4_watchdog.sh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

[[ $# -ge 4 ]] || {
  echo "Usage: stage4_validate_outputs.sh PHASE RUN_ID ARRAY_MANIFEST OWNER_SELECTION [ARGS...]" >&2
  exit 2
}
PHASE="$1"
RUN_ID="$2"
ROOT_MANIFEST="$3"
OWNER_SELECTION="$4"
shift 4

case "${PHASE}" in
  preparation|annotation|merge) ;;
  *) echo "ERROR: unsupported Stage 4 output-validation phase: ${PHASE}" >&2; exit 2 ;;
esac
ecoda_validate_run_id "${RUN_ID}" || exit 1
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: Stage 4 run root is missing: ${RUN_ROOT}" >&2; exit 1; }
ecoda_validate_run_owned_path "${ROOT_MANIFEST}" "${RUN_ROOT}" || exit 1
ecoda_validate_manifest "${ROOT_MANIFEST}" 3 || exit 1
ecoda_validate_run_owned_path "${OWNER_SELECTION}" "${RUN_ROOT}" || exit 1
ecoda_validate_manifest "${OWNER_SELECTION}" 2 || exit 1

ARTIFACT_RECORD_USED=0
validate_artifact() {
  local path="$1" producer="${2:-${RUN_ID}}" record=""
  ARTIFACT_RECORD_USED=0
  record="$(ecoda_artifact_record_path "${path}" "${RUN_ID}" 2>/dev/null || true)"
  if [[ -n "${record}" && -f "${record}" && ! -L "${record}" ]]; then
    ecoda_validate_artifact_record "${path}" "${producer}" "${RUN_ID}" || return 1
    ARTIFACT_RECORD_USED=1
    ECODA_CHECKSUM_MD5="${ECODA_ARTIFACT_RECORD_MD5}"
    ECODA_CHECKSUM_SIZE="${ECODA_ARTIFACT_RECORD_SIZE}"
  else
    ecoda_validate_checksum "${path}"
  fi
}

validate_preparation() {
  local ds views run_root expected_root_real run_root_real union_path
  local prepare_script="${PROJECT_ROOT}/src/4_cell_type_annotation/1.1_prepare_chunks.py"
  [[ -r "${prepare_script}" ]] || return 1
  expected_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
  [[ -n "${expected_root_real}" ]] || return 1
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && -n "${views}" && -n "${run_root}" ]] || return 1
    run_root_real="$(ecoda_realpath_existing "${run_root}" 2>/dev/null || true)"
    [[ -n "${run_root_real}" && "${expected_root_real}" == "${run_root_real}" ]] || return 1
    export DS_NAME="${ds}" ANNOTATION_VIEWS="${views}" \
      ANNOTATION_RUN_ROOT="${run_root}" ANNOTATION_RUN_ID="${RUN_ID}"
    "${PYTHON_BIN}" "${prepare_script}" --views "${views}" \
      --run-root "${run_root}" --validate-only || return 1
    union_path="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    validate_artifact "${union_path}" "${RUN_ID}" || return 1
  done < "${ROOT_MANIFEST}"
}

validate_annotation() {
  local ds chunk feather_dir chunk_num feather union_path sample expected_union expected_chunk_dir
  local feather_sidecar_flag expected_args validated_unions=""
  while IFS=$'\t' read -r ds chunk feather_dir; do
    [[ -n "${ds}" && -n "${chunk}" && -n "${feather_dir}" ]] || return 1
    ecoda_validate_run_owned_path "${chunk}" "${RUN_ROOT}" || return 1
    expected_chunk_dir="${RUN_ROOT}/datasets/${ds}/chunks"
    expected_union="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    [[ "${chunk}" == "${expected_chunk_dir}/chunk_"*.txt ]] || return 1
    [[ "${feather_dir}" == "${RUN_ROOT}/datasets/${ds}/annotations" ]] || return 1
    chunk_num="${chunk##*/chunk_}"
    chunk_num="${chunk_num%.txt}"
    [[ "${chunk_num}" =~ ^[1-9][0-9]*$ ]] || return 1
    feather="${feather_dir}/annotations_chunk_${chunk_num}.feather"
    [[ -s "${feather}" ]] || return 1
    union_path="$(sed -n '1p' "${chunk}")"
    [[ "${union_path}" == "${expected_union}" ]] || return 1
    case " ${validated_unions} " in
      *" ${union_path} "*) ;;
      *)
        validate_artifact "${union_path}" "${RUN_ID}" || return 1
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
    [[ ${ARTIFACT_RECORD_USED} -eq 1 ]] && feather_sidecar_flag="--sidecar-validated"
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${feather}" "${feather_sidecar_flag}" "${expected_args[@]}" \
      >/dev/null 2>&1 || return 1
  done < "${ROOT_MANIFEST}"
}

marker_field() {
  local marker="$1" field="$2"
  sed -n "s/^${field}=//p" "${marker}" | head -1
}

validate_merge() {
  local ds views run_root marker view name path expected_sources expected_records
  local union_path union_md5 union_size source_md5 source_size record
  while IFS=$'\t' read -r ds views run_root; do
    [[ -n "${ds}" && "${run_root}" == "${RUN_ROOT}" ]] || return 1
    marker="${RUN_ROOT}/datasets/${ds}/merge.ok"
    [[ -s "${marker}" ]] || return 1
    [[ "$(marker_field "${marker}" STATE)" == "OK" &&
       "$(marker_field "${marker}" DATASET)" == "${ds}" &&
       "$(marker_field "${marker}" VIEWS)" == "${views}" ]] || return 1
    union_path="${RUN_ROOT}/datasets/${ds}/union/union.h5ad"
    [[ "$(marker_field "${marker}" UNION_PATH)" == "${union_path}" ]] || return 1
    validate_artifact "${union_path}" "${RUN_ID}" || return 1
    union_md5="${ECODA_CHECKSUM_MD5}"
    union_size="${ECODA_CHECKSUM_SIZE}"
    [[ "$(marker_field "${marker}" UNION_MD5)" == "${union_md5}" &&
       "$(marker_field "${marker}" UNION_SIZE)" == "${union_size}" ]] || return 1
    expected_sources=""
    expected_records=""
    IFS=',' read -r -a view_list <<< "${views}"
    for view in "${view_list[@]}"; do
      name="$(jq -r --arg ds "${ds}" --arg view "${view}" \
        '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' \
        "${DATASETS_JSON_FILE}")"
      [[ -n "${name}" ]] || return 1
      path="${HPC_SCRATCH_DIR}/${ds}/output/${name}"
      [[ -s "${path}" ]] || return 1
      validate_artifact "${path}" stage4_merge || return 1
      source_md5="${ECODA_CHECKSUM_MD5}"
      source_size="${ECODA_CHECKSUM_SIZE}"
      "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
        --h5ad "${path}" --sidecar-validated >/dev/null 2>&1 || return 1
      record="${path}|${source_md5}|${source_size}"
      [[ -z "${expected_sources}" ]] && expected_sources="${path}" ||
        expected_sources="${expected_sources};${path}"
      [[ -z "${expected_records}" ]] && expected_records="${record}" ||
        expected_records="${expected_records};${record}"
    done
    [[ "$(marker_field "${marker}" SOURCE_H5ADS)" == "${expected_sources}" &&
       "$(marker_field "${marker}" SOURCE_RECORDS)" == "${expected_records}" ]] || return 1
  done < "${ROOT_MANIFEST}"
}

case "${PHASE}" in
  preparation) validate_preparation ;;
  annotation) validate_annotation ;;
  merge) validate_merge ;;
esac
