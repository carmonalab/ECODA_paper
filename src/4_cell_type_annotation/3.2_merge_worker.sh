#!/bin/bash
#SBATCH --job-name=annotation_merge
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == "1" ]] || {
  echo "ERROR: Stage 4 workers require an immutable source snapshot." >&2
  exit 1
}
SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-}"
RUN_ID="${ECODA_RUN_ID:-}"
ANNOTATION_ID="${ANNOTATION_RUN_ID:-}"
[[ -z "${ANNOTATION_ID}" || "${ANNOTATION_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] ||
  { echo "ERROR: ANNOTATION_RUN_ID is invalid." >&2; exit 1; }
[[ -z "${ANNOTATION_ID}" || "${ANNOTATION_ID}" == "${RUN_ID}" ]] ||
  { echo "ERROR: ANNOTATION_RUN_ID does not match ECODA_RUN_ID." >&2; exit 1; }
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
  [[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* &&
     "${RUNTIME_IMAGE}" = /* && "${RUNTIME_MANIFEST}" = /* ]] ||
    { echo "ERROR: immutable snapshot/runtime paths are missing in the container." >&2; exit 1; }
else
  [[ "${SOURCE_ROOT}" = /* && -d "${SOURCE_ROOT}" ]] ||
    { echo "ERROR: immutable source root is missing or invalid." >&2; exit 1; }
  [[ "${SOURCE_MANIFEST}" = /* && -f "${SOURCE_MANIFEST}" &&
     ! -L "${SOURCE_MANIFEST}" && -r "${SOURCE_MANIFEST}" ]] ||
    { echo "ERROR: immutable source manifest is missing or invalid." >&2; exit 1; }
  [[ "${RUNTIME_IMAGE}" = /* && -f "${RUNTIME_IMAGE}" &&
     ! -L "${RUNTIME_IMAGE}" && -r "${RUNTIME_IMAGE}" ]] ||
    { echo "ERROR: immutable runtime image is missing or invalid." >&2; exit 1; }
  [[ "${RUNTIME_MANIFEST}" = /* && -f "${RUNTIME_MANIFEST}" &&
     ! -L "${RUNTIME_MANIFEST}" && -r "${RUNTIME_MANIFEST}" ]] ||
    { echo "ERROR: immutable runtime manifest is missing or invalid." >&2; exit 1; }
  SOURCE_ROOT="$(realpath "${SOURCE_ROOT}" 2>/dev/null || true)"
  [[ -n "${SOURCE_ROOT}" ]] ||
    { echo "ERROR: immutable source root cannot be canonicalized." >&2; exit 1; }
  SOURCE_MANIFEST="$(realpath "${SOURCE_MANIFEST}" 2>/dev/null || true)"
  [[ -n "${SOURCE_MANIFEST}" ]] ||
    { echo "ERROR: immutable source manifest cannot be canonicalized." >&2; exit 1; }
fi
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
SNAPSHOT_ID="${SNAPSHOT_ROOT##*/}"
[[ "${SOURCE_ROOT##*/}" == "tree" &&
   "${SNAPSHOT_ID}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] ||
  { echo "ERROR: source paths are not a commit-keyed snapshot identity." >&2; exit 1; }
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
RUN_ROOT="${ECODA_RUN_ROOT:-}"
if [[ -n "${RUN_ROOT}" ]]; then
  [[ "${RUN_ROOT}" = /* ]] ||
    { echo "ERROR: ECODA_RUN_ROOT must be absolute." >&2; exit 1; }
else
  RUNS_ROOT="${ECODA_RUNS_ROOT:-}"
  [[ "${RUNS_ROOT}" = /* ]] ||
    { echo "ERROR: ECODA_RUNS_ROOT is required when ECODA_RUN_ROOT is unset." >&2; exit 1; }
  if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
    [[ -d "${RUNS_ROOT}" ]] ||
      { echo "ERROR: ECODA_RUNS_ROOT does not exist." >&2; exit 1; }
  fi
  RUN_ROOT="${RUNS_ROOT}/${RUN_ID}"
fi
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" ]]; then
  [[ "${RUN_ROOT##*/}" == "${RUN_ID}" ]] ||
    { echo "ERROR: container run root does not match ECODA_RUN_ID." >&2; exit 1; }
else
  RUN_ROOT="$(realpath "${RUN_ROOT}" 2>/dev/null || true)"
  [[ -n "${RUN_ROOT}" && -d "${RUN_ROOT}" && "${RUN_ROOT##*/}" == "${RUN_ID}" ]] ||
    { echo "ERROR: run root is missing or does not match ECODA_RUN_ID." >&2; exit 1; }
  SOURCE_COPY="${RUN_ROOT}/manifests/source.manifest"
  RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
  [[ -f "${SOURCE_COPY}" && ! -L "${SOURCE_COPY}" && -r "${SOURCE_COPY}" &&
     -f "${RUNTIME_IDENTITY}" && ! -L "${RUNTIME_IDENTITY}" && -r "${RUNTIME_IDENTITY}" ]] ||
    { echo "ERROR: run-bound source/runtime manifests are missing or invalid." >&2; exit 1; }
  cmp -s "${SOURCE_MANIFEST}" "${SOURCE_COPY}" ||
    { echo "ERROR: run-bound source manifest does not match ECODA_SOURCE_MANIFEST." >&2; exit 1; }
  identity_image="$(sed -n 's/^RUNTIME_IMAGE=//p' "${RUNTIME_IDENTITY}" | head -1)"
  identity_manifest="$(sed -n 's/^RUNTIME_MANIFEST=//p' "${RUNTIME_IDENTITY}" | head -1)"
  [[ "${identity_image}" == "${RUNTIME_IMAGE}" &&
     "${identity_manifest}" == "${RUNTIME_MANIFEST}" ]] ||
    { echo "ERROR: run-bound runtime identity does not match runtime environment." >&2; exit 1; }
fi
BOUND_RUN_ROOT="${RUN_ROOT}"
export ECODA_RUN_ID="${RUN_ID}" ECODA_RUN_ROOT="${RUN_ROOT}"
export ECODA_RUNTIME_PROFILE=stage4
[[ -r "${SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh" ]] ||
  { echo "ERROR: snapshot runtime helper is missing." >&2; exit 1; }
source "${SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 \
  "${SOURCE_ROOT}/src/4_cell_type_annotation/3.2_merge_worker.sh" || exit 1
SCRIPT_DIR="${SOURCE_ROOT}/src/4_cell_type_annotation"
[[ -r "${SCRIPT_DIR}/../slurm_config.sh" &&
   -r "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh" &&
   -r "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh" ]] ||
  { echo "ERROR: snapshot runtime helpers are missing." >&2; exit 1; }
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  ecoda_runtime_validate_bound_run || exit 1
fi
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"
ECODA_RUNS_ROOT="${BOUND_RUN_ROOT%/*}"

MANIFEST="${ANNOTATION_MERGE_MANIFEST:-}"
[[ -r "${MANIFEST}" ]] || { echo "ERROR: merge manifest is unreadable: ${MANIFEST}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
RUN_ROOT="${BOUND_RUN_ROOT}"
ecoda_validate_run_owned_path "${MANIFEST}" "${RUN_ROOT}" ||
  { echo "ERROR: merge manifest is outside the Stage 4 run root." >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
[[ -n "${line}" ]] || { echo "ERROR: no merge row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME VIEWS MANIFEST_RUN_ROOT <<< "${line}"
[[ -n "${DS_NAME}" && -n "${VIEWS}" && -n "${MANIFEST_RUN_ROOT}" ]] || { echo "ERROR: malformed merge manifest row" >&2; exit 1; }
run_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
manifest_root_real="$(ecoda_realpath_existing "${MANIFEST_RUN_ROOT}" 2>/dev/null || true)"
[[ -n "${run_root_real}" && "${run_root_real}" == "${manifest_root_real}" ]] ||
  { echo "ERROR: merge manifest row run root is not canonical." >&2; exit 1; }
RUN_ROOT="${run_root_real}"
[[ -r "${RUN_ROOT}/metadata" ]] || { echo "ERROR: Stage 4 run metadata is missing." >&2; exit 1; }
metadata_stage="$(sed -n 's/^STAGE=//p' "${RUN_ROOT}/metadata" | head -1)"
metadata_run="$(sed -n 's/^RUN_ID=//p' "${RUN_ROOT}/metadata" | head -1)"
[[ "${metadata_stage}" == "stage4" && "${metadata_run}" == "${RUN_ID}" ]] ||
  { echo "ERROR: Stage 4 run metadata does not match merge run." >&2; exit 1; }
OWNER_DIR="$(ecoda_owner_dir stage4 "${DS_NAME}")"
owner_state="$(ecoda_owner_state "${OWNER_DIR}" 2>/dev/null || true)"
owner_run="$(ecoda_owner_run "${OWNER_DIR}" 2>/dev/null || true)"
owner_stage="$(ecoda_owner_field "${OWNER_DIR}" STAGE 2>/dev/null || true)"
owner_key="$(ecoda_owner_field "${OWNER_DIR}" KEY 2>/dev/null || true)"
[[ "${owner_run}" == "${RUN_ID}" && "${owner_stage}" == "stage4" &&
   "${owner_key}" == "${DS_NAME}" &&
   ( "${owner_state}" == "ACTIVE" || "${owner_state}" == "OK" ) ]] ||
  { echo "ERROR: Stage 4 dataset owner is missing or foreign." >&2; exit 1; }
ANNOT_DIR="${RUN_ROOT}/datasets/${DS_NAME}/annotations"
[[ -d "${ANNOT_DIR}" ]] || { echo "ERROR: annotation feather directory missing: ${ANNOT_DIR}" >&2; exit 1; }
IFS=',' read -r -a VIEW_LIST <<< "${VIEWS}"
[[ ${#VIEW_LIST[@]} -gt 0 ]] || { echo "ERROR: merge view list is empty" >&2; exit 1; }
UNION_PATH="${RUN_ROOT}/datasets/${DS_NAME}/union/union.h5ad"
ecoda_validate_checksum "${UNION_PATH}" || { echo "ERROR: annotation union checksum failed: ${UNION_PATH}" >&2; exit 1; }
SOURCE_PATHS=""
SOURCE_RECORDS=""
for view in "${VIEW_LIST[@]}"; do
  output_name="$(jq -r --arg ds "${DS_NAME}" --arg view "${view}" '.[$ds].views[$view].output_file_name // .[$ds].views[$view].output_file // empty' "${DATASETS_JSON_FILE}")"
  [[ -n "${output_name}" ]] || { echo "ERROR: missing output contract for ${DS_NAME}/${view}" >&2; exit 1; }
  h5ad_path="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${output_name}"
  [[ -s "${h5ad_path}" ]] || { echo "ERROR: h5ad missing before merge: ${h5ad_path}" >&2; exit 1; }
  ecoda_validate_checksum "${h5ad_path}" || { echo "ERROR: h5ad checksum failed before merge: ${h5ad_path}" >&2; exit 1; }
  "${PYTHON_BIN}" "${SCRIPT_DIR}/3.1_merge_annotations.py" \
    --h5ad-path "${h5ad_path}" --annot-dir "${ANNOT_DIR}" --union-path "${UNION_PATH}" \
    --output-path "${h5ad_path}"
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
    --h5ad "${h5ad_path}" --require-sidecar >/dev/null
  ecoda_write_checksum "${h5ad_path}"
  ecoda_write_artifact_record "${h5ad_path}" "stage4_merge" "${RUN_ID}" >/dev/null || {
    echo "ERROR: merged H5AD artifact record publication failed: ${h5ad_path}" >&2
    exit 1
  }
  source_record="${h5ad_path}|${ECODA_CHECKSUM_MD5}|${ECODA_CHECKSUM_SIZE}"
  [[ -z "${SOURCE_PATHS}" ]] && SOURCE_PATHS="${h5ad_path}" ||
    SOURCE_PATHS="${SOURCE_PATHS};${h5ad_path}"
  [[ -z "${SOURCE_RECORDS}" ]] && SOURCE_RECORDS="${source_record}" ||
    SOURCE_RECORDS="${SOURCE_RECORDS};${source_record}"
done
union_md5="$(ecoda_md5_file "${UNION_PATH}")"
union_size="$(wc -c < "${UNION_PATH}" | tr -d '[:space:]')"
marker="${RUN_ROOT}/datasets/${DS_NAME}/merge.ok"
tmp="${marker}.tmp.$$"
printf 'STATE=OK\nDATASET=%s\nVIEWS=%s\nSOURCE_H5ADS=%s\nSOURCE_RECORDS=%s\nUNION_PATH=%s\nUNION_MD5=%s\nUNION_SIZE=%s\n' \
  "${DS_NAME}" "${VIEWS}" "${SOURCE_PATHS}" "${SOURCE_RECORDS}" \
  "${UNION_PATH}" "${union_md5}" "${union_size}" > "${tmp}"
mv -f "${tmp}" "${marker}"
