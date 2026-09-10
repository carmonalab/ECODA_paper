#!/bin/bash
#SBATCH --job-name=scrna_worker
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
  "${SOURCE_ROOT}/src/4_cell_type_annotation/2.1_run_worker.sh" || exit 1
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

MANIFEST_PATH="${CHUNKS_MANIFEST:-}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: CHUNKS_MANIFEST is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
RUN_ROOT="${BOUND_RUN_ROOT}"
ecoda_validate_run_owned_path "${MANIFEST_PATH}" "${RUN_ROOT}" ||
  { echo "ERROR: annotation chunk manifest is outside the Stage 4 run root." >&2; exit 1; }
MANIFEST_LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
[[ -n "${MANIFEST_LINE}" ]] || { echo "ERROR: no chunk manifest row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME CHUNK_FILE FEATHER_DIR <<< "${MANIFEST_LINE}"
[[ -n "${DS_NAME}" && -n "${CHUNK_FILE}" ]] || { echo "ERROR: malformed chunk manifest row" >&2; exit 1; }
export DS_NAME
EXPECTED_FEATHER_DIR="${RUN_ROOT}/datasets/${DS_NAME}/annotations"
[[ "${FEATHER_DIR:-${EXPECTED_FEATHER_DIR}}" == "${EXPECTED_FEATHER_DIR}" ]] ||
  { echo "ERROR: annotation feather directory is not run-owned." >&2; exit 1; }
ecoda_validate_run_owned_path "${CHUNK_FILE}" "${RUN_ROOT}" ||
  { echo "ERROR: annotation chunk is outside the Stage 4 run root." >&2; exit 1; }
FEATHER_DIR="${EXPECTED_FEATHER_DIR}"
export ANNOTATION_OUTPUT_DIR="${FEATHER_DIR}"
NORMAL_TISSUE="$(jq -r --arg ds "${DS_NAME}" '
  if has($ds) and (.[$ds] | has("normal_tissue")) then
    if .[$ds].normal_tissue == true then "true"
    elif .[$ds].normal_tissue == false then "false"
    else empty end
  else empty end
' "${DATASETS_JSON_FILE}")"
[[ "${NORMAL_TISSUE}" == "true" || "${NORMAL_TISSUE}" == "false" ]] ||
  { echo "ERROR: configured normal_tissue is missing or invalid for ${DS_NAME}." >&2; exit 1; }
export NORMAL_TISSUE
[[ -f "${CHUNK_FILE}" ]] || { echo "ERROR: chunk file not found: ${CHUNK_FILE}" >&2; exit 1; }

CHUNK_NUM="$(basename "${CHUNK_FILE}")"
CHUNK_NUM="${CHUNK_NUM#chunk_}"
CHUNK_NUM="${CHUNK_NUM%.txt}"
FEATHER_FILE="${FEATHER_DIR}/annotations_chunk_${CHUNK_NUM}.feather"
remove_feather_record() {
  local record
  record="$(ecoda_artifact_record_path "${FEATHER_FILE}" "${RUN_ID}" 2>/dev/null || true)"
  [[ -z "${record}" ]] || rm -f "${record}"
}
mkdir -p "${FEATHER_DIR}"
if [[ -s "${FEATHER_FILE}" ]]; then
  if "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${FEATHER_FILE}" --require-sidecar >/dev/null 2>&1; then
    if ! ecoda_validate_artifact_record "${FEATHER_FILE}" "stage4_annotation" "${RUN_ID}" >/dev/null 2>&1; then
      ecoda_write_artifact_record "${FEATHER_FILE}" "stage4_annotation" "${RUN_ID}" >/dev/null || {
        echo "ERROR: annotation feather artifact record is invalid: ${FEATHER_FILE}" >&2
        exit 1
      }
    fi
    echo "Annotation feather already exists and passed schema/checksum validation: ${FEATHER_FILE}"
    exit 0
  fi
  echo "Existing annotation feather failed schema/checksum validation; rebuilding: ${FEATHER_FILE}" >&2
  remove_feather_record
  rm -f "${FEATHER_FILE}" "${FEATHER_FILE}.md5"
fi

source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
export_worker_thread_env
set +e
${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.1.1_process_chunk.R" "${CHUNK_FILE}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  [[ -s "${FEATHER_FILE}" ]] || { echo "ERROR: annotation worker exited without a feather: ${FEATHER_FILE}" >&2; exit 1; }
  # The canonical annotation contract performs the strict sidecar check once.
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
    --path "${FEATHER_FILE}" --require-sidecar >/dev/null 2>&1 || {
    echo "ERROR: annotation worker produced an invalid feather: ${FEATHER_FILE}" >&2
    exit 1
  }
  ecoda_write_artifact_record "${FEATHER_FILE}" "stage4_annotation" "${RUN_ID}" >/dev/null || {
    echo "ERROR: annotation feather artifact record publication failed: ${FEATHER_FILE}" >&2
    exit 1
  }
  exit 0
fi
ERR_FILE="${ANNOTATION_ERROR_PREFIX:-${LOGS_DIR}/4_cell_type_annotation}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  rm -f "${FEATHER_FILE}"
  remove_feather_record
  exit 0
fi
exit ${RC}
