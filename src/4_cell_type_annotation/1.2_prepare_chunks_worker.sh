#!/bin/bash
#SBATCH --job-name=annotation_prepare
#SBATCH --time=01:00:00
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
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}" ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
SNAPSHOT_ID="${SNAPSHOT_ROOT##*/}"
[[ "${SOURCE_ROOT##*/}" == "tree" &&
   "${SNAPSHOT_ID}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] ||
  { echo "ERROR: source paths are not a commit-keyed snapshot identity." >&2; exit 1; }
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
  "${SOURCE_ROOT}/src/4_cell_type_annotation/1.2_prepare_chunks_worker.sh" || exit 1
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

MANIFEST="${ANNOTATION_PREP_MANIFEST:-}"
[[ -r "${MANIFEST}" ]] || { echo "ERROR: annotation prep manifest is unreadable: ${MANIFEST}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
EXPECTED_ROOT="${BOUND_RUN_ROOT}"
ecoda_validate_run_owned_path "${MANIFEST}" "${EXPECTED_ROOT}" ||
  { echo "ERROR: annotation prep manifest is outside the Stage 4 run root." >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
[[ -n "${line}" ]] || { echo "ERROR: no annotation prep row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME VIEWS RUN_ROOT <<< "${line}"
[[ -n "${DS_NAME}" && -n "${VIEWS}" && -n "${RUN_ROOT}" ]] || { echo "ERROR: malformed annotation prep row" >&2; exit 1; }
expected_root_real="$(ecoda_realpath_existing "${EXPECTED_ROOT}" 2>/dev/null || true)"
run_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
[[ -n "${expected_root_real}" && "${expected_root_real}" == "${run_root_real}" ]] ||
  { echo "ERROR: annotation prep row run root is not canonical." >&2; exit 1; }
export DS_NAME ANNOTATION_VIEWS="${VIEWS}" ANNOTATION_RUN_ROOT="${RUN_ROOT}" ANNOTATION_RUN_ID="${RUN_ID}"
FORCE_FLAG=()
[[ "${FORCE_ANNOTATION:-0}" == "1" ]] && FORCE_FLAG=(--force)

source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" \
  --views "${VIEWS}" --run-root "${RUN_ROOT}" "${FORCE_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  exit 0
fi
ERR_FILE="${ANNOTATION_PREP_ERROR_PREFIX:-${LOGS_DIR}/4_annotation_prepare}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
