#!/bin/bash
#SBATCH --job-name=alzheimer_donor_assay
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
SOURCE_SNAPSHOT_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
RUN_ROOT_FROM_STAGE="${STAGE2_RUN_ROOT:-}"
RUN_ROOT_FROM_ENV="${ECODA_RUN_ROOT:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUNTIME_IMAGE_ENV="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST_ENV="${ECODA_RUNTIME_MANIFEST:-}"
SCRIPT_RELATIVE="src/2_dataset_specific_preprocessing/1.7_submit_alzheimer_donor_assay.sh"
WORKER_RELATIVE="src/2_dataset_specific_preprocessing/1.7.1_create_alzheimer_donor_assay.py"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ $# -eq 0 ]] || fail "Alzheimer donor-assay hook does not accept positional arguments"
[[ "${SOURCE_SNAPSHOT_REQUIRED}" == "1" ]] ||
  fail "Alzheimer donor-assay hook requires an immutable source snapshot"
[[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* ]] ||
  fail "Alzheimer donor-assay hook requires absolute source root and manifest"
[[ "${SOURCE_ROOT##*/}" == "tree" ]] ||
  fail "Stage 2 source root is not a snapshot tree: ${SOURCE_ROOT}"
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
[[ "${SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] ||
  fail "Stage 2 source manifest is not bound to the source snapshot"
[[ -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" ]] ||
  fail "Stage 2 immutable source root is missing or unsafe: ${SOURCE_ROOT}"
command -v realpath >/dev/null 2>&1 || fail "realpath is required for Stage 2 source isolation"
SOURCE_ROOT_REAL="$(realpath -e "${SOURCE_ROOT}" 2>/dev/null || realpath "${SOURCE_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize Stage 2 source root: ${SOURCE_ROOT}"
[[ "${SOURCE_ROOT_REAL}" == "${SOURCE_ROOT}" ]] ||
  fail "Stage 2 source root is not canonical"
WORKER_SCRIPT="${SOURCE_ROOT}/${WORKER_RELATIVE}"
[[ -f "${WORKER_SCRIPT}" && ! -L "${WORKER_SCRIPT}" && -r "${WORKER_SCRIPT}" ]] ||
  fail "Alzheimer donor-assay worker is missing from the source snapshot: ${WORKER_SCRIPT}"

if [[ -n "${RUN_ROOT_FROM_STAGE}" ]]; then
  [[ -z "${RUN_ROOT_FROM_ENV}" || "${RUN_ROOT_FROM_STAGE}" == "${RUN_ROOT_FROM_ENV}" ]] ||
    fail "Stage 2 run-root exports disagree"
  RUN_ROOT="${RUN_ROOT_FROM_STAGE}"
else
  RUN_ROOT="${RUN_ROOT_FROM_ENV}"
fi
[[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] ||
  fail "Alzheimer donor-assay hook requires an existing run root"
RUN_ROOT_REAL="$(realpath -e "${RUN_ROOT}" 2>/dev/null || realpath "${RUN_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize Stage 2 run root: ${RUN_ROOT}"
[[ "${RUN_ROOT_REAL}" == "${RUN_ROOT}" ]] || fail "Stage 2 run root is not canonical"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
   "${RUN_ROOT##*/}" == "${RUN_ID}" ]] ||
  fail "Stage 2 run root and run ID do not match"
[[ -n "${RUNTIME_IMAGE_ENV}" && "${RUNTIME_IMAGE_ENV}" = /* &&
   -n "${RUNTIME_MANIFEST_ENV}" && "${RUNTIME_MANIFEST_ENV}" = /* ]] ||
  fail "Alzheimer donor-assay hook requires bound runtime image and manifest paths"
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" ]] ||
  fail "Stage 2 run is missing its source/runtime identity manifests"

if [[ -n "${SLURM_JOB_ID:-}" && "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  command -v scontrol >/dev/null 2>&1 || fail "scontrol is required to verify the immutable worker command"
  SCHEDULER_SCRIPT="$(scontrol show job "${SLURM_JOB_ID}" -o |
    grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)" ||
    fail "could not recover Stage 2 worker script from Slurm"
  [[ "${SCHEDULER_SCRIPT}" == "${SOURCE_ROOT}/${SCRIPT_RELATIVE}" ]] ||
    fail "Slurm worker command is not the immutable Alzheimer donor-assay hook"
fi

source "${SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  [[ -f "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" && -r "${SOURCE_MANIFEST}" ]] ||
    fail "Stage 2 source manifest is missing or unsafe: ${SOURCE_MANIFEST}"
  cmp -s "${RUN_SOURCE_MANIFEST}" "${SOURCE_MANIFEST}" ||
    fail "Stage 2 run source manifest does not match the bound snapshot"
fi
source "${SOURCE_ROOT}/src/slurm_config.sh"
export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID="${RUN_ID}"

IDENTITY_IMAGE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || exit 1
IDENTITY_MANIFEST="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || exit 1
[[ "${IDENTITY_IMAGE}" == "${RUNTIME_IMAGE_ENV}" &&
   "${IDENTITY_MANIFEST}" == "${RUNTIME_MANIFEST_ENV}" &&
   "${ECODA_RUNTIME_IMAGE}" == "${RUNTIME_IMAGE_ENV}" &&
   "${ECODA_RUNTIME_MANIFEST}" == "${RUNTIME_MANIFEST_ENV}" ]] ||
  fail "Stage 2 run-bound runtime paths do not match the exported runtime"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  IDENTITY_IMAGE_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || exit 1
  IDENTITY_MANIFEST_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || exit 1
  IDENTITY_IMAGE_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || exit 1
  IDENTITY_MANIFEST_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || exit 1
  [[ "${ECODA_RUNTIME_IMAGE_SHA256:-}" == "${IDENTITY_IMAGE_SHA}" &&
     "${ECODA_RUNTIME_MANIFEST_SHA256:-}" == "${IDENTITY_MANIFEST_SHA}" &&
     "${ECODA_RUNTIME_IMAGE_SIZE:-}" == "${IDENTITY_IMAGE_SIZE}" &&
     "${ECODA_RUNTIME_MANIFEST_SIZE:-}" == "${IDENTITY_MANIFEST_SIZE}" ]] ||
    fail "Stage 2 run-bound runtime identity does not match exported runtime metadata"
  ecoda_runtime_validate_bound_run || exit 1
fi

ecoda_runtime_reexec_worker stage2 "${SOURCE_ROOT}/${SCRIPT_RELATIVE}" || exit 1
cd "${PROJECT_ROOT}"

RAW_INPUT="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer.h5ad"
OUTPUT_PATH="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
[[ -f "${RAW_INPUT}" && ! -L "${RAW_INPUT}" && -r "${RAW_INPUT}" ]] ||
  fail "authoritative Alzheimer raw H5AD is missing or unsafe: ${RAW_INPUT}"
[[ "${RAW_INPUT}" != "${OUTPUT_PATH}" ]] || fail "Alzheimer derivative would overwrite the raw input"
[[ -d "$(dirname "${OUTPUT_PATH}")" && ! -L "$(dirname "${OUTPUT_PATH}")" ]] ||
  fail "Alzheimer output directory is missing or unsafe"
[[ -x "${PYTHON_BIN}" ]] || fail "pinned Stage 2 Python interpreter is not executable: ${PYTHON_BIN}"

WORKER_ARGS=(
  --input-file "${RAW_INPUT}"
  --output-file "${OUTPUT_PATH}"
  --expected-samples 104
  --expected-donors 83
  --expected-assay-counts "10x3v3=83,10xmultiome=21"
  --expected-sex-counts "female=59,male=45"
  --require-example-ids
)
case "${FORCE_PREPROCESS:-0}" in
  0) ;;
  1) WORKER_ARGS+=(--force) ;;
  *) fail "FORCE_PREPROCESS must be 0 or 1" ;;
esac

"${PYTHON_BIN}" "${WORKER_SCRIPT}" "${WORKER_ARGS[@]}"
"${PYTHON_BIN}" "${WORKER_SCRIPT}" \
  --input-file "${RAW_INPUT}" \
  --output-file "${OUTPUT_PATH}" \
  --expected-samples 104 \
  --expected-donors 83 \
  --expected-assay-counts "10x3v3=83,10xmultiome=21" \
  --expected-sex-counts "female=59,male=45" \
  --require-example-ids \
  --validate-only
