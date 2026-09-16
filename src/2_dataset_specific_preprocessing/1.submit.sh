#!/bin/bash
# Generic immutable Stage 2 worker boundary. The gate owns the step table;
# this file validates the run once, enters the runtime, and dispatches only
# the scientific worker selected by that table.
set -euo pipefail

SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
SOURCE_SNAPSHOT_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
RUN_ROOT_FROM_STAGE="${STAGE2_RUN_ROOT:-}"
RUN_ROOT_FROM_ENV="${ECODA_RUN_ROOT:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUNTIME_IMAGE_ENV="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST_ENV="${ECODA_RUNTIME_MANIFEST:-}"
GENERIC_RELATIVE="src/2_dataset_specific_preprocessing/1.submit.sh"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ $# -eq 2 && "$1" == "--step" && -n "$2" ]] ||
  fail "usage: 1.submit.sh --step STEP"
STEP="$2"

[[ "${SOURCE_SNAPSHOT_REQUIRED}" == "1" ]] ||
  fail "Stage 2 worker requires an immutable source snapshot"
[[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* ]] ||
  fail "Stage 2 worker requires absolute source root and manifest"
[[ "${SOURCE_ROOT##*/}" == "tree" ]] ||
  fail "Stage 2 source root is not a snapshot tree: ${SOURCE_ROOT}"
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
[[ "${SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] ||
  fail "Stage 2 source manifest is not bound to the source snapshot"
[[ -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" ]] ||
  fail "Stage 2 immutable source root is missing or unsafe: ${SOURCE_ROOT}"
command -v realpath >/dev/null 2>&1 ||
  fail "realpath is required for Stage 2 source isolation"
SOURCE_ROOT_REAL="$(realpath -e "${SOURCE_ROOT}" 2>/dev/null || realpath "${SOURCE_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize Stage 2 source root: ${SOURCE_ROOT}"
[[ "${SOURCE_ROOT_REAL}" == "${SOURCE_ROOT}" ]] ||
  fail "Stage 2 source root is not canonical"
GENERIC_SCRIPT="${SOURCE_ROOT}/${GENERIC_RELATIVE}"
[[ -f "${GENERIC_SCRIPT}" && ! -L "${GENERIC_SCRIPT}" && -r "${GENERIC_SCRIPT}" ]] ||
  fail "generic Stage 2 worker boundary is missing from the source snapshot: ${GENERIC_SCRIPT}"

if [[ -n "${RUN_ROOT_FROM_STAGE}" ]]; then
  [[ -z "${RUN_ROOT_FROM_ENV}" || "${RUN_ROOT_FROM_STAGE}" == "${RUN_ROOT_FROM_ENV}" ]] ||
    fail "Stage 2 run-root exports disagree"
  RUN_ROOT="${RUN_ROOT_FROM_STAGE}"
else
  RUN_ROOT="${RUN_ROOT_FROM_ENV}"
fi
[[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] ||
  fail "Stage 2 worker requires an existing run root"
RUN_ROOT_REAL="$(realpath -e "${RUN_ROOT}" 2>/dev/null || realpath "${RUN_ROOT}" 2>/dev/null)" ||
  fail "could not canonicalize Stage 2 run root: ${RUN_ROOT}"
[[ "${RUN_ROOT_REAL}" == "${RUN_ROOT}" ]] ||
  fail "Stage 2 run root is not canonical"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
   "${RUN_ROOT##*/}" == "${RUN_ID}" ]] ||
  fail "Stage 2 run root and run ID do not match"
[[ -n "${RUNTIME_IMAGE_ENV}" && "${RUNTIME_IMAGE_ENV}" = /* &&
   -n "${RUNTIME_MANIFEST_ENV}" && "${RUNTIME_MANIFEST_ENV}" = /* ]] ||
  fail "Stage 2 worker requires bound runtime image and manifest paths"
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" ]] ||
  fail "Stage 2 run is missing its source/runtime identity manifests"

source "${SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh"
source "${SOURCE_ROOT}/src/slurm_config.sh"
cd "${PROJECT_ROOT}"

if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  [[ -f "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" && -r "${SOURCE_MANIFEST}" ]] ||
    fail "Stage 2 source manifest is missing or unsafe: ${SOURCE_MANIFEST}"
  cmp -s "${RUN_SOURCE_MANIFEST}" "${SOURCE_MANIFEST}" ||
    fail "Stage 2 run source manifest does not match the bound snapshot"
fi

STEP_ROW="$(BASH_ENV=/dev/null bash "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
  --step-config "${STEP}")" ||
  fail "unknown or malformed Stage 2 step: ${STEP}"
IFS='|' read -r ROW_STEP JOB_NAME STEP_TIME STEP_CPUS STEP_MEMORY EXECUTOR \
  WORKER_RELATIVE OUTPUTS DEPENDENCY <<< "${STEP_ROW}"
[[ "${ROW_STEP}" == "${STEP}" && -n "${JOB_NAME}" && -n "${STEP_TIME}" &&
   -n "${STEP_CPUS}" && -n "${STEP_MEMORY}" && -n "${EXECUTOR}" &&
   -n "${WORKER_RELATIVE}" && -n "${OUTPUTS}" && -n "${DEPENDENCY}" ]] ||
  fail "Stage 2 step table row is malformed: ${STEP}"
case "${WORKER_RELATIVE}" in
  ""|/*|*"/"*|*..*|*$'\n'*|*$'\t'*)
    fail "Stage 2 worker path is not a source-relative scientific worker: ${WORKER_RELATIVE}"
    ;;
esac
WORKER_SCRIPT="${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/${WORKER_RELATIVE}"
[[ -f "${WORKER_SCRIPT}" && ! -L "${WORKER_SCRIPT}" && -r "${WORKER_SCRIPT}" ]] ||
  fail "scientific Stage 2 worker is missing from the source snapshot: ${WORKER_SCRIPT}"

if [[ -n "${SLURM_JOB_ID:-}" && "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  command -v scontrol >/dev/null 2>&1 ||
    fail "scontrol is required to verify the immutable Stage 2 worker boundary"
  SCHEDULER_SCRIPT="$(scontrol show job "${SLURM_JOB_ID}" -o |
    grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)" ||
    fail "could not recover Stage 2 worker script from Slurm"
  [[ "${SCHEDULER_SCRIPT}" == "${GENERIC_SCRIPT}" ]] ||
    fail "Slurm worker command is not the immutable generic Stage 2 boundary"
fi

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
  export ECODA_RUNTIME_PROFILE=stage2
  ecoda_runtime_validate_bound_run || exit 1
fi

ecoda_runtime_reexec_worker stage2 "${GENERIC_SCRIPT}" --step "${STEP}" || exit 1

run_alzheimer() {
  local raw_input="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer.h5ad"
  local output_path="${HPC_SCRATCH_DIR}/Alzheimer/data/SEAAD_Alzheimer_donor_assay.h5ad"
  local -a worker_args
  [[ -f "${raw_input}" && ! -L "${raw_input}" && -r "${raw_input}" ]] ||
    fail "authoritative Alzheimer raw H5AD is missing or unsafe: ${raw_input}"
  [[ "${raw_input}" != "${output_path}" ]] ||
    fail "Alzheimer derivative would overwrite the raw input"
  [[ -d "$(dirname "${output_path}")" && ! -L "$(dirname "${output_path}")" ]] ||
    fail "Alzheimer output directory is missing or unsafe"
  [[ -x "${PYTHON_BIN}" ]] ||
    fail "pinned Stage 2 Python interpreter is not executable: ${PYTHON_BIN}"
  worker_args=(
    --input-file "${raw_input}"
    --output-file "${output_path}"
    --expected-samples 104
    --expected-donors 83
    --expected-assay-counts "10x3v3=83,10xmultiome=21"
    --expected-sex-counts "female=59,male=45"
    --require-example-ids
  )
  case "${FORCE_PREPROCESS:-0}" in
    0) ;;
    1) worker_args+=(--force) ;;
    *) fail "FORCE_PREPROCESS must be 0 or 1" ;;
  esac
  "${PYTHON_BIN}" "${WORKER_SCRIPT}" "${worker_args[@]}"
  "${PYTHON_BIN}" "${WORKER_SCRIPT}" \
    --input-file "${raw_input}" \
    --output-file "${output_path}" \
    --expected-samples 104 \
    --expected-donors 83 \
    --expected-assay-counts "10x3v3=83,10xmultiome=21" \
    --expected-sex-counts "female=59,male=45" \
    --require-example-ids \
    --validate-only
}

case "${EXECUTOR}" in
  python)
    "${PYTHON_BIN}" "${WORKER_SCRIPT}"
    ;;
  python_module)
    if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
      module load GCCcore/12.2.0
    fi
    "${PYTHON_BIN}" "${WORKER_SCRIPT}"
    ;;
  python_force)
    FORCE_FLAG=()
    case "${FORCE_PREPROCESS:-0}" in
      0) ;;
      1) FORCE_FLAG=(--force) ;;
      *) fail "FORCE_PREPROCESS must be 0 or 1" ;;
    esac
    "${PYTHON_BIN}" "${WORKER_SCRIPT}" "${FORCE_FLAG[@]}"
    ;;
  python_alzheimer)
    run_alzheimer
    ;;
  r)
    ${PIXI_RSCRIPT} "${WORKER_SCRIPT}"
    ;;
  *)
    fail "unsupported Stage 2 worker executor: ${EXECUTOR}"
    ;;
esac
