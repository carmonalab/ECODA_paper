#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
RUN_ROOT_FROM_ENV="${ECODA_RUN_ROOT:-}"
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
SOURCE_MANIFEST_RUN_ENV="${ECODA_SOURCE_MANIFEST_RUN:-}"
SOURCE_SNAPSHOT_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
RUN_ROOT_FROM_STAGE="${PREPROCESS_RUN_ROOT:-}"
RUN_ROOT_FROM_ENV="${RUN_ROOT_FROM_ENV:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUNTIME_IMAGE_ENV="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST_ENV="${ECODA_RUNTIME_MANIFEST:-}"
RUNTIME_IDENTITY_ENV="${ECODA_RUNTIME_IDENTITY:-}"

[[ "${SOURCE_SNAPSHOT_REQUIRED}" == "1" ]] || {
  echo "ERROR: Stage 3 worker requires an immutable source snapshot." >&2
  exit 1
}
[[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* ]] || {
  echo "ERROR: Stage 3 worker requires absolute source root and manifest." >&2
  exit 1
}
[[ "${SOURCE_ROOT##*/}" == "tree" ]] || {
  echo "ERROR: Stage 3 source root is not a snapshot tree: ${SOURCE_ROOT}" >&2
  exit 1
}
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
[[ "${SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] || {
  echo "ERROR: Stage 3 source manifest is not bound to the source snapshot." >&2
  exit 1
}
[[ -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" ]] || {
  echo "ERROR: Stage 3 immutable source root is missing or unsafe: ${SOURCE_ROOT}" >&2
  exit 1
}
WORKER_SCRIPT="${SOURCE_ROOT}/${SCRIPT_RELATIVE}"
[[ -f "${WORKER_SCRIPT}" && ! -L "${WORKER_SCRIPT}" && -r "${WORKER_SCRIPT}" ]] || {
  echo "ERROR: Stage 3 worker script is missing from the source snapshot: ${WORKER_SCRIPT}" >&2
  exit 1
}
WORKER_SCRIPT="$(ecoda_require_source_script_path "${WORKER_SCRIPT}" "${SOURCE_ROOT}")" || exit 1

if [[ -n "${RUN_ROOT_FROM_STAGE}" ]]; then
  [[ -z "${RUN_ROOT_FROM_ENV}" || "${RUN_ROOT_FROM_STAGE}" == "${RUN_ROOT_FROM_ENV}" ]] || {
    echo "ERROR: Stage 3 run-root exports disagree." >&2
    exit 1
  }
  RUN_ROOT="${RUN_ROOT_FROM_STAGE}"
else
  RUN_ROOT="${RUN_ROOT_FROM_ENV}"
fi
[[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 3 worker requires an existing run root." >&2
  exit 1
}
RUN_ROOT_REAL="$(realpath -e "${RUN_ROOT}" 2>/dev/null || realpath "${RUN_ROOT}" 2>/dev/null)" || {
  echo "ERROR: could not canonicalize Stage 3 run root: ${RUN_ROOT}" >&2
  exit 1
}
[[ "${RUN_ROOT_REAL}" == "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 3 run root is not canonical." >&2
  exit 1
}
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
   "${RUN_ROOT##*/}" == "${RUN_ID}" ]] || {
  echo "ERROR: Stage 3 run root and run ID do not match." >&2
  exit 1
}
[[ -n "${RUNTIME_IMAGE_ENV}" && "${RUNTIME_IMAGE_ENV}" = /* &&
   -n "${RUNTIME_MANIFEST_ENV}" && "${RUNTIME_MANIFEST_ENV}" = /* ]] || {
  echo "ERROR: Stage 3 worker requires bound runtime image and manifest paths." >&2
  exit 1
}
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" ]] || {
  echo "ERROR: Stage 3 run is missing its source/runtime identity manifests." >&2
  exit 1
}
ecoda_validate_run_owned_path "${RUN_SOURCE_MANIFEST}" "${RUN_ROOT}" || exit 1
ecoda_validate_run_owned_path "${RUN_RUNTIME_IDENTITY}" "${RUN_ROOT}" || exit 1
[[ -z "${SOURCE_MANIFEST_RUN_ENV}" || "${SOURCE_MANIFEST_RUN_ENV}" == "${RUN_SOURCE_MANIFEST}" ]] || {
  echo "ERROR: Stage 3 run source manifest export does not match the run root." >&2
  exit 1
}
[[ -z "${RUNTIME_IDENTITY_ENV}" || "${RUNTIME_IDENTITY_ENV}" == "${RUN_RUNTIME_IDENTITY}" ]] || {
  echo "ERROR: Stage 3 runtime identity export does not match the run root." >&2
  exit 1
}
cmp -s "${RUN_SOURCE_MANIFEST}" "${SOURCE_MANIFEST}" || {
  echo "ERROR: Stage 3 run source manifest does not match the immutable snapshot." >&2
  exit 1
}
export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID="${RUN_ID}"
IDENTITY_IMAGE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || exit 1
IDENTITY_MANIFEST="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || exit 1
[[ "${IDENTITY_IMAGE}" == "${RUNTIME_IMAGE_ENV}" &&
   "${IDENTITY_MANIFEST}" == "${RUNTIME_MANIFEST_ENV}" ]] || {
  echo "ERROR: Stage 3 run-bound runtime paths do not match the exported runtime." >&2
  exit 1
}
IDENTITY_IMAGE_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || exit 1
IDENTITY_MANIFEST_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || exit 1
IDENTITY_IMAGE_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || exit 1
IDENTITY_MANIFEST_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || exit 1
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  [[ "${ECODA_RUNTIME_IMAGE_SHA256:-}" == "${IDENTITY_IMAGE_SHA}" &&
     "${ECODA_RUNTIME_MANIFEST_SHA256:-}" == "${IDENTITY_MANIFEST_SHA}" &&
     "${ECODA_RUNTIME_IMAGE_SIZE:-}" == "${IDENTITY_IMAGE_SIZE}" &&
     "${ECODA_RUNTIME_MANIFEST_SIZE:-}" == "${IDENTITY_MANIFEST_SIZE}" ]] || {
    echo "ERROR: Stage 3 run-bound runtime identity does not match exported runtime metadata." >&2
    exit 1
  }
fi
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" &&
      "${ECODA_RUNTIME_MODE:-apptainer}" == "host" ]]; then
  ecoda_runtime_validate_bound_run || exit 1
fi
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_AUX_ROOT="${SOURCE_ROOT%/}/aux"
export ECODA_RUNTIME_IMAGE="${IDENTITY_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${IDENTITY_MANIFEST}"
export ECODA_RUNTIME_IDENTITY="${RUN_RUNTIME_IDENTITY}"
export ECODA_RUNTIME_PROFILE=stage3
ecoda_runtime_reexec_worker stage3 "${WORKER_SCRIPT}" || exit 1
cd "${PROJECT_ROOT}"

if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCHEDULER_SCRIPT="$(scontrol show job "${SLURM_JOB_ID}" -o |
    grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)" || {
    echo "ERROR: could not recover Stage 3 worker script from Slurm." >&2
    exit 1
  }
  [[ "${SCHEDULER_SCRIPT}" == "${WORKER_SCRIPT}" ]] || {
    echo "ERROR: Slurm worker command is not the immutable snapshot script." >&2
    exit 1
  }
fi

MANIFEST_PATH="${PREPROCESS_SELECTION_FILE:-${PREPROCESS_DATASETS_FILE:-}}"
DS_NAME=""
VIEW_NAME=""
if [[ -n "${MANIFEST_PATH}" ]]; then
  [[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: preprocessing manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
  ecoda_validate_run_owned_path "${MANIFEST_PATH}" "${RUN_ROOT}" || {
    echo "ERROR: Stage 3 preprocessing manifest is not run-owned." >&2
    exit 1
  }
  ecoda_validate_manifest "${MANIFEST_PATH}" 2 || {
    echo "ERROR: Stage 3 preprocessing manifest is malformed." >&2
    exit 1
  }
  manifest_line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
  [[ -n "${manifest_line}" ]] || { echo "ERROR: no manifest row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
  IFS=$'\t' read -r DS_NAME VIEW_NAME <<< "${manifest_line}"
else
  command -v jq >/dev/null 2>&1 || { echo "ERROR: jq unavailable and no preprocessing manifest" >&2; exit 1; }
  DS_NAME="$(jq -r 'keys[]' "${DATASETS_JSON_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
  VIEW_NAME="${PREPROCESS_VIEW:-}"
fi
[[ -n "${DS_NAME}" ]] || { echo "ERROR: empty dataset for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
export DS_NAME
export PREPROCESS_VIEW="${VIEW_NAME}"
echo "Processing dataset/view: ${DS_NAME}/${VIEW_NAME:-all} (array task ${SLURM_ARRAY_TASK_ID})"

DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"
OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"
FORCE_FLAG=()
[[ "${FORCE_PREPROCESS:-0}" == "1" ]] && FORCE_FLAG=(--force)
VIEW_FLAG=()
[[ -n "${VIEW_NAME}" ]] && VIEW_FLAG=(--view "${VIEW_NAME}")

source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_preprocess.py" \
  --config_path "${DATASETS_JSON_FILE}" \
  --input_dir "${DATA_DIR}" \
  --output_dir "${OUTPUT_DIR}" \
  --ds_name "${DS_NAME}" \
  "${FORCE_FLAG[@]}" \
  "${VIEW_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  if [[ -n "${VIEW_NAME}" ]]; then
    output_name="$(ecoda_view_output_name "${DS_NAME}" "${VIEW_NAME}")"
    [[ -n "${output_name}" ]] || {
      echo "ERROR: missing Stage 3 output contract for ${DS_NAME}/${VIEW_NAME}" >&2
      exit 1
    }
    output_path="${OUTPUT_DIR}/${output_name}"
    "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/benchmark_h5ad_contract.py" \
      --path "${output_path}" --view "${VIEW_NAME}" --method "Stage 3 worker publication" >/dev/null 2>&1 || {
      echo "ERROR: Stage 3 worker H5AD publication contract failed: ${output_path}" >&2
      exit 1
    }
    ecoda_write_artifact_record "${output_path}" stage3 "${RUN_ID}" >/dev/null || {
      echo "ERROR: Stage 3 worker artifact record publication failed: ${output_path}" >&2
      exit 1
    }
  fi
  worker_clear_retry_count
  echo "Preprocessing complete for ${DS_NAME}/${VIEW_NAME:-all}"
  exit 0
fi
ERR_PREFIX="${PREPROCESS_ERROR_PREFIX:-${LOGS_DIR}/3_scrnaseq_preprocessing}"
ERR_FILE="${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  if [[ -n "${VIEW_NAME}" ]] && command -v jq >/dev/null 2>&1; then
    output_name="$(jq -r --arg ds "${DS_NAME}" --arg view "${VIEW_NAME}" '.[$ds].views[$view].output_file_name // empty' "${DATASETS_JSON_FILE}")"
    [[ -n "${output_name}" ]] && rm -f "${OUTPUT_DIR}/${output_name}" "${OUTPUT_DIR}/${output_name}.md5"
  fi
  exit 0
fi
exit ${RC}
