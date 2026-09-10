#!/bin/bash
#SBATCH --job-name=myocardial_prep
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

# ---------------------------------------------------------------------------
# Myocardial Infarction raw count reconstruction step.
#
# Runs 1.5.1_reconstruct_myocardial_counts.py, which inverts the log1p-normalized
# expression in Myocardial_Infarc_2.h5ad into exact raw UMI integer counts via
# cell-wise minimum step inversion and vaults them into adata.layers["counts"].
#
# Runs 1.5.1_reconstruct_myocardial_counts.py; FORCE_PREPROCESS=1 is translated
# to --force so intentional recomputation survives the scheduler hook.
set -euo pipefail

SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
SOURCE_SNAPSHOT_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
RUN_ROOT_FROM_STAGE="${STAGE2_RUN_ROOT:-}"
RUN_ROOT_FROM_ENV="${ECODA_RUN_ROOT:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUNTIME_IMAGE_ENV="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST_ENV="${ECODA_RUNTIME_MANIFEST:-}"
SCRIPT_RELATIVE="src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh"

[[ "${SOURCE_SNAPSHOT_REQUIRED}" == "1" ]] || {
  echo "ERROR: Stage 2 worker requires an immutable source snapshot." >&2
  exit 1
}
[[ "${SOURCE_ROOT}" = /* && "${SOURCE_MANIFEST}" = /* ]] || {
  echo "ERROR: Stage 2 worker requires absolute source root and manifest." >&2
  exit 1
}
[[ "${SOURCE_ROOT##*/}" == "tree" ]] || {
  echo "ERROR: Stage 2 source root is not a snapshot tree: ${SOURCE_ROOT}" >&2
  exit 1
}
SNAPSHOT_ROOT="${SOURCE_ROOT%/tree}"
[[ "${SNAPSHOT_ROOT##*/}" =~ ^[[:xdigit:]]{40}$ &&
   "${SOURCE_MANIFEST}" == "${SNAPSHOT_ROOT}/identity/source.manifest" ]] || {
  echo "ERROR: Stage 2 source manifest is not bound to the source snapshot." >&2
  exit 1
}
[[ -d "${SOURCE_ROOT}" && ! -L "${SOURCE_ROOT}" ]] || {
  echo "ERROR: Stage 2 immutable source root is missing or unsafe: ${SOURCE_ROOT}" >&2
  exit 1
}
command -v realpath >/dev/null 2>&1 || {
  echo "ERROR: realpath is required for Stage 2 source isolation." >&2
  exit 1
}
SOURCE_ROOT_REAL="$(realpath -e "${SOURCE_ROOT}" 2>/dev/null || realpath "${SOURCE_ROOT}" 2>/dev/null)" || {
  echo "ERROR: could not canonicalize Stage 2 source root: ${SOURCE_ROOT}" >&2
  exit 1
}
[[ "${SOURCE_ROOT_REAL}" == "${SOURCE_ROOT}" ]] || {
  echo "ERROR: Stage 2 source root is not canonical." >&2
  exit 1
}
WORKER_SCRIPT="${SOURCE_ROOT}/${SCRIPT_RELATIVE}"
[[ -f "${WORKER_SCRIPT}" && ! -L "${WORKER_SCRIPT}" && -r "${WORKER_SCRIPT}" ]] || {
  echo "ERROR: Stage 2 worker script is missing from the source snapshot: ${WORKER_SCRIPT}" >&2
  exit 1
}

if [[ -n "${RUN_ROOT_FROM_STAGE}" ]]; then
  [[ -z "${RUN_ROOT_FROM_ENV}" || "${RUN_ROOT_FROM_STAGE}" == "${RUN_ROOT_FROM_ENV}" ]] || {
    echo "ERROR: Stage 2 run-root exports disagree." >&2
    exit 1
  }
  RUN_ROOT="${RUN_ROOT_FROM_STAGE}"
else
  RUN_ROOT="${RUN_ROOT_FROM_ENV}"
fi
[[ "${RUN_ROOT}" = /* && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 2 worker requires an existing run root." >&2
  exit 1
}
RUN_ROOT_REAL="$(realpath -e "${RUN_ROOT}" 2>/dev/null || realpath "${RUN_ROOT}" 2>/dev/null)" || {
  echo "ERROR: could not canonicalize Stage 2 run root: ${RUN_ROOT}" >&2
  exit 1
}
[[ "${RUN_ROOT_REAL}" == "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 2 run root is not canonical." >&2
  exit 1
}
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ &&
   "${RUN_ROOT##*/}" == "${RUN_ID}" ]] || {
  echo "ERROR: Stage 2 run root and run ID do not match." >&2
  exit 1
}
[[ -n "${RUNTIME_IMAGE_ENV}" && "${RUNTIME_IMAGE_ENV}" = /* &&
   -n "${RUNTIME_MANIFEST_ENV}" && "${RUNTIME_MANIFEST_ENV}" = /* ]] || {
  echo "ERROR: Stage 2 worker requires bound runtime image and manifest paths." >&2
  exit 1
}
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" ]] || {
  echo "ERROR: Stage 2 run is missing its source/runtime identity manifests." >&2
  exit 1
}

if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCHEDULER_SCRIPT="$(scontrol show job "${SLURM_JOB_ID}" -o |
    grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)" || {
    echo "ERROR: could not recover Stage 2 worker script from Slurm." >&2
    exit 1
  }
  [[ "${SCHEDULER_SCRIPT}" == "${WORKER_SCRIPT}" ]] || {
    echo "ERROR: Slurm worker command is not the immutable snapshot script." >&2
    exit 1
  }
fi
source "${SOURCE_ROOT}/src/utils/bash/ecoda_runtime.sh"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  [[ -f "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" && -r "${SOURCE_MANIFEST}" ]] || {
    echo "ERROR: Stage 2 source manifest is missing or unsafe: ${SOURCE_MANIFEST}" >&2
    exit 1
  }
  cmp -s "${RUN_SOURCE_MANIFEST}" "${SOURCE_MANIFEST}" || {
    echo "ERROR: Stage 2 run source manifest does not match the bound snapshot." >&2
    exit 1
  }
fi
source "${SOURCE_ROOT}/src/slurm_config.sh"
export ECODA_RUN_ROOT="${RUN_ROOT}" ECODA_RUN_ID="${RUN_ID}"

IDENTITY_IMAGE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE)" || exit 1
IDENTITY_MANIFEST="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST)" || exit 1
[[ "${IDENTITY_IMAGE}" == "${RUNTIME_IMAGE_ENV}" &&
   "${IDENTITY_MANIFEST}" == "${RUNTIME_MANIFEST_ENV}" &&
   "${ECODA_RUNTIME_IMAGE}" == "${RUNTIME_IMAGE_ENV}" &&
   "${ECODA_RUNTIME_MANIFEST}" == "${RUNTIME_MANIFEST_ENV}" ]] || {
  echo "ERROR: Stage 2 run-bound runtime paths do not match the exported runtime." >&2
  exit 1
}
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  IDENTITY_IMAGE_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SHA256)" || exit 1
  IDENTITY_MANIFEST_SHA="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SHA256)" || exit 1
  IDENTITY_IMAGE_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_IMAGE_SIZE)" || exit 1
  IDENTITY_MANIFEST_SIZE="$(_ecoda_runtime_require_identity_value "${RUN_RUNTIME_IDENTITY}" RUNTIME_MANIFEST_SIZE)" || exit 1
  [[ "${ECODA_RUNTIME_IMAGE_SHA256:-}" == "${IDENTITY_IMAGE_SHA}" &&
     "${ECODA_RUNTIME_MANIFEST_SHA256:-}" == "${IDENTITY_MANIFEST_SHA}" &&
     "${ECODA_RUNTIME_IMAGE_SIZE:-}" == "${IDENTITY_IMAGE_SIZE}" &&
     "${ECODA_RUNTIME_MANIFEST_SIZE:-}" == "${IDENTITY_MANIFEST_SIZE}" ]] || {
    echo "ERROR: Stage 2 run-bound runtime identity does not match exported runtime metadata." >&2
    exit 1
  }
  ecoda_runtime_validate_bound_run || exit 1
fi
ecoda_runtime_reexec_worker stage2 "${WORKER_SCRIPT}" "$@" || exit 1
cd "${PROJECT_ROOT}"


FORCE_FLAG=()
[[ "${FORCE_PREPROCESS:-0}" == "1" ]] && FORCE_FLAG=(--force)
"${PYTHON_BIN}" "${SOURCE_ROOT}/src/2_dataset_specific_preprocessing/1.5.1_reconstruct_myocardial_counts.py" "${FORCE_FLAG[@]}" "$@"
