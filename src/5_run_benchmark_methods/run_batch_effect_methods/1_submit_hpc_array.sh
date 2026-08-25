#!/bin/bash
set -euo pipefail

# Submit the explicit two-pass batch-effect method suite.
#
# Usage:
#   ./1_submit_hpc_array.sh --pass uncorrected
#   ./1_submit_hpc_array.sh --pass uncorrected --ds_name _debug \
#       --methods prepare_pseudobulk,composition,mrvi
#   ./1_submit_hpc_array.sh --pass corrected --ds_name Joanito
#   ./1_submit_hpc_array.sh --pass uncorrected --sync-only 12345,12346
#
# Only the declared batch-effect views and method allow-list are accepted.
# Every pass owns a separate scratch/NAS root, manifest set, execution-log
# prefix, cache namespace, and result stem. This submitter intentionally uses
# the strict terminal gate rather than legacy benchmark watchdog paths: the
# shared watchdog still assumes the ordinary benchmark cache naming contract.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/../benchmark_submit_common.sh"
cd "${PROJECT_ROOT}"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

PASS=""
DS_NAME_ARG=""
METHODS_ARG=""
PARTITION_ARG=""
FORCE_ARG=0
SYNC_ONLY_IDS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pass)
      PASS="${2:-}"
      shift 2
      ;;
    --pass=*)
      PASS="${1#*=}"
      shift
      ;;
    --ds_name)
      DS_NAME_ARG="${2:-}"
      shift 2
      ;;
    --ds_name=*)
      DS_NAME_ARG="${1#*=}"
      shift
      ;;
    --methods)
      METHODS_ARG="${2:-}"
      shift 2
      ;;
    --methods=*)
      METHODS_ARG="${1#*=}"
      shift
      ;;
    --partition)
      PARTITION_ARG="${2:-}"
      shift 2
      ;;
    --partition=*)
      PARTITION_ARG="${1#*=}"
      shift
      ;;
    --force)
      FORCE_ARG=1
      shift
      ;;
    --sync-only)
      SYNC_ONLY_IDS="${2:-}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires at least one job id (comma-separated)."
        exit 1
      fi
      shift 2
      ;;
    --sync-only=*)
      SYNC_ONLY_IDS="${1#*=}"
      if [[ -z "${SYNC_ONLY_IDS}" ]]; then
        echo "ERROR: --sync-only requires at least one job id (comma-separated)."
        exit 1
      fi
      shift
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ "${PASS}" != "uncorrected" && "${PASS}" != "corrected" ]]; then
  echo "ERROR: --pass must be exactly uncorrected or corrected."
  exit 1
fi
if [[ -n "${SYNC_ONLY_IDS}" && ${FORCE_ARG} -eq 1 ]]; then
  echo "ERROR: --sync-only cannot be combined with --force."
  exit 1
fi

ANALYSIS_VIEW="batch_effect_${PASS}"
ANALYSIS_ROOT="${HPC_SCRATCH_DIR}/batch_effect/${PASS}"
ANALYSIS_NAS_ROOT="${NAS_TARGET_DIR}/batch_effect/${PASS}"
ANALYSIS_LOG_PREFIX="execution_times_batch_effect_${PASS}_"
export ANALYSIS_VIEW ANALYSIS_PASS="${PASS}" ANALYSIS_ROOT ANALYSIS_NAS_ROOT
export ANALYSIS_LOG_PREFIX ANALYSIS_HIGH_RES_ONLY=1
export FORCE_BENCHMARK="${FORCE_ARG}"
# The batch workers use ANALYSIS_MANIFEST. Never populate the legacy manifest
# variable, even when the caller's shell inherited one.
unset BENCHMARK_MANIFEST
mkdir -p "${ANALYSIS_ROOT}/manifests" "${ANALYSIS_ROOT}/embeddings" \
  "${ANALYSIS_ROOT}/results" "${ANALYSIS_ROOT}/pseudobulks" \
  "${ANALYSIS_ROOT}/gloscope_dists" "${LOGS_DIR}"

DATASET_NAMES=()
if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" --arg view "${ANALYSIS_VIEW}" \
      '.[$ds] != null and .[$ds].use_for_batch_effect == true and .[$ds].views[$view] != null' \
      "${DATASETS_JSON_FILE}" >/dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not enabled for ${ANALYSIS_VIEW} in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  DATASET_NAMES=("${DS_NAME_ARG}")
else
  while IFS= read -r name; do
    DATASET_NAMES+=("${name}")
  done < <(jq -r --arg view "${ANALYSIS_VIEW}" \
    'to_entries[] |
     select(.value.use_for_batch_effect == true) |
     select(.value.views[$view] != null) |
     .key | select(startswith("_") | not)' "${DATASETS_JSON_FILE}")
fi

if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: No datasets expose ${ANALYSIS_VIEW} with use_for_batch_effect=true."
  exit 1
fi

if [[ "${PASS}" == "corrected" ]]; then
  for DS in "${DATASET_NAMES[@]}"; do
    BATCH_COL="$(jq -r --arg ds "${DS}" '.[$ds].columns.batch // empty' "${DATASETS_JSON_FILE}")"
    if [[ -z "${BATCH_COL}" ]]; then
      echo "ERROR: corrected batch-effect view requires a confirmed columns.batch (dataset ${DS})."
      exit 1
    fi
  done
fi

echo "Found ${#DATASET_NAMES[@]} ${ANALYSIS_VIEW} datasets."

METHODS=(prepare_pseudobulk pseudobulk gloscope composition mrvi pilot pilotgm qot)
if [[ -n "${METHODS_ARG}" ]]; then
  IFS=',' read -r -a METHODS <<< "${METHODS_ARG}"
fi
VALID_METHODS="prepare_pseudobulk|pseudobulk|gloscope|composition|mrvi|pilot|pilotgm|qot"
UNIQUE_METHODS=()
for M in "${METHODS[@]}"; do
  if ! [[ "${M}" =~ ^(${VALID_METHODS})$ ]]; then
    echo "ERROR: Unknown batch-effect method '${M}'."
    exit 1
  fi
  if [[ ! " ${UNIQUE_METHODS[*]} " =~ " ${M} " ]]; then
    UNIQUE_METHODS+=("${M}")
  fi
done
METHODS=("${UNIQUE_METHODS[@]}")

NEEDS_PREP=0
for M in "${METHODS[@]}"; do
  if [[ "${M}" == "pseudobulk" || "${M}" == "composition" ]]; then
    NEEDS_PREP=1
  fi
done
if [[ ${NEEDS_PREP} -eq 1 ]]; then
  HAS_PREP=0
  for M in "${METHODS[@]}"; do
    [[ "${M}" == "prepare_pseudobulk" ]] && HAS_PREP=1
  done
  if [[ ${HAS_PREP} -eq 0 ]]; then
    echo "NOTE: adding prepare_pseudobulk before pseudobulk/composition."
    METHODS=(prepare_pseudobulk "${METHODS[@]}")
  fi
fi

ALL_DATASET_NAMES=("${DATASET_NAMES[@]}")

method_spec() {
  local METHOD="$1"
  METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
  METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
  METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}"
                --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
  if [[ "${METHOD}" == "mrvi" || "${METHOD}" == "pilotgm" ]]; then
    METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_GPU}"
    METHOD_THROTTLE="${BENCHMARK_GPU_ARRAY_THROTTLE}"
    METHOD_FLAGS=(--gpus="${BENCHMARK_GPU_COUNT}"
                  --constraint="${BENCHMARK_GPU_CONSTRAINT}"
                  --cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}")
  elif [[ "${METHOD}" == "pilot" || "${METHOD}" == "qot" ]]; then
    METHOD_PARTITION="${SLURM_PARTITION_BENCHMARK_CPU}"
    METHOD_THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
    METHOD_FLAGS=(--constraint="${BENCHMARK_CPU_CONSTRAINT}"
                  --cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}")
  fi
  if [[ -n "${PARTITION_ARG}" ]]; then
    METHOD_PARTITION="${PARTITION_ARG}"
    FILTERED_FLAGS=()
    for FLAG in "${METHOD_FLAGS[@]}"; do
      case "${FLAG}" in
        --constraint|--constraint=*) ;;
        *) FILTERED_FLAGS+=("${FLAG}") ;;
      esac
    done
    METHOD_FLAGS=("${FILTERED_FLAGS[@]}")
  fi
}

submit_array() {
  local METHOD="$1"
  local MANIFEST="${ANALYSIS_ROOT}/manifests/${PASS}_batch_effect_manifest_${METHOD}_$$.txt"
  : > "${MANIFEST}"
  for DS in "${DATASET_NAMES[@]}"; do
    echo "${DS}" >> "${MANIFEST}"
  done
  export ANALYSIS_MANIFEST="${MANIFEST}" METHOD="${METHOD}"
  method_spec "${METHOD}"

  if [[ "${METHOD}" == "mrvi" || "${METHOD}" == "pilot" ||
        "${METHOD}" == "pilotgm" || "${METHOD}" == "qot" ]]; then
    WORKER_SCRIPT="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh"
    JOB_PREFIX="5_batch_effect_${PASS}_${METHOD}"
  else
    WORKER_SCRIPT="${PROJECT_ROOT}/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
    JOB_PREFIX="5_batch_effect_${PASS}_${METHOD}"
  fi
  export JOB_LOG_PREFIX="${JOB_PREFIX}"

  echo "Submitting ${METHOD} array (${#DATASET_NAMES[@]} datasets, pass=${PASS}, partition=${METHOD_PARTITION})" >&2
  SUBMIT_MSG=$(sbatch \
    --array="1-${#DATASET_NAMES[@]}%${METHOD_THROTTLE}" \
    --partition="${METHOD_PARTITION}" \
    "${METHOD_FLAGS[@]}" \
    --mem="${BENCHMARK_MEM}" \
    --output="${LOGS_DIR}/${JOB_PREFIX}_%A_%a.log" \
    --error="${LOGS_DIR}/${JOB_PREFIX}_%A_%a.err" \
    --mail-user="${USER_EMAIL}" \
    "${WORKER_SCRIPT}")
  ARRAY_JOB_ID="$(echo "${SUBMIT_MSG}" | grep -oE '[0-9]+')"
  echo "  ${METHOD} array job ID: ${ARRAY_JOB_ID}" >&2
  echo "${ARRAY_JOB_ID}"
}

run_one_dataset() {
  local DATASET="$1"
  local PREP_JOB_ID NEW_ARRAY_ID
  local -a ARRAY_JOB_IDS=()
  local -a ARRAY_JOB_METHODS=()
  DATASET_NAMES=("${DATASET}")
  echo "=== Running ${ANALYSIS_VIEW} for ${DATASET} ==="

  PREP_JOB_ID=""
  for M in "${METHODS[@]}"; do
    if [[ "${M}" == "prepare_pseudobulk" ]]; then
      PREP_JOB_ID="$(submit_array "${M}")"
      benchmark_wait_for_array "${PREP_JOB_ID}" "${ANALYSIS_VIEW}:${DATASET}:${M}"
      break
    fi
  done
  for M in "${METHODS[@]}"; do
    [[ "${M}" == "prepare_pseudobulk" ]] && continue
    NEW_ARRAY_ID="$(submit_array "${M}")"
    ARRAY_JOB_IDS+=("${NEW_ARRAY_ID}")
    ARRAY_JOB_METHODS+=("${M}")
  done
  for i in "${!ARRAY_JOB_IDS[@]}"; do
    benchmark_wait_for_array \
      "${ARRAY_JOB_IDS[$i]}" "${ANALYSIS_VIEW}:${DATASET}:${ARRAY_JOB_METHODS[$i]}"
  done

  # Sync/checksum/cleanup completes for this cohort before the next one starts.
  analysis_merge_sync_cleanup "${METHODS[@]}"
}

if [[ -n "${SYNC_ONLY_IDS}" ]]; then
  if [[ -z "${DS_NAME_ARG}" ]]; then
    echo "ERROR: --sync-only requires --ds_name so its logs can be scoped safely."
    exit 1
  fi
  DATASET_NAMES=("${DS_NAME_ARG}")
  IFS=',' read -r -a ARRAY_JOB_IDS <<< "${SYNC_ONLY_IDS}"
  echo "=== Sync-only resume mode for ${DS_NAME_ARG}: ${SYNC_ONLY_IDS} ==="
  for JOB_ID in "${ARRAY_JOB_IDS[@]}"; do
    benchmark_wait_for_array "${JOB_ID}" "${ANALYSIS_VIEW}:${DS_NAME_ARG}:resume:${JOB_ID}"
  done
  analysis_merge_sync_cleanup "${METHODS[@]}"
else
  for DS in "${ALL_DATASET_NAMES[@]}"; do
    run_one_dataset "${DS}"
  done
fi
