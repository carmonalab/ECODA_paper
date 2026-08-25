#!/bin/bash
# Durable-gate entrypoint for produced annotation roles. Author-annotated
# cohorts do not need a redundant annotation pass; the current onboarding
# derived role is Lupus HiTME layer1/layer2. Independent datasets may be passed
# as a comma-separated list when another produced role is added.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG="Lupus_PBMC"
MEMORY="32G"
MAX_MEMORY="128G"
PARTITION="${SLURM_PARTITION_BENCHMARK_CPU:-shared-cpu}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets)
      DATASETS_ARG="${2:-}"
      shift 2
      ;;
    --datasets=*)
      DATASETS_ARG="${1#*=}"
      shift
      ;;
    --mem)
      MEMORY="${2:-}"
      shift 2
      ;;
    --max-mem)
      MAX_MEMORY="${2:-}"
      shift 2
      ;;
    --partition)
      PARTITION="${2:-}"
      shift 2
      ;;
    --throttle)
      THROTTLE="${2:-}"
      shift 2
      ;;
    -h|--help)
      echo "Usage: 1_submit_onboarding_stage.sh [--datasets LIST] [--mem VALUE] [--max-mem VALUE] [--partition NAME] [--throttle N]"
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required for annotation stage validation." >&2
  exit 1
fi

DATASET_NAMES=()
OLD_IFS="${IFS}"
IFS=','
read -r -a DATASET_NAMES <<< "${DATASETS_ARG}"
IFS="${OLD_IFS}"
if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: no annotation datasets selected." >&2
  exit 1
fi
for DS_NAME in "${DATASET_NAMES[@]}"; do
  jq -e --arg ds "${DS_NAME}" \
    '.[$ds] and .[$ds].use_for_batch_effect == true' \
    "${DATASETS_JSON_FILE}" >/dev/null || {
      echo "ERROR: ${DS_NAME} is not a batch-effect dataset." >&2
      exit 1
    }
done

# Chunk preparation uses srun synchronously on the compute partition. It is
# intentionally before the array so every selected dataset contributes to one
# immutable manifest and the watchdog can retry exact chunk rows.
for DS_NAME in "${DATASET_NAMES[@]}"; do
  "${SCRIPT_DIR}/1_prepare_chunks.sh" production "${DS_NAME}"
done

MANIFEST_DIR="${HPC_SCRATCH_DIR}/_annotation_stage_manifests"
mkdir -p "${MANIFEST_DIR}" "${LOGS_DIR}"
MANIFEST="${MANIFEST_DIR}/onboarding_$(date +%Y%m%d%H%M%S)_$$.txt"
: > "${MANIFEST}"
TOTAL_CHUNKS=0
for DS_NAME in "${DATASET_NAMES[@]}"; do
  CHUNKS_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks"
  for CHUNK_FILE in "${CHUNKS_DIR}"/chunk_*.txt; do
    [[ -f "${CHUNK_FILE}" ]] || continue
    printf '%s\t%s\n' "${DS_NAME}" "${CHUNK_FILE}" >> "${MANIFEST}"
    TOTAL_CHUNKS=$((TOTAL_CHUNKS + 1))
  done
done
if [[ ${TOTAL_CHUNKS} -eq 0 ]]; then
  echo "ERROR: no annotation chunks found for ${DATASET_NAMES[*]}." >&2
  exit 1
fi

ARRAY_ID="$(sbatch --parsable \
  --array="1-${TOTAL_CHUNKS}%${THROTTLE}" \
  --mem="${MEMORY}" \
  --partition="${PARTITION}" \
  --output="${LOGS_DIR}/4_cell_type_annotation_%A_%a.log" \
  --error="${LOGS_DIR}/4_cell_type_annotation_%A_%a.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,CHUNKS_MANIFEST=${MANIFEST},ANNOTATION_ERROR_PREFIX=${LOGS_DIR}/4_cell_type_annotation" \
  "${SCRIPT_DIR}/2.1_run_worker.sh")"
ARRAY_ID="${ARRAY_ID%%;*}"
[[ "${ARRAY_ID}" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid annotation array id: ${ARRAY_ID}" >&2; exit 1; }

WATCHDOG_ID="$(sbatch --parsable --wait \
  --dependency="afterany:${ARRAY_ID}" \
  --partition="${PARTITION}" \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=2G \
  --time="${ANNOTATION_WATCHDOG_TIME_LIMIT:-12:00:00}" \
  --output="${LOGS_DIR}/4_cell_type_annotation_watchdog_%j.log" \
  --error="${LOGS_DIR}/4_cell_type_annotation_watchdog_%j.err" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL" \
  "${SCRIPT_DIR}/1.2_annotation_watchdog.sh" \
  "${ARRAY_ID}" "${MANIFEST}" "${MEMORY}" "${MAX_MEMORY}" \
  "${PARTITION}" "${THROTTLE}")"
WATCHDOG_ID="${WATCHDOG_ID%%;*}"
[[ "${WATCHDOG_ID}" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid annotation watchdog id: ${WATCHDOG_ID}" >&2; exit 1; }

printf 'ANNOTATION_ARRAY_JOB_ID=%s\n' "${ARRAY_ID}"
printf 'ANNOTATION_WATCHDOG_JOB_ID=%s\n' "${WATCHDOG_ID}"
printf 'ANNOTATION_MANIFEST=%s\n' "${MANIFEST}"
printf 'ANNOTATION_DATASETS=%s\n' "${DATASET_NAMES[*]}"

# Merge is login-side because compute nodes cannot reach NAS. It runs only
# after the watchdog's terminal scratch validation and is part of this gate.
for DS_NAME in "${DATASET_NAMES[@]}"; do
  "${SCRIPT_DIR}/3_submit_merge.sh" "${DS_NAME}"
done
printf 'ANNOTATION_MERGE=OK\n'
