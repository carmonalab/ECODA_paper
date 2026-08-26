#!/bin/bash
# Durable-gate entrypoint for cell type annotation.
# Runs dual automated annotation (HiTME + scATOMIC) across all benchmark and
# suitable batch-effect cohorts (skipping datasets flagged
# not_suitable_for_auto_annotation in datasets.json).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

DATASETS_ARG=""
VIEW_ARG=""
FORCE_ARG=0
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
    --view)
      VIEW_ARG="${2:-}"
      shift 2
      ;;
    --view=*)
      VIEW_ARG="${1#*=}"
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
    --force)
      FORCE_ARG=1
      shift
      ;;
    -h|--help)
      echo "Usage: 1_submit_onboarding_stage.sh [--datasets LIST] [--view NAME] [--mem VALUE] [--max-mem VALUE] [--partition NAME] [--throttle N] [--force]"
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
if [[ -n "${DATASETS_ARG}" ]]; then
  OLD_IFS="${IFS}"
  IFS=','
  read -r -a DATASET_NAMES <<< "${DATASETS_ARG}"
  IFS="${OLD_IFS}"
else
  # Default: all datasets in datasets.json (skipping keys starting with "_" unless explicitly passed)
  while IFS= read -r name; do
    [[ "${name}" == _* ]] && continue
    DATASET_NAMES+=("${name}")
  done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")
fi

if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: no annotation datasets selected." >&2
  exit 1
fi

# Filter datasets: validate existence and filter out unsuitable cohorts
ACTIVE_DATASET_NAMES=()
for DS_NAME in "${DATASET_NAMES[@]}"; do
  if ! jq -e --arg ds "${DS_NAME}" 'has($ds)' "${DATASETS_JSON_FILE}" >/dev/null 2>&1; then
    echo "ERROR: '${DS_NAME}' is not defined in ${DATASETS_JSON_FILE}." >&2
    exit 1
  fi
  if jq -e --arg ds "${DS_NAME}" '.[$ds].not_suitable_for_auto_annotation // [] | (index("hitme") and index("scatomic"))' "${DATASETS_JSON_FILE}" >/dev/null 2>&1; then
    echo "NOTE: ${DS_NAME} is flagged not_suitable_for_auto_annotation in ${DATASETS_JSON_FILE}; skipping from annotation stage."
    continue
  fi
  ACTIVE_DATASET_NAMES+=("${DS_NAME}")
done

if [[ ${#ACTIVE_DATASET_NAMES[@]} -eq 0 ]]; then
  echo "NOTE: all selected datasets are marked not_suitable_for_auto_annotation. Nothing to annotate."
  exit 0
fi

# Chunk preparation: run sequentially before the array so every selected dataset
# contributes to one immutable manifest and the watchdog can retry exact chunk rows.
PREPARE_EXTRA_ARGS=()
if [[ -n "${VIEW_ARG}" ]]; then
  PREPARE_EXTRA_ARGS+=(--view "${VIEW_ARG}")
fi
if [[ ${FORCE_ARG} -eq 1 ]]; then
  PREPARE_EXTRA_ARGS+=(--force)
fi

for DS_NAME in "${ACTIVE_DATASET_NAMES[@]}"; do
  "${SCRIPT_DIR}/1_prepare_chunks.sh" production "${DS_NAME}" "${PREPARE_EXTRA_ARGS[@]}"
done

MANIFEST_DIR="${HPC_SCRATCH_DIR}/_annotation_stage_manifests"
mkdir -p "${MANIFEST_DIR}" "${LOGS_DIR}"
MANIFEST="${MANIFEST_DIR}/annotation_$(date +%Y%m%d%H%M%S)_$$.txt"
: > "${MANIFEST}"
TOTAL_CHUNKS=0
for DS_NAME in "${ACTIVE_DATASET_NAMES[@]}"; do
  CHUNKS_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks"
  shopt -s nullglob
  CHUNK_FILES=("${CHUNKS_DIR}"/chunk_*.txt)
  shopt -u nullglob
  for CHUNK_FILE in "${CHUNK_FILES[@]}"; do
    [[ -f "${CHUNK_FILE}" ]] || continue
    printf '%s\t%s\n' "${DS_NAME}" "${CHUNK_FILE}" >> "${MANIFEST}"
    TOTAL_CHUNKS=$((TOTAL_CHUNKS + 1))
  done
done

if [[ ${TOTAL_CHUNKS} -eq 0 ]]; then
  echo "NOTE: No pending chunks found for ${ACTIVE_DATASET_NAMES[*]}. Datasets appear already annotated."
  echo "Proceeding directly to merge step."
  for DS_NAME in "${ACTIVE_DATASET_NAMES[@]}"; do
    "${SCRIPT_DIR}/3_submit_merge.sh" "${DS_NAME}"
  done
  printf 'ANNOTATION_MERGE=OK\n'
  exit 0
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
printf 'ANNOTATION_DATASETS=%s\n' "${ACTIVE_DATASET_NAMES[*]}"

# Merge is login-side because compute nodes cannot reach NAS. It runs only
# after the watchdog's terminal scratch validation and is part of this gate.
for DS_NAME in "${ACTIVE_DATASET_NAMES[@]}"; do
  "${SCRIPT_DIR}/3_submit_merge.sh" "${DS_NAME}"
done
printf 'ANNOTATION_MERGE=OK\n'
