#!/bin/bash
# Submit the uncorrected onboarding preprocessing stage as one durable SLURM
# array plus a compute-node watchdog. The submitter is intended to run inside
# durable-hpc-gate; it blocks on the watchdog via scheduler-native --wait.
#
# The submitter never implements a login-node squeue/sacct polling loop.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

VIEW="batch_effect_uncorrected"
MEMORY="128G"
MAX_MEMORY="${PREPROCESS_MEM_MAX:-500G}"
PARTITION="${SLURM_PARTITION_BENCHMARK_CPU:-shared-cpu}"
THROTTLE="${MAX_NUM_CHUNKS_PARALLEL}"
FORCE_PREPROCESS=0
DATASETS_ARG=""

usage() {
  cat <<'EOF'
Usage: 1_submit_batch_effect_stage.sh [options]

Submit all onboarding cohorts for one preprocessing view as one SLURM array.
The command submits the worker array and compute-node OOM watchdog, then
blocks on the watchdog with scheduler-native --wait for durable-gate use.

Options:
  --view NAME       preprocessing view (only batch_effect_uncorrected here)
  --datasets LIST   comma-separated dataset IDs; default is the onboarding set
  --mem VALUE       initial worker memory (default: 128G)
  --max-mem VALUE   OOM retry ceiling (default: PREPROCESS_MEM_MAX or 500G)
  --partition NAME  SLURM partition (default: shared-cpu)
  --throttle N      maximum concurrent array tasks
  --force           recompute existing outputs
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --view)
      VIEW="${2:-}"
      shift 2
      ;;
    --view=*)
      VIEW="${1#*=}"
      shift
      ;;
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
    --mem=*)
      MEMORY="${1#*=}"
      shift
      ;;
    --max-mem)
      MAX_MEMORY="${2:-}"
      shift 2
      ;;
    --max-mem=*)
      MAX_MEMORY="${1#*=}"
      shift
      ;;
    --partition)
      PARTITION="${2:-}"
      shift 2
      ;;
    --partition=*)
      PARTITION="${1#*=}"
      shift
      ;;
    --throttle)
      THROTTLE="${2:-}"
      shift 2
      ;;
    --throttle=*)
      THROTTLE="${1#*=}"
      shift
      ;;
    --force)
      FORCE_PREPROCESS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ "${VIEW}" != "batch_effect_uncorrected" ]]; then
  echo "ERROR: this stage submitter only accepts --view batch_effect_uncorrected." >&2
  exit 1
fi
if [[ -z "${MEMORY}" || -z "${MAX_MEMORY}" || -z "${PARTITION}" || -z "${THROTTLE}" ]]; then
  echo "ERROR: memory, partition, and throttle values must be non-empty." >&2
  exit 1
fi
if ! [[ "${THROTTLE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --throttle must be a positive integer." >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required to resolve onboarding datasets." >&2
  exit 1
fi

DATASET_NAMES=()
if [[ -n "${DATASETS_ARG}" ]]; then
  OLD_IFS="${IFS}"
  IFS=','
  read -r -a DATASET_NAMES <<< "${DATASETS_ARG}"
  IFS="${OLD_IFS}"
else
  # The nine newly onboarded cohorts are the use_for_batch_effect-only entries
  # other than the legacy CombinedPBMC entry. Joanito and Stephenson are the
  # two established reference cohorts with confirmed technical batches.
  while IFS= read -r DS_NAME; do
    [[ -n "${DS_NAME}" ]] || continue
    DATASET_NAMES+=("${DS_NAME}")
  done < <(jq -r 'to_entries[]
    | select(.key != "CombinedPBMC")
    | select(.key != "Joanito")
    | select(.key != "Stephenson")
    | select(.value.use_for_batch_effect == true)
    | select(.value.use_for_benchmark == false)
    | .key' "${DATASETS_JSON_FILE}")
  DATASET_NAMES+=("Joanito" "Stephenson")
fi

if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: no onboarding datasets selected." >&2
  exit 1
fi

for DS_NAME in "${DATASET_NAMES[@]}"; do
  if ! jq -e --arg ds "${DS_NAME}" \
      '.[$ds] and .[$ds].use_for_batch_effect == true' \
      "${DATASETS_JSON_FILE}" >/dev/null; then
    echo "ERROR: ${DS_NAME} is not a batch-effect dataset in datasets.json." >&2
    exit 1
  fi
done

MANIFEST_DIR="${HPC_SCRATCH_DIR}/_preprocessing_stage_manifests"
WATCHDOG_STATUS_DIR="${HPC_SCRATCH_DIR}/_preprocessing_watchdog"
mkdir -p "${MANIFEST_DIR}" "${WATCHDOG_STATUS_DIR}" "${LOGS_DIR}"
MANIFEST="${MANIFEST_DIR}/${VIEW}_$(date +%Y%m%d%H%M%S)_$$.txt"
printf '%s\n' "${DATASET_NAMES[@]}" > "${MANIFEST}"

sync_outputs_to_nas() {
  local DS_NAME OUTPUT_NAME OUTPUT_FILE DEST_DIR
  if [[ ! -d "${NAS_TARGET_DIR}" ]]; then
    echo "ERROR: NAS path is unreachable: ${NAS_TARGET_DIR}" >&2
    return 1
  fi
  for DS_NAME in "${DATASET_NAMES[@]}"; do
    OUTPUT_NAME="$(jq -er --arg ds "${DS_NAME}" --arg view "${VIEW}" \
      '.[$ds].views[$view].output_file_name' "${DATASETS_JSON_FILE}")"
    OUTPUT_FILE="${HPC_SCRATCH_DIR}/${DS_NAME}/output/${OUTPUT_NAME}"
    if [[ ! -s "${OUTPUT_FILE}" ]]; then
      echo "ERROR: expected preprocessing artifact is missing or empty: ${OUTPUT_FILE}" >&2
      return 1
    fi
    DEST_DIR="${NAS_TARGET_DIR}/${DS_NAME}/output"
    mkdir -p "${DEST_DIR}"
    rsync -rlptD "${HPC_SCRATCH_DIR}/${DS_NAME}/output/" "${DEST_DIR}/"
  done
}


ARRAY_OUTPUT="${LOGS_DIR}/3_scrnaseq_batch_effect_%A_%a.log"
ARRAY_ERROR="${LOGS_DIR}/3_scrnaseq_batch_effect_%A_%a.err"
WATCHDOG_OUTPUT="${LOGS_DIR}/3_scrnaseq_batch_effect_watchdog_%j.log"
WATCHDOG_ERROR="${LOGS_DIR}/3_scrnaseq_batch_effect_watchdog_%j.err"
EXPORTS="ALL,PREPROCESS_DATASETS_FILE=${MANIFEST},PREPROCESS_VIEW=${VIEW},FORCE_PREPROCESS=${FORCE_PREPROCESS}"

ARRAY_ID="$(sbatch --parsable \
  --array="1-${#DATASET_NAMES[@]}%${THROTTLE}" \
  --mem="${MEMORY}" \
  --partition="${PARTITION}" \
  --output="${ARRAY_OUTPUT}" \
  --error="${ARRAY_ERROR}" \
  --mail-user="${USER_EMAIL}" \
  --export="${EXPORTS}" \
  "${SCRIPT_DIR}/1.1_run_worker.sh")"
ARRAY_ID="${ARRAY_ID%%;*}"
if ! [[ "${ARRAY_ID}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: sbatch returned an invalid preprocessing array ID: ${ARRAY_ID}" >&2
  exit 1
fi

WATCHDOG_ID="$(sbatch --parsable --wait \
  --dependency="afterany:${ARRAY_ID}" \
  --partition="${PARTITION}" \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=2G \
  --time="${PREPROCESS_WATCHDOG_TIME_LIMIT:-12:00:00}" \
  --output="${WATCHDOG_OUTPUT}" \
  --error="${WATCHDOG_ERROR}" \
  --mail-user="${USER_EMAIL}" \
  --export="ALL,PREPROCESS_VIEW=${VIEW},FORCE_PREPROCESS=${FORCE_PREPROCESS}" \
  "${SCRIPT_DIR}/1.2_preprocess_watchdog.sh" \
  "${ARRAY_ID}" "${MANIFEST}" "${MEMORY}" "${MAX_MEMORY}" \
  "${PARTITION}" "${THROTTLE}")"
WATCHDOG_ID="${WATCHDOG_ID%%;*}"
if ! [[ "${WATCHDOG_ID}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: sbatch returned an invalid preprocessing watchdog ID: ${WATCHDOG_ID}" >&2
  exit 1
fi

printf 'PREPROCESS_ARRAY_JOB_ID=%s\n' "${ARRAY_ID}"
printf 'PREPROCESS_WATCHDOG_JOB_ID=%s\n' "${WATCHDOG_ID}"
printf 'PREPROCESS_DATASET_MANIFEST=%s\n' "${MANIFEST}"
printf 'PREPROCESS_VIEW=%s\n' "${VIEW}"
printf 'PREPROCESS_DATASETS=%s\n' "${DATASET_NAMES[*]}"
printf 'PREPROCESS_MEMORY=%s\n' "${MEMORY}"
printf 'PREPROCESS_MAX_MEMORY=%s\n' "${MAX_MEMORY}"
printf 'PREPROCESS_PARTITION=%s\n' "${PARTITION}"
if [[ "${PREPROCESS_SUBMITTER_TEST:-0}" == "1" ]]; then
  printf 'PREPROCESS_NAS_SYNC=TEST_SKIPPED\n'
elif ! sync_outputs_to_nas; then
  exit 1
else
  printf 'PREPROCESS_NAS_SYNC=OK\n'
fi
