#!/bin/bash
#SBATCH --job-name=scrna_worker
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 \
  "${SCRIPT_DIR}/2.1_run_worker.sh" || exit 1
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

MANIFEST_PATH="${CHUNKS_MANIFEST:-}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: CHUNKS_MANIFEST is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
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
mkdir -p "${FEATHER_DIR}"
if [[ -s "${FEATHER_FILE}" ]]; then
  if "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
      --path "${FEATHER_FILE}" --require-sidecar >/dev/null 2>&1 &&
     ecoda_validate_checksum "${FEATHER_FILE}"; then
    echo "Annotation feather already exists and passed schema/checksum validation: ${FEATHER_FILE}"
    exit 0
  fi
  echo "Existing annotation feather failed schema/checksum validation; rebuilding: ${FEATHER_FILE}" >&2
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
  "${PYTHON_BIN}" "${PROJECT_ROOT}/src/utils/py/annotation_contract.py" \
    --path "${FEATHER_FILE}" --require-sidecar >/dev/null 2>&1 || {
    echo "ERROR: annotation worker produced an invalid feather: ${FEATHER_FILE}" >&2
    exit 1
  }
  ecoda_validate_checksum "${FEATHER_FILE}" || {
    echo "ERROR: annotation worker feather checksum failed: ${FEATHER_FILE}" >&2
    exit 1
  }
  exit 0
fi
ERR_FILE="${ANNOTATION_ERROR_PREFIX:-${LOGS_DIR}/4_cell_type_annotation}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  rm -f "${FEATHER_FILE}"
  exit 0
fi
exit ${RC}
