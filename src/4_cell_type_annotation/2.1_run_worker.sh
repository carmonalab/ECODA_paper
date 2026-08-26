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
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

MANIFEST_PATH="${CHUNKS_MANIFEST:-}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: CHUNKS_MANIFEST is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
MANIFEST_LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
[[ -n "${MANIFEST_LINE}" ]] || { echo "ERROR: no chunk manifest row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME CHUNK_FILE FEATHER_DIR <<< "${MANIFEST_LINE}"
[[ -n "${DS_NAME}" && -n "${CHUNK_FILE}" ]] || { echo "ERROR: malformed chunk manifest row" >&2; exit 1; }
FEATHER_DIR="${FEATHER_DIR:-${HPC_SCRATCH_DIR}/${DS_NAME}/output}"
export DS_NAME
export ANNOTATION_OUTPUT_DIR="${FEATHER_DIR}"
export TISSUE_TYPE="$(jq -r --arg ds "${DS_NAME}" '.[$ds].tissue // empty' "${DATASETS_JSON_FILE}")"
export NORMAL_TISSUE="$(jq -r --arg ds "${DS_NAME}" 'if .[$ds].normal_tissue == null then empty else .[$ds].normal_tissue end' "${DATASETS_JSON_FILE}")"
[[ -f "${CHUNK_FILE}" ]] || { echo "ERROR: chunk file not found: ${CHUNK_FILE}" >&2; exit 1; }

CHUNK_NUM="$(basename "${CHUNK_FILE}")"
CHUNK_NUM="${CHUNK_NUM#chunk_}"
CHUNK_NUM="${CHUNK_NUM%.txt}"
FEATHER_FILE="${FEATHER_DIR}/annotations_chunk_${CHUNK_NUM}.feather"
mkdir -p "${FEATHER_DIR}"
if [[ -s "${FEATHER_FILE}" ]]; then
  echo "Annotation feather already exists and is nonempty: ${FEATHER_FILE}"
  exit 0
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
  exit 0
fi
ERR_FILE="${ANNOTATION_ERROR_PREFIX:-${LOGS_DIR}/4_cell_type_annotation}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  rm -f "${FEATHER_FILE}"
  exit 0
fi
exit ${RC}
