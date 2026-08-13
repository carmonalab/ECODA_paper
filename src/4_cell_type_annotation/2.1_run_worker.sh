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
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available on worker node; cannot read manifest/datasets.json."
  exit 1
fi

if [[ -z "${CHUNKS_MANIFEST:-}" ]]; then
  echo "ERROR: CHUNKS_MANIFEST is not set. Export it before submitting the array."
  exit 1
fi

# Read this task's line from the global manifest (written by
# 2_submit_hpc_array.sh):  DS_NAME<TAB>absolute_chunk_path
MANIFEST_LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHUNKS_MANIFEST}")"
if [[ -z "${MANIFEST_LINE}" ]]; then
  echo "ERROR: No manifest entry for array task ${SLURM_ARRAY_TASK_ID} in ${CHUNKS_MANIFEST}."
  exit 1
fi
IFS=$'\t' read -r DS_NAME CHUNK_FILE <<< "${MANIFEST_LINE}"

# Per-dataset env for 2.1.1_process_chunk.R (read via Sys.getenv()). If jq
# finds no tissue key for the dataset, the R defaults apply.
export DS_NAME
export TISSUE_TYPE="$(jq -r --arg ds "${DS_NAME}" '.[$ds].tissue // empty' "${DATASETS_JSON_FILE}")"
export NORMAL_TISSUE="$(jq -r --arg ds "${DS_NAME}" '.[$ds].normal_tissue // empty' "${DATASETS_JSON_FILE}")"
echo "Task ${SLURM_ARRAY_TASK_ID}: DS_NAME=${DS_NAME}, chunk=${CHUNK_FILE}"
echo "Exported TISSUE_TYPE=${TISSUE_TYPE}, NORMAL_TISSUE=${NORMAL_TISSUE} for ${DS_NAME}"

if [[ ! -f "${CHUNK_FILE}" ]]; then
  echo "Task ${SLURM_ARRAY_TASK_ID}: ERROR: Chunk file ${CHUNK_FILE} not found."
  exit 1
fi

# Skip chunks whose annotation feather already exists (re-runs of
# 2_submit_hpc_array.sh on an annotated-but-not-yet-merged dataset; same
# chunk_N.txt -> annotations_chunk_N.feather mapping as the merge). The chunk
# manifest stays unfiltered, so the 3_submit_merge.sh coverage gate (chunk
# files vs feather count) is unchanged; 1.1_prepare_chunks.py
# deletes stale feathers on every chunk rebuild (production), so an existing
# feather always matches the current chunk set.
CHUNK_NUM="$(basename "${CHUNK_FILE}")"
CHUNK_NUM="${CHUNK_NUM#chunk_}"
CHUNK_NUM="${CHUNK_NUM%.txt}"
FEATHER_FILE="${HPC_SCRATCH_DIR}/${DS_NAME}/output/annotations_chunk_${CHUNK_NUM}.feather"
if [[ -f "${FEATHER_FILE}" ]]; then
  echo "Task ${SLURM_ARRAY_TASK_ID}: ${FEATHER_FILE} already exists — annotation already done, skipping."
  exit 0
fi

# Unified retry handling: pin BLAS/OMP threads so CPU time ~= wall time, then
# run the chunk. R stderr lands in the Slurm .err file, so one
# transient-signature grep covers both. The retry counter is cleared on
# success and only bumped on requeue.
source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
export_worker_thread_env

set +e
${PIXI_RSCRIPT} \
  "${SCRIPT_DIR}/2.1.1_process_chunk.R" \
  "${CHUNK_FILE}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Task ${SLURM_ARRAY_TASK_ID}: chunk processing complete."
  exit 0
fi
ERR_FILE="${LOGS_DIR}/4_cell_type_annotation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  exit 0   # requeued; the script restarts, likely on another node
fi
exit ${RC}
