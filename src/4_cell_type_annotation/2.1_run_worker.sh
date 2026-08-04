#!/bin/bash
#SBATCH --job-name=scrna_worker
#SBATCH --time=02:00:00                    
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                  
#SBATCH --mem=16G                          
#SBATCH --mail-type=END,FAIL
# NOTE: 16G should be sufficient for 5 samples/chunk (each subset from h5ad).
# If OOM occurs, bump to 32G and resubmit the failed chunk only.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

# jq module loads in the submit script do not propagate to workers; load here.
module load jq/1.6 >/dev/null 2>&1 || true
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

echo "Task ${SLURM_ARRAY_TASK_ID}: running annotation (pixi run Rscript --vanilla)"
${PIXI_RSCRIPT} \
  "${SCRIPT_DIR}/2.1.1_process_chunk.R" \
  "${CHUNK_FILE}"
echo "Task ${SLURM_ARRAY_TASK_ID}: chunk processing complete."
