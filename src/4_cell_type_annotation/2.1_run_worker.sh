#!/bin/bash
#SBATCH --job-name=scrna_worker
#SBATCH --partition=shared-cpu             
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

# HOME_CHUNKS_DIR is exported by 2_submit_hpc_array.sh; set a fallback if missing
: "${HOME_CHUNKS_DIR:=${SCRATCH_OUTPUT_DIR}/chunks}"

CHUNK_FILE="${HOME_CHUNKS_DIR}/chunk_${SLURM_ARRAY_TASK_ID}.txt"

# PROJECT_ROOT is sourced (and exported) by 2.1.1_process_chunk.sh via slurm_config.sh
bash "${SCRIPT_DIR}/2.1.1_process_chunk.sh" "${SLURM_ARRAY_TASK_ID}" "${CHUNK_FILE}"