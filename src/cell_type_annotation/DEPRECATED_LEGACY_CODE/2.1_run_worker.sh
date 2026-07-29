#!/bin/bash
#SBATCH --job-name=scrna_worker
#SBATCH --partition=shared-cpu             
#SBATCH --time=02:00:00                    
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                  
#SBATCH --mem=8G                           
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source "./config.env"
cd "${PROJECT_ROOT}"

CHUNK_FILE="${HOME_CHUNKS_DIR}/chunk_${SLURM_ARRAY_TASK_ID}.txt"

bash "./2.2_process_chunk.sh" "${SLURM_ARRAY_TASK_ID}" "${PROJECT_ROOT}" "${CHUNK_FILE}"