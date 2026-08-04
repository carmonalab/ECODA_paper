#!/bin/bash
#SBATCH --job-name=combine_pbmc
#SBATCH --partition=shared-cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
# NOTE: 64G baseline — GongSharma is huge and may need more.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load GCCcore/12.2.0

python "${SCRIPT_DIR}/_create_combinedpbmc_dataset.py"
