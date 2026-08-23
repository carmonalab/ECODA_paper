#!/bin/bash
#SBATCH --job-name=bassez_fill_subtype
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/ {print $2}' | xargs dirname)"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

${PIXI_RSCRIPT} "${SCRIPT_DIR}/1.6.1_fill_bassez_cellsubtype.R"
