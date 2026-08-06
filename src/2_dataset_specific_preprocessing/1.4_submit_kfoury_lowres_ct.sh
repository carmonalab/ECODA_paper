#!/bin/bash
#SBATCH --job-name=kfoury_lowres_ct
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

${PIXI_RSCRIPT} "${SCRIPT_DIR}/1.4.1_create_kfoury_lowres_ct.R"
