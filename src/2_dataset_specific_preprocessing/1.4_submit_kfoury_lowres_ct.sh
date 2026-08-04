#!/bin/bash
#SBATCH --job-name=kfoury_lowres_ct
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
# NOTE: 32G baseline — the whole Kfoury .rds is read via readRDS + saveRDS
# in a single process.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

${PIXI_RSCRIPT} "${SCRIPT_DIR}/1.4.1_create_kfoury_lowres_ct.R"
