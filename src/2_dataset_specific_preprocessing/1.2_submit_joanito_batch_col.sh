#!/bin/bash
#SBATCH --job-name=joanito_batch_col
#SBATCH --partition=shared-cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
# NOTE: 32G baseline — the whole Joanito .rds is read via readRDS + saveRDS
# in a single process; 8G was too tight for the full object.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

${PIXI_RSCRIPT} "${SCRIPT_DIR}/_create_joanito_batch_col.R"
