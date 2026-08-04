#!/bin/bash
#SBATCH --job-name=joanito_batch_col
#SBATCH --partition=shared-cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

"${HOME}/.pixi/bin/pixi" run Rscript --vanilla "${SCRIPT_DIR}/_create_joanito_batch_col.R"
