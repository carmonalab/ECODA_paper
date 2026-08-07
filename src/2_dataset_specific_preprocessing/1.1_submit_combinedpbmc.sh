#!/bin/bash
#SBATCH --job-name=combine_pbmc
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
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

module load GCCcore/12.2.0

"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_create_combinedpbmc_dataset.py"
