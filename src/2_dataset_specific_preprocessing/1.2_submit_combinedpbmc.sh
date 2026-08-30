#!/bin/bash
#SBATCH --job-name=combine_pbmc
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage2 \
  "${SCRIPT_DIR}/1.2_submit_combinedpbmc.sh" || exit 1
cd "${PROJECT_ROOT}"

if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  module load GCCcore/12.2.0
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/1.2.1_create_combinedpbmc_dataset.py"
