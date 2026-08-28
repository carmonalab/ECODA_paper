#!/bin/bash
#SBATCH --job-name=myocardial_prep
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

# ---------------------------------------------------------------------------
# Myocardial Infarction raw count reconstruction step.
#
# Runs 1.5.1_reconstruct_myocardial_counts.py, which inverts the log1p-normalized
# expression in Myocardial_Infarc_2.h5ad into exact raw UMI integer counts via
# cell-wise minimum step inversion and vaults them into adata.layers["counts"].
#
# Runs 1.5.1_reconstruct_myocardial_counts.py; FORCE_PREPROCESS=1 is translated
# to --force so intentional recomputation survives the scheduler hook.
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

FORCE_FLAG=()
[[ "${FORCE_PREPROCESS:-0}" == "1" ]] && FORCE_FLAG=(--force)
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.5.1_reconstruct_myocardial_counts.py" "${FORCE_FLAG[@]}"
