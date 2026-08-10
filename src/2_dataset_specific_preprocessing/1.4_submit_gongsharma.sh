#!/bin/bash
#SBATCH --job-name=gongsharma_cap
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL

# ---------------------------------------------------------------------------
# GongSharma per-sample 5000-cell cap step.
#
# Runs 1.4.1_subset_gongsharma.py, which caps each sample
# (specimen.specimenGuid) at 5000 cells in the two STAGED SoundLife h5ads and
# OVERWRITES THEM IN PLACE (atomic temp-file write + os.replace):
#   ${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data/
#     SoundLife_YoungAdult_Male_CMVneg.h5ad  1,712,244 -> 531,291 cells
#     SoundLife_YoungAdult_Male_CMVpos.h5ad  1,206,761 -> 365,000 cells
# Total 896,291 cells / 180 samples / max 5000 per sample (seed 42; identical
# cell set to the historical NAS artifact Gongsharma_cmv_young_males.h5ad).
# Fixes the preprocess-array OOM (2.92M-cell sc.concat + densified
# sc.pp.scale exceeded the 128G worker).
#
# RE-STAGING CAVEAT: src/1_stage_data/1_stage_data.sh rsyncs the NAS
# originals back over the capped files — re-run this step (via 1_submit_hpc.sh)
# before the preprocess array after any re-stage.
#
# ORDERING: 1_submit_hpc.sh submits this step FIRST and gates the CombinedPBMC
# step (1.1_submit_combinedpbmc.sh) behind it via --dependency=afterok: that
# step reads the SAME staged files in backed mode, so an in-place overwrite
# racing its read would nondeterminize the CombinedPBMC dataset.
# ---------------------------------------------------------------------------

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

"${PYTHON_BIN}" "${SCRIPT_DIR}/1.4.1_subset_gongsharma.py"
