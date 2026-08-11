#!/bin/bash
#SBATCH --job-name=transzeroimp_worker
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
# This script sits two levels below src/ (slurm_config.sh lives at src/).
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"

# ANALYSIS and BENCHMARK_MANIFEST are exported by 1_submit_hpc_array.sh
# (sbatch propagates the submit script's environment).
if [[ -z "${ANALYSIS:-}" ]]; then
  echo "ERROR: ANALYSIS is not set. Export it before submitting the array."
  exit 1
fi
if [[ "${ANALYSIS}" != "trans" && "${ANALYSIS}" != "zeroimp" ]]; then
  echo "ERROR: ANALYSIS must be 'trans' or 'zeroimp' (got '${ANALYSIS}')."
  exit 1
fi
if [[ -z "${BENCHMARK_MANIFEST:-}" ]]; then
  echo "ERROR: BENCHMARK_MANIFEST is not set. Export it before submitting the array."
  exit 1
fi

# Read this task's dataset from the per-analysis manifest (one DS per line,
# written by 1_submit_hpc_array.sh): sed -n matches SLURM_ARRAY_TASK_ID 1:1.
DS_NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BENCHMARK_MANIFEST}")"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No manifest entry for array task ${SLURM_ARRAY_TASK_ID} in ${BENCHMARK_MANIFEST}."
  exit 1
fi
export DS_NAME
echo "Task ${SLURM_ARRAY_TASK_ID}: ANALYSIS=${ANALYSIS}, DS_NAME=${DS_NAME}"

OUT_DIR="${HPC_SCRATCH_DIR}/benchmark/embeddings"
mkdir -p "${OUT_DIR}"
# Log name is deterministic per (ANALYSIS, DS_NAME): the trans and zeroimp
# arrays each have a distinct ANALYSIS, so there is no read-modify-write
# collision, and re-runs overwrite the same file. The merge script globs
# the (analysis x dataset) cross product.
LOG_FILE="${OUT_DIR}/execution_times_${ANALYSIS}_${DS_NAME}.feather"

# FORCE_BENCHMARK is exported by 1_submit_hpc_array.sh (--force); forward it
# to the R scripts so existing outputs are recomputed.
FORCE_FLAG=""
if [[ "${FORCE_BENCHMARK:-0}" == "1" ]]; then
  FORCE_FLAG="--force"
fi

if [[ "${ANALYSIS}" == "trans" ]]; then
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_transformation_analysis.R"
else
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_zeroimp_analysis.R"
fi

# Staging + unified retry handling: stage the pixi R library to node-local
# /scratch (immune to stale BeeGFS client-cache views), then run the analysis.
# No thread pinning (benchmark hardware is pinned for runtime comparability).
# Both staging and R stderr land in the Slurm .err file, so one
# transient-signature grep covers both. Analysis outputs are overwritten on
# re-run (idempotent per-dataset files).
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"

echo "Task ${SLURM_ARRAY_TASK_ID}: running ${ANALYSIS} analysis on ${DS_NAME}"
# PIXI_RSCRIPT word-splits into `pixi run -e py-cuda13 Rscript --vanilla`;
# it must stay unquoted (established convention, see 2.1_run_worker.sh).
set +e
stage_env_rlib "benchmark" && ${PIXI_RSCRIPT} "${R_SCRIPT}" \
    --config_path "${DATASETS_JSON_FILE}" \
    --ds_name "${DS_NAME}" \
    --view benchmark_analysis \
    --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" \
    --output_dir "${HPC_SCRATCH_DIR}/benchmark/results" \
    --log_file "${LOG_FILE}" \
    ${FORCE_FLAG}
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Task ${SLURM_ARRAY_TASK_ID}: ${ANALYSIS} on ${DS_NAME} complete."
  exit 0
fi
ERR_FILE="${LOGS_DIR}/5_transzeroimp_${ANALYSIS}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  exit 0   # requeued; the script restarts, likely on another node
fi
exit ${RC}
