#!/bin/bash
#SBATCH --job-name=benchmark_r_worker
#SBATCH --time=12:00:00
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
# This script sits two levels below src/ (slurm_config.sh lives at src/).
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"

# METHOD and the generic analysis manifest are exported by the submitter.
# Ordinary benchmark runs still use BENCHMARK_MANIFEST as the fallback; batch
# runs use ANALYSIS_MANIFEST and never populate the legacy variable.
if [[ -z "${METHOD:-}" ]]; then
  echo "ERROR: METHOD is not set. Export it before submitting the array."
  exit 1
fi
MANIFEST_PATH="${ANALYSIS_MANIFEST:-${BENCHMARK_MANIFEST:-}}"
if [[ -z "${MANIFEST_PATH}" ]]; then
  echo "ERROR: ANALYSIS_MANIFEST (or legacy BENCHMARK_MANIFEST) is not set."
  exit 1
fi

DS_NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No manifest entry for array task ${SLURM_ARRAY_TASK_ID} in ${MANIFEST_PATH}."
  exit 1
fi
export DS_NAME
echo "Task ${SLURM_ARRAY_TASK_ID}: METHOD=${METHOD}, DS_NAME=${DS_NAME}"

ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
ANALYSIS_VIEW="${ANALYSIS_VIEW:-benchmark_analysis}"
RESULTS_DIR="${ANALYSIS_ROOT}/results"
PSEUDOBULK_DIR="${ANALYSIS_ROOT}/pseudobulks"
GLOSCOPE_DIR="${ANALYSIS_ROOT}/gloscope_dists"
EMBED_DIR="${ANALYSIS_ROOT}/embeddings"
mkdir -p "${RESULTS_DIR}" "${PSEUDOBULK_DIR}" "${GLOSCOPE_DIR}" "${EMBED_DIR}"
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  LOG_FILE="${EMBED_DIR}/execution_times_batch_effect_${ANALYSIS_PASS}_${METHOD}_${DS_NAME}.feather"
else
  LOG_FILE="${EMBED_DIR}/execution_times_${METHOD}_${DS_NAME}.feather"
fi

# FORCE_BENCHMARK is exported by 1_submit_hpc_array.sh (--force); forward it
# to the R scripts so existing per-combo caches are recomputed.
FORCE_FLAG=""
if [[ "${FORCE_BENCHMARK:-0}" == "1" ]]; then
  FORCE_FLAG="--force"
fi
ANALYSIS_PASS_FLAG=()
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  ANALYSIS_PASS_FLAG=(--analysis_pass "${ANALYSIS_PASS}")
fi

if [[ "${METHOD}" == "prepare_pseudobulk" ]]; then
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_prepare_pseudobulk.R"
else
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_benchmark_methods_r.R"
fi

# Unified retry handling: run the method, then grep the Slurm .err for
# transient signatures on failure. No thread pinning (benchmark hardware is
# pinned for runtime comparability). R stderr lands in the Slurm .err file,
# so one transient-signature grep covers it. Method outputs are overwritten
# on re-run (idempotent per-dataset files).
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"

echo "Task ${SLURM_ARRAY_TASK_ID}: running ${METHOD} on ${DS_NAME}"
# PIXI_RSCRIPT word-splits into `pixi run --as-is -e py-cuda13 Rscript
# --vanilla` (--as-is: no lockfile/env mutation from workers); it must stay
# unquoted (established convention, see 2.1_run_worker.sh).
set +e
${PIXI_RSCRIPT} "${R_SCRIPT}" \
    --config_path "${DATASETS_JSON_FILE}" \
    --ds_name "${DS_NAME}" \
    --view "${ANALYSIS_VIEW}" \
    --method "${METHOD}" \
    --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" \
    --results_dir "${RESULTS_DIR}" \
    --pseudobulk_dir "${PSEUDOBULK_DIR}" \
    --gloscope_cache_dir "${GLOSCOPE_DIR}" \
    --log_file "${LOG_FILE}" \
    ${FORCE_FLAG} \
    "${ANALYSIS_PASS_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Task ${SLURM_ARRAY_TASK_ID}: ${METHOD} on ${DS_NAME} complete."
  exit 0
fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_r_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  exit 0   # requeued; the script restarts, likely on another node
fi
exit ${RC}
