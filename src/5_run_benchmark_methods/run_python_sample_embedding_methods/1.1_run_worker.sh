#!/bin/bash
#SBATCH --job-name=benchmark_worker
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
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

# Read this task's dataset from the per-method manifest.
DS_NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No manifest entry for array task ${SLURM_ARRAY_TASK_ID} in ${MANIFEST_PATH}."
  exit 1
fi
export DS_NAME
echo "Task ${SLURM_ARRAY_TASK_ID}: METHOD=${METHOD}, DS_NAME=${DS_NAME}"

ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
ANALYSIS_VIEW="${ANALYSIS_VIEW:-benchmark_analysis}"
OUT_DIR="${ANALYSIS_ROOT}/embeddings"
mkdir -p "${OUT_DIR}"
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  LOG_FILE="${OUT_DIR}/execution_times_batch_effect_${ANALYSIS_PASS}_${METHOD}_${DS_NAME}.feather"
else
  LOG_FILE="${OUT_DIR}/execution_times_${METHOD}_${DS_NAME}.feather"
fi

# FORCE_BENCHMARK is exported by 1_submit_hpc_array.sh (--force); forward it
# to 1.1.1_benchmark_methods_py.py so existing feathers are recomputed.
FORCE_FLAG=""
if [[ "${FORCE_BENCHMARK:-0}" == "1" ]]; then
  FORCE_FLAG="--force"
fi
ANALYSIS_PASS_FLAG=()
HIGH_RES_FLAG=()
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  ANALYSIS_PASS_FLAG=(--analysis_pass "${ANALYSIS_PASS}")
fi
if [[ "${ANALYSIS_HIGH_RES_ONLY:-0}" == "1" ]]; then
  HIGH_RES_FLAG=(--high_resolution_only)
fi

# Unified retry handling: transient-failure signatures (stale BeeGFS
# client-cache views, missing imports) trigger a self-requeue, capped by
# WORKER_MAX_RETRIES. No R library staging (python env too large), no thread
# pinning (benchmark hardware is pinned for runtime comparability). Python
# stderr lands in the Slurm .err file, so one transient-signature grep covers
# it. Method outputs are overwritten on re-run (idempotent per-dataset files).
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"

echo "Task ${SLURM_ARRAY_TASK_ID}: running ${METHOD} on ${DS_NAME}"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_benchmark_methods_py.py" \
    --config_path "${DATASETS_JSON_FILE}" \
    --ds_name "${DS_NAME}" \
    --view "${ANALYSIS_VIEW}" \
    --method "${METHOD}" \
    --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" \
    --output_dir "${OUT_DIR}" \
    --log_file "${LOG_FILE}" \
    ${FORCE_FLAG} \
    "${ANALYSIS_PASS_FLAG[@]}" \
    "${HIGH_RES_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Task ${SLURM_ARRAY_TASK_ID}: ${METHOD} on ${DS_NAME} complete."
  exit 0
fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  # pyarrow to_feather is NOT atomic: a partial feather from the failed
  # attempt would be skipped as "already processed" by the combo skip-check
  # (out_path.exists() in 1.1.1_benchmark_methods_py.py) on the re-run.
  # Delete this task's per-dataset outputs (embeddings/dists + exec log all
  # embed ${DS_NAME}) so the re-run recomputes them fresh. Filenames are
  # per (method, dataset), so other tasks' files are untouched.
  rm -f "${OUT_DIR}"/*"${DS_NAME}"*.feather
  exit 0   # requeued; the script restarts, likely on another node
fi
exit ${RC}
