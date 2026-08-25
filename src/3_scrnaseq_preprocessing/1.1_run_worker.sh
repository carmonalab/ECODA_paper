#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
# NOTE: 16G is the baseline. Datasets with >100k cells may need 32-64G.
# If a dataset OOMs, increase --mem for that specific dataset's worker.

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

if [[ -z "${PREPROCESS_DATASETS_FILE:-}" ]] && ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available on worker node and PREPROCESS_DATASETS_FILE is unset; cannot derive DS_NAME."
  exit 1
fi

# Durable stage wrappers provide a one-dataset-per-line manifest so selected
# onboarding cohorts can share one array and OOM retries can resubmit only the
# failed rows. Preserve the datasets.json mapping for legacy submitter calls.
if [[ -n "${PREPROCESS_DATASETS_FILE:-}" ]]; then
  if [[ ! -r "${PREPROCESS_DATASETS_FILE}" ]]; then
    echo "ERROR: preprocessing dataset manifest is unreadable: ${PREPROCESS_DATASETS_FILE}"
    exit 1
  fi
  DS_NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PREPROCESS_DATASETS_FILE}")"
else
  # Bash arrays do not propagate through sbatch; derive DS_NAME from datasets.json
  # directly (jq 'keys[]' is sorted, matching the array indices in 1_submit_hpc_array.sh).
  DS_NAME="$(jq -r 'keys[]' "${DATASETS_JSON_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
fi
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No dataset for array task ${SLURM_ARRAY_TASK_ID}."
  exit 1
fi
echo "Processing dataset: ${DS_NAME} (array task ${SLURM_ARRAY_TASK_ID})"

DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"
OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"

# FORCE_PREPROCESS is exported by 1_submit_hpc_array.sh (--force); forward it
# to 1.1.1_preprocess.py so existing outputs are recomputed (debug re-runs).
FORCE_FLAG=""
if [[ "${FORCE_PREPROCESS:-0}" == "1" ]]; then
  FORCE_FLAG="--force"
fi
VIEW_FLAG=()
if [[ -n "${PREPROCESS_VIEW:-}" ]]; then
  VIEW_FLAG=(--view "${PREPROCESS_VIEW}")
fi

# Unified retry handling: transient-failure signatures (stale BeeGFS
# client-cache views, missing imports) trigger a self-requeue, capped by
# WORKER_MAX_RETRIES. No R library staging (python env too large) and no
# thread pinning (preprocessing is multi-threaded by design). Python stderr
# lands in the Slurm .err file, so one transient-signature grep covers it.
source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"

set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_preprocess.py" \
    --config_path "${DATASETS_JSON_FILE}" \
    --input_dir "${DATA_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --ds_name "${DS_NAME}" \
    ${FORCE_FLAG} \
    "${VIEW_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Preprocessing complete for ${DS_NAME}"
  exit 0
fi
if [[ -n "${PREPROCESS_ERROR_PREFIX:-}" ]]; then
  ERR_FILE="${PREPROCESS_ERROR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
else
  ERR_FILE="${LOGS_DIR}/3_scrnaseq_preprocessing_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
fi
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  exit 0   # requeued; the script restarts, likely on another node
fi
exit ${RC}
