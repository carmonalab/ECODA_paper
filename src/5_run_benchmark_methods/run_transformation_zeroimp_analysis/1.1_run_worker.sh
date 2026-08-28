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
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"
[[ "${ANALYSIS}" == trans || "${ANALYSIS}" == zeroimp ]] || { echo "ERROR: ANALYSIS must be trans or zeroimp." >&2; exit 1; }
MANIFEST_PATH="${ANALYSIS_MANIFEST:-${BENCHMARK_MANIFEST:-}}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: analysis manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
[[ -n "${line}" ]] || { echo "ERROR: no analysis row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME ROW_VIEW ROW_LABEL <<< "${line}"
export DS_NAME
ANALYSIS_VIEW="${ROW_VIEW:-benchmark_analysis}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
OUT_DIR="${ANALYSIS_ROOT}/embeddings"
EXECUTION_LOG_DIR="${EXECUTION_LOG_DIR:-${OUT_DIR}}"
mkdir -p "${OUT_DIR}" "${EXECUTION_LOG_DIR}"
LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${ANALYSIS}_${DS_NAME}.feather"
FORCE_FLAG=()
[[ "${FORCE_BENCHMARK:-0}" == 1 ]] && FORCE_FLAG=(--force)
if [[ "${ANALYSIS}" == trans ]]; then R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_transformation_analysis.R"; else R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_zeroimp_analysis.R"; fi
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
set +e
${PIXI_RSCRIPT} "${R_SCRIPT}" \
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}" --view "${ANALYSIS_VIEW}" \
  --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" --output_dir "${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}/results" \
  --log_file "${LOG_FILE}" "${FORCE_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then worker_clear_retry_count; exit 0; fi
ERR_FILE="${LOGS_DIR}/5_transzeroimp_${ANALYSIS}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
