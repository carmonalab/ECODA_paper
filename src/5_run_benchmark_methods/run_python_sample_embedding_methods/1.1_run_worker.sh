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
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/../../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage5 \
  "${SCRIPT_DIR}/1.1_run_worker.sh" || exit 1
cd "${PROJECT_ROOT}"
[[ -n "${METHOD:-}" ]] || { echo "ERROR: METHOD is not set." >&2; exit 1; }
MANIFEST_PATH="${ANALYSIS_MANIFEST:-${BENCHMARK_MANIFEST:-}}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: benchmark manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
[[ -n "${line}" ]] || { echo "ERROR: no benchmark row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME ROW_VIEW ROW_LABEL ROW_COMBO ROW_EXTRA <<< "${line}"
[[ -z "${ROW_EXTRA:-}" ]] || {
  echo "ERROR: benchmark row has more than four tab-separated fields" >&2
  exit 1
}
export DS_NAME
ANALYSIS_VIEW="${ROW_VIEW:-${ANALYSIS_VIEW:-benchmark_analysis}}"
export ANALYSIS_VIEW
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
OUT_DIR="${ANALYSIS_ROOT}/embeddings"
EXECUTION_LOG_DIR="${EXECUTION_LOG_DIR:-${OUT_DIR}}"
mkdir -p "${OUT_DIR}" "${EXECUTION_LOG_DIR}"
if [[ -n "${ROW_COMBO:-}" ]]; then
  LOG_SUFFIX="_${ROW_COMBO}"
else
  LOG_SUFFIX=""
fi
if [[ -n "${ANALYSIS_PASS:-}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_batch_effect_${ANALYSIS_PASS}_${METHOD}_${DS_NAME}${LOG_SUFFIX}.feather"
else
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${METHOD}_${DS_NAME}${LOG_SUFFIX}.feather"
fi
FORCE_FLAG=()
[[ "${FORCE_BENCHMARK:-0}" == 1 ]] && FORCE_FLAG=(--force)
ANALYSIS_PASS_FLAG=()
[[ -n "${ANALYSIS_PASS:-}" ]] && ANALYSIS_PASS_FLAG=(--analysis_pass "${ANALYSIS_PASS}")
HIGH_RES_FLAG=()
[[ "${ANALYSIS_HIGH_RES_ONLY:-0}" == 1 ]] && HIGH_RES_FLAG=(--high_resolution_only)
COMBO_FLAG=()
[[ -n "${ROW_COMBO:-}" ]] && COMBO_FLAG=(--combo "${ROW_COMBO}")
# GPU-backed methods must receive an explicit CUDA device.  The Python entry
# point rejects the default "auto" value so a missing runtime export cannot
# silently move a scheduled GPU job onto CPU.
GPU_DEVICE_ARGS=()
case "${ECODA_APPTAINER_NV:-0}" in
  0)
    if [[ "${METHOD}" == "mrvi" && -n "${ROW_COMBO:-}" ]]; then
      GPU_DEVICE_ARGS=(--device cpu)
    fi
    ;;
  1) GPU_DEVICE_ARGS=(--device cuda) ;;
  *) echo "ERROR: ECODA_APPTAINER_NV must be 0 or 1." >&2; exit 1 ;;
esac
PYTHON_ARGS=(
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}" --view "${ANALYSIS_VIEW}"
  --method "${METHOD}" --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output"
  --output_dir "${OUT_DIR}" --log_file "${LOG_FILE}"
)
if [[ ${#GPU_DEVICE_ARGS[@]} -gt 0 ]]; then
  PYTHON_ARGS+=("${GPU_DEVICE_ARGS[@]}")
fi
if [[ ${#FORCE_FLAG[@]} -gt 0 ]]; then
  PYTHON_ARGS+=("${FORCE_FLAG[@]}")
fi
if [[ ${#ANALYSIS_PASS_FLAG[@]} -gt 0 ]]; then
  PYTHON_ARGS+=("${ANALYSIS_PASS_FLAG[@]}")
fi
if [[ ${#HIGH_RES_FLAG[@]} -gt 0 ]]; then
  PYTHON_ARGS+=("${HIGH_RES_FLAG[@]}")
fi
if [[ ${#COMBO_FLAG[@]} -gt 0 ]]; then
  PYTHON_ARGS+=("${COMBO_FLAG[@]}")
fi
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_benchmark_methods_py.py" "${PYTHON_ARGS[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then worker_clear_retry_count; exit 0; fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  rm -f "${LOG_FILE}" "${LOG_FILE}.md5"
  exit 0
fi
exit ${RC}
