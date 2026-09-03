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
[[ -n "${DS_NAME}" ]] || { echo "ERROR: malformed benchmark row" >&2; exit 1; }
[[ -z "${ROW_EXTRA:-}" ]] || {
  echo "ERROR: benchmark row has more than four tab-separated fields" >&2
  exit 1
}
if [[ -n "${ROW_COMBO:-}" ]]; then
  [[ "${METHOD}" == gloscope ]] || {
    echo "ERROR: combo token is only supported for method gloscope" >&2
    exit 1
  }
  [[ -z "${ANALYSIS_PASS:-}" ]] || {
    echo "ERROR: GloScope combo shards are only supported for ordinary runs" >&2
    exit 1
  }
  case "${ROW_COMBO}" in
    hvg2000_pcadims10|hvg2000_pcadims30|hvg2000_pcadims50|\
      hvg1000_pcadims30|hvg3000_pcadims30) ;;
    *)
      echo "ERROR: invalid GloScope combo token: ${ROW_COMBO}" >&2
      exit 1
      ;;
  esac
fi
export DS_NAME
ANALYSIS_VIEW="${ROW_VIEW:-${ANALYSIS_VIEW:-benchmark_analysis}}"
export ANALYSIS_VIEW
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${HPC_SCRATCH_DIR}/benchmark}"
RESULTS_DIR="${ANALYSIS_ROOT}/results"
PSEUDOBULK_DIR="${ANALYSIS_ROOT}/pseudobulks"
GLOSCOPE_DIR="${ANALYSIS_ROOT}/gloscope_dists"
EMBED_DIR="${ANALYSIS_ROOT}/embeddings"
EXECUTION_LOG_DIR="${EXECUTION_LOG_DIR:-${EMBED_DIR}}"
mkdir -p "${RESULTS_DIR}" "${PSEUDOBULK_DIR}" "${GLOSCOPE_DIR}" "${EMBED_DIR}" "${EXECUTION_LOG_DIR}"
if [[ -n "${ROW_COMBO:-}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${METHOD}_${DS_NAME}_${ROW_COMBO}.feather"
elif [[ -n "${ANALYSIS_PASS:-}" ]]; then
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_batch_effect_${ANALYSIS_PASS}_${METHOD}_${DS_NAME}.feather"
else
  LOG_FILE="${EXECUTION_LOG_DIR}/execution_times_${METHOD}_${DS_NAME}.feather"
fi
FORCE_FLAG=()
[[ "${FORCE_BENCHMARK:-0}" == 1 ]] && FORCE_FLAG=(--force)
ANALYSIS_PASS_FLAG=()
[[ -n "${ANALYSIS_PASS:-}" ]] && ANALYSIS_PASS_FLAG=(--analysis_pass "${ANALYSIS_PASS}")
COMBO_FLAG=()
[[ -n "${ROW_COMBO:-}" ]] && COMBO_FLAG=(--combo "${ROW_COMBO}")
if [[ "${METHOD}" == prepare_pseudobulk ]]; then
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_prepare_pseudobulk.R"
else
  R_SCRIPT="${SCRIPT_DIR}/1.1.1_run_benchmark_methods_r.R"
fi
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
set +e
${PIXI_RSCRIPT} "${R_SCRIPT}" \
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}" --view "${ANALYSIS_VIEW}" \
  --method "${METHOD}" --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" \
  --results_dir "${RESULTS_DIR}" --pseudobulk_dir "${PSEUDOBULK_DIR}" \
  --gloscope_cache_dir "${GLOSCOPE_DIR}" --log_file "${LOG_FILE}" \
  "${FORCE_FLAG[@]}" "${ANALYSIS_PASS_FLAG[@]}" "${COMBO_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then worker_clear_retry_count; exit 0; fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_r_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
