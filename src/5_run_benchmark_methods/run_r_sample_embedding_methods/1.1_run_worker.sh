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
WORKER_SOURCE_RELATIVE="src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh"
SOURCE_REQUIRED="${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}"
SOURCE_ROOT="${ECODA_SOURCE_ROOT:-}"
SOURCE_MANIFEST="${ECODA_SOURCE_MANIFEST:-}"
RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-}"
RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-}"
RUN_ID="${ECODA_RUN_ID:-}"
RUN_ROOT="${ECODA_RUN_ROOT:-}"
[[ "${SOURCE_REQUIRED}" == "1" && "${SOURCE_ROOT}" = /* &&
   "${SOURCE_MANIFEST}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires an immutable source snapshot." >&2
  exit 1
}
[[ -n "${RUNTIME_IMAGE}" && "${RUNTIME_IMAGE}" = /* &&
   -n "${RUNTIME_MANIFEST}" && "${RUNTIME_MANIFEST}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires a bound runtime image and manifest." >&2
  exit 1
}
[[ -n "${RUN_ID}" && "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ &&
   -n "${RUN_ROOT}" && "${RUN_ROOT}" = /* ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 worker requires a valid bound run root and ID." >&2
  exit 1
}
[[ -d "${SOURCE_ROOT}" && -f "${SOURCE_MANIFEST}" &&
   -r "${SOURCE_MANIFEST}" && ! -L "${SOURCE_MANIFEST}" ]] || {
  echo "ERROR: immutable source manifest is missing or unreadable: ${SOURCE_MANIFEST}" >&2
  exit 1
}
SNAPSHOT_WORKER="${SOURCE_ROOT%/}/${WORKER_SOURCE_RELATIVE}"
[[ -f "${SNAPSHOT_WORKER}" && -r "${SNAPSHOT_WORKER}" && ! -L "${SNAPSHOT_WORKER}" ]] || {
  echo "ERROR: immutable Stage 5 worker is missing or unreadable: ${SNAPSHOT_WORKER}" >&2
  exit 1
}
command -v realpath >/dev/null 2>&1 || {
  echo "ERROR: realpath is required for immutable Stage 5 workers." >&2
  exit 1
}
SNAPSHOT_WORKER_REAL="$(realpath "${SNAPSHOT_WORKER}" 2>/dev/null || true)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  command -v scontrol >/dev/null 2>&1 || {
    echo "ERROR: scontrol is required to verify the immutable worker command." >&2
    exit 1
  }
  SCHEDULED_COMMAND="$(scontrol show job "${SLURM_JOB_ID}" -o 2>/dev/null || true)"
  SCHEDULED_SCRIPT="${SCHEDULED_COMMAND#*Command=}"
  SCHEDULED_SCRIPT="${SCHEDULED_SCRIPT%% *}"
  SCHEDULED_SCRIPT_REAL="$(realpath "${SCHEDULED_SCRIPT}" 2>/dev/null || true)"
  [[ -n "${SCHEDULED_SCRIPT_REAL}" && "${SCHEDULED_SCRIPT_REAL}" == "${SNAPSHOT_WORKER_REAL}" ]] || {
    echo "ERROR: mutable or non-snapshot Stage 5 worker command rejected." >&2
    exit 1
  }
else
  CURRENT_SCRIPT_REAL="$(realpath "${BASH_SOURCE[0]}" 2>/dev/null || true)"
  [[ -n "${CURRENT_SCRIPT_REAL}" && "${CURRENT_SCRIPT_REAL}" == "${SNAPSHOT_WORKER_REAL}" ]] || {
    echo "ERROR: mutable or non-snapshot Stage 5 worker path rejected." >&2
    exit 1
  }
fi
source "${SOURCE_ROOT%/}/src/slurm_config.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUNTIME_IMAGE="${RUNTIME_IMAGE}"
export ECODA_RUNTIME_MANIFEST="${RUNTIME_MANIFEST}"
export ECODA_RUN_ID="${RUN_ID}"
export ECODA_RUN_ROOT="${RUN_ROOT}"
EXPECTED_RUN_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
[[ "${RUN_ROOT}" == "${EXPECTED_RUN_ROOT}" && -d "${RUN_ROOT}" && ! -L "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 5 worker run root is not the exact bound run root: ${RUN_ROOT}" >&2
  exit 1
}
RUN_SOURCE_MANIFEST="${RUN_ROOT}/manifests/source.manifest"
RUN_RUNTIME_IDENTITY="${RUN_ROOT}/manifests/runtime.identity"
[[ -f "${RUN_SOURCE_MANIFEST}" && -r "${RUN_SOURCE_MANIFEST}" && ! -L "${RUN_SOURCE_MANIFEST}" &&
   -f "${RUN_RUNTIME_IDENTITY}" && -r "${RUN_RUNTIME_IDENTITY}" && ! -L "${RUN_RUNTIME_IDENTITY}" ]] || {
  echo "ERROR: legacy_source_unpinned: Stage 5 run-bound source/runtime manifests are missing." >&2
  exit 1
}
cmp -s "${SOURCE_MANIFEST}" "${RUN_SOURCE_MANIFEST}" || {
  echo "ERROR: Stage 5 run source manifest differs from the immutable source manifest." >&2
  exit 1
}
ecoda_runtime_reexec_worker stage5 \
  "${SNAPSHOT_WORKER}" || exit 1
SCRIPT_DIR="${SOURCE_ROOT%/}/src/5_run_benchmark_methods/run_r_sample_embedding_methods"
source "${SOURCE_ROOT%/}/src/slurm_config.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_runtime.sh"
source "${SOURCE_ROOT%/}/src/utils/bash/ecoda_run_common.sh"
export ECODA_SOURCE_ROOT="${SOURCE_ROOT}"
export ECODA_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export ECODA_RUN_ID="${RUN_ID}"
export ECODA_RUN_ROOT="${RUN_ROOT}"
if [[ "${ECODA_RUNTIME_MODE:-host}" == "host" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  ecoda_runtime_validate_bound_run || {
    echo "ERROR: Stage 5 run-bound runtime validation failed." >&2
    exit 1
  }
fi
ecoda_validate_run_id "${ECODA_RUN_ID}" || exit 1
[[ "${ECODA_RUN_ROOT}" == "${RUN_ROOT}" ]] || {
  echo "ERROR: Stage 5 run root changed during immutable worker reexec." >&2
  exit 1
}
cd "${PROJECT_ROOT}"
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
[[ -n "${METHOD:-}" ]] || { echo "ERROR: METHOD is not set." >&2; exit 1; }
case "${METHOD}" in
  gloscope|mofa|pseudobulk|scitd|composition|trans|zeroimp)
    export ECODA_ARTIFACT_PRODUCER="stage5_${METHOD}"
    ;;
esac
MANIFEST_PATH="${MATRIX_RETRY_MANIFEST:-${ANALYSIS_MANIFEST:-${BENCHMARK_MANIFEST:-}}}"
[[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: benchmark manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
ecoda_validate_run_owned_path "${MANIFEST_PATH}" "${ECODA_RUN_ROOT}" || {
  echo "ERROR: benchmark manifest is outside the bound Stage 5 run root: ${MANIFEST_PATH}" >&2
  exit 1
}
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
set +e
R_ARGS=(
  "${R_SCRIPT}"
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}" --view "${ANALYSIS_VIEW}"
  --method "${METHOD}" --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output"
  --results_dir "${RESULTS_DIR}" --pseudobulk_dir "${PSEUDOBULK_DIR}"
  --gloscope_cache_dir "${GLOSCOPE_DIR}" --log_file "${LOG_FILE}"
)
if [[ ${#FORCE_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${FORCE_FLAG[@]}"); fi
if [[ ${#ANALYSIS_PASS_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${ANALYSIS_PASS_FLAG[@]}"); fi
if [[ ${#COMBO_FLAG[@]} -gt 0 ]]; then R_ARGS+=("${COMBO_FLAG[@]}"); fi
${PIXI_RSCRIPT} "${R_ARGS[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then worker_clear_retry_count; exit 0; fi
ERR_PREFIX="${JOB_LOG_PREFIX:-5_benchmark_r_${METHOD}}"
ERR_FILE="${LOGS_DIR}/${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
