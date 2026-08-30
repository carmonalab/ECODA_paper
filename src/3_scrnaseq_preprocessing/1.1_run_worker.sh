#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --time=04:00:00
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
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage3 \
  "${SCRIPT_DIR}/1.1_run_worker.sh" || exit 1
cd "${PROJECT_ROOT}"

MANIFEST_PATH="${PREPROCESS_SELECTION_FILE:-${PREPROCESS_DATASETS_FILE:-}}"
DS_NAME=""
VIEW_NAME=""
if [[ -n "${MANIFEST_PATH}" ]]; then
  [[ -r "${MANIFEST_PATH}" ]] || { echo "ERROR: preprocessing manifest is unreadable: ${MANIFEST_PATH}" >&2; exit 1; }
  manifest_line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST_PATH}")"
  [[ -n "${manifest_line}" ]] || { echo "ERROR: no manifest row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
  IFS=$'\t' read -r DS_NAME VIEW_NAME <<< "${manifest_line}"
else
  command -v jq >/dev/null 2>&1 || { echo "ERROR: jq unavailable and no preprocessing manifest" >&2; exit 1; }
  DS_NAME="$(jq -r 'keys[]' "${DATASETS_JSON_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
  VIEW_NAME="${PREPROCESS_VIEW:-}"
fi
[[ -n "${DS_NAME}" ]] || { echo "ERROR: empty dataset for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
export DS_NAME
export PREPROCESS_VIEW="${VIEW_NAME}"
echo "Processing dataset/view: ${DS_NAME}/${VIEW_NAME:-all} (array task ${SLURM_ARRAY_TASK_ID})"

DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"
OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"
FORCE_FLAG=()
[[ "${FORCE_PREPROCESS:-0}" == "1" ]] && FORCE_FLAG=(--force)
VIEW_FLAG=()
[[ -n "${VIEW_NAME}" ]] && VIEW_FLAG=(--view "${VIEW_NAME}")

source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_preprocess.py" \
  --config_path "${DATASETS_JSON_FILE}" \
  --input_dir "${DATA_DIR}" \
  --output_dir "${OUTPUT_DIR}" \
  --ds_name "${DS_NAME}" \
  "${FORCE_FLAG[@]}" \
  "${VIEW_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  echo "Preprocessing complete for ${DS_NAME}/${VIEW_NAME:-all}"
  exit 0
fi
ERR_PREFIX="${PREPROCESS_ERROR_PREFIX:-${LOGS_DIR}/3_scrnaseq_preprocessing}"
ERR_FILE="${ERR_PREFIX}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
  if [[ -n "${VIEW_NAME}" ]] && command -v jq >/dev/null 2>&1; then
    output_name="$(jq -r --arg ds "${DS_NAME}" --arg view "${VIEW_NAME}" '.[$ds].views[$view].output_file_name // empty' "${DATASETS_JSON_FILE}")"
    [[ -n "${output_name}" ]] && rm -f "${OUTPUT_DIR}/${output_name}" "${OUTPUT_DIR}/${output_name}.md5"
  fi
  exit 0
fi
exit ${RC}
