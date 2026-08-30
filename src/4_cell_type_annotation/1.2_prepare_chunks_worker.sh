#!/bin/bash
#SBATCH --job-name=annotation_prepare
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SLURM_JOB_ID:-}" &&
      "${ECODA_RUNTIME_IN_CONTAINER:-0}" != "1" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
source "${SCRIPT_DIR}/../slurm_config.sh"
source "${SCRIPT_DIR}/../utils/bash/ecoda_runtime.sh"
ecoda_runtime_reexec_worker stage4 \
  "${SCRIPT_DIR}/1.2_prepare_chunks_worker.sh" || exit 1
source "${SCRIPT_DIR}/../utils/bash/ecoda_run_common.sh"
cd "${PROJECT_ROOT}"

MANIFEST="${ANNOTATION_PREP_MANIFEST:-}"
[[ -r "${MANIFEST}" ]] || { echo "ERROR: annotation prep manifest is unreadable: ${MANIFEST}" >&2; exit 1; }
RUN_ID="${ANNOTATION_RUN_ID:-}"
ecoda_validate_run_id "${RUN_ID}" || { echo "ERROR: ANNOTATION_RUN_ID is invalid or missing." >&2; exit 1; }
EXPECTED_ROOT="${HPC_SCRATCH_DIR}/_ecoda_runs/${RUN_ID}"
ecoda_validate_run_owned_path "${MANIFEST}" "${EXPECTED_ROOT}" ||
  { echo "ERROR: annotation prep manifest is outside the Stage 4 run root." >&2; exit 1; }
line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
[[ -n "${line}" ]] || { echo "ERROR: no annotation prep row for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IFS=$'\t' read -r DS_NAME VIEWS RUN_ROOT <<< "${line}"
[[ -n "${DS_NAME}" && -n "${VIEWS}" && -n "${RUN_ROOT}" ]] || { echo "ERROR: malformed annotation prep row" >&2; exit 1; }
expected_root_real="$(ecoda_realpath_existing "${EXPECTED_ROOT}" 2>/dev/null || true)"
run_root_real="$(ecoda_realpath_existing "${RUN_ROOT}" 2>/dev/null || true)"
[[ -n "${expected_root_real}" && "${expected_root_real}" == "${run_root_real}" ]] ||
  { echo "ERROR: annotation prep row run root is not canonical." >&2; exit 1; }
export DS_NAME ANNOTATION_VIEWS="${VIEWS}" ANNOTATION_RUN_ROOT="${RUN_ROOT}" ANNOTATION_RUN_ID="${RUN_ID}"
FORCE_FLAG=()
[[ "${FORCE_ANNOTATION:-0}" == "1" ]] && FORCE_FLAG=(--force)

source "${SCRIPT_DIR}/../utils/bash/worker_retry.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" \
  --views "${VIEWS}" --run-root "${RUN_ROOT}" "${FORCE_FLAG[@]}"
RC=$?
set -e
if [[ ${RC} -eq 0 ]]; then
  worker_clear_retry_count
  exit 0
fi
ERR_FILE="${ANNOTATION_PREP_ERROR_PREFIX:-${LOGS_DIR}/4_annotation_prepare}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then exit 0; fi
exit ${RC}
