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
# This script sits two levels below src/ (slurm_config.sh lives at src/).
source "${SCRIPT_DIR}/../../slurm_config.sh"
cd "${PROJECT_ROOT}"

# METHOD and BENCHMARK_MANIFEST are exported by 1_submit_hpc_array.sh
# (sbatch propagates the submit script's environment).
if [[ -z "${METHOD:-}" ]]; then
  echo "ERROR: METHOD is not set. Export it before submitting the array."
  exit 1
fi
if [[ -z "${BENCHMARK_MANIFEST:-}" ]]; then
  echo "ERROR: BENCHMARK_MANIFEST is not set. Export it before submitting the array."
  exit 1
fi

# Read this task's dataset from the per-method manifest (one DS per line,
# written by 1_submit_hpc_array.sh): sed -n matches SLURM_ARRAY_TASK_ID 1:1.
DS_NAME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BENCHMARK_MANIFEST}")"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No manifest entry for array task ${SLURM_ARRAY_TASK_ID} in ${BENCHMARK_MANIFEST}."
  exit 1
fi
export DS_NAME
echo "Task ${SLURM_ARRAY_TASK_ID}: METHOD=${METHOD}, DS_NAME=${DS_NAME}"

OUT_DIR="${HPC_SCRATCH_DIR}/benchmark/embeddings"
mkdir -p "${OUT_DIR}"
# Log name carries the array JOB id: the mrvi/scpoli/pilot arrays run
# concurrently and all use task ids 1..N — without the job id they would
# read-modify-write the same file. The merge script globs per job id.
LOG_FILE="${OUT_DIR}/execution_times_task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.feather"

# FORCE_BENCHMARK is exported by 1_submit_hpc_array.sh (--force); forward it
# to 1.1.1_benchmark_methods_py.py so existing feathers are recomputed.
FORCE_FLAG=""
if [[ "${FORCE_BENCHMARK:-0}" == "1" ]]; then
  FORCE_FLAG="--force"
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: running ${METHOD} on ${DS_NAME}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_benchmark_methods_py.py" \
    --config_path "${DATASETS_JSON_FILE}" \
    --ds_name "${DS_NAME}" \
    --view benchmark_analysis \
    --method "${METHOD}" \
    --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output" \
    --output_dir "${OUT_DIR}" \
    --log_file "${LOG_FILE}" \
    ${FORCE_FLAG}

echo "Task ${SLURM_ARRAY_TASK_ID}: ${METHOD} on ${DS_NAME} complete."
