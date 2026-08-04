#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=shared-cpu
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
# NOTE: 16G is the baseline. Datasets with >100k cells may need 32-64G.
# If a dataset OOMs, increase --mem for that specific dataset's worker.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load GCCcore/12.2.0
module load jq/1.6
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available on worker node; cannot derive DS_NAME from ${DATASETS_JSON_FILE}."
  exit 1
fi

# Bash arrays do not propagate through sbatch; derive DS_NAME from datasets.json
# directly (jq 'keys[]' is sorted, matching the array indices in 1_submit_hpc_array.sh).
DS_NAME="$(jq -r 'keys[]' "${DATASETS_JSON_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: No dataset for array task ${SLURM_ARRAY_TASK_ID} in ${DATASETS_JSON_FILE}."
  exit 1
fi
echo "Processing dataset: ${DS_NAME} (array task ${SLURM_ARRAY_TASK_ID})"

DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"
OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"

python "${SCRIPT_DIR}/1.1.1_preprocess.py" \
    --config_path "${DATASETS_JSON_FILE}" \
    --input_dir "${DATA_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --ds_name "${DS_NAME}"

echo "Preprocessing complete for ${DS_NAME}"
