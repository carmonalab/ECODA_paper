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

DATASET_NAMES=()
while IFS= read -r name; do
  DATASET_NAMES+=("$name")
done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")

IDX=$((SLURM_ARRAY_TASK_ID - 1))
DS_NAME="${DATASET_NAMES[$IDX]}"
echo "Processing dataset: ${DS_NAME} (array task ${SLURM_ARRAY_TASK_ID})"

DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"
OUTPUT_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}"
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"

python "${SCRIPT_DIR}/1.1.1_preprocess.py" \
    --config_path "${DATASETS_JSON_FILE}" \
    --base_path "${DATA_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --ds_name "${DS_NAME}"

echo "Preprocessing complete for ${DS_NAME}"
