#!/bin/bash
set -euo pipefail

# Stage raw data from NAS to HPC scratch (login node has NAS access).
# Run this BEFORE src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py
# and src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh.

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load GCCcore/12.2.0
module load jq/1.6

echo "=== Staging data from NAS to HPC scratch ==="

jq -r '
  to_entries[] |
  .key as $key |
  .value.folder_name as $folder |
  .value.views |
  to_entries[] |
  .value.input_file_name |
  if type == "array" then .[] else . end |
  "\($key) \($folder) \(.)"
' "${DATASETS_JSON_FILE}" | sort -u | \
while read -r KEY FOLDER_NAME RAW_FILE_NAME; do
    if [ "$FOLDER_NAME" == "null" ] || [ -z "$FOLDER_NAME" ]; then
        continue
    fi

    NAS_FILE_PATH="${NAS_SC_DIR}/${FOLDER_NAME}/output/${RAW_FILE_NAME}"
    DEST_DIR="${HPC_SCRATCH_DIR}/${KEY}/data"
    echo "Dataset folder: ${FOLDER_NAME}, file: ${RAW_FILE_NAME}"

    if [ -f "$NAS_FILE_PATH" ]; then
        mkdir -p "${DEST_DIR}"
        rsync -ah --progress "$NAS_FILE_PATH" "${DEST_DIR}/"
    else
        echo "  -> [WARNING] Source not found: ${NAS_FILE_PATH}"
    fi
done

echo "Data staging complete."
