#!/bin/bash

# ==============================================================================
# Script: copy_data_from_nas_to_hpc_scratch.sh
# Purpose: Parse datasets.json and copy corresponding dataset files from NAS to HPC
# ==============================================================================

source "$(dirname "${BASH_SOURCE[0]}")/config.env"

# Check if jq is installed (required for JSON parsing)
module load GCCcore/12.2.0
module load jq/1.6
if ! command -v jq &> /dev/null; then
    echo "Error: 'jq' is not installed or not in PATH."
    echo "Please load it (e.g., 'module load jq' or 'apt-get install jq') and try again."
    exit 1
fi

# Check if datasets.json exists
if [ ! -f "${DATASETS_JSON_FILE}" ]; then
    echo "Error: Cannot find datasets.json at ${DATASETS_JSON_FILE}."
    exit 1
fi

# Create the target scratch directory if it doesn't exist
mkdir -p "${HPC_SCRATCH_DIR}"
echo "Target directory ready: ${HPC_SCRATCH_DIR}"
echo "Starting data transfer..."
echo "------------------------------------------------------"

# Parse JSON and Transfer Files
# jq extracts folder_name and file_name as a space-separated string per dataset
jq -r '.datasets | to_entries[] | "\(.value.folder_name) \(.value.file_name)"' "${DATASETS_JSON_FILE}" | \
while read -r DS_NAME RAW_FILE_NAME; do
    
    # Handle empty/null edge cases safely
    if [ "$DS_NAME" == "null" ] || [ -z "$DS_NAME" ]; then
        continue
    fi

    NAS_FILE_PATH="${NAS_SC_DIR}/${DS_NAME}/output/${RAW_FILE_NAME}"

    echo "Processing dataset: ${DS_NAME}"
    
    # Check if the source file actually exists on the NAS
    if [ -f "$NAS_FILE_PATH" ]; then
        echo "  -> Copying ${RAW_FILE_NAME} to scratch..."
        
        # Using rsync for robust copying with a progress bar. 
        rsync -ah --progress "$NAS_FILE_PATH" "$HPC_SCRATCH_DIR"
        
        echo "  -> Done."
    else
        echo "  -> [WARNING] Source file not found: ${NAS_FILE_PATH}"
    fi
    echo "------------------------------------------------------"
done

echo "All dataset transfers completed successfully!"