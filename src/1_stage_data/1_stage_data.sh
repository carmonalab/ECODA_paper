#!/bin/bash
set -euo pipefail

# Stage raw data from NAS to HPC scratch (login node has NAS access).
# Run this BEFORE src/2_dataset_specific_preprocessing/1_submit_hpc.sh
# and src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh.
#
# Usage:
#   ./1_stage_data.sh                  # stage all datasets (skips keys starting with _)
#   ./1_stage_data.sh --ds_name _debug # stage only the _debug dataset
#
# Convention: default-all loops skip datasets whose key starts with "_" (e.g.
# the _debug entry) unless explicitly requested via --ds_name. The _debug raw
# files must exist on the NAS under
# Standardized_SingleCell_Datasets/debug/output/ (user places them there after
# running src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R),
# or be copied to ${HPC_SCRATCH_DIR}/_debug/data/ manually.

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"

module load GCCcore/12.2.0
module load jq/1.6

DS_NAME_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ds_name)
      DS_NAME_ARG="${2:-}"
      shift 2
      ;;
    --ds_name=*)
      DS_NAME_ARG="${1#*=}"
      shift
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  echo "=== Staging data from NAS to HPC scratch (dataset: ${DS_NAME_ARG}) ==="
else
  echo "=== Staging data from NAS to HPC scratch (all datasets; '_' keys skipped) ==="
fi

# Datasets with views emit each view's input_file_name; view-less datasets
# (e.g. Zhu) fall back to the dataset-level file_names (string or array).
# The later `sort -u` dedups files shared across views.
# When --ds_name is given, only that dataset's files are emitted; otherwise
# keys starting with "_" are skipped.
jq -r --arg only_ds "${DS_NAME_ARG}" '
  to_entries[] |
  (if $only_ds != "" then select(.key == $only_ds) else select(.key | startswith("_") | not) end) |
  .key as $key |
  .value.folder_name as $folder |
  (if (.value.views | length) > 0 then .value.views | to_entries[] | .value.input_file_name else .value.file_names end) |
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
