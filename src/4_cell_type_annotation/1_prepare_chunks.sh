#!/bin/bash
#
# 1_prepare_chunks.sh — Build dataset chunks (Supports production or test modes)
#
# Usage:
#   ./1_prepare_chunks.sh                        # production, all datasets (chunk_size = 2)
#   ./1_prepare_chunks.sh test                   # test mode, all datasets (chunk_size = 1)
#   ./1_prepare_chunks.sh production <DS_NAME>   # production, single dataset
#   ./1_prepare_chunks.sh test <DS_NAME>         # test mode, single dataset
#

set -euo pipefail

# 1. Load central config
source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${LOGS_DIR}"

# 2. Parse arguments: MODE (production/test) + optional single dataset
MODE_ARG="${1:-production}"
DS_NAME_ARG="${2:-}"
PY_ARGS=""

if [ "$MODE_ARG" = "test" ]; then
  echo ">>> CONFIGURING PIPELINE IN TEST MODE <<<"
  PY_ARGS="--test"
elif [ "$MODE_ARG" = "production" ]; then
  echo ">>> CONFIGURING PIPELINE IN PRODUCTION MODE <<<"
else
  echo "ERROR: Unknown mode '${MODE_ARG}'. Use 'production' or 'test'."
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi

# 3. Build the list of datasets (all keys of datasets.json, or a single
#    validated key when DS_NAME_ARG is given)
DATASET_NAMES=()
if [[ -n "${DS_NAME_ARG}" ]]; then
  if ! jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}" >/dev/null 2>&1; then
    echo "ERROR: '${DS_NAME_ARG}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
  DATASET_NAMES+=("${DS_NAME_ARG}")
else
  while IFS= read -r name; do
    DATASET_NAMES+=("$name")
  done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")
fi

if [[ ${#DATASET_NAMES[@]} -eq 0 ]]; then
  echo "ERROR: No datasets found in ${DATASETS_JSON_FILE}."
  exit 1
fi
echo "Datasets to process: ${DATASET_NAMES[*]}"

# -------------------------------------------------------------------------
# STAGE REFERENCE MAPS: Copy ref files from UNIGE NAS to Cluster Scratch
# -------------------------------------------------------------------------

# Reference maps for cell type annotation (light reference atlases in blood
# and tumors, Garnica et al. 2024; carmonalab/Reference_maps):
# https://doi.org/10.6084/m9.figshare.26310994
# Staged NAS-first with a Figshare download fallback. Every file is MD5-
# verified against the manifest regardless of source (NAS, Figshare, or
# pre-existing); fails closed if any of the 4 files is missing or corrupt
# after both attempts.
REF_MAP_NAMES=(
  "sketched_CD8T_human_ref_v1.rds"
  "sketched_CD4T_human_ref_v2.rds"
  "sketched_DC_human_ref_v2.rds"
  "sketched_MoMac_human_v1.rds"
)
declare -A REF_MAP_FILE_IDS=(
  ["sketched_CD8T_human_ref_v1.rds"]="47714158"
  ["sketched_CD4T_human_ref_v2.rds"]="47714155"
  ["sketched_DC_human_ref_v2.rds"]="47714161"
  ["sketched_MoMac_human_v1.rds"]="47714164"
)
declare -A REF_MAP_MD5S=(
  ["sketched_CD8T_human_ref_v1.rds"]="be86058ddafdd0154faf0485286b86e7"
  ["sketched_CD4T_human_ref_v2.rds"]="5540a0ee287e291528c96d476794b194"
  ["sketched_DC_human_ref_v2.rds"]="033d491ba7ca9bbf0badcae828e55b2c"
  ["sketched_MoMac_human_v1.rds"]="3043cd9058a8746d972c7be195b18e36"
)

mkdir -p "${HOME_REF_DIR}"

# Sweep temp files orphaned by interrupted runs (login-node SSH drops; the
# PID-suffixed temp names are never re-targeted by later runs).
rm -f "${HOME_REF_DIR}"/.*.tmp.*

ref_map_md5_ok() {
  local f="$1" path="$2"
  [[ "$(md5sum "${path}" | awk '{print $1}')" == "${REF_MAP_MD5S[$f]}" ]]
}

ref_maps_staged() {
  local f
  for f in "${REF_MAP_NAMES[@]}"; do
    [[ -f "${HOME_REF_DIR}/${f}" ]] || return 1
    ref_map_md5_ok "$f" "${HOME_REF_DIR}/${f}" || return 1
  done
  return 0
}

drop_bad_ref_files() {
  local f
  for f in "${REF_MAP_NAMES[@]}"; do
    if [[ -e "${HOME_REF_DIR}/${f}" ]] && ! ref_map_md5_ok "$f" "${HOME_REF_DIR}/${f}"; then
      echo "WARNING: ${HOME_REF_DIR}/${f} failed MD5 verification; removing for re-staging."
      rm -f "${HOME_REF_DIR}/${f}"
    fi
  done
}

if ref_maps_staged; then
  echo ">>> Reference maps already staged and MD5-verified in ${HOME_REF_DIR}. Skipping staging. <<<"
else
  drop_bad_ref_files
  echo "Staging reference maps from NAS to home directory..."
  if rsync -a "${NAS_REF_DIR}" "${HOME_REF_DIR}/"; then
    echo "Reference maps staged from NAS."
  else
    echo "WARNING: NAS rsync of reference maps failed (${NAS_REF_DIR} unavailable?); falling back to Figshare download."
  fi
  drop_bad_ref_files

  for f in "${REF_MAP_NAMES[@]}"; do
    if [[ -f "${HOME_REF_DIR}/${f}" ]]; then
      continue
    fi
    url="https://ndownloader.figshare.com/files/${REF_MAP_FILE_IDS[$f]}"
    tmp="${HOME_REF_DIR}/.${f}.tmp.$$"
    echo "Downloading ${f} from Figshare (${url})..."
    if ! curl -f -L --retry 3 -o "${tmp}" "${url}"; then
      rm -f "${tmp}"
      echo "ERROR: Download of ${f} from Figshare failed."
      exit 1
    fi
    if ! ref_map_md5_ok "$f" "${tmp}"; then
      rm -f "${tmp}"
      echo "ERROR: MD5 checksum mismatch for ${f}; expected ${REF_MAP_MD5S[$f]}."
      exit 1
    fi
    mv "${tmp}" "${HOME_REF_DIR}/${f}"
  done

  if ! ref_maps_staged; then
    echo "ERROR: Reference maps incomplete or failing MD5 checks in ${HOME_REF_DIR}; need all of: ${REF_MAP_NAMES[*]}."
    exit 1
  fi
fi

# 4. Build chunks per dataset (sequential, one short-lived compute session each;
#    --time=00:30:00 covers the full dataset loop, per-dataset runs are usually
#    much shorter)
FAILED_DATASETS=()
SKIPPED_DATASETS=()

for DS_NAME in "${DATASET_NAMES[@]}"; do
  echo ""
  echo ">>> Building chunks for dataset: ${DS_NAME} <<<"

  # Skip datasets without preprocessed .h5ad input (e.g. Zhu has no views)
  if ! ls "${HPC_SCRATCH_DIR}/${DS_NAME}/output"/*.h5ad >/dev/null 2>&1; then
    echo "WARNING: No preprocessed .h5ad files in ${HPC_SCRATCH_DIR}/${DS_NAME}/output; skipping ${DS_NAME}."
    SKIPPED_DATASETS+=("${DS_NAME}")
    continue
  fi

  LOG_FILE="${LOGS_DIR}/prepare_chunks_${MODE_ARG}_${DS_NAME}.log"
  export DS_NAME

  echo "Allocating short-lived compute session to build ${DS_NAME} chunks..."
  if ! srun --partition="${SLURM_PARTITION}" \
       --time=00:30:00 \
       --ntasks=1 \
       --cpus-per-task=1 \
       --mem=4G \
       --output="${LOG_FILE}" \
       --error="${LOG_FILE}" \
       "${PYTHON_BIN}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" ${PY_ARGS}; then
    echo "ERROR: Chunk generation failed for ${DS_NAME}. See ${LOG_FILE}."
    FAILED_DATASETS+=("${DS_NAME}")
    continue
  fi

  echo "✓ Chunks generated for ${DS_NAME}. Log saved to: ${LOG_FILE}"
done

echo ""
echo "=== Chunk preparation summary ==="
echo "Processed: ${#DATASET_NAMES[@]} datasets"
echo "Skipped (no preprocessed .h5ad): ${#SKIPPED_DATASETS[@]}"
if [[ ${#SKIPPED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${SKIPPED_DATASETS[*]}"
fi
echo "Failed: ${#FAILED_DATASETS[@]}"
if [[ ${#FAILED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${FAILED_DATASETS[*]}"
  echo "Resolve the failures, then rerun this script before 2_submit_hpc_array.sh."
  exit 1
fi
