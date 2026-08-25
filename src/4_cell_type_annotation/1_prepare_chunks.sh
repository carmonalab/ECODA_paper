#!/bin/bash
#
# 1_prepare_chunks.sh — Build dataset chunks (Supports production or test modes)
# Usage:
#   ./1_prepare_chunks.sh                        # production, all datasets (chunk_size = 2)
#   ./1_prepare_chunks.sh test                   # test mode, all datasets (chunk_size = 1)
#   ./1_prepare_chunks.sh production <DS_NAME>   # production, single dataset
#   [--view NAME]                                # only require this preprocessed view
#   [--force]                                    # recompute chunks even if the dataset is
#                                                # already annotated (accepted in any position)
# Datasets that are already annotated are skipped unless --force is given —
# but only after a clean-entry check (1.1_prepare_chunks.py --check-clean)
# confirms the views carry no legacy annotation columns; flagged datasets are
# rebuilt for re-annotation (see the per-dataset loop below).
#

set -euo pipefail

# 1. Load central config
source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${LOGS_DIR}"

# 2. Parse arguments: MODE (production/test) + optional single dataset,
#    optional view, and optional --force flag (accepted in any position).
FORCE_ARG=0
VIEW_ARG=""
POS_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)
      FORCE_ARG=1
      shift
      ;;
    --view)
      VIEW_ARG="${2:-}"
      if [[ -z "${VIEW_ARG}" ]]; then
        echo "ERROR: --view requires a view name."
        exit 1
      fi
      shift 2
      ;;
    --view=*)
      VIEW_ARG="${1#*=}"
      if [[ -z "${VIEW_ARG}" ]]; then
        echo "ERROR: --view requires a view name."
        exit 1
      fi
      shift
      ;;
    -*)
      echo "ERROR: Unknown argument: $1"
      exit 1
      ;;
    *)
      POS_ARGS+=("$1")
      shift
      ;;
  esac
done

MODE_ARG="${POS_ARGS[0]:-production}"
DS_NAME_ARG="${POS_ARGS[1]:-}"
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
SKIPPED_ANNOTATED=()
SKIPPED_INCOMPLETE=()
FLAGGED_LEGACY=()

for DS_NAME in "${DATASET_NAMES[@]}"; do
  echo ""
  echo ">>> Building chunks for dataset: ${DS_NAME} <<<"

  # Skip datasets that are already annotated, unless --force. Two done states:
  #   Branch 1 — annotated, not yet merged: >=1 chunk_*.txt in output/chunks/
  #     and every chunk_N.txt has its matching annotations_chunk_N.feather in
  #     output/ (mapping: chunk_3.txt -> annotations_chunk_3.feather).
  #   Branch 2 — already merged (3_submit_merge.sh is the only step that
  #     deletes these dirs): output/chunks/ absent AND annotation_union/ absent
  #     AND >=1 annotations_chunk_*.feather exists in output/.
  # Partial feather coverage (some chunks missing their feather) is the "not
  # done" state and falls through to the rebuild below.
  shopt -s nullglob
  CHUNK_FILES=("${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks"/chunk_*.txt)
  ANNOT_FILES=("${HPC_SCRATCH_DIR}/${DS_NAME}/output"/annotations_chunk_*.feather)
  shopt -u nullglob

  ANNOTATED=0
  if [[ ${#CHUNK_FILES[@]} -gt 0 ]]; then
    ANNOTATED=1
    for chunk_file in "${CHUNK_FILES[@]}"; do
      chunk_name=$(basename "${chunk_file}")
      chunk_num="${chunk_name#chunk_}"
      chunk_num="${chunk_num%.txt}"
      if [[ ! -f "${HPC_SCRATCH_DIR}/${DS_NAME}/output/annotations_chunk_${chunk_num}.feather" ]]; then
        ANNOTATED=0
        break
      fi
    done
  elif [[ ! -d "${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union" && ${#ANNOT_FILES[@]} -gt 0 ]]; then
    ANNOTATED=1
  fi

  if [[ ${ANNOTATED} -eq 1 ]]; then
    if [[ ${FORCE_ARG} -eq 1 ]]; then
      echo "NOTE: ${DS_NAME} is already annotated; --force given, rebuilding chunks."
    else
      # Clean-entry check: an annotated dataset may only be skipped if every
      # tier-matching obs column of its views was produced by THIS pipeline
      # (present in its own annotation feathers). Legacy columns mean a
      # previous annotation leaked into the views -> NOT clean -> rebuild
      # chunks and re-annotate (the worker wipe + merge tiered drop then
      # scrub them). Never use --force semantics to circumvent.
      # DS_NAME is exported here too (not only at the chunk-generation srun
      # below): 1.1_prepare_chunks.py requires it from the environment.
      export DS_NAME
      echo "Already annotated: ${DS_NAME} — running clean-entry check..."
      CHECK_LOG_FILE="${LOGS_DIR}/prepare_chunks_${MODE_ARG}_${DS_NAME}_check.log"
      set +e
      srun --partition="${SLURM_PARTITION}" \
           --time=00:30:00 \
           --ntasks=1 \
           --cpus-per-task=1 \
           --mem=4G \
           --output="${CHECK_LOG_FILE}" \
           --error="${CHECK_LOG_FILE}" \
           "${PYTHON_BIN}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" --check-clean
      CHECK_EXIT=$?
      set -e
      if [[ ${CHECK_EXIT} -eq 2 ]]; then
        echo "WARNING: ${DS_NAME} carries legacy annotation columns (see ${CHECK_LOG_FILE}); rebuilding chunks for re-annotation."
        FLAGGED_LEGACY+=("${DS_NAME}")
        ANNOTATED=0
      elif [[ ${CHECK_EXIT} -eq 0 ]]; then
        echo "Clean: ${DS_NAME} — skipping chunk generation (use --force to recompute)"
        SKIPPED_ANNOTATED+=("${DS_NAME}")
        continue
      else
        echo "ERROR: clean-entry check failed for ${DS_NAME} (exit ${CHECK_EXIT}); see ${CHECK_LOG_FILE}."
        FAILED_DATASETS+=("${DS_NAME}")
        continue
      fi
    fi
  fi

  # Preprocessing-completeness guard: for datasets with expected views
  # (datasets.json `views.*` carrying input/output file names, mirroring the
  # 1.1.1_preprocess.py skip semantics), ALL expected output .h5ad files must
  # already exist in output/ before chunking. A still-running preprocess array
  # can have written only some of a multi-view dataset's views; chunking on a
  # partial view set would mark the dataset "already annotated" and silently
  # stay incomplete forever.
  EXPECTED_VIEW_FILES=()
  while IFS= read -r view_file; do
    EXPECTED_VIEW_FILES+=("${view_file}")
  done < <(jq -r --arg ds "${DS_NAME}" --arg view "${VIEW_ARG}" \
    '.[$ds].views | to_entries[] | select($view == "" or .key == $view) | select(.value.input_file_name != null) | select(.value.output_file_name != null) | .value.output_file_name' \
    "${DATASETS_JSON_FILE}")

  MISSING_VIEW_FILES=()
  for view_file in "${EXPECTED_VIEW_FILES[@]}"; do
    if [[ ! -f "${HPC_SCRATCH_DIR}/${DS_NAME}/output/${view_file}" ]]; then
      MISSING_VIEW_FILES+=("${view_file}")
    fi
  done

  if [[ ${#MISSING_VIEW_FILES[@]} -gt 0 ]]; then
    echo "WARNING: preprocessing incomplete for ${DS_NAME}: missing view file(s): ${MISSING_VIEW_FILES[*]} — re-run this script after the preprocess array finishes."
    SKIPPED_INCOMPLETE+=("${DS_NAME}")
    continue
  fi

  # Datasets with no expected views (e.g. Zhu, "views": {}) fall back to the
  # plain "any preprocessed .h5ad" check (Zhu path unchanged).
  if [[ ${#EXPECTED_VIEW_FILES[@]} -eq 0 ]] && ! ls "${HPC_SCRATCH_DIR}/${DS_NAME}/output"/*.h5ad >/dev/null 2>&1; then
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
echo "Skipped (already annotated): ${#SKIPPED_ANNOTATED[@]}"
if [[ ${#SKIPPED_ANNOTATED[@]} -gt 0 ]]; then
  echo "  ${SKIPPED_ANNOTATED[*]}"
fi
echo "Flagged (legacy annotation columns, rebuilt for re-annotation): ${#FLAGGED_LEGACY[@]}"
if [[ ${#FLAGGED_LEGACY[@]} -gt 0 ]]; then
  echo "  ${FLAGGED_LEGACY[*]}"
fi
echo "Skipped (no preprocessed .h5ad): ${#SKIPPED_DATASETS[@]}"
if [[ ${#SKIPPED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${SKIPPED_DATASETS[*]}"
fi
echo "Skipped (preprocessing incomplete): ${#SKIPPED_INCOMPLETE[@]}"
if [[ ${#SKIPPED_INCOMPLETE[@]} -gt 0 ]]; then
  echo "  ${SKIPPED_INCOMPLETE[*]}"
fi
echo "Failed: ${#FAILED_DATASETS[@]}"
if [[ ${#FAILED_DATASETS[@]} -gt 0 ]]; then
  echo "  ${FAILED_DATASETS[*]}"
  echo "Resolve the failures, then rerun this script before 2_submit_hpc_array.sh."
  exit 1
fi
