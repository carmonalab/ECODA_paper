#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 3_submit_merge.sh — Merge per-chunk annotation feathers into every view
# h5ad of one dataset, then sync the annotated files back to NAS.
#
# Run AFTER 2_submit_hpc_array.sh <DS_NAME> completed for the same dataset.
#
# Usage:
#   ./3_submit_merge.sh <DS_NAME>
#
# Pipeline stage: annotation runs once per DATASET on the per-dataset union
# h5ad (see 1.1_prepare_chunks.py); the resulting annotations_chunk_*.feather
# files are merged here into EACH view h5ad of the dataset on the (sample,
# barcode) join key (3.1_merge_annotations.py). Afterwards the union file
# (${HPC_SCRATCH_DIR}/${DS}/annotation_union/) and the now-stale chunk files
# (${HPC_SCRATCH_DIR}/${DS}/output/chunks, which pin the deleted union path)
# are removed, and the annotated view h5ads are rsynced to
# ${NAS_TARGET_DIR}/${DS}/output/.
# ---------------------------------------------------------------------------

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils/bash/sync_status_email.sh"
mkdir -p "${LOGS_DIR}"

DS_NAME="${1:-}"
if [[ -z "${DS_NAME}" ]]; then
  echo "ERROR: Usage: ./3_submit_merge.sh <DS_NAME>"
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi
if ! jq -e --arg ds "${DS_NAME}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
  echo "ERROR: '${DS_NAME}' is not a dataset in ${DATASETS_JSON_FILE}."
  exit 1
fi

OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
UNION_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union"
# Memory per merge srun; overridable for the largest views (e.g. Stephenson's
# batch-effect view may OOM on the 64G default).
MERGE_MEM="${MERGE_MEM:-64G}"

# View h5ads to annotate (exclude rds->h5ad conversion caches). Annotation
# feathers live directly in OUTPUT_DIR (written there by 2.1.1_process_chunk.R).
shopt -s nullglob
VIEW_FILES=("${OUTPUT_DIR}"/*.h5ad)
VIEW_FILES_CLEAN=()
for f in "${VIEW_FILES[@]}"; do
  [[ "$(basename "${f}")" != *_raw.h5ad ]] && VIEW_FILES_CLEAN+=("${f}")
done
if [[ ${#VIEW_FILES_CLEAN[@]} -eq 0 ]]; then
  echo "ERROR: No preprocessed view .h5ad files found in ${OUTPUT_DIR}."
  exit 1
fi

ANNOT_FILES=("${OUTPUT_DIR}"/annotations_chunk_*.feather)
if [[ ${#ANNOT_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No annotations_chunk_*.feather files found in ${OUTPUT_DIR}."
  echo "       Run 2_submit_hpc_array.sh ${DS_NAME} first."
  exit 1
fi

# Coverage gate: every chunk submitted by the last 2_submit_hpc_array.sh run
# must have produced its feather (2.1.1_process_chunk.R writes exactly one
# feather per chunk file). A partial annotation array (some chunks failed)
# would otherwise be merged here and rsynced over good NAS data. The manifest
# survives merge re-runs (it is not deleted by this script), so this check
# also passes on legitimate re-runs.
CHUNKS_MANIFEST="${HPC_SCRATCH_DIR}/chunks_manifest.txt"
EXPECTED_CHUNKS=0
if [[ -f "${CHUNKS_MANIFEST}" ]]; then
  EXPECTED_CHUNKS=$(grep -c "^${DS_NAME}[[:space:]]" "${CHUNKS_MANIFEST}" 2>/dev/null || true)
fi
if [[ ${EXPECTED_CHUNKS} -ne ${#ANNOT_FILES[@]} ]]; then
  echo "ERROR: Annotation coverage mismatch for ${DS_NAME}:"
  echo "       ${#ANNOT_FILES[@]} feathers found, but ${EXPECTED_CHUNKS} chunks were submitted"
  echo "       (see ${CHUNKS_MANIFEST}). Re-run 2_submit_hpc_array.sh ${DS_NAME} first."
  notify_sync_status \
    "ECODA: annotation merge NOT synced (${DS_NAME})" \
    "Annotation merge + NAS sync failed for ${DS_NAME}: coverage mismatch (${#ANNOT_FILES[@]} feathers found, ${EXPECTED_CHUNKS} chunks submitted). Re-run 2_submit_hpc_array.sh ${DS_NAME} first."
  exit 1
fi
echo "Found ${#ANNOT_FILES[@]} annotation feather files (matching ${EXPECTED_CHUNKS} chunks) and ${#VIEW_FILES_CLEAN[@]} view h5ads."

# NAS reachability check BEFORE any destructive step or expensive merge: an
# unreachable NAS (e.g. no VPN) must not destroy scratch artifacts.
if ! ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
  echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
  notify_sync_status \
    "ECODA: annotation merge NOT synced (${DS_NAME})" \
    "Annotation merge + NAS sync failed for ${DS_NAME}: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
  exit 1
fi

export DS_NAME

for VIEW_FILE in "${VIEW_FILES_CLEAN[@]}"; do
  VIEW_NAME="$(basename "${VIEW_FILE}")"
  LOG_FILE="${LOGS_DIR}/merge_annotations_${DS_NAME}_${VIEW_NAME%.h5ad}.log"
  echo "Merging annotations into ${VIEW_NAME}..."
  if ! srun --partition="${SLURM_PARTITION}" \
       --time=02:00:00 \
       --ntasks=1 \
       --cpus-per-task=1 \
       --mem="${MERGE_MEM}" \
       --output="${LOG_FILE}" \
       --error="${LOG_FILE}" \
       "${PYTHON_BIN}" "${SCRIPT_DIR}/3.1_merge_annotations.py" \
         --h5ad-path "${VIEW_FILE}" \
         --annot-dir "${OUTPUT_DIR}"; then
    echo "ERROR: Merge failed for ${VIEW_NAME}. See ${LOG_FILE}."
    notify_sync_status \
      "ECODA: annotation merge NOT synced (${DS_NAME})" \
      "Annotation merge + NAS sync failed for ${DS_NAME}: merge failed for ${VIEW_NAME} (see ${LOG_FILE})."
    exit 1
  fi
  echo "✓ Merged ${VIEW_NAME}. Log saved to: ${LOG_FILE}"
done

# Clean up: the union h5ad and the chunk files pinning it are stale now (they
# are regenerated by 1_prepare_chunks.sh on the next annotation run).
# NAS reachability was verified above, before the merges.
rm -rf "${UNION_DIR}"
rm -rf "${OUTPUT_DIR}/chunks"
echo "Removed ${UNION_DIR} and ${OUTPUT_DIR}/chunks (regenerated on next annotation run)."

# Sync annotated view h5ads back to NAS (login node has NAS access).
# --exclude='*.tmp' skips leftover atomic-write temp files (a crash between
# the python write and os.replace would otherwise sync a stale .tmp).
echo "Syncing annotated h5ads to NAS..."
mkdir -p "${NAS_TARGET_DIR}/${DS_NAME}/output"
rsync -rlptDv --exclude='*.tmp' "${OUTPUT_DIR}/" "${NAS_TARGET_DIR}/${DS_NAME}/output/"

# Stale intermediates: 2_submit_hpc_array.sh already synced output/ to NAS
# before the merge, so annotations_chunk_*.feather and chunks/ persist there
# and accumulate across runs. Remove them from NAS (keep the LOCAL feathers —
# the coverage gate above counts them on re-runs).
shopt -s nullglob
NAS_ANNOT_FEATHERS=("${NAS_TARGET_DIR}/${DS_NAME}/output"/annotations_chunk_*.feather)
shopt -u nullglob
if (( ${#NAS_ANNOT_FEATHERS[@]} > 0 )); then
  rm -f "${NAS_ANNOT_FEATHERS[@]}"
fi
rm -rf "${NAS_TARGET_DIR}/${DS_NAME}/output/chunks"
echo "Success: annotated h5ads synchronized to ${NAS_TARGET_DIR}/${DS_NAME}/output/"
notify_sync_status \
  "ECODA: annotations merged + synced (${DS_NAME})" \
  "Annotations for ${DS_NAME} merged into all view h5ads and synced to ${NAS_TARGET_DIR}/${DS_NAME}/output/."
