#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 3_submit_merge.sh — Merge per-chunk annotation feathers into every view
# h5ad of a dataset, then sync the annotated files back to NAS.
#
# Run AFTER 2_submit_hpc_array.sh completed for the same dataset(s).
#
# Usage:
#   ./3_submit_merge.sh                    # default-all: every dataset in
#                                          #   datasets.json (INCLUDING _debug).
#                                          #   Already-merged datasets skip,
#                                          #   no-feather datasets WARNING-skip,
#                                          #   failures are collected; one
#                                          #   summary email; exit 1 if any
#                                          #   dataset failed.
#   ./3_submit_merge.sh --force            # default-all, force re-merge. An
#                                          #   already-merged dataset has no
#                                          #   chunks, so it fails the coverage
#                                          #   gate per dataset — expected.
#   ./3_submit_merge.sh <DS_NAME>          # merge + NAS sync (skips if already merged)
#   ./3_submit_merge.sh <DS_NAME> --force  # force re-merge (coverage gate still applies)
#
# Pipeline stage: annotation runs once per DATASET on the per-dataset union
# h5ad (see 1.1_prepare_chunks.py); the resulting annotations_chunk_*.feather
# files are merged here into EACH view h5ad of the dataset on the (sample,
# barcode) join key (3.1_merge_annotations.py). Afterwards the union file
# (${HPC_SCRATCH_DIR}/${DS}/annotation_union/) and the now-stale chunk files
# (${HPC_SCRATCH_DIR}/${DS}/output/chunks, which pin the deleted union path)
# are removed, and the annotated view h5ads are rsynced to
# ${NAS_TARGET_DIR}/${DS}/output/.
#
# Single-dataset mode sends per-event sync-status emails (each success email
# carries per-view merge durations); default-all mode sends exactly one summary
# email (processed/failed lists + per-dataset merge durations) after the loop.
# ---------------------------------------------------------------------------

source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
cd "${PROJECT_ROOT}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils/bash/sync_status_email.sh"
mkdir -p "${LOGS_DIR}"

# Merge durations for the default-all summary email: one "<DS_NAME>|<seconds>"
# entry per processed dataset, appended by the caller loop below. Per-view
# durations for single-dataset mode live in merge_one_ds's VIEW_DURATIONS.
DS_DURATIONS=()

fmt_seconds() {
  local s="$1"
  printf '%02d:%02d:%02d' "$(( s / 3600 ))" "$(( (s % 3600) / 60 ))" "$(( s % 60 ))"
}

FORCE_ARG=0
POS_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)
      FORCE_ARG=1
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

DS_NAME="${POS_ARGS[0]:-}"
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi
if [[ -n "${DS_NAME}" ]]; then
  if ! jq -e --arg ds "${DS_NAME}" 'has($ds)' "${DATASETS_JSON_FILE}" > /dev/null 2>&1; then
    echo "ERROR: '${DS_NAME}' is not a dataset in ${DATASETS_JSON_FILE}."
    exit 1
  fi
fi

# Per-dataset merge: skip check -> view/feather globs -> coverage gate -> NAS
# check -> per-view srun merge -> cleanup -> rsync. Single-dataset mode calls
# this once and uses its return code as the exit code (per-event emails);
# default-all mode calls it once per dataset and collects failures (one summary
# email from the caller). Hard `exit`s would kill a default-all loop, so every
# failure path `return 1` instead; per-event notify_sync_status calls are
# guarded by EMAIL_MODE so single-dataset behavior is preserved exactly.
merge_one_ds() {
  local DS_NAME="$1"
  local DS_START="$(date +%s)"
  local OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"
  local UNION_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union"
  # Memory per merge srun; overridable for the largest views (e.g. Stephenson's
  # batch-effect view may OOM on the 64G default).
  local MERGE_MEM="${MERGE_MEM:-64G}"
  # Per-view srun durations for the per-event email: "<VIEW_NAME>|<seconds>".
  local VIEW_DURATIONS=()

  # Post-merge skip: this script is the only step that deletes output/chunks/ and
  # annotation_union/, so a dataset with both absent and >=1 feather still in
  # output/ was already merged (feathers are kept deliberately for the coverage
  # gate on re-runs). Nothing to merge or sync — return 0 (no srun, no NAS sync,
  # no sync-status email) unless --force.
  shopt -s nullglob
  local ANNOT_FILES=("${OUTPUT_DIR}"/annotations_chunk_*.feather)
  shopt -u nullglob
  if [[ ${FORCE_ARG} -eq 0 && ! -d "${OUTPUT_DIR}/chunks" && ! -d "${UNION_DIR}" && ${#ANNOT_FILES[@]} -gt 0 ]]; then
    echo "Already merged: ${DS_NAME} — skipping merge + NAS sync (use --force to re-merge)"
    return 0
  fi

  # View h5ads to annotate (exclude rds->h5ad conversion caches). Annotation
  # feathers live directly in OUTPUT_DIR (written there by 2.1.1_process_chunk.R).
  shopt -s nullglob
  local VIEW_FILES=("${OUTPUT_DIR}"/*.h5ad)
  shopt -u nullglob
  local VIEW_FILES_CLEAN=()
  if (( ${#VIEW_FILES[@]} > 0 )); then
    for f in "${VIEW_FILES[@]}"; do
      [[ "$(basename "${f}")" != *_raw.h5ad ]] && VIEW_FILES_CLEAN+=("${f}")
    done
  fi
  if [[ ${#VIEW_FILES_CLEAN[@]} -eq 0 ]]; then
    if [[ "${EMAIL_MODE}" == "per-event" ]]; then
      echo "ERROR: No preprocessed view .h5ad files found in ${OUTPUT_DIR}."
      return 1
    fi
    echo "WARNING: No preprocessed view .h5ad files found in ${OUTPUT_DIR}; skipping ${DS_NAME}."
    return 0
  fi

  shopt -s nullglob
  ANNOT_FILES=("${OUTPUT_DIR}"/annotations_chunk_*.feather)
  shopt -u nullglob
  if [[ ${#ANNOT_FILES[@]} -eq 0 ]]; then
    if [[ "${EMAIL_MODE}" == "per-event" ]]; then
      echo "ERROR: No annotations_chunk_*.feather files found in ${OUTPUT_DIR}."
      echo "       Run 2_submit_hpc_array.sh ${DS_NAME} first."
      return 1
    fi
    echo "WARNING: No annotations_chunk_*.feather files found in ${OUTPUT_DIR}; skipping ${DS_NAME}."
    return 0
  fi

  # Coverage gate: every chunk submitted by the last 2_submit_hpc_array.sh run
  # must have produced its feather (2.1.1_process_chunk.R writes exactly one
  # feather per chunk file). A partial annotation array (some chunks failed)
  # would otherwise be merged here and rsynced over good NAS data. The manifest
  # survives merge re-runs (it is not deleted by this script), so this check
  # also passes on legitimate re-runs.
  local CHUNKS_MANIFEST="${HPC_SCRATCH_DIR}/chunks_manifest.txt"
  local EXPECTED_CHUNKS=0
  if [[ -f "${CHUNKS_MANIFEST}" ]]; then
    EXPECTED_CHUNKS=$(grep -c "^${DS_NAME}[[:space:]]" "${CHUNKS_MANIFEST}" 2>/dev/null || true)
  fi
  if [[ ${EXPECTED_CHUNKS} -ne ${#ANNOT_FILES[@]} ]]; then
    echo "ERROR: Annotation coverage mismatch for ${DS_NAME}:"
    echo "       ${#ANNOT_FILES[@]} feathers found, but ${EXPECTED_CHUNKS} chunks were submitted"
    echo "       (see ${CHUNKS_MANIFEST}). Re-run 2_submit_hpc_array.sh ${DS_NAME} first."
    local GATE_HINT=""
    if [[ ! -d "${OUTPUT_DIR}/chunks" ]]; then
      GATE_HINT="Dataset has no chunks (already merged or never prepared); run 1_prepare_chunks.sh ${DS_NAME} --force, then 2_submit_hpc_array.sh ${DS_NAME} first."
      echo "       ${GATE_HINT}"
    fi
    if [[ "${EMAIL_MODE}" == "per-event" ]]; then
      notify_sync_status \
        "ECODA: annotation merge NOT synced (${DS_NAME})" \
        "Annotation merge + NAS sync failed for ${DS_NAME}: coverage mismatch (${#ANNOT_FILES[@]} feathers found, ${EXPECTED_CHUNKS} chunks submitted). Re-run 2_submit_hpc_array.sh ${DS_NAME} first.${GATE_HINT:+ ${GATE_HINT}}"
    fi
    return 1
  fi
  echo "Found ${#ANNOT_FILES[@]} annotation feather files (matching ${EXPECTED_CHUNKS} chunks) and ${#VIEW_FILES_CLEAN[@]} view h5ads."

  # NAS reachability check BEFORE any destructive step or expensive merge: an
  # unreachable NAS (e.g. no VPN) must not destroy scratch artifacts.
  if ! ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
    if [[ "${EMAIL_MODE}" == "per-event" ]]; then
      notify_sync_status \
        "ECODA: annotation merge NOT synced (${DS_NAME})" \
        "Annotation merge + NAS sync failed for ${DS_NAME}: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
    fi
    return 1
  fi

  export DS_NAME

  for VIEW_FILE in "${VIEW_FILES_CLEAN[@]}"; do
    local VIEW_NAME="$(basename "${VIEW_FILE}")"
    local LOG_FILE="${LOGS_DIR}/merge_annotations_${DS_NAME}_${VIEW_NAME%.h5ad}.log"
    local VIEW_START="$(date +%s)"
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
      if [[ "${EMAIL_MODE}" == "per-event" ]]; then
        notify_sync_status \
          "ECODA: annotation merge NOT synced (${DS_NAME})" \
          "Annotation merge + NAS sync failed for ${DS_NAME}: merge failed for ${VIEW_NAME} (see ${LOG_FILE})."
      fi
      return 1
    fi
    VIEW_DURATIONS+=("${VIEW_NAME}|$(( $(date +%s) - VIEW_START ))")
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
  local NAS_ANNOT_FEATHERS=("${NAS_TARGET_DIR}/${DS_NAME}/output"/annotations_chunk_*.feather)
  shopt -u nullglob
  if (( ${#NAS_ANNOT_FEATHERS[@]} > 0 )); then
    rm -f "${NAS_ANNOT_FEATHERS[@]}"
  fi
  rm -rf "${NAS_TARGET_DIR}/${DS_NAME}/output/chunks"
  echo "Success: annotated h5ads synchronized to ${NAS_TARGET_DIR}/${DS_NAME}/output/"
  if [[ "${EMAIL_MODE}" == "per-event" ]]; then
    # Per-view merge durations (wall clock around each srun; login-node driven).
    local DUR_BLOCK="Merge durations (view, elapsed):
"
    local line view_name view_secs
    for line in "${VIEW_DURATIONS[@]}"; do
      view_name="${line%%|*}"
      view_secs="${line##*|}"
      DUR_BLOCK+="  ${view_name}  $(fmt_seconds "${view_secs}")"$'\n'
    done
    DUR_BLOCK+="Total merge duration: $(fmt_seconds "$(( $(date +%s) - DS_START ))")"
    notify_sync_status \
      "ECODA: annotations merged + synced (${DS_NAME})" \
      "Annotations for ${DS_NAME} merged into all view h5ads and synced to ${NAS_TARGET_DIR}/${DS_NAME}/output/.
${DUR_BLOCK}"
  fi
  return 0
}

if [[ -n "${DS_NAME}" ]]; then
  # Single-dataset mode: per-event sync-status emails; exit code = merge result.
  EMAIL_MODE="per-event"
  RC=0
  merge_one_ds "${DS_NAME}" || RC=$?
  exit ${RC}
fi

# Default-all mode: bare invocation (optionally with --force) = every dataset
# key in datasets.json, INCLUDING _debug (mirrors the annotation default-all
# convention in 1_prepare_chunks.sh / 2_submit_hpc_array.sh — a fresh _debug is
# WARNING-skipped like any no-feather dataset). Already-merged datasets skip,
# no-feather datasets WARNING-skip, merge failures are collected; one summary
# email after the loop; exit 1 if any dataset failed.
EMAIL_MODE="summary"

# NAS is needed for every dataset — check once up front (fail fast with one
# email + exit 1) instead of per dataset inside the loop.
if ! ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
  echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
  notify_sync_status \
    "ECODA: annotation merge NOT synced (all datasets)" \
    "Default-all annotation merge aborted: NAS path ${NAS_TARGET_DIR} is unreachable (check VPN/NAS mount)."
  exit 1
fi

FAILED_DATASETS=()
OK_DATASETS=()
while IFS= read -r DS; do
  echo ""
  echo ">>> Merging annotations for dataset: ${DS} <<<"
  DS_LOOP_START="$(date +%s)"
  if merge_one_ds "${DS}"; then
    OK_DATASETS+=("${DS}")
  else
    FAILED_DATASETS+=("${DS}")
  fi
  DS_DURATIONS+=("${DS}|$(( $(date +%s) - DS_LOOP_START ))")
done < <(jq -r 'keys[]' "${DATASETS_JSON_FILE}")

DURATIONS_BLOCK="Merge durations (dataset, elapsed):
"
for line in "${DS_DURATIONS[@]}"; do
  DURATIONS_BLOCK+="  ${line%%|*}  $(fmt_seconds "${line##*|}")"$'\n'
done

if (( ${#FAILED_DATASETS[@]} > 0 )); then
  notify_sync_status \
    "ECODA: annotation merge NOT synced (${#FAILED_DATASETS[@]} dataset(s) failed)" \
    "Default-all annotation merge finished with failures. Failed: ${FAILED_DATASETS[*]}. Processed OK: ${OK_DATASETS[*]}. See the run stdout for details.
${DURATIONS_BLOCK}"
  exit 1
fi
notify_sync_status \
  "ECODA: annotations merged + synced (all datasets)" \
  "Default-all annotation merge finished: ${#OK_DATASETS[@]} datasets processed (merged or skipped) and synced to ${NAS_TARGET_DIR}/<DS_NAME>/output/ (${OK_DATASETS[*]}).
${DURATIONS_BLOCK}"
exit 0
