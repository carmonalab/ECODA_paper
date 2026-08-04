# Plan: Per-dataset HPC scratch/NAS output layout

## Goal

Restructure the HPC folder layout so each dataset is self-contained. Preprocessed outputs move from the global `${HPC_SCRATCH_DIR}/output/<DS_NAME>/` to `${HPC_SCRATCH_DIR}/<DS_NAME>/output/`, mirroring the existing per-dataset raw-input convention (`<DS_NAME>/data/`). The NAS mirror becomes per-dataset too (`${NAS_TARGET_DIR}/<DS_NAME>/output/`). `SCRATCH_OUTPUT_DIR` env var is removed; all consumers derive the per-dataset output dir as `${HPC_SCRATCH_DIR}/${DS_NAME}/output`.

## Decisions (confirmed with user)

- Scratch: `<DS_NAME>/data/` (raw inputs, unchanged) + `<DS_NAME>/output/` (preprocessed .h5ad, `chunks/`, `annotations_chunk_*.feather`, annotated .h5ad). `CombinedPBMC/data/` stays (combine output + rds→h5ad cache input); `CombinedPBMC/output/` holds its preprocessed output (CombinedPBMC is a `datasets.json` key, picked up by the preprocess array like any other).
- NAS: mirror per-dataset → `Projects/ECODA_paper/<DS_NAME>/output/`. No global `output/` subdir on NAS anymore.
- `SCRATCH_OUTPUT_DIR` env var: removed entirely (no legacy/reused name). Every consumer either has `DS_NAME` in scope (workers) or iterates over datasets (login-node scripts).
- No migration of existing scratch data (pipelines not run yet).
- `chunks_manifest.txt` moves to `${HPC_SCRATCH_DIR}/chunks_manifest.txt` (no global output dir left).

## Target layout (docs)

```
$HOME/scratch/ECODA_paper               # HPC_SCRATCH_DIR
├── <DS_NAME>/data/                     # staged raw inputs per dataset (1_stage_data.sh)
└── <DS_NAME>/output/                   # preprocessed .h5ad per view, chunks/, annotations, annotated .h5ad
$HOME/reference_atlases/sketched_200ct/ # HOME_REF_DIR (HiTME reference maps)

NAS Projects/ECODA_paper/               # NAS_TARGET_DIR
└── <DS_NAME>/output/                   # rsynced per-dataset from ${HPC_SCRATCH_DIR}/<DS_NAME>/output
```

## Changes

### 1. `src/slurm_config.sh`
- Delete line 25: `export SCRATCH_OUTPUT_DIR="${HPC_SCRATCH_DIR}/output"`. Keep `HPC_SCRATCH_DIR`.

### 2. `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`
- Line 36: `OUTPUT_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}"` → `OUTPUT_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/output"`.

### 3. `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`
- Replace the single rsync (lines 48-55) with a per-dataset loop over existing output dirs:
  ```bash
  echo "Array Job ${ARRAY_JOB_ID} finished. Syncing results to NAS..."
  SYNCED_COUNT=0
  if ls "${NAS_TARGET_DIR}/.." > /dev/null 2>&1; then
    for DS_DIR in "${HPC_SCRATCH_DIR}"/*/output; do
      [[ -d "${DS_DIR}" ]] || continue
      DS_NAME="$(basename "$(dirname "${DS_DIR}")")"
      mkdir -p "${NAS_TARGET_DIR}/${DS_NAME}/output"
      rsync -rlptDv "${DS_DIR}/" "${NAS_TARGET_DIR}/${DS_NAME}/output/"
      SYNCED_COUNT=$((SYNCED_COUNT + 1))
    done
    echo "Results synchronized to ${NAS_TARGET_DIR}/<DS_NAME>/output/ (${SYNCED_COUNT} datasets)"
  else
    echo "ERROR: NAS path ${NAS_TARGET_DIR} is unreachable."
    exit 1
  fi
  ```
  (`set -euo pipefail` is on; the `[[ -d ]]` guard handles the unmatched-glob case.)

### 4. `src/4_cell_type_annotation/1_prepare_chunks.sh`
- Lines 87-88: `"${SCRATCH_OUTPUT_DIR}/${DS_NAME}"/*.h5ad` → `"${HPC_SCRATCH_DIR}/${DS_NAME}/output"/*.h5ad` (both the check and the warning message).

### 5. `src/4_cell_type_annotation/1.1_prepare_chunks.py`
- Docstring line 9: `${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks` → `${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks`.
- Remove the `SCRATCH_OUTPUT_DIR` env read + fallback (lines 43-47); keep the `HPC_SCRATCH_DIR` guard. Set `path_data = Path(hpc_scratch_dir) / ds_name / "output"` (line 49).

### 6. `src/4_cell_type_annotation/2_submit_hpc_array.sh`
- Line 41: `CHUNKS_MANIFEST="${SCRATCH_OUTPUT_DIR}/chunks_manifest.txt"` → `"${HPC_SCRATCH_DIR}/chunks_manifest.txt"`.
- Line 50: `CHUNKS_DIR="${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks"` → `"${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks"`.
- Lines 100-108: replace the single rsync with the same per-dataset loop as in change 3 (identical code).

### 7. `config_helper.R`
- Remove the `SCRATCH_OUTPUT_DIR` env read + fallback (lines 15-18). Keep the `HPC_SCRATCH_DIR` and `DS_NAME` guards.
- `scratch_output_dir <- file.path(hpc_scratch_dir, ds_name, "output")`.
- Update the comment at lines 30-31: paths live under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (still matches `2_submit_hpc_array.sh` chunk dir).

### 8. `docs/ARCHITECTURE.md`
- Lines 36-40: HPC layout block → target layout above.
- Lines 43-49: NAS block → `Projects/ECODA_paper/<DS_NAME>/output/` (drop the `output/` line).
- Lines 52-61 path table: replace the `SCRATCH_OUTPUT_DIR` row with a `${HPC_SCRATCH_DIR}/<DS_NAME>/output/` row (`HPC_SCRATCH_DIR`); update NAS row; update the HPC_SCRATCH_DIR row description.
- Line 91 (NAS ↔ Scratch data flow): `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` → `${HPC_SCRATCH_DIR}/${DS_NAME}/output`.
- Line 126 (`2_submit_hpc_array.sh` role): `${SCRATCH_OUTPUT_DIR}/chunks_manifest.txt` → `${HPC_SCRATCH_DIR}/chunks_manifest.txt`.
- Line 131 (`config_helper.R` role): drop `SCRATCH_OUTPUT_DIR` from the env-var list; `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` → `${HPC_SCRATCH_DIR}/${DS_NAME}/output`.
- Line 137 (data flow): `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` → `${HPC_SCRATCH_DIR}/${DS_NAME}/output`.
- Line 138 (environment propagation): remove `SCRATCH_OUTPUT_DIR` from the exported-vars list.

### 9. `AGENTS.md`
- Line 131: `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` → `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (2 occurrences in the sentence).
- Line 136: HPC folder layout bullet → per-dataset layout (data/ + output/).
- Line 138: NAS bullet → `${NAS_TARGET_DIR}/<DS_NAME>/output/` rsync target.

### 10. `TODO.md`
- Add a completed entry noting the per-dataset output layout change (matches repo convention of logging completed work in TODO.md).

## Explicitly unchanged

- `src/1_stage_data/1_stage_data.sh` (`${HPC_SCRATCH_DIR}/${KEY}/data`), `_create_combinedpbmc_dataset.py` (`${HPC_SCRATCH_DIR}/CombinedPBMC/data`), `_create_joanito_batch_col.R` — already per-dataset input paths.
- `2.1_run_worker.sh` / `2.1.1_process_chunk.sh` / `2.1.1.1_process_chunk.R` / `3_merge_annotations.py` — manifest/chunk-file driven or explicit-arg; code unchanged (paths they consume move with the new layout).
- Notebooks (`benchmark_analysis.rmd`, `batch_effect_analysis.rmd`) — read local `data/` dir; no scratch/NAS references.

## Validation

- No pipeline runs (AGENTS.md rule).
- `grep -rn "SCRATCH_OUTPUT_DIR"` across the repo → zero matches (including docs).
- `bash -n` on the four edited shell scripts.
- `python -m py_compile src/4_cell_type_annotation/1.1_prepare_chunks.py`.
- Parse `config_helper.R` (`pixi run Rscript -e 'parse("config_helper.R")'` or equivalent syntax check).
- Grep ARCHITECTURE.md/AGENTS.md for stale `output/` layout snippets.
