# Simplify Joanito debug dataset: build it on HPC inside the existing Joanito step

## Goal

Replace the "create debug subset locally + manual rsync to HPC scratch" workflow with a
single HPC dataset-specific step. The debug subset is derived from the already-staged
Joanito raw `.rds` on HPC scratch, inside the existing Joanito step
(`1.2.1_create_joanito_batch_col.R`), which already reads the full Joanito object into
memory. Building the subset from the same in-memory object is nearly free and removes the
local run, the manual rsync, the gitignored `data/debug/` folder, and the NAS `debug/output/`
fallback.

## Decisions (user-approved)

- **Merge debug creation into the existing Joanito step** (no new submit script, no
  dispatcher change). This also eliminates the read/write race that a separate parallel
  step would have with the in-place `saveRDS` of the Joanito `.rds`.
- **Rename the step pair** (user suggested; shorter name chosen):
  - `1.2_submit_joanito_batch_col.sh` → `1.2_submit_joanito.sh`
  - `1.2.1_create_joanito_batch_col.R` → `1.2.1_prepare_joanito.R`
- **`_debug` outputs stay synced to NAS** (`${NAS_TARGET_DIR}/_debug/output/` via the
  existing preprocess-array and merge syncs) — desired for traceability/reproducibility.
  No sync logic changes anywhere.
- **`_debug` raw input never lives on NAS**: `datasets.json` `_debug.folder_name` →
  `null` (CombinedPBMC precedent), so `1_stage_data.sh` cleanly skips it.

## Current workflow (being replaced)

1. Local: `pixi run Rscript src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R`
   → `data/debug/JoaI_2022_35773407_debug_5samples.{rds,h5ad}` (gitignored).
2. Manual: `rsync data/debug/ ${HPC_SCRATCH_DIR}/_debug/data/`
   (or NAS `debug/output/` + `1_stage_data.sh --ds_name _debug`).
3. HPC: preprocess array / annotation chunks+array+merge with `--ds_name _debug`.

## New workflow

1. `src/1_stage_data/1_stage_data.sh` (stages Joanito raw `.rds` to
   `${HPC_SCRATCH_DIR}/Joanito/data/`; `_debug` skipped).
2. `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — the (renamed) Joanito step
   computes `seqtec`, saves it back in place, and additionally writes
   `${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad`.
3. Preprocess: `1_submit_hpc_array.sh --ds_name _debug`; annotation:
   `1_prepare_chunks.sh test _debug`, `2_submit_hpc_array.sh _debug`,
   `3.1_submit_merge.sh _debug` — all unchanged.

## Implementation tasks

1. **Rename + extend `1.2.1_create_joanito_batch_col.R` → `1.2.1_prepare_joanito.R`**
   - Keep: env-driven (`HPC_SCRATCH_DIR` required, `stop()` if unset), read staged Joanito
     raw via `DATASETS_JSON_FILE` env var
     (`config[["Joanito"]][["views"]][["batch_effect_analysis"]][["input_file_name"]]`),
     compute `seqtec`, in-place idempotent `saveRDS`.
   - Add (after the `saveRDS`, reusing the in-memory `seurat`):
     - Seeded (`set.seed(321)`) selection of 5 samples covering
       `(sample.origin × seqtec × Site)` combos, candidates with ≥500 cells (port
       `select_debug_samples()` from the current `1.3.1` script), 500 cells/sample.
     - Keep minimal obs columns (current `cols_keep` list: `sample.ID`, `sample.origin`,
       `patient.ID`, `iCMS`, `dataset`, `Gender`, `Site`, `cell.type`, `seqtec`).
     - Rebuild clean v5 object (`CreateSeuratObject(counts=..., meta.data=...)`).
     - `write_h5ad()` (anndataR; X=None + `layers["counts"]` — handled downstream by the
       X-promotion in `1.1.1_preprocess.py`) to
       `${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad`
       (hardcoded stem must match datasets.json `file_names`/view `input_file_name`;
       `dir.create(recursive=TRUE)`).
     - Keep the summary print (samples, seqtec, Site, n_cells) for verification.
   - Drop: local/NAS input fallbacks, CLI args (`--input`/`--output-dir`/…), the `.rds`
     debug output, and the duplicated `seqtec` mapping comment (now a single source of
     truth — `seqtec` is computed once in this script).
2. **Rename + update `1.2_submit_joanito_batch_col.sh` → `1.2_submit_joanito.sh`**
   - Update `--job-name` (e.g. `joanito_prep`) and the R script path. Keep 32G/1h/cpu
     profile (same full-object read). No dispatcher change needed (glob
     `1.*_submit_*.sh` picks it up).
3. **Delete `src/2_dataset_specific_preprocessing/1.3.1_create_debug_dataset.R`**.
4. **`datasets.json`**: `_debug.folder_name` → `null`. Everything else in the `_debug`
   entry unchanged. (User-approved change; AGENTS.md rule noted.)
5. **`.gitignore`**: remove the `data/debug/` line.
6. **Docs**:
   - `AGENTS.md`: Stage-2 bullet and the `_debug` bullet in "datasets.json" — the Joanito
     step now also builds the debug subset into `${HPC_SCRATCH_DIR}/_debug/data/`; remove
     local-script/manual-rsync/NAS-fallback wording; note `_debug` raw input has
     `folder_name: null` (staging skips it) while `_debug` *outputs* sync to NAS as usual.
   - `TODO.md`: Phase 1 §1.1 (rewrite to the merged-step design), Phase 2 prereqs (drop
     local build + rsync; debug subset comes from the Joanito step), Human-managed tasks,
     Changelog.
   - `docs/ARCHITECTURE.md`: stage-2 table row (renamed step + debug creation) and the
     "Debug subset build" section (line ~100) — now part of `1.2.1_prepare_joanito.R` on
     HPC scratch.

## Notes / constraints

- Do **not** copy the raw debug h5ad into `${HPC_SCRATCH_DIR}/_debug/output/`: the merge
  script treats every `output/*.h5ad` (except `*_raw.h5ad`) as a view to annotate/merge
  into. Raw subset stays in `_debug/data/` only.
- Debug subset creation in every full pipeline run is accepted (cheap: object already in
  memory; keeps NAS artifacts reproducible).
- `1_stage_data.sh --ds_name _debug` becomes a clean no-op (`folder_name: null` → skip
  loop), which is the intended behavior; no code change needed there.
- Debug re-runs downstream still use `--force` on the preprocess array (already
  implemented).
- anndataR + Seurat are available in the pixi env (`PIXI_RSCRIPT`), already used by the
  preprocess rds→h5ad path on workers.

## Validation (no pipeline runs; per AGENTS.md)

- `bash -n` on `1_submit_joanito.sh` and `1_submit_hpc.sh`; R parse-check
  `1.2.1_prepare_joanito.R` (e.g. `pixi run Rscript --vanilla -e "invisible(parse(file='src/2_dataset_specific_preprocessing/1.2.1_prepare_joanito.R'))"`).
- `jq empty datasets.json`; jq smoke: default-all staging emits no `_debug` line,
  `--ds_name _debug` emits a line that the null-folder guard skips.
- Full validation deferred to the user's HPC Phase-2 debug run:
  stage (`1_stage_data.sh`), `1_submit_hpc.sh`, then
  `1_submit_hpc_array.sh --ds_name _debug` → `1_prepare_chunks.sh test _debug` →
  `2_submit_hpc_array.sh _debug` → `3.1_submit_merge.sh _debug`; verify
  `${HPC_SCRATCH_DIR}/_debug/data/*.h5ad` exists after the Joanito step and the debug
  preprocess output lands in `${NAS_TARGET_DIR}/_debug/output/`.
