# HPC Pipeline Review — Fixes & Cleanup (src/2, 3, 4)

Reviewed all scripts in `src/2_dataset_specific_preprocessing/`, `src/3_scrnaseq_preprocessing/`,
`src/4_cell_type_annotation/`, plus `src/slurm_config.sh`, `config_helper.R`, `src/datasets_io.py`,
`src/utils/preprocess_utils.py`, `src/utils/seurat_utils.R`, `docs/ARCHITECTURE.md`.

## Findings (verified against code)

### Critical bugs
1. **`1.1_prepare_chunks.py` globs `*.h5ad` and picks up the `*_raw.h5ad` rds→h5ad conversion caches**
   written into `${HPC_SCRATCH_DIR}/${DS_NAME}/output` by `preprocess_utils.load_single_input()`
   (every `.rds`-input dataset: Adams, Bassez, Joanito, Kfoury, Kim, Lee, Pelka, Smillie, Stephenson, Wu, Zhang).
   The cache file lacks the standardized `Sample` column → `KeyError` crash for most datasets; for
   datasets whose original sample col *is* "Sample" (e.g. Stephenson) it silently creates chunks for the
   **raw, un-subsetted** file → wrong annotation input.
2. **`2.1.1_process_chunk.R` feeds log-normalized `adata.X` to scATOMIC/HiTME as counts**
   (`get_sample_seurat_obj`, line 36: `subset_py$X$astype("float64")$tocsc()`). Preprocessing vaults raw
   counts into `layers["counts"]`; X is `normalize_total + log1p` data. Annotations would be computed on
   wrong input. The repo already has a correct helper: `get_seurat_obj_from_h5ad()` in
   `src/utils/seurat_utils.R:266` (takes `counts_layer="counts"`).
3. **`3_merge_annotations.py` joins on `cell_barcode` only** — barcodes repeat across samples →
   pandas duplicate-index join explodes `obs` rows. Must join on (sample, barcode). Also this script is
   an orphan: no bash wrapper calls it (ARCHITECTURE.md presents it as a pipeline stage, but it is manual).
4. **Stale feather hazard**: `annotations_chunk_<SLURM_ARRAY_TASK_ID>.feather` (task IDs global across
   datasets). Single-dataset reruns or chunk-size changes renumber task IDs → stale feathers from earlier
   runs get merged → duplicate/corrupt annotations.

### Cleanup / inconsistencies
5. **Dead env vars / params**: `HOME_CHUNKS_DIR` exported in `2.1_run_worker.sh:43` but never read;
   `config_helper.R` `test_mode`/`chunk_size`/`max_test_array_jobs` unused by `2.1.1_process_chunk.R`.
6. **Python invocation inconsistency**: bare `python` in `1.1_submit_combinedpbmc.sh` and preprocess
   `1.1_run_worker.sh`; `${PROJECT_ROOT}/.pixi/envs/default/bin/python` in `1_prepare_chunks.sh`;
   `${HOME}/.pixi/bin/pixi run Rscript` elsewhere. Bare `python` on HPC worker nodes likely picks a
   python without scanpy/anndata → failure.
7. **`SLURM_PARTITION` / `SLURM_ACCOUNT` in slurm_config.sh unused** — partition hardcoded per script.
8. **No `#SBATCH --output/--error`** in `1.1_submit_combinedpbmc.sh` / `1.2_submit_joanito_batch_col.sh`
   → `slurm-<id>.out` files dumped into PROJECT_ROOT on HPC.
9. **`mkdir -p "${LOGS_DIR}"` missing** in `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` and
   `4_cell_type_annotation/2_submit_hpc_array.sh` (only `1_prepare_chunks.sh:17` creates it; `logs/` is
   untracked → likely absent on HPC → `sbatch --output` fails).
10. **No failure detection before NAS sync**: both array submit scripts sync and exit 0 even if tasks
    FAILED/TIMEOUT/OOM (dispatcher `1_submit_hpc.sh` does check sacct states — mirror that).
11. **`2.1.1_process_chunk.R` `get_scGateDB(branch="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4", force_update=TRUE)` per chunk** → up to 500
    concurrent downloads; use a cached copy.
12. **Sample-name standardization runs BEFORE `apply_subset_vars`** in `1.1.1_preprocess.py:213-238` —
    latent mismatch if a subset value contains `-` or a leading digit (current datasets unaffected).
13. `2_submit_hpc_array.sh` syncs **all** datasets even in single-dataset mode.
14. Local artifact `src/2_dataset_specific_preprocessing/__pycache__/*.pyc` (gitignored; delete).
15. TODO.md:31-33 references `TODO_STUMP_1_submit_hpc_array.sh`/`TODO_STUMP_1.1_run_worker.sh` in
    `src/5_run_benchmark_methods/run_python_sample_embedding_methods/` — these files don't exist.
16. **Verify items** (no code change without user confirmation): `1.2_submit_joanito_batch_col.sh` 8G
    may be tight for whole-Joanito `readRDS`; annotation worker `--time=02:00:00` tight given 5×2
    retry timeouts; preprocess worker 16G baseline for GongSharma.
17. **Verified OK**: dispatcher glob `1.*_submit_*.sh` does NOT match `1_submit_hpc.sh` (no recursion);
    docs (ARCHITECTURE.md) match code; array task↔dataset index mapping consistent (`jq keys[]` sorted +
    `sed -n`); `create_clean_seuratv5_object` keeps all meta.data columns (Zhang `Treatment` subset is fine).

### Debug-mode status
- **Debug dataset NOT implemented**: no `_create_debug_dataset.R`, no `data/debug/`, no `_debug` entry
  in datasets.json; TODO.md Step 2 marked "Stale".
- Partial test support that DOES exist: `1_prepare_chunks.sh test` + `1.1_prepare_chunks.py --test`
  (chunk size 1); preprocess `--ds_name` filter; single-dataset args for chunk prep + annotation submit.
- Per user instruction, the full debug/test execution mode goes to TODO.md (Section "TODO.md additions").

---

## Phase 1 — Direct fixes (ordered; each step independently safe)

### 1. `src/4_cell_type_annotation/1.1_prepare_chunks.py`
- In the h5ad loop, skip cache files: `if f.name.endswith("_raw.h5ad"): continue` (comment: rds→h5ad
  conversion cache from `preprocess_utils.load_single_input()`).
- Additionally delete stale `annotations_chunk_*.feather` from `path_data` at start (next to the chunks-dir
  wipe) to fix stale-feather merges.

### 2. `src/4_cell_type_annotation/2.1.1_process_chunk.R`
- Replace `get_sample_seurat_obj()` with the existing `get_seurat_obj_from_h5ad()` helper
  (`src/utils/seurat_utils.R:266`): `counts_layer="counts"`, with fallback to `X` + warning if the layer
  is missing. Requires sourcing `src/utils/seurat_utils.R` (load_all_functions.R already sourced via
  preprocess_utils.py on the Python side — in R, add `source(file.path(project_root, "src/utils/load_all_functions.R"))`
  or source seurat_utils.R directly; verify the R script's existing imports).
- Stable feather naming: derive from the chunk file instead of `SLURM_ARRAY_TASK_ID`:
  `annot_file <- file.path(paths$path_output, paste0("annotations_", sub("\\.txt$", ".feather", basename(args$chunk_file))))`
  (chunk_1.txt → annotations_chunk_1.feather; still matches the merge glob).
- Cache scGateDB: run `get_scGateDB(branch="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4", force_update=TRUE)` once and save into a file (e.g. `aux/scGateDB.rds` and document how it was created) and pass via env.

### 3. `src/4_cell_type_annotation/3_merge_annotations.py`
- Join key = (sample, barcode): in the annotations df build `_key = sample + "_" + cell_barcode`;
  in `adata.obs` build `_key = obs[SAMPLE_COLNAME].astype(str) + "_" + obs_names` (read SAMPLE_COLNAME
  from env, default "Sample"); join on `_key`, then drop the `_key` columns. Keep the column whitelist logic.

### 4. `src/4_cell_type_annotation/2.1_run_worker.sh`
- Remove `export HOME_CHUNKS_DIR=...` (dead). Keep TISSUE_TYPE/NORMAL_TISSUE exports.

### 5. `config_helper.R`
- Remove unused `test_mode` / `chunk_size` / `max_test_array_jobs` entries (or mark clearly if the future
  debug mode will re-use them — see TODO.md item below).

### 6. `src/slurm_config.sh` + python invocation standardization
- Add: `export PYTHON_BIN="${PROJECT_ROOT}/.pixi/envs/default/bin/python"` (mirror `1_prepare_chunks.sh`).
- Use it in: `1.1_submit_combinedpbmc.sh`, `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`,
  `src/4_cell_type_annotation/1_prepare_chunks.sh`.
- Optionally also add `PIXI_RSCRIPT="${HOME}/.pixi/bin/pixi run Rscript --vanilla"` and use in
  `1.2_submit_joanito_batch_col.sh` + `2.1_run_worker.sh` (consistency only).

### 7. `src/2_dataset_specific_preprocessing/` step scripts
- `1.1_submit_combinedpbmc.sh` and `1.2_submit_joanito_batch_col.sh`: add
  `#SBATCH --output=${LOGS_DIR}/...` / `#SBATCH --error=${LOGS_DIR}/...` (sbatch does not expand vars in
  directives — use `--output="${LOGS_DIR}/2_dataset_specific_preprocessing_%j.log"` passed on the sbatch
  command line from `1_submit_hpc.sh`, or plain paths). Simplest: pass `--output`/`--error` flags in
  `1_submit_hpc.sh` when submitting each step.
- `1.2_submit_joanito_batch_col.sh`: bump `--mem` (e.g. 32G) — whole-Joanito `readRDS`+`saveRDS`; verify on cluster.

### 8. Array submit scripts — logging dir + failure detection
- `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` and
  `src/4_cell_type_annotation/2_submit_hpc_array.sh`: add `mkdir -p "${LOGS_DIR}"` before `sbatch`.
- After the monitoring loop, before rsync: check states via
  `sacct -j "${ARRAY_JOB_ID}" --format=State -n | grep -qE 'FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL'`
  → if hit, exit 1 without syncing (mirror the dispatcher's philosophy).
- `2_submit_hpc_array.sh`: when `DS_NAME_ARG` is given, sync only that dataset's output dir.

### 9. `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`
- Move `apply_subset_vars()` before the sample-col standardization block (subset on original values), and
  add a guard: if the subset is empty, raise an error with the view name.

### 10. Housekeeping
- Delete `src/2_dataset_specific_preprocessing/__pycache__/`.
- Fix stale TODO.md:31-33 `TODO_STUMP_*` references (files don't exist under
  `run_python_sample_embedding_methods/` — only `1.2_benchmark_methods_py.qmd`).
  - Add `TODO.md` entries for the new files.

### 11. Docs (after code changes)
- Update `docs/ARCHITECTURE.md` + `AGENTS.md` + `TODO.md` only where behavior changed (feather naming,
  counts-layer use, merge join key, PYTHON_BIN var, per-step log files). Keep existing "not yet polished"
  notes for `2.1.1_process_chunk.R`.

## Validation (per AGENTS.md: do NOT run pipeline scripts; syntax-level only)
- `bash -n` on every touched `.sh`; `python -m py_compile` on touched `.py`; `Rscript -e 'parse(file=...)'`
  on touched `.R`.
- Grep: no remaining `HOME_CHUNKS_DIR`, bare `python` invocations in the three
  pipeline dirs.
- On HPC (when user is on UNIGE network): `sbatch --test-only` for the step scripts + array workers;
  optionally a real run of the combinedpbmc/joanito steps.

---

## Phase 2 — TODO.md additions (large items, out of scope for direct fixes)

### Item A — Debug/test execution mode with subsetted Joanito debug dataset (revive TODO.md Step 2)
Current status to record: Step 2 is stale; debug dataset NOT implemented (no `_create_debug_dataset.R`,
no `data/debug/`, no `_debug` entry in datasets.json). Already-supported pieces: `1_prepare_chunks.sh test`
→ `1.1_prepare_chunks.py --test` (1 sample/chunk); preprocess `--ds_name`; single-dataset args in chunk
prep + annotation submit. Plan:
- Implement `src/2_dataset_specific_preprocessing/_create_debug_dataset.R` per old Step 2b (5 Joanito
  samples, 500 cells/sample, minimal obs cols) → `data/debug/*.rds` + `*.h5ad`.
- Register as `_debug` in datasets.json (or separate `debug_datasets.json`; decide) with
  `batch_effect_analysis` + `benchmark_analysis` views.
- Wire test mode end-to-end: preprocess array `--ds_name _debug` (+ view filter), chunk prep `test`,
  annotation array on the debug dataset, validation checklist (h5ad loads, keys `X_pca`/`X_pca_harmony`,
  runtime < 30s locally).
- Coverage matrix: which scripts support test/single-dataset mode today vs. which need adaptions
  (`1_submit_hpc_array.sh`, `1.1_run_worker.sh`, `2_submit_hpc_array.sh`, `2.1_run_worker.sh` — the last
  two already support a single DS_NAME arg).

### Item B — Centralize SLURM partition (extends existing TODO Step 6)
- Decide: pass `--partition="${SLURM_PARTITION}"` at submit time (sbatch directives don't expand vars);
  remove `#SBATCH --partition` lines from `1.1_submit_combinedpbmc.sh`, `1.2_submit_joanito_batch_col.sh`,
  preprocess + annotation workers, and the `srun --partition=shared-cpu` in `1_prepare_chunks.sh`.
- Note per-pipeline partition needs (shared-cpu for preprocess/annotation; GPU for some benchmark methods).

### Item C — Integrate `3_merge_annotations.py` into the pipeline
- Decide: add a `3.1_submit_merge.sh`/`3_submit_merge.sh` wrapper (per-dataset, after array completion) or
  document it as an explicit manual step in ARCHITECTURE.md (remove the "stage" implication).
  - Could this be implemented directly in a previous step? (per-dataset, after array completion)
  - Could this be implemented in `2_submit_hpc_array.sh` (after all workers complete) or is it cleaner and safer to keep this separate?
