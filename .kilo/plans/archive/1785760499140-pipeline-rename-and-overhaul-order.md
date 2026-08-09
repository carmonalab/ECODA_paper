# Pipeline Folder Rename + Overhaul — Implementation Plan

## Goal

Implement TODO.md items in the best logical order: first the folder rename (TODO.md:3-9, previously missing from the plan), then the open reviewer points (TODO.md:11-27) in their existing priority order, all executed in the final (renamed) locations.

## Decisions (confirmed with user)

1. **Rename first.** Pure mechanical `git mv` (history preserved, no logic change); prerequisite for the pipeline sequence (staging must be split out so `_create_combinedpbmc_dataset.py` runs between staging and preprocessing). Subsequent work lands in final paths; docs updated once.
2. **`_preprocess_utils.py` → `src/utils/preprocess_utils.py`.** Hard constraint: numbered dirs (e.g. `src/3_scrnaseq_preprocessing/`) cannot be Python packages — `from src.3_scrnaseq_preprocessing._preprocess_utils import ...` is a SyntaxError. Both importers update to `from src.utils.preprocess_utils import ...`.
3. **SAMPLE_COLNAME: keep the env var**, and additionally read it in `1.1.1_preprocess.py` (currently hardcodes `"Sample"`) so the constant is centralized in `slurm_config.sh` for both R and Python consumers.

## Overall order

- **Phase A** — Folder rename + staging extraction + reference/docs update (do FIRST)
- **Phase B** — Reviewer points in TODO.md's priority order: B1 `NAS_TARGET_DIR`, B2 `SAMPLE_COLNAME`, B3 rewrite `1.1.1_preprocess.py`, B4 rewrite `_create_combinedpbmc_dataset.py`, B5 remove `getwd()` fallbacks, B6 `config_helper.R`/Step 1 simplification, B7 centralisation sweep + final docs pass.

---

## Phase A — Folder rename (git mv only, no logic changes)

### A1. Create new directories

```
src/1_stage_data/
src/2_dataset_specific_preprocessing/
src/3_scrnaseq_preprocessing/
src/4_cell_type_annotation/
src/5_run_benchmark_methods/
```

### A2. Move map (`git mv`, preserves history; working-tree modifications on `1.1.1_preprocess.py`, `_create_combinedpbmc_dataset.py`, `TODO.md` are carried along untouched)

| Current | Target |
|---|---|
| `src/preprocess/1_submit_hpc_array.sh` | `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (Phase 1 removed — see A3) |
| `src/preprocess/1.1_run_worker.sh` | `src/3_scrnaseq_preprocessing/1.1_run_worker.sh` |
| `src/preprocess/1.1.1_preprocess.py` | `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` |
| `src/preprocess/_preprocess_utils.py` | `src/utils/preprocess_utils.py` |
| `src/preprocess/_create_combinedpbmc_dataset.py` | `src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py` |
| `src/preprocess/_create_joanito_batch_col.R` | `src/2_dataset_specific_preprocessing/_create_joanito_batch_col.R` |
| `src/preprocess/TODO_STUMP_preprocess_sikkema.qmd`, `preprocess_gongsharma.qmd` | `src/3_scrnaseq_preprocessing/` |
| `src/cell_type_annotation/*` (incl. `DEPRECATED_LEGACY_CODE/`, `3_merge_annotations.py`) | `src/4_cell_type_annotation/` |
| `src/benchmark/benchmark_methods_r.R`, `benchmark_pipeline.R`, `run_python_sample_embedding_methods/` | `src/5_run_benchmark_methods/` |

Delete now-empty `src/preprocess/`, `src/cell_type_annotation/`, `src/benchmark/`. Ignore `__pycache__`/`.DS_Store`.

### A3. Extract staging into `src/1_stage_data/1_stage_data.sh`

- New standalone login-node script: `source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"`, `cd "${PROJECT_ROOT}"`, `module load GCCcore/12.2.0` + `jq/1.6`, then the Phase 1 loop from `1_submit_hpc_array.sh:10-42` verbatim (jq over `datasets.json` skipping `folder_name: null`, rsync `NAS_SC_DIR/...` → `${HPC_SCRATCH_DIR}/${KEY}/data`).
- `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` keeps Phase 2 (array submit) + Phase 3 (monitor + rsync to NAS). Optionally rename log files (`logs/preprocess_%A_%a.log` → `logs/3_scrnaseq_preprocessing_%A_%a.log`).
- Pipeline sequence becomes explicit: `1_stage_data` → `2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py` → `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` → `4_cell_type_annotation/` → `5_run_benchmark_methods/`.

Note: `../slurm_config.sh` depth is unchanged for all new dirs (all are one level under `src/`, same as before; `run_python_sample_embedding_methods/` stays depth-3 with `../../slurm_config.sh` — verify).

### A4. Reference updates

- `src/utils/load_all_functions.R:17-18`: `src/benchmark/...` → `src/5_run_benchmark_methods/...`
- `src/slurm_config.sh:34` comment: `src/preprocess/1.1.1_preprocess.py` → `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:12` and `src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py:45`: `from src.preprocess._preprocess_utils import ...` → `from src.utils.preprocess_utils import ...`
- `_create_combinedpbmc_dataset.py` header docstring (lines 19-27): update `src/preprocess/1_submit_hpc_array.sh` mention and R-interop note (relative to new dir, `parents[2]` still resolves to project root — verify `sys.path` insert still points at repo root)
- Docs: `AGENTS.md`, `README.md`, `docs/ARCHITECTURE.md`, `TODO.md` (mark rename section done; update all path mentions incl. Completed entries)
- `.kilo/plans/*.md`: leave untouched (historical)

### A5. Rename-phase validation (no pipeline runs, per AGENTS.md)

- `bash -n` on all modified `.sh`; `python3 -m py_compile` on modified `.py`
- Grep sweeps return nothing outside docs/plans: `src/preprocess`, `src/cell_type_annotation`, `src/benchmark`, `from src.preprocess`
- `git status` shows renames (R) + intended reference edits; `git log` history preserved

---

## Phase B — Reviewer points (in TODO.md priority order, final locations)

### B1. `NAS_TARGET_DIR` (config-only, unblocks result syncing)

- `src/slurm_config.sh:20`: `NAS_TARGET_DIR="${NAS_PREFIX}/Projects/ECODA_paper"` (= `/Volumes/Shared/Projects/ECODA_paper` locally)
- rsync steps in `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` + `4_cell_type_annotation/2_submit_hpc_array.sh` need no change (they use the var)

### B2. `SAMPLE_COLNAME` (keep env var, extend to Python)

- Keep `export SAMPLE_COLNAME="Sample"` in `slurm_config.sh:41`; keep R guards in `1.1_prepare_chunks.r:66-67` / `2.1.1.1_process_chunk.R:51`
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:213-216`: replace hardcoded `"Sample"` with `sample_col_out = os.environ.get("SAMPLE_COLNAME", "Sample")` and write `adata_full.obs[sample_col_out]`

### B3. Rewrite `1.1.1_preprocess.py` CLI (unblocks HPC preprocessing)

- Verify argparse (`--config_path`/`--base_path`/`--output_dir`) is honored by `main()` (partially done per TODO Completed); ensure relative args are joined against `PROJECT_ROOT` only, absolute args kept as-is; missing env → clear `stop()`-style errors

### B4. `_create_combinedpbmc_dataset.py` (HPC-capable; run before preprocess array)

- Verify per-dataset layout default when `HPC_SCRATCH_DIR` set, local flat fallback (mostly done per TODO Completed)
- Optional: `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — single `sbatch` job running the combine script (for GongSharma OOM on login node)
- Human step (documented, not executed): verify Zhu raw file has a `Sample` obs column + raw counts; run on HPC from `${PROJECT_ROOT}` with `module load GCCcore/12.2.0`, after `1_stage_data` and before the preprocess array

### B5. Remove `getwd()` fallbacks (env vars now guaranteed)

- `4_cell_type_annotation/1.1_prepare_chunks.r:1-2` and `2.1.1.1_process_chunk.R:8-10`: `if (project_root == "") stop("PROJECT_ROOT not set — source slurm_config.sh")`
- `config_helper.R:36-37`: replace `home_ref_dir`/`gene_ref` fallbacks (`$HOME` path, `getwd()` gene file) with `stop()` guards

### B6. Step 1 — `config_helper.R` + cell type annotation simplification

- **Keep `config_helper.R`** at project root: R cannot source `slurm_config.sh` natively; the env-var-based `get_pipeline_config()` is the R-facing equivalent — cleaner than bash-side replacement
- Verify `1.1_prepare_chunks.r` / `2.1.1.1_process_chunk.R` use only `get_pipeline_config()` paths (already per-dataset under `${SCRATCH_OUTPUT_DIR}/${DS_NAME}`)
- Check `1_prepare_chunks.sh` / `2.1.1_process_chunk.sh` for further path centralisation into `slurm_config.sh`; `HOME_CHUNKS_DIR` (empty at `slurm_config.sh:55-56`, set in `2_submit_hpc_array.sh:17`) — leave as-is unless a clean lazy derivation exists (`DS_NAME` is not set when slurm_config.sh is sourced)
- Remove `gene_ref` remnants if confirmed unused (GENE_REF_FILE comment in slurm_config.sh:34-36)

### B7. Centralisation sweep + final docs pass

- Sweep ALL scripts, notebooks, `src/utils/` for hardcoded paths (`data/`, `logs/`, `aux/`, `EnsemblGenes...`, NAS paths) → route through `slurm_config.sh` env vars
- `src/benchmark/` handling: **pending decision** (local vs HPC, per TODO.md:31) — gates only benchmark-path handling; if decided later, needs a SLURM wrapper + figure-sync strategy. Note in plan/TODO, do not block on it
- Final docs update: AGENTS.md, README.md, docs/ARCHITECTURE.md, TODO.md (Completed section with all decisions)

---

## Out of scope (explicit)

- Unified single `1_submit_hpc_array.sh` across pipelines (TODO "Other major goals":44-48) — follow-up item
- New datasets, new methods (PILOT-GM-VAE/QOT/PULSAR), batch-effect method implementation
- `.kilo/plans/*.md` historical plan files

## Validation (final, after B7)

- `bash -n` all `.sh`; `py_compile` all modified `.py`; `Rscript -e 'parse(file="...")'` on modified `.R`/`.r` (syntax only)
- Grep sweeps for old paths and hardcoded dirs → clean
- No pipeline scripts executed (per AGENTS.md rule; validated on HPC later with the Joanito debug dataset)

## Risks

| Risk | Mitigation |
|---|---|
| Missed path references after rename | A5 + B7 grep sweeps; `git status`/`git log` review |
| Python import breakage from numbered dirs | Resolved by D2 (`src/utils/preprocess_utils.py`); verify `sys.path` insert still resolves in both importers |
| Dirty working tree during `git mv` | `git mv` handles modified files; re-verify diffs after move |
| Docs drift | All doc updates happen in A4 + B7, then final grep sweep |
