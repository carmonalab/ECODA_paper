# HPC Layout Docs + TODO Verification & Fixes + Dataset-Specific Preprocessing Dispatcher

## Goal

1. Document the HPC folder structure (what files go where) in `docs/ARCHITECTURE.md` + summary in `AGENTS.md`.
2. Verify the TODO.md "DONE" items against the code; fix the gaps found (preprocess array bug, TODO 3e, Joanito HPC step, dead dirs, legacy code, jq module).
3. Restructure `src/2_dataset_specific_preprocessing/` for HPC execution: per-step sbatch jobs + a parallel dispatcher.

## Decisions (confirmed with user)

- **HPC layout** (documented; assumed repo at `~/ECODA_paper` on HPC — only affects docs, `PROJECT_ROOT` derives automatically):
  - `$HOME/scratch/ECODA_paper` (`HPC_SCRATCH_DIR`): `<DS_NAME>/data/` (staged raw inputs), `CombinedPBMC/data/` (combine output + rds→h5ad cache), `output/<DS_NAME>/` (`SCRATCH_OUTPUT_DIR`: preprocessed .h5ad per view, `chunks/` during annotation, `annotations_chunk_*.feather`, annotated .h5ad).
  - `$HOME/reference_atlases/sketched_200ct/` (`HOME_REF_DIR`); `$PROJECT_ROOT/logs`, `aux/`, `.pixi/`.
  - NAS: `NAS_SC_DIR` (raw source), `NAS_REF_DIR`; `NAS_TARGET_DIR` = `Projects/ECODA_paper/` with `output/` (rsynced from `SCRATCH_OUTPUT_DIR`), `benchmark/{embeddings,plots}/`, `batch_effect_analysis/{embeddings,plots}/` (targets for method .feathers + notebook plots; filled once the `5_run_benchmark_methods` decision is made — TODO).
- **Remove dead dirs**: `path_plots`, `path_output_samples`, `path_output_ecoda` in `config_helper.R` + `1.1_prepare_chunks.r` (created but never written to; annotation feathers go to `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` directly).
- **Preprocess array fix**: add `--ds_name` filter to `1.1.1_preprocess.py`; `1.1_run_worker.sh` passes its `DS_NAME` (currently every worker loops over ALL datasets with its own base/output dirs → guaranteed missing-file failure).
- **TODO 3e**: unified per-view pipeline — HVG `batch_key="Sample"` (benchmark) / `batch_col` (batch effect); `X_pca` stored for every HVG size; Harmony + unsupervised clustering (neighbors + Leiden, on both `X_pca` and `X_pca_harmony`) only at the 2000-HVG pass (saves significant compute; reviewer point only needs the 2000-HVG cluster annotations). Batch views also get Leiden on `X_pca_harmony` (`BATCH_VIEW_RES = []` removed). Key names may be renamed safely (pipeline never run; no external consumers of `leiden_res_*`/`nobatch` names found).
- **Dataset-specific preprocessing**: per-step sbatch scripts (`1.1_submit_combinedpbmc.sh`, `1.2_submit_joanito_batch_col.sh`) + dispatcher `1_submit_hpc.sh` that submits all steps **in parallel** (`sbatch --parsable`), waits for all, reports per-job state via `sacct`. Steps must be mutually independent (documented; current two are).
- **Delete** `src/4_cell_type_annotation/DEPRECATED_LEGACY_CODE/` (stale `getwd()` fallbacks contradicting current guards).
- **`5_run_benchmark_methods` decision deferred** → record in TODO.md (user's plan: python methods in parallel on HPC; `benchmark_methods_r.R`/`benchmark_pipeline.R` run on HPC sequentially via a single-worker bash script e.g. `5_run_benchmark_methods/1_run_worker.sh`, subfolder TBD; data processing moved out of `notebooks/benchmark_analysis.rmd`).

---

## Tasks

### Task 1 — Preprocess rewrite (`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`)

Unified per-view pipeline — **both view types get the same treatment**; only the HVG `batch_key` differs (benchmark: `"Sample"` from `SAMPLE_COLNAME` env, batch-effect: `batch_col`). `X_pca` is stored for every HVG size; Harmony (`X_pca_harmony`) and unsupervised clustering (neighbors + Leiden) run **only at the 2000-HVG pass** (constant `CLUSTER_N_HVG = 2000`), on both the uncorrected and harmony-corrected embeddings.

1. **Split `run_downstream_for_gene_set` (lines 69-102) into three small functions** — the subset/scale/PCA part is shared and runs once per gene set, whether harmony is computed or not:
   - `compute_pca_and_store(adata, genes, key_suffix, n_pcs=50)` → `sub = adata[:, adata.var_names.isin(genes)].copy()`; `sc.pp.scale(sub, max_value=10)`; `sc.pp.pca(sub, n_comps=min(n_pcs, sub.n_vars-1, sub.n_obs-1), svd_solver="arpack")`; `adata.obsm[f"X_pca_{key_suffix}"] = sub.obsm["X_pca"]`; returns `sub`.
   - `compute_harmony_and_store(adata, sub, batch_key, key_suffix)` → `sc.external.pp.harmony_integrate(sub, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony")`; `adata.obsm[f"X_pca_harmony_{key_suffix}"] = sub.obsm["X_pca_harmony"]`.
   - `run_clustering(adata, rep_key, key_suffix, resolutions)` → `neighbors_key = f"neighbors_{key_suffix}"`; `sc.pp.neighbors(adata, n_pcs=adata.obsm[rep_key].shape[1], use_rep=rep_key, key_added=neighbors_key)`; loop `sc.tl.leiden(adata, resolution=r, key_added=f"leiden_res_{r}_{key_suffix}", neighbors_key=neighbors_key)`.
2. **`select_hvgs_ranked` gains a `batch_key` parameter** (passed through to `sc.pp.highly_variable_genes`); delete now-redundant `compute_hvgs` (lines 18-33).
3. **`process_view` becomes unified** (drop `use_harmony` flag / `BATCH_VIEW_RES`):
   - `select_hvgs_ranked(adata, n_top_genes=max(n_hvg_sizes), batch_key=batch_key)`.
   - Per HVG size `n` (key_suffix = `f"{view_name}_hvg{n}"`):
     - `sub = compute_pca_and_store(adata, genes, key_suffix)` — always.
     - Only when `n == CLUSTER_N_HVG` (2000):
       - `run_clustering(adata, f"X_pca_{key_suffix}", key_suffix, resolutions)`.
       - `compute_harmony_and_store(adata, sub, batch_key, key_suffix)`.
       - `run_clustering(adata, f"X_pca_harmony_{key_suffix}", f"{key_suffix}_harmony", resolutions)`.
   - Resulting keys (example, benchmark view): `X_pca_benchmark_analysis_hvg3000`, `X_pca_benchmark_analysis_hvg2000`, `X_pca_benchmark_analysis_hvg1000` (PCA only); at 2000 additionally `X_pca_harmony_benchmark_analysis_hvg2000`, `leiden_res_{r}_benchmark_analysis_hvg2000`, `leiden_res_{r}_benchmark_analysis_hvg2000_harmony` (batch views analogously with `batch_effect_analysis_hvg2000`). Old `_nobatch_hvg{n}` / `_batch_hvg2000` suffixes replaced.
4. **`main()` changes**:
   - `--ds_name` filter: new argparse arg + `main()` parameter; when set, skip all other datasets in the `config.items()` loop (view-less datasets like Zhu still exit via "No views defined.").
   - `batch_key` per view type: benchmark views → `os.environ.get("SAMPLE_COLNAME", "Sample")`; batch-effect views → `batch_col` (keep the existing `batch_col not in obs` guard, lines 225-229).
   - `n_hvg_sizes` per view type: benchmark `(3000, 2000, 1000)` unchanged; batch-effect `[2000]` (keep `BATCH_VIEW_N_HVG`, single size — easy to extend). Resolutions `(0.1, 0.4, 2, 5, 20, 50)` for both.
   - Harmony + clustering now always computed at the 2000-HVG pass for every view (no `use_harmony` argument).
5. Update module docstring/header comments to document the per-view output keys and the clustering-size rule.

### Task 2 — Worker passes dataset (`src/3_scrnaseq_preprocessing/1.1_run_worker.sh`)

- Add `--ds_name "${DS_NAME}"` to the `1.1.1_preprocess.py` invocation (line ~35).
- Note: view-less datasets (Zhu) exit cleanly via existing "Skipping ... No views defined." path — resolves TODO.md's optional cleanup item.

### Task 3 — Dataset-specific preprocessing HPC steps (`src/2_dataset_specific_preprocessing/`)

1. `git mv 1_submit_hpc.sh 1.1_submit_combinedpbmc.sh`; keep content (sources `slurm_config.sh`, `cd ${PROJECT_ROOT}`, `module load GCCcore/12.2.0`, runs `_create_combinedpbmc_dataset.py`).
2. New `1.2_submit_joanito_batch_col.sh`:
   - `#SBATCH`: `--job-name=joanito_batch_col`, shared-cpu, ~`--time=01:00:00`, `--mem=8G`, `--ntasks=1`.
   - Sources `slurm_config.sh`, `cd "${PROJECT_ROOT}"`, runs `"${HOME}/.pixi/bin/pixi" run Rscript --vanilla "${SCRIPT_DIR}/_create_joanito_batch_col.R"` (mirrors `2.1.1_process_chunk.sh` pixi invocation).
3. Adapt `_create_joanito_batch_col.R` (HPC-capable):
   - Input path from env: `input <- file.path(Sys.getenv("HPC_SCRATCH_DIR"), "Joanito", "data", "JoaI_2022_35773407_Nofilt_whole.rds")`; `stop()` guard if `HPC_SCRATCH_DIR` unset.
   - Keep in-place `saveRDS` (idempotent — recomputes `seqtec`); must run **after** `1_stage_data.sh`, **before** the preprocess array.
4. New `1_submit_hpc.sh` dispatcher (replaces the old single-step script):
   - Sources `slurm_config.sh`, `cd "${PROJECT_ROOT}"`; no module loads (steps own their env).
   - Submits all `1.*_submit_*.sh` step scripts in the folder via `sbatch --parsable` (array `STEP_SCRIPTS=(...)`), collects job IDs.
   - Waits for all (squeue loop on collected IDs), then prints `sacct -j <ids> --format=JobID,JobName,State,ExitCode` summary; exits non-zero if any step failed.
   - Header comment: steps run in parallel and MUST be independent; a future dependent step must be submitted after the wait (or via `--dependency`).

### Task 4 — Dead dirs removal

- `config_helper.R`: remove `path_plots`, `path_output_samples`, `path_output_ecoda` (keep `path_data`, `path_output`, `path_output_chunks`, `path_ref`).
- `src/4_cell_type_annotation/1.1_prepare_chunks.r`: remove the corresponding `dir.create(...)` calls (lines ~39, ~45).

### Task 5 — Cleanup + minor fixes

- `git rm -r src/4_cell_type_annotation/DEPRECATED_LEGACY_CODE/`.
- `src/4_cell_type_annotation/2_submit_hpc_array.sh`: add `module load jq/1.6` (with `command -v jq` fallback retained) before the `TISSUE_TYPE`/`NORMAL_TISSUE` export block — currently jq is likely absent on the HPC login node without the module, silently skipping the export.

### Task 6 — Documentation updates

- `docs/ARCHITECTURE.md`:
  - New section "HPC Folder Layout" (tree + per-folder "what goes here" table incl. `CombinedPBMC/data/` dual role, `output/<DS_NAME>/` contents, NAS target structure).
  - Update preprocessing table: `1.1.1_preprocess.py` row (view output keys: `X_pca_{view}_hvg{n}` for every HVG size; `X_pca_harmony_{view}_hvg2000`, `neighbors_{...}`, `leiden_res_{r}_{view}_hvg2000` + `leiden_res_{r}_{view}_hvg2000_harmony` only at the 2000-HVG pass; HVG batch_key = Sample/batch_col), `1.1_run_worker.sh` row (`--ds_name`).
  - Update `src/2_dataset_specific_preprocessing/` row + CombinedPBMC bullet: dispatcher + per-step scripts + Joanito step; remove `1_submit_hpc.sh` mention.
- `AGENTS.md`: brief HPC layout summary in "HPC general information"; update Stage 2 bullet (dataset-specific steps via `1_submit_hpc.sh` dispatcher; Joanito step required before preprocess); note annotation dirs (`samples/`, `ecoda/`, `plots/` no longer created).
- `README.md`: usage section — dataset-specific preprocessing via `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh`; preprocess array unchanged.
- `TODO.md`:
  - Mark as done: TODO 3e (benchmark views batch-aware HVG + harmony-based Leiden; batch views Leiden on harmony), `--ds_name` array fix, `_create_joanito_batch_col.R` HPC-capable + `1.2_submit_joanito_batch_col.sh`, dispatcher `1_submit_hpc.sh`, dead-dir removal, `DEPRECATED_LEGACY_CODE` deletion, jq module fix, HPC layout docs.
  - Add "later" point: `src/5_run_benchmark_methods/` — `run_python_sample_embedding_methods/` run in parallel on HPC; `benchmark_methods_r.R` + `benchmark_pipeline.R` run on HPC sequentially via a simple single-worker bash script (e.g. `5_run_benchmark_methods/1_run_worker.sh`, subfolder structure TBD); move this data-processing step out of `notebooks/benchmark_analysis.rmd`; NAS targets `benchmark/{embeddings,plots}/` + `batch_effect_analysis/{embeddings,plots}/`.

## Validation (no pipeline runs, per AGENTS.md)

- `bash -n` on all modified/new `.sh`.
- `python3 -m py_compile` on modified `.py`.
- `Rscript -e 'parse(file="...")'` on modified `.R`/`.r` (syntax only).
- Grep sweeps clean: `nobatch`, `compute_hvgs`, `use_harmony`, `BATCH_VIEW_RES`, `path_output_samples`, `path_output_ecoda`, `path_plots` (outside `DEPRECATED_LEGACY_CODE` which is deleted), `1_submit_hpc.sh` (now the dispatcher), `DEPRECATED_LEGACY_CODE`.
- Verify `git status` shows the renames/deletions as expected.

## Assumptions / Open

- Repo lives at `~/ECODA_paper` on the HPC (docs only; `PROJECT_ROOT` derives automatically via `slurm_config.sh`).
- `5_run_benchmark_methods` HPC/local decision deferred (recorded in TODO.md).
- New obsm/obs key names are free to change (pipeline never run; no external consumers found).

## Risks

| Risk | Mitigation |
|---|---|
| Renamed obsm/obs keys break later R/python consumers | Pipeline not yet run; grep sweep for old names in Task 6; naming documented in ARCHITECTURE.md |
| Parallel dispatcher steps not independent in the future | Header comment + `sacct` failure report; dependent steps must be sequenced |
| Joanito rds modified in place while preprocess array reads it | Order enforced: `1_stage_data` → `1_submit_hpc.sh` → `3_scrnaseq_preprocessing` (documented); in-place edit is idempotent |
| jq still missing on login node | `module load jq/1.6` with fallback retained; other scripts already load the module |
