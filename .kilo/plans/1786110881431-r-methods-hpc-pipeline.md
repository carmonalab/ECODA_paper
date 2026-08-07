# R methods HPC pipelines + benchmark_analysis split

## Goal

1. Move the heavy R benchmark methods (GloScope, MOFA, Pseudobulk, scITD) off the desktop onto an HPC pipeline, mirroring `run_python_sample_embedding_methods/`.
2. Split the old `run_analyses()` (in `src/5_run_benchmark_methods/benchmark_pipeline.R`) into separate scripts:
   - `run_benchmark_analysis` → HPC pipeline (methods above).
   - `run_transformation_analysis` + `run_zeroimp_analysis` → one fast pipeline with **two workers** (one per analysis).
3. GloScope: keep only the sqrtmat variant, drop the non-sqrtmat variant, rename results `GloScope_hvg2000_pcadims30_sqrtmat` → `GloScope_hvg2000_pcadims30`.
4. Precompute (and save) the shared DESeq2 pseudobulks *before* MOFA/Pseudobulk run; method timing must include pseudobulk creation + method runtime (as the old code did: `exec_time(method) + exec_time_pb_norm`).
5. Minimal notebook adaptation so results are consumable end-to-end (full h5ad notebook rework stays in TODO Phase 3.4).

## Confirmed decisions (user)

- **Scope**: HPC pipeline A runs `gloscope`, `mofa`, `pseudobulk`, `scitd` (+ `prepare_pseudobulk` prep step). ECODA variants, GloProp, EPIC deconv, Avg_PCA_embedding, Freq_highres stay in the notebook (fast, composition-based).
- **Structure**: separate pipeline dir (not merged with the Python one); shared exec-log schema + merge script + NAS `benchmark/` target reused.
- **Pseudobulk**: keep the existing R calculation (`get_pb_deseq2` / DESeq2.normalize) — do NOT switch to scanpy.
- **Notebook**: minimal adaptation included now.

## Key design

- New dir `src/5_run_benchmark_methods/run_r_sample_embedding_methods/` (Pipeline A) and `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/` (Pipeline B).
- One SLURM array per method (Pipeline A: `prepare_pseudobulk` first, monitored to completion, then `gloscope`/`mofa`/`pseudobulk`/`scitd`) — mirrors the Python pipeline (per-method manifests, `sed -n ${SLURM_ARRAY_TASK_ID}p`, sacct fail-closed gate, NAS check, exec-log merge via the existing `run_python_sample_embedding_methods/1.1.2_merge_execution_times.py`, rsync, log cleanup).
- Pipeline B: two arrays (`trans`, `zeroimp`) submitted by one submitter → two workers per dataset.
- R workers run via `${PIXI_RSCRIPT}` (already exported by `slurm_config.sh`), source `src/utils/load_all_functions.R` + a new `benchmark_hpc_utils.R`, read the preprocessed benchmark-view h5ad via reticulate, save RDS bundles + per-task exec-log feathers (same schema as Python: `dataset, method, time_secs, mem_GB` with `mem_GB = NA`).
- Outputs (all under `${HPC_SCRATCH_DIR}/benchmark/`, rsynced to `${NAS_TARGET_DIR}/benchmark/`):
  - `results/<ds>_<method>.rds` — method-level named list of result bundles (bundle = `create_result_bundle()` + `exec_time` in seconds).
  - `results/<ds>_<method>_<combo>.rds` — per-combo cache files (skip-if-exists unless `--force`, failure-resume for 12h tasks).
  - `results/<ds>_trans.rds`, `results/<ds>_zeroimp.rds`.
  - `pseudobulks/<ds>_pseudobulk_<variant>.rds` — `list(pb = <samples x genes DESeq2-normalized matrix>, time_secs = <creation time>)`.
  - `gloscope_dists/<ds>_gloscope_hvg<n>_pcadims<d>_dists.rds` — raw GloScope distance cache (sqrt applied at processing).
  - `embeddings/execution_times_task_<ARRAY_JOB_ID>_<ARRAY_TASK_ID>.feather` — per-task exec logs (merged by the shared merge script).
- Hardware: Pipeline A pinned to the CPU benchmark class (`SLURM_PARTITION_BENCHMARK_CPU`, `BENCHMARK_CPU_CONSTRAINT` EPYC-7742, 16 cpus, 128G — same as PILOT, so cross-method runtime comparison stays valid; `--partition` override drops the constraint pin, as in the Python submitter). Pipeline B: `shared-cpu`, 4 cpus, 32G, no pin (not part of runtime comparisons).

## Task list

### T1. pixi.toml — add missing R packages (worker needs the full `load_all_functions.R` set + DESeq2)

Verified against `pixi.lock` (Aug 2026):

- **Missing from the lock file — must add to `[dependencies]`**: `r-vegan`, `r-robcompositions`, `bioconductor-deseq2` (pulls SummarizedExperiment/MatrixGenerics transitively; leave unversioned or `>=1.46,<2` and let bioconda resolve the r-base 4.5.2-compatible `r45` build rather than pinning `1.46.*` by hand), `r-doparallel`, `r-plotly`, `r-pheatmap`.
- **Already present transitively (no action needed, though explicit pins are harmless)**: r-foreach, bioconductor-biocparallel, r-cluster, r-parallelly, r-jsonlite, r-proc, r-zoo, r-patchwork, r-ggpubr, r-ggrepel, r-rstatix, r-scales, bioconductor-matrixgenerics.
- NOTE: `DESeq2`, `vegan`, `robCompositions` are NOT in `imports.R`'s `require()` list (used via `::` and the foreach `.packages` arg), so the T8 `source("load_all_functions.R")` check does NOT validate them — T8 must additionally run an explicit `requireNamespace(c("DESeq2", "vegan", "robCompositions"))` check.
- Verify locally with `pixi install` + `pixi run -e py-cpu Rscript --vanilla -e 'source("src/utils/load_all_functions.R")'` + the requireNamespace check above. Note: HPC env update runs later via the standard README pixi setup (py-cuda13; GitHub-installed pkgs — Seurat, anndataR, MOFA2, scITD, GloScope, scECODA, STACAS, thisutils, ... — come from `[tasks.setup]`, not `pixi install`).

### T2. R module refactors (shared by notebook + HPC)

- `src/5_run_benchmark_methods/benchmark_methods_r.R`:
  - Delete `process_gloscope_fig()` (non-sqrtmat). Rename `process_gloscope_sqrtmat_fig()` → `process_gloscope_fig()` with new signature `(embedding_matrix, sample_ids, metadata, label_col, gloscope_dist_file, n_pca_dims = 30, dens = "KNN", dist_metric = c("KL"), k = 25)` (no longer takes a Seurat object): computes `GloScope::gloscope()` on `embedding_matrix[, 1:n_pca_dims]` (cached as raw `gloscope_dist_file`), `sqrt()` + NA→0, standardizes sample names, returns bundle. Robustness:
    - `k <- min(k, n_samples - 1)` with a warning when adjusted (needed for the 5-sample `_debug` dataset; real datasets have ≥23 samples so k=25 unchanged).
    - `n_pca_dims <- min(n_pca_dims, ncol(embedding_matrix))` with a warning: the obsm stores `min(50, n_vars-1, n_obs-1)` PCs (scanpy `compute_pca_and_store`), so `_debug` (5 samples) has only 4 PCs and pcadims 10/30/50 would otherwise crash.
    - Timing semantics (matches legacy): on gloscope-cache miss the combo time includes dist computation; on cache hits it is sqrt+read only — same behavior as the old notebook's `path_data` dist-cache.
  - Keep `process_mofa_bulk_fig`, `process_pseudobulk_fig`, `process_pseudobulk_ct_fig`, `process_scitd_fig` unchanged (HPC workers call them).
- `src/5_run_benchmark_methods/benchmark_pipeline.R`:
  - `run_benchmark_analysis()`: remove the Pseudobulk, MOFA, GloScope and scITD blocks (keep everything else: ECODA variants, GloProp, Avg_PCA, Freq_highres, python-feather ingest + missing-file check, `seurat@misc` usage). Specifically:
    - **Drop the `ds_filename <- "GongSharma_all"` remap** (lines 237–241): both the HPC Python pipeline and the legacy qmd write feathers named with the datasets.json key (`Gongsharma_cmv_young_males_hvg{n}_mrvi_dists.feather` etc.), so the remap would fail the missing-file check for GongSharma. All feather lookups use the datasets.json key directly.
    - **Keep the local `pb_norm <- get_pb_deseq2(seurat, sample_col, hvg = NULL, n_hvg = 2000)` computation**: still required by `ECODA_deconv` (`process_deconv_fig(t(pb_norm), labels)`). Do NOT wrap it in `exec_time` (legacy ECODA_deconv exec_time excluded pb creation).
    - In the `factors_test` loop keep ONLY `ECODA_authors_HR_*_PCA_dims`; `Pseudobulk_*_PCA_dims` move to HPC.
    - Remove the now-unneeded `seurat@version <- package_version("3.1.5")` hack (MOFA-only) and the `hvg <- get_current_hvgs(seurat)` variable if nothing retained uses it.
  - Add HPC driver functions (ported from the removed blocks, preserving exact combo lists and result names):
    - `prepare_pseudobulks_hpc(seurat, sample_col, hvg_rank_genes)` → variants `schvg2000` (`get_pb_deseq2(hvg = top-2000 hvg_rank genes)`), `hvg2000`, `hvg500`, `hvg2000_bl` (STACAS `data("default_black_list")`, `black_list = "default_without_sex_genes"` — **preserve current behavior including the existing `%in% black_list` no-op typo at pseudobulk.R:84**, do not silently fix), `hvg1000`, `hvg3000`. Returns per-variant `list(pb, time_secs)`.
    - `run_gloscope_hpc(seurat, metadata, label_col, embedding_key = "pca_benchmark_analysis_hvg{n}", cache_dir, ...)` → combos hvg2000 × pcadims {10,30,50}; hvg1000, hvg3000 × pcadims 30.
    - `run_mofa_hpc(seurat, metadata, labels, pb_variants, ...)` → `MOFA_hvg2000_factors{2,3,5,10,15}`, `MOFA_hvg{1000,3000}_factors15` (each time = pb creation time + MOFA runtime). Skip combos with `num_factors >= n_samples` (warning; makes `_debug` usable).
    - `run_pseudobulk_hpc(seurat, labels, pb_variants, ...)` → `Pseudobulk_schvg2000`, `Pseudobulk_hvg2000`, `Pseudobulk_hvg500`, `Pseudobulk_hvg2000_bl`, `Pseudobulk_CT_LR_hvg2000`, `Pseudobulk_CT_HR_hvg2000`, `Pseudobulk_CT_LR_hvg500`, `Pseudobulk_CT_HR_hvg500` (gate each CT variant on its own ct col being non-null — note: legacy gates `CT_LR_hvg500` on `cell_type_high_res` at benchmark_pipeline.R:403-405; the per-col gating is cleaner and behaviorally identical since datasets.json always defines both columns), `Pseudobulk_{2,3,5,10,15}_PCA_dims`, `Pseudobulk_hvg1000`, `Pseudobulk_hvg3000`.
    - `run_scitd_hpc(seurat, label_col, hvg_sets, ...)` → `scITD_hvg2000_factors{2,3,5,10,15}`, `scITD_hvg{1000,3000}_factors5` (uses `cell_type_low_res` ct col, as today). Skip combos where `num_factors + 5 >= n_samples` (tucker rank `c(f, f+5)` must be < n_samples; `_debug` with 5 samples cannot run scITD at all — validate scITD on a real dataset, e.g. Kim).
    - `load_hpc_benchmark_results(result_list, ds, path_results_nas, methods = c("gloscope","mofa","pseudobulk","scitd"))` — notebook-side: preserves legacy rerun semantics (entries already present in `result_list$bmark[[ds]]` are kept; only missing methods read from NAS), reads `<ds>_<method>.rds` bundles into `result_list$bmark[[ds]]` (warning on missing file), reads `<ds>_trans.rds`/`<ds>_zeroimp.rds` into `result_list$trans[[ds]]`/`$zeroimp[[ds]]`, saves `result_list.rds`, returns the list. Replaces `run_analyses()` (delete it; `run_transformation_analysis`/`run_zeroimp_analysis` stay, now only called by Pipeline B workers).
- `src/utils/constants.R`: `method_label_map_main` entry `GloScope_hvg2000_pcadims30_sqrtmat` → `GloScope_hvg2000_pcadims30`.
- New `src/5_run_benchmark_methods/benchmark_hpc_utils.R` (NOT added to `load_all_functions.R`; sourced by HPC scripts): tiny `--flag value` arg parser, `get_h5ad_path(config, ds, view)` (reuse `read_datasets_json` from `datasets_io.R`, which already maps `columns.label` → `label_col`, `output_file_name` → `output_file`, etc.), obs-only backed read + full-counts Seurat build via reticulate/anndata (`get_seurat_obj_from_h5ad()` from `seurat_utils.R` with `counts_layer = "counts"` + `fetch_embedding = c("X_pca_benchmark_analysis_hvg1000","X_pca_benchmark_analysis_hvg2000","X_pca_benchmark_analysis_hvg3000")`; embeddings land in reductions named `pca_benchmark_analysis_hvg{n}`), `VariableFeatures(seurat) <- top-n hvg_rank genes` (from `adata$var["hvg_rank"]`, mirroring the python `top_n_hvg_genes` — also needed by the `prepare_pseudobulk` worker for `schvg2000`), Sample-name standardization (`standardize_sample_names`) before computing labels/metadata, and `log_exec_row(dataset, method, time_secs, log_file)` writing the shared exec-log feather schema with read-modify-write append semantics (mirror `log_execution_time` in `1.1.1_benchmark_methods_py.py`: read existing per-task feather if present, append, dedup on (dataset, method) keep=last).

### T3. Pipeline A — `src/5_run_benchmark_methods/run_r_sample_embedding_methods/`

- `1_submit_hpc_array.sh`: `[--ds_name <DS>] [--methods gloscope,mofa,pseudobulk,scitd,prepare_pseudobulk] [--partition <P>] [--force]`. Dataset resolution identical to the Python submitter (jq: `use_for_benchmark == true` + `benchmark_analysis` view, skip `_*` unless `--ds_name`). Validate `--methods` against the known set (error on unknown, like the Python submitter). If `mofa` or `pseudobulk` requested and `prepare_pseudobulk` not explicitly listed → auto-prepend it. Submit order: `prepare_pseudobulk` array → poll `squeue` until done → sacct gate (fail-closed) → submit the other arrays → poll all → sacct gate → NAS reachability check → merge exec logs via `${PYTHON_BIN} run_python_sample_embedding_methods/1.1.2_merge_execution_times.py --output_dir ${HPC_SCRATCH_DIR}/benchmark/embeddings --no-cleanup --job_ids <ids> [--existing-log NAS log]` → `rsync ${HPC_SCRATCH_DIR}/benchmark/ ${NAS_TARGET_DIR}/benchmark/` → delete per-task logs. Per-method manifest `${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_$$.txt`; export `METHOD`, `BENCHMARK_MANIFEST`, `FORCE_BENCHMARK`. All methods pinned to the CPU class (constraint dropped on `--partition` override); throttle `MAX_NUM_CHUNKS_PARALLEL`; logs `${LOGS_DIR}/5_benchmark_r_<method>_%A_%a.log/.err`.
- `1.1_run_worker.sh`: same SBATCH header (12h, 16 cpus, 128G defaults)/`scontrol show job` SCRIPT_DIR recovery/`slurm_config.sh` + `cd ${PROJECT_ROOT}` boilerplate as `run_python_sample_embedding_methods/1.1_run_worker.sh`; requires `METHOD` + `BENCHMARK_MANIFEST`; DS_NAME from manifest; calls `${PIXI_RSCRIPT} 1.1.1_run_benchmark_methods_r.R` (or `1.1.1_prepare_pseudobulk.R` for `prepare_pseudobulk`) with `--config_path --ds_name --view benchmark_analysis --input_dir ${HPC_SCRATCH_DIR}/${DS_NAME}/output --results_dir ${HPC_SCRATCH_DIR}/benchmark/results --pseudobulk_dir .../benchmark/pseudobulks --gloscope_cache_dir .../benchmark/gloscope_dists --log_file ${HPC_SCRATCH_DIR}/benchmark/embeddings/execution_times_task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.feather [--force]`.
- `1.1.1_run_benchmark_methods_r.R`: sources `load_all_functions.R` + `benchmark_hpc_utils.R`; loads the h5ad → Seurat + obsm embeddings (reticulate path above); sets `seurat@misc$cell_type_low_res`/`label_col` from datasets.json; dispatches on `--method` to the T2 driver; writes per-combo bundle files + method-level RDS; per-combo exec-log rows; skip-if-exists semantics (method RDS exists → skip all unless `--force`; else per-combo files reused).
- `1.1.1_prepare_pseudobulk.R`: loads Seurat (counts + `var["hvg_rank"]` only; no embeddings), runs `prepare_pseudobulks_hpc()`, writes `pseudobulks/<ds>_pseudobulk_<variant>.rds` (atomic tmp+rename) and exec-log rows `prepare_pseudobulk_<variant>`.
- MOFA/Pseudobulk workers load pb variants from `pseudobulks/`; if a variant RDS is missing (e.g. `--methods mofa` without prep on a stale scratch), compute it on the fly (save via tmp+rename) and still report time = pb time + method time.

### T4. Pipeline B — `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/`

- `1_submit_hpc_array.sh`: `[--ds_name <DS>] [--analysis trans,zeroimp] [--partition <P>] [--force]`; two arrays (one per analysis; manifests `benchmark_manifest_<analysis>_$$.txt`; export `ANALYSIS`, `BENCHMARK_MANIFEST`, `FORCE_BENCHMARK`); same monitor/gate/merge/sync/cleanup tail as T3. Partition `shared-cpu`, 4 cpus, 32G, `--time=02:00:00`, throttle `MAX_NUM_CHUNKS_PARALLEL`.
- `1.1_run_worker.sh`: boilerplate as above; branches on `ANALYSIS` (`trans` → `1.1.1_run_transformation_analysis.R`, `zeroimp` → `1.1.1_run_zeroimp_analysis.R`); args `--config_path --ds_name --view benchmark_analysis --input_dir --output_dir ${HPC_SCRATCH_DIR}/benchmark/results --log_file ...`.
- `1.1.1_run_transformation_analysis.R` / `1.1.1_run_zeroimp_analysis.R`: obs-only backed read (reticulate, no counts matrix); `ct_comps = table(Sample_std, obs[[cell_type_high_res]])` cast to `as.data.frame.matrix` (rowSums != 0 filter, as `get_ct_comp_df_seurat`; data.frame so the dplyr verbs in `run_zeroimp_analysis` work), `labels = get_labels`-equivalent (per-sample slice(1) of `label_col`, names = Sample); call `run_transformation_analysis(ct_comps, labels)` / `run_zeroimp_analysis(ct_comps, labels)`; save `<ds>_trans.rds` / `<ds>_zeroimp.rds`; one exec-log row per dataset (`trans_analysis`/`zeroimp_analysis`).

### T5. slurm_config.sh

- Comment-only update: note that the R benchmark pipeline reuses the CPU benchmark class (EPYC-7742, 16 cores, 128G) for runtime comparability with PILOT.

### T6. Minimal notebook adaptation (`notebooks/benchmark_analysis.rmd`)

- Replace the `run_analyses(result_list, ds, seurat, path_data, path_plots)` call in the per-dataset loop with: `result_list <- load_hpc_benchmark_results(result_list, ds, path_nas_benchmark)` (new `path_nas_benchmark` variable, default `/Volumes/Shared/Projects/ECODA_paper/benchmark/results`, documented + easy to adapt to the user's mount) followed by the (now light) `run_benchmark_analysis(...)`. Stephenson `Site == "Ncl"` subset, `replace_HiTMElayer3_annot()`, `seurat@misc` setup stay unchanged.
- **Exec-times block: filter `p_exec_times` to python-only method rows** (e.g. `grepl("^(MrVI|scPoli|PILOT)", method)` before the existing recodes): the merged `execution_times.feather` will now also contain R-method rows (gloscope/mofa/pseudobulk/scitd/prepare_pseudobulk_*, trans_analysis, zeroimp_analysis). Without the filter, R methods appear twice in `exec_times.rds` (bundle-derived `r_exec_times` + feather rows), double-counting points and skewing the Supp-fig-14 lm fits. (Bundle `exec_time` stays the single source for R-method timing.)
- GloScope rename everywhere: `GloScope_hvg2000_pcadims30_sqrtmat` → `GloScope_hvg2000_pcadims30` (filter_methods lists, recodes, default_methods), simplify the Supp-fig-2 `grepl("sqrtmat", method)` filter (remove the `| grepl("sqrtmat", method)` clause and the `_sqrtmat` str_replace).

### T7. Docs

- `docs/ARCHITECTURE.md`: update the "6 sbatch workers" list (→ 8, incl. the two new workers in the `scontrol show job` recovery note); add the two new pipeline sections (workflow diagrams, file-role tables, usage) under "Benchmark, ECODA Transformation and ECODA Zero Imputation Analyses"; mark the Layers 1–5 note as superseded for the R side (the existing placeholder at lines ~187–190 already anticipates this); update the benchmark HPC folder layout (`results/`, `pseudobulks/`, `gloscope_dists/` under `benchmark/`); update the call-flow diagram (run_analyses split).
- `TODO.md`: mark Phase 3.2 code-complete (HPC debug validation pending), add changelog entry.
- `README.md`: usage for the two new pipelines.
- `AGENTS.md`: Stage 3 bullet + worker-list update (pointers only).

### T8. Validation (no pipeline runs — per AGENTS.md)

- `bash -n` on all new/modified shell scripts.
- R parse check on all new/modified R files: `Rscript --vanilla -e 'parse(file = "<file>")'` (local R or the pixi env).
- `pixi install` locally (validates T1 resolution), then `pixi run -e py-cpu Rscript --vanilla -e 'source("src/utils/load_all_functions.R")'` to confirm the full imports.R set loads (do this on the Mac env; HPC py-cuda13 env is updated by the user via the standard setup), **plus** `requireNamespace(c("DESeq2","vegan","robCompositions"))` (not covered by load_all_functions — see T1).
- No changes to `datasets.json`.
- HPC debug smoke tests (user-executed, after Phase 2 prerequisites): Pipeline B first (`--ds_name _debug --analysis trans,zeroimp`), then Pipeline A `--ds_name _debug --methods prepare_pseudobulk,pseudobulk` (gloscope works via the k-adjust + n_pca_dims clamp, but note `_debug` has only 4 PCs — validate pcadims 10/30/50 on a real dataset; MOFA runs only factors 2,3 — 5+ skipped by the n_samples guard; scITD cannot run on 5 samples (tucker ranks) — validate on a real dataset, e.g. Kim). Check `benchmark/results/`, `pseudobulks/`, `execution_times.feather` on NAS.

## Notes / assumptions

- The preprocessed view h5ad is the single input (counts layer + obs + `X_pca_benchmark_analysis_hvg{n}` obsm); the legacy R-side PCA/FindVariableFeatures recomputation is replaced by stored obsm + `var["hvg_rank"]` (consistent with how PILOT/MrVI/scPoli already consume preprocessing; gene sets may differ slightly from the legacy Seurat `var.features` — mainly affects `schvg2000` and scITD hvg1000/hvg3000 — which is expected; the new run is internally consistent).
- Result bundles keep the exact legacy names (minus the GloScope `_sqrtmat` suffix); the notebook's exec-time section (reads bundle `exec_time`, numeric seconds) and all plot code work unchanged on the loaded bundles.
- `Pseudobulk_hvg2000_bl` behavior (incl. the pre-existing black-list typo in `get_pb_deseq2`) is preserved for legacy-result equivalence; flag to the user, do not fix silently.
- The RDS files are written with the default R serialization; loaded on the user's Mac (R 4.5.x both sides — compatible).
- `Gongsharma_cmv_young_males` files use the datasets.json key everywhere (T2 drops the legacy `GongSharma_all` feather-name remap so HPC-produced feathers are found); the notebook's existing `GongSharma` rename covers plotting.
- `datrans` (`run_transformation_analysis`) hardcodes `n_cores = 1`; scECODA/STACAS/thisutils etc. come from the `[tasks.setup]` installs (not pixi.toml) — the py-cuda13 env needs a `pixi run -e py-cuda13 setup` after `pixi install` for the HPC workers (standard README workflow).
