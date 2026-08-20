# Implement PILOT-GM-VAE + QOT benchmark methods (with _debug NAS testing); defer PULSAR/UCE

## Context & decisions (verified against code/env)

- **PILOT-GM-VAE** = PyPI package `pilotgm` 0.1.1 (CostaLab/PILOT-GM-VAE, BIB 2025). API used: `train_gmvae(adata, dataset_name, pca_key, epochs, num_classes, cuda, ...)` then `gmmvae_wasserstein_distance(adata, emb_matrix, sample_col, status, wass_dis=True)` → `adata.uns['EMD_df']` (sample×sample dist matrix, index = sample IDs) + `adata.uns['EMD']`. All deps already in the lockfile — only `pilotgm` itself is missing: torch (2.12.0 in the local py-cpu env; lockfile holds per-platform/per-feature variants in the 2.9–2.13 range — no version claim needed), numba 0.65.1, joblib 1.5.3, tqdm, pilotpy 2.0.15, scanpy, sklearn 1.9.0.
- **QOT** = PennShenLab/QOT `qot_utils_re.py` (commit `28cd529880c1`, 2024-12-04). Only module-level dep missing from our env: **`phate`** (only used by `trajectory_analysis`, which we never call, but it is a module-level import — must be installed). Everything else present and version-compatible: POT 0.9.6.post1 (incl. `ot.gaussian.bures_wasserstein_distance`, `ot.emd`), statsmodels 0.14.6, adjustText 1.4.0, gprofiler 1.0.0, scikit-learn 1.9.0 (`rand_score` still importable from `sklearn.metrics.cluster`), seaborn 0.13.2 (`ci=None` only warns), numpy 2.4.6 (`pd.DataFrame.from_dict(ndarray)` OK). API used: `Run_QOT(adata, gene_matrix, type_cell, id_col, progession, dataset_type='rna', num_components_list=[1], random_state=2, min_samples_for_gmm=0, qot_method='cosine', normalized_set=False)` → `adata.uns['QOT_Distance']`.
- **The 3–5-sample bug** (user-confirmed recollection): in `Gaussian_Mixture_Representation` the chain `if num_samples >= num_components + min_samples_gmm: GMM … elif num_samples == 1: … elif num_samples == 2: …` has NO else — (sample, cell-type) groups with `3 <= num_samples < num_components + min_samples_gmm` are silently dropped from `params`, leaving 0-distance rows/cols in `QOT_Distance`. Fix = add one `else` branch reusing the single-Gaussian representation (mean + zero covariance + weight 1, same as the `== 2` branch). Not triggered at our settings (num_components=1, min_samples=0 → threshold 1), but the user wants it fixed (harmless, one branch).
- **Duplicate-key rename bug in BOTH upstream APIs** (found during plan review, verified by simulation with the repo's py-cpu pandas): QOT's `Extract_Info` builds `rename(columns={type_cell: 'Cell_type', id: 'sampleID', progession: 'status'})` and pilotgm's `gmmvae_wasserstein_distance` builds `rename_dict = {current_columns[-3]: 'sampleID', current_columns[-2]: 'cell_type', current_columns[-1]: 'status'}`. When the id/sample column equals the progession/status column (both `"Sample"` in our call), the dict literals collapse the duplicate key to `{'Sample': 'status'}`, BOTH `Sample` columns get renamed to `status`, no `sampleID` column survives, and the GMM step (`groupby(['sampleID', 'Cell_type'])` / `df.drop(['sampleID'])`) raises KeyError. **Fix lives in OUR wrappers, not in vendored code**: create a distinct temp obs column (e.g. `adata.obs["_bench_prog"] = adata.obs["Sample"]`) and pass THAT as `progession` (QOT) / `status` (pilotgm). Also avoids passing any bio label (no-leakage).
- **Index-dropping step** (`Extract_Info` does `.reset_index(drop=True)` on PCA/obs columns before concat): verified positional-order-aligned → works for our use. Per user: do NOT fix.
- **PULSAR (+UCE)**: foundation-model scale (UCE 1280-dim embeddings, multi-GB pretrained weights, PBMC-specific pretrained checkpoints, unclear GPU/VRAM on our cluster; out-of-domain for multi-tissue benchmark) → **NOT included**; record as a later step in TODO.md.
- **Combos** (user decision): both new methods mirror PILOT exactly — `_highres` × hvg{1000,2000,3000} + `_lowres` × hvg2000; default combo = hvg2000_highres. PILOT-GM-VAE: `num_classes = max(2, n_unique_cell_types_in_ct_col)`, `epochs=50` (function default). QOT: `num_components_list=[1]`, `min_samples_for_gmm=0`, `qot_method='cosine'`.
- **Figure placement** (user decision): both methods are added to Figure 2A (main benchmark figure) AND all extended lists (Supp fig 15, presentation, exec-time, funky heatmaps).
- **Testing**: `_debug` benchmark-view h5ad exists on NAS at `/Volumes/Shared/Projects/ECODA_paper/_debug/output/JoaI_2022_35773407_debug_5samples_benchmark_analysis_ECODAprocessed.h5ad` (2500 cells, 5 samples, 10 cell types, obsm `X_pca_benchmark_analysis_hvg{1000,2000,3000}`, `layers["counts"]`, `hvg_rank`) → run the py script directly against it from the Mac (py-cpu env), then HPC-validate on `_debug`.

## Method naming (exact strings; R/notebook depend on them)

- Feather: `{ds}_hvg{n}{res}_qot_dists.feather`, exec-method string `QOT_hvg{n}{res}`
- Feather: `{ds}_hvg{n}{res}_pilotgm_dists.feather`, exec-method string `PILOT-GM-VAE_hvg{n}{res}`
- R bundle keys: `QOT_hvg{n}` (+ `QOT_hvg{n}_lowres`), `PILOT-GM-VAE_hvg{n}` (+ `_lowres`)
- Labels (constants.R / figures): `QOT`, `PILOT-GM-VAE`
- Submit/manifest labels (merge + log filenames): `qot`, `pilotgm`
- Feather layout = plain `DataFrame.to_feather()` with pandas index (sample names) — identical to PILOT; R ingest via `column_to_rownames(ncol)`.

## Task list

### 1. Dependencies (`pixi.toml`)
- Add to `[pypi-dependencies]`: `pilotgm = "==0.1.1"` and `phate = ">=1.0,<2"` (both pure-python; resolve on osx-arm64 and linux-64). Do not touch any existing pins.
- Local (Mac): `pixi install` (plain; no HPC guard needed), then verify: `.pixi/envs/py-cpu/bin/python -c "import pilotgm, phate"`.
- HPC: user refreshes via the guarded entry points only (`sbatch src/utils/bash/setup_env_sbatch.sh` preferred; `src/utils/bash/refresh_env.sh` for a small login-node re-run, inside tmux, no active jobs) and waits ~1h for BeeGFS stale-cache views before submitting arrays (documented HPC convention; worker self-retry covers residual flakes).

### 2. Vendor QOT script + the single bug fix
- Copy `qot_utils_re.py` (PennShenLab/QOT @ `28cd529880c1`) to `src/5_run_benchmark_methods/run_python_sample_embedding_methods/qot_utils_re.py`, prepend a header comment: provenance (repo, commit SHA, date), MIT license note, "modified: one bug fix in Gaussian_Mixture_Representation (see below) — no other changes".
- Fix in `Gaussian_Mixture_Representation`: change the `elif num_samples == 2:` tail into an `else` branch covering `3 <= num_samples < num_components + min_samples_gmm` with the identical single-Gaussian dict (mean of group, `np.zeros((1, d, d))` covariances, weight `[1]`, same `proportion` formula). Keep the `== 1` and `== 2` branches as-is.
- Do NOT touch anything else (including `Extract_Info` index dropping — verified fine).
- Known latent issues left unfixed (document in a comment/plan note, they are on dead paths for us): `sns.barplot(ci=None)` deprecation (only in `plot_hor_vs_vert`, subgroup-DE path), hardcoded `Results_PILOT` paths (`compute_diff_expressions`), `compute_shapley_values`/`trajectory_analysis` unused by us.

### 3. Extend `1.1.1_benchmark_methods_py.py`
- `--method` choices: add `qot`, `pilotgm`.
- Lazy imports (mirror `get_scpoli()` pattern): `import qot_utils_re` inside `run_qot` (avoids phate import for other methods); `import pilotgm` inside `run_pilotgm`. Note: `qot_utils_re.py` sits next to the py script and is importable only because the script's directory is `sys.path[0]` when run as `python <path>/1.1.1_benchmark_methods_py.py` — do NOT convert to a package-relative import.
- Combo rules (mirror PILOT; reuse `run_pilot_for()` logic — rename/alias as `run_extra_for` or add per-method predicates):
  - `qot`: highres × {1000,2000,3000}, lowres × 2000; `run_qot(sub, ct_col, view, n, out_path)`; skip lowres combos when `cell_type_low_res` is None (same guard as pilot/scpoli).
  - `pilotgm`: same combos; `run_pilotgm(sub, ct_col, view, n, out_path)`.
- `is_default_combo`: `qot` → n==2000 && highres; `pilotgm` → same. Defaults-first stable sort unchanged.
- `run_qot(...)`:
  - emb key `X_pca_{view}_hvg{n}`; pass ndarray as-is (script's `Extract_Info` builds its own DataFrame).
  - NaN fill in ct column → `"Unknown"` (scPoli pattern; script itself filters `Cell_type != 'Unknown'`, NaN would otherwise leak into groupby keys).
  - Distinct temp obs column for `progession`: `adata.obs["_bench_prog"] = adata.obs["Sample"]`, then `Run_QOT(adata, gene_matrix=emb_key, type_cell=ct_col, id_col="Sample", progession="_bench_prog", dataset_type="rna", num_components_list=[1], random_state=2, min_samples_for_gmm=0, qot_method="cosine", normalized_set=False)`. Do NOT pass `progession="Sample"` (duplicate-key rename bug, see context). The `status` column is only carried into `Datafame_for_use` and dropped in the GMM step — never used in the distance path.
  - Build `df_dists = pd.DataFrame(adata.uns["QOT_Distance"], index=samples, columns=samples)` with `samples = adata.uns["Datafame_for_use"]["sampleID"].unique()`; `df_dists.to_feather(out_path)`.
  - Exec-time + mem logging identical to `run_pilot`.
- `run_pilotgm(...)`:
  - emb key `X_pca_{view}_hvg{n}`; wrap in named-columns DataFrame and store back in `adata.obsm[emb_key]` (required: `gmmvae_wasserstein_distance` → `extract_data_anno_scRNA_from_h5ad` accesses `.columns`, same as PILOT).
  - NaN fill in ct column → `"Unknown"`; `num_classes = max(2, int(adata.obs[ct_col].nunique()))`.
  - Distinct temp obs column for `status`: `adata.obs["_bench_status"] = adata.obs["Sample"]`, pass `status="_bench_status"` (duplicate-key rename bug, see context).
  - **cwd protection**: `train_gmvae` hardcodes `./trained_models/<dataset_name>/` — `os.chdir()` into a node-local tempdir (`tempfile.mkdtemp(prefix="pilotgm_")`) around the pilotgm block, `os.chdir()` back in a `finally` (worker cwd is `${PROJECT_ROOT}` on HPC; prevents repo pollution). NOT `${output_dir}`: that lives under `${HPC_SCRATCH_DIR}/benchmark/`, which the submit tail rsyncs wholesale to the NAS — saved weights would pollute the NAS sync. Weights are ephemeral by design anyway (`load_weights=False`; retries re-train from scratch), so the tempdir is sufficient; keep `save_model=True` (small files, harmless).
  - `train_gmvae(adata, dataset_name=ds_name, pca_key=emb_key, labels_column=None, epochs=50, num_classes=num_classes, cuda=use_cuda, gpuID=0, load_weights=False, save_model=True, seed=1)` with `use_cuda = args.device == "cuda" or (args.device == "auto" and torch.cuda.is_available())` — plain `device == "cuda"` would silently run CPU under the default `--device auto` on GPU nodes.
  - `gmmvae_wasserstein_distance(adata, emb_matrix=emb_key, clusters_col="component_assignment", sample_col="Sample", status="_bench_status", wass_dis=True, covariance_type="full", num_components=None)` — `num_components` is recomputed inside from `component_assignment` when `wass_dis=True`; keep defaults for the rest.
  - Save `adata.uns["EMD_df"]` (`sampleID`-indexed, already the required layout) to feather; exec-time + mem logging as above.

### 4. HPC plumbing
- `1_submit_hpc_array.sh`: `METHODS=(mrvi scpoli pilot qot pilotgm)`; case mapping: `qot` → CPU partition (same branch as `pilot`); `pilotgm` → GPU partition (same branch as `mrvi|scpoli`); update the `*)` error message. Everything else (manifest, throttles, `--sync-only`, merge labels via `benchmark_merge_sync_cleanup "${METHODS[@]}"`, per-task log naming `execution_times_<method>_<ds>.feather`) works unchanged.
- `1.1_run_worker.sh`: no change (METHOD env passthrough, worker_retry, feather deletion on requeue all generic). Verify nothing hardcodes the three old method names.

### 5. R ingest + notebook
- `benchmark_methods_r.R`: add `process_qot_fig()` and `process_pilotgm_fig()` — copies of `process_pilot_fig()` (read_feather → `column_to_rownames(ncol)` → `standardize_sample_names` → `create_result_bundle(..., dist_mat = as.dist(feat_mat))`). No `load_all_functions`/worker-loader changes needed (these files are already in both loaders).
- `benchmark_pipeline.R::run_benchmark_analysis`: add the qot + pilotgm ingest blocks mirroring the PILOT block exactly (per-HVG highres file + lowres@2000 block), keys `QOT_hvg{i}`/`QOT_hvg{i}_lowres`, `PILOT-GM-VAE_hvg{i}`/`PILOT-GM-VAE_hvg{i}_lowres`.
- `src/utils/constants.R`: `method_label_map_main` += `"PILOT-GM-VAE_hvg2000" = "PILOT-GM-VAE"`, `"QOT_hvg2000" = "QOT"`.
- `notebooks/benchmark_analysis.rmd` (every list where `PILOT_hvg2000` appears):
  - `default_methods` (line ~87): add `"PILOT-GM-VAE_hvg2000_HR\\b"`, `"QOT_hvg2000_HR\\b"` — the regex matches the RECODED screening-figure names (existing entry is `PILOT_hvg2000_HR\b`); a bare `QOT_hvg2000\b` would NOT match `QOT_hvg2000_HR` (`_` is a word char, so no word boundary sits between the trailing `0` and `_`).
  - Supp fig 2 screening chunk (~line 1272, user-confirmed: INCLUDE the new methods): add `grepl("^QOT_", method) ~ "QOT"` and `grepl("^PILOT-GM-VAE_", method) ~ "PILOT-GM-VAE"` branches to the `method_plot` case_when (the existing `^PILOT_` branch does NOT match the dashed name — the char after `PILOT` is `-`, not `_`) and extend the `factor(levels = ...)` with `"QOT"`, `"PILOT-GM-VAE"`. Uncovered rows turn NA in the facet column and are silently dropped from the figure. Also add the analogous `str_replace("QOT_hvg(\\d+)", "QOT_hvg\\1_HR")` / `"PILOT-GM-VAE_hvg(\\d+)"` rules next to the existing PILOT one (line ~1297) so the x-labels + default-highlight work. Lowres rows render as `QOT_hvg2000_HR_LR` etc. — same cosmetic quirk as the existing `PILOT_hvg2000_HR_LR`, keep for consistency.
  - All method vectors: exhaustive grep-driven rule — every `filter_methods <- c(...)` definition (lines 338, 438, 1380, 1447, 1956, 2115, 2199, 2367, 2464, 2524; 1026 is ECODA-only — skip) and every `"PILOT_hvg2000"` occurrence (25 sites: vectors at 358, 450, 1398, 1462, 1968, 2123, 2219, 2380, 2475, 2532 + the PCA-plot `methods <- list(...)` recodes at 658, 765, 835 + recode-map entries) — add both methods next to `PILOT_hvg2000` (Figure 2A funky heatmap at 2118/2123, Supp fig 15 at 2375/2380, presentation at 2470/2475 + 2524/2532, exec-time plot at 2199/2219).
  - All ggplot recode maps: `"PILOT_hvg2000" = "PILOT"` blocks → also `"PILOT-GM-VAE_hvg2000" = "PILOT-GM-VAE"`, `"QOT_hvg2000" = "QOT"` (main + extended figures).
  - Line ~1297 method-key normalization regex: add analogues for `PILOT-GM-VAE_hvg(\d+)` and `QOT_hvg(\d+)` (the existing `PILOT_hvg(\d+)` pattern does NOT match the dashed name — that's correct, add explicit rules).
  - `normalize_exec_keys` (~line 1927): add `gsub("^(PILOT-GM-VAE_hvg[0-9]+)_highres$", "\\1", method)` and same for `QOT_hvg[0-9]+` (feather rows `_highres` ↔ bundle keys without suffix).
  - Check `method_label_map_main`-based figures consume the new keys through the shared recodes (follow how PILOT flows).
- `batch_effect_analysis.rmd`: out of scope (Phase 4); PILOT-GM-VAE on `X_pca_harmony` for batch-effect views stays in TODO.md Phase 4.

### 6. Testing (user-requested: implementation + validation with `_debug`, direct from NAS)
Local (Mac, py-cpu env; h5ad read from NAS read-only, outputs to a scratch dir):
1. After `pixi install`: import check `pilotgm`, `phate`, vendored `qot_utils_re`.
2. **QOT bug-fix unit test**: tiny synthetic AnnData (2 samples × 2 cell types, one (sample,ct) group with exactly 4 cells, `min_samples_for_gmm` high enough to trigger the old fall-through) → assert the group is no longer dropped (all sample pairs have non-zero distance; `params` contains the (source, ct) key).
3. **End-to-end on `_debug`** (from repo root):
   - `.pixi/envs/py-cpu/bin/python src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py --config_path datasets.json --ds_name _debug --view benchmark_analysis --method qot --input_dir /Volumes/Shared/Projects/ECODA_paper/_debug/output --output_dir <scratch> --hvg 2000`
   - same with `--method pilotgm` (2500 cells × 50 epochs on Mac CPU is a few minutes; keep `--hvg 2000` only).
    - Assert: `{ds}_hvg2000_highres_qot_dists.feather` / `..._pilotgm_dists.feather` are 5×5 with sample IDs as index; `execution_times_qot__debug.feather` has rows `QOT_hvg2000_highres`, `PILOT-GM-VAE_hvg2000_highres` — plus the `_lowres` rows (with `--hvg 2000`, each method runs 2 combos: highres + lowres). Note `_debug` has `cell_type_low_res == cell_type_high_res == "cell.type"` (flattened from `columns` by `read_datasets_json`), so the lowres combo runs on identical input — fine for a smoke test, but don't infer highres/lowres differences from `_debug` numbers.
4. **R ingest smoke test**: Rscript (pixi R) loading `src/utils/load_all_functions.R` + `benchmark_methods_r.R`, call `process_qot_fig`/`process_pilotgm_fig` on the feathers with a minimal labels df → bundle has `feat_mat` 5×5 + `dist_mat`, no error.
5. Optional: render `benchmark_analysis.rmd` on `_debug` if the user wants the figure path exercised before HPC.

HPC validation (user runs, per AGENTS.md conventions — full copy&paste commands in the final summary):
- After env refresh + BeeGFS wait: `source src/slurm_config.sh && cd "${PROJECT_ROOT}" && ./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods qot,pilotgm` (qot → CPU array, pilotgm → GPU array; blocking tail syncs `benchmark/embeddings` to NAS).
- Then full-dataset rollout: same submitter without `--ds_name` (all benchmark datasets), or per-method for parallelism across terminals (`--methods qot` and `--methods pilotgm` can run concurrently).

### 7. Docs + housekeeping
- `TODO.md`: mark 3.2 PILOT-GM-VAE + QOT DONE (both benchmark-view; batch-effect-view Harmony input stays Phase 4); replace the PULSAR bullet with a "later step" note (foundation-model/UCE scale, PBMC-specific weights — feasibility check needed; candidate for a dedicated plan).
- `.gitignore`: add `trained_models/` (belt-and-braces; the script's hardcoded relative model dir).
- README.md / docs/ARCHITECTURE.md: add the two methods to the python-methods pipeline description + feather-naming/table sections (keep brief).

### 8. Risks / notes
- `train_gmvae` trains on the whole dataset per (hvg, res) combo — GPU array wall time will be dominated by PILOT-GM-VAE (expect hours on large datasets; `--time 12:00:00` in the worker header may need raising for the largest datasets — verify on `_debug` first and adjust if needed; keep within partition limits).
- `gmmvae_wasserstein_distance` uses `Parallel(n_jobs=-1)` + `ot.emd2(..., numThreads=50)` — fine with the pinned CPU counts; `compute_emd` is called on the GPU worker (pilotgm array) — CPU threads do the OT while torch trains.
- Requeue path: `1.1_run_worker.sh` deletes `*${DS_NAME}*.feather` on transient failure → pilotgm re-trains from scratch on retry (idempotent; expected cost).
- `phate` adds a small dependency closure (scprep, graphtools, igraph + new transitive `kneed`/`tasklogger`) — already largely present (leidenalg/igraph/scikit-network/elpigraph/pydiffmap in lock); the resolver adds the rest automatically, all pure-python.
- Do not run HPC pipeline scripts as part of implementation validation except the explicitly requested `_debug` test; all sbatch-submitted work happens on HPC by the user.
