# Plan: Adapt `notebooks/benchmark_analysis.rmd` to the new pipeline

## Goal

Make `notebooks/benchmark_analysis.rmd` run end-to-end against the new pipeline
outputs (per-view h5ads on NAS, unified `execution_times.feather`, HPC result
bundles + pseudobulks), add the RAM-usage plot, fix plot correctness issues,
remove broken/orphaned code, refactor duplicated blocks — **without committing
or pushing anything** (user reviews first; see Commit policy).

The notebook runs locally on a macOS machine (64 GB RAM): data loading must NOT
fully load h5ads (backed reads only: obs/metadata, one obsm array).

## Decisions (user-confirmed)

1. **Loading**: adopt new-pipeline naming and reuse pipeline helpers. Tap h5ads
   (backed): obs only for composition methods; single obsm read for Avg_PCA.
   The counts → DESeq2-pseudobulk step (needed only for ECODA_deconv) is
   already on HPC — the notebook reuses the synced NAS pseudobulk variant
   `hvg2000` (verified identical computation: `get_pb_deseq2(hvg = NULL,
   n_hvg = 2000, black_list = "none")`, same call as the legacy notebook's
   `pb_norm`). → No counts access in the notebook at all.
2. **RAM plot**: standalone, info-only (no numbers/regression stats), placed
   directly below the exec-time figure section; Python-method rows only
   (MrVI/scPoli/PILOT; R rows have `mem_GB = NA`); log-log axes.
   R-method RAM (sacct MaxRSS backfill) → TODO.md "Ideas for later".
3. **Paper figures are sacred**: every plot saved with a filename starting
   `Figure` or `Supp_fig` is kept and FIXED if broken — never removed
   (add this rule to AGENTS.md). This includes Figure 4 CE/DF (they work —
   `ECODA_authors_HR_top_varexp0` yields the top-2 CTs, verified via
   `get_hvcs` min-2 fallback).
4. **ECODA+PB legacy chunks**: KEPT but commented out / flagged as legacy
   (internal testing, decided NOT to show in publication). AGENTS.md note.
5. **Funkyheatmap**: refactor the 3 near-identical blocks into one
   parameterized function. Metrics become a notebook-level parameter; default
   set = ANOSIM (`anosim_score`), Modularity (`mod_knn3_score`), ARI
   (`cluster_score`); any combination selectable (e.g. only silhouette +
   `mod_knn9_score`).
6. **GongSharma naming**: internal key stays `Gongsharma_cmv_young_males`;
   displayed as "GongSharma" in figures (current line-144 rename kept); the
   exec-time join key is normalized to match.
7. **pixi**: add `r-ggthemes` + `r-writexl` (pinned) to pixi.toml (approved);
   local `pixi install` by the implementer; HPC env refresh is a user to-do.
8. **Zero-imp method naming**: the underscore-less names (`counts_all0.5`)
   come from the HPC refactor of `run_zeroimp_analysis` (paste0 without
   separator); the legacy notebook referenced `counts_all_0.5` / `counts_all_1`
   (verified in git history). Fix PIPELINE-side: rename to underscore style
   with string-formatted values (`counts_all_0.5`, `counts_all_2/3`,
   `percentage_all_0.1%`, `multLN_0.1%`, ...). Notebook refs then match
   legacy as-is. Requires a zeroimp array re-run (--force) for already-computed
   datasets → user to-do.
9. **`standardize_sample_names()` is NOT needed in the notebook** — verified
   the preprocessing already standardizes the Sample column
   (`1.1.1_preprocess.py:285-291`: leading digit → `g` prefix, `-` → `_`, same
   logic as the R function).
10. **`replace_HiTMElayer3_annot`**: comment out the call in the notebook
    (legacy test, keep as documentation); the function stays in
    `seurat_utils.R` untouched.

## Verified findings (checked against current code)

Hard breaks (notebook errors today):
- `datasets[[ds]][["ds_name"]]` → `NULL` (field no longer in
  `read_datasets_json()`); legacy `data/<ds>_ECODAprocessed.rds` files no
  longer exist — lines 98, 3732.
- `data/execution_times.feather` no longer produced locally; unified feather is
  at `${NAS_TARGET_DIR}/benchmark/embeddings/execution_times.feather`
  (schema: dataset, method, time_secs, mem_GB; `mem_GB = NA_real_` for R rows).
- Zero-imp flattening errors on the new nested `multLN*`/`multRepl*` bundles
  (`Can't combine ..1$value <double> and ..2$value <list>`, verified) —
  lines 3378-3391.
- `RNA_snn_res.*` columns gone → `ECODA_seuratres_*` fails; new columns:
  `leiden_res_{r}_benchmark_analysis_hvg2000` (r = 0.1, 0.4, 2, 5, 20, 50) —
  `benchmark_pipeline.R:551-557`.
- `theme_classic2` (ggthemes) and `writexl` not in pixi env (verified absent
  from pixi.toml and installed R libraries).
- Chunk-order bugs: `default_methods_regex` used before defined (line 393 vs
  1116-1144); `feat_mat_recalc` used before defined (line 1605 vs 1627).
- `_debug` dataset is returned by `read_datasets_json(view =
  "benchmark_analysis")` → must be skipped (`_*` convention).
- PILOT exec-time naming mismatch: feather rows `PILOT_hvg{n}_highres` vs
  bundle names `PILOT_hvg{n}` → `merge_exec_times()` silently misses all PILOT
  rows. GongSharma key mismatch breaks the join for that dataset.
- ggpubr `stat_compare_means(ref.group = "All cell types")` errors
  ("missing value where TRUE/FALSE needed", verified) — HVC chunks lines
  1498, 1543, 1566: no such level; correct ref = `"100"`.

Orphaned/broken (to remove per decision 3 — none of these save a
Figure/Supp_fig output):
- CLR-of-two-parts demo chunks (lines 1600-1620) — undefined
  `feat_mat_recalc`, dataset-specific cell types; unsaved exploratory output.
  (Figure 4 CE/DF at 1622-1708 is a separate, working section — KEEP.)
- Dead legacy dataset renames (lines 183-184, 1802) and commented-out blocks
  (lines 755-757, 2180-2184, 2609-2652, 2933-2976, 3288-3297,
  benchmark_pipeline.R:658-671).
- Empty `mod_score` zeroimp boxplot (lines 3628-3651; `mod_score` rows are
  dropped at line 3394 → always-empty PDF).

Kept but commented/guarded (decisions 3, 4):
- ECODA+PB sections (lines 1229-1431): `ECODA_PB_combo_*` is disabled in
  `run_benchmark_analysis`; comment the 2 plot chunks + distance-distribution
  demo as legacy (or guard with a "legacy, methods not in pipeline" message);
  TODO.md "Ideas for later" entry.

Correctness fixes (figures preserved):
- funkyheatmap blocks include 8 metric columns (3 raw-name leftovers; dplyr
  1.2 `recode` keeps unmatched values — verified) → restrict to the
  parameterized default metric set (decision 5).
- `Supp_fig_18_Kfoury_PCA_authors_vs_auto_Kfoury.pdf` duplicated-ds filename
  (lines 478, 514) — cosmetic; flag at review.
- Supp fig 3-13 header vs actual idx range (3..12) — cosmetic; flag at review.
- Supp fig 1 depends on `clean_data` from a previous chunk → make
  self-contained.
- Header "Inputs" text (lines 22-23) names `Preprocess_datasets.Rmd` /
  `Process_data.ipynb` → update to new pipeline.

## Implementation tasks (ordered)

### 0. pixi env
- Add `r-ggthemes` + `r-writexl` (pinned) to `pixi.toml`; run local
  `pixi install`. (HPC refresh: user to-do.)

### 1. Shared helpers (`src/utils/`)
- `seurat_utils.R`:
  - `get_ct_comp_df(obs, sample_col, ct_col)` — data.frame variant of
    `get_ct_comp_df_seurat` (same table → output shape).
  - `rename_leiden_cols(obs, view, resolutions)` — map
    `leiden_res_{r}_{view}_hvg2000` → `RNA_snn_res.{r}`.
- `datasets_io.R`: `get_view_h5ad_path(datasets, ds, view)` — wraps the
  `views[[view]]$output_file` lookup (mirrors `get_h5ad_path` in
  `benchmark_hpc_utils.R`).
- DRY: `run_ct_comps_analysis_worker` (`benchmark_hpc_utils.R:297`) reuses
  `get_ct_comp_df` instead of its inline table().

### 2. `benchmark_pipeline.R`
- `run_zeroimp_analysis`: rename method keys to underscore style with
  string-formatted values — loop `i` over strings `c("0.5", "2/3", "1", "5",
  "10", "20", "50", "100", "200")` (`as.numeric(i)` for computation) →
  `counts_zeros_0.5`, `counts_all_0.5`, `counts_all_2/3`, ...; same for
  `percentage_all_` / `percentage_zeros_` (→ `percentage_all_0.1%`),
  `multLN_` / `multRepl_` (→ `multLN_0.1%`). (Fixes the current
  `counts_all0.666666666666667` artifacts too.)
- `run_benchmark_analysis`: add an obs/adata input path (Seurat path kept as
  fallback): labels/metadata from `obs`; composition via `get_ct_comp_df`;
  `ECODA_seuratres_*` via `rename_leiden_cols(obs, view)` (new `view`
  argument, default `"benchmark_analysis"`); `pb_norm` passed in (from NAS
  pseudobulk `hvg2000`, decision 1) — drop the internal
  `get_pb_deseq2(seurat, ...)` call on the new path; Avg_PCA from a passed
  PCA-embedding matrix (obsm read, per-sample means).
- `process_coda_fig` / `process_gloprop_fig` / `process_avg_pca_embedding_fig`:
  add obs/embedding input variants (Seurat paths unchanged).

### 3. Notebook setup chunk (lines 30-55)
- Paths: `path_nas_project <- /Volumes/Shared/Projects/ECODA_paper`;
  `path_nas_benchmark` → `.../benchmark/results` (kept);
  `path_nas_embeddings` → `.../benchmark/embeddings`; pseudobulk =
  `.../benchmark/pseudobulks/<ds>_pseudobulk_hvg2000.rds` (list(pb,
  time_secs) → `$pb`); h5ad = `file.path(path_nas_project, ds, "output",
  <view output_file_name>)`; keep local `path_data` for caches/xlsx;
  `path_plots` unchanged.
- `datasets[[ds]][["ds_name"]]` usages → `ds` (the key).
- Keep the GongSharma display rename (line 144); drop legacy renames
  (183-184).
- New parameter: `benchmark_metrics <- c("anosim_score", "mod_knn3_score",
  "cluster_score")`, validated against `names(score_label_map)`, with a
  comment showing how to change it (e.g. add `"sil_score"`, `"mod_knn9_score"`).
- Move `default_methods` / `default_methods_regex` definitions up here (fixes
  line-393 ordering bug); remove the later duplicate (line 1132-1144).
- Extend fail-fast guard: `benchmark/results`, `benchmark/embeddings`,
  `benchmark/pseudobulks`, h5ad dir.
- Update "Inputs" description (lines 20-28).

### 4. Notebook dataset loop (lines 75-125) + later re-reads
- Skip `_*` keys.
- Open h5ad backed (reticulate `anndata read_h5ad(backed = "r")`, as in the R
  workers); read `obs` only; no `standardize_sample_names` (already done
  upstream, decision 9); `replace_HiTMElayer3_annot` call commented out
  (decision 10); set `misc` fields; read obsm
  `X_pca_benchmark_analysis_hvg2000` once; load NAS pseudobulk `hvg2000`;
  call the extended `run_benchmark_analysis(..., obs = obs, adata = adata,
  pca_emb = ..., pb_norm = ...)`; close/GC.
- Drop the Stephenson `subset(Site == "Ncl")` (view already applies it).
- Collect `n_cells <- nrow(obs)` per dataset during the loop (replaces the
  full re-load at lines 1850-1870).
- Supp table 1 (lines 3724-3738): compute from obs (n samples, n cells,
  cells/sample mean±sd, n cell types).

### 5. Unified execution times (lines 1766-1982)
- Replace the `r_exec_times` + `p_exec_times` chunks with one source:
  1. NAS merged feather (`path_nas_embeddings/execution_times.feather`) —
     HPC rows (R heavy + Python + trans/zeroimp + prepare_pseudobulk_*).
  2. Bundle-derived rows for local-only methods missing from the feather
     (ECODA_* variants, GloProp, Freq_highres, Avg_PCA_embedding,
     ECODA_deconv, ECODA_seuratres_*, ECODA_HiTME_*, ECODA_scATOMIC_HR):
     union by (dataset, method); `mem_GB = NA`.
- Normalize: dataset `Gongsharma_cmv_young_males` → `GongSharma`; method
  `PILOT_hvg{n}_highres` → `PILOT_hvg{n}`.
- Drop `n_samples`/`n_features` from the table (unused there; already in
  `df_results`); keep `n_cells` (task 4).
- Cache to `data/exec_times.rds`; delete stale caches on schema change.
- Supp fig 14 kept byte-identical (paper figure).

### 6. RAM plot (new, directly below Supp fig 14)
- Info-only chunk: exec_times rows with non-NA `mem_GB` (MrVI/scPoli/PILOT),
  same method recoding as Supp fig 14; `aes(n_cells, mem_GB, color = method)`,
  `scale_x_log10` + `scale_y_log10` + `annotation_logticks` both axes, points
  only (shape 4, alpha 0.5), method colors, `theme_minimal`, legend right.
  Save e.g. `Supp_fig_14B_benchmark_mem_GB.pdf`.

### 7. Zero-imp section (lines 3375-3679)
- Fix flattening: handle one-level-deeper `multLN_*`/`multRepl_*` entries
  (rows for scalar scores only; drop list-valued rows with a message).
- Notebook refs already use legacy names (`counts_all_0.5`, `counts_all_1`)
  — they match the renamed pipeline (task 2); no notebook-side renames.
- Remove the empty `mod_score` boxplot block (lines 3628-3651).

### 8. funkyheatmap refactor + metric parameter
- Extract `build_funky_heatmap(...)` (place in `src/utils/plotting.R`):
  score filter (`score_set = benchmark_metrics`) → `score_label_map` recode →
  dataset display recode → ranks → min-max scale → aggregates → column_info
  (metrics from `score_set` only) → palettes/legends → plot + ggsave.
- Rewire Figure 2A, Supp fig 15, "For presentation" chunks; outputs
  byte-identical except the metric-column restriction (decision 5).
- Supp fig 1: build `clean_data` locally (df_results + exec_times, same
  filter/recode) — no cross-chunk dependency.

### 9. Remaining correctness fixes
- HVC chunks: `ref.group = "All cell types"` → `"100"` (lines 1498, 1543,
  1566).
- Figure 4 CE/DF: verify unchanged output (top_varexp0 = top-2 CTs ✓).
- Remove approved orphaned code (task 3 list); comment out ECODA+PB sections
  with a legacy note (decision 4).
- Fix duplicate-ds filename in Supp fig 18 ggsave calls (flag for review).

### 10. AGENTS.md additions
- Rule: plots saved with filenames starting `Figure`/`Supp_fig` are paper
  figures — keep and fix, never remove.
- ECODA_PB legacy: `ECODA_PB_combo_*` kept (commented) but NOT shown /
  implemented — internal testing only, not in publication figures.
- Benchmark figure hierarchy: main benchmark figure = results across all
  datasets with only the default/main parameter setting per method; extended
  figure = additional non-standard methods; parameter screening = main methods
  across parameter ranges (robustness check that default parameters are
  faithful — no cherry-picking).

### 11. TODO.md + docs
- TODO.md: mark 3.3 notebook adaptation done; "Ideas for later":
  (a) R-method RAM via sacct MaxRSS backfill into `execution_times.feather`;
  (b) ECODA+Pseudobulk combos (legacy, not shown);
  (c) optional future HPC worker dumping per-dataset `obs.rds` + PCA-mean
  matrix if the notebook should ever avoid h5ad reads entirely (currently not
  needed — obs/obsm reads are light).
- ARCHITECTURE.md: small update if it documents the old exec-times reading.

## Validation (no HPC/pipeline runs; AGENTS.md)
- `Rscript` parse-checks: extract all R chunks from the Rmd and `parse()`.
- Unit-test new helpers on synthetic data.frames: `get_ct_comp_df`,
  `rename_leiden_cols`, zeroimp flattening (nested multLN), `build_funky_heatmap`.
- Optional: render the notebook restricted to `_debug` (temporary `ds_filter`
  param in setup) if the NAS is mounted; final full render is done by the user.

## Commit policy (explicitly overrides AGENTS.md task-completion workflow)
- The implementing agent must NOT commit, push, or archive anything.
- After implementation, the user reviews the changes. Only after explicit user
  confirmation: stage the selected subset, move this plan to
  `.kilo/plans/archive/`, commit (including the archived plan), and push.

## User to-dos after implementation
1. Locally: `pixi install` (pulls r-ggthemes + r-writexl) — needed before the
   notebook renders.
2. HPC: refresh the pixi env when no jobs are active
   (`sbatch src/utils/bash/setup_env_sbatch.sh`, or
   `src/utils/bash/refresh_env.sh` inside tmux) — pixi.toml changed.
3. Re-run the zeroimp array with `--force` for any dataset whose zeroimp
   bundles are already on NAS (method names changed, task 2):
   `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh`
   (with the zeroimp analysis flag), then sync.
4. Render `notebooks/benchmark_analysis.rmd` locally (NAS mounted) and check
   the regenerated Figure/Supp_fig PDFs (incl. the new RAM figure).
5. Review the changes, then confirm commit + push + plan archiving.
