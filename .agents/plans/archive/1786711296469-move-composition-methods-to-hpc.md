# Move composition-based benchmark methods to the HPC pipeline (+ fix Adams pseudobulk naming)

## Goal

1. **Phase 0 (bug fix):** unblock `notebooks/benchmark_analysis.rmd`, which crashes on the first dataset (Adams) with
   `Error in if (!all(x == round(x))) ... : missing value where TRUE/FALSE needed` inside `cluster::silhouette`.
2. **Phase 1 (architecture):** move the remaining (fast, composition-based) benchmark methods from the notebook to the HPC
   R benchmark pipeline (`run_r_sample_embedding_methods/`), so the notebook reads **zero h5ad files** and all exec times
   are uniformly HPC-measured (currently Supp fig 14 mixes laptop times for local methods with HPC times).

## Root cause (verified)

- `create_result_bundle()` (`src/5_run_benchmark_methods/benchmark_methods_r.R:493`, NA check :504) does
  `labels <- labels[rownames(feat_mat)]`.
  For Adams, feat_mat rownames come from the HPC pseudobulk bundle, whose sample names are `g022C_a`-style (underscores),
  while the current Adams benchmark h5ad obs has `g022C-a`-style names (hyphens) → 58/100 samples get `NA` labels →
  `calc_sil()` → `silhouette()` → the crash (reproduced: NA labels give exactly this error; NA in dist gives a different
  one). Since the Phase-0.4-style hardening is already in place (`stop()` on NA labels, :504-514), the notebook today
  crashes with that explicit stop message instead of the silhouette error — the fix below is still required either way.
- Why Adams only: `1.1.1_prepare_pseudobulk.R:59` (and the R method worker, `1.1.1_run_benchmark_methods_r.R:111,120,137`)
  re-apply `standardize_sample_names()` (`gsub("-", "_")`, `src/utils/seurat_utils.R:291`) to an obs column already
  standardized upstream (`1.1.1_preprocess.py:285-290`). The Adams h5ad predates the `-`→`_` change in the python
  preprocessor and keeps hyphens. Verified: all other 10 benchmark datasets have pseudobulk sample names ⊆ obs names.
- HPC bundle labels for Adams are also underscore-based, but scores are grouping-based (unaffected) and the notebook
  figures use the local-method labels (`ECODA_authors_HR$labels`, notebook lines 972/1728/1823) — no downstream cross-reference breaks.

## Phase 0 — Adams pseudobulk fix (quick)

1. `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R:59` — remove the
   `standardize_sample_names()` re-application (keep the function for the legacy Seurat path at `benchmark_pipeline.R:435`).
2. `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:111,120,137` — remove the same
   redundant calls (otherwise future Adams R-method runs reintroduce the naming split).
3. Regenerate the Adams bundles on HPC:
   ```bash
   source src/slurm_config.sh && cd "${PROJECT_ROOT}" && \
   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
     --ds_name Adams --methods prepare_pseudobulk
   ```
4. Hardening: **ALREADY DONE — verify only, no code change.** `create_result_bundle` (`benchmark_methods_r.R:504-514`)
   already `stop()`s on NA labels listing the missing-sample count + first 5 names (the plan's original :473 reference is
   stale; the function moved to :493 when the check was added). Confirm the message is present, nothing else to do.
5. Re-knit `notebooks/benchmark_analysis.rmd` → must pass Adams and all datasets. **Save a reference copy of the local
   results** (`cp result_list.rds data/result_list_local_reference.rds`) for the Phase-1 parity check.

## Phase 0.5 — Peak-RAM recording for R benchmark methods

R workers currently log `mem_GB = NA` for every row (`log_exec_row`'s default; the notebook's RAM figure
`Supp_fig_14B_benchmark_mem_GB.pdf` therefore shows python methods only). Mirror the python worker
(`peak_rss_gb()` via `getrusage().ru_maxrss`, `1.1.1_benchmark_methods_py.py:125-133`) so R rows carry real peaks.

### Design

- New helper `peak_rss_gb()` in `benchmark_hpc_utils.R` (next to `log_exec_row`): on Linux read `VmHWM` from
  `/proc/self/status` (kB → GB, the process peak-RSS equivalent of `ru_maxrss`); off-Linux return `NA_real_` (workers
  run Linux; base-R only, no new packages). Same monotonic-cumulative semantics as python: log the value at each
  combo's completion; combos running earlier report the least bloated peak.
- Store the measurement in each per-combo bundle rds as `res[["mem_GB"]]` so the re-emit paths replay the *original*
  value on cache reuse (logging the live cumulative peak on a resume would overstate the combo's RAM).
- `log_exec_row` call sites to pass `mem_gb = peak_rss_gb()` (compute path) / `mem_gb = cached$mem_GB` (re-emit path):
  - `benchmark_pipeline.R`: `run_gloscope_hpc` (compute :964, re-emit :937), `run_mofa_hpc` (:1031/:1015),
    `run_pseudobulk_hpc` (:1107/:1097 and :1139/:1125), `run_scitd_hpc` (:1172/:1157 and :1239/:1223) — all combos.
  - `1.1.1_prepare_pseudobulk.R`: compute :77-78 and re-emit :88-89. Store `mem_GB` per variant in
    `prepare_pseudobulks_hpc` (`benchmark_pipeline.R:889-892`, additive 3rd list element — backward compatible,
    old bundles without it re-emit NA).
  - `run_ct_comps_analysis_worker` (`benchmark_hpc_utils.R:308`, trans/zeroimp single row per dataset).
  - Worker method-RDS re-emit path `1.1.1_run_benchmark_methods_r.R:75-80`: pass `cached[[nm]]$mem_GB`.
  - The new `composition` driver (Phase 1) is written with this from the start.
- Update `log_exec_row`'s docstring (`benchmark_hpc_utils.R:190-193`: no longer "mem_GB = NA_real_ for R rows").
- Notebook (`benchmark_analysis.rmd` exec-times chunk ~:1873-1956): bundle-derived rows (:1913-1920) carry
  `mem_GB` from the bundle when present (`if ("mem_GB" %in% names(bundle)) ...`); update the stale comments at
  :1875-1883 (":1877" "NA for R rows" is no longer true). **Delete `data/exec_times.rds` before the first re-knit**
  (the plan's previously-optional cache-drop becomes mandatory: an existing cache would keep the NA mem rows and
  suppress the new HPC values — the schema check at :1889 won't catch it).
- Datasets whose R bundles predate this change keep `NA` until their next `--force` run (the in-progress Phase 3.1
  rollout covers them); the sacct-`MaxRSS` backfill idea in TODO.md stays as an optional backfill for legacy bundles.

## Phase 1 — Move composition methods to the HPC pipeline

### Design decisions (verified against code)

- One new method name: `composition` → one array task per dataset, one worker pass computing all ~52 combos + the metadata bundle.
  `EPIC` (slowest combo) is minutes; the 12 h / 128 G default worker limits are fine.
- New driver lives in `benchmark_pipeline.R` (HPC DRIVERS section, after `run_scitd_hpc`) — automatically available to
  workers via `load_worker_functions.R` (sources `benchmark_pipeline.R` + `seurat_utils.R` → `rename_leiden_cols` available).
- `composition` consumes the hvg2000 pseudobulk variant (ECODA_deconv) → `1_submit_hpc_array.sh` MUST auto-prepend
  `prepare_pseudobulk` (extend the existing `NEEDS_PREP` condition, lines 149-166).
- Obs-only path in the worker: **no Seurat materialization** (backed h5ad obs + obsm `X_pca_benchmark_analysis_hvg2000` +
  pb variant via `load_pb_variants()` from `benchmark_hpc_utils.R:130`). Labels/metadata straight from obs, **no**
  `standardize_sample_names` (Phase 0 removes it anyway).
- `rename_leiden_cols(obs, view="benchmark_analysis")` must be applied on the worker for the `ECODA_seuratres_*` combos
  (defaults: resolutions `c(0.1, 0.4, 2, 5, 20)` — matches `run_benchmark_analysis`'s `seurat_res`).
- Seed: `set.seed(123)` at driver start (one worker process per dataset → deterministic across HPC re-runs). Note:
  `ECODA_authors_HR_NULL` (shuffle_labels) will NOT be bit-identical to the old notebook run (sequential stream) — it is a
  null control; document this in the notebook comment.
- Checksums + exec-log merge are generic: `benchmark_merge_sync_cleanup` hashes all `*.rds` under `results/` and merges
  per-(label × dataset) exec logs → new bundles/rows are covered with no submit-script changes beyond the method list.
- Python-feather methods (MrVI/PILOT/QOT/PILOT-GM-VAE/scPoli) **stay in the notebook** (feather reads + scoring,
  seconds, no h5ad I/O). Trans/zeroimp already on HPC — unchanged.
- **Drop `result_list.rds` checkpointing (Option A).** The notebook starts fresh each knit:
  `load_hpc_benchmark_results` re-reads all bundles (fast `readRDS`), the feather methods recompute (seconds),
  and stats come from `<ds>_metadata.rds`. This eliminates the stale-in-memory clobber hazard observed in
  practice (a session holding an old list silently overwrote the good file via the per-dataset `saveRDS`
  inside `load_hpc_benchmark_results`) and makes knits fully deterministic. Verified: `result_list.rds` and
  `load_hpc_benchmark_results` are used ONLY by `benchmark_analysis.rmd` (grep across `notebooks/*.rmd`) —
  the batch-effect notebook does not use them, so the change is isolated.

### Tasks

1. **Driver** `run_composition_methods_hpc(...)` in `benchmark_pipeline.R` — follow the `run_gloscope_hpc` pattern exactly
   (per-combo cache file `<ds>_<combo>.rds`, skip-if-exists with exec-time re-emit, `exec_time` appended via `exec_time()`,
   `save_rds_atomic`, `log_exec_row`). Per Phase 0.5, also store `res[["mem_GB"]]` in each cache file and re-emit it
   (`mem_gb = cached$mem_GB`) on cache reuse. Reuse existing functions verbatim:
   - `Avg_PCA_embedding` → `process_avg_pca_embedding_fig(seurat=NULL, labels, pca_emb=…, obs=…)` (needs hvg2000 obsm)
   - `ECODA_deconv` → `process_deconv_fig(t(pb_hvg2000), labels)` (needs pb variant)
   - `ECODA_authors_LR`, `ECODA_authors_HR`, `ECODA_authors_HR_NULL` → `process_coda_fig(obs=…)`
   - `ECODA_authors_HR_top_varexp{0,0.1,…,0.9}` (10) → `process_coda_fig(ECODA_top_varexp_hvct=…)`
   - `ECODA_authors_HR_3most_varcts`, `_2least_varcts`, `_3least_varcts` → `process_coda_fig(ECODA_top_n_hvct=…, var_ct_desc=…)`
   - `ECODA_authors_HR_{2,3,5,10,15}_PCA_dims` → `process_coda_fig(pca_dims=…)`
   - `ECODA_seuratres_{0.1,0.4,2,5,20}` → `process_coda_fig(ct_col="RNA_snn_res.<r>")`
   - `GloProp` → `process_gloprop_fig(metadata, ct_col=ct_high_res, label_col, obs=…)`
   - `Freq_highres` → `process_coda_fig(calc_clr=FALSE)`
   - `ECODA_HiTME_HR_layer2`/`layer3` (+10 varexp each), `ECODA_scATOMIC_HR` → `process_coda_fig(ct_col="layer2"/"layer3"/"scATOMIC_pred")`
     **guarded**: skip combo with a warning when the ct column is absent from obs (annotation produces these columns —
     `src/4_cell_type_annotation/2.1.1_process_chunk.R:103-105` — but availability varies; the old notebook would crash
     on missing columns, so the guard is strictly better).
   - Defaults must mirror `run_benchmark_analysis` (factors_test `c(2,3,5,10,15)`, seurat_res `c(0.1,0.4,2,5,20)`,
     ECODA_top_varexp_hvct `seq(0,0.9,0.1)`).
   - Also emit `<ds>_metadata.rds` = `list(labels = named factor, n_cells, n_samples, cells_per_sample = named int,
     n_cell_types_high_res)` — replaces the notebook's obs reads (stats, Supp table 1, exec-times n_cells).
2. **Worker** `1.1.1_run_benchmark_methods_r.R`: add `composition` to the allowed-method check (line 47); new obs-only
   branch: read backed h5ad obs, fetch hvg2000 obsm, `load_pb_variants(..., variants="hvg2000")`, `rename_leiden_cols(obs)`,
   build labels/metadata from obs; call the driver; `save_rds_atomic` to `<ds>_composition.rds` (and `<ds>_metadata.rds`).
   Method-level skip-if-exists + exec-row re-emit as for the other methods.
3. **Submit script** `1_submit_hpc_array.sh`: add `composition` to `VALID_METHODS` (line 134) and to `NEEDS_PREP` (line 149-166).
   No other changes needed (worker args already include `--pseudobulk_dir`; `1.1_run_worker.sh` else-branch dispatches to
   the right R script; watchdog manifest is built from the method list).
4. **Notebook** `benchmark_analysis.rmd` (loop chunk lines ~119-251):
   - Drop the h5ad-backed read, `pb_norm` read, `pca_emb` read and the `run_benchmark_analysis(...)` call (lines ~139-231).
   - Call `load_hpc_benchmark_results(result_list, ds, path_nas_benchmark,
     methods = c("gloscope","mofa","pseudobulk","scitd","composition"))`.
   - Read `<ds>_metadata.rds` → `labels` for the feather methods + `n_cells_by_ds`/`ds_stats_list` (replace obs-based computation).
   - Keep the fail-fast NAS checks; replace the per-dataset h5ad/pb checks with a `<ds>_metadata.rds` existence check.
   - Keep the python-feather block (MrVI/PILOT/QOT/PILOT-GM-VAE/scPoli, currently inside `run_benchmark_analysis` lines
     680-792) — call `process_*_fig` directly with labels from the metadata bundle.
    - Exec-times chunk (lines 1873-1956): `n_cells_by_ds` source changes to the metadata bundles; HPC-logged
      rows for the new combos (with `mem_GB` per Phase 0.5) flow in automatically via the merged
      `execution_times.feather`; bundle-derived rows (:1913-1920) pick up bundle `mem_GB` when present; delete
      `data/exec_times.rds` before re-knitting.
    - **Option A — remove `result_list.rds` checkpointing:**
      - Replace the result-list chunk (lines ~107-113, `readRDS(result_file)` branch; drop any visibility
        `message()` the user may have added interactively) with a fresh `result_list <- list()`.
      - Remove the `saveRDS(result_list, file = "result_list.rds")` from `load_hpc_benchmark_results`
        (`benchmark_pipeline.R:307`) and update its docstring: bundles are always loaded fresh, no persistence,
        no "entries already present are kept" rerun semantics.
      - Update the loop comment (lines 213-216, "legacy rerun semantics") and the fail-fast chunk comment
        (lines 125-127: the "AND overwrite result_list.rds with the partial list" hazard is gone).
      - Delete the local `data/exec_times.rds` cache (previously "optional", now **mandatory** per Phase 0.5:
        a stale cache pins the old NA `mem_GB` rows and the :1889 schema check does not detect it) —
        regenerate from the merged feather + bundle timings + metadata stats each knit.
    - Update the stale comment at lines 215-216 ("fast composition-based benchmark methods run locally").
   - `run_benchmark_analysis` itself: leave in place, add a deprecation comment (notebook-only caller; no deletion).
5. **Docs**: update AGENTS.md (Stage 3: notebook no longer computes composition methods; R workers log peak RAM),
   `docs/ARCHITECTURE.md` (worker dispatch, new method, metadata bundle, `peak_rss_gb`), TODO.md:
   - Phase 3.2: note the PILOT NaN fix (commit `c547613`, `fill_unknown_ct` in `1.1.1_benchmark_methods_py.py:326`) and
     that the Lee/Zhang PILOT bundles on NAS are stale (all-zero EMD distances) — re-run pending (see follow-ups below);
     update the QOT/PILOT-GM-VAE `_debug` HPC validation status.
   - Remove the "R-method peak RAM: backfill sacct MaxRSS" entry from "Ideas for later" (superseded by Phase 0.5
     in-worker VmHWM logging; reword to "optional sacct-MaxRSS backfill for bundles computed before 2026-08-14").

## Validation

1. `--ds_name _debug --methods composition` on HPC (Joanito 5-sample subset) — smoke test:
   ```bash
   source src/slurm_config.sh && cd "${PROJECT_ROOT}" && \
   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
     --ds_name _debug --methods composition
   ```
2. Full run for all benchmark datasets (`--methods composition`; auto-prepends prepare_pseudobulk). Verify the Adams
   bundle sample names now match obs (hyphens).
3. Parity: compare `<ds>_composition.rds` bundles against `data/result_list_local_reference.rds` (Phase 0 step 5) for the
   deterministic combos — `scores`, `feat_mat`, `dist_mat`, `labels` must be identical (allow `exec_time` + `ECODA_authors_HR_NULL`
   differences). Add a one-off comparison script; do not commit it.
4. Re-knit the notebook: must produce identical Figure 2A/3A/B, Supp figs 2/14/15; exec-time figures now show only HPC times.
5. Verify no dataset silently loses HiTME/scATOMIC combos that it previously had (check `layer2`/`layer3`/`scATOMIC_pred`
   presence per dataset in the validation run; if a dataset lacks them, confirm the old notebook never completed that dataset).
6. **Option A determinism check**: knit twice in a row → identical figures; confirm no `result_list.rds` is written
   anywhere during the knit (the stale-overwrite failure class is gone; a stale file must not reappear).
7. **RAM check (Phase 0.5)**: confirm R rows in the merged `execution_times.feather` have non-NA `mem_GB` for
   datasets run after the change; `Supp_fig_14B_benchmark_mem_GB.pdf` now includes GloScope/MOFA/Pseudobulk/scITD/
   composition + prepare_pseudobulk rows. Verify `data/exec_times.rds` was regenerated (deleted before re-knit).
8. Cleanup after the parity check (step 3) passes: delete `result_list.rds`, `result_list.rds.bak` (repo root) and
   `data/result_list_local_reference.rds` — the reference copy is only needed until the bundles are parity-verified.

## Post-implementation user follow-ups (reminders — the implementing agent must surface these to the user)

1. **Lee/Zhang PILOT re-run** (PILOT NaN bug, fix `c547613` confirmed in the code): HiTME `layer2` annotations on
   Lee/Zhang contain NaN for unclassified cells, which silently produced all-zero PILOT EMD distances; the fix
   (`fill_unknown_ct`) is committed but the NAS bundles are stale. Ask the user to run (parallel-friendly, one
   terminal each; needs `source src/slurm_config.sh && cd "${PROJECT_ROOT}"` first):
   ```bash
   ./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods pilot --force
   ./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Zhang --methods pilot --force
   ```
2. **QOT + PILOT-GM-VAE `_debug` validation** (implementation commit `19ac675`): the `_debug` validation run is
   currently on the HPC; once it passes, the user should roll both methods out to all benchmark datasets
   (`--methods qot,pilotgm` on the python submitter; qot → CPU, pilotgm → GPU arrays). The notebook already ingests
   their feathers (skip-with-message while absent) and TODO.md Phase 3.2 documents the status.
3. Record both items in TODO.md Phase 3.2 (done as part of task 5 above; update once the runs complete).

## Risks / guardrails

- **Paper figures are sacred**: bundle fields (`scores`, `feat_mat`, `dist_mat`, `labels`, `exec_time`) and combo names must
  match the legacy keys exactly — figures read `result_list$bmark[[ds]][[<key>]]` unchanged.
- **No-leakage**: labels remain score-only; no design/batch covariates added anywhere.
- Determinism: per-dataset `set.seed(123)`; HR_NULL differs from the old local run — null control, acceptable.
- Stale-overwrite hazard: eliminated by dropping `result_list.rds` (Option A); no shared use with
  `batch_effect_analysis.rmd` (verified by grep).
- `checksums.md5` regenerated by the sync (hashes all `results/*.rds`) — new bundles covered automatically.
- If a dataset lacks `layer2`/`layer3`/`scATOMIC_pred`, those combos are skipped with a warning (better than the old
  notebook, which would crash) — verify against validation step 5.

## Out of scope

- Moving python-feather methods (MrVI/PILOT/scPoli/QOT/PILOT-GM-VAE) to HPC (kept in the notebook; low value).
- Re-running Adams preprocessing to the current underscore standard (would require re-running annotation + all benchmark
  arrays for Adams; only needed if Adams is ever re-preprocessed — revisit then).
- Trans/zeroimp pipelines (already on HPC).

## Completion workflow (per AGENTS.md)

Archive this plan file to `.kilo/plans/archive/`, stage related changes, commit, push.
