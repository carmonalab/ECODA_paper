# ECODA benchmark additions plan

## Context

Add three narrowly scoped benchmark analyses without rerunning any completed pipeline, dataset, or existing method: Harmony-based unsupervised Leiden composition, ECODA_authors_HR cell-depth subsetting, and MOFAcellulaR low/high-resolution variants. Existing Pipeline 1–5 outputs, checksums, manifests, and gates are treated as immutable. Harmony and cell-subsetting are local derived analyses; MOFAcellulaR is conditional on a successful `_debug` feasibility run before any benchmark-dataset execution.

Dataset selection is the config-derived union, not the existing Stage 5 intersection:

```text
union = entry.use_for_benchmark == true
     OR entry.views.benchmark_analysis exists
```

Exclude names beginning with `_` except the explicit `_debug` probe. Every union member intended to run must have a declared, readable, checksum-valid `benchmark_analysis` H5AD. If a union member lacks that view/artifact, emit a blocked-dataset report and stop before producing partial output; never substitute a batch-effect view and never silently use the current submitter selector at `src/5_run_benchmark_methods/1_submit_hpc_array.sh:718-721`.

## Approach

### 1. Isolate new artifacts and source manifests

1. Add a local derived-analysis runner at `src/5_run_benchmark_methods/run_local_ecoda_derived.R`; no equivalent existing runner handles both fast derived analyses. Give it explicit arguments:
   - `--config_path`
   - `--input_dir`
   - `--output_dir` (required, unique run-owned directory under the benchmark-results area)
   - `--analysis harmony|cell_subsetting`
   - `--scope benchmark_union|debug`
   - `--seeds` only for `cell_subsetting`, defaulting to the fixed sequence `101:120`
2. Make the runner build and write a source manifest before processing. Each row contains dataset, `benchmark_analysis` view, H5AD path, source checksum, and requested analysis. Validate that every row is in the config-derived union and that no batch-effect view is present.
3. Reuse the existing backed-H5AD observation loading path used by the benchmark worker (`load_h5ad_counts_free()` in `src/5_run_benchmark_methods/benchmark_hpc_utils.R`) and the existing atomic/checksum helpers (`save_rds_atomic()`, `artifact_checksum_ok()`). Do not materialize or rewrite input H5AD files.
4. New outputs must be written atomically to the explicit run directory. Existing canonical `<DS>_composition.rds`, `<DS>_mofa.rds`, pseudobulks, `.md5` files, gates, and manifests are never opened for writing.
5. Do not add either analysis to the default Stage 5 method list, shared worker dispatch, or generic submitter. Do not invoke `1_submit_hpc_array.sh` for these analyses.

### 2. Implement Harmony Leiden composition locally

1. For each manifest row, require this exact observation column in the existing benchmark H5AD:
   `leiden_res_2_benchmark_analysis_hvg2000_harmony`.
2. Do not recompute PCA, Harmony, neighbors, or Leiden. Pipeline 3 already creates the Harmony representation and Leiden columns in `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:474-513`.
3. Build sample labels using `Sample` and the configured `entry$label_col`; labels are used only for scoring.
4. Reuse `process_coda_fig()` from `src/5_run_benchmark_methods/benchmark_methods_r.R:23-100` with:
   - `seurat = NULL`
   - `ct_col = "leiden_res_2_benchmark_analysis_hvg2000_harmony"`
   - `obs =` the backed-H5AD observation table
   - default ECODA CLR/imputation behavior
5. Write exactly one standalone bundle per dataset:

```text
<DS>_ECODA_seuratres_2_harmony.rds
<DS>_ECODA_seuratres_2_harmony.rds.md5
```

The bundle must contain `feat_mat`, `labels`, `dist_mat`, `scores`, and `counts`, using the existing `create_result_bundle()` contract. The internal method key is `ECODA_seuratres_2_harmony`; the display label is `ECODA_Leiden_res_2_harmony`.
6. Do not map the Harmony column to `RNA_snn_res.2` and do not add it to the existing composition RDS. The ordinary `ECODA_seuratres_2` path in `src/5_run_benchmark_methods/benchmark_pipeline.R:1650-1677` remains unchanged.

### 3. Implement local ECODA_authors_HR cell subsetting

1. Use the same local runner with `--analysis cell_subsetting`. Read only `Sample`, the configured `cell_type_high_res`, and the configured biological label column from each existing benchmark H5AD.
2. Do not reapply Pipeline 3’s 500-cell sample filter. That filtering already occurred upstream.
3. Use these exact ordered targets:

```text
all cells, 2000, 1000, 500, 400, 300, 200, 150, 100, 50
```

For every sample and target, retain `min(original_cell_count, target)` cells. Samples below a target pass through unchanged. Keep all samples in the subset experiment.
4. For targets below `all cells`, draw without replacement for each sample using exactly the fixed seeds `101:120`. The all-cells baseline is one unmodified replicate (`replicate = 0`, `seed = NA`), not 20 duplicated rows.
5. For every target/seed, call the existing ECODA composition implementation with `ct_col = entry$cell_type_high_res`; use biological labels only for ANOSIM scoring. Do not add this experiment to the benchmark method matrix.
6. Write a tidy, atomic run-owned artifact and a tabular companion:

```text
ECODA_authors_HR_cell_subsetting.rds
ECODA_authors_HR_cell_subsetting.csv
```

Store dataset, target, replicate, seed, sample IDs, effective cells per sample, total cells, and ANOSIM.
7. Add a new final section in `notebooks/benchmark_analysis.rmd`, labelled `Supp fig x` until the next unused publication number is assigned. Follow the existing Figure 4A plotting style at `benchmark_analysis.rmd:1947-2035`:
   - x-axis in the exact target order above;
   - y-axis ANOSIM;
   - points are each dataset’s mean across the 20 subsamples;
   - lines connect points within each dataset;
   - bars are means across dataset-level means;
   - whiskers are standard error across dataset-level means;
   - retain replicate-level values in the artifact for diagnostics.

### 4. Add a strict MOFAcellulaR `_debug` feasibility gate

1. Add `src/5_run_benchmark_methods/run_mofacellular.R` as a dedicated wrapper; the current `mofa` implementation is single-view MOFA2 (`benchmark_methods_r.R:329-357`) and has no MOFAcellulaR equivalent.
2. The wrapper must support `--scope debug|benchmark_union`, but production scope must refuse to run unless a checksum-verified debug-pass record is supplied. The debug scope is the only MOFAcellulaR execution allowed before feasibility is established.
3. Run the debug probe only on `_debug`, with two factors because the five-sample debug cohort cannot fit 15 factors. Test both configured cell-type resolutions. The probe must verify package installation, pseudobulk/view construction, MOFA2 execution, factor extraction, finite features, and non-modification of the source H5AD.
4. Install MOFAcellulaR through the repository’s guarded R environment setup at an exact Git commit resolved during the debug probe. Never install an unpinned mutable branch. Record the tested commit in the debug-pass record and every production artifact.
5. Build cell-type-specific pseudobulk views from existing `layers["counts"]` using `Sample` plus `entry$cell_type_low_res` or `entry$cell_type_high_res`. Biological labels are scoring metadata only and must not be supplied as model covariates.
6. Restrict input genes to the repository’s existing `var["hvg_rank"]` top 2000 genes. Any MOFAcellulaR per-view HVG selection occurs in memory only and must not write to or replace H5AD `var` fields.
7. Use the documented MOFAcellulaR preparation settings:
   - `filt_profiles(..., ncells = 0)`;
   - `filt_gex_byexpr(min.count = 5, min.prop = 0.25)`;
   - `filt_views_bysamples(nsamples = 2)`;
   - `filt_views_bygenes(ngenes = 15)`;
   - `filt_samples_bycov(prop_coverage = 0.9)`;
   - `tmm_trns(scale_factor = 1000000)`;
   - `filt_gex_byhvg(prior_hvg = NULL, var.threshold = 0)`;
   - final view filtering at the documented 15-gene threshold;
   - no biological marker/covariate filtering.
8. MOFAcellulaR is allowed to drop samples. Do not pad missing factor rows or fabricate latent embeddings. Record per-view and final dropped sample IDs. Retain the dataset result whenever the model produces a valid factor matrix with at least two samples; fail that variant/dataset closed if the model produces zero/one samples, missing IDs, duplicate IDs, nonfinite factors, or an invalid bundle. Do not drop an entire dataset merely because some samples were removed.
9. Use MOFA2 with `num_factors = 2` for `_debug` and `num_factors = 15` for production, deterministic seed `42`, fast convergence, and `spikeslab_weights = FALSE` as documented by MOFAcellulaR’s workflow. If production `num_factors >= available_samples`, fail that dataset/variant closed rather than changing the requested factor count.
10. Extract the actual MOFA factor rows, subset labels to those exact row IDs, and call `create_result_bundle()` without adding fabricated rows. Add a dedicated standalone validator for these artifacts that explicitly allows dropped samples only for MOFAcellulaR and reports the dropped IDs. Do not broaden the generic `matrix_artifact_validator.py` method contract or existing R method set.
11. Only after the `_debug` pass, permit benchmark scope and write exactly:

```text
<DS>_MOFAcellulaR_hvg2000_factors15_lowres.rds
<DS>_MOFAcellulaR_hvg2000_factors15_highres.rds
```

Each bundle records annotation column, package commit, factor count, filtering parameters, actual sample IDs, dropped sample IDs, source H5AD checksum, `feat_mat`, `labels`, `dist_mat`, and `scores`.
12. If the debug probe fails, stop without adding MOFAcellulaR to the benchmark pipeline or notebook method lists. If it passes and a full benchmark run is later required, use a target-only durable HPC gate with exactly two method rows per eligible dataset; never route through the broad Stage 5 matrix.

### 5. Load standalone artifacts into the notebook

1. Add a small strict loader near the existing benchmark result-loading block in `notebooks/benchmark_analysis.rmd:203-233`. Its signature is `load_derived_bundle(path, dataset, method, expected_samples = NULL, allow_dropped = FALSE)`.
2. The loader must require a nonempty checksum-valid RDS and validate `feat_mat`, `labels`, `dist_mat`, and `scores`. Harmony requires the full expected sample set; MOFAcellulaR uses actual factor sample IDs and reports dropped IDs.
3. After `load_hpc_benchmark_results()`, insert Harmony bundles under `ECODA_seuratres_2_harmony`. After the MOFA debug/benchmark gate succeeds, insert low/high MOFAcellulaR bundles under their exact method keys. Do not change the existing HPC method loader’s method list.
4. In `src/utils/constants.R:25-34`, add:

```r
"ECODA_seuratres_2_harmony" = "ECODA_Leiden_res_2_harmony"
```

5. Use a dedicated Figure 3A/Supp fig 18 annotation vector, adding only `ECODA_seuratres_2_harmony`. Keep author, HiTME, scATOMIC, and ordinary Leiden methods unchanged. The existing Figure 3A/Supp fig 18 plots are at `benchmark_analysis.rmd:1275-1401`.
6. Leave Figure 2A’s default-method list unchanged. Add both MOFAcellulaR variants only to Supp fig 2. Update the Supp fig 2 grouping predicate at `benchmark_analysis.rmd:1655-1676` to group both `MOFA_` and `MOFAcellulaR_` method keys. Do not select a low/high Figure 2A default yet.
7. Add a final `NOTES.md` entry recording the standalone artifact paths, source column, no-rerun guarantee, 20-seed subsetting contract, MOFAcellulaR debug gate, exact package commit, dropped-sample policy, and deferred Figure 2A default.

## Critical files & anchors

- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`, `process_view()` around lines 474-513 — existing Harmony embedding and Leiden columns; read-only source.
- `src/5_run_benchmark_methods/benchmark_methods_r.R`, `process_coda_fig()` and `create_result_bundle()` around lines 23-100 and 660-684 — reusable composition/result contracts.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh`, selection and method construction around lines 718-866 — explicitly do not add new methods here; its current selector is not the required union.
- `src/utils/constants.R`, annotation map around lines 25-34 — add the Harmony display label only.
- `notebooks/benchmark_analysis.rmd`, result loading and figures around lines 203-233, 1275-1401, 1614-1728, and 1947-2035 — standalone loading and requested plot integrations.

## Verification

No verification is performed during plan approval. After implementation, run only focused local checks until MOFAcellulaR’s debug gate passes.

1. **Selection/source preflight**
   - Working directory: repository root.
   - Environment: repository Pixi default environment; explicit input and output directories.
   - Use the local derived runner’s validation mode or an equivalent manifest-only path.
   - Expected output: union membership is listed, every runnable row has a valid `benchmark_analysis` H5AD, batch-effect views are rejected, and no existing result or gate path appears in the write set.

2. **Harmony smoke test**
   - Run the local runner on `_debug` only if `_debug` declares a valid `benchmark_analysis` view; otherwise stop with the missing-prerequisite report.
   - Expected output: `<debug>_ECODA_seuratres_2_harmony.rds` with the exact standalone key/contract, finite scores, and complete sample IDs.
   - Confirm the input H5AD checksum and modification time are unchanged.

3. **Cell-subsetting smoke test**
   - Run the local runner on `_debug` with seeds `101:120`.
   - Expected output: all ten targets, one baseline row, 20 rows per non-baseline target/dataset, all sample IDs retained, and effective counts equal to `min(original, target)`.
   - Re-run into a second temporary output directory and compare the tidy output byte-for-byte or by digest.

4. **MOFAcellulaR debug gate**
   - Run only `_debug` with two factors and both low/high annotations.
   - Expected output: valid low/high factor bundles, exact tested package commit, actual sample IDs plus any dropped IDs, finite factor values, at least two retained samples, and unchanged H5AD checksum/mtime.
   - A package/API error, empty view, zero/one retained samples, invalid factor matrix, or source modification is a hard failure; do not proceed to benchmark integration.

5. **Production MOFAcellulaR scope check**
   - Only after the debug-pass record exists, validate the benchmark-union manifest before any full-cohort launch.
   - Expected row count: exactly two MOFAcellulaR variants per runnable union dataset.
   - Any existing method, batch view, unrelated dataset, or extra row aborts the run before worker submission. If full-cohort execution is required, submit through `durable-hpc-gate-ecoda` only.

6. **Notebook contract check**
   - Execute only the affected `benchmark_analysis.rmd` chunks using the repository’s established non-knitr workflow.
   - Expected outputs: Figure 3A and Supp fig 18 include `ECODA_Leiden_res_2_harmony`; the new `Supp fig x` contains the cell-subsetting ANOSIM plot; Supp fig 2 includes both MOFAcellulaR variants after the debug/benchmark gate; Figure 2A remains unchanged.

## Assumptions & contingencies

- Existing completed Pipeline 1–5 artifacts remain valid and are not recomputed. If an input H5AD is missing or invalid, stop the affected new analysis instead of rerunning preprocessing.
- The config-derived union is authoritative for selection. A union member without a runnable `benchmark_analysis` view is reported and blocked; no other view is substituted.
- MOFAcellulaR sample dropping is permitted and is recorded explicitly. No latent rows are fabricated and no sample-retention exception is hidden in a generic validator.
- The `_debug` MOFA probe uses two factors solely because the debug cohort is too small for 15 factors. Production remains fixed at 15 factors.
- Low/high MOFAcellulaR variants remain separate artifacts. Figure 2A’s default is selected only after their benchmark results are reviewed.

## Implementation checkpoint — 2026-09-10

### Completed in the working tree

- Added the root `AGENTS.md` override documenting that Pipeline 1–5
  benchmark completion is authoritative, stale gates do not trigger reruns,
  local Harmony/res50/subsetting work must not launch existing HPC pipelines,
  and full H5ADs must remain on HPC scratch/NAS.
- Added explicit standalone selectors to
  `run_local_ecoda_derived.R`: `res50`, `harmony`, and `cell_subsetting`.
  Res50 and Harmony no longer emit each other’s outputs.
- Added run-owned source/output/status manifests, checksums, owner markers,
  path containment, per-sample label consistency, missing-sentinel handling,
  exact target/seed contracts, and bounded one-dataset-at-a-time processing.
- Added `h5ad_obs_free.py` so cell-subsetting metadata access does not
  materialize expression/count values or PCA/Harmony embeddings.
- Added `extract_derived_composition.R` for a one-time HPC metadata read per
  dataset and exact 181-key cell-subsetting composition cache generation.
- Added snapshot-backed local consumption to
  `run_local_ecoda_derived.R`; it consumes compact counts/labels without
  reopening H5ADs and preserves remote H5AD provenance.
- Added strict snapshot-aware notebook loading, independent res50/Harmony
  roots, the Figure 3A/X3B res50/Harmony integration, Supp fig 2
  `MOFAcellulaR_` grouping, the `Supp_fig_X` cell-depth plot, and notes.

### Verification evidence

- Pixi R parsing passes for the standalone runner and snapshot extractor.
- Pixi Python compilation passes for `h5ad_obs_free.py`.
- `knitr::purl()` plus R parsing passes for the modified benchmark notebook.
- Pure helper/manifest/selector/target contract smoke passes in a temporary
  R script.
- The local subset mirrors were tested only as diagnostics and correctly
  failed the full benchmark H5AD contract; they were not accepted as sources.
- A read-only remote scratch inspection confirmed all 11 authoritative
  `benchmark_analysis` H5ADs and adjacent checksums exist, with all required
  source obs columns and persisted schema nodes.
- Temporary local H5AD copies were limited to the small `_debug` diagnostic;
  no full benchmark H5AD was copied to the workstation.

### Remote standalone run findings

- One explicitly scoped standalone temporary job, `4398127`, ran only the
  new derived wrapper against `/home/users/h/halterc/scratch/ECODA_paper` and
  wrote under `_user_temp/derived_local_20260910165057`; it did not invoke
  Pipeline 1–5.
- Res50 and Harmony completed with 11 datasets and exact standalone
  artifacts/manifests.
- Cell-subsetting initially blocked four datasets because missing author
  high-resolution annotations were rejected. The later correction now uses a
  metadata-only reader, but the retry was stopped before completion while
  changing missing-cell semantics.
- The raw `_debug` input contains five sample IDs, while the available
  processed debug H5AD contains 12 different sample IDs. Its temporary
  Harmony/res50 output is therefore invalid and must not be consumed.
- No production MOFAcellulaR work, package installation, Pixi mutation, or
  existing Pipeline 1–5 launch occurred.

### Current blockers and remaining work

1. The compact snapshot extractor job `4398845` was stopped before
   publication because the immediate figure request only required a faster
   metadata-only matrix. The exact run-owned snapshot workflow remains
   available for the later full subsetting analysis.
2. Full cell-subsetting output has not been produced. The source-level
   implementation now drops missing high-resolution cells before grouping and
   sampling and records total versus annotated counts.
3. MOFAcellulaR has not been executed. The existing Bamboo environment has
   `MOFA2` and `reticulate`, but not `MOFAcellulaR`; no package installation or
   Pixi mutation was performed. The wrapper remains conditional and
   source-only.
4. Production MOFAcellulaR remains deferred until a package strategy and
   reviewed debug result exist. Existing MOFA2 `1.20.2` and `mofapy2 0.7.3`
   pins remain unchanged.
5. The user clarified that `_debug` may contain an expanded sample universe.
   The wrapper therefore records the current unique Sample IDs and requires
   at least two valid samples instead of assuming five fixed IDs.

### Immediate figure result

The fast Bamboo obs-only job `4398864` read the 11 configured
`benchmark_analysis` H5ADs and wrote compact local artifacts under
`data/benchmark/results/derived/fast_obs_20260910/`. The updated
`Figure_3_A_annotationmethods_barplot_anosim.pdf`,
`Supp_fig_18_annotationmethods_barplot_mod_ari.pdf`, and
`Figure_X3_B_number_of_celltypes.pdf` include both
`ECODA_Leiden_res_50` and `ECODA_Leiden_res_2_harmony`; Figure X3B starts its
log10 y-axis at 2 cell types, labels it `Number of cell types (log10 scale)`,
uses intermediate 2–3–5 ticks through 1000, and remains 5×5 inches so the
legend is not clipped.

### Resource and execution boundary

- Do not stage full benchmark H5ADs on the 64-GB-RAM workstation with less
  than 100 GB free disk.
- Do not use local subset mirrors as full-cohort sources.
- Do not invoke any existing `src/` Pipeline 1–5 submitter or repair path.
- Future full-source extraction must remain standalone, explicitly scoped,
  run-owned, and metadata-only for these derived analyses.
