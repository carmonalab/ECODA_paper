# Batch uncorrected analysis implementation plan

## Context

Create the first batch-effect analysis notebook at `notebooks/batch_effect_analysis_uncorrected.rmd`. It must analyze the twelve canonical batch-effect datasets using the already persisted uncorrected method distance artifacts under `data/batch_effect/uncorrected`; it must not rerun preprocessing, recompute method embeddings, read the legacy notebook inputs, or use corrected artifacts. The later corrected notebook remains deferred until this notebook’s evidence has been reviewed and method-specific correction covariates have been selected.

The analysis must compare exactly these seven methods: `ECODA_authors_HR`, `ECODA_seuratres_2`, `Pseudobulk_hvg2000`, `GloScope_hvg2000_pcadims30`, `MrVI_hvg2000`, `PILOT_hvg2000`, and `QOT_hvg2000`. It must produce method-specific MDS/ANOSIM, univariate PERMANOVA, biology-informed PERMANOVA, full-candidate joint PERMANOVA plus global decomposition, grouped score tables, and one metadata-only NMI collinearity heatmap per dataset.

## Approach

### 1. Add a shared, pass-scoped R analysis contract

Add `src/utils/batch_effect_analysis.R`. No equivalent shared R helper currently owns this exact pass-qualified artifact inventory, strict sidecar/sample-order checks, method-specific metric tables, joint decomposition, and NMI output.

Implement these exact helpers:

- `batch_uncorrected_method_specs()` returns the immutable seven-method named specification list. Each entry contains its artifact kind and exact path/key contract:
  - `ECODA_authors_HR`: `results/<ds>_batch_effect_uncorrected_composition.rds`, bundle key `ECODA_authors_HR`;
  - `ECODA_seuratres_2`: the same composition bundle, key `ECODA_seuratres_2`;
  - `Pseudobulk_hvg2000`: `results/<ds>_batch_effect_uncorrected_pseudobulk.rds`, bundle key `Pseudobulk_hvg2000`;
  - `GloScope_hvg2000_pcadims30`: `results/<ds>_batch_effect_uncorrected_gloscope.rds`, bundle key `GloScope_hvg2000_pcadims30`;
  - `MrVI_hvg2000`: `embeddings/<ds>_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather`;
  - `PILOT_hvg2000`: `embeddings/<ds>_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather`;
  - `QOT_hvg2000`: `embeddings/<ds>_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather`.
  The shuffled `ECODA_authors_HR_NULL`, HiTME/scATOMIC entries, scPoli, PILOT-GM-VAE, MOFA, and scITD are not plotted or loaded by this notebook.

- `validate_batch_artifact(path, kind = "artifact")` enforces the existing `MD5=`, `SIZE=`, `PATH=` sidecar contract used by `src/5_run_benchmark_methods/benchmark_hpc_utils.R:545-555`: regular nonempty file, required `.md5`, exact recorded path, exact byte size, and exact MD5. It stops on every mismatch.

- `write_batch_metadata_sidecar(metadata, path, expected_sample_ids = NULL)` writes a temporary Feather file, atomically renames it, writes a strict `.md5` sidecar, and revalidates it. It requires a unique nonblank `Sample` column; when `expected_sample_ids` is supplied, order and membership must match exactly. This is the only writer for the new sample metadata sidecars.

- `read_batch_metadata_sidecar(path, expected_sample_ids, required_label)` validates the sidecar and returns one row per standardized `Sample`. It requires `Sample` and the dataset’s primary biological label. Registered secondary or technical candidates absent from the sidecar are retained as unavailable candidates with explicit status; the loader never substitutes another column.

- `batch_candidate_registry(config, dataset_order)` imports `notebooks/dataset_onboarding/dataset_specs.py` through `reticulate` and combines it with `read_datasets_json(view = "batch_effect_uncorrected")`. Use `BATCH_EFFECT_DATASET_ORDER` as the twelve-dataset order, `BATCH_EFFECT_SPECS` as the technical-candidate registry, and `DATASET_SPECS[[ds]]$sample_stable_cols` minus `bio_col` as secondary biology for the nine cohorts that have that specification. Joanito, Stephenson, and CombinedPBMC have no invented secondary-biology fields. Preserve original metadata names for display, but assign safe internal aliases for model formulas.

- `read_batch_method_dist(input_root, dataset, method_spec, expected_sample_ids)` loads one persisted distance matrix. RDS bundles must contain finite square `dist_mat` values with nonblank unique row/column IDs; Feather matrices must use the final Feather column as the sample ID, have finite numeric square values, and have column names exactly equal to that ordered ID vector. The returned matrix must match `expected_sample_ids` exactly in both membership and order. Do not call `align_result_samples()` to reorder persisted matrices; fail closed on order mismatch.

- `load_batch_uncorrected_dataset(input_root, metadata_root, dataset, registry)` validates the metadata sidecar, the corresponding `<ds>_batch_effect_uncorrected_metadata.rds` summary bundle, and all seven method artifacts. It compares the RDS metadata label names/values and sample count to the sidecar and every method matrix. It returns metadata, covariate classes, method distance objects, and artifact paths. It must not search for corrected or legacy filenames and must not fall back to H5AD or onboarding subset files.

- `compute_batch_anosim(dist_mat, grouping, permutations = 999L)` returns the ANOSIM statistic and permutation p-value for one covariate after the caller has validated the grouping. Use `vegan::anosim()` on the supplied persisted distance object; never build a new feature or CLR distance.

- `compute_batch_permanova(dist_mat, metadata, candidate, biology_col = NULL, permutations = 999L)` returns R², pseudo-F, raw p-value, status, and warning for either a candidate-only or biology-informed model. Candidate-only uses `dist ~ candidate`; biology-informed uses `dist ~ bio_col + candidate` and reports the candidate term. Build formulas from safe temporary metadata aliases with `stats::reformulate()`; never interpolate raw names such as `Cognitive status` or `CoVID-19 severity`. Use 999 permutations and `vegan::adonis2(..., by = "margin")`.

- `compute_batch_joint_permanova(dist_mat, metadata, candidate_registry, permutations = 999L)` builds one complete-case model containing the primary biology, every valid secondary biological candidate, and every valid technical candidate. It returns:
  - candidate-level marginal R², pseudo-F, raw p-value, and status for every estimable full-model term;
  - Holm-adjusted p-values across valid candidate terms for each dataset/method;
  - global `Unique Biological`, `Unique Technical`, `Shared / Confounded`, and `Residual / Unexplained` R² components.
  The global decomposition must use the same method-specific `dist_mat`, complete-case sample rows, safe design aliases, and full-model terms as the candidate bars. Compute `Unique Biological` as the sum of full-model marginal R² terms classified biological, `Unique Technical` as the sum of technical terms, `Shared / Confounded` as `max(0, full_model_R2 - unique_biological - unique_technical)`, and residual as `max(0, 1 - full_model_R2)`. If the full design is rank-deficient or a term is aliased, preserve the affected candidate as `NON_ESTIMABLE` with a warning and do not fabricate zero values; mark the global decomposition unavailable for that method rather than producing a misleading stack.

- `compute_batch_nmi(metadata, candidate_registry)` reproduces the onboarding NMI contract from `build_batch_candidate_evidence.R:157-165` / `onboarding_utils.py:2001-2065`: sample-level rows, exclude `Sample`, drop blank/sentinel values, require at least two levels, use entropy-normalized arithmetic NMI in `[0, 1]`, and return `NA` for insufficient pairs. `plot_batch_nmi_heatmap()` uses the ordered candidate matrix, 0–1 scale, annotated cells, gray missing cells, and visibly documents NMI `> 0.70` as severe collinearity.

- `make_batch_metric_table()` and `make_batch_joint_table()` return normalized long-form data frames with original dataset/method/covariate names, covariate class, sample count, statistics, adjusted p-values, availability status, artifact path, and warnings. `make_batch_nmi_table()` returns dataset/covariate-pair/NMI/severity rows. All output rows retain invalid candidates as status records even when they are excluded from plotted bars.

### 2. Provision and validate the uncorrected metadata sidecars

Before notebook analysis, create one sidecar per canonical dataset at:

```text
data/batch_effect/uncorrected/metadata/<dataset>_sample_metadata.feather
```

Use the exact processed `batch_effect_uncorrected` H5AD view resolved through `datasets.json`, the existing counts-free `load_h5ad_sample_metadata()` reader in `src/5_run_benchmark_methods/benchmark_hpc_utils.R:237-274`, and the existing first-row-per-sample ordering contract in `src/utils/seurat_utils.R:333-362`. Include `Sample`, the primary biological label, and every available registered secondary/technical candidate under its original metadata name. Do not use `data/new_dataset_checks/subsets`, raw source files, or legacy RDS metadata for this export.

Write the sidecars through `write_batch_metadata_sidecar()`. Validate each sidecar against the ordered labels in `<ds>_batch_effect_uncorrected_metadata.rds` and all seven persisted method matrices before any plots are made. Corrected sidecars use the same schema and are deferred with the corrected notebook.

### 3. Create the uncorrected analysis notebook

Add `notebooks/batch_effect_analysis_uncorrected.rmd` using the setup conventions from `notebooks/benchmark_analysis.rmd:38-149` and the existing `rmdformats::downcute` format. The first setup chunk must:

- normalize the working directory when launched from `notebooks/`;
- source `src/utils/load_all_functions.R` and `src/utils/batch_effect_analysis.R`;
- read `datasets.json` through `read_datasets_json(view = "batch_effect_uncorrected")`;
- import the onboarding registry and reorder to the exact twelve-dataset order;
- set `input_root = file.path(getwd(), "data", "batch_effect", "uncorrected")`, `metadata_root = file.path(input_root, "metadata")`, `plot_root = file.path(getwd(), "plots", "batch_effect_uncorrected")`, and `analysis_output_root = file.path(input_root, "analysis")`;
- set `set.seed(42)` and `n_permutations <- 999L`;
- fail before looping if any dataset, sidecar, or one of the seven method artifacts fails validation.

Organize the notebook in this order:

1. Libraries, registry, paths, and strict artifact preflight.
2. Per-dataset metadata/covariate summary and sample-order assertions.
3. One seven-panel MDS figure per dataset. Reuse `plot_mds()` from `src/utils/plotting.R:164-250` with `cluster_score = FALSE`, `mod_score = FALSE`, `sil_score = FALSE`, `anosim_score = TRUE`, and the primary biological label. Do not recompute any feature matrix or distance matrix.
4. One grouped ANOSIM figure per dataset. Use covariate on the x-axis, method as the fill/group, primary biology first, secondary biology second, technical candidates last, and separate biological/technical facets where necessary. Plot only valid scores; retain unavailable candidates in the table.
5. One grouped univariate PERMANOVA figure per dataset using candidate-only R².
6. One grouped biology-informed PERMANOVA figure per dataset using candidate R² from `bio_col + candidate`; omit the primary biology as a candidate in this panel.
7. One two-panel joint PERMANOVA figure per dataset. The upper panel is the candidate-level marginal R² from the full model across all valid candidates and methods. The lower panel is one stacked global decomposition bar per method using the four qmd-style components. Both panels use the same persisted method distance matrix and complete-case model.
8. One NMI heatmap per dataset after all MDS and metric figures. It is metadata-only and has no method facet.
9. Write normalized score, decomposition, and NMI tables atomically under `data/batch_effect/uncorrected/analysis/`, each with an MD5/SIZE/PATH sidecar. Do not modify the source artifact checksum manifest.

Use `method_label_map_main` from `src/utils/constants.R` for display names. Use patchwork/ggplot2 patterns already present in `notebooks/benchmark_analysis.rmd:550-668`; do not introduce publication `Figure*` or `Supp_fig*` names for these exploratory batch outputs.

Use these exact output names per dataset:

```text
plots/batch_effect_uncorrected/<dataset>_mds.pdf
plots/batch_effect_uncorrected/<dataset>_anosim.pdf
plots/batch_effect_uncorrected/<dataset>_permanova_univariate.pdf
plots/batch_effect_uncorrected/<dataset>_permanova_bio_informed.pdf
plots/batch_effect_uncorrected/<dataset>_permanova_joint.pdf
plots/batch_effect_uncorrected/<dataset>_nmi.pdf
```

Use these exact analysis tables:

```text
data/batch_effect/uncorrected/analysis/batch_effect_uncorrected_scores.csv
 data/batch_effect/uncorrected/analysis/batch_effect_uncorrected_decomposition.csv
 data/batch_effect/uncorrected/analysis/batch_effect_uncorrected_nmi.csv
```

The score table must contain one row per dataset × selected method × registered candidate, including ANOSIM, univariate PERMANOVA, biology-informed PERMANOVA, joint candidate R²/p-values, status, and warnings. The decomposition table must contain one row per dataset × selected method × component with sample count and status. The NMI table must contain one row per dataset × ordered covariate pair with NMI and the `>0.70` severity flag.

### 4. Add focused synthetic contract coverage

Add `tests/test_batch_effect_analysis.R` as a standalone `Rscript` test. Use a temporary artifact root and deterministic synthetic data; do not read a production cohort. The fixture must:

- create a six-or-more-sample metadata sidecar with a non-syntactic biological column such as `CoVID-19 severity`, one secondary biological field, and one technical field;
- create checksummed minimal RDS/Feather artifacts under the exact seven uncorrected method naming contracts;
- assert that the loader discovers exactly the seven selected methods, rejects corrected/legacy fallbacks, validates exact sample IDs/order, and fails on missing or stale sidecars;
- assert that `plot_mds()` produces an MDS object whose title contains the ANOSIM score;
- assert finite univariate ANOSIM and candidate-only PERMANOVA values;
- assert finite biology-informed candidate values and that the primary biology is not emitted as a biology-informed candidate bar;
- assert full-model joint candidate rows include all estimable registered terms, Holm-adjusted p-values are present, and the global decomposition has the four exact component names and sums to one within tolerance;
- assert that a rank-deficient/aliased synthetic candidate is marked `NON_ESTIMABLE` and is never converted to zero;
- assert the NMI matrix is symmetric, has diagonal one, preserves candidate labels, and flags a deliberately confounded pair above 0.70;
- assert table row completeness and status/warning propagation.

Do not add tests that modify `datasets.json`, onboarding README, corrected artifacts, or scheduler state.

## Critical files & anchors

- `notebooks/batch_effect_analysis_uncorrected.rmd` — new orchestration notebook; setup, twelve-dataset loop, six per-dataset plot outputs, and atomic table outputs.
- `src/utils/batch_effect_analysis.R` — new strict loader, metadata-sidecar contract, method-specific metrics, joint decomposition, NMI, and plotting/table helpers.
- `tests/test_batch_effect_analysis.R` — new deterministic synthetic contract test for the observable analysis behavior.
- `notebooks/dataset_onboarding/dataset_specs.py:330-377` — authoritative twelve-dataset order and registered technical candidates; import, do not duplicate or edit.
- `src/utils/plotting.R:164-250` — existing `plot_mds()`/ANOSIM implementation to reuse without changing its contract.

## Verification

Run all commands from the repository root with the Pixi default environment.

1. Run the new focused contract test:

   ```bash
   pixi run -e default Rscript --vanilla tests/test_batch_effect_analysis.R
   ```

   Expected observable result: the test prints its success line; all seven synthetic method artifacts load; non-syntactic formulas, all requested metric tables, the four-component decomposition, and NMI matrix pass; corrupted/misaligned fixtures fail closed.

2. Run the existing focused artifact/routing contracts that guard the inputs reused by the notebook:

   ```bash
   pixi run -e default python tests/test_batch_effect_registry_and_modes.py
   pixi run -e default python tests/test_benchmark_matrix_validator.py
   pixi run -e default Rscript --vanilla tests/test_batch_rds_contract.R
   pixi run -e default python tests/test_h5ad_pseudobulk.py
   bash tests/test_ecoda_run_common.sh
   bash tests/test_checksum_reuse.sh
   bash tests/test_atomic_artifact_writers.sh
   ```

   Expected observable result: all commands pass without submitting workers, changing scheduler state, or modifying source artifacts.

3. Before rendering against the real local artifact mirror, run the existing exact batch matrix/RDS validators against the run-owned uncorrected selection and source identity. Use the run’s validated artifact root, selection file, source H5AD input root, and source-identity record; the commands must use:

   ```text
   matrix_artifact_validator.py --root ... --selection ... --labels prepare_pseudobulk pseudobulk gloscope composition mrvi pilot qot --batch --batch-pass uncorrected --exact --input-root ... --config datasets.json --source-identity ...
   validate_benchmark_rds_contract.R --root ... --selection ... --labels gloscope,pseudobulk,composition --batch-pass uncorrected --input-root ... --config datasets.json
   ```

   Expected observable result: all twelve datasets and all seven method artifacts have valid checksums, finite square matrices, identical ordered sample IDs, and valid metadata before notebook execution. A missing/stale sidecar or corrected/legacy artifact must stop the gate.

4. After the uncorrected metadata sidecars have been provisioned and validated, render the actual notebook with the repository environment:

   ```bash
   pixi run -e default Rscript --vanilla -e 'rmarkdown::render("notebooks/batch_effect_analysis_uncorrected.rmd", output_dir="plots/batch_effect_uncorrected/rendered")'
   ```

   Expected observable result: rendering completes without warning-level data fallback; twelve MDS PDFs, twelve ANOSIM PDFs, twelve univariate PERMANOVA PDFs, twelve biology-informed PERMANOVA PDFs, twelve joint/decomposition PDFs, twelve NMI PDFs, and the three checksummed analysis tables exist. Every table contains exactly the twelve canonical datasets and seven selected methods where the method-level metric is estimable.

5. Inspect one rendered dataset (Joanito) and one non-syntactic-label dataset (Covid19_PBMC) to verify visually that each MDS panel shows its ANOSIM title, each grouped bar legend contains the seven selected methods, the joint figure has candidate bars above method-specific decomposition stacks, and the NMI heatmap uses the same candidate labels as the tables. This is analysis-output inspection only; no full-cohort preprocessing or benchmark launch is allowed.

## Assumptions & contingencies

- The uncorrected method artifact inventory is fixed to the seven names above. If any selected artifact is absent or fails validation, stop before plotting; do not silently replace it with `ECODA_authors_HR_NULL`, an annotation-specific method, a corrected artifact, or a legacy file.
- Metadata sidecars must come from the exact processed `batch_effect_uncorrected` H5AD views and use the standardized `Sample` universe. If a sidecar is missing, stale, incomplete, duplicated, or out of order, stop; never read onboarding subset H5ADs or infer sample covariates from a method matrix.
- Missing, constant, incomplete, sample-unique, or perfectly biology-confounded candidates are emitted as status/warning rows and omitted from bars. They are never represented by zero scores.
- If a biology-informed or full joint model is non-estimable because of rank deficiency, preserve the candidate/status record and omit the affected bar. Do not assign a zero. If the full joint design is non-estimable, leave the global decomposition unavailable for that dataset/method and retain the NMI evidence.
- `Cognitive status`, `CoVID-19 severity`, `Single cell sequencing platform`, and other non-syntactic names remain unchanged in metadata and displayed output; only internal formula aliases are safe names.
- The notebook reads and validates existing method distances; it never recalculates CLR, pseudobulk, embeddings, PCA, or method distances. The QMD implementation supplies the four decomposition labels and candidate-screening behavior only; all numerical method metrics in this notebook come from persisted `dist_mat` objects.
- The corrected notebook and corrected sidecar provisioning are not part of this implementation. The helper’s sidecar schema and method registry are pass-aware so the later corrected task can reuse them without changing this uncorrected contract.
