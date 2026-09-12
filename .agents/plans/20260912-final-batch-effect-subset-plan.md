# ECODA final batch-effect subset processing

## Context

Implement the approved final batch-effect analysis without rerunning completed
cohorts or broad historical selections. The production source of truth is the
current `datasets.json` plus the authoritative full-cohort data on Bamboo/HPC;
local subset mirrors are diagnostic only. The final run must correct the Covid
same-column subset rule, regenerate the remaining changed batch views, recover
only the missing Kidney uncorrected Stage 5 rows, synchronize the new results to
the workstation, and run the final uncorrected analysis notebook without
writing into the legacy analysis lane.

The repository already contains user-approved changes to `AGENTS.md` and the
concise onboarding README. Preserve unrelated working-tree changes. Do not
modify or delete existing H5ADs, RDS bundles, pseudobulks, Feather files,
plots, manifests, gates, or checksums. Existing legacy artifacts remain
available as read-only sources for the mixed final analysis.
## Implementation orchestration

Use subagents for the substantial independent implementation units to reduce
context load and keep long-file edits isolated. Each implementation subagent
must receive the relevant contract from this plan, must not edit outside its
assigned files, and must skip formatters, linters, and project-wide tests.
The parent agent owns phase-level verification and integration.

After the subset-rule and final-variant contracts are fixed as written here,
the implementation proceeds in dependency waves. Substantial disjoint edits
use subagents; the parent agent owns phase-level verification and integration.

Wave 1 is the subset/preflight contract unit. It owns
`src/utils/py/preprocess_utils.py`,
`src/3_scrnaseq_preprocessing/1.0_audit_input_views.py`,
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`,
`src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, and the new
`src/utils/bash/h5ad_obs_audit_worker.sh`. A separate focused regression
subtask owns `tests/test_subset_vars.py`. These changes define the mask,
sample-consistency audit, strict direct-H5AD obs-only preflight, and the
eight-row Stage 3 release gate.

After Wave 1 passes the focused test, the Stage 5 variant unit and final
analysis unit run in parallel because they own disjoint files. The Stage 5
unit owns `src/5_run_benchmark_methods/1_submit_hpc_array.sh`,
`src/utils/bash/ecoda_run_common.sh`,
`src/5_run_benchmark_methods/benchmark_submit_common.sh`,
`src/5_run_benchmark_methods/validate_benchmark_rds_contract.R`, the
Python/R Stage 5 workers and wrappers, and
`src/utils/py/export_h5ad_sample_metadata.py`. It owns final roots, stems,
provenance, direct-path validation, metadata export, worker propagation, and
synchronization.
The final analysis unit owns `src/utils/batch_effect_analysis.R` and
`notebooks/batch_effect_analysis_uncorrected.rmd`; it owns the final-nine
scope text, two mixed-source manifests, manifest-aware loading, bundle-key
semantics, and final-only analysis outputs.

A separate Stage 5/manifest contract test subtask owns the existing focused
tests `tests/test_benchmark_selection_file.sh`,
`tests/test_benchmark_sync.sh`, `tests/test_benchmark_rds_contract.R`, and
`tests/test_batch_effect_analysis.R`. It adds assertions for legacy-path
compatibility, final variant roots/stems, direct pseudobulk result paths,
composition/null shared bundles, explicit empty Feather bundle keys, and
per-artifact-kind sample-ID rules. These tests are written in parallel but
run only after the owning source units are integrated.

After the subset test passes, the registry/configuration unit owns only
`datasets.json`, `AGENTS.md`, and `tests/test_batch_effect_registry_and_modes.py`.
It changes the Covid rule, the eight target H5AD names, stale scope paragraphs,
and focused assertions; it must preserve unrelated user edits and must not
touch frozen rows.
The parent then runs the focused configuration/routing checks and integrates
all source/test units before any HPC launch.

Stage 2/3/5 launch preparation is serialized after all code/config units are
integrated: Stage 2 Joanito first, Stage 3 eight-row array second, Stage 5
changed-dataset wave third, and targeted Kidney Stage 5 recovery last. The
changed-dataset Stage 5 wave and Kidney recovery share the final analysis root
and therefore do not launch concurrently.


## Decisions and traceability

| Decision or finding | Evidence | Consequence |
|---|---|---|
| Covid filtering is not implemented | `src/utils/py/preprocess_utils.py:356-366` only calls `Series.isin()` and treats every operator other than `in` as complement; `datasets.json:513-522` uses scalar `30` with `op: "<="` | Fix the shared evaluator and the Covid rules before any real target row runs. |
| Covid day values are categorical strings | `data/new_dataset_checks/subsets/Covid19_PBMC_meta.json` is a historical local diagnostic only: it reports a categorical sampling-day field, numeric-looking levels, and `control`, but has no current HPC source identity or predicate result | Before any Stage 3 compute, run the separate read-only HPC `obs` preflight described in step 4; record the current input identity, dtype/levels, and predicate counts there. The preflight, not this local JSON, authorizes the target run. |
| Final Covid predicate | User decision | Retain sampling-day `control` or finite numeric `<= 30`; exclude blank, unknown, malformed, and missing values. Apply to both corrected and uncorrected views. |
| Subset evaluation scope | Stage 3 applies subset rows before sample standardization, sample-count filtering, normalization, HVG/PCA, and corrected embeddings (`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:717-743`) | Do not slice existing processed H5ADs after the fact; regenerate changed views from authoritative HPC inputs. |
| Standard versus obs-only preflight | `src/utils/bash/h5ad_preflight_worker.sh:253-267` can publish an artifact record even when validating an existing H5AD | Covid source inspection and final metadata export use the separate no-record `h5ad_obs_audit_worker.sh`; they never attach a record to the immutable input. |
| Stage 5 metadata source | `benchmark_hpc_utils.R:782-818` reads sample metadata, but Stage 5 workers only publish method RDS/embedding and execution-time artifacts | Add an explicit obs-only exporter and checksum its Feather output before final manifest creation; never assume a method worker created it. |
| Pseudobulk path contract | `1.1.1_run_benchmark_methods_r.R:337-341,419-423` names the result as `<stem>_<method>.rds`; the cache is separately named by `1.1.1_prepare_pseudobulk.R:293-296,411-470` | Manifest `Pseudobulk_hvg2000` points to the final result bundle; the pseudobulk cache is a dependency-only path. |
| Final registry size | `batch_candidate_registry():221-240` and `dataset_specs.py:486-502` enforce the historical twelve-row order | Generalize only the R registry caller to accept an explicit final nine-row subset while retaining the historical twelve-row assertion. |
| Corrected Stage 5 boundary | User decision and final lane contract | `corrected_final` is Stage 3-only in this run; no corrected Stage 5 work or artifact root is created. |
| Batch annotation scope | User-approved `AGENTS.md` exception | Batch-effect views do not run Pipeline 4. Preserve configured source/author cell-type columns. |
| Frozen cohorts | User decision | `Alzheimer`, `Breast_cancer`, `Lupus_PBMC`, and `Stephenson` are final, must not occur in any new job, validator selection, or compute manifest, and are included in final notebook plots only through their approved legacy result artifacts; no legacy dataset H5AD is read for the final analysis. |
| Changed final-view targets | User decision and current registry | `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`, both uncorrected and corrected Stage 3 views. |
| Kidney recovery | User decision | `Kidney_KPMP_full` receives only missing uncorrected Stage 5 method rows; do not run Stage 2/3/4 for Kidney. |
| Disabled cohorts | Current registry flags | `CombinedPBMC`, `Kidney_KPMP`, `Myocardial_infarction`, and `Parkinson` receive no new work. `_debug` is not a production target. |
| Stage 2 coverage | `src/2_dataset_specific_preprocessing/1_submit_hpc.sh:233-259,394-444` | Only `Joanito` has a target-specific Stage 2 hook. Covid, Diabetes, and Lung use already staged direct H5AD inputs; do not invent Stage 2 jobs for them. |
| Corrected batch variables | `datasets.json` target contracts | The corrected Stage 3 views use `Covid19_PBMC=datasets`, `Diabetes=dataset`, `Joanito=seqtec`, and `Lung=dataset`; biological labels remain evaluation-only. |
| Final naming | User decision | Target H5ADs and uncorrected Stage 5 artifacts use `_final`; `corrected_final` is reserved for corrected Stage 3 H5AD/view outputs and is not a Stage 5 root in this run. |
| Notebook boundary | User decision and current notebook code | `batch_effect_analysis_uncorrected.rmd` becomes final-only; the contingency notebook remains legacy-only and untouched; `batch_effect_analysis_legacy.rmd` remains out of scope. |
| Current gate policy | `AGENTS.md`, durable profile, submitters | Existing durable-gate/snapshot/accounting/checksum policy remains unchanged for new target work. Avoid redundant rows, but do not bypass required target-run contracts. |

## Approach

### 1. Establish the immutable target manifests

Create run-owned, headerless selection manifests before any scheduler launch.
Never derive these manifests from the default historical twelve-row list.

Stage 2 logical selector (passed to the Stage 2 submitter, not a
headerless dataset/step scheduler manifest):

```text
src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Joanito \
  --steps joanito
```

The durable wrapper must materialize and checksum its own run-scoped
selection/step manifest from that selector and record exactly one `Joanito`
`joanito` row. Do not hand this selector text to a manifest parser expecting
`DATASET<TAB>VIEW` rows.

Stage 3 selection, exactly eight rows:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected
Covid19_PBMC<TAB>batch_effect_corrected
Diabetes<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_corrected
Joanito<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_corrected
Lung<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_corrected
```

Stage 5 final changed-dataset selection, exactly four rows:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
```

Kidney Stage 5 selection, exactly one dataset row:

```text
Kidney_KPMP_full<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
```

The four frozen cohorts, `_debug`, and all disabled cohorts must be absent from
every new scheduler selection and validator input. Do not pass
`--exact-batch-selection`, because that mode requires the obsolete historical
twelve-row matrix.

### 2. Replace the shared subset evaluator

Update `src/utils/py/preprocess_utils.py` while preserving the public signature
`apply_subset_vars(adata, subset_vars, copy=True)` and its `copy=True/False`
behavior.

Add a small shared evaluator, such as
`evaluate_subset_mask(adata, subset_vars) -> pandas.Series`, and a sample-level
consistency helper, such as
`assert_subset_sample_consistency(adata, mask, sample_col, context)`. The new
wrapper must use the evaluator rather than maintaining a second implementation.

Implement this exact rule contract:

- `in` retains exact membership and `notin` excludes exact membership.
- Scalar `values` are normalized to one-item sequences, so existing scalar
  Joanito/Lung rules remain valid.
- `<=`, `<`, `>=`, and `>` parse trimmed categorical/string values with
  `pandas.to_numeric(errors="coerce")`; malformed, blank, missing, and
  non-finite values are false for comparison masks.
- Comparison thresholds require one finite numeric `values` value.
- `include_values` is optional, normalized to a sequence, and ORed into the
  positive comparison mask after trimmed, case-insensitive string matching.
- The Covid rule therefore evaluates as `numeric_day <= 30 OR day ==
  "control"` on the sampling-day column only.
- Multiple declared columns remain combined with boolean `AND`; the
  `include_values` exception is local to its own column rule.
- Missing columns, malformed rule objects, missing operators, invalid
  thresholds, and unknown operators raise explicit `KeyError`/`ValueError`
  failures. Never silently treat an unknown operator as exclusion.
- The returned mask is indexed exactly by `adata.obs_names`.

Before slicing a final view, call the sample-level helper with the configured
raw sample column. It must fail if one sample has both retained and dropped
cell rows under the final mask. This check is row-level evaluation followed by
sample-level consistency auditing; it must use configured `sampleID` for Covid,
`sample.ID` for Joanito, `sample` for Lung, and `donor_id` for Diabetes. Do not
collapse repeated `PatientID` values into one Covid sample.

Update the Stage 3 preprocessing caller
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:717-743` to evaluate and
audit the mask before slicing, then retain the existing downstream order:
subset, standardize `Sample`, remove samples below 500 cells, and compute
view-specific representations. Update the input-view audit caller
`src/3_scrnaseq_preprocessing/1.0_audit_input_views.py:124-131` to use the same
mask/audit helper. Keep the CombinedPBMC caller
`src/2_dataset_specific_preprocessing/1.2.1_create_combinedpbmc_dataset.py:108-120`
on the shared wrapper even though CombinedPBMC is disabled.

Do not make diagnostic onboarding specs/QMDs into a second production
configuration. They remain diagnostic and must not override `datasets.json`.
If their local helper is exercised by a new test, it must delegate to the same
shared evaluator rather than copy comparison logic.

### 3. Encode the approved Covid rule and final H5AD names

After the evaluator is implemented and tested, update only the necessary
`datasets.json` entries:

1. In both `Covid19_PBMC.views.batch_effect_uncorrected.subset_vars` and
   `.batch_effect_corrected.subset_vars`, use:

   ```json
   {
     "Sampling day (Days after symptom onset)": {
       "values": 30,
       "op": "<=",
       "include_values": ["control"]
     }
   }
   ```

2. Change only the output names for the four regenerated datasets, in both
   batch views, by inserting `_final` before `_ECODAprocessed.h5ad`:

   ```text
   Covid19_Ren2021_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad
   Covid19_Ren2021_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad
   diabetes_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad
   diabetes_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad
   JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad
   JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad
   lungatlas_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad
   lungatlas_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad
   ```

Leave all existing legacy output names, including the frozen cohorts and the
old changed-dataset outputs, in place on HPC. The new names are the canonical
paths for the final target views after this change; old files are not deleted
or overwritten.
Before implementation launches, replace the stale batch-effect baseline
paragraphs in `AGENTS.md:177-219` with the approved current scope: final
Stage 3 targets are `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`; the
final Stage 5 changed wave uses those four; `Kidney_KPMP_full` receives only
targeted missing uncorrected Stage 5 rows; `Alzheimer`, `Breast_cancer`,
`Lupus_PBMC`, and `Stephenson` are frozen and absent from all new jobs and
validator selections; disabled cohorts and `_debug` are absent from production
selection. Preserve the existing durable-gate and snapshot requirements while
updating only the stale scope statements.
The registry/configuration unit must update the focused registry contract
test (currently `tests/test_batch_effect_registry_and_modes.py`) to cover the
final nine-dataset order and the eight manifest method keys, while retaining
the historical twelve-dataset assertion for the legacy contingency path. It
must assert that only the four regenerated target view names gain `_final`;
frozen, disabled, and legacy changed-dataset names remain unchanged.


### 4. Add the synthetic subset regression before HPC data

Add `tests/test_subset_vars.py` following the standalone style of
`tests/test_preprocessing_sample_filter.py`.

Use two tiny AnnData fixtures rather than one ambiguous sample table.

The Covid predicate fixture has observation IDs `c0` through `c5` and:

```text
Sample:       A, A, B, C, D, E
sampleID:     sA, sA, sB, sC, sD, sE
PatientID:    pA, pA, pB, pC, pD, pE
sampling_day: 29, 30, 30.5, control, unknown, malformed
```

Store `sampling_day` as a pandas categorical/string column. Apply the Covid
rule and assert that exactly `c0`, `c1`, and `c3` are retained (`A/29`,
`A/30`, and `C/control`); `c2`, `c4`, and `c5` are dropped. This explicitly
covers the inclusive `30` boundary, numeric-looking categorical strings,
literal `control`, and malformed/unknown values.

Use a separate split-sample fixture:

```text
Sample:       S, S
sampleID:     sS, sS
PatientID:    pS, pS
sampling_day: 30, 30.5
```

The row mask is `[True, False]`; `assert_subset_sample_consistency()` must
raise because one configured `sampleID` is split. A separate repeated-identity
case uses `PatientID=[p-repeat,p-repeat]`, `sampleID=[s1,s2]`, and days
`[30,31]`; it must not raise merely because `PatientID` repeats.

Also include explicit assertions for scalar `in`/`notin` values, missing/blank
and non-finite comparison values, unknown-operator rejection, and preservation
of `copy=True/False` behavior.

Run only the focused test before real data:

```bash
pixi run -e default python tests/test_subset_vars.py
```

Before any Stage 3 compute, perform a distinct, read-only HPC metadata
preflight for the current Covid direct input. Extend
`src/3_scrnaseq_preprocessing/1.0_audit_input_views.py` with an explicit
`--obs-only` mode and a direct `--input-file PATH` option. In this mode,
`--input-file` is required, must resolve to an existing `.h5ad`, must match
the configured Covid `input_file_name`, and is read only through the existing
`h5py`/`read_dataframe` observation reader. The mode must reject
`--input-root`-only invocation, RDS inputs, and the `load_input` fallback.
It writes only a run-owned JSON report plus checksum. Implement the separate
worker as `src/utils/bash/h5ad_obs_audit_worker.sh`; it must validate the
immutable source/runtime binding but must never call
`ecoda_write_artifact_record` on the inspected H5AD.

The worker invokes the immutable snapshot’s Python interpreter with the
following two view-specific reports, writing only below the run root:

```text
${PYTHON_BIN} src/3_scrnaseq_preprocessing/1.0_audit_input_views.py \
  --config ${ECODA_SOURCE_ROOT}/datasets.json \
  --input-file ${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad \
  --output-root ${ECODA_RUN_ROOT}/preflight \
  --output ${ECODA_RUN_ROOT}/preflight/Covid19_PBMC_batch_effect_uncorrected.json \
  --view batch_effect_uncorrected --ds-name Covid19_PBMC --obs-only

${PYTHON_BIN} src/3_scrnaseq_preprocessing/1.0_audit_input_views.py \
  --config ${ECODA_SOURCE_ROOT}/datasets.json \
  --input-file ${HPC_SCRATCH_DIR}/Covid19_PBMC/data/Covid19_Ren2021.h5ad \
  --output-root ${ECODA_RUN_ROOT}/preflight \
  --output ${ECODA_RUN_ROOT}/preflight/Covid19_PBMC_batch_effect_corrected.json \
  --view batch_effect_corrected --ds-name Covid19_PBMC --obs-only
```

The durable wrapper records these reports and their checksums as preflight
inputs before releasing the eight Stage 3 rows; it does not treat them as
canonical H5AD outputs. Run the preflight against the current HPC-resolved
direct input before the Stage 3 worker array is released. The report must
include the immutable source/runtime identity, absolute input path,
H5AD checksum/size, `obs` row count, sampling-day dtype, raw unique values
(including `control`, unknown/blank/malformed values), missing/non-finite
counts, configured `sampleID`/`PatientID` cardinalities, retained/dropped
cells and samples for the exact predicate, and split-sample count. It must
fail closed on a missing column, source/checksum mismatch, malformed
metadata, or any split sample. The Stage 3 selection cannot proceed unless
both reports are present, checksummed, and record the exact configured rule.
The local diagnostic JSON is evidence of expected schema only and must not
authorize the production run.


A failure blocks all Stage 2/3/5 launches. Do not use the local Covid subset
H5AD as a production test source; the synthetic fixture and the HPC
obs-only report are the required deterministic/authoritative checks.

### 5. Run the only required Stage 2 hook

Run Stage 2 after the focused test and before Stage 3, using the immutable
snapshot-backed durable workflow and the exact selector `Joanito` plus step
`joanito`:

```text
src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Joanito \
  --steps joanito
```

This selector is not a headerless `DATASET<TAB>VIEW` manifest. The durable
wrapper must materialize and checksum its own run-scoped step manifest from
the selector and record exactly one `Joanito`/`joanito` row. Do not include
Covid19_PBMC, Diabetes, Lung, Kidney_KPMP_full, the four frozen cohorts,
`_debug`, or disabled cohorts. Do not use `--force` unless the selected
Joanito output is specifically identified as invalid; valid `seqtec` and
`cell.type_new` metadata is a no-op and must not be reprocessed.

Covid19_PBMC, Diabetes, and Lung have direct staged H5AD inputs and no Stage 2
hook in the current submitter. Their absence from Stage 2 is intentional. No
Pipeline 4 annotation work occurs for any batch-effect view.

Use the repository’s snapshot/runtime preparation, exact run ID, and
run-scoped selection. After launch, arm the one required unbounded durable
wait, inspect every emitted Stage 2/watchdog ID once, and obtain the required
terminal review before Stage 3. Do not poll repeatedly or rerun an ambiguous
wrapper.

### 6. Regenerate the four final Stage 3 views

Run Stage 3 only with the eight-row manifest from step 1:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected
Covid19_PBMC<TAB>batch_effect_corrected
Diabetes<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_corrected
Joanito<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_corrected
Lung<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_corrected
```

Use the canonical Stage 3 submitter with `--selection-file`; do not use the
historical exact-selection mode and do not include Kidney or any
frozen/disabled cohort. Before the worker array is released, require the
Covid read-only HPC `obs` preflight from step 4 for both declared views. The
preflight report is run-owned evidence of the current source; it is not an
H5AD artifact record and must not write beside the immutable direct input.

The Stage 3 worker must resolve the final output names from the updated
`datasets.json` and write:

- Covid19_PBMC: final uncorrected and corrected views;
- Diabetes: final uncorrected and corrected views using existing
  `cell_type`/`cell_type_reannotatedIntegrated` columns;
- Joanito: final uncorrected and corrected views after the targeted Stage 2
  metadata repair, using `cell.type`/`cell.type_new`;
- Lung: final uncorrected and corrected views using `ann_coarse`/`ann_fine`.

The uncorrected representation remains `Sample`-keyed raw PCA/neighbors/Leiden
without Harmony. The corrected representation uses the configured technical
batch variables and Harmony, with the biological label excluded from all
processing covariates. Both views use the same final subset for each dataset.
The corrected view is a Stage 3 output only; it does not authorize a
corrected Stage 5 run.

Preserve source/author cell-type metadata in every output. Do not call Pipeline
4, prepare annotation chunks, run annotation workers, merge annotation
Feathers, or add HiTME/scATOMIC columns.

The Stage 3 run must be snapshot-backed and use the current durable gate. Do
not set `--force` for a broad target list. If a final target row has a valid
existing final artifact, skip that row; if it is absent or invalid, recompute
only that explicitly selected row.

After terminal review, record the eight output paths and their subset audit
summaries in the run-owned manifest. Each audit must include dataset, view,
configured raw sample column, total/retained/dropped cell counts,
total/retained/dropped sample counts, the exact HPC input identity, and
split-sample count (which must be zero).

### 7. Add a final variant to Stage 5 without changing legacy mode

Extend the canonical Stage 5 wrapper and shared artifact-path helpers with an
explicit `--analysis-variant final` option. The default with no variant must
remain byte-for-byte compatible with the existing legacy roots and stems.
`final` is valid only for the approved uncorrected batch-effect selection
files; reject it with `benchmark_analysis`, a corrected Stage 5 pass, a broad
default selection, or a missing explicit selection file.

Set `ANALYSIS_VARIANT`, `ANALYSIS_ROOT`, `ANALYSIS_NAS_ROOT`,
`ANALYSIS_PASS`, and `ANALYSIS_LOG_PREFIX` before constructing pending
selection state or writing run metadata. Final run metadata must contain:

```text
ANALYSIS_VARIANT=final
ANALYSIS_ROOT=${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final
ANALYSIS_NAS_ROOT=${NAS_TARGET_DIR}/batch_effect/uncorrected_final
ANALYSIS_PASS=uncorrected
PASS=uncorrected
ROOT=${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final
```

Move or refactor the current root assignment before the `RUN_METADATA` block
around `src/5_run_benchmark_methods/1_submit_hpc_array.sh:2398-2425`, so the
metadata cannot record a legacy root for a final run. `ecoda_run_audit.sh`,
artifact ownership records, worker environments, validator path
reconstruction, watchdogs, and synchronization must consume the same
variant-qualified root. The legacy mode must continue to emit
`batch_effect/uncorrected` or `benchmark` roots exactly as before.

Update these exact layers together:

- `src/5_run_benchmark_methods/1_submit_hpc_array.sh`: parse/validate the
  option, establish final roots, export the variant, and record the fields
  above in both active and `NOOP_VALIDATED` metadata.
- `src/utils/bash/ecoda_run_common.sh:_ecoda_stage5_artifacts_for`: derive
  variant-qualified stems for ownership and sync expansion.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh:benchmark_artifacts_for`:
  use the same centralized stem rule for every Stage 5 method.
- `src/5_run_benchmark_methods/benchmark_submit_common.sh:878-945`:
  enumerate the same variant-qualified paths during sync and validation.
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:1857-1863`:
  add `_final` to batch output names when `ANALYSIS_VARIANT=final` while
  keeping embedding keys semantic-view based.
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:337-342,419-423`:
  add `_final` to batch cache/result stems when the final variant is set.
  The result stem is the cache stem followed by the method name; it is not
  the pseudobulk cache filename.
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R`:
  accept the variant-qualified cache root and preserve legacy cache behavior.
- Both worker wrappers must forward `ANALYSIS_VARIANT=final` and use it in
  execution-log filenames.
- `src/5_run_benchmark_methods/validate_benchmark_rds_contract.R` and
  `src/utils/bash/ecoda_run_common.sh` must validate direct final paths and
  never reconstruct a final artifact under the legacy pass root.

Final Stage 5 roots:

```text
${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final
${NAS_TARGET_DIR}/batch_effect/uncorrected_final
```

The name `corrected_final` is reserved for corrected Stage 3 H5AD/view
outputs only. This run creates no corrected Stage 5 root, selection,
metadata, method artifact, sync list, or notebook lane.

Final Stage 5 paths must distinguish the pseudobulk cache from the result
bundle:

```text
${ANALYSIS_ROOT}/pseudobulks/<DS>_batch_effect_uncorrected_final_pseudobulk_hvg2000.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_uncorrected_final_pseudobulk.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_uncorrected_final_gloscope.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_uncorrected_final_composition.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_uncorrected_final_metadata.rds
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_pilot_dists.feather
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_qot_dists.feather
```

The final manifest key `Pseudobulk_hvg2000` maps to
`results/<DS>_batch_effect_uncorrected_final_pseudobulk.rds` and its
`Pseudobulk_hvg2000` bundle key. The `pseudobulks/..._pseudobulk_hvg2000.rds`
path is only the dependency cache. The final composition RDS contains
`ECODA_authors_HR`, `ECODA_seuratres_2`, and
`ECODA_authors_HR_NULL`; the null manifest row intentionally shares the
composition path and selects the null bundle key.

The final Stage 5 suite is exactly:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

Do not select `mofa`, `scitd`, `scpoli`, `pilotgm`, `trans`, `zeroimp`, or
any post-baseline method. The four changed datasets form one explicit final
Stage 5 wave: four dataset rows and seven method rows per dataset (28 method
rows, with `prepare_pseudobulk` as the declared dependency/method row).
Before either Stage 5 gate, record the exact remote wrapper command and
expected row count in the durable-gate manifest. The changed-dataset gate’s
`--exact-command` must be the snapshot `ecoda_source_snapshot.sh exec`
wrapper required by `AGENTS.md`, with this script/argument tail:

```text
--script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- \
  --selection-file <run-root>/manifests/stage5_changed_final.tsv \
  --pass uncorrected \
  --analysis-variant final \
  --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

The recorded selection file is exactly the four changed-dataset rows from
step 1. The wrapper must emit exactly four dataset rows times seven methods,
28 method rows including the declared `prepare_pseudobulk` dependency. No
`--analyses`, `--exact-batch-selection`, broad dataset list, corrected pass,
or `--force` is permitted.

After the changed gate reaches terminal review, construct the Kidney
run-owned missing-method manifest from the current validated legacy inventory
and record the resulting exact method list before submission. Its
`--exact-command` uses the same snapshot wrapper and this tail:

```text
--script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- \
  --selection-file <run-root>/manifests/kidney_missing_final.tsv \
  --pass uncorrected \
  --analysis-variant final \
  --target-methods <comma-separated-missing-methods>
```

That manifest contains only
`Kidney_KPMP_full<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected`.
The recorded `<comma-separated-missing-methods>` is the missing subset of
`prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`, and the
expected method-row count equals its length. If its length is zero, record
`NOOP_VALIDATED` and submit no Kidney method array. These two gates are
serialized because both own `batch_effect/uncorrected_final`.


Add an explicit, idempotent metadata-only export for the four regenerated
H5ADs and `Kidney_KPMP_full`; current Stage 5 workers do not emit the
notebook’s sample-metadata Feather. Add
`src/utils/py/export_h5ad_sample_metadata.py`, using the existing h5py-only
reader in `src/utils/py/h5ad_pseudobulk.py`, and add the separate
read-only worker `src/utils/bash/h5ad_obs_audit_worker.sh`. The worker must
never call `ecoda_write_artifact_record` against the inspected H5AD.

For final Stage 5 rows, the worker reads only `obs` from the selected
uncorrected H5AD and atomically writes:

```text
${ANALYSIS_ROOT}/metadata/<DS>_sample_metadata.feather
${ANALYSIS_ROOT}/metadata/<DS>_sample_metadata.feather.md5
```

The export includes `Sample`, the configured primary biological label,
configured batch keys, configured cell-type columns, and every candidate
column required by `dataset_specs.py` for the final registry. It validates
non-empty unique sample IDs, preserves source sample order, and never opens
`X`, `raw`, or `layers["counts"]`. The final Stage 5 wrapper records the
exporter’s run-owned manifest, checksum, and terminal status before method
no-op selection; an already valid final metadata Feather is skipped
individually. This covers the four changed H5ADs and the existing Kidney
H5AD without copying any H5AD to the workstation.

Run a separate targeted Kidney Stage 5 recovery against the same final
analysis root. Before submission, inspect only the existing
`Kidney_KPMP_full/batch_effect_uncorrected` method artifacts and emit a
run-owned method manifest containing the missing subset of the seven-method
suite. Existing valid Kidney rows remain outside the recovery selection. If
no method is missing, write a no-op report and submit no Kidney method job;
the required obs-only metadata export still runs or validates independently.
If methods are missing, submit only those method rows with `--target-methods`
and a specific dependency reason; never force or rerun the full suite.

Both Stage 5 runs use the existing snapshot-backed durable workflow, exact
selection files, and current gate policy. They must be serialized when they
share the final analysis root. The four frozen cohorts and all disabled
cohorts are absent from both runs. Stage 5 runs only the uncorrected pass; the
corrected outputs required here are Stage 3 inputs for future use, not a
request to run corrected Stage 5 methods now.

### 8. Build the mixed-source final analysis manifests

Create two local, final-lane manifests only after Stage 5 artifacts and the
obs-only metadata export are available. They are analysis inputs, not
scheduler selections.

`data/batch_effect/uncorrected_final/final_analysis_metadata.tsv`:

```text
dataset<TAB>summary_lane<TAB>feather_lane<TAB>metadata_summary_path<TAB>metadata_feather_path
```

Write exactly one row for each final dataset in this order:

```text
Alzheimer
Breast_cancer
Covid19_PBMC
Kidney_KPMP_full
Diabetes
Lupus_PBMC
Lung
Joanito
Stephenson
```

For `Alzheimer`, `Breast_cancer`, `Lupus_PBMC`, and `Stephenson`, both paths
are the existing local legacy files:

```text
data/batch_effect/uncorrected/results/<DS>_batch_effect_uncorrected_metadata.rds
data/batch_effect/uncorrected/metadata/<DS>_sample_metadata.feather
```

For `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`, both paths use the
final root and exact stems:

```text
data/batch_effect/uncorrected_final/results/<DS>_batch_effect_uncorrected_final_metadata.rds
data/batch_effect/uncorrected_final/metadata/<DS>_sample_metadata.feather
```

For `Kidney_KPMP_full`, the Feather path is always the final obs-only export:

```text
data/batch_effect/uncorrected_final/metadata/Kidney_KPMP_full_sample_metadata.feather
```

Its summary path is the final
`<DS>_batch_effect_uncorrected_final_metadata.rds` if a final composition
artifact was produced and validated; otherwise it is the valid existing
legacy `data/batch_effect/uncorrected/results/<DS>_batch_effect_uncorrected_metadata.rds`.
Record `summary_lane` and `feather_lane` independently so this mixed source
is explicit rather than inferred.

`data/batch_effect/uncorrected_final/final_analysis_artifacts.tsv`:

```text
dataset<TAB>method<TAB>lane<TAB>artifact_kind<TAB>artifact_path<TAB>bundle_key
```

Each dataset has exactly these eight method keys:

```text
ECODA_authors_HR
ECODA_seuratres_2
Pseudobulk_hvg2000
GloScope_hvg2000_pcadims30
MrVI_hvg2000
PILOT_hvg2000
QOT_hvg2000
ECODA_authors_HR_NULL
```

For every legacy-lane row, use the existing root
`data/batch_effect/uncorrected`. For each `<DS>`, the exact
`method<TAB>artifact_kind<TAB>relative_path<TAB>bundle_key` mapping is:

```text
ECODA_authors_HR<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_composition.rds<TAB>ECODA_authors_HR
ECODA_seuratres_2<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_composition.rds<TAB>ECODA_seuratres_2
Pseudobulk_hvg2000<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_pseudobulk.rds<TAB>Pseudobulk_hvg2000
GloScope_hvg2000_pcadims30<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_gloscope.rds<TAB>GloScope_hvg2000_pcadims30
MrVI_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather<TAB>
PILOT_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather<TAB>
QOT_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_hvg2000_highres_qot_dists.feather<TAB>
ECODA_authors_HR_NULL<TAB>standalone_scores<TAB>results/<DS>_batch_effect_uncorrected_ECODA_authors_HR_NULL.rds<TAB>scores
```

The four frozen cohorts use those legacy rows read-only. A valid existing
Kidney method row also uses the legacy mapping.

For every final-lane row, use `data/batch_effect/uncorrected_final` and the
same method-relative layout:

```text
ECODA_authors_HR<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_final_composition.rds<TAB>ECODA_authors_HR
ECODA_seuratres_2<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_final_composition.rds<TAB>ECODA_seuratres_2
Pseudobulk_hvg2000<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_final_pseudobulk.rds<TAB>Pseudobulk_hvg2000
GloScope_hvg2000_pcadims30<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_final_gloscope.rds<TAB>GloScope_hvg2000_pcadims30
MrVI_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather<TAB>
PILOT_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_pilot_dists.feather<TAB>
QOT_hvg2000<TAB>distance_feather<TAB>embeddings/<DS>_batch_effect_uncorrected_final_hvg2000_highres_qot_dists.feather<TAB>
ECODA_authors_HR_NULL<TAB>rds_bundle<TAB>results/<DS>_batch_effect_uncorrected_final_composition.rds<TAB>ECODA_authors_HR_NULL
```
The six-column TSV always includes the physical sixth field. The
`bundle_key` field is explicitly empty (a trailing `<TAB>`) for
`distance_feather` rows, and the loader rejects any non-empty key for those
rows. This prevents physical-file counting from being confused with the
eight logical method keys.

For `Kidney_KPMP_full`, choose `lane=legacy` and the legacy path for each
method whose current artifact contract is valid; choose `lane=final` and the
final path for each method emitted by the targeted recovery. The final
composition bundle supplies the final null key; a legacy null row uses the
standalone `scores` file above. Never infer a path from a single common root.

Modify `src/utils/batch_effect_analysis.R` with a manifest-aware loader while
leaving `load_batch_uncorrected_dataset()` available for the legacy
contingency notebook. Add helpers such as:

```r
read_batch_final_manifest(path, expected_datasets, expected_methods)
load_batch_uncorrected_dataset_from_manifest(
  metadata_manifest, artifact_manifest, dataset, registry
)
```

For `lane=final` rows, validate headers/schema, exact nine-dataset order, one
metadata row per dataset, exactly eight method keys per dataset, valid
lane/artifact-kind values, and path containment under the repository root.
Validate final metadata-summary labels/sample counts against the final
metadata Feather sample IDs. For final `rds_bundle` artifacts, deserialize
the explicit `bundle_key` and require the bundle’s sample IDs to match the
dataset metadata. For final `distance_feather` artifacts, validate the
Feather schema/checksum and require its distance sample IDs to match. For
final `standalone_scores` artifacts, require the explicit `scores` key and
required score fields but do not require sample IDs: the legacy null artifact
is an aggregate score bundle without a sample-ID axis. The final composition
null row is instead an `rds_bundle` and must use its explicit
`ECODA_authors_HR_NULL` key.

Support the explicit `rds_bundle`/`standalone_scores`/`distance_feather`
mappings above, including the shared final composition/null path.
`lane=legacy` rows are approved read-only inputs for the four frozen cohorts
and valid Kidney rows: parse their explicit paths and bundle keys and read
them for the notebook, but do not run a separate preflight, checksum audit,
artifact validator, or eligibility/reuse check on those rows. An unreadable
legacy file is a notebook input error, not a reason to schedule a job. Load
every artifact from its manifest path and preserve existing RDS/Feather
distance and score semantics.

Generalize `batch_candidate_registry()` so an explicit non-empty unique
dataset order is accepted. Preserve the exact twelve-dataset validation when
the historical order is supplied, but for the final nine-dataset order require
that every selected dataset is covered by `DATASET_SPECS` and
`BATCH_EFFECT_SPECS`, preserving the caller’s order. Do not alter the
historical `dataset_specs.py` order or its diagnostic twelve-row contract.
The final notebook uses the nine-row registry; the contingency notebook keeps
its legacy twelve-row call.

### 9. Synchronize only final Stage 5 outputs to the workstation

After each final Stage 5 gate reaches terminal completion and review, sync
only the explicit final result artifacts, final metadata exports, checksums,
and final manifests:

```text
data/batch_effect/uncorrected_final/
  results/
  embeddings/
  pseudobulks/
  metadata/
  final_analysis_metadata.tsv
  final_analysis_artifacts.tsv
```

Generate the HPC `rsync --files-from` list from the final Stage 5 method
manifest plus the separate metadata-export manifest. Include final outputs
for the four changed datasets, any newly produced Kidney method/cache/result
rows, and the four regenerated-target/Kidney sample-metadata Feathers. Include
the exact `.md5` sidecars. Do not copy full H5ADs, raw counts, annotation
unions, or legacy frozen result files. The frozen rows in the local mixed
manifest point to their existing legacy paths.

Resolve `BAMBOO_HOME` with `ssh bamboo 'printf %s "$HOME"'` and use an
explicit `rsync --files-from` list rooted at the HPC scratch tree. Preserve
relative paths and do not use a recursive all-dataset sync. Verify the local
manifest paths, checksums, and file sizes against the terminal HPC manifest
before running the notebook.

### 10. Run the final uncorrected notebook only against the final lane

Update `notebooks/batch_effect_analysis_uncorrected.rmd` to:
- replace the opening “twelve canonical batch-effect cohorts” scope text with
  the nine-dataset mixed-source final-lane contract;
- replace every hard-coded twelve-dataset order/cardinality assertion,
  including the `Kidney_KPMP` legacy key, with the manifest’s exact nine-row
  order and `Kidney_KPMP_full`; leave those assertions in the untouched
  contingency notebook unchanged;

- use the exact nine-dataset order from `final_analysis_metadata.tsv`;
- read `final_analysis_metadata.tsv` and `final_analysis_artifacts.tsv`;
- include all nine final plot rows: frozen cohorts use only their approved
  legacy result/metadata artifacts read-only, while changed/Kidney rows use
  final artifacts;
- never read legacy dataset H5ADs for final plots;
- never reconstruct artifact paths from a common root;
- write only under `data/batch_effect/uncorrected_final/analysis/` and
  `plots/batch_effect_uncorrected_final/`;
- never read or write the legacy contingency output root for final outputs;
- retain the existing seven displayed methods plus the persisted shuffled
  baseline, using the explicit composition/null bundle mapping.

Expected minimum final outputs:

```text
plots/batch_effect_uncorrected_final/batch_effect_uncorrected_final_funkyheatmap.pdf
plots/batch_effect_uncorrected_final/<dataset>_mds.pdf
plots/batch_effect_uncorrected_final/<dataset>_anosim.pdf
data/batch_effect/uncorrected_final/analysis/batch_effect_uncorrected_final_scores.csv
data/batch_effect/uncorrected_final/analysis/batch_effect_uncorrected_final_decomposition.csv
data/batch_effect/uncorrected_final/analysis/batch_effect_uncorrected_final_nmi.csv
```

Retain the notebook’s existing per-dataset PERMANOVA/decomposition/NMI
outputs under the final plot root and use final-qualified aggregate table
names. Do not edit or run
`notebooks/batch_effect_analysis_uncorrected_batchconfounding_contingency.rmd`;
it remains legacy-only. Do not edit or run
`notebooks/batch_effect_analysis_legacy.rmd`.

Execute only the affected notebook chunks in one persistent Pixi R session,
not a full knitr render:

```text
hub start name=final-batch-r application=pixi \
  args=["run","-e","default","R"] cwd=<repository-root> pty=true ready.log=">"
```

Load the notebook setup, manifest-loading, analysis, and output chunks once in
dependency order using only the synced final result/metadata files. Do not
run unrelated notebook chunks. Verify the final plots/tables and that legacy
plot/analysis files retain their pre-run modification state.
## Critical files & anchors

- `src/utils/py/preprocess_utils.py:353-383` — shared subset evaluator and
  comparison/include-values/sample-consistency helpers.
- `src/3_scrnaseq_preprocessing/1.0_audit_input_views.py`,
  `src/utils/bash/h5ad_obs_audit_worker.sh`, and
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:618-859` — direct-H5AD
  obs-only Covid preflight, no-record audit worker, pre-subset processing
  order, final H5AD output, and embedding keys.
- `datasets.json:113-145,494-525,554-644,674-729` — Joanito, Covid, Kidney,
  Diabetes, and Lung source/view/filter/output contracts.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh`,
  `benchmark_submit_common.sh`, `ecoda_run_common.sh`, and the R/Python
  worker/preparation scripts — Stage 5 variant roots, provenance, paths,
  metadata export, worker dispatch, and synchronization.
- `src/utils/batch_effect_analysis.R:23-413` and
  `notebooks/batch_effect_analysis_uncorrected.rmd` — legacy method specs/
  loader to preserve plus the manifest-aware final loader and outputs.

## Verification

Run verification in this order and stop before the next stage on failure:

1. **Subset contract:** run
   `pixi run -e default python tests/test_subset_vars.py`. Observable proof:
   categorical strings `29`, `30`, `30.5`, `control`, `unknown`, and
   malformed input produce exactly `c0,c1,c3`; the separate mixed sample
   raises; repeated `PatientID` values with distinct `sampleID` values do not
   raise; scalar membership, invalid operators, non-finite values, and
   `copy=True/False` behave as specified.
2. **Configuration and routing contract:** parse `datasets.json`; assert both
   Covid views contain `op="<="`, scalar `values=30`, and
   `include_values=["control"]`; assert the eight final H5AD names, unchanged
   frozen/disabled names, final-nine registry order, and historical twelve-row
   registry behavior. Assert the Stage 2 selector, eight Stage 3 rows, four
   changed Stage 5 rows, and one Kidney Stage 5 dataset row contain no frozen,
   disabled, or `_debug` dataset. Run the focused registry/routing test.
3. **Focused source-contract tests:** run these before any gate:

   ```bash
   bash tests/test_benchmark_selection_file.sh
   bash tests/test_benchmark_sync.sh
   pixi run -e default Rscript --vanilla tests/test_benchmark_rds_contract.R
   pixi run -e default Rscript --vanilla tests/test_batch_effect_analysis.R
   ```

   These tests must prove legacy root/stem behavior remains unchanged, final
   variant roots and direct pseudobulk result paths are selected, the
   composition and null bundle keys are explicit, distance rows carry an
   empty `bundle_key`, and standalone null scores do not require sample IDs.
   Do not replace these focused checks with a project-wide suite.
4. **HPC Covid obs preflight:** before any Stage 3 worker compute, inspect the
   run-owned report from `h5ad_obs_audit_worker.sh`. Confirm it references the
   current HPC H5AD and source/runtime manifests, has valid checksum/size,
   records actual sampling-day dtype and raw levels, includes `control` and
   unknown/blank/malformed counts as observed, evaluates the exact predicate,
   and reports retained/dropped cells/samples with zero split samples. Confirm
   no H5AD-side artifact record or source mutation was written. A missing or
   mismatched report blocks Stage 3.
5. **Stage 2:** terminal evidence must show only the explicit Joanito hook,
   with no accidental other dataset/step rows. If valid `seqtec` and
   `cell.type_new` already exist, the run must be `NOOP_VALIDATED`, not a
   reprocessing job.
6. **Stage 3:** for each of the eight target views, confirm the final H5AD
   exists on HPC with valid checksum and required configured sample/label/
   cell-type columns; subset audits have zero split samples; uncorrected and
   corrected embedding keys match the semantic view; and the Covid output
   counts agree with the accepted HPC preflight. Confirm no frozen/disabled
   H5AD was selected or modified and no Pipeline 4 artifacts were created.
7. **Stage 5 changed wave:** expected selection is four rows and the fixed
   seven-method suite; expected method-row count is 28. Confirm final run
   metadata contains `ANALYSIS_VARIANT=final`, `PASS=uncorrected`,
   `ANALYSIS_PASS=uncorrected`, `ANALYSIS_ROOT`/`ROOT` ending in
   `batch_effect/uncorrected_final`, and no legacy root. Confirm the separate
   sample-metadata export exists and every final cache/result/embedding path
   has the exact `_final` stem and checksum.
8. **Stage 5 Kidney recovery:** inspect the current legacy artifact inventory
   immediately before submission. The emitted method manifest contains only
   missing Kidney rows; valid legacy rows are not selected. If the missing set
   is empty, proof is a no-op report with no Kidney method job. Otherwise
   every submitted row is one of the seven methods, the final metadata export
   is present, and no valid old Kidney row appears.
9. **Local sync:** compare the explicit sync list with terminal final method
   and metadata-export manifests. Confirm result bundles map to the actual
   pseudobulk result path (not its cache), all required `.md5` sidecars are
   present, no H5AD/raw-count/legacy frozen file is copied, and the four
   frozen rows remain references to existing local legacy artifacts.
10. **Final notebook:** execute only the final notebook chunks and confirm
    non-empty final funky heatmap, per-dataset MDS and ANOSIM PDFs for all nine
    rows, and final scores/decomposition/NMI tables. Confirm no legacy dataset
    H5AD was read, all paths came from manifests, `corrected_final` was not
    created as a Stage 5 lane, and legacy plot/analysis files retain their
    pre-run modification state.

All full-cohort Stage 2/3/5 launches and their required preflight workers must
use the checked-in snapshot-backed `durable-hpc-gate-ecoda` workflow. After
each launch, arm exactly one unbounded durable wait, perform one terminal
inspect with every emitted scheduler/watchdog ID, and obtain the required
reviewer approval before the dependent stage. These checks are scoped to the
selected rows; never add frozen or disabled cohorts to satisfy a historical
matrix.

## Assumptions and contingencies

- HPC scratch/NAS is authoritative. A local subset mirror, including
  `Covid19_PBMC_meta.json`, never qualifies as a production source or as
  evidence that a target H5AD is complete. The current HPC obs-only Covid
  report is the required pre-compute evidence.
- The four frozen cohorts are complete by user declaration. Do not schedule,
  preflight, validate, or recompute them. If a final notebook cannot find a
  declared local legacy result, stop with a missing-input report rather than
  adding that cohort to a job.
- The Covid obs preflight is a distinct read-only path. Never substitute the
  standard `h5ad_preflight_worker.sh`, which may publish run-owned artifact
  records. A preflight failure stops the run before Stage 3 processing.
- `Joanito` is the only required Stage 2 target hook. If valid `seqtec` and
  `cell.type_new` metadata already exist, Stage 2 emits no compute; otherwise
  it runs only the Joanito hook.
- Covid, Diabetes, and Lung direct H5AD inputs are already staged on HPC. If a
  direct input is missing, stop before submission; do not invent a Stage 2
  conversion hook.
- The existing Kidney uncorrected H5AD is present and is reused read-only. If
  it is absent, do not silently add Stage 3/4; stop and report the missing
  prerequisite because the approved scope is Stage 5-only.
- Kidney’s missing method set is determined from the current legacy artifact
  inventory immediately before its targeted Stage 5 gate. Valid rows are
  skipped individually; only missing rows are submitted. Its obs-only sample
  metadata export is independent of the method-row decision.
- Existing final target paths are skipped only when their current contract is
  valid. An invalid selected final path is recomputed as that same targeted
  row, never through a broad `--force` selection.
- `corrected_final` is a reserved Stage 3 namespace only. No corrected Stage 5
  root or artifact may be created by this run.
- If any emitted selection contains `Alzheimer`, `Breast_cancer`,
  `Lupus_PBMC`, or `Stephenson`, cancel the emitted scheduler IDs and durable
  runner, preserve evidence, and mark the run failed; do not let the
  unintended wave finish.
- If the submitted wrapper emits any dataset, view, method, or output path
  outside the written manifests, stop immediately and inspect the run as a
  scope mismatch.
- `PILOT-GM-VAE`, MOFA, scITD, scPoli, ordinary benchmark views, corrected
  Stage 5, and Pipeline 4 annotation work are not part of this plan.
- Legacy analysis outputs remain untouched. Final notebook plots/tables and
  final Stage 5 artifacts are separate even when a frozen cohort contributes
  a legacy result path to the mixed-source final manifest.

## Implementation status and findings

Updated during execution on 2026-09-12.

### Completed implementation

- Wave 1 added the shared fail-closed subset evaluator and sample-consistency audit, wired preprocessing and input-view audit callers, added strict direct-H5AD Covid obs-only auditing, added the immutable no-record obs audit worker, and added exact final Stage 3 target selection/release handling.
- Added `tests/test_subset_vars.py` with the categorical Covid boundary/control/malformed-value fixture, split-sample failure, repeated-PatientID non-failure, scalar membership, invalid-rule, nonfinite, and copy-semantics assertions.
- Configuration now contains the approved Covid `<= 30` plus case-insensitive `control` rule in both batch views and only the eight approved regenerated target H5AD names with `_final`. Stale batch scope policy and registry assertions were updated; historical/frozen/disabled contracts remain explicit.
- Stage 5 gained the explicit final uncorrected variant, final roots/stems, run metadata provenance, obs-only sample-metadata Feather export, worker propagation, final RDS/Feather validation, and final synchronization path expansion while retaining legacy mode.
- Final analysis gained explicit mixed-source nine-dataset manifests/loaders, artifact-kind and bundle-key validation, final composition/null sharing, legacy read-only lanes, registry generalization, and a final-only uncorrected notebook.
- Implementation commit: `e9ee50add76c2e7826980d7333e6f9440d5c647b`; pushed to `origin/master`; unrelated pre-existing files remain unstaged.

### Verification completed

- `pixi run -e default python tests/test_subset_vars.py` passed.
- `pixi run -e default python tests/test_batch_effect_registry_and_modes.py` passed after pinning the corrected-null guard fixture to a deliberately missing batch column.
- `bash tests/test_benchmark_selection_file.sh` passed.
- `bash tests/test_benchmark_sync.sh` passed after restoring `SYNC_ARTIFACTS=()` at the shared sync helper boundary.
- `pixi run -e default Rscript --vanilla tests/test_benchmark_rds_contract.R` passed.
- `pixi run -e default Rscript --vanilla tests/test_batch_effect_analysis.R` passed after hardening manifest path traversal rejection and the physical-field test fixture.
- `bash tests/test_preprocessing_stage_submitter.sh` passed with its expected stubbed ownership-error output.
- Integrated shell syntax, Python compilation, R parsing, and `datasets.json` JSON validation passed.

### Findings and repairs

- Wave 1 audit edits required structural repair of `read_obs_only()`, restoration of repository-root `sys.path` insertion, MD5 fields in Covid report identities, run-owned report/sidecar validation, target-only preflight gating, row-isolated worker dispatch, and scheduler-ID persistence.
- Final analysis path containment now rejects explicit `..` path components even when the target does not yet exist; this closes `normalizePath(mustWork=FALSE)` lexical traversal ambiguity.
- First durable Stage 2 attempt failed before scheduler submission because the snapshot executor correctly rejected the textual symlink `${BAMBOO_HOME}/scratch/ECODA_paper`; the failure evidence was preserved and the no-ID terminal inspect failed closed as required. The retried gate uses canonical `/srv/beegfs/scratch/users/h/halterc/ECODA_paper` scratch/log paths while retaining allowed gate manifest/status paths.

### Current gate status

- Immutable source snapshots exist for implementation commit e9ee50add76c2e7826980d7333e6f9440d5c647b, gate commit fe380b880bedad958d0e1929a2565a0bc7e3e2fe, and the spool-recovery commit 70c81fdfb33c651be5dcd51789211c21650810e5 under `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/`.
- Reviewed relocated runtime identity is `ecoda-py-cuda13-6bbf70b-relocated`; its Pixi TOML/lock hashes match the implementation source.
- Corrected durable Stage 2 gate `stage2_joanito_final_20260912b` completed with scheduler IDs 4403663/4403664; its exact run-scoped audit, terminal inspect, and Luna Max reviewer approval passed.
- No Stage 3 gate is currently active or release-eligible. Gate `stage3_batch_final_20260912c` failed in its Covid obs-only preflight with ID `4403668`; terminal inspection completed with failed accounting/release, reviewer approval was not performed, and retry remains paused.
- Frozen cohorts are not selected for new jobs; disabled cohorts and `_debug` remain outside production selections. No H5AD was copied to the workstation and no legacy artifact/output was intentionally overwritten.


### Live gate update

- The first Stage 2 gate (`stage2_joanito_final_20260912`) failed before scheduler submission because the snapshot executor rejected the textual symlink scratch root; its evidence remains preserved and no scheduler IDs were emitted.
- The corrected Stage 2 gate (`stage2_joanito_final_20260912b`) completed with scheduler IDs 4403663/4403664; the exact run-scoped audit, one terminal inspect, and Luna Max reviewer approval passed.
- The first Stage 3 prepared manifest was abandoned before launch after detecting a duplicated `scratch` path. The corrected gate `stage3_batch_final_20260912b` then failed before preprocessing because its Slurm spool copy could not resolve `../../slurm_config.sh`; the targeted `stage3_batch_final_20260912c` retry also failed before preprocessing in the container/spool bootstrap path. No Stage 3 worker rows were released.
- Stage 3 c terminal accounting/inspection completed as failed for preflight ID `4403668`; reviewer approval was not performed. The required container guard repair, four-row uncorrected retry, independent corrected wave, Stage 5 recoveries, final artifact sync, and final notebook execution remain pending.

### Live gate update 2

- Stage 3 gate `stage3_batch_final_20260912b` failed before preprocessing because the Slurm spool copy could not resolve `../../slurm_config.sh`; its one preflight array ID 4403666 was terminally inspected as failed and no Stage 3 worker rows were released.
- The bootstrap fix is committed as `70c81fdfb33c651be5dcd51789211c21650810e5` and has a verified immutable snapshot at `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/70c81fdfb33c651be5dcd51789211c21650810e5`.
- The Stage 3 recovery gate `stage3_batch_final_20260912c` is terminally failed in the Covid obs-only preflight with ID `4403668`; its terminal inspection completed with failed accounting/release, reviewer approval was not performed, and no downstream Stage 5 gate is authorized.
### User clarification superseding execution scope

The initial approved plan above is retained as implementation history. The
following clarification from the user on 2026-09-12 is now authoritative for
the remaining execution work wherever it conflicts with the earlier
frozen-cohort/four-plus-four Stage 3 wording.

1. Do not launch more gates, run more commands, or rebuild the HPC SIF while
   source contracts are still being repaired. Finish local code, configuration,
   manifest, worker, validator, and notebook stabilization first.
2. After those local/source/runtime contracts are proven, publish and validate
   one versioned final `py-cuda13` SIF/runtime identity for all subsequent
   waves. Do not repeatedly rebuild it after each failed worker bootstrap.
   Reuse the immutable SIF unless its identity is invalid or does not bind the
   finalized source/runtime contract.
3. The missing **uncorrected** Pipeline 3 subset is only
   `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`. Its next selection is
   four `batch_effect_uncorrected` rows, not the historical twelve-row matrix.
   No Alzheimer, Breast_cancer, Lupus_PBMC, Stephenson, or other completed
   cohort is selected for that uncorrected subset.
4. `Kidney_KPMP_full` is added separately in Pipeline 5 only. Its recovery
   selects only missing rows from
   `prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`;
   valid existing rows remain outside the selection. Kidney is not added to
   the uncorrected Pipeline 3 subset or Pipeline 4 annotation.
5. Corrected mode is independent of the uncorrected subset and may run in
   parallel after the local/runtime freeze. Its explicit Pipeline 3 selection
   must be generated from the current `datasets.json` at launch time:
   every non-underscore entry with `use_for_batch_effect == true`, using
   `batch_effect_corrected`. Do not carry forward the earlier frozen-cohort
   exclusion into this corrected wave, and do not infer membership from stale
   gate history. `_debug` remains diagnostic-only.
6. Corrected mode is a first full-cohort exercise of the recently added
   multi-batch integration. Treat failures as expected engineering risk:
   fail closed, preserve all evidence, and repair only failed dataset/view
   rows with a recorded dependency reason. A corrected failure must not trigger
   a broad rerun of successful corrected rows, the uncorrected subset, or
   frozen artifacts. Corrected Pipeline 3 output is separate from the final
   uncorrected Pipeline 5 analysis; no corrected Pipeline 5 lane is authorized
   by this clarification.

### Feasibility assessment

This clarified workflow is feasible. Pipeline 3 already accepts explicit
selection manifests, corrected and uncorrected views have distinct output
contracts, and the corrected wave can be isolated in its own durable run.
Pipeline 5 can remain serialized on the shared final uncorrected root while
the independent corrected Pipeline 3 run uses its own view/output namespace.
The parallelism constraint is explicit: the current durable profile uses one
`ecoda-benchmark` serialization group, and the custom Stage 3 Covid obs-only
preflight is intentionally triggered only by the exact four-dataset
uncorrected target selection. Therefore “parallel corrected mode” cannot be
implemented by inventing a second serialization group or by silently
bypassing the Covid release evidence. The safest supported design is one
validated Stage 3 scheduler manifest/wave containing the four uncorrected
target rows plus the dynamically generated corrected rows, while retaining
the exact target-row Covid preflight, or a separately implemented independent
gate with an explicit corrected-mode preflight/release contract. This choice
must be fixed in source and manifests before resuming; no ad hoc concurrent
gate is safe under the current shared lock.
The main feasibility risk is not the decomposition; it is first-run behavior
of the corrected multi-batch path on full cohorts: missing/constant batch
levels, unexpected metadata encodings, Harmony/resource failures, and
dataset-specific source columns may require targeted code or configuration
repairs.

Before the one final SIF publication and any resumed gate, complete this
stabilization checklist:

- Fix `h5ad_obs_audit_worker.sh` bootstrap recovery to use the established
  `SLURM_JOB_ID && ECODA_RUNTIME_IN_CONTAINER != 1` spool condition. Inside
  Apptainer, retain the inherited immutable source path and do not require
  `scontrol`. Gate `stage3_batch_final_20260912c` failed in its Covid
  preflight with ID `4403668`; terminal inspection completed as failed,
  reviewer approval was not performed, and retry remains paused.
- Make direct obs-only output-scope validation check absolute ancestry,
  run-root containment, and symlink safety before creating an external
  directory. Keep reports/checksums strictly under the bound run root.
- Route diagnostic onboarding subset evaluation through the shared evaluator
  (or remove the duplicate comparison implementation) so the scalar Covid
  threshold and comparison/include-values semantics cannot diverge.
- Bind obs-worker run roots to the global `${ECODA_RUNS_ROOT}/${RUN_ID}`
  layout before trusting manifests/reports. Revalidate final selection
  membership with parsed rows rather than newline counts.
- Enforce the final root identity in shared Stage 5 artifact/sync helpers,
  not only in the CLI, and verify final no-op reuse against the prior
  terminal global owner/producer record so valid final artifacts are skipped
  individually across runs.
- Make the final analysis loader bind `lane=final` to
  `data/batch_effect/uncorrected_final` and `lane=legacy` to the approved
  legacy lane, with frozen/Kidney lane policy explicit. Ensure the metadata
  exporter includes `BATCH_EFFECT_SPECS` candidates such as Joanito/Stephenson
  `Site` in addition to active `DATASET_SPECS` fields.
- Exercise the corrected path with validator-only metadata/source checks and
  deterministic small fixtures before full cohorts. Confirm biological labels
  remain evaluation-only covariates and that all configured corrected batch
  keys are present before allocating full jobs.

### Pause status for compaction

No further shell, SSH, test, durable-gate, scheduler, or SIF-build command is
authorized until the user completes clarification and compacts the main-agent
context. The initial implementation commits, failed-gate evidence, and local
plan history remain preserved. After compaction, resume with the checklist
above, then publish one final runtime identity, then launch the four-row
uncorrected Pipeline 3 subset, the independent all-configured corrected
Pipeline 3 wave, and the serialized targeted Pipeline 5 recoveries only after
their exact manifests and predecessor reviews are complete.
### Plan maintenance

All subsequent implementation, verification, gate, SIF, failure, repair, and
scope updates MUST be recorded in this plan. Keep each status entry concise:
date, status, evidence/path, and next blocker or action. Do not duplicate the
full rationale when a short evidence-linked update is sufficient.

### Resume update

- 2026-09-12, resumed after the user-authorized context compaction; local
  stabilization work is dispatched in disjoint units for the obs worker,
  obs-only/onboarding paths, Stage 5 root ownership, final analysis/export,
  and the combined Stage 3 wave. No new gate or SIF command is authorized
  until those units are integrated and the focused contracts pass.

- 2026-09-12, corrected-contract scout confirmed exactly nine dynamic corrected
  datasets and their configured sample/label/batch keys; evidence is the
  current `datasets.json` plus `src/utils/py/batch_contract.py`. The
  pre-allocation contract must inspect authoritative source `obs` for required
  batch columns, at least two levels, within-sample constancy, missing/sentinel
  values, near-unique levels, and full-rank composite design. Biological
  labels remain evaluation-only. No source-level counts can be inferred from
  static configuration.

- 2026-09-12, local verification caught a red combined Stage 3 fixture:
  `tests/test_preprocessing_stage_submitter.sh` reaches the validator-only
  combined run but its stubbed Covid preflight evidence is rejected. A
  focused test repair is dispatched; no HPC/SIF action is allowed until it
  passes.
- 2026-09-12, advisory-driven Stage 3 repairs now classify exact four-row
  uncorrected manifests independently (with the required Covid preflight),
  reserve combined mode for four-plus-dynamic-corrected rows, and fail closed
  when corrected RDS inputs lack a bound raw H5AD metadata cache. The
  corrected dynamic set remains the nine current config rows; missing RDS
  source metadata is now an explicit blocker rather than `CONFIG_ONLY_RDS`.

- 2026-09-12, reconciled `AGENTS.md` with the clarification: frozen cohorts
  remain excluded from uncorrected/final-lane work but are explicitly allowed
  in the independent corrected Stage 3 wave when selected by current
  `datasets.json`; disabled cohorts and `_debug` remain excluded everywhere.
  The corrected RDS prerequisite/source audit remains pending before runtime
  publication.

- 2026-09-12, stabilization contracts and the combined Stage 3 regression
  fixture pass the focused Python/R/shell checks, including the new 13-row
  combined selection test; commit `986c6c7` is pushed to `origin/master`.
  The next blocker is authoritative read-only corrected-source auditing on
  Bamboo before publishing/reusing the final runtime identity.

- 2026-09-12, incorporated the external read-only architecture review as
  evidence, not as authorization for a broad Ponytail refactor. The proposed
  deletion/DRY rewrite of Pipelines 2--5 is outside this execution scope and
  would risk the repository’s scientific and ownership invariants; preserve
  trust-boundary validation, atomic writes, checksums, and targeted recovery.
  The review did identify a real missing `read_datasets_json` import, now
  dispatched for correction, plus the required snapshot-bound RDS/source
  validation gaps.

### Pause status after structural review

- 2026-09-12, the user requested a pause before further plan execution:
  broader Pipeline 2--5 structural simplification/DRY/deletion work is
  explicitly out of scope for this session and belongs in a separate session.
  No runtime, gate, scheduler, HPC, SIF-publication, structural-refactor,
  deletion, or downstream Stage 5/analysis action is authorized here.
- 2026-09-12, the in-flight stabilization agents were cancelled after the
  read-only review. `986c6c7` remains the last pushed implementation commit.
  The attempted snapshot of that commit failed because the command used a
  short hash; no final `986c6c7` source snapshot or final runtime publication
  exists.
- 2026-09-12, the uncommitted worktree contains the plan update, the
  clarified `AGENTS.md` scope text, the missing `read_datasets_json` import,
  and prior obs-worker changes; unrelated pre-existing modifications remain
  untouched. The import correction and any partial cancelled-agent edits are
  not yet revalidated or committed.
- 2026-09-12, gate state is unchanged: Stage 2
  `stage2_joanito_final_20260912b` remains reviewed/completed with IDs
  `4403663/4403664`; Stage 3
  `stage3_batch_final_20260912c` remains failed in Covid preflight with ID
  `4403668`; no downstream Stage 5 gate is authorized. The focused suite was
  green before the latest unverified/cancelled edits; rerun it only after a
  future resume.
- 2026-09-12, next-session blockers are the pre-mkdir obs-worker path
  validation, snapshot-owned RDS `@meta.data` corrected-source validation,
  immutable-source binding for corrected preflight imports, strict
  corrected-only/partial-selection guards, and a fresh full-hash source
  snapshot/runtime freeze. Reassess the whole architecture before applying
  further patches.

### Structural review disposition

- 2026-09-12, the reviewer feedback is accepted as architecture triage, not
  as a pre-run refactor specification. Concrete trust-boundary blockers
  remain in this plan: verify/commit the `read_datasets_json` import,
  complete authoritative RDS `@meta.data` validation, finish obs-worker
  pre-mkdir ordering and strict Stage 3 selection guards, bind corrected
  preflight imports to the immutable snapshot, and reconcile the clarified
  `AGENTS.md` scope.
- 2026-09-12, deletion of stubs/shims, Stage 2 wrapper consolidation,
  submitter DRY rewrites, MD5 unification, and retirement of the R batch
  contract are explicitly deferred. They require a separate call-graph,
  schema-compatibility, migration, and focused-regression review. The
  separate follow-up plan is
  `.agents/plans/1789248891506-ecoda-pipeline-structure-plan.md`.
- 2026-09-12, reviewer line counts, commit identifiers, and checksum claims
  are not treated as authoritative while the worktree is moving; verify them
  independently before any structural change. The current execution remains
  paused and no runtime, gate, scheduler, HPC, SIF, or downstream analysis
  action is authorized in this session.

### Resume update after structural triage

- 2026-09-12, the user resumed implementation of this plan while keeping
  broad structural cleanup in the separate follow-up plan. Concrete fixes now
  include the authoritative `read_datasets_json` import, obs-worker
  pre-creation path validation, strict four-row/combined/corrected-only
  selection classification, immutable-source-bound corrected preflight, and
  snapshot-owned RDS `@meta.data` auditing.
- 2026-09-12, focused verification passed after these fixes:
  `test_subset_vars.py`, `test_batch_effect_registry_and_modes.py`,
  `test_multibatch_contracts.py`, `test_preprocessing_stage_submitter.sh`,
  `test_benchmark_selection_file.sh`, `test_benchmark_sync.sh`,
  `test_benchmark_rds_contract.R`, `test_batch_effect_analysis.R`, and
  `test_corrected_source_metadata.R`; shell syntax, Python compilation,
  `datasets.json` parsing, and R auditor parsing also passed. No scheduler,
  full-cohort worker, SIF, or durable-gate action has been taken since the
  pause; the subsequent read-only SSH/source inventory and failed RDS audit
  are recorded below.

- 2026-09-12, the H5AD corrected-subset audit patch was paused/cancelled
  while the user considered a Snakemake-based structural replacement; no
  canceled-agent output is accepted as verified. The production corrected
  wave remains blocked on the authoritative RDS audit resource failure.

- 2026-09-12, the resumed trust-boundary implementation and focused
  regression are committed and pushed as `332f7c4` (following `986c6c7`).
  The new snapshot-owned RDS auditor is covered by
  `tests/test_corrected_source_metadata.R`. The next action is read-only
  authoritative source auditing and then one full-hash snapshot/runtime
  freeze; no broad structural refactor is part of this execution.

- 2026-09-12, authoritative corrected-source audit was started against the
  run-owned audit root. The direct Joanito RDS metadata audit was killed with
  exit 137 after about 173 seconds while deserializing the 3.88-GB RDS on
  Bamboo; no report was produced and the sequential Stephenson audit did not
  start. This is a no-compute preflight resource failure, not validation
  evidence. The corrected wave remains unauthorized until metadata auditing
  moves to a memory-sufficient/bounded path without bypassing the RDS
  contract.

### H5AD corrected-source audit update

- 2026-09-12, the user explicitly selected continuation of this plan and
  dropped Snakemake from the current decision. A source-bound H5AD obs-only
  corrected-source auditor now mirrors the RDS subset/sample/batch contract;
  its deterministic regression passes. The earlier remote H5AD audit attempts
  correctly stopped because the existing `332f7c4` snapshot predated this
  helper; no source data was treated as validated.
- 2026-09-12, the existing Joanito/Stephenson raw H5AD caches were inspected
  read-only and their global producer records were absent. They remain
  diagnostic only, not authoritative replacements for the killed direct RDS
  audit. The direct Joanito RDS audit remains blocked by exit 137/OOM.
