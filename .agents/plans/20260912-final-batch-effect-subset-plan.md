# ECODA final batch-effect subset processing

## Context

Implement the approved final batch-effect analysis without rerunning completed
cohorts or broad historical selections. The production source of truth is the
current `datasets.json` plus the authoritative full-cohort data on Bamboo/HPC;
local subset mirrors are diagnostic only. Pipeline 3 uses separate parallel
uncorrected and corrected selections. The uncorrected lane feeds the final
uncorrected Stage 5 analysis, while the corrected lane feeds its own
`corrected_final` Stage 5 lane; both lanes reuse valid outputs idempotently and
write outside the legacy analysis lane.

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
sample-consistency audit, strict direct-H5AD obs-only preflight, and two
separate Stage 3 release gates: exactly four uncorrected rows and exactly nine
configured corrected rows.

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

Stage 2/3/5 launch preparation is serialized only across code/config
integration and documented data dependencies: Stage 2 Joanito precedes both
Stage 3 selections; each Stage 5 lane follows its own reviewed Stage 3
predecessor. The uncorrected Stage 5 final lane is one five-dataset wave:
its four regenerated datasets and `Kidney_KPMP_full` are selected together.
The reviewed uncorrected Stage 3 predecessor and the current four-row
`NOOP_VALIDATED` check are complete, so the uncorrected Stage 5 gate may be
prepared and launched now while corrected Stage 3 is still running. Once
corrected Stage 3 reaches terminal `COMPLETED` and Luna Max review, the
corrected Stage 5 gate may be prepared and launched without waiting for the
uncorrected Stage 5 gate to finish. The corrected Stage 5 lane is a distinct
nine-dataset root and serialization group; the two Stage 5 gates are
independent and may run concurrently. Valid Kidney artifacts are skipped
individually inside the uncorrected five-dataset wave; there is no separate
uncorrected Kidney gate or same-root gate serialization.


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
| Corrected Stage 5 lane | Latest user clarification | Run `--pass corrected --analysis-variant corrected_final` against the nine current configured corrected datasets at `batch_effect/corrected_final`, with seven methods, concurrent dataset rows, and an explicit eight-key artifact manifest; reuse valid rows idempotently. |
| Batch annotation scope | User-approved `AGENTS.md` exception | Batch-effect views do not run Pipeline 4. Preserve configured source/author cell-type columns. |
| Frozen cohorts | User decision | `Alzheimer`, `Breast_cancer`, `Lupus_PBMC`, and `Stephenson` are frozen and absent from new uncorrected/final-lane jobs, validator selections, and compute manifests; they remain eligible in the independently configured corrected Stage 3 and corrected Stage 5 nine-dataset lanes. Final notebook plots use their approved legacy result artifacts only; no legacy dataset H5AD is read for the final analysis. |
| Changed final-view targets | User decision and current registry | `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`, both uncorrected and corrected Stage 3 views. |
| Kidney coverage | User decision | `Kidney_KPMP_full` is included in the one uncorrected five-dataset Stage 5 wave; a run-owned legacy inventory skips valid rows individually and leaves only missing/invalid rows pending in that same wave. No separate Kidney gate is created. |
| Disabled cohorts | Current registry flags | `CombinedPBMC`, `Kidney_KPMP`, `Myocardial_infarction`, and `Parkinson` receive no new work. `_debug` is not a production target. |
| Stage 2 coverage | `src/2_dataset_specific_preprocessing/1_submit_hpc.sh:233-259,394-444` | Only `Joanito` has a target-specific Stage 2 hook. Covid, Diabetes, and Lung use already staged direct H5AD inputs; do not invent Stage 2 jobs for them. |
| Corrected batch variables | `datasets.json` target contracts | Every corrected Stage 3 row uses its configured technical batch variables directly; the current corrected examples include `Covid19_PBMC=datasets`, `Diabetes=dataset`, `Joanito=seqtec`, and `Lung=dataset`. Biological labels remain evaluation-only. |
| Final naming and roots | Latest user clarification | Uncorrected Stage 5 uses `_final` stems under `batch_effect/uncorrected_final`; corrected Stage 5 uses `_corrected_final` stems under `batch_effect/corrected_final`. Corrected Stage 3 H5ADs remain distinct inputs to the corrected Stage 5 lane. |
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

Pipeline 3 uses two separate selection manifests and gates. Dataset rows
dispatch concurrently within each selection; manifest order is deterministic
ordering only.

Uncorrected Stage 3 selection, exactly four rows:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_uncorrected
```

Corrected Stage 3 selection, exactly every current non-underscore
`datasets.json` entry with `use_for_batch_effect=true`, in config order:

```text
Joanito<TAB>batch_effect_corrected
Stephenson<TAB>batch_effect_corrected
Alzheimer<TAB>batch_effect_corrected
Breast_cancer<TAB>batch_effect_corrected
Covid19_PBMC<TAB>batch_effect_corrected
Kidney_KPMP_full<TAB>batch_effect_corrected
Diabetes<TAB>batch_effect_corrected
Lupus_PBMC<TAB>batch_effect_corrected
Lung<TAB>batch_effect_corrected
```

Uncorrected Stage 5 final selection, exactly five datasets:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
Kidney_KPMP_full<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected
```

Corrected Stage 5 final selection, exactly the nine corrected datasets above,
in the same order:

```text
Joanito<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Stephenson<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Alzheimer<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Breast_cancer<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Covid19_PBMC<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Kidney_KPMP_full<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Diabetes<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Lupus_PBMC<TAB>batch_effect_corrected<TAB>batch_effect_corrected
Lung<TAB>batch_effect_corrected<TAB>batch_effect_corrected
```

The four frozen cohorts, `_debug`, and all disabled cohorts must be absent from
new uncorrected jobs and validator inputs. Frozen cohorts may occur in the
independent corrected selections above only when present in the current
configuration rule. Do not pass `--exact-batch-selection`, because that mode
requires the obsolete historical twelve-row matrix.

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
Before full-cohort launch, `AGENTS.md` records the approved current scope:
Pipeline 3 uses separate parallel selections, exactly four uncorrected target
rows and exactly nine current configured corrected rows; dataset rows dispatch
concurrently and manifest order is deterministic only. The corrected Pipeline 3
path uses original cell-level technical metadata directly, with no
sample-level constancy/majority check and no retired corrected-source RDS
preflight; configuration/content/provenance checks and the Covid obs-only
subset preflight remain. Stage 5 has separate uncorrected and corrected lanes:
the uncorrected five-dataset selection
(`Covid19_PBMC`, `Diabetes`, `Joanito`, `Lung`, and `Kidney_KPMP_full`) uses
`batch_effect/uncorrected_final`, while the corrected nine-dataset lane uses
`--pass corrected --analysis-variant corrected_final` and
`batch_effect/corrected_final`, the seven methods, and an explicit eight-key
artifact manifest. Valid rows are reused idempotently. Preserve the existing
durable-gate, snapshot, and separate structural-cleanup-plan boundaries while
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

### 6. Regenerate the separate Stage 3 selections

Run two separate Stage 3 manifests and durable gates. The uncorrected manifest
contains exactly the four rows in step 1; the corrected manifest contains
exactly the nine current configured rows in step 1. Dataset rows dispatch
concurrently within each gate. Do not combine the selections or serialize the
four and nine rows behind one matrix.

Uncorrected Stage 3:

```text
Covid19_PBMC<TAB>batch_effect_uncorrected
Diabetes<TAB>batch_effect_uncorrected
Joanito<TAB>batch_effect_uncorrected
Lung<TAB>batch_effect_uncorrected
```

Corrected Stage 3:

```text
Joanito<TAB>batch_effect_corrected
Stephenson<TAB>batch_effect_corrected
Alzheimer<TAB>batch_effect_corrected
Breast_cancer<TAB>batch_effect_corrected
Covid19_PBMC<TAB>batch_effect_corrected
Kidney_KPMP_full<TAB>batch_effect_corrected
Diabetes<TAB>batch_effect_corrected
Lupus_PBMC<TAB>batch_effect_corrected
Lung<TAB>batch_effect_corrected
```

Use the canonical Stage 3 submitter with each explicit `--selection-file`; do
not use the historical exact-selection mode. The uncorrected gate excludes
Kidney, frozen cohorts, disabled cohorts, and `_debug`; the corrected gate
includes the nine config-selected rows above. Before the Covid row in either
gate is released, require the corresponding read-only HPC `obs` preflight.
Those reports are run-owned evidence of the current source, not H5AD artifact
records, and must not write beside the immutable direct input. The retired
corrected-source metadata/RDS release preflight is not part of either gate.

The Stage 3 worker must resolve output names from the updated `datasets.json`.
The four uncorrected targets use their final-qualified names. Corrected rows
use their configured corrected-view names, including the corrected outputs for
the nine current config-selected datasets.

The uncorrected representation remains `Sample`-keyed raw
PCA/neighbors/Leiden without Harmony. The corrected representation passes each
cell's original configured technical batch metadata directly to Harmony/HVG,
with the biological label excluded from all processing covariates. Do not
perform a sample-level batch constancy check, within-`Sample` check, or
cell-level majority rewrite in Pipeline 3. Configuration, content, provenance,
and the Covid subset preflight remain required.

Preserve source/author cell-type metadata in every output. Do not call
Pipeline 4, prepare annotation chunks, run annotation workers, merge
annotation Feathers, or add HiTME/scATOMIC columns.

Each selected row is idempotent: validate non-empty output, schema, checksum,
and ownership before dispatch; reconcile any prior owner-state discrepancy
validator-only before reuse. A valid output from the prior reviewed gate is
skipped as `NOOP_VALIDATED`; only missing or invalid rows are recomputed, and
never through a broad `--force` selection.

After terminal review, record all selected output paths and their subset audit
summaries in the run-owned manifest. Each audit must include dataset, view,
configured raw sample column, total/retained/dropped cell counts,
total/retained/dropped sample counts, the exact HPC input identity, and
split-sample count (which must be zero).

### 7. Add a final variant to Stage 5 without changing legacy mode

Extend the canonical Stage 5 wrapper and shared artifact-path helpers with
explicit `--analysis-variant final` and
`--analysis-variant corrected_final` options. The default with no variant must
remain byte-for-byte compatible with the existing legacy roots and stems.
`final` is valid only for the explicit five-dataset uncorrected batch-effect
selection; `corrected_final` is valid only for the explicit nine-dataset
corrected selection. Reject either variant with `benchmark_analysis`, the
wrong pass, a broad default selection, or a missing explicit selection file.

Set `ANALYSIS_VARIANT`, `ANALYSIS_ROOT`, `ANALYSIS_NAS_ROOT`,
`ANALYSIS_PASS`, and `ANALYSIS_LOG_PREFIX` before constructing pending
selection state or writing run metadata. The two final run metadata contracts
must contain:

```text
ANALYSIS_VARIANT=final
ANALYSIS_ROOT=${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final
ANALYSIS_NAS_ROOT=${NAS_TARGET_DIR}/batch_effect/uncorrected_final
ANALYSIS_PASS=uncorrected
PASS=uncorrected
ROOT=${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final

ANALYSIS_VARIANT=corrected_final
ANALYSIS_ROOT=${HPC_SCRATCH_DIR}/batch_effect/corrected_final
ANALYSIS_NAS_ROOT=${NAS_TARGET_DIR}/batch_effect/corrected_final
ANALYSIS_PASS=corrected
PASS=corrected
ROOT=${HPC_SCRATCH_DIR}/batch_effect/corrected_final
```

Move or refactor the current root assignment before the `RUN_METADATA` block
around `src/5_run_benchmark_methods/1_submit_hpc_array.sh:2398-2425`, so the
metadata cannot record a legacy root for either final run.
`ecoda_run_audit.sh`, artifact ownership records, worker environments,
validator path reconstruction, watchdogs, and synchronization must consume the
same variant-qualified root. The legacy mode must continue to emit
`batch_effect/uncorrected` or `benchmark` roots exactly as before.

Update these exact layers together:

- `src/5_run_benchmark_methods/1_submit_hpc_array.sh`: parse/validate both
  options, establish the selected final root, export the variant, and record
  the fields above in both active and `NOOP_VALIDATED` metadata.
- `src/utils/bash/ecoda_run_common.sh:_ecoda_stage5_artifacts_for`: derive
  variant-qualified stems for ownership and sync expansion.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh:benchmark_artifacts_for`:
  use the same centralized stem rule for every Stage 5 method.
- `src/5_run_benchmark_methods/benchmark_submit_common.sh:878-945`:
  enumerate the same variant-qualified paths during sync and validation.
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:1857-1863`:
  add `_final` or `_corrected_final` to batch output names for the selected
  variant while keeping embedding keys semantic-view based.
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:337-342,419-423`:
  add the selected final suffix to batch cache/result stems. The result stem
  is the cache stem followed by the method name; it is not the pseudobulk
  cache filename.
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R`:
  accept the variant-qualified cache root and preserve legacy cache behavior.
- Both worker wrappers must forward the selected variant and use it in
  execution-log filenames.
- `src/5_run_benchmark_methods/validate_benchmark_rds_contract.R` and
  `src/utils/bash/ecoda_run_common.sh` must validate direct final paths and
  never reconstruct a final artifact under the legacy pass root.

Final Stage 5 roots:

```text
${HPC_SCRATCH_DIR}/batch_effect/uncorrected_final
${NAS_TARGET_DIR}/batch_effect/uncorrected_final
${HPC_SCRATCH_DIR}/batch_effect/corrected_final
${NAS_TARGET_DIR}/batch_effect/corrected_final
```

The two roots are independent lanes. The corrected lane is a Stage 5
benchmark lane, not a Stage 3-only namespace; its outputs consume corrected
Stage 3 H5ADs and remain separate from the uncorrected analysis lane.

Each final Stage 5 lane must distinguish the pseudobulk cache from the result
bundle. Uncorrected paths are:

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

Corrected paths are analogous, with the corrected suffix:

```text
${ANALYSIS_ROOT}/pseudobulks/<DS>_batch_effect_corrected_final_pseudobulk_hvg2000.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_corrected_final_pseudobulk.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_corrected_final_gloscope.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_corrected_final_composition.rds
${ANALYSIS_ROOT}/results/<DS>_batch_effect_corrected_final_metadata.rds
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather
${ANALYSIS_ROOT}/embeddings/<DS>_batch_effect_corrected_final_hvg2000_highres_qot_dists.feather
```

For either lane, the logical manifest has exactly these eight keys:
`ECODA_authors_HR`, `ECODA_seuratres_2`, `Pseudobulk_hvg2000`,
`GloScope_hvg2000_pcadims30`, `MrVI_hvg2000`, `PILOT_hvg2000`,
`QOT_hvg2000`, and `ECODA_authors_HR_NULL`. The composition path is shared by
`ECODA_authors_HR`/`ECODA_seuratres_2` and the null key, with explicit bundle
keys `ECODA_authors_HR`, `ECODA_seuratres_2`, and
`ECODA_authors_HR_NULL`. `Pseudobulk_hvg2000` maps to the pseudobulk result
bundle, not the cache. Distance rows use the corresponding Feather paths and
an explicitly empty `bundle_key`; the physical manifest always retains its
sixth field.

The final Stage 5 suite is exactly:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

Do not select `mofa`, `scitd`, `scpoli`, `pilotgm`, `trans`, `zeroimp`, or
any post-baseline method. The uncorrected final lane is one explicit
five-dataset selection in the order from step 1:
`Covid19_PBMC`, `Diabetes`, `Joanito`, `Lung`, and `Kidney_KPMP_full`. It
declares a maximum of five dataset rows times seven methods, 35 method rows;
the run-owned Kidney legacy inventory skips valid legacy artifacts
individually, so only missing or invalid Kidney methods remain pending in this
same wave. The corrected wave contains nine dataset rows and at most 63
method rows. Dataset rows dispatch concurrently within each lane, and valid
rows are skipped individually.
Before each Stage 5 gate, record the exact remote wrapper command and expected
selection/method-row ceiling in the durable-gate manifest. The uncorrected
final gate's `--exact-command` must be the snapshot
`ecoda_source_snapshot.sh exec` wrapper required by `AGENTS.md`, with this
script/argument tail:

```text
--script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- \
  --selection-file <run-root>/manifests/stage5_uncorrected_final.tsv \
  --pass uncorrected \
  --analysis-variant final \
  --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

The recorded selection file is exactly the five uncorrected dataset rows from
step 1. The wrapper's declared method-row ceiling is 35, including the
`prepare_pseudobulk` dependency row for each selected dataset; valid legacy
Kidney rows are omitted from pending/owner/validation/sync expansion by the
run-owned inventory. No `--target-methods`, `--analyses`,
`--exact-batch-selection`, broad dataset list, corrected pass, or `--force` is
permitted for this five-dataset wave.

Before the uncorrected gate is released, the same run-owned validator-only
legacy inventory records all seven Kidney methods. Valid legacy artifacts stay
outside the pending selection, while missing/invalid methods are submitted in
the same five-dataset wave. Do not construct a second `kidney_missing_final`
manifest or launch a separate Kidney recovery gate.

The corrected gate uses the same snapshot wrapper and this tail:

```text
--script src/5_run_benchmark_methods/1_submit_hpc_array.sh -- \
  --selection-file <run-root>/manifests/stage5_corrected_final.tsv \
  --pass corrected \
  --analysis-variant corrected_final \
  --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot
```

Its recorded selection is exactly the nine corrected rows from step 1 and its
maximum method-row count is 63. No broad dataset list, wrong pass, or `--force`
is permitted.

Both final Stage 5 lanes use the existing snapshot-backed durable workflow,
exact selection files, and current gate policy. The uncorrected five-dataset
gate may start now from its reviewed uncorrected Stage 3 predecessor. The
corrected nine-dataset gate may start as soon as its corrected Stage 3
predecessor reaches terminal `COMPLETED` and Luna Max review; it does not wait
for the uncorrected gate to finish. Their dataset rows and method rows dispatch
concurrently within each distinct root, and the two roots may be gated in
parallel. The uncorrected five-dataset wave owns the only uncorrected final
root; the corrected lane uses exactly the nine current config-selected
corrected datasets. Existing valid rows remain reusable and are never forced.


Add explicit, idempotent metadata-only exports for both final lanes; current
Stage 5 workers do not emit the notebook's sample-metadata Feather. Add
`src/utils/py/export_h5ad_sample_metadata.py`, using the existing h5py-only
reader in `src/utils/py/h5ad_pseudobulk.py`, and add the separate read-only
worker `src/utils/bash/h5ad_obs_audit_worker.sh`. The worker must never call
`ecoda_write_artifact_record` against the inspected H5AD.

For uncorrected final rows, the worker reads only `obs` from the four
regenerated final H5ADs and existing `Kidney_KPMP_full` uncorrected H5AD and
atomically writes:

```text
${ANALYSIS_ROOT}/metadata/<DS>_sample_metadata.feather
${ANALYSIS_ROOT}/metadata/<DS>_sample_metadata.feather.md5
```

For corrected final rows it reads only `obs` from the nine corrected Stage 3
H5ADs and writes the same metadata/checksum pair below the corrected root.
Corrected sample-level consumers read this explicit Feather path, verify its
MD5 sidecar and Sample order against the selected H5AD, and do not recompute
votes. The exporter applies majority only to
`Alzheimer/assay`, `Breast_cancer/suspension_dissociation_time` (literal
`unknown` is an ordinary configured class), and `Lupus_PBMC/batch_cov`.
Uncorrected exports perform no majority assignment; cell-level corrected
methods retain source cell values.

Each export includes `Sample`, the configured primary biological label,
configured batch keys, configured cell-type columns, and every candidate
column required by `dataset_specs.py` for the applicable final registry. It
validates non-empty unique sample IDs, preserves source sample order, and
never opens `X`, `raw`, or `layers["counts"]`. Each final Stage 5 wrapper
records the exporter's run-owned manifest, checksum, and terminal status
before method no-op selection; an already valid metadata Feather is skipped
individually.

For the uncorrected five-dataset wave, the run-owned Kidney legacy inventory
is computed and validated before pending method selection. Valid legacy
Kidney rows remain outside the pending/owner/validation/sync sets; missing or
invalid methods join the same uncorrected wave. The required uncorrected
obs-only metadata export still runs or validates independently. There is no
separate targeted Kidney Stage 5 recovery.

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

### 9. Synchronize only approved final Stage 5 outputs to the workstation

After each final Stage 5 gate reaches terminal completion and review, sync only
the explicit lane-specific result artifacts, metadata exports, checksums, and
run-owned manifests. The uncorrected lane is:

```text
data/batch_effect/uncorrected_final/
  results/
  embeddings/
  pseudobulks/
  metadata/
  final_analysis_metadata.tsv
  final_analysis_artifacts.tsv
```

The corrected lane is a separate synchronized root:

```text
data/batch_effect/corrected_final/
  results/
  embeddings/
  pseudobulks/
  metadata/
  final_analysis_artifacts.tsv
```

Generate each HPC `rsync --files-from` list from that lane's final Stage 5
method manifest plus its separate metadata-export manifest. Include uncorrected
outputs for the four changed datasets, any newly produced Kidney
method/cache/result rows, and the four regenerated-target/Kidney
sample-metadata Feathers. Include corrected outputs for the nine configured
corrected datasets and their metadata Feathers. Include the exact `.md5`
sidecars. Do not copy full H5ADs, raw counts, annotation unions, or legacy
frozen result files. The frozen rows in the uncorrected local mixed manifest
point to their existing legacy paths; corrected artifacts remain in their
separate lane and are not inferred from that manifest.

Resolve `BAMBOO_HOME` with `ssh bamboo 'printf %s "$HOME"'` and use an
explicit `rsync --files-from` list rooted at the HPC scratch tree. Preserve
relative paths and do not use a recursive all-dataset sync. Verify each local
lane's manifest paths, checksums, and file sizes against its terminal HPC
manifest before running the uncorrected notebook.

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
6. **Stage 3:** validate each of the four uncorrected and nine corrected
   selected outputs independently. Confirm non-empty files, required schema,
   checksum, ownership, configured sample/label/cell-type columns, and
   zero-split subset audits. Confirm uncorrected and corrected embedding keys
   match their semantic views, the Covid counts agree with the accepted
   obs-only preflight, and no Pipeline 4 artifacts were created. Frozen rows
   are excluded from the uncorrected gate but appear in corrected validation
   only when selected by the current configuration rule.
7. **Stage 5 final lanes:** the uncorrected final lane is one five-dataset
   selection and has a maximum of 35 declared method rows. Confirm its
   metadata contains `ANALYSIS_VARIANT=final`, `PASS=uncorrected`,
   `ANALYSIS_PASS=uncorrected`, and a root ending in
   `batch_effect/uncorrected_final`. Confirm the run-owned Kidney inventory
   proves valid legacy rows were skipped individually and only missing/invalid
   rows entered the same wave. The corrected wave has nine dataset rows and at
   most 63 method rows; confirm `ANALYSIS_VARIANT=corrected_final`,
   `PASS=corrected`, `ANALYSIS_PASS=corrected`, and a root ending in
   `batch_effect/corrected_final`. Both lanes must use exact `_final` or
   `_corrected_final` stems, explicit eight-key manifests, checksums, and
   individual valid-row reuse.
8. **Stage 5 Kidney coverage:** inspect the current legacy artifact inventory
   immediately before the uncorrected five-dataset gate. Valid legacy rows
   must remain outside the pending selection; missing/invalid rows must be
   submitted in that same wave. The metadata export is independent of the
   method-row decision. There is no separate Kidney method job or
   `--target-methods` recovery gate.
9. **Local sync:** compare the explicit sync list with terminal final method
   and metadata-export manifests. Confirm result bundles map to the actual
   pseudobulk result path (not its cache), all required `.md5` sidecars are
   present, no H5AD/raw-count/legacy frozen file is copied, and the four
   frozen rows remain references to existing local legacy artifacts.
10. **Final notebook:** execute only the uncorrected final notebook chunks and
    confirm non-empty final funky heatmap, per-dataset MDS and ANOSIM PDFs for
    all nine rows, and final scores/decomposition/NMI tables. Confirm no legacy
    dataset H5AD was read, all paths came from manifests, corrected Stage 5
    artifacts remain in their separate lane, and legacy plot/analysis files
    retain their pre-run modification state.

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
- The four frozen cohorts are complete by user declaration for the
  uncorrected/final analysis lane. Do not schedule, preflight, validate, or
  recompute them in that lane. They remain eligible for the independently
  configured corrected Pipeline 3 and corrected Pipeline 5 selections. If the
  final notebook cannot find a declared local legacy result, stop with a
  missing-input report rather than adding that cohort to an uncorrected job.
- The Covid obs preflight is a distinct read-only path. Never substitute the
  standard `h5ad_preflight_worker.sh`, which may publish run-owned artifact
  records. A preflight failure stops the run before Stage 3 processing.
- `Joanito` is the only required Stage 2 target hook. If valid `seqtec` and
  `cell.type_new` metadata already exist, Stage 2 emits no compute; otherwise
  it runs only the Joanito hook.
- Covid, Diabetes, and Lung direct H5AD inputs are already staged on HPC. If a
  direct input is missing, stop before submission; do not invent a Stage 2
  conversion hook.
- The existing Kidney uncorrected H5AD is present and is reused read-only for
  the uncorrected/final Stage 5 lane. If it is absent, do not silently add
  uncorrected Stage 3/4; stop and report the missing prerequisite because the
  approved uncorrected scope is Stage 5-only. The independent corrected
  nine-dataset lane follows its own corrected Stage 3 input contract.
- Kidney's legacy artifact inventory is determined immediately before the
  uncorrected five-dataset Stage 5 gate. Valid rows are skipped individually
  inside that wave; only missing or invalid methods are submitted there. Its
  obs-only sample metadata export is independent of the method-row decision.
- Existing final target paths are skipped only when their current contract is
  valid. An invalid selected final path is recomputed as that same targeted
  row, never through a broad `--force` selection. This idempotent contract
  applies independently to both final Stage 5 roots.
- `corrected_final` is the corrected Stage 5 root and variant. It is distinct
  from corrected Stage 3 H5AD/view paths and must not be reconstructed under
  the legacy pass root.
- If an emitted **uncorrected** selection contains `Alzheimer`,
  `Breast_cancer`, `Lupus_PBMC`, or `Stephenson`, cancel the emitted scheduler
  IDs and durable runner, preserve evidence, and mark the run failed; do not
  let the unintended wave finish. The corrected nine-dataset selection is
  governed by the current non-underscore, batch-enabled configuration rule.
- If the submitted wrapper emits any dataset, view, method, or output path
  outside the written manifests, stop immediately and inspect the run as a
  scope mismatch.
- `PILOT-GM-VAE`, MOFA, scITD, scPoli, ordinary benchmark views, and Pipeline 4
  annotation work are not part of these final Stage 5 lanes.
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
   frozen artifacts. Corrected Pipeline 3 output feeds the independent
   corrected Stage 5 lane as well as remaining downstream contracts; valid
   Stage 5 rows remain reusable.

### Feasibility assessment

This clarified workflow is feasible. Pipeline 3 accepts two explicit
selection manifests and distinct output contracts: exactly four uncorrected
rows and exactly nine corrected rows. Each gate dispatches its dataset rows in
parallel; manifest order is deterministic only. Pipeline 5 has distinct
uncorrected and corrected roots, so their independent five- and nine-dataset
gates may run in parallel after their respective Pipeline 3 predecessors are
reviewed. The uncorrected lane is one five-dataset wave; its Kidney inventory
filters valid legacy rows within that wave, so no same-root recovery gate is
serialized afterward.

The corrected Pipeline 3 contract passes each cell's original technical batch
metadata directly to Harmony/HVG. It performs no sample-level batch constancy
or majority check and does not invoke the retired corrected-source metadata/RDS
release preflight. Configuration, content, provenance, and the Covid obs-only
subset preflight remain required. The main first-run risks are malformed
metadata encodings, missing configured columns, Harmony/resource failures, and
dataset-specific source columns; any repair must remain row-targeted.

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
above, then publish one final runtime identity, then launch the separate
four-row uncorrected and nine-row corrected Pipeline 3 gates, followed by the
uncorrected five-dataset and corrected nine-dataset Pipeline 5 lanes after
their exact manifests and predecessor reviews are complete. Dataset rows in
each selection dispatch concurrently; the uncorrected five-dataset inventory
filters valid Kidney rows inside that single wave, and no separate uncorrected
recovery dependency is serialized.
### Plan maintenance

All subsequent implementation, verification, gate, SIF, failure, repair, and
scope updates MUST be recorded in this plan. Keep each status entry concise:
date, status, evidence/path, and next blocker or action. Do not duplicate the
full rationale when a short evidence-linked update is sufficient.

### Resume update

- 2026-09-12, resumed after the user-authorized context compaction; local
  stabilization work is dispatched in disjoint units for the obs worker,
  obs-only/onboarding paths, Stage 5 root ownership, final analysis/export,
  and the separate Stage 3 gates. No new gate or SIF command is authorized
  until those units are integrated and the focused contracts pass.

- 2026-09-12, corrected-contract scout confirmed exactly nine dynamic
  corrected datasets and their configured sample/label/batch keys; evidence
  was the current `datasets.json` plus `src/utils/py/batch_contract.py`.
  This was exploratory evidence, not a release requirement. The later
  approved contract passes each cell's original technical metadata directly
  to Harmony/HVG and does not require sample-level constancy, majority
  rewriting, or the retired corrected-source RDS preflight. Configuration,
  content, provenance, and Covid obs-only subset checks remain required.

- 2026-09-12, local verification caught a red combined Stage 3 fixture:
  `tests/test_preprocessing_stage_submitter.sh` reaches the validator-only
  combined run but its stubbed Covid preflight evidence is rejected. A
  focused test repair is dispatched; no HPC/SIF action is allowed until it
  passes.
- 2026-09-12, advisory-driven Stage 3 repairs classified exact four-row
  uncorrected and nine-row corrected manifests independently. The earlier
  combined-mode and corrected-source RDS preflight proposals are historical
  implementation evidence only and are superseded by the separate-gate,
  direct-cell-metadata contract above. Corrected rows must not be rejected for
  lack of the retired RDS release evidence.

- 2026-09-12, reconciled `AGENTS.md` with the clarification: frozen cohorts
  remain excluded from uncorrected/final-lane work but are explicitly allowed
  in the independent corrected Stage 3 and corrected Stage 5 lanes when
  selected by current `datasets.json`; disabled cohorts and `_debug` remain
  excluded everywhere. The later clarification retires the corrected RDS
  prerequisite as a release condition; historical audit evidence remains
  preserved below.

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

- 2026-09-12, the lightweight subset boundary and source-bound H5AD auditor
  were committed/pushed as `b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d`; a
  verified snapshot and fresh run-owned audit root were created at
  `_ecoda_source_snapshots/b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d` and
  `_ecoda_runs/corrected_h5ad_audit_b8f7aec_20260912b`.
- 2026-09-12, authoritative H5AD corrected-source auditing passed for
  `Covid19_PBMC`, `Kidney_KPMP_full`, `Diabetes`, and `Lung`, with zero
  split samples and estimable batch contracts; their checksummed reports are
  retained below the fresh audit root. It failed closed for `Alzheimer`
  because `assay` disagrees within 21 configured `donor_id` samples, for
  `Breast_cancer` because `suspension_dissociation_time` contains 65,359
  `unknown` sentinel cells, and for `Lupus_PBMC` because `batch_cov`
  disagrees within configured `sampleID` values. These are source/config
  contract failures requiring targeted scientific review, not reasons for a
  broad rerun. No Stage 3 worker or downstream gate was launched.

### Uncorrected Stage 3 launch

- 2026-09-12, the exact four-row uncorrected selection was staged under the
  canonical Bamboo gate input root and prepared/reconciled through
  `durable-hpc-gate-ecoda` as
  `stage3_uncorrected_final_20260912`, with reviewed Stage 2
  `stage2_joanito_final_20260912b` as the only predecessor. The wrapper uses
  the full-hash `b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d` source snapshot, the
  reviewed relocated runtime, canonical scratch/log roots, and no
  `--force`.
- 2026-09-12, the gate launched once and exactly one unbounded durable waiter
  was armed. Scheduler IDs and terminal status are intentionally not claimed
  until waiter completion, run-scoped audit, terminal inspect, and Luna Max
  review. The corrected wave remains blocked by the recorded source metadata
  failures/OOM and has not been launched.

### Uncorrected Stage 3 failed preflight

- 2026-09-12, `stage3_uncorrected_final_20260912` failed before any Stage 3
  worker array or watchdog was released. Its Covid obs-only preflight array
  ID was `4403790`; the remote error was the Slurm-spool bootstrap resolving
  the worker’s relative `slurm_config.sh` path instead of the inherited
  immutable snapshot source root.
- 2026-09-12, the required exact run-scoped audit was attempted and failed on
  the terminal `FAIL` status, then one terminal durable inspect ran with
  `4403790`. Accounting was queried once, reported `4403790|FAILED|1:0`,
  and the gate is `FAILED`/not release-eligible with no reviewer approval.
  The failed run and logs remain preserved; no retry or downstream Stage 5
  action is authorized until the container source-root bootstrap fix is in a
  fresh full-hash snapshot.

### Latest findings and pause boundary

- 2026-09-12, the authoritative corrected H5AD audit on snapshot
  `b8f7aec15dfc91bc493caeb04f3d2cc066db8c5d` passed
  `Covid19_PBMC`, `Kidney_KPMP_full`, `Diabetes`, and `Lung`, but failed
  `Alzheimer` (`assay` disagrees within 21 `donor_id` samples),
  `Breast_cancer` (`suspension_dissociation_time` has 65,359 `unknown`
  sentinel cells), and `Lupus_PBMC` (`batch_cov` disagrees within
  `sampleID`). The direct Joanito RDS audit remains exit-137/OOM; its raw
  H5AD cache has no producer record and is not authoritative.
- 2026-09-12, the container-spool bootstrap repair is committed/pushed as
  full commit `734b174a0b0b2a9c4e07edbf1e11d03c9fbf8206`; local syntax,
  subset, H5AD-auditor, Stage 3, Stage 5, and R contract checks passed before
  this commit. A fresh snapshot for this commit has not yet been created.
- 2026-09-12, terminal inspect of failed gate
  `stage3_uncorrected_final_20260912` queried accounting once for preflight
  ID `4403790` and observed an additional `4403791|RUNNING|0:0` row beside
  `4403790|FAILED|1:0`. The extra scheduler job is unresolved; do not
  relaunch, cancel, mark the failed gate settled, or treat the four-row
  selection as validated until its identity and terminal state are resolved.
  Preserve the existing wait/inspect evidence.
- 2026-09-12, the user requested context compaction and a pause after this
  status update. No further snapshot, scheduler, gate, Stage 5, sync, or
  analysis action is authorized until the next explicit resume signal.

### Explicit resume and repaired gate

- 2026-09-12, the user’s later instruction to “forget about Snakemake, and
  continue with the plan implementation” is the explicit resume signal that
  supersedes the earlier pause boundary. The prior extra job
  `4403791|RUNNING|0:0` was identified as `h5ad_obs_audit_worker.sh` and
  subsequently settled `COMPLETED|0:0`; it required no cancellation.
- 2026-09-12, the container source-root repair was snapshotted at
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/734b174a0b0b2a9c4e07edbf1e11d03c9fbf8206`.
  A new exact four-row uncorrected gate
  `stage3_uncorrected_final_20260912b` was prepared, reconciled, launched
  once, and placed under its single unbounded durable waiter. The gate uses
  no `--force`; terminal scheduler IDs/status remain pending waiter
  completion and the required single inspect/review sequence.
- 2026-09-13, the user clarified that `majority_v1` is the approved
  sample-level technical metadata policy for Alzheimer, Breast_cancer, and
  Lupus_PBMC. No minimum winner fraction is required; exact ties and actual
  missing/blank/non-finite values remain hard failures. Breast treats only the
  literal `unknown` value in `suspension_dissociation_time` as an ordinary
  configured class.
- The authoritative Bamboo obs-only audit measured majority fractions of
  66.96--100% for Alzheimer `assay` (median 100%), 100% for every configured
  Breast technical key (including 65,359 `unknown` dissociation-time cells
  across 11 samples), and 30.95--100% for Lupus `batch_cov` (median 100%).
  It found no ties. The audit read only configured `obs` metadata through the
  pinned container and found no source `.md5` sidecars.
- Majority implementation is Python-canonical: the shared Python validator
  extracts H5AD/RDS metadata, selects unique per-sample technical winners,
  and records winner values, counts, and fractions. Corrected R consumers
  consume that serialized result through `reticulate`; they do not recompute
  votes. Original cell-level batch values remain unchanged for Harmony/HVG,
  while ECODA/limma/pseudobulk receive the sample-level technical winners.
  Biological labels and sample IDs are never voted.
- New majority summaries use schema 2; strict schema-1 summaries and legacy
  first-observation paths remain readable. The policy is carried in
  `datasets.json` only for the three approved datasets. Focused Python,
  R/Python handoff, H5AD/RDS, Stage 3 submitter, Stage 5 selection/sync,
  batch-correction, and multibatch contract checks pass. No pipeline rerun,
  scheduler submission, or existing-artifact invalidation occurred for this
  policy change.
- The completed gate `stage3_uncorrected_final_20260912b` used the intended
  four-row uncorrected scope. Scheduler roots `4403794` and `4403795`
  completed successfully; terminal inspect passed and Luna Max approved the
  gate. Its run-scoped audit still reports a NAS owner-state discrepancy:
  scratch H5AD ownership is terminal `OK`, while the synchronized NAS owner
  remains `ACTIVE`. Do not advance to further pipeline work until that
  ownership issue and the explicit user approval for any rerun are resolved.
- 2026-09-13, the user narrowed the majority policy to downstream
  sample-level metadata only. The implementation must not apply majority to
  Pipeline 3, Harmony, HVG, or any cell-level batch column. A separate
  Pipeline 3 cell-level cutover is pending verification: it removes only
  sample-level batch checks while retaining ordinary H5AD content and
  configuration/provenance checks.
- The remaining policy is one obs-only Python Feather exporter. The scoped
  `datasets.json` policies use `majority_keys=["assay"]` for Alzheimer,
  `["suspension_dissociation_time"]` for Breast_cancer, and `["batch_cov"]`
  for Lupus_PBMC. Breast accepts literal `unknown` only for its dissociation
  key. Other technical keys, biological labels, and Sample IDs retain
  first-observation behavior.
- The exporter streams per-sample counts, selects the unique highest-count
  class without a minimum fraction threshold, and writes only sample-level
  metadata plus its checksum. Focused exporter, configuration, source
  storage, subset, Stage 5 selection/sync, RDS, batch-correction, and
  analysis contract checks pass. `NOTES.md` now records the audit ranges,
  Lupus lower-tail distribution, policy, and limitation.
- The user authorized remaining HPC processing after this implementation and
  documentation step. No new pipeline launch, rerun, artifact rewrite, or
  existing-artifact invalidation has occurred. Before the next launch, resolve
  the existing `stage3_uncorrected_final_20260912b` NAS owner-state audit
  discrepancy and record the exact approved selection/runtime command.
- 2026-09-13, the user superseded the earlier corrected-Stage-5 exclusion
  boundary: the end goal now includes a separate corrected Stage 5 benchmark
  lane. This does not change the cell-level Stage 3 contract:
  Harmony/HVG continue to use each cell's original batch metadata, while
  sample-level corrected Stage 5 consumers use the approved majority
  assignments only for affected keys.
- Corrected Stage 5 must use a separate non-legacy analysis root and explicit
  run-owned selection/manifests. It is a separate submission from uncorrected
  Stage 5; independent corrected/uncorrected dataset and method rows may run
  in parallel after their required Stage 3 inputs are reviewed. Existing
  valid rows remain reusable; only missing/invalid rows are selected.
- The corrected Stage 5 contract is now frozen in section 7: nine configured
  datasets, seven methods, `_corrected_final` stems, the
  `batch_effect/corrected_final` root, and an explicit eight-key artifact
  manifest. No corrected Stage 5 scheduler job has been submitted.
- 2026-09-13, the user clarified that uncorrected and corrected processing
  are separate jobs/gates and may run in parallel. `CombinedPBMC` is removed
  from the active analysis registry (`use_for_batch_effect=false`) and must
  not appear in any new selection. The current corrected registry therefore
  contains the nine non-underscore, batch-enabled datasets only.
- Pipeline 3 corrected processing has no sample-level batch check. It keeps
  each cell's original technical metadata for Harmony/HVG; no cell-level
  majority assignment is performed. The old corrected-source metadata
  release/evidence preflight is no longer invoked. The separate Covid
  obs-only subset preflight remains for both declared views because it audits
  the approved subset, not batch constancy.
- Corrected Stage 5 is now an explicit end-goal lane. It consumes the
  corrected Stage 3 H5ADs, writes to a separate `batch_effect/corrected_final`
  root, and uses majority-derived sample metadata only for the affected
  `Alzheimer/assay`, `Breast_cancer/suspension_dissociation_time` (with
  literal `unknown` as a class), and `Lupus_PBMC/batch_cov` covariates.
  Uncorrected Stage 5 performs no batch correction and no majority assignment.
  Corrected and uncorrected Stage 5 gates are separate and may run in
  parallel after their respective Stage 3 predecessors are terminally
  reviewed; no corrected Stage 5 scheduler job has yet been submitted.
- The uncorrected final Stage 5 selection is one five-dataset wave:
  `Covid19_PBMC`, `Diabetes`, `Joanito`, `Lung`, and `Kidney_KPMP_full`.
  It uses the seven-method batch suite and reuses valid Kidney rows, selecting
  only missing/invalid methods. The corrected Stage 5 selection is one
  nine-dataset wave with the same seven methods, concurrent dataset rows,
  explicit `_corrected_final` stems, and the separate corrected root described
  in section 7.
- 2026-09-13, the user clarified that Pipeline 3 has no sample-level batch
  constancy check: each cell's original technical metadata is passed directly
  to Harmony/HVG, with no within-`Sample` constancy check and no cell-level
  majority assignment. The majority assignment is exclusively a corrected
  Stage 5 sample-level metadata operation.
- The user confirmed the execution topology: uncorrected and corrected
  processing use separate jobs/gates and may run in parallel. `CombinedPBMC`
  is excluded from all new selections because its current
  `use_for_batch_effect` flag is false. The corrected Stage 3 selection is
  the nine current non-underscore, batch-enabled datasets.
- The corrected Stage 5 lane is part of the end goal. It is a separate
  submission/root from uncorrected Stage 5; it consumes corrected Stage 3
  H5ADs, applies majority only to corrected sample-level technical
  covariates for the three affected datasets/keys, and otherwise preserves
  native cell-level handling. Uncorrected Stage 5 performs no batch
  correction and no majority assignment.
- The uncorrected Stage 5 wave contains the four changed datasets plus
  `Kidney_KPMP_full` in one explicit selection. The corrected Stage 5 wave
  contains the nine active corrected datasets in one explicit selection.
  Both use seven baseline batch methods; existing valid artifacts remain
  reusable and valid Kidney rows are not forced.
- Implementation and HPC launch are paused at the user's request for memory
  compaction. No pipeline job is to be submitted until the user explicitly
  resumes execution. The next resume must verify the pending code cutover,
  update manifests/runtime identity, reconcile the prior NAS owner-state
  discrepancy, and then launch the separate reviewed gates in this order
  (or concurrently where the documented dependencies permit).

### Documentation-only clarification status

- 2026-09-13, the latest user clarification was reconciled into `AGENTS.md`
  and this plan: Pipeline 3 now has separate parallel four-row uncorrected
  and nine-row corrected selections; corrected Pipeline 3 uses original
  cell-level batch metadata without sample-level constancy/majority checks or
  the retired source-RDS release preflight; and both Stage 5 lanes are
  explicit, idempotent, and root-separated.
- The corrected Stage 5 contract is
  `--pass corrected --analysis-variant corrected_final`, root
  `batch_effect/corrected_final`, nine datasets, seven methods, concurrent
  dataset rows, and an explicit eight-key artifact manifest. No scheduler
  launch, test, or validation was performed by this documentation-only update.
  Historical gate evidence and the separate structural-cleanup follow-up plan
  remain preserved.

### Interrupted-agent audit

- **(a) Canceled partial work superseded and restored.** `FixRPolicyBoundaries`
  (R majority/validator plumbing), `InferCompositionPolicy` (R composition
  fallback), `FixPythonPolicyAndMrvi` (corrected Stage 5 Python worker),
  and `StreamMetadataExporter` (overbroad exporter cache/streaming changes)
  were cancelled while their broad implementations were in progress.
  `RestoreRCoreBaseline`, `RestoreStage5Baseline`, and
  `RestorePythonCoreBaseline` replaced their assigned source files with exact
  `734b174a0b0b2a9c4e07edbf1e11d03c9fbf8206` blobs; each restoration reported
  a clean baseline diff for its assigned paths. `WriteMetadataDecisionNotes`
  was cancelled before its documentation result; the parent later wrote and
  re-read the final `NOTES.md` section.
- **(b) Completed edits retained but stopped for compaction.**
  `CellLevelPreprocess` owns the pending Pipeline 3 cell-level cutover in
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`;
  `AdjustH5ADContentContract` and `FixStage3H5ADCli` own the corresponding
  summary-optional Pipeline 3 H5AD validation in
  `src/utils/py/benchmark_h5ad_contract.py`;
  `RemoveCorrectedPreflights` owns the removal of corrected sample-level
  source-preflight invocations in
  `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (and its focused test).
  These agents were stopped after reporting completed edits when the user
  requested documentation-only work for compaction. Their agent-local
  `py_compile`, H5AD smoke, and `bash -n` checks passed, but the parent must
  re-read and run the final union of focused checks before launch; they are
  not marked as a completed clean cutover yet.
- **(c) V2 shell/test work superseded by Stage 3 restoration.**
  `HardenV2Evidence`, `ScopeConfigOnlyEvidence`, and `AddPolicySchemaTest`
  modified the temporary schema-2 Stage 3 shell/test path. Those changes were
  superseded and removed by `RestoreStage3Baseline`; no V2 shell gate or
  V2-only test edits are retained. The restoration agent checked the assigned
  Stage 3 files against the exact baseline before the retained cell-level
  cutover edits above.
- The cancelled exporter left a malformed intermediate block, but
  `MinimalMetadataExporter` restored the exporter from the baseline and
  implemented the final narrow obs-only Feather reducer; `FixAffectedKeyLoop`
  made the final affected-key-only correction. The focused
  `tests/test_batch_majority_contract.py` regression passed, including
  source-storage preservation, tie/missing handling, Breast `unknown`, and
  first-row preservation for unaffected fields.
- No interrupted agent launched a scheduler job or intentionally rewrote an
  existing production artifact. The prior completed uncorrected gate and
  its preserved NAS owner-state discrepancy remain separate evidence.
- On resume, first inspect the retained Pipeline 3 edits and current
  `datasets.json`/plan scope, run the focused cutover checks, verify the
  separate Stage 3 and Stage 5 manifests/runtime identities, and only then
  submit the explicitly reviewed HPC gates.
### Resume implementation status

- 2026-09-13, the user resumed implementation after compaction. The local
  implementation wave is integrated across the clarified Pipeline 3 split,
  corrected-final Stage 5 lane, majority metadata handoff, Kidney legacy
  inventory reuse, run-audit variant binding, and corrected-final watchdog
  retry identity. `AGENTS.md` and this plan now describe the same separate
  parallel scopes and roots.
- The corrected Stage 5 lane is frozen as
  `--pass corrected --analysis-variant corrected_final`, with the nine current
  configured datasets, seven methods, `_batch_effect_corrected_final_` stems,
  root `batch_effect/corrected_final`, and explicit eight-key artifact mapping.
  The uncorrected final lane is one five-dataset wave including
  `Kidney_KPMP_full`; valid legacy Kidney methods are recorded in a
  run-owned inventory and skipped rather than recomputed.
- Pipeline 3 retains original cell-level batch metadata for Harmony/HVG and
  performs no sample-level batch constancy or majority operation. Pending
  Covid work in either approved Stage 3 selection requires both run-owned
  obs-only reports; the retired corrected-source metadata/RDS release
  preflight is not invoked.
- Parent verification is green for the subset, majority/exporter, Stage 3
  submitter, Stage 5 selection (including Kidney legacy reuse), Stage 5 sync,
  RDS, batch-analysis, registry, multibatch, corrected-source, batch-correction,
  and corrected-final watchdog regressions. Shell syntax, Python compilation,
  R parsing, and `datasets.json` validation also pass. The new
  `tests/test_ecoda_run_audit.sh` fixture has received several test-only
  repairs; its latest remaining failure was a literal-tab encoding in the
  synthetic metadata-export row, corrected by `FixAuditManifestTabs`, but the
  parent has not rerun that test after the edit.
- The no-variant corrected R path remains strict and legacy-compatible after
  `FixLegacyCorrectedMode`. Corrected-final Breast cell composites accept only
  the configured literal `unknown` for
  `suspension_dissociation_time`; all other sentinel paths remain strict.
  The prior Joanito RDS audit exit 137 remains historical evidence and is not
  a Pipeline 3 prerequisite under the clarified contract.
- `FixAuditManifestTabs` was canceled when the user requested this
  documentation update after reporting its test-only edit; that edit remains
  unverified. No implementation agent launched HPC work, submitted a
  scheduler job, rewrote an existing production artifact, or invalidated a
  gate. The prior reviewed Stage 2/Stage 3 evidence and the Stage 3 NAS owner
  discrepancy remain preserved.
- On the next resume, first rerun `tests/test_ecoda_run_audit.sh` and the full
  focused union after the last test-fixture edit. Then perform the
  validator-only NAS owner reconciliation, validate the reviewed Joanito
  Stage 2 predecessor contract, freeze one full-hash source snapshot/runtime,
  and generate exact idempotent Stage 3/Stage 5 manifests. Only after those
  checks pass may the separate durable Stage 3 and Stage 5 gates be prepared
  or launched, with one durable wait, terminal inspect, and reviewer approval
  per gate.
- The Stage 3 artifact-owner lifecycle blocker is resolved. The watchdog now
  reconstructs both scratch and NAS owners from the bound root selection,
  persists those artifact owners beside the stage-owner manifest, keeps them
  `ACTIVE` until the submitter verifies scratch-to-NAS sync, and transitions
  every tracked owner to `FAIL` on failure without deleting owner directories.
  `bash tests/test_preprocessing_stage_submitter.sh` passed with explicit
  scratch/NAS `ACTIVE` → `OK` and `ACTIVE` → `FAIL` assertions.
- Validator-only reconciliation of the reviewed
  `stage3_uncorrected_final_20260912b` gate verified the stale NAS owner PID
  was absent, all four NAS H5ADs were regular/nonempty, their sizes and
  sidecar/checksum identities matched the run records, and their owners were
  `RUN_ID`/`STAGE`-bound before changing only those four owner records to
  terminal `OK`. No H5AD, checksum, or artifact record was rewritten.
- The reviewed Joanito Stage 2 predecessor remains release-eligible and its
  run-owned records are `PUBLISHED`. The terminal watchdog log records the
  exact `joanito` hook, 373,058 cells across 189 samples, current
  `seqtec`/`cell.type_new`, and the 2,500-cell five-sample debug artifact.
  A separate read-only direct RDS recheck was attempted but terminated with
  remote exit 255 during startup; it did not mutate the RDS, so the reviewed
  watchdog semantic evidence remains authoritative.
- Implementation commit `59b781fa31e6d9fb015e1d7911aa283f5131ea6c` was
  pushed to `origin/master`. Its verified source snapshot is
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots/59b781fa31e6d9fb015e1d7911aa283f5131ea6c`
  with `COMPLETE`, a full source archive, and a matching
  `identity/source.manifest`. The first create attempt was rejected before
  writing because the textual `$HOME/scratch` path is a symlink; retrying with
  the canonical `/srv/beegfs/scratch` parent succeeded.
- Versioned runtime
  `ecoda-py-cuda13-6bbf70b-relocated` was reused after direct SHA-256, size,
  regular-file, and manifest checks. The image identity is
  `8fcd00b7a02592f82e8ad96ec9dffdcfb40214b3c47218f770dfb8e828e05e34` and
  the runtime-manifest identity is
  `1c4c2cc340353fe86f1b89235f1f3e5c7b3545f7f303d2b85612cfa8cdedfd50`;
  its Pixi TOML/lock hashes match the source snapshot. No SIF build or
  scheduler work was performed.
- The exact Stage 3 manifests are now run-owned under the canonical scratch
  gate tree: `stage3_uncorrected_final_20260913a/selection.tsv` has four
  rows, MD5 `b06c87a68cabce5d6ce792aa03ee80e4`, and size 135; 
  `stage3_corrected_final_20260913a/selection.tsv` has the nine rows in the
  declared config order, MD5 `3c43044510b35ab9b88a86fdd3b8d467`, and size
  305. No frozen or disabled cohort appears in the uncorrected manifest;
  the corrected manifest is exactly the approved nine-dataset lane.
- The snapshot-backed uncorrected Stage 3 selection was validator-only:
  `stage3_uncorrected_final_20260913a` validated all four existing final
  H5ADs, synchronized their checksums, wrote `NOOP_VALIDATED`, and emitted no
  scheduler IDs. Its reviewed predecessor therefore releases the independent
  uncorrected Stage 5 lane.
- Corrected Stage 3 gate `stage3_corrected_final_20260913a` was prepared,
  reconciled as absent, and launched once through `durable-hpc-gate-ecoda`
  with the 59b781f snapshot, validated relocated runtime, reviewed Joanito
  predecessor, and nine-row corrected selection. The initial durable wait
  lost its completion transport and entered `PRELAUNCH_STOP`; recovery
  `status` later found the remote runner terminal `FAILED` with exit 1 at
  `2026-09-13T11:05:07Z`. The single terminal inspect used every emitted
  scheduler ID: preflight `4403877`, initial array `4403887`, OOM retry
  arrays `4403903` and `4403989`, and watchdog `4403888`. Accounting failed
  closed because `4403903` was `OUT_OF_MEMORY|0:125` while its later retry
  `4403989` completed; no Luna Max reviewer approval was requested.
- All nine corrected scratch H5ADs and run-owned artifact records passed
  content, layer, sidecar, and checksum validation. The wrapper failed during
  NAS synchronization because the corrected batch-contract validator tried to
  treat its long inline JSON identity as a filesystem path and hit the
  platform filename-length limit; this was a validator parser defect, not
  invalid H5AD content. A later unsafe sync-only probe was stopped before
  synchronization and overwrote only the remote terminal reason; the original
  durable inspect evidence remains preserved locally, and the failed gate is
  not release-eligible. Corrected Stage 5 remains blocked pending a distinct
  validator-only repair path, terminal evidence, and required review.
- User clarification on 2026-09-13 supersedes the stale split wording in
  earlier Stage 5 subsections: all five uncorrected datasets run together in
  one explicit `final` wave. The canonical selection is
  `Covid19_PBMC`, `Diabetes`, `Joanito`, `Lung`, and `Kidney_KPMP_full`, with a
  35-row method ceiling and validator-only per-method Kidney reuse. Missing
  or invalid Kidney methods remain part of that same wave; no separate
  `--target-methods` recovery gate or same-root serialization follows it.
- The first uncorrected Stage 5 gate attempt,
  `stage5_uncorrected_final_20260913a`, failed closed before emitting any
  scheduler ID because the active corrected Stage 3 executor held the shared
  source-snapshot parent lock. Its one terminal inspect attempt stopped before
  accounting because no scheduler ID existed; this was not a dataset or
  artifact result, and no worker or production artifact was created.
- To preserve the requested overlap, a separate same-commit read-only source
  snapshot was created under
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots_stage5_uncorrected`.
  Recovery gate `stage5_uncorrected_final_20260913b` uses that snapshot,
  retains the exact five-row selection and reviewed uncorrected Stage 3
  dependency, was prepared, reconciled, and launched successfully, and has
  one unbounded durable wait active. Its terminal scheduler IDs, accounting,
  run-scoped audit, and Luna Max review remain pending.
- The recovered uncorrected Stage 5 gate
  `stage5_uncorrected_final_20260913b` reached terminal `FAILED` at the
  aggregate gate. Its one terminal inspect queried all recorded metadata,
  preflight, method-array, watchdog, and aggregate IDs exactly once; the
  audit is not release-eligible and no reviewer approval was requested.
  Metadata export `4403898`, H5AD preflight `4403904`, R-environment
  preflight `4403909`, preparation/pseudobulk/pilot/qot watchdogs and their
  arrays completed. The failed method classes are GloScope
  (`4403920`/`4403925`), composition (`4403931`/`4403940`), and MRVI
  (`4403946`/`4403947`); aggregate gate `4403988` failed accordingly.
- The terminal worker evidence is concrete and row-scoped: MRVI rejected
  implicit `auto` execution because the worker omitted its already-built
  `--device cuda` argument; GloScope fell through to the count-backed
  `load_benchmark_seurat` path even though the selected H5ADs contain
  `layers=['counts']`; composition found no visible final `hvg2000`
  pseudobulk cache at worker time. The four final hvg2000 caches are now
  present, checksummed, and run-owned. The valid Kidney legacy inventory
  skipped six methods; only Kidney GloScope was missing and remains a
  targeted affected row.
- Successful prepare, pseudobulk, pilot, and qot artifacts remain immutable
  and reusable. The four final hvg2000 pseudobulk caches were revalidated
  against their checksums and run-owned publication records. The valid Kidney
  legacy inventory skipped six methods; only Kidney GloScope was missing.
- Targeted repair is implemented: the Python worker appends its explicit
  `GPU_DEVICE_ARGS`, and the R worker has a standalone GloScope counts-free
  path plus a separate post-cache no-op branch before scITD. The focused
  worker regressions pass (`test_benchmark_matrix_submitter.sh` and
  `test_benchmark_worker_dispatch.R`), as do the Stage 5 selection,
  synchronization, RDS, and batch-analysis contracts. No production artifact
  was changed by the repair.
- The affected-row manifest
  `stage5_uncorrected_repair_20260913a/affected_methods.tsv` records exactly
  13 failed method/dataset rows: five GloScope, four composition, and four
  MRVI. Its durable gate was launched but stopped through a verified
  process-tree termination before `pending_selection.tsv` or any scheduler
  ID was emitted, after validation showed that the failed producer had marked
  successful prepare/pseudobulk/pilot/qot owners `FAIL`. No recovery compute
  or production artifact was created by that attempt.
- The 16 successful method rows were revalidated validator-only against their
  terminal watchdog states, exact run-owned matrix manifests, scratch
  checksums, records, producer identities, and stage/global owner metadata.
  The run-owned `owner_promotion_16.tsv` report records one validated scratch
  payload per prepare/pseudobulk row and payload plus runtime metadata for each
  pilot/qot row. The historical stage and scratch payload owners were promoted
  to `OK`; the missing NAS destinations were explicitly left unsettled.
- A separately scoped no-compute selected sync then transferred exactly those
  16 method rows (30 payload/runtime files plus the five metadata outputs and
  merged execution log), normalized NAS sidecar `PATH` fields atomically, and
  verified every transferred digest. Its run-owned
  `selected_sync_16_post_audit.tsv` and status report record 16 `NAS_OWNER_OK`
  rows and 16 finalized NAS payload owners. No worker or scheduler job was
  created; the original aggregate `FAIL` evidence remains preserved.
- The targeted final selector now permits either the one-row Kidney exception
  or an exact five-row selection with explicit method classes. A targeted
  composition repair requires a validated hvg2000 prepare cache and never
  broadens the pending method scope. The five-row
  `gloscope,composition,mrvi` selector produces exactly 13 pending rows
  (5 + 4 + 4), verified by `tests/test_benchmark_selection_file.sh`. No
  corrected Stage 5 gate is permitted before corrected Stage 3 review.
- The corrected-batch H5AD contract loader now recognizes long inline JSON
  identities before filesystem probing, while preserving readable-path and
  malformed/object validation. `tests/test_benchmark_h5ad_contract.py` passes
  for long inline JSON, file-backed JSON, malformed JSON, and non-object JSON.
- Stage 5 failure finalization now preserves a run-owned stage or global
  artifact owner only when every matching run-owned matrix manifest has a
  correctly labeled terminal watchdog `STATE=OK`; missing, malformed, failed,
  or foreign status/owner records remain fail-closed without deleting owner
  directories. An empty acquisition list is treated as “nothing to finalize,”
  not as a malformed owner. Sync cleanup is subshell-wrapped so failures
  return through `stage5_abort`. `tests/test_benchmark_matrix_submitter.sh`
  passes its aggregate-failure and true sync-boundary owner-state regressions.
- The focused repair union is green after these fixes: Stage 5 matrix
  ownership, H5AD contract parsing, exact selection, synchronization, worker
  dispatch, RDS contracts, and batch-effect analysis. The expected negative
  runtime/source-escape diagnostics remain covered by the matrix test.
- The selected-sync post-audit is run-owned under
  `stage5_uncorrected_final_20260913b/manifests/selected_sync_16_post_audit.tsv`
  and its `STATE=OK` status. It covers 30 method payload/runtime files plus
  five metadata outputs and the merged execution log; no scheduler work was
  emitted by either validator-only command.
- Commit `837eafb80c5b204ba6999f3c51ec8e4d95b09ddb` contains the exact
  five-row targeted selector contract and its validated hvg2000 dependency
  reuse. Its verified snapshot is
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_source_snapshots_stage5_targeted_recovery/837eafb80c5b204ba6999f3c51ec8e4d95b09ddb`;
  source archive SHA-256 is
  `d0034d6730f5496bb8041d60c5edcbdb4f029e5b4097d2df833f3d72a320a36d`.
  The reused runtime remains the versioned
  `ecoda-py-cuda13-6bbf70b-relocated` identity recorded above.
- A first targeted Stage 5 recovery manifest was prepared with the obsolete
  `ecoda-stage5-uncorrected` serialization group, then canceled before any
  remote launch; its local `launch_intent` and absence of remote status,
  runner, and scheduler IDs are preserved. It was superseded by
  `stage5_uncorrected_targeted_recovery_20260913b`, prepared and launched
  with `ecoda-benchmark`, the exact five-row selection, explicit target
  methods `gloscope,composition,mrvi`, and expected pending count 13. Its
  single unbounded durable wait is armed; terminal accounting, run-scoped
  audit, and Luna Max review remain pending.
- A validator-only corrected Stage 3 sync repair is now implemented without
  changing `stage3_load_bound_run`: the new report utility validates all
  corrected H5AD content with the new validator while recording both the
  failed run's immutable 59b source identity and the repair snapshot identity.
  `--validated-sync-report` makes the existing sync-only path consume that
  report while retaining H5AD artifact-record/owner checks and selected sync.
  The focused report generator and Stage 3 submitter tests pass. The repair
  source changes are not yet committed or snapshotted; corrected Stage 5
  remains blocked until the new sync repair is terminally inspected and
  reviewed.
