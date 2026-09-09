# Onboarding and established batch-effect uncorrected run and batch-column gate

## Status

Proposed plan only. No scheduler submission, source publication, commit, or push is authorized by this plan. The existing authoritative corrective plan remains the implementation baseline.

## Scientific and publication boundaries

Scope the first pass to the twelve cohorts that are batch-effect datasets in `datasets.json`: the nine onboarding cohorts (`Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Diabetes`, `Kidney_KPMP`, `Lupus_PBMC`, `Lung`, `Myocardial_infarction`, and `Parkinson`) plus `Joanito`, `Stephenson`, and `CombinedPBMC`.
- Run only uncorrected-compatible views in the first preprocessing/annotation/benchmark pass. Never let default all-view resolution submit a corrected view.

- `Joanito` and `Stephenson` use their declared `batch_effect_uncorrected` views.
- `CombinedPBMC` must be migrated from its legacy `batch_effect_analysis` view
  to the canonical `batch_effect_uncorrected` view before execution. Do not
  retain a preprocessor compatibility alias or submit the legacy name. Its
  corrected view remains blocked unless a separate explicit view-contract
  decision adds one.
- The Stephenson row means the full declared batch-effect cohort (`Status` Healthy/Covid with the two configured sample exclusions), not the narrower `benchmark_analysis` view restricted to `Site == Ncl`.
- Stage 4 skips exactly `Alzheimer`, `Diabetes`, and `Parkinson`; the other nine cohorts require both automated annotation methods.
- Keep the nine onboarding `columns.batch` values `null` until the candidate-evidence review. Re-evaluate the existing `Joanito`, `Stephenson`, and `CombinedPBMC` values rather than treating them as permanently accepted.
- Biological labels remain evaluation metadata. They are not preprocessing, HVG, normalization, Harmony, embedding, model, or correction covariates.
- **No `git commit`, `git push`, or equivalent publication command may run before explicit user confirmation authorizing publication.** The implementation turn remains uncommitted and unpushed. The existing durable-gate plan additionally requires the later exact `go` boundary before scheduler/artifact reconciliation and publication.
- Do not alter or delete preserved `.gate`, `.kilo`, archived-plan, figure, or other user-owned artifacts while implementing this plan. Reconcile stopped/active scheduler state only after the user-controlled execution boundary.


## User decisions and implementation rationale

- This is a two-pass batch-effect analysis. The first pass is strictly
  uncorrected and is the evidence gate for technical batch selection. The
  second pass is separate and must not run while any selected
  `columns.batch` value remains unconfirmed.
- The biological label is evaluation-only. It must never enter sample/cell
  filtering, HVG selection, normalization, PCA, Harmony, model covariates,
  annotation inputs, or correction formulas. Technical candidates are
  evaluated after the uncorrected run; they are not selected by convenience or
  p-value alone.
- The raw-observation invariant is strict: for each configured sample unit,
  fewer than 500 cells means the entire sample unit is removed; 500 or more
  cells means the sample unit is retained. This is why existing pre-filter
  h5ads cannot be reused as final inputs. The raw/staged source, selected
  `subset_vars`, and configured sample column determine the units.
- All suitable cohorts must receive both HiTME and scATOMIC annotation even
  when author annotations already exist; author annotations remain baseline
  metadata, while the two automated methods provide standardized outputs.
  Only datasets explicitly carrying both method names in
  `not_suitable_for_auto_annotation` are skipped.
- Parallelism is required within every stage: independent Stage 2 hooks,
  Stage 3 dataset/view rows, Stage 4 preparation/annotation/merge rows, and
  Stage 5 dataset/method rows use arrays. Dataset-by-dataset serial loops are
  not an acceptable substitute. Stage dependencies still impose barriers:
  complete and review one stage before starting the next.
- `meta_cols_keep` was a dead legacy field. It was never used to subset
  observation columns; it has been removed from the registry, adapters,
  audits, reports, and active documentation. Do not recreate it under another
  name or use it as a preprocessing contract.
- The current local worktree contains an unfinished redesign of the run
  ownership, manifest, checksum, and array orchestration. The implementation
  agent must first finish and locally verify that redesign, then establish one
  clean immutable Bamboo revision. No HPC execution may consume the mixed
  pre-redesign artifacts as if they were final.

## Findings from the current worktree

### Stage 2 — correct dispatcher, one force-propagation gap

`src/2_dataset_specific_preprocessing/1_submit_hpc.sh` is the correct production dispatcher. It narrows `--datasets` to the selected hook, keeps GongSharma -> CombinedPBMC as the only dependency edge, and maps `myocardial_counts` to `1.5_submit_myocardial.sh`.

The planned Myocardial rerun must therefore use the dispatcher, not call the hook directly. One repair is still required before a forced rerun: the dispatcher exports `FORCE_PREPROCESS`, but `1.5_submit_myocardial.sh` does not pass `--force` to `1.5.1_reconstruct_myocardial_counts.py`. Because the reconstruction helper skips an existing valid `layers["counts"]`, dispatcher `--force` currently invalidates only the sidecar and cannot force recomputation.

Before Stage 3, the same dispatcher must also validate or build the `Joanito` preparation (`seqtec` and `_debug`) and the `CombinedPBMC` input (`gongsharma_cap` -> `combinedpbmc`) when their checksums/artifacts are absent or invalid.

### Stage 3 — legacy batch submitter is gone; canonical entrypoint is correct

Only `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` remains under the Stage 3 submitter pattern. It builds an immutable `DATASET<TAB>VIEW` selection, validates existing h5ads/checksums, submits one pending array, and gates it with the watchdog.

Use an explicit selection manifest. All twelve rows use
`batch_effect_uncorrected` after the CombinedPBMC view migration. Do not use
the deleted `1_submit_batch_effect_stage.sh`, and do not omit the view because
omission resolves every declared view, including corrected.

The legacy token must also be removed from the allowed-view sets and h5ad
contracts (`1.1.1_preprocess.py`, `benchmark_h5ad_contract.py`, loaders, tests,
and evidence code). This is a schema migration, not a compatibility alias:
rename the CombinedPBMC registry view and its output to
`batch_effect_uncorrected`, then regenerate any dependent manifests.

`NOTES.md` records the raw audit counts that motivate the refresh: Breast
cancer drops 2 sample units, COVID-19 PBMC 8, Kidney 2, Lupus 1, Lung 19,
and Parkinson 1; Myocardial drops none but still requires Stage 2 raw-count
reconstruction; Alzheimer and Diabetes drop none. Joanito also requires a
fresh Stage 3 run because its prior uncorrected preprocessing did not apply
the `<500`-cell sample-unit filter. For clean traceability, the execution
decision is to run all twelve batch-effect datasets through Stage 3 after the
Myocardial Stage 2 prerequisite, even when the audit predicts no removals.
Do not reuse any pre-filter h5ad.

The unfinished local Stage 4 redesign introduces run-owned roots, preparation,
annotation, and merge arrays, per-dataset ownership, atomic manifests, and
parallel synchronization. Retain that concurrency model, but do not consume
its outputs until the strict dual-method contract below passes. Its
`--skip-prepare --reuse-run` path must validate the immutable preparation,
chunk, feather, and checksum records rather than trusting file presence.

### Stage 4 — explicit view selection works, dual-method contract is too weak


Observed failures that motivate these repairs:

- The prior suitable-cohort run produced Breast feathers with HiTME but no
  scATOMIC. Existing nonempty feathers were skipped without checking method
  completeness, and the worker did not pass scATOMIC's `breast_mode=TRUE` for
  breast tissue. The final worker must pass breast-specific mode to both
  `run_scATOMIC()` and `create_summary_matrix()`, and must preserve
  `normal_tissue=false` (a jq `false // empty` expression incorrectly erased
  false values).
- A merge then failed while writing mixed object-valued HiTME labels
  (`layer2`/`layer3`) to HDF5. Merge code must normalize nonnumeric annotation
  labels to nullable strings while retaining numeric score columns, before the
  atomic h5ad write.
- Parallel R workers encountered missing lazy-load files in shared packages
  (`Biostrings`, `AnnotationDbi`, `arrow`, `SummarizedExperiment`) after a
  partial shared environment mutation. The guarded environment refresh must
  complete successfully, include the pinned `robCompositions` dependency, and
  pass the full package-integrity/import smoke check before Stage 4.

`src/4_cell_type_annotation/1_submit_onboarding_stage.sh` is the canonical production entrypoint and its explicit selection-file path is suitable for the first pass. Its default path still resolves all declared views, so the first run must use a selection manifest containing only the twelve uncorrected-compatible rows.

The current `src/utils/py/annotation_contract.py` accepts any one known annotation column. The merge path likewise accepts any nonempty annotation column, and `2.1.1_process_chunk.R` deliberately converts all-NA/missing-method output into a warning/checkpoint. This permits a suitable cohort to pass with HiTME but no scATOMIC (or vice versa), matching the reported Breast failure.

The contract must instead require, for every runnable cohort and every selected view:

- all HiTME layers `layer1`, `layer2`, `layer3`;
- all scATOMIC layers `layer_1` through `layer_6`;
- `scATOMIC_pred`, `classification_confidence`, `S.Score`, `G2M.Score`, and `Phase`;
- unique `(Sample, cell_barcode)` keys and complete cell-key coverage; and
- nonblank/non-NA output from each method for every sample, while permitting per-cell unclassified values where the method explicitly leaves a cell unresolved and recording coverage rates.

The three declared exclusions remain clean skips, not weakly annotated outputs. A stale artifact that lacks either method must fail validation and be recomputed or fail closed; it must never be accepted as a completed Stage 4 run. The same strict contract applies to Joanito, Stephenson, and CombinedPBMC.

### Stage 5 — pass-scope violations remain

The canonical `src/5_run_benchmark_methods/1_submit_hpc_array.sh` currently
emits `BENCHMARK_*` runtime markers even when `--pass uncorrected|corrected`
is active. Its `worker_env` also includes `BENCHMARK_MANIFEST` alongside
`ANALYSIS_MANIFEST`, and `matrix_watchdog.sh` exports both variables for retry
arrays. These violate the pass-scoped batch contract.

The fix is to use a single pass-scoped manifest variable for canonical matrix
workers and retries (`ANALYSIS_MANIFEST`), reject an accidental
`BENCHMARK_MANIFEST` in pass mode, and emit pass-scoped runtime markers (for
example `BATCH_EFFECT_*` for batch passes and `BENCHMARK_*` only for ordinary
benchmark mode). Preserve the ordinary `watchdog_main.sh` contract only for
callers that explicitly use ordinary benchmark mode; it must not leak into the
matrix batch path. The uncorrected batch input resolver must use the canonical
`batch_effect_uncorrected` view for every selected cohort, including
CombinedPBMC.

The matrix wrapper currently sources shared helpers before assigning
`ANALYSIS_PASS`, so a conditional watchdog-status path can be initialized to
the ordinary benchmark directory. Assign/export the pass namespace before
helper initialization or make the helper resolve it at submission time.
Batch-mode status files, retry manifests, scheduler markers, and log prefixes
must all remain pass-scoped; only the ordinary compatibility path may use
benchmark-named runtime fields.

### Documentation and regression drift

`notebooks/dataset_onboarding/README.md` still names the deleted
`1_submit_batch_effect_stage.sh`. It must point to the canonical Stage 3
dispatcher and document the twelve-cohort batch-effect scope, nine
dual-annotation cohorts, three auto-annotation exclusions, CombinedPBMC's
canonical uncorrected view, and explicit selection-manifest usage.

The current focused submitter tests pass, but they do not fully defend the new
requirements: Myocardial force propagation, exact twelve-cohort selection,
CombinedPBMC view migration, full Stephenson scope, strict dual
HiTME/scATOMIC completeness, and absence of `BENCHMARK_MANIFEST`/`BENCHMARK_*`
markers in batch mode need targeted regressions.

## Implementation phases

### 1. Repair contracts and pass scoping locally

1. Add shared HiTME/scATOMIC required-column and coverage checks to `annotation_contract.py` for Feather and h5ad artifacts.
2. Make `3.1_merge_annotations.py`, `3.2_merge_worker.sh`, and the Stage 4 final/watchdog validation call the same strict contract after merge and before sync.
3. Change the annotation worker so a runnable cohort cannot turn a failed method into a successful all-NA checkpoint. Keep per-cell NA coverage statistics for genuinely unclassified cells.
4. Validate existing run-owned feathers before any worker skip; an invalid old feather is not a valid checkpoint.
5. Repair Stage 2 force propagation for `myocardial_counts`; add a focused stub
   test proving `--force` reaches the reconstruction helper.
6. Remove the legacy `batch_effect_analysis` token from the canonical
   preprocessor allowed-view sets, h5ad contracts, loaders, and tests. Rename
   CombinedPBMC's registry view/output to `batch_effect_uncorrected`; do not
   preserve an alias.
7. Remove `BENCHMARK_MANIFEST` and all `BENCHMARK_*` scheduler/status markers
   from canonical batch worker/retry exports, enforce the pass-scoped variable
   contract in workers, and switch batch runtime/status/log paths to the pass
   namespace. Ensure `ANALYSIS_PASS` is available before shared-helper
   initialization or resolve pass paths dynamically; ordinary compatibility
   callers retain their old contract only outside pass mode.
8. Extend the batch-candidate evidence registry/loader to cover Joanito,
   Stephenson, and CombinedPBMC, using CombinedPBMC's canonical
   `batch_effect_uncorrected` view and the full Stephenson batch-effect view.
9. Update focused tests for strict dual annotation, exact exclusions, exact
   twelve-row selection, CombinedPBMC legacy-name rejection and migration,
   pass-scoped exports/markers, and preserved ordinary compatibility behavior.
10. Update `datasets.json`, onboarding README commands, and scope documentation
    for the CombinedPBMC view rename. Do not change technical batch-column
    values in this phase; the view-schema migration is required and is not a
    scientific batch-column decision.

### 2. Local and `_debug` verification before HPC

Run only local/stub checks first:

- shell syntax, Python compilation, R parsing, profile JSON parsing, and focused Stage 2–5 tests;
- a deterministic contract fixture with complete dual annotations, missing HiTME, missing scATOMIC, all-NA method output, duplicate keys, and incomplete key coverage;
- an `_debug` Joanito five-sample annotation/preprocessing smoke path, without launching a full onboarding cohort;
- confirm no test or source path reintroduces the deleted Stage 3 submitter or the two-variable batch manifest export.

No full-cohort scheduler work occurs during this implementation/verification phase.

### 3. Reconcile the stopped work only after explicit execution approval

At the later user-controlled execution boundary, reconcile existing durable-gate
manifests, remote status, tmux/process identity, and scheduler state. Preserve
all discrepancy evidence. Do not relaunch an ambiguous gate. Establish one
clean Bamboo source revision after the local fixes are accepted. Freeze that
revision before scheduler submission: the technical batch-column values remain
unchanged until the scientific decision phase, while the required
CombinedPBMC view-schema migration and any dependency-lock update needed by the
validated annotation environment must already be complete.

### 4. Stage 2 prerequisite gates

Before Stage 3, validate or produce the three dataset-specific prerequisites through the canonical Stage 2 dispatcher:

```bash
./src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Myocardial_infarction \
  --steps myocardial_counts \
  --force

./src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets CombinedPBMC,Joanito
```

The first command is used when an intentional Myocardial recomputation is required. The second command selects `gongsharma_cap` and `combinedpbmc` automatically for CombinedPBMC and the `joanito` hook for Joanito; valid artifacts are skipped by checksum. Each command is a separate durable Stage 2 gate if both are needed. Do not call dataset-specific hooks directly.

Every Stage 2 gate must source `src/slurm_config.sh`, run from `$HOME/ECODA_paper`, and complete one durable wait, one terminal inspect with every emitted scheduler ID, and required reviewer approval before Stage 3 consumes its outputs.

### 5. Stage 3 uncorrected preprocessing gate

Create one immutable twelve-row selection TSV on Bamboo. Use these exact rows:

```text
Alzheimer	batch_effect_uncorrected
Breast_cancer	batch_effect_uncorrected
Covid19_PBMC	batch_effect_uncorrected
Kidney_KPMP	batch_effect_uncorrected
Myocardial_infarction	batch_effect_uncorrected
Diabetes	batch_effect_uncorrected
Lupus_PBMC	batch_effect_uncorrected
Lung	batch_effect_uncorrected
Parkinson	batch_effect_uncorrected
Joanito	batch_effect_uncorrected
Stephenson	batch_effect_uncorrected
CombinedPBMC	batch_effect_uncorrected
```

Use the canonical dispatcher through a new durable gate:

```bash
./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --selection-file "$BATCH_UNCORRECTED_SELECTION"
```

Acceptance: exactly twelve selection rows, one array over pending rows, no corrected output submission, all selected h5ads pass the repaired h5ad contract and checksum, Stephenson uses its full declared batch-effect view rather than `benchmark_analysis`, and the terminal audit includes the original array, every OOM retry array, and watchdog IDs.

### 6. Stage 4 uncorrected dual annotation gate

Use the same twelve-row selection TSV with the canonical Stage 4 submitter:

```bash
./src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  --selection-file "$BATCH_UNCORRECTED_SELECTION"
```

Expected result: nine runnable rows (`Breast_cancer`, `Covid19_PBMC`, `Kidney_KPMP`, `Lupus_PBMC`, `Lung`, `Myocardial_infarction`, `Joanito`, `Stephenson`, and `CombinedPBMC`) and three explicit `SKIP_NOT_SUITABLE` records (`Alzheimer`, `Diabetes`, and `Parkinson`). Preparation, annotation, and merge remain separate arrays; independent datasets run concurrently; no selected row is corrected.

Acceptance: every runnable dataset has complete HiTME and scATOMIC required columns and per-sample method coverage, complete `(Sample, barcode)` keys, an atomic merged h5ad, a checksum, and a passing final contract. Breast must not be accepted on the basis of HiTME-only output. The three exclusions have no fabricated annotation artifact.

### 7. Stage 5 uncorrected matrix and evidence gate

After reviewed Stage 4 completion, run only the twelve batch-effect datasets in the uncorrected-compatible pass with the fixed batch method suite:

```text
prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,pilotgm,qot
```

Use the canonical Stage 5 wrapper with `--pass uncorrected`, an explicit
twelve-dataset list, and no ordinary benchmark selection. Every selected row
resolves `batch_effect_uncorrected`, including CombinedPBMC after the required
view migration. Submit one array per method across all selected datasets, then
gate all method watchdogs together; only pseudobulk and composition depend on
the pseudobulk preparation watchdog. Canonical workers and matrix retries
receive only `ANALYSIS_MANIFEST` in pass mode.

Run `notebooks/dataset_onboarding/build_batch_candidate_evidence.R` as a
separate durable evidence gate against the uncorrected h5ads and pass-qualified
method root. Update it to read all twelve candidate registries and the
canonical CombinedPBMC uncorrected view, and produce one per-cohort CSV plus
`batch_candidate_review.csv`; malformed artifacts and sample-order mismatches
fail closed. Evidence includes candidate completeness, levels/samples per
level, biological NMI, marginal/joint PERMANOVA with Holm adjustment, and
constant/sample-unique/perfect-confounding warnings.

Excluded cohorts have no annotation-dependent method rows; methods that do not
require automated cell-type labels still run when their input contracts pass.
Any unavailable method applicability is recorded explicitly, never silently
imputed or treated as a pass.

### 8. Scientific batch-column decision checkpoint

Review the same candidate classes for all twelve cohorts. The nine onboarding candidate lists remain those declared in `notebooks/dataset_onboarding/dataset_specs.py`:

- `Alzheimer`: `assay`, `tissue_type`, `PMI`;
- `Breast_cancer`: `assay`, `sequencing_platform`, `sample_preservation_method`, `suspension_dissociation_time`, `suspension_dissociation_reagent`;
- `Covid19_PBMC`: `Single cell sequencing platform`, `City`, `datasets`, `Sample type`;
- `Diabetes`: `batch_integration`, `dataset`, `design`, `assay`;
- `Kidney_KPMP`: `experiment`, `library`, `tissue_type`, `region.l1`, `region.l2`, `assay`;
- `Lupus_PBMC`: `batch_cov`, `Processing_Cohort`, `ind_cov_batch_cov`;
- `Lung`: `dataset`, `study`, `platform`, `assay`, `origin_fine`;
- `Myocardial_infarction`: `batch`, `sampleType`, `dissociation_s1`;
- `Parkinson`: `Brain_bank`, `assay`, `tissue_type`, `PMI`, `RIN`.

Add the established batch-effect cohorts to the same evidence table:

- `Joanito`: `seqtec`, `Site` (with `seqtec` derived by the Stage 2 preparation hook);
- `Stephenson`: `Site` from the full declared batch-effect view;
- `CombinedPBMC`: `batch` from the combined multi-source view.

Reject candidates that are missing/incomplete, constant, sample-unique, or perfectly confounded with biology. Prefer candidates with replication across samples and a defensible technical meaning; do not choose by p-value alone. The existing `Joanito.seqtec`, `Stephenson.Site`, and `CombinedPBMC.batch` values must be re-evaluated through this same checkpoint. A scientific owner explicitly confirms one technical column per cohort, or records that the cohort cannot support corrected analysis.

Only after that explicit decision may `datasets.json` be edited for any changed value. Existing values may remain unchanged if re-confirmed. Every change must preserve all other dataset/view fields and pass the registry regression. No corrected preprocessing or correction model may run before the confirmed batch column is present in the h5ad and the configuration.

### 9. Corrected pass (later, separately authorized)

For only cohorts with an explicitly confirmed batch column:

1. Run Stage 3 with an explicit corrected selection manifest. Joanito and Stephenson have declared `batch_effect_corrected` views; the nine onboarding cohorts use their declared corrected views. CombinedPBMC currently has no `batch_effect_corrected` view, so corrected CombinedPBMC processing is blocked unless a separate, explicit `datasets.json` view-contract decision adds one. Do not fabricate that view.
2. Run Stage 4 for the corrected view for the nine auto-annotatable cohorts that have a corrected view, or use a validated two-view run when paired annotation provenance requires it. Verify paired `(Sample, barcode)` identity against the uncorrected view.
3. Run Stage 5 with `--pass corrected` and the same explicitly confirmed datasets/methods. Verify Harmony/PCA keys, batch-only model arguments, corrected pseudobulk (`blind=FALSE`, confirmed batch, `correct_batch=TRUE`, `~ 1`), row alignment, checksums, and NAS destinations.
4. Keep corrected artifacts in the corrected pass root; never overwrite uncorrected artifacts.

## Acceptance checklist

- Stage 2 runs the canonical Myocardial, Joanito, and CombinedPBMC prerequisite hooks with real force semantics and the required CombinedPBMC dependency.
- Stage 3 selection contains exactly twelve batch-effect datasets and only
  `batch_effect_uncorrected` views.
- Stage 4 skips exactly the three declared dual-annotation exclusions and rejects any of the other nine runnable datasets missing either HiTME or scATOMIC output.
- Stage 5 batch workers/retries use the pass-scoped manifest contract and do not emit ordinary benchmark markers in batch mode.
- Candidate evidence is reproducible, sample-order aligned, and retained as the decision record for all twelve cohorts.
- No biological label enters preprocessing or correction.
- No corrected run starts while any selected `columns.batch` remains unconfirmed; CombinedPBMC remains blocked until a corrected-view contract exists.
- No commit or push occurs before explicit user publication confirmation; all implementation changes remain uncommitted and unpushed until then.
