# Dataset onboarding — uncorrected batch-effect evidence gate

The nine Joodaki et al. 2025 cohorts plus Joanito, Stephenson, and the
derived CombinedPBMC cohort are registered from full-file metadata audits.
The authoritative sources are:

- `datasets.json` for inputs, subsets, views, output names, and activation flags;
- `notebooks/dataset_onboarding/dataset_specs.py` for user-confirmed roles,
  candidate technical columns, decision notes, and annotation provenance;
- regenerated `<key>_meta.json` audit files for observed counts, hierarchy,
  conflict warnings, and gate evidence.

Provisional root drafts are archived only for provenance. They are not registry
inputs.

## Canonical user-confirmed roles

| Key | Sample | Label | Low tier | High tier | Annotation source |
|---|---|---|---|---|---|
| Alzheimer | `donor_id` | `Cognitive status` | `Subclass` | `Supertype` | author / author |
| Breast_cancer | `sample_id` | `disease` | `broad_cell_type` | `author_cell_type` | author / author |
| Covid19_PBMC | `sampleID` | `CoVID-19 severity` | `majorType` | `celltype` | author / author |
| Diabetes | `donor_id` | `disease` | `cell_type` | `cell_type_reannotatedIntegrated` | author / author |
| Kidney_KPMP | `specimen` | `condition.l1` | `subclass.l1` | `subclass.l3` | author / author |
| Lung | `sample` | `disease` | `ann_coarse` | `ann_fine` | author / author |
| Lupus_PBMC | `sampleID` | `Status` | `layer1` | `layer2` | HiTME / HiTME |
| Myocardial_infarction | `orig_ident` | `patient_group` | `cell_type` | `cell_subtype` | author / author |
| Parkinson | `donor_id` | `disease` | `cell_type` | Leiden res-5 | author / Leiden |
| Joanito | `sample.ID` | `sample.origin` | `cell.type` | `cell.type_new` | author / derived |
| Stephenson | `Sample` | `Status` | `initial_clustering` | `full_clustering` | author / author |
| CombinedPBMC | `Sample` | `cond` | `layer1` | `layer2` | HiTME / HiTME |

The previous heuristic choice, stable-field conflicts, and aggregation warnings
remain in each audit. They explain the decision; they do not silently replace
the declared role. Missing IDs, standardized-ID collisions, missing labels, and
failed declared author hierarchies remain hard failures. HiTME and Leiden
columns are produced-output roles and remain pending until processed h5ad
evidence validates them.

## Full-file audit workflow

Source files are staged to `${HPC_SCRATCH_DIR}/_downloads/` on Bamboo and audit
metadata/subsets are written to its `subsets/` directory:

```bash
cd "${HOME}/ECODA_paper"
./notebooks/dataset_onboarding/run_subset_hpc.sh
```

The worker applies a spec's `subset_vars` before sample and annotation audits.
It records pre-filter and post-filter cell/sample counts plus the exact filter
expression. Diagnostic subsets are for reports only; they cannot promote a
failed gate.

Pull evidence and compare it with the registry:

```bash
mkdir -p data/new_dataset_checks/subsets
rsync -avP bamboo:scratch/ECODA_paper/_downloads/subsets/ \
  data/new_dataset_checks/subsets/
pixi run python notebooks/dataset_onboarding/_debug_validation.py \
  --metadata-only \
  --registry-audit-dir data/new_dataset_checks/subsets \
  --config datasets.json
```

For derived annotations, pass the processed output root. The validator writes
and then consumes `<key>_postprocess_gate.json`; it never treats raw-audit
absence as processed annotation evidence:

```bash
pixi run python notebooks/dataset_onboarding/_debug_validation.py \
  --registry-audit-dir data/new_dataset_checks/subsets \
  --processed-artifact-dir "$HPC_SCRATCH_DIR/batch_effect/uncorrected" \
  --config datasets.json
```

## Lung 10x registry subset

Lung is the exact categorical subset
`{"platform": {"values": ["10x"], "op": "in"}}`, applied before all role and
hierarchy audits. The platform mask is authoritative if it differs from the
assay-name mask containing `10x`; the audit records both masks row-by-row.
Observed filtered cell/sample totals replace paper expectations only after the
full-file run verifies them. The `ann_coarse → ann_fine` hierarchy must pass.
No copied stage-2 dataset is created.
The primary biological label is `disease`, with observed values including
`normal`, `chronic obstructive pulmonary disease`, `lung adenocarcinoma`,
`squamous cell lung carcinoma`, and `non-small cell lung carcinoma`.
`origin` remains available as a secondary tumor-versus-normal metadata field.

## Two-pass registry views

Every selected dataset declares `views.batch_effect_uncorrected` and
`views.batch_effect_corrected`. The logical view key is never
`batch_effect_analysis`. Existing non-Combined output basenames retain their
historical `batch_effect_analysis_uncorrected` filename component; this is only
a filename component, not an accepted view.

CombinedPBMC is the explicit basename migration:

```text
raw:       CombinedPBMC/data/combined_pbmc.h5ad
uncorrected: combined_pbmc_batch_effect_uncorrected_ECODAprocessed.h5ad
corrected:   combined_pbmc_batch_effect_corrected_ECODAprocessed.h5ad
```

Its low/high roles are `layer1`/`layer2`. Joanito's roles are
`cell.type`/`cell.type_new`; the latter is derived from `cell.type` and `iCMS`.
The batch-effect analysis consumes only each registry entry's high-resolution
role. Stephenson's batch subset is its full declared batch-effect view; its
only candidate is `Site`.

The Stage 3 and Stage 4 canonical dispatchers consume one immutable,
headerless two-column selection file. Its adjacent `MD5/SIZE/PATH` sidecar is
part of the contract. Exact mode validates these twelve rows in order and
rejects legacy or corrected views before any array submission:

```bash
./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --selection-file "$BATCH_UNCORRECTED_SELECTION" \
  --exact-batch-selection --force
```

The Stage 3 dispatcher supports two non-submitting recovery forms:
`--sync-only RUN_ID` revalidates the run-owned manifests, scheduler records,
watchdog state, h5ads, and checksums before completing an interrupted tail;
numeric/CSV `--sync-only` IDs require the original `--datasets` and
`--view/--views` selection and never infer a broader scope.

The uncorrected view always uses `batch_key=Sample`; corrected execution
remains gated on a confirmed technical column. All nine new cohorts retain
`columns.batch: null`; Joanito (`seqtec`) and Stephenson (`Site`) remain the
existing confirmed values. The Stage 4 exact run records three a-priori
`SKIP_NOT_SUITABLE` rows (`Alzheimer`, `Diabetes`, `Parkinson`) and runs nine
runnable datasets.

## Durable stage execution

Stage 2 derives the Myocardial counts, Joanito metadata/debug artifact, and
CombinedPBMC raw input. The latter accepts the old raw basename only through a
checksum/content-validated one-time rename to `combined_pbmc.h5ad`.

Stage 3 is the canonical
`src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` array. Stage 4 is
`src/4_cell_type_annotation/1_submit_onboarding_stage.sh`; it records the
three explicit unsuitable-cohort skips and runs all nine eligible datasets in
parallel. Both stages use OOM-only watchdog retries and fail closed on stale
checksums, malformed manifests, or invalid content.

Launch each full-cohort stage only through the checked-in
`durable-hpc-gate-ecoda` profile. Arm one unbounded durable `wait`, run one
terminal `inspect` with every emitted array/retry/watchdog ID, and obtain Luna
Max reviewer approval before starting the dependent stage. Do not poll
`squeue` or `sacct` from the agent session.

## Pass-specific preprocessing

`batch_effect_uncorrected` runs one hvg2000 pass with `Sample`, raw PCA,
neighbors, and Leiden, with no Harmony. `batch_effect_corrected` requires a
confirmed batch, selects HVGs by that technical column, and computes raw PCA
plus Harmony neighbors/Leiden. Exact keys are:

```text
X_pca_batch_effect_uncorrected_hvg2000
leiden_res_<r>_batch_effect_uncorrected_hvg2000
X_pca_batch_effect_corrected_hvg2000
X_pca_harmony_batch_effect_corrected_hvg2000
leiden_res_<r>_batch_effect_corrected_hvg2000_harmony
```

The fixed suite uses resolutions `0.1, 0.4, 2, 5, 20, 50`; reported ECODA
uses res-2, except Parkinson's configured res-5 tier.

Stage 5 matrix selection rows are `DATASET<TAB>VIEW<TAB>SCOPE`; in batch
mode both `VIEW` and `SCOPE` must equal the selected
`batch_effect_<pass>` view. The third field is pass scope, not a method;
methods come only from the fixed `--methods` list below.

## Fixed batch method suite

The fixed Stage 5 pass list is:

`prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`.

PILOT-GM-VAE is intentionally excluded from batch-effect passes because its
training cost is disproportionate for the large batch-effect cohorts.

It consumes only the configured high-resolution role for author ECODA
composition. `layer2` and `scATOMIC_pred` remain separate standardized
annotation-derived features. For the three a-priori exclusions,
`ECODA_HiTME_HR_layer2` and `ECODA_scATOMIC_HR` are recorded as unavailable
with reason `not_suitable_for_auto_annotation`; all other applicable outputs
remain required.

`Avg_PCA`, MOFA, scITD, scPoli, GloProp, cell-frequency-only baselines,
LR ECODA, top-variable-cell-type variants, zero-imputation screens, and
parameter screens are excluded. ECODA defaults are exactly
`clr_zero_impute_method="counts_all"` and `clr_zero_impute_num=0.5`: add 0.5
to every count before CLR. The shuffled baseline shares features and uses
deterministic label shuffling; labels remain evaluation-only.

Corrected execution is deferred until the evidence decision. It uses batch-only
LMM correction and corrected pseudobulk/Harmony/MrVI settings only after
technical columns are explicitly confirmed.

## Evidence checkpoint

The uncorrected method bundles feed
`notebooks/dataset_onboarding/submit_batch_candidate_evidence.sh`, which runs
the strict twelve-row evidence builder and synchronizes only validated outputs:
one CSV per cohort, `batch_candidate_review.csv`, their MD5 sidecars, and the
checksum manifest. The builder records candidate completeness, levels/samples
per level, NMI with biology, marginal/joint PERMANOVA $R^2$/Holm-adjusted
p-values, explicit method availability/reasons, and
constant/sample-unique/perfect-confounding warnings. It uses 999 permutations
and strict sample-order checks.

At this checkpoint all nine new `columns.batch` values remain `null`. Stop for
one explicit user-confirmed technical column per cohort. Only then run the
corrected pass.

Render local reports with:

```bash
for report in notebooks/dataset_onboarding/dataset_check_*.qmd; do
  quarto render "${report}"
done
```