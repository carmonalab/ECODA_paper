# Dataset onboarding — PILOT-GM-VAE study (Phase 5)

The nine Joodaki et al. 2025 cohorts are registered from full-file metadata
audits. The authoritative sources are:

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
| Lung | `sample` | `origin` | `ann_coarse` | `ann_fine` | author / author |
| Lupus_PBMC | `sampleID` | `Status` | `layer1` | `layer2` | HiTME / HiTME |
| Myocardial_infarction | `orig_ident` | `patient_group` | `cell_type` | `cell_subtype` | author / author |
| Parkinson | `donor_id` | `disease` | `cell_type` | Leiden res-5 | author / Leiden |

The previous heuristic choice, stable-field conflicts, and aggregation warnings
remain in each audit. They explain the decision; they do not silently replace
the declared role. Missing IDs, standardized-ID collisions, missing labels,
missing retained metadata, and failed declared author hierarchies remain hard
failures. HiTME and Leiden columns are produced-output roles and remain pending
until processed h5ad evidence validates them.

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

## Two-pass registry views

Every new cohort, Joanito, and Stephenson exposes exactly:

- `views.batch_effect_uncorrected`;
- `views.batch_effect_corrected`.

Both views use the same input and subset and distinct output names:

```text
<stem>_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad
<stem>_batch_effect_analysis_corrected_ECODAprocessed.h5ad
```

All nine new cohorts use `use_for_benchmark: false`,
`use_for_batch_effect: true`, and `columns.batch: null`. The uncorrected view
is runnable and always uses `batch_key=Sample`. The corrected view fails closed
until a confirmed technical column is written. Joanito (`seqtec`) and
Stephenson (`Site`) are already confirmed. Parkinson's high tier is
view-specific: raw res-5 in the uncorrected view and Harmony-qualified res-5
in the corrected view. Python and R loaders merge optional view-level column
overrides over the dataset-level columns.

No active batch-effect artifact, key, manifest, cache, log, or result variable
uses a benchmark-named identifier.

## Durable stage execution

The uncorrected onboarding preprocessing stage is one parallel array over the
nine new cohorts plus Joanito and Stephenson. Launch
`src/3_scrnaseq_preprocessing/1_submit_batch_effect_stage.sh` through the
checked-in `durable-hpc-gate-ecoda` profile; it submits the worker array and a
compute-node watchdog, then returns scheduler IDs without waiting. The
watchdog retries only `OUT_OF_MEMORY` rows with doubled memory up to its
configured ceiling, validates every pass-qualified h5ad, and synchronizes
completed output to the canonical NAS root.

Use one durable gate per pipeline stage. Arm one unbounded durable `wait`, run
one terminal `inspect` with the array, every watchdog retry array, and watchdog
IDs emitted in the wrapper log, then obtain Luna Max reviewer approval before
starting annotation. Do not poll `squeue` or `sacct` from the agent session or
launch the next stage before the reviewed gate is terminal.

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

## Fixed batch method suite

Configured high-resolution methods are:

`ECODA_authors_HR`, applicable `ECODA_HiTME_HR_layer2`, applicable
`ECODA_scATOMIC_HR`, `ECODA_seuratres_2`, `ECODA_authors_HR_NULL`,
Pseudobulk, GloScope, PILOT, PILOT-GM-VAE, MrVI, and QOT.

Avg_PCA, MOFA, scITD, scPoli, GloProp, cell-frequency-only baselines,
LR ECODA, top-variable-cell-type variants, zero-imputation screens, and
parameter screens are excluded. ECODA defaults are exactly
`clr_zero_impute_method="counts_all"` and `clr_zero_impute_num=0.5`: add 0.5
to every count before CLR. The shuffled baseline shares features and uses
deterministic label shuffling; labels remain evaluation-only.

Corrected ECODA uses batch-only LMM correction and exact row-zero-sum
recentring. Corrected pseudobulk uses `blind=FALSE`, the confirmed batch,
`correct_batch=TRUE`, and `~ 1`. GloScope/PILOT/PILOT-GM-VAE/QOT use Harmony
PCA; MrVI receives only the confirmed technical `batch_key`.

## Evidence checkpoint

The uncorrected method bundles feed
`notebooks/dataset_onboarding/build_batch_candidate_evidence.R`, which emits
one CSV per cohort and `batch_candidate_review.csv` with completeness,
levels/samples per level, NMI with biology, marginal/joint PERMANOVA
$R^2$/Holm-adjusted p-values, and constant/sample-unique/perfect-confounding
warnings. It uses 999 deterministic permutations and strict sample-order
checks.

At this checkpoint all nine new `columns.batch` values remain `null`. Stop for
one explicit user-confirmed technical column per cohort. Only then run the
corrected pass, verify paired identities, pass-specific checksums/NAS sync,
exact keys, zero-sum CLR rows, batch-only pseudobulk settings, and native MrVI
batch arguments.

Render local reports with:

```bash
for report in notebooks/dataset_onboarding/dataset_check_*.qmd; do
  quarto render "${report}"
done
```