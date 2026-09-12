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

## Active registry roles (`datasets.json`)

| Key | Sample | Label | Low tier | High tier | Usable high-tier categories | Annotation source |
|---|---|---|---|---|---:|---|
| Alzheimer | `donor_id` | `Cognitive status` | `Subclass` | `Supertype` | 131 | author / author |
| Breast_cancer | `sample_id` | `disease` | `broad_cell_type` | `author_cell_type` | 58 | author / author |
| Covid19_PBMC | `sampleID` | `CoVID-19 severity` | `majorType` | `celltype` | 58 | author / author |
| Diabetes | `donor_id` | `disease` | `cell_type` | `cell_type_reannotatedIntegrated` | 20 | author / author |
| Kidney_KPMP_full | `specimen` | `condition.l1` | `subclass.l1` | `subclass.l3` | pending | author / author |
| Lung | `sample` | `disease` | `ann_coarse` | `ann_fine` | 44 | author / author |
| Lupus_PBMC | `sampleID` | `Status` | `layer1` | `louvain` | 24 | HiTME / derived |
| Myocardial_infarction | `orig_ident` | `patient_group` | `cell_type` | `cell_subtype` | 33 | author / author |
| Parkinson | `donor_id` | `disease` | `cell_type` | `cell_type` | 11 | author / author |
| Joanito | `sample.ID` | `sample.origin` | `cell.type` | `cell.type_new` | 12 | author / derived |
| Stephenson | `Sample` | `Status` | `initial_clustering` | `full_clustering` | 50 | author / author |
| CombinedPBMC | `Sample` | `cond` | `layer1` | `layer2` | 40 | HiTME / HiTME |

> Counts are distinct non-missing values in the configured High tier source
> column, not method-output feature-column counts. `CombinedPBMC` has 40
> valid `layer2` labels; its composition metadata contains 41 unique values
> because one missing/sentinel level accounts for 77,813 cells and is excluded. `Kidney_KPMP_full`
> remains pending because no active full-cohort artifact is available; the
> legacy `Kidney_KPMP` count of 67 is not carried over. Historical Lupus
> `layer2`/Parkinson Leiden counts are not used; current roles are `louvain`
> (24) and `cell_type` (11).

> **Legacy compatibility note — not active:** `Kidney_KPMP` is retained only
> for the preserved historical onboarding QMD and reproducibility of its
> existing artifacts. Its registry flags remain disabled; it is excluded from
> active onboarding, batch-effect order, selection files, and worker dispatch.

The previous heuristic choice, stable-field conflicts, and aggregation warnings
remain in each audit. They explain the decision; they do not silently replace
the declared role. Missing IDs, standardized-ID collisions, missing labels, and
failed declared author hierarchies remain hard failures.

The table follows the active `datasets.json` High tier roles. Historical
Lupus `layer2` and Parkinson Leiden selections do not replace the declared
roles; produced-output roles are published only when corresponding processed
evidence is available.

## Kidney_KPMP_full combined cohort

`Kidney_KPMP_full` is the active Kidney cohort and represents a combined
single-cell/single-nucleus (sc/sn) source. Its canonical staged input is
`Kidney_KPMP_full.h5ad`; source identity and verification details come from
the authoritative source catalog and registry rather than from this
documentation. No source URL or paper/legacy count is assumed here. Expected
cell and independent-unit counts remain unset until the full-file audit
establishes and the user accepts the source contract.

The declared roles are `specimen` (sample), `condition.l1` (biological
condition), `subclass.l1` (low-resolution annotation), and `subclass.l3`
(high-resolution annotation), with author/author provenance. The observed
`suspension_type` field is an explicit batch-effect candidate and must be
reported separately for its single-cell and single-nucleus coverage. Missing
declared roles or required metadata are hard audit findings; the report does
not silently substitute another column.

The active batch-effect registry and its authoritative twelve-row order use
`Kidney_KPMP_full` in place of the disabled legacy key. The active uncorrected
view is `batch_effect_uncorrected`: preprocessing uses `Sample` and no
technical correction covariate. Its Stage 5 batch run uses the fixed suite
`prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot`.

The active order is:

```text
Alzheimer, Breast_cancer, Covid19_PBMC, Kidney_KPMP_full,
Myocardial_infarction, Diabetes, Lupus_PBMC, Lung, Parkinson,
Joanito, Stephenson, CombinedPBMC
```

The Kidney onboarding report
`dataset_check_Kidney_KPMP_full.qmd` is diagnostic information only. It reads
the source or a diagnostic subset, displays audit findings, and cannot edit
`datasets.json`, selection files, or submit jobs. A source, modality,
metadata, hierarchy, or count discrepancy stops the onboarding decision for
explicit user direction; a report rendering never authorizes automated
registry or compute changes.


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
remains gated on a confirmed technical batch definition: either one scalar
key or an ordered, nonempty list of independent keys. A corrected
`columns.batch: null` or malformed key list is rejected before a run root,
manifest, scheduler ID, or worker state exists. All nine new cohorts retain
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
neighbors, and Leiden, with no Harmony. `batch_effect_corrected` accepts either
one scalar batch key or an ordered, nonempty list of independent additive
batch factors (a one-element list has the same execution semantics as its
scalar key), then computes raw PCA plus Harmony neighbors/Leiden. Corrected
`columns.batch: null` is rejected; the raw reader shape remains unchanged
(`string`, list, or `null`) and corrected execution-boundary normalization
rejects non-string, empty, duplicated, or whitespace-only keys, including a
key equal to the biological label or standardized `Sample`. After view
subsetting, `Sample` standardization, and the 500-cell filter, every
configured key is validated over the full selected cell metadata before HVG,
PCA, or Harmony: it must exist, contain no missing/blank value or
case-insensitive `NA`, `nan`, `None`, `<NA>`, `n/a`, `null`, or `Unknown`
sentinel, be constant within `Sample`, and have at least two levels.
Disconnected, rank-deficient, near-unique, or otherwise non-estimable
correction designs fail closed; no `Unknown` imputation is allowed. Biological
labels are evaluation-only. Exact keys are:

```text
X_pca_batch_effect_uncorrected_hvg2000
leiden_res_<r>_batch_effect_uncorrected_hvg2000
X_pca_batch_effect_corrected_hvg2000
X_pca_harmony_batch_effect_corrected_hvg2000
leiden_res_<r>_batch_effect_corrected_hvg2000_harmony
```

The fixed suite uses resolutions `0.1, 0.4, 2, 5, 20, 50`; reported ECODA
uses res-2. Parkinson's uncorrected High tier is the configured `cell_type`;
its corrected view uses the explicit res-5 Harmony column.

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

Corrected execution is deferred until the evidence decision and starts only
after the confirmed scalar or ordered key list passes the full-cell and
estimability checks above. Corrected ECODA uses additive per-key random
intercepts; corrected Pseudobulk uses a scalar composite for batch-only
removal with DESeq2 `~ 1`, `blind=FALSE`, and `correct_batch=TRUE`. No
biological label enters correction or design.

## Corrected-mode batch contract

The corrected method matrix is exactly:

```text
ECODA_authors_HR
ECODA_seuratres_2
ECODA_authors_HR_NULL
Pseudobulk
GloScope
PILOT
MrVI
QOT
```

No other method is in the corrected matrix. In particular,
`pilot-gm-vae` (`PILOT-GM-VAE`) is explicitly excluded from corrected work
and must not be selected, validated, or used as a fallback. The historical
fixed Stage 5 list above and ordinary benchmark method references remain
unchanged.

### Ordered factors, validation, and identity

Readers preserve the raw `columns.batch` JSON shape (`string`, ordered list,
or `null`). At the corrected execution boundary, a scalar is one ordered key;
a nonempty list is an ordered set of independent additive batch factors, and a
one-element list has the same direct-column semantics as its scalar key.
Corrected `null` is rejected, as are non-string, empty, duplicated, or
whitespace-only key names and keys that equal the biological label or
standardized `Sample`. These failures occur before any run root, manifest,
scheduler ID, or worker state is created.

After view subsetting, `Sample` standardization, and the 500-cell filter, every
configured key is validated across the full selected cell metadata before
first-observation metadata collapse or any HVG, PCA, Harmony, CLR, DESeq2,
pseudobulk, or MRVI setup. Each key must exist, contain no actual missing or
blank values and no case-insensitive `NA`, `nan`, `None`, `<NA>`, `n/a`,
`null`, or `Unknown` sentinel, remain constant within standardized `Sample`,
and have at least two observed levels. Disconnected, rank-deficient,
near-unique, or otherwise non-estimable correction designs fail closed; no
`Unknown` imputation is allowed. (`Unknown` in a cell-type annotation column
is outside this batch-column rule.) Biological labels remain evaluation-only
and never enter preprocessing, correction, design, or MRVI setup.

For two or more keys, the exact cross-language composite encoding is
`ecoda_batch_composite_v1`. Composite construction is used only for two or
more keys; a scalar or one-key list remains a direct column. Each nonmissing
scalar category value is canonicalized to typed Unicode text: character or
factor values use `s:` followed by the exact label; booleans use `b:true` or
`b:false`; signed integers use `i:` followed by base-10 decimal; and finite
IEEE-754 binary64 values use `f64:` followed by their 16-lowercase-hex-digit
big-endian bit pattern. Other, list, date, and object values are rejected.

The UTF-8 bytes of each key name and canonical value are encoded as lowercase
hexadecimal. Each token is exactly
`ecoda_batch_composite_v1|<key-count>|<pair-1>;<pair-2>;...`, with each pair
exactly `<key-byte-length>:<key-hex>,<value-byte-length>:<value-hex>`.
Lengths count raw UTF-8 bytes rather than characters, every hex field has two
characters per byte, and delimiters never occur inside hex. Configured key
order is preserved in every token; categorical levels may be sorted by raw
UTF-8 token bytes only when category metadata is constructed. Missing values
are never coerced into a category.

The reserved in-memory observation name is
`__ecoda_batch_combined_v1`; assert that it is absent before construction and
absent again before writing any corrected H5AD. Record the encoding version,
ordered source keys, per-key levels, composite-level count, sample-constancy
result, correction formula/mode, and key-set/configuration fingerprint in
run-owned metadata.

The fingerprint input is exact bytes, with no JSON/object serialization,
locale dependence, or map iteration:
`ecoda_batch_contract_v1\0` followed in this fixed order by
`field("encoding","ecoda_batch_composite_v1")`,
`field("keys",ordered_key_vector_v1)`,
`field("scalarization",scalarization_id)`,
`field("method",method_id)`, and `field("model",model_id)`.
`field(name,value)` is
`<name-byte-length>:<lowercase-hex-UTF-8(name)>,<value-byte-length>:<lowercase-hex-UTF-8(value)>;`,
with lengths counting raw UTF-8 bytes. `ordered_key_vector_v1` is
`key-count|<key-byte-length>:<lowercase-hex-UTF-8(key)>;...` in configured
order. `scalarization_id` is exactly `direct_v1` for a scalar or one-key
direct column, or `composite_v1` for two or more keys. `method_id` is the
exact canonical corrected method or `preprocess`; `model_id` is one of
`hvg_composite_v1`, `harmony_native_list_v1`,
`ecoda_additive_random_intercepts_v1`, `pseudobulk_composite_v1`,
`mrvi_composite_v1`, or `embedding_consumer_harmony_v1`. Hash the resulting
bytes with SHA-256 and record lowercase hexadecimal output; source H5AD and
configuration checksums remain separate.

For example, keys `site`, `tech` with values `A`, `x` encode as
`ecoda_batch_composite_v1|2|4:73697465,3:733a41;4:74656368,3:733a78`.

### Corrected method representations

| Method/path | Corrected multi-key representation |
|---|---|
| `ECODA_authors_HR`, `ECODA_seuratres_2`, `ECODA_authors_HR_NULL` | Use additive per-key random intercepts in the CLR correction model. The null variant shuffles only evaluation labels after the same corrected features are produced. |
| Scanpy HVG | Build `__ecoda_batch_combined_v1` temporarily and pass that one scalar name as `batch_key`; delete it immediately after HVG ranking, including retry or error paths. |
| Harmony | Pass the original ordered batch-key list natively (`vars_use=list` or the equivalent supported API); do not scalarize it. Preserve the exact pass-qualified embedding and orientation handling. |
| `Pseudobulk` | Build the composite at sample-metadata level and pass one scalar composite through the DESeq2/limma batch-removal path. Corrected DESeq2 uses exactly `~ 1`, `blind=FALSE`, and `correct_batch=TRUE`, with batch-only removal and no biological design protection. |
| `GloScope`, `PILOT`, `QOT` | No method-local batch argument. Resolve only `X_pca_harmony_batch_effect_corrected_hvg2000`; a missing exact key is a hard error. |
| `MrVI` | Load all original batch columns, recreate the equivalent composite in memory, and pass only its one scalar name to `MRVI.setup_anndata(batch_key=...)`. Never pass a list and never persist the temporary column. |

The HVG, MRVI, and R pseudobulk builders may be separate implementations,
but they must produce identical composite golden-vector tokens and the same
key-set fingerprint for equivalent scalar, two-key, and three-key fixtures.
New corrected artifacts are reusable only when this ordered key identity,
encoding, scalarization policy, model policy, and fingerprint all match; a
checksum alone is not semantic proof of the configured batch factors.

These rules apply only to `batch_effect_corrected`. Historical
`batch_effect_uncorrected` rows and order, its `Sample`/no-Harmony behavior,
and ordinary benchmark paths and method references remain unchanged.

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
one explicit user-confirmed technical batch definition per cohort: either one
scalar key or an ordered, nonempty list. Only then run the corrected pass.

Render local reports with:

```bash
for report in notebooks/dataset_onboarding/dataset_check_*.qmd; do
  quarto render "${report}"
done
```