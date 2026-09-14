# Dataset onboarding

## Purpose

The onboarding notebooks investigate newly added cohorts before analysis. They
describe dataset structure, metadata, cell-type columns, candidate technical
variables, biological/technical confounding, and possible inclusion or
subsetting decisions. They support the user's scientific decision; they are
not a production job manifest.

## Sources and authority

- `datasets.json` is the production registry for inputs, views, subsets,
  activation flags, and configured cell-type columns.
- The full-cohort data on Bamboo/HPC is authoritative. Local subset mirrors
  under `data/` are diagnostic and may be stale.
- `dataset_specs.py`, rendered QMDs, and `<key>_meta.json` files provide
  onboarding evidence and notes. They do not override `datasets.json`.

## Cohorts examined

The onboarding set comprises nine newly added Joodaki et al. cohorts, plus the
existing `Joanito` and `Stephenson` cohorts. `CombinedPBMC` is documented
separately below as disabled legacy context, not as an active cohort.

The source/paper-reported major cell-type counts below preserve the useful
cohort-level annotation context while remaining distinct from configured
high-resolution labels.

| Cohort | General description | Configured batch variable(s) | Source/paper-reported major cell types | Current batch-analysis decision |
|---|---|---|---:|---|
| `Alzheimer` | Brain multiome cohort with limited assay diversity. | `assay + sex` | 18 | Included |
| `Breast_cancer` | Breast-tumor cohort with assay and sequencing imbalances. | `assay + sequencing_platform + suspension_dissociation_time` | 10 | Included |
| `Covid19_PBMC` | PBMC cohort with healthy controls and COVID samples across time after symptom onset. | `datasets` | 10 | Included with an approved acute/control subset |
| `Diabetes` | Pancreatic cohort spanning multiple studies and disease groups. | `dataset` | 13 | Included after confounded disease groups are excluded |
| `Kidney_KPMP_full` | Combined single-cell/single-nucleus KPMP cohort. | `suspension_type + sex` | 14 | Included |
| `Lupus_PBMC` | PBMC cohort with case/healthy metadata and technical batch candidates. | `batch_cov` | 11 | Included |
| `Lung` | Lung atlas with platform, tissue, disease, and assay metadata. | `dataset` | 12 | Included after the configured COPD exclusion |
| `Joanito` | Cohort with sample-origin, sequencing, and site metadata. | `seqtec` | Unavailable (not documented in onboarding specs) | Included after the configured LymphNode exclusion |
| `Stephenson` | Cohort with biological status and site/batch metadata. | `Site` | Unavailable (not documented in onboarding specs) | Included |
| `Myocardial_infarction` | Small myocardial-infarction cohort with biological/technical confounding concerns. | — | 11 | Excluded |
| `Parkinson` | Brain multi-tissue cohort with weak overall batch-analysis signal. | — | 11 | Excluded |
| `Kidney_KPMP` | Preserved legacy Kidney key. | — | — | Disabled; not interchangeable with `Kidney_KPMP_full` |
| `_debug` | Small non-production Joanito verification fixture. | — | 10 | Diagnostic only |

The counts in the `Source/paper-reported major cell types` column are
source/paper-reported major cell-type counts. They are not configured
high-resolution unique-label counts: configured high-resolution columns (for
example `Supertype`, `author_cell_type`, or `celltype`) can contain many more
unique labels. `Joanito` and `Stephenson` are marked unavailable because no
verified source/paper-reported count is documented in the onboarding specs;
no count is inferred.

### Disabled legacy context: `CombinedPBMC`

`CombinedPBMC` was artificially composed from existing cohorts. Its biological
groups were confounded with batch/cohort, so those groups were not biologically
identifiable. It is retained only as stale legacy code/artefact context; it is
disabled, excluded from the active/current-final cohorts, and must not be
re-enabled.

The current final batch-analysis set is the nine included cohorts above; it
excludes `CombinedPBMC` and all other disabled or diagnostic entries:
`Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Diabetes`,
`Kidney_KPMP_full`, `Lupus_PBMC`, `Lung`, `Joanito`, and `Stephenson`.
The exact processing rows and artifact handling belong in the implementation
plan, not in this README.

## Cell-type columns used for batch analysis

Batch-effect views use the configured source/author cell-type columns already
present in the cohort data. The user-approved batch-analysis decision is not
to run HiTME or scATOMIC Pipeline 4 for these views.

| Cohort | Low-resolution column | High-resolution column |
|---|---|---|
| `Alzheimer` | `Subclass` | `Supertype` |
| `Breast_cancer` | `broad_cell_type` | `author_cell_type` |
| `Covid19_PBMC` | `majorType` | `celltype` |
| `Diabetes` | `cell_type` | `cell_type_reannotatedIntegrated` |
| `Kidney_KPMP_full` | `subclass.l1` | `subclass.l3` |
| `Lupus_PBMC` | `layer1` | `louvain` |
| `Lung` | `ann_coarse` | `ann_fine` |
| `Joanito` | `cell.type` | `cell.type_new` |
| `Stephenson` | `initial_clustering` | `full_clustering` |

## Confounding and subset decisions

Onboarding showed that some biological conditions are confounded with batch
or study structure. The registry records the resulting decisions:

- `Diabetes`: exclude `endocrine pancreas disorder` and `type 1 diabetes
  mellitus`, whose biological groups are strongly tied to study structure.
- `Joanito`: exclude `sample.origin == "LymphNode"`.
- `Lung`: retain `platform == "10x"` and `tissue == "lung"`, and exclude
  `disease == "chronic obstructive pulmonary disease"`.
- `Covid19_PBMC`: the approved final rule retains the sampling-day literal
  `control` or a finite numeric value `<= 30`; unknown, malformed, blank, and
  missing values are excluded. The current categorical/string encoding
  requires an explicit registry/helper correction before this is final.
- `Myocardial_infarction`, `Parkinson`, and `CombinedPBMC` are excluded from
  the current batch-analysis decision for the documented confounding, sample,
  or signal concerns.

Wholly batch-collinear biological groups were excluded rather than treated as
identifiable biological signal.
When a biological group was perfectly collinear with a technical batch/cohort, that group was excluded rather than protected or corrected; this is a scientific subset decision.

These are scientific registry decisions. Implementation-specific pipeline
selection, targeted rows, output paths, and artifact reuse belong in the
implementation plan.