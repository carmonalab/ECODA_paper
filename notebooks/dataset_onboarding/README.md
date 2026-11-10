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
- Production participation requires both the appropriate activation flag and
  configured view; view existence alone does not enable a cohort.
- Bamboo is the active HPC host after maintenance. The user reports that the
  latest data was synchronized back; Yggdrasil is stale/unmaintained. Verify
  selected Bamboo sources and provenance before reuse. Local subset mirrors
  under `data/` remain diagnostic and may be stale.
- `dataset_specs.py`, rendered QMDs, and `<key>_meta.json` files provide
  onboarding evidence and notes. They do not override `datasets.json`.
- The [active reprocessing plan](../../.agents/plans/20261108-confounding-subsets-benchmark-reprocessing-plan.md)
  owns approved decisions, exact raw-ID selection, execution scope, and applied
  status. Applied progress is summarized in the dated checkpoint below;
  registry edits or historical checksums alone do not establish release approval.
- Only benchmark views supply comparative timings and retain their configured
  CPU/GPU model constraints. Both batch views permit any compatible hardware
  model. Batch MrVI defaults to CPU; future GPU work requires explicit
  `--gpu-policy any` or `default`, with pinned CUDA/runtime validation and any
  explicit VRAM request preserved. This does not alter current jobs or require
  recomputing valid results.

## Cohorts examined

The onboarding set comprises nine newly added Joodaki et al. cohorts, plus the
existing `Joanito` and `Stephenson` cohorts. `CombinedPBMC` is documented
separately below as disabled legacy context, not as an active cohort.

The source/paper-reported major cell-type counts below preserve the useful
cohort-level annotation context while remaining distinct from configured
high-resolution labels.

| Cohort | General description | Configured batch variable(s) | Source/paper-reported major cell types | Current batch-analysis decision |
|---|---|---|---:|---|
| `Alzheimer` | Brain multiome cohort with limited assay diversity. | `assay` | 18 | Included after Reference exclusion |
| `Breast_cancer` | Breast-tumor cohort with assay and sequencing imbalances. | `assay + suspension_dissociation_time` | 10 | Included with two retained v3/dissociation batches |
| `Covid19_PBMC` | PBMC cohort with healthy controls and COVID samples across time after symptom onset. | `datasets` | 10 | Included with acute/control and d02/d09/d13 restrictions |
| `Diabetes` | Pancreatic cohort spanning multiple studies and disease groups. | `dataset` | 13 | Excluded because most batches have samples from only one condition. After removal of confounded batches, too few samples remain and the remaining two batches are both heavily imbalanced (3 normal and 12 type 2 diabetes left). |
| `Kidney_KPMP_full` | Combined single-cell/single-nucleus KPMP cohort. | `suspension_type` | 14 | Included after Papilla exclusion; no sex correction |
| `Lupus_PBMC` | PBMC cohort with case/healthy metadata and technical batch candidates. | `batch_cov` | 11 | Included with AbFlare-4 and selected pools 20–22 |
| `Lung` | Lung atlas with platform, tissue, disease, and assay metadata. | `dataset` | 12 | Included with six retained batches and the existing COPD exclusion |
| `Joanito` | Cohort with sample-origin, sequencing, and site metadata. | `seqtec` | Unavailable (not documented in onboarding specs) | Included after LymphNode and Tumor-2 exclusions |
| `Stephenson` | Cohort with biological status and site/batch metadata. | `Site` | Unavailable (not documented in onboarding specs) | Included after Sanger and sample-quality exclusions |
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

The approved batch-analysis set is the eight included cohorts above; it
excludes `CombinedPBMC` and all other disabled or diagnostic entries:
`Alzheimer`, `Breast_cancer`, `Covid19_PBMC`,
`Kidney_KPMP_full`, `Lupus_PBMC`, `Lung`, `Joanito`, and `Stephenson`.
The benchmark set has 11 enabled production datasets; only `Stephenson`
overlaps these batch cohorts. Alzheimer, Joanito and Kidney full are batch-only.
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

The sample-level matrices in `plots/batch_effect_bio_vs_batch/` motivated the
approved restrictions below. Apply fixed selections, then validate balance;
do not substitute an automatic pruning algorithm. Each batch must represent
at least two biological groups. Every represented batch/group combination
needs at least two samples and at least 20% of its batch for two represented
groups, or 10% for three. Equality passes. Absent global groups are reported
structural zeros. More than three represented groups needs an explicit decision.

| Cohort | Approved restriction | Expected samples per batch view |
|---|---|---:|
| Alzheimer | Remove all three `Reference` donor-by-assay samples; retain Dementia/No dementia | 101 |
| Breast_cancer | Keep only `10x 3' v3` with `12 hour` or `15 minute` dissociation | 33 |
| Covid19_PBMC | Keep d02 (22), d09 (12), d13 (10), intersected with literal `control` or finite numeric sampling day <=30; d07 excluded | 44 |
| Joanito | Remove `Tumor-2` and retain the LymphNode exclusion | 124 |
| Kidney_KPMP_full | Remove `region.l1 == "Papilla"`; retain unknown region and correct only suspension type | 80 |
| Lupus_PBMC | Keep `dmx_AbFlare-4` (11), `dmx_YS-JY-20_pool3/4`, `dmx_YS-JY-21_pool1/2`, `dmx_YS-JY-22_pool5/6` (83 combined) | 94 |
| Lung | Keep Goveia 2020, Lambrechts 6149v2/6653, and Leader beads_cite/digest-deadcell_cite/5p_v1_CD2 batches; preserve 10x/lung/COPD restrictions | 80 |
| Stephenson | Remove Sanger and the four excluded sample IDs documented below | 110 |

The active plan records exact Lung literals and all expected counts. These
are expectations to reconcile against Bamboo and >=500-cell/QC eligibility,
not proof of regenerated artifacts. Preserve whole samples: resolve approved
sample membership from audited sample-level batch assignments, then keep all
cells of those raw IDs. Cell-level batch predicates can otherwise split samples.

- **COVID verification:** the existing matrix's checksum-bound sample metadata
  already satisfies the control/day<=30 rule. Intersecting the three retained
  batches gives 44 samples; each retained batch includes controls. Unknown,
  malformed, blank, missing and nonfinite day values fail the temporal rule.
- **Stephenson exception:** after removing `MH8919226` and `MH8919227`
  (Healthy/Ncl samples with no annotated T cells), Ncl has 43 Covid / 10 Healthy
  (18.87% Healthy). The user explicitly retains this below-20% batch. Preserve
  the earlier exclusions `BGCV10_CV0198` and `MH8919230`. The Ncl-only benchmark
  target has 53 samples and requires all affected method/screen results refreshed.
- **Papilla rationale:** the four Papilla specimens were nuclei from deceased
  donors. Single-cell biopsies have literal `region.l1 == "NA"` and
  `tissue == "kidney"`; do not infer Cortex/Medulla for them. The lack of
  expected papillary tissue in biopsies motivates the exclusion, without
  treating unknown biopsy anatomy as a known compartment.
- **Diabetes:** disable the cohort entirely. Removing confounded batches leaves
  only three normal and 12 type 2 diabetes samples; the earlier two-group
  exclusion rule is superseded. Myocardial infarction, Parkinson and
  CombinedPBMC remain excluded.

These are prespecified scientific cohort decisions. Biological labels remain
ground truth for evaluation and approved cohort selection; they never enter
normalization, HVG selection, embeddings or correction-model covariates.

## Approved evaluation and repair work

- Kidney reports `condition.l1` and literal `tissue` separately on the same 80
  samples/distances: kidney (45), cortex of kidney (23), renal medulla (12).
  The scalar primary label remains supported; dataset-level
  `columns.evaluation_labels` declares both axes without a second physical
  dataset, H5AD, feature/distance calculation or method fit. Initially both
  evaluations count equally as independent summary entries.
  Tissue recovery is confounded with suspension and must be described as
  recorded tissue/sampling-category recovery. Region is diagnostic only.
- Default benchmark and batch summaries use five metrics: ANOSIM, adjusted
  kNN-3 modularity, mean Ward/PAM ARI, silhouette and LISI separation. Min–max
  normalize within dataset/view/evaluation x metric, using neutral 0.5 for
  constant ranges. Rank by the descending arithmetic mean of the five equally
  weighted scores; include the two added right-hand metric columns. LISI is
  already higher-is-better separation and must not be inverted again.
  Missing/nonfinite scores fail release.
- Corrected batch comparisons will distinguish unsupervised
  `ecoda_unsupervised_harmony` (Harmony -> Leiden -> CLR) from
  `ecoda_unsupervised_harmony_limma` (the same CLR followed once by technical-only
  limma). Both use HVG2000/Leiden resolution 2 and share the same clusters;
  uncorrected comparisons retain non-Harmony PCA -> Leiden -> CLR.
- Bassez's HR metadata derivative is complete: valid `cellSubType` was retained
  and 50,809 missing values were filled from `cellType` across 74,948 cells and
  26 samples, yielding 40 categories. Source inputs and protected counts,
  embeddings and other metadata are unchanged. The approved HR-dependent
  method work is complete, with selected releases reviewed and NAS-validated.
- The whole-sample blacklist pseudobulk variant is complete across all 11
  benchmark datasets. The exclusion source is the unchanged frozen
  human aux groups G1/S, G2/M, Mito, Ribo, TCR. Filter before fitting and select
  2,000 eligible genes afterward. Add no sex/chromosome mask or exemption.

Report intake now requires an explicit reviewed control and its checksum:
benchmark uses `ECODA_BENCHMARK_RELEASE_RDS`/`ECODA_BENCHMARK_RELEASE_MD5`, and
batch uses `ECODA_BATCH_RELEASE_CONTROL`/`ECODA_BATCH_RELEASE_CONTROL_MD5`.
The finite ledger selects exact aggregate keys or new standalone results;
missing required new results cannot fall back to old artifacts. Figure 2A's
12-method default pool and Supp fig 15's 21-method extended pool include NULL;
MDS means use the corresponding pool. Preserve stored primary NULL labels,
while Kidney tissue NULL records seed 124. Batch grids contain 360 uncorrected
or 405 corrected values (nine evaluation entries, five metrics).

Operational checkpoint (**2026-11-10**): original P0–P5 processing and selected method
releases are complete, scientifically reviewed and NAS-validated; no method
GPU work remains queued. Saved-label and zero-imputation screens, PCA rendering
and cell-depth diagnostics are complete and reviewed. The thin projection and
both batch reports have passed review. The parameter screen is being completed
with supported PILOT-GM-VAE HR-1000/2000/3000 and LR-2000 results, reusing valid
outputs and fitting only missing or current-input-invalid settings. Final
screen/workbook QA and the finite NAS backup remain pending. Filenames and numbering are
preserved. The [active plan](../../.agents/plans/20261108-confounding-subsets-benchmark-reprocessing-plan.md)
owns exact execution/release status and preserved failed incidents.
