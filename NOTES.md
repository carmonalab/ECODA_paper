# ECODA Paper — Architectural Decisions & Methodology Notes

> **Document Purpose:**  
> This document serves as the central engineering and methodology log for the *ECODA* manuscript revisions and benchmarking extensions. It records the statistical rationale, architectural designs, diagnostic tool implementations, and open discussion points for the project.

---


### Breast cancer corrected batch policy

- **Decision:** Breast's production configuration uses the dataset-level
  `columns.batch` list exactly
  `["assay", "suspension_dissociation_time"]`. Corrected Stage 3 and corrected
  Stage 5 consume this same list; no `views.*.columns` override is active and
  Stage 5 cannot select a separate correction key.
- `sequencing_platform` remains in H5AD `obs` and exported sample metadata
  when present, but is not a configured Breast correction column.
- `disease` is a biological label and remains evaluation-only; it never enters
  filtering, preprocessing, Harmony, or corrected Stage 5 models.
- The previous Breast corrected Stage 3 output was generated under the
  superseded column contract. The regenerated output
  `BreastCncr_processed_batch_effect_analysis_corrected_assay_dissociation_ECODAprocessed.h5ad`
  completed under the global pair. Its Stage 3 gate was terminally inspected
  and reviewed; a separate read-only semantic contract check remains before
  corrected Stage 5 release.
- **Historical rank evidence (retained, non-authoritative):** The 165-sample
  validator report showed the full additive design at rank `8/10` with
  residual degrees of freedom `157`. The observed dependencies are exact: all
  141 `10x 3' v3` samples use NovaSeq; the 24 `10x 3' v2` samples use HiSeq
  3000/4000; and all 11 `unknown` dissociation-time samples use HiSeq 4000.
  The pairwise ranks were `assay + sequencing_platform = 3/4`,
  `sequencing_platform + suspension_dissociation_time = 8/9`, and
  `assay + suspension_dissociation_time = 8/8`. This historical evidence
  supports the selected dataset-level pair; it does not establish a separate
  Stage 5 policy.
- **Historical implementation interpretation (retained, non-authoritative):**
  The prior report described an intercept-preserving design equivalent to
  `~ 1 + batch_key_1 + batch_key_2` after internal aliasing and noted that no
  `assay__sequencing_platform` composite was introduced. These statements
  document the superseded analysis interpretation only; they do not establish
  a separate Stage 5 key or replace the dataset-level configuration.
- The historical durable validator evidence is
  `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/_ecoda_logs/stage5_eight_corrected_ygg_parallel_20260915T131349Z/BREAST_BATCH_RANK_DIAGNOSIS.tsv`.
  The old eight-dataset gate validated seven dataset-level consumer rows and
  failed only Breast under the obsolete three-key design; no target method
  array was emitted. That gate is historical evidence only. The current
  regenerated Breast Stage 3 gate is the reviewed predecessor for the
  confirmed seven-method Breast Stage 5 scope.


> **Historical-note boundary.** The donor-only Alzheimer snapshot, majority-vote
> policy, `lme4` correction description, and pre-change final-lane descriptions
> below are retained for provenance and are explicitly labelled historical;
> they do not override this status.


## 1. Dataset Onboarding & Quality Control Pipeline

### Executive Summary
- Standardized onboarding check notebooks (`dataset_check_<Name>.qmd`) were developed across all target cohorts (9 new cohorts + `_debug` subset).
- Diagnostic workflow operates strictly on unintegrated raw counts to establish ground-truth data characteristics before any pipeline processing.
- Implemented **Normalized Mutual Information (NMI)** to detect metadata collinearity and confounding upfront.
- Implemented **Global Compositional Distance PERMANOVA** to partition inter-sample distance variance into Unique Biological, Unique Technical, Shared/Confounded, and Residual components.

### 1.1 Diagnostic Architecture & Verification Pillars
Each onboarding notebook adheres to the following sequence:
1. **Count Layer & Integer Sanity:** Verifies that `.X` or `layers['counts']` contains un-normalized, non-log-transformed integer counts with expected sparsity.
2. **Cohort Demographics & Harmonization:** Checks sample-level representation, cell-type label hierarchies (evaluating automated annotation suitability), and matches against paper metadata.
3. **Confounding & Metadata Collinearity Matrix (NMI):** Computes sample-level pairwise Normalized Mutual Information across all candidate biological and technical covariates.
4. **Unintegrated Baseline Embeddings:** Generates fresh PCA and UMAP embeddings directly from raw counts (2,000 HVGs, 50 PCs) without any pre-existing batch correction.
5. **Variance Partition & Distance PERMANOVA:** 
   - Feature-level variance partitioning (Sikkema et al. 2023 / HLCA) per cell type.
   - Global distance-based PERMANOVA on sample-by-sample Aitchison distance ($D = \|\text{clr}(x_i) - \text{clr}(x_j)\|_2$).
6. **Cell-Level Separation (LISI) & Automated Verdict:** Evaluates unintegrated cell-level LISI (mixing across batch keys vs. separation across biological labels) and generates an automated onboarding verdict.

### 1.2 Metadata Collinearity Diagnostic (Normalized Mutual Information — NMI)
- **Motivation:** Single-cell cohorts frequently exhibit partial or complete confounding between biological condition and technical variables (e.g. disease severity confounded with sequencing batch or tissue source). Simple regression over-attributes variance if collinearity is ignored.
- **Mathematical Definition:**
  $$\text{NMI}(X, Y) = \frac{2 \cdot I(X; Y)}{H(X) + H(Y)}$$
  where $I(X; Y)$ is mutual information and $H(X)$ is Shannon entropy. NMI evaluates non-linear statistical dependence at the sample level without parametric assumptions.
- **Thresholding:**
  - $\text{NMI} \approx 0.0$: Covariates are statistically independent.
  - $\text{NMI} \in [0.3, 0.7]$: Moderate overlap / partial confounding (requires multi-variable joint modeling).
  - $\text{NMI} > 0.70$: Severe collinearity (e.g. `Cognitive status` vs. `APOE4 status` in Alzheimer); batch correction against such covariates risks removing biological signal.

### 1.3 Inter-Sample Distance PERMANOVA & Variance Decomposition
- **Motivation:** While feature-level $R^2$ indicates which individual cell types fluctuate, patient stratification methods evaluate global pairwise sample distances.
- **Formulation (Anderson 2001 Marginal PERMANOVA):**
  - Given sample CLR composition matrix $X \in \mathbb{R}^{N \times P}$, compute pairwise Euclidean distance matrix $D$ (Aitchison distance).
  - Compute Gower's centered matrix $G = H (-0.5 D^2) H$, where $H = I - \frac{1}{N}\mathbf{1}\mathbf{1}^T$.
  - Fit full multi-variable model $D \sim \text{cov}_1 + \text{cov}_2 + \dots + \text{cov}_k$ and compute marginal sums of squares $\text{SS}_{\text{marginal}}(k)$ by dropping covariate $k$ from the joint projector.
- **Variance Partitioning:**
  $$\text{Total Distance Variance} = R^2_{\text{Unique Bio}} + R^2_{\text{Unique Tech}} + R^2_{\text{Shared (Confounded)}} + R^2_{\text{Residual}}$$
  - Evaluated via permutation Pseudo-$F$ tests ($B=999$ permutations) with Benjamini-Hochberg FDR correction.
  - Visualized via a dual-panel plot: (1) Global 100% variance decomposition stacked bar, and (2) Marginal $R^2$ bar chart with permutation $p$-values and significance codes.

### 1.4 Historical onboarding cohort sample-count comparison with PILOT-GM-VAE

The historical counts below come from the full-file onboarding audits in
`data/new_dataset_checks/subsets/*_meta.json`, using the sample column registered
in [`datasets.json`](datasets.json). The cell threshold is strict: samples with
fewer than 500 cells are dropped, while samples with exactly 500 cells are
retained.

| Dataset | Historical sample column | Historical samples | Dropped (<500) | Retained | PILOT-GM-VAE reported samples |
| :--- | :--- | ---: | ---: | ---: | ---: |
| Alzheimer | `donor_id` | 83 | 0 | 83 | 83 |
| Breast cancer | `sample_id` | 167$ | 2 | 165 | 126$ |
| Covid-19 PBMC | `sampleID` | 172 | 8 | 164 | 151 |
| Diabetes | `donor_id` | 56 | 0 | 56 | 52* |
| Kidney (KPMP) | `specimen` | 47 | 2 | 45 | 45 |
| Kidney (KPMP) full (sc/sn) | `specimen` | 88 | 4 | 84 | not reported |
| Myocardial infarction (MI-2) | `orig_ident` | 24 | 0 | 24 | 23 |
| Lung atlas | `sample` | 304† | 19 | 285 | 165‡ |
| Lupus PBMC | `sampleID` | 261 | 1 | 260 | 261 |
| Parkinson | `donor_id` | 97 | 1 | 96 | 97 |

The PILOT-GM-VAE column is the reported count from [Table 1 of the
study](https://academic.oup.com/bib/article/26/5/bbaf547/8287234#536377145);
the same dataset descriptions are preserved in the
[author-provided PDF](notebooks/dataset_onboarding/datasets.pdf). These
reported units are not necessarily identical to the current registry units:
PILOT uses donor/patient/sample terminology by cohort, whereas the current
table follows the configured column in `datasets.json`.

* `$` For Breast cancer, the onboarding audit finds 126 `donor_id` values and 167
  `sample_id` values. `sample_id` nests `donor_id` (up to four sample IDs per
  donor), while `accSample` is equivalent to `donor_id`. The PILOT-GM-VAE
  report says 126 donors, so the most likely explanation is that the authors
  used `donor_id` (or the equivalent `accSample`) rather than `sample_id`.
  This explains the 126-versus-167 discrepancy before the 500-cell filter.
  We found that 25 donors (mostly labelled as `disease` = `breast cancer`)
  had samples taken from the left and right breast. However, it was not clear
  whether both samples from the same donor were cancer samples, respectively,
  or whether one is `cancer` and the other `normal`. 4 `normal` donors had two
  samples, each, prepared with with different dissociation times. 4 `normal`
  donors had 3-4 samples, each, from left and right breast (1-2 samples from each side)
  prepared with with different `suspension_uuid`.
* `*` The PILOT-GM-VAE Table 1 lists 52 Diabetes samples. The author-provided
  PDF additionally says that four embryo samples were excluded as outliers, but
  the paper does not reconcile that statement with the 52-sample Table 1
  entry. If the exclusion is applied to those 52 samples, the effective count
  would be 48. The current onboarding file contains 56 `donor_id` values and
  has no `subset_vars` exclusion for those samples.
* `†` The lung dataset has 304 `sample` IDs,
  and 165 `donor_id` values. It has 304 nonzero donor/sample pairs; every
  sample maps to exactly one donor, 79 donors have multiple sample IDs, and
  the maximum is 16 sample IDs for one donor. The per-donor distribution is:
  86 donors with one sample, 56 with two, 4 with three, 17 with four, 1 with
  ten, and 1 with sixteen. The 19 sample IDs below 500 cells leave 285
  sample-level units. For donors with multiple samples, different samples come
  from different tissues, different conditions, storage (fresh or frozen),
  dissociation protocols, suspension types (cell or nuclei), or sequencing (3' vs 5').
* `‡` The PILOT-GM-VAE report gives 941,504 cells and 165 donors for Lung.
  Those values now match cell and `donor_id` totals
  exactly. The remaining 304-versus-165 discrepancy is therefore a unit
  choice: ECODA's configured column counts `sample` IDs, while PILOT reports
  donors.

### 1.5 Historical HPC-backed donor-to-sample audit for the batch-effect onboarding view

This is a historical, read-only audit snapshot. Its donor-only/configured-column
descriptions predate the strict `donor_id_assay` mapping in the current status
and remain for provenance; they are not the current Alzheimer sample identity.
The audit targeted the canonical `batch_effect_uncorrected` output on Bamboo. It used the configured `py-cuda13` runtime for HDF5/AnnData-backed metadata reads only; the larger cohort reads were dispatched through `srun` on Bamboo. Only `obs` metadata and categorical codes were read, and `.X`/count matrices were never materialized. Canonical HPC sample/donor counts below are post-500-cell-filter output units, whereas raw/local audit counts are pre-filter source metadata; the PILOT values remain publication-reported units.

| Cohort | Configured / sample-like column | Donor/patient column | Canonical HPC sample units | Donor/patient units | Exact donor→sample frequency | Observed explanation |
| :--- | :--- | :--- | ---: | ---: | :--- | :--- |
| [Alzheimer](https://pubmed.ncbi.nlm.nih.gov/39402379/) | `donor_id` (configured); `Specimen ID` (audit) | `donor_id` | 202 Specimen IDs; 1,395,601 cells | 83 | 53×2, 24×3, 6×4 | Every specimen maps to one donor; the configured registry sample remains `donor_id`. |
| [Breast cancer](https://pubmed.ncbi.nlm.nih.gov/37380767/) | `sample_id` (configured); `Sample` (canonical) | `donor_id` | 165 Samples; 713,691 cells | 124 | 90×1, 29×2, 3×3, 2×4 | 34 repeated donors contribute 75 samples; `Sample` maps to at most one donor. |
| [Covid-19 PBMC](https://pubmed.ncbi.nlm.nih.gov/33657410/) | `sampleID` (configured); `Sample` (canonical) | `PatientID` | 164 Samples; 991,227 cells | 149 | 138×1, 8×2, 2×3, 1×4 | 11 repeated patients contribute 26 retained samples; `Sample` maps to at most one patient. |
| [Diabetes](https://pubmed.ncbi.nlm.nih.gov/37697055/) | `donor_id` (configured); `dataset__design__sample` | `donor_id` | 56 Samples; 301,796 cells | 56 | 56×1 | One-to-one donor/sample view; no repeated donor/sample IDs. |
| [Kidney (KPMP)](https://pubmed.ncbi.nlm.nih.gov/37468583/) | `specimen` (configured/canonical) | `donor_id` | 45 specimens; 103,642 cells | 43 | 42×1, PRE019×3 | `specimen` maps to at most one donor; PRE019 is the sole repeated donor. |
| [Lung atlas](https://pubmed.ncbi.nlm.nih.gov/37291214/) | `sample` (configured); `Sample` (canonical; registry `platform=10x AND tissue=lung`) | `donor_id` | 285 Samples; 936,636 cells | 156 | 82×1, 53×2, 5×3, 14×4, 1×10, 1×16 | 74 repeated donors contribute 203 samples; `Sample` maps to at most one donor. |
| [Lupus PBMC](https://pubmed.ncbi.nlm.nih.gov/35389781/) | `sampleID` (configured); `Sample` (canonical) | unavailable | 260 Samples; 1,263,220 cells | unavailable | unavailable | No donor/patient/subject column exists, so repeated-donor status cannot be inferred. |
| [Myocardial infarction](https://pubmed.ncbi.nlm.nih.gov/35948637/) | `orig_ident` (configured); `Sample` (canonical) | `patient` | 24 Samples; 132,888 cells | 19 | 15×1, P2/P3/P15×2, P9×3 | Four repeated patients contribute 9 samples; current 24 `orig_ident` values differ from the QMD/PILOT 23-sample report. |
| [Parkinson](https://pubmed.ncbi.nlm.nih.gov/35513515/) ([PMID 39580497](https://pubmed.ncbi.nlm.nih.gov/39580497/)) | `donor_id` (configured); `Sample` (canonical) | `donor_id` | 96 Samples; 2,095,732 cells | 96 | 96×1 | One Sample per donor, but donor-level aggregation does not remove within-donor tissue heterogeneity. |

The publication comparison is the [PILOT-GM-VAE Table 1
report](https://academic.oup.com/bib/article/26/5/bbaf547/8287234#536377145).
Its donor/patient/sample terminology is a publication unit choice and is not
automatically identical to the configured column or to the canonical
`batch_effect_uncorrected` `Sample` field.

**Per-cohort audit details.** The existing Breast and Lung explanatory claims
and all Section 1.4 footnotes above remain verbatim; the details below add
the raw-versus-post-filter and canonical-unit clarification.

- **Alzheimer (historical donor-only/configured-column snapshot).** The audit unit is `Specimen ID` (202) against `donor_id` (83):
  53 donors have 2 specimens, 24 have 3, and 6 have 4, for 202 specimens
  total, and every specimen maps to one donor. All specimens have stable
  `Cognitive status`, `disease`, `sex`, and `tissue`/`tissue_type`; 21 donors
  have assay values spanning `10x 3' v3` and `10x multiome`. Metadata columns
  are `Cognitive status`, `Specimen ID`, `donor_id`, `assay`, `PMI`, and
  `tissue_type`.

- **Breast cancer.** The pre-filter audit has 167 `sample_id` values and 126
  `donor_id` values; two low-cell sample units are removed in the canonical
  output, leaving 165 `Sample` units and 124 donors across 713,691 cells.
  The exact donor-to-Sample frequencies are 90 donors with 1, 29 with 2, 3
  with 3, and 2 with 4; 34 donors repeat and account for 75 retained samples,
  while Sample-to-donor is at most one. Within those repeated donors,
  `tissue_location` differs for 25 donors (left/right patterns),
  `suspension_dissociation_time` differs for 5 (12 hour versus 15 minute for
  P09, P11, P20, and P23; 12 hour versus 6 hour for P24),
  `suspension_uuid` differs for all 34, and `sample_uuid` differs for 31.
  Disease is constant within donor and all preserved samples are fresh.
  Retained sample-level totals are disease: breast cancer 43, normal 122;
  `sample_source`: Baylor 103, MD Anderson 37, UCI 25; assay: `3'v2` 24,
  `3'v3` 141; and dissociation time: 3h 68, 2h 27, 15min 27, 12h 19, 4h 8,
  6h 5, unknown 11.

- **Covid-19 PBMC.** The raw/local audit has 172 `sampleID` values and 151
  `PatientID` values; 8 sample units below 500 cells are dropped, leaving 164
  canonical Samples and 149 patients across 991,227 cells. The retained
  donor-to-Sample frequencies are 138 patients with 1, 8 with 2, 2 with 3,
  and 1 with 4; 11 patients repeat and account for 26 retained samples, while
  Sample-to-Patient is at most one. The exact repeated-patient
  `Sample time`/sampling-day patterns are:

  | PatientID | Retained samples | `Sample time` / `Sampling day` |
  | :--- | ---: | :--- |
  | `P-M009` | 2 | convalescence day36/day36 |
  | `P-M026` | 3 | progression days7,11; convalescence22 |
  | `P-M041` | 2 | progression4; convalescence21 |
  | `P-M042` | 2 | progression10; convalescence27 |
  | `P-M043` | 2 | progression5; convalescence9 |
  | `P-M044` | 2 | progression6; convalescence19 |
  | `P-S032` | 2 | progression16/22 |
  | `P-S035` | 4 | progression9/9; convalescence14/14 |
  | `P-S036` | 3 | progression14; convalescence19/21 |
  | `P-S069` | 2 | progression3/15 |
  | `P-S070` | 2 | progression7/14 |

  `Sample time` differs for 7 repeated patients (`P-M026`, `P-M041`,
  `P-M042`, `P-M043`, `P-M044`, `P-S035`, `P-S036`); sampling day differs for
  10 (all above except `P-M009`). Within each repeated patient,
  `CoVID-19 severity`, `Outcome`, `Sex`, `City`, `datasets`, `Single cell
  sequencing platform`, and `Sample type` stay constant; no repeated-patient
  group changes them. Retained sample-level totals are severity: control 20,
  mild-moderate 67, severe-critical 77; `Sample time`: control 20,
  convalescence 84, progression 60; `Sample type`: fresh PBMC 83,
  frozen PBMC 81; platform: `10X5'` 136, `10X3'` 28. Metadata columns are
  `sampleID`/`Sample`, `PatientID`, `CoVID-19 severity`, `Sample time`,
  `Sampling day (Days after symptom onset)`, `Outcome`, `Sex`, `City`,
  `datasets`, `Single cell sequencing platform`, and `Sample type`.

- **Diabetes.** The canonical output has 56 `donor_id` values and 56
  `dataset__design__sample` values in a one-to-one mapping (56 donors × 1
  sample), with no repeated donor/sample IDs, across 301,796 cells. Ten donor
  IDs mix sex metadata in the audit; this warning remains. Metadata includes
  `disease` (endocrine pancreas disorder 12, normal 26, T1D 6, T2D 12),
  `dataset` (9 source-study labels), `design`, `batch_integration`, `assay`,
  `strain`, `diabetes_model`, and `age`. The PILOT Table 1 value of 52 versus
  the current 56 remains unexplained, including the author-provided
  embryo-exclusion note.

- **Kidney (KPMP).** The raw audit has 47 specimens and 45 donors (and 49
  library values); two raw specimen units are removed by the `<500` filter,
  leaving 45 canonical specimens and 43 donors across 103,642 cells. The
  exact donor-to-specimen frequency is 42 donors × 1 and the repeated donor
  `PRE019` × 3; the three records have experiment values `PREMIERE31`,
  `PREMIERE36`, and `PREMIERE37`, and library values `PRE019-4`,
  `PRE19-025`, and `PRE19-05`. Specimen-to-donor is at most one. Conditions
  are stable within donor; sample totals are `condition.l1`: Ref 20, CKD 14,
  AKI 11; and `condition.long`: Normal Reference 20, DKD 11, AKI 11, H-CKD 3.
  `experiment` and `library` are sample-specific technical identifiers;
  `10x 3' v3` assay and biopsy tissue occur for all samples.

- **Kidney_KPMP_full (combined sc/sn).** The user-confirmed CellxGene download
  passes the full-file audit with 304,652 cells and 88 `specimen` units.
  The audit also observes 67 `donor_id` values and 93 `library` values; the
  declared `specimen` role is the passing canonical unit, while donor/library
  alternatives remain audit evidence. `subclass.l1` and `subclass.l3` both
  pass the declared author-annotation gates, `condition.l1` is the biological
  label, and count sanity passes on sparse integer `raw.X`. The
  `suspension_type` modality split is 104,314 single-cell (`cell`) and 200,338
  single-nucleus (`nucleus`) cells. The uncorrected view uses `Sample` and no
  technical correction covariate; `suspension_type` and `sex` remain explicit
  batch candidates. This is a separate active cohort; the legacy
  `Kidney_KPMP` record and its historical counts remain unchanged.

- **Lupus PBMC.** The canonical output has 1,263,220 cells and 260 Sample IDs;
  the raw/local audit has 261 `sampleID` values and one is dropped by the
  `<500` filter. No donor, patient, or subject column exists in canonical
  `obs`, so repeated-donor status cannot be inferred and no repeated-donor
  conclusion is possible. Metadata includes `Status` (Case 162, Healthy 98),
  `SLE_status` (SLE 162, Healthy 98), `batch_cov` (23 levels),
  `Processing_Cohort` (1=38, 2=120, 3=26, 4=76), `sampleID`/`Sample`, `Age`,
  `Sex`, `pop_cov`, and `ind_cov_batch_cov`. The existing raw audit's
  sampleID stable-field conflict for exactly `1130_1130` and `1772_1772` is
  preserved.

- **Lung.** The canonical HPC view is the registry intersection
  `platform=10x AND tissue=lung`, not the local platform-only audit. It has
  285 Samples, 156 donors, and 936,636 cells; 74 donors repeat and account
  for 203 Samples, with exact donor-to-Sample frequencies 82×1, 53×2, 5×3,
  14×4, 1×10, and 1×16, and Sample-to-donor at most one. The 304 samples and
  165 donors in the existing source-subset note, with its exact distribution
  86×1, 56×2, 4×3, 17×4, 1×10, and 1×16 and 79 repeated donors, are
  pre-filter counts; 19 sample units below 500 cells leave the post-filter HPC
  285/156 view. In the current output, `origin` differs for 57 of 74
  repeated donors (mostly `normal_adjacent` plus `tumor_primary`);
  `origin_fine` has multiple non-missing values for 30 (including
  `tumor_edge`/`tumor_middle`); `assay` differs for 3 (`10x 3' v2` versus
  `10x 5' v1`); and `dataset` differs for 4. Disease is constant within
  donors; `platform=10x`, `tissue=lung`, and `suspension_type=cell`, with
  sex/study/tissue otherwise stable at sample level. Post-filter sample-level
  disease totals are lung adenocarcinoma 133, normal 71, squamous cell lung
  carcinoma 41, NSCLC 19, COPD 21; `origin` totals are `tumor_primary` 114,
  `normal` 92, `normal_adjacent` 76, `tumor_metastasis` 3; and assay totals
  are `3'v2` 256, `5'v1` 16, `3'v3` 8, `3'v1` 5. The metadata columns are
  `sample`/`Sample`, `donor_id`, `disease`, `origin`, `origin_fine`, `tissue`,
  `assay`, `dataset`, `study`, `platform`, `suspension_type`, and `sex`.
  The current local `Lung_meta.json` platform-only audit reports 336 samples
  and 193 donors because it omits registry `tissue=lung`; those values must
  not be conflated with the canonical batch view.

- **Myocardial infarction.** The canonical output has 24 `Sample`/`orig_ident`
  values and 19 patients across 132,888 cells. Patient-to-Sample frequency is
  15×1, 3×2 (`P2`, `P3`, `P15`), and 1×3 (`P9`): four patients repeat and
  account for 9 samples, while Sample-to-patient is at most one. `P2` has
  `patient_region_id` `RZ_BZ_P2`/`RZ_GT_P2`, `sampleType` border/remote,
  batch A/B, and cell types Fib/CM; `P3` has
  `RZ_BZ_P3`/`RZ_P3`, border/remote, and CM/Myeloid; `P15` has
  `GT_IZ_P15`/`IZ_P15`, sample types not mixed in the stored summary, and
  Fib/Myeloid; `P9` has `GT_IZ_P9`/`RZ_P9` plus a third sample,
  ischemic/remote, `patient_group` ischemic/myogenic, and CM/Endo.
  `batch` differs only for P2; `patient_region_id` differs for all four;
  `sampleType` differs for P2, P3, and P9; and `dissociation_s1` differs for
  all four. Sample-level `patient_group` totals are myogenic 13, fibrotic 5,
  ischemic 6. Metadata columns are `orig_ident`/`Sample`, `patient`,
  `patient_region_id`, `patient_group`, `sampleType`, `batch`,
  `dissociation_s1`, and `cell_type`. The existing QMD/PILOT 23-sample,
  19-patient statement conflicts with the current 24 `orig_ident` count and
  is retained as a flagged discrepancy.

- **Parkinson.** The canonical output has 96 `Sample` values and 96
  `donor_id` values in a one-to-one mapping (96 donors × 1 Sample) across
  2,095,732 cells; one raw donor below 500 cells is removed. There is one
  Sample per donor, but each donor has 2–5 tissue values, so donor-level
  aggregation is not the same as no within-donor heterogeneity: donor-to-tissue
  frequency is 4 donors × 2 tissues, 4 × 3, 17 × 4, and 71 × 5, and all 96
  donors span multiple tissue values. The five tissue totals are dorsal motor
  nucleus vagus 20, medial globus pallidus 22, prefrontal cortex 20,
  primary motor cortex 17, and primary visual cortex 17; disease is
  Parkinson 73 and normal 23; `Brain_bank` is Harvard 25, MSSM 23, UD 21,
  UM 27; and all assays are `10x 3' v3`. No separate sample/patient ID exists,
  so configured `donor_id` is the aggregation unit. Metadata columns are
  `donor_id`/`Sample`, `tissue`, `tissue_type`, `disease`, `Brain_bank`, `PMI`,
  `RIN`, and `assay`; the existing 96-donor/multi-tissue warning is retained.


### 1.6 Inclusion or exclusion of cohorts

Some cohorts were either completely or partially excluded from the analysis for the following reasons:
- No high resolution cell type annotations provided by the authors
- Too few samples
- Biological group completely confounded with batch (biological signal cannot be separated from batch signal)
- No clear biological signal


Detailed overview of cohorts. Can be used or biased by batch?

Included:
- [X] Alzheimer -> very few multiome
- [X] Breast -> imbalanced assay and sequencing_platform, but suspension_dissociation_time very mixed
- [X] Covid-19 -> use only accute COVID, i.e. <=30 days (above, no separation) (ECODA still shows separation on mds1/2 dims)
- [X] Diabetes -> heavily biased (only “normal” strongly mixed, with 9 datasets. “type 1 diabetes” only one dataset, “endocrine pancreas disorder” only one dataset, “type 2 diabetes” only two datasets) -> drop confounded batch/bio groups (endo + type 1)
- [X] Joanito -> seqtec mostly balanced (except a few lymphnode all 5’), tissue site imbalanced but very mixed for all -> drop lymphnode
- [X] Kidney_KPMP_full -> very well balanced
- [X] Lung -> copd (chronic obstructive pulmonary disease) only in one dataset (adams, which also makes up half of all normal samples. adams seems to cluster away the most) -> drop adams? or at least drop copd?
- [X] Lupus -> well mixed
- [X] Stephenson -> well mixed

Excluded:
- [  ] CombinedPBMC -> biased
- [  ] Kidney Cancer -> too few samples, no high resolution cell type annotation, two biological groups defined by presence/absence of one specific annotated cell type ("tumor cells")
- [  ] Myocardial -> very low batch effect but two bio conds mostly from same batch, third bio cond balanced -> ischemic confounded -> after exclusion too few samples
- [  ] Pancreas (PDAC) -> o high resolution cell type annotation, biological groups defined by presence/absence of one specific annotated cell type ("ductal cells") -> Separates by "Ductal 2" cells (which might be cancer cells)
- [  ] Parkinson -> quite well mixed but overall almost no signal -> show in appendix, not for ranking. Also no high resolution cell type annotation
- [  ] Kidney Cancer -> too few samples, no high resolution cell type annotation, biological groups defined by presence/absence of one specific annotated cell type


---

## 2. Batch Effects in Single-Cell Cohorts: Expression vs. Composition

### Executive Summary
- Single-cell technical noise operates through distinct physical mechanisms across modalities: **gene expression artifacts** (sequencing chemistry, depth, ambient RNA) vs. **compositional artifacts** (dissociation protocols, cryopreservation, cell-type destruction).
- Benchmark evaluations distinguish between expression-level batch keys and composition-level batch keys.
- Detailed cohort profiles and candidate metadata columns are documented in [`notebooks/dataset_onboarding/README.md`](notebooks/dataset_onboarding/README.md).

### 2.1 Modality-Specific Noise Mechanisms

| Modality | Primary Physical Batch Mechanisms | Example Candidate Metadata Columns |
| :--- | :--- | :--- |
| **Gene Expression** | Sequencing platform, library prep chemistry (10x 3' v2 vs v3, 5' vs 3'), flowcell/lane, sequencing depth, ambient RNA contamination | `assay`, `seqtec`, `sequencing_run`, `10x_chemistry`, `library_id` |
| **Cell Composition** | Enzymatic dissociation protocol, tissue digestion time, cold ischemia time, fresh vs. frozen / cryopreservation, FACS gating / cell enrichment | `tissue_type`, `PMI` (post-mortem interval), `dissociation_protocol`, `sample_preservation`, `enrichment` |
| **Cohort / Demographic** | Clinical collection site, hospital center, donor sex, age, ethnic background | `Site`, `Center`, `sex`, `Age`, `self_reported_ethnicity` |

---

## 3. Batch Effect Correction & Benchmark Strategy (`batch_effect_analysis.rmd`)

### Historical pre-change final-lane summary

> This subsection records the pre-change uncorrected mixed-source final-lane
> plan. Keep it for provenance; it is not current corrected-final
> compute-completion evidence.

Batch-effect analysis uses the canonical logical views
`batch_effect_uncorrected` and `batch_effect_corrected`. The final analysis
uses an artifact variant with `_final` H5AD/result names and separate
`uncorrected_final`/`corrected_final` roots; the historical legacy roots remain
unchanged.

The final mixed-source scope includes the existing uncorrected results for
`Alzheimer`, `Breast_cancer`, `Lupus_PBMC`, and `Stephenson`, newly regenerated
final views for `Covid19_PBMC`, `Diabetes`, `Joanito`, and `Lung`, and only the
missing targeted uncorrected Stage 5 rows for `Kidney_KPMP_full`. `CombinedPBMC`,
`Kidney_KPMP`, `Myocardial_infarction`, and `Parkinson` receive no new work.

The biological label is evaluation-only. It never enters filtering, HVG
selection, normalization, PCA, Harmony, CLR correction, pseudobulk design,
MrVI covariates, or any other model input.

Batch-effect views do not run Pipeline 4. They use the configured source/author
cell-type columns:

| Dataset | Low-resolution column | High-resolution column |
| :--- | :--- | :--- |
| Covid19_PBMC | `majorType` | `celltype` |
| Diabetes | `cell_type` | `cell_type_reannotatedIntegrated` |
| Joanito | `cell.type` | `cell.type_new` |
| Lung | `ann_coarse` | `ann_fine` |

The frozen cohorts use their already-produced uncorrected result artifacts
read-only and are absent from all new jobs and validator selections. Diabetes
also remains covered by its explicit automatic-annotation exemption, but the
batch workflow does not invoke annotation work for any target.

The approved Covid final subset is same-column filtering on
`Sampling day (Days after symptom onset)`: retain the literal `control` level
or a finite numeric value `<= 30`; exclude blank, unknown, malformed, and
missing values. The source column is categorical with string-valued levels, so
numeric-looking values require explicit parsing.

### 3.1 Fixed method contract

The batch suite contains exactly these configured high-resolution methods:

- ECODA configured author high tier (`ECODA_authors_HR`);
- ECODA Leiden resolution 2 (`ECODA_seuratres_2`);
- deterministic shuffled-label ECODA baseline (`ECODA_authors_HR_NULL`);
- Pseudobulk;
- GloScope;
- PILOT;
- MrVI;
- QOT.

The suite excludes `Avg_PCA`, MOFA, scITD, scPoli, GloProp,
cell-frequency-only baselines, LR ECODA, top-variable-cell-type variants,
zero-imputation screens, and all parameter screens. QOT belongs to the
Harmony-corrected expression group alongside GloScope, PILOT, and
PILOT-GM-VAE.

ECODA uses the exact default `clr_zero_impute_method="counts_all"` and
`clr_zero_impute_num=0.5`: add 0.5 to **every** count before the CLR
transformation, not only entries equal to zero. The shuffled baseline uses the
same features and deterministic label shuffling; labels remain evaluation-only.

### 3.2 Pass-specific preprocessing

`batch_effect_uncorrected` performs one hvg2000 pass with `batch_key=Sample`,
raw PCA, neighbors, and Leiden only. It emits no Harmony representation.

`batch_effect_corrected` requires a confirmed non-null top-level dataset
`columns.batch`; view-level column objects cannot replace or narrow that
dataset-level source. It performs one hvg2000 pass with HVGs selected by the
configured technical columns, computes raw PCA, then Harmony and
neighbors/Leiden on Harmony. For `Breast_cancer`, the exact correction list is
`["assay", "suspension_dissociation_time"]`; `sequencing_platform` remains
metadata when present but is not a correction column. The biological label is
never protected in correction.

Pass-qualified keys are literal and never fall back:

- `X_pca_batch_effect_uncorrected_hvg2000`;
- `leiden_res_<r>_batch_effect_uncorrected_hvg2000`;
- `X_pca_batch_effect_corrected_hvg2000`;
- `X_pca_harmony_batch_effect_corrected_hvg2000`;
- `leiden_res_<r>_batch_effect_corrected_hvg2000_harmony`.

All configured batch-effect datasets use the cell-type columns documented in
the final-scope table above. Disabled cohorts and the four frozen cohorts are
not part of new final computation.

### Historical modality-specific corrected-input descriptions (pre-limma)

The following correction paths are retained as historical implementation
notes. Current composition and pseudobulk correction uses separate-covariate
limma fixed effects; combined/artificial keys and `lme4` are prohibited.

- **Historical ECODA composition:** each CLR cell-type feature is fit with
  `lme4::lmer(y ~ 1 + (1 | batch), REML=TRUE)`. Subtract only the fitted batch
  random effect, then recenter every corrected row to an exact zero sum.
  Missing IDs, fewer than two batch levels, nonconvergence, and sample-order
  mismatches fail closed. The biological label is absent from the formula.
- **Pseudobulk expression:** uncorrected uses
  `blind=TRUE`, `batch_col=NULL`, `correct_batch=FALSE`, design `~ 1`.
  Corrected uses `blind=FALSE`, the confirmed technical `batch_col`,
  `correct_batch=TRUE`, design `~ 1`.
- **GloScope, PILOT, PILOT-GM-VAE, and QOT:** uncorrected resolves
  `X_pca_batch_effect_uncorrected_hvg2000`; corrected resolves
  `X_pca_harmony_batch_effect_corrected_hvg2000`. Missing exact keys are
  errors.
- **MrVI:** uncorrected receives no technical covariate; corrected receives
  only the confirmed technical column as native `batch_key`.

### 3.4 Historical pre-change final-lane artifact and evidence contract

The paths and mixed-source artifact descriptions in this subsection are
historical pre-change notes; they are not completion evidence for the current
corrected-final recovery.

Existing legacy pass artifacts remain under
`${HPC_SCRATCH_DIR}/batch_effect/<pass>/` and
`${NAS_TARGET_DIR}/batch_effect/<pass>/` at their current paths. New final
method bundles, distances, pseudobulks, execution logs, manifests, and plots
use separate `_final` analysis/plot roots and variant-qualified names. The
final analysis manifest records which rows reuse legacy results and which rows
use newly generated final artifacts.

The uncorrected evidence report records completeness, levels and samples per
candidate, NMI with biology, marginal and joint PERMANOVA $R^2$ and
Holm-adjusted p-values, and constant/sample-unique/perfect-confounding
warnings. It uses 999 deterministic permutations and strict sample-order
checks. It emits one CSV per cohort plus `batch_candidate_review.csv`.

The final run is targeted rather than a historical twelve-row rerun. The four
frozen cohorts are read-only analysis inputs and are absent from all new
selection, validation, and scheduler manifests. `Diabetes` receives no
annotation work; the final batch H5ADs retain its configured source columns.

### 3.5 No-leakage invariant

Technical correction is never allowed to protect or model a biological label.
Confounded technical variables remain documented warnings, not reasons to
silently change the confirmed sample or label roles. Missing IDs, collisions,
missing labels, empty derived annotations, invalid hierarchies, and missing
exact pass keys remain hard failures.


### Historical majority-vote technical batch-covariate policy (2026-09-13; superseded)

This is a historical policy snapshot retained for provenance. It is superseded
by the current strict `donor_id_assay`/limma contract and must not be read as
the current correction policy.

On 2026-09-13, the user approved `majority_v1` per configured sample for
technical batch covariates. No minimum winner fraction threshold is applied;
exact ties still fail because they have no majority, and real NA/blank/non-finite
values still fail.

Breast treats literal `unknown` as an ordinary
`suspension_dissociation_time` class only; it does not ignore other
dissociation times or impute/filter them. Unknown sentinels are not globally
accepted.

| Dataset | Configured sample column | Cells / samples | Read-only Bamboo audit |
| :--- | :--- | :--- | :--- |
| Alzheimer | `donor_id` | 1,395,601 cells/83 samples | `assay` mixed in 21 samples; winner fraction 66.96–100% (median 100%); `10x 3' v3` winner in all 83; `sex` 100% in every sample. |
| Breast_cancer | `sample_id` | 714,331 cells/167 samples | `assay`, `sequencing_platform`, and `suspension_dissociation_time` each 100% within every sample; 65,359 cells/11 samples have literal `unknown` dissociation time. |
| Lupus_PBMC | `sampleID` | 1,263,676 cells/261 samples | `batch_cov` winner fraction 30.95–100% (median 100%); 70 samples mixed, 23 levels, no ties. |
For Lupus, the all-sample winner-fraction quantiles (minimum, 1%, 5%, 10%,
25%, median, 75%, 90%, 95%, 99%, maximum) were 30.95%, 33.40%, 42.66%,
51.95%, 65.81%, 100%, 100%, 100%, 100%, 100%, and 100%. Among the 70 mixed
samples, 10 were in the 30--40% bin, 9 in 40--50%, 33 in 50--60%, 17 in
60--70%, none in 70--80% or 80--90%, and one in 90--<100%. The lowest
sample was `IGTB195_IGTB195` at 3,952/12,768 (30.95%), followed by
`IGTB514_IGTB514` at 4,022/12,491 (32.20%) and `IGTB469_IGTB469` at
4,374/13,543 (32.30%).


**Historical implementation boundary.** The majority policy is downstream-only. The
obs-only Python exporter reads H5AD `obs` metadata and writes the
sample-level Feather table plus its checksum; it does not alter H5ADs or
open `X`, `raw`, or `layers`. It votes only the explicitly affected keys:
`assay` for Alzheimer, `suspension_dissociation_time` for Breast_cancer, and
`batch_cov` for Lupus_PBMC. All other technical covariates remain at their
first-observation values, and labels and Sample IDs are never voted.

The exporter implementation is
`src/utils/py/export_h5ad_sample_metadata.py`, with the scoped policy in
`datasets.json`. Pipeline 3 corrected processing keeps the original
cell-level batch metadata for Harmony/HVG and performs no sample-level batch
check or majority assignment; ordinary H5AD content, checksum, and
configuration/provenance gates remain. The majority assignments are used only
by corrected Stage 5 sample-level consumers such as ECODA composition and
pseudobulk/limma. Uncorrected Stage 5 performs no batch correction and does
not use the majority assignments.

**Interpretation and limitation.** The mode is a declared pragmatic metadata
policy, not proof that minority cells are mislabeled. Exact ties and actual
missing/blank/non-finite values remain errors, and no minimum winner fraction
is imposed. Future source changes require rerunning the obs-only metadata
export/audit. Existing H5AD/RDS/Feather artifacts remain untouched.

No pipeline rerun or existing artifact invalidation occurred while implementing
this policy.

#### Lupus_PBMC: `batch_cov` and replicate-aware sample metadata

##### Summary of findings

1. The authors explicitly used `batch_cov`, representing the 23 multiplexed
   library pools, to correct single-cell batch effects with ComBat.
2. `sampleID` is the donor/patient identifier (261 unique donors). The
   experiment deliberately included replicates across pools and processing
   batches: 355 total sample runs across 23 pools. Consequently, 70 donors
   were sequenced in multiple pools, so their cells naturally split across
   multiple `batch_cov` levels.
3. `Processing_Cohort` is not a better substitute. Donors were also replicated
   across the four processing cohorts, while that variable is much coarser
   than the 23-pool `batch_cov` variable.

##### 1. What the original study did (Perez et al. 2022, *Science*)

The main text and supplementary methods (`science.abf1970_sm.v2 (1.pdf)`)
describe the following:

- **Experimental design (Fig. S1A–B and Supplementary Methods pp. 2–3):**
  - 355 total sample runs were profiled across 23 multiplexed pools over four
    processing batches.
  - The 355 runs came from 264 individuals (162 SLE cases, 49 CLUES healthy
    controls, 50 ImmVar healthy controls, and 19 flare cases), including 94
    replicates and 10 longitudinal samples.
  - Each pool contained 7–19 multiplexed donors processed in one 10x Chromium
    channel and demultiplexed with freemuxlet.
- **Single-cell batch correction (Supplementary Methods p. 3):**

  > “In total, 1,263,676 cells remained in the final dataset. The data was
  > then adjusting for pool using COMBAT (60). The most variable 1,999 were
  > retained and the count matrix was rescaled.”

  In the released AnnData/Seurat object, the covariate encoding these 23 pools
  is named `batch_cov`, with values such as
  `dmx_YS-JY-21_pool2` and `dmx_YE_7-19`.
- **Visualization (Fig. S2C–D):**
  - Fig. S2C displays “UMAP projection colored by processing pool
    (batch_cov).”
  - Fig. S2D displays “UMAP projection colored by processing batch
    (Processing_Cohort).”
- **Pseudobulk/differential expression (Supplementary Methods p. 5):**
  - For cell-type-specific pseudobulk EdgeR differential expression, the
    authors aggregated counts per individual donor and included processing
    batch (that is, `Processing_Cohort`) and age as covariates.

##### 2. Relationship between the metadata columns

In the canonical dataset:

| Column | Unique levels | Meaning |
| :--- | ---: | :--- |
| `sampleID` | 261 | Individual donors/patients after excluding ImmVar case-control samples and samples with fewer than 100 cells. |
| `batch_cov` | 23 | Multiplexed 10x library pools: the 23 physical capture pools, for example `dmx_*`. |
| `Processing_Cohort` | 4 | Processing batches: Batch 1 = 10 pools, Batch 2 = 3 pools, Batch 3 = 4 pools, Batch 4 = 6 pools. |
| `ind_cov_batch_cov` | 355 | Donor × pool runs (`<sampleID>:<batch_cov>`), the physical sequencing units. |

Because `sampleID` is defined at the donor level (261 donors), rather than at
the run level (355 runs):

- 191 donors were sequenced in exactly one pool, giving 100% constant
  `batch_cov`.
- 70 donors were intentionally sequenced as replicates across multiple pools
  and batches, giving mixed `batch_cov` values; the lowest observed majority
  was 30.95%.

##### 3. Would `Processing_Cohort` be better?

No. Switching to `Processing_Cohort` does not resolve the issue:

1. Replicates also span processing cohorts. Fig. S1B specifically notes
   “Processing Batch 1: Repeated samples” and “Processing Batch 3: Repeated
   samples”; Fig. S2E assesses correlation between cell-type percentages for
   biological replicates from different batches. In the test subset, 9 of 10
   samples are split across multiple `Processing_Cohort` levels.
2. `Processing_Cohort` is coarser: it collapses 10 pools into batch 1, 3 into
   batch 2, and so on. The primary technical batch drivers in droplet
   scRNA-seq (cell capture efficiency, ambient RNA contamination, GEM
   emulsion, and sequencing lane) occurred at the individual-pool level
   represented by `batch_cov`, which is why the authors applied ComBat to the
   pools.

##### 4. Does `batch_cov` make sense to correct for?

- **Cell-level workflows (Harmony, Scanpy, or ComBat):** Yes. Every cell
  belongs unambiguously to one library pool, so correcting for `batch_cov`
  directly mirrors the original study’s ComBat strategy.
- **Sample-level workflows (pseudobulk, CLR composition, or limma/DESeq2):**
  - With `sampleID` (261 donors) as the sample unit, assigning one batch label
    requires the documented majority vote. This is a pragmatic heuristic; 70
    replicated donors contain cells from other pools, and 19 have a majority
    below 50%.
  - `ind_cov_batch_cov` is the only column that is 100% unmixed per sample
    (355 samples). Treating it as the sample column would create replicate
    samples for the same patient and introduce pseudoreplication into patient
    classification (Case vs. Healthy).
  - Keeping `batch_cov` as the configured batch covariate is therefore
    scientifically justified and matches the original paper’s batch
    definition.


### 3.6 Historical pre-change notebook and publication boundary

The notebook and publication-path statements below are historical pre-change
final-lane notes retained for provenance; they do not establish completion of
the current corrected-final recovery.

`notebooks/batch_effect_analysis_uncorrected_batchconfounding_contingency.rmd`
remains legacy-only and is not updated. It continues to read and write its
existing legacy paths.

`notebooks/batch_effect_analysis_uncorrected.rmd` is the final-only analysis
notebook. It reads the mixed-source final manifest and writes only to the
`uncorrected_final` analysis and plot roots. `notebooks/batch_effect_analysis_legacy.rmd`
is out of scope and remains untouched.

No publication figure is renamed, removed, or rewritten as part of this
batch-analysis processing work.

The `Supp_fig_22` publication-figure name and placement remain deferred. No
publication figure is renamed, removed, or rewritten while the naming and
figure-hierarchy decision is pending.

## 4. Benchmark regeneration and annotation invariants

- Stage 3 counts observations per raw/staged view before Scanpy's per-cell and
  per-gene filters. Samples with fewer than 500 observations are removed;
  exactly 500 are retained. The authoritative raw/view audit identifies
  `BIOKEY_8_Pre` (365) and `BIOKEY_25_Pre` (296) in Bassez, `LB4180T` (496) in
  Lee, six Smillie samples (`N58.LPB2=498`, `N19.LPB=485`, `N8.LPB=482`,
  `N12.LPB=441`, `N14.LPA=432`, `N12.LPA=243`), and Zhang
  `Pre_P010_t=8`/`Pre_P018_b=437`. Current processed mirrors are evidence only.
- Missing high-resolution Bassez annotations are filled from the configured
  low-resolution/broad annotation. This is an accepted annotation-contract
  change, not a biological-label substitute. Bassez and Smillie still run the
  supported HiTME/scATOMIC Stage 4 path after Stage 3 changes so derived
  annotation artifacts used by Figure 3 and Supp fig 19 are fresh.
- scPoli's label encoder can receive mixed strings and `NaN` values and fail
  under NumPy 2.x. PILOT's cost matrix can treat `NaN` as a pseudo-cell type
  with a zero-cell centroid, yielding `NaN`/invalid EMD distances. The worker
  therefore replaces missing cell-type annotations with an explicit `Unknown`
  category, preserving every cell and sample; complete annotation columns are
  unchanged. Missing cells are never dropped and the biological label is never
  used in this handling.
- PILOT-GM uses the configured cell-type annotation only for the historical
  number-of-components choice; its distance routine consumes the model-generated
  component assignments rather than the biological annotations. The upstream
  full-covariance path can emit undefined `np.cov` values for singleton
  component assignments, so the worker applies a finite positive-semidefinite
  covariance repair before EMD and rejects any remaining nonfinite output.
- MrVI is stochastic. CPU/GPU choice and preprocessing/input changes can alter
  learned distances despite `scvi.settings.seed = 0`.
- Gene-expression methods can differ slightly because the old Seurat
  preprocessing and current Scanpy preprocessing are not bit-identical.
- Zhang old/new differences additionally include the two low-cell samples,
  explicit `Unknown` cell-type handling, stored-HVG/raw-count inputs, and
  PILOT's switch from recomputed PCA to the stored preprocessing embedding.

## 5. Modularity graph update in the new benchmark pipeline

The modularity implementation changed during the repository migration from
the March legacy pipeline to the current shared scoring code. The change is
intentional and affects the modularity score, not the underlying biological
labels or pseudobulk features.

### 5.1 Legacy graph

The March implementation in `functions.R` iterated only over each sample's
directed k-nearest-neighbor list. It calculated the number of shared neighbors
for those candidate pairs and symmetrized the result. This produced a
kNN-restricted SNN graph: pairs that shared neighbors but were not direct kNN
neighbors were omitted.

### 5.2 Current graph

`src/utils/scoring_metrics.R::compute_snn_graph()` builds a sparse binary
sample-to-neighbor incidence matrix `A` and computes:

$$
S = A A^T
$$

After removing the diagonal, `S_{ij}` is the number of nearest neighbors shared
by samples `i` and `j`, for all sample pairs. This is the standard full SNN
edge set and matches the edge construction used by Seurat's `ComputeSNN`;
the current ECODA implementation intentionally retains raw shared-neighbor
weights and does not apply Seurat's Jaccard normalization or pruning threshold.

`Matrix::tcrossprod()` is used rather than bare `t()` so sparse-matrix
operations remain valid in the transformation workers. `knn_k` is clamped to
`n_samples - 1`, which prevents invalid neighbor indices on small datasets such
as `_debug`.

### 5.3 Score naming and comparability

The old `mod_score` was the default modularity using
`k = max(3, round(sqrt(n_samples)))`. The current equivalent is
`mod_knnsqrtn_score`; fixed-neighbor scores are explicitly reported as
`mod_knn3_score`, `mod_knn6_score`, and `mod_knn9_score`. The old and current
values are not numerically interchangeable because the graph edge set changed.
The current notebook therefore drops the obsolete `mod_score` field rather than
mixing it with the new score names.

`igraph::modularity()` remains the underlying weighted modularity calculation.
ECODA additionally divides it by `1 - 1 / n_groups` as a project-specific
group-count adjustment. This adjustment is unchanged by the graph migration.
For exact March reproducibility, the legacy kNN-restricted graph must be used;
for the current method definition, report the full unpruned raw-overlap SNN
score and do not describe it as Seurat's Jaccard-pruned score.

### 5.4 Supp fig X2 modularity diagnostics

Supp fig X2A and X2B are information-only diagnostics for the four current
modularity scores: `mod_knn3_score`, `mod_knn6_score`, `mod_knn9_score`, and
`mod_knnsqrtn_score` (the square-root-of-sample-size variant). Panel X2A shows
their pairwise distributions and correlations. Panel X2B plots each variant
against the same shared `Mean_Score` used in panel 1B:

$$
\text{Mean\_Score} =
\operatorname{mean}(\text{ANOSIM},\ \text{Modularity}_{k=3},\
\text{ARI},\ \text{Silhouette},\ \text{LISI})
$$

The mean is not recomputed separately for each modularity variant; this keeps
the X2B comparisons on the same reference scale as 1B. In the current
benchmark-analysis scope (110 dataset/method rows), the pairwise Pearson
correlations are 0.93--0.98, confirming that the modularity choices are very
similar. The `mod_knn3_score` variant remains the default because it is the
established benchmark metric and is the `Modularity` entry in
`score_label_map`; panels X2A and X2B do not introduce an alternative default.

## 6. New methods

- PILOT-GM-VAE: Implemented
- MOFAcellulaR might be added later as an optional R benchmark method
- PULSAR: The PULSAR model card explicitly warns that the checkpoint may not work for tissues other than PBMC

## 7. PERMANOVA assumptions and limitations

- Centroid vs. Dispersion Confounding: PERMANOVA is sensitive to differences in multivariate dispersion (spread/variance) between groups. If Group A has higher variance than Group B, PERMANOVA can return a significant $p$-value even if the centroids are identical.
  - Unbalanced designs exacerbate bias:
    - If the group with a larger sample size has higher dispersion, PERMANOVA becomes overly conservative (loss of power).
    - If the group with a smaller sample size has higher dispersion, PERMANOVA becomes overly liberal (elevated Type I error rate).
  - We did not run a dispersion test prior to PERMANOVA (e.g., betadisper in R / PERMDISP) alongside it to test dispersion per group.
  - Possible mitigation strategy: use Wd test (Welch-type PERMANOVA):
    - Proposed by Hamidi et al. (2019) for distance-based tests. It modifies the test statistic using a Welch-Satterthwaite-like adjustment to remain robust against heteroscedasticity, even with unbalanced sample sizes.
- Exchangeability & Strata: If data have nested or repeated-measures structures (e.g., subjects sampled across multiple timepoints), permutations must be restricted within blocks/strata to preserve the underlying dependency structure.
  - We did not assess detailed nesting and repeated-measures structures for the tested batch datasets, due to time constraints and batch effect analysis/mitigation was not the main scope of this project.

### ANOSIM vs PERMANOVA

Would They Give the Same Result for a Single Covariate? - No, they will generally not yield identical results.
- Nonlinear monotonic shifts: Monotonically stretching or compressing dissimilarities will change the PERMANOVA pseudo-$F$ statistic and potentially its $p$-value, but will leave the ANOSIM $R$ and $p$-value completely unchanged.
- Sensitivity to magnitude: PERMANOVA accounts for the actual magnitude of distances, giving larger absolute differences greater weight, whereas ANOSIM flattens these differences into rank orders.
- Agreement: They usually point in the same qualitative direction (both tend to reject $H_0$ under strong separation), but their $p$-values and effect sizes will differ.

## 8. Reduced joint batch-covariate designs

The original batch-effect joint PERMANOVA attempted to include every valid
registered biological and technical covariate as a factor in one full model.
Nine of twelve datasets were correctly rejected as `NON_ESTIMABLE` because the
design matrix contained aliased terms or more factor columns than the sample
universe could support. This was not primarily a missing-value problem:
near-unique identifiers, nested metadata partitions, exact duplicate
partitions, and high-cardinality categorical fields were the dominant causes.

The uncorrected analysis now uses an explicit reduced-design contract while
retaining every original candidate in the score table:

1. The primary biological label is never dropped or imputed.
2. Non-primary candidates with more than 50% unique levels are excluded from
   the joint model as `DROPPED_NEAR_UNIQUE`. This removes subject/sample-like
   identifiers that would consume most model degrees of freedom.
3. Exact duplicate metadata partitions are represented by one composite model
   term whose name concatenates the original fields with `__`, for example
   `field_a__field_b`. The original fields remain in the table with
   `MERGED_EXACT_PARTITION` status and are not entered separately. A duplicate
   of the primary biological label maps to the primary term instead of
   replacing or concatenating the ground-truth label.
4. Nested reductions are explicit, not inferred from correlation alone:
   `Lung: origin_fine -> origin` and
   `Kidney_KPMP: condition.l2/condition.long -> condition.l1`. The less
   granular field is retained when it is present and estimable. No generic
   rule drops technical fields merely because one field appears nested in
   another; fields such as site and dissociation protocol remain candidates
   unless the observed design matrix makes them algebraically aliased.
5. Remaining algebraic aliases are resolved deterministically in registry
   order: primary biology first, registered technical candidates next, and
   secondary biology last. When an unrelated term must be removed to obtain an
   estimable design, the latest aliased retained term is marked
   `DROPPED_ALIASED` with the retained-term list in its warning. This is an
   estimability decision, not a claim that the dropped field is scientifically
   irrelevant.

Reduced models are labelled `ESTIMABLE_REDUCED`; their four decomposition
components still sum to one. Candidate-level rows retain `joint_model_term`,
`joint_model_class`, `reduction_status`, `reduction_warning`, and the global
design warning so downstream interpretation can distinguish a reduced
estimable model from the original unreduced specification.

Missing values remain missing. Literal `NA`, `nan`, `None`, and `<NA>` are
sentinels, not levels. Technical missingness could be represented as an
explicit `Unknown` level only in a separately documented sensitivity analysis;
the main batch-effect pass does not impute metadata because that would create
an artificial covariate level and could manufacture batch signal. Biological
labels remain ground-truth evaluation fields and are never imputed.

## 9. Standalone ECODA derived analyses

The first approved derived-analysis phase is read-only with respect to the
completed Stage 3--5 outputs.  It does not edit `datasets.json`, source H5ADs,
canonical composition/pseudobulk/MOFA RDS files, or generic Stage 5 dispatch
and validators.  No existing benchmark method is rerun or overwritten.

### 9.1 Independent res50 and Harmony bundles

The local derived runner has three explicit, mutually exclusive selectors:
`--analysis res50`, `--analysis harmony`, and `--analysis cell_subsetting`.
The first two selectors are independent standalone analyses; neither selector
reruns or overwrites an ordinary Stage 5 method, and neither selector emits the
other selector's output.

The `res50` run emits exactly one RDS per selected dataset,
`<DS>_ECODA_seuratres_50.rds` (+ `.rds.md5`).  Its source contract is the
existing `leiden_res_50_benchmark_analysis_hvg2000` observation column and
`X_pca_benchmark_analysis_hvg2000` embedding; its manifest method is
`ECODA_seuratres_50`, displayed as `ECODA_Leiden_res_50`.

The `harmony` run emits exactly one RDS per selected dataset,
`<DS>_ECODA_seuratres_2_harmony.rds` (+ `.rds.md5`).  Its source contract is
the existing `leiden_res_2_benchmark_analysis_hvg2000_harmony` observation
column and `X_pca_harmony_benchmark_analysis_hvg2000` embedding; its manifest
method is `ECODA_seuratres_2_harmony`, displayed as
`ECODA_Leiden_res_2_harmony`.

Each invocation is a separate run identity: use a distinct run-owned output
directory and `run_id` for res50 versus Harmony.  The manifest `analysis`
field, not the filename alone, is authoritative for selecting the notebook
root.

The runner's top-level output root is `ECODA_DERIVED_ROOT`, defaulting to
`data/benchmark/results/derived`; every invocation must use a fresh strict
descendant run directory.  Keep selector identities separate with
analysis-specific paths such as `data/benchmark/results/derived/res50/<run_id>`
and `data/benchmark/results/derived/harmony/<run_id>`.  The notebook roots
point directly to those completed run directories through
`ECODA_DERIVED_RES50_DIR` and `ECODA_DERIVED_HARMONY_DIR`.  When either
variable is unset, that derived analysis is not loaded and baseline notebook
behavior is unchanged.  When set, the root must name an existing run
directory and the loader fails closed on non-completed status, wrong analysis,
mismatched run/source provenance, escaped paths, checksums, source H5AD
identities, or invalid bundle fields.  The runner and notebook both require
complete, ordered benchmark sample bundles.
Run and source manifests record the exact selector analysis.  A completed
res50 or Harmony run manifest must record exactly one expected output method
and filename per source dataset; the cell-subsetting manifest retains its exact
RDS/CSV pair.  Raw-config union selection, source checksum/schema/root/owner
checks, and all existing fail-closed validation remain in force.

The standalone res50 method remains in the shared Figure 3A/X3B annotation
method block when its root is opted in.  It is intentionally not added to
Figure 2A or Supp fig 2, and the ordinary Stage 5 `ECODA_seuratres_2` method
is unchanged.  Neither derived path recomputes PCA, Harmony, neighbors, or
Leiden, and Harmony is never mapped to an `RNA_snn_res.*` alias.

### 9.2 Author cell-depth subsetting

The cell-subsetting run is a separate explicit run-owned directory selected by
`ECODA_DERIVED_CELL_SUBSETTING_DIR`; it is opt-in and has no notebook default.
The local-composition invocation is:
`run_local_ecoda_derived.R --analysis cell_subsetting
--scope benchmark_union --composition_dir data/benchmark/results --cores 8`.
Selection comes from the `datasets.json` benchmark union, excludes underscore-
prefixed keys in production scope, and requires exactly one local
`<DS>_ECODA_authors_HR.rds` plus `<DS>_metadata.rds` pair per selected dataset.

The local source bundles persist `counts_all, 0.5` values rather than raw
integer counts.  The runner first verifies the integer-plus-0.5 contract,
recovers counts by subtracting 0.5, and preserves the metadata
`cells_per_sample` totals as diagnostics.  When those totals exceed recovered
annotated counts (Bassez, Kim, Lee, and Zhang), the recovered annotated totals
define the subsetting population; metadata totals are never forced into the HR
composition.  Source paths, size, mtime, and computed MD5 are recorded and
rechecked; canonical input RDS files and their missing sidecars are untouched.

The exact target order is `all cells, 2000, 1000, 500, 400, 300, 200, 150,
100, 50`.  The `all cells` baseline is one unmodified row
(`replicate = 0`, `seed = NA`).  Every other target has exactly 20
without-replacement subsamples with seeds `101:120`; each sample keeps
`min(recovered_annotated_cells, target)` cells and all samples remain present.
Because cell IDs are absent, the runner uses category-level
multivariate-hypergeometric sampling.  This has the same distribution as
uniform cell-level sampling but is explicitly an aggregate approximation.

Independent datasets run in parallel PSOCK workers (`--cores 8`); each worker
reads one source pair and returns its 181 validated result rows.  Replicate-
level ANOSIM values, sample IDs, effective per-sample counts, metadata totals,
recovered annotated counts, target, replicate, seed, and the approximation
method remain in the RDS/CSV diagnostics.  The completed run is
`data/benchmark/results/derived/local_composition_subsetting_parallel_20260910_release/run/`.

The final notebook plot uses dataset means across the 20 subsamples, connects
dataset-level points with lines, and draws aggregate dataset-mean bars with
standard-error whiskers.  Error bars are drawn first so bars cover their
central segment; bar width is 0.9 and the output is 4 x 3.2 inches:
`Supp_fig_X_ECODA_authors_HR_cell_subsetting.pdf`.

### 9.3 One-time HPC composition snapshot and local reuse
The full processed `benchmark_analysis` H5ADs remain on Bamboo scratch at
`/home/users/h/halterc/scratch/ECODA_paper/<dataset>/output/<benchmark_analysis output_file_name>`
with their adjacent `.md5` files.  The composition extractor is a one-time,
explicitly scoped HPC job: it reads each full H5AD once, pulls only the
cell-level metadata needed to form the compact composition tables, and
publishes a fresh run-owned `derived_composition_snapshot.rds` with its
checksum and `composition_snapshot_manifest.json`.  The cache records the
verified remote H5AD path, MD5, size, mtime, and configured label/high-
resolution column names as provenance; those source files are not copied to
macOS.

After the compact cache is synchronized locally, each standalone selector
(`res50`, `harmony`, and `cell_subsetting`) may reuse that same snapshot with
`--composition_snapshot`.  Snapshot-backed derived run/source manifests record
`source_mode=composition_snapshot`, the local snapshot path, and its MD5.
Notebook loading verifies the local snapshot manifest, RDS sidecar checksum,
size, config identity, source provenance, and every derived-output checksum.
Remote H5AD paths are retained for audit only and are never reopened in this
mode; ordinary H5AD-mode runs continue to require local source files and
matching sidecars.

This cache workflow does not invoke any Pipeline 1--5 script and does not
change the Figure 2A/default method set.  When opted in, res50 and Harmony
remain in the shared Figure 3A/X3B scope, while cell-subsetting remains the
independent `Supp_fig_X_ECODA_authors_HR_cell_subsetting.pdf` analysis.  No
expression/count matrix, embedding, or full H5AD is transferred to or
retained on macOS.

### 9.4 Deferred MOFAcellulaR decision

MOFAcellulaR remains conditional and is not loaded, installed, rerun, or added
to any method list in this phase.  Its package source/verified Git SHA and
runtime/packaging strategy are still a user-reviewed decision because the
current Pixi lock cannot represent the verified SHA natively.  A separate
`_debug` feasibility pass (and review of its result) is required before any
production MOFAcellulaR artifact or durable HPC launch.  Until that decision
and debug review are complete, Figure 2A defaults and all existing MOFA
artifacts remain unchanged.

### 9.5 Fast obs-only figure refresh — 2026-09-10

For the immediate Figure 3A/Supp fig 18 refresh, Bamboo job `4398864` read
only `obs` metadata from the 11 configured benchmark-analysis H5ADs:
`Sample`, the configured biological-label column, and the persisted res50 and
Harmony Leiden columns. It did not read `X`, `layers["counts"]`, `obsm`, or
transfer any H5AD. The compact local artifacts are under
`data/benchmark/results/derived/fast_obs_20260910/`, including the long
composition table, labels, source metadata, and res50/Harmony sample-by-cell
type count matrices.

The refreshed outputs are
`plots/Figure_3_A_annotationmethods_barplot_anosim.pdf`,
`plots/Supp_fig_18_annotationmethods_barplot_mod_ari.pdf`, and
`plots/Figure_X3_B_number_of_celltypes.pdf`. Figure 3A and Supp fig 18 include
`ECODA_Leiden_res_50` and `ECODA_Leiden_res_2_harmony`; Figure X3B includes
both methods, starts its log10 y-axis at 2 cell types, labels the axis
`Number of cell types (log10 scale)`, and uses intermediate 2–3–5 ticks
through 1000. It remains 5×5 inches so the legend is not clipped. Figure 2A
defaults are unchanged. The compact-input refresh is a fast figure-generation
artifact, not a replacement for the run-owned snapshot workflow used by the
standalone runner.

The fast MOFA sanity check found `MOFA2` and `reticulate`, but
`MOFAcellulaR` is not installed in the existing Bamboo environment
(`MOFAcellulaR=FALSE`), so no `_debug` model could run. The wrapper remains
source-only and fail-closed. Its debug source check uses the current source's
actual unique Sample IDs (at least two), rather than assuming a fixed
five-sample universe.

### 9.6 One-off cell-type filtering diagnostic — 2026-09-10

One read-only diagnostic used the same `datasets.json` benchmark-union
selection and recovered HR counts.  It retained all-zero factor levels for
reporting and used metadata total cells as the low-abundance denominator;
annotated totals were reported alongside them.  No diagnostic artifact was
written.

The number of cell types zero in more than 50/60/70/80/90% of samples was:

| Dataset | >50% | >60% | >70% | >80% | >90% |
| --- | ---: | ---: | ---: | ---: | ---: |
| Adams | 3 | 2 | 2 | 1 | 0 |
| Bassez | 4 | 2 | 1 | 1 | 1 |
| GongSharma | 2 | 1 | 1 | 0 | 0 |
| Kfoury | 3 | 2 | 2 | 0 | 0 |
| Kim | 28 | 22 | 13 | 11 | 3 |
| Lee | 21 | 19 | 18 | 14 | 9 |
| Pelka | 22 | 9 | 5 | 3 | 1 |
| Smillie | 4 | 1 | 1 | 1 | 1 |
| Stephenson | 3 | 1 | 1 | 0 | 0 |
| Wu | 5 | 3 | 2 | 1 | 1 |
| Zhang | 17 | 16 | 9 | 4 | 1 |

The corresponding all-cell ANOSIM results are shown below as
`score (change from the all-cell baseline)` after removing cell types zero in
more than the indicated fraction of samples:

| Dataset | Baseline | >50% | >60% | >70% | >80% | >90% |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Adams | 0.4605 | 0.4311 (-0.0294) | 0.4472 (-0.0133) | 0.4472 (-0.0133) | 0.4592 (-0.0013) | 0.4605 (0) |
| Bassez | 0.4046 | 0.3643 (-0.0403) | 0.3874 (-0.0172) | 0.3880 (-0.0166) | 0.3880 (-0.0166) | 0.3880 (-0.0166) |
| GongSharma | 0.6953 | 0.6952 (-0.0000) | 0.6950 (-0.0003) | 0.6950 (-0.0003) | 0.6953 (0) | 0.6953 (0) |
| Kfoury | 0.4252 | 0.4368 (+0.0116) | 0.4376 (+0.0124) | 0.4376 (+0.0124) | 0.4252 (0) | 0.4252 (0) |
| Kim | 0.8757 | 0.8421 (-0.0335) | 0.8419 (-0.0338) | 0.8685 (-0.0072) | 0.8693 (-0.0064) | 0.8804 (+0.0047) |
| Lee | 0.2423 | 0.2839 (+0.0416) | 0.2821 (+0.0399) | 0.2655 (+0.0232) | 0.2468 (+0.0046) | 0.2327 (-0.0096) |
| Pelka | 0.8448 | 0.8017 (-0.0432) | 0.8083 (-0.0366) | 0.8524 (+0.0076) | 0.8502 (+0.0053) | 0.8448 (-0.0000) |
| Smillie | 0.1227 | 0.1280 (+0.0053) | 0.1251 (+0.0024) | 0.1251 (+0.0024) | 0.1251 (+0.0024) | 0.1251 (+0.0024) |
| Stephenson | 0.3953 | 0.3795 (-0.0158) | 0.3806 (-0.0147) | 0.3806 (-0.0147) | 0.3953 (0) | 0.3953 (0) |
| Wu | 0.1799 | 0.2361 (+0.0562) | 0.1991 (+0.0192) | 0.1935 (+0.0136) | 0.1855 (+0.0057) | 0.1855 (+0.0057) |
| Zhang | 0.7556 | 0.6815 (-0.0741) | 0.6741 (-0.0815) | 0.7887 (+0.0331) | 0.7666 (+0.0110) | 0.7566 (+0.0010) |

The effect is dataset-specific rather than uniformly inflationary: removal
increases ANOSIM for Kfoury, Lee, and Wu at the 50% threshold, while it
decreases ANOSIM for Adams, Bassez, GongSharma, Kim, Pelka, and Stephenson.
Smillie changes only slightly.  Zhang decreases at 50--60% but increases
after the stricter filters.  Bassez and Kim baseline values include their
pre-existing all-zero `NA` factor level, which is retained here to reproduce
the persisted all-cell composition.

Bassez and Kim each contain an all-zero `NA` factor level; it is included in
the table but is not a biological annotation.  For Wu, the five >50% zero
types are `B cells Naive`, `Mature Luminal`,
`Myeloid_c5_Macrophage_3_SIGLEC1`, `Myoepithelial`, and
`T_cells_c5_CD8+_GZMK`.  Removing those from the all-cell composition changes
ANOSIM from 0.1799 to 0.2361.  At >60%, >70%, >80%, and >90%, the corresponding
scores are 0.1991, 0.1935, 0.1855, and 0.1855.

At the 50-cell target, Wu has at least one all-zero category in 19/20
replicates, with 2.05 all-zero categories on average (maximum 4).  The
repeatedly affected types include `Myeloid_c5_Macrophage_3_SIGLEC1` (zero in
17/20 replicates), `Cycling PVL` (8/20), `Myeloid_c7_Monocyte_3_FCGR3A`
(6/20), and `Myeloid_c0_DC_LAMP3` (5/20).  The 50-cell Wu ANOSIM mean is
0.2486 (range 0.1678–0.3201), so the upward shift is sampling variability
consistent with rare-category loss rather than a monotone depth effect.

Low-abundance all-cell filtering also changes scores in both directions.  For
Wu, removing types below 1% of metadata-total cells removes 17 types and
changes ANOSIM to 0.2072 (+0.0273); below 0.1% removes 3 types and changes it
to 0.1887 (+0.0088).  The diagnostic was informational only and is not used
by the production subsetting run.