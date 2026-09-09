# ECODA Paper — Architectural Decisions & Methodology Notes

> **Document Purpose:**  
> This document serves as the central engineering and methodology log for the *ECODA* manuscript revisions and benchmarking extensions. It records the statistical rationale, architectural designs, diagnostic tool implementations, and open discussion points for the project.

---

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

### 1.4 New onboarding cohorts: sample-count comparison with PILOT-GM-VAE

The current counts below come from the full-file onboarding audits in
`data/new_dataset_checks/subsets/*_meta.json`, using the sample column registered
in [`datasets.json`](datasets.json). The cell threshold is strict: samples with
fewer than 500 cells are dropped, while samples with exactly 500 cells are
retained.

| Dataset | Current sample column | Current samples | Dropped (<500) | Retained | PILOT-GM-VAE reported samples |
| :--- | :--- | ---: | ---: | ---: | ---: |
| Alzheimer | `donor_id` | 83 | 0 | 83 | 83 |
| Breast cancer | `sample_id` | 167$ | 2 | 165 | 126$ |
| Covid-19 PBMC | `sampleID` | 172 | 8 | 164 | 151 |
| Kidney (KPMP) | `specimen` | 47 | 2 | 45 | 45 |
| Myocardial infarction (MI-2) | `orig_ident` | 24 | 0 | 24 | 23 |
| Diabetes | `donor_id` | 56 | 0 | 56 | 52* |
| Lupus PBMC | `sampleID` | 261 | 1 | 260 | 261 |
| Lung atlas | `sample` | 304† | 19 | 285 | 165‡ |
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

### 1.5 HPC-backed donor-to-sample audit for the batch-effect onboarding view

To reconcile the configured registry units with the units in the batch-effect onboarding view, a read-only audit targeted the canonical `batch_effect_uncorrected` output on Bamboo. It used the configured `py-cuda13` runtime for HDF5/AnnData-backed metadata reads only; the larger cohort reads were dispatched through `srun` on Bamboo. Only `obs` metadata and categorical codes were read, and `.X`/count matrices were never materialized. Canonical HPC sample/donor counts below are post-500-cell-filter output units, whereas raw/local audit counts are pre-filter source metadata; the PILOT values remain publication-reported units.

| Cohort | Configured / sample-like column | Donor/patient column | Canonical HPC sample units | Donor/patient units | Exact donor→sample frequency | Observed explanation |
| :--- | :--- | :--- | ---: | ---: | :--- | :--- |
| [Alzheimer](https://pubmed.ncbi.nlm.nih.gov/39402379/) | `donor_id` (configured); `Specimen ID` (audit) | `donor_id` | 202 Specimen IDs; 1,395,601 cells | 83 | 53×2, 24×3, 6×4 | Every specimen maps to one donor; the configured registry sample remains `donor_id`. |
| [Breast cancer](https://pubmed.ncbi.nlm.nih.gov/37380767/) | `sample_id` (configured); `Sample` (canonical) | `donor_id` | 165 Samples; 713,691 cells | 124 | 90×1, 29×2, 3×3, 2×4 | 34 repeated donors contribute 75 samples; `Sample` maps to at most one donor. |
| [Covid-19 PBMC](https://pubmed.ncbi.nlm.nih.gov/33657410/) | `sampleID` (configured); `Sample` (canonical) | `PatientID` | 164 Samples; 991,227 cells | 149 | 138×1, 8×2, 2×3, 1×4 | 11 repeated patients contribute 26 retained samples; `Sample` maps to at most one patient. |
| [Diabetes](https://pubmed.ncbi.nlm.nih.gov/37697055/) | `donor_id` (configured); `dataset__design__sample` | `donor_id` | 56 Samples; 301,796 cells | 56 | 56×1 | One-to-one donor/sample view; no repeated donor/sample IDs. |
| [Kidney (KPMP)](https://pubmed.ncbi.nlm.nih.gov/37468583/) | `specimen` (configured/canonical) | `donor_id` | 45 specimens; 103,642 cells | 43 | 42×1, PRE019×3 | `specimen` maps to at most one donor; PRE019 is the sole repeated donor. |
| [Lupus PBMC](https://pubmed.ncbi.nlm.nih.gov/35389781/) | `sampleID` (configured); `Sample` (canonical) | unavailable | 260 Samples; 1,263,220 cells | unavailable | unavailable | No donor/patient/subject column exists, so repeated-donor status cannot be inferred. |
| [Lung atlas](https://pubmed.ncbi.nlm.nih.gov/37291214/) | `sample` (configured); `Sample` (canonical; registry `platform=10x AND tissue=lung`) | `donor_id` | 285 Samples; 936,636 cells | 156 | 82×1, 53×2, 5×3, 14×4, 1×10, 1×16 | 74 repeated donors contribute 203 samples; `Sample` maps to at most one donor. |
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

- **Alzheimer.** The audit unit is `Specimen ID` (202) against `donor_id` (83):
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

### Executive Summary

Batch-effect analysis is a validated two-pass workflow over the explicit
registry views `batch_effect_uncorrected` and `batch_effect_corrected`. The
uncorrected pass is the evidence gate: it always preprocesses with
`batch_key=Sample`, runs the method suite without technical correction, and
selects no batch column automatically. Only after reviewing its evidence may
one confirmed technical column per cohort be written to `datasets.json`.
Corrected execution is fail-closed while that column is `null`.

The biological label is evaluation-only. It never enters filtering, HVG
selection, normalization, PCA, Harmony, CLR correction, pseudobulk design,
MrVI covariates, or any other model input.

The exact batch roles for the onboarding pass are `Joanito` low/high
`cell.type`/`cell.type_new` and `CombinedPBMC` low/high `layer1`/`layer2`.
Batch composition consumes only the configured high-resolution role; low
resolution remains a registry consistency field. `not_suitable_for_auto_annotation`
is an a-priori skip for `Alzheimer`, `Diabetes`, and `Parkinson`, not a failed
annotation result.

The annotation worker intentionally leaves scATOMIC `breast_mode` at its
default `FALSE` for cross-cohort comparability. No caller may pass that option.

Stephenson's batch-effect view uses the full declared subset and candidate
`Site`. CombinedPBMC's raw input is `combined_pbmc.h5ad`, with explicit
uncorrected/corrected outputs `combined_pbmc_batch_effect_uncorrected_ECODAprocessed.h5ad`
and `combined_pbmc_batch_effect_corrected_ECODAprocessed.h5ad`. The old raw
basename is accepted only for guarded one-time migration.

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

`batch_effect_corrected` requires a confirmed non-null `columns.batch`. It
performs one hvg2000 pass with HVGs selected by that technical column, computes
raw PCA, then Harmony and neighbors/Leiden on Harmony. The biological label is
never protected in correction.

Pass-qualified keys are literal and never fall back:

- `X_pca_batch_effect_uncorrected_hvg2000`;
- `leiden_res_<r>_batch_effect_uncorrected_hvg2000`;
- `X_pca_batch_effect_corrected_hvg2000`;
- `X_pca_harmony_batch_effect_corrected_hvg2000`;
- `leiden_res_<r>_batch_effect_corrected_hvg2000_harmony`.

All cohorts retain Leiden resolutions `(0.1, 0.4, 2, 5, 20, 50)`. The fixed
suite consumes resolution 2, except Parkinson's configured high tier, which
uses res-5 from the corresponding pass.

### 3.3 Modality-specific corrected inputs

- **ECODA composition:** each CLR cell-type feature is fit with
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

### 3.4 Artifact and evidence contract

Every pass is isolated under
`${HPC_SCRATCH_DIR}/batch_effect/<pass>/` and
`${NAS_TARGET_DIR}/batch_effect/<pass>/`. Method bundles, distances,
pseudobulks, execution logs, manifests, watchdog status, and checksums are
pass-scoped. Active filenames and runtime identifiers use
`batch_effect_uncorrected` or `batch_effect_corrected`; no pass artifact uses a
benchmark-named identifier.

The uncorrected evidence report records completeness, levels and samples per
candidate, NMI with biology, marginal and joint PERMANOVA $R^2$ and
Holm-adjusted p-values, and constant/sample-unique/perfect-confounding
warnings. It uses 999 deterministic permutations and strict sample-order
checks. It emits one CSV per cohort plus `batch_candidate_review.csv`.

The evidence checkpoint keeps all nine new `columns.batch` values `null`.
After explicit user confirmation, each corrected run verifies paired
cell/sample identities, pass-specific checksums, NAS synchronization, CLR
zero-sum recentering, batch-only pseudobulk mode, exact Harmony keys, and
native MrVI batch arguments.

### 3.5 No-leakage invariant

Technical correction is never allowed to protect or model a biological label.
Confounded technical variables remain documented warnings, not reasons to
silently change the confirmed sample or label roles. Missing IDs, collisions,
missing labels, empty derived annotations, invalid hierarchies, and missing
exact pass keys remain hard failures.


### 3.6 Historical notebook and publication boundary

`notebooks/batch_effect_analysis.rmd` remains intentionally unchanged in this
remediation. Its historical analysis surface is not part of the uncorrected
pipeline cutover; a future extension must consume the explicit
`batch_effect_uncorrected`/`batch_effect_corrected` contracts without changing
the biological-label or `columns.batch` invariants. The future handoff is to
add any new pass-specific analysis in a separate notebook or an explicitly
reviewed revision, after the scientific decision checkpoint.

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