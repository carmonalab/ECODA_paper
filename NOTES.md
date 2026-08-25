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

### 3.1 Fixed method contract

The batch suite contains exactly these configured high-resolution methods:

- ECODA configured author high tier (`ECODA_authors_HR`);
- applicable HiTME ECODA (`ECODA_HiTME_HR_layer2`);
- applicable scATOMIC ECODA (`ECODA_scATOMIC_HR`);
- ECODA Leiden resolution 2 (`ECODA_seuratres_2`);
- deterministic shuffled-label ECODA baseline (`ECODA_authors_HR_NULL`);
- Pseudobulk;
- GloScope;
- PILOT;
- PILOT-GM-VAE;
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
