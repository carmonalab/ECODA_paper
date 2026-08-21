# ECODA Paper — Architectural Decisions & Methodology Notes

> **Document Purpose:**  
> This document serves as the central engineering and methodology log for the *ECODA* manuscript revisions and benchmarking extensions. It records the statistical rationale, architectural designs, diagnostic tool implementations, and open discussion points for the project.

---

## 1. Dataset Onboarding & Quality Control Pipeline

### Executive Summary
- Standardized onboarding check notebooks (`dataset_check_<Name>.qmd`) were developed across all target cohorts (9 new cohorts + `_debug` subset).
- Replaced manual inspection with an automated 6-step diagnostic suite operating strictly on unintegrated raw counts.
- Implemented **Normalized Mutual Information (NMI)** to detect metadata collinearity/confounding upfront.
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
  - $\text{NMI} > 0.70$: Severe collinearity (e.g. `Cognitive status` vs. `APOE4 status` in Alzheimer, or `Status` vs. `tissue_source` in Myocardial Infarction); batch correction against such covariates risks removing biological signal.

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
- Single-cell technical noise operates through distinct physical mechanisms across modalities: **gene expression artifacts** (sequencing chemistry, depth, ambient RNA) vs. **compositional artifacts** (dissociation protocols, cryopreservation, cell fragile destruction).
- Benchmark evaluation must distinguish between expression-level batch keys and composition-level batch keys.
- Detailed cohort profiles and candidate variables are documented in [`notebooks/dataset_onboarding/README.md`](notebooks/dataset_onboarding/README.md).

### 2.1 Modality-Specific Noise Mechanisms

| Modality | Primary Physical Batch Mechanisms | Example Candidate Metadata Columns |
| :--- | :--- | :--- |
| **Gene Expression** | Sequencing platform, library prep chemistry (10x 3' v2 vs v3, 5' vs 3'), flowcell/lane, sequencing depth, ambient RNA contamination | `assay`, `seqtec`, `sequencing_run`, `10x_chemistry`, `library_id` |
| **Cell Composition** | Enzymatic dissociation protocol, tissue digestion time, cold ischemia time, fresh vs. frozen / cryopreservation, FACS gating / cell enrichment | `tissue_type`, `PMI` (post-mortem interval), `dissociation_protocol`, `sample_preservation`, `enrichment` |
| **Cohort / Demographic** | Clinical collection site, hospital center, donor sex, age, ethnic background | `Site`, `Center`, `sex`, `Age`, `self_reported_ethnicity` |

### 2.2 Target Cohort Batch Structures Summary

- **Alzheimer (SEA-AD):** Primary Bio = `Cognitive status`; Secondary Bio = `ADNC`, `Braak stage`, `CERAD score`, `APOE4 status`; Tech = `assay` (Expression), `tissue_type`, `PMI` (Composition).
- **Breast Cancer:** Primary Bio = `clinical_subtype` (ER+, HER2+, TNBC); Secondary Bio = `tumor_grade`, `stage`; Tech = `library_id` (Expression), `patient_treatment_status` (Composition/Clinical).
- **Covid-19 PBMC:** Primary Bio = `Disease_Identity` (Control, Mild, Severe, Critical); Tech = `Processing_Site` / `Center` (Cohort), `Batch` (Expression).
- **Diabetes (Mouse Islets):** Primary Bio = `disease` / `genotype`; Secondary Bio = `age`, `diet`; Tech = `assay`, `batch` (Expression/Sample Prep).
- **Kidney KPMP:** Primary Bio = `Disease_Identity` (Reference, CKD, AKI); Tech = `Institution` / `Site` (Cohort), `Assay` (Expression).
- **Human Lung Cell Atlas (HLCA / Sikkema):** Primary Bio = `disease` (COPD, IPF, COVID-19); Tech = `sequencing_platform`, `assay` (Expression), `tissue_dissociation_protocol`, `anatomical_region` (Composition).
- **Lupus PBMC:** Primary Bio = `disease_status` (SLE vs. Healthy Control); Tech = `batch_cov` (Expression / Processing Pool).
- **Myocardial Infarction:** Primary Bio = `cell_type_annotation_level` / `disease_zone` (Myocardial Infarct, Border Zone, Remote); Tech = `sample_source`, `assay` (Expression / Prep).
- **Parkinson's Disease:** Primary Bio = `disease` (PD vs Control); Secondary Bio = `Braak_stage`; Tech = `assay` (Expression), `PMI` (Composition).

---

## 3. Batch Effect Correction & Benchmark Strategy (`batch_effect_analysis.rmd`)

### Executive Summary
- The multi-batch benchmark evaluates patient stratification methods under two conditions: **Uncorrected (Raw)** and **Corrected**.
- Batch correction is applied modality-appropriately while strictly enforcing the **No-Leakage Principle** (`design ~ 1`, no biological label guidance).
- Method outputs will be evaluated using PERMANOVA on resulting sample distance matrices.

### 3.1 Modality-Appropriate Correction Implementations

1. **Cell-Type Composition Methods (ECODA, Pseudobulk):**
   - Corrected via `limma::removeBatchEffect()` directly on the sample $\times$ cell-type CLR matrix (or DESeq2 normalized pseudobulk counts).
   - **Guardrail:** Strictly unsupervised (`design = matrix(1, nrow=N, ncol=1)`), preventing biological group label leakage.
   - **Simplex Recentering:** Post-correction row recentering ($\text{clr}^* - \text{mean}(\text{clr}^*)$) restores exact zero-sum simplex constraints.
2. **Cell-Level Expression Embedding Methods (PILOT, PILOT-GM, GloScope):**
   - Corrected using `Harmony` (`X_pca_harmony`) computed on cell-level PCA prior to Earth Mover's Distance (EMD) or GMM fitting.
3. **Deep Generative / Integrated Models (MrVI):**
   - Corrected using native `batch_key` parameter support in the model architecture.

### 3.2 Evaluation via Distance Matrix PERMANOVA
- **Why PERMANOVA is the Gold Standard for Distance Matrices:**
  - Operates **directly on pairwise distance matrices $D$** without intermediate projection or coordinate distortion.
  - Accommodates non-Euclidean metrics (e.g. PILOT Earth Mover's Distance, pseudobulk correlation distances).
  - Uses non-parametric permutation testing without multivariate normality assumptions.
  - Marginal $R^2$ (`by = "margin"`) accurately partitions unique variance between partially collinear covariates.
- **Expected Benchmark Behavior:**
  - **Ideal Batch Correction:** $\Delta R^2_{\text{Batch}} = R^2_{\text{Batch}}(\text{Raw}) - R^2_{\text{Batch}}(\text{Corr}) \gg 0$ (batch variance suppressed to near zero), while $R^2_{\text{Bio}}(\text{Corr}) \ge R^2_{\text{Bio}}(\text{Raw})$ (biological separation preserved or enhanced).

### 3.3 Open Decisions & Discussion Points (Tracked in `TODO.md`)

1. **Visualization Strategy for Multi-Batch Benchmark:**
   - *Option A: Grouped / Stacked PERMANOVA Bar Plot:* A grouped bar chart displaying $R^2_{\text{Bio}}$, $R^2_{\text{Batch}}$, and $R^2_{\text{Residual}}$ for each method under Raw vs. Corrected conditions. Clean, statistically exact, and intuitive.
   - *Option B: FunkyHeatmap:* Compact tabular heatmap summarizing multiple metrics (PERMANOVA $R^2$, Silhouette, LISI, ANOSIM, ARI) across methods. Matches the style of Figure 2A / Supp Fig 15, but may become dense with multiple batch covariates.
2. **Biological Signal vs. Batch Effect Metric ($\text{Bio} / \text{Batch}$ Ratio):**
   - *Formulation 1 (PERMANOVA Ratio):* $\text{Ratio}_{\text{PERMANOVA}} = \frac{R^2_{\text{Bio}}}{R^2_{\text{Batch}} + \epsilon}$
   - *Formulation 2 (Aggregate Metric Ratio):* Ratio of biological cluster separation (e.g. $\text{Silhouette}_{\text{Bio}}$, $\text{ANOSIM}_{\text{Bio}}$) to batch mixing ($\text{LISI}_{\text{Batch}}$).
   - Needs empirical comparison across benchmark cohorts to decide the most robust and interpretable index.
