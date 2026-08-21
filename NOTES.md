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
- The multi-batch benchmark evaluates patient stratification methods across a **Two-Pass Workflow**: (1) **Uncorrected (Raw)**, followed by (2) **Corrected**.
- **Uncorrected Pass as the Decision Gate:** Running all methods without batch correction first and inspecting distance-level PERMANOVA and NMI collinearity across modalities provides the complete picture needed to decide which specific batch variables to target per method/modality.
- **Cross-Modality Confounding Attribution:** When biological conditions and technical batches are partially confounded (e.g. `Site` with unbalanced disease/control allocations), single-modality evaluations cannot separate true biology from batch noise. Cross-modality comparison (expression vs. composition) clarifies whether separation is modality-specific or driven by shared technical confounders.
- **Evaluation & Ratio-of-Ratios Metric:** Evaluates the shift in $R^2_{\text{Bio}}$ and $R^2_{\text{Batch}}$, the $\text{Bio} / \text{Batch}$ signal-to-noise ratio, and the **Ratio-of-Ratios** ($\text{Ratio}_{\text{Corr}} / \text{Ratio}_{\text{Raw}}$).

### 3.1 Two-Pass Benchmark Workflow

```
                        Processed Cohort Data
                                  |
              +-------------------+-------------------+
              |                                       |
    [Pass 1: Uncorrected (Raw)]             [Pass 2: Corrected]
    - All methods run raw                   - Modality-appropriate correction
    - Metadata Collinearity (NMI)           - limma for ECODA/Pseudobulk
    - Distance Matrix PERMANOVA             - Harmony on PCA for PILOT/GloScope
              |                             - Native batch_key for MrVI
              v                                       |
    [Attribution Decision Gate]                       v
    - Identify active batch drivers         [Distance Matrix PERMANOVA]
    - Cross-modality confounding check      - Measure batch suppression & bio retention
              |                                       |
              +-------------------+-------------------+
                                  |
                        [Comparative Synthesis]
                        - Bio / Batch Ratios (Raw vs. Corr)
                        - Delta R² (Bio vs. Batch)
                        - Ratio-of-Ratios (Correction Benefit Index)
```

1. **Pass 1: Uncorrected (Raw) Benchmark Evaluation:**
   - Execute all benchmark methods (ECODA family, Pseudobulk, Avg_PCA, GloScope, PILOT, PILOT-GM, MrVI, QOT) on raw, uncorrected embeddings/counts.
   - For each method output distance matrix ($D_{\text{raw}}$), evaluate:
     - **Metadata Collinearity (NMI):** Maps confounding structure across metadata columns.
     - **Marginal PERMANOVA:** Quantifies baseline $R^2_{\text{Bio}}$, $R^2_{\text{Batch}}$, $R^2_{\text{Shared}}$, and $R^2_{\text{Residual}}$.
   - **Cross-Modality Confounding Rationale:** In datasets where `Site` or `Institution` has unbalanced disease distributions, it is mathematically impossible within one modality to prove whether patient clustering reflects disease biology or site-specific sample prep. Comparing composition-based methods (ECODA) with expression-based methods (Pseudobulk, PILOT, MrVI) indicates whether one modality recovers the biological signal more cleanly than the other despite the technical confounder.
   - **Decision Gate:** Finalize which batch variables to correct for each method and modality based on these empirical uncorrected PERMANOVA results.

2. **Pass 2: Modality-Appropriate Batch Correction:**
   - **Cell-Type Composition Methods (ECODA, Pseudobulk):** `limma::removeBatchEffect(x, batch = batch_col)` applied to the sample $\times$ cell-type CLR matrix or DESeq2-normalized pseudobulk expression matrix.
     - **No-Leakage Invariant:** The correction model removes only the technical `batch_col`. The biological label (ground truth) is strictly excluded (`design = NULL` or intercept-only `~ 1`).
     - **Simplex Recentering (for CLR):** Post-correction row recentering ($\text{clr}^* - \text{mean}(\text{clr}^*)$) restores exact zero-sum simplex constraints.
   - **Cell-Level Expression Embeddings (PILOT, PILOT-GM, GloScope):** `Harmony` (`X_pca_harmony`) computed on cell-level PCA prior to sample-level distribution distance calculation (e.g. Earth Mover's Distance or GMM fitting).
   - **Deep Generative Models (MrVI):** Native multi-covariate `batch_key` integration within the autoencoder architecture.

### 3.2 Quantitative Evaluation & Comparative Metrics

1. **Marginal Distance Variance Explained (PERMANOVA):**
   $$\text{Total Distance Variance} = R^2_{\text{Bio}} + R^2_{\text{Batch}} + R^2_{\text{Shared}} + R^2_{\text{Residual}}$$
2. **Biological Signal vs. Batch Effect Ratio ($\text{Bio} / \text{Batch}$ Ratio):**
   $$\text{Ratio}_{\text{Raw}} = \frac{R^2_{\text{Bio}}(\text{Raw})}{R^2_{\text{Batch}}(\text{Raw}) + \epsilon}, \quad \text{Ratio}_{\text{Corr}} = \frac{R^2_{\text{Bio}}(\text{Corr})}{R^2_{\text{Batch}}(\text{Corr}) + \epsilon}$$
   *(Also evaluated using summarized benchmark metrics: $\text{Silhouette}_{\text{Bio}} / \text{LISI}_{\text{Batch}}$).*
3. **Signal Shifts ($\Delta R^2$):**
   - Batch Suppression: $\Delta R^2_{\text{Batch}} = R^2_{\text{Batch}}(\text{Raw}) - R^2_{\text{Batch}}(\text{Corr}) \quad (\text{Target: } \gg 0)$
   - Biological Retention: $\Delta R^2_{\text{Bio}} = R^2_{\text{Bio}}(\text{Corr}) - R^2_{\text{Bio}}(\text{Raw}) \quad (\text{Target: } \ge 0)$
4. **Ratio-of-Ratios (Correction Benefit Index):**
   $$\text{RoR} = \frac{\text{Ratio}_{\text{Corr}}}{\text{Ratio}_{\text{Raw}}}$$
   An $\text{RoR} > 1.0$ quantitatively proves that batch correction improved the biological signal-to-noise ratio rather than removing biological variation.

### 3.3 Open Decisions & Visualization Options (Tracked in `TODO.md`)

1. **Visualization Format for Multi-Batch Benchmark:**
   - *Option A: Grouped / Stacked PERMANOVA Bar Plot:* Method-by-method grouped bar chart showing $R^2_{\text{Bio}}$, $R^2_{\text{Batch}}$, and $R^2_{\text{Residual}}$ under Raw vs. Corrected conditions. Highly interpretable for variance shifts.
   - *Option B: FunkyHeatmap Overview:* Compact tabular heatmap summarizing PERMANOVA $R^2$, Silhouette, LISI, ANOSIM, and ARI across all methods. Provides visual consistency with Main Figure 2A / Supp Fig 15.
2. **Metric Ratio Selection:** Decide between PERMANOVA $R^2$ ratio vs. aggregate benchmark metric ratio ($\text{Silhouette} / \text{LISI}$) for summary reporting across cohorts.
