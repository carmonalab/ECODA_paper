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
| Breast cancer | `sample_id` | 167 | 2 | 165 | 126 |
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

* `*` The PILOT-GM-VAE Table 1 lists 52 Diabetes samples. The author-provided
  PDF additionally says that four embryo samples were excluded as outliers, but
  the paper does not reconcile that statement with the 52-sample Table 1
  entry. If the exclusion is applied to those 52 samples, the effective count
  would be 48. The current onboarding file contains 56 `donor_id` values and
  has no `subset_vars` exclusion for those samples.
* For Breast cancer, the onboarding audit finds 126 `donor_id` values and 167
  `sample_id` values. `sample_id` nests `donor_id` (up to four sample IDs per
  donor), while `accSample` is equivalent to `donor_id`. The PILOT-GM-VAE
  report says 126 donors, so the most likely explanation is that the authors
  used `donor_id` (or the equivalent `accSample`) rather than `sample_id`.
  This explains the 126-versus-167 discrepancy before the 500-cell filter.
* `†` The lung dataset has 304 `sample` IDs,
  and 165 `donor_id` values. It has 304 nonzero donor/sample pairs; every
  sample maps to exactly one donor, 79 donors have multiple sample IDs, and
  the maximum is 16 sample IDs for one donor. The per-donor distribution is:
  86 donors with one sample, 56 with two, 4 with three, 17 with four, 1 with
  ten, and 1 with sixteen. The 19 sample IDs below 500 cells leave 285
  sample-level units.
* `‡` The PILOT-GM-VAE report gives 941,504 cells and 165 donors for Lung.
  Those values now match cell and `donor_id` totals
  exactly. The remaining 304-versus-165 discrepancy is therefore a unit
  choice: ECODA's configured column counts `sample` IDs, while PILOT reports
  donors.

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
