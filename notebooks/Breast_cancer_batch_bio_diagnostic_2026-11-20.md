# `Breast_cancer` Batch-vs-Biology Cell-Type Diagnostic (2026-11-20)

- **Cohort & Source:** `Breast_cancer` (Kumar et al. 2023 / HBCA + breast tumor; PILOT-GM-VAE batch dataset)
- **Active Subset (`datasets.json`):** `assay == "10x 3' v3"` (constant across all 33 samples) and `suspension_dissociation_time` in `["12 hour", "15 minute"]` -> **33 samples, 114,579 cells** (`3,472.09 ± 1,442.53` cells/sample, min `657`, max `6,190`), **58 `author_cell_type` (`Authors_HR`)** subsets and **10 `broad_cell_type` (`Authors_LR`)** subsets.
- **Data Provenance:** Derived from `data/derived/ECODA_becm_comp_20261116_v1/outputs/Breast_cancer.rds`, `data/derived/publication_figures_20261110/p6_batch_majority_current_20261110_v3/batch_effect_corrected/correction_arrow_scores.csv`, and Bamboo Stage 5 composition bundles (`Breast_cancer_batch_effect_{uncorrected,corrected}_final_composition.rds`). Full 58-cell-type table: `data/derived/phase1_metadata_audit_20261120/breast_cancer_bio_vs_batch_celltypes.csv`.

---

## 1. Cohort-Level Contingency & Multi-Method Behavior

| Dissociation Batch | `breast cancer` | `normal` | Total Samples | Cancer Share |
|---|---:|---:|---:|---:|
| **`12 hour`** | 12 | 5 | 17 | **70.6%** |
| **`15 minute`** | 7 | 9 | 16 | **43.8%** |
| **Total** | **19** | **14** | **33** | 57.6% |

- **Residual bio–batch imbalance:** Even after excluding `unknown` and non-`10x 3' v3` samples, `12 hour` samples remain ~27 percentage points enriched for `breast cancer` (`OR = 3.09`, Fisher's exact `p = 0.166`), and `N = 33` is the smallest sample count in the batch benchmark.
- **Uncorrected batch effect is >2× stronger than biological signal in ECODA:** `ECODA_authors` uncorrected gives **`bio ANOSIM = 0.183`** vs **`batch ANOSIM = 0.379`** (and Leiden `ECODA_seuratres_2` uncorrected gives `bio = 0.251` vs `batch = 0.863`).
- **Simultaneous collapse of bio & batch separation upon linear/global correction:**

| Method | `bio ANOSIM` (Unc → Corr) | `batch ANOSIM` (Unc → Corr) | `bio Silhouette` (Unc → Corr) | Behavior |
|---|---:|---:|---:|---|
| `ECODA_authors` (`limma`) | `0.183 → 0.094` | `0.379 → -0.055` | `0.068 → 0.036` | Batch removed; bio halved |
| `ECODA_Harmony_Leiden_CLR` | `0.251 → 0.235` | `0.863 → 0.805` | `0.111 → 0.101` | Harmony does **not** remove `12h/15m` compositional batch |
| `ECODA_Harmony_Leiden_CLR_limma` | `0.251 → 0.024` | `0.863 → -0.008` | `0.111 → 0.035` | `limma` collapses bio completely |
| `Pseudobulk` (`limma`) | `0.150 → 0.073` | `0.771 → -0.026` | `0.070 → 0.025` | Batch removed; bio halved |
| `MrVI` | `0.172 → 0.190` | `0.839 → 0.035` | `0.091 → 0.060` | Bio ANOSIM slightly up, Silhouette drops |
| `GloScope` | `0.180 → 0.169` | `0.820 → 0.646` | `0.076 → 0.069` | Batch largely persists (`0.646`) |
| `QOT` | `0.165 → 0.163` | `0.799 → 0.635` | `0.091 → 0.084` | Batch largely persists (`0.635`) |
| `PILOT` | `0.106 → 0.104` | `0.428 → 0.409` | `0.065 → 0.064` | Low bio & persistent batch |

---

## 2. Which Cell Types Drive Batch vs. Biology (`58` `Authors_HR` Subsets)

### A. The 11 Dissociation-Associated Cell Types (`p_bh_batch < 0.05`, `12h` vs `15m`)

| Cell Type | Total Cells | Zero Frac | `Δ_batch` (`12h-15m`) | `p_bh_batch` | `Δ_bio` (`ca-no`) | `p_bio` | `Δ_bio` (in `12h`) | `Δ_bio` (in `15m`) | `Δ_batch` (in cancer) | `Δ_batch` (in normal) | Dropped at `t=2`? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **`pericytes`** | 7,176 | 0.030 | **+3.486** | 0.00037 | +1.792 | 0.0607 | **-0.498** | +1.173 | +2.367 | **+4.037** | **Yes** |
| **`Fibro-SFRP4`** | 643 | 0.242 | **+3.254** | 0.00037 | +2.094 | 0.0558 | **-0.555** | +1.915 | +1.871 | **+4.341** | **Yes** |
| **`vsmc`** | 7,113 | 0.000 | **+2.390** | 0.00037 | +1.077 | 0.3722 | **-0.794** | **-1.103** | +3.142 | +2.833 | **Yes** |
| `Vas-capillary` | 6,917 | 0.000 | +1.906 | 0.00037 | +1.562 | 0.0114 | +0.377 | +1.278 | +1.071 | +1.972 | No |
| `Fibro-major` | 25,784 | 0.000 | +1.849 | 0.01520 | +1.722 | 0.0060 | +1.716 | -0.185 | +1.971 | +0.070 | No |
| `Vas-arterial` | 3,140 | 0.000 | +1.513 | 0.00037 | +1.123 | 0.0275 | **-1.029** | +1.523 | +0.452 | **+3.004** | No |
| `CD4-Tem` | 539 | 0.152 | -1.511 | 0.01222 | -0.283 | 0.2666 | -0.364 | +1.513 | -2.539 | -0.662 | No |
| `Lumsec-basal` | 1,841 | 0.091 | **-1.422** | 0.02407 | **-1.133** | 0.0362 | -0.919 | -1.135 | -1.249 | -1.464 | No |
| `Lumsec-prol` | 83 | 0.515 | **-1.200** | 0.02719 | -0.325 | 0.2992 | +0.070 | -0.057 | -1.101 | -1.228 | No (`droplow50pct` drops) |
| `cDC2` | 708 | 0.121 | -1.186 | 0.02146 | -0.672 | 0.0836 | +0.118 | -0.389 | -0.778 | -1.284 | No |
| `bmem_unswitched` | 49 | 0.636 | +0.458 | 0.03446 | +0.433 | 0.0275 | +0.752 | +0.816 | +0.942 | +1.006 | No (`droplow50pct` drops) |

### B. True Biological Signal vs. Batch-Confounded Variance

1. **Strongest biological (`breast cancer` vs `normal`) differences are sparse immune & luminal-HR shifts — NOT stromal mural cells:**
   - **Only 1 cell type out of 58 passes `p_bh_bio < 0.05`**: **`CD4-activated`** (`317` cells, `zero_frac = 0.515`, `Δ_bio = -1.76`, `p = 1.59e-5`, `p_bh = 0.00092`) — depleted in cancer within *both* `12h` (`-1.41`, `p = 0.013`) and `15m` (`-2.44`, `p = 0.0021`) and surviving `limma` (`-1.96`, `p = 8e-5`).
   - Next strongest non-confounded signals (`p_bh_bio = 0.055–0.076`): `Lymph-major` (`+1.89`, `p = 0.0030`), `Mast` (`-1.29`, `p = 0.0030`, `zero_frac = 0.667`), `CD8-activated` (`-0.98`, `p = 0.0038`, `zero_frac = 0.636`), `Macro-lipo` (`-1.32`, `p = 0.0067`), and **`LummHR-SCGB`** (`-2.04`, `p = 0.0092`, consistent in `12h = -1.68`, `15m = -2.15`, `limma = -2.07`).
2. **Epithelial compartments (`LummHR-*`, `Lumsec-*`, `basal`):**
   - Author annotations contain **no separate `"cancer"` label**; epithelial cells are grouped under the 11 HBCA luminal/basal phenotypes.
   - Only **`Lumsec-basal`** (`batch -1.42`, `bio -1.13`, `p = 0.036`) and **`Lumsec-prol`** (`batch -1.20`, `bio -0.33`, `p = 0.30`) show significant dissociation sensitivity (`15m > 12h`), both slightly depleted in cancer in the same direction as the `12h`-cancer enrichment. `basal` (`6,464` cells, `var = 3.74`) and `LummHR-major` (`2,898` cells, `var = 2.55`) show high inter-sample variance unrelated to either batch (`p_bh > 0.53`) or cancer-vs-normal (`p > 0.35`).
3. **Why `limma` drops `ECODA_authors` `bio ANOSIM` (`0.183 → 0.094`):**
   - `12 hour` digestion massively increases **stromal/vascular recovery** (`pericytes` +3.49, `Fibro-SFRP4` +3.25, `vsmc` +2.39, `Vas-capillary` +1.91, `Fibro-major` +1.85, `Vas-arterial` +1.51). Because `12h` is 70.6% cancer, **7 of the 11 batch-significant types appear positively or negatively shifted in marginal cancer-vs-normal comparisons** even though within-batch `Fibro-SFRP4` (`-0.56`), `pericytes` (`-0.50`), and `Vas-arterial` (`-1.03`) **disappear or reverse inside `12h`**.
   - Without protecting `disease` (strictly unsupervised correction), `limma` strips the dissociation-associated variance from `pericytes` (`var 4.07 → 1.65`), `Fibro-SFRP4` (`3.87 → 1.72`), `vsmc` (`3.23 → 1.31`), `Fibro-major` (`2.70 → 1.96`), `Vas-capillary` (`1.77 → 0.84`), and `Vas-arterial` (`1.57 → 0.81`). Once that batch-inflated stromal/vascular separation is removed, the true marginal biological separation across all 58 cell types in Euclidean CLR space drops from `0.183` to `0.094` (and `Pseudobulk` similarly drops `0.150 → 0.073`).
4. **Why `sigdiff_ctdrop` preserves more signal than `limma`, EXCEPT when `droplow50pct` is added:**
   - `ECODA_authors_sigdiff_ctdrop_0` (drops all 11 BH-significant batch types) lowers `batch ANOSIM` (`0.379 → 0.088`) while keeping `bio ANOSIM = 0.156` (still above `limma = 0.094` and `Pseudobulk = 0.073`).
   - `ECODA_authors_sigdiff_ctdrop_2` (drops only the 3 extreme stromal/mural types: `Fibro-SFRP4`, `vsmc`, `pericytes`) retains `bio ANOSIM = 0.180` and reduces `batch ANOSIM` (`0.379 → 0.218`).
   - Adding `droplow50pct` (`sigdiff_ctdrop_2_droplow50pct`, 16 types dropped) drops 13 sparse types — inadvertently deleting **`CD4-activated`** (`#1` bio hit, `zero_frac = 0.515`), **`Mast`** (`0.667`), and **`CD8-activated`** (`0.636`) — causing `bio ANOSIM` to fall to **`0.116`**.

---

## 3. Recommendation for User Decision Gate 1

**Recommend DROPPING `Breast_cancer` from the main batch-effect benchmark (leaving 6 clean batch cohorts: `Alzheimer`, `Covid19_PBMC`, `Joanito`, `Kidney_KPMP_full`, `Lupus_PBMC`, `Stephenson`), and documenting the exclusion transparently in `NOTES.md` + `Supp Table 5`**, for four converging scientific reasons:
1. **Small $N$ + residual confounding:** After removing confounded sites/assays/unknowns, only 33 samples remain with 70.6% vs 43.8% cancer imbalance across `12h` vs `15m`.
2. **Weak intrinsic biological signal dominated by dissociation:** Baseline biological ANOSIM is only `0.106–0.251` across *all* 8 methods, whereas dissociation batch ANOSIM is `0.379–0.863` (>2× larger), and high-variance epithelial subsets (`basal`, `LummHR-major`) barely distinguish cancer from normal.
3. **Bio/batch collinearity makes linear batch correction destroy signal across methods:** Unsupervised batch removal (`limma` in ECODA & Pseudobulk) drops biological ANOSIM in half (`0.183 → 0.094`, `0.251 → 0.024`, `0.150 → 0.073`) because `Fibro-major`, `Fibro-SFRP4`, `Vas-capillary`, `Vas-arterial`, `pericytes`, and `Lumsec-basal` vary jointly with dissociation and cancer proportion.
4. **Revision commitment (`revision-plan.docx` Point 10) remains solidly satisfied without it:** 11 benchmark + 5 new non-Stephenson batch cohorts + `Yan` cancer cohort = **17 unique cohorts** (≥16 required; 6 new cohorts ≥5 required).
*(Alternatively, if retained as a 7th cohort, keep `Breast_cancer` explicitly framed as a small-$N$, dissociation-confounded negative/cautionary edge case where `sigdiff_ctdrop_0` outperforms `limma`.)*
