# Implementation Plan — Enhanced Dataset Onboarding & Diagnostic Analytics

**Status:** COMPLETED  
**Scope:** `notebooks/dataset_onboarding/`, `TODO.md`

---

## Overview & Context

This document details the completed implementation for the single-cell dataset onboarding diagnostic suite (`notebooks/dataset_onboarding/`) across all candidate cohorts. Every dataset onboarded into ECODA undergoes thorough, dual evaluation of:
1. **Gene Expression Batch Effects** (Unintegrated PCA + UMAP embeddings, cell-level LISI separation).
2. **Cell-Type Compositional Variance Partitioning** (Quantifying the fraction of total inter-sample variance explained by biological vs. technical covariates across cell types, following Sikkema et al. 2023 *Nature Medicine*, Fig. 4a).
3. **Cross-Study Cell Type Harmonization & Granularity** (Sharing across batches/studies).
4. **Metadata Categorization & Attribution** (Executive concise overview + detailed categorization across 7 standardized roles).

---

## 1. Cell Type Harmonization & Granularity Evaluation (Section 1)

ECODA requires fine-grained, biologically meaningful cell type annotations that are **consistently harmonized across all batches/studies** in multi-cohort datasets.

### Implemented Architecture:
- Implemented `ou.cell_type_harmonization_check(obs, ct_cols, batch_col, sample_col)` in `onboarding_utils.py`.
- For each candidate annotation column (e.g. `cell_type`, `broad_cell_type`, `ann_fine`, `author_cell_type`, `subclass.l1`, etc.):
  - Counts total unique cell types (`n_unique`).
  - Cross-tabulates with the primary study/batch variable (`study`, `dataset`, `Single cell sequencing platform`, `batch_cov`, `batch`, `assay`).
  - Calculates sharing statistics across batches (present in $\ge 80\%$ vs single-batch-only).
  - Classifies annotation status:
    - **`Harmonized (Atlas-Wide)`**: Labels shared across all batches/studies (optimal for ECODA).
    - **`Partially Harmonized`**: Core cell lineages shared, but fine subtypes are batch-restricted.
    - **`Study-Specific / Unharmonized`**: Disjoint label sets per study.
- Displays an untruncated comparison table with automated recommendations for the finest harmonized annotation level.

---

## 2. Metadata Column Categorization Schema

In each notebook, all `obs` columns are parsed and classified into 7 distinct roles:

1. **`Main Biological Condition (PILOT-GM-VAE / Primary Contrast)`**: Primary biological label for patient stratification (e.g., `Cognitive status` for Alzheimer, `disease` for Breast Cancer/Diabetes/Parkinson, `CoVID-19 severity` for Covid-19 PBMC, `condition.l1` for Kidney KPMP, `disease`/`origin` for Lung, `Status` for Lupus PBMC, `patient_group` for Myocardial Infarction).
2. **`Secondary / Demographic Biological Covariates`**: Tissue of origin, sub-phenotypes, pathology stages, anatomical location, and donor demographics (`Age`, `Sex`, `BMI`, `Braak stage`, `Smoking status`, `harmonized_ethnicity`, `anatomical_region` / CCF score).
3. **`Batch Effect Candidates (Technical Covariates)`**: Technical factors introducing non-biological variance:
   - Sequencing chemistry & platform (`assay`, `sequencing_platform`, `3' or 5'`).
   - Library & processing batches (`batch_cov`, `batch`, `Processing_Cohort`, `City`, `Brain_bank`, `library`).
   - Tissue handling & sampling protocols (`tissue_dissociation_protocol`, `tissue_sampling_method`, `donor_status`, `sample_preservation_method`, `tissue_type`, `PMI`).
4. **`Cell Type Annotations`**: Major cell lineages, sub-clusters, fine cell types, and ontology terms (`cell_type`, `broad_cell_type`, `author_cell_type`, `ann_fine`, `cell_type_ontology_term_id`).
5. **`Sample & Donor Identifiers`**: Sample, patient, donor, and specimen IDs (`sample_id`, `donor_id`, `patient`, `sampleID`, `specimen`).
6. **`Technical QC Metrics & Single-Cell Artifacts`**: Continuous cell-level QC parameters (`n_counts`, `n_features`, `percent_mito`, `percent_rb`, `doublet_score`, `cell_cycle_phase`).
7. **`Uninformative / Constant Columns`**: Constants (`organism = Homo sapiens`), join keys (`observation_joinid`), or binary study flags (`is_primary_data`).

### Implemented Visualization:
- Added `ou.summarize_obs_categories()` producing a concise executive summary table of assigned metadata categories, column counts, and comma-separated column lists.
- Global pandas display settings configured (`pd.set_option('display.max_colwidth', None)`, `pd.set_option('display.max_columns', None)`, etc.) to prevent cell truncation.

---

## 3. Variance Partitioning Across Covariates (Section 5)

Evaluates sample-level centered log-ratio (CLR) abundances ($+0.5$ pseudocount zero-imputation) and partitions variance across covariates following Sikkema et al. 2023 (*Nature Medicine* 29:1563–1577, Fig. 4a).

### Implemented Methods:
1. **`compute_variance_partition(clr_df, meta_df, sample_col, bio_cols, tech_cols)`:**
   - Fits linear models $\text{CLR}_{ct} \sim C(X)$ for each covariate $X$ across each cell type $ct$ and the whole cohort average (`Whole atlas`).
   - Returns `(var_tech_df, var_bio_df)` containing the fraction of total variance explained ($R^2$).
   - Masks non-evaluated / single-level covariates as `np.nan` (rendered in grey).
2. **`plot_variance_partition_heatmap(var_tech, var_bio, title=None, vmax=0.40)`:**
   - Dual-panel heatmap displaying:
     - **Left panel:** `Covariate (technical)` with 90° rotated labels.
     - **Right panel:** `Covariate (biological)` with 90° rotated labels.
     - **Y-axis:** `Cell type` (with `Whole atlas` at the top).
     - **Colorbar:** `Fraction of total variance` using `Reds` colormap (`0.0` to `0.40`, `extend='max'`).
     - **Grey fill:** Inapplicable or constant covariates.
3. **Multi-Variable Joint OLS Linear Regression (`compute_compositional_joint_lm`):**
   - Fits all covariates simultaneously ($\text{CLR} \sim \text{cov}_1 + \text{cov}_2 + \dots$), computing Type-II ANOVA sum of squares and FDR-adjusted $p$-values.

---

## 4. Unintegrated Embeddings & UMAP Panel Enhancement (Section 4)

- Fresh unintegrated PCA+UMAP computed directly from raw integer counts (`10k` normalization, `log1p`, 2000 HVGs, 50 PCs, Scanpy default UMAP).
- Integrated `SAMPLE_COL` coloring into the UMAP panels alongside biological conditions, batch candidates, and cell types.

---

## 5. Separation Heatmap Layout Overhaul (Section 6)

- Dynamic proportional sizing in `plot_separation_heatmap()`:
  - `col_width = 1.8` inches per score column.
  - `label_margin = max(4.0, max_label_len * 0.14)` inches for y-axis labels.
  - Figure width = `label_margin + n_cols * col_width + 1.8` inches.
  - Figure height = `max(5.0, n_rows * 0.38 + 1.5)` inches.
- Full formatting with explicit `fmt=".2f"` and contrast-aware bold typography.

---

## 6. Subsetting Strategy: Balanced $\min(N_{\text{samples}}, 20)$ Sample Allocation

- Implemented `ou.select_balanced_samples()`:
  - Targets $\min(N_{\text{samples}}, 20)$ samples with 500 cells per sample.
  - Stratifies across biological conditions and round-robins across batch candidates to guarantee multi-batch representation per condition.

---

## 7. Master Summary Table in `README.md`

Overhauled `notebooks/dataset_onboarding/README.md` with the 17-column master summary table:
- **`Bio Condition (PILOT-GM-VAE)`** vs **`Bio Condition (ECODA)`**
- **`Batch Condition (PILOT-GM-VAE)`** vs **`Batch Sequencing (ECODA)`** vs **`Batch Sample Prep (ECODA)`**
- **`Suitable for Auto-Annotation`** (`Yes` vs `No`)
- **`Count Check`** (`PASS (raw.X: raw integer counts)` vs `NOTE (X: log-normalized)`)
- **`Recommended Use`** (`benchmark` vs `batch-effect` vs `reference`)

---

## 8. Verification Results

All 10 Quarto onboarding notebooks rendered with **100% success** (exit code 0):
- `dataset_check__debug.qmd`
- `dataset_check_Alzheimer.qmd`
- `dataset_check_Breast_cancer.qmd`
- `dataset_check_Covid19_PBMC.qmd`
- `dataset_check_Diabetes.qmd`
- `dataset_check_Kidney_KPMP.qmd`
- `dataset_check_Lung.qmd`
- `dataset_check_Lupus_PBMC.qmd`
- `dataset_check_Myocardial_infarction.qmd`
- `dataset_check_Parkinson.qmd`
