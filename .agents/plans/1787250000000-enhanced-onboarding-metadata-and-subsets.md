# Implementation Plan — Enhanced Dataset Onboarding & Diagnostic Analytics

## Overview & Context

This implementation plan details the full technical design and execution strategy for upgrading the single-cell dataset onboarding diagnostic suite (`notebooks/dataset_onboarding/`) across all candidate cohorts. It ensures that every new dataset onboarded into ECODA undergoes thorough, unbiased evaluation of data integrity, biological signal, technical batch effects, and cross-study cell type harmonization.

---

## 1. Cell Type Harmonization & Granularity Evaluation (Section 1)

ECODA requires fine-grained, biologically meaningful cell type annotations that are **consistently harmonized across all batches/studies** in multi-cohort datasets. If a dataset is a compilation of multiple studies with disjoint or study-specific cell labels, fine annotations cannot be used for cross-study patient stratification.

### Implementation:
- Add `ou.cell_type_harmonization_check(obs, ct_cols, batch_col, sample_col)` to `onboarding_utils.py`.
- For each candidate annotation column (e.g. `cell_type`, `broad_cell_type`, `ann_fine`, `author_cell_type`, `subclass.l1`, etc.):
  - Count total unique cell types (`n_unique`).
  - Cross-tabulate with the primary study/batch variable (e.g., `study`, `dataset`, `Single cell sequencing platform`, `batch_cov`, `batch`, `assay`).
  - Calculate sharing statistics:
    - Number of cell types present across $\ge 80\%$ of batches vs single-batch-only.
    - Distribution of batches per cell type.
  - Classify annotation status:
    - **`Harmonized (Atlas-Wide)`**: Labels shared across all batches/studies (optimal for ECODA).
    - **`Partially Harmonized`**: Core cell lineages shared, but fine subtypes are batch-restricted.
    - **`Study-Specific / Unharmonized`**: Disjoint label sets per study (cannot be used for cross-study ECODA stratification).
- Present an executive comparison table and automated recommendation for the finest harmonized annotation level.

---

## 2. Metadata Column Categorization Schema

In each notebook, all `obs` columns are parsed and classified into 7 distinct roles:

1. **`Main Biological Condition (PILOT-GM-VAE / Primary Contrast)`**: Primary biological label for patient stratification (e.g., `Cognitive status` for Alzheimer, `disease` for Breast Cancer/Diabetes, `CoVID-19 severity` for Covid-19 PBMC, `condition.l1` for Kidney KPMP, `origin` for Lung, `Status` for Lupus PBMC, `patient_group` for Myocardial Infarction, `disease` for Parkinson).
2. **`Secondary / Demographic Biological Covariates`**: Tissue of origin, sub-phenotypes, pathology stages, anatomical location, and donor demographics (`Age`, `Sex`, `BMI`, `Braak stage`, `Smoking status`, `harmonized_ethnicity`, `anatomical_region` / CCF score).
3. **`Batch Effect Candidates (Technical Covariates)`**: Technical factors introducing non-biological variance:
   - Sequencing chemistry & platform (`assay`, `sequencing_platform`, `3' or 5'`).
   - Library & processing batches (`batch_cov`, `batch`, `Processing_Cohort`, `hospital`, `City`, `Brain_bank`).
   - Tissue handling & sampling protocols (`tissue_dissociation_protocol`, `tissue_sampling_method`, `donor_status`, `sample_preservation_method`, `tissue_type`).
4. **`Cell Type Annotations`**: Major cell lineages, sub-clusters, fine cell types, and ontology terms (`cell_type`, `broad_cell_type`, `author_cell_type`, `ann_fine`, `cell_type_ontology_term_id`).
5. **`Sample & Donor Identifiers`**: Sample, patient, donor, and specimen IDs (`sample_id`, `donor_id`, `patient`, `sampleID`, `specimen`).
6. **`Technical QC Metrics & Single-Cell Artifacts`**: Continuous cell-level QC parameters (`n_counts`, `n_features`, `percent_mito`, `percent_rb`, `doublet_score`, `cell_cycle_phase`).
7. **`Uninformative / Constant Columns`**: Constants (e.g., `organism = Homo sapiens`), join keys (`observation_joinid`), or binary study flags (`is_primary_data`).

---

## 3. UMAP Panel Enhancement (Sample-Colored Embedding)

- Update `ou.embed_and_umap_workflow()` to always generate a dedicated **UMAP panel colored by Sample ID** (`SAMPLE_COL`).
- Visual grid displays:
  - **Panel A:** Biological condition (`BIO_COL`)
  - **Panel B:** Primary batch candidate(s) (`BATCH_COLS`)
  - **Panel C:** Cell type annotation (`CT_COL`)
  - **Panel D:** Sample ID (`SAMPLE_COL`) — to visually detect donor-level clustering vs proper within-condition mixing.

---

## 4. Separation Heatmap Layout Overhaul (Section 5 Fix)

### Problem Identified:
- In datasets with many cell types and long labels (e.g., Lung, Alzheimer), the heatmap columns get squished horizontally into thin vertical bars, cell values overlap (`0.99 0.99 0.99`), and the colorbar compresses the plot.

### Redesign in `plot_separation_heatmap`:
- **Dynamic Proportional Sizing:**
  - `col_width = 1.8` inches per score column (`bio_separation`, `batch_<var>_separation`).
  - `label_margin = max(4.0, max_label_len * 0.14)` inches for y-axis labels.
  - Figure width = `label_margin + n_cols * col_width + 1.8` inches.
  - Row height = `0.38` inches per cell type row; Figure height = `max(5.0, n_rows * 0.38 + 1.5)` inches.
- **Visual Typography & Formatting:**
  - Explicit numeric formatting `fmt=".2f"`.
  - Contrast-aware annotation font (`fontsize=9`, `fontweight="bold"`).
  - Standalone colorbar formatting with `cbar_kws={"shrink": 0.5, "aspect": 20, "label": "LISI Separation (1=Separated, 0=Mixed)"}` so it does not compress heatmap cells.

---

## 5. Subsetting Strategy: Balanced $\min(N_{\text{samples}}, 20)$ Sample Allocation

### Budget:
- **Sample Count:** Target $\min(N_{\text{samples}}, 20)$ samples per dataset (e.g., 20 samples if cohort has $\ge 20$, or all available samples if $< 20$).
- **Cell Count:** 500 cells per sample (total $2,500 - 10,000$ cells), matching the validated `_debug` strategy.

### Joint Balancing Algorithm (`select_balanced_samples`):
- Filter candidate samples with $\ge 500$ cells (or $\ge 200$ cells if sample counts are small).
- Cross-tabulate candidate samples across the primary biological condition $\times$ suspected batch candidates (e.g. `assay`, `tissue`/`region`, `batch`).
- Round-robin sample selection ensuring that for each biological condition, samples from multiple batch levels are selected (e.g., Normal in Assay A + Normal in Assay B; Disease in Assay A + Disease in Assay B).
- This ensures cell-level LISI can score batch separation within shared cell types without collinearity dropping all tests.

---

## 6. Dataset-Specific Covariate Configurations

### 1. Alzheimer (SEA-AD)
- **Primary Bio:** `Cognitive status` (Dementia, Mild Cognitive Impairment, No dementia)
- **Secondary Bio:** `ADNC`, `Braak stage`, `CERAD score`, `Age at death`, `APOE4 status`, `sex`
- **Batch Candidates:** `assay` (10x v2 vs v3), `tissue_type` / `tissue`, `PMI`
- **Sample Col:** `donor_id` | **Cell Type:** `cell_type` / `Subclass`

### 2. Breast Cancer (Wu et al. 2021 / HBCA)
- **Primary Bio:** `disease` (`normal` vs `breast cancer`)
- **Secondary Bio:** `tissue_location`, `bmi_group`, `procedure_group`, `breast_density`, `age_group`
- **Batch Candidates:** `assay` (10x 3' v2 vs v3), `sequencing_platform`, `sample_source`, `suspension_dissociation_time`
- **Sample Col:** `donor_id` (or `sample_id`) | **Cell Type:** `cell_type` / `broad_cell_type`

### 3. Covid-19 PBMC (Ren et al. 2021)
- **Primary Bio:** `CoVID-19 severity` (`Control`, `Mild`, `Severe`, `Critical`)
- **Secondary Bio:** `Age`, `Sex`, `Outcome`, `Comorbidities`, `Sampling day`
- **Batch Candidates:** `Single cell sequencing platform`, `datasets`, `City`, `Sample type`
- **Sample Col:** `sampleID` (or `PatientID`) | **Cell Type:** `celltype` / `majorType`

### 4. Diabetes (Mouse Pancreatic Islets)
- **Primary Bio:** `disease` (`T2D` vs `ND`)
- **Secondary Bio:** `diabetes_model`, `age`, `cell_cycle_phase`
- **Batch Candidates:** `dataset`, `design`, `assay`
- **Sample Col:** `donor_id` (or `dataset__design__sample`) | **Cell Type:** `cell_type`

### 5. Kidney KPMP
- **Primary Bio:** `condition.l1` (`CKD`, `AKI`, `Healthy Reference`)
- **Secondary Bio:** `condition.l2`, `eGFR`, `hypertension`, `diabetes_history`, `sex`
- **Batch Candidates:** `assay` (10x 3' v3 vs v2 vs single-nucleus), `tissue_type`, `region.l1` (Cortex vs Medulla)
- **Sample Col:** `donor_id` | **Cell Type:** `subclass.l1` / `cell_type`

### 6. Lung Atlas (Sikkema et al. 2023 Nature Medicine, Fig 4)
- **Primary Bio:** `origin` (`tumor` vs `normal` / non-malignant) or `disease`
- **Secondary Bio:** `anatomical_region` (proximal vs distal / CCF score), `harmonized_ethnicity`, `BMI`, `Sex`, `Smoking status`, `age`, `tumor_stage`, `EGFR_mutation`
- **Batch Candidates:** `tissue_dissociation_protocol`, `tissue_sampling_method`, `donor_status`, `single_cell_platform` / `platform`, `cell_ranger_version`, `sequencing_platform`, `3' or 5'`, `dataset` / `study`
- **Sample Col:** `donor_id` (or `sample`) | **Cell Type:** `cell_type` / `ann_fine` / `cell_type_major`

### 7. Lupus PBMC (Perez et al. 2022)
- **Primary Bio:** `Status` (`SLE` vs `Control`)
- **Secondary Bio:** `SLE_status`, `Age`, `Sex`, `pop_cov` (ancestry)
- **Batch Candidates:** `batch_cov`, `Processing_Cohort`
- **Sample Col:** `sampleID` | **Cell Type:** `cell_types` / `ct_cov`

### 8. Myocardial Infarction (Kuppe et al. 2022)
- **Primary Bio:** `patient_group` (`Control`, `Myogenic`, `Fibrotic`, `Ischemic`)
- **Secondary Bio:** `sampleType`, `cell_subtype`
- **Batch Candidates:** `batch` (distinct library batches 'A' vs 'B')
- **Sample Col:** `patient` | **Cell Type:** `cell_type`

### 9. Parkinson (Kamath et al. 2022)
- **Primary Bio:** `disease` (`PD` vs `Control`)
- **Secondary Bio:** `path_braak_lb`, `PMI`, `RIN`, `sex`
- **Batch Candidates:** `Brain_bank`, `assay`, `tissue_type`
- **Sample Col:** `donor_id` | **Cell Type:** `cell_type`

---

## 7. Documentation Updates (`notebooks/dataset_onboarding/README.md` & `TODO.md`)

- Document concise summary findings per dataset in `README.md`:
  - Identified batch effect variables.
  - Degree of collinearity with biological condition.
  - Within-cell-type LISI mixing/separation.
  - Recommended cell type annotation column based on cross-batch harmonization.
- Track post-onboarding extension in `TODO.md`:
  - Cross-metadata collinearity matrices (Cramér's V) and atlas-wide LISI heatmaps on full HPC cohorts.

---

## Proposed Changes by Component

### Component A: Onboarding Utilities & Visualization
- [MODIFY] `notebooks/dataset_onboarding/onboarding_utils.py`:
  - Add `cell_type_harmonization_check(obs, ct_cols, batch_col, sample_col)`.
  - Add `categorize_obs_columns(obs, config)`.
  - Overhaul `plot_separation_heatmap(sep_df, name)` with dynamic width/height and non-compressing colorbars.
  - Update `embed_and_umap_workflow()` to add sample-colored UMAP panels.
  - Implement `select_balanced_samples()` for $\min(N_{\text{samples}}, 20)$ joint multi-condition multi-batch sample allocation.

### Component B: Candidate Dataset Notebooks
- [MODIFY] `notebooks/dataset_onboarding/dataset_check_*.qmd`:
  - Add Section 1 cell type harmonization comparison table.
  - Add metadata categorization table section (`ou.categorize_obs_columns`).
  - Wire dataset-specific primary/secondary bio columns and technical batch candidates (including Sikkema Lung Fig 4 covariates).
  - Update UMAP panel rendering to display sample-colored plots.

### Component C: Documentation
- [MODIFY] `notebooks/dataset_onboarding/README.md`:
  - Document Myocardial Infarction batch effect, Sikkema Lung covariates, harmonization criteria, and balanced subsetting findings.
- [MODIFY] `TODO.md`:
  - Track post-onboarding full-cohort collinearity synthesis.

---

## Verification Plan

1. **Local Render Verification:**
   - Execute `pixi run -e default quarto render notebooks/dataset_onboarding/dataset_check_<name>.qmd` for each dataset.
   - Verify that all 10 HTML reports generate with exit code 0, complete metadata categorization tables, cell type harmonization tables, sample-colored UMAPs, and readable, un-squeezed LISI separation heatmaps.
2. **Collinearity & Batch Evaluation Verification:**
   - Verify that balanced subsets allow LISI to evaluate batch candidates without dropping all cell types as "confounded".
