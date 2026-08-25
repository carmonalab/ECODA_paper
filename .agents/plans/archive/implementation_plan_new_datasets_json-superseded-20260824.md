# SUPERSEDED — 2026-08-25

This provisional plan is retained for provenance only. It is superseded by
the checked-in `datasets.json`, the explicit onboarding roles in
`notebooks/dataset_onboarding/dataset_specs.py`, and regenerated audit JSON.
Technical batch columns remain candidates until the uncorrected evidence gate
and explicit user confirmation.

## Final user-confirmed registry decisions

| Dataset | Sample | Label | Low-resolution role | High-resolution role | Annotation source |
| --- | --- | --- | --- | --- | --- |
| Alzheimer | `donor_id` | `Cognitive status` | `Subclass` | `Supertype` | author / author |
| Breast cancer | `sample_id` | `disease` | `broad_cell_type` | `author_cell_type` | author / author |
| Covid-19 PBMC | `sampleID` | `CoVID-19 severity` | `majorType` | `celltype` | author / author |
| Diabetes | `donor_id` | `disease` | `cell_type` | `cell_type_reannotatedIntegrated` | author / author |
| Kidney KPMP | `specimen` | `condition.l1` | `subclass.l1` | `subclass.l3` | author / author |
| Lung | `sample` | `origin` | `ann_coarse` | `ann_fine` | author / author |
| Lupus PBMC | `sampleID` | `Status` | `layer1` | `layer2` | HiTME / HiTME |
| Myocardial infarction | `orig_ident` | `patient_group` | `cell_type` | `cell_subtype` | author / author |
| Parkinson | `donor_id` | `disease` | `cell_type` | Leiden res-5 (view-specific) | author / Leiden |

Lung is the `platform == "10x"` registry subset. Lupus remains tissue
`"Blood"`. No new cohort has a selected batch column at this stage; the
corrected view is declared but fail-closed until confirmation.

# Implementation Plan: Onboarding New Datasets into `datasets.json`

This plan proposes the exact configuration for integrating the 9 new datasets (from the Joodaki et al. 2025 / PILOT-GM-VAE cohort, PMID 41097818) into `datasets.json`, answers all specific technical questions regarding cell type annotations and batch effect columns, and defines the phased rollout strategy.

---

## 1. Batch Effect Architecture: Expression vs. Composition Batch Keys

> [!TIP]
> **Assessment of the User's Proposal:**
> Your proposed **Standardized Dataset Registry in `datasets.json`** with `columns.batch_expression` and `columns.batch_composition` (with fallback to `columns.batch`) is **conceptually sound and strongly recommended**.

### Why Separate Expression & Composition Batch Keys?
1. **Gene Expression Methods (DESeq2, limma, Harmony, Scanorama, MrVI):**
   - Sensitive to molecular and technical variations across the $\sim 20,000$ gene transcriptome (sequencing platform chemistry, library preparation kit, read depth, capture efficiency, e.g., `10x 3' v2` vs `10x 3' v3` or `NovaSeq` vs `HiSeq`).
2. **Cell Type Compositional Methods (ECODA, GloProp, propeller, scCODA):**
   - Sensitive to macroscopic cell isolation and sampling biases (tissue dissociation protocol, enzymatic digestion time, sampling region/site, sorting gates, tissue preservation method, e.g., `Cortex` vs `Medulla`, `Frozen` vs `Fresh`, `Biopsy` vs `Resection`).

### Proposed Phased Rollout
1. **Phase 1 (This Step):** Onboard all 9 datasets into `datasets.json` in clean baseline mode (`use_for_benchmark: true`, uncorrected benchmark views) with primary candidate batch columns.
2. **Phase 2:** Execute uncorrected preprocessing and baseline benchmarking across the new cohorts on HPC.
3. **Phase 3:** Using empirical LISI separation scores and Sikkema-style variance partitioning computed on full cohorts, finalize specific `columns.batch_expression` and `columns.batch_composition` keys for dedicated batch-effect views.

---

## 2. Dataset-by-Dataset Technical Answers & Proposed Columns

Below is the detailed analysis addressing every question raised for each dataset:

---

### A. Alzheimer (`SEAAD_Alzheimer.h5ad`)
*Gabitto et al. 2024 Nat Neurosci (PMID 42486312) — 1,395,601 cells, 83 donors*

- **Cell Type Resolution:**
  - `Class` (3 types): Coarse lineages (`Glutamatergic`, `GABAergic`, `Non-neuronal`).
  - `Supertype` (18 types) & `Subclass` (18 types): In SEA-AD, `Supertype` and `Subclass` have a **1-to-1 correspondence** with the 18 official Allen Brain Institute cell subclasses (`Astrocyte`, `Microglia-PVM`, `Oligo`, `OPC`, `Endothelial`, `VLMC`, `L2/3 IT`, `L4 IT`, `L5 IT`, `L6 IT`, `L5 ET`, `L5/6 NP`, `L6 CT`, `L6b`, `Pvalb`, `Sst`, `Lamp5`, `Sncg`, `Vip`, `Chandelier`).
  - **Proposed:** `cell_type_low_res: "Class"`, `cell_type_high_res: "Supertype"`.
- **Sample & Bio Label:**
  - `sample: "donor_id"`, `label: "Cognitive status"` (`Dementia`, `No dementia`, `Reference`).
- **Auto-Annotation:** Flag `not_suitable_for_auto_annotation: ["hitme", "scatomic"]` (brain tissue lacks tumor/immune lineages).

---

### B. Breast Cancer (`BreastCncr_processed.h5ad`)
*Wu et al. 2021 Nat Genet (PMID 34493872) — 714,331 cells, 126 donors / 167 samples*

- **Cell Type Resolution & Harmonization:**
  - **Why did the `.qmd` report note *"Sharing Status: Partially Harmonized"*?
    The diagnostic subset sampled 20 donors. In Wu et al., the fine-grained `author_cell_type` (58 subtypes, e.g. `Fibro-prematrix`, `LummHR-SCGB`, `CD8-Trm`) and `cell_type` (39 ontology subtypes) contain rare cell states that appear only in subsets of patients or specific tumor histologies.
    However, **`broad_cell_type` contains 10 major lineages** (`fibroblasts`, `basal`, `vascular`, `tcells`, `lumhr`, `lumsec`, `pericytes`, `myeloid`, `bcells`, `lymphatic`) that are **100% shared across normal and cancer samples**.
  - **Proposed:** `cell_type_low_res: "broad_cell_type"` (10 lineages), `cell_type_high_res: "cell_type"` (39 harmonized ontology types).
- **Sample & Bio Label:**
  - `sample: "donor_id"`, `label: "disease"` (`normal` vs `breast cancer`).
- **Batch Candidates:** `sequencing_platform` (expression), `sample_preservation_method` / `suspension_dissociation_time` (composition).

---

### C. Covid-19 PBMC (`Covid19_Ren2021.h5ad`)
*Ren et al. 2021 Cell (PMID 33657410) — 993,171 cells, 151 donors*

- **Cell Type Resolution:**
  - `majorType` (6 types): `T/NK`, `B`, `Mono`, `DC`, `Plasma`, `Other`.
  - `celltype` (10 types): `T`, `NK`, `B`, `CD14_Mono`, `CD16_Mono`, `cDC`, `pDC`, `Plasma`, `Platelet`, `Megakaryocyte`.
  - **Proposed:** `cell_type_low_res: "majorType"`, `cell_type_high_res: "celltype"`.
- **Sample & Bio Label:**
  - `sample: "sampleID"`, `label: "CoVID-19 severity"` (`Healthy`, `Mild`, `Moderate`, `Severe`, `Critical`).
- **Batch Candidates:** `Single cell sequencing platform` (expression), `datasets` / `City` / `Sample type` (composition).

---

### D. Diabetes Mouse Islet (`diabetes.h5ad`)
*Hrovatin et al. 2023 Nat Metab (PMID 37697055) — 264,235 cells, 52 samples*

- **Cell Type Resolution:**
  - `cell_type` (13 types): Standard islet cell ontology (`pancreatic A cell`, `type B pancreatic cell`, `pancreatic D cell`, `pancreatic PP cell`, `endothelial cell`, `pancreatic stellate cell`, `pancreatic ductal cell`, `pancreatic acinar cell`, `hematopoietic cell`, etc.).
  - `cell_type_reannotatedIntegrated` (20 types): Hrovatin et al.'s integrated consensus cell types across all diabetes models (`alpha`, `beta`, `delta`, `gamma`, `endothelial`, `stellate a.`, `stellate q.`, `immune`, `E endo.`, `E non-endo.`, etc.).
  - **Proposed:** `cell_type_low_res: "cell_type"`, `cell_type_high_res: "cell_type_reannotatedIntegrated"`.
- **Sample & Bio Label:**
  - `sample: "dataset__design__sample"`, `label: "disease"` (`normal`, `type 1 diabetes mellitus`, `type 2 diabetes mellitus`, `endocrine pancreas disorder`).
- **Auto-Annotation:** Flag `not_suitable_for_auto_annotation: ["hitme", "scatomic"]` (mouse gene nomenclature).

---

### E. Kidney KPMP (`Kidney_KPMP.h5ad`)
*Lake et al. 2023 Nature (PMID 41648348) — 104,314 cells, 45 specimens*

- **What is the difference between `subclass.l3` and `subclass.full`?**
  - They represent the **exact same 63–64 cell subtypes**:
    - `subclass.l3` uses the **official KPMP short 3-to-5 letter codes** (`aPT`, `dFIB`, `PL`, `dPT`, `aTAL2`, `NKC/T`, `ncMON`, `PT-S1/2`, `T`, `dOMCD-PC`, `B`, `EC-AEA`, `C-TAL`, `dCNT`).
    - `subclass.full` uses the **expanded descriptive clinical names** (`Adaptive / Maladaptive / Repairing Proximal Tubule Epithelial Cell`, `Monocyte-derived Cell`, `Degenerative Fibroblast`, `Plasma Cell`, `Natural Killer Cell / Natural Killer T Cell`, etc.).
  - For coarse resolution, `subclass.l1` defines 15 major nephron/immune compartments (`PT`, `IMM`, `FIB`, `TAL`, `PC`, `EC`, `CNT`, `IC`, `PEC`, `VSM/P`, `DCT`, `PapE`, `DTL`, `POD`, `ATL`).
  - **Proposed:** `cell_type_low_res: "subclass.l1"`, `cell_type_high_res: "subclass.l3"` (clean short codes matching KPMP publications).
- **Sample & Bio Label:**
  - `sample: "specimen"`, `label: "condition.l1"` (`Healthy Reference`, `Chronic Kidney Disease`, `Acute Kidney Injury`).

---

### F. Lupus PBMC (`Lupus_Perez2022.h5ad`)
*Perez et al. 2022 Science (PMID 35389779) — 1,263,676 cells, 261 donors*

- **Cell Type Resolution:**
  - `cell_types` (11 types): Broad PBMC lineages (`T4`, `T8`, `B`, `NK`, `cM`, `ncM`, `cDC`, `pDC`, `Prolif`, `PB`, `Progen`).
  - `ct_cov` (14 types): Author fine-grained states (`T4_naive`, `T4_em`, `T8_naive`, `T8_em`, `CytoT_GZMK+`, `CytoT_GZMH+`, `B_naive`, `B_mem`, `B_plasma`, `NK_dim`, `NK_bright`, `cM`, `ncM`, `Progen`).
  - **Proposed:** `cell_type_low_res: "cell_types"`, `cell_type_high_res: "ct_cov"`. (HiTME will also run as extended auto-annotation).
- **Sample & Bio Label:**
  - `sample: "sampleID"`, `label: "Status"` (`Healthy` vs `Case`).

---

### G. Lung Atlas (`lungatlas.h5ad`)
*Sikkema et al. 2023 Nat Med (PMID 37291214) — 941,504 cells, 556 samples / 318 donors*

- **Harmonization & Preventing Disease Leakage:**
  > [!IMPORTANT]
  > **Crucial Finding:** Your concern about annotation bias is 100% correct!
  > - In `cell_type_tumor` (51 types), cancer cells are labeled by condition (`Tumor cells LUAD`, `Tumor cells LUSC`). If evaluated against disease status, this causes complete artificial circular separation because those cell types exist only in one disease!
  > - In contrast, **`ann_coarse` (12 lineages)** and **`cell_type_major` (24 lineages)** unbiasly unify all malignant cells under a single harmonized label: `"Tumor cells"`. All other lineages (`T cell CD4`, `T cell CD8`, `Macrophage alveolar`, `Monocyte`, `B cell`, `Endothelial cell`, `Stromal`, etc.) are present across all conditions and tissue origins (`normal`, `normal_adjacent`, `tumor_primary`).
  - **Proposed:** `cell_type_low_res: "ann_coarse"`, `cell_type_high_res: "cell_type_major"`.
- **Sample & Bio Label:**
  - `sample: "sample"`, `label: "disease"` (`normal`, `chronic obstructive pulmonary disease`, `lung adenocarcinoma`, `squamous cell lung carcinoma`, `non-small cell lung carcinoma`).

---

### H. Myocardial Infarction (`Myocardial_Infarc_2.h5ad`)
*Kuppe et al. 2022 Nature (PMID 35948637) — 132,888 cells, 23 patients*

- **What is the difference between `cell_subtype` and `cell_subtype2`?**
  - Both contain 33 fine subtypes. The only difference is naming conventions:
    - `cell_subtype2` prefixes the cluster index and gene markers (e.g. `Fib1_SCARA5`, `Fib2_Myofib`, `Fib3_C7`, `Fib4_COL15A1`, `CD_8`, `CD_4`).
    - `cell_subtype` uses standard clean anatomical names (e.g. `Fib_SCARA5`, `Myofib`, `Fib_0`, `Fib_3`, `CD8`, `CD4`).
  - `cell_type` provides the 11 major cardiac lineages (`Endo`, `Fib`, `CM`, `Neuronal`, `Myeloid`, `PC`, `prolif`, `vSMCs`, `Mast`, `Lymphoid`, `Adipo`).
  - **Proposed:** `cell_type_low_res: "cell_type"`, `cell_type_high_res: "cell_subtype"`.
- **Sample & Bio Label:**
  - `sample: "patient"`, `label: "patient_group"` (`fibrotic`, `myogenic`, `ischemic`).

---

### I. Parkinson (`Parkinson.h5ad`)
*Kamath et al. 2022 Nat Neurosci (PMID 35513515) — 2,096,155 cells, 97 donors*

- **What is the difference between `cell_type` and `derived_class2`?**
  - They are a **direct 1-to-1 mapping**:
    - `derived_class2` (11 types): Short author acronyms (`Oligo`, `EN`, `Astro`, `IN`, `OPC`, `Micro_PVM`, `Endo`, `Ependymal`, `Mural`, `DMNX_Neu`, `Adaptive`).
    - `cell_type` (11 types): Full Cell Ontology terms (`oligodendrocyte`, `glutamatergic neuron`, `astrocyte`, `GABAergic neuron`, `oligodendrocyte precursor cell`, `central nervous system macrophage`, `endothelial cell`, `ependymal cell`, `mural cell`, `central nervous system neuron`, `leukocyte`).
  - **Proposed:** `cell_type_low_res: "derived_class2"`, `cell_type_high_res: "cell_type"`.
- **Sample & Bio Label:**
  - `sample: "donor_id"`, `label: "disease"` (`normal` vs `Parkinson disease`).
- **Auto-Annotation:** Flag `not_suitable_for_auto_annotation: ["hitme", "scatomic"]` (brain tissue).

---

## 3. Proposed Additions to `datasets.json`

```json
  "Alzheimer": {
    "display_name": "Alzheimer (SEA-AD)",
    "file_names": "SEAAD_Alzheimer.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Brain",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": false,
    "not_suitable_for_auto_annotation": ["hitme", "scatomic"],
    "columns": {
      "sample": "donor_id",
      "label": "Cognitive status",
      "batch": "assay",
      "cell_type_low_res": "Class",
      "cell_type_high_res": "Supertype"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "SEAAD_Alzheimer.h5ad",
        "output_file_name": "SEAAD_Alzheimer_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Breast_cancer": {
    "display_name": "Breast cancer (Wu 2021)",
    "file_names": "BreastCncr_processed.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Breast tumor",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "donor_id",
      "label": "disease",
      "batch": "sequencing_platform",
      "cell_type_low_res": "broad_cell_type",
      "cell_type_high_res": "cell_type"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "BreastCncr_processed.h5ad",
        "output_file_name": "BreastCncr_processed_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Covid19_PBMC": {
    "display_name": "Covid-19 PBMC (Ren 2021)",
    "file_names": "Covid19_Ren2021.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Blood",
    "normal_tissue": true,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "sampleID",
      "label": "CoVID-19 severity",
      "batch": "Single cell sequencing platform",
      "cell_type_low_res": "majorType",
      "cell_type_high_res": "celltype"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Covid19_Ren2021.h5ad",
        "output_file_name": "Covid19_Ren2021_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Diabetes": {
    "display_name": "Diabetes mouse islet (Hrovatin 2023)",
    "file_names": "diabetes.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Pancreas",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "not_suitable_for_auto_annotation": ["hitme", "scatomic"],
    "columns": {
      "sample": "dataset__design__sample",
      "label": "disease",
      "batch": "dataset",
      "cell_type_low_res": "cell_type",
      "cell_type_high_res": "cell_type_reannotatedIntegrated"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "diabetes.h5ad",
        "output_file_name": "diabetes_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Kidney_KPMP": {
    "display_name": "Kidney KPMP (Lake 2023)",
    "file_names": "Kidney_KPMP.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Kidney",
    "normal_tissue": true,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "specimen",
      "label": "condition.l1",
      "batch": "library",
      "cell_type_low_res": "subclass.l1",
      "cell_type_high_res": "subclass.l3"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Kidney_KPMP.h5ad",
        "output_file_name": "Kidney_KPMP_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Lupus_PBMC": {
    "display_name": "Lupus PBMC (Perez 2022)",
    "file_names": "Lupus_Perez2022.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Blood",
    "normal_tissue": true,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "sampleID",
      "label": "Status",
      "batch": "batch_cov",
      "cell_type_low_res": "cell_types",
      "cell_type_high_res": "ct_cov"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Lupus_Perez2022.h5ad",
        "output_file_name": "Lupus_Perez2022_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Lung": {
    "display_name": "Lung Atlas (Sikkema 2023)",
    "file_names": "lungatlas.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Lung",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "sample",
      "label": "disease",
      "batch": "dataset",
      "cell_type_low_res": "ann_coarse",
      "cell_type_high_res": "cell_type_major"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "lungatlas.h5ad",
        "output_file_name": "lungatlas_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Myocardial_infarction": {
    "display_name": "Myocardial infarction (Kuppe 2022)",
    "file_names": "Myocardial_Infarc_2.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Heart",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": true,
    "columns": {
      "sample": "patient",
      "label": "patient_group",
      "batch": "batch",
      "cell_type_low_res": "cell_type",
      "cell_type_high_res": "cell_subtype"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Myocardial_Infarc_2.h5ad",
        "output_file_name": "Myocardial_Infarc_2_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  },
  "Parkinson": {
    "display_name": "Parkinson (Kamath 2022)",
    "file_names": "Parkinson.h5ad",
    "folder_name": "JooM_2025_41097818",
    "tissue": "Brain",
    "normal_tissue": false,
    "use_for_benchmark": true,
    "use_for_batch_effect": false,
    "not_suitable_for_auto_annotation": ["hitme", "scatomic"],
    "columns": {
      "sample": "donor_id",
      "label": "disease",
      "batch": "Brain_bank",
      "cell_type_low_res": "derived_class2",
      "cell_type_high_res": "cell_type"
    },
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Parkinson.h5ad",
        "output_file_name": "Parkinson_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  }
```

---

## 4. Verification Plan

1. Validate JSON syntax and structure with `jq . datasets.json`.
2. Verify Python parser: `pixi run python -c "from src.utils.py.datasets_io import read_datasets_json; ds = read_datasets_json(); print('Loaded', len(ds), 'datasets')"`
3. Verify R parser: `pixi run Rscript -e "source('src/utils/datasets_io.R'); ds <- read_datasets_json(view='benchmark_analysis'); cat('Benchmark datasets:', length(ds), '\n')"`
4. Verify data staging resolution test: test `1_stage_data.sh` staging dry-run over the new keys.
