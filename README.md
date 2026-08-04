[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the results and figures from the
paper: **"Cell type composition drives patient stratification in single-cell
RNA-seq cohorts"**.

## **Overview**

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization
of cellular heterogeneity, but summarizing this data for cohort-level analysis
remains a challenge. Using **11 scRNA-seq cohorts** (697 samples) across
different biological conditions, we benchmarked seven state-of-the-art sample
representation methods—**MOFA+, scITD, GloScope, GloProp, MrVI, PILOT,
scPoli**—plus two baselines (**Pseudobulk** and cell-type composition via
**ECODA**). The benchmark evaluates their ability to recover known biological
groupings in a fully unsupervised setting.

### **Key Findings**

-   **Performance:** Centered log-ratio (CLR)-transformed cell-type proportions
    (ECODA) consistently match or outperform more complex methods in recovering
    known biological groupings in an unsupervised setting.
-   **Efficiency:** ECODA requires orders of magnitude fewer computational
    resources and produces embeddings in seconds.
-   **Robustness:** The approach is highly robust to technical batch effects and
    various cell-type annotation strategies. On the Joanito dataset, ECODA
    achieved a batch ANOSIM of 0.041 (effectively no batch separation) while
    preserving biological signal (biological ANOSIM = 0.640); in contrast,
    pseudobulk showed strong batch separation (batch ANOSIM = 0.706).
-   **Interpretability:** Biological stratification is often driven by a small
    subset of highly variable cell types (HVCs), providing direct mechanistic
    insights.

The **scECODA** R package for scalable cohort-level analysis is available at
[github.com/carmonalab/scECODA](https://github.com/carmonalab/scECODA).


## **Usage**

This section explains how to reproduce this paper's results and figures and the general workflow of the analysis.

### **Repository Contents**

- `docs/ARCHITECTURE.md`: Full pipeline architecture, call flow, and module documentation.

- `datasets.json`: Centralized dataset metadata (sample/label columns, subsetting rules, batch information).

- `notebooks/`: .rmd notebooks for:
  - Benchmark analysis (`notebooks/benchmark_analysis.rmd`)
  - Batch effect analysis (`notebooks/batch_effect_analysis.rmd`)
  - Per-dataset QC filtering (`notebooks/QC_filtering/`)

- `src`: Source code for:
  - Raw-data staging from NAS to HPC scratch (`src/1_stage_data/`)
  - Dataset-specific preprocessing scripts (`src/2_dataset_specific_preprocessing/`)
  - Standardized scRNA-seq preprocessing pipeline (Python/Scanpy) (`src/3_scrnaseq_preprocessing/`)
  - Cell type annotation with scATOMIC + HiTME (`src/4_cell_type_annotation/`)
  - Scripts to run the benchmarked methods (`src/5_run_benchmark_methods/`)
  - Various utility functions (`src/utils/`), mainly for the benchmark pipeline analysis
  - HPC SLURM configuration, e.g. paths (`src/slurm_config.sh`)

### **Reference data**

- `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`: Ensembl 105 human gene reference
  (GRCh38.p13) used for gene-name standardization in the preprocessing pipeline
  (`src/gene_utils.py`). Originally retrieved on 14.02.2022 from the `aux/` folder
  of the carmonalab `scRNAseq_data_processing` repository:
  https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz


### **Installation**

1. Install [Pixi](https://pixi.sh) — the project's package and environment manager.
2. Clone the repository:
   ```bash
   git clone <repo-url> && cd ECODA_paper
   ```
3. Install dependencies:
   ```bash
   pixi install
   ```
   This creates the default `py-cpu` environment (macOS / development). For HPC
   with CUDA, use the `py-cuda13` environment instead:
   ```bash
   pixi install --environment py-cuda13
   ```
4. Install R packages (Seurat, anndataR, SignatuR, scATOMIC deps, scITD deps,
   HiTME deps, and all benchmark method packages) via the chained setup task:
   ```bash
   pixi run setup
   ```


### **Workflow**

The analysis proceeds through four stages:

- **Stage 1 — QC Filtering:** Done manually in per-dataset .rmd notebooks in
  `notebooks/QC_filtering/`.
- **Stage 2 — Preprocessing + Cell Type Annotation:**
  - **Raw-data staging** (`src/1_stage_data/`): login-node script that stages raw
    inputs from NAS to HPC scratch (`./src/1_stage_data/1_stage_data.sh`).
  - **Dataset-specific preprocessing** (`src/2_dataset_specific_preprocessing/`):
    per-step sbatch jobs (e.g. `_create_combinedpbmc_dataset.py`,
    `_create_joanito_batch_col.R`) submitted in parallel via the `1_submit_hpc.sh`
    dispatcher, run after staging and before the preprocess array.
  - **Preprocessing** (`src/3_scrnaseq_preprocessing/`): Standardized preprocessing pipeline (Python/Scanpy):
    - Filter cells (min_genes=100) and genes (min_cells=3)
    — Sample/gene name standardization
    - Normalize counts and log1p
    - HVG selection
    - Scale and run PCA
    - Harmony integration
    - Leiden clustering
  - **Cell Type Annotation** (`src/4_cell_type_annotation/`): HPC-parallelized scATOMIC + HiTME
    annotation via SLURM array jobs.
- **Stage 3 — Benchmark Analysis:** Render `notebooks/benchmark_analysis.rmd` in RStudio.
  Python benchmark methods (MrVI, PILOT, scPoli) are invoked automatically via
  rpy2. The R pipeline orchestrates all method processors, scoring metrics, and
  figure generation.
- **Stage 4 — Batch Effect Analysis:** Render `notebooks/batch_effect_analysis.rmd` in
  RStudio.

**HPC execution:** Submit SLURM array jobs for:
- **Raw-data staging** via `src/1_stage_data/1_stage_data.sh` (login node, NAS → scratch):
  ```bash
  ./src/1_stage_data/1_stage_data.sh
  ```
- **Dataset-specific preprocessing** via `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
  (submits all per-step sbatch jobs in parallel, waits, reports via sacct):
  ```bash
  sbatch src/2_dataset_specific_preprocessing/1_submit_hpc.sh
  ```
- **Preprocessing** via `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (submits array + syncs results):
  ```bash
  sbatch src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh
  ```
- **Cell type annotation** via `src/4_cell_type_annotation/` (all datasets by default; an optional dataset name restricts to one dataset):
  ```bash
  ./src/4_cell_type_annotation/1_prepare_chunks.sh             # production: 5 samples/chunk
  ./src/4_cell_type_annotation/1_prepare_chunks.sh test        # test: 1 sample/chunk
  ./src/4_cell_type_annotation/2_submit_hpc_array.sh
  ```
- **Benchmark methods** via `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh` (stages data + submits array + syncs results):

See the [Architecture documentation](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b)
for more details on workflow and usage.

#### **Expected Outputs**

- `.feather` files — cross-language distance matrices and embeddings produced
  by Python benchmark methods and consumed by R processors.
- Publication figures — MDS plots, PCA biplots, benchmark bar charts,
  separation metric heatmaps, and transformation analysis panels.
- Execution time logs documenting per-method and per-dataset runtime.


## Reference

If you use ECODA or this benchmark code in your research, please cite our
preprint:

**Cell type composition drives patient stratification in single-cell RNA-seq
cohorts.** Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*. doi:
[10.64898/2026.03.27.714811v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1)
