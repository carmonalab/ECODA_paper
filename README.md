[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the results and figures from the
paper: **"Cell type composition drives patient stratification in single-cell
RNA-seq cohorts"**.

## **Overview**

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization
of cellular heterogeneity, but summarizing this data for cohort-level analysis
remains a challenge. Using **14 scRNA-seq datasets** (11 benchmark +
batch-effect cohorts; see [datasets.json](datasets.json) for the authoritative
list), we benchmarked seven state-of-the-art sample representation
methods—**MOFA+, scITD, GloScope, GloProp, MrVI, PILOT, scPoli**—plus two
baselines (**Pseudobulk** and cell-type composition via **ECODA**). The
benchmark evaluates their ability to recover known biological groupings in a
fully unsupervised setting.

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

### **Repository Contents**

- `docs/ARCHITECTURE.md`: Full pipeline architecture, call flow, and module documentation.
- `datasets.json`: Centralized dataset metadata (sample/label columns, subsetting rules, batch information).
- `notebooks/`: benchmark analysis, batch effect analysis, and per-dataset QC filtering notebooks (`.rmd`).
- `src`: Source code for the pipeline stages — staging (`src/1_stage_data/`), dataset-specific preprocessing (`src/2_dataset_specific_preprocessing/`), standardized preprocessing (`src/3_scrnaseq_preprocessing/`), cell type annotation (`src/4_cell_type_annotation/`), benchmark methods (`src/5_run_benchmark_methods/`), utilities (`src/utils/`), and HPC SLURM config (`src/slurm_config.sh`).

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
   with CUDA, use `pixi install --environment py-cuda13` instead.
4. Install R packages (Seurat, anndataR, SignatuR, scATOMIC/scITD/HiTME deps,
   and all benchmark method packages) via the chained setup task:
   ```bash
   pixi run setup
   ```


### **Workflow**

The analysis proceeds through four stages. See
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md#preprocessing-pipeline) for the
full pipeline call flow, file-role tables, and HPC folder layout.

- **Stage 1 — QC Filtering:** Per-dataset .rmd notebooks in
  `notebooks/QC_filtering/`.
- **Stage 2 — Preprocessing + Cell Type Annotation:** Stage raw data from NAS
  to HPC scratch, run dataset-specific preprocessing steps, then the
  standardized Python/Scanpy preprocessing pipeline and HPC-parallelized
  scATOMIC + HiTME annotation (see
  [ARCHITECTURE.md](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-src4_cell_type_annotation)).
- **Stage 3 — Benchmark Analysis:** Render `notebooks/benchmark_analysis.rmd`
  in RStudio (current local flow, to be superseded by the planned HPC benchmark
  pipeline — see TODO.md). R-native methods run in the notebook; Python methods
  (MrVI, PILOT, scPoli) run separately and exchange data via `.feather` files.
- **Stage 4 — Batch Effect Analysis:** Render
  `notebooks/batch_effect_analysis.rmd` in RStudio (under expansion — see
  TODO.md).

**HPC execution** (SLURM; see
[ARCHITECTURE.md](docs/ARCHITECTURE.md#hpc-folder-layout) for the folder layout):

```bash
./src/1_stage_data/1_stage_data.sh                            # stage raw data (NAS → scratch, login node)
sbatch src/2_dataset_specific_preprocessing/1_submit_hpc.sh   # dataset-specific preprocessing steps
sbatch src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh     # standardized preprocessing + sync to NAS
./src/4_cell_type_annotation/1_prepare_chunks.sh              # production: 5 samples/chunk (test: 1)
./src/4_cell_type_annotation/2_submit_hpc_array.sh            # scATOMIC + HiTME annotation array
./src/4_cell_type_annotation/3.1_submit_merge.sh <DS>         # merge annotations into every view + sync to NAS
```

- **Benchmark methods** (`src/5_run_benchmark_methods/run_python_sample_embedding_methods/`) — PLANNED, not yet implemented: no SLURM submit script exists yet; the folder contains only the method notebook `1.2_benchmark_methods_py.qmd` (to be converted to `.py`). See TODO.md.

#### **Expected Outputs**

- `.feather` files — cross-language distance matrices and embeddings produced
  by Python benchmark methods and consumed by R processors.
- Publication figures (MDS plots, PCA biplots, benchmark bar charts, separation
  metric heatmaps, transformation analysis panels) and per-method/per-dataset
  execution-time logs.


## Reference

If you use ECODA or this benchmark code in your research, please cite our
preprint:

**Cell type composition drives patient stratification in single-cell RNA-seq
cohorts.** Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*. doi:
[10.64898/2026.03.27.714811v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1)
