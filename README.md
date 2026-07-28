[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the results and figures from the
paper: **"Cell type composition drives patient stratification in single-cell
RNA-seq cohorts"**.

### **Overview**

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

### **Repository Contents**

-   `datasets.json`: Centralized dataset metadata (sample/label columns, subsetting rules, batch information).
-   `notebooks/QC_filtering/`: Per-dataset R Markdown notebooks for standard scRNA-seq QC.
-   `src/preprocess/1.2_preprocess.py`: Standardized preprocessing pipeline (Python/Scanpy) — sample/gene name standardization, HVG selection, unsupervised clustering, Harmony integration.
-   `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`: Python benchmark methods (MrVI, PILOT, scPoli).
-   `src/utils/` + `src/benchmark/`: Modular R functions — benchmark pipeline orchestration, scoring metrics, pseudobulk processing, HVC selection, math utilities, plotting.
-   `notebooks/benchmark_analysis.rmd`: Core analysis script orchestrating the benchmark pipeline and generating paper figures.
-   `notebooks/batch_effect_analysis.rmd`: Batch effect analysis and correction evaluation.
-   `src/`: SLURM submission scripts for HPC parallel execution (preprocessing, cell type annotation, benchmark methods).
-   `docs/ARCHITECTURE.md`: Full pipeline architecture, call flow, and module documentation.

The **scECODA** R package for scalable cohort-level analysis is available at
[github.com/carmonalab/scECODA](https://github.com/carmonalab/scECODA).

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

### **Usage / Workflow**

The analysis proceeds through four stages:

- **Stage 1 — QC Filtering:** Open the per-dataset `.Rmd` notebooks in
  `notebooks/QC_filtering/` and render in RStudio.
- **Stage 2 — Preprocessing:** Standardized sample/gene name standardization,
  HVG selection, clustering, and Harmony integration:
  ```bash
  pixi run -e py-cpu python src/preprocess/1.2_preprocess.py
  ```
  Dataset metadata (sample columns, subsetting rules, batch info) is driven by
  `datasets.json`.
- **Stage 3 — Benchmark Analysis:** Render `notebooks/benchmark_analysis.rmd` in RStudio.
  Python benchmark methods (MrVI, PILOT, scPoli) are invoked automatically via
  rpy2. The R pipeline orchestrates all method processors, scoring metrics, and
  figure generation.
- **Stage 4 — Batch Effect Analysis:** Render `notebooks/batch_effect_analysis.rmd` in
  RStudio.

**HPC execution:** Submit SLURM array jobs for cell-type annotation via
`src/cell_type_annotation/`. Preprocessing via `src/preprocess/1_submit_hpc_array.sh` (stages data + submits array + syncs results):
```bash
sbatch src/preprocess/1_submit_hpc_array.sh
```

### **Expected Outputs**

- `.feather` files — cross-language distance matrices and embeddings produced
  by Python benchmark methods and consumed by R processors.
- Publication figures — MDS plots, PCA biplots, benchmark bar charts,
  separation metric heatmaps, and transformation analysis panels.
- Execution time logs documenting per-method and per-dataset runtime.

---

## Reference

If you use ECODA or this benchmark code in your research, please cite our
preprint:

**Cell type composition drives patient stratification in single-cell RNA-seq
cohorts.** Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*. doi:
[10.64898/2026.03.27.714811v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1)
