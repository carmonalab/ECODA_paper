[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the results and figures from the
paper: **"Cell type composition drives patient stratification in single-cell
RNA-seq cohorts"**.

### **Overview** Single-cell RNA

sequencing (scRNA-seq) enables high-resolution characterization of cellular
heterogeneity, but summarizing this data for cohort-level analysis remains a
challenge. We benchmarked several state-of-the-art sample representation
methods—including deep generative models and factor decomposition—against a
simple baseline: **ECODA (Exploratory Compositional Data Analysis)**.

### **Key Findings**

-   **Performance:** Centered log-ratio (CLR)-transformed cell-type proportions
    (ECODA) consistently match or outperform more complex methods in recovering
    known biological groupings in an unsupervised setting.
-   **Efficiency:** ECODA requires orders of magnitude fewer computational
    resources and produces embeddings in seconds.
-   **Robustness:** The approach is highly robust to technical batch effects and
    various cell-type annotation strategies.
-   **Interpretability:** Biological stratification is often driven by a small
    subset of highly variable cell types (HVCs), providing direct mechanistic
    insights.

### **Repository Contents**

-   `datasets.json`: Centralized dataset metadata (sample/label columns, subsetting rules, batch information).
-   `QC_filtering/`: Per-dataset R Markdown notebooks for standard scRNA-seq QC.
-   `src/py/preprocess.py`: Standardized preprocessing pipeline (Python/Scanpy) — sample/gene name standardization, HVG selection, unsupervised clustering, Harmony integration.
-   `src/py/benchmark_methods_py.qmd`: Python benchmark methods (MrVI, PILOT, scPoli).
-   `src/R/`: Modular R functions (12 files) — benchmark pipeline orchestration, scoring metrics, pseudobulk processing, HVC selection, math utilities, plotting.
-   `benchmark_analysis.rmd`: Core analysis script orchestrating the benchmark pipeline and generating paper figures.
-   `batch_effect_analysis.rmd`: Batch effect analysis and correction evaluation.
-   `src/bash/`: SLURM submission scripts for HPC parallel execution (preprocessing, cell type annotation, benchmark methods).
-   `docs/ARCHITECTURE.md`: Full pipeline architecture, call flow, and module documentation.

The **scECODA** R package for scalable cohort-level analysis is available at
[github.com/carmonalab/scECODA](https://github.com/carmonalab/scECODA). ---

## Reference

If you use ECODA or this benchmark code in your research, please cite our
preprint:

**Cell type composition drives patient stratification in single-cell RNA-seq
cohorts.** Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*. doi:
[10.64898/2026.03.27.714811v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1)
