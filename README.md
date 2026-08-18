[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the analyses, benchmark, and figures from the paper:  
**"Cell type composition drives patient stratification in single-cell RNA-seq cohorts"**.

---

## Overview

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization of cellular heterogeneity, but summarizing single-cell data for cohort-level patient stratification remains a challenge. We benchmarked nine state-of-the-art sample representation methods—**MOFA+, scITD, GloScope, GloProp, MrVI, PILOT, PILOT-GM-VAE, QOT, and scPoli**—alongside baseline approaches (**Pseudobulk** and cell-type composition via **ECODA**) across standardized clinical cohorts (see [`datasets.json`](datasets.json)).

### Key Findings
- **Performance:** Centered log-ratio (CLR)-transformed cell-type proportions (**ECODA**) consistently match or outperform complex gene-expression-based embeddings in recovering known biological groupings in an unsupervised setting.
- **Efficiency:** ECODA produces patient embeddings in seconds with minimal computational overhead compared to deep-learning or optimal-transport methods.
- **Robustness:** Highly robust to technical batch effects and across various cell-type annotation granularities.
- **Interpretability:** Stratification is often driven by a small subset of **Highly Variable Cell Types (HVCs)**, providing direct mechanistic insights.

The standalone **scECODA** R package for cohort-level compositional analysis is available at [github.com/carmonalab/scECODA](https://github.com/carmonalab/scECODA).

---

## Installation & Environment Setup

This project uses [Pixi](https://pixi.sh) for reproducible multi-language (R & Python) environment management.

### 1. Install Pixi
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Set Up the Environment

- **Local macOS / Workstation (Lightweight analysis & figure generation):**
  ```bash
  pixi install && pixi run setup
  ```
  *Note:* Only the interactive analysis notebooks (`notebooks/benchmark_analysis.rmd` and `notebooks/batch_effect_analysis.rmd`) are designed to run locally on precomputed results.

- **HPC Cluster (SLURM worker build for full pipeline execution):**
  ```bash
  sbatch src/utils/bash/setup_env_sbatch.sh
  ```

---

## Workflow Overview

The reproducible pipeline is structured into four main stages. Detailed execution instructions, SLURM array commands, and folder architectures are documented in [**`docs/ARCHITECTURE.md`**](docs/ARCHITECTURE.md).

```
Stage 1: QC Filtering           -> notebooks/QC_filtering/
Stage 2: Preprocessing & Annot  -> src/1_stage_data/ -> src/2_dataset_specific_preprocessing/
                                   -> src/3_scrnaseq_preprocessing/ -> src/4_cell_type_annotation/
Stage 3: Benchmark Methods      -> src/5_run_benchmark_methods/ (Python & R arrays)
                                   -> notebooks/benchmark_analysis.rmd (Figures & Evaluation)
Stage 4: Batch Effect Analysis  -> notebooks/batch_effect_analysis.rmd
```

---

## Repository Structure

```
├── datasets.json               # Authoritative dataset metadata, columns, and views
├── docs/                       # Architecture documentation and hotfix notes
│   └── ARCHITECTURE.md         # Full pipeline architecture, HPC layout, and call flow
├── notebooks/                  # R Markdown analysis notebooks & dataset onboarding
│   ├── benchmark_analysis.rmd  # Main benchmark figures, separation metrics, and tables
│   ├── batch_effect_analysis.rmd # Batch effect quantification & mitigation
│   ├── QC_filtering/           # Per-dataset raw QC filtering notebooks
│   └── dataset_onboarding/     # Onboarding scripts for new cohorts
├── src/                        # Pipeline source code
│   ├── 1_stage_data/           # Raw data staging from storage to compute scratch
│   ├── 2_dataset_specific_preprocessing/ # Dataset-specific cohort adaptations
│   ├── 3_scrnaseq_preprocessing/         # Standardized Scanpy preprocessing pipeline
│   ├── 4_cell_type_annotation/           # Parallel scATOMIC + HiTME annotation
│   ├── 5_run_benchmark_methods/          # Python & R benchmark method arrays
│   ├── slurm_config.sh         # Centralized HPC environment and path configuration
│   └── utils/                  # Core math, scoring, plotting, and I/O utilities
└── aux/                        # Reference maps, gene blocklists, and gene annotation tables
```

---

## Citation

If you use ECODA or this benchmark suite in your research, please cite:

> **Cell type composition drives patient stratification in single-cell RNA-seq cohorts.**  
> Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*.  
> doi: [10.64898/2026.03.27.714811v1](https://doi.org/10.64898/2026.03.27.714811v1)

---

## Reference Data

- **Ensembl 105 (GRCh38.p13):** Human gene reference for name standardization (`aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`).
- **Carmona Lab Reference Atlases:** ScRNA-seq reference atlases for cell-type mapping (Garnica et al., 2024; Figshare DOI: [10.6084/m9.figshare.26310994](https://doi.org/10.6084/m9.figshare.26310994); repo: [`carmonalab/Reference_maps`](https://github.com/carmonalab/Reference_maps)).