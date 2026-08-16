[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the results and figures from the
paper: **"Cell type composition drives patient stratification in single-cell
RNA-seq cohorts"**.

## Overview

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization
of cellular heterogeneity, but summarizing this data for cohort-level analysis
remains a challenge. Using **14 scRNA-seq datasets** (11 benchmark +
batch-effect cohorts; see [datasets.json](datasets.json) for the authoritative
list), we benchmarked nine state-of-the-art sample representation
methods—**MOFA+, scITD, GloScope, GloProp, MrVI, PILOT, PILOT-GM-VAE, QOT,
scPoli**—plus two baselines (**Pseudobulk** and cell-type composition via
**ECODA**). The benchmark evaluates their ability to recover known biological
groupings in a fully unsupervised setting.

### Key Findings

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


## Usage

### Repository Contents

- `docs/ARCHITECTURE.md`: Full pipeline architecture, call flow, and module documentation.
- `datasets.json`: Centralized dataset metadata (sample/label columns, subsetting rules, batch information).
- `notebooks/`: benchmark analysis, batch effect analysis, and per-dataset QC filtering notebooks (`.rmd`).
- `src`: Source code for the pipeline stages — staging (`src/1_stage_data/`), dataset-specific preprocessing (`src/2_dataset_specific_preprocessing/`), standardized preprocessing (`src/3_scrnaseq_preprocessing/`), cell type annotation (`src/4_cell_type_annotation/`), benchmark methods (`src/5_run_benchmark_methods/`), utilities (`src/utils/`), and HPC SLURM config (`src/slurm_config.sh`).


### Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/carmonalab/ECODA_paper
    ```
2. Install [Pixi](https://pixi.sh) — the project's package and environment manager.
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    echo 'export PATH="${HOME}/.pixi/bin:${PATH}"' >> ~/.bashrc
    echo '[ -f ~/.bashrc ] && source ~/.bashrc' >> ~/.bash_profile
    ```
3. Setup your environment:
    ```bash
    cd ~/ECODA_paper
    ```
    **HPC (recommended — main focus).** Build the `py-cuda13` environment on a
    worker node via sbatch (16 cpus / 64 GB / 2 h defaults; overridable, e.g.
    `sbatch --cpus-per-task=32 --mem=128G ...`). R source packages compile in
    parallel (`MAKEFLAGS=-jN`), so it is far faster and more robust than a
    login-node install:
    ```bash
    sbatch src/utils/bash/setup_env_sbatch.sh
    ```
    Never submit it while other jobs are active or a manual `refresh_env.sh`
    run is in progress — both entry points serialize on the shared
    `logs/env_refresh.lock` and refuse to start otherwise. Watch the build
    with `tail -f logs/setup_env_<jobid>.out`; success prints
    "R library integrity check OK" then "All packages in src/utils/imports.R + worker subsets + env_check.R load OK".
    For small re-runs on the login node: `tmux new -s env-refresh` →
    `src/utils/bash/refresh_env.sh`. After pulling changes to
    `pixi.toml`/`pixi.lock` onto the HPC clone, re-sync the env before
    submitting jobs — concurrent `pixi run` re-syncs from parallel jobs can
    race (observed: `failed to remove directory ... os error 2`); a serial
    re-sync avoids this. `pixi install` does NOT detect corrupt package files;
    if a job reports a missing R package, repair with:
    `rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup`.

    **Local macOS (optional — lightweight notebooks only):**
    ```bash
    pixi install && pixi run setup
    ```
    The first `setup` also compiles the R packages (Seurat, anndataR, ...) from
    source and takes a while. Only `notebooks/benchmark_analysis.rmd` and
    `notebooks/batch_effect_analysis.rmd` (which run on precomputed
    distance-matrix results) are intended to run locally; the heavy pipeline
    stages are HPC-only.


### Workflow

The analysis proceeds through four stages. See
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md#preprocessing-pipeline) for the
full pipeline call flow, file-role tables, and HPC folder layout.

- **Stage 1 — QC Filtering:** Per-dataset .rmd notebooks in
  `notebooks/QC_filtering/`.
- **Stage 2 — Data Crunching:** HPC pipelines for:
  - Dataset-specific preprocessing
  - Standard scRNA-seq preprocessing with Python/Scanpy.
  - Cell type annotation with scATOMIC + HiTME.
  - Running benchmark methods (Python/R)
  - ECODA transformation and zero-imputation analyses
- **Stage 3 — Benchmark Analysis:** Load crunched data and summarize results 
  in `notebooks/benchmark_analysis.rmd`.
- **Stage 4 — Batch Effect Analysis:** In `notebooks/batch_effect_analysis.rmd`.

#### HPC execution
Uses SLURM; see [ARCHITECTURE.md](docs/ARCHITECTURE.md#hpc-folder-layout) for the folder layout.

- **Important notes**
  ```
  `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh` creates a small "_debug" dataset for testing and debugging.
  Most pipelines have additional parameters, e.g.:
  - --ds_name: to run only for a specific dataset
  - --force: to overwrite existing outputs (otherwise they are typically skipped and not re-processed to safe ressources)

  Example:
  `./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods mrvi`
  ```

Running the pipeline:

- **Staging** (`src/1_stage_data/`):
  ```bash
  ./src/1_stage_data/1_stage_data.sh                         # stage raw data from cold storage (e.g. NAS → HPC scratch)
  ```

- **Dataset-specific preprocessing** (`src/2_dataset_specific_preprocessing/`):
  ```bash
  ./src/2_dataset_specific_preprocessing/1_submit_hpc.sh     # dataset-specific preprocessing steps -> also creates a small "_debug" dataset for testing and debugging
  ```

- **scRNA-seq preprocessing** (`src/3_scrnaseq_preprocessing/`):
  ```bash
  ./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh       # standardized preprocessing + sync to NAS
  ```

- **Cell type annotation** (`src/4_cell_type_annotation/`):
  ```bash
  ./src/4_cell_type_annotation/1_prepare_chunks.sh           # prepares chunk files required for the next step. See docs/ARCHITECTURE.md for details.
  ./src/4_cell_type_annotation/2_submit_hpc_array.sh         # scATOMIC + HiTME annotation array
  ./src/4_cell_type_annotation/3_submit_merge.sh             # merge annotations back to h5ad + sync to NAS
  ```

- **Benchmark methods** (`src/5_run_benchmark_methods/`):
  - **Python methods** (`run_python_sample_embedding_methods/`): MrVI, scPoli,
    PILOT-GM-VAE (GPU), PILOT, QOT (CPU).
    ```bash
    ./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh
    ```
  - **R methods** (`run_r_sample_embedding_methods/`): GloScope, MOFA,
    Pseudobulk, scITD + `composition` (the ECODA_* family, GloProp, EPIC
    deconv, Avg_PCA_embedding, Freq_highres — obs-only worker, also emits
    `<ds>_metadata.rds`) on the CPU benchmark class (pinned like PILOT for
    runtime comparability); `mofa`/`pseudobulk`/`composition` auto-prepend
    the `prepare_pseudobulk` prep step.
    ```bash
    ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh
    ```
  - **Transformation + zero-imputation analyses** (`run_transformation_zeroimp_analysis/`):
    two arrays (`trans`, `zeroimp`).
    ```bash
    ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh
    ```
  All benchmark submitters monitor their arrays, verify every task via `sacct`
  (fail-closed), merge the per-task execution-time logs (incl. per-worker peak
  RAM) into `benchmark/execution_times.feather` and rsync results to
  `Projects/ECODA_paper/benchmark/` on the NAS; the notebook
  (`notebooks/benchmark_analysis.rmd`) loads the R result bundles via
  `load_hpc_benchmark_results()` (path: `path_nas_benchmark` in the setup
  chunk) and reads zero h5ad files (labels/stats from the `<ds>_metadata.rds`
  bundles).

#### Expected Outputs

- `.feather` files — cross-language distance matrices and embeddings produced
  by Python benchmark methods and consumed by R processors.
- `.rds` result bundles — per-method/per-combo benchmark results (GloScope,
  MOFA, Pseudobulk, scITD, composition incl. `<ds>_metadata.rds`) +
  transformation/zero-imputation results, computed
  on HPC and loaded by the notebook (`benchmark/{results,pseudobulks,
  gloscope_dists}/` on the NAS).
- Publication figures (MDS plots, PCA biplots, benchmark bar charts, separation
  metric heatmaps, transformation analysis panels) and per-method/per-dataset
  execution-time logs (`benchmark/execution_times.feather`).


## How to Cite

If you use ECODA or this benchmark code in your research, please cite our
preprint:

**Cell type composition drives patient stratification in single-cell RNA-seq
cohorts.** Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*. doi:
[10.64898/2026.03.27.714811v1](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1)


## Reference data

- Ensembl 105 human gene reference (GRCh38.p13) used for gene-name standardization in the preprocessing pipeline
    - Retrieved from: https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz
- Reference maps for cell type annotation were obtained from the Carmona Lab
  light reference atlases (Garnica et al., 2024)
    - Hosted on Figshare: [https://doi.org/10.6084/m9.figshare.26310994](https://doi.org/10.6084/m9.figshare.26310994).
    - Associated GitHub Code Repository: [carmonalab/Reference_maps](https://github.com/carmonalab/Reference_maps)
    - The pipeline auto-stages the 4 sketched maps (`sketched_CD8T_human_ref_v1.rds`, `sketched_CD4T_human_ref_v2.rds`, `sketched_DC_human_ref_v2.rds`, `sketched_MoMac_human_v1.rds`) into `HOME_REF_DIR` (`${HOME}/reference_atlases/sketched_200ct/`), from the NAS when available, otherwise downloaded from the Figshare article above; the local repo no longer carries a mirror.