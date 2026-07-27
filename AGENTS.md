# Paper/repo review and update strategy
This repo is about: ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts
Link to paper: https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full

> Kilo Code indexes code structure and function signatures automatically.
> AGENTS.md focuses on domain concepts, pipeline logic, and project
> conventions that indexing cannot infer.

## Major reviewer points to be addressed
- Extend batch effect analysis
    - More datasets have to be added
        - A requirement for this: Batch effect analysis pipeline has to be implemented in a more structured and standardized way
        - Preprocessing has to be harmonized to fit with the benchmark analysis
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are defined in the TODO.md
        - A draft strategy for specific pipeline and code implementation are defined in the TODO.md
- Extend benchmark analysis
    - More datasets have to be added
        - Benchmark analysis has mainly to be adapted to be run on the HPC cluster and to be cleaned up with minor adaptions
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are PILOT-GM-VAE (very similar to PILOT which is already implemented, trivial) and possibly PULSAR (needs to be tested if it can be run at all, as it is a foundation model, requiring substantial hardware only available on the cluster, including GPU. Input for PULSAR is Universal Cell Embedding (UCE)).
            - PILOT-GM-VAE can be added by agent
            - PULSAR needs to be tested for requirements


# datasets.json
This acts as ground truth for the datasets evaluated in this study. See datasets.json for most up-to-date list of datasets used and conditions.
Do not change this file without asking.

# Batch effect analysis dataset info
- Most datasets are monolithic h5ad files with a batch_col, e.g.:
    - Stephenson has batch effect by batch_col "Site" (both, in terms of gene expression (major) and cell type composition (minor, just one or a few monocyte subtypes))
- A "combined PBMC" dataset was created from multiple other available datasets (included for method benchmark analysis and or batch effect analysis):
    - Combined PBMC (Stephenson, GongSharma, Zhu) (see batch_effect_analysis.rmd, see also TODO.md)


# HPC
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/bash/copy_data_from_nas_to_hpc_scratch.sh)
- results must be copied back to NAS after processing
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start


# Pipeline Overview

Four-stage end-to-end pipeline:

- **Stage 1 — QC Filtering** (`./QC_filtering/`): Per-dataset R Markdown notebooks (12 cohorts). Standard scRNA-seq QC: mitochondrial genes, gene/transcript count thresholds, doublet removal. Produces QC-filtered input for Stage 2.
- **Stage 2 — Preprocessing** (`src/py/preprocess.py` + `src/bash/cell_type_annotation/`): Python/Scanpy pipeline. Sample/gene name standardization, sample subsetting (driven by `datasets.json`), low cell-count removal, HVG selection (batch-aware and non-batch modes, multiple sizes), unsupervised clustering (multiple Leiden resolutions), Harmony integration. Cell type annotation via scATOMIC + HiTME (HPC-parallelized, borrowed from another project, not yet polished).
- **Stage 3 — Benchmark Analysis** (`benchmark_analysis.rmd` + `src/R/` + `src/py/benchmark_methods_py.qmd`): `run_analyses()` orchestrates three sub-pipelines: (3.1) benchmark methods (R-native: ECODA, Pseudobulk, MOFA, scITD, GloScope, GloProp, Avg PCA, Deconvolution; Python: MrVI, PILOT, scPoli — output `.feather` files ingested by R processors), (3.2) transformation analysis comparing ECODA transformations via `datrans()` parallel engine, (3.3) zero imputation analysis (4 strategies + multiLN/multiRepl).
- **Stage 4 — Batch Effect Analysis** (`batch_effect_analysis.rmd`): Under expansion per reviewer comments. Methods: ECODA (batch-associated CT removal), Pseudobulk (DESeq2 + limma), MrVI (native batch handling), GloScope (Harmony space), PILOT-GM-VAE (Harmony space).
- **HPC wrappers** (`src/bash/`): SLURM scripts for parallel per-dataset execution. `config.env`, `copy_data_from_nas_to_hpc_scratch.sh`, `cell_type_annotation/` (array job pipeline). Planned: `submit_preprocess_py.sh`, `submit_benchmark_methods_py.sh`.

## R Modules (`src/R/`)

12 files loaded in dependency order by `src/R/load_all_functions.R`:

| Module | Description |
|---|---|
| `imports.R` | Package loading — 47 packages via `load_my_packages()` |
| `constants.R` | Centralized label maps: method names, score labels, dataset display names |
| `helpers.R` | `apply_method_labels()` (label recoding), `merge_exec_times()` |
| `math_utils.R` | `clr()`, `impute_zeros()` (4 strategies), `cv()`, `calc_perc_df()`, `min_max()`, `zscore_transform()`, `global_quantile_norm_gaussian()` |
| `scoring_metrics.R` | `calc_sep_score()` (ANOSIM + Silhouette + Modularity + ARI + LISI), `calc_sil()`, `calc_modularity()` (KNN: sqrt(n), 3, 6, 9), `compute_KNN_from_dist()`, `compute_snn_graph()`, `clust_eval()`, `calc_lisi()` |
| `pseudobulk.R` | Pseudobulk extraction, DESeq2 normalization with batch correction, gene blacklist filtering |
| `hvcs.R` | Highly variable cell type selection and mean-variance plotting |
| `seurat_utils.R` | Seurat object creation, h5ad loading, metadata extraction (some functions deprecated, moved to Python) |
| `plotting.R` | PCA, MDS, PCA contribution plots |
| `benchmark_methods_r.R` | R-side method processors: `process_coda_fig()`, `process_pseudobulk_fig()`, `process_mofa_bulk_fig()`, `process_scitd_fig()`, `process_mrvi_fig()`, `process_pilot_fig()`, `process_scpoli_fig()`, `process_gloscope_fig()`, `process_gloprop_fig()`, `process_deconv_fig()`, `process_avg_pca_embedding_fig()`, `create_result_bundle()` |
| `benchmark_pipeline.R` | Orchestration: `run_analyses()`, `run_benchmark_analysis()`, `datrans()` (parallel transformation engine), `run_transformation_analysis()`, `run_zeroimp_analysis()`, `exec_time()` |

See `docs/ARCHITECTURE.md` for the full call flow, function-level documentation, and cross-language dependencies.


# Domain Terminology

- **ECODA** (Exploratory Compositional Data Analysis): Uses CLR-transformed cell-type proportions for cohort-level patient stratification in an unsupervised setting.
- **CLR** (Centered Log-Ratio): Transformation for compositional data: `log(x_i / geometric_mean(x))`. Requires zero-imputation beforehand. Implemented in `src/R/math_utils.R:6`.
- **HVCs** (Highly Variable Cell Types): Cell types with highest variance across samples, selected for stratification. Implemented in `src/R/hvcs.R`.
- **Zero imputation**: Four strategies for handling zero cell-type counts before CLR: `counts_zeros` (replace zeros with fixed count), `counts_all` (add fixed count to all), `percentage_zeros` (replace zeros with percentage of row total), `percentage_all` (add percentage to all). Implemented in `src/R/math_utils.R:30`. Also uses `multiLN` and `multiRepl` from `zCompositions`.
- **Pseudobulk**: Aggregating single-cell counts per sample before DESeq2 normalization. Implemented in `src/R/pseudobulk.R`.
- **Separation metrics**: Evaluate how well methods recover known biological groupings. All in `src/R/scoring_metrics.R`:
  - **ANOSIM**: Analysis of Similarities (`calc_sep_score()`)
  - **ARI**: Adjusted Rand Index (`clust_eval()`)
  - **Silhouette**: `calc_sil()`
  - **Modularity**: `calc_modularity()` with multiple KNN variants (sqrt(n), 3, 6, 9)
  - **LISI**: Local Inverse Simpson's Index (`calc_lisi()`, :159)
- **Harmony integration**: Batch correction by integrating PCA embeddings across samples/batches. Computed in `src/py/preprocess.py`.
- **`.feather` files**: Apache Arrow format for cross-language data exchange. Python methods in `src/py/benchmark_methods_py.qmd` output distance matrices/embeddings as `.feather`; R method processors in `src/R/benchmark_methods_r.R` ingest them.