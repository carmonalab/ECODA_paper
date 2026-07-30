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


# Documentation
- README.md
- ARCHITECTURE.md
- AGENTS.md
- Todo list in TODO.md
- Should be kept up-to-date upon any changes


# datasets.json
This acts as ground truth for the datasets evaluated in this study. See datasets.json for most up-to-date list of datasets used and conditions.
Do not change this file without asking.

# Batch effect analysis dataset info
- Most datasets are monolithic h5ad files with a batch_col, e.g.:
    - Stephenson has batch effect by batch_col "Site" (both, in terms of gene expression (major) and cell type composition (minor, just one or a few monocyte subtypes))
- A "combined PBMC" dataset was created from multiple other available datasets (included for method benchmark analysis and or batch effect analysis):
    - Combined PBMC (Stephenson, GongSharma, Zhu) (see batch_effect_analysis.rmd, see also TODO.md)


# HPC
- Current repo lives on a local MacOS computer.
- If you need to test HPC cluster bash scripts:
    - The HPC cluster is only available if logged in to the UNIGE network (user might work from home (needs to connect to VPN) or from the office (has access to UNIGE network)).
        - If in the UNIGE network, you can log in with `ssh halterc@login1.bamboo.hpc.unige.ch` (user might need to enter a password).
- Heavy scripts are run on the HPC cluster, specifically located in these folders:
    - `src/preprocess`
    - `src/cell_type_annotation`
- These HPC pipelines are not finished yet and need to be updated to be finalized:
    - `src/preprocess`
    - Not implemented yet: HPC pipeline bash scripts for `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` (needs to be adapted, e.g. converted to .py and run on the cluster by calling from a bash script)
    - `src/cell_type_annotation` — HPC-parallelized scATOMIC + HiTME cell type annotation.
        - renv remnants removed; R environment fully managed by pixi (`pixi run Rscript`).
        - R code extracted from bash heredoc into standalone `2.1.1.1_process_chunk.R`.
        - `config_helper.R` moved from `DEPRECATED_LEGACY_CODE/` to project root (env-var based).
        - See [ARCHITECTURE.md](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b) for full pipeline documentation.
- `slurm_config.sh` is the HPC config file, used by all bash scripts, containing paths to the HPC cluster and other settings.
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/preprocess/1_submit_hpc_array.sh)
- results must be copied back to NAS after processing (typically implemented in `*_submit_hpc_array.sh` scripts upon completion of the HPC jobs)
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start


# Pipeline Overview

Four-stage end-to-end pipeline:

- **Stage 1 — QC Filtering** (`notebooks/QC_filtering/`): Per-dataset R Markdown notebooks (12 cohorts). Standard scRNA-seq QC: mitochondrial genes, gene/transcript count thresholds, doublet removal. Produces QC-filtered input for Stage 2.
- **Stage 2 — Preprocessing** (`src/preprocess/1.1.1_preprocess.py` + `src/cell_type_annotation/`): Python/Scanpy pipeline. Sample/gene name standardization, sample subsetting (driven by `datasets.json`), low cell-count removal, HVG selection (batch-aware and non-batch modes, multiple sizes), unsupervised clustering (multiple Leiden resolutions), Harmony integration. Cell type annotation via scATOMIC + HiTME (HPC-parallelized, borrowed from another project, not yet polished).
- **Stage 3 — Benchmark Analysis** (`notebooks/benchmark_analysis.rmd` + `src/utils/` + `src/benchmark/` + `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`): `run_analyses()` orchestrates three sub-pipelines: (3.1) benchmark methods (R-native: ECODA, Pseudobulk, MOFA, scITD, GloScope, GloProp, Avg PCA, Deconvolution; Python: MrVI, PILOT, scPoli — output `.feather` files ingested by R processors), (3.2) transformation analysis comparing ECODA transformations via `datrans()` parallel engine, (3.3) zero imputation analysis (4 strategies + multiLN/multiRepl).
- **Stage 4 — Batch Effect Analysis** (`notebooks/batch_effect_analysis.rmd`): Under expansion per reviewer comments. Methods: ECODA (batch-associated CT removal), Pseudobulk (DESeq2 + limma), MrVI (native batch handling), GloScope (Harmony space), PILOT-GM-VAE (Harmony space).
- **HPC wrappers** (`src/`): SLURM scripts for parallel per-dataset execution. `slurm_config.sh`, `preprocess/1_submit_hpc_array.sh` (stages data + submits array), `cell_type_annotation/` (array job pipeline). Planned: `submit_benchmark_methods_py.sh`.

## R Modules (`src/utils/` + `src/benchmark/`)

11 utility files loaded by `src/utils/load_all_functions.R`, plus 2 benchmark-specific files in `src/benchmark/`:

| Module | Location | Description |
|---|---|---|
| `imports.R` | `src/utils/` | Package loading — 47 packages via `load_my_packages()` |
| `constants.R` | `src/utils/` | Centralized label maps: method names, score labels, dataset display names |
| `helpers.R` | `src/utils/` | `apply_method_labels()` (label recoding), `merge_exec_times()` |
| `math_utils.R` | `src/utils/` | `clr()`, `impute_zeros()` (4 strategies), `cv()`, `calc_perc_df()`, `min_max()`, `zscore_transform()`, `global_quantile_norm_gaussian()` |
| `scoring_metrics.R` | `src/utils/` | `calc_sep_score()` (ANOSIM + Silhouette + Modularity + ARI + LISI), `calc_sil()`, `calc_modularity()` (KNN: sqrt(n), 3, 6, 9), `compute_KNN_from_dist()`, `compute_snn_graph()`, `clust_eval()`, `calc_lisi()` |
| `pseudobulk.R` | `src/utils/` | Pseudobulk extraction, DESeq2 normalization with batch correction, gene blacklist filtering |
| `hvcs.R` | `src/utils/` | Highly variable cell type selection and mean-variance plotting |
| `seurat_utils.R` | `src/utils/` | Seurat object creation, h5ad loading, metadata extraction (some functions deprecated, moved to Python) |
| `plotting.R` | `src/utils/` | PCA, MDS, PCA contribution plots |
| `benchmark_methods_r.R` | `src/benchmark/` | R-side method processors: `process_coda_fig()`, `process_pseudobulk_fig()`, `process_mofa_bulk_fig()`, `process_scitd_fig()`, `process_mrvi_fig()`, `process_pilot_fig()`, `process_scpoli_fig()`, `process_gloscope_fig()`, `process_gloprop_fig()`, `process_deconv_fig()`, `process_avg_pca_embedding_fig()`, `create_result_bundle()` |
| `benchmark_pipeline.R` | `src/benchmark/` | Orchestration: `run_analyses()`, `run_benchmark_analysis()`, `datrans()` (parallel transformation engine), `run_transformation_analysis()`, `run_zeroimp_analysis()`, `exec_time()` |

See `docs/ARCHITECTURE.md` for the full call flow, function-level documentation, and cross-language dependencies.


# Domain Terminology

- **ECODA** (Exploratory Compositional Data Analysis): Uses CLR-transformed cell-type proportions for cohort-level patient stratification in an unsupervised setting.
- **CLR** (Centered Log-Ratio): Transformation for compositional data: `log(x_i / geometric_mean(x))`. Requires zero-imputation beforehand. Implemented in `src/utils/math_utils.R:6`.
- **HVCs** (Highly Variable Cell Types): Cell types with highest variance across samples, selected for stratification. Implemented in `src/utils/hvcs.R`.
- **Zero imputation**: Four strategies for handling zero cell-type counts before CLR: `counts_zeros` (replace zeros with fixed count), `counts_all` (add fixed count to all), `percentage_zeros` (replace zeros with percentage of row total), `percentage_all` (add percentage to all). Implemented in `src/utils/math_utils.R:30`. Also uses `multiLN` and `multiRepl` from `zCompositions`.
- **Pseudobulk**: Aggregating single-cell counts per sample before DESeq2 normalization. Implemented in `src/utils/pseudobulk.R`.
- **Separation metrics**: Evaluate how well methods recover known biological groupings. All in `src/utils/scoring_metrics.R`:
  - **ANOSIM**: Analysis of Similarities (`calc_sep_score()`)
  - **ARI**: Adjusted Rand Index (`clust_eval()`)
  - **Silhouette**: `calc_sil()`
  - **Modularity**: `calc_modularity()` with multiple KNN variants (sqrt(n), 3, 6, 9)
  - **LISI**: Local Inverse Simpson's Index (`calc_lisi()`, :159)
- **Harmony integration**: Batch correction by integrating PCA embeddings across samples/batches. Computed in `src/preprocess/1.1.1_preprocess.py`.
- **`.feather` files**: Apache Arrow format for cross-language data exchange. Python methods in `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` output distance matrices/embeddings as `.feather`; R method processors in `src/benchmark/benchmark_methods_r.R` ingest them.