# Plan: General TODO — Simplest Tasks (lines 1–37)

## Scope

The "General" section of `TODO.md` (lines 1–37) contains several meta-tasks. This plan covers only the **simplest, most mechanical tasks** that require no external research, no paper reading, and no design deliberation.

## Identified Simplest Tasks

### Task 0: Rewrite `docs/ARCHITECTURE.md` to match current codebase and pipeline

**Problem:** ARCHITECTURE.md is heavily outdated. It was written for the preprint-era codebase (monolithic `functions.R`, `MAIN_Analysis.Rmd`, R-only preprocessing) and does not reflect the post-reviewer-refactoring (modular `src/R/`, Python preprocessing, extended batch effect analysis, HPC infrastructure required to manage all the heavy computing required for preprocessing and running the benchmark methods).

#### Historical context to capture

The project evolved in three phases:
1. **Preprint (submitted 2026-03-31)**: All code in `functions.R` + `MAIN_Analysis.Rmd`, preprocessing in R (`Preprocess_datasets.Rmd`).
2. **Cleanup and streamlining (2026-03-31 → present)**: Repo documentation, file renames (`MAIN_Analysis.Rmd` → `benchmark_analysis.rmd`), dataset metadata centralized to `datasets.json`. Full preprocessing migrated to Python (`src/py/preprocess.py`), R functions modularized (`src/R/*`), HPC workflow needs to added (`src/bash/`) to speed up pipeline to incorporate new datasets, methods and streamlined preprocessing pipeline to also be ready for batch effect analysis downstream.
3. **Reviewer feedback (2026-04-29 → present)**: Batch effect analysis needs to be expanded (`batch_effect_analysis.rmd`)

#### Specific inaccuracies to fix

| Current (Wrong) | Should Be |
|---|---|
| Title: "batch effect correction methods" | "sample representation benchmarking for patient stratification" |
| All references to `functions.R` | Correct modular file in `src/R/` |
| `MAIN_Analysis.Rmd` / `Preprocess_datasets.Rmd` as entry points | `benchmark_analysis.rmd`, `batch_effect_analysis.rmd` |
| Scoring: ANOSIM, ARI, Silhouette, Modularity | Add LISI score (`calc_lisi`, `scoring_metrics.R:159`) |
| `calc_modularity` single KNN | Multiple KNN variants: sqrt(n), 3, 6, 9 |
| No Python pipeline | Add `src/py/preprocess.py` and `benchmark_methods_py.qmd` |
| No batch effect analysis | Add `batch_effect_analysis.rmd` workflow |
| No HPC infrastructure | Add `src/bash/` SLURM scripts and cell type annotation pipeline |
| `datrans` as "transformation dispatcher" | Parallel analysis engine with `foreach`/`doParallel` (`benchmark_pipeline.R:5`) to test transformation of cell type compositional data (from counts and/or relative abundances in percent), as described in the literature on CoDA (compositional data analysis) |
| "42 named functions in `functions.R`" | Functions distributed across 12 modules |

#### Pipeline flow to document (4 stages)

The ARCHITECTURE.md should describe the end-to-end pipeline as a 4-stage flow:

**Stage 1 — QC Filtering** (`./QC_filtering/`)
- 12 per-dataset R Markdown notebooks (one per cohort (for most but not all datasets, some were taken from the authors as-is): Adams, Bassez, Kfoury, Kim, Lee, Pelka, Smillie, Stephenson, Wu, Zhang, Zhu, JoaI).
- Standard scRNA-seq QC: mitochondrial gene expression, gene count, transcript count thresholds, simple doublet removal based on upper threshold cutoff of number of transcripts per droplet (cell).
- Not changed since preprint; produces QC-filtered input files for Stage 2.

**Stage 2 — Preprocessing** (`src/py/preprocess.py` for general preprocessing and `src/bash/cell_type_annotation` for cell type annotation)
- Migrated from R (`Preprocess_datasets.Rmd`) to Python for performance and streamlining.
- Reads `.h5ad` or `.rds` files (previously → converted and saved to both formats (`.h5ad` for Python benchmark methods, `.rds` Seurat for R benchmark methods), now safes only .h5ad files but has a function to convert to Seurat object when needed for R benchmark methods).
- Standardization: sample names (no leading digits, hyphens → underscores), gene names via `bionty` (previously `STACAS` in R).
- Sample subsetting via `subset_vars` defined in `datasets.json`.
- Low cell-count sample removal.
- HVG selection (excluding blacklisted genes (TBD, possibly using blacklist from STACAS/SignatuR)): batch-aware (`batch_key`) and non-batch modes (batch_key="Sample" current common practice in the field), multiple HVG sizes (3000, 2000, 1000).
- Unsupervised clustering: ScaleData → RunPCA → FindNeighbors → FindClusters (multiple Leiden resolutions).
- Harmony integration for batch views (and separate test for unsupervised clustering-based cell type annotation in non-batch-effect datasets (for paper appendix)).
- Config driven by `datasets.json` (centralized dataset metadata, replaced `Processed_dataset_metadata.R` and various other re-definitions of dataset metadata across multiple scripts).
- Cell type annotation: scATOMIC + HiTME (parallelization issues → borrowed cell type annotation pipeline from another project, scripts in `src/bash/cell_type_annotation/` — not yet polished/adopted).

**Stage 3 — Benchmark Analysis** (`benchmark_analysis.rmd` + `src/R/` + `src/py/benchmark_methods_py.qmd`)
- Main entry: `benchmark_analysis.rmd` orchestrates `run_analyses()` → three sub-pipelines:
  - **3.1 `run_benchmark_analysis()`**: Runs all benchmarked sample embedding methods.
    - R-native methods: ECODA (CoDA), Pseudobulk (DESeq2), MOFA, scITD, GloScope, GloProp, Avg PCA Embedding, Deconvolution (EPIC). Implemented in `src/R/benchmark_methods_r.R`.
    - Python methods: MrVI (scvi-tools), PILOT (pilotpy), scPoli (scarches). Run in `src/py/benchmark_methods_py.qmd`, output Euclidean inter-sample distances as `*_method_dists.feather` or embeddings as `*_embs.feather`.
  - **3.2 `run_transformation_analysis()`**: Compares ECODA transformation methods (clr, alr, arcsine_sqrt, freq, counts).
  - **3.3 `run_zeroimp_analysis()`**: Compares zero imputation strategies (counts_zeros, counts_all, percentage_zeros, percentage_all, multLN, multRepl).
  - **3.4 Plotting/summarization**: `benchmark_analysis.rmd` generates heatmaps, MDS plots, performance comparisons.
- R module structure (`src/R/`):
  - `benchmark_pipeline.R` — orchestration (`run_analyses`, `run_benchmark_analysis`, `datrans`, `run_transformation_analysis`, `run_zeroimp_analysis`, `exec_time`)
  - `benchmark_methods_r.R` — R-side method processors (`process_coda_fig`, `process_pseudobulk_fig`, `process_mofa_bulk_fig`, `process_scitd_fig`, `process_mrvi_fig`, `process_pilot_fig`, `process_scpoli_fig`, `process_gloscope_fig`, `process_gloprop_fig`, `process_deconv_fig`, `process_avg_pca_embedding_fig`, `create_result_bundle`) and ingests output from `src/py/benchmark_methods_py.qmd` to incorporate the methods run in python
  - `scoring_metrics.R` — evaluation core (`calc_sep_score`, `calc_sil`, `calc_modularity`, `clust_eval`, `calc_lisi`, `compute_KNN_from_dist`, `compute_snn_graph`)
  - `math_utils.R` — CLR, zero imputation (4 strategies), min-max, z-score, quantile normalization, CV
  - `hvcs.R` — highly variable cell type selection and mean-variance plotting
  - `pseudobulk.R` — pseudobulk extraction, DESeq2 normalization with batch correction, gene blacklist filtering (needs to be moved to preprocessing HVG selection)
  - `seurat_utils.R` — Seurat object creation, h5ad loading, metadata/label extraction, sample name standardization (deprecated, moved to `src/py/preprocess.py`), multi-resolution clustering (deprecated legacy code), HiTME annotation replacement (deprecated legacy function, replace_HiTMElayer3_annot should be removed and thresh = 0.1 should be set directly during cell type annotation with HiTME in `src/bash/cell_type_annotation/2.2_process_chunk.sh`)
  - `plotting.R` — PCA, MDS, PCA contribution plots (legacy function)
  - `helpers.R` — method label recoding (`apply_method_labels`), exec time merging (crutch function, should be properly handled by harmonizing `src/py/benchmark_methods_py.qmd` and `src/R/benchmark_methods_r.R`)
  - `constants.R` — centralized label maps for easy adaptation in `benchmark_analysis.rmd` plots (method display names, score labels, dataset display names, GongSharma facet labels)
  - `imports.R` — package loading for R scripts (47 packages via `load_my_packages()`)
  - `load_all_functions.R` — dependency-ordered module loader

**Stage 4 — Batch Effect Analysis** (`batch_effect_analysis.rmd`)
- Initially rudimentary (preprint: only Joanito + Combined PBMC datasets).
- Needs to be expanded per reviewer comments: more datasets, more methods, batch-aware preprocessing.
- Methods: ECODA (remove batch-associated cell types (TBD)), Pseudobulk (DESeq2 + limma batch correction (TBD)), MrVI (native batch handling), GloScope (harmony integrated space (calculated in `src/py/preprocess.py`)), PILOT-GM-VAE (harmony integrated space (calculated in `src/py/preprocess.py`)).

**Stage 2 and 3 wrappers: SLURM submission scripts for parallel execution per dataset on HPC Infrastructure** (`src/bash/`)
- `config.env` — general environment configuration (file paths to NAS, repo subfolders, pixi, etc.).
- `copy_data_from_nas_to_hpc_scratch.sh` — data transfer (NAS → HPC local scratch).
- Stage 2-associated: `submit_preprocess_py.sh` — parallelized preprocessing pipeline (HPC array job submission of `src/py/preprocess.py` per dataset) (TBD)
- Stage 2-associated: `cell_type_annotation/` — parallelized cell type annotation pipeline (chunk preparation, HPC array job submission, worker/process scripts). Borrowed from another project, not yet adapted and needs to be polished.
- Stage 3-associated: `submit_benchmark_methods_py.sh` — parallelized benchmark analysis pipeline (HPC array job submission of `src/py/benchmark_methods_py.qmd` per dataset and per method) (TBD)

#### Steps

1. Rewrite title and overview paragraph with correct project description and historical context
2. Replace the layered architecture model with the 4-stage pipeline flow described above and bash wrappers for HPC job submission
3. Replace all `functions.R` references with correct `src/R/` module paths and function locations
4. Update scoring metrics section to include LISI and multiple KNN modularity variants
5. Rewrite call flow diagram to reflect the 4-stage pipeline and modular structure
6. Update dependency summary with correct file counts, function locations, and cross-language dependencies

### Task 1: Fix outdated file references in README.md

**Problem:** README.md references files that no longer exist under those names. This is an independent task — the stale references are verifiable by listing files on disk, without needing ARCHITECTURE.md.

| README Reference | Actual File | Status |
|---|---|---|
| `Process_data.ipynb` | Python methods in `src/py/benchmark_methods_py.qmd` + `src/py/preprocess.py` | Stale |
| `MAIN_Analysis.Rmd` | `benchmark_analysis.rmd` | Renamed |
| `functions.R` | `src/R/` (12 modular files) | Restructured |

**Also missing from README:**
- `batch_effect_analysis.rmd`
- `src/R/` modular structure
- `src/py/` Python modules
- `src/bash/` HPC scripts

**Steps:**
1. Replace `Process_data.ipynb` → `src/py/benchmark_methods_py.qmd` (Python benchmark methods) + `src/py/preprocess.py` (preprocessing)
2. Replace `MAIN_Analysis.Rmd` → `benchmark_analysis.rmd`
3. Replace `functions.R` → modular `src/R/` structure (key modules: `benchmark_pipeline.R`, `scoring_metrics.R`, `hvcs.R`, `pseudobulk.R`, etc.)
4. Add `batch_effect_analysis.rmd` to Repository Contents
5. Add brief mention of `src/bash/` for HPC workflows

### Task 2: Add concise workflow summary to AGENTS.md

**Problem:** AGENTS.md lacks any description of the project's data processing pipeline and analysis workflow (TODO.md line 9: "Add this in very concise form also to AGENTS.md").

**Dependency:** Requires Task 0 (updated ARCHITECTURE.md) for easier implementation (e.g. to source accurate information).

**Steps:**
1. Add a "Pipeline Overview" section to AGENTS.md with a concise summary derived from the rewritten `docs/ARCHITECTURE.md`, covering the 4-stage pipeline:
   - **Stage 1** — QC filtering: per-dataset notebooks in `./QC_filtering/`
   - **Stage 2** — Preprocessing: `src/py/preprocess.py` (Python) + `src/bash/cell_type_annotation/` (HPC cell type annotation, TBD). HVG selection, sample/gene standardization, unsupervised clustering, Harmony integration, config driven by `datasets.json`.
   - **Stage 3** — Benchmark Analysis: `benchmark_analysis.rmd` orchestrates `run_analyses()` → three sub-pipelines: benchmark methods (R-native + Python via `.feather` exchange), transformation analysis, zero imputation analysis. R functions modularized in `src/R/` (12 files).
   - **Stage 4** — Batch Effect Analysis: `batch_effect_analysis.rmd` (rudimentary, under expansion per reviewer comments).
   - **HPC wrappers**: `src/bash/` SLURM scripts for parallel execution per dataset.
2. List `src/R/` modules with 1-line descriptions (derived from Task 0)
3. Add reference to `docs/ARCHITECTURE.md` for full call graph and detailed function-level documentation

### Task 3: Add domain terminology glossary to AGENTS.md

**Dependency:** Requires Task 2 (updated AGENTS.md) to not collide with code changes in AGENTS.md.

**Problem:** TODO.md (line 15–16) asks to "note any domain-specific terminology that needs to be defined for the agent" and "highlight hard-to-understand mathematical transformations (e.g., clr, impute_zeros)."

**Steps:**
1. Add a "Domain Terminology" section to AGENTS.md defining:
   - **ECODA**: Exploratory Compositional Data Analysis — uses CLR-transformed cell-type proportions for cohort-level patient stratification
   - **CLR (Centered Log-Ratio)**: Transformation for compositional data; `log(x_i / geometric_mean(x))`; requires zero-imputation beforehand (see `src/R/math_utils.R:6`)
   - **HVCs (Highly Variable Cell Types)**: Cell types with highest variance across samples, selected for stratification (see `src/R/hvcs.R`)
   - **Zero imputation**: Four strategies for handling zero cell-type counts before CLR: `counts_zeros`, `counts_all`, `percentage_zeros`, `percentage_all` (see `src/R/math_utils.R:30`)
   - **Pseudobulk**: Aggregating single-cell counts per sample before DESeq2 normalization (see `src/R/pseudobulk.R`)
   - **Separation metrics**: ANOSIM, ARI, Silhouette, Modularity (multiple KNN), LISI — evaluate how well methods recover known biological groupings (see `src/R/scoring_metrics.R`)

## Tasks Explicitly EXCLUDED (not "simplest")

The following General-section tasks are excluded because they require paper reading, external research, design deliberation, or iterative human-in-the-loop work:

- **Adding install/run instructions to README.md** (line 4) — requires understanding the pixi/pixi.toml dependency setup and testing that instructions actually work in a fresh environment.
- **Adding data sources and preprocessing step descriptions to README.md** (line 6–7) — requires paper reading to accurately describe the biological context and dataset provenance.
- **Reading the bioRxiv paper** and extracting methodology/findings (line 10–12) — requires external research.
- **Researching AGENTS.md best practices online** (line 26–27) — requires external research.
- **Proposing folder structure reorganization** (line 23) — requires design deliberation and human approval.
- **Standardizing file naming conventions** (line 24) — requires design deliberation and `datasets.json` coordination.
- **Consolidating duplicate scripts** (line 35) — requires code analysis and refactoring.
- **Iterating on suggestions until complete** (line 37) — requires human review cycles.

## Execution Order

1. **Task 0** (ARCHITECTURE.md rewrite) — foundational prerequisite for Tasks 2 and 3 which source from it.
2. **Task 1** (README.md fixes) — independent; can run in parallel with Task 0 (stale references are verifiable by listing files on disk).
3. **Task 2** (AGENTS.md workflow) — depends on Task 0 (sources pipeline summary from rewritten ARCHITECTURE.md).
4. **Task 3** (AGENTS.md terminology) — depends on Task 2 to avoid edit collisions on the same file (AGENTS.md).

Tasks 0 and 1 can execute in parallel. Tasks 2 and 3 must run sequentially after Task 0.

## Validation

- After Task 0: Verify every file/function reference in ARCHITECTURE.md exists in the codebase; verify no references to `functions.R` or `MAIN_Analysis.Rmd` remain; verify all 4 pipeline stages are documented; verify Python data exchange (`.feather` files) is described; verify HPC infrastructure section exists
- After Task 1: Verify every file referenced in README.md exists on disk; verify no stale references (`Process_data.ipynb`, `MAIN_Analysis.Rmd`, `functions.R`) remain
- After Task 2: Verify AGENTS.md pipeline section matches the actual code flow in `src/R/`, `src/py/`, and `src/bash/`
- After Task 3: Verify all defined terms map to actual code locations (grep for `clr`, `hvcs`, `anosim`, `lisi`, `impute_zeros`)
