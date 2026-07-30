# ECODA_paper Architecture & Call Graph

## Overview
This is an **R-based bioinformatics pipeline** for benchmarking batch effect correction methods in single-cell RNA-seq data. The code is organized as a **multi-layered, top-down architecture** with clear separation between orchestration, method processing, evaluation, and utilities.


## Datasets and Analysis Definitions (Benchmark and/or Batch Effect Analysis)
```
datasets.json (read by src/utils/datasets_io.R::read_datasets_json())
└── per-view entries
    ├── "Adams" → benchmark_analysis view
    │   ├── output_file: AdamsT_2020_32832599_benchmark_analysis_ECODAprocessed.h5ad
    │   ├── label_col, low_res_ct_col, hi_res_ct_col (from dataset-level columns)
    ├── ...
```
Each view declares input_file_name (raw input) and output_file_name (preprocessed output).
Dataset-level metadata (columns, display_name, tissue) is shared across views.


## Preprocessing Pipeline

### Standard scRNA-seq preprocessing pipeline (`src/preprocess/`)



### Cell Type Annotation Pipeline (`src/cell_type_annotation/`)

Cell type annotation runs as a **separate HPC-parallelized pipeline** on SLURM, independent of the benchmark analysis call graph above. It takes monolithic `.h5ad` files from preprocessing and annotates them cell-by-cell using **scATOMIC + HiTME** in parallel array jobs.

#### Workflow

```
[ Monolithic h5ad Files on NAS ]
               │
                ▼ (1_prepare_chunks.sh + 1.1_prepare_chunks.r)
     [ chunk_1.txt ]  [ chunk_2.txt ]  ...  [ chunk_N.txt ]  (5 samples each)
               │               │
               ▼ (SLURM Array) ▼ (SLURM Array)
          [ Worker 1 ]     [ Worker 2 ]  ...              (2.1_run_worker.sh)
               │               │
               ▼               ▼
     [ chunk_1.feather ] ... [ chunk_N.feather ]           (annotations on scratch)
               │               │
               ▼  (3_merge_annotations.py)  ▼
     [ annotated .h5ad on NAS ]
```

#### Files

| File | Role |
|---|---|
| `1_prepare_chunks.sh` | Thin bash wrapper: sources `slurm_config.sh`, stages data + ref maps from NAS → scratch, calls `1.1_prepare_chunks.r` via \`srun\` (4G, 10min). Supports `test` mode (`./1_prepare_chunks.sh test` → 1 sample/chunk). |
| `1.1_prepare_chunks.r` | Reads each .h5ad in backed mode (reticulate + anndata), extracts unique sample IDs from `sample_col` (env var `SAMPLE_COLNAME`), groups them into chunks of 5 (or 1 in test mode), writes `chunk_N.txt` files (1st line = h5ad path, subsequent lines = sample IDs). |
| `2_submit_hpc_array.sh` | Reads chunk count from scratch, submits a SLURM array job (`--array=1-N`, `MAX_NUM_CHUNKS_PARALLEL` concurrency), monitors for completion, then syncs results back to NAS via rsync. |
| `2.1_run_worker.sh` | `#SBATCH` worker (shared-cpu, 16G, 2h). Sources `slurm_config.sh`, reads `CHUNK_FILE` from `SLURM_ARRAY_TASK_ID`, calls `2.1.1_process_chunk.sh`. |
| `2.1.1_process_chunk.sh` | Thin wrapper that calls `2.1.1.1_process_chunk.R` via `pixi run Rscript --vanilla`. |
| `2.1.1.1_process_chunk.R` | Core annotation logic: reads the chunk's h5ad file, iterates per sample, extracts sample-level Seurat objects, runs **scATOMIC** (5 attempts with timeout) then **HiTME** (5 attempts with timeout). Writes per-chunk `.feather` file with annotation columns. Dual annotation: scATOMIC provides layer-1..6 predictions + confidence; HiTME provides layer1/2/3 UCell signatures + scGate/ProjecTILs refinement. Only annotation columns are kept (not full Seurat objects). |
| `3_merge_annotations.py` | Reads all `annotations_chunk_*.feather` files, joins them into the input `.h5ad`'s `obs` by cell barcode, keeps only the whitelisted annotation columns, writes annotated `.h5ad`. |
| `config_helper.R` | (Project root) Builds path config from env vars exported by `slurm_config.sh`. Called by both `1.1_prepare_chunks.r` and `2.1.1.1_process_chunk.R`. |

#### Key design details

- **scATOMIC + HiTME dual annotation**: scATOMIC provides hierarchical cell-type predictions (layer_1..6) with confidence scores; HiTME annotates using scGate models + ProjecTILs reference maps, producing layer1/2/3 labels. Both are run on each sample independently (i.e. independent from the rest of the dataset).
- **Retry loops**: Both annotation methods have up to 5 retry attempts with dynamic timeouts (max(60s, n_cells/10000 × 600s)) to handle HPC node variability.
- **NAS ↔ Scratch data flow**: Login node copies data from NAS to scratch before array starts. Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS.
- **Output format**: Per-chunk `.feather` files (Apache Arrow, cross-language) → merged into original `.h5ad` by `3_merge_annotations.py`.

#### Usage

```bash
# 1. Prepare chunks (stages data + generates chunk files)
export DS_NAME="Stephenson"
./1_prepare_chunks.sh                     # production: 5 samples/chunk
./1_prepare_chunks.sh test                # test: 1 sample/chunk

# 2. Submit SLURM array (monitors + syncs to NAS after completion)
./2_submit_hpc_array.sh

# 3. (Optional) Merge annotations into .h5ad if not done automatically
python 3_merge_annotations.py <path/to/input.h5ad> <path/to/annot_dir>
```

#### Test mode

`1_prepare_chunks.sh test` sets chunk_size = 1 (vs 5 for production). This means each chunk contains only 1 sample, producing more but smaller array jobs. Useful for quick validation. In the future, this will be replaced by the Joanito 5-sample debug dataset (see TODO.md).


---

## Benchmark, ECODA Transformation and ECODA Zero Imputation Analyses

### Data Processing

#### Complete Call Flow (Simplified)

```
notebooks/benchmark_analysis.rmd
    │
    ▼
run_analyses(ds, seurat)
    │
    ├──► run_benchmark_analysis()
    │       │
    │       ├──► get_pb_deseq2() ──► DESeq2.normalize()
    │       ├──► process_coda_fig() ──► datrans() ──► clr()
    │       ├──► process_pseudobulk_fig()
    │       ├──► process_mofa_bulk_fig()
    │       ├──► process_scitd_fig()
    │       ├──► process_gloscope_fig()
    │       └──► ... (15+ process_* functions)
    │                │
    │                └──► create_result_bundle()
    │                    ├──► calc_sep_score()
    │                    │   ├──► calc_modularity() ──► compute_KNN_from_dist()
    │                    │   │                              compute_snn_graph()
    │                    │   ├──► calc_sil()
    │                    │   ├──► clust_eval()
    │                    │   └──► anosim()
    │                    └──► plot_mds() ──► calc_modularity(), clust_eval()
    │
    ├──► run_transformation_analysis()
    │       └──► datrans() ──► clr(), impute_zeros()
    │                          └──► create_result_bundle()
    │
    └──► run_zeroimp_analysis()
            └──► impute_zeros() ──► create_result_bundle()
```


##### Key Design Patterns

1. **Factory Pattern**: `create_result_bundle()` is called by every processing function, producing a standardized result structure with `{scores, feat_mat, dist_mat, labels}`.

2. **Strategy Pattern**: `datrans()` dispatches to different transformation strategies (CLR, log, etc.) based on a method parameter.

3. **Pipeline Pattern**: Data flows linearly: AnnData object → Benchmark method data processing (Py/R) → Feature Matrix → Distance Matrix → Scores.
Note:
- Some methods are run in python before resulting in a distance matrix that is then read in notebooks/benchmark_analysis.rmd
- R methods are directly run in notebooks/benchmark_analysis.rmd
- Since not all methods provide a feature matrix, some directly output a distance matrix.

4. **Caching/Skip Logic**: `if (!method_name %in% names(res_list))` checks prevent re-running already-computed methods.

---

#### LAYER 1: Entry Point / Orchestration
```
run_analyses(ds, seurat, ...)
├── run_benchmark_analysis(...)
├── run_transformation_analysis(ct_comps, labels)
└── run_zeroimp_analysis(ct_comps, labels)
```
**`run_analyses`** is the main entry point. For each dataset it dispatches three separate analysis pipelines:
- **Benchmark analysis**: Runs the benchmarking method to calculate separation scores
- **Transformation analysis**: Tests different data transformations for ECODA (like CLR)
- **Zero imputation analysis**: Tests different zero-handling strategies for ECODA

---

#### LAYER 2: Benchmark Analysis Dispatcher
```
run_benchmark_analysis(...)
├── get_pb_deseq2(...)           → Pseudobulk preprocessing
├── process_pseudobulk_fig(...)  → Pseudobulk
├── process_coda_fig(...)        → ECODA variants (using different cell type annotations: LR, HR, shuffled, HiTME, scATOMIC, Leiden)
├── process_gloprop_fig(...)     → GloProp
├── process_deconv_fig(...)      → EPIC deconvolution
├── process_mofa_bulk_fig(...)   → MOFA
├── process_scitd_fig(...)       → scITD
├── process_gloscope_fig(...)    → GloScope
├── process_scpoli_fig(...)      → scPoli
├── process_mrvi_fig(...)        → MrVI
├── process_pilot_fig(...)       → PILOT
└── process_avg_pca_embedding(...) → Average PCA embedding as used by the authors of MrVI for their "Pseudobulk baseline"
```
This is the **largest function** (~500+ lines). It iteratively runs **dozens of method variants** by looping over parameters (HVG counts, factor numbers, resolution values, PCA dimensions). Each method is wrapped in `exec_time()` to track execution time.

---

#### LAYER 3: Method Processing Functions
Each `process_*_fig` function follows the same pattern:
```
process_*_fig → compute feature matrix → create_result_bundle from distance matrix → return scores
```

---

#### LAYER 4: Evaluation & Separation Scoring Core
```
create_result_bundle(feat_mat, labels, dist_mat)
└── calc_sep_score(dist_mat, labels)
    ├── calc_modularity(dist_mat, labels, knn_k)
    │   ├── compute_KNN_from_dist(dist_mat, knn_k)
    │   └── compute_snn_graph(knn)
    ├── calc_sil(dist_mat, labels)
    ├── clust_eval(dist_mat, labels)
    └── vegan::anosim(dist_mat, labels)
```

**This is the heart of the evaluation:**

| Function | Purpose |
|---|---|
| `calc_sep_score` | Aggregates all separation metrics into one score bundle |
| `calc_modularity` | Graph-based modularity (builds KNN → SNN graph, uses igraph) |
| `calc_sil` | Silhouette width (cluster::silhouette) |
| `clust_eval` | ARI from hierarchical + PAM clustering vs ground truth |
| `compute_KNN_from_dist` | Build K-nearest-neighbor matrix from distance matrix |
| `compute_snn_graph` | Build Shared Nearest Neighbor weighted graph (sparse matrix multiplication) |

---

#### LAYER 5: Utility Functions
```
Data extraction:
├── get_metadata(seurat)
├── get_labels(seurat, label_col)
├── get_ct_comp_df_seurat(seurat, sample_col, ct_col)
├── get_pb(seurat, sample_col)
└── get_current_hvgs(seurat)

ECODA zero imputation, data transformation and Pseudobulk preprocessing:
├── clr(df)                        → Centered Log-Ratio transformation
├── impute_zeros(df, method, num)  → Zero imputation (4 strategies)
├── calc_perc_df(df)               → Row-wise percentages
├── datrans(feat_mat, method, ...) → Main transformation dispatcher
├── DESeq2.normalize(matrix, ...)  → DESeq2-like normalization
└── standardize_sample_names()     → Fix digit-prefixed names

Visualization:
├── plot_mds(dist_mat, labels, ...)     → MDS plot with scores in title
├── plot_pca(...)                        → PCA visualization
└── varmeanplot(...)                     → Variance-mean relationship

Feature selection (highly variable cell types (HVCs)):
├── get_ct_var(df)          → Cell type variance ranking
└── get_hvcs(df_var, ...)   → Select top HVCs by variance threshold

Benchmarking:
├── prepare_benchmark_data(...)    → Prepare heatmap data
├── merge_exec_times(...)          → Join results with timing
├── apply_method_labels(...)       → Recode method names for display
├── min_max(values)                → Min-max normalization to [0,1]
└── exec_time(fun)                 → Time a code block
```

### Data Analysis

Data analysis is performed in a single notebook, `notebooks/benchmark_analysis.rmd`.

---

## Batch Effect Analysis # TODO