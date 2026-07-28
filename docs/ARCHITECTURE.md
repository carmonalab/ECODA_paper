# ECODA_paper Architecture & Call Graph

## Overview
This is an **R-based bioinformatics pipeline** for benchmarking batch effect correction methods in single-cell RNA-seq data. The code is organized as a **multi-layered, top-down architecture** with clear separation between orchestration, method processing, evaluation, and utilities.

---

## Architecture Layers (Top-Down)

### LAYER 1: Entry Points / Orchestration
```
run_analyses(ds, seurat, ...)
├── run_benchmark_analysis(...)
├── run_transformation_analysis(ct_comps, labels)
└── run_zeroimp_analysis(ct_comps, labels)
```
**`run_analyses`** is the main entry point. For each dataset it dispatches three parallel analysis pipelines:
- **Benchmark analysis**: Runs all batch correction method variants
- **Transformation analysis**: Tests different data transformations (like CLR)
- **Zero imputation analysis**: Tests different zero-handling strategies

---

### LAYER 2: Benchmark Analysis Dispatcher
```
run_benchmark_analysis(...)
├── get_pb_deseq2(...)          → Pseudobulk normalization
├── process_pseudobulk_fig(...)  → Multiple variants (HVG500, HVG2000, PCA dims)
├── process_coda_fig(...)        → ECODA variants (LR, HR, shuffled, HiTME, scATOMIC, Leiden)
├── process_gloprop_fig(...)     → GloProp method
├── process_deconv_fig(...)      → EPIC deconvolution
├── process_mofa_bulk_fig(...)   → MOFA (multiple factor counts)
├── process_scitd_fig(...)       → scITD (multiple factor counts)
├── process_gloscope_fig(...)    → GloScope (multiple PCA dims)
├── process_scpoli_fig(...)      → scPoli
├── process_mrvi_fig(...)        → MrVI
├── process_pilot_fig(...)       → PILOT
└── process_avg_pca_embedding(...) → PCA embedding baseline
```
This is the **largest function** (~500+ lines). It iteratively runs **dozens of method variants** by looping over parameters (HVG counts, factor numbers, resolution values, PCA dimensions). Each method is wrapped in `exec_time()` to track execution time.

---

### LAYER 3: Method Processing Functions (Leaf Workers)
Each `process_*_fig` function follows the same pattern:
```
process_*_fig → compute feature matrix → create_result_bundle → return scores
```

---

### LAYER 4: Evaluation & Scoring Core
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

### LAYER 5: Utility Functions
```
Data extraction:
├── get_metadata(seurat)
├── get_labels(seurat, label_col)
├── get_ct_comp_df_seurat(seurat, sample_col, ct_col)
├── get_pb(seurat, sample_col)
└── get_current_hvgs(seurat)

Data transformation:
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

Feature selection:
├── get_ct_var(df)          → Cell type variance ranking
└── get_hvcs(df_var, ...)   → Select top HVCs by variance threshold

Benchmarking:
├── prepare_benchmark_data(...)    → Prepare heatmap data
├── merge_exec_times(...)          → Join results with timing
├── apply_method_labels(...)       → Recode method names for display
├── min_max(values)                → Min-max normalization to [0,1]
└── exec_time(fun)                 → Time a code block
```

---

### LAYER 6: Configuration & Data
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

---

## Complete Call Flow (Simplified)

```
MAIN_Analysis.Rmd / Preprocess_datasets.Rmd
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

---

## Key Design Patterns

1. **Factory Pattern**: `create_result_bundle()` is called by every processing function, producing a standardized result structure with `{scores, feat_mat, dist_mat, labels}`.

2. **Strategy Pattern**: `datrans()` dispatches to different transformation strategies (CLR, log, etc.) based on a method parameter.

3. **Pipeline Pattern**: Data flows linearly: Raw Seurat Object → Pseudobulk → Feature Matrix → Distance Matrix → Scores.

4. **Centralized Label Mapping**: Method/score/dataset label mappings are defined once at the top of `functions.R` and reused via `apply_method_labels()`.

5. **Caching/Skip Logic**: `if (!method_name %in% names(res_list))` checks prevent re-running already-computed methods.

---

## Dependency Summary

- **42 named functions** in `functions.R`
- **`create_result_bundle`** has the most callers (12) — it's the central evaluation hub
- **`calc_sep_score`** has 5 callers
- **`process_coda_fig`** is the most frequently invoked method processor (called ~20+ times with different parameters)
- **No unit tests** exist for any function (all marked as "no covering tests found")
- External dependencies: `igraph`, `cluster`, `mclust`, `vegan`, `Seurat`, `Matrix`, `ggplot2`, `dplyr`, `tidyr`