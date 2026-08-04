# ECODA_paper Architecture & Call Graph

## Overview
This is an **R-based bioinformatics pipeline** for benchmarking batch effect correction methods in single-cell RNA-seq data. The code is organized as a **multi-layered, top-down architecture** with clear separation between orchestration, method processing, evaluation, and utilities.


## Datasets and Analysis Definitions (Benchmark and/or Batch Effect Analysis)
```
datasets.json (read by src/utils/datasets_io.R::read_datasets_json() and
               src/datasets_io.py::read_datasets_json())
└── per-dataset entries
    ├── "Adams" → benchmark_analysis view
    │   ├── output_file: AdamsT_2020_32832599_benchmark_analysis_ECODAprocessed.h5ad
    │   ├── label_col, low_res_ct_col, hi_res_ct_col (from dataset-level columns)
    ├── ...
```
Each entry carries dataset-level metadata (columns, display_name, tissue,
file_names) plus one or more views. `file_names` holds the dataset's raw
input file(s) (string or list), independent of views; view-level
`input_file_name` / `output_file_name` declare raw input and preprocessed
output per view. Both readers return all dataset-level fields plus the
matching views. R keeps skipping view-less datasets (callers always pass a
view filter); the Python reader includes them (e.g. `Zhu`, a view-less raw
source for the CombinedPBMC dataset that is never preprocessed standalone).


## HPC Folder Layout

Repo lives at `~/ECODA_paper` on the HPC (only affects docs; `PROJECT_ROOT`
derives automatically via `slurm_config.sh`). Worker nodes only access local
scratch; the NAS is reachable from the login node only.

```
# HPC home (Bamboo cluster)
$HOME/ECODA_paper                       # repo clone (PROJECT_ROOT): logs/, aux/, .pixi/
$HOME/scratch/ECODA_paper               # HPC_SCRATCH_DIR
├── <DS_NAME>/data/                     # staged raw inputs per dataset (1_stage_data.sh)
├── <DS_NAME>/output/                   # preprocessed .h5ad per view, chunks/, annotations, annotated .h5ad
├── CombinedPBMC/data/                  # combine output + rds→h5ad cache (CombinedPBMC/output/ holds its preprocessed output)
└── chunks_manifest.txt                 # global chunk manifest (2_submit_hpc_array.sh)
$HOME/reference_atlases/sketched_200ct/ # HOME_REF_DIR (HiTME reference maps)

# NAS (carmona_smb; login node only)
DataCollections/Standardized_SingleCell_Datasets/   # NAS_SC_DIR — raw source datasets
DataCollections/reference_atlases/sketched_200ct/   # NAS_REF_DIR
Projects/ECODA_paper/                               # NAS_TARGET_DIR — results back-sync target
└── <DS_NAME>/output/                   # rsynced per-dataset from ${HPC_SCRATCH_DIR}/<DS_NAME>/output
benchmark/{embeddings,plots}/           # method .feathers + notebook plots (TODO: 5_run_benchmark_methods decision)
batch_effect_analysis/{embeddings,plots}/           # same, for batch effect analysis
```

| Path | Env var | What goes here |
|---|---|---|
| `${HPC_SCRATCH_DIR}/<DS_NAME>/data/` | `HPC_SCRATCH_DIR` | Staged raw inputs per dataset (`1_stage_data.sh` rsyncs from `NAS_SC_DIR`) |
| `${HPC_SCRATCH_DIR}/<DS_NAME>/output/` | `HPC_SCRATCH_DIR` | Preprocessed .h5ad per view (keys `X_pca_{view}_hvg{n}`, ...); during annotation additionally `chunks/chunk_*.txt`, `annotations_chunk_*.feather`, merged annotated .h5ad |
| `${HPC_SCRATCH_DIR}/CombinedPBMC/data/` | `HPC_SCRATCH_DIR` | Dual role: output of the CombinedPBMC combine step and rds→h5ad cache; input to the preprocess array (`CombinedPBMC/output/` holds its preprocessed output) |
| `${HPC_SCRATCH_DIR}/chunks_manifest.txt` | `HPC_SCRATCH_DIR` | Global chunk manifest (one `DS_NAME<TAB>chunk_path` line per chunk), rebuilt by `2_submit_hpc_array.sh` on every run |
| `${HOME_REF_DIR}` | `HOME_REF_DIR` | HiTME reference maps (staged from `NAS_REF_DIR` by `1_prepare_chunks.sh`) |
| `${PROJECT_ROOT}/logs`, `aux/`, `.pixi/` | `PROJECT_ROOT` | SLURM logs (`LOGS_DIR`), auxiliary files (gene maps), pixi env |
| `${NAS_SC_DIR}` | `NAS_SC_DIR` | Raw source datasets (`folder_name` in `datasets.json`) |
| `${NAS_REF_DIR}` | `NAS_REF_DIR` | Reference atlas source (HiTME) |
| `${NAS_TARGET_DIR}/<DS_NAME>/output/` | `NAS_TARGET_DIR` | Results back-sync target: rsynced per-dataset from `${HPC_SCRATCH_DIR}/<DS_NAME>/output`; `benchmark/{embeddings,plots}/`, `batch_effect_analysis/{embeddings,plots}/` (targets for method `.feather`s + notebook plots; filled once the `5_run_benchmark_methods` decision is made — TODO) |


## Preprocessing Pipeline

### Standard scRNA-seq preprocessing pipeline (`src/3_scrnaseq_preprocessing/`)

The preprocessing stage is split across three `src/` folders run in sequence:

| Folder | Role |
|---|---|
| `src/1_stage_data/` | `1_stage_data.sh` — login-node script: sources `slurm_config.sh`, stages raw inputs from NAS to HPC scratch (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`) via jq over `datasets.json` + rsync. |
| `src/2_dataset_specific_preprocessing/` | `1_submit_hpc.sh` dispatcher submits all per-step sbatch jobs in this folder in parallel (`1.1_submit_combinedpbmc.sh`, `1.2_submit_joanito_batch_col.sh`), passing per-step `--output/--error` into `${LOGS_DIR}` (SLURM directives cannot expand env vars, so the flags are given on the sbatch command line), waits for all, reports per-job state via `sacct` and exits non-zero on any failure. Steps must be mutually independent. Run after staging, before the preprocess array. |
| `src/3_scrnaseq_preprocessing/` | `1_submit_hpc_array.sh` (array submit + monitor + rsync back to NAS), `1.1_run_worker.sh`, `1.1.1_preprocess.py`. |

#### Files # TODO

| File | Role |
|---|---|
| `1_submit_hpc_array.sh` | Thin bash wrapper: sources `slurm_config.sh`, submits the preprocess array (`1.1_run_worker.sh`), monitors completion, verifies via `sacct` that every task state is `COMPLETED` (fail-closed: aborts without syncing on any non-COMPLETED state or empty sacct output), rsyncs results back to NAS (login node). Raw-data staging no longer lives here — see `src/1_stage_data/1_stage_data.sh`. |
| `1.1_run_worker.sh` | `#SBATCH` worker: sources `slurm_config.sh`, resolves its dataset from `SLURM_ARRAY_TASK_ID`, calls `1.1.1_preprocess.py` (via `${PYTHON_BIN}`) with `--config_path/--input_dir/--output_dir/--ds_name`. |
| `1.1.1_preprocess.py` | Standardized preprocessing: filtering, gene/sample name standardization, sample subsetting (from `datasets.json`; subsetting runs on original values BEFORE the sample column is standardized, and errors out if the subset is empty), batch-aware HVG selection, PCA, Harmony integration, Leiden clustering. Writes one .h5ad per view. **CSR-on-disk by construction**: `base_preprocessing()` forces `tocsr()` on `X` and `layers["counts"]` unconditionally (not only dense inputs); scanpy ops afterwards (normalize_total, log1p, HVG selection, subsetting) preserve CSR and `write_h5ad()` preserves the in-memory format, so the written files are always CSR. |

- **View output keys** (unified per-view pipeline, `1.1.1_preprocess.py`): both view types get the same treatment; only the HVG `batch_key` differs (benchmark views: standardized sample column from `SAMPLE_COLNAME`; batch-effect views: the dataset's `batch_col`).
  - `X_pca_{view}_hvg{n}` — stored for **every** HVG size (`benchmark_analysis`: 3000/2000/1000; `batch_effect_analysis`: 2000).
  - Harmony + unsupervised clustering (neighbors + Leiden) run **only at the 2000-HVG pass** (`CLUSTER_N_HVG`), on both embeddings:
    - `X_pca_harmony_{view}_hvg2000`, `neighbors_{view}_hvg2000`, `leiden_res_{r}_{view}_hvg2000`
    - `neighbors_{view}_hvg2000_harmony`, `leiden_res_{r}_{view}_hvg2000_harmony` (on the harmony-corrected embedding)
  - Example (benchmark view): `X_pca_benchmark_analysis_hvg3000`, `X_pca_benchmark_analysis_hvg2000`, `X_pca_benchmark_analysis_hvg1000`, plus `X_pca_harmony_benchmark_analysis_hvg2000`, `leiden_res_0.1_benchmark_analysis_hvg2000`, `leiden_res_0.1_benchmark_analysis_hvg2000_harmony`, ... (batch views analogously with `batch_effect_analysis_hvg2000`).

- **NAS ↔ Scratch data flow**: Raw-data staging from NAS to scratch happens in `src/1_stage_data/1_stage_data.sh` (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`); Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS (per-dataset `${NAS_TARGET_DIR}/${DS_NAME}/output`).
- **Working-directory convention**: every HPC bash script sources `src/slurm_config.sh` and then `cd "${PROJECT_ROOT}"` (all existing scripts do this; keep it for any new script). `preprocess_utils.py` pins the embedded R working directory to `${PROJECT_ROOT}` at import time (its rpy2 interop sources `src/utils/load_all_functions.R` + 11 module files via repo-relative paths), so Python callers are CWD-independent — the `cd` remains belt-and-braces and required by convention.
- **Dataset-specific preprocessing steps** (all submitted via the `1_submit_hpc.sh` dispatcher, in parallel):
  - **CombinedPBMC combine** (`1.1_submit_combinedpbmc.sh` → `1.1.1_create_combinedpbmc_dataset.py`): HPC-capable. With `HPC_SCRATCH_DIR` set it defaults to `--layout per-dataset` (sources at `${HPC_SCRATCH_DIR}/<ds>/data`, output `${HPC_SCRATCH_DIR}/CombinedPBMC/data`); locally it defaults to the flat `PROJECT_ROOT/data` layout. 64G baseline — GongSharma is huge and may need more. Script is CWD-independent; still submitted via the `1_submit_hpc.sh` dispatcher.
  - **Joanito batch column** (`1.2_submit_joanito_batch_col.sh` → `1.2.1_create_joanito_batch_col.R`): computes the `seqtec` batch column in place (idempotent) from `${HPC_SCRATCH_DIR}/Joanito/data/JoaI_2022_35773407_Nofilt_whole.rds`. 32G mem baseline — the whole .rds is read and re-saved in a single process. Required before the preprocess array (the `batch_effect_analysis` view uses it as `batch_col`).
  - Both must run **after** `1_stage_data.sh` and **before** `1_submit_hpc_array.sh` (staging skips `folder_name: null` datasets, and the preprocess array task reads the combined file). Processing scripts follow the decimal depth convention (`N.N.N_<action>.<ext>`, mirroring `1.1.1_preprocess.py`/`2.1.1_process_chunk.R`).


### Cell Type Annotation Pipeline (`src/4_cell_type_annotation/`)

Cell type annotation runs as a **separate HPC-parallelized pipeline** on SLURM, independent of the benchmark analysis call graph above. It takes monolithic `.h5ad` files from preprocessing and annotates them cell-by-cell using **scATOMIC + HiTME** in parallel array jobs.

#### Workflow

```
[ Monolithic h5ad Files on NAS ]
               │
                ▼ (1_prepare_chunks.sh → 1.1_prepare_chunks.py)
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
| `1_prepare_chunks.sh` | Thin bash wrapper: sources `slurm_config.sh`, stages ref maps from NAS → scratch (raw data is staged by `src/1_stage_data/1_stage_data.sh`), then iterates over all datasets in `datasets.json` (or a single dataset passed as 2nd positional arg) calling `1.1_prepare_chunks.py` (`${PYTHON_BIN}`) via `srun` (4G, 30min) per dataset. Datasets without preprocessed `.h5ad` input are skipped with a warning. Supports `test` mode (`./1_prepare_chunks.sh test` → 1 sample/chunk). |
| `1.1_prepare_chunks.py` | Native Python (anndata, no reticulate): reads each .h5ad in backed mode (skipping `*_raw.h5ad` rds→h5ad conversion caches left by `preprocess_utils.load_single_input()`), extracts unique sample IDs from `sample_col` (env var `SAMPLE_COLNAME`), groups them into chunks of 5 (or 1 in test mode), writes `chunk_N.txt` files (1st line = h5ad path, subsequent lines = sample IDs). In production mode only, deletes stale `annotations_chunk_*.feather` from the output dir AFTER all chunks were generated successfully, so reruns cannot merge leftover annotations from an earlier chunk numbering (test mode leaves production feathers untouched). |
| `2.0_create_scgate_db.R` | One-time download of the scGate model DB (`get_scGateDB(..., force_update=TRUE)`), persisted to `${SCGATE_DB_PATH}` (`aux/scGateDB.rds`). Run from `2_submit_hpc_array.sh` before the annotation array so workers do not download in parallel. Idempotent (exits early if the cache exists). |
| `2_submit_hpc_array.sh` | Creates the scGate model DB cache once via `srun` (`2.0_create_scgate_db.R` → `${SCGATE_DB_PATH}`; failure is non-fatal, workers download + persist themselves), builds a global chunk manifest (`${HPC_SCRATCH_DIR}/chunks_manifest.txt`, one tab-separated `DS_NAME<TAB>chunk_path` line per chunk across all datasets or a single dataset passed as positional arg), submits a single SLURM array job (`--array=1-TOTAL_CHUNKS`, `MAX_NUM_CHUNKS_PARALLEL` concurrency), monitors for completion, verifies via `sacct` that EVERY task state is `COMPLETED` (fail-closed: aborts without syncing on any non-COMPLETED state or empty sacct output), then syncs results back to NAS via rsync (per-dataset `${NAS_TARGET_DIR}/${DS_NAME}/output`; in single-dataset mode only that dataset is synced). |
| `2.1_run_worker.sh` | `#SBATCH` worker (shared-cpu, 16G, 2h). Sources `slurm_config.sh`, loads jq (module loads do not propagate from the submit script), reads its `DS_NAME`/`CHUNK_FILE` from the global manifest line matching `SLURM_ARRAY_TASK_ID` (`sed -n`), auto-exports per-dataset `TISSUE_TYPE`/`NORMAL_TISSUE` from `datasets.json`, checks the chunk file exists, then calls `2.1.1_process_chunk.R` directly via `${PIXI_RSCRIPT}`. |
| `2.1.1_process_chunk.R` | Core annotation logic: reads the chunk's h5ad file **in backed mode** (`read_h5ad(..., backed="r")`; `obs` is metadata-only — read once via `py_to_r(adata$obs)` and reused for every sample), warns (not stops) if the on-disk X format is not `csr` (`adata$X$format`; anndata only overrides selective row-indexing for backed CSR — CSC would fall back to a full in-memory read per subset, i.e. silent OOM), then iterates per sample, extracts sample-level Seurat objects from the raw counts layer (`get_seurat_obj_from_h5ad()`, `counts_layer="counts"` with `X` fallback + warning; sourced from `src/utils/seurat_utils.R`), runs **scATOMIC** (5 attempts with timeout) then **HiTME** (5 attempts with timeout). scGate models load from the shared `${SCGATE_DB_PATH}` cache (download + persist fallback). Writes per-chunk `.feather` named after the chunk file (`chunk_1.txt` → `annotations_chunk_1.feather`), so reruns are stable regardless of array task renumbering. Dual annotation: scATOMIC provides layer-1..6 predictions + confidence; HiTME provides layer1/2/3 UCell signatures + scGate/ProjecTILs refinement. Only annotation columns are kept (not full Seurat objects). |
| `3_merge_annotations.py` | Reads all `annotations_chunk_*.feather` files, joins them into the input `.h5ad`'s `obs` on a `(sample, barcode)` composite key (barcodes repeat across samples and views; duplicate keys are dropped with a warning; `SAMPLE_COLNAME` read from env, default "Sample"), keeps only the whitelisted annotation columns, writes annotated `.h5ad`. |
| `config_helper.R` | (Project root) Builds path config from env vars exported by `slurm_config.sh` (`DS_NAME`, `HPC_SCRATCH_DIR`). All dataset-specific paths are per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (`path_data`, `path_output`, `path_output_chunks`); `path_ref` is a home resource (`HOME_REF_DIR`). Annotation feathers go to `${HPC_SCRATCH_DIR}/${DS_NAME}/output` directly (`path_output`) — the old `samples/`, `ecoda/`, `plots/` dirs are no longer created. Missing env vars raise `stop()` errors (no silent fallbacks). Called by `2.1.1_process_chunk.R`. |

#### Key design details

- **scATOMIC + HiTME dual annotation**: scATOMIC provides hierarchical cell-type predictions (layer_1..6) with confidence scores; HiTME annotates using scGate models + ProjecTILs reference maps, producing layer1/2/3 labels. Both are run on each sample independently (i.e. independent from the rest of the dataset).
- **Retry loops**: Both annotation methods have up to 5 retry attempts with dynamic timeouts (max(60s, n_cells/10000 × 600s)) to handle HPC node variability.
- **NAS ↔ Scratch data flow**: Raw-data staging from NAS to scratch happens only in `src/1_stage_data/1_stage_data.sh` (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`); cell type annotation consumes the preprocessed output of that pipeline (`${HPC_SCRATCH_DIR}/${DS_NAME}/output` per dataset, matching `config_helper.R`). `1_prepare_chunks.sh` stages reference maps; `2_submit_hpc_array.sh` creates the scGate model DB cache (gene standardization moved into `1.1.1_preprocess.py`; the `GENE_REF_FILE` staging block was removed). The gene reference is committed to the repo at `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` (Ensembl 105, GRCh38.p13), originally downloaded 14.02.2022 from https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; it is consumed by `src/gene_utils.py`. Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS.
- **Environment propagation**: `slurm_config.sh` `export`s all core vars (`PROJECT_ROOT`, `DATASETS_JSON_FILE`, `HPC_SCRATCH_DIR`, `SAMPLE_COLNAME`, ref/gene files, `PYTHON_BIN`, `PIXI_RSCRIPT`, `SCGATE_DB_PATH`, ...) so they reach R via `Sys.getenv()` and Python via `os.environ` through both `srun` (`1_prepare_chunks.sh`, `2_submit_hpc_array.sh`) and `sbatch` (`2_submit_hpc_array.sh`). Bash arrays do not propagate through `sbatch`, so workers derive `DS_NAME` from `datasets.json` (via jq, loaded on the worker — module loads do not propagate either); `2.1_run_worker.sh` auto-exports per-task `TISSUE_TYPE`/`NORMAL_TISSUE` from `datasets.json`.
- **Counts input**: scATOMIC/HiTME receive the raw counts from `layers["counts"]` (vaulted by `base_preprocessing`), not the log-normalized `X`; if the layer is missing, `X` is used with a warning.
- **Backed per-sample reads are selective**: preprocessed `.h5ad` files are CSR-on-disk for both `X` and `layers["counts"]` (see `1.1.1_preprocess.py` above), so `get_seurat_obj_from_h5ad()`'s per-sample row subset (`adata[cell_indices]`) reads only the selected rows' segments (`backed_csr_matrix` selective indexing); `obs` is metadata-only and never triggers matrix I/O. On a CSC-on-disk file the same subset would materialize the full matrix in memory per sample (anndata has no selective row-indexing override for CSC) — `2.1.1_process_chunk.R` warns on non-CSR input.
- **Output format**: Per-chunk `.feather` files (Apache Arrow, cross-language), named after the chunk file (`annotations_chunk_<N>.feather`) → merged into original `.h5ad` by `3_merge_annotations.py` on a `(sample, barcode)` composite key.

#### Usage

```bash
# 1. Prepare chunks (stages data + generates chunk files)
./1_prepare_chunks.sh                     # production, all datasets: 5 samples/chunk
./1_prepare_chunks.sh test                # test, all datasets: 1 sample/chunk
./1_prepare_chunks.sh production Stephenson   # production, single dataset

# 2. Submit SLURM array (auto-creates scGate DB cache if missing, monitors + syncs to NAS after completion)
./2_submit_hpc_array.sh                   # all datasets with chunks
./2_submit_hpc_array.sh Stephenson        # single dataset

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