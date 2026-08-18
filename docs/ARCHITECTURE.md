# ECODA Pipeline Architecture & System Design

This document serves as the single technical source of truth for the **ECODA** reproducible benchmark suite: pipeline stages, HPC execution layout, data flow, and module reference.

---

## 1. System Overview & Core Design Principles

The ECODA benchmarking pipeline is designed for **cohort-level exploratory analysis of scRNA-seq datasets**. It compares compositional representations against gene-expression embeddings and baselines across standardized cohorts in a strictly unsupervised setting.

```
+---------------------------------------------------------------------------------------------------+
|                                       RAW DATA (NAS Storage)                                      |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v (1_stage_data.sh)
+---------------------------------------------------------------------------------------------------+
|                               STAGE 2: PREPROCESSING & ANNOTATION                                 |
|                                                                                                   |
|  [2_dataset_specific_preprocessing] -> [3_scrnaseq_preprocessing] -> [4_cell_type_annotation]       |
|    - GongSharma sample cap               - Standardized HVG, PCA       - Chunked scATOMIC & HiTME |
|    - CombinedPBMC combine                - Harmony batch correction    - Backed CSR union reads   |
|    - Joanito seqtec & _debug             - Leiden clusterings          - Atomic obs metadata join |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+---------------------------------------------------------------------------------------------------+
|                         STAGE 3: BENCHMARK & TRANSFORMATION ARRAYS                                |
|                                                                                                   |
|  [Python Benchmark Array]           [R Benchmark Array]            [Transformation & Zero-Imp]    |
|   - MrVI, scPoli, PILOT-GM-VAE       - GloScope, MOFA, Pseudobulk   - 7 CLR transformations       |
|   - PILOT (EMD), QOT                 - scITD, Composition (ECODA)   - 6 Zero-imputation methods   |
|   -> Emits .feather dists/embs       -> Emits .rds result bundles   -> Emits .rds bundles         |
+---------------------------------------------------------------------------------------------------+
                                                  |
                                                  v
+---------------------------------------------------------------------------------------------------+
|                             EVALUATION & PUBLICATION NOTEBOOKS                                    |
|                                                                                                   |
|  notebooks/benchmark_analysis.rmd                      notebooks/batch_effect_analysis.rmd        |
|  - Ingests precomputed .rds & .feather results          - Evaluates bio vs batch separation       |
|  - Computes ANOSIM, ARI, Silhouette, Modularity, LISI   - Tests batch-correction robustness       |
|  - Renders Figure 2A, Supp Fig 2, Supp Fig 14-21       - Generates batch evaluation panels        |
+---------------------------------------------------------------------------------------------------+
```

### Key Technical Contracts
1. **Backed HDF5 / CSR on Disk:** All preprocessed `.h5ad` files are stored in Compressed Sparse Row (CSR) format. This allows backed-mode chunk workers (`read_h5ad(..., backed="r")`) to slice sample subsets efficiently without loading entire dense count matrices into RAM.
2. **Minimal Annotation Union:** To prevent eager HDF5 group loading overhead during cell-type annotation, multi-view datasets use a minimal union layout (`X` = raw counts CSR, `obs` + `var` only, no `layers`/`obsp`/`obsm`).
3. **Cross-Language Data Exchange (`.feather`):** Python methods serialize embeddings and pairwise distance matrices to Apache Arrow `.feather` files, which R consumers read with zero format drift.
4. **Atomic Writes & MD5 Checksums:** All R result bundles are written atomically (temp file + rename) and cataloged in `checksums.md5` before being synchronized to persistent NAS storage.

---

## 2. Dataset Metadata & Schema (`datasets.json`)

All evaluated cohorts and views are configured in [`datasets.json`](../datasets.json).

### Schema Structure
```json
{
  "<DS_NAME>": {
    "display_name": "Human-readable label",
    "file_names": "Raw source file(s) in NAS storage",
    "folder_name": "Source subfolder in NAS_SC_DIR",
    "tissue": "Tissue of origin",
    "normal_tissue": true,
    "use_for_benchmark": true,
    "use_for_batch_effect": false,
    "columns": {
      "sample": "Metadata column holding sample/patient ID",
      "label": "Ground truth biological condition column",
      "batch": "Technical batch variable column (if applicable)",
      "cell_type_low_res": "Broad cell type annotation column",
      "cell_type_high_res": "High-granularity cell type annotation column"
    },
    "meta_cols_keep": ["List", "of", "obs", "columns", "to", "retain"],
    "views": {
      "benchmark_analysis": {
        "input_file_name": "Input raw file",
        "output_file_name": "<DS>_benchmark_analysis_ECODAprocessed.h5ad",
        "subset_vars": { "column": { "values": ["exclude_val"], "op": "notin" } }
      },
      "batch_effect_analysis": {
        "input_file_name": "Input raw file",
        "output_file_name": "<DS>_batch_effect_analysis_ECODAprocessed.h5ad",
        "subset_vars": {}
      }
    }
  }
}
```

### The `_debug` Dataset
- Built by `src/2_dataset_specific_preprocessing/1.3.1_prepare_joanito.R` into `${HPC_SCRATCH_DIR}/_debug/data/`.
- A 5-sample subset (500 cells/sample) covering combinations of biological origin and technical batches.
- Used for fast end-to-end integration tests of stages 2, 3, and 4. Default execution scripts skip `_*` entries unless `--ds_name _debug` is explicitly specified.

---

## 3. HPC & Storage Directory Layout

Execution is split between compute scratch (fast local I/O on worker nodes) and the NAS archive (persistent storage accessible via login nodes).

```
# HPC Scratch ($HOME/scratch/ECODA_paper -> ${HPC_SCRATCH_DIR})
${HPC_SCRATCH_DIR}/
├── <DS_NAME>/data/               # Staged raw inputs per dataset (1_stage_data.sh)
├── <DS_NAME>/output/             # Preprocessed .h5ad, annotation chunks, and merged outputs
├── <DS_NAME>/annotation_union/   # Minimal counts union for chunked annotation (temporary)
├── benchmark/
│   ├── embeddings/               # Python & R method distance matrices (.feather)
│   ├── results/                  # Evaluated R result bundles (<ds>_<method>.rds)
│   ├── pseudobulks/              # Shared DESeq2 normalized pseudobulk matrices
│   ├── gloscope_dists/           # Cached GloScope distance matrices
│   └── checksums.md5             # MD5 verification sidecar
├── chunks_manifest_<pid>.txt     # Per-submission chunk manifests
├── _worker_retries/              # Worker transient retry counters
└── _benchmark_watchdog/          # SLURM watchdog status logs

# NAS Storage (carmona_smb / Shared Collections; login node only)
${NAS_SC_DIR}/                    # Standardized raw source datasets
${NAS_REF_DIR}/                   # Carmona Lab reference maps (HiTME)
${NAS_TARGET_DIR}/                # Synced project results
├── <DS_NAME>/output/             # Rsynced preprocessed and annotated .h5ad files
└── benchmark/                    # Rsynced embeddings, results, and execution times
```

### Centralized Environment Variables (`src/slurm_config.sh`)

| Variable | Default Path / Value | Description |
|---|---|---|
| `PROJECT_ROOT` | `$HOME/ECODA_paper` | Root directory of the git clone |
| `HPC_SCRATCH_DIR` | `$HOME/scratch/ECODA_paper` | High-speed scratch directory on compute cluster |
| `NAS_SC_DIR` | `/srv/.../Standardized_SingleCell_Datasets` | Cold storage path for raw datasets |
| `NAS_TARGET_DIR` | `/srv/.../Projects/ECODA_paper` | Target sync directory on persistent storage |
| `HOME_REF_DIR` | `$HOME/reference_atlases/sketched_200ct` | Local sketched reference maps for HiTME |
| `SCGATE_DB_PATH` | `${PROJECT_ROOT}/aux/scGateDB.rds` | Pre-downloaded scGate model database cache |
| `PYTHON_BIN` | `${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python` | Target Python binary on HPC |
| `PIXI_RSCRIPT` | `pixi run --as-is -e py-cuda13 Rscript --vanilla` | Target Rscript invocation command |
| `SAMPLE_COLNAME` | `"Sample"` | Standardized sample metadata column name |

---

## 4. Pipeline Stages & Execution

### Stage 1 — QC Filtering (`notebooks/QC_filtering/`)
Per-dataset R Markdown notebooks performing study-specific initial quality control (filtering low-quality cells, doublets, and uninformative genes) before dataset standardization.

---

### Stage 2 — Preprocessing & Cell Type Annotation

```
+------------------+      +--------------------------------+      +-------------------------------+      +--------------------------+
| 1_stage_data.sh  | ---> | 2_dataset_specific_preproc     | ---> | 3_scrnaseq_preprocessing      | ---> | 4_cell_type_annotation   |
| (NAS -> Scratch) |      | (GongSharma, Combined, Joanito)|      | (Standardized Scanpy array)   |      | (scATOMIC + HiTME array) |
+------------------+      +--------------------------------+      +-------------------------------+      +--------------------------+
```

#### 1. Data Staging (`src/1_stage_data/`)
- `1_stage_data.sh`: Login-node utility that queries `datasets.json` and rsyncs required raw input files from `${NAS_SC_DIR}` to `${HPC_SCRATCH_DIR}/<DS_NAME>/data/`. Supports `--ds_name <DS>`.

#### 2. Dataset-Specific Preprocessing (`src/2_dataset_specific_preprocessing/`)
- `1_submit_hpc.sh`: Orchestration dispatcher running cohort-specific adjustments in parallel:
  - `1.1_submit_gongsharma.sh` (`1.1.1_subset_gongsharma.py`): Caps samples at ≤5,000 cells to prevent worker OOM while preserving deterministic sampling.
  - `1.2_submit_combinedpbmc.sh` (`1.2.1_create_combinedpbmc_dataset.py`): Assembles the multi-study CombinedPBMC cohort (Stephenson + GongSharma + Zhu) with lazy backed chunk loading.
  - `1.3_submit_joanito.sh` (`1.3.1_prepare_joanito.R`): Annotates `seqtec` batch column and generates the `_debug` 5-sample testing dataset.
  - `1.4_submit_kfoury_lowres_ct.sh` (`1.4.1_create_kfoury_lowres_ct.R`): Derives the low-resolution cell type column (`cells_lowres`).

#### 3. Standardized Preprocessing (`src/3_scrnaseq_preprocessing/`)
- `1_submit_hpc_array.sh` -> `1.1_run_worker.sh` -> `1.1.1_preprocess.py`:
  - Standardizes gene names using Ensembl 105 (`aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`).
  - Standardizes sample columns and applies view-specific `subset_vars`.
  - Performs library normalization, $\log(1+p)$ scaling, and batch-aware HVG selection (1000, 2000, 3000 HVGs).
  - Computes PCA (`X_pca_{view}_hvg{n}`), Harmony batch correction (`X_pca_harmony_{view}_hvg2000`), and multi-resolution Leiden clusterings.
  - Writes outputs in **CSR format on disk**.

#### 4. Parallel Cell Type Annotation (`src/4_cell_type_annotation/`)
- `1_prepare_chunks.sh` (`1.1_prepare_chunks.py`): Builds a minimal union `.h5ad` (`X` = counts CSR) and splits samples into chunk manifests (2 samples/chunk).
- `2_submit_hpc_array.sh` (`2.1_run_worker.sh` -> `2.1.1_process_chunk.R`): SLURM array executing dual annotation cell-by-cell on each sample:
  - **scATOMIC:** Hierarchical classification (layer 1–6) with bounded EM optimization.
  - **HiTME:** Reference mapping via scGate and ProjecTILs sketched reference atlases.
  - Generates per-chunk `.feather` annotations with per-sample checkpointing.
- `3_submit_merge.sh` (`3.1_merge_annotations.py`): Atomically joins feather annotations into every preprocessed view `.h5ad` on `(Sample, Barcode)` composite keys and synchronizes final annotated files to NAS.

---

### Stage 3 — Benchmark & Method Analyses (`src/5_run_benchmark_methods/`)

```
                                  +---------------------------------------+
                                  | 1_submit_hpc_array.sh (Submitter)     |
                                  +---------------------------------------+
                                         /            |            \
                                        /             |             \
                                       v              v              v
                        +------------------+  +------------------+  +------------------+
                        | Python Methods   |  | R Methods        |  | Transformation & |
                        | (GPU / CPU)      |  | (CPU Class)      |  | Zero-Imputation  |
                        | MrVI, scPoli,    |  | GloScope, MOFA,  |  | 7 transforms,    |
                        | PILOT, QOT,      |  | Pseudobulk,      |  | 6 zero-imp       |
                        | PILOT-GM-VAE     |  | scITD,           |  | strategies       |
                        |                  |  | Composition      |  |                  |
                        +------------------+  +------------------+  +------------------+
                                        \             |             /
                                         \            |            /
                                          v           v           v
                                  +---------------------------------------+
                                  | watchdog_main.sh (OOM Auto-Escalation)|
                                  +---------------------------------------+
                                                      |
                                                      v
                                  +---------------------------------------+
                                  | NAS Sync & checksums.md5 Verification |
                                  +---------------------------------------+
```

#### 1. Python Benchmark Array (`run_python_sample_embedding_methods/`)
- `1.1.1_benchmark_methods_py.py`: Implements MrVI (sample latent space), scPoli (sample embeddings), PILOT (Wasserstein distance), PILOT-GM-VAE (Gaussian Mixture VAE), and QOT (Quantum Optimal Transport).
- Outputs distance matrices and embeddings as `.feather` files to `${HPC_SCRATCH_DIR}/benchmark/embeddings/`.

#### 2. R Benchmark Array (`run_r_sample_embedding_methods/`)
- `1.1.1_benchmark_methods_r.R`: Evaluates GloScope, MOFA+, Pseudobulk (DESeq2 normalized counts + Euclidean distance), scITD (tensor decomposition), and the **`composition`** suite:
  - **ECODA variants:** Standard CLR, Highly Variable Cell Types (`ECODA_HVCs_*`), and unsupervised Leiden resolution scans (`ECODA_seuratres_*`).
  - **Baselines:** GloProp, EPIC deconvolution, Avg_PCA_embedding, and raw cell-type frequencies.
- Evaluates separation metrics and saves `.rds` result bundles alongside `<ds>_metadata.rds` to `${HPC_SCRATCH_DIR}/benchmark/results/`.

#### 3. Transformation & Zero-Imputation Arrays (`run_transformation_zeroimp_analysis/`)
- `1.1.1_transformation_zeroimp.R`: Systematically benchmarks 7 compositional transformations (CLR, ILR, ALR, Log, Arcsin, Identity, Log-ratio) and 6 zero-handling strategies (`counts_zeros`, `counts_all`, `percentage_zeros`, `percentage_all`, `multiLN`, `multiRepl`).

#### 4. Benchmark Watchdog & Memory Auto-Escalation (`watchdog_main.sh`)
- Submitted alongside array jobs as a lightweight monitoring process on a compute node.
- Survives SSH client disconnects.
- Detects `OUT_OF_MEMORY` task failures via `sacct` and automatically resubmits failed tasks with doubled RAM allocations (`128G -> 256G -> 500G`, clamped to `BENCHMARK_MEM_MAX`).
- On successful completion of all tasks, merges execution-time and peak RSS logs into `execution_times.feather` and triggers the NAS sync.

#### 5. Benchmark Analysis Notebook (`notebooks/benchmark_analysis.rmd`)
- Loads precomputed `.rds` result bundles and `.feather` matrices via `load_hpc_benchmark_results()`.
- Does **not** load monolithic raw `.h5ad` files into memory, enabling rapid local knitting on macOS.
- Generates publication figures (Figure 2A, Supp Fig 2, Supp Figs 14–21) and summary rank heatmaps.

---

### Stage 4 — Batch Effect Analysis (`notebooks/batch_effect_analysis.rmd`)

Evaluates method resilience against technical batch confounding across datasets with defined batch keys (e.g., Stephenson center differences, Joanito sequencing technology, CombinedPBMC multi-cohort merge):
- Tests ECODA batch-associated cell-type removal (testing each cell type's association with batch covariates).
- Evaluates Pseudobulk with batch-only correction (`blind=FALSE, correct_batch=TRUE`).
- Benchmarks sample representations derived from uncorrected vs. Harmony-integrated spaces.

---

### Dataset Onboarding Architecture (`notebooks/dataset_onboarding/`)

Standardized protocol for onboarding new external cohorts:
1. **HPC Download Worker (`download_datasets_hpc.sh` / `run_download_worker.sh`):** Array-based resumable `curl` downloading directly into BeeGFS scratch, followed by MD5 verification and NAS sync.
2. **Sample-First Subsetting (`onboarding_utils.py`):** Fast diagnostic inspection reading ~10–20 stratified samples directly from backed `.h5ad` files into memory (<2 GB RAM) for immediate exploratory analysis.
3. **Count Integrity Check:** Identifies raw integer count layers (`layers["counts"]` vs `raw.X` vs `X`).
4. **Standalone LISI Scoring (`onboarding_metrics.R`):** Calculates cell-level biological and batch LISI on unintegrated PCA space.
5. **Onboarding Check Notebooks (`dataset_check_<Name>.qmd`):** Diagnostic Quarto notebooks for review before registering cohorts into `datasets.json`.

---

## 5. File & Module Reference Tables

### Core Utilities (`src/utils/`)

| Script / Module | Primary Functions & Purpose |
|---|---|
| `src/utils/math_utils.R` | `clr_transform()`, `zero_imputation()`, ALR/ILR transformations |
| `src/utils/hvcs.R` | `select_hvcs()` — variance-based selection of Highly Variable Cell Types |
| `src/utils/pseudobulk.R` | `get_pb_deseq2()`, `DESeq2.normalize()` — pseudobulk aggregation & normalization |
| `src/utils/scoring_metrics.R` | `calc_sep_score()` (ANOSIM), `clust_eval()` (ARI), `calc_sil()`, `calc_lisi()`, `calc_modularity()` |
| `src/utils/datasets_io.R` | `read_datasets_json()`, `get_dataset_view_info()` (R parser) |
| `src/utils/py/datasets_io.py` | `read_datasets_json()` (Python parser matching R semantics) |
| `src/utils/py/gene_utils.py` | `standardize_gene_symbols()` using Ensembl 105 reference dictionary |
| `src/utils/bash/worker_retry.sh` | Sourced by SLURM workers for automated self-requeue on transient I/O faults |
| `src/utils/bash/sync_status_email.sh` | Best-effort email notification utility with per-task duration reports |
| `src/utils/env_check.R` | Validates package attachment and namespace loads across Pixi environments |

---

## 6. Resilience & Fault-Tolerance Mechanisms

1. **Transient Fault Recovery:** Array workers source `worker_retry.sh`. On non-zero exit codes caused by transient BeeGFS cache misses, workers grep the task `.err` log for known signatures and self-requeue via `scontrol requeue` (capped at `WORKER_MAX_RETRIES=3`).
2. **Deterministic Locking:** Environment refreshes lock on `logs/env_refresh.lock` to prevent parallel array workers from racing on the shared R library directory.
3. **Fail-Closed Verification:** All submitter sync tails verify that `sacct` reports `COMPLETED` for 100% of tasks before initiating NAS synchronization. If any task failed, synchronization is aborted and an alert email is dispatched.
4. **Resume via `--sync-only <job-id>`:** If an interactive SSH session disconnects while array jobs are running, re-running the submitter with `--sync-only <job-id>` reconnects to the monitoring tail without resubmitting duplicate jobs.