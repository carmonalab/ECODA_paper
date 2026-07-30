# Plan: Document Cell Type Annotation Pipeline

## Background

The deprecated `src/cell_type_annotation/DEPRECATED_LEGACY_CODE/README.md` describes an older version of the cell type annotation pipeline. The pipeline has been cleaned up (renv → pixi, config.env → slurm_config.sh + config_helper.R, R code extracted from heredoc, merge step changed from `3_summarize_ecoda.R` to `3_merge_annotations.py`). The **workflow concept is still accurate** but the implementation details differ. Currently, no documentation file describes the cell type annotation pipeline at the architecture level.

## Goal

Add a dedicated **Cell Type Annotation Pipeline** section to `docs/ARCHITECTURE.md`, then update `README.md` and `AGENTS.md` with cross-references.

---

## Task 1 — Add Cell Type Annotation Section to `docs/ARCHITECTURE.md`

Insert after the "Architecture Layers" block and before "Complete Call Flow" (after line 124).

### 1a. Workflow diagram

Add a section titled `### Cell Type Annotation Pipeline (Stage 2b)` with this text diagram:

```
[ Monolithic h5ad Files on NAS ]
               │
               ▼ (1_prepare_chunks.sh + 1_prepare_chunks.r)
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

### 1b. File-by-file description

| File | Role |
|---|---|
| `1_prepare_chunks.sh` | Thin bash wrapper: sources `slurm_config.sh`, stages data + ref maps from NAS → scratch, calls `1_prepare_chunks.r` via `srun` (4G, 10min). Supports `test` mode (`./1_prepare_chunks.sh test` → 1 sample/chunk). |
| `1_prepare_chunks.r` | Reads each .h5ad in backed mode (reticulate + anndata), extracts unique sample IDs from `sample_col` (env var `SAMPLE_COLNAME`), groups them into chunks of 5 (or 1 in test mode), writes `chunk_N.txt` files (1st line = h5ad path, subsequent lines = sample IDs). |
| `2_submit_hpc_array.sh` | Reads chunk count from scratch, submits a SLURM array job (`--array=1-N`, `MAX_NUM_CHUNKS_PARALLEL` concurrency), monitors for completion, then syncs results back to NAS via rsync. |
| `2.1_run_worker.sh` | `#SBATCH` worker (shared-cpu, 16G, 2h). Sources `slurm_config.sh`, reads `CHUNK_FILE` from `SLURM_ARRAY_TASK_ID`, calls `2.2_process_chunk.sh`. |
| `2.2_process_chunk.sh` | Thin wrapper that calls `2.2_process_chunk.R` via `pixi run Rscript --vanilla`. |
| `2.2_process_chunk.R` | Core annotation logic: reads the chunk's h5ad file, iterates per sample, extracts sample-level Seurat objects, runs **scATOMIC** (5 attempts with timeout) then **HiTME** (5 attempts with timeout). Writes per-chunk `.feather` file with annotation columns. Dual annotation: scATOMIC provides layer-1..6 predictions + confidence; HiTME provides layer1/2/3 UCell signatures + scGate/ProjecTILs refinement. Only annotation columns are kept (not full Seurat objects). |
| `3_merge_annotations.py` | Reads all `annotations_chunk_*.feather` files, joins them into the input `.h5ad`'s `obs` by cell barcode, keeps only the whitelisted annotation columns, writes annotated `.h5ad`. |
| `config_helper.R` | (Project root) Builds path config from env vars exported by `slurm_config.sh`. Called by both `1_prepare_chunks.r` and `2.2_process_chunk.R`. |

### 1c. Key design details

- **scATOMIC + HiTME dual annotation**: scATOMIC provides hierarchical cell-type predictions (layer_1..6) with confidence scores; HiTME refines these using scGate models + ProjecTILs reference maps, producing layer1/2/3 labels. Both are run on every sample.
- **Retry loops**: Both annotation methods have up to 5 retry attempts with dynamic timeouts (max(60s, n_cells/10000 × 600s)) to handle HPC node variability.
- **NAS ↔ Scratch data flow**: Login node copies data from NAS to scratch before array starts. Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS.
- **Output format**: Per-chunk `.feather` files (Apache Arrow, cross-language) → merged into original `.h5ad` by `3_merge_annotations.py`. No RDS files are produced by the current pipeline.

### 1d. Usage

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

### 1e. Test mode

`1_prepare_chunks.sh test` sets chunk_size = 1 (vs 5 for production). This means each chunk contains only 1 sample, producing more but smaller array jobs. Useful for quick validation. In the future, this will be replaced by the Joanito 5-sample debug dataset (see TODO.md).

---

## Task 2 — Update `AGENTS.md`

Replace lines 58-63 which say `cell_type_annotation` is "not finished yet" with:

```
- `src/cell_type_annotation` — HPC-parallelized scATOMIC + HiTME cell type annotation.
  - renv remnants removed; R environment fully managed by pixi (`pixi run Rscript`).
  - R code extracted from bash heredoc into standalone `2.2_process_chunk.R`.
  - `config_helper.R` moved from `DEPRECATED_LEGACY_CODE/` to project root (env-var based).
  - See [ARCHITECTURE.md](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b) for full pipeline documentation.
```

---

## Task 3 — Update `README.md`

Expand the Stage 2 bullet (lines 79-85) to mention cell type annotation clearly:

```
- **Stage 2 — Preprocessing + Cell Type Annotation:**
  - **Preprocessing** (`src/preprocess/1.2_preprocess.py`): Standardized sample/gene name
    standardization, HVG selection, clustering, and Harmony integration.
  - **Cell Type Annotation** (`src/cell_type_annotation/`): HPC-parallelized scATOMIC + HiTME
    annotation via SLURM array jobs. See the [Architecture
    documentation](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b) for workflow and
    usage.
```

Update the HPC execution paragraph (lines 93-94) to be more explicit:

```
**HPC execution:** Submit SLURM array jobs for:
- **Cell type annotation** via `src/cell_type_annotation/` (see [Architecture](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b) for step-by-step usage):
  ```bash
  export DS_NAME="Stephenson"
  ./src/cell_type_annotation/1_prepare_chunks.sh
  ./src/cell_type_annotation/2_submit_hpc_array.sh
  ```
- **Preprocessing** via `src/preprocess/1_submit_hpc_array.sh` (stages data + submits array + syncs results):
  ```bash
  sbatch src/preprocess/1_submit_hpc_array.sh
  ```
```

---

## Files to Modify

1. `docs/ARCHITECTURE.md` — Insert new section (major content addition)
2. `AGENTS.md` — Update cell_type_annotation bullet (minor)
3. `README.md` — Expand Stage 2 and HPC execution paragraphs (minor)

## Not in Scope

- Creating a standalone README for `src/cell_type_annotation/` (the info belongs in ARCHITECTURE.md)
- Updating the deprecated `DEPRECATED_LEGACY_CODE/README.md` (leave as-is for reference)
- Any code changes
