# Pipeline Overhaul Execution Plan

## Rename folders to represent sequence of pipelines to be run — DONE
- [x] src/1_stage_data (extracted staging code from `src/preprocess/1_submit_hpc_array.sh` into `src/1_stage_data/1_stage_data.sh`)
- [x] src/2_dataset_specific_preprocessing (moved `src/preprocess/_create_combinedpbmc_dataset.py` and `src/preprocess/_create_joanito_batch_col.R` there)
- [x] src/3_scrnaseq_preprocessing
- [x] src/4_cell_type_annotation
- [x] src/5_run_benchmark_methods
- [x] Update documentation (ARCHITECTURE.md, AGENTS.md, README.md, TODO.md)

## Open reviewer points (to be addressed) — DONE except noted
- [x] **Centralise all hardcoded paths in the pipeline to `slurm_config.sh`**: sweep ALL scripts for hardcoded paths (incl. `src/5_run_benchmark_methods/`, notebooks, `src/utils/`); preprocess + cell type annotation are done, `src/5_run_benchmark_methods/` decision pending (see `## TODO:`) — pipeline scripts done; notebooks stay `getwd()`-based (local analysis); `src/5_run_benchmark_methods/` decision still pending
- [x] **`SAMPLE_COLNAME`**: decided — keep the env var export and additionally read it in `1.1.1_preprocess.py` (`os.environ.get("SAMPLE_COLNAME", "Sample")`), so the constant lives only in `slurm_config.sh` for both R and Python
- [x] **Can `config_helper.R` be replaced?** → see `## Step 1` — decided: keep `config_helper.R` at project root (R cannot source `slurm_config.sh` natively; the env-var-based `get_pipeline_config()` is the R-facing equivalent)
- [x] **Remove fallback code**: `getwd()` fallbacks dropped in `1.1_prepare_chunks.r` / `2.1.1.1_process_chunk.R` (now `stop()` guards on `PROJECT_ROOT`); `home_ref_dir`/`gene_ref` fallbacks dropped in `config_helper.R` (now `stop()` guards); `GENE_REF_FILE`/`GENE_REF_URL` removed from `slurm_config.sh` and the gene-annotation staging block removed from `1_prepare_chunks.sh` (unused — gene standardization lives in `gene_utils.py` via `aux/EnsemblGenes...`)
- [x] **Reticulate chunk preparation removed**: 1.1_prepare_chunks.r replaced by native Python 1.1_prepare_chunks.py (anndata backed read, --test flag, same chunk format); 1_prepare_chunks.sh now calls the pixi python binary; docs updated.
- [x] **NAS back-sync target**: `NAS_TARGET_DIR` now `/Volumes/Shared/Projects/ECODA_paper` (local) = `/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper` (HPC); rsync steps use the var, no script edits needed
- [x] **Re-write `src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py`** → see `## TODO:` — done (HPC-capable per-dataset layout default when `HPC_SCRATCH_DIR` set, explicit CLI); per-step sbatch wrapper added: `1.1_submit_combinedpbmc.sh` (later superseded by the `1_submit_hpc.sh` dispatcher)
- [x] **Re-write `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`** → see `## TODO:` — done (argparse `--config_path`/`--base_path`/`--output_dir` honored, absolute-vs-relative path resolution)

### Suggested order / priority — DONE
1. `NAS_TARGET_DIR` → `Projects/ECODA_paper` (config-only, unblocks correct result syncing)
2. `SAMPLE_COLNAME` decision (env var vs hardcoded `"Sample"`)
3. Re-write `1.1.1_preprocess.py` (explicit CLI parsing + env-driven path resolution) — unblocks HPC preprocessing
4. Re-write `_create_combinedpbmc_dataset.py` (HPC-capable; run before the preprocess array)
5. Remove `getwd()` fallbacks (once env vars are guaranteed)
6. `config_helper.R` replacement decision + `## Step 1` simplification
7. Centralisation sweep of remaining hardcoded paths (incl. `src/5_run_benchmark_methods/` decision) + docs update


## TODO:
- **`src/5_run_benchmark_methods/` (later)**: decision made — run on the HPC cluster:
  - `run_python_sample_embedding_methods/` in parallel across datasets (SLURM array, per-method+dataset workers)
  - `benchmark_methods_r.R` + `benchmark_pipeline.R` run on HPC sequentially via a simple single-worker bash script (e.g. `5_run_benchmark_methods/1_run_worker.sh`, subfolder structure TBD)
  - Move this data-processing step out of `notebooks/benchmark_analysis.rmd`
  - NAS targets: `benchmark/{embeddings,plots}/` and `batch_effect_analysis/{embeddings,plots}/` (see `NAS_TARGET_DIR`, HPC folder layout in ARCHITECTURE.md)

## Completed (HPC layout + TODO verification plan)
- [x] **TODO 3e — unified per-view pipeline** (`1.1.1_preprocess.py`): benchmark views now use batch-aware HVG (`batch_key="Sample"` from `SAMPLE_COLNAME`) and harmony-based Leiden; batch views get Leiden on `X_pca_harmony` too (`BATCH_VIEW_RES` removed); `X_pca` stored for every HVG size; Harmony + neighbors + Leiden (6 resolutions) only at the 2000-HVG pass
- [x] **Preprocess array fix**: `--ds_name` filter added to `1.1.1_preprocess.py`; `1.1_run_worker.sh` passes its `DS_NAME` (previously every worker looped over all datasets with its own dirs → guaranteed missing-file failure). View-less datasets (Zhu) still exit cleanly via "No views defined."
- [x] **Joanito HPC step**: `_create_joanito_batch_col.R` is HPC-capable (input from `HPC_SCRATCH_DIR`, in-place idempotent `seqtec` recompute); new `1.2_submit_joanito_batch_col.sh` (`pixi run Rscript`); must run after `1_stage_data.sh`, before the preprocess array
- [x] **Dataset-specific preprocessing dispatcher**: `1_submit_hpc.sh` submits all `1.*_submit_*.sh` steps in parallel (`sbatch --parsable`), waits for all, reports per-job state via `sacct`, exits non-zero on failure; old single-step script renamed `1.1_submit_combinedpbmc.sh` (steps must be mutually independent)
- [x] **Dead dirs removed**: `path_plots`, `path_output_samples`, `path_output_ecoda` deleted from `config_helper.R` + `1.1_prepare_chunks.r` (created but never written to; annotation feathers go to `${HPC_SCRATCH_DIR}/${DS_NAME}/output` directly)
- [x] **`DEPRECATED_LEGACY_CODE/` deleted** (`src/4_cell_type_annotation/`)
- [x] **jq module fix**: `2_submit_hpc_array.sh` loads `module load jq/1.6` before the `TISSUE_TYPE`/`NORMAL_TISSUE` export block (`command -v jq` fallback retained)
- [x] **HPC layout documented**: new "HPC Folder Layout" section in `docs/ARCHITECTURE.md` (tree + per-folder table incl. `CombinedPBMC/data/` dual role, `output/<DS_NAME>/` contents, NAS targets) + summary in `AGENTS.md`

## Step 1 - Check if cell type annotation pipeline can be simplified
- [x] **Wrapper chain collapsed to 3 levels** (mirrors `src/3_scrnaseq_preprocessing/`): `2.1.1_process_chunk.sh` merged into `2.1_run_worker.sh`; `2.1.1.1_process_chunk.R` renamed to `2.1.1_process_chunk.R`; legacy `chunk_file__` arg encoding removed (R script now takes the chunk path as its single positional arg)
- Is `config_helper.R` still needed or should this better be handled in `slurm_config.sh`? (open reviewer point — see `## Open reviewer points`)
  - If `config_helper.R` is easier and cleaner, then check `1.1_prepare_chunks.r` and `2.1.1_process_chunk.R` to see if they can be simplified
- Check if any `*_prepare_chunks.sh` and `*_process_chunk.sh` can be simplified by centralizing any more paths into `slurm_config.sh`



## Other major goals
- Preprocessing and subsequently cell type annotation should be run calling a single script for each.
  - Check if a single unified `1_submit_hpc_array.sh` can be used for `src/3_scrnaseq_preprocessing/` and `src/4_cell_type_annotation/`
    - The script should simply read `slurm_config.sh` and pass the correct arguments to run a pipeline across all (or the specified) datasets to `1.1_run_worker.sh` and `2.1.1_run_worker.sh`
    - Can also be used subsequently for other pipelines, e.g. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/` (TBD)
- After the pipelines are finished, update the documentation and write up the execution plan for the new pipelines.



## Step 2 — Debug Dataset (Joanito 5-sample)

---
Stale:
- Should also be able to test the preprocessing pipeline
  - This should start from the Joanito input dataset (see datasets.json, "input_file_name": "JoaI_2022_35773407_Nofilt_whole.rds")
  - Should be done after running `src/2_dataset_specific_preprocessing/_create_joanito_batch_col.R`, as the batch metadata column is needed for debugging.

### 2a. Locate source
- Existing file: `data/JoaI_2022_35773407_Nofilt_whole_ECODAprocessed.rds` (old-pipeline Seurat output, already QC-filtered)

### 2b. Create subset script
- Write `src/2_dataset_specific_preprocessing/_create_debug_dataset.R` that:
  - Reads the .rds file
  - Selects 5 samples covering both biological conditions and batches
  - Subsets to 500 cells per sample (random)
  - Saves as both:
    - `data/debug/JoaI_2022_35773407_debug_5samples.rds` (RDS for R workflow)
    - `data/debug/JoaI_2022_35773407_debug_5samples.h5ad` (h5ad for Python workflow)
  - Creates minimal obs columns: `sample.ID`, `sample.origin`, `patient.ID`, `Site`, `cell.type`

### 2c. Add debug dataset reference
- Temporarily register in `datasets.json` a `"_debug"` entry with a `benchmark_analysis` and a `batch_effect_analysis` view pointing to the debug .h5ad file
- OR keep the debug entry in a separate `debug_datasets.json` to avoid cluttering the real config?

### Validation
- Output .h5ad loads in `sc.read_h5ad()` without error
- Shape matches: ~2500 cells × genes, with all 5 samples present
- Runs in under 30s on a laptop

---

## Step 3 — Preprocessing Core (preprocess.py + Blacklist + Batch HVGs + Harmony + SLURM Wrappers)

### 3a. Fix datasets.json reading in preprocess.py

Current blocker at line 247-251: `json.load(f)["datasets"]` and `ds_info.get("file_name")`.

Changes:
- Remove `["datasets"]` wrapper — iterate `config.items()` directly
- Read `input_file_name` and `output_file_name` from each `view_info` (per-view config)
- Handle `input_file_name` as either `str` (single file) or `list` (array for Gongsharma):
  - If list, concatenate multiple h5ad files with `sc.concat()` after loading
- Replace `file_name_no_ext` with `stem` derived from the dataset key (not the filename), since output names are now explicit in `output_file_name`
- `output_file_name` is taken from `view_info["output_file_name"]` directly

### 3b. Handle array input_file_name (Gongsharma)

- Write a helper `_load_input(adata_full, input_file_name, base_path)`:
  - If `str`: load single file (existing logic)
  - If `list`: check all exist, load each, `sc.concat()` keep all obs from each partial file
- Support `.rds` (via R interop `convert_rds_to_raw_h5ad`) and `.h5ad` directly

### 3c. Gene blacklist filtering

The blacklist RDS already exists at `aux/genes.blocklist.rds` (STACAS default_black_list).

Approach (no R dependency for this in Python):
1. Write a one-time R extraction script to dump the blacklist as a text file: `aux/genes_blocklist.txt`
   - One gene symbol per line
   - Two variants: `full` (all genes) and `no_sex` (excluding X/Y chromosome genes)
2. In preprocess.py, add a `load_blacklist(path, exclude_sex=True)` function
3. Apply before HVG selection: `adata = adata[:, ~adata.var_names.isin(blacklist)].copy()`

Also:
- For non-batch views, ensure `sc.pp.highly_variable_genes(batch_key="Sample")` is used
- For batch views, `batch_key=batch_col` (from datasets.json)
- Create harmony embeddings for both view types (integrated by sample for benchmark, by batch_col for batch effect)
- Check strategy to handle Gongsharma dataset (huge, provided in chunks)
- Cleanup preprocess_gongsharma.qmd if necessary
- Determine which Gongsharma files need pre-processing (per datasets.json)

### 3d. Batch column support

Currently `preprocess.py:260` looks for `columns.batch` in datasets.json, which doesn't exist for any dataset.

Changes:
- Add `"batch"` key to `datasets.json` for relevant datasets:
  - Stephenson: `"Site"`
  - Joanito: `"Site"` (needs to be constructed from metadata — was done manually in batch_effect_analysis.rmd)
  - Zhu: `null` or leave absent (single batch)
- In `process_view()`: pass `batch_key` correctly based on view type:
  - `batch_effect_analysis` views → use dataset-level `columns.batch`
  - Keep current logic where `batch_key = None` for `benchmark_analysis` views (will be changed to `"Sample"` in 3e)

### 3e. HVG & Harmony improvements — DONE (see "Completed (HPC layout + TODO verification plan)")

- Non-batch views (benchmark_analysis): `batch_key="Sample"` for HVG selection
- Batch views (batch_effect_analysis): `batch_key=batch_col` (from datasets.json)
- Harmony embeddings always computed for batch views (`X_pca_harmony`)
- For benchmark views, also compute harmony embeddings for cell type annotation (not used by benchmark methods directly)
- Change `process_view()` call for `benchmark_analysis` views to use `batch_key="Sample"` (from `columns.sample`)
- Remove the old distinction where only batch views use `batch_key`
- In `run_downstream_for_gene_set()`: always compute `X_pca_harmony` when `batch_key` is set (harmony integration by sample/batch)
- Save harmony results to `obsm` with appropriate key

### 3f. Write SLURM wrappers for preprocess

The two wrapper files are currently empty:

`src/3_scrnaseq_preprocessing/2_submit_hpc_array.sh`:
- Source `src/slurm_config.sh`
- Iterate over datasets from datasets.json (using jq or a Python helper)
- For each dataset+view pair, submit `2.1_run_worker.sh` as an array job
- Each array task = one view of one dataset
- After all array jobs complete, sync back to NAS

`src/3_scrnaseq_preprocessing/2.1_run_worker.sh`:
- `#SBATCH` headers (shared-cpu, 8G, 2h, etc.)
- Source `src/slurm_config.sh`
- Copy data from NAS to scratch (`1_copy_data_from_nas_to_hpc_scratch.sh` or inline)
- Run `2.2_preprocess.py` with the dataset/view as argument
- After completion, copy result back to NAS


### 3g. Fix copy_data_from_nas_to_hpc_scratch.sh

Current script (line 33) uses `.datasets` and `.file_name`:
```
jq -r '.datasets | to_entries[] | "\(.value.folder_name) \(.value.file_name)"'
```

Rewrite to:
```
jq -r 'to_entries[] | .key as $ds | .value.folder_name as $folder | .value.views | to_entries[] | "\($folder) \(.value.input_file_name)"'
```

Handle array `input_file_name` (for Gongsharma, copy all files in the array).

### Validation
- `preprocess.py` runs on the debug dataset without crashing
- Output .h5ad contains expected keys: `X_pca`, `X_pca_harmony` (if batch_key set), HVG masks, cluster annotations
- `process_view` runs for both `benchmark_analysis` and `batch_effect_analysis` views
- Harmony integration produces embeddings with shape `(n_cells, 50)`
- Gene blacklist removes known noisy genes (MT-, MALAT1, etc.)
- SLURM wrappers parse correctly (`sbatch --test-only` or dry-run)

---

## Step 4 — R Downstream Data Loading & Path Adaptations

### 4a. Choose h5ad loading approach in R

Options:
| Approach | Pros | Cons |
|---|---|---|
| **anndataR** | Native R package, no Python runtime | Newer package, may lack features; Seurat v5 support uncertain |
| **reticulate** | Battle-tested, access to full Python anndata API | Requires Python env; slower on large objects |
| **SeuratDisk (h5seurat)** | Direct Seurat conversion | Deprecated-ish, may have issues with v5 objects |
| **Seurat's built-in ReadH5AD** | Available in Seurat v5 | Experimental, may not handle all .h5ad variants |

Decision criteria: the preprocessed .h5ad files need to be loaded into Seurat objects for benchmark methods that require Seurat (ECODA, Pseudobulk, etc.). **Recommended: use `SeuratData::ReadH5AD()` or `reticulate` for initial loading, benchmark on the debug dataset.**

### 4b. Update file-path construction in .rmd files

The old code uses `datasets[[ds]][["ds_name"]]` (file stem) to construct RDS input paths. With the new pipeline, input is `.h5ad`, and the path is constructed from `read_datasets_json()` output.

Current references to update:
- `benchmark_analysis.rmd`: lines ~73, ~1832, ~3749
- `batch_effect_analysis.rmd`: lines ~435, ~475

The new pattern:
```r
ds_config <- read_datasets_json("datasets.json", "benchmark_analysis")
file_path <- file.path("data", ds_config$output_file_name[ds_config$dataset == ds])
```

Or more precisely, the view-specific output filename from datasets.json.

### 4c. Test loading path

- Load the debug dataset's preprocessed .h5ad in R
- Convert to Seurat
- Run `run_analyses()` on the debug dataset (or a subset of methods)
- Verify `process_coda_fig()`, `process_pseudobulk_fig()`, etc. produce expected output

### Validation
- `read_datasets_json("datasets.json", "benchmark_analysis")` returns a data.frame with expected columns
- R can load the preprocessed .h5ad into a Seurat object
- `run_analyses()` runs without file-not-found errors on the debug dataset
- All benchmark method processors receive the correct data

---

## Step 5 — Downstream Python Benchmarks & Batch Mitigation

Re-purpose `benchmark_methods_py.qmd` (and rename) to process py methods for benchmark and batch effect analyses, respectively.
Methods need to be adapted whether they run on benchmark view or batch effect view.

### 5a. Convert benchmark_methods_py.qmd to .py

The .qmd file (Quarto notebook) needs to be converted to a standalone .py script suitable for SLURM batch submission.

Changes:
- Strip Quarto markdown/chunk headers
- Pass dataset name and view as CLI arguments (instead of hardcoded loop)
- Strip redundant preprocessing steps now done in preprocess.py:
  - `sc.pp.normalize_total()`
  - `sc.pp.log1p()`
  - `sc.pp.highly_variable_genes()`
  - `sc.pp.scale()`
  - `sc.pp.pca()`
  → These are already applied in the preprocessed .h5ad; just load and use `adata.obsm["X_pca"]` etc.
- Make sure every method gets the correct input (counts, or specific embedding)
- Double-check compatibility with current preprocess.py workflow (scanpy overwrites adata.X at every step)

### 5b. Add batch correction parameters per method (only for batch effect analysis views)

- **MrVI**: pass `batch_key` to `MRVI.setup_anndata(adata, sample_key="Sample", batch_key=batch_col)` (for batch effect analysis, use current setup for benchmark view)
- **PILOT-GM-VAE**: Needs to be implemented. Read `X_pca_harmony` from preprocessed .h5ad as input embeddings (for batch effect analysis, use current setup from PILOT (use uncorrected PCA) for benchmark view)
- **scPoli**: no batch correction (used for benchmark only, not batch effect)
- **PILOT**: keep as-is (uses uncorrected PCA)
- Add QOT, PULSAR

### 5c. Write SLURM wrappers

`src/5_run_benchmark_methods/run_python_sample_embedding_methods/TODO_STUMP_1_submit_hpc_array.sh`:
- Source `src/slurm_config.sh`
- For each dataset+view combination that needs Python methods, submit an array job
- Job concurrency controlled by `MAX_NUM_CHUNKS_PARALLEL`

`src/5_run_benchmark_methods/run_python_sample_embedding_methods/TODO_STUMP_1.1_run_worker.sh`:
- `#SBATCH` headers (shared-cpu, 16G, 4h or GPU if needed for PILOT-GM-VAE)
- Copy preprocessed .h5ad to scratch
- Run converted .py script with dataset/view/method arguments
- Output `.feather` files back to NAS

### 5d. Harmonize execution time & memory logging

- Check if Python's `log_execution_time_and_memory()` output format matches what `benchmark_analysis.rmd` expects at `p_exec_times`
- If mismatch, align Python output to match R's `exec_time()` format
- Peak memory may be unreliable; consider dropping it or check if it works on HPC

### Validation
- Python benchmark script runs on debug dataset output .h5ad for each method (MrVI, PILOT, scPoli, PILOT-GM-VAE)
- Output `.feather` files are produced and R can read them via `arrow::read_feather()`
- Execution time format matches R's expected input
- SLURM wrappers parse correctly (needs access to HPC cluster. See `AGENTS.md` for details)
  - Debug dataset (subsetted Joanito) dataset can be used for validation


## Step 6 — Update SLURM config

### Step 6a. Update `slurm_config.sh`
Needs to be adjusted:
src/slurm_config.sh:39-39
```
SLURM_PARTITION="shared-cpu"
```
- shared-cpu for:
  - cell_type_annotation
  - preprocess
- specified infrastructure (including gpu for some methods) for:
  - benchmark_analysis (includes a runtime check, so needs clearly defined infrastructure used across all methods), used for:
    - `src/5_run_benchmark_methods/run_python_sample_embedding_methods/`
    - `src/5_run_benchmark_methods/benchmark_pipeline.R`

## Step 7 Update README.md section **Usage / Workflow**
- Update `README.md` with new steps from `cell_type_annotation`

---



## New Datasets to Be Added

- **Batch effect analysis:**
  - Whole Stephenson by batch/center (n = 143)
  - From PILOT-GM-VAE paper:
    - KPMP Kidney subset (n = 45) — needs batch effect check first
    - Breast cancer (n = 126)
    - Covid-19 PBMC (n = 151)
    - Diabetes (n = 52)
    - Possibly: Sikkema Lung (n = 165) — `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd`
- **Benchmark analysis:**
  - From PILOT-GM-VAE paper:
    - Alzheimer (n = 83)
    - Lupus PBMC (n = 261)
    - Myocardial infarction (n = 23)
    - Possibly: Kidney KPMP subset (n = 45)
  - Appendix: GongSharma for other subsetting conditions (author annotations only, no preprocessing needed)

## New Methods to Be Added

- **PILOT-GM-VAE**: very similar to PILOT; add to `benchmark_methods_py.qmd`, `benchmark_methods_r.R`, and `constants.R`
- **QOT**: check https://github.com/PennShenLab/QOT (no package, `qot_utils_re.py` script; read dependencies from `QOT_PDAC_Example.ipynb`)
- **PULSAR**: check https://github.com/snap-stanford/PULSAR/ — does it require UCE? Needs HPC GPU, may exceed resources for larger datasets

---

## Batch Effect Analysis

### Preprocessing
- Handle HVG calculation using correct `batch_key`
- `datasets.json`: add `"columns"` `"batch"` to `batch_effect_analysis` views
- Create low-res cell types for Kfoury dataset (see `Preprocess_datasets.Rmd`)
- For Joanito, create batch column BEFORE preprocessing (from sequencing technology/metadata) so `preprocess.py` can use it as `batch_key` — DONE (`_create_joanito_batch_col.R` + `1.2_submit_joanito_batch_col.sh`)
- Update datasets with batch column mapping in `datasets.json`
- Handle case where multiple datasets are combined (e.g. Combined PBMC in `batch_effect_analysis.rmd`)

### Downstream
- `DESeq2.normalize()`: check `batch_col` is correctly implemented
- `get_pb_deseq2()`: add `batch_col` argument, implement batch correction

### Analysis methods
- **ECODA**: remove cell types significantly different across batches (t-test/Wilcoxon for 2 batches, ANOVA/Kruskal-Wallis for >2 batches)
- **Pseudobulk** (DESeq2 + limma): re-use `process_pseudobulk_fig()` with added `batch_col` argument
- **MrVI**: native batch effect correction
- **GloScope**: run on Harmony-integrated space
- **PILOT-GM-VAE**: run on Harmony-integrated space

---

## Risk Register

| Risk | Mitigation |
|---|---|
| preprocess.py changes touch many interleaved concerns | Implement 3a → 3b → 3c → 3d → 3e → 3f sequentially, testing each |
| Debug dataset may not represent all edge cases (array files, batch columns) | Step 3 validated against debug, then run on one real dataset (e.g., Kim) before full rollout |
| R h5ad loading may fail on complex .h5ad files | Test multiple loading strategies in Step 4a; fall back to reticulate |
| .rmd file-path updates may miss references | Grep for `ds_name` and `processed_file_stem` patterns across both .rmd files after update |
| SLURM wrappers depend on HPC environment (paths, modules) not available locally | Write and validate syntax; actual testing requires HPC access; document test mode |

---

## Ideas for Later

- Possibly run GloScope on the HPC cluster
- Implement MOFAcellular? Needs requirements testing
- How to get final cell/sample/annotation counts from .h5ad (without loading full object into memory)?

---

# Completed

- [x] **Step 1 Foundation**: SLURM config centralized (`src/slurm_config.sh`), files restructured (notebooks/, src/utils/, src/preprocess/, src/benchmark/), all internal paths updated, docs (ARCHITECTURE.md, AGENTS.md, README.md) updated
- [x] **Cell type annotation pipeline**: Adopted from another project (HPC-parallelized), path-adjusted for new src/ structure, HiTME/scATOMIC column whitelisting from Preprocess_datasets.Rmd, scATOMIC added to `2.1.1_process_chunk.sh` and `3_merge_annotations.py`, retry loops with timeout for R worker
- [x] **SLURM config migration**: `src/bash/config.env` → `src/slurm_config.sh`, all bash scripts updated to source centralized file, `config.env` deleted
- [x] **HPC wrappers**: `src/preprocess/1_submit_hpc_array.sh` / `1.1_run_worker.sh` implemented (now `src/3_scrnaseq_preprocessing/`), data copy from NAS to scratch absorbed into submit script, `copy_data_from_nas_to_hpc_scratch.sh` deleted
- [x] **File migration**: QC_filtering/ → notebooks/, .rmd files moved, `src/R/` → `src/utils/` (11 files), `src/bash/cell_type_annotation/` → `src/cell_type_annotation/`, stale files deleted
- [x] **Step 3d Follow-up**: Plan verified, DRY `_preprocess_utils.py` extracted (now `src/utils/preprocess_utils.py`), RAM optimization (in-place gene subsetting, early obs trimming, `del`/`gc`), NAS syncing confirmed for preprocessing and cell type annotation, `config_helper.R` fixed, `run_worker.sh` memory bumped to 16G baseline
- [x] **Shared datasets.json reader**: Python `src/datasets_io.py::read_datasets_json()` (stdlib-only) + enriched `src/utils/datasets_io.R::read_datasets_json()` (harmonized key structure, backward compatible), `file_names` added to all datasets.json entries, Zhu re-added (view-less raw source for CombinedPBMC), preprocess + combine scripts refactored to use the shared reader, Stephenson combine source switched to `benchmark_analysis` view (Site = Ncl)
- [x] Raw-data staging consolidated in `src/preprocess/1_submit_hpc_array.sh` (per-dataset dirs, jq module loads moved to top); redundant staging removed from `src/cell_type_annotation/1_prepare_chunks.sh`
- [x] **Folder rename to pipeline sequence** (`src/1_stage_data`, `src/2_dataset_specific_preprocessing`, `src/3_scrnaseq_preprocessing`, `src/4_cell_type_annotation`, `src/5_run_benchmark_methods`): `git mv` (history preserved), staging extracted from `1_submit_hpc_array.sh` into `src/1_stage_data/1_stage_data.sh`, `_preprocess_utils.py` → `src/utils/preprocess_utils.py`, all imports/references/docs updated (AGENTS.md, README.md, docs/ARCHITECTURE.md, TODO.md)
- [x] **Reviewer points (priority order, in final locations)**:
    - `NAS_TARGET_DIR` → `${NAS_PREFIX}/Projects/ECODA_paper` (`slurm_config.sh`); rsync steps in `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` + `4_cell_type_annotation/2_submit_hpc_array.sh` unchanged (use the var)
    - `SAMPLE_COLNAME` kept as env var in `slurm_config.sh`; `1.1.1_preprocess.py` now writes `os.environ.get("SAMPLE_COLNAME", "Sample")` (was hardcoded `"Sample"`)
    - `1.1.1_preprocess.py` CLI verified (argparse `--config_path`/`--base_path`/`--output_dir` honored by `main()`, relative args joined against `PROJECT_ROOT`, absolute kept)
    - `_create_combinedpbmc_dataset.py` verified HPC-capable (per-dataset layout default when `HPC_SCRATCH_DIR` set); optional single-`sbatch` wrapper added: `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
    - `getwd()` fallbacks removed: `1.1_prepare_chunks.r` + `2.1.1.1_process_chunk.R` use `stop()` guards on `PROJECT_ROOT`; `config_helper.R` `home_ref_dir`/`gene_ref` fallbacks → `stop()` guards; `GENE_REF_FILE`/`GENE_REF_URL` removed from `slurm_config.sh` (unused; gene standardization is in `src/gene_utils.py` via `aux/EnsemblGenes...`), gene-annotation staging block removed from `1_prepare_chunks.sh`
    - `config_helper.R` kept at project root (R-facing env-based config; R cannot source `slurm_config.sh`); `LOGS_DIR` added to `slurm_config.sh`, all `logs/` references routed through it; cell type annotation array logs renamed to `4_cell_type_annotation_%A_%a.{log,err}`
    - Fixed pre-existing relative-script-reference bugs in `4_cell_type_annotation/2_submit_hpc_array.sh` (`sbatch` worker path) and `2.1_run_worker.sh` (`bash` chunk-processor path) — now resolved via `${SCRIPT_DIR}`/`dirname ${BASH_SOURCE[0]}`
- [x] **Replace base_path / project_root with slurm_config.sh env vars** (per-dataset layout):
    - `src/slurm_config.sh` now exports all core env vars (`DATASETS_JSON_FILE`, `NAS_SC_DIR`, `NAS_TARGET_DIR`, `NAS_REF_DIR`, `HOME_REF_DIR`, `GENE_REF_FILE`, `PIXI_R_LIB`, `HPC_SCRATCH_DIR`) + new `SAMPLE_COLNAME="Sample"` so they reach R (`Sys.getenv()`) and Python (`os.environ`) through `srun`/`sbatch`
    - Decision: the TODO example's flat `SCRATCH_DATA_DIR=${HPC_SCRATCH_DIR}/data` / `SCRATCH_CHUNKS_DIR=${HPC_SCRATCH_DIR}/chunks` were NOT added — they conflict with the per-dataset convention; dataset-specific dirs are derived per script (`${HPC_SCRATCH_DIR}/${DS_NAME}/data`, `${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks`)
    - `config_helper.R` paths are now per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (annotation pipeline functional end-to-end, matches `2_submit_hpc_array.sh` CHUNKS_DIR)
    - R scripts (`1.1_prepare_chunks.r`, `2.1.1.1_process_chunk.R`) resolve `project_root` from `PROJECT_ROOT` env (getwd() fallback), source `config_helper.R` explicitly, SAMPLE_COLNAME guard added
    - `1.1.1_preprocess.py`: explicit absolute-vs-relative path resolution + argparse for `--config_path/--base_path/--output_dir` (was silently ignoring CLI args)
    - `_create_combinedpbmc_dataset.py` is now HPC-capable (`--layout per-dataset` default when `HPC_SCRATCH_DIR` set; must run before the preprocess array, from `${PROJECT_ROOT}` with `module load GCCcore/12.2.0`)
    - `2_submit_hpc_array.sh` auto-exports per-dataset `TISSUE_TYPE`/`NORMAL_TISSUE` from `datasets.json`; `2.1_run_worker.sh`/`2.1.1_process_chunk.sh` no longer pass PROJECT_ROOT as an arg (sourced from slurm_config.sh)
- [x] **Per-dataset HPC output layout**: preprocessed outputs moved from the global `${HPC_SCRATCH_DIR}/output/<DS_NAME>/` to `${HPC_SCRATCH_DIR}/<DS_NAME>/output/`; NAS mirror is per-dataset too (`${NAS_TARGET_DIR}/<DS_NAME>/output/`); `SCRATCH_OUTPUT_DIR` env var removed entirely — all consumers derive the per-dataset output dir from `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (preprocess worker, `1_prepare_chunks.sh`/`1.1_prepare_chunks.py`, `config_helper.R`); both submit scripts now rsync per-dataset (glob loop over `${HPC_SCRATCH_DIR}/*/output` with `[[ -d ]]` guard, `set -euo pipefail` safe); `chunks_manifest.txt` moved to `${HPC_SCRATCH_DIR}/chunks_manifest.txt`; docs updated (ARCHITECTURE.md, AGENTS.md)
