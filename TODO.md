# Pipeline Overhaul Execution Plan

## Overview

Five ordered steps to modernize the ECODA_paper pipeline: foundation restructure → debug dataset → preprocess.py rewrite → R downstream adaptations → Python benchmark wrappers. Each step is validated before the next begins.

---

## Step 1 — Foundation: SLURM Config & Directory Restructure

### 1a. Centralize SLURM environment variables
- Create `src/slurm_bash_config.env` from `src/bash/config.env`
- Add: `SLURM_ACCOUNT`, `SLURM_PARTITION`, `MAX_NUM_CHUNKS_PARALLEL`, `USER_EMAIL`, `NAS_TARGET_DIR`, `SCRATCH_OUTPUT_DIR`, `HOME_CHUNKS_DIR`
- Source this file from all bash scripts that need it

### 1b. Move files to target structure

Target layout (from TODO.md):

```
ECODA_paper/
├── notebooks/
│   ├── QC_filtering/               # entire folder + contents
│   ├── batch_effect_analysis.rmd
│   └── benchmark_analysis.rmd
└── src/
    ├── slurm_bash_config.env
    ├── utils/                      # R files: constants.R, helpers.R, hvcs.R, ...
    ├── preprocess/
    │   ├── 1_copy_data_from_nas_to_hpc_scratch.sh
    │   ├── TODO_STUMP_2_submit_hpc_array.sh → 2_submit_hpc_array.sh
    │   ├── TODO_STUMP_2.1_run_worker.sh → 2.1_run_worker.sh
    │   ├── 2.2_preprocess.py       # moved from src/py/preprocess.py
    │   ├── TODO_STUMP_preprocess_sikkema.qmd
    │   └── preprocess_gongsharma.qmd
    ├── cell_type_annotation/       # stays as-is (already in src/bash/)
    └── benchmark/
        ├── benchmark_methods_r.R
        ├── benchmark_pipeline.R
        └── run_python_sample_embedding_methods/
            ├── TODO_STUMP_1_submit_hpc_array.sh
            ├── TODO_STUMP_1.1_run_worker.sh
            └── 1.2_benchmark_methods_py.qmd
```

Move operations (git mv):
1. `QC_filtering/` → `notebooks/QC_filtering/`
2. `batch_effect_analysis.rmd` → `notebooks/batch_effect_analysis.rmd`
3. `benchmark_analysis.rmd` → `notebooks/benchmark_analysis.rmd`
4. `src/R/*.R` → `src/utils/*.R` (each of the 13 .R files)
5. `src/py/preprocess.py` → `src/preprocess/2.2_preprocess.py`
6. `src/py/DRAFT_BARE_preprocess_sikkema.qmd` → `src/preprocess/TODO_STUMP_preprocess_sikkema.qmd`
7. `src/py/preprocess_gongsharma.qmd` → `src/preprocess/preprocess_gongsharma.qmd`
8. `src/bash/copy_data_from_nas_to_hpc_scratch.sh` → `src/preprocess/1_copy_data_from_nas_to_hpc_scratch.sh`
9. `src/bash/preprocess/submit_job.sh` → `src/preprocess/2_submit_hpc_array.sh` (rename + implement)
10. `src/bash/preprocess/run_worker_preprocess.sh` → `src/preprocess/2.1_run_worker.sh` (rename + implement)
11. `src/bash/config.env` → `src/slurm_bash_config.env` (rename + extend)
12. `src/bash/cell_type_annotation/*` → `src/cell_type_annotation/*` (move up one level)
13. `src/R/load_all_functions.R` → `src/utils/load_all_functions.R`
14. `src/R/benchmark_pipeline.R` → `src/benchmark/benchmark_pipeline.R`
15. `src/R/benchmark_methods_r.R` → `src/benchmark/benchmark_methods_r.R`
16. `src/py/benchmark_methods_py.qmd` → `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`

### 1c. Update all internal path references
- `src/R/load_all_functions.R` sources → `src/utils/`
- `preprocess.py` line 14: `ro.r('source("src/R/load_all_functions.R")')` → `src/utils/load_all_functions.R`
- `src/bash/cell_type_annotation/*` hardcoded paths → use `src/slurm_bash_config.env`
- `notebooks/batch_effect_analysis.rmd` and `benchmark_analysis.rmd` source calls → `src/utils/`
- `copy_data_from_nas_to_hpc_scratch.sh` dataset iteration → use new datasets.json schema

### 1d. Update documentation
- `docs/ARCHITECTURE.md`: update file paths, module table, call graph
- `AGENTS.md`: update "Pipeline Overview" directory paths, "R Modules" table
- `README.md` if needed

### Validation
- `source("src/utils/load_all_functions.R")` executes without error in R
- `ls` confirms all target directories exist with correct files
- All stale files at old locations are removed

---

## Step 2 — Debug Dataset (Joanito 5-sample)

### 2a. Locate source
- Existing file: `data/JoaI_2022_35773407_Nofilt_whole_ECODAprocessed.rds` (old-pipeline Seurat output, already QC-filtered)

### 2b. Create subset script
- Write `src/preprocess/_create_debug_dataset.R` that:
  - Reads the .rds file
  - Selects 5 samples covering both biological conditions and batches
  - Subsets to 500 cells per sample (random)
  - Saves as both:
    - `data/debug/JoaI_2022_35773407_debug_5samples.rds` (RDS for R workflow)
    - `data/debug/JoaI_2022_35773407_debug_5samples.h5ad` (h5ad for Python workflow)
  - Creates minimal obs columns: `sample.ID`, `sample.origin`, `patient.ID`, `Site`, `cell.type`

### 2c. Add debug dataset reference
- Temporarily register in `datasets.json` a `"_debug"` entry with a single `benchmark_analysis` view pointing to the debug .h5ad file
- OR keep the debug entry in a separate `debug_datasets.json` to avoid cluttering the real config

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

### 3d. Batch column support

Currently `preprocess.py:260` looks for `columns.batch` in datasets.json, which doesn't exist for any dataset.

Changes:
- Add `"batch"` key to `datasets.json` for relevant datasets:
  - Stephenson: `"Site"` (line 10)
  - Joanito: `"Site"` (needs to be constructed from metadata — was done manually in batch_effect_analysis.rmd)
  - Zhu: `null` or leave absent (single batch)
- In `process_view()`: pass `batch_key` correctly based on view type:
  - `batch_effect_analysis` views → use dataset-level `columns.batch`
  - Keep current logic where `batch_key = None` for `benchmark_analysis` views (will be changed to `"Sample"` in 3e)

### 3e. HVG & Harmony improvements

The TODO specifies:
- Non-batch views (benchmark_analysis): `batch_key="Sample"` for HVG selection
- Batch views (batch_effect_analysis): `batch_key=batch_col` (from datasets.json)
- Harmony embeddings always computed for batch views (`X_pca_harmony`)
- For benchmark views, also compute harmony embeddings for cell type annotation (not used by benchmark methods directly)

Changes:
- Change `process_view()` call for `benchmark_analysis` views to use `batch_key="Sample"` (from `columns.sample`)
- Remove the old distinction where only batch views use `batch_key`
- In `run_downstream_for_gene_set()`: always compute `X_pca_harmony` when `batch_key` is set (harmony integration by sample/batch)
- Save harmony results to `obsm` with appropriate key

### 3f. Write SLURM wrappers for preprocess

The two wrapper files are currently empty:

`src/preprocess/2_submit_hpc_array.sh`:
- Source `src/slurm_bash_config.env`
- Iterate over datasets from datasets.json (using jq or a Python helper)
- For each dataset+view pair, submit `2.1_run_worker.sh` as an array job
- Each array task = one view of one dataset
- After all array jobs complete, sync back to NAS

`src/preprocess/2.1_run_worker.sh`:
- `#SBATCH` headers (shared-cpu, 8G, 2h, etc.)
- Source `src/slurm_bash_config.env`
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

### 5b. Add batch correction parameters per method

- **MrVI**: pass `batch_key` to `MRVI.setup_anndata(adata, sample_key="Sample", batch_key=batch_col)`
- **PILOT-GM-VAE**: read `X_pca_harmony` from preprocessed .h5ad as input embeddings
- **scPoli**: no batch correction (used for benchmark only, not batch effect)
- **PILOT**: keep as-is (uses uncorrected PCA)

### 5c. Write SLURM wrappers

`src/benchmark/run_python_sample_embedding_methods/TODO_STUMP_1_submit_hpc_array.sh`:
- Source `src/slurm_bash_config.env`
- For each dataset+view combination that needs Python methods, submit an array job
- Job concurrency controlled by `MAX_NUM_CHUNKS_PARALLEL`

`src/benchmark/run_python_sample_embedding_methods/TODO_STUMP_1.1_run_worker.sh`:
- `#SBATCH` headers (shared-cpu, 16G, 4h or GPU if needed for PILOT-GM-VAE)
- Copy preprocessed .h5ad to scratch
- Run converted .py script with dataset/view/method arguments
- Output `.feather` files back to NAS

### 5d. Harmonize execution time & memory logging

- Check if Python's `log_execution_time_and_memory()` output format matches what `benchmark_analysis.rmd` expects at `p_exec_times`
- If mismatch, align Python output to match R's `exec_time()` format
- Peak memory may be unreliable; consider dropping it entirely

### Validation
- Python benchmark script runs on debug dataset output .h5ad for each method (MrVI, PILOT, scPoli, PILOT-GM-VAE)
- Output `.feather` files are produced and R can read them via `arrow::read_feather()`
- Execution time format matches R's expected input
- SLURM wrappers parse correctly

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

## Quick Reference: Files to Create vs Modify

| File | Action |
|---|---|
| `src/slurm_bash_config.env` | **Create** (from config.env + new vars) |
| `src/preprocess/2.2_preprocess.py` | **Move + rewrite** (from src/py/preprocess.py) |
| `src/preprocess/2_submit_hpc_array.sh` | **Create** (empty → full implementation) |
| `src/preprocess/2.1_run_worker.sh` | **Create** (empty → full implementation) |
| `src/preprocess/1_copy_data_from_nas_to_hpc_scratch.sh` | **Move + rewrite** (from src/bash/) |
| `src/preprocess/_create_debug_dataset.R` | **Create** |
| `aux/genes_blocklist.txt` | **Create** (extract from genes.blocklist.rds) |
| `datasets.json` | **Modify** (add `columns.batch`, add debug entry) |
| `notebooks/benchmark_analysis.rmd` | **Move + modify** |
| `notebooks/batch_effect_analysis.rmd` | **Move + modify** |
| `src/benchmark/benchmark_methods_r.R` | **Move** (from src/R/) |
| `src/benchmark/benchmark_pipeline.R` | **Move** (from src/R/) |
| `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` | **Move + convert to .py** |
| Remaining `src/R/*.R` → `src/utils/*.R` | **Move** (13 files) |
| `docs/ARCHITECTURE.md` | **Modify** |
| `AGENTS.md` | **Modify** |




---------------------------------------------------


# TBD (HPC cluster overhaul — out of scope for initial per-view migration)
- Update `preprocess.py` to read flat datasets.json (no `["datasets"]` wrapper), use `input_file_name`/`output_file_name` per view, and handle array `input_file_name` for multi-file datasets (e.g. Gongsharma).
    - Currently preprocess.py reads `json.load(f)["datasets"]` (KeyError) and `ds_info.get("file_name")` (None), so it crashes on the new datasets.json.
    - The new bash wrappers (`run_worker_preprocess.sh`, `submit_job.sh`) need to be implemented as part of this overhaul (see "Update preprocessing" below).
- Update `copy_data_from_nas_to_hpc_scratch.sh` to use the new `datasets.json` structure.

## Integrate multi-file preprocessing (e.g. Gongsharma h5ad files) into preprocess.py
- Move the downsample_by_group() logic from preprocess_gongsharma.qmd into preprocess.py
- Support array input_file_name in datasets.json views
- This may require additional fields in datasets.json (out of scope for the initial per-view migration)



# Implement h5ad loading in R (for preprocessed .h5ad files produced by preprocess.py)
- Evaluate and implement the best approach to load h5ad files in benchmark_analysis.rmd and batch_effect_analysis.rmd:
  - anndataR (R package for native AnnData loading)
  - reticulate (Python anndata via rpy2)
  - Convert h5ad to Seurat object (since some benchmarked R methods may need Seurat objects)
- Update the loading code in both .rmd files once the approach is chosen
- Additionally, the old `datasets[[ds]][["ds_name"]]` (file stem) field no longer exists in the new `read_datasets_json()` output. The RDS file-path construction in both .rmd files (5 references: benchmark_analysis.rmd lines 73, 1832, 3749; batch_effect_analysis.rmd lines 435, 475) needs to be updated as part of the migration to preprocessed .h5ad files.





# How to get the final number of cells and samples and cell type annotations?
- Might be easiest to simply get them from .h5ad files (is this possible without loading the full object into memory?) (instead of adding them to the datasets.json file, which would cause it to bload and outdate)


# Update repo structure
Final structure:
ECODA_paper/
├── notebooks/                              
│   ├── QC_filtering/                       # Move whole QC_filtering folder here (including all files it contains)
│   ├── batch_effect_analysis.rmd
│   └── benchmark_analysis.rmd
└── src/
    ├── slurm_bash_config.env               # Renamed from config.env, contains only the SLURM environment variables for all the bash scripts (e.g. paths). -> Pull out paths from src/bash/cell_type_annotation .sh bash scripts into this centralized file (later to also be used for the other .sh bash scripts)
    ├── utils/
    │   ├── constants.R
    │   ├── helpers.R
    │   ├── hvcs.R
    │   ├── load_all_functions.R
    │   ├── math_utils.R
    │   ├── plotting.R
    │   ├── pseudobulk.R
    │   ├── scoring_metrics.R
    │   └── seurat_utils.R
    ├── preprocess/                         # Self-contained preprocessing
    │   ├── 1_copy_data_from_nas_to_hpc_scratch.sh
    │   ├── TODO_STUMP_2_submit_hpc_array.sh
    │   ├── TODO_STUMP_2.1_run_worker.sh
    │   ├── 2.2_preprocess.py
    │   ├── TODO_STUMP_preprocess_sikkema.qmd
    │   └── preprocess_gongsharma.qmd
    ├── cell_type_annotation/               # Self-contained annotation
    │   ├── 1.1_prepare_chunks.r
    │   ├── 1.2_prepare_chunks.sh
    │   ├── 2_submit_hpc_array.sh
    │   ├── 2.1_run_worker.sh
    │   └── 2.2_process_chunk.sh
    └── benchmark/                          # Self-contained benchmark
        ├── benchmark_methods_r.R
        ├── benchmark_pipeline.R
        └── run_python_sample_embedding_methods/
                ├── TODO_STUMP_1_submit_hpc_array.sh # Needs to be created (only empty file, implementation follows later)
                ├── TODO_STUMP_1.1_run_worker.sh # Needs to be created  (only empty file, implementation follows later)
                └── 1.2_benchmark_methods_py.qmd # (for now, only move this file here)(-> Needs to be converted to .py file and adapted to run on the HPC cluster to be run per dataset per method)

## IMPORTANT
- For every file or folder that needs to be moved:
    - check if it is used in any of the other files, then move it and update the code (path/dependencies) accordingly
    - Update docs/ARCHITECTURE.md, AGENTS.md and README.md if necessary
    - After moving the file/folder, double-check again that the code (path/dependencies) is updated correctly everywhere



# HPC cluster implementation
- Centralize SLURM environment variables in `src/slurm_bash_config.env`

## Update cell type annotation
- Cell type annotation bash scripts need to be updated to accomodate:
    - the new preprocessing workflow
    - pixi environment
    - updated paths

## Update preprocessing
- `Preprocess_datasets.Rmd` is superseded by `src/py/preprocess.py`, which still needs to be implemented with bash scripts to be run on the HPC cluster
    - `src/preprocess/TODO_STUMP_2_submit_hpc_array.sh` needs to be adapted to run on the HPC cluster
    - `src/preprocess/TODO_STUMP_2.1_run_worker.sh` needs to be adapted to run on the HPC cluster

## Update benchmark methods that need to be run in python (does not affect benchmark_methods_r.R)

- `src/benchmark/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` Needs to be converted to .py file and adapted to run on the HPC cluster to be run per dataset per method
    - `src/benchmark/run_python_sample_embedding_methods/TODO_STUMP_1_submit_hpc_array.sh` needs to be adapted to run on the HPC cluster
    - `src/benchmark/run_python_sample_embedding_methods/TODO_STUMP_1.1_run_worker.sh` needs to be adapted to run on the HPC cluster

- datasets (including batching strategy) now defined in datasets.json
- remove redundant preprocessing step (e.g. sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat")
        sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata, n_comps=50), etc.)
    - noew implemented already in preprocess.py
    - Make sure that every method gets the correct input (counts, or specific embedding)
    - Double-check whether and how this works with the current preprocess.py workflow (because scanpy overwrites adata.X at every pre-processing step!)
- Make it compatible with the benchmark pipeline (current script is only built for this)
    - Either adapt to also incorporate the batch effect analysis datasets and files or create completely separate script)
        - Batch effect analysis will only be run for MrVI and PILOT-GM-VAE python methods (and some other methods in R)
    - Add batch effect mitigation strategies (e.g. MrVI provides this with a parameter, PILOT_GM_VAE needs X_pca_harmony as input)
        - MrVI has native "nuisance variable" that handles batch effect (see here: https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/MrVI_tutorial.html)
        - Possibly just add batch_col argument to MrVI constructor: MRVI.setup_anndata(adata, sample_key="Sample", batch_key=batch_col) ?
        - scPoli will not be used for batch effect analysis
        - PILOT and PILOT-GM-VAE take as input the X_pca_harmony embeddings
- double-check def log_execution_time_and_memory() output format and if it conforms with the expected input in benchmark_analysis.rmd at p_exec_times (which needs to be combined also with r_exec_times)
    - Suggest plan to harmonize, if necessary
- peak_memory was previously implemented but not running correctly in wsl2 on windows desktop
    - check online for reasons why it might not have worked correctly
    - might work if run on hpc cluster
        - because of different memory allocation and file system handling
        - if each dataset (or dataset-method combination) is run in separate instance
    - otherwise comment out or drop completely (probably cleaner and still backed up in git history)

## Explain new cell type annotation pipeline
Was adopted from another project because previous workflow was in r but parallelization constantly kept freezing workers and no approach was found that could prevent it.
Thus, a new cell type annotation was adopted that can be run on the HPC cluster in parallel for any number of datasets and any number of samples for scalability.
- Moved from Preprocess_datasets.Rmd to ./src/bash/cell_type_annotation/ (added from another project, so needs to be adapted and polished)
    - Add documentation explaining the new pipeline structure and usage
    - Ensure compatibility with current preprocessing outputs (preprocessed h5ad files)
    - Update any hardcoded paths or dataset references to use standardized sample names from preprocess step

## Add small test dataset for debugging `preprocess.py` and other HPC pipelines
- Subset the Joanito dataset to a small number of samples (e.g., 5, covering both biological conditions and batches (sequencing technologies)) for testing and subset to 500 cells per sample (using random subsetting)
- Make it so that the output creates a non-batch as well as a batch view (think about how to best implement this. possibly it should not be added to data.json, but instead only be used in the pipeline)
- check resulting file size of.h5ad file and decide whether to store in the repo or on the nas (and where on the nas)



# New datasets to be added:
Implement new datasets (mainly check QC filtering, column names and check for possible batch effects on UMAP, the rest will run in pipeline)
- batch effect analysis:
    - Whole Stephenson by batch/center (n = 143)
    - From PILOT-GM-VAE paper:
        - KPMP Kidney (subset used in PILOT-GM-VAE paper)(needs to be checked for batch effects first) (n = 45)(alternatively with full dataset)
        - Breast cancer (n = 126)
        - Covid-19 PBMC (n = 151)
        - Diabetes (n = 52)
        - Possibly: Sikkema Lung (n = 165)
    - `src/py/DRAFT_BARE_preprocess_sikkema.qmd` draft/exploratory
- benchmark analysis:
    - From PILOT-GM-VAE paper:
        - Alzheimer (n = 83)
        - Lupus PBMC (n = 261)
        - Myocardial infarction (n = 23)
        - Possibly: Kidney KPMP (subset used in PILOT-GM-VAE paper)(needs to be checked for batch effects first) (n = 45)
    - For appendix:
        - GongSharma dataset for other subsetting conditions for ECODA based on author annotations only (does not need pre-processing)

# New methods to be added:
- PILOT-GM-VAE (very similar to PILOT, needs to be added by benchmark_methods_py.qmd and benchmark_methods_r.R and constants.R)
- Check if possible:
    - QOT (https://github.com/PennShenLab/QOT, no package, just qot_utils_re.py script file. Dependencies need to be read out from e.g. QOT_PDAC_Example.ipynb)
    - PULSAR (https://github.com/snap-stanford/PULSAR/, does it require UCE (https://github.com/snap-stanford/UniversalCellEmbedding)?)
        - Needs to be run on the HPC cluster, possibly exceeding GPU requirements for larger datasets




# preprocess.py
- Add blacklist as default for filtering genes
    - Load gene blacklist from file (e.g. from STACAS, see call to `default_black_list` in get_pb_deseq2 in src/R/pseudobulk.R)
        - Maybe save blacklist file to this repo for clarity and add explanation
    - Filter out blacklisted genes before HVG calculation
Possible solutions (but using the blacklist from STACAS):
# Example: Identify mitochondrial genes
is_mito = adata.var_names.str.startswith('MT-')
# Invert the mask to keep non-blacklisted features
adata = adata[:, ~is_mito].copy()
--------------------------------------------------
# Define your blacklisted features
blacklisted_features = ['MALAT1', 'HBB-BS', 'HBZ']  # example genes

# Create a boolean mask for features NOT in the blacklist
keep_features = ~adata.var_names.isin(blacklisted_features)

# Subset the AnnData object
adata = adata[:, keep_features].copy()
--------------------------------------------------
# Mark blacklisted genes in var dataframe
adata.var['blacklisted'] = adata.var_names.isin(blacklisted_features)

# Subset to keep only those where 'blacklisted' is False
adata = adata[:, ~adata.var['blacklisted']].copy()
--------------------------------------------------

- hvgs: for non-batch views, make sure that sc.pp.highly_variable_genes(batch_key="Sample") is used
- need to run and create new harmony embeddings (integrated by samlpe) based on pca embeddings for the "benchmark_analysis" views and create cell type annotations based on unsupervised clustering based on pca and additionally also based on harmony embeddings
- for Batch_effect.Rmd, ensure it uses the updated preprocessing pipeline with batch-aware HVG calculation and new harmony embeddings
- check strategy to handle gongsharma dataset. it's huge, that's why the authors provided the datasets in smaller chunks of .h5ad files subsetted by gender, cmv and age
    - cleanup preprocess_gongsharma.qmd if necessary (still contains legacy code for other conditions that are not used in the current analysis (see datasets.json for most up-to-date list of datasets used and conditions))
    - it does not make sense to combine into one file, so keep it separate
    - check which files need to be actually pre-processed (in datasets.json, do not change this file)



# Batch effect
- Double-check and updated if necessary:
    - "views": "benchmark_analysis" also requires HVGs to be calculated using batch_key="Sample"
    - "views": "benchmark_analysis" also requires to calculate harmony embeddings (but benchmark methods will use uncorrected embeddings or counts. "X_pca_harmony" is only required for a new additional analysis for unsupervised clustering cell type annotation (to be added to benchmark_analysis.rmd. The goal is to compare unsupervised clustering methods on both corrected and uncorrected embeddings).)
    - "views": "batch_effect" requires HVGs to be calculated using batch_key="batch"

## Background information
- For the preprint, the batch effect analysis was minimal (including only the Joanito dataset and the "Combined PBMC (Stephenson, GongSharma, Zhu)").
- The code implementation was drafty, partly because it was just two datasets and partly because we did not make it a major point. However, now the code and batch effect analysis needs to be expanded and handled more cleanly.
- The reviewers raised this as a major point (probably the most important point) to be improve added after reviewing the preprint.

# Batch effect analysis dataset info

Should it be done once without batch correction, and once with?
- It is more important to do WITH batch correction
- Show biological of non-batch corrected results for paper appendix?

The final analysis for batch effect correction needs to be run on the following methods and batch effect mitigation strategies:
- ECODA: remove cell types that are significantly different across batches
    - Possibly re-use process_coda_fig() with added argument batch_col = NULL?
    - Use metadata column batch_col to test for batch-associated cell types (after clr-transformation)
        - Possible methods:
            - Simplest method: remove them and re-calculate clr-transformed cell type composition
            - More complex method: batch correct applying something like a (mixed) linear model or similar method. Could limma be used?
        - Which statistical method to use and which threshold to use for significant difference between batches? -> Print warning naming cell types and Significance. -> Possibly depends on the number of batches:
            - Test every cell type separately or, possibly better, check global variance of cell type composition across batches and see if it is different.
            - If 2 batches: use t-test or Wilcoxon rank-sum test, p-value < 0.05
            - If >2 batches: use ANOVA or Kruskal-Wallis test, p-value < 0.05
- Pseudobulk (DESeq2 + limma)
    - Possibly re-use `process_pseudobulk_fig()` with added argument `batch_col`?
- MrVI (native batch effect correction)
- GloScope (run on harmony integrated space)
- PILOT-GM-VAE (run on harmony integrated space)


## Preprocessing
- handle hvg calculation using correct batch_key
    - datasets.json: add "columns" "batch" to batch_effect_analysis views
    - Create low res cell types for Kfoury dataset (see Preprocess_datasets.Rmd for details)
    - For Joanito, the batch column must be created BEFORE preprocessing in preprocess.py (currently defined manually at end of batch_effect_analysis.rmd). The column needs to exist in the input data so preprocess.py can use it as batch_key for HVG selection and harmony integration.
        - Approach: create batch column (e.g. combining sequencing technology/metadata fields) as part of data preparation, before preprocess.py runs.
    - update datasets with batch column mapping in datasets.json
    - handle case where multiple datasets are combined (e.g. in batch_effect_analysis.rmd at "# Combined PBMC (Stephenson, GongSharma, Zhu)")

## Down-stream
- DESeq2.normalize():
    - check batch_col is correctly implemented
    - check that it does not get "Sample" as a batch column
- get_pb_deseq2(): add argument batch_col, implement batch_col

## Analysis


# Ideas for later:
- Possibly run GloScope on the HPC cluster? Leave it for now
- Imp,ement MOFAcellular?
    - MOFAcellular needs to be tested for requirements