# -- RESOLVED: Gene name standardization --
# bionty was replaced with a flat-file Ensembl 105 reference
# (aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz), the same reference used by
# STACAS. No compatibility concerns with scATOMIC/HiTME.
# See _load_ensembl105_map() + standardize_gene_symbols() in src/py/preprocess.py.


# TBD (HPC cluster overhaul — out of scope for initial per-view migration)
- Update `preprocess.py` to read flat datasets.json (no `["datasets"]` wrapper), use `input_file_name`/`output_file_name` per view, and handle array `input_file_name` for multi-file datasets (e.g. Gongsharma).
    - Currently preprocess.py reads `json.load(f)["datasets"]` (KeyError) and `ds_info.get("file_name")` (None), so it crashes on the new datasets.json.
    - The new bash wrappers (`run_worker_preprocess.sh`, `submit_job.sh`) need to be implemented as part of this overhaul (see "Update preprocessing" below).
- Update `copy_data_from_nas_to_hpc_scratch.sh` to use the new `datasets.json` structure.

# Implement h5ad loading in R (for preprocessed .h5ad files produced by preprocess.py)
- Evaluate and implement the best approach to load h5ad files in benchmark_analysis.rmd and batch_effect_analysis.rmd:
  - anndataR (R package for native AnnData loading)
  - reticulate (Python anndata via rpy2)
  - Convert h5ad to Seurat object (since some benchmarked R methods may need Seurat objects)
- Update the loading code in both .rmd files once the approach is chosen
- Additionally, the old `datasets[[ds]][["ds_name"]]` (file stem) field no longer exists in the new `read_datasets_json()` output. The RDS file-path construction in both .rmd files (5 references: benchmark_analysis.rmd lines 73, 1832, 3749; batch_effect_analysis.rmd lines 435, 475) needs to be updated as part of the migration to preprocessed .h5ad files.


# Integrate multi-file preprocessing (e.g. Gongsharma h5ad files) into preprocess.py
- Move the downsample_by_group() logic from preprocess_gongsharma.qmd into preprocess.py
- Support array input_file_name in datasets.json views
- This may require additional fields in datasets.json (out of scope for the initial per-view migration)




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