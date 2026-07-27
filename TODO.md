# General

> All items in this section have been addressed. See commit history for details.

- [x] README.md updated: Installation, Usage/Workflow, Expected Outputs, paper context (Overview, Key Findings)
- [x] Paper read and incorporated into documentation
- [x] AGENTS.md updated: Pipeline Overview, R Module Table, Domain Terminology, reviewer priorities, HPC notes, kilo code indexing note
- [x] File naming conventions audited; critical issues fixed (spaces removed from 2 QC_filtering filenames)
- [x] Information partitioning strategy defined: README (human), AGENTS (agent), TODO (actionable)
- [x] Script consolidation audited:
  - No duplicate R functions across 12 modules
  - No duplicate Python scripts (preprocess_gongsharma.qmd is dataset-specific; benchmark_methods_py.qmd is general)
  - `Preprocess_datasets.Rmd` orphaned (superseded by `src/py/preprocess.py`) — safe to remove
  - `Process_data.ipynb` does not exist (replaced by `benchmark_methods_py.qmd`)
  - `src/py/DRAFT_BARE_preprocess_sikkema.qmd` draft/exploratory — clean up when finalized
  - `Figure_workflow_schematic.Rmd` orphaned (figure gen, unreferenced) — document or move
  - `Processed_dataset_metadata.R` overlaps `datasets.json` — consolidate if needed
  - 22 `plots/DELETE_*` files are stale — safe to remove

## Explain new cell type annotation pipeline
Was adopted from another project because previous workflow was in r but parallelization constantly kept freezing workers and no approach was found that could prevent it.
Thus, a new cell type annotation was adopted that can be run on the HPC cluster in parallel for any number of datasets and any number of samples for scalability.
- Moved from Preprocess_datasets.Rmd to ./src/bash/cell_type_annotation/ (added from another project, so needs to be adapted and polished)
    - Add documentation explaining the new pipeline structure and usage
    - Ensure compatibility with current preprocessing outputs (preprocessed h5ad files)
    - Update any hardcoded paths or dataset references to use standardized sample names from preprocess step

## Explain migration from R pre-processing in Preprocess_datasets.Rmd to Python in src/py/preprocess.py
- Compare the R-based preprocessing workflow (remove_low_cellcount_samples, standardize_sample_names, create_clean_seuratv5_object, gene standardization) with the Python implementation using Scanpy
- Document the mapping of R Seurat functions to Scanpy equivalents (e.g., sc.pp.filter_cells -> remove_low_cellcount_samples, sc.pp.highly_variable_genes -> standardize variable genes)
- Explain how the gene standardization step (STACAS in R) is handled in Python (bionty's Gene.standardize)
- Detail the data format conversions (Seurat v5 to AnnData) and how metadata columns are preserved (previously metadata had to be subset as the obj.list <- SplitObject(seurat, split.by = "Sample") was used to split by sample to process each sample separately and then re-merge everything and JoinLayers (which every additional metadata column added degrades speed exponentially the more samples there are (with too many columns this never finished for some datasets)) -> problem solved by using .h5ad anndata objects that can be queried (already implemented))
- Note any functional differences or improvements made during the migration



## New datasets to be added:
- batch effect analysis:
    - Whole Stephenson by batch/center (n = 143)
    - From PILOT-GM-VAE paper:
        - KPMP Kidney (subset used in PILOT-GM-VAE paper)(needs to be checked for batch effects first) (n = 45)(alternatively with full dataset)
        - Breast cancer (n = 126)
        - Covid-19 PBMC (n = 151)
        - Diabetes (n = 52)
        - Possibly: Sikkema Lung (n = 165)
- benchmark analysis:
    - From PILOT-GM-VAE paper:
        - Alzheimer (n = 83)
        - Lupus PBMC (n = 261)
        - Myocardial infarction (n = 23)
        - Possibly: Kidney KPMP (subset used in PILOT-GM-VAE paper)(needs to be checked for batch effects first) (n = 45)

## New methods to be added:
- PILOT-GM-VAE (very similar to PILOT, needs to be added by agent to Process_data.ipynb)
- MOFAcellular
    - MOFAcellular needs to be tested for requirements


# Add small test dataset for debugging `preprocess.py` and HPC pipeline
- Subset the Joanito dataset to a small number of samples (e.g., 5, covering both biological conditions and batches (sequencing technologies)) for testing and subset to 500 cells per sample (using random subsetting)
- Make it so that the output creates a non-batch as well as a batch view (think about how to best implement this. possibly it should not be added to data.json, but instead only be used in the pipeline)
- check resulting file size of.h5ad file and decide whether to store in the repo or on the nas (and where on the nas)


# preprocess.py
- Add blacklist as default for filtering genes
    - Load gene blacklist from file (e.g. from STACAS, see call to `default_black_list` in get_pb_deseq2 in src/R/pseudobulk.R)
        - Maybe save blacklist file to this repo for clarity and add explanation
    - Filter out blacklisted genes before HVG calculation
- hvgs: for non-batch views, make sure that sc.pp.highly_variable_genes(batch_key="Sample") is used
- need to run and create new harmony embeddings (integrated by samlpe) based on pca embeddings for the "benchmark_analysis" views and create cell type annotations based on unsupervised clustering based on pca and additionally also based on harmony embeddings
- for Batch_effect.Rmd, ensure it uses the updated preprocessing pipeline with batch-aware HVG calculation and new harmony embeddings
- check strategy to handle gongsharma dataset. it's huge, that's why the authors provided the datasets in smaller chunks of .h5ad files subsetted by gender, cmv and age
    - convert Preprocess_gongsharma.ipynb to qmd
    - cleanup Preprocess_gongsharma if necessary (still contains legacy code for other conditions that are not used in the current analysis (see datasets.json for most up-to-date list of datasets used and conditions))
    - it does not make sense to combine into one file, so keep it separate
    - check which files need to be actually pre-processed (in datasets.json, do not change this file)

# Process_data.ipynb
Requires complete overhaul
- First step: Possibly convert ipynb to quarto for better reproducibility and agentic workflow
- Major point: Process_data should be able to be run on the hpc cluster
    - Add bash script to submit jobs to cluster scheduler
        - separate scripts can be saved to src/bash/run_python_sample_embedding_methods/ folder
        - Create SLURM submission script with appropriate resource requests (memory, time, nodes)
        - Use array jobs for each dataset or dataset-method combination
- Keep Process_data.qmd at first
- datasets (including batching strategy) now defined in datasets.json
- remove redundant preprocessing step (e.g. sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor="seurat")
        sc.pp.scale(adata, max_value=10)
        sc.pp.pca(adata, n_comps=50), etc.)
    - see preprocess.py for the correct implementation
    - Make sure that every method gets the correct input (counts, or specific embedding)
    - Double-check whether this works with the current preprocess.py workflow (because scanpy overwrites adata.X at every pre-processing step!)
- Add batch effect mitigation strategies
    - MrVI has native "nuisance variable" that handles batch effect (see here: https://docs.scvi-tools.org/en/1.3.3/tutorials/notebooks/scrna/MrVI_tutorial.html)
        - Possibly just add batch_col argument to MrVI constructor: MRVI.setup_anndata(adata, sample_key="Sample", batch_key=batch_col) ?
    - scPoli will not be used for batch effect analysis
    - PILOT (and PILOT-GM-VAE) takes either the PCA embedding or the Harmony integrated space
- peak_memory was previously implemented but not running correctly in wsl2 on windows desktop
    - check online for reasons why it might not have worked correctly
    - might work if run on hpc cluster
        - because of different memory allocation and file system handling
        - if each dataset (or dataset-method combination) is run in separate instance
    - otherwise comment out or drop completely (probably cleaner and still backed up in git history)
- double-check def log_execution_time_and_memory() output format and if it conforms with the expected input in benchmark_analysis.rmd at p_exec_times (which needs to be combined also with r_exec_times)
    - Suggest plan to harmonize, if necessary

# Batch effect

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
    - datasets.json: add "columns" "batch"
         - Joanito batch was manually defined at the end of batch_effect_analysis.rmd
    - update datasets with batch column mapping
    - handle case where multiple datasets are combined (e.g. in batch_effect_analysis.rmd at "# Combined PBMC (Stephenson, GongSharma, Zhu)")

## Down-stream
- DESeq2.normalize():
    - check batch_col is correctly implemented
    - check that it does not get "Sample" as a batch column
- get_pb_deseq2(): add argument batch_col, implement batch_col

## Analysis