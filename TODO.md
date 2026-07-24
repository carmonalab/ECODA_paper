# General
- check README.md and update if necessary:
    - README.md is intended mainly for human readers, e.g.:
    - how to install and run the repo
    - General workflow (data processing pipeline, analysis steps, output generation)
        - Overview of data sources and preprocessing steps
        - Description of analysis methods and batch effect correction strategies (see also below on initial minimal batch effect analysis implemented in the preprint and major point to be addressed as per reviewers)
        - Expected outputs and figure generation process
        - Add this in very concise form also to AGENTS.md
- read the paper related to this repository (https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)
    - get the purpose of this analysis, hypothesis, goal and methodology
    - get key findings and understand the methodology
- based on the above, check AGENTS.md and suggest a plan how to improve it (AGENTS.md was just created with minimal information, mainly only with key info on datasets.json and info on the specific HPC cluster it will be run on)
    - it should give necessary context to an agentic coding model but not bloat it, so keep it concise
    - as you read the paper, note any domain-specific terminology that needs to be defined for the agent
        - Highlight any hard-to-understand mathematical transformations (e.g., clr, impute_zeros)
    - as kilo code uses code base indexing, the following points might not be necessary but double-check if they are actually needed or not or how minimal additional information added to AGENTS.md might improve the agent's performance, given it has a code base indexing from kilo code:
        - Ensure it explains the project folder structure and covers all functions in src/
        - Include key function signatures and expected inputs/outputs
        - Reference the modular structure and dependencies between files
        - Add brief explanations of core concepts (ECODA, CLR, batch correction strategies)
- suggest improvements, e.g:
    - Clearer folder structure and/or organization
    - Standardize file naming conventions for consistency
    - Propose best strategy on how to structure README.md, AGENTS.md and TODO.md:
        - Check online for current reccomendations and guidelines on how to structure the AGENTS.md file, e.g.:
            - Structure
            - Critical information and what it should contain
            - How to keep it concise to prevent bloat filling up the agents context window
        - Which information has to be in both or can redundancy be reduced? Redundant info might still be necessary, as README.md is for humans and might need more or other information than AGENTS.md which possibly could be more concise.
        - Some information in the AGENTS.md might need to be moved to TODO.md and vice versa. Check and suggest best habits.
    - Update README.md
    - Update AGENTS.md
    - Update TODO.md
    - Consolidate duplicate or redundant scripts into single, well-documented modules (explained in further details below)
- Iterate on the above suggestions until README.md and AGENTS.md are complete and a general repo structure was found that is clear and organized

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

Should it be done once without batch correction, and once with? -> probably more important to only do WITH batch correction.

The final analysis for batch effect correction needs to be run on the following methods and batch effect mitigation strategies:
- ECODA: remove cell types that are significantly different across batches
    - Possibly re-use process_coda_fig() with added argument batch_col = NULL?
    - Use metadata column batch_col to test for batch-associated cell types (after clr-transformation), remove them and re-calculate clr-transformed cell type composition
        - Which statistical method to use and which threshold to use for significant difference between batches? -> Print warning naming cell types and Significance. -> Possibly depends on the number of batches:
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