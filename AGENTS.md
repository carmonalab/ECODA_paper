> Kilo Code indexes code structure and function signatures automatically.
> AGENTS.md focuses on domain concepts, pipeline logic, and project conventions that indexing cannot infer.

# Paper/repo review and update strategy
This repo is about: ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts
Link to paper: https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full

## Major reviewer points to be addressed
- Extend batch effect analysis
    - More datasets have to be added
        - A requirement for this: Batch effect analysis pipeline has to be implemented in a more structured and standardized way
        - Preprocessing has to be harmonized to fit with the benchmark analysis
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are defined in the TODO.md
        - A draft strategy for specific pipeline and code implementation are defined in the TODO.md
- Extend benchmark analysis
    - More datasets have to be added
        - Benchmark analysis has mainly to be adapted to be run on the HPC cluster and to be cleaned up with minor adaptions
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are PILOT-GM-VAE (very similar to PILOT which is already implemented, trivial) and possibly PULSAR (needs to be tested if it can be run at all, as it is a foundation model, requiring substantial hardware only available on the cluster, including GPU. Input for PULSAR is Universal Cell Embedding (UCE)).
            - PILOT-GM-VAE can be added by agent
            - PULSAR needs to be tested for requirements

# General rules
- Do not run pipeline scripts (e.g. .R, .py or .sh) for validation checks after implementing new code, unless the user asks for.
    - Validation of HPC pipeline scripts (e.g. .R, .py or .sh) will be run once the pipeline has been fully implemented, using a small debugging dataset (e.g. derived from the Joanito dataset)


# Repo structure

## Documentation
- `AGENTS.md`: specific instructions, rules and information for agents
- `README.md`: simple overview of the repo and how to use it
- `docs/ARCHITECTURE.md`: overview of the full workflow, function-level documentation, and cross-language dependencies.
- `TODO.md`: Todo list
- Documentation files should be kept up-to-date upon any changes

## datasets.json
- This acts as ground truth for the datasets evaluated in this study
    - See datasets.json for most up-to-date list of datasets used and conditions.
- Do not change this file without asking
- Files defined in datasets.json are stored on the NAS
    - The user and agents exclusively work on the user's computer, so the NAS is only accessed by the user from the computer
    - The user and agents can work on the HPC but always connecting from the user's computer
    - NAS dataset path called from user computer (needs user to connect to NAS server first, ask him to connect if needed): `/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets`
    - NAS dataset path called from HPC (needs user to connect to NAS server first, ask him to connect if needed; no need to connect to additionally connect NAS server but can only be called from login node, see section "# HPC general information" below): `/srv/smednas515.unige.ch/carmona_smb/DataCollections/Standardized_SingleCell_Datasets`


## Pipeline Overview
Four-stage end-to-end pipeline:

- **Stage 1 — QC Filtering** (`notebooks/QC_filtering/`): Per-dataset R Markdown notebooks (12 cohorts). Standard scRNA-seq QC: mitochondrial genes, gene/transcript count thresholds. Produces QC-filtered input for Stage 2.
- **Stage 2 — Preprocessing**:
    - Needs to be run first: Dataset-specific preprocessing steps in `src/2_dataset_specific_preprocessing/`, submitted via the `1_submit_hpc.sh` dispatcher (parallel per-step sbatch jobs: `1.1_submit_combinedpbmc.sh` → `_create_combinedpbmc_dataset.py`, `1.2_submit_joanito_batch_col.sh` → `_create_joanito_batch_col.R`; the Joanito step computes the `seqtec` batch column and must run before the preprocess array)
    - Raw-data staging from NAS to HPC scratch: `src/1_stage_data/1_stage_data.sh` (login-node script, run before dataset-specific preprocessing and the preprocess array)
    - HPC pipelines for `src/3_scrnaseq_preprocessing/` + `src/4_cell_type_annotation/`, both run arrays across datasets (spawning one worker node per dataset):
        - `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`: bash script to run standardized scRNA-seq preprocessing Python/Scanpy pipeline on the HPC cluster:
            - Filter cells (min_genes=100) and genes (min_cells=3) + low cell-count sample removal
            - Sample/gene name standardization, sample subsetting (driven by `datasets.json`)
            - HVG selection (batch-aware and non-batch modes, multiple sizes)
            - PCA and Harmony integration
            - Unsupervised cell clustering (multiple Leiden resolutions)
        - `src/4_cell_type_annotation/`: HPC-parallelized scATOMIC + HiTME cell type annotation via SLURM array jobs
            - `src/4_cell_type_annotation/1_prepare_chunks.sh`: bash script that reads the preprocessed .h5ad files and creates multiple chunk files, one chunk file per 5 samples.
                - Each chunk `chunk_*.txt` file contains the information to run per cpu core (file and samples).
            - `src/4_cell_type_annotation/2_submit_hpc_array.sh`: bash script to run cell type annotation via scATOMIC + HiTME (HPC-parallelized, borrowed from another project, not yet polished).
                - Spawns workers, one worker per chunk file, to process samples in parallel.
- **Stage 3 — Benchmark Analysis** (`notebooks/benchmark_analysis.rmd` + `src/utils/` + `src/5_run_benchmark_methods/` + `src/5_run_benchmark_methods/run_python_sample_embedding_methods/`):
    - TODO:
        - `run_python_sample_embedding_methods/`: SLURM bash scripts to run Python sample embedding methods on the HPC cluster. To bo done:
            - Adapt and convert `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` to .py and run on the cluster by calling from a bash script
                - Rename `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` to `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`
            - Add bash scripts to `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh` and `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`
                - `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh` calls `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`
            - Needs to be checked if the following methods can be run (QOT might fail due to bugs and PULSAR might fail due to hardware requirements):
                Add methods QOT (https://github.com/PennShenLab/QOT) and PULSAR (https://github.com/snap-stanford/PULSAR) to `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`
                    - QOT is not a package but only. provided in a .py script. Dependencies need to be extracted from one of the .ipynb files they provide.
                    - PULSAR is a package, but I think it requires "Universel Cell Embeddings" (UCE, see https://github.com/snap-stanford/UCE) as input and needs to be run with HPC GPU with possibly large VRAM requirement.
    - `notebooks/benchmark_analysis.rmd`:`run_analyses()` orchestrates three sub-pipelines:
        - 3.1 runs R-native benchmark methods and incorporates output from `run_python_sample_embedding_methods/` scripts (.feather files containing embeddings or distance matrices)
        - 3.2 transformation analysis comparing ECODA transformations via `datrans()` parallel engine (in R)
        - 3.3 zero imputation analysis (4 strategies + multiLN/multiRepl)
- **Stage 4 — Batch Effect Analysis** TODO: create batch effect analysis. (`notebooks/batch_effect_analysis.rmd`): Under expansion per reviewer comments. Methods (and batch mitigation strategy per method) to implement:
    - ECODA (batch-associated CT removal)
    - Pseudobulk (DESeq2 + limma)
    - MrVI (provides native batch handling, simply add parameter and pass batch_col)
    - GloScope (Harmony space in adata.obsm["X_pca_harmony"] (created in `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`))
    - PILOT-GM-VAE (Harmony space in adata.obsm["X_pca_harmony"] (created in `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`))
- **HPC wrappers**:
    - `slurm_config.sh` provides basic configuration (e.g. paths)
    - SLURM scripts:
        - `*_submit_hpc_array.sh` submits job arrays, spawning one worker node with `*_run_worker.sh` per chunk file or dataset or (method and dataset). Planned:
            - `1_stage_data/1_stage_data.sh` (stages raw data from NAS to HPC scratch folder, login node)
            - `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (submits array to preprocess datasets in parallel)
            - `4_cell_type_annotation/1_prepare_chunks.sh` (reads preprocessed .h5ad files and creates multiple chunk files, one chunk file per 5 samples)
            - `4_cell_type_annotation/2_submit_hpc_array.sh` (submits array to process chunk files in parallel (one chunk file contains ~5 samples) to annotate cells per sample)
            - `5_run_benchmark_methods/run_python_sample_embedding_methods/`, planned to run python benchmark methods in parallel across all datasets.

## R Modules for benchmark analysis (`src/5_run_benchmark_methods/` and `src/utils/`)

11 utility files loaded by `src/utils/load_all_functions.R`, plus 2 benchmark-specific files in `src/5_run_benchmark_methods/`

# HPC general information

## IMPORTANT:
- Do not use the login node to execute any code
    - If you do, you are disturbing all other users and this is unacceptable. When this happens we will most likely kill your process.
    - The login node should only be used to compile your code and submit a Slurm job. You must even use Slurm to run your tests. The debug-cpu and debug-gpu partitions are dedicated for small tests.

## Additional important points:
- Current repo lives on a local MacOS computer
- The user and agents can work on the HPC but always connecting from the user's computer
- If you need to test HPC cluster bash scripts:
    - The HPC cluster is only available if logged in to the UNIGE network (user might work from home (needs to connect to VPN) or from the office (has access to UNIGE network)).
        - If in the UNIGE network, you can log in with `ssh [REDACTED_HOST]` (user needs to enter the password).
- Heavy scripts are run on the HPC cluster, specifically located in these folders:
    - `src/1_stage_data`
    - `src/2_dataset_specific_preprocessing`
    - `src/3_scrnaseq_preprocessing`
    - `src/4_cell_type_annotation`
    - `src/5_run_benchmark_methods/run_python_sample_embedding_methods`
- These HPC pipelines are not finished yet and need to be updated to be finalized:
    - `src/3_scrnaseq_preprocessing`
    - Not implemented yet: HPC pipeline bash scripts for `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` (needs to be adapted, e.g. converted to .py and run on the cluster by calling from a bash script)
    - `src/4_cell_type_annotation` — HPC-parallelized scATOMIC + HiTME cell type annotation.
        - renv remnants removed; R environment fully managed by pixi (`pixi run Rscript`).
        - R code extracted from bash heredoc into standalone `2.1.1_process_chunk.R`.
        - `config_helper.R` moved from `DEPRECATED_LEGACY_CODE/` (deleted) to project root (env-var based).
        - Annotation paths are per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (see `config_helper.R`); annotation feathers go to `${HPC_SCRATCH_DIR}/${DS_NAME}/output` directly — `samples/`, `ecoda/`, `plots/` dirs are no longer created. `SAMPLE_COLNAME="Sample"` is exported by `slurm_config.sh` (preprocess standardizes the sample column), and `TISSUE_TYPE`/`NORMAL_TISSUE` are auto-exported per array task from `datasets.json` by `2.1_run_worker.sh` (via jq, `module load jq/1.6` on the worker).
        - See [ARCHITECTURE.md](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-stage-2b) for full pipeline documentation.
- `slurm_config.sh` is the HPC config file, used by all bash scripts, containing paths to the HPC cluster and other settings.

### HPC folder layout
- `$HOME/scratch/ECODA_paper` (`HPC_SCRATCH_DIR`): `<DS_NAME>/data/` (staged raw inputs), `<DS_NAME>/output/` (preprocessed .h5ad per view, `chunks/` during annotation, `annotations_chunk_*.feather`, annotated .h5ad), `CombinedPBMC/data/` (combine output + rds→h5ad cache), `chunks_manifest.txt` (global chunk manifest)
- `$HOME/reference_atlases/sketched_200ct/` (`HOME_REF_DIR`); `$PROJECT_ROOT/logs`, `aux/`, `.pixi/`
- NAS: `NAS_SC_DIR` (raw source), `NAS_REF_DIR`; `NAS_TARGET_DIR` = `Projects/ECODA_paper/` with `<DS_NAME>/output/` (rsynced per-dataset from `${HPC_SCRATCH_DIR}/<DS_NAME>/output`), `benchmark/{embeddings,plots}/` and `batch_effect_analysis/{embeddings,plots}/` (targets for method .feathers + notebook plots; filled once the `5_run_benchmark_methods` decision is made — TODO)
- See [ARCHITECTURE.md](docs/ARCHITECTURE.md#hpc-folder-layout) for the full layout
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/1_stage_data/1_stage_data.sh)
- results must be copied back to NAS after processing (typically implemented in `*_submit_hpc_array.sh` scripts upon completion of the HPC jobs)
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start


# Batch effect analysis dataset info
- Most datasets are monolithic h5ad files with a batch_col, e.g.:
    - Stephenson has batch effect by batch_col "Site" (both, in terms of gene expression (major) and cell type composition (minor, just one or a few monocyte subtypes))
- A "combined PBMC" dataset was created from multiple other available datasets (included for method benchmark analysis and or batch effect analysis):
    - Combined PBMC (Stephenson, GongSharma, Zhu) (see batch_effect_analysis.rmd, see also TODO.md)


# Domain Terminology

- **ECODA** (Exploratory Compositional Data Analysis): Uses CLR-transformed cell-type proportions for cohort-level patient stratification in an unsupervised setting.
- **CLR** (Centered Log-Ratio): Transformation for compositional data: `log(x_i / geometric_mean(x))`. Requires zero-imputation beforehand. Implemented in `src/utils/math_utils.R:6`.
- **HVCs** (Highly Variable Cell Types): Cell types with highest variance across samples, selected for stratification. Implemented in `src/utils/hvcs.R`.
- **Zero imputation**: Four strategies for handling zero cell-type counts before CLR: `counts_zeros` (replace zeros with fixed count), `counts_all` (add fixed count to all), `percentage_zeros` (replace zeros with percentage of row total), `percentage_all` (add percentage to all). Implemented in `src/utils/math_utils.R:30`. Also uses `multiLN` and `multiRepl` from `zCompositions`.
- **Pseudobulk**: Aggregating single-cell counts per sample before DESeq2 normalization. Implemented in `src/utils/pseudobulk.R`.
- **Separation metrics**: Evaluate how well methods recover known biological groupings. All in `src/utils/scoring_metrics.R`:
  - **ANOSIM**: Analysis of Similarities (`calc_sep_score()`)
  - **ARI**: Adjusted Rand Index (`clust_eval()`)
  - **Silhouette**: `calc_sil()`
  - **Modularity**: `calc_modularity()` with multiple KNN variants (sqrt(n), 3, 6, 9)
  - **LISI**: Local Inverse Simpson's Index (`calc_lisi()`, :159)
- **Harmony integration**: Batch correction by integrating PCA embeddings across samples/batches. Computed in `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`.
- **`.feather` files**: Apache Arrow format for cross-language data exchange. Python methods in `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd` output distance matrices/embeddings as `.feather`; R method processors in `src/5_run_benchmark_methods/benchmark_methods_r.R` ingest them.