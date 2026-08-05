> Kilo Code indexes code structure and function signatures automatically.
> AGENTS.md focuses on domain concepts, pipeline logic, and project conventions that indexing cannot infer.

# Paper/repo review and update strategy
This repo is about: ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts
Link to paper: https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full

## Open reviewer work
- Extend batch effect analysis and benchmark analysis (more datasets + more methods).
    - Additional datasets will be added by the user (human) — not by agents.
    - Pipeline/code work, method feasibility (PILOT-GM-VAE, QOT, PULSAR) and implementation drafts: see TODO.md (Phase 3 + Phase 4).

# General rules
- Do not run pipeline scripts (e.g. .R, .py or .sh) for validation checks after implementing new code, unless the user asks for.
    - Validation of HPC pipeline scripts (e.g. .R, .py or .sh) will be run once the pipeline has been fully implemented, using a small debugging dataset (e.g. derived from the Joanito dataset)
- All HPC bash scripts must run with the working directory set to ${PROJECT_ROOT}:
  source `src/slurm_config.sh`, then `cd "${PROJECT_ROOT}"`. This is the established
  convention in every existing script — keep it for any new script (Python/R interop
  resolves repo-relative paths; see docs/ARCHITECTURE.md).


# Repo structure

## Documentation
- README.md (overview + usage), `docs/ARCHITECTURE.md` (pipeline call flow, file-role tables, HPC layout — the authority), TODO.md (task tracking). Keep documentation files up-to-date upon any changes.

## datasets.json
- This acts as ground truth for the datasets evaluated in this study
    - See datasets.json for most up-to-date list of datasets used and conditions.
- Do not change this file without asking
- The `_debug` entry (Joanito 5-sample subset, built by the Joanito step `1.2.1_prepare_joanito.R` into `${HPC_SCRATCH_DIR}/_debug/data/`) is registered here with both views. Convention: default-all script loops (`1_stage_data.sh`, `1_submit_hpc_array.sh`) skip `_*` keys unless explicitly requested via `--ds_name _debug`; `_debug.folder_name` is `null`, so staging skips it — the raw subset never lives on the NAS, only `_debug` *outputs* sync to NAS (`${NAS_TARGET_DIR}/_debug/output/`). Used for debugging `3_scrnaseq_preprocessing/`, `4_cell_type_annotation/`, `5_benchmark_methods/` (not the simple `1_stage_data/` / `2_dataset_specific_preprocessing/`). Details: see ARCHITECTURE.md.
- Files defined in datasets.json are stored on the NAS
    - The user and agents exclusively work on the user's computer, so the NAS is only accessed by the user from the computer
    - NAS dataset path from user computer (connect to NAS server first, ask user if needed): `/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets`
    - NAS dataset path from HPC (login node only): `/srv/smednas515.unige.ch/carmona_smb/DataCollections/Standardized_SingleCell_Datasets`

## data/
- `data/ARCHIVE_LEGACY_DATA/`: legacy data from previous workflow (local seurat objects in .rds files) — do not use; will most likely not work with the current pipeline.


## Pipeline Overview (Stage 1–4)
Four-stage end-to-end pipeline; file-level details live in docs/ARCHITECTURE.md.

- **Stage 1 — QC Filtering** (`notebooks/QC_filtering/`): per-dataset R Markdown notebooks.
- **Stage 2 — Preprocessing + Cell Type Annotation** (see [ARCHITECTURE.md](docs/ARCHITECTURE.md#preprocessing-pipeline) and [annotation section](docs/ARCHITECTURE.md#cell-type-annotation-pipeline-src4_cell_type_annotation)):
    - Staging (`src/1_stage_data/`) → dataset-specific steps (`src/2_dataset_specific_preprocessing/`, e.g. `1.2.1_prepare_joanito.R` builds the `_debug` subset + `seqtec` batch column; `1.4.1_create_kfoury_lowres_ct.R` creates `cells_lowres`) → preprocess array (`src/3_scrnaseq_preprocessing/`) → annotation chunks + array + merge (`src/4_cell_type_annotation/`; `3.1_submit_merge.sh <DS>` merges `annotations_chunk_*.feather` into every view h5ad and syncs to NAS).
    - Preprocessed .h5ad files are **CSR-on-disk by construction** — required for selective backed-mode per-sample reads in annotation (details in ARCHITECTURE.md).
    - **Drafts (keep, not dead code)**: `preprocess_gongsharma.qmd` (GongSharma other-subsetting conditions) and `TODO_STUMP_preprocess_sikkema.qmd` (future Sikkema Lung dataset) in `src/3_scrnaseq_preprocessing/` are intentional drafts for future implementation — do NOT delete.
- **Stage 3 — Benchmark Analysis** (see [ARCHITECTURE.md](docs/ARCHITECTURE.md#benchmark-ecoda-transformation-and-ecoda-zero-imputation-analyses)): `notebooks/benchmark_analysis.rmd`'s `run_analyses()` orchestrates 3.1 benchmark methods (R-native + Python via `.feather`), 3.2 transformation analysis (`datrans()`), 3.3 zero imputation. Pending pipeline work (Python methods on HPC, PILOT-GM-VAE, QOT/PULSAR): see TODO.md Phase 3.
- **Stage 4 — Batch Effect Analysis** (see [ARCHITECTURE.md](docs/ARCHITECTURE.md#batch-effect-analysis)): `notebooks/batch_effect_analysis.rmd`, under expansion (methods: ECODA batch-associated CT removal, Pseudobulk DESeq2+limma, MrVI, GloScope, PILOT-GM-VAE): see TODO.md Phase 4.

## R Modules for benchmark analysis (`src/5_run_benchmark_methods/` and `src/utils/`)

11 utility files loaded by `src/utils/load_all_functions.R`, plus 2 benchmark-specific files in `src/5_run_benchmark_methods/` (details in ARCHITECTURE.md).

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
- `slurm_config.sh` is the HPC config file, used by all bash scripts, containing paths to the HPC cluster and other settings.
- **Worker environment invariants** (details in ARCHITECTURE.md):
    - Python is invoked via `PYTHON_BIN` and R via `PIXI_RSCRIPT` from `slurm_config.sh` — never bare `python`/`Rscript` (worker nodes may not have scanpy/anndata); `RETICULATE_PYTHON` is also exported so R workers' reticulate always uses the pixi python (mirrors the project-root `.Rprofile`, which only applies to non-vanilla sessions).
    - Annotation (`2.1.1_process_chunk.R`) builds Seurat objects from the raw counts layer via `get_seurat_obj_from_h5ad()` (`layers["counts"]`, X fallback with warning) — NOT from log-normalized `X`; feather names derive from the chunk file (`chunk_<N>.txt` → `annotations_chunk_<N>.feather`), not `SLURM_ARRAY_TASK_ID`; scGate models load from the shared `${SCGATE_DB_PATH}` cache (`aux/scGateDB.rds`) created by `2.0_create_scgate_db.R`.
    - Annotation paths are per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (see `config_helper.R`); `SAMPLE_COLNAME="Sample"` is exported by `slurm_config.sh`; `TISSUE_TYPE`/`NORMAL_TISSUE` are auto-exported per array task from `datasets.json` by `2.1_run_worker.sh` (via jq, `module load jq/1.6` on the worker).

### HPC folder layout
- Full layout (scratch, reference atlases, NAS targets, env-var table): see [ARCHITECTURE.md](docs/ARCHITECTURE.md#hpc-folder-layout)
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/1_stage_data/1_stage_data.sh); results copied back to NAS after processing (typically in `*_submit_hpc_array.sh` upon completion)
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
