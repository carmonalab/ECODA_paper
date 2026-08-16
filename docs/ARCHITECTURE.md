# ECODA_paper Architecture & Call Graph

## Overview
This is a **multi-layered, top-down architecture** for exploratory compositional data analysis of scRNA-seq cohorts: it benchmarks sample-representation methods (MOFA+, scITD, GloScope, GloProp, MrVI, PILOT, scPoli, Pseudobulk, ECODA) for recovering known biological groupings, with batch effect analysis as one stage. The code is organized with clear separation between orchestration, method processing, evaluation, and utilities.


## Datasets and Analysis Definitions (Benchmark and/or Batch Effect Analysis)
```
datasets.json (read by src/utils/datasets_io.R::read_datasets_json() and
               src/utils/py/datasets_io.py::read_datasets_json())
└── per-dataset entries
    ├── "Adams" → benchmark_analysis view
    │   ├── output_file: AdamsT_2020_32832599_benchmark_analysis_ECODAprocessed.h5ad
    │   ├── label_col, cell_type_low_res, cell_type_high_res (from dataset-level columns)
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
$HOME/ECODA_paper                       # repo clone (PROJECT_ROOT): logs/, aux/, .pixi/ (HPC envs at .pixi/envs/py-cuda13/)
$HOME/scratch/ECODA_paper               # HPC_SCRATCH_DIR
├── <DS_NAME>/data/                     # staged raw inputs per dataset (1_stage_data.sh)
├── <DS_NAME>/output/                   # preprocessed .h5ad per view, chunks/, annotations (incl. annotation_tmp/ checkpoints), annotated .h5ad
├── <DS_NAME>/annotation_union/         # per-dataset union h5ad for annotation (outside synced output dir)
├── CombinedPBMC/data/                  # combine output + rds→h5ad cache (CombinedPBMC/output/ holds its preprocessed output)
├── benchmark/embeddings/               # method .feathers + per-task exec logs (Python + R benchmark arrays)
├── benchmark/results/                  # per-method/per-combo R result bundles (<ds>_<method>.rds, <ds>_<combo>.rds, <ds>_trans.rds, <ds>_zeroimp.rds)
├── benchmark/pseudobulks/              # shared DESeq2 pseudobulks (<ds>_pseudobulk_<variant>.rds, prepare_pseudobulk worker)
├── benchmark/gloscope_dists/           # raw GloScope distance cache (<ds>_gloscope_hvg<n>_pcadims<d>_dists.rds)
├── benchmark/checksums.md5             # md5 sidecar over the RDS bundles (written by benchmark_merge_sync_cleanup; verified by the notebook before readRDS)
├── benchmark_manifest_<method>.txt     # per-method benchmark manifests (1_submit_hpc_array.sh, rebuilt every run)
├── chunks_manifest_<pid>.txt            # per-submission chunk manifest (2_submit_hpc_array.sh)
├── _worker_retries/                    # per-(job,task) self-requeue counters (worker_retry.sh; harmless when stale)
└── _benchmark_watchdog/                # benchmark watchdog status files (<watchdog_id>.status; never synced)
$HOME/reference_atlases/sketched_200ct/ # HOME_REF_DIR (HiTME reference maps)

# NAS (carmona_smb; login node only)
DataCollections/Standardized_SingleCell_Datasets/   # NAS_SC_DIR — raw source datasets
DataCollections/reference_atlases/sketched_200ct/   # NAS_REF_DIR
Projects/ECODA_paper/                               # NAS_TARGET_DIR — results back-sync target
└── <DS_NAME>/output/                   # rsynced per-dataset from ${HPC_SCRATCH_DIR}/<DS_NAME>/output
benchmark/{embeddings,results,pseudobulks,gloscope_dists,plots}/   # method .feathers (embeddings/) + R result bundles (results/, pseudobulks/, gloscope_dists/) + merged execution_times.feather, all filled by the benchmark pipelines; plots/ by the notebook
batch_effect_analysis/{embeddings,plots}/           # same, for batch effect analysis
```

| Path | Env var | What goes here |
|---|---|---|
| `${HPC_SCRATCH_DIR}/<DS_NAME>/data/` | `HPC_SCRATCH_DIR` | Staged raw inputs per dataset (`1_stage_data.sh` rsyncs from `NAS_SC_DIR`) |
| `${HPC_SCRATCH_DIR}/<DS_NAME>/output/` | `HPC_SCRATCH_DIR` | Preprocessed .h5ad per view (keys `X_pca_{view}_hvg{n}`, ...); during annotation additionally `chunks/chunk_*.txt`, `annotations_chunk_*.feather`, `annotation_tmp/` per-sample checkpoints (never synced to NAS), merged annotated .h5ad |
| `${HPC_SCRATCH_DIR}/<DS_NAME>/annotation_union/` | `HPC_SCRATCH_DIR` | Per-dataset union h5ad (`union.h5ad`, dedup on `(sample, barcode)`) in the **minimal annotation layout**: `X` = raw counts (float64 CSR), `obs`+`var` only, no `layers`/`obsp`/`obsm`/`uns`/`varm`/`varp` — anndata's `read_h5ad(backed="r")` eagerly loads the whole `layers` group into RAM at open, so a full-layout union would add a multi-GB memory floor per annotation worker (Kfoury: +1.31 GB; a GongSharma-class union: 26–40 GB); chunked + annotated ONCE per dataset; OUTSIDE the synced `output/` dir so NAS-sync globs stay clean; deleted by `3_submit_merge.sh` after merging |
| `${HPC_SCRATCH_DIR}/CombinedPBMC/data/` | `HPC_SCRATCH_DIR` | Dual role: output of the CombinedPBMC combine step and rds→h5ad cache; input to the preprocess array (`CombinedPBMC/output/` holds its preprocessed output) |
| `${HPC_SCRATCH_DIR}/chunks_manifest_<pid>.txt` | `HPC_SCRATCH_DIR` | Per-submission chunk manifest (one `DS_NAME<TAB>chunk_path` line per chunk), written by `2_submit_hpc_array.sh`, read by its own array workers via the exported `CHUNKS_MANIFEST`; PID suffix keeps parallel submissions from clobbering each other |
| `${HPC_SCRATCH_DIR}/_worker_retries/` | `HPC_SCRATCH_DIR` | Per-(job,task) self-requeue counter files (`<jobid>_<taskid>.count`, `src/utils/bash/worker_retry.sh`); cleared on task success, stale leftovers harmless (ids never collide across submissions) |
| `${HPC_SCRATCH_DIR}/_benchmark_watchdog/` | `HPC_SCRATCH_DIR` | Benchmark watchdog status files (`<watchdog_id>.status` — `STATE=OK|FAIL` + `JOB_REPORT=` lines; written by `watchdog_main.sh`, read by the login tail's `benchmark_wait_watchdog`); job-id named, stale leftovers harmless, never rsync'd (outside `benchmark/`) |
| `${HOME_REF_DIR}` | `HOME_REF_DIR` | HiTME reference maps (staged by `1_prepare_chunks.sh`: NAS first via `NAS_REF_DIR`, Figshare download fallback, DOI 10.6084/m9.figshare.26310994) |
| `${PROJECT_ROOT}/logs`, `aux/`, `.pixi/` | `PROJECT_ROOT` | SLURM logs (`LOGS_DIR`), auxiliary files (gene maps), pixi envs (HPC: `.pixi/envs/py-cuda13/`) |
| `${NAS_SC_DIR}` | `NAS_SC_DIR` | Raw source datasets (`folder_name` in `datasets.json`) |
| `${NAS_REF_DIR}` | `NAS_REF_DIR` | Reference atlas source (HiTME) |
| `${NAS_TARGET_DIR}/<DS_NAME>/output/` | `NAS_TARGET_DIR` | Results back-sync target: rsynced per-dataset from `${HPC_SCRATCH_DIR}/<DS_NAME>/output`; `benchmark/{embeddings,results,pseudobulks,gloscope_dists,plots}/` (method `.feather`s + R result bundles + merged `execution_times.feather` from the benchmark pipelines; `plots/` filled by the notebook), `batch_effect_analysis/{embeddings,plots}/` (Phase 4) |


## Preprocessing Pipeline

### Standard scRNA-seq preprocessing pipeline (`src/3_scrnaseq_preprocessing/`)

The preprocessing stage is split across three `src/` folders run in sequence:

| Folder | Role |
|---|---|
| `src/1_stage_data/` | `1_stage_data.sh` — login-node script: sources `slurm_config.sh`, stages raw inputs from NAS to HPC scratch (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`) via jq over `datasets.json` + rsync. Datasets with views emit each view's `input_file_name`; view-less datasets (`"views": {}`, e.g. Zhu) fall back to the dataset-level `file_names`; `sort -u` dedups files shared across views and datasets with `folder_name: null` (CombinedPBMC) are skipped. Supports `--ds_name <DS>` single-dataset mode; default-all runs skip keys starting with `_` (e.g. the `_debug` entry) unless explicitly requested. Sends a best-effort completion email (files staged, missing-source warnings, wall-clock duration; `--ds_name`-scoped subject in single-dataset mode) via `src/utils/bash/sync_status_email.sh` — no sbatch involved, so no sacct/job reports. |
| `src/2_dataset_specific_preprocessing/` | `1_submit_hpc.sh` dispatcher submits all per-step sbatch jobs in this folder in parallel (`1.1_submit_gongsharma.sh`, `1.2_submit_combinedpbmc.sh`, `1.3_submit_joanito.sh`, `1.4_submit_kfoury_lowres_ct.sh`), passing per-step `--partition`/`--output`/`--error`/`--mail-user` on the sbatch command line (SLURM directives cannot expand env vars, so the flags are given on the sbatch command line; partition comes from the `${SLURM_PARTITION}` comma list `shared-cpu,shared-gpu,private-carmona-gpu` in `slurm_config.sh` — stages 2–4 jobs may land on any of the three, whichever frees resources first; `--mail-user="${USER_EMAIL}"` overrides the step scripts' bare `#SBATCH --mail-type=END,FAIL`, which would otherwise email the cluster default user), waits for all, reports per-job state via `sacct` and exits non-zero on any failure. Supports `--sync-only <id1,id2,...>` resume mode (skip submission, re-run the wait + sacct gate + summary for previously submitted jobs; no sync happens here — the preprocess array syncs afterwards). Sends ONE best-effort completion email at the end ("COMPLETED (N steps)" / "FAILED (k/N)" + per-step JobName/State/Elapsed/ExitCode lines from sacct; no NAS sync by design). All steps run in parallel EXCEPT the one explicit dependency: the GongSharma cap step (`1.1`) is submitted first (and skipped by name in the parallel loop, so it is never double-submitted) and the CombinedPBMC step (`1.2`) is gated behind it via `--dependency=afterok` (1.1 overwrites the staged SoundLife h5ads in place, which 1.2 reads in backed mode — fail-closed, so a failed cap never submits CombinedPBMC). Run after staging, before the preprocess array. The Joanito step (`1.3.1_prepare_joanito.R`) additionally builds the `_debug` 5-sample subset into `${HPC_SCRATCH_DIR}/_debug/data/`. |
| `src/3_scrnaseq_preprocessing/` | `1_submit_hpc_array.sh` (array submit + monitor + rsync back to NAS), `1.1_run_worker.sh`, `1.1.1_preprocess.py`. |

#### Files

| File | Role |
|---|---|
| `1_submit_hpc_array.sh` | Thin bash wrapper: sources `slurm_config.sh`, submits the preprocess array (`1.1_run_worker.sh`; `--partition="${SLURM_PARTITION}"` passed on the sbatch command line — SLURM directives do not expand env vars), monitors completion, verifies via `sacct` that every task state is `COMPLETED` (fail-closed: aborts without syncing on any non-COMPLETED state or empty sacct output), rsyncs results back to NAS (login node; single-dataset mode syncs only the requested dataset). Supports `--ds_name <DS>` single-dataset mode (1-task array at the dataset's position in the sorted jq key list, so the worker's `sed -n ${SLURM_ARRAY_TASK_ID}p` mapping needs no change) and `--force` (recomputes existing outputs; forwarded to the workers via the `FORCE_PREPROCESS` env var); default-all runs skip `_*` keys — this relies on `_` sorting after all uppercase-initial dataset keys, making the non-`_` keys a prefix of the sorted keys. Supports `--sync-only <job-id>` resume mode (skip submission, re-run the monitor/gate/sync tail for a previously submitted job — reconnect after an SSH drop; repeat the original `--ds_name`; unknown/purged ids fail closed without syncing). The tail sends a best-effort sync-status email via `src/utils/bash/sync_status_email.sh` (`notify_sync_status`; login-node mail CLI, skipped silently if absent): "preprocess synced to NAS" on success, "NOT synced — reason" (with a mapped per-task report on gate failure, replacing the raw sacct dump) before every fail-closed `exit 1`. Both emails include a per-task report (task → dataset from the FULL sorted jq key list — `_` keys included, so single-dataset `_debug` tasks map correctly — state, elapsed, exit code) + array wall time (sacct `Elapsed`; `n/a` when sacct is unavailable/purged). Raw-data staging no longer lives here — see `src/1_stage_data/1_stage_data.sh`. |
| `1.1_run_worker.sh` | `#SBATCH` worker: sources `slurm_config.sh`, resolves its dataset from `SLURM_ARRAY_TASK_ID`, calls `1.1.1_preprocess.py` (via `${PYTHON_BIN}`) with `--config_path/--input_dir/--output_dir/--ds_name` (+ `--force` when `FORCE_PREPROCESS=1`), under the worker self-healing wrapper (`worker_retry.sh`: transient-signature grep on the Slurm `.err` + counter-capped self-requeue). No R library staging (python env too large), no thread pinning (multi-threaded by design). Retry is safe: `write_h5ad` is atomic (tmp + `os.replace`), so a mid-write crash leaves no target file and the "Already processed" skip never reuses partial outputs. |
| `1.1.1_preprocess.py` | Standardized preprocessing: filtering, gene/sample name standardization, sample subsetting (from `datasets.json`; subsetting runs on original values BEFORE the sample column is standardized, and errors out if the subset is empty), batch-aware HVG selection, PCA, Harmony integration, Leiden clustering. Writes one .h5ad per view. Promotes an `X=None` input with a `counts` layer to `X` (anndataR-written files, e.g. the `_debug` h5ad). `--force` bypasses the "Already processed" skip (needed for debug re-runs). **CSR-on-disk by construction**: `base_preprocessing()` forces `tocsr()` on `X` and `layers["counts"]` unconditionally (not only dense inputs); scanpy ops afterwards (normalize_total, log1p, HVG selection, subsetting) preserve CSR and `write_h5ad()` preserves the in-memory format, so the written files are always CSR. |

- **View output keys** (unified per-view pipeline, `1.1.1_preprocess.py`): both view types get the same treatment; only the HVG `batch_key` differs (benchmark views: standardized sample column from `SAMPLE_COLNAME`; batch-effect views: the dataset's `batch_col`).
  - `X_pca_{view}_hvg{n}` — stored for **every** HVG size (`benchmark_analysis`: 3000/2000/1000; `batch_effect_analysis`: 2000).
  - Harmony + unsupervised clustering (neighbors + Leiden) run **only at the 2000-HVG pass** (`CLUSTER_N_HVG`), on both embeddings:
    - `X_pca_harmony_{view}_hvg2000`, `neighbors_{view}_hvg2000`, `leiden_res_{r}_{view}_hvg2000`
    - `neighbors_{view}_hvg2000_harmony`, `leiden_res_{r}_{view}_hvg2000_harmony` (on the harmony-corrected embedding)
  - Example (benchmark view): `X_pca_benchmark_analysis_hvg3000`, `X_pca_benchmark_analysis_hvg2000`, `X_pca_benchmark_analysis_hvg1000`, plus `X_pca_harmony_benchmark_analysis_hvg2000`, `leiden_res_0.1_benchmark_analysis_hvg2000`, `leiden_res_0.1_benchmark_analysis_hvg2000_harmony`, ... (batch views analogously with `batch_effect_analysis_hvg2000`).

- **NAS ↔ Scratch data flow**: Raw-data staging from NAS to scratch happens in `src/1_stage_data/1_stage_data.sh` (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`); Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS (per-dataset `${NAS_TARGET_DIR}/${DS_NAME}/output`).
- **Working-directory convention**: every HPC bash script sources `src/slurm_config.sh` and then `cd "${PROJECT_ROOT}"` (all existing scripts do this; keep it for any new script). `src/utils/py/preprocess_utils.py` pins the embedded R working directory to `${PROJECT_ROOT}` at import time (its rpy2 interop sources `src/utils/seurat_utils.R` + attaches only Seurat/anndataR — a deliberately minimal import, not the 20-package `load_all_functions.R`; preprocessing only needs `readRDS → create_clean_seuratv5_object → write_h5ad`, mirroring the annotation worker's lighter pattern), so Python callers are CWD-independent — the `cd` remains belt-and-braces and required by convention. **sbatch-submitted scripts must not resolve `slurm_config.sh` from `BASH_SOURCE`**: Slurm copies submitted scripts to `/var/spool/slurmd/job<id>/slurm_script`, so `BASH_SOURCE` points at the spool dir. The 9 sbatch workers (`1.1_submit_gongsharma.sh`, `1.2_submit_combinedpbmc.sh`, `1.3_submit_joanito.sh`, `1.4_submit_kfoury_lowres_ct.sh`, `1.1_run_worker.sh`, `4_cell_type_annotation/2.1_run_worker.sh`, `5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`, `5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh`, `5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh`) recover `SCRIPT_DIR` from the job record via `scontrol show job` (`Command=` field) when `SLURM_JOB_ID` is set, with the `BASH_SOURCE` fallback for login-node execution (`srun`-based scripts such as `1_prepare_chunks.sh` run without `SLURM_JOB_ID` and are unaffected).
- **Dataset-specific preprocessing steps** (all submitted via the `1_submit_hpc.sh` dispatcher; all run in parallel EXCEPT the GongSharma cap step `1.1`, which the dispatcher submits first (and skips by name in the parallel loop, so it is never double-submitted) and gates the CombinedPBMC step `1.2` behind via `--dependency=afterok` — see the 1.1 bullet):
  - **CombinedPBMC combine** (`1.2_submit_combinedpbmc.sh` → `1.2.1_create_combinedpbmc_dataset.py`): HPC-capable. With `HPC_SCRATCH_DIR` set it defaults to `--layout per-dataset` (sources at `${HPC_SCRATCH_DIR}/<ds>/data`, output `${HPC_SCRATCH_DIR}/CombinedPBMC/data`); locally it defaults to the flat `PROJECT_ROOT/data` layout. Parallel in-job: the three sources load concurrently in 3 fork workers (`ProcessPoolExecutor`); each worker writes a trimmed per-source intermediate to `output_dir/_intermediates/<source>_subset.h5ad` (overwritten on rerun), and the main process reloads those small intermediates for the common-gene intersection (with `<5000` union fallback), Sample prefixing, `sc.concat(join="outer")` and write. rpy2 is imported lazily inside the two R workers (Stephenson, Zhu) only — the parent never initializes R before forking. GongSharma is read in backed mode (`sc.read_h5ad(..., backed='r')`): non-counts layers dropped on disk, only the HDF5 chunks covering the rng(123)-picked ~15 samples materialized via `to_memory()` (full in-memory load fallback on failure). The GongSharma component reads the **capped** staged files (≤5000 cells/sample, seed-42 deterministic subset — see the 1.1 bullet); the CombinedPBMC dataset is rebuilt automatically on every `1_submit_hpc.sh` run, so its GongSharma content updates with the cap. 16 cpus / 256G sbatch baseline (rds→h5ad cache `{stem}_raw.h5ad` makes reruns faster). Script is CWD-independent; still submitted via the `1_submit_hpc.sh` dispatcher.
  - **Joanito prepare** (`1.3_submit_joanito.sh` → `1.3.1_prepare_joanito.R`): computes the `seqtec` batch column in place (idempotent; single source of truth for the batch mapping) from `${HPC_SCRATCH_DIR}/Joanito/data/JoaI_2022_35773407_Nofilt_whole.rds` (path from `datasets.json` via `DATASETS_JSON_FILE`), then builds the `_debug` 5-sample subset from the same in-memory object (see "Debug subset build" below). 32G mem baseline — the whole .rds is read in a single process. On re-runs, if `seqtec` already exists in the staged `.rds`, both the recompute and the in-place `saveRDS` are skipped (fast path). Required before the preprocess array (the `batch_effect_analysis` view uses it as `batch_col`).
  - **Kfoury low-res cell types** (`1.4_submit_kfoury_lowres_ct.sh` → `1.4.1_create_kfoury_lowres_ct.R`): ported from `Preprocess_datasets.Rmd` — collapses the author `cells` annotations into `cells_lowres` (Tcells/NKcells/Bcells/MoMac/DCcells, other labels kept) in place (idempotent) from `${HPC_SCRATCH_DIR}/Kfoury/data/Kfoury_2021_34719426.rds`. 32G mem baseline. Required before the preprocess array (`columns.cell_type_low_res = "cells_lowres"` for Kfoury).
  - **GongSharma per-sample cap** (`1.1_submit_gongsharma.sh` → `1.1.1_subset_gongsharma.py`): caps each sample (`columns.sample` = `specimen.specimenGuid`) at 5000 cells in the two staged SoundLife h5ads, **in place** (atomic temp-file write + `os.replace`), at `${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data/` — fixes the preprocess-array OOM (the worker held both files + `sc.concat` = 2.92M cells and densified `sc.pp.scale` matrices per HVG pass ≈ 70 GB for 3000 top genes, exceeding the 128G worker). Sampling mirrors the historical `downsample_by_group` (git 3a4711e, `src/py/preprocess_gongsharma.qmd`): a seed-42 `RandomState` per file, groups in pandas `unique()` order of first appearance, `RandomState.choice` per group, applied as a boolean keep-mask (order-preserving; identical selected cell set) → 531,291 + 365,000 = **896,291 cells / 180 samples / max 5000/sample**, matching the historical NAS artifact `Gongsharma_cmv_young_males.h5ad` (which is the validation target for the sampling). Pure Python (`${PYTHON_BIN}`; scanpy/anndata/numpy only, no R interop); backed `to_memory()` materializes only the needed HDF5 chunks (peak memory ≈ one capped matrix). Idempotent (re-capping an already-capped file keeps every cell). **Re-staging caveat**: `1_stage_data.sh` restores the NAS originals — re-run `1_submit_hpc.sh` before the preprocess array after any re-stage. **Ordering**: the dispatcher submits this step first (and skips it by name in the parallel loop, so it is never double-submitted) and gates the CombinedPBMC step (1.2) behind it via `--dependency=afterok` (1.2 reads the same staged files in backed mode; fail-closed).
  - **Debug subset build** (inside `1.3.1_prepare_joanito.R`, HPC scratch): NOT a separate step — the Joanito step derives it from the in-memory object after the `seqtec` computation (and in-place save, on first run). Seeded (seed 321) selection of 5 samples covering `(sample.origin × seqtec × Site)` combos (candidates with ≥500 cells), 500 cells/sample, minimal obs cols (incl. `seqtec`, `Site`, sample/patient cols), rebuilt as a clean Seurat v5 object and written via `write_h5ad()` (anndataR: X=None + `layers["counts"]`, handled by the preprocess X-promotion) to `${HPC_SCRATCH_DIR}/_debug/data/JoaI_2022_35773407_debug_5samples.h5ad` (stem must match `datasets.json`). Raw subset only ever lives in `_debug/data/` on scratch — never copy it into `_debug/output/` (the merge script treats every `output/*.h5ad` except `*_raw.h5ad` as a view). `_debug.folder_name` is `null`, so staging skips it; `_debug` *outputs* sync to NAS as usual (`${NAS_TARGET_DIR}/_debug/output/`).
  - Steps must run **after** `1_stage_data.sh` and **before** `1_submit_hpc_array.sh` (staging skips `folder_name: null` datasets, and the preprocess array task reads the combined file). Processing scripts follow the decimal depth convention (`N.N.N_<action>.<ext>`, mirroring `1.1.1_preprocess.py`/`2.1.1_process_chunk.R`).


### Cell Type Annotation Pipeline (`src/4_cell_type_annotation/`)

Cell type annotation runs as a **separate HPC-parallelized pipeline** on SLURM, independent of the benchmark analysis call graph above. It takes monolithic `.h5ad` files from preprocessing and annotates them cell-by-cell using **scATOMIC + HiTME** in parallel array jobs.

#### Workflow

```
[ Monolithic h5ad Files on NAS ]
               │
               ▼ (1_prepare_chunks.sh → 1.1_prepare_chunks.py)
   [ per-dataset union.h5ad (annotation_union/, dedup (sample, barcode)) ]
               │
               ▼
      [ chunk_1.txt ]  [ chunk_2.txt ]  ...  [ chunk_N.txt ]  (2 samples each)
               │               │
               ▼ (SLURM Array) ▼ (SLURM Array)
          [ Worker 1 ]     [ Worker 2 ]  ...              (2.1_run_worker.sh)
               │               │
               ▼               ▼
      [ chunk_1.feather ] ... [ chunk_N.feather ]           (annotations on scratch)
               │               │
               ▼  (3_submit_merge.sh → 3.1_merge_annotations.py, per view h5ad) ▼
      [ annotated .h5ad on NAS ]   (union.h5ad + chunks deleted after merge)
```

Annotations run ONCE per dataset (on the union); the merge step joins the same
feather set into EVERY view h5ad of the dataset on the `(sample, barcode)` key.
The union is always written as a minimal h5ad (raw counts in `X`, no
`layers`/`obsm`/`obsp`) — including for single-view datasets, whose
preprocessed views keep their log-norm `X` + `layers["counts"]` layout and are
never chunked directly.

#### Files

| File | Role |
|---|---|
| `1_prepare_chunks.sh` | Thin bash wrapper: sources `slurm_config.sh`, stages ref maps NAS → `$HOME` (skips if all 4 present and MD5-verified; Figshare download fallback via DOI 10.6084/m9.figshare.26310994) (raw data is staged by `src/1_stage_data/1_stage_data.sh`), then iterates over all datasets in `datasets.json` (or a single dataset passed as 2nd positional arg) calling `1.1_prepare_chunks.py` (`${PYTHON_BIN}`) via `srun` (4G, 30min, `--partition="${SLURM_PARTITION}"` at submit time) per dataset. **Preprocessing-completeness guard** (before chunking): every expected view output from `datasets.json` (`views.*` entries with `input_file_name` + `output_file_name`, mirroring `1.1.1_preprocess.py`'s skip semantics — all views present after the `_raw.h5ad` exclusions) must exist in `${HPC_SCRATCH_DIR}/${DS_NAME}/output/`; a still-running preprocess array can have written only some of a multi-view dataset's views (Stephenson, `_debug`), and chunking on a partial view set would mark the dataset "already annotated" and silently stay incomplete forever — datasets with missing view files get a loud WARNING listing them ("preprocessing incomplete ... re-run this script after the preprocess array finishes"), land in a `SKIPPED_INCOMPLETE` summary bucket, exit 0, and stay "not annotated" so they are picked up on re-run. Datasets with no expected views (e.g. Zhu, `"views": {}`) keep the plain "any preprocessed `.h5ad`" skip check as before. Datasets without preprocessed `.h5ad` input are skipped with a warning. **Already-annotated datasets are skipped** unless `--force` (accepted in any position): done = either (branch 1) `output/chunks/` holds ≥1 `chunk_*.txt` and every `chunk_N.txt` has its matching `annotations_chunk_N.feather` in `output/`, or (branch 2, the post-merge state) `output/chunks/` and `annotation_union/` are both absent but ≥1 feather remains — partial feather coverage falls through to the rebuild; `--force` restores the old unconditional rebuild + (production) feather deletion. Skip applies in both `production` and `test` modes (test-mode leftovers with complete feathers are a valid done state; annotation values are identical regardless of chunk grouping — `--force` is the escape hatch). **Clean-entry check** before any skip: `1.1_prepare_chunks.py --check-clean` (same short srun shape) reads each view's obs column names (h5py-only) and the feather column sets (pyarrow schema-only) and flags a dataset NOT clean when any Tier-1/Tier-2 legacy-matching obs column is absent from this run's feathers — flagged datasets (`FLAGGED_LEGACY` summary bucket, exit 2) are rebuilt and re-annotated instead of skipped (worker wipe + merge tiered drop then scrub them); clean datasets are skipped as before; check errors land in `FAILED_DATASETS`. Supports `test` mode (`./1_prepare_chunks.sh test <DS>` → 1 sample/chunk). |
| `1.1_prepare_chunks.py` | Native Python (anndata + h5py, no reticulate): builds ONE chunk set per DATASET. Always constructs the per-dataset **minimal annotation union** at `${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad` (dedup on `(sample, barcode)`): `X` = raw counts (float64 CSR — the same dtype as the views' `layers["counts"]`), `obs`+`var` only, no `layers`/`obsp`/`obsm` (anndata's `read_h5ad(backed="r")` eagerly loads the entire `layers` group into RAM at open, so a full-layout union would add a multi-GB memory floor per worker — Kfoury: +1.31 GB; a GongSharma-class union: 26–40 GB). Single-view datasets and multi-view datasets whose largest view already equals the union (e.g. Stephenson, benchmark ⊂ batch-effect view) get a streaming h5py copy of that view (h5py-only key-set reads: obs `_index` + sample column with categorical decoding + counts `indptr` row-count guard — no anndata open, so no eager `layers` materialization at the 4G srun; no in-memory concat); only truly partial-overlap views fall back to an in-memory `sc.concat` + dedup (warning printed; needs a raised srun `--mem`). Every written union is self-verified backed-mode (lazy `X`, nnz == source counts nnz, no layers/obsp/obsm, sample column present). Reads the union in backed mode (skipping `*_raw.h5ad` rds→h5ad conversion caches), extracts unique sample IDs from `sample_col` (env var `SAMPLE_COLNAME`), groups them into chunks of 2 (or 1 in test mode), writes `chunk_N.txt` files (1st line = union h5ad path, subsequent lines = sample IDs). **Defensive fail-closed completeness check**: after the "no h5ad files" exit, reads `DATASETS_JSON_FILE` from env and requires every expected view output (`views.*` with `input_file_name` + `output_file_name`) to be present in `output/` — a missing one exits with a CRITICAL error listing the files (dataset lands in `FAILED_DATASETS`), catching bypasses/drift of the bash-level guard; if `DATASETS_JSON_FILE` is unset it warns and skips the check (never on HPC, where `slurm_config.sh` exports it). **`--check-clean` mode** (no chunk/union writes; exit 0 = clean, 2 = legacy found, 1 = error): for each view h5ad (excluding `*_raw.h5ad`) reads the obs column names (h5py-only, `read_obs_colnames_h5py` — mirrors `read_obs_keys_h5py`'s fail-closed dataframe-encoding guard) and the annotation feather column sets (pyarrow schema-only `pa.ipc.open_file`); a dataset is clean iff every Tier-1/Tier-2 legacy-matching obs column is also a column of its own feathers; per-view output prints matched names/counts, the Tier-1 (previous-annotation-like) vs Tier-2-only (possibly author metadata) split, and the legacy subset. In production mode only, deletes stale `annotations_chunk_*.feather` from the output dir AFTER all chunks were generated successfully, so reruns cannot merge leftover annotations from an earlier chunk numbering (test mode leaves production feathers untouched). |
| `2.0_create_scgate_db.R` | One-time download of the scGate model DB (`get_scGateDB(..., force_update=TRUE)`), persisted to `${SCGATE_DB_PATH}` (`aux/scGateDB.rds`). Run from `2_submit_hpc_array.sh` before the annotation array so workers do not download in parallel. Idempotent (exits early if the cache exists). |
| `2_submit_hpc_array.sh` | Creates the scGate model DB cache once via `srun` (`2.0_create_scgate_db.R` → `${SCGATE_DB_PATH}`; `--partition="${SLURM_PARTITION}"`; failure is non-fatal, workers download + persist themselves), builds a **per-submission chunk manifest** (`${HPC_SCRATCH_DIR}/chunks_manifest_<pid>.txt`, one tab-separated `DS_NAME<TAB>chunk_path` line per chunk across all datasets or a single dataset passed as positional arg; PID suffix = benchmark-submitter pattern — the manifest is a shared mutable resource and a parallel second submission (e.g. another dataset) must not truncate it under a running array, since workers read `${CHUNKS_MANIFEST}` at task start and each sbatch submission propagates its OWN path through the environment), submits a single SLURM array job (`--array=1-TOTAL_CHUNKS`, `MAX_NUM_CHUNKS_PARALLEL` concurrency, `--partition="${SLURM_PARTITION}"`), monitors for completion, verifies via `sacct` that EVERY task state is `COMPLETED` (fail-closed: aborts without syncing on any non-COMPLETED state or empty sacct output), then syncs results back to NAS via rsync (per-dataset `${NAS_TARGET_DIR}/${DS_NAME}/output`; in single-dataset mode only that dataset is synced). Workers skip chunks whose annotation feather already exists (see `2.1_run_worker.sh`) — re-runs of the array on an annotated-but-not-yet-merged dataset do no redundant annotation; the chunk manifest stays unfiltered, so the `3_submit_merge.sh` coverage gate (chunk files vs feather count) is unchanged. The sync-tail rsync excludes `annotations_chunk_*.feather`, `chunks/` and `annotation_tmp/` (intermediates never consumed from NAS — the merge reads local scratch only; previously cleaned up post-sync by `3_submit_merge.sh`, which no longer fires for already-merged datasets). Supports `--sync-only <job-id>` resume mode (skip the scGate cache, the chunk-manifest build and the submission; re-run the monitor/gate/sync tail for a previously submitted job). The tail sends a best-effort sync-status email via `src/utils/bash/sync_status_email.sh` (success / "NOT synced — reason" with a per-dataset aggregated report on gate failure, replacing the raw sacct dump; skipped silently if no mail CLI). Because chunk arrays can be thousands of tasks, the report aggregates per dataset (chunks COMPLETED/FAILED + summed task elapsed, task → dataset from the submission's `${CHUNKS_MANIFEST}` line number; manifest missing/unknown — e.g. `--sync-only` — → task-id-only fallback) and bounds the failed-chunk task-id list to 20 + "+ N more", plus the array wall time. |
| `2.1_run_worker.sh` | `#SBATCH` worker (32G, 4 cpus, 2h; partition from `${SLURM_PARTITION}` at submit time). Sources `slurm_config.sh`, gets jq from it (guarded `module load GCCcore/12.2.0` + `jq/1.6` — module loads do not propagate through sbatch, hence the central location), reads its `DS_NAME`/`CHUNK_FILE` from the global manifest line matching `SLURM_ARRAY_TASK_ID` (`sed -n`), auto-exports per-dataset `TISSUE_TYPE`/`NORMAL_TISSUE` from `datasets.json`, checks the chunk file exists, then **skips with exit 0** (task counts as COMPLETED in sacct) when the chunk's annotation feather already exists — `${HPC_SCRATCH_DIR}/${DS_NAME}/output/annotations_chunk_<N>.feather` for `chunk_<N>.txt` (chunk number derived from the chunk file's basename, the same `chunk_N.txt → annotations_chunk_N.feather` mapping used by the merge; safe against stale feathers because `1.1_prepare_chunks.py` deletes them on every production chunk rebuild). After the feather skip (which stays fast), sources `src/utils/bash/worker_retry.sh`, pins BLAS/OMP threads to 1 (`export_worker_thread_env` — annotation-only, so CPU time ≈ wall time), and runs `2.1.1_process_chunk.R` via `${PIXI_RSCRIPT}` under a unified retry wrapper: on non-zero exit it greps the task's Slurm `.err` for transient signatures (stale BeeGFS client-cache views, missing files/imports) and self-requeues via `scontrol requeue` (counter-capped, `WORKER_MAX_RETRIES`, default 3); success clears the counter. R packages are read directly from the env library (no per-task staging since 2026-08-13); the retry wrapper covers residual stale-view flakes. |
| `2.1.1_process_chunk.R` | Core annotation logic: reads the chunk's h5ad file (the per-dataset union) **in backed mode** (`read_h5ad(..., backed="r")`; `obs` is metadata-only — read once via `py_to_r(adata$obs)` and reused for every sample), warns (not stops) if the on-disk X format is not `csr` (`adata$X$format`; anndata only overrides selective row-indexing for backed CSR — CSC would fall back to a full in-memory read per subset, i.e. silent OOM), **wipes ALL pre-existing obs columns except the sample column after the sample-column check** (`obs[, sample_col, drop=FALSE]`, `message()` naming the dropped columns) — annotation starts from a minimal Seurat object whose `meta.data` holds only `Sample`, so legacy annotation columns (`scGate_*`, `functional.cluster`, `*_UCell`, `layer_1..6`, ...) can never leak into the annotation object (they would silently skip scATOMIC's `layer_1` guard or confuse `Run.HiTME`); no pattern lists needed. Then iterates per sample, extracts sample-level Seurat objects from the raw counts (`get_seurat_obj_from_h5ad()`, sourced from `src/utils/seurat_utils.R`; `layer_keys`/`counts_layer` computed ONCE per chunk outside the loop — the layers group never changes between samples): the annotation union carries counts in `X` by design (no `layers` group; `counts_layer="X"` is the designed primary path, announced with a `message()`), while preprocessed view files keep `layers["counts"]` (`counts_layer="counts"`). Runs **scATOMIC** (up to 5 attempts with budget-aware timeouts) then **HiTME** (up to 5 attempts with budget-aware timeouts). scGate models load from the shared `${SCGATE_DB_PATH}` cache (download + persist fallback). **Wall-clock policy**: per-attempt timeouts are capped at 1800 s and only taken when the remaining wall time (from `SLURM_TIME_LIMIT` × 60, i.e. minutes→seconds, minus `proc.time()[3]` elapsed — `SLURM_JOB_START_TIME` is not set by Slurm; the under-count by staging time is the safe direction) exceeds the attempt by a 300 s margin; `R.utils::withTimeout` is best-effort only (it fires at R evaluation points, not inside blocking reticulate python calls), so the real budget is enforced at R level and Slurm's wall limit remains the backstop. **Per-sample checkpoints**: each sample is written atomically (tmp + `file.rename`) to `output/annotation_tmp/chunk_<N>/sample_<NN>.feather` as it completes; a sample-level `tryCatch` records failures and continues; if ANY sample failed the script exits 1 WITHOUT the final feather (intermediates persist → a re-run resumes only the failed samples), otherwise all intermediates are re-read in index order, `rbind`ed, and the final `annotations_chunk_<N>.feather` is written atomically, then `annotation_tmp/` is removed. The `annotation_tmp/` subdirectory never matches the merge/coverage glob `annotations_chunk_*.feather`, and `1.1_prepare_chunks.py` deletes it on every production chunk rebuild. Writes per-chunk `.feather` named after the chunk file (`chunk_1.txt` → `annotations_chunk_1.feather`), so reruns are stable regardless of array task renumbering. Dual annotation: scATOMIC provides layer-1..6 predictions + confidence; HiTME provides layer1/2/3 UCell signatures + scGate/ProjecTILs refinement. Only the fresh annotation output columns (`ANNOT_OUTPUT_COLS` — the sole columns the pipeline emits; a fresh-output extraction set, NOT a keep-whitelist) are kept. |
| `3_submit_merge.sh` | Per-dataset merge wrapper (`./3_submit_merge.sh <DS_NAME> [--force]`): **skips with exit 0** (no srun, no NAS sync, no sync-status email) when the dataset is already merged — `output/chunks/` and `annotation_union/` absent with ≥1 `annotations_chunk_*.feather` still in `output/` (this script is the only step that deletes those dirs; the feathers are kept deliberately for the coverage gate on re-runs); `--force` falls through to the full merge, still behind the coverage gate. Otherwise: fails if no `annotations_chunk_*.feather` exist or if the feather count does not match the chunk-file count (`${OUTPUT_DIR}/chunks/chunk_*.txt` — race-free: chunk files are per-dataset and only deleted by this script after a successful merge, whereas the shared `chunks_manifest.txt` is clobbered by parallel `2_submit_hpc_array.sh` runs and caused a false "0 chunks submitted" mismatch; guards against merging a partially failed annotation array; NOT bypassed by `--force`), verifies NAS reachability BEFORE any destructive step, loops over the dataset's view h5ads in `${HPC_SCRATCH_DIR}/${DS_NAME}/output`, runs `3.1_merge_annotations.py` per view via `srun` (64G baseline — big views like Stephenson's batch-effect view may need more; overridable via `MERGE_MEM` env var, `--partition` at submit time), then deletes the stale `annotation_union/` dir and `output/chunks/` (both are regenerated by the next `1_prepare_chunks.sh` run; the chunk files pin the deleted union path and must not be reused), removes the stale `annotations_chunk_*.feather` and `chunks/` from NAS after the sync (they accumulate otherwise — `2_submit_hpc_array.sh` synced them upstream; local feathers are kept for the coverage gate), and rsyncs the annotated h5ads to `${NAS_TARGET_DIR}/${DS_NAME}/output/` (skipping `*.tmp` atomic-write leftovers and `annotation_tmp/` — checkpoints must never reach NAS). Re-runs are safe (idempotent merge). Sends a best-effort sync-status email on success and before the coverage-gate/NAS-unreachable/merge-failure exits (`notify_sync_status`; no `--sync-only` — re-running the merge is already safe). **Default-all mode** (bare `./3_submit_merge.sh [--force]`, no positional): loops over all `datasets.json` keys (`_debug` included, matching the annotation default-all convention; a fresh `_debug` is WARNING-skipped like any no-feather dataset), skips already-merged datasets, WARNING-skips no-view/no-feather datasets, checks NAS once up front (fail fast, one email + exit 1), collects failures and sends ONE summary email (processed/failed lists + per-dataset merge durations measured with local `date +%s` around each `merge_one_ds` call — merges are login-node-driven `srun` calls, not a tracked array, so no sacct), exiting 1 if any dataset failed; single-dataset mode is unchanged (per-event emails, each success email carrying per-view merge durations + total). In default-all `--force`, already-merged datasets fail the coverage gate per dataset (no chunk files — `EXPECTED_CHUNKS=0`) — documented, fail-closed. |
| `3.1_merge_annotations.py` | Reads all `annotations_chunk_*.feather` files (argparse CLI: `--h5ad-path` required, `--annot-dir` defaults to the .h5ad's parent, `--output-path` defaults to overwriting the input), joins them into the input `.h5ad`'s `obs` on a `(sample, barcode)` composite key (barcodes repeat across samples and views; duplicate keys are dropped with a warning; `SAMPLE_COLNAME` read from env, default "Sample", on BOTH sides — the feather annotations column and the h5ad obs sample column — keeping the join key consistent), exits non-zero if no feathers exist, the sample column is missing, the feathers contain no annotation output columns (`ANNOT_OUTPUT_COLS`, mirror of the worker — annotation produced no results), or the merge matches 0 cells (key drift). **Tiered pre-existing-column drop** (backstop for the worker's obs wipe; `ANNOT_OUTPUT_COLS` is NOT a keep-whitelist — pre-existing versions are always dropped): Tier 1 unconditional regex patterns (`^scGate`, `^functional\.cluster`, `_UCell$`, `^scATOMIC`, `^layer_?\d`) + Tier 2 exact names (`S.Score`, `G2M.Score`, `Phase`, `classification_confidence`, `cellCycle.G1S_UCell`, `cellCycle.G2M_UCell`) dropped with a loud warning flagging them as "possibly author metadata (exact-name match)". Idempotent re-merge (existing drops before the join). **Post-merge invariant (fail loudly)**: every final obs column matching the tier patterns must be present in this run's feather columns, else `sys.exit(1)` (a missed drop pattern = stale data leaked in). Writes `adata.uns["ecoda_annotation_version"] = "1"` (provenance marker), then the annotated `.h5ad` atomically (temp file + `os.replace`). |
| `config_helper.R` | (Project root) Builds path config from env vars exported by `slurm_config.sh` (`DS_NAME`, `HPC_SCRATCH_DIR`). All dataset-specific paths are per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (`path_data`, `path_output`); `path_ref` is a home resource (`HOME_REF_DIR`). Annotation feathers go to `${HPC_SCRATCH_DIR}/${DS_NAME}/output` directly (`path_output`) — the old `samples/`, `ecoda/`, `plots/` dirs are no longer created. Missing env vars raise `stop()` errors (no silent fallbacks). Called by `2.1.1_process_chunk.R`. |

#### Key design details

- **scATOMIC + HiTME dual annotation**: scATOMIC provides hierarchical cell-type predictions (layer_1..6) with confidence scores; HiTME annotates using scGate models + ProjecTILs reference maps, producing layer1/2/3 labels. Both are run on each sample independently (i.e. independent from the rest of the dataset).
- **Retry loops**: Both annotation methods have up to 5 retry attempts with dynamic timeouts (max(60s, n_cells/10000 × 600s), per-attempt cap 1800 s, gated on remaining wall time) to handle HPC node variability.
- **Worker self-healing** (`src/utils/bash/worker_retry.sh`, sourced by all 5 array workers after `slurm_config.sh`): on a non-zero task exit the worker greps its Slurm `.err` file for transient-failure signatures (`TRANSIENT_REQEX`: `No such file or directory`, `cannot open file`, `package or namespace load failed`, `No module named`, ... — the stale BeeGFS client-cache failure mode) and self-requeues via `scontrol requeue <array_job>_<task>` from the RUNNING task, capped by a per-(job,task) counter under `${HPC_SCRATCH_DIR}/_worker_retries/` (`WORKER_MAX_RETRIES`, default 3; cleared on success). R workers read packages **directly from the env library** — the per-task `stage_env_rlib` staging (2026-08-13 removal) had been copying the whole 1.3 GB library per task, which amplified the metadata storm it was meant to fix; direct-env execution is exactly the proven `WORKER_STAGE_R_LIB=0` behavior (Python workers have always run unstaged), and the transient retry still covers residual stale-view flakes. The annotation worker alone also pins BLAS/OMP threads to 1 (`export_worker_thread_env`) so CPU time ≈ wall time; benchmark workers do NOT (hardware-pinned for runtime comparability) and preprocess does NOT (intentionally multi-threaded). Retry safety: annotation/preprocess are safe (final feathers and `write_h5ad` are atomic — partial files never pass the skip checks); the python benchmark worker deletes its per-(method,dataset) feathers on the requeue path because pyarrow `to_feather` is non-atomic and its combo skip-check would otherwise reuse a partial file; R benchmark writes are all `save_rds_atomic` (tmp+rename).
- **NAS ↔ Scratch data flow**: Raw-data staging from NAS to scratch happens only in `src/1_stage_data/1_stage_data.sh` (per-dataset dirs `${HPC_SCRATCH_DIR}/${DS_NAME}/data`); cell type annotation consumes the preprocessed output of that pipeline (`${HPC_SCRATCH_DIR}/${DS_NAME}/output` per dataset, matching `config_helper.R`). `1_prepare_chunks.sh` stages reference maps; `2_submit_hpc_array.sh` creates the scGate model DB cache (gene standardization moved into `1.1.1_preprocess.py`; the `GENE_REF_FILE` staging block was removed). The gene reference is committed to the repo at `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`, consumed by `src/utils/py/gene_utils.py` (see README "Reference data" for provenance). Worker nodes only access scratch. After array completes, login node rsyncs results back to NAS. Because the sync tails are login-bound (NAS is login-node-only), an SSH drop never kills a running pipeline: on reconnect, re-run the submitter with `--sync-only <job-id>` (same other flags) to re-check and sync; the tails also send best-effort sync-status emails (`src/utils/bash/sync_status_email.sh`) so the "synced to NAS?" answer reaches you even if the terminal session is gone — since the job-duration work, those emails include a per-task/per-dataset (or per-method) completion breakdown + job durations (sacct `Elapsed`, array wall time).
- **Environment propagation**: `slurm_config.sh` `export`s all core vars (`PROJECT_ROOT`, `DATASETS_JSON_FILE`, `HPC_SCRATCH_DIR`, `SAMPLE_COLNAME`, ref/gene files, `PYTHON_BIN`, `PIXI_RSCRIPT`, `RETICULATE_PYTHON`, `SCGATE_DB_PATH`, ...) so they reach R via `Sys.getenv()` and Python via `os.environ` through both `srun` (`1_prepare_chunks.sh`, `2_submit_hpc_array.sh`) and `sbatch` (`2_submit_hpc_array.sh`). On HPC, `PYTHON_BIN`/`PIXI_RSCRIPT`/`RETICULATE_PYTHON` all point at the `py-cuda13` pixi env (`.pixi/envs/py-cuda13/`); `PIXI_RSCRIPT` expands to `pixi run --as-is -e py-cuda13 Rscript --vanilla` (word-splits; `--as-is` = `--no-install --frozen`: workers never update the lockfile or install/repair the env at runtime — refreshes only via the guarded `setup_env_sbatch.sh`/`refresh_env.sh`; `--frozen` alone would still install the prefix when it mismatches the lockfile). `slurm_config.sh` also prepends `.pixi/envs/py-cuda13/bin` to `PATH` so rpy2 (imported at module level by `src/utils/py/preprocess_utils.py`) finds R/Rscript when workers invoke `${PYTHON_BIN}` directly, and exports `LD_LIBRARY_PATH` with the env lib dir (`.pixi/envs/py-cuda13/lib`) FIRST so R's `dyn.load()` resolves package `.so` dependencies against the conda toolchain they were built with instead of module (GCCcore) or node-system libs — pixi's activation sets this for `PIXI_RSCRIPT` workers automatically, but bare-`${PYTHON_BIN}` rpy2 workers would otherwise resolve against module/node libs and fail to attach source-built packages on node images with older libstdc++/GLIBCXX (preprocessing array 4294806); the export sits AFTER the module loads so the env lib dir wins. `RETICULATE_PYTHON` pins R's reticulate to the pixi python (reticulate's own discovery may otherwise pick a stray `~/.virtualenvs/r-reticulate` on the worker); the project-root `.Rprofile` mirrors this but only applies to non-vanilla R sessions (`PIXI_RSCRIPT` uses `--vanilla`) and points at `.pixi/envs/default` on macOS only (`py-cuda13` is linux-64-target-scoped). Bash arrays do not propagate through `sbatch`, so workers derive `DS_NAME` from `datasets.json` (via jq); `slurm_config.sh` now also performs the guarded module loads for jq centrally (toolchain + version pinned — `GCCcore/12.2.0` + `jq/1.6`, `|| true` — because Lmod's hierarchical tree hides jq/1.6 behind its GCCcore prerequisite and module loads do not propagate through sbatch either); `2.1_run_worker.sh` auto-exports per-task `TISSUE_TYPE`/`NORMAL_TISSUE` from `datasets.json`.
- **Counts input**: scATOMIC/HiTME receive the raw counts, never the log-normalized data. For preprocessed view files the counts come from `layers["counts"]` (vaulted by `base_preprocessing`); the annotation union carries the counts in `X` by design (minimal layout, no `layers` group — see `1.1_prepare_chunks.py`), so the worker's `X` fallback is the designed primary path for union files and announces itself with a `message()` rather than a warning.
- **Backed per-sample reads are selective**: the annotation union is CSR-on-disk with counts in `X` (written by `1.1_prepare_chunks.py`), so `get_seurat_obj_from_h5ad()`'s per-sample row subset (`adata[cell_indices]`) reads only the selected rows' segments (`backed_csr_matrix` selective indexing) — the minimal layout's lazy `X` makes the open cost ~10 MB instead of the multi-GB eager `layers` load; `obs` is metadata-only and never triggers matrix I/O. On a CSC-on-disk file the same subset would materialize the full matrix in memory per sample (anndata has no selective row-indexing override for CSC) — `2.1.1_process_chunk.R` warns on non-CSR input.
- **Output format**: Per-chunk `.feather` files (Apache Arrow, cross-language), named after the chunk file (`annotations_chunk_<N>.feather`) → merged into original `.h5ad` by `3.1_merge_annotations.py` on a `(sample, barcode)` composite key.
- **Bounded `cutoff.scATOMIC::em()` worker monkey-patch**: `2.1.1_process_chunk.R` (right after `library(cutoff.scATOMIC)`) replaces `em` in the `cutoff.scATOMIC` namespace with a bounded version (default `t = 1e-6` — flexmix's tolerance, empirically identical to the bit-stable fit on real Pelka scores — `max_iter = 200`, `message()` on cap). Upstream bug: `em()`'s default `t = 1e-64` is below machine epsilon and the loop has no iteration cap; `automatic_threshold()`'s `runif` jitter makes `mean(lambda)` never bit-stabilize on ill-conditioned score mixtures → infinite EM loop (each iteration a full Nelder-Mead fit), swallowing `withTimeout` interrupts until the 2h wall limit kills the task (observed on Pelka `C171_TA`, 10,709 cells). The patch is applied via `unlockBinding` + `assignInNamespace` (required: `automatic_threshold` calls `cutoff.scATOMIC::em` via `::`) with `environment()` set to the package namespace so internal helpers (`hash`, `startval`, `mLL`) resolve. See abelson-lab/scATOMIC issue <NUMBER>.
- **Column-aligned chunk assembly**: the final `annotations_chunk_<N>.feather` is assembled from the per-sample checkpoints with NA-filling of missing columns (plain `rbind` crashed on heterogeneous checkpoint schemas — "numbers of columns of arguments do not match" — when a sample's scATOMIC/HiTME produced no annotation columns, e.g. CombinedPBMC chunk_28 / Zhang chunk_13 resumes).

#### Usage

```bash
# 1. Prepare chunks (stages data + generates chunk files on the per-dataset union;
#    already-annotated datasets are skipped unless --force is given)
./1_prepare_chunks.sh                     # production, all datasets: 2 samples/chunk
./1_prepare_chunks.sh test                # test, all datasets: 1 sample/chunk
./1_prepare_chunks.sh production Stephenson   # production, single dataset
./1_prepare_chunks.sh production Stephenson --force   # force rebuild of an annotated dataset

# 2. Submit SLURM array (auto-creates scGate DB cache if missing, monitors + syncs to NAS after completion)
./2_submit_hpc_array.sh                   # all datasets with chunks
./2_submit_hpc_array.sh Stephenson        # single dataset

# 3. Merge annotations into every view h5ad of the dataset + sync to NAS
#    (already-merged datasets skip; --force re-merges through the coverage gate;
#    bare invocation = all datasets, no-feather datasets WARNING-skip, one summary email)
./3_submit_merge.sh                     # all datasets (incl. _debug)
./3_submit_merge.sh Stephenson          # single dataset
./3_submit_merge.sh Stephenson --force  # force re-merge
```

#### Test mode

`1_prepare_chunks.sh test <DS>` sets chunk_size = 1 (vs 5 for production). This
means each chunk contains only 1 sample, producing more but smaller array jobs.
Useful for quick validation together with the Joanito 5-sample `_debug`
dataset (see datasets.json + TODO.md).


---

## Benchmark, ECODA Transformation and ECODA Zero Imputation Analyses

Note: the Layers 1–5 call flow below documents the notebook-based pipeline
and is now SUPERSEDED for the R side: ALL R benchmark methods — heavy
(GloScope, MOFA, Pseudobulk, scITD) AND composition-based (ECODA variants,
GloProp, EPIC deconv, Avg_PCA_embedding, Freq_highres, via the
`composition` method) — plus the transformation/zero-imputation analyses run
on HPC (see "R benchmark methods on HPC" and "Transformation and
zero-imputation analyses on HPC" below); the notebook reads ZERO h5ad files
(bundles + `<ds>_metadata.rds` via `load_hpc_benchmark_results()`, python
method feathers via the `process_*_fig` ingest functions). Phase 3.1 (Python
methods) is already HPC-based — see "Python benchmark methods on HPC" below.

### Python benchmark methods on HPC (`src/5_run_benchmark_methods/run_python_sample_embedding_methods/`)

The Python benchmark methods (MrVI, scPoli, PILOT-GM-VAE on GPU; PILOT, QOT
on CPU) run as SLURM arrays on the preprocessed benchmark view h5ad (inputs
under `${HPC_SCRATCH_DIR}/<DS_NAME>/output`, per `datasets.json` view output
files). The legacy notebook `1.2_benchmark_methods_py.qmd` is kept as an
archive (do NOT delete) — `1.1.1_benchmark_methods_py.py` preserves its
feather naming, method-string format and data layout exactly (R ingest
`process_mrvi_fig`/`process_scpoli_fig`/`process_pilot_fig`/
`process_qot_fig`/`process_pilotgm_fig`, `constants.R` label map and
notebook recodes depend on them). QOT runs the vendored `qot_utils_re.py`
(PennShenLab/QOT @ `28cd529880c1`, two hotfixes — traceability in
`docs/qot_hotfixes.md`); PILOT-GM-VAE runs the `pilotgm` PyPI package
(CostaLab/PILOT-GM-VAE).

#### Workflow

```
1_submit_hpc_array.sh (login; per method: manifest + sbatch array;
  GPU: shared-gpu/H200/8 cores, CPU: shared-cpu/EPYC-7742/16 cores — pinned
  for runtime comparability, see slurm_config.sh benchmark vars; an explicit
  --partition <P> override drops the constraint pin)
   ├─ watchdog_main.sh (compute node, one per method array: terminal wait +
   │    per-task gate + OOM escalation, status file
   │    ${HPC_SCRATCH_DIR}/_benchmark_watchdog/<watchdog_id>.status — survives
   │    SSH drops of the login tail)
   │     └─ (OOM only) re-submits the OOM'd datasets with doubled --mem,
   │        clamped to BENCHMARK_MEM_MAX, via its own resubmit closure
   └─ 1.1_run_worker.sh (worker; METHOD + benchmark_manifest_<method>_<pid>.txt -> DS_NAME)
        └─ 1.1.1_benchmark_methods_py.py (CLI; one (method, dataset) task)
             input:  ${HPC_SCRATCH_DIR}/<DS_NAME>/output/<view output_file> h5ad
                     (obsm X_pca_{view}_hvg{n} for PILOT/QOT/PILOT-GM-VAE;
                      var["hvg_rank"] + layers["counts"] raw counts + Sample
                      obs col + ct annotation cols for MrVI/scPoli)
             output: ${HPC_SCRATCH_DIR}/benchmark/embeddings/
                     {ds}_hvg{n}[_lowres|_highres]_{mrvi_dists,scpoli_dims<d>_embs,pilot_dists,qot_dists,pilotgm_dists}.feather
                     execution_times_<METHOD>_<DS>.feather  (per-task exec log)
   └─ login tail: wait per method on the WATCHDOG id (benchmark_wait_watchdog:
      reads STATE from the status file; OK -> merge its JOB_REPORT lines; FAIL
      -> email + fail closed) ->
      NAS reachability check ->
      merge (--no-cleanup --labels <methods> --datasets <DS> [--existing-log NAS log]; login)
        -> execution_times.feather (keeps NAS log continuity across partial runs)
   └─ rsync -> ${NAS_TARGET_DIR}/benchmark/  (embeddings/ + execution_times.feather)
   └─ only then: delete this run's per-task logs (label x dataset cross product,
      plus a legacy execution_times_task_* sweep; never the merged log)
```

#### Files

| File | Role |
|---|---|
| `1_submit_hpc_array.sh` | Login-node submitter: `[--ds_name <DS>] [--methods mrvi,scpoli,pilot,qot,pilotgm] [--partition <P>] [--force] [--sync-only <id1,id2,...>]`. Resolves benchmark datasets from `datasets.json` (jq: `use_for_benchmark == true` AND a `benchmark_analysis` view; skips `_*` keys unless `--ds_name` is given). Submits ONE array per method with hardware pinned on the **sbatch command line** (SLURM directives do not expand env vars): mrvi/scpoli/pilotgm → `${SLURM_PARTITION_BENCHMARK_GPU}` + `--gpus=${BENCHMARK_GPU_COUNT}` + `--constraint=${BENCHMARK_GPU_CONSTRAINT}` + `${BENCHMARK_GPU_CPUS_PER_TASK}`; pilot/qot → `${SLURM_PARTITION_BENCHMARK_CPU}` + `--constraint=${BENCHMARK_CPU_CONSTRAINT}` + `${BENCHMARK_CPU_CPUS_PER_TASK}`; all with `--mem=${BENCHMARK_MEM}`. `--partition <P>` overrides the per-method partition AND drops the method's `--constraint` pin (debug-only partitions — `_debug` smoke tests on `private-carmona-gpu` or an ad-hoc `shared-gpu` override; the pinned constraint would never match non-pinned hardware and jobs would hang PENDING; `--gpus`/`--cpus-per-task`/`--mem` are kept). Writes `${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_$$.txt` (one dataset per line, rebuilt every run; PID suffix prevents an overlapping second submission from clobbering queued arrays' manifests) and exports `METHOD` + `BENCHMARK_MANIFEST` + `FORCE_BENCHMARK` per submission (sbatch propagates only the exported environment). Array throttle: `${BENCHMARK_GPU_ARRAY_THROTTLE}` (4 = the 4 H200s on gpu006) for GPU methods, `${MAX_NUM_CHUNKS_PARALLEL}` for CPU. Logs to `${LOGS_DIR}/5_benchmark_<method>_%A_%a.log/.err`. Each method array gets a **compute-node watchdog** (`benchmark_submit_watchdog` → `watchdog_main.sh`, 1 cpu/2G/`WATCHDOG_TIME_LIMIT` default 12h, per-method partition WITHOUT the pinned constraint class; logs `${LOGS_DIR}/5_benchmark_watchdog_<method>_<id>.log/.err`) that owns the terminal wait (shared `squeue` poll, exact-id match via `-o %A`, 60s interval + bounded `sacct` poll-until-terminal) + the **OOM-aware fail-closed per-task gate** and writes its status file (`${HPC_SCRATCH_DIR}/_benchmark_watchdog/<id>.status`) — an SSH drop of this login tail can no longer interrupt an escalation chain (the watchdog survives on a compute node; recovery after a tail death is `--sync-only <watchdog_id>` via the status-file branch, or an idempotent re-run). The login tail waits per watchdog id via `benchmark_wait_watchdog` (terminal wait → status-file grace ≤2 min → `STATE=OK`: merge the `JOB_REPORT=` lines into the final email's job-durations block; `STATE=FAIL`: email "NOT synced — watchdog failed" with the report + exit 1; watchdog non-COMPLETED or COMPLETED without a status file: fail closed with its sacct `State,ExitCode` + a pointer to its logs), then **NAS reachability check FIRST** (fail before doing any destructive merge work), then merges exec logs via `1.1.2_merge_execution_times.py` on the login node (`${PYTHON_BIN}`) scoped to this run's methods x datasets (`--labels`/`--datasets`) with `--no-cleanup` and `--existing-log` pointing at the NAS log (keeps full-log continuity across `--ds_name` partial runs), then rsyncs `${HPC_SCRATCH_DIR}/benchmark/` → `${NAS_TARGET_DIR}/benchmark/`, and ONLY after the rsync succeeds deletes this run's per-task logs (scoped to the run's label x dataset cross product, plus a separate legacy `execution_times_task_*` sweep). `--sync-only <id1,id2,...>` (comma-separated, one id per submitted array) skips the submission loop and re-runs the gate + merge/sync/cleanup tail for the given ids. |
| `1.1_run_worker.sh` | `#SBATCH` worker (defaults overridden at submit time): 12h (`shared-*` partition max), 8 cores, 128G. Sources `slurm_config.sh`, `cd "${PROJECT_ROOT}"`; requires `METHOD` + `BENCHMARK_MANIFEST` env (error if unset); reads its `DS_NAME` from the manifest via `sed -n "${SLURM_ARRAY_TASK_ID}p"`. Output dir `${HPC_SCRATCH_DIR}/benchmark/embeddings` (created), per-task log `execution_times_<METHOD>_<DS_NAME>.feather` (deterministic name: each concurrent array has a distinct METHOD, and re-runs overwrite the same file). Forwards `--force` when `FORCE_BENCHMARK=1`. Calls `1.1.1_benchmark_methods_py.py` via `${PYTHON_BIN}` under the worker self-healing wrapper (`worker_retry.sh`: transient-signature grep on the Slurm `.err` + counter-capped self-requeue; the requeue path deletes this task's `*${DS_NAME}*.feather` outputs because pyarrow `to_feather` is non-atomic and the combo skip-check would reuse a partial file). No R library staging (python env too large), no thread pinning. |
| `1.1.1_benchmark_methods_py.py` | CLI py script (replaces the qmd logic): `--config_path --ds_name --view benchmark_analysis --method {mrvi,scpoli,pilot,qot,pilotgm} --input_dir --output_dir --log_file --hvg {1000,2000,3000} --force --device {auto,cpu,cuda}`. Reads datasets.json via `src/utils/py.datasets_io.read_datasets_json` (repo-root import, like `1.1.1_preprocess.py`). **Input resolution**: view output h5ad via `sc.read_h5ad` + `.to_memory()`; ct columns from datasets.json (`cell_type_low_res`/`cell_type_high_res`); sample column `obs["Sample"]`; HVG subset via stored `var["hvg_rank"]` (`top_n_hvg_genes`, sorted ranks — no recomputation), **subset FIRST, then point X at the raw counts layer** (`layers["counts"]`; X is log-normalized — MrVI keeps the sparse counts, scPoli densifies only the small HVG subset via `layers["counts"].toarray().astype("float32", copy=False)`; warning + log-normalized X fallback if the layer is missing). **Combos** (legacy rules; skipped if the output feather exists unless `--force`): MrVI on the lowres resolution for every HVG size → `{ds}_hvg{n}_mrvi_dists.feather`; scPoli highres (n=2000 → dims 2,3,5,10,15; n=1000/3000 → dim 15) + lowres (n=2000 → dim 15) → `{ds}_hvg{n}{_lowres|_highres}_scpoli_dims{d}_embs.feather`; PILOT/QOT/PILOT-GM-VAE highres every HVG size + lowres n=2000 → `{ds}_hvg{n}{_lowres|_highres}_{pilot,qot,pilotgm}_dists.feather`. h5ad loaded ONCE per task. Combos run **defaults-first** (MrVI_hvg2000, scPoli_hvg2000_dims15_highres, PILOT_hvg2000_highres, QOT_hvg2000_highres, PILOT-GM-VAE_hvg2000_highres — the `constants.R` `method_label_map_main` rows) because `ru_maxrss` peak RSS is monotonic within a process: earlier combos report the least bloated `mem_GB`; non-default combos keep their relative order (stable sort; no-op if `--hvg` excludes the default size). **Method bodies** (qmd semantics preserved): MrVI `setup_anndata(sample_key="Sample")` + `train(max_epochs=50, accelerator=<device>)` (default `auto` → GPU on shared-gpu nodes; qmd hard-coded CPU) + `get_local_sample_distances(keep_cell=False, groupby="dummy_col", batch_size=32)`; scPoli with qmd `early_stopping_kwargs` + `recon_loss="nb"`, `get_conditional_embeddings()` → `DataFrame` indexed by sample names; PILOT consumes the **preprocessed obsm** (`emb_matrix=f"X_pca_{view}_hvg{n}"`, qmd recomputed PCA) → `adata.uns["EMD_df"]`. QOT runs the **vendored** `qot_utils_re.py` (`Run_QOT`, `num_components_list=[1]`, `min_samples_for_gmm=0`, `qot_method="cosine"`; two hotfixes — see `docs/qot_hotfixes.md`) → `adata.uns["QOT_Distance"]`, wrapper passes a distinct temp obs col (`_bench_prog`) as `progession` to dodge the upstream duplicate-key rename bug. PILOT-GM-VAE (`pilotgm` PyPI) runs `train_gmvae` (50 epochs, `num_classes = max(2, n_unique_ct)`, `use_cuda = device=="cuda" or (device=="auto" and torch.cuda.is_available())`) + `gmmvae_wasserstein_distance` (`sample_col="Sample"`, `status="_bench_status"` temp col, `wass_dis=True`) inside a node-local tempdir (`train_gmvae` hardcodes `./trained_models/<ds>/`; keeps weights off the repo and off the NAS sync) → `adata.uns["EMD_df"]`; the obsm entry is a plain ndarray during training (torch 2.x rejects DataFrames) and a named-columns DataFrame for the distance step (extract needs `.columns`); `get_pilotgm()` shim keeps the pilotgm package dir on `sys.path` (its module-level `from networks.Networks import *` — loky workers inherit `sys.path` at spawn, removing it breaks their `import pilotgm`). Feathers are plain `DataFrame.to_feather()` with the pandas index (sample names) kept — index becomes the last feather column, matching R's `column_to_rownames(ncol)`. **Exec-time/memory logging**: per combo, `time.time()` around the method body only (excludes h5ad loading) + `ru_maxrss` (Linux KB / macOS bytes → GB); exact legacy method strings `MrVI_hvg{n}` / `scPoli_hvg{n}_dims{d}{res}` / `PILOT_hvg{n}{res}` / `QOT_hvg{n}{res}` / `PILOT-GM-VAE_hvg{n}{res}`; row `{dataset, method, time_secs, mem_GB}` appended (read-modify-write) to the per-task `--log_file` (one process per file). `scvi.settings.seed = 0`; prints `scvi.__version__` + `torch.cuda.is_available()`; `jax` dropped (MrVI is pytorch now). |
| `1.1.2_merge_execution_times.py` | Login-node merge: `--output_dir` (default `benchmark/embeddings`) + `--labels <names> --datasets <ds>` (scope the task-log glob `execution_times_<label>_<ds>.feather` to this run's method/analysis × dataset cross product — stale logs from previous failed runs never leak in) + `--existing-log <path>` (e.g. the NAS `execution_times.feather`: this run's rows win via `drop_duplicates(keep="last")`, untouched historical rows are preserved, so partial `--ds_name` runs extend the full log instead of overwriting it) + `--cleanup`/`--no-cleanup` (default on; the submit script passes `--no-cleanup` and deletes the logs itself after the rsync). Dedup on (dataset, method) keep=last matches qmd overwrite-on-rerun semantics. The notebook's unified exec-time section reads this feather directly and unions bundle-derived rows (mem_GB from the bundle when present) — regenerated fresh every knit, no cache. The extra `mem_GB` column carries each worker's process peak RSS (python `ru_maxrss`, R `VmHWM`); legacy R rows computed before 2026-08-16 keep NA (optional sacct-`MaxRSS` backfill, see TODO.md). |

#### Key design details

- **Hardware pinned for runtime comparability**: same GPU model (H200,
  `nvidia_h200_nvl`, gpu[005-006]) for all GPU tasks, same CPU model (EPYC-7742,
  cpu[001-052]) + core count + RAM for all CPU tasks — see the benchmark block
  in `slurm_config.sh` (all env-overridable; formal partition strategy is
  TODO Phase 3.5). An explicit `--partition <P>` override drops the
  `--constraint` pin (explicit partition choice = user accepts non-pinned
  hardware) — used for `_debug` smoke tests on the private node
  (`${SLURM_PARTITION_PRIVATE}`), which is deliberately excluded from the
  pinned benchmark partitions (different CPU model would flaw runtime
  comparisons; H100 ≠ H200 constraint).
- **Feather layout is contractual**: R ingest (`process_*_fig` + `constants.R`
  label map + notebook recodes) requires the qmd's exact file names, method
  strings and index-as-last-column layout — do not change.
- **PILOT input change**: the qmd recomputed normalize/log1p/scale/PCA from
  scratch; the HPC version consumes the preprocessed obsm
  `X_pca_{view}_hvg{n}` (TODO 3.1 "consume preprocessed obsm" requirement).
- **Counts input for MrVI/scPoli**: both models need the raw counts (not the
  log-normalized X) — HVG-subset first (`var["hvg_rank"]`), then point X at
  `layers["counts"]`; scPoli densifies only the small subset
  (`toarray().astype("float32", copy=False)`) to keep memory bounded.
- **Per-task exec logs**: one deterministic feather per (method, dataset) —
  `execution_times_<METHOD>_<DS>.feather`; no job id needed since the
  concurrent arrays each have a distinct method/analysis, and re-runs
  overwrite the same file. Merged on the login node after the sacct gate,
  scoped to this run's (method × dataset) cross product via `--labels` +
  `--datasets` (stale logs from failed runs never leak in), with
  `--existing-log` preserving the NAS log across partial runs; per-task logs
  are deleted only after the NAS rsync succeeded (scoped to this run's label x
  dataset cross product, so an overlapping submission's logs survive; legacy
  task-job-id logs get a separate always-safe sweep).
  Concurrent submissions of the two benchmark submitters are unsupported:
  the shared merged feather + full-dir rsync already interleave (the
  cleanup is scoped to each run's own logs, so overlapping runs do not
  delete each other's per-task logs, but the merged feather/rsync race
  remains).
- **`--device auto`**: MrVI uses the allocated GPU on any GPU node — the
  `shared-gpu` nodes, or the private node for `_debug` runs via
  `--partition ${SLURM_PARTITION_PRIVATE}` (`torch.cuda.is_available()` is
  printed for verification).
- **Skip semantics**: combos whose output feather exists are skipped unless
  `--force` (qmd skip logic); ct columns missing from datasets.json skip the
  affected resolution with a warning (MrVI without `cell_type_low_res` skips
  entirely).

#### Usage

```bash
# 1. (prereq) preprocess + annotate the benchmark datasets (stages 2-4)

# 2. Run Python benchmark methods (all benchmark datasets; arrays monitored
#    + sacct-gated + exec logs merged + NAS-synced by the script)
./1_submit_hpc_array.sh                                  # mrvi, scpoli, pilot, qot, pilotgm
./1_submit_hpc_array.sh --ds_name _debug --methods qot,pilotgm  # debug single
./1_submit_hpc_array.sh --methods mrvi,scpoli --force    # recompute existing feathers
./1_submit_hpc_array.sh --partition debug-cpu            # override partitions
./1_submit_hpc_array.sh --ds_name _debug --methods pilot \
  --partition private-carmona-gpu                        # _debug on the private node
                                                         # (constraint pin dropped)
```

#### Test mode

The Phase 3.1 debug smoke test runs `./1_submit_hpc_array.sh --ds_name _debug
--methods mrvi` (then scpoli, pilot, qot, pilotgm) against the preprocessed
`_debug` benchmark h5ad; check the feathers + exec log; R-ingest
compatibility is exercised by Phase 3.4. Heavy Python deps
(scvi-tools/scarches/pilotpy/pilotgm/torch/phate) live in the `py-cuda13` pixi
env (`PYTHON_BIN`).

### R benchmark methods on HPC (`src/5_run_benchmark_methods/run_r_sample_embedding_methods/`)

All R benchmark methods (heavy: GloScope, MOFA, Pseudobulk, scITD;
composition-based: `composition` — the ECODA_* family, GloProp, EPIC deconv,
Avg_PCA_embedding, Freq_highres) run as SLURM
arrays on the preprocessed benchmark view h5ad (inputs under
`${HPC_SCRATCH_DIR}/<DS_NAME>/output`, per `datasets.json` view output files),
mirroring the Python pipeline (per-method manifests, `sed -n
${SLURM_ARRAY_TASK_ID}p`, fail-closed sacct gates, shared exec-log schema +
merge script + NAS `benchmark/` target). Workers run R via `${PIXI_RSCRIPT}`
(exported by `slurm_config.sh`), source `src/utils/imports_worker_core.R` +
`src/utils/load_worker_functions.R` + `benchmark_hpc_utils.R` (MOFA2/scITD
attached conditionally on `--method`, EPIC/GloScope for `composition`; see
the R-module notes in AGENTS.md),
and read the preprocessed h5ad through reticulate
(obs-only backed read; raw counts from `layers["counts"]` + stored obsm
`X_pca_benchmark_analysis_hvg{n}` embeddings + `var["hvg_rank"]` gene ranks —
no PCA/FindVariableFeatures recomputation, consistent with how
PILOT/MrVI/scPoli consume preprocessing; the `composition` worker is
obs-only — backed obs + hvg2000 obsm + precomputed hvg2000 pseudobulk, no
Seurat materialization). All methods are pinned to the CPU
benchmark class (`${SLURM_PARTITION_BENCHMARK_CPU}` EPYC-7742,
`${BENCHMARK_CPU_CPUS_PER_TASK}` cpus, `${BENCHMARK_MEM}`) — the same class as
PILOT, so cross-method runtime comparisons stay valid; an explicit
`--partition <P>` override drops the constraint pin.

#### Workflow

```
1_submit_hpc_array.sh (login; per method: manifest + sbatch array on the
  pinned CPU class shared-cpu/EPYC-7742/16 cores/128G; --partition override
  drops the constraint pin; mofa/pseudobulk/composition auto-prepend
  prepare_pseudobulk)
   ├─ prepare_pseudobulk array FIRST + its own soft-gate WATCHDOG
   │    (watchdog_main.sh; artifact gate: all PB_VARIANT_NAMES variants
   │    present per dataset in benchmark/pseudobulks/ -> STATE=OK without the
   │    strict task-state gate; the strict OOM-aware gate applies under
   │    --force or when variant files are missing — a stale-node prep failure
   │    with the variants already on disk must not block the method arrays,
   │    while an OOM'd prep task auto-escalates its memory)
   │    └─ 1.1_run_worker.sh -> 1.1.1_prepare_pseudobulk.R
   │         input:  h5ad (counts + var["hvg_rank"] only; no embeddings)
   │         output: ${HPC_SCRATCH_DIR}/benchmark/pseudobulks/<ds>_pseudobulk_<variant>.rds
   │                 (list(pb, time_secs, mem_GB) per variant; exec-log rows
   │                 prepare_pseudobulk_<variant>)
   └─ gloscope/mofa/pseudobulk/scitd/composition arrays (after the prep
        watchdog) — each with its own strict watchdog (terminal wait +
        per-task gate + OOM escalation: OUT_OF_MEMORY tasks re-submit only
        their datasets with doubled --mem, clamped to BENCHMARK_MEM_MAX;
        non-OOM failures fail closed; status file
        ${HPC_SCRATCH_DIR}/_benchmark_watchdog/<id>.status; survives SSH
        drops of the login tail)
        └─ 1.1_run_worker.sh -> 1.1.1_run_benchmark_methods_r.R
             input:  h5ad -> Seurat (counts + obsm pca_benchmark_analysis_hvg{n}
                     reductions + hvg_rank VariableFeatures); mofa/pseudobulk
                     reuse the precomputed pseudobulks/ (on-the-fly fallback);
                     mofa skips the counts/embeddings materialization unless
                     the fallback triggers (metadata/labels from obs);
                     composition is obs-only (backed obs + hvg2000 obsm +
                     precomputed hvg2000 pb; EPIC/GloScope attached)
             output: ${HPC_SCRATCH_DIR}/benchmark/results/
                     <ds>_<combo>.rds (per-combo bundles; combo names are
                     method-prefixed, so no method infix; each bundle carries
                     exec_time + mem_GB = peak_rss_gb() VmHWM)
                     <ds>_<method>.rds (named list of bundles)
                     <ds>_metadata.rds (composition worker: labels, n_cells,
                     n_samples, cells_per_sample, n_cell_types_high_res —
                     replaces the notebook's obs reads)
                     gloscope_dists/<ds>_gloscope_hvg<n>_pcadims<d>_dists.rds
                     (raw GloScope distance cache; sqrt applied at processing)
                     execution_times_<METHOD>_<DS>.feather (per-task exec log)
   └─ login tail: wait per method on the WATCHDOG id (benchmark_wait_watchdog:
      reads STATE from the status file; OK -> merge JOB_REPORT lines; FAIL ->
      email + fail closed) ->
      NAS reachability check -> checksums.md5 RDS integrity sidecar ->
      merge (shared 1.1.2_merge_execution_times.py, --no-cleanup
      --labels <methods> --datasets <DS> [--existing-log NAS log])
      -> execution_times.feather
   └─ rsync -> ${NAS_TARGET_DIR}/benchmark/
   └─ only then: delete this run's per-task logs (label x dataset cross product)
```

#### Files

| File | Role |
|---|---|
| `1_submit_hpc_array.sh` | Login-node submitter: `[--ds_name <DS>] [--methods prepare_pseudobulk,gloscope,mofa,pseudobulk,scitd,composition] [--partition <P>] [--force] [--sync-only <id1,id2,...>]`. Same dataset resolution as the Python submitter (jq: `use_for_benchmark == true` + `benchmark_analysis` view; `_*` keys skipped unless `--ds_name`) — via the shared `benchmark_submit_common.sh` helpers. Validates `--methods` (error on unknown), dedupes, and auto-prepends `prepare_pseudobulk` when `mofa`/`pseudobulk`/`composition` is requested without it (composition consumes the hvg2000 variant for ECODA_deconv). Submit order: `prepare_pseudobulk` array first with its own soft-gate WATCHDOG (artifact-completeness pass — all `PB_VARIANT_NAMES` variants present in `benchmark/pseudobulks/` per dataset — with the strict OOM-aware task-state gate only under `--force` or when variant files are missing), waited via its watchdog id, then the other arrays, each with its own strict watchdog (`benchmark_submit_watchdog` → `watchdog_main.sh`; `OUT_OF_MEMORY` tasks re-submit only their datasets with doubled `--mem` via the watchdog's own resubmit closure, clamped to `BENCHMARK_MEM_MAX`, default 500G; non-OOM failures fail closed in the watchdog's status file). All methods pinned to the CPU benchmark class on the sbatch command line (`--constraint=${BENCHMARK_CPU_CONSTRAINT}` + `--cpus-per-task=${BENCHMARK_CPU_CPUS_PER_TASK}` + `--mem=${BENCHMARK_MEM}`); `--partition <P>` override drops the constraint pin. Per-method manifest `${HPC_SCRATCH_DIR}/benchmark_manifest_<method>_$$.txt` (PID suffix; rebuilt every run); exports `METHOD` + `BENCHMARK_MANIFEST` + `FORCE_BENCHMARK`. `--sync-only` skips all submission (including the blocking `prepare_pseudobulk` wait); watchdog ids are gated via their status files (`benchmark_wait_watchdog`), other ids keep the strict `benchmark_wait_for_array` gate. The shared merge/sync/cleanup tail (NAS check → `checksums.md5` RDS integrity sidecar → sacct gates → exec-log merge → rsync → per-task log cleanup) runs via `benchmark_merge_sync_cleanup`. |
| `1.1_run_worker.sh` | `#SBATCH` worker (defaults overridden at submit time: 12h, 16 cpus, 128G). Same boilerplate as the Python worker (scontrol `Command=` SCRIPT_DIR recovery, `slurm_config.sh` + `cd ${PROJECT_ROOT}`); requires `METHOD` + `BENCHMARK_MANIFEST`; `DS_NAME` from the manifest via `sed -n ${SLURM_ARRAY_TASK_ID}p`. Per-task exec log `execution_times_<METHOD>_<DS_NAME>.feather` (deterministic name: each concurrent array has a distinct METHOD, re-runs overwrite the same file). Calls `1.1.1_prepare_pseudobulk.R` (for `prepare_pseudobulk`) or `1.1.1_run_benchmark_methods_r.R` via `${PIXI_RSCRIPT}` with `--config_path --ds_name --view benchmark_analysis --method --input_dir --results_dir --pseudobulk_dir --gloscope_cache_dir --log_file` (+ `--force` when `FORCE_BENCHMARK=1`), under the worker self-healing wrapper: transient-signature self-requeue only — R packages are read directly from the env library (slim `imports_worker_core.R` import set; no per-task staging since 2026-08-13), no thread pinning (hardware pinned for runtime comparability). |
| `1.1.1_run_benchmark_methods_r.R` | Method worker: sources `imports_worker_core.R` (Seurat/reticulate/dplyr) + `load_worker_functions.R` + `benchmark_hpc_utils.R`, attaching `MOFA2`/`scITD`/`EPIC`/`GloScope` conditionally on `--method` (mofa/scitd/composition only; gloscope stays namespaced); loads the h5ad → Seurat via `load_benchmark_seurat()` (reticulate: raw counts layer with X fallback + `X_pca_benchmark_analysis_hvg{1000,2000,3000}` obsm → reductions named `pca_benchmark_analysis_hvg{n}`), sets `VariableFeatures(seurat)` from the top-2000 `hvg_rank` genes and `seurat@misc$cell_type_low_res`/`label_col` from datasets.json. Memory: mofa skips the counts/embeddings materialization (metadata/labels from obs) and builds the Seurat only when the on-the-fly pb fallback triggers; gloscope fetches the embeddings, pseudobulk/scitd the counts; composition is OBS-ONLY (no Seurat at all: backed obs with `rename_leiden_cols(obs, view="benchmark_analysis")` + hvg2000 obsm + precomputed hvg2000 pb via `load_pb_variants()`). Dispatches on `--method` to the T2 drivers in `benchmark_pipeline.R` (`run_gloscope_hpc`/`run_mofa_hpc`/`run_pseudobulk_hpc`/`run_scitd_hpc`/`run_composition_methods_hpc`); mofa/pseudobulk/composition load the precomputed pb variants via `load_pb_variants()` (on-the-fly computation of ONLY the missing variants + atomic save if a variant RDS is missing — composition never triggers it: the submitter auto-prepends the prep array). Skip-if-exists: the method RDS exists → re-emits its exec-log rows (failure-resume must not lose timing from an aborted run) and skips all unless `--force`; else per-combo cache files are reused (also re-logged). Writes per-combo bundles + the method-level RDS + per-combo exec-log rows (composition additionally emits `<ds>_metadata.rds`). |
| `1.1.1_prepare_pseudobulk.R` | Prep worker: loads the Seurat (counts + `var["hvg_rank"]` only; no embeddings), runs `prepare_pseudobulks_hpc()` (variants `schvg2000`/`hvg2000`/`hvg500`/`hvg2000_bl`/`hvg1000`/`hvg3000`), writes `pseudobulks/<ds>_pseudobulk_<variant>.rds` atomically (tmp+rename) and one exec-log row per variant (`prepare_pseudobulk_<variant>`; `mem_GB` from the bundle). Sample names come from the obs column as-is (already standardized upstream — no `standardize_sample_names` re-application, which diverged hyphenated names for h5ads predating the python change, e.g. Adams). Skip-if-exists per variant unless `--force`. |
| `benchmark_hpc_utils.R` | (Not in `load_all_functions.R`; sourced explicitly by the HPC scripts.) Tiny `--flag value` arg parser (`parse_flags`), `get_h5ad_path()` (reuses `read_datasets_json`, which already maps `columns.label` → `label_col`/`output_file_name` → `output_file`), hvg_rank gene helpers (`get_hvg_rank_genes`, `make_hvg_sets`), the single source of truth for the pseudobulk variant names (`PB_VARIANT_NAMES`, `pb_variants_missing`), `load_benchmark_seurat()` (full Seurat build via `get_seurat_obj_from_h5ad`, counts layer + obsm embeddings), `load_pb_variants()` (read, or on-the-fly-compute only the missing variants, of the shared pseudobulks), `save_rds_atomic()` (tmp+rename), `peak_rss_gb()` (process peak RSS in GB from `/proc/self/status` VmHWM on Linux — the R-side equivalent of the python worker's `ru_maxrss`; NA_real_ off-Linux; monotonic-cumulative, logged at each combo's completion), `log_exec_row()` (shared exec-log feather schema `dataset, method, time_secs, mem_GB` with `mem_GB = NA_real_` when no peak measurement is available (nullable double, matching the Python writer's float64), read-modify-write append + dedup on (dataset, method) keep=last — mirrors `log_execution_time()` in `1.1.1_benchmark_methods_py.py`) and `run_ct_comps_analysis_worker()` (shared Pipeline B worker driver). |
| `watchdog_main.sh` | Compute-node watchdog job (submitted by `benchmark_submit_watchdog`, one per method array): 1 cpu / 2G / `WATCHDOG_TIME_LIMIT` (default 12h = shared-* MaxTime), method partition without the pinned constraint class, logs `${LOGS_DIR}/5_benchmark_watchdog_<label>_<id>.log/.err`, `#SBATCH --mail-type=END,FAIL`. Standard worker boilerplate (scontrol `Command=` SCRIPT_DIR recovery + `slurm_config.sh` + `cd ${PROJECT_ROOT}` + common source). Args after `--`: partition/throttle/log-prefix/worker-script + per-method flags. Owns the terminal wait + per-task gate + OOM escalation for its array via `benchmark_wait_oom_retry` with its own `watchdog_resubmit` closure (writes the reduced retry manifest, exports `METHOD`/`BENCHMARK_MANIFEST`, sbatch's the retry with the forwarded spec), and writes its status file `${WATCHDOG_STATUS_DIR}/<SLURM_JOB_ID>.status` (self-named; the id is unknowable at submit time) before exiting — survives SSH drops of the login tail. `soft-gate` mode = prepare_pseudobulk artifact gate (all `PB_VARIANT_NAMES` variants present per manifest-listed dataset → `STATE=OK` without the task-state gate; `--force` or missing variants → strict OOM-aware gate). No emailing / NAS access / exec-log merging (login tail only). Pure bash + slurm CLI — no pixi/R/Python. |
| `src/utils/imports_worker_core.R` + `src/utils/imports_worker_transzeroimp.R` + `src/utils/load_worker_functions.R` | Slim worker loaders (`src/utils/`, subset of the notebook loader): `imports_worker_core.R` attaches Seurat/reticulate/dplyr (R benchmark + prepare_pseudobulk workers; MOFA2/scITD/EPIC/GloScope attach conditionally per-method in `1.1.1_run_benchmark_methods_r.R`), `imports_worker_transzeroimp.R` attaches doParallel/foreach/reticulate/dplyr (trans/zeroimp workers; obs-only reads, no Seurat, no scECODA — `datrans` is a local function), and `load_worker_functions.R` sources the `load_all_functions.R` util files minus `imports.R` (notebook attach list, ~20 pkgs; notebooks only) and `plotting.R` (notebook-only); repo-wide env verification (attach ∪ namespaced-only ∪ worker/annotation packages) lives in `src/utils/env_check.R`. Both import subsets + env_check.R are smoke-checked by the guarded env-refresh scripts (`setup_env_sbatch.sh` [4/4], `refresh_env.sh` [3/3]). |
| `benchmark_submit_common.sh` | Shared helper sourced by the three benchmark submitters (python, R, trans/zeroimp) after `slurm_config.sh`: `benchmark_resolve_datasets <ds_name_arg>` (jq dataset resolution + `_*` skip convention), `benchmark_wait_array_terminal <job_id> <label>` (shared `squeue` exact-id poll + bounded `sacct` poll-until-terminal), `benchmark_wait_for_array <job_id> <label>` (terminal wait + fail-closed sacct gate over ALL rows), `benchmark_wait_oom_retry <job_id> <label> <resubmit_fn> <manifest> [status_file]` (OOM-auto-escalating variant for the benchmark submitters' own arrays, and the compute-node watchdog's engine: per-TASK states only — `.batch`/`.extern`/master rows excluded; all-COMPLETED records a `JOB_REPORTS` entry (only the FINAL successful retry) and passes; any non-COMPLETED, non-OOM state fails closed (per-task report + email); `OUT_OF_MEMORY` tasks map to datasets via `sed -n <task>p <manifest>` and re-submit ONLY those datasets with doubled `--mem` via `${resubmit_fn} <label> <ds_csv> <new_mem> <new_manifest>` (which echoes the new array job id), with the doubled value clamped to `BENCHMARK_MEM_MAX` (default 500G), looping until the ceiling or a 4-attempt cap — at the ceiling it fails closed with a per-task MaxRSS report). With the optional 5th arg `status_file` (watchdog mode) every terminal path writes the status file (`STATE=OK|FAIL`, `LABEL=`, one `JOB_REPORT=<label>|<id>|<wall>` line per gated array incl. the final retry id, `FAIL_REASON=`/`REPORT=` on FAIL) instead of emailing (compute nodes have no mail CLI) and instead of appending to `JOB_REPORTS`), `benchmark_submit_watchdog <array_id> <label> <manifest> <mode> <partition> <throttle> <log_prefix> <worker_script> <flags...>` (submits one `watchdog_main.sh` job per method array — 1 cpu/2G/`WATCHDOG_TIME_LIMIT`, method partition without the constraint pin; `strict` mode for method arrays, `soft-gate` for prepare_pseudobulk; forwards partition/throttle/log-prefix/worker/flags for its retry-array submissions; echoes only the watchdog job id; logs `${LOGS_DIR}/5_benchmark_watchdog_<label>_<id>.log/.err`), `benchmark_wait_watchdog <watchdog_id> <label>` (login-tail counterpart: terminal wait → status-file grace ≤2 min → `STATE=OK` merges the `JOB_REPORT=` lines into `JOB_REPORTS` and returns 0; `STATE=FAIL` emails "NOT synced — watchdog failed" with the report and exits 1; watchdog non-COMPLETED or COMPLETED without a status file fails closed with its sacct `State,ExitCode` + a pointer to its logs), `benchmark_pb_variant_names <benchmark_hpc_utils.R path>` (sed/grep parse of `PB_VARIANT_NAMES` — single source of truth, used by the soft-gate watchdog), `benchmark_bump_mem <mem>`/`benchmark_mem_ge <a> <b>` (mem-string helpers: `<N>G`/`<N>T` → 2N, suffix-aware compare), `benchmark_merge_sync_cleanup <labels...>` (NAS reachability check → writes the `checksums.md5` RDS integrity sidecar → merges per-task exec logs via the shared `1.1.2_merge_execution_times.py` (`--labels` × `DATASET_NAMES` cross product) → rsyncs `${HPC_SCRATCH_DIR}/benchmark/` → `${NAS_TARGET_DIR}/benchmark/` → deletes this run's per-task logs (scoped to the run's label x dataset cross product) plus a separate legacy `execution_times_task_*` sweep, so an overlapping submission's logs are never deleted). Sources `src/utils/bash/sync_status_email.sh` and sends best-effort sync-status emails (`notify_sync_status`): on gate failure ("NOT synced — task states", with a per-task report mapping task → `DATASET_NAMES[i-1]` (or the run's manifest, when a retry is gating) (state, elapsed, exit code) + array wall time instead of the raw sacct dump), on NAS-unreachable, and after a successful rsync — the final success/NAS-unreachable emails also carry a "Job durations" block (label, job id, array wall time per gated array, accumulated in the global `JOB_REPORTS` by `benchmark_wait_for_array`/`benchmark_wait_oom_retry` in submission order). The submitters' `--sync-only <id1,id2,...>` resume mode reuses these helpers as-is (skip submission; ids with a watchdog status file are gated via `benchmark_wait_watchdog`, other ids via the strict `benchmark_wait_for_array`, then merge/sync/cleanup). |

#### Key design details

- **Combo lists and result names are the legacy ones** (minus the GloScope
  `_sqrtmat` suffix, which is merged into `GloScope_hvg2000_pcadims30`):
  GloScope `hvg2000 × pcadims {10,30,50}` + `hvg{1000,3000} × pcadims30`;
  MOFA `hvg2000_factors{2,3,5,10,15}` + `hvg{1000,3000}_factors15`; Pseudobulk
  `schvg2000/hvg2000/hvg500/hvg2000_bl/hvg1000/hvg3000` + `CT_LR/HR_hvg{2000,500}`
  (each CT variant gated on its own ct col) + `{2,3,5,10,15}_PCA_dims`;
  scITD `hvg2000_factors{2,3,5,10,15}` + `hvg{1000,3000}_factors5`;
  composition `Avg_PCA_embedding` + `ECODA_deconv` + `ECODA_authors_LR` +
  `ECODA_authors_HR{,_NULL}` + `ECODA_authors_HR_top_varexp{0,0.1,…,0.9}` +
  `ECODA_authors_HR_{3most,2least,3least}_varcts` +
  `ECODA_authors_HR_{2,3,5,10,15}_PCA_dims` +
  `ECODA_seuratres_{0.1,0.4,2,5,20}` + `GloProp` + `Freq_highres` +
  `ECODA_HiTME_HR_layer{2,3}{,_top_varexp{0,0.1,…,0.9}}` + `ECODA_scATOMIC_HR`
  (HiTME/scATOMIC combos guarded: skipped with a warning when the ct column
  is absent from obs — the old notebook crashed). Composition defaults mirror
  `run_benchmark_analysis` (`factors_test`/`seurat_res`/
  `ECODA_top_varexp_hvct`) and run under `set.seed(123)` per dataset
  (ECODA_authors_HR_NULL is a null control, so its RNG stream difference
  from the old notebook runs is inconsequential).
- **Method timing includes the shared pseudobulk creation**: MOFA/Pseudobulk
  combo times = `pb_variant$time_secs` + method runtime (as the legacy code
  did with `exec_time(method) + exec_time_pb_norm`); GloScope combo times
  include the distance computation on a cache miss and are sqrt+read only on
  a hit (legacy `path_data` dist-cache semantics).
- **Sample-count guards** (make `_debug` usable): MOFA combos with
  `num_factors >= n_samples` are skipped (warning); scITD combos with
  `num_factors + 5 >= n_samples` are skipped (the tucker rank `c(f, f+5)`
  must be `< n_samples` — the 5-sample `_debug` dataset cannot run scITD at
  all); GloScope clamps `k` to `min(k, n_samples - 1)` and `n_pca_dims` to
  `min(n_pca_dims, ncol(embedding_matrix))` (the obsm stores
  `min(50, n_vars-1, n_obs-1)` PCs, so `_debug`'s 4 PCs would crash pcadims
  10/30/50 without the clamp).
- **The stored obsm + `var["hvg_rank"]` replace the legacy R-side
  recomputation**: `schvg2000` uses the top-2000 `hvg_rank` genes and scITD
  hvg1000/hvg3000 use the ranked sets; gene sets may differ slightly from the
  legacy Seurat `var.features` (expected; the new run is internally
  consistent).
- **`Pseudobulk_hvg2000_bl` behavior is preserved for legacy-result
  equivalence**, including the pre-existing `%in% black_list` no-op typo in
  `get_pb_deseq2` (pseudobulk.R:84) — flagged to the user, not silently
  fixed.
- **Result bundles keep the exact legacy names**; each bundle additionally
  carries `exec_time` (numeric seconds) and `mem_GB` (`peak_rss_gb()` at combo
  completion; replayed on cache reuse — the live cumulative peak would
  overstate a resumed combo's RAM). The notebook's exec-time section combines
  the NAS merged `execution_times.feather` (HPC rows) with bundle-derived rows
  for methods the feather does not cover (union by dataset+method,
  key-normalized: `Gongsharma_cmv_young_males` → `GongSharma`,
  `PILOT_hvg{n}_highres` → `PILOT_hvg{n}`); regenerated fresh on every knit
  (no `data/exec_times.rds` cache — a stale cache would pin old rows, e.g. NA
  mem_GB for R methods). RDS files use the default R
  serialization (R 4.5.x on both HPC and the user's Mac — compatible).
- **Bundle integrity sidecar**: the submit tail writes `benchmark/checksums.md5`
  (GNU md5sum over the RDS bundles); the notebook's `load_hpc_benchmark_results()`
  verifies listed files before `readRDS()` (RDS deserialization can execute
  code) and fails on mismatch — files not listed (pre-sidecar results) are
  read unverified.
- **Failure-resume re-logs timings**: cache hits (method RDS, per-combo
  bundles, prep variants) re-emit their stored `exec_time` + `mem_GB` rows
  into the current run's per-task log, so timing computed in an aborted run
  is not lost (the merge is scoped to the current run's labels × datasets).
- **Skip semantics**: per-combo cache files + method-level RDS; `--force`
  recomputes everything.

#### Usage

```bash
# 1. (prereq) preprocess + annotate the benchmark datasets (stages 2-4)

# 2. Run the R benchmark methods (all benchmark datasets; prep array gated
#    first; arrays monitored + sacct-gated + exec logs merged + NAS-synced
#    by the script). mofa/pseudobulk/composition auto-prepend
#    prepare_pseudobulk.
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name _debug --methods prepare_pseudobulk,pseudobulk   # debug smoke test
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --methods mofa --force                                     # recompute MOFA
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name _debug --methods composition                     # smoke test (obs-only)
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name _debug --methods gloscope --partition private-carmona-gpu  # constraint pin dropped
```

#### Test mode

Debug smoke tests (`_debug`, 5 samples): Pipeline B first
(`--ds_name _debug --analysis trans,zeroimp`), then Pipeline A
`--ds_name _debug --methods prepare_pseudobulk,pseudobulk`. GloScope works via
the k-adjust + `n_pca_dims` clamp, but `_debug` has only 4 PCs — validate
pcadims 10/30/50 on a real dataset; MOFA runs only factors 2,3 (5+ skipped by
the `n_samples` guard); scITD cannot run on 5 samples (tucker ranks) —
validate on a real dataset (e.g. Kim). Check `benchmark/results/`,
`benchmark/pseudobulks/` and the merged `execution_times.feather` on the NAS.

### Transformation and zero-imputation analyses on HPC (`src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/`)

The (fast) ECODA transformation analysis (`run_transformation_analysis` →
`datrans`) and zero-imputation analysis (`run_zeroimp_analysis`) run as two
SLURM arrays (one per analysis) submitted by one script — two workers per
dataset. Not part of the pinned runtime comparisons: partition
`${SLURM_PARTITION}` (default shared-cpu), 4 cpus, 32G, `--time=02:00:00`, no
constraint pin.

#### Workflow

```
1_submit_hpc_array.sh (login; [--ds_name <DS>] [--analysis trans,zeroimp]
  [--partition <P>] [--force]; two arrays, manifests
  benchmark_manifest_<analysis>_<pid>.txt, exports ANALYSIS +
  BENCHMARK_MANIFEST + FORCE_BENCHMARK)
   └─ 1.1_run_worker.sh (worker; branches on ANALYSIS)
        ├─ 1.1.1_run_transformation_analysis.R  (trans array)
        └─ 1.1.1_run_zeroimp_analysis.R         (zeroimp array)
             input:  h5ad obs-only backed read (reticulate, no counts matrix)
             output: ${HPC_SCRATCH_DIR}/benchmark/results/<ds>_trans.rds
                     ${HPC_SCRATCH_DIR}/benchmark/results/<ds>_zeroimp.rds
                     execution_times_<ANALYSIS>_<DS>.feather
                     (one row per dataset: trans_analysis / zeroimp_analysis)
   └─ NAS reachability check -> sacct gate (fail-closed) ->
      merge (shared 1.1.2_merge_execution_times.py) -> rsync -> delete logs
```

#### Files

| File | Role |
|---|---|
| `1_submit_hpc_array.sh` | Login-node submitter: `[--ds_name <DS>] [--analysis trans,zeroimp] [--partition <P>] [--force] [--sync-only <id1,id2,...>]`. Same dataset resolution as the other benchmark submitters; validates + dedupes the analyses; submits one array per analysis (`--partition=${SLURM_PARTITION}` default shared-cpu, `--cpus-per-task=4`, `--mem=32G`, `--time=02:00:00`, `MAX_NUM_CHUNKS_PARALLEL` throttle); `--sync-only` skips submission and gates the given job ids; same monitor/gate/merge/sync/cleanup tail as the other benchmark submitters. |
| `1.1_run_worker.sh` | `#SBATCH` worker (2h, 4 cpus, 32G): same boilerplate (scontrol `Command=` SCRIPT_DIR recovery, `slurm_config.sh` + `cd ${PROJECT_ROOT}`); requires `ANALYSIS` (trans or zeroimp) + `BENCHMARK_MANIFEST`; `DS_NAME` via `sed -n ${SLURM_ARRAY_TASK_ID}p`; per-task exec log `execution_times_<ANALYSIS>_<DS_NAME>.feather` (deterministic name: each concurrent array has a distinct ANALYSIS, re-runs overwrite the same file); calls the per-analysis R script via `${PIXI_RSCRIPT}` with `--config_path --ds_name --view benchmark_analysis --input_dir --output_dir --log_file` (+ `--force`), under the worker self-healing wrapper: transient-signature self-requeue only — R packages are read directly from the env library (slim `imports_worker_transzeroimp.R` import set, no Seurat; no per-task staging since 2026-08-13), no thread pinning (hardware pinned for runtime comparability). |
| `1.1.1_run_transformation_analysis.R` | obs-only backed read; `ct_comps = as.data.frame.matrix(table(Sample, obs[[cell_type_high_res]]))` with the `rowSums != 0` filter (as `get_ct_comp_df_seurat`; a plain data.frame so the dplyr verbs in `run_transformation_analysis` work); `labels` = `get_labels`-equivalent (per-sample slice(1) of `label_col`, names = Sample); calls `run_transformation_analysis()`; saves `<ds>_trans.rds` atomically; one exec-log row (`trans_analysis`). Skip-if-exists unless `--force`. |
| `1.1.1_run_zeroimp_analysis.R` | Same shape, calling `run_zeroimp_analysis()`; saves `<ds>_zeroimp.rds`; one exec-log row (`zeroimp_analysis`). |

#### Usage

```bash
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh            # both analyses, all datasets
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --ds_name _debug --analysis trans,zeroimp    # debug smoke test
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --analysis trans --force                     # recompute the trans analysis
```

### Data Processing

#### Complete Call Flow (Simplified)

```
notebooks/benchmark_analysis.rmd
    │
    ▼
load_hpc_benchmark_results(result_list, ds, path_nas_benchmark,
                           methods = c(gloscope, mofa, pseudobulk, scitd,
                                       composition))
    │   reads NAS benchmark/results/<ds>_{gloscope,mofa,pseudobulk,scitd,
    │   composition}.rds (HPC Pipeline A) + <ds>_{trans,zeroimp}.rds (HPC
    │   Pipeline B); bundles loaded FRESH every knit (no result_list.rds
    │   persistence, no rerun semantics)
    │   + <ds>_metadata.rds (labels, n_cells, n_samples,
    │     cells_per_sample, n_cell_types_high_res)
    ▼
python-feather methods (in the notebook loop, with labels from the metadata
bundle):
    ├──► process_mrvi_fig() / process_pilot_fig() /
    │    process_scpoli_fig() / process_qot_fig() /
    │    process_pilotgm_fig()  ◄── .feather files (HPC Python pipeline)
    └──► ... (create_result_bundle() -> calc_sep_score())

HPC Pipeline A (run_r_sample_embedding_methods/, SLURM arrays on
${HPC_SCRATCH_DIR}/benchmark/, synced to NAS benchmark/):
  prepare_pseudobulk -> prepare_pseudobulks_hpc() -> pseudobulks/<ds>_*.rds
  gloscope     -> run_gloscope_hpc()      -> process_gloscope_fig()  -> GloScope_*
  mofa         -> run_mofa_hpc()          -> process_mofa_bulk_fig() -> MOFA_*
  pseudobulk   -> run_pseudobulk_hpc()    -> process_pseudobulk_fig()/…ct_fig() -> Pseudobulk_*
  scitd        -> run_scitd_hpc()         -> process_scitd_fig()     -> scITD_*
  composition  -> run_composition_methods_hpc() -> process_coda_fig()/
                     process_gloprop_fig()/process_deconv_fig()/
                     process_avg_pca_embedding_fig()  -> ECODA_*/GloProp/
                     Freq_highres/Avg_PCA_embedding/ECODA_deconv
                     (+ <ds>_metadata.rds)

HPC Pipeline B (run_transformation_zeroimp_analysis/, 2 arrays):
  trans   -> run_transformation_analysis() -> datrans() -> <ds>_trans.rds
  zeroimp -> run_zeroimp_analysis()                       -> <ds>_zeroimp.rds

Shared scoring core (used by both notebook and HPC bundles):
  create_result_bundle()
      ├──► calc_sep_score()
      │   ├──► calc_modularity() ──► compute_KNN_from_dist()
      │   │                              compute_snn_graph()
      │   ├──► calc_sil()
      │   ├──► clust_eval()
      │   └──► anosim()
      └──► plot_mds() ──► calc_modularity(), clust_eval()
```


##### Key Design Patterns

1. **Factory Pattern**: `create_result_bundle()` is called by every processing function, producing a standardized result structure with `{scores, feat_mat, dist_mat, labels}`.

2. **Strategy Pattern**: `datrans()` dispatches to different transformation strategies (CLR, log, etc.) based on a method parameter.

3. **Pipeline Pattern**: Data flows linearly: AnnData object → Benchmark method data processing (Py/R) → Feature Matrix → Distance Matrix → Scores.
Note:
- Python methods (MrVI, scPoli, PILOT, QOT, PILOT-GM-VAE) and ALL R
  benchmark methods (GloScope, MOFA, Pseudobulk, scITD, composition:
  ECODA variants, GloProp, EPIC deconv, Avg_PCA_embedding, Freq_highres) +
  the transformation/zero-imputation analyses run on HPC and are read back
  into notebooks/benchmark_analysis.rmd via `load_hpc_benchmark_results()`
  (R bundles; labels/stats from `<ds>_metadata.rds`) or the `.feather`
  ingest functions (`process_mrvi_fig` etc., called directly in the dataset
  loop with the metadata-bundle labels)
- The notebook reads ZERO h5ad files
- Since not all methods provide a feature matrix, some directly output a distance matrix.

4. **Caching/Skip Logic**: `if (!method_name %in% names(res_list))` checks prevent re-running already-computed methods; HPC workers additionally skip per-combo cache files / method-level RDS unless `--force`.

---

#### LAYER 1: Entry Point / Orchestration
```
load_hpc_benchmark_results(result_list, ds, path_nas_benchmark,
                           methods = c(gloscope, mofa, pseudobulk, scitd,
                                       composition))            [notebook]
├── reads HPC Pipeline A bundles (<ds>_{gloscope,mofa,pseudobulk,scitd,
│   composition}.rds; <ds>_metadata.rds read separately in the loop)
├── reads HPC Pipeline B results (<ds>_{trans,zeroimp}.rds)
└── loads everything fresh every knit (no result_list.rds persistence,
    no "entries already present are kept" rerun semantics)
```
**`load_hpc_benchmark_results`** is the notebook entry point. For each dataset
it loads the analysis outputs (benchmark separation scores, transformation
analysis, zero-imputation analysis) from the NAS `benchmark/results/`
directory; the composition bundle + `<ds>_metadata.rds` are produced by
`run_composition_methods_hpc()` (HPC Pipeline A), and
`run_transformation_analysis()` / `run_zeroimp_analysis()` (now only called by
the HPC Pipeline B workers) provide the computation behind them.

---

#### LAYER 2: Benchmark Analysis Dispatcher
```
[notebook loop] python-feather methods only (labels from <ds>_metadata.rds):
├── process_mrvi_fig(...)        → MrVI
├── process_pilot_fig(...)       → PILOT
├── process_qot_fig(...)         → QOT
├── process_pilotgm_fig(...)     → PILOT-GM-VAE
└── process_scpoli_fig(...)      → scPoli

(run_benchmark_analysis(...) — DEPRECATED, kept for reference only:
 old notebook dispatcher; all R methods below moved to HPC)
```
HPC Pipeline A drivers (in `benchmark_pipeline.R`, called by the
`run_r_sample_embedding_methods/` workers; result names preserved):
```
prepare_pseudobulks_hpc(...) → shared DESeq2 pseudobulks (schvg2000/hvg2000/
                                hvg500/hvg2000_bl/hvg1000/hvg3000), per-variant
                                list(pb, time_secs, mem_GB)
run_gloscope_hpc(...)        → GloScope_hvg2000_pcadims{10,30,50},
                                GloScope_hvg{1000,3000}_pcadims30
run_mofa_hpc(...)            → MOFA_hvg2000_factors{2,3,5,10,15},
                                MOFA_hvg{1000,3000}_factors15
run_pseudobulk_hpc(...)      → Pseudobulk_schvg2000/hvg2000/hvg500/hvg2000_bl/
                                hvg1000/hvg3000, CT_LR/HR_hvg{2000,500},
                                {2,3,5,10,15}_PCA_dims
run_scitd_hpc(...)           → scITD_hvg2000_factors{2,3,5,10,15},
                                scITD_hvg{1000,3000}_factors5
run_composition_methods_hpc(...) → Avg_PCA_embedding, ECODA_deconv,
                                ECODA_authors_LR/HR/HR_NULL,
                                ECODA_authors_HR_top_varexp{0,0.1,…,0.9},
                                ECODA_authors_HR_{3most,2least,3least}_varcts,
                                ECODA_authors_HR_{2,3,5,10,15}_PCA_dims,
                                ECODA_seuratres_{0.1,0.4,2,5,20}, GloProp,
                                Freq_highres, ECODA_HiTME_HR_layer{2,3}*
                                (guarded), ECODA_scATOMIC_HR (guarded);
                                + <ds>_metadata.rds
```
Each driver times every combo with `exec_time()` (MOFA/Pseudobulk include the
shared pseudobulk creation time) and returns bundles carrying a numeric
`exec_time` + `mem_GB` (peak RSS); the workers write per-combo cache files,
method-level RDS and exec-log rows.

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
├── DESeq2.normalize(matrix, metadata, n_hvg, batch_col, blind, correct_batch)
│                       → DESeq2-like normalization (defaults `blind=TRUE`, no
│                         batch; batch-effect mode: `batch_col` + `blind=FALSE` +
│                         `correct_batch=TRUE` for batch-only limma removal; vst
│                         sparsity fallback chain for <1000 genes > mean 5)
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

## Batch Effect Analysis

Analysis lives in `notebooks/batch_effect_analysis.rmd` (notebook-based, under
expansion per reviewer comments). Methods to implement — with per-method batch
mitigation strategy — are tracked in TODO.md Phase 4: ECODA batch-associated
cell-type removal, Pseudobulk (DESeq2 + limma), MrVI (native `batch_key`), and
GloScope / PILOT-GM-VAE on `adata.obsm["X_pca_harmony"]` (the Harmony space
created in `1.1.1_preprocess.py`).

## Legacy pipeline notes

- `Preprocess_datasets.Rmd` (repo root, 905 lines) was the superseded legacy
  Seurat preprocessing pipeline; it was deleted after the audit (git history
  preserves it). Its steps are ported as follows:
  - scGate models + ProjecTILs ref maps → `2.0_create_scgate_db.R` /
    `${SCGATE_DB_PATH}` + `2.1.1_process_chunk.R` (`HOME_REF_DIR`).
  - Kfoury `cells_lowres` creation → `src/2_dataset_specific_preprocessing/1.4.1_create_kfoury_lowres_ct.R`.
  - Gene symbol standardization + preprocess loop → `src/utils/py/gene_utils.py` /
    `1.1.1_preprocess.py`; HiTME/scATOMIC column whitelist →
    `3.1_merge_annotations.py`.
  - GongSharma "clearcut age" split → dropped; deferred to the kept draft
    `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` (other-subsetting
    conditions).
  - "Export datasets without author annotation" (Lee/Zhang for scPoli) →
    dropped as a pipeline step. NOTE for benchmark interpretation: in the
    legacy pipeline Lee and Zhang were exported WITHOUT author annotation
    (scPoli was trained on them); this context explains certain decisions in
    `notebooks/benchmark_analysis.rmd` and in the paper.