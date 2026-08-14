> Kilo Code indexes code structure and function signatures automatically.
> AGENTS.md focuses on domain concepts, pipeline logic, and project conventions that indexing cannot infer.

# Paper/repo review and update strategy
This repo is about: ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts
Link to paper: https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full

## Open reviewer work
- Extend batch effect analysis and benchmark analysis (more datasets + more methods).
    - Additional datasets will be added by the user (human) — not by agents.
    - Pipeline/code work, method feasibility (PILOT-GM-VAE, QOT, PULSAR) and implementation drafts: see TODO.md (Phase 3-6).

# Agent Guardrails & Domain Terms

## Paper figures are sacred
- Every plot saved with a filename starting `Figure` or `Supp_fig` is a paper
  figure: KEEP and FIX if broken — never remove it. (Figure 4 CE/DF, Supp fig
  14, Supp fig 15, Supp fig 16AB, Supp fig 17-21, Figure 2A/B, Figure 3A/B.)
- Benchmark figure hierarchy:
  - Main benchmark figure (Figure 2A) = results across all datasets with only
    the default/main parameter setting per method.
  - Extended figure (Supp fig 15 / "For presentation") = additional
    non-standard methods (HiTME/scATOMIC annotations, cell-type pseudobulk,
    frequency-based composition, ...).
  - Parameter screening (Supp fig 2) = main methods across parameter ranges
    (robustness check that default parameters are faithful — no
    cherry-picking).
- `ECODA_PB_combo_*` (ECODA+Pseudobulk distance combos): legacy, kept
  commented-out in `benchmark_analysis.rmd` for internal testing only — NOT
  implemented in `run_benchmark_analysis` and NOT shown in publication
  figures (see TODO.md "Ideas for later").

## No-leakage (central premise)
- Bio labels (Status, sample.origin, cond, …) are ground truth **only**: they
  must NEVER be passed as a design covariate, batch key, or any other input to
  preprocessing, DESeq2 normalization, batch correction, HVG selection, or
  embedding steps. Batch correction is batch-only (no `design` argument in
  `removeBatchEffect`). Violation = supervised analysis, invalidates the
  paper's premise.
- `DESeq2.normalize()` semantics: defaults `blind=TRUE`, `batch_col=NULL`,
  `correct_batch=FALSE` → benchmark mode (design `~ 1`, legacy-equivalent, no
  correction). Batch-effect analysis uses `batch_col` + `blind=FALSE` +
  `correct_batch=TRUE` (batch-only `limma::removeBatchEffect`, no design
  protection). `get_pb_deseq2()` passes the same three params through.

# General rules
- Never drop defined versions in pixi.toml! This breaks reproducibility. If unsure (e.g. when removing a package or changing how it's imported, ask the user).
- Do not run pipeline scripts (e.g. .R, .py or .sh) for validation checks after implementing new code, unless the user asks for.
    - Validation of HPC pipeline scripts (e.g. .R, .py or .sh) will be run once the pipeline has been fully implemented, using a small debugging dataset (e.g. derived from the Joanito dataset)
- All HPC bash scripts must run with the working directory set to ${PROJECT_ROOT}:
  source `src/slurm_config.sh`, then `cd "${PROJECT_ROOT}"`. This is the established
  convention in every existing script — keep it for any new script (Python/R interop
  resolves repo-relative paths; see docs/ARCHITECTURE.md).
- sbatch-submitted scripts must NOT resolve `src/slurm_config.sh` from
  `BASH_SOURCE`: Slurm copies submitted scripts to
  `/var/spool/slurmd/job<id>/slurm_script`, so `BASH_SOURCE` points at the spool
  dir. The 10 sbatch workers (`1.1_submit_gongsharma.sh`,
  `1.2_submit_combinedpbmc.sh`, `1.3_submit_joanito.sh`,
  `1.4_submit_kfoury_lowres_ct.sh`,
  `3_scrnaseq_preprocessing/1.1_run_worker.sh`,
  `4_cell_type_annotation/2.1_run_worker.sh`,
  `5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`,
  `5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh`,
  `5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh`, and
  `utils/bash/setup_env_sbatch.sh`)
  recover `SCRIPT_DIR` from the job record via `scontrol show job`
  (`Command=` field) when `SLURM_JOB_ID` is set, with the `BASH_SOURCE` fallback
  for login-node execution. Keep this block in any new sbatch-submitted script.
- `slurm_config.sh` prepends the py-cuda13 env bin
  (`.pixi/envs/py-cuda13/bin`) to `PATH` so rpy2 (`src/utils/py/preprocess_utils.py`)
  finds R/Rscript when workers invoke `${PYTHON_BIN}` directly, and exports
  `LD_LIBRARY_PATH` with the env lib dir first (after the module loads) so R
  package `.so` deps resolve against the conda toolchain, mirroring what
  `pixi run` does automatically for `PIXI_RSCRIPT` workers.
- Search code with the built-in Grep/semantic-search tools or `git grep` (tracked files
  only). Never run raw `grep -rn "..." .` — it scans the gitignored 97 GB `data/` and
  `.pixi/` and will time out. If plain `grep` must be used, scope the path and exclude
  heavy dirs, e.g. `grep -rn "pattern" --exclude-dir={data,.pixi} src notebooks docs`.
  Note: `--exclude-dir` only prunes dirs found during traversal — the positional roots
  (`src`, `notebooks`, `docs`) are always searched.
- When giving code or suggestions to the user to check or run something, always provide the full copy&paste-ready code step-by-step (including e.g. if he needs to source slurm_config.sh first for paths on the hpc to be recognized).
- The user can run several terminals, each one accessing a different HPC login node, so he can run HPC bash scripts from any one of them, also in parallel. If possible, suggest which commands on the HPC can be run in parallel.

## Documentation organization
- Project Summary & Citation → README.md (primary home) (keep it short and concise! no details! ask user if unsure)
- Pipeline Call Graphs & HPC Layout → docs/ARCHITECTURE.md (keep here; remove detailed file-by-file logic trees from AGENTS.md)
- Agent Guardrails & Domain Terms → AGENTS.md (high-level rules; reference ARCHITECTURE.md for step details)
- Pending Tasks & Method Extensions → TODO.md (centralize planned conversions, script additions, reviewer extensions)


## Task Completion Workflow

Whenever you finish implementing a plan located in `.kilo/plans/`:
1. Move the completed plan file (but not other unrelated plan files!) into `.kilo/plans/archive/` (create the directory if it doesn't exist).
2. Stage all changed files related to the plan (not others), including the (specific) archived plan (`git add .`).
3. Create a git commit summarizing the work done.
4. Push the commit to the remote repository.


# Repo structure

## Documentation
- README.md (overview + usage), `docs/ARCHITECTURE.md` (pipeline call flow, file-role tables, HPC layout — the authority), TODO.md (task tracking). Keep documentation files up-to-date upon any changes.

## datasets.json
- This acts as ground truth for the datasets evaluated in this study
    - See datasets.json for most up-to-date list of datasets used and conditions.
- Do not change this file without asking
- The `_debug` entry (Joanito 5-sample subset, built by the Joanito step `1.3.1_prepare_joanito.R` into `${HPC_SCRATCH_DIR}/_debug/data/`) is registered here with both views. Convention: default-all script loops (`1_stage_data.sh`, `1_submit_hpc_array.sh`) skip `_*` keys unless explicitly requested via `--ds_name _debug`; `_debug.folder_name` is `null`, so staging skips it — the raw subset never lives on the NAS, only `_debug` *outputs* sync to NAS (`${NAS_TARGET_DIR}/_debug/output/`). Used for debugging `3_scrnaseq_preprocessing/`, `4_cell_type_annotation/`, `5_benchmark_methods/` (not the simple `1_stage_data/` / `2_dataset_specific_preprocessing/`). Details: see ARCHITECTURE.md.
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
    - Staging (`src/1_stage_data/`) → dataset-specific steps (`src/2_dataset_specific_preprocessing/`, e.g. `1.3.1_prepare_joanito.R` builds the `_debug` subset + `seqtec` batch column; `1.4.1_create_kfoury_lowres_ct.R` creates `cells_lowres`) → preprocess array (`src/3_scrnaseq_preprocessing/`) → annotation chunks + array + merge (`src/4_cell_type_annotation/`; `3_submit_merge.sh <DS>` merges `annotations_chunk_*.feather` into every view h5ad and syncs to NAS).
    - Preprocessed .h5ad files are **CSR-on-disk by construction** — required for selective backed-mode per-sample reads in annotation (details in ARCHITECTURE.md).
    - **Drafts (keep, not dead code)**: `preprocess_gongsharma.qmd` (GongSharma other-subsetting conditions) and `TODO_STUMP_preprocess_sikkema.qmd` (future Sikkema Lung dataset) in `src/3_scrnaseq_preprocessing/` are intentional drafts for future implementation — do NOT delete.
- **Stage 3 — Benchmark Analysis** (see [ARCHITECTURE.md](docs/ARCHITECTURE.md#benchmark-ecoda-transformation-and-ecoda-zero-imputation-analyses)): heavy R methods (GloScope, MOFA, Pseudobulk, scITD) on HPC via `run_r_sample_embedding_methods/` (pinned CPU class, `prepare_pseudobulk` prep array gated first); transformation/zero-imputation analyses on HPC via `run_transformation_zeroimp_analysis/` (two arrays: `trans`, `zeroimp`); Python methods (MrVI/scPoli/PILOT) on HPC via `run_python_sample_embedding_methods/`. `notebooks/benchmark_analysis.rmd` loads the HPC bundles via `load_hpc_benchmark_results()` and runs the fast composition-based methods (ECODA variants, GloProp, EPIC deconv, Avg_PCA, Freq_highres). Pending pipeline work (PILOT-GM-VAE, QOT/PULSAR): see TODO.md Phase 3.
- **Stage 4 — Batch Effect Analysis** (see [ARCHITECTURE.md](docs/ARCHITECTURE.md#batch-effect-analysis)): `notebooks/batch_effect_analysis.rmd`, under expansion (methods: ECODA batch-associated CT removal, Pseudobulk DESeq2+limma, MrVI, GloScope, PILOT-GM-VAE): see TODO.md Phase 4.

## R Modules for benchmark analysis (`src/5_run_benchmark_methods/` and `src/utils/`)

11 utility files loaded by `src/utils/load_all_functions.R` (notebook loader; also sources `src/utils/imports.R` — the canonical 42-package env-verification list — and `src/utils/plotting.R`, notebook-only). Workers use the slim counterparts instead: `src/utils/load_worker_functions.R` (the same 11 util files minus `imports.R`/`plotting.R`) + one of `src/utils/imports_worker_core.R` (R benchmark workers: Seurat/reticulate/dplyr) or `src/utils/imports_worker_transzeroimp.R` (trans/zeroimp: doParallel/foreach/reticulate/dplyr, no Seurat) — both smoke-checked by the guarded env-refresh scripts. Plus 3 benchmark-specific files in `src/5_run_benchmark_methods/` (`benchmark_methods_r.R`, `benchmark_pipeline.R` — both in `load_all_functions.R` and `load_worker_functions.R` — and `benchmark_hpc_utils.R`, which is HPC-only and sourced explicitly by the HPC worker scripts; details in ARCHITECTURE.md).

# HPC general information

## IMPORTANT:
- Do not use the login node to execute any code
    - If you do, you are disturbing all other users and this is unacceptable. When this happens we will most likely kill your process.
    - The login node should only be used to compile your code and submit a Slurm job. You must even use Slurm to run your tests. The debug-cpu and debug-gpu partitions are dedicated for small tests.
    - Defer from using tmux on the login node, except for the pixi install scripts (critical for HPC setup). Later pipeline steps are not critical (i.e. syncs can be done later at any time manually).

## Additional important points:
- Current repo lives on a local MacOS computer
- The user and agents can work on the HPC but always connecting from the user's computer
- If you need to test HPC cluster bash scripts:
    - The HPC cluster is only available if logged in to the UNIGE network (user might work from home (needs to connect to VPN) or from the office (has access to UNIGE network)).
        - If in the UNIGE network, you can log in with `ssh [REDACTED_HOST]` (user needs to enter the password).
    - By default, some scripts might need to be given the rights to be executed first, e.g. `chmod +x src/1_stage_data/1_stage_data.sh` or, for the whole repo: `find src -name "*.sh" -exec chmod +x {} \;1_stage_data.sh`
    - All `src/*.sh` are committed with the executable bit set (`100755`), so fresh clones are ready to run. If your HPC working copy predates that commit, set `git config core.fileMode false` (already done on the current HPC clone) instead of chmod-ing — otherwise mode-change noise blocks `git pull`.
- Heavy scripts are run on the HPC cluster, specifically located in these folders:
    - `src/1_stage_data`
    - `src/2_dataset_specific_preprocessing`
    - `src/3_scrnaseq_preprocessing`
    - `src/4_cell_type_annotation`
    - `src/5_run_benchmark_methods/run_python_sample_embedding_methods`
    - `src/5_run_benchmark_methods/run_r_sample_embedding_methods`
    - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis`
- The blocking submit-script tails (monitor → sacct gate → NAS rsync) are login-bound by design (NAS is login-node-only): an SSH drop never kills a running pipeline. Disconnect recovery = re-run the same command with `--sync-only <job-id>` (repeat the original `--ds_name`/`--methods`/`--analysis` flags) — never a bare re-submit; unknown/purged ids fail closed (exit 1, no sync). Sync-status emails are best-effort via the login node's mail CLI (`src/utils/bash/sync_status_email.sh`, `notify_sync_status`): "synced to NAS" / "NOT synced — reason", skipped silently if no mailx/mail/sendmail exists. Since the job-duration work, emails include a per-task/per-dataset (or per-method) completion breakdown + job durations: per-task lines (task → dataset, state, elapsed, exit code) + array wall time from sacct `Elapsed` (helpers `array_task_report`/`array_wall_time` in `sync_status_email.sh`; `n/a` when sacct is unavailable/purged), per-dataset aggregation for annotation chunk arrays (bounded failed-chunk list), per-step lines for the stage-2 dispatcher, per-view/per-dataset merge durations (`date +%s` wall clock, no sacct) for `3_submit_merge.sh`, and a per-method "Job durations" block in the benchmark final emails; `1_stage_data.sh` sends its own completion email (files staged, warnings, duration). `--sync-only` is supported by `1_submit_hpc.sh`, `3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, `4_cell_type_annotation/2_submit_hpc_array.sh`, and the three benchmark submitters (one comma-separated id per submitted array). To receive the emails, `export USER_EMAIL="you@example.com"` in the login-node shell profile (`~/.bashrc`); sbatch propagates it — without it, Slurm falls back to `${USER}@unige.ch`, which may be non-deliverable.
- `slurm_config.sh` is the HPC config file, used by all bash scripts, containing paths to the HPC cluster and other settings.
- pixi is a user-space binary on HPC at `~/.pixi/bin/pixi` (no module exists); the `py-cuda13` env lives at `.pixi/envs/py-cuda13/`. First-time env setup runs on a **worker node** via `sbatch src/utils/bash/setup_env_sbatch.sh` (submitted from the login node; 16 cpus/64G/2h defaults, parallel `MAKEFLAGS=-jN` via `R_SETUP_JOBS="${SLURM_CPUS_PER_TASK}"`; conda `compilers`/`make` deps make the R source build node-independent, GCCcore-module fallback); small login-node re-runs use `src/utils/bash/refresh_env.sh` (installs are okay to run on the login node, but never while jobs are active — both entry points serialize on `logs/env_refresh.lock`, see below).
- **`slurm_config.sh` path vars are NOT set in interactive login shells.** `HPC_SCRATCH_DIR`, `NAS_TARGET_DIR`, etc. only exist in scripts that `source src/slurm_config.sh` — typing `ls ${HPC_SCRATCH_DIR}/Wu/` in a fresh SSH session expands to an empty prefix and fails. In instructions/commands for the user, always give the concrete path (`${HOME}/scratch/ECODA_paper`, i.e. `[REDACTED_PATH]/scratch/ECODA_paper` on the current HPC clone) or prefix with `source src/slurm_config.sh &&`. The user can also export the vars in `~/.bashrc` (e.g. `echo 'export HPC_SCRATCH_DIR="${HOME}/scratch/ECODA_paper"' >> ~/.bashrc`).
- **Env refresh is a single serialized step — never run while jobs are active.** `pixi install` (conda) and `pixi run -e py-cuda13 setup` (remotes/GitHub source packages) both write into the same R library (`.pixi/envs/py-cuda13/lib/R/library`); a mutation racing with running array tasks corrupts it (observed: `digest/Meta/package.rds` missing → "there is no package called 'digest'"; `mime` lazyload DB missing; ABI warnings from R/Matrix bumps without source-package rebuilds). Two entry points, both guarded: `src/utils/bash/setup_env_sbatch.sh` (worker node, preferred — 16 cpus/64G parallel `make`; refuses to start while any OTHER job is active, own `SLURM_JOB_ID` excluded) and `src/utils/bash/refresh_env.sh` (login-node quick path for small re-runs; refuses while `squeue -u $USER` is non-empty; run it inside `tmux` so an SSH drop cannot interrupt an install). They share the lockfile `logs/env_refresh.lock` (PID + timestamp; stale if the PID is dead or > 24 h old), so a manual refresh and an sbatch build can never run concurrently. Note `pixi install` does NOT verify installed package files — a corrupt package dir survives re-installs; the `setup` task now runs an integrity check (Meta/package.rds presence, skipping R's own `translations` component, + full loader package check via `src/utils/imports.R`) and fails loudly. Definitive repair: `rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup` (rebuilds everything, incl. source packages against the current R/Matrix). HPC (2026-08-11): after env mutations (pixi install/prefix update, install.packages, setup), worker nodes can serve stale BeeGFS client-cache views ($HOME is BeeGFS, not NFS — see docs/hpc_docs/storage_on_hpc.md) — arrays may fail ~1 random task at the R/python import stage with ENOENT/empty reads on files that exist (files fine on the login node). Wait up to ~1h (typically faster for most nodes) after env changes before submitting arrays, or re-run (idempotent). Details: 
`.kilo/plans/1786440267018-setup-lazyload-integrity-check.md`. Since 2026-08-12 `PIXI_RSCRIPT` runs with `pixi run --as-is`,
so workers can never trigger an implicit lockfile update or env install at runtime — env drift now fails as stale-env errors
(missing package / `package or namespace load failed` signatures, transient-retry requeues then fail) instead of
auto-installing; the refresh-before-arrays workflow is mechanically enforced. **Worker self-healing** (since 2026-08-11) mitigates the remaining staleness risk at task level: all 5 array workers (`4_cell_type_annotation/2.1_run_worker.sh`, `3_scrnaseq_preprocessing/1.1_run_worker.sh`, and the three benchmark workers) source `src/utils/bash/worker_retry.sh` and, on a non-zero exit, grep the task's Slurm `.err` for transient signatures (`TRANSIENT_REQEX`: `No such file or directory`, `cannot open file`, `package or namespace load failed`, `No module named`, `missing from the pixi environment`, `.onLoad failed in loadNamespace`, `attempt to apply non-function` — the R workers' custom `load_my_packages` stop message, added 2026-08-13 after stale-view Seurat misses on 2 nodes never requeued, and the arrow `.onLoad` stale-view class, added 2026-08-14 after the `'arrow' ... attempt to apply non-function` observation, ...) and self-requeue via `scontrol requeue <array_job>_<task>` from the RUNNING task, capped by a per-(job,task) counter in `${HPC_SCRATCH_DIR}/_worker_retries/` (`WORKER_MAX_RETRIES`, default 3; cleared on success). R workers read packages **directly from the env library** (no per-task staging — the 2026-08-13 `stage_env_rlib` removal) with slim per-class imports: `imports_worker_core.R` (Seurat/reticulate/dplyr) for the R benchmark workers, `imports_worker_transzeroimp.R` (doParallel/foreach/reticulate/dplyr, no Seurat) for trans/zeroimp, both sourcing the utils via `load_worker_functions.R`; `imports.R` remains the env-verification list only. This keeps per-task startup I/O small, and the retry mechanism covers the residual stale-view flakes. OOM kills cannot self-requeue (the task is dead), so benchmark OOM auto-escalation runs as a **compute-node watchdog job** (`src/5_run_benchmark_methods/watchdog_main.sh`, one per method array, submitted right after each array by `benchmark_submit_watchdog`): it owns the terminal wait + per-task gate + `OUT_OF_MEMORY` detection and re-submits only those datasets with doubled `--mem` (128G → 256G → 500G, clamped to `BENCHMARK_MEM_MAX`, default 500G — the doubled value is never submitted above the ceiling, which fits the 512000 MB shared-cpu nodes; `BENCHMARK_MEM` itself is env-overridable) before failing closed — see `benchmark_wait_oom_retry` in `benchmark_submit_common.sh`. The watchdog survives SSH drops of the login tail (observed 2026-08-14: `4313942_3` scitd OOM'd at 128G, the login tail died, no retry was ever submitted): the tail only waits for the watchdog job id (`benchmark_wait_watchdog`, reading its status file `${HPC_SCRATCH_DIR}/_benchmark_watchdog/<watchdog_id>.status` — `STATE=OK|FAIL` + one `JOB_REPORT=` line per gated array) and then syncs; recovery after a tail death is `--sync-only <watchdog_id>` (status-file branch) or an idempotent re-run. `prepare_pseudobulk` gets a `soft-gate` watchdog (artifact-completeness pass first, strict OOM-aware gate fallback — the old login-side `benchmark_wait_prep_array` logic). Watchdog specs are modest (1 cpu/2G; `WATCHDOG_TIME_LIMIT` default 12h = shared-* partition MaxTime — never set it above the target partition's MaxTime, sbatch rejects at submit time). Annotation tasks: `--time` 02:00:00, per-attempt timeout cap 1800 s, budget-aware attempts vs remaining wall time (`SLURM_TIME_LIMIT` is minutes → ×60; `SLURM_JOB_START_TIME` is not set by Slurm → `proc.time()[3]`, conservative), per-sample checkpoints under `output/annotation_tmp/chunk_<N>/sample_<NN>.feather` (resume-only-failed-samples; final feather written atomically only when all samples succeed; never synced to NAS — excluded from both rsyncs; deleted by `1.1_prepare_chunks.py` on production rebuilds). Thread pinning (`export_worker_thread_env`, BLAS/OMP=1) is **annotation-only**: benchmark workers are hardware-pinned for runtime comparability and preprocess is intentionally multi-threaded. Retry safety: annotation/preprocess writes are atomic (partial files never pass the skip checks); the python benchmark worker deletes its per-dataset feathers on the requeue path (pyarrow `to_feather` non-atomic); R benchmark writes are `save_rds_atomic`. Node-shared `/srv/share/users/...` staging is a documented follow-up (metadata load from up-to-1000 concurrent copies), not implemented.
- **Worker environment invariants** (details in ARCHITECTURE.md):
    - Python is invoked via `PYTHON_BIN` and R via `PIXI_RSCRIPT` from `slurm_config.sh` — never bare `python`/`Rscript` (worker nodes may not have scanpy/anndata); on HPC all three interpreter pointers (`PYTHON_BIN`, `PIXI_RSCRIPT`, `RETICULATE_PYTHON`) come from the `py-cuda13` pixi env; any `pixi run`/`pixi run setup` on HPC must use `-e py-cuda13`; `PIXI_RSCRIPT` runs with `--as-is`
(`pixi run --as-is -e py-cuda13 Rscript --vanilla`) so workers can never update the lockfile or install/repair the env at
runtime (`--frozen` alone is NOT sufficient — it still installs the prefix when it mismatches the lockfile); env refreshes
happen only via the guarded `setup_env_sbatch.sh`/`refresh_env.sh`. `RETICULATE_PYTHON` is also exported so R workers' reticulate always uses the pixi python (mirrors the project-root `.Rprofile`, which only applies to non-vanilla sessions; the `.Rprofile` fallback targets `.pixi/envs/default` on macOS only — `py-cuda13` is linux-64-scoped and does not exist on osx-arm64).
    - Annotation (`2.1.1_process_chunk.R`) builds Seurat objects from the raw counts via `get_seurat_obj_from_h5ad()` — for preprocessed views from `layers["counts"]`, NOT from log-normalized `X`; the annotation union (`1.1_prepare_chunks.py`) carries counts in `X` by design (minimal layout: `X` = counts, no `layers`/`obsp`/`obsm` — anndata eagerly loads `layers` at backed open, which would add a multi-GB RAM floor per worker), so the `X` fallback is the designed primary path for union files (a `message()`, not a warning); feather names derive from the chunk file (`chunk_<N>.txt` → `annotations_chunk_<N>.feather`), not `SLURM_ARRAY_TASK_ID`; scGate models load from the shared `${SCGATE_DB_PATH}` cache (`aux/scGateDB.rds`) created by `2.0_create_scgate_db.R`; `2.1.1_process_chunk.R` also monkey-patches `cutoff.scATOMIC::em` right after `library(cutoff.scATOMIC)` with a bounded EM loop (upstream hang bug, see abelson-lab/scATOMIC issue <NUMBER>) — keep this patch on any edit touching the library block.
    - Annotation paths are per-dataset under `${HPC_SCRATCH_DIR}/${DS_NAME}/output` (see `config_helper.R`); `SAMPLE_COLNAME="Sample"` is exported by `slurm_config.sh`; `TISSUE_TYPE`/`NORMAL_TISSUE` are auto-exported per array task from `datasets.json` by `2.1_run_worker.sh` (via jq — `slurm_config.sh` performs the guarded `module load GCCcore/12.2.0` + `jq/1.6` centrally; consumers' `command -v jq` guards fail closed).

### HPC folder layout
- Full layout (scratch, reference atlases, NAS targets, env-var table): see [ARCHITECTURE.md](docs/ARCHITECTURE.md#hpc-folder-layout)
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/1_stage_data/1_stage_data.sh); results copied back to NAS after processing (typically in `*_submit_hpc_array.sh` upon completion)
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start


## HPC knowledge base (offline, indexed)
- Local snapshots of the UNIGE HPC documentation and community forum are committed
  here: `docs/hpc_docs/` (official wiki, /hpc/ namespace) and `docs/hpc_forum/`
  (forum topics from 2024-01-01, plus pinned; full topic list in
  `docs/hpc_forum/forum_index.md`). Snapshot date: 2026-08-11 — content may drift
  from the live sites. Refresh: `python3 src/utils/py/fetch_hpc_docs.py docs` and
  `... forum`, then commit.
- Retrieval protocol: before guessing cluster specifics (partitions, modules,
  quotas, storage semantics, known issues), run a `semantic_search` query targeted
  at `docs/hpc_docs` and `docs/hpc_forum` (e.g. "BeeGFS attribute cache
  staleness", "current issues on HPC cluster", "debug-cpu partition limits").
  Fall back to the live sites only if results are ambiguous:
  https://doc.eresearch.unige.ch/hpc/start and https://hpc-community.unige.ch/
  (2026 "current issues" thread: https://hpc-community.unige.ch/t/4185).
- The forum is Discourse; its JSON API (/latest.json, /t/<id>.json) is used by the
  scraper. The anonymous API does not expose post markdown (`raw`); the scraper
  converts the rendered `cooked` HTML to Markdown instead.

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
