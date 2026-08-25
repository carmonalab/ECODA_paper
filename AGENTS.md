# Repository Guidelines

## Project Overview

ECODA (Exploratory Compositional Data Analysis) is a reproducible R/Python workflow for unsupervised patient stratification from single-cell cohorts. It compares CLR-based cell-type composition, pseudobulk, and sample-embedding methods, then scores recovery of known biological groups.

### Non-negotiable scientific and repository rules

- **No label leakage.** Biological labels such as `Status`, `sample.origin`, `cond`, and `Disease_Identity` are ground truth only. Never pass them to preprocessing, HVG selection, normalization, batch correction, embeddings, or model covariates.
- `DESeq2.normalize()` benchmark defaults are `blind=TRUE`, `batch_col=NULL`, `correct_batch=FALSE` (`~ 1`). Batch-effect mode is batch-only: `blind=FALSE`, `batch_col=<batch>`, `correct_batch=TRUE`; never protect biological labels in `removeBatchEffect`.
- `datasets.json` is the dataset/view ground truth. **Do not modify it without explicit user confirmation.**
- Files beginning with `Figure` or `Supp_fig` are publication figures: fix them, never remove them. Figure hierarchy: `Figure 2A` uses default/main settings; `Supp fig 15` contains extended methods; `Supp fig 2` is parameter screening. Exclude legacy `ECODA_PB_combo_*` from publication figures.
- Preserve all version constraints in `pixi.toml` and the resolved `pixi.lock`.
- Use the `_debug` Joanito five-sample subset for routine verification. Do not run full cohorts for minor checks.

## Architecture & Data Flow

1. **Configuration:** `datasets.json` defines datasets, metadata columns, views, and filenames. `src/utils/datasets_io.R` and `src/utils/py/datasets_io.py` are the language-specific access layer.
2. **Data staging:** `src/1_stage_data/1_stage_data.sh` copies raw data from NAS to `$HOME/scratch/ECODA_paper`; `src/2_dataset_specific_preprocessing/` converts cohort-specific inputs.
3. **Canonical preprocessing:** `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` filters data, preserves raw counts in `layers["counts"]`, normalizes/log-transforms `X`, ranks HVGs, computes PCA, and creates Harmony/neighbors/Leiden outputs. RDS conversion and subset validation live in `src/utils/py/preprocess_utils.py`.
4. **Cell-type annotation:** `src/4_cell_type_annotation/` prepares sample chunks, runs annotation workers, checkpoints per-sample Feather output, and merges chunks.
5. **Benchmarking:** `src/5_run_benchmark_methods/` runs R and Python methods through SLURM arrays. Methods converge on sample feature matrices, distance matrices, or `create_result_bundle(feat_mat, labels, dist_mat)` bundles.
6. **Scoring and persistence:** `src/utils/scoring_metrics.R` computes silhouette, modularity, ANOSIM, ARI, and LISI. Results are saved atomically as `.rds` bundles with `checksums.md5`; Feather carries cross-language embeddings, distances, and execution logs.
7. **Analysis:** local notebooks consume precomputed results and generate publication figures; they do not rerun cohort preprocessing.

Operational concurrency is explicit rather than application-async: R uses `foreach`/`doParallel`, Python/R workers run in SLURM arrays, and shell watchdogs gate synchronization on `sacct`/`squeue`. Missing status, checksum mismatch, worker failure, or exhausted OOM retry must fail closed.

### Durable HPC execution

- Every full-cohort preprocessing, annotation, benchmark, evidence, or
  correction run MUST be launched through the checked-in
  `durable-hpc-gate-ecoda` profile. Direct SSH-launched long-running wrappers
  are not an acceptable substitute.
- Independent datasets MAY and SHOULD run in one SLURM array for a pipeline
  stage. The durable gate owns the array's terminal wait, accounting
  inspection, checksum/NAS audit, and Luna Max review; the next pipeline stage
  starts only after that gate is reviewed `COMPLETED`.
- Any array that can OOM MUST use a compute-node watchdog with automatic
  OOM-only resubmission of the affected manifest rows, bounded memory
  escalation, and fail-closed handling of non-OOM failures or an exhausted
  ceiling.
- After `launch`, arm exactly one unbounded durable `wait`; do not repeatedly
  poll `squeue`/`sacct` from the agent session. Perform one terminal `inspect`
  with every scheduler ID emitted by the wrapper, then the reviewer approval
  inspect. Use `status` only for non-mutating recovery checks.
- Short SSH commands for staging, code synchronization, and reading terminal
  evidence remain allowed. If the durable gate or required watchdog cannot
  represent a new wrapper, stop before launching full-cohort work and
  generalize the shared wrapper rather than bypassing the gate.

## Key Directories

- `src/1_stage_data/` — NAS-to-scratch staging.
- `src/2_dataset_specific_preprocessing/` — cohort-specific conversion and harmonization.
- `src/3_scrnaseq_preprocessing/` — shared Scanpy preprocessing and h5ad production.
- `src/4_cell_type_annotation/` — chunk preparation, annotation workers, and merge.
- `src/5_run_benchmark_methods/` — R/Python benchmark workers, submitters, watchdogs, and result synchronization.
- `src/utils/` — dataset I/O, scoring, imports, environment checks, preprocessing utilities, and shell environment setup.
- `notebooks/` — local analysis, publication figures, and dataset-onboarding reports.
- `tests/` — focused standalone R regressions. Shell watchdog tests remain beside benchmark code.
- `data/` — large/gitignored data; never scan recursively or delete recursively without explicit confirmation.
- `$HOME/scratch/ECODA_paper` on `bamboo` — data storage only, not a git clone. The HPC repository is `$HOME/ECODA_paper`.

## Development Commands

### Local macOS setup and analysis

```bash
pixi install
pixi run setup
pixi run Rscript tests/test_bassez_and_benchmark_regressions.R
bash src/5_run_benchmark_methods/test_oom_retry.sh
pixi run check-r-deps
```

Render analysis notebooks locally on macOS, never on HPC:

```bash
pixi run Rscript -e 'rmarkdown::render("notebooks/benchmark_analysis.rmd")'
pixi run Rscript -e 'rmarkdown::render("notebooks/batch_effect_analysis.rmd")'
```

### HPC setup

From the repository clone on `bamboo`:

```bash
cd "$HOME/ECODA_paper"
sbatch src/utils/bash/setup_env_sbatch.sh
```

For a guarded login-node refresh, first ensure no array jobs are active and use a persistent session:

```bash
tmux new -s env-refresh
cd "$HOME/ECODA_paper"
source src/slurm_config.sh
src/utils/bash/refresh_env.sh
```

### Representative pipeline commands

```bash
cd "$HOME/ECODA_paper"
source src/slurm_config.sh
src/1_stage_data/1_stage_data.sh --ds_name _debug
src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name _debug
src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods mrvi
```

Submit heavy work to SLURM; never preprocess or benchmark on a login node. Use `debug-cpu`/`debug-gpu` for short interactive checks (they have a max of 15 minutes). For non-blocking monitoring:

```bash
squeue -u "$USER"
sacct -j <job-id> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```


## Code Conventions & Common Patterns

- **Configuration over duplication:** resolve datasets and columns through `datasets.json` helpers; source `src/slurm_config.sh` in every HPC shell entry point.
- **Clean cross-language contracts:** h5ad stores raw counts in `layers["counts"]`; Feather stores tabular matrices with sample identity; RDS result bundles carry features, labels, and distances. Preserve row names/sample order and reject NA or mismatched identifiers.
- **Fail closed:** R uses `stop()`/`stopifnot()` and benchmark parallelism uses `.errorhandling="stop"`; shell scripts use strict status gates; Python validates required observation columns. Warn-and-skip is reserved for explicitly optional/missing artifacts.
- **Artifact safety:** use temporary files plus atomic rename, per-sample checkpoints, and MD5 verification before `readRDS`. Do not weaken checksum checks.
- **Naming:** numbered scripts encode pipeline order (`1_submit...`, `2_process...`, `3_merge...`); dataset-specific code stays under stage 2; shared helpers belong in `src/utils/`; shell environment variables are uppercase.
- **Dependencies/state:** there is no framework dependency-injection container or centralized in-memory state manager. Pass data/config explicitly through functions and CLI arguments; persistent state is files, manifests, status files, and checksums.
- **Performance:** retain sparse matrices and subset before densifying. Do not add avoidable full-cohort copies. `scPoli` intentionally densifies only the selected HVG subset.
- **R/Python imports:** use namespaced R package calls in HPC workers where established. Keep lazy imports for optional heavyweight Python methods.
- **Shell compatibility:** shared benchmark shell helpers and tests support Bash 3.2; avoid unsupported newer Bash syntax unless the target script explicitly requires it.
- **Comments/docs:** document scientific invariants and non-obvious scheduler behavior, not line-by-line mechanics. Update commands when entry points change.
- **SLURM spool recovery:** submitted scripts may execute from `/var/spool/slurmd/`. Recover their source directory before sourcing configuration:

  ```bash
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
      SCRIPT_DIR="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/ {print $2}' | xargs dirname)"
  else
      SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  fi
  ```

- **Completed implementation plans:** move plans from `.agents/plans/` to `.agents/plans/archive/`, stage the implementation and archived plan, then commit and push.

## Important Files

- `datasets.json` — central dataset, metadata, and view contract; confirmation required before edits.
- `pixi.toml`, `pixi.lock` — runtime environments, pinned dependencies, and Pixi tasks.
- `src/slurm_config.sh` — canonical HPC paths, interpreters, modules, resources, and retry ceilings.
- `src/utils/datasets_io.R`, `src/utils/py/datasets_io.py` — shared dataset configuration access.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` — canonical single-cell preprocessing entry point.
- `src/5_run_benchmark_methods/benchmark_pipeline.R` — transformations, zero imputation, distances, and parallel scoring.
- `src/5_run_benchmark_methods/benchmark_methods_r.R` — R method wrappers and result-bundle contract.
- `src/5_run_benchmark_methods/benchmark_submit_common.sh` — submission, retry, synchronization, and checksum logic.
- `src/utils/scoring_metrics.R` — benchmark metrics.
- `src/utils/bash/setup_env_sbatch.sh`, `src/utils/bash/refresh_env.sh` — serialized environment mutation and smoke checks.
- `README.md`, `docs/ARCHITECTURE.md` — operator workflow and pipeline map.

## Runtime/Tooling Preferences

- **Package/environment manager:** Pixi. Do not introduce Conda, renv, pip-only, npm, or a second lockfile for project dependencies.
- **Required versions:** R `4.5.2`; Python `3.13.*`. Supported base platforms are `osx-arm64` and `linux-64`; HPC workers use the `py-cuda13` Pixi environment.
- **Worker invocation:** never use bare `python`, `Rscript`, or ordinary `pixi run` inside jobs. Source the config and use its immutable commands:

  ```bash
  source src/slurm_config.sh
  cd "${PROJECT_ROOT}"
  "${PYTHON_BIN}" path/to/worker.py
  ${PIXI_RSCRIPT} path/to/worker.R
  ```

  `PIXI_RSCRIPT` includes `pixi run --as-is -e py-cuda13 Rscript --vanilla`, preventing runtime lock/environment mutation.
- Environment setup and refresh serialize on `logs/env_refresh.lock` and must not run while arrays are active.
- `bamboo` is the HPC cluster. Login nodes are for editing, compilation, staging, NAS sync, and SLURM submission only.
- Never run `rm -rf` against `$HOME/scratch` or `data/` without explicit user confirmation.
- No project-wide formatter, linter, Makefile, or CI workflow is configured. Match surrounding R, Python, and shell style rather than inventing a second convention.

## Testing & QA

QA is focused and script-based; there is no aggregate test runner, CI gate, coverage threshold, `pytest`, or `testthat` suite.

- `tests/test_bassez_and_benchmark_regressions.R` uses base-R assertions and temporary fixtures. It covers Bassez metadata fallback/validation and RDS checksum loading.
- `src/5_run_benchmark_methods/test_oom_retry.sh` uses deterministic shell stubs. It covers memory escalation, scheduler states, status files, notifications, retry ceilings, and fail-closed watchdog behavior.
- `src/utils/env_check.R` plus import smoke checks run during environment setup/refresh. Run manually with `pixi run check-r-deps` locally or `pixi run -e py-cuda13 check-r-deps` on Linux HPC.
- `src/4_cell_type_annotation/1_prepare_chunks.sh test [<DS_NAME>]` is a small pipeline mode, not the unit-test suite.

For changes, run the narrow contract test and exercise the `_debug` path. Add tests only for new observable behavior, boundary conditions, failure handling, or scientific invariants. Keep fixtures temporary and deterministic; stub Slurm/NAS effects rather than requiring live infrastructure.