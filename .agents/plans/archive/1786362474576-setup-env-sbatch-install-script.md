# Plan: HPC env install via sbatch worker node (setup_env_sbatch.sh) + install docs

## Context / problem

- `pixi run -e py-cuda13 setup` compiles R packages from source with `R CMD INSTALL`, which
  is **single-threaded per package** (Seurat's Rcpp code is the bulk) → 15–30+ min on the
  2-core/8 GB login node, blocking it and vulnerable to SSH-drop-interrupted installs
  (the known corruption history: `digest/Meta/package.rds`, `mime` lazyload DB).
- Decision made with the user: dedicated **sbatch install script** that builds the env on a
  worker node (16 cpus / 64G / 2h defaults, overridable at submit), with HPC-first install
  docs; local (macOS) install stays possible but is explicitly secondary (lightweight
  notebooks only).
- Worker nodes have internet (scGate DB download evidence in `2.0_create_scgate_db.R`) and
  shared `$HOME`/scratch; they may lack system gcc/make → add conda compilers to the env so
  `R CMD INSTALL` is node-independent.
- Speedup comes from parallelizing compilation (`MAKEFLAGS=-jN`), not from more cores alone.

## Decisions (locked)

- **Resources**: `#SBATCH --cpus-per-task=16 --mem=64G --time=02:00:00 --partition=shared-cpu
  --mail-type=END,FAIL`; overridable via `sbatch --cpus-per-task=N --mem=M ...`.
- **New file**: `src/utils/bash/setup_env_sbatch.sh` — a plain sbatch worker script
  (repo-consistent; submitted from the login node with `sbatch src/utils/bash/setup_env_sbatch.sh`).
- **Logs**: `#SBATCH --output=[REDACTED_PATH]/ECODA_paper/logs/setup_env_%j.out` and
  `#SBATCH --error=...%j.err` (hardcoded HPC repo path — single-user repo; slurm does NOT
  expand `${HOME}` in `#SBATCH` lines; script does `mkdir -p "${LOGS_DIR}"` at start).
- **`R_SETUP_JOBS="${SLURM_CPUS_PER_TASK:-16}"`** exported by the script; the pixi.toml
  setup task reads it for `MAKEFLAGS`.
- **Lockfile guard** shared by `refresh_env.sh` AND `setup_env_sbatch.sh`:
  `${PROJECT_ROOT}/logs/env_refresh.lock` (PID + timestamp; stale if PID dead or age > 24 h);
  both scripts refuse to run while it is held. `setup_env_sbatch.sh` additionally refuses if
  any OTHER slurm job of the user is active (`squeue -u "$USER"` excluding own
  `SLURM_JOB_ID`) — same guard `refresh_env.sh` already has.
- **Toolchain preflight** in the job, after `pixi install`: verify `gcc`/`make` are on PATH
  inside a `pixi run -e py-cuda13` session (conda compilers activation). If missing: try
  `module load GCCcore/12.2.0` (already loaded via `slurm_config.sh` for jq) and re-check;
  else fail with a clear message (run `refresh_env.sh` once on the login node, or verify the
  `compilers`/`make` pixi deps were installed).

## Tasks (ordered)

### 1. pixi.toml — hermetic toolchain (cuda13 feature only, linux-64)
Under `[feature.cuda13.target.linux-64.dependencies]` add:
```toml
# Hermetic C/C++/Fortran build toolchain for the R source-package compilation
# in [tasks.setup]: worker nodes may not have system gcc/make, and conda
# compilers make R CMD INSTALL node-independent (login node or worker).
compilers = "*"
make = "*"
```
(osx-arm64 unaffected — the feature is linux-64-scoped.)

### 2. pixi.toml — parallel compile in the setup task
Insert at the top of the `[tasks.setup]` R code (right after the `options(...)` line, before
`Sys.setenv(R_REMOTES_UPGRADE = 'never')`):
```r
  # Parallel source compilation: R CMD INSTALL builds with make; MAKEFLAGS=-jN
  # compiles translation units in parallel. Without it every Rcpp-heavy package
  # (Seurat, MOFA2, GloScope, ...) compiles single-threaded and setup takes
  # 15-30+ min regardless of node cores. Cap at 8 on the shared login node;
  # override with the R_SETUP_JOBS env var (e.g. 16 on a dedicated worker).
  setup_jobs <- suppressWarnings(as.integer(Sys.getenv('R_SETUP_JOBS', unset = '')))
  if (is.na(setup_jobs) || setup_jobs < 1) {
    setup_jobs <- min(8, max(2, parallel::detectCores() - 1))
  }
  Sys.setenv(MAKEFLAGS = paste0('-j', setup_jobs))
  cat('Using', setup_jobs, 'parallel compile jobs (MAKEFLAGS)')
```
**Hard constraints on this block** (it lives inside TOML `"""..."""` AND inside a shell
`Rscript -e "..."` wrapper): no `"`, no `$`, no backticks, and no `\n` escape inside R
string literals (TOML converts `\n` to a real newline, breaking the R string — validated
during the earlier session). Single quotes only.

### 3. New `src/utils/bash/setup_env_sbatch.sh`
- `#SBATCH` header per Decisions; `set -euo pipefail`.
- SCRIPT_DIR recovery via `scontrol show job` (AGENTS.md sbatch convention; `BASH_SOURCE`
  fallback for login-node execution), then `source "${SCRIPT_DIR}/../../slurm_config.sh"`,
  `cd "${PROJECT_ROOT}"`, `mkdir -p "${LOGS_DIR}"`.
- Lockfile guard (acquire at start, release on EXIT) + own-job-excluding squeue guard.
- Toolchain preflight helper as described (runs after `pixi install`).
- `export R_SETUP_JOBS="${SLURM_CPUS_PER_TASK:-16}"`.
- Steps with `echo` separators + `SECONDS`/`date +%s` timings:
  1. `"${HOME}/.pixi/bin/pixi" install --environment py-cuda13`
  2. toolchain preflight
  3. `"${HOME}/.pixi/bin/pixi" run -e py-cuda13 setup`
  4. smoke check: `pixi run -e py-cuda13 Rscript --vanilla -e 'invisible(lapply(c("digest","SeuratObject","Seurat","anndataR"), library, character.only=TRUE)); cat("All critical packages load OK\n")'`
- Executable bit `100755` (repo convention; `chmod +x` before commit).

### 4. Update `src/utils/bash/refresh_env.sh`
- Add the same lockfile acquire/release (so a manual login-node refresh and an sbatch build
  can never run concurrently). Keep the existing squeue guard. No other behavior change.

### 5. Docs
- **README.md** — restructure "Setup your environment" (step 3):
  - *HPC (recommended, main focus)*: clone → pixi binary → `sbatch src/utils/bash/setup_env_sbatch.sh`
    (worker node, 16 cpus/64G, parallel make; never while other jobs or a manual
    `refresh_env.sh` run are active; watch `logs/setup_env_<jobid>.out`; success =
    "R library integrity check OK" + "All critical packages load OK"). Login-node quick path
    for small re-runs: `tmux new -s env-refresh` → `src/utils/bash/refresh_env.sh`. Keep the
    existing serial-re-sync warning (concurrent `pixi run` re-syncs race).
  - *Local macOS (optional — lightweight notebooks only)*: `pixi install && pixi run setup`
    (py-cpu); note the first setup also compiles the R packages (Seurat, anndataR, ...) and
    takes a while; the notebooks `notebooks/benchmark_analysis.rmd` /
    `notebooks/batch_effect_analysis.rmd` run on precomputed distance-matrix results and are
    the intended local use; the heavy pipeline stages are HPC-only.
- **AGENTS.md** — fix the stale line "First-time env setup ... runs via the sbatch job in
  README.md": point to `src/utils/bash/setup_env_sbatch.sh` (submit from login node) and
  refresh_env.sh; extend the env-refresh bullet with the shared lockfile
  (`logs/env_refresh.lock`) and the new worker-node path.

### 6. Validation (no pipeline runs — per AGENTS.md)
- `bash -n` on `setup_env_sbatch.sh` and `refresh_env.sh`.
- `python3 -c "import tomllib"` parse of pixi.toml; extract the setup `cmd` (strip the
  `Rscript -e "` wrapper), assert no `"`/`$`/backtick and no `\n` inside R string literals,
  and `Rscript --vanilla -e 'parse(...)'` the extracted R code.
- HPC run (user-driven, when no other jobs are active): `git pull` → `sbatch
  src/utils/bash/setup_env_sbatch.sh`; expect the job to print the integrity check OK + smoke
  check OK and exit COMPLETED; confirm `logs/setup_env_<jobid>.out` timings.
- Follow the repo's Task Completion Workflow: archive this plan to `.kilo/plans/archive/`,
  `git add .`, commit (repo message style), push; then HPC `git pull` and the sbatch run.

## Risks / notes

- **conda compilers under pixi**: if pixi does not run the compilers' activation scripts,
  `gcc`/`make` may not appear on PATH in `pixi run` sessions — the toolchain preflight
  catches this; fallback is the GCCcore module (via `slurm_config.sh`) or one login-node
  `refresh_env.sh` run; do not silently proceed without a compiler.
- First run after adding `compilers`/`make` re-runs `pixi install` (adds ~1–2 GB to the env;
  no R-package recompiles are triggered by that).
- `-j16`/64G is ample for Rcpp compiles; if the partition's default mem is lower, `--mem=64G`
  is explicit in the script.
- A concurrent **manual** (non-slurm) setup on the login node is only stopped by the shared
  lockfile — document that both entry points use it.
- The `#SBATCH --output` path is hardcoded to the current HPC repo location (single-user
  repo); comment this in the script.

## Out of scope

- Running the heavy pipeline stages locally; local install beyond the lightweight notebooks.
- `ccache` for future rebuilds; speeding up `pixi install` downloads (already parallel).
