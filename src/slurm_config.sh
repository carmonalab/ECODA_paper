#!/bin/bash
#
# Centralized SLURM and project environment variables
#
# Source this from any bash script that needs project paths or SLURM config:
#   source "$(dirname "${BASH_SOURCE[0]}")/../slurm_config.sh"
# (adjust relative path based on script depth)
#

# --- Project Paths ---
# NOTE: All path/env vars below are exported so that srun/sbatch children and
# Sys.getenv()/os.environ consumers (R config_helper.R, Python scripts) see them.
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export DATASETS_JSON_FILE="${PROJECT_ROOT}/datasets.json"
export LOGS_DIR="${PROJECT_ROOT}/logs"

# --- NAS Paths ---
export NAS_PREFIX="/srv/smednas515.unige.ch/carmona_smb"
NAS_BASE_DIR="${NAS_PREFIX}/DataCollections"
export NAS_SC_DIR="${NAS_BASE_DIR}/Standardized_SingleCell_Datasets"
export NAS_TARGET_DIR="${NAS_PREFIX}/Projects/ECODA_paper"

# --- HPC Scratch Paths ---
export HPC_SCRATCH_DIR="${HOME}/scratch/ECODA_paper"

# --- Pixi-managed interpreter commands ---
# HPC interpreters come from the py-cuda13 pixi env; plain `python` on worker
# nodes may resolve to a bare system python without scanpy/anndata etc.
# PYTHON_BIN: py-cuda13 python (has scanpy/anndata/torch etc.).
# PIXI_RSCRIPT: pixi Rscript entry point for the py-cuda13 env (word-splits
# into `pixi run --as-is -e py-cuda13 Rscript --vanilla`; the -e flag must stay
# between `run` and `Rscript` — pixi run without -e uses the default env instead).
# `--as-is` (= `--no-install --frozen`) is mandatory: a plain `pixi run` will
# update the lockfile and install/repair the env when the installed prefix does
# not match the lockfile, i.e. any worker could silently trigger a solve + env
# install at runtime — the "prefix update via `pixi run`" mutation class behind
# the 2026-08-11 corrupted library/ENOENT array failures. Env changes happen
# ONLY via the guarded entry points (setup_env_sbatch.sh / refresh_env.sh).
# Note `--frozen` alone is NOT sufficient: it freezes the lockfile but still
# installs the prefix when it mismatches the lockfile.
export PYTHON_BIN="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"
export PIXI_RSCRIPT="${HOME}/.pixi/bin/pixi run --as-is -e py-cuda13 Rscript --vanilla"

# rpy2 (imported by src/utils/py/preprocess_utils.py) needs R/Rscript on PATH;
# workers run PYTHON_BIN directly, so prepend the env bin (keeps python/R
# consistent across login node and workers).
export PATH="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin:${PATH}"

# --- jq (JSON parsing) ---
# Lmod's hierarchical tree on UNIGE clusters hides jq/1.6 behind its toolchain
# prerequisite GCCcore/12.2.0 (verified via `module spider jq` on bamboo, 2026-08-07).
# Both loads are guarded: a failing module load must never abort a `set -e` script
# nor print noise — the `command -v jq` guards in each consumer script are the
# fail-closed backstop. If the module tree updates (jq 1.6 no longer builds on newer
# GCCcore; EasyBuild pairs jq/1.7.1-1.8.1 with GCCcore/13.x and jq/1.8.1 with
# GCCcore/14.x), re-run `module spider jq` and bump both lines — or add jq to the
# py-cuda13 pixi env to drop the module dependency entirely.
module load GCCcore/12.2.0 >/dev/null 2>&1 || true
module load jq/1.6 >/dev/null 2>&1 || true

# --- R dynamic library resolution (rpy2 workers) ---
# R package .so files built by `pixi run setup` resolve their dependencies
# via LD_LIBRARY_PATH, which pixi's activation sets automatically but bare
# ${PYTHON_BIN} execution (rpy2 workers, e.g. 1.1.1_preprocess.py) does not.
# Without the env lib dir first, dyn.load() resolves against module (GCCcore)
# or node-system libs, and packages built with newer conda toolchains fail to
# attach on node images with older libstdc++/GLIBCXX (preprocessing array
# 4294806). Must come AFTER the module loads so the env lib dir wins.
export LD_LIBRARY_PATH="${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib:${LD_LIBRARY_PATH:-}"

# --- reticulate python (R workers) ---
# Pinned explicitly so R (2.1.1_process_chunk.R, imports.R) always uses the
# py-cuda13 python: reticulate's own discovery may otherwise pick a stray
# ~/.virtualenvs/r-reticulate or system python on the worker. Exported so it
# propagates through sbatch/srun. Mirrors .Rprofile (project root), which only
# applies to non-vanilla R sessions (.Rprofile is not read with --vanilla).
# Note: .Rprofile points at .pixi/envs/default on macOS only (py-cuda13 is
# linux-64-target-scoped and does not exist on osx-arm64).
export RETICULATE_PYTHON="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"

# --- Reference atlas paths (cell type annotation) ---
export NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
export HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"

# --- scGate model DB cache (cell type annotation) ---
# Created once by 2.0_create_scgate_db.R (via 2_submit_hpc_array.sh); loaded by
# 2.1.1_process_chunk.R workers so they do not download in parallel.
# SCGATE_DB_BRANCH is the single source of truth for the model DB version.
export SCGATE_DB_PATH="${PROJECT_ROOT}/aux/scGateDB.rds"
export SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"

# --- Sample column name (cell type annotation) ---
# 1.1.1_preprocess.py standardizes every dataset's sample column to "Sample".
export SAMPLE_COLNAME="Sample"

# --- Benchmark methods (src/5_run_benchmark_methods) ---
# Hardware PINNED for runtime comparability (same model/cores/RAM within each
# resource class). GPU methods (MrVI, scPoli): shared-gpu, H200
# (nvidia_h200_nvl, gpu[005-006]), 8 cores, 128G. CPU methods (PILOT, and the
# R benchmark pipeline — GloScope, MOFA, Pseudobulk, scITD via
# run_r_sample_embedding_methods/, which reuses this CPU class so cross-method
# runtime comparisons stay valid): shared-cpu, EPYC-7742 (V8, cpu[001-052]),
# 16 cores, 128G. All env-overridable.
# The private-carmona-gpu node is DELIBERATELY excluded here (different CPU
# model would flaw cross-method runtime comparisons; its GPU is H100, not H200,
# so BENCHMARK_GPU_CONSTRAINT would never match). _debug runs may target it via
# --partition ${SLURM_PARTITION_PRIVATE} — the benchmark submitter drops the
# --constraint pin on any --partition override. shared-gpu is the pinned GPU
# benchmark partition (real runs); ad-hoc --partition shared-gpu overrides are
# for debugging only (non-pinned, constraint dropped).
SLURM_PARTITION_BENCHMARK_GPU="shared-gpu"
SLURM_PARTITION_BENCHMARK_CPU="shared-cpu"
BENCHMARK_GPU_CONSTRAINT="nvidia_h200_nvl"
BENCHMARK_CPU_CONSTRAINT="EPYC-7742"
BENCHMARK_GPU_COUNT=1
BENCHMARK_GPU_CPUS_PER_TASK=8
BENCHMARK_CPU_CPUS_PER_TASK=16
BENCHMARK_MEM="${BENCHMARK_MEM:-128G}"
# Doubling ceiling for the benchmark submitters' OOM auto-escalation: an
# OUT_OF_MEMORY task is re-submitted with doubled --mem (128G -> 256G -> 500G),
# with the doubled value CLAMPED to this ceiling — a retry never requests
# more memory than the nodes fit (500G = 512000 MB = exactly the pinned
# shared-cpu node RAM, so it only schedules on a fully idle big node; 512G =
# 524288 MB would never fit and would hang the retry PENDING forever, which
# the submitter's squeue poll has no timeout for). OOM at the ceiling fails
# closed with a per-task MaxRSS report. Env-overridable per command, e.g.
# BENCHMARK_MEM_MAX=256G ./1_submit_hpc_array.sh.
BENCHMARK_MEM_MAX="${BENCHMARK_MEM_MAX:-500G}"
BENCHMARK_GPU_ARRAY_THROTTLE=4   # 4 H200s on gpu006
# Private node for _debug benchmark runs (not part of the pinned benchmark
# hardware): pass --partition "${SLURM_PARTITION_PRIVATE}" to the submitter.
SLURM_PARTITION_PRIVATE="private-carmona-gpu"

# --- SLURM Configuration ---
# Passed at submit time via `--partition="${SLURM_PARTITION}"` (sbatch
# directives do not expand variables). Comma list: jobs may land on any of
# the partitions (whichever frees resources first; order is not a preference).
# Used by the stages 2-4 submit scripts (CPU + GPU shared nodes and the
# private node). The benchmark submitter uses its own pinned vars above and
# is NOT part of this list; there shared-gpu appears only as the pinned GPU
# benchmark partition or a debug-only --partition override.
SLURM_PARTITION="shared-cpu,shared-gpu,private-carmona-gpu" # TODO: Adapt for specific pipelines

# --- User Info ---
# USER_EMAIL is the recipient for Slurm --mail-user and sync-status emails
# (notify_sync_status). It must be set by the user in their HPC shell profile
# (e.g. ~/.bashrc) — personal addresses must NOT be hardcoded in this repo.
# The default below is a best-effort guess and may be non-deliverable.
if [[ -z "${USER_EMAIL:-}" ]]; then
  export USER_EMAIL="${USER}@unige.ch"
  echo "WARNING: USER_EMAIL unset — falling back to ${USER_EMAIL}. Set USER_EMAIL in your HPC ~/.bashrc to receive Slurm + sync-status emails." >&2
else
  export USER_EMAIL
fi

# --- Parallelism ---
MAX_NUM_CHUNKS_PARALLEL=1000

# --- Worker self-healing (src/utils/bash/worker_retry.sh) ---
# WORKER_MAX_RETRIES: self-requeue cap on transient-failure signatures (grep
# on the Slurm .err file + scontrol requeue from the RUNNING task; counter
# file per (job, task) under ${HPC_SCRATCH_DIR}/_worker_retries/). R workers
# read packages directly from the env library (slim per-class imports, no
# per-task staging); the retry mechanism covers residual stale BeeGFS
# client-cache flakes.
WORKER_MAX_RETRIES=6
