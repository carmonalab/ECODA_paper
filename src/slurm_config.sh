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
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" &&
      -n "${HPC_SCRATCH_DIR:-}" ]]; then
  export HPC_SCRATCH_DIR
else
  export HPC_SCRATCH_DIR="${HOME}/scratch/ECODA_paper"
fi
# --- Immutable worker runtime contract ---
export ECODA_RUNTIME_MODE="${ECODA_RUNTIME_MODE:-host}"
export ECODA_RUNTIME_IMAGE="${ECODA_RUNTIME_IMAGE:-${HPC_SCRATCH_DIR}/_ecoda_runtime/ecoda-py-cuda13.sif}"
export ECODA_RUNTIME_MANIFEST="${ECODA_RUNTIME_MANIFEST:-${ECODA_RUNTIME_IMAGE}.manifest}"
export ECODA_RUNTIME_IN_CONTAINER="${ECODA_RUNTIME_IN_CONTAINER:-0}"
export ECODA_RUNTIME_PROFILE="${ECODA_RUNTIME_PROFILE:-default}"
export ECODA_APPTAINER_NV="${ECODA_APPTAINER_NV:-0}"
export ECODA_RUNTIME_PREFIX="${ECODA_RUNTIME_PREFIX:-}"
export ECODA_HOST_ENV_PREFIX="${PROJECT_ROOT}/.pixi/envs/py-cuda13"
export ECODA_HOST_PYTHON_BIN="${ECODA_HOST_ENV_PREFIX}/bin/python"
export ECODA_HOST_PIXI_RSCRIPT="${ECODA_HOST_ENV_PREFIX}/bin/Rscript --vanilla"


# --- Pixi-managed interpreter commands ---
# HPC interpreters come from the py-cuda13 pixi env; plain `python` on worker
# nodes may resolve to a bare system python without scanpy/anndata etc.
# PYTHON_BIN: py-cuda13 python (has scanpy/anndata/torch etc.).
# PIXI_RSCRIPT: direct Rscript binary from the py-cuda13 env. Do not launch
# worker R through `pixi run`: even `pixi run --as-is` executes r-base's
# activation-r-base.sh, which runs `R CMD javareconf` and rewrites the shared
# lib/R/etc/Makeconf and ldpaths on every activation. Concurrent workers then
# perform shared-prefix writes. The direct binary is safe because PATH,
# LD_LIBRARY_PATH, and RETICULATE_PYTHON are set below; only the guarded
# setup_env_sbatch.sh / refresh_env.sh entry points may invoke pixi for
# environment mutation.
export PYTHON_BIN="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"
export PIXI_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla"

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
# R package .so files built by `src/utils/setup_r_packages.R` resolve their dependencies
# via LD_LIBRARY_PATH, which pixi's activation sets automatically but bare
# ${PYTHON_BIN} execution (rpy2 workers, e.g. 1.1.1_preprocess.py) does not.
# Without the env lib dir first, dyn.load() resolves against module (GCCcore)
# or node-system libs, and packages built with newer conda toolchains fail to
# attach on node images with older libstdc++/GLIBCXX (preprocessing array
# 4294806). Must come AFTER the module loads so the env lib dir wins.
export LD_LIBRARY_PATH="${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib:${LD_LIBRARY_PATH:-}"

# --- reticulate python (R workers) ---
# Pinned explicitly so R always uses the py-cuda13 python: the annotation
# worker (2.1.1_process_chunk.R) and the R benchmark/transzeroimp workers
# (imports_worker_core.R / imports_worker_transzeroimp.R) attach reticulate;
# its own discovery may otherwise pick a stray ~/.virtualenvs/r-reticulate or
# system python on the worker. Exported so it propagates through sbatch/srun. Mirrors .Rprofile (project root), which only
# applies to non-vanilla R sessions (.Rprofile is not read with --vanilla).
# Note: .Rprofile points at .pixi/envs/default on macOS only (py-cuda13 is
# linux-64-target-scoped and does not exist on osx-arm64).
export RETICULATE_PYTHON="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"

# Container workers enter this branch only after ecoda_runtime_reexec_worker
# crosses the Apptainer boundary.  ECODA_RUNTIME_MODE alone never selects
# image-internal paths: submitters and host-side validators keep the host
# Pixi paths until a worker is inside the image.
if [[ "${ECODA_RUNTIME_IN_CONTAINER}" == "1" ]]; then
  if [[ -z "${ECODA_RUNTIME_PREFIX}" ]]; then
    echo "ERROR: ECODA_RUNTIME_PREFIX is required inside the immutable runtime." >&2
    return 1
  fi
  if [[ ! -x "${ECODA_RUNTIME_PREFIX}/bin/python" ||
        ! -x "${ECODA_RUNTIME_PREFIX}/bin/Rscript" ]]; then
    echo "ERROR: immutable runtime prefix is missing Python/Rscript: ${ECODA_RUNTIME_PREFIX}" >&2
    return 1
  fi

  _ecoda_strip_host_runtime_path() {
    local input_path="${1:-}"
    local item filtered=""
    local path_items=()
    local IFS=:
    read -r -a path_items <<< "${input_path}"
    for item in "${path_items[@]}"; do
      case "${item}" in
        "${ECODA_HOST_ENV_PREFIX}"|"${ECODA_HOST_ENV_PREFIX}"/*) ;;
        "") ;;
        *) if [[ -n "${filtered}" ]]; then
             filtered="${filtered}:${item}"
           else
             filtered="${item}"
           fi ;;
      esac
    done
    printf '%s' "${filtered}"
  }

  ECODA_RUNTIME_HOST_PATH_REST="$(_ecoda_strip_host_runtime_path "${PATH}")"
  ECODA_RUNTIME_HOST_LD_PATH_REST="$(_ecoda_strip_host_runtime_path "${LD_LIBRARY_PATH:-}")"
  export PYTHON_BIN="${ECODA_RUNTIME_PREFIX}/bin/python"
  export PIXI_RSCRIPT="${ECODA_RUNTIME_PREFIX}/bin/Rscript --vanilla"
  export PATH="${ECODA_RUNTIME_PREFIX}/bin${ECODA_RUNTIME_HOST_PATH_REST:+:${ECODA_RUNTIME_HOST_PATH_REST}}"
  export LD_LIBRARY_PATH="${ECODA_RUNTIME_PREFIX}/lib${ECODA_RUNTIME_HOST_LD_PATH_REST:+:${ECODA_RUNTIME_HOST_LD_PATH_REST}}"
  export R_HOME="${ECODA_RUNTIME_PREFIX}/lib/R"
  export RETICULATE_PYTHON="${ECODA_RUNTIME_PREFIX}/bin/python"
  export PYTHONNOUSERSITE=1
  # Read-only SIF layers cannot host Numba/Matplotlib user caches.  Keep
  # generated caches on the explicitly bound node-local temporary filesystem.
  export NUMBA_DISABLE_CACHING=1
  export XDG_CACHE_HOME="${TMPDIR:-/tmp}/ecoda-xdg-cache"
  export MPLCONFIGDIR="${TMPDIR:-/tmp}/ecoda-matplotlib"
  mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}"
  unset PYTHONHOME PYTHONPATH R_LIBS_USER R_LIBS_SITE R_ENVIRON_USER R_PROFILE_USER
  unset ECODA_RUNTIME_HOST_PATH_REST ECODA_RUNTIME_HOST_LD_PATH_REST
fi

# --- Reference atlas paths (cell type annotation) ---
export NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
if [[ "${ECODA_RUNTIME_IN_CONTAINER:-0}" == "1" &&
      -n "${HOME_REF_DIR:-}" ]]; then
  export HOME_REF_DIR
else
  export HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"
fi

# --- scGate model DB cache (cell type annotation) ---
# Created once by the canonical Stage 4 submitter via 2.0_create_scgate_db.R;
# loaded by annotation workers so they do not download in parallel.
# SCGATE_DB_BRANCH is the single source of truth for the model DB version.
export SCGATE_DB_PATH="${PROJECT_ROOT}/aux/scGateDB.rds"
export SCGATE_DB_BRANCH="41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4"
# The HiTME CellOntology dictionary is staged once beside the reference maps
# and symlinked into each R worker's tempdir; workers must never download it.
export SCGATE_MODEL_CACHE_DIR="${HOME_REF_DIR}/scGate_models"
export SCGATE_ONTOLOGY_BRANCH="master"

# --- Sample column name (cell type annotation) ---
# 1.1.1_preprocess.py standardizes every dataset's sample column to "Sample".
export SAMPLE_COLNAME="Sample"

# --- Benchmark methods (src/5_run_benchmark_methods) ---
# Standard benchmark MrVI/scPoli jobs remain pinned to H200 for comparable
# headline results. Flexible GPU jobs (debug, non-default, and batch-effect
# workloads) use the configured any-GPU partition list without a model pin.
SLURM_PARTITION_PRIVATE="${SLURM_PARTITION_PRIVATE:-private-carmona-gpu}"
SLURM_PARTITION_BENCHMARK_GPU="${SLURM_PARTITION_BENCHMARK_GPU:-shared-gpu}"
SLURM_PARTITION_BENCHMARK_CPU="${SLURM_PARTITION_BENCHMARK_CPU:-shared-cpu}"
BENCHMARK_GPU_CONSTRAINT="${BENCHMARK_GPU_CONSTRAINT:-nvidia_h200_nvl}"
BENCHMARK_CPU_CONSTRAINT="${BENCHMARK_CPU_CONSTRAINT:-EPYC-7742}"
BENCHMARK_GPU_COUNT="${BENCHMARK_GPU_COUNT:-1}"
BENCHMARK_GPU_CPUS_PER_TASK="${BENCHMARK_GPU_CPUS_PER_TASK:-8}"
BENCHMARK_CPU_CPUS_PER_TASK="${BENCHMARK_CPU_CPUS_PER_TASK:-16}"
BENCHMARK_GPU_DEFAULT_METHODS="${BENCHMARK_GPU_DEFAULT_METHODS:-mrvi,scpoli}"
BENCHMARK_GPU_DEFAULT_PARTITION="${BENCHMARK_GPU_DEFAULT_PARTITION:-${SLURM_PARTITION_BENCHMARK_GPU}}"
BENCHMARK_GPU_ANY_PARTITION="${BENCHMARK_GPU_ANY_PARTITION:-${SLURM_PARTITION_BENCHMARK_GPU},${SLURM_PARTITION_PRIVATE}}"
BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE="${BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE:-4}"
BENCHMARK_GPU_ANY_ARRAY_THROTTLE="${BENCHMARK_GPU_ANY_ARRAY_THROTTLE:-4}"
BENCHMARK_GPU_DEFAULT_TIME_LIMIT="${BENCHMARK_GPU_DEFAULT_TIME_LIMIT:-12:00:00}"
BENCHMARK_GPU_ANY_TIME_LIMIT="${BENCHMARK_GPU_ANY_TIME_LIMIT:-03:00:00}"
BENCHMARK_CPU_TIME_LIMIT="${BENCHMARK_CPU_TIME_LIMIT:-12:00:00}"
BENCHMARK_GPU_ANY_VRAM_PER_GPU="${BENCHMARK_GPU_ANY_VRAM_PER_GPU:-}"
# Backward-compatible name used by older callers/tests for the default class.
BENCHMARK_GPU_ARRAY_THROTTLE="${BENCHMARK_GPU_ARRAY_THROTTLE:-${BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE}}"
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
# Compute-node watchdog jobs (watchdog_main.sh, one per method array): own the
# terminal wait + OOM escalation so an SSH drop of the login tail cannot
# interrupt escalation. 1 cpu / 2G jobs on the method's partition (no
# constraint pin). WATCHDOG_TIME_LIMIT bounds them — the 12h default is the
# shared-* partition MaxTime (the workers' #SBATCH --time=12:00:00 is the
# partition max); a higher value is rejected by sbatch at submit time, so it
# must never be set above the target partition's MaxTime.
WATCHDOG_TIME_LIMIT="${WATCHDOG_TIME_LIMIT:-12:00:00}"

# --- SLURM Configuration ---
# Passed at submit time via `--partition="${SLURM_PARTITION}"` (sbatch
# directives do not expand variables). Comma list: jobs may land on any of
# the partitions (whichever frees resources first; order is not a preference).
# Used by stages 2-4 submit scripts (CPU + GPU shared nodes and the private
# node). Stage 5 uses the explicit default/flexible GPU classes above rather
# than inheriting this mixed pipeline partition list.
SLURM_PARTITION="shared-cpu,shared-gpu,private-carmona-gpu"

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
