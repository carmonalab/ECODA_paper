# Fix HPC sbatch path resolution + rpy2 R_HOME

## Context — root causes (established by live HPC debugging)

1. **`BASH_SOURCE` breaks under `sbatch`.** Slurm copies submitted scripts to
   `/var/spool/slurmd/job<id>/slurm_script` and runs them from there, so
   `SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"` resolves to the
   spool dir and `source "${SCRIPT_DIR}/../slurm_config.sh"` fails → exit 1.
   Affects every sbatch-submitted script (all 6 below). Validated fix on HPC:
   recover the real script path from the job record via `scontrol show job`
   (kfoury job COMPLETED after the patch).
2. **rpy2 cannot find R.** `src/utils/preprocess_utils.py:4` imports rpy2 at
   module level; jobs invoke Python as the raw `${PYTHON_BIN}` with no pixi env
   `bin/` on `PATH` → `RuntimeError: Unable to determine R_HOME`. Hits the
   stage-2 CombinedPBMC step and will hit stage-3 preprocess workers (same
   import). Fix once in `slurm_config.sh` (sourced by every script).
3. **Pixi env staleness + parallel re-sync race (operational).**
   `pixi.toml`/`pixi.lock` changed recently (commits `58201b9`, `dfcefcf`,
   `0ff2dc6`); if the HPC clone's env is behind the lock, each `pixi run`
   re-syncs — the parallel Joanito+Kfoury jobs raced and corrupted the update
   (`failed to remove directory ... No such file or directory (os error 2)`).
   Fix: serial env re-sync on the login node + README note.

## Decisions

- **Inline scontrol block, no refactor.** Scripts need `SCRIPT_DIR` *before*
  sourcing `slurm_config.sh`, so a shared helper cannot help. Use the exact
  block already validated on bamboo; keep the `BASH_SOURCE` fallback for
  login-node execution.
- **Unconditional PATH export** in `slurm_config.sh` (mirrors the
  unconditional `PYTHON_BIN`/`PIXI_RSCRIPT`). `py-cuda13` is linux-only and
  does not exist on osx-arm64 — a nonexistent PATH entry is harmless on macOS.
- **Keep `module load GCCcore/12.2.0`** in `1.1_submit_combinedpbmc.sh`
  (verified working on bamboo).
- **Keep default `pixi run`** (no `--no-install`): env staleness is solved by
  the documented login-node re-sync; `--no-install` would silently run stale
  envs.

## Tasks

### 1. Patch the 6 sbatch-submitted scripts

Files (all share the same `SCRIPT_DIR=` line; their `slurm_config.sh` source
line keeps its existing depth — `../` for 5 files, `../../` for the benchmark
worker):

- `src/2_dataset_specific_preprocessing/1.1_submit_combinedpbmc.sh`
- `src/2_dataset_specific_preprocessing/1.2_submit_joanito.sh`
- `src/2_dataset_specific_preprocessing/1.3_submit_kfoury_lowres_ct.sh`
- `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`
- `src/4_cell_type_annotation/2.1_run_worker.sh`
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`

Insert this block directly after the existing `SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"` line:

```bash
# Under sbatch the script is copied to /var/spool/slurmd/job<id>/slurm_script,
# so BASH_SOURCE no longer points at the repo. Recover the real path from the
# job record (Command= field); BASH_SOURCE fallback for login-node execution.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  SCRIPT_DIR="$(dirname "$(scontrol show job "${SLURM_JOB_ID}" -o | grep -o 'Command=[^ ]*' | head -1 | cut -d= -f2)")"
fi
```

No other changes to these files.

### 2. `src/slurm_config.sh`: prepend py-cuda13 bin to PATH

After the `PYTHON_BIN`/`PIXI_RSCRIPT` block (after line 34), add:

```bash
# rpy2 (imported by src/utils/preprocess_utils.py) needs R/Rscript on PATH;
# workers run PYTHON_BIN directly, so prepend the env bin (keeps python/R
# consistent across login node and workers).
export PATH="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin:${PATH}"
```

### 3. Documentation

- **README.md** (Installation, after the step-3 pixi setup code block, ~line 71):
  add a note that after pulling changes to `pixi.toml`/`pixi.lock` onto the HPC
  clone, re-run `pixi install -e py-cuda13` on the login node *before*
  submitting jobs — concurrent `pixi run` re-syncs from parallel jobs can race
  (observed: `failed to remove directory ... os error 2`).
- **AGENTS.md** (HPC conventions + worker env invariants): update the
  "All HPC bash scripts must run with the working directory set to
  ${PROJECT_ROOT}: source `src/slurm_config.sh`..." convention to state that
  sbatch-submitted scripts must NOT resolve `slurm_config.sh` from
  `BASH_SOURCE` (Slurm copies them to `/var/spool/slurmd/job<id>/slurm_script`);
  they recover `SCRIPT_DIR` via `scontrol show job` (reference the 6 patched
  scripts). Note that `slurm_config.sh` prepends the py-cuda13 env bin to
  `PATH` so rpy2 finds R.
- **docs/ARCHITECTURE.md** (light touch): mention the PATH export in the
  `slurm_config.sh` file-role entry and the sourcing requirement for
  sbatch-submitted scripts.

### 4. Local validation (no pipeline runs — per AGENTS.md)

- `bash -n` on all 7 modified bash files.
- `git diff` review: 6 identical blocks, 1 PATH line, doc edits only.

### 5. Commit (needed for HPC git-based sync)

- Stage only the intended files; message in repo style, e.g.
  `Fix sbatch SCRIPT_DIR resolution and rpy2 R_HOME on HPC`.

### 6. HPC rollout (user-run; HPC clone currently has the scontrol patch applied manually, uncommitted)

1. On the HPC clone, discard the manually applied patch so `git pull` applies
   cleanly (or `git stash` if other local changes exist):
   `git checkout -- <the 6 scripts>` (or `git stash`).
2. `git pull` on the HPC clone.
3. Serial env sync on the login node (never parallel pixi while stale):
   `pixi install -e py-cuda13 && pixi run -e py-cuda13 setup`
   (fallback if the remove-directory error recurs:
   `rm -rf .pixi/envs/py-cuda13` then repeat).
4. Re-run stage 2: `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
   — expect all 3 steps COMPLETED.
5. Debug pipeline: stage 3
   `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name _debug`,
   then stage 4 (`1_prepare_chunks.sh test _debug` →
   `2_submit_hpc_array.sh _debug` → `3_submit_merge.sh _debug`).

## Risks / notes

- `scontrol` on compute nodes assumed (standard Slurm client); validated on
  bamboo. `Command=` path is space-free (repo path), so the `grep -o` cut is safe.
- PATH prepend shadows system `python`/`R` on the login node — intended for
  interpreter consistency.
- If the Joanito step later fails with a *different* error (e.g. missing staged
  `.rds`), that is a staging issue — out of scope here.

## Out of scope

- Re-staging data; `datasets.json` changes.
- Refactoring the `SCRIPT_DIR` pattern into a shared helper.
- HPC-side execution by the implementing agent (per AGENTS.md, pipeline
  scripts are only run when the user asks).
