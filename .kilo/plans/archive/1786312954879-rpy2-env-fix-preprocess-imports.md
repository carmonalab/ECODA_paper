# Fix rpy2 R-package load failures in preprocessing workers (array 4294806)

## Context

Array 4294806 (preprocessing, 2026-08-09 23:39): tasks 1–3 (gpu004) COMPLETED,
tasks 4–14 (cpu007–010) FAILED with exit 1:

```
Error in load_my_packages(my_packages):
The following packages are missing from the pixi environment:
HiTME, MOFA2, scITD, GloScope, scECODA   (task 6: scITD absent from list; task 7: scECODA absent; task 8: + tidyverse)
```

Trigger chain: `1.1.1_preprocess.py` imports `src/utils/py/preprocess_utils.py`,
whose module level runs `ro.r('source("src/utils/load_all_functions.R")')`
(preprocess_utils.py:16). `load_all_functions.R` → `imports.R` →
`load_my_packages()` fails closed (`stop()`) if any of the 42 packages fails
`require()`.

### Evidence (all collected on HPC)

- All 5 packages + their deps (basilisk, rhdf5, SignatuR, SingleCellExperiment)
  ARE installed in `.pixi/envs/py-cuda13/lib/R/library/`; newest DESCRIPTION
  mtime = Aug 7 00:05 (scATOMIC) → env unchanged since; no `pixi install`
  touched the R library.
- Earlier R-heavy jobs succeeded: benchmark R workers (load_all_functions.R
  sourced unconditionally before any caching) COMPLETED on cpu007 at 17:14;
  annotation + Kfoury preprocess COMPLETED on gpu004.
- rcheck (pixi run Rscript, LD_LIBRARY_PATH = env libs): all packages load on
  cpu007 (00:20).
- A/B test (bare `${PYTHON_BIN}` + rpy2, LD_LIBRARY_PATH **unset**): all 42
  load on cpu007 (00:25). With env libs prepended: all 42 load.
- NodeList: successes = gpu004 only; failures = cpu007–010 only.

### Diagnosis

Not "missing packages", not "rpy2 vs Rscript", not "node permanently
incompatible". The untested difference between the failing array workers and
every passing test: the worker sources `slurm_config.sh`, which runs
`module load GCCcore/12.2.0` + `jq/1.6` (slurm_config.sh:50–51), putting the
GCC-12 module lib dirs (and the login-shell environment) into
`LD_LIBRARY_PATH`. R's `dyn.load()` then resolves the source-built packages'
`.so` dependencies against module/system libstdc++/libgomp instead of the
newer conda toolchain libs the packages were built against → packages needing
newer GLIBCXX/GLIBC symbols fail `require()` on node images that don't
provide them (per-node varying subset: cpu007–010 differ, gpu004 works).
`pixi run` (used by PIXI_RSCRIPT workers) sets LD_LIBRARY_PATH to the env
libs automatically; bare `${PYTHON_BIN}` execution (rpy2 workers) does not —
the PATH-prepend mirror in slurm_config.sh:39 is incomplete.

Caveat: A (LD_LIBRARY_PATH unset) also passed, so a residual trigger
(14 concurrent NFS package loads) is not fully excluded; the array rerun is
the decisive validation. The fix is strictly an improvement either way
(option B is proven to pass).

## Decisions

1. **Primary fix (general)**: export `LD_LIBRARY_PATH` with the py-cuda13 env
   lib dir FIRST in `src/slurm_config.sh`, placed AFTER the module loads.
   Completes the existing activation mirror (PATH already prepended) so any
   bare-binary rpy2/R usage resolves package deps from the build-time env on
   every node class. Fixes the failure class, not just these 5 packages.
2. **Structural fix (defense-in-depth)**: shrink the preprocessing rpy2 R
   import from the full 42-package `load_all_functions.R` to what
   `convert_rds_to_raw_h5ad` actually needs (Seurat + anndataR +
   `src/utils/seurat_utils.R`) — the same lighter pattern the annotation
   worker already uses (2.1.1_process_chunk.R). Removes benchmark-only
   packages (HiTME, MOFA2, scITD, GloScope, scECODA) from the preprocessing
   path entirely.
3. **Rejected alternatives**: `pixi run` for python workers (PYTHON_BIN is
   invoked quoted everywhere — word-splitting breaks; per-worker pixi
   startup cost; rattler NFS-cache warning noise in every log); rebuilding R
   packages with RPATH (brittle, non-relocatable, ~1h rebuild); sourcing
   `pixi shell-hook` in slurm_config.sh (overrides PATH/module conventions,
   pixi call on every source).

## Implementation tasks

1. **`src/slurm_config.sh`** — after the `module load` block (line ~51), add:

   ```bash
   # R package .so files built by `pixi run setup` resolve their dependencies
   # via LD_LIBRARY_PATH, which pixi's activation sets automatically but bare
   # ${PYTHON_BIN} execution (rpy2 workers, e.g. 1.1.1_preprocess.py) does not.
   # Without the env lib dir first, dyn.load() resolves against module (GCCcore)
   # or node-system libs, and packages built with newer conda toolchains fail to
   # attach on node images with older libstdc++/GLIBCXX (preprocessing array
   # 4294806). Must come AFTER the module loads so the env lib dir wins.
   export LD_LIBRARY_PATH="${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib:${LD_LIBRARY_PATH:-}"
   ```

2. **`src/utils/py/preprocess_utils.py`** — replace line 16:
   `ro.r('source("src/utils/load_all_functions.R")')` with a minimal load:

   ```python
   ro.r('''
   # Lighter than load_all_functions.R: preprocessing only converts .rds -> h5ad
   # (readRDS -> create_clean_seuratv5_object -> write_h5ad); the benchmark-only
   # packages (HiTME, MOFA2, scITD, GloScope, scECODA, ...) must not be pulled in
   # (same pattern as the annotation worker 2.1.1_process_chunk.R).
   source("src/utils/seurat_utils.R")
   library(Seurat)
   library(anndataR)
   ''')
   ```

   Update the adjacent comment block (lines 8–14) accordingly. No other R
   interop in this file needs imports.R symbols (verify: `convert_rds_to_raw_h5ad`
   uses readRDS / create_clean_seuratv5_object / write_h5ad only).

3. **Docs** — keep in sync with the repo convention:
   - `docs/ARCHITECTURE.md` line ~161 (environment propagation paragraph):
     document the new `LD_LIBRARY_PATH` export and why (rpy2 dyn.load;
     mirror of pixi activation; complements the PATH prepend).
   - `AGENTS.md` (the `slurm_config.sh` prepends … `PATH` bullet): mention
     `LD_LIBRARY_PATH` alongside PATH.
   - `docs/ARCHITECTURE.md` (preprocessing section, 1.1_run_worker.sh /
     preprocess_utils.py roles): note the minimal R import in the preprocessing
     path.

## Validation

1. Local syntax check: `python3 -m py_compile src/utils/py/preprocess_utils.py`
   (py_compile does not import rpy2/scanpy — safe on macOS).
2. HPC repro before/after (optional but recommended, pins the trigger):
   sbatch on `--partition=shared-cpu --nodelist=cpu007` that sources
   `src/slurm_config.sh` then runs `${PYTHON_BIN} -c
   'import rpy2.robjects as ro; ro.r("source(\"src/utils/load_all_functions.R\")"); print("OK")'`.
   Expect FAIL pre-fix, OK post-fix. (Post-fix, the repro should use the new
   minimal preprocess import as well.)
3. HPC: `git pull`; re-run
   `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (idempotent —
   tasks 1–3 skip via "Already processed"; failed tasks left no partial
   output files since `write_h5ad` happens only at the end of each view).
   Success = all 14 tasks COMPLETED + NAS rsync + sync-status email.
4. If CPU tasks still fail post-fix: residual trigger (concurrency/NFS) —
   next step is throttling the array (`%N` in the sbatch `--array` spec) and
   comparing node images; do NOT revert the fix (option B is proven correct).

## Risks / notes

- Env-libs-first LD_LIBRARY_PATH mirrors what `pixi run` already does for
  the benchmark R workers (PIXI_RSCRIPT); no behavior change expected there.
- `1.1.1_create_combinedpbmc_dataset.py` (stage 2) lazily imports
  preprocess_utils in its workers → automatically benefits from both changes.
- After implementation, follow the repo's Task Completion Workflow: move this
  plan to `.kilo/plans/archive/`, `git add .`, commit, push.
