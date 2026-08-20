# CWD Convention: Document + Harden

## Context / Assessment

`src/utils/preprocess_utils.py:11` runs `ro.r('source("src/utils/load_all_functions.R")')`, and `load_all_functions.R:7-18` sources 11 more files via repo-relative paths. Any Python script importing `preprocess_utils` (`_create_combinedpbmc_dataset.py`, `1.1.1_preprocess.py`) therefore required CWD = `${PROJECT_ROOT}` for the rpy2 interop.

Current state: all 9 HPC bash scripts already source `src/slurm_config.sh` and then `cd "${PROJECT_ROOT}"` (verified in `1_stage_data.sh:9`, `1_submit_hpc.sh:22`, `1.1_submit_combinedpbmc.sh:16`, `1.2_submit_joanito_batch_col.sh:15`, `1_submit_hpc_array.sh:5`, `1.1_run_worker.sh:17`, `1_prepare_chunks.sh:16`, `2_submit_hpc_array.sh:9`, `2.1_run_worker.sh:17`). So the convention works today but is implicit, enforced by nothing, and documented nowhere except the script docstring and historical `.kilo/plans`.

Decision (user confirmed): document the convention AND harden the code so Python callers become CWD-independent.

## Task 1 — Harden `src/utils/preprocess_utils.py`

Before the existing `ro.r('source(...)')` line, pin the embedded R working directory to the repo root (rpy2's R has its own wd; this does NOT change the Python process CWD):

```python
# load_all_functions.R sources its module files via repo-relative paths, so pin
# the embedded R working directory to PROJECT_ROOT at import time. All R
# interop calls below use absolute paths; this makes callers CWD-independent.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
ro.r('setwd')(str(PROJECT_ROOT))
ro.r('source("src/utils/load_all_functions.R")')
```

Notes:
- `Path(__file__).resolve().parents[2]` already resolves repo root correctly for `src/utils/preprocess_utils.py` (same idiom as `_create_combinedpbmc_dataset.py:44`).
- Use `ro.r('setwd')(str(...))` (rpy2 function call), not string interpolation, to avoid quoting issues.
- Safe: the only R functions used downstream (`convert_rds_to_raw_h5ad`, `create_clean_seuratv5_object`, `write_h5ad`) are called with absolute paths.

## Task 2 — Update stale docstring in `_create_combinedpbmc_dataset.py:25-27`

Replace the "Must run from ${PROJECT_ROOT}" caveat (now resolved by Task 1):

```
- Requires `module load GCCcore/12.2.0` for the R interop (rds->h5ad conversion).
  CWD-independent: preprocess_utils.py pins the embedded R working directory to
  ${PROJECT_ROOT} at import time.
```

Keep the heavy-load sbatch note (line 28-29) unchanged.

## Task 3 — Document in `docs/ARCHITECTURE.md`

In the "Preprocessing Pipeline" section, add a bullet right after the "NAS ↔ Scratch data flow" bullet (line 92):

```
- **Working-directory convention**: every HPC bash script sources `src/slurm_config.sh` and then `cd "${PROJECT_ROOT}"` (all existing scripts do this; keep it for any new script). `preprocess_utils.py` pins the embedded R working directory to `${PROJECT_ROOT}` at import time (its rpy2 interop sources `src/utils/load_all_functions.R` + 11 module files via repo-relative paths), so Python callers are CWD-independent — the `cd` remains belt-and-braces and required by convention.
```

Optionally extend the CombinedPBMC bullet (line 94) with "(script is CWD-independent; still submitted via the `1_submit_hpc.sh` dispatcher)".

## Task 4 — Add rule to `AGENTS.md` ("General rules" section, after line 26)

```
- All HPC bash scripts must run with the working directory set to ${PROJECT_ROOT}:
  source `src/slurm_config.sh`, then `cd "${PROJECT_ROOT}"`. This is the established
  convention in every existing script — keep it for any new script (Python/R interop
  resolves repo-relative paths; see docs/ARCHITECTURE.md).
```

## Out of scope

- README.md: usage examples already run from repo root; no change (optional one-liner if desired).
- Keep `cd "${PROJECT_ROOT}"` in all bash scripts (including future `src/5_run_benchmark_methods/...` scripts) — do not remove.
- `module load GCCcore/12.2.0` requirement is orthogonal (toolchain for rpy2 on HPC) — untouched.

## Validation

- `python -m py_compile src/utils/preprocess_utils.py src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py` (local, safe).
- Full R-interop smoke test (per AGENTS.md rule, ask the user or run on HPC debug partition): from a CWD outside the repo, `python -c 'import sys; sys.path.insert(0, "<repo>"); from src.utils.preprocess_utils import convert_rds_to_raw_h5ad_r'` — the import must succeed (exercises the R source chain via the pinned wd).
- No behavior change for existing bash-launched runs (they already cd'd to PROJECT_ROOT).

## Risks

- rpy2's `setwd` only affects the embedded R, not Python's `os.getcwd()`; no code in the two callers relies on R's wd differing from repo root.
- If a future caller of `preprocess_utils` relies on R's wd being the launch directory, the pinned wd would change behavior — unlikely (all R calls take absolute paths).
