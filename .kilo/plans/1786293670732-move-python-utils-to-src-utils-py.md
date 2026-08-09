# Move Python utils to src/utils/py/

## Decision

Move all 3 shared Python utility modules into a dedicated `src/utils/py/` folder
(user-confirmed). Rationale:

- Interop is one-directional (Python → R via rpy2 in `preprocess_utils.py`); no R
  file calls Python. The rpy2 interop sources R files via a **repo-relative** path
  pinned at import time (`ro.r('source("src/utils/load_all_functions.R")')`), which
  is independent of the py file's own location — the move is safe.
- `src/utils/` becomes purely R (matching its documented role as the R utility
  suite), `src/` root becomes purely stage folders + `slurm_config.sh`, and all
  Python shared code lives in one place.
- Both internal `Path(__file__).parents[N]` computations must be fixed; if missed
  they fail **loudly at import** (FileNotFoundError / R source error), not silently.

## Files to move (use `git mv` to preserve history)

1. `src/datasets_io.py` → `src/utils/py/datasets_io.py`
2. `src/gene_utils.py` → `src/utils/py/gene_utils.py`
3. `src/utils/preprocess_utils.py` → `src/utils/py/preprocess_utils.py`

No `__init__.py` files needed (repo uses namespace packages; none exist in
`src/` or `src/utils/` today).

## Task 1 — Fix internal path computations in the moved files

- `src/utils/py/preprocess_utils.py:14`: `PROJECT_ROOT = Path(__file__).resolve().parents[2]`
  → `parents[3]`. The `ro.r('source("src/utils/load_all_functions.R")')` call and the
  embedded R function definition stay unchanged (repo-root-relative).
- `src/utils/py/gene_utils.py:11`: `project_root = Path(__file__).resolve().parents[1]`
  → `parents[3]` (resolves `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` from repo root).
- `src/utils/py/datasets_io.py`: no path logic — no change.

## Task 2 — Update imports in the 3 caller files (5 import statements)

Import path `src.X` → `src.utils.py.X` in:

1. `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:11-13`
   - `from src.gene_utils import standardize_gene_symbols` → `from src.utils.py.gene_utils import standardize_gene_symbols`
   - `from src.datasets_io import read_datasets_json` → `from src.utils.py.datasets_io import read_datasets_json`
   - `from src.utils.preprocess_utils import load_input, apply_subset_vars` → `from src.utils.py.preprocess_utils import load_input, apply_subset_vars`
2. `src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py:62-63` (same renames) and `:93` (lazy import `from src.utils.preprocess_utils import load_input, apply_subset_vars` → `src.utils.py.preprocess_utils`)
3. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:55`
   - `from src.datasets_io import read_datasets_json` → `from src.utils.py.datasets_io import read_datasets_json`

The `sys.path.insert(repo_root)` lines in all 3 files are unaffected (repo root
is still inserted; `src.utils.py` resolves as a namespace package).

## Task 3 — Update comment-only references (module-path strings)

- `src/slurm_config.sh:36`: comment `src/utils/preprocess_utils.py` → `src/utils/py/preprocess_utils.py`
- `src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py`: comments at
  lines 45 (`preprocess_utils.py pins the embedded R working directory...`) and 250
  (`...never imports src.utils.preprocess_utils...`) → new path
- `src/4_cell_type_annotation/1.1_prepare_chunks.py:147`: comment mentions
  `preprocess_utils.load_single_input()` — function-name only, no path; update wording
  to `src.utils.py.preprocess_utils` for consistency (or leave; it is a bare function ref)

## Task 4 — Update documentation

- `docs/ARCHITECTURE.md` (7 spots — exact path strings to change):
  - Lines 9-10: `src/datasets_io.py::read_datasets_json()` → `src/utils/py/datasets_io.py::read_datasets_json()`
  - Line 101: `preprocess_utils.py pins the embedded R working directory...` → add
    `src/utils/py/` prefix; keep the `src/utils/load_all_functions.R` source path as-is
    (it stays correct)
  - Line 157: `consumed by src/gene_utils.py` → `consumed by src/utils/py/gene_utils.py`
  - Line 158: `imported at module level by src/utils/preprocess_utils.py` → `src/utils/py/preprocess_utils.py`
  - Line 242: `via src/datasets_io.read_datasets_json` → `via src/utils/py.datasets_io.read_datasets_json`
    (also the parenthetical "(repo-root import, like `1.1.1_preprocess.py`)" stays valid)
  - Line 714: `→ src/gene_utils.py /` → `→ src/utils/py/gene_utils.py /`
- `AGENTS.md:34`: `rpy2 (src/utils/preprocess_utils.py)` → `rpy2 (src/utils/py/preprocess_utils.py)`
- `README.md`: **no update needed** (zero references; keep concise per user request)

## Task 5 — Validation

1. `python -m py_compile` the 3 moved files (no deps required; catches syntax
   breakage from the `parents[N]` edits).
2. `git grep -n "src/datasets_io\|src/gene_utils\|src.utils.preprocess_utils\|src/utils/preprocess_utils"` —
   must only hit `src/utils/py/` paths, docs, and (unmodified) `.kilo/plans/` history.
3. `git grep -n "from src\.\|import src\." -- '*.py'` — all shared-module imports
   resolve to `src.utils.py.*`.
4. Optional (per AGENTS.md, pipeline runs are user-validated): from outside the repo,
   `python -c 'import sys; sys.path.insert(0, "<repo>"); from src.utils.py.preprocess_utils import load_input'`
   in the pixi env (`pixi run -e default python ...`) exercises the pinned-wd R source
   chain end-to-end (needs rpy2 + R on macOS).

## Out of scope / do not touch

- `.kilo/plans/*.md` — historical plans keep old paths.
- `notebooks/*.rmd` — only source `src/utils/load_all_functions.R` (R side, unchanged).
- `src/utils/load_all_functions.R` and all R modules — untouched.
- `datasets.json`, `slurm_config.sh` logic (comment only).
