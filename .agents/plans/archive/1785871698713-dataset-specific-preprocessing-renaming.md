# Rename `src/2_dataset_specific_preprocessing/` scripts to standard pipeline naming

## Goal

Make the file names in `src/2_dataset_specific_preprocessing/` follow the repo's
standardised pipeline script naming convention (decimal call-depth numbering,
established in `.kilo/plans/1785446558424-pipeline-script-renaming.md`):

- Depth `N`: top-level entry point, `N_action.sh`
- Depth `N.N`: directly called script, `N.N_action.sh`
- Depth `N.N.N`: grandchild (actual processing script), `N.N.N_action.py|R`
- No leading underscore prefix; the `_` prefix is reserved for non-call-chain
  helpers (e.g. `src/utils/preprocess_utils.py`), which these two files are not
  — they are pipeline processing steps invoked by the step scripts.

Reference implementations: `1.1.1_preprocess.py`, `2.1.1_process_chunk.R`.

## Current state

| File | Status |
|---|---|
| `1_submit_hpc.sh` | Conforms (depth-1 dispatcher, non-array → stays `1_submit_hpc.sh`, NOT `_hpc_array.sh`) |
| `1.1_submit_combinedpbmc.sh` | Conforms (depth 1.1; matches dispatcher glob `1.*_submit_*.sh`) |
| `1.2_submit_joanito_batch_col.sh` | Conforms (depth 1.2) |
| `_create_combinedpbmc_dataset.py` | **Violates** (underscore prefix, no depth number) |
| `_create_joanito_batch_col.R` | **Violates** (underscore prefix, no depth number) |

## Task 1 — Rename files (use `git mv` to preserve history)

| Old | New |
|---|---|
| `src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py` | `src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py` |
| `src/2_dataset_specific_preprocessing/_create_joanito_batch_col.R` | `src/2_dataset_specific_preprocessing/1.2.1_create_joanito_batch_col.R` |

Call chain after rename: `1_submit_hpc.sh` → `1.1_submit_combinedpbmc.sh` → `1.1.1_create_combinedpbmc_dataset.py`; `1_submit_hpc.sh` → `1.2_submit_joanito_batch_col.sh` → `1.2.1_create_joanito_batch_col.R`.

## Task 2 — Update internal references (source files)

| File | Change |
|---|---|
| `src/2_dataset_specific_preprocessing/1.1_submit_combinedpbmc.sh:20` | `"${SCRIPT_DIR}/_create_combinedpbmc_dataset.py"` → `"${SCRIPT_DIR}/1.1.1_create_combinedpbmc_dataset.py"` |
| `src/2_dataset_specific_preprocessing/1.2_submit_joanito_batch_col.sh:19` | `"${SCRIPT_DIR}/_create_joanito_batch_col.R"` → `"${SCRIPT_DIR}/1.2.1_create_joanito_batch_col.R"` |
| `src/1_stage_data/1_stage_data.sh:5` (comment) | `src/2_dataset_specific_preprocessing/_create_combinedpbmc_dataset.py` → `.../1.1.1_create_combinedpbmc_dataset.py` |
| `1.1.1_create_combinedpbmc_dataset.py` docstring (2 occurrences, former lines 5–6) | `python _create_combinedpbmc_dataset.py` → `python 1.1.1_create_combinedpbmc_dataset.py` (both `Local` and `HPC` usage lines) |

No logic changes. No other source references exist (repo-wide grep verified:
only files above, the docs below, and historical `.kilo/plans/*`).

## Task 3 — Update documentation

- `docs/ARCHITECTURE.md`:
  - Line 95: `` `1.1_submit_combinedpbmc.sh` → `_create_combinedpbmc_dataset.py` `` → `` `1.1_submit_combinedpbmc.sh` → `1.1.1_create_combinedpbmc_dataset.py` ``
  - Line 96: `` `1.2_submit_joanito_batch_col.sh` → `_create_joanito_batch_col.R` `` → `` `1.2_submit_joanito_batch_col.sh` → `1.2.1_create_joanito_batch_col.R` ``
  - Optionally append to the folder-2 row (line 74) or bullets a note that processing scripts follow the `N.N.N_<action>` depth convention (mirrors `1.1.1_preprocess.py`/`2.1.1_process_chunk.R`).
- `README.md` lines 106–107: `_create_combinedpbmc_dataset.py`, `_create_joanito_batch_col.R` → new names.
- `AGENTS.md` line 60: `` `1.1_submit_combinedpbmc.sh` → `_create_combinedpbmc_dataset.py`, `1.2_submit_joanito_batch_col.sh` → `_create_joanito_batch_col.R` `` → new names.
- `TODO.md`: update all `_create_combinedpbmc_dataset.py` / `_create_joanito_batch_col.R` references to the new names (lines 5, 18, 25, 69, 99, 389, 441, 451; line 44 references `1.1_submit_combinedpbmc.sh` which is unchanged). Optionally add a checked `[x]` entry noting the rename under the "Completed" section.

Leave `.kilo/plans/*` historical plan files untouched (established precedent).

## Validation

Per AGENTS.md, do NOT run the pipeline scripts. Verify with grep only:

1. `grep -rn '_create_combinedpbmc_dataset\|_create_joanito_batch_col' src/ docs/ README.md AGENTS.md TODO.md` → must return zero matches (only historical hits in `.kilo/plans/` remain, which is expected).
2. `grep -rn '1.1.1_create_combinedpbmc_dataset\|1.2.1_create_joanito_batch_col' src/ docs/ README.md AGENTS.md TODO.md` → matches in the expected files only.
3. `git status` shows the two renames as renames (git detects via `git mv`).

## Out of scope

- No logic/content changes to the renamed scripts (docstring usage lines only).
- No renaming of `1_submit_hpc.sh`, `1.1_submit_combinedpbmc.sh`, `1.2_submit_joanito_batch_col.sh` — they already conform.
- `.kilo/plans/*` historical documents.
