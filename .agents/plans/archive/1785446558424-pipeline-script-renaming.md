# Pipeline Script Renaming Plan

## Goal

Rename pipeline call-chain scripts in `src/preprocess/` and `src/cell_type_annotation/` to use a hierarchical decimal numbering scheme where the depth of the number reflects the call depth.

## Proposed Convention

- Top-level entry point scripts keep a single number: `N_action.sh`
- Directly called scripts get one more decimal: `N.1_action.ext`
- Grandchildren get another decimal: `N.1.1_action.ext`
- And so on.

Helper/utility scripts (underscore-prefixed, e.g. `_preprocess_utils.py`) and standalone notebooks (`.qmd`) are not part of the call chain and keep their names. `DEPRECATED_LEGACY_CODE/` is frozen and left untouched.

## Files to Rename (4 files)

| Directory | Current name | New name | Reason |
|---|---|---|---|
| `src/cell_type_annotation/` | `1_prepare_chunks.r` | `1.1_prepare_chunks.r` | Called by `1_prepare_chunks.sh` (depth 1.1) |
| `src/cell_type_annotation/` | `2.2_process_chunk.sh` | `2.1.1_process_chunk.sh` | Called by `2.1_run_worker.sh`, not a sibling (depth 1.1.1) |
| `src/cell_type_annotation/` | `2.2_process_chunk.R` | `2.1.1.1_process_chunk.R` | Called by `2.1.1_process_chunk.sh` (depth 1.1.1.1) |
| `src/preprocess/` | `1.2_preprocess.py` | `1.1.1_preprocess.py` | Called by `1.1_run_worker.sh`, not a sibling (depth 1.1.1) |

## Files Unchanged

- `1_prepare_chunks.sh`, `2_submit_hpc_array.sh`, `3_merge_annotations.py` — top-level entry points
- `2.1_run_worker.sh` — already correctly numbered as 2.1
- `1_submit_hpc_array.sh`, `1.1_run_worker.sh` — already correctly numbered
- All underscore-prefixed helpers: `_preprocess_utils.py`, `_create_combinedpbmc_dataset.py`, `_create_joanito_batch_col.R`
- All `.qmd` notebooks in `src/preprocess/`
- `DEPRECATED_LEGACY_CODE/` — frozen

## Internal File Reference Updates

After renaming, update references within these 7 files (8 edits):

| File | Line | Old reference | New reference |
|---|---|---|---|
| `src/cell_type_annotation/1_prepare_chunks.sh` | 86 | `1_prepare_chunks.r` | `1.1_prepare_chunks.r` |
| `src/cell_type_annotation/2.1_run_worker.sh` | 24 | `2.2_process_chunk.sh` | `2.1.1_process_chunk.sh` |
| `src/cell_type_annotation/2.1.1_process_chunk.sh` | 3 | `2.2_process_chunk.sh` | `2.1.1_process_chunk.sh` |
| `src/cell_type_annotation/2.1.1_process_chunk.sh` | 33 | `2.2_process_chunk.R` | `2.1.1.1_process_chunk.R` |
| `src/cell_type_annotation/2.1.1.1_process_chunk.R` | 2 | `2.2_process_chunk.R` | `2.1.1.1_process_chunk.R` |
| `src/cell_type_annotation/2.1.1.1_process_chunk.R` | 4 | `2.2_process_chunk.sh` | `2.1.1_process_chunk.sh` |
| `src/preprocess/1.1_run_worker.sh` | 35 | `1.2_preprocess.py` | `1.1.1_preprocess.py` |
| `src/slurm_config.sh` | 32 | `src/preprocess/1.2_preprocess.py` | `src/preprocess/1.1.1_preprocess.py` |

## Documentation Updates

### `docs/ARCHITECTURE.md`
- Line 35: `1_prepare_chunks.r` → `1.1_prepare_chunks.r`
- Line 52-53: `1_prepare_chunks.r` → `1.1_prepare_chunks.r`
- Line 55: `2.2_process_chunk.sh` → `2.1.1_process_chunk.sh`
- Line 56: `2.2_process_chunk.sh` → `2.1.1_process_chunk.sh`; `2.2_process_chunk.R` → `2.1.1.1_process_chunk.R`
- Line 57-58: Full row — rename both file references in the table
- Line 59: `1_prepare_chunks.r` → `1.1_prepare_chunks.r`; `2.2_process_chunk.R` → `2.1.1.1_process_chunk.R`
- Line 80: `3_merge_annotations.py` — unchanged (top-level)
- Line 85: `1_prepare_chunks.sh` — unchanged (top-level)

### `AGENTS.md`
- Line 60: `2.2_process_chunk.R` → `2.1.1.1_process_chunk.R`
- Line 77: `src/preprocess/1.2_preprocess.py` → `src/preprocess/1.1.1_preprocess.py`
- Line 116: `src/preprocess/1.2_preprocess.py` → `src/preprocess/1.1.1_preprocess.py`

### `TODO.md`
- Line 10: `2.2_process_chunk.R` → `2.1.1.1_process_chunk.R`
- Line 343: `2.2_process_chunk.sh` → `2.1.1_process_chunk.sh`; `3_merge_annotations.py` — unchanged

### `README.md` — No changes needed (only top-level script names are referenced, none of which change)

## Execution Order

1. Rename 4 files using `git mv`
2. Update internal references in the 7 source files listed above
3. Update 3 documentation files listed above
4. Verify: `grep -rn '1_prepare_chunks\.r\|2\.2_process_chunk\.sh\|2\.2_process_chunk\.R\|1\.2_preprocess\.py' src/ docs/ AGENTS.md README.md TODO.md` — should return only DEPRECATED_LEGACY_CODE/ and .kilo/plans/ hits (which are intentionally excluded)
