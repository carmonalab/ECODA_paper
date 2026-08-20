# Collapse cell-type-annotation wrapper chain (3 levels, mirroring stage 3)

## Context

`src/4_cell_type_annotation/` currently has a 4-level execution chain:

```
2_submit_hpc_array.sh  →  2.1_run_worker.sh  →  2.1.1_process_chunk.sh  →  2.1.1.1_process_chunk.R
(login node: submit/     (compute node:        (thin bash wrapper:      (core annotation logic)
 monitor/rsync)          manifest lookup,       logging, chunk check,
                         DS_NAME + TISSUE_      pixi Rscript call)
                         TYPE/NORMAL_TISSUE
                         exports)
```

A reviewer suggested collapsing to 2 scripts. **2 is not achievable**: the submit
script must stay on the login node (sbatch submission, completion monitoring,
NAS rsync — compute nodes have no NAS access), so the minimum is 3 levels.
The correct target is the 3-level pattern already established by the sibling
pipeline `src/3_scrnaseq_preprocessing/` (`1_submit_hpc_array.sh` →
`1.1_run_worker.sh` → `1.1.1_preprocess.py`). The 4th level here is inherited
legacy (TODO.md line ~400) and is genuinely redundant:

- `2.1.1_process_chunk.sh` only adds: a timestamped log helper, a chunk-file
  existence check, `cd "${PROJECT_ROOT}"`, and the pixi invocation.
- `2.1_run_worker.sh` already sources `slurm_config.sh`, `cd`s to
  `PROJECT_ROOT`, and has `SCRIPT_DIR` (2.1_run_worker.sh:15-17).

Naming deviation from the reviewer: keep `2_submit_hpc_array.sh` (renumbering
to `1_` would collide with `1_prepare_chunks.sh`/`1.1_prepare_chunks.py` in the
same directory) and keep the `X.X_run_worker → X.X.X_core` convention.

## Changes

### 1. Merge `2.1.1_process_chunk.sh` into `2.1_run_worker.sh`
- Keep all existing worker logic (SBATCH header, slurm_config source, jq load,
  manifest lookup via `sed -n ${SLURM_ARRAY_TASK_ID}p`, `DS_NAME`/`CHUNK_FILE`
  split, `TISSUE_TYPE`/`NORMAL_TISSUE` exports).
- Add after the exports:
  - chunk-file existence check (`[[ -f "${CHUNK_FILE}" ]]` → error exit, as in
    the old wrapper lines 24-27),
  - timestamped `log_msg()` helper (optional; plain `echo` with task id is fine),
  - pixi invocation:
    `"${HOME}/.pixi/bin/pixi" run Rscript --vanilla "${SCRIPT_DIR}/2.1.1_process_chunk.R" "${CHUNK_FILE}"`.
- Remove the final `bash ... 2.1.1_process_chunk.sh` call and the stale comment
  at line 49 ("PROJECT_ROOT is sourced (and exported) by 2.1.1_process_chunk.sh
  ..." — the worker already sources slurm_config.sh itself).

### 2. Rename `2.1.1.1_process_chunk.R` → `2.1.1_process_chunk.R`
- `git mv`; update header comment ("Called by 2.1_run_worker.sh (pixi run
  Rscript --vanilla)").

### 3. Simplify R arg passing (drop legacy `chunk_file__` encoding)
In `2.1.1_process_chunk.R`, replace the `key__value` parsing block
(lines ~70-85: `raw_args`/`parsed_args_list`/`modifyList` machinery) with:
```r
raw_args <- commandArgs(trailingOnly = TRUE)
args <- defaults
if (length(raw_args) > 0) args$chunk_file <- raw_args[1]
```
Keep the existing `file.exists()` guard (lines 83-85) unchanged. All other
params already come from env vars (`defaults` block, lines 51-62). Only the
worker calls this script, and HPC chunk paths never contain `__`, so this is
safe. (If a minimal-diff refactor is preferred, this step is optional — the
merge works without it.)

### 4. Delete `2.1.1_process_chunk.sh`
After step 1, nothing references it (verified: only `2.1_run_worker.sh` calls it).

### 5. Update docs
- `docs/ARCHITECTURE.md`:
  - workflow diagram (line ~112) — worker label unchanged (`2.1_run_worker.sh`),
  - file table (lines ~127-132): remove the `2.1.1_process_chunk.sh` row; update
    `2.1_run_worker.sh` role ("...then calls `2.1.1_process_chunk.R` directly via
    `pixi run Rscript --vanilla`"); rename `2.1.1.1_process_chunk.R` →
    `2.1.1_process_chunk.R`;
  - `config_helper.R` row (line ~132): "Called by `2.1.1_process_chunk.R`";
  - env-propagation paragraph (line ~139) — no structural change, but remove any
    `2.1.1_process_chunk.sh` mention.
- `AGENTS.md` line ~129: `2.1.1.1_process_chunk.R` → `2.1.1_process_chunk.R`.
- `TODO.md`:
  - open item (lines 48-51, "Check if cell type annotation pipeline can be
    simplified" / "*_process_chunk.sh" check): mark addressed — wrapper collapsed
    to 3 levels, arg encoding removed;
  - live filename references to `2.1.1.1_process_chunk.R` (line 50) → new name.
  - Historical changelog entries (lines ~15, ~400, ~413, ~420, ~423) stay as-is.
- `README.md`: no change (only references `2_submit_hpc_array.sh`).

## Validation
- `bash -n` on the modified `2.1_run_worker.sh`.
- Grep repo for stale references: `2.1.1_process_chunk.sh`, `2.1.1.1_process_chunk` —
  only expected hits: git history / historical TODO.md entries.
- Do **not** run the pipeline (per AGENTS.md rule).

## Risks
- Low. Single caller of the removed wrapper; `SLURM_ARRAY_TASK_ID` env propagates
  through `pixi run` unchanged (feather naming in R depends on it).
- Doc drift: grep for old filenames after edits (covered in validation).

## Out of scope
- `1_prepare_chunks.sh`/`1.1_prepare_chunks.py` (prepare step stays a separate
  entry point; it runs before submission, not on workers).
- Unifying `2_submit_hpc_array.sh` with stage 3's submit script (TODO.md
  "Other major goals") — separate open goal.
- `config_helper.R` vs `slurm_config.sh` question (TODO.md line 49) — unchanged.
