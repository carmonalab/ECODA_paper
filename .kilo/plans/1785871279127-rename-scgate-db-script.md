# Rename 0.1_create_scgate_db.R → 2.0_create_scgate_db.R and move its invocation

## Goal
Rename `src/4_cell_type_annotation/0.1_create_scgate_db.R` to `2.0_create_scgate_db.R` and call it from `2_submit_hpc_array.sh` instead of `1_prepare_chunks.sh`. Update all documentation/comments accordingly.

## Why
The scGate DB cache is only consumed by the annotation array (`2.1.1_process_chunk.R` workers), not by chunk preparation. Running it in `2_submit_hpc_array.sh` moves cache creation closer to where it is used. It also lets `1_prepare_chunks.sh` and `2_submit_hpc_array.sh` run independently.

## Steps

### 1. Rename the R script (git mv)
- `git mv src/4_cell_type_annotation/0.1_create_scgate_db.R src/4_cell_type_annotation/2.0_create_scgate_db.R`
- Update the header comment of the file:
  - `0.1_create_scgate_db.R` → `2.0_create_scgate_db.R`
  - "Called by 1_prepare_chunks.sh (srun compute session, pixi run Rscript --vanilla) BEFORE the annotation array" → "Called by 2_submit_hpc_array.sh (srun compute session, pixi run Rscript --vanilla) BEFORE the annotation array"
  - Keep all logic unchanged (idempotent, `SCGATE_DB_PATH`/`SCGATE_DB_BRANCH` env vars).

### 2. Remove the scGate DB block from `1_prepare_chunks.sh`
- Delete the section "STAGE scGate MODEL DB CACHE ..." (currently lines ~73–96): the comment block plus the `SCRIPT_DIR=...` line (only used by that block), the `if [[ -f "${SCGATE_DB_PATH}" ]] ... else ... fi` block, and its `echo` statements.
- `SCRIPT_DIR` is used later in the file for `1.1_prepare_chunks.py` (line 126) — keep the second definition at the top of the chunk loop. Verify it still exists after removal.
- Note: `1_prepare_chunks.sh` no longer needs `PIXI_RSCRIPT` for this step; nothing else in the file uses it — do NOT remove any `slurm_config.sh` sourcing.

### 3. Add the scGate DB block to `2_submit_hpc_array.sh`
Insert after the dataset-list validation block (after line 31, before the manifest build):
- Add near top (after `cd "${PROJECT_ROOT}"`): `SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"` and `mkdir -p "${LOGS_DIR}"` (the current `mkdir -p "${LOGS_DIR}"` at line 81 can then be dropped or left — prefer moving it up).
- Copy the block from `1_prepare_chunks.sh` (same `srun` pattern on login node: `--partition=shared-cpu`, `--time=00:30:00`, `--ntasks=1`, `--cpus-per-task=1`, `--mem=4G`, log to `${LOGS_DIR}/prepare_scgatedb.log`):
  - Idempotent check `[[ -f "${SCGATE_DB_PATH}" ]]` → skip with message.
  - Run `${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.0_create_scgate_db.R"`.
  - On failure: WARNING + continue (workers still have download+persist fallback in `2.1.1_process_chunk.R`).
- Keep log filename `prepare_scgatedb.log` unchanged (minimal churn).

### 4. Update code comments
- `src/slurm_config.sh` (~line 41): `0.1_create_scgate_db.R (via 1_prepare_chunks.sh)` → `2.0_create_scgate_db.R (via 2_submit_hpc_array.sh)`.
- `src/4_cell_type_annotation/2.1.1_process_chunk.R` (~line 68): comment "(created by 0.1_create_scgate_db.R via ...)" → "(created by 2.0_create_scgate_db.R via 2_submit_hpc_array.sh)".

### 5. Update documentation
- `docs/ARCHITECTURE.md`:
  - Row `1_prepare_chunks.sh` (line 126): remove "creates the scGate model DB cache once (`0.1_create_scgate_db.R` via `srun` → `${SCGATE_DB_PATH}`)".
  - Row `0.1_create_scgate_db.R` (line 128): rename to `2.0_create_scgate_db.R`; "Run from `1_prepare_chunks.sh`" → "Run from `2_submit_hpc_array.sh`".
  - Row `2_submit_hpc_array.sh` (line 129): add that it creates the scGate DB cache once via `srun` before submitting.
  - "NAS ↔ Scratch data flow" bullet (line 139): "`1_prepare_chunks.sh` stages reference maps and the scGate model DB cache" → "`1_prepare_chunks.sh` stages reference maps; `2_submit_hpc_array.sh` creates the scGate model DB cache".
  - "Environment propagation" bullet (line 140): update srun mention to also include `2_submit_hpc_array.sh` (e.g. "through both `srun` (`1_prepare_chunks.sh`, `2_submit_hpc_array.sh`) and `sbatch` (`2_submit_hpc_array.sh`)").
  - Usage section (~line 152): optional note that `2_submit_hpc_array.sh` auto-creates the cache if missing.
- `AGENTS.md`:
  - Line 136: "created by `0.1_create_scgate_db.R` (run via `1_prepare_chunks.sh`)" → "created by `2.0_create_scgate_db.R` (run via `2_submit_hpc_array.sh`)".
  - Lines 72 and 102 (`2_submit_hpc_array.sh` descriptions): mention it creates the scGate model DB cache (if missing) before submitting.
- `TODO.md` line 41 (checked item): `0.1_create_scgate_db.R` (run via srun in `1_prepare_chunks.sh`) → `2.0_create_scgate_db.R` (run via srun in `2_submit_hpc_array.sh`).

## Validation (no pipeline execution per AGENTS.md)
- `git status` / `git diff` to confirm the rename and edits.
- `grep -rn "0\.1_create_scgate_db\|1_prepare_chunks.*scgate\|scgate.*1_prepare" src/ docs/ AGENTS.md TODO.md` → no stale references.
- `grep -rn "2\.0_create_scgate_db" src/ docs/ AGENTS.md TODO.md` → references: `2_submit_hpc_array.sh`, the R file itself, docs.
- `bash -n` on `1_prepare_chunks.sh` and `2_submit_hpc_array.sh` for syntax check (no HPC execution).

## Out of scope
- No change to `datasets.json`, `2.1.1_process_chunk.R` logic, or the worker fallback download.
- No pipeline scripts are run (per AGENTS.md; full validation happens later on HPC with the debugging dataset).
