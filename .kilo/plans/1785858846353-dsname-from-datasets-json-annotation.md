# Derive DS_NAME from datasets.json in the cell type annotation pipeline

## Goal

Remove the manual `export DS_NAME="..."` requirement from the cell type annotation pipeline (`src/4_cell_type_annotation/`). `DS_NAME` must be derived from `datasets.json` keys, mirroring `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, so `1_prepare_chunks.sh` and `2_submit_hpc_array.sh` process all datasets automatically (with an optional single-dataset argument for targeted reruns).

## Critical finding (drives the design)

Bash arrays do **not** propagate through `sbatch` (verified locally: even `declare -ax` arrays never appear in the child's environment; `env` is empty). Therefore:

- `src/3_scrnaseq_preprocessing/1.1_run_worker.sh:28` (`DS_NAME="${DATASET_NAMES[$IDX]}"`) is a latent bug — the array is empty in preprocessing workers.
- The annotation pipeline must NOT replicate this. Workers must derive `DS_NAME` themselves from `datasets.json` (exported via `slurm_config.sh` as `DATASETS_JSON_FILE`, readable from workers since the repo lives on the shared filesystem).

## Decisions (resolved)

1. **Annotation parallelization**: single global SLURM array over ALL chunks of ALL datasets, mapped via a manifest file. One submission, one monitor loop, one rsync.
2. **Targeted runs**: both scripts accept an optional positional `DS_NAME` arg (validated against `datasets.json` keys); default = all datasets.
3. **Preprocessing fix**: include the minimal 2-line fix in `src/3_scrnaseq_preprocessing/1.1_run_worker.sh` (same jq-based derivation), since it is the reference pattern.

## Changes

### 1. `src/4_cell_type_annotation/1_prepare_chunks.sh`

- Replace the `DS_NAME` env requirement (lines 28–32) with:
  - Build `DATASET_NAMES` from `jq -r 'keys[]' "${DATASETS_JSON_FILE}"` (keep the array only for iteration in this script — it is not passed to workers).
  - Optional positional args: `MODE="${1:-production}"` (values `production`/`test`, backward compatible), optional `DS_NAME_ARG="${2:-}"`.
  - If `DS_NAME_ARG` given: validate with `jq -e --arg ds "${DS_NAME_ARG}" 'has($ds)' "${DATASETS_JSON_FILE}"`; error out if not a key. Otherwise iterate over all keys.
- For each dataset: `export DS_NAME=<ds>` (R script reads it via `Sys.getenv`), then run the existing `srun ... 1.1_prepare_chunks.r "${R_PASS_ARG}"` step unchanged (lines 46–60).
- Per-dataset robustness: before srun, check for preprocessed input — `if ! ls "${SCRATCH_OUTPUT_DIR}/${DS_NAME}"/*.h5ad >/dev/null 2>&1` → warn "no preprocessed .h5ad, skipping" and continue (e.g. Zhu has no views). Wrap srun in `if ! srun ...; then` to record failures; continue; at the end print summary and `exit 1` if any dataset failed.
- Keep ref-map staging (lines 38–44) as-is, before the loop.

### 2. `src/4_cell_type_annotation/2_submit_hpc_array.sh`

- Replace the `DS_NAME` env requirement (lines 4–15) with optional positional `DS_NAME_ARG="${1:-}"`, validated against `datasets.json` keys as above. Default = all datasets.
- Build a global manifest:
  - For each dataset (arg or all keys): `NUM_CHUNKS=$(ls -1 "${SCRATCH_OUTPUT_DIR}/${DS_NAME}"/chunks/chunk_*.txt 2>/dev/null | wc -l)`; if 0 → warn and skip (explicit single-dataset arg with 0 chunks = hard error, as today); if 0 across all → error "Run 1_prepare_chunks.sh first".
  - Write `export CHUNKS_MANIFEST="${SCRATCH_OUTPUT_DIR}/chunks_manifest.txt"` (regenerate each run): one line per chunk, tab-separated `DS_NAME<TAB>absolute_chunk_path`.
  - `TOTAL_CHUNKS=$(wc -l < "${CHUNKS_MANIFEST}")`.
- Submit `sbatch --array=1-${TOTAL_CHUNKS}%${MAX_NUM_CHUNKS_PARALLEL} ... 2.1_run_worker.sh` (remove `HOME_CHUNKS_DIR`/`NUM_CHUNKS` logic and the `TISSUE_TYPE`/`NORMAL_TISSUE` export block — moves to the worker, since it is per-task now).
- Keep monitor loop + final rsync (lines 48–66) unchanged.

### 3. `src/4_cell_type_annotation/2.1_run_worker.sh`

- Load jq: `module load jq/1.6` with `command -v jq` fallback error (module loads in the submit script do not propagate to workers).
- Read manifest line: `IFS=$'\t' read -r DS_NAME CHUNK_FILE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHUNKS_MANIFEST}")`; error if empty.
- Export `DS_NAME`, `HOME_CHUNKS_DIR="$(dirname "${CHUNK_FILE}")"`, and per-dataset `TISSUE_TYPE`/`NORMAL_TISSUE` via the same jq one-liners currently in `2_submit_hpc_array.sh` (lines 30–31).
- Remove the `${HOME_CHUNKS_DIR:=${SCRATCH_OUTPUT_DIR}/chunks}` fallback (line 20).
- Call `2.1.1_process_chunk.sh "${SLURM_ARRAY_TASK_ID}" "${CHUNK_FILE}"` as before.

### 4. `src/3_scrnaseq_preprocessing/1.1_run_worker.sh` (minimal fix)

- Replace `DS_NAME="${DATASET_NAMES[$IDX]}"` (line 28) with `DS_NAME="$(jq -r 'keys[]' "${DATASETS_JSON_FILE}" | sed -n "${SLURM_ARRAY_TASK_ID}p")"` (`sed` is 1-based, matching the array indices). Add `module load jq/1.6` + `command -v jq` guard.

### 5. Docs (keep in sync)

- `README.md` lines 132–137: replace `export DS_NAME="Stephenson"` usage with the new interface; mention optional `[DS_NAME]` arg.
- `docs/ARCHITECTURE.md`:
  - Usage block (lines 143–154): drop `export DS_NAME`, document new args.
  - File table rows for `1_prepare_chunks.sh`, `2_submit_hpc_array.sh`, `2.1_run_worker.sh` (lines 124–127).
  - "Environment propagation" note (line 138): `TISSUE_TYPE`/`NORMAL_TISSUE` now exported per-task by the worker.
- `AGENTS.md`: update the sentence about `TISSUE_TYPE`/`NORMAL_TISSUE` being exported by `2_submit_hpc_array.sh` (now by `2.1_run_worker.sh`).
- `src/slurm_config.sh` comment (lines 48–50): `HOME_CHUNKS_DIR` is now worker-internal, not set by callers.

## Unchanged

- `1.1_prepare_chunks.r`, `2.1.1_process_chunk.sh`, `2.1.1.1_process_chunk.R`, `config_helper.R` — all read env vars, which the wrappers now set per dataset.
- `datasets.json` (not touched).
- `AUTHOR_ANNOT_COLNAMES` remains manually exported (documented behavior).

## Risks / notes

- Feather naming stays unique: `annotations_chunk_<SLURM_ARRAY_TASK_ID>.feather` written to per-dataset `path_output`; with a global manifest task IDs are globally unique.
- A stale `chunks_manifest.txt` is impossible to misread because it is regenerated on every run; chunks dirs are recreated fresh by `1.1_prepare_chunks.r`.
- Worker derivation via `sed -n "${SLURM_ARRAY_TASK_ID}p"` requires jq on workers — add the module load + guard.
- If `1_prepare_chunks.sh` is run while `2_submit_hpc_array.sh` is running, the manifest could be regenerated mid-run — acceptable (documented run order: prepare chunks, then submit).

## Validation

Per AGENTS.md rules: no pipeline runs after implementation. Validation is limited to:
- `bash -n` syntax check on all modified bash scripts (local).
- Manual review of arg parsing and manifest format.
- Full HPC validation deferred: run `1_prepare_chunks.sh test` then `2_submit_hpc_array.sh` on the HPC with the small debug dataset once the pipeline is finalized.
