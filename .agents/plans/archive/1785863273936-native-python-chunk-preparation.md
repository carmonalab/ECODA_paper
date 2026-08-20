# Reviewer point 2: Replace reticulate-based chunk preparation with native Python

## Goal

Replace `src/4_cell_type_annotation/1.1_prepare_chunks.r` (R script that uses reticulate to read `.h5ad` files via Python's anndata) with a native Python script. The chunk files it produces are consumed downstream by `2_submit_hpc_array.sh` → `2.1_run_worker.sh` → `2.1.1.1_process_chunk.R`, so the on-disk chunk format **must stay byte-compatible**:

- `chunk_N.txt` per file: line 1 = absolute `.h5ad` path, subsequent lines = sample IDs.
- Chunks dir: `${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks`, deleted + recreated fresh on each run.
- Chunk size: 5 (production) / 1 (`test` mode).
- One globally increasing `chunk_N` counter across all `.h5ad` files of the dataset.

## Current flow (context)

- `1_prepare_chunks.sh` loops over datasets, per dataset `export DS_NAME` then `srun --partition=shared-cpu --time=00:30:00 --ntasks=1 --cpus-per-task=1 --mem=4G "${PROJECT_ROOT}/.pixi/envs/default/bin/Rscript" --vanilla 1.1_prepare_chunks.r "test__false|test__true"`.
- `1.1_prepare_chunks.r` reads `DS_NAME`/`SAMPLE_COLNAME` env vars, sources `config_helper.R`, derives `path_data = ${SCRATCH_OUTPUT_DIR}/${DS_NAME}` and `path_output_chunks = ${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks`, reads each `*.h5ad` in backed mode via reticulate, `as.character(unique(obs[[SAMPLE_COLNAME]]))`, chunks consecutively, writes files.
- `config_helper.R` remains R-only and is still needed by `2.1.1.1_process_chunk.R` — it is NOT touched by this change (except doc references).

## Changes

### 1. Create `src/4_cell_type_annotation/1.1_prepare_chunks.py`

Native Python replacement. Follow conventions of `1.1.1_preprocess.py` (argparse CLI; env vars via `os.environ`).

- **CLI**: `argparse` with `--test` flag (production default). Replaces the R packed-arg convention `test__false`/`test__true`.
- **Env guards** (mirror the R script's `stop()` guards, exit non-zero with clear message):
  - `PROJECT_ROOT` (sanity check that `slurm_config.sh` was sourced)
  - `DS_NAME`
  - `SAMPLE_COLNAME` (defaults to nothing; must be set — do NOT default to `"Sample"`, keep R's strict guard)
  - `SCRATCH_OUTPUT_DIR`, with fallback `os.path.join(HPC_SCRATCH_DIR, "output")` (same fallback as `config_helper.R`)
- **Paths**: `path_data = SCRATCH_OUTPUT_DIR/DS_NAME`; `path_output_chunks = path_data/chunks`.
- **Clean start**: `shutil.rmtree(path_output_chunks, ignore_errors=True)` then `mkdir`.
- **h5ad scan**: `sorted(Path(path_data).glob("*.h5ad"))` (R `list.files` default `sorted=TRUE` — keep ordering deterministic).
- **Per file** (behavior parity with R):
  - `ad.read_h5ad(f, backed="r")` (anndata 0.12.19 is already in pixi deps).
  - Sample IDs: `adata.obs[sample_col].astype(str).unique()` — pandas `unique()` preserves first-appearance order like R's `unique()`; `astype(str)` mirrors R's `as.character()`.
  - If sample column missing → raise with clear error (parity with R's failure mode).
  - Chunk consecutive samples: `[samples[i:i+chunk_size] for i in range(0, len(samples), chunk_size)]` (parity with R `split(..., ceiling(seq_along()/n))`).
  - Write each group as `chunk_<global_counter>.txt`, first line = absolute h5ad path (as `str(Path.resolve())` or `os.path.abspath`), then sample IDs; `chunk_size = 1` if `--test` else `5`.
- **Messages**: print the same progress info as the R script (dataset path, files found, per-file processing, total chunk count) so logs remain comparable.
- **Exit codes**: unhandled exception → non-zero, caught by `set -euo pipefail` in the bash wrapper.

### 2. Update `src/4_cell_type_annotation/1_prepare_chunks.sh`

- Replace `ENV_RSCRIPT="${PROJECT_ROOT}/.pixi/envs/default/bin/Rscript"` (line 73) with `ENV_PYTHON="${PROJECT_ROOT}/.pixi/envs/default/bin/python"` (same pixi-env pattern; anndata already a pixi dependency).
- Remove `export R_LIBS_SITE="${PIXI_R_LIB}:${R_LIBS_SITE:-}"` (line 76) — R-specific, no R invoked by this script anymore (`PIXI_R_LIB` stays exported in `slurm_config.sh` for the annotation workers).
- Replace `R_PASS_ARG` (lines 22–26) with a Python-args variable: `PY_ARGS=""` (production), `PY_ARGS="--test"` (test mode).
- Update the `srun` call (line 106) to:
  `"${ENV_PYTHON}" "${SCRIPT_DIR}/1.1_prepare_chunks.py" ${PY_ARGS}` (unquoted expansion is safe here: `PY_ARGS` is always set, either empty or `--test`).
- Keep `srun` resources (4G / 30min / 1 cpu) — native Python backed read is lighter than R+reticulate; no resource change needed.

### 3. Delete `src/4_cell_type_annotation/1.1_prepare_chunks.r`

Fully superseded; grep confirms only `1_prepare_chunks.sh`, `2_submit_hpc_array.sh` (comment), `docs/ARCHITECTURE.md`, and `TODO.md` reference it.

### 4. Update `src/4_cell_type_annotation/2_submit_hpc_array.sh`

Line 39 comment: `1.1_prepare_chunks.r` → `1.1_prepare_chunks.py`.

### 5. Update `docs/ARCHITECTURE.md`

- Line 107 diagram: `(1_prepare_chunks.sh + 1.1_prepare_chunks.r)` → `1.1_prepare_chunks.py`.
- Line 124 (`1_prepare_chunks.sh` row): "calling `1.1_prepare_chunks.r` via `srun`" → "calling `1.1_prepare_chunks.py` (pixi python) via `srun` (4G, 30min) per dataset".
- Line 125: replace the `1.1_prepare_chunks.r` row with `1.1_prepare_chunks.py`: "Native Python (anndata, no reticulate): reads each .h5ad in backed mode, extracts unique sample IDs from `sample_col` (env var `SAMPLE_COLNAME`), groups them into chunks of 5 (or 1 in test mode), writes `chunk_N.txt` files (1st line = h5ad path, subsequent lines = sample IDs)."
- Line 131 (`config_helper.R` row): "Called by both `1.1_prepare_chunks.r` and `2.1.1.1_process_chunk.R`" → "Called by `2.1.1.1_process_chunk.R`".
- Line 138 (environment propagation): adjust "so they reach R via `Sys.getenv()`" to mention both R (`Sys.getenv()`) and Python (`os.environ`) consumers.

### 6. Update `TODO.md`

Add a completed entry under `## Open reviewer points (to be addressed)` (after line 15's entry), e.g.:
`- [x] **Reticulate chunk preparation removed**: 1.1_prepare_chunks.r replaced by native Python 1.1_prepare_chunks.py (anndata backed read, --test flag, same chunk format); 1_prepare_chunks.sh now calls the pixi python binary; docs updated.`

## Validation

Per AGENTS.md general rules: do NOT run pipeline scripts for validation; full HPC validation happens later with the debug (Joanito-derived) dataset. Only:

1. `python -m py_compile src/4_cell_type_annotation/1.1_prepare_chunks.py` (local, syntax check).
2. `bash -n src/4_cell_type_annotation/1_prepare_chunks.sh` (syntax check, no execution).
3. `bash -n src/4_cell_type_annotation/2_submit_hpc_array.sh` (syntax check).
4. Grep to confirm no remaining references to `1.1_prepare_chunks.r`.

## Risks / notes

- **Chunk-format compatibility** is the main risk; the new script must reproduce the R output exactly (sorted file order, first-appearance sample order, line 1 = h5ad path, consecutive sample grouping). A stale manifest cannot be misread because `2_submit_hpc_array.sh` regenerates it and the chunks dir is recreated fresh.
- **NA sample values**: R `as.character(NA)` yields `"NA"`; pandas `astype(str)` yields `"nan"`. Sample columns are expected clean post-preprocessing (`1.1.1_preprocess.py` standardizes them), so no special handling — behavior parity on valid data.
- `config_helper.R`, `slurm_config.sh`, `2.1_run_worker.sh`, `2.1.1_process_chunk.sh`, `2.1.1.1_process_chunk.R`, `3_merge_annotations.py` are unchanged (the last four still need R).
