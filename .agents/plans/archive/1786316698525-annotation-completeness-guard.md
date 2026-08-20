# Phase 5 — Annotation completeness guard

Fix the partial-preprocessing race in the cell-type annotation pipeline (TODO.md Phase 5, items 1 + 2).

## Problem

`1_prepare_chunks.sh` decides "preprocessing done?" only via `ls output/*.h5ad`. For multi-view datasets (Stephenson, `_debug`; CombinedPBMC is single-view) a still-running preprocess task can have written only 1 of 2 views, so the annotation union/chunks are built from an incomplete view set. The dataset is then marked "already annotated" and skipped forever by later runs — silent incomplete annotations; only `--force` repairs. Additionally, re-running `2_submit_hpc_array.sh` on annotated-but-not-yet-merged datasets redundantly re-annotates every chunk (no skip logic).

## Decisions (confirmed with user)

- **Incomplete preprocessing → WARNING-skip** (loud message listing missing files, listed in the summary, exit 0; dataset stays "not annotated" and is picked up on re-run).
- **Completeness signal = view-completeness check** (all expected `views.*.output_file_name` from datasets.json exist in `output/`), NOT a preprocess completion marker (stale-marker risk on `--force` re-runs; would touch preprocess workers).
- **Include the worker-side per-chunk skip** in `2.1_run_worker.sh`; the chunk manifest stays unfiltered so the `3_submit_merge.sh` coverage gate is unchanged.

## Changes

### 1. `src/4_cell_type_annotation/1_prepare_chunks.sh`

Replace the readiness check at lines 235–240 (inside the per-dataset loop, after the ANNOTATED skip):

- Compute expected view outputs, mirroring `1.1.1_preprocess.py:226-235` skip semantics:
  `jq -r --arg ds "${DS_NAME}" '.[$ds].views | to_entries[] | select(.value.input_file_name != null) | select(.value.output_file_name != null) | .value.output_file_name' "${DATASETS_JSON_FILE}"` (while-read loop; no filenames contain whitespace).
- **Expected set empty** (e.g. Zhu, `"views": {}`) → keep the existing `ls output/*.h5ad` check and its "No preprocessed .h5ad files ... skipping" WARNING as-is (Zhu path unchanged).
- **Some expected files missing** → loud WARNING: "preprocessing incomplete for ${DS_NAME}: missing view file(s): <list> — re-run this script after the preprocess array finishes", append to a new `SKIPPED_INCOMPLETE` array, `continue` (no srun slot wasted).
- **All expected files present** → proceed as today.
- Summary block: add "Skipped (preprocessing incomplete): ${#SKIPPED_INCOMPLETE[@]}" with the names; exit code unchanged (0 unless FAILED_DATASETS non-empty).

### 2. `src/4_cell_type_annotation/1.1_prepare_chunks.py`

Defensive fail-closed check in `main()` right after the existing "No preprocessed .h5ad files" exit (lines 252–254), before union building:

- Read `DATASETS_JSON_FILE` from env (exported by `slurm_config.sh`); if unset, warn and skip the check (never hit on HPC).
- `expected = {v["output_file_name"] for v in ds_entry["views"].values() if v.get("input_file_name") is not None and v.get("output_file_name")}`.
- `missing = expected - {f.name for f in h5ad_files}` (glob already excludes `*_raw.h5ad`); if non-empty → `sys.exit("CRITICAL Error: preprocessing incomplete for <ds>: missing expected view file(s) <sorted(missing)> in <path_data> (run the preprocess array first).")` — catches bypasses/drift of the bash check (dataset lands in FAILED_DATASETS, fail-closed).

### 3. `src/4_cell_type_annotation/2.1_run_worker.sh`

After the chunk-file existence check (lines 49–52), add a pure-bash skip (no jq):

- Derive `CHUNK_NUM` from `basename "${CHUNK_FILE}"` (`chunk_<N>.txt` → N), feather path `${HPC_SCRATCH_DIR}/${DS_NAME}/output/annotations_chunk_${CHUNK_NUM}.feather` (same `chunk_N.txt → annotations_chunk_N.feather` mapping used by the merge).
- If it exists → `echo "Task ${SLURM_ARRAY_TASK_ID}: ${FEATHER_FILE} already exists — annotation already done, skipping."` and `exit 0` (task counts as COMPLETED in sacct).
- Safe against stale feathers: `1.1_prepare_chunks.py` deletes feathers on every chunk rebuild (production), so an existing feather always matches the current chunk set.
- Manifest untouched → `3_submit_merge.sh` coverage gate (manifest lines vs feather count) unchanged; partial-failure re-runs still converge (failed chunks re-annotate, done chunks skip).

### 4. Docs + wrap-up

- `docs/ARCHITECTURE.md`:
  - Row `1_prepare_chunks.sh` (line 146): document the preprocessing-completeness predicate (all expected view h5ads from datasets.json must exist; missing → WARNING-skip listing the files, picked up on re-run; no-views datasets still skip).
  - Row `1.1_prepare_chunks.py` (line 147): document the defensive fail-closed expected-view check.
  - Row `2_submit_hpc_array.sh` (line 149): document that workers skip chunks whose feather already exists (re-runs before merge do no redundant annotation; manifest unfiltered, coverage gate unchanged).
- `TODO.md`: tick the two Phase 5 checkboxes.
- Per AGENTS.md Task Completion Workflow: move this plan to `.kilo/plans/archive/`, `git add .`, commit, push.

## Validation (no pipeline runs; AGENTS.md rule)

- `bash -n` on the two modified bash scripts.
- `py_compile` on `1.1_prepare_chunks.py`.
- jq sanity on datasets.json: Stephenson → 2 expected files, Zhu → 0.
- HPC checks for the user (when convenient): (a) 1_prepare_chunks default-all while a multi-view dataset is mid-preprocess → WARNING-skip with missing files, exit 0; (b) re-run 2_submit on an annotated-but-unmerged dataset → tasks exit fast with the skip message, sacct COMPLETED; (c) merge coverage gate passes unchanged.

## Failure modes / edge cases

- Dataset mid-preprocess with 0 views written → expected non-empty, all missing → WARNING-skip (message lists all).
- View file corrupt/mid-write → python read fails → FAILED_DATASETS (fail-closed, same as today).
- `--force` preprocess re-run in progress → same guard applies (no stale-marker problem — why the marker approach was rejected).
- 2_submit on a fully-annotated dataset before merge → array slots still allocated (manifest size unchanged) but all tasks exit 0 instantly; no redundant compute.

## Out of scope

- Preprocess worker changes (no marker), `datasets.json` changes, `3_submit_merge.sh` changes.
- The negligible NAS-sync race between the preprocess tail and the merge rsync (documented in session chat; not part of Phase 5).
