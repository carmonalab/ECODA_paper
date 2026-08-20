# Safeguard: skip re-annotation of already-annotated datasets in `src/4_cell_type_annotation`

## Goal

Prevent the annotation pipeline from unnecessarily re-running for datasets that already have complete annotations, and add a `--force` flag (repo convention: `--force` on the submitter, bypasses the skip) to force recomputation when needed. Scope decided with user: **`1_prepare_chunks.sh` + `3_submit_merge.sh` only** — the annotation array (`2_submit_hpc_array.sh`) stays unchanged (its workers already skip chunks whose feather exists, and its manifest/coverage-gate logic must not be disturbed).

## Background / root cause

- `1.1_prepare_chunks.py` (called via srun from `1_prepare_chunks.sh`, its only caller) **unconditionally** rebuilds `chunks/` + `annotation_union/` and, in production mode, **deletes all `annotations_chunk_*.feather`** (lines 256–305). Deleting the feathers is what forces a full re-annotation by the array.
- `2.1.1_process_chunk.R:145` already skips chunks whose feather exists — the array is safe; the damage happens upstream in chunk prep.
- `3_submit_merge.sh` re-merges feathers into every view h5ad (atomic temp+`os.replace` rewrite of multi-GB files) and re-syncs to NAS every time it is invoked — wasteful when the dataset is already merged.
- No `--force` exists anywhere in `4_cell_type_annotation` (verified by grep). Convention from other pipelines (`FORCE_PREPROCESS`, `FORCE_BENCHMARK`): parse `--force`, guard `--force` + `--sync-only` combination where relevant.

## "Already annotated" predicate

A dataset is considered done (skip, unless `--force`) if either:

- **Branch 1 — annotated, not yet merged**: `${HPC_SCRATCH_DIR}/<DS>/output/chunks/` contains ≥1 `chunk_*.txt` **and** every `chunk_N.txt` has its matching `annotations_chunk_N.feather` in `${HPC_SCRATCH_DIR}/<DS>/output/` (mapping: `chunk_3.txt` → `annotations_chunk_3.feather`).
- **Branch 2 — already merged** (merge is the only step that deletes these dirs): `output/chunks/` absent **and** `annotation_union/` absent **and** ≥1 `annotations_chunk_*.feather` exists in `output/`.

Anything else (no chunks/no feathers → fresh; partial feathers → incomplete) proceeds with the current rebuild behavior. Partial coverage correctly triggers rebuild + re-annotation (stale annotation artifacts are the "not done" state).

## Changes

### Task 1 — `src/4_cell_type_annotation/1_prepare_chunks.sh`

1. **Argument parsing**: add optional `--force` flag accepted in any position alongside the existing positional `MODE` (`production`/`test`) and optional `DS_NAME` (`./1_prepare_chunks.sh production Kfoury --force`, `./1_prepare_chunks.sh --force`, etc.). Keep current defaults (`production`, all datasets). Error on unknown flags. There is no `--sync-only` here, so no combination guard needed.
2. **Per-dataset skip check** (inside the existing `for DS_NAME in ...` loop, *before* the h5ad-existence check and the srun):
   - Evaluate the predicate above for `${HPC_SCRATCH_DIR}/${DS_NAME}` (use `shopt -s nullglob` for the globs; mind the existing `shopt` usage).
   - If annotated and no `--force`: print `Already annotated: ${DS_NAME} — skipping chunk generation (use --force to recompute)` and `continue` (add to a new `SKIPPED_ANNOTATED` list).
   - If `--force`: print a note that the run is forced.
3. **Summary**: extend the end summary to report `Skipped (already annotated):` counts/names separately from the existing `Skipped (no preprocessed .h5ad)` category.
4. Skip applies in **both** `production` and `test` modes (test mode never deletes feathers anyway; `--force` restores the old behavior).

No changes to `1.1_prepare_chunks.py` (single caller, keeps the destructive step only when the skip did not trigger or `--force` was given).

### Task 2 — `src/4_cell_type_annotation/3_submit_merge.sh`

1. **Argument parsing**: add optional `--force` flag (`./3_submit_merge.sh <DS_NAME> [--force]`); `DS_NAME` stays the first positional.
2. **Post-merge skip** (immediately after the `datasets.json` DS validation, *before* the view/feather globs and the coverage gate):
   - If `output/chunks/` absent **and** `annotation_union/` absent **and** ≥1 `annotations_chunk_*.feather` in `output/` and no `--force`: print `Already merged: ${DS_NAME} — skipping merge + NAS sync (use --force to re-merge)` and `exit 0` (no srun, no NAS sync, no sync-status email).
   - `--force` falls through to the existing logic unchanged — including the fail-closed coverage gate (do NOT bypass it; the gate protects NAS from partial merges).
3. Update the usage header comment.

### Task 3 — Docs

- Update `docs/ARCHITECTURE.md` (annotation pipeline table rows for `1_prepare_chunks.sh` and `3_submit_merge.sh`) to document the skip predicate and `--force`.
- Update the usage header comments of both scripts (Task 1/2 include this).
- No changes needed in AGENTS.md / README.md (high-level, unaffected).

## Behavior matrix (after implementation)

| State | Command | Outcome |
|---|---|---|
| Fresh (no chunks, no feathers) | `1_prepare_chunks.sh` | rebuilds (unchanged) |
| Annotated, not merged (chunks + complete feathers) | `1_prepare_chunks.sh` | **skip** |
| Merged (chunks/union gone, feathers kept) | `1_prepare_chunks.sh` | **skip** (branch 2) |
| Any done state | `1_prepare_chunks.sh --force` | rebuild + (production) delete feathers → re-annotation |
| Merged | `3_submit_merge.sh <DS>` | **skip** |
| Merged | `3_submit_merge.sh <DS> --force` | re-merge (gate still applies) |
| Annotated, not merged | `3_submit_merge.sh <DS>` | merges (unchanged) |
| Partial feathers | `1_prepare_chunks.sh` | rebuild + re-annotate (unchanged) |

## Edge cases / caveats

- **Test-mode leftover chunks** (chunk_size=1, 2-line chunk files) with complete feathers trigger branch 1 and skip a later production run. Not harmful (annotation results identical; array would skip anyway); document in the ARCHITECTURE.md row that `--force` is the escape hatch after test-mode runs.
- **Manual deletion of `chunks/` only** (union still present): branch 2 requires the union to be absent → rebuild proceeds → correct.
- **`--force` merge on a post-merge dataset** whose chunks were never re-added to `chunks_manifest.txt` (e.g., an array run after the merge excluded the DS): the coverage gate still fails closed — expected; message instructs to re-run prepare + array first.
- **`_debug`**: default-all runs of `1_prepare_chunks.sh` still include it (unchanged); once annotated it is now skipped like any other dataset. The `_*` default-all filtering convention mismatch with preprocessing is out of scope.

## Out of scope (explicitly)

- Changes to `2_submit_hpc_array.sh` / `2.1.1_process_chunk.R` (array-level skip + worker `--force`; would require adapting the merge coverage gate).
- Changes to `1.1_prepare_chunks.py` (`--force` arg / internal skip) — bash-only check since it has a single caller.
- `_*` key filtering in annotation default-all loops.

## Validation

Per AGENTS.md: do not run pipeline scripts for validation until the pipeline is fully implemented; then validate with the `_debug` debugging dataset. Intermediate checks that require no compute:

1. `bash -n` on both modified scripts.
2. Login node (no srun triggered): `1_prepare_chunks.sh production Kfoury` → expect `Already annotated: Kfoury — skipping` and no srun allocation (works whether Kfoury is in the annotated or merged state).
3. Login node: `1_prepare_chunks.sh production Kfoury --force` → expect the forced-rebuild note (defer the actual srun to `_debug` validation).
4. Login node: `3_submit_merge.sh Kfoury` → expect `Already merged: Kfoury — skipping` and exit 0.
5. Simulate partial coverage on the login node (temporarily rename one Kfoury feather): `1_prepare_chunks.sh production Kfoury` → expect the rebuild path (no skip); restore the feather.
6. Full validation on `_debug` once implemented (per AGENTS.md): prepare → skip as no-op; `--force` flow → rebuild → array → merge → verify annotation columns in the merged h5ad.
