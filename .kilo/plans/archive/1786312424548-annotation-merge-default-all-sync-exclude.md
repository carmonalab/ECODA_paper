# Annotation follow-up: default-all merge mode + keep NAS sync free of annotation intermediates

## Goal

Follow-up to the annotation skip safeguard (archived plan `1786310705787`). Three fixes + docs:

1. Add a **default-all mode to `3_submit_merge.sh`** (bare invocation = all datasets), making README.md's "running the pipeline" block correct (today a bare `3_submit_merge.sh` errors: `ERROR: Usage: ./3_submit_merge.sh <DS_NAME> [--force]`). Scope decided with user: **fix**, not document-only.
2. **Exclude annotation intermediates from the array's NAS sync** (`2_submit_hpc_array.sh` default-all sync tail rsyncs every `*/output` dir, re-uploading the kept local feathers of already-merged datasets; the merge skip no longer cleans them, so they accumulate on NAS). Scope decided with user: **exclude from array sync** (root cause), not clean-in-skip / document-only.
3. Docs: README.md + ARCHITECTURE.md (usage block + affected table rows).

Out of scope (unchanged): the annotation array (`2_submit_hpc_array.sh` submission/manifest/coverage-gate logic), `2.1.1_process_chunk.R`, `1.1_prepare_chunks.py`, `_debug` default-all filtering.

## Task 1 — `src/4_cell_type_annotation/3_submit_merge.sh`: default-all mode

**Behavior contract** (decided with user — "continue + summary email"):

- `./3_submit_merge.sh <DS_NAME> [--force]` — **single-dataset mode, unchanged byte-for-byte**: skip exit 0, per-event `notify_sync_status` emails, exit codes.
- `./3_submit_merge.sh [--force]` (no positional) — **default-all mode** over all `jq -r 'keys[]' datasets.json` keys (**includes `_debug`**, matching the annotation default-all convention; a fresh `_debug` is WARNING-skipped like any no-feather dataset):
  - already merged → print `Already merged: <DS> — skipping merge + NAS sync (use --force to re-merge)` + continue;
  - no `annotations_chunk_*.feather` in `output/` → `WARNING: No annotations_chunk_*.feather files ... skipping <DS>` + continue (mirrors `1_prepare_chunks.sh`'s no-h5ad skip);
  - merge failure (coverage gate, view merge error) → print error, add to `FAILED_DATASETS`, continue;
  - **NAS reachability is hoisted before the loop** (one `ls "${NAS_TARGET_DIR}/.."` check, fail fast with one email + exit 1 — every dataset needs NAS);
  - after the loop: **one summary email** (all-merged list / failed list) and `exit 1` if `FAILED_DATASETS` non-empty, else `exit 0`.

**Implementation notes**

- Restructure the per-dataset flow into a bash function (e.g. `merge_one_ds()`) so it can `return` per dataset instead of `exit`; the current `exit 1`/`exit 0` calls become returns in default-all mode. Guard the per-event `notify_sync_status` calls inside it with a mode flag (e.g. `EMAIL_MODE=per-event` for single-dataset, `summary` for default-all) so single-dataset behavior is preserved exactly; the summary email is sent only by the default-all caller.
- Keep the existing `while [[ $# -gt 0 ]]` case-loop parsing (`--force`, unknown-flag error) as-is; `DS_NAME="${POS_ARGS[0]:-}"` empty → default-all mode.
- Keep the existing `shopt -s/-u nullglob` pairs intact around each glob.
- `--force` in default-all mode: non-merged datasets re-merge; already-merged datasets have no chunks, so the coverage gate fails closed per dataset (`EXPECTED_CHUNKS=0`, `ANNOT_FILES≥0` → collected into `FAILED_DATASETS` with the existing "Re-run 2_submit_hpc_array.sh ... first" error). This is the documented expected behavior — do not special-case it.
- Update the header usage comment: bare invocation = all datasets; note skip behavior and `--force` gate semantics for merged datasets.

## Task 2 — `src/4_cell_type_annotation/2_submit_hpc_array.sh`: exclude intermediates from NAS sync

- In the post-pipeline sync loop (`SYNC_DIRS=...; for DS_DIR in ...; rsync -rlptDv "${DS_DIR}/" "${NAS_TARGET_DIR}/${DS_NAME}/output/"`, ~lines 229-235), add to the rsync:
  `--exclude='annotations_chunk_*.feather' --exclude='chunks/'`
- Rationale (document in the script comment): nothing consumes feathers/chunks from NAS — the merge reads local scratch only; the merge script's post-sync NAS cleanup (`3_submit_merge.sh` lines ~146-151) proves they are unwanted there. With the merge skip, that cleanup no longer fires for merged datasets, so the exclude is the root-cause fix.
- Do NOT touch the submission/manifest/sacct-gate logic. The rsync has no `--delete`, so pre-existing NAS artifacts are not removed by this change.
- `3_submit_merge.sh`'s NAS feather/chunks cleanup becomes a guarded no-op — leave it as-is (its glob-count guard already handles absence); no churn.

## Task 3 — optional: clearer coverage-gate message for the merged+`--force` case

- In the coverage-mismatch error block of `3_submit_merge.sh` (the `Re-run 2_submit_hpc_array.sh ${DS_NAME} first.` message), when `output/chunks/` is absent, extend the message to also point at chunk preparation: e.g. `Dataset has no chunks (already merged or never prepared); run 1_prepare_chunks.sh ${DS_NAME} --force, then 2_submit_hpc_array.sh ${DS_NAME} first.`
- Gate logic stays exactly as-is (still fail closed). Marked optional — skip if churn concerns; the ARCHITECTURE row already documents the behavior.

## Task 4 — Docs

- **README.md:115**: keep the bare `./src/4_cell_type_annotation/3_submit_merge.sh` line (now correct) and extend the comment: `# merge annotations back to h5ad + sync to NAS (all datasets; per-dataset: ./3_submit_merge.sh <DS> [--force])`.
- **docs/ARCHITECTURE.md**:
  - `3_submit_merge.sh` table row: add default-all mode sentence (bare invocation = all datasets, `_debug` included, no-feather datasets WARNING-skipped, failures collected with one summary email + exit 1, NAS checked once up front; single-dataset mode unchanged; `--force` default-all fails the gate per already-merged dataset — expected).
  - `2_submit_hpc_array.sh` table row: add that the sync tail rsyncs exclude `annotations_chunk_*.feather` and `chunks/` (intermediates never consumed from NAS; previously cleaned up by the merge script after its sync).
  - Usage block (annotation section): update step 3 to show bare (all) + per-dataset + `--force` lines with a one-line explanation.
- No AGENTS.md changes (its `3_submit_merge.sh <DS>` mention stays valid).

## Edge cases / caveats

- **Default-all merge after a partial cycle**: datasets skipped by prepare (already annotated) are WARNING-skipped by the array (no chunks → not in manifest) and merged normally by the merge loop (chunks exist → gate matches feathers) — unchanged behavior, now one command.
- **Merged datasets + default-all `--force`**: gate fail-closed per dataset, collected as failure with the (improved) message — expected, not a bug.
- **`_debug`**: included in default-all merge (annotation convention, no `_` filter); fresh `_debug` → WARNING-skip; processed `_debug` → merges like any dataset.
- **NAS-unreachable in default-all**: hoisted single check → one email + exit 1 before any merge attempt.

## Validation

Per AGENTS.md: no pipeline script runs until the pipeline is fully implemented; then validate with `_debug` (TODO.md:27 flow). Intermediate checks (no compute):

1. `bash -n` on both modified scripts.
2. Login node, no srun: `./3_submit_merge.sh Kfoury` → `Already merged: Kfoury — skipping ...` + exit 0 (unchanged single-dataset path).
3. Login node: `./3_submit_merge.sh` (bare) → default-all loop: Kfoury skipped, no-feather datasets (Zhu, fresh `_debug`) WARNING-skipped, remaining datasets merge; verify the summary email/exit code shape on the login-node stdout (emails best-effort, may be skipped without a mail CLI).
4. Login node: `./3_submit_merge.sh Kfoury --force` → coverage gate fail-closed with the (improved) message; exit 1 (expected).
5. `_debug` full flow once the pipeline is validated: prepare → array → `./3_submit_merge.sh` (bare, or `_debug`) → verify annotation columns in the merged h5ad (TODO.md:27); then re-run bare merge → all datasets skip.
