# Rename merge scripts per naming convention + fix reviewed bugs

## Context

`src/4_cell_type_annotation/` pairs are named by call depth (see
`docs/ARCHITECTURE.md:101` and the earlier renaming plan):
top-level entry point = `N_action.sh`, directly-called script = `N.N_action.ext`.

- `1_prepare_chunks.sh` (top-level) → `1.1_prepare_chunks.py` (called) ✓
- `2_submit_hpc_array.sh` (top-level) → `2.1_run_worker.sh` (called) → `2.1.1_process_chunk.R` ✓
- `3.1_submit_merge.sh` (top-level) → `3_merge_annotations.py` (called) ✗ **inverted**

The merge pair must be swapped to match. Scope per user: rename + fix ALL
findings from the review.

## Task 1 — Rename (git mv, preserves history)

1. `git mv src/4_cell_type_annotation/3.1_submit_merge.sh src/4_cell_type_annotation/3_submit_merge.sh`
2. `git mv src/4_cell_type_annotation/3_merge_annotations.py src/4_cell_type_annotation/3.1_merge_annotations.py`

## Task 2 — Update references

- `src/4_cell_type_annotation/3_submit_merge.sh` (renamed): header comment (line 5), usage (lines 11, 30), prose comment (line 16), and the `"${SCRIPT_DIR}/3_merge_annotations.py"` call (line 102) → `3.1_merge_annotations.py`; self-name refs → `3_submit_merge.sh`.
- `src/4_cell_type_annotation/3.1_merge_annotations.py` (renamed): no self-name refs inside; nothing to change.
- `README.md:107`, `AGENTS.md:46`, `TODO.md:14`, `TODO.md:36`: `3.1_submit_merge.sh` → `3_submit_merge.sh`.
- `docs/ARCHITECTURE.md`: lines 57, 125, 143, 144, 155, 170, 369 — replace name tokens both ways (`3.1_submit_merge.sh` → `3_submit_merge.sh`, `3_merge_annotations.py` → `3.1_merge_annotations.py`); keep the descriptive prose unchanged.
- `.kilo/plans/*.md`: leave untouched (historical records, established practice from the earlier renaming plan).

## Task 3 — Bug fixes in `3.1_merge_annotations.py` (renamed)

Review findings (verified against `2.1.1_process_chunk.R` feather writer and
`1.1_prepare_chunks.py` union builder):

1. **Silent exit-0 on failure paths (HIGH)** — `merge_annotations()` returns
   normally (exit 0) on "No annotation feather files found" (line 33-35) and
   "'Sample' column missing in annotation feathers" (line 46-48). The shell
   wrapper then proceeds and rsyncs an unannotated h5ad over good NAS data.
   Fix: replace both `return` paths with `sys.exit(1)` (add `import sys`).
2. **No match-rate validation (HIGH)** — a key-format drift (e.g. obs_names /
   sample-dtype mismatch between view h5ad and union feather) yields 0 joined
   rows with no error. Fix: after computing `merged_count` (line 88), if
   `merged_count == 0` → print clear error + `sys.exit(1)`; else log the
   match rate (`merged_count / original_n_obs`).
3. **`_key` column clobbering (LOW)** — if `obs` already contains a `_key`
   column it is silently overwritten and lost from the output. Fix: use a
   collision-resistant internal name (e.g. `__annot_key`) for both sides and
   the post-join drop.

## Task 4 — Bug fixes in `3_submit_merge.sh` (renamed)

1. **Stale intermediate files synced to NAS (MED)** —
   `2_submit_hpc_array.sh` rsyncs `output/` to NAS *before* the merge, so
   `annotations_chunk_*.feather` and `chunks/` land on NAS; the merge's own
   rsync (line 121) has no `--delete`, so they persist and accumulate on NAS.
   Fix: after the successful rsync, remove them from the NAS output dir:
   `rm -f "${NAS_TARGET_DIR}/${DS_NAME}/output"/annotations_chunk_*.feather`
   and `rm -rf "${NAS_TARGET_DIR}/${DS_NAME}/output/chunks"`. Keep the LOCAL
   feathers — the coverage gate (lines 69-79) counts them on re-runs.
2. **`.tmp` leftovers rsynced to NAS (LOW)** — a crash between the atomic
   write and `os.replace` in the python leaves `*.h5ad.tmp` in `OUTPUT_DIR`.
   Fix: add `--exclude='*.tmp'` to the rsync on line 121.
3. **Fixed `--mem=64G` (LOW)** — may OOM on the largest views (Stephenson
   batch-effect view, already flagged in ARCHITECTURE.md). Fix: make it
   overridable, e.g. `MERGE_MEM="${MERGE_MEM:-64G}"` and use `--mem="${MERGE_MEM}"`.

## Validation (no pipeline runs — AGENTS.md rule)

1. `bash -n src/4_cell_type_annotation/3_submit_merge.sh`
2. `python -m py_compile src/4_cell_type_annotation/3.1_merge_annotations.py`
3. `grep -rn '3\.1_submit_merge\|3_merge_annotations' src/ docs/ README.md AGENTS.md TODO.md` → only `.kilo/plans/` hits remain.
4. `git status` shows both renames as detected (similarity rename).

## Out of scope

- `2_submit_hpc_array.sh` pre-merge NAS sync (feathers landing on NAS
  upstream is mitigated by Task 4.1); the `_raw.h5ad` files excluded from the
  merge loop are deliberate; the srun-per-view merge strategy stays as-is.
