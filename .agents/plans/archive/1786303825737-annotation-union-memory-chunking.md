# Annotation pipeline: memory-lean union + chunk size 2 + worker tuning

## Goal

Three measured improvements to `src/4_cell_type_annotation/`:

1. **Kill the per-dataset memory floor**: anndata's `read_h5ad(backed="r")` eagerly loads the entire `layers` group (verified in anndata 0.12.10 and the HPC-pinned 0.12.19), so every worker loads the full union `counts` layer into RAM at open (Kfoury: +1.31 GB; a GongSharma-class union: 26–40 GB — OOM at any chunk size). Fix: write the annotation union with **counts as `X`** and no `layers`/`obsp`/`obsm`/`uns`/`varm`/`varp` (measured open cost: +10 MB; `X` stays lazy, per-sample subsetting stays selective).
2. **Production chunk size 5 → 2** (2 samples per chunk; `--test` stays 1).
3. **Worker sbatch**: `--cpus-per-task=2`, `--mem=32G` (keep — safe side), `--time=01:00:00` (measured: ~25 s/sample at 4k cells, ~40 s at 23k cells, ~15–30 s startup → tasks stay ~1–2 min; 5-retry worst case ≪ 1 h).

Measured anchors (from profiling session): scATOMIC peak RSS 6.0 GB @ 2.5k cells → 17.6 GB @ 23.2k cells (~20–23 GB for 20k-cell samples with imputation on → 32G keeps margin; 16G would not). Startup: R packages 7.7 s / 1.06 GB, ref maps ~245 MB total (~5–15 s), scGate DB < 1 s.

## Design decisions

- **Union layout**: `X` = counts (float64 CSR — same dtype as today's `layers["counts"]`, so the worker's `astype("float64")` path is a no-op), keep `obs` + `var` only, everything else dropped. Always write the union (even for single-view datasets — the current "chunk the view directly" fast path must go, because the view carries `layers` and log-norm `X`).
- **Union write method**: streaming h5py copy (copy `obs`/`var` groups with attrs; copy `layers/counts/{data,indices,indptr}` into a new `X` group; set `encoding-type=csr_matrix` / `encoding-version=0.1.0` on `X` and root attrs). Memory-lean → the existing 4G prepare `srun` stays. Fallback if anndata rejects the file: in-memory anndata rewrite for small files only.
- **Worker**: the `counts`-absent → `X` fallback becomes the designed primary path for unions; downgrade the `warning(...)` to a `message(...)`.
- **Chunking/merge/manifest/feathers**: unchanged (all chunk-size-agnostic; coverage gate in `3_submit_merge.sh` untouched).
- **Preprocessed views unchanged**: views keep `X` (log-norm) + `layers["counts"]`; only the annotation union changes.

## Tasks

### 1. `src/4_cell_type_annotation/1.1_prepare_chunks.py`
- Add `write_annotation_union(source_path, union_path)`:
  - h5py streaming copy of `obs` + `var` groups (incl. their attrs and `_index`),
  - copy `layers/counts/{data,indices,indptr}` datasets into `X` (preserve float64/int32 dtypes),
  - set `X` and root `encoding-type`/`encoding-version` attrs,
  - self-verify: `ad.read_h5ad(union_path, backed="r")` → `X` lazy (`_CSRDataset`), `X.nnz == source counts nnz`, no `layers`/`obsp`/`obsm`, sample column present.
- Rework `build_union()`:
  - 1 view → `write_annotation_union(view, union_path)` (drop the "chunk view directly" branch; the key-set comparison is unnecessary for a single view).
  - ≥ 2 views → keep the existing backed obs-only key-set comparison: largest view equals union → `write_annotation_union(largest_view, ...)`; partial overlap → existing in-memory `sc.concat` + dedup, then write with `X = counts` (CSR), obs/var only, no layers/obsp/obsm.
  - Always return `union_path`.
- `chunk_size = 1 if args.test else 2`.
- Update the module docstring: union layout rationale (anndata eager-`layers` load), chunks of 2.

### 2. `src/4_cell_type_annotation/2.1.1_process_chunk.R`
- When `"counts"` is not in `layer_keys` (i.e., annotation-union file), replace the `warning(...)` with `message(...)` stating that the union carries counts in `X` by design (fallback logic unchanged).
- Keep everything else (`ncores = 1`, timeouts, feather naming, CSR-format warning).

### 3. `src/4_cell_type_annotation/2.1_run_worker.sh`
- `#SBATCH --cpus-per-task=4` → `2`
- `#SBATCH --mem=32G` (unchanged)
- `#SBATCH --time=02:00:00` → `01:00:00`

### 4. `src/4_cell_type_annotation/1_prepare_chunks.sh`
- Usage header comment: `chunk_size = 5` → `2`.

### 5. Docs
- `docs/ARCHITECTURE.md`:
  - line ~123: `(5 samples each)` → `(2 samples each)`
  - line ~144 (`1.1_prepare_chunks.py`): union is now always written as a minimal h5ad (`X` = counts, no layers/obsp/obsm — anndata eager-`layers` rationale), chunks of 2
  - line ~147 (`2.1_run_worker.sh`): `(16G, 2h` → `(32G, 2 cpus, 1h` (also fixes the stale 16G text)
  - line ~148 (`2.1.1_process_chunk.R`): counts come from `X` for the annotation union (fallback is the designed path)
  - line ~63 (annotation_union dir) and line ~160 (backed-reads note): document the minimal-union layout
  - line ~167: `5 samples/chunk` → `2 samples/chunk`
- `AGENTS.md` line ~115 (worker invariant): note that the annotation union carries counts in `X` (the `X` fallback is by design, not a warning case).
- README: no chunk-size text — no change.

### 6. Validation
- **Local smoke (Mac, NAS-mounted)**: run the new union writer on the Kfoury view; verify backed-open RSS (+10 MB class, not +1.3 GB), `X.nnz` == source counts nnz, sample listing + chunk files correct; run the worker's R extraction path on the written union (X fallback) and confirm byte-identical matrix (nnz 6,056,210 for GSM4274683).
- **HPC (per repo convention, once the pipeline is fully implemented, using `_debug`)**: `./1_prepare_chunks.sh test _debug` → `./2_submit_hpc_array.sh _debug` → `./3_submit_merge.sh _debug`; verify feathers merge, 100% match rate, `sacct` MaxRSS per task ≪ 32G. Then re-run Kfoury production (chunk=2) and compare against the existing NAS feathers/annotated h5ad.
- Also verify a `chunk_*.txt` file contains exactly 2 samples in production mode.

### 7. Commit (AGENTS.md workflow)
- Archive this plan to `.kilo/plans/archive/`, `git add .`, commit, push.

## Risks / notes

- **Scratch disk**: the union copy now exists for every dataset (single-view included) — Kfoury +~1 GB; GongSharma-class ~26–40 GB. Acceptable; the prepare `srun` (4G, 30 min) is fine since the copy is streaming I/O; bump only if a giant dataset times out.
- **h5py attrs/encoding correctness**: guarded by the self-verify step in task 1 and the R extraction check in validation; in-memory fallback for small files if anndata rejects the layout.
- **Memory**: 32G keeps margin for 20k-cell samples (~20–23 GB with imputation on). If a future dataset has samples ≫ 40k cells, revisit.
- `MAX_NUM_CHUNKS_PARALLEL=500` unchanged.

## Out of scope

- HiTME memory profiling (not measured locally; Kfoury HPC timestamps show it is fast).
- Changing the preprocessed view layout (views keep log-norm `X` + `layers["counts"]`).
- Any change to `2_submit_hpc_array.sh` / `3_submit_merge.sh` (chunk-size-agnostic).
