# Speed up `1.1.1_create_combinedpbmc_dataset.py` (backed read + in-job parallelism)

## Why it is slow (diagnosis, verified against code)

1. **GongSharma 40 GB fully loaded before downsampling** — `sc.read_h5ad` at
   `1.1.1_create_combinedpbmc_dataset.py:158` reads both files (28 GB + 12 GB) into RAM and
   single-threadedly decompresses them; the 15-sample pick happens *after* (lines 161-165).
2. **Stephenson (3.6 GB .rds) + Zhu (1 GB .rds) via R interop** — `convert_rds_to_raw_h5ad`
   (`src/utils/preprocess_utils.py:18`): `readRDS` + `CreateSeuratObject` re-copy +
   `data <- counts` duplicate + anndataR `write_h5ad`, all single-threaded. Cached to
   `{stem}_raw.h5ad` after the first run (reruns skip this).
3. **Strictly serial**: Stephenson → GongSharma → Zhu → concat → write. Wall time is the sum.
4. **CPUs are irrelevant**: the sbatch already requests 16 (`1.1_submit_combinedpbmc.sh:6`), but
   h5ad reads, R conversion and anndata writes are effectively single-threaded. The previous
   answer's "raise cpus-per-task from 4 to 16" advice was stale.

## Goal

Same output semantics (concat order Stephenson→GongSharma→Zhu, common-gene intersection with
`<5000` union fallback, `Sample`/`batch`/`cond` obs trimming, `rng(123)` 15-sample pick, output
`combined_pbmc_batch_effect_analysis.h5ad`), but wall time ≈ max(source) instead of sum, and
GongSharma read cost cut ~3-5x via backed-mode slicing.

## Design

Single sbatch job (unchanged: 16 cpus, 128G, `module load GCCcore/12.2.0`). The Python script:

1. Runs the three sources in **3 worker processes** (`concurrent.futures.ProcessPoolExecutor`,
   `max_workers=3`, fork context on Linux HPC). Each worker writes its own small intermediate
   h5ad into `output_dir / "_intermediates"` (`<source>_subset.h5ad`) and returns `(n_obs, path)`.
2. Main process (after workers complete) loads the 3 small intermediates, computes common genes,
   prefixes Sample IDs, `sc.concat(join="outer")`, writes the combined file — exactly as today.

### Per-source workers

- **Stephenson worker** (needs R): lazy `from src.utils.preprocess_utils import load_input,
  apply_subset_vars` *inside* the worker function so rpy2/R init happens only in worker processes
  (parent never imports rpy2 → clean fork). Same pipeline as today: `load_input` (rds→h5ad
  cache) → `apply_subset_vars("benchmark_analysis")` → `standardize_gene_symbols` +
  `var_names_make_unique` → keep cols → set `batch`/`cond`/`Sample` → write intermediate.
- **Zhu worker** (needs R): same, no view/subsetting.
- **GongSharma worker** (no R needed):
  1. `sc.read_h5ad(path, backed='r')` for each of the two files.
  2. Union of unique `sample_col` values across both files; pick `min(15, n)` with
     `np.random.default_rng(123)` — identical seed/logic to today's lines 161-163.
  3. On each backed object: drop all layers except `counts` (backed deletion does not load X;
     avoids materializing heavy scaled layers), then `subset = adata[mask].to_memory()` (loads
     only the HDF5 chunks containing the chosen samples).
  4. `standardize_gene_symbols` + `var_names_make_unique` → keep cols → set
     `batch`/`cond`/`Sample` → write intermediate.
  5. **Fallback**: wrap the backed path in try/except; on failure fall back to today's full
     `sc.read_h5ad` in-memory load with a warning (correctness over speed, one-shot script).

### Robustness / details

- Keep the existing `--config-path` / `--base-path` / `--output-dir` / `--layout` CLI and the
  per-dataset vs flat path resolution exactly as-is.
- Move the top-level `from src.utils.preprocess_utils import ...` import (line 47) into the
  worker functions to keep the parent rpy2-free.
- Preserve `standardize_gene_symbols` + `var_names_make_unique` ordering, obs column trimming,
  `keep_base = ["Sample", "batch", "cond"]`, and the `<5000` common-gene → union fallback in main.
- Intermediates live in `output_dir/_intermediates/` (underscore dir, no glob collisions with
  the `_raw.h5ad` cache or any consumer). Overwritten unconditionally on rerun.
- Print per-worker progress (`Loading Stephenson...`, GongSharma sample pick counts, etc.) as
  today so the log stays recognizable.

## Tasks

1. **Refactor `src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py`**
   - Extract per-source load+prepare logic into three worker functions (Stephenson, GongSharma,
     Zhu) that write intermediates; lazy rpy2 imports inside the two R workers.
   - Implement the backed-mode GongSharma path with the counts-layer-only + to_memory slicing
     and the full-load fallback.
   - Replace the serial `main()` body with `ProcessPoolExecutor(3)` + intermediate loading +
     shared final steps (common genes, prefix, concat, write).
   - Update the module docstring (parallel workers, backed read, `_intermediates/`).
2. **`src/2_dataset_specific_preprocessing/1.1_submit_combinedpbmc.sh`** — no change required
   (16 cpus / 128G already set; verify comment mentioning memory if any becomes stale).
3. **`docs/ARCHITECTURE.md`** — update the "CombinedPBMC combine" bullet (line ~99): note backed
   GongSharma read, 3 parallel in-job workers, `_intermediates/` dir.
4. **`TODO.md`** — update the cluster-verify line (~37) to note the CombinedPBMC step runtime
   after optimization (was ~40-60 min serial; expect ~20-25 min first run, faster on reruns due
   to the `_raw.h5ad` cache).

## Validation

- Local, safe: `python -m py_compile src/2_dataset_specific_preprocessing/1.1.1_create_combinedpbmc_dataset.py`.
- HPC (user runs; do not run scripts from this agent): resubmit
  `1.1_submit_combinedpbmc.sh` and check:
  - Log lines match the previous run's expectations: 3 batches; GongSharma picked exactly 15
    samples; `Combined shape` printed; output written to
    `${HPC_SCRATCH_DIR}/CombinedPBMC/data/combined_pbmc_batch_effect_analysis.h5ad`.
  - If the previous run's combined h5ad still exists on scratch: compare `n_obs`, `n_vars`,
    `obs` columns (Sample/batch/cond) and sample counts per batch.
  - `squeue`/`sacct` wall time vs. previous run.

## Risks / notes

- Backed-mode slicing reads HDF5 chunks containing the selected cells: with 15 of ~50+ samples
  this reads ~30-50% of the GongSharma data instead of 100% — the win scales with sample count;
  the full-load fallback guarantees correctness regardless.
- Fork + rpy2 in workers is the standard pattern on Linux (only the two R workers init R; the
  parent stays rpy2-free, so no fork-after-R-init concerns).
- If the current HPC job (4284619) completes first, its `_raw.h5ad` caches for Stephenson/Zhu
  make the optimized run's R conversion nearly free.
- Peak memory ≈ GongSharma subset (~5-15 GB) + Stephenson raw (~5-10 GB) + Zhu (~1 GB) +
  join (~10-20 GB) across processes — well under 128G; no sbatch mem change.
