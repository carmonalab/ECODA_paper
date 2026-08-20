# Fix rds→h5ad conversion crash on Seurat v5 (Assay5) objects

## Symptom

Job `4288836` (`1.1_submit_combinedpbmc.sh`) FAILED. Both the Stephenson and
Zhu workers crashed in the embedded-R rds→h5ad conversion:

```
Error in (function (input_path, output_path)  :
  no slot of name "counts" for this object of class "Assay5"
```

(same error twice = Stephenson + Zhu, the two `.rds` sources of CombinedPBMC).

## Root cause

`src/utils/preprocess_utils.py:17-28` embeds R code in `convert_rds_to_raw_h5ad`
that uses Seurat **v4-only slot syntax**:

```r
seurat@assays$RNA@data <- seurat@assays$RNA@counts   # line 23 — crashes
...
seurat@assays$RNA@data <- NULL                       # line 25 — also v4-only
```

- In Seurat 5.2.0, `CreateSeuratObject()` (used inside
  `create_clean_seuratv5_object()`) produces **Assay5** objects, which have no
  `counts`/`data` slots (data lives in `layers`). `@counts` on Assay5 raises
  exactly the reported error.
- Reproduced locally (Seurat 5.2.0 + anndataR 1.3.0; SeuratObject 5.2.0 is
  pinned in pixi.lock for both linux-64/HPC and osx-arm64/local).
- Empirical anndataR behavior: `write_h5ad()` on an Assay5 object writes
  **X=None + layers["counts"]** (even if a `data` layer is populated), so the
  old hack would not have produced X=counts for v5 anyway — X must be promoted
  in Python (mirror of `base_preprocessing()` at
  `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:108-110`).

### Why other jobs "complete"

Joanito/Kfoury jobs completed only because their `{stem}_raw.h5ad` caches
exist in HPC scratch from earlier (Seurat-4-era) successful conversions; the
`if (!file.exists(output_path))` guard skips conversion. Stephenson/Zhu never
had caches → crash on every run. The same latent bug affects **all 12
datasets with `.rds` inputs** (Adams, Bassez, Joanito, Kfoury, Kim, Lee,
Pelka, Smillie, Stephenson, Wu, Zhang, Zhu) in the preprocess array
(`load_input` at `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:231`),
triggered whenever a cache is missing or `--force` is used.

### Secondary bug

`create_clean_seuratv5_object()` (`src/utils/seurat_utils.R:6-19`):
- else-if branch passes `counts = seurat@assays$RNA$counts` (always NULL there)
  → `CreateSeuratObject(counts = NULL)` → broken object.
- `@layers` slot access in that branch errors on classic v4 `Assay` objects.

## Fix

### 1. `src/utils/preprocess_utils.py`

- Embedded R `convert_rds_to_raw_h5ad`: delete the two `@data`/`@counts` slot
  lines (23-25); write the cleaned object directly:

```r
convert_rds_to_raw_h5ad <- function(input_path, output_path) {
  seurat <- readRDS(input_path)
  seurat <- create_clean_seuratv5_object(seurat)
  if (!file.exists(output_path)) {
    write_h5ad(seurat, output_path)
  }
}
```

- `load_single_input()` rds branch: after `sc.read_h5ad(raw_h5ad_path)`,
  promote the counts layer to X (exact mirror of `base_preprocessing`):

```python
adata = sc.read_h5ad(raw_h5ad_path)
if adata.X is None and "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()
return adata
```

This keeps the invariant "X is the raw counts" for all downstream code
(combinedpbmc workers, subsetting, intermediates, concat).

### 2. `src/utils/seurat_utils.R` — `create_clean_seuratv5_object()`

Robust counts extraction, v4- and v5-safe:

```r
create_clean_seuratv5_object <- function(seurat) {
  rna <- seurat@assays$RNA
  counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")
  if (is.null(counts) && inherits(rna, "Assay5")) {
    counts <- rna@layers[["X"]]
  }
  if (is.null(counts)) {
    counts <- GetAssayData(seurat, assay = "RNA", layer = "data")
  }
  if (!is.null(counts)) {
    seurat <- CreateSeuratObject(counts = counts, meta.data = seurat@meta.data)
  }
  return(seurat)
}
```

Verified semantics: `GetAssayData(layer="counts"/"data")` works for both
classic `Assay` (v4) and `Assay5` in Seurat 5.2.0; layer named `"X"` (produced
by SeuratIO `read_h5ad`) is only reachable on Assay5.

## Scope / non-changes

- No changes to `datasets.json`, NAS, staging, or the HPC submit scripts.
- `1.2.1_prepare_joanito.R:133` uses `$counts` (layer access, v5-safe) — no
  change needed.
- `benchmark_pipeline.R:705` caller of `create_clean_seuratv5_object()`
  benefits from the fix (h5ad-derived objects have layer `"X"`).

## Validation

### Done (local, inline repro — Seurat 5.2.0 / anndataR 1.3.0)

1. Mock "Stephenson-like" Assay5 RDS (counts layer + Sample/Status/Site meta)
   → fixed `create_clean_seuratv5_object` + `convert_rds_to_raw_h5ad` →
   `write_h5ad` → h5ad has X=None + layers["counts"] → Python promotion gives
   X = counts (shape/dtype/meta preserved). ✓
2. v4 classic `Assay` path: `GetAssayData`/`$counts` semantics verified — no
   regression. ✓

### To run after implementation

- Local smoke test: re-run the inline experiment with the actual edited files
  (source `src/utils/load_all_functions.R`, call the real
  `convert_rds_to_raw_h5ad` via rpy2 or Rscript on a mock Assay5 RDS).
- HPC (user): `git pull` in `~/ECODA_paper`, rerun
  `src/2_dataset_specific_preprocessing/1.1.1_run_worker.sh`-equivalent
  submission (the combinedpbmc sbatch worker) — expect COMPLETED; verify
  `${HPC_SCRATCH_DIR}/CombinedPBMC/data/` now contains the two `*_raw.h5ad`
  caches + `combined_pbmc_batch_effect_analysis.h5ad`.
- Spot-check an uncached rds dataset through the preprocess array
  (`--force` on one view, debug partition) to confirm the fix covers the
  3_scrnaseq_preprocessing path too.
