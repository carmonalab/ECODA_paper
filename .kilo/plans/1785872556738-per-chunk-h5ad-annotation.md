# Guarantee CSR-on-disk for preprocessed h5ads (memory-light backed annotation)

## Answer to the design question

**Can obs be sliced with the adata in backed mode?** Yes. obs is stored as plain metadata
(separate HDF5 group), so `adata$obs` reads metadata only (few columns x cells; tens of MB even at
1M cells) and `adata[rows]` creates a lazy view whose `.obs`/`.layers` are sliced on access. obs
never triggers matrix I/O.

The memory-heavy operation is not obs but the sparse-matrix row subset `adata[cell_indices]`
inside `get_seurat_obj_from_h5ad()` (src/4_cell_type_annotation/2.1.1_process_chunk.R:147). Its
behavior depends only on the on-disk sparse format (anndata 0.12.10, `.pixi/`):

- **CSR on disk** -> `backed_csr_matrix._get_arrayXslice` -> `get_compressed_vectors`
  (`anndata/_core/sparse_dataset.py:208-216, 274-289`): reads ONLY the selected rows' segments.
  Empirically verified: subsetting 2000 cells of a 3M-nnz CSR file ~ 0 MB extra RSS.
- **CSC on disk** -> not overridden -> `self.to_memory()[rows]` (`sparse_dataset.py:382-389, 466-503,
  659-667`): full `data`/`indices` materialized per call. Verified: ~full matrix size extra RSS.

## Key finding: the CSR guarantee is incomplete

`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:116-119` only converts **dense** inputs to CSR:

```python
if not sp.issparse(adata.X):                       # sparse (e.g. CSC) inputs pass through unchanged
    adata.X = sp.csr_matrix(adata.X)
if not sp.issparse(adata.layers["counts"]):
    adata.layers["counts"] = sp.csr_matrix(adata.layers["counts"])
```

Sparse inputs (CSC is the h5ad default and typical of R-written raw files) are written as CSC ->
per-sample subsets in annotation full-load the dataset -> the reported RAM problem. With the fix
below, per-sample subsets are selective and the current backed design is correct as-is:
**no per-chunk h5ad files, no prepare_chunks / slurm changes, no changes to
`get_seurat_obj_from_h5ad()`.**

All preprocessed files are created fresh by the pipeline (no legacy preprocessed files exist), so
the fix fully determines the on-disk format; the raw input files are never rewritten.

## Tasks

1. **`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` — force CSR unconditionally** in
   `base_preprocessing()` (replace lines 116-119):
   ```python
   adata.X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
   adata.layers["counts"] = adata.layers["counts"].tocsr() if sp.issparse(adata.layers["counts"]) else sp.csr_matrix(adata.layers["counts"])
   ```
   (`tocsr()` on an already-CSR matrix returns itself — no copy.) scanpy ops after this
   (normalize_total, log1p, HVG selection, subsetting) preserve CSR, so the written file stays CSR.

2. **Optional hardening:** in `2.1.1_process_chunk.R`, after `read_h5ad(h5ad_file, backed="r")`,
   warn (not stop) when the on-disk X format is `csc` (`adata$X$format` via reticulate), so a
   non-CSR file fails loudly instead of silently OOMing.

3. **Docs (`AGENTS.md`, `docs/ARCHITECTURE.md`):** note that preprocessed h5ads are CSR-on-disk
   (X AND `layers/counts`) by construction; backed-mode per-sample subsets in annotation are then
   selective; obs is metadata-only (small). Update any statements from the previous per-chunk-h5ad
   draft (not implemented).

## Explicitly unchanged

- `src/4_cell_type_annotation/2.1.1_process_chunk.R` (keeps backed mode, per-sample
  `get_seurat_obj_from_h5ad()` calls)
- `src/4_cell_type_annotation/1.1_prepare_chunks.py`, `1_prepare_chunks.sh`, `2_submit_hpc_array.sh`
- `src/utils/seurat_utils.R` (`get_seurat_obj_from_h5ad()`, incl. its `astype("float64")` cast —
  per-sample cost only; leave for scATOMIC/HiTME compatibility)

## Risks

- Raw input sparse format (CSC vs dense) on the NAS datasets is unknown; the unconditional
  `tocsr()` handles both, so this is mitigated by construction.
- The change is memory/I-O behavior only; annotation results are unaffected.

## Out of scope

- Validation runs (need NAS/HPC access; to be done by the user in the HPC debug run).
- Per-chunk h5ad files (dropped — unnecessary once files are CSR), QOT/PULSAR benchmark methods,
  batch-effect analysis, `datasets.json`.
