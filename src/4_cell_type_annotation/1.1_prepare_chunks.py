#!/usr/bin/env python3
#
# 1.1_prepare_chunks.py — Build per-dataset sample chunk files (production or test mode)
#
# Native Python replacement for the former reticulate-based 1.1_prepare_chunks.r.
#
# One chunk set per DATASET (not per view): all preprocessed view .h5ad files in
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output are merged into a single union h5ad
# (deduplicated on (sample, barcode)), written OUTSIDE the output dir to
# ${HPC_SCRATCH_DIR}/${DS_NAME}/annotation_union/union.h5ad so the NAS sync
# glob over output/*.h5ad stays clean. Samples are grouped into consecutive
# chunks of 2 (or 1 in --test mode) and chunk_<N>.txt files are written under
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks. Line 1 of each chunk file is the
# absolute union .h5ad path, subsequent lines are sample IDs.
#
# The union is ALWAYS written as a minimal h5ad (even for single-view datasets,
# whose preprocessed view keeps log-norm X + layers["counts"] and must not be
# chunked directly): X = raw counts (float64 CSR), obs + var only, no
# layers/obsp/obsm/uns/varm/varp. anndata's read_h5ad(backed="r") eagerly
# materializes the entire `layers` group into RAM at open (verified in anndata
# 0.12.10 and the HPC-pinned 0.12.19), so a union that kept counts in layers
# would add a full-matrix memory floor per annotation worker (Kfoury:
# +1.31 GB; a GongSharma-class union: 26-40 GB — OOM at any chunk size). With
# counts in X the backed open costs ~10 MB; X stays lazy and per-sample backed
# subsetting stays selective (CSR-on-disk). The annotation worker
# (2.1.1_process_chunk.R) treats the counts-absent -> X fallback as the
# designed primary path for union files (message, not warning).
#
# Union construction is memory-lean for the common case:
#   - exactly one view  -> streaming h5py copy (obs/var groups + counts as X)
#   - multiple views, one of which already contains the full union (e.g.
#     Stephenson: benchmark view is a cell subset of the batch-effect view)
#     -> same streaming copy of that view (backed obs-only reads)
#   - otherwise (partial overlaps between views) -> full in-memory concat +
#     dedup; this path needs a large srun allocation and prints a warning.

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
import scanpy as sc


def cell_keys(adata, sample_col):
    """(sample, barcode) composite keys for the obs of `adata`."""
    samples = adata.obs[sample_col].astype(str).values
    return np.char.add(np.char.add(samples, "_"), adata.obs_names.astype(str).values)


def verify_annotation_union(union_path, expected_nnz, sample_col):
    """
    Backed-mode self-check of a written annotation union: X must be a lazy
    backed CSR dataset (not eagerly materialized), its nnz must equal the
    source counts nnz, the file must carry no layers/obsp/obsm data, and the
    sample column must be present. Raises on any mismatch.
    """
    adata = ad.read_h5ad(union_path, backed="r")
    try:
        if type(adata.X).__name__ != "_CSRDataset":
            raise RuntimeError(
                f"Union X is '{type(adata.X).__name__}' (expected lazy backed "
                "CSR) — on-disk layout is not the minimal union."
            )
        if len(adata.layers) or len(adata.obsp) or len(adata.obsm):
            raise RuntimeError("Union unexpectedly contains layers/obsp/obsm data.")
        if sample_col not in adata.obs.columns:
            raise KeyError(f"Sample column '{sample_col}' not found in union obs.")
        with h5py.File(union_path, "r") as f:
            nnz = int(f["X/indptr"][-1])
        if nnz != expected_nnz:
            raise RuntimeError(
                f"Union X nnz ({nnz}) != source counts nnz ({expected_nnz})."
            )
    finally:
        adata.file.close()
    print(f"  Union self-check OK: lazy X ({nnz} nnz), no layers/obsp/obsm, "
          f"'{sample_col}' present.")


def write_annotation_union(source_path, union_path, sample_col):
    """
    Streaming, memory-lean annotation-union writer: copies the obs/var groups
    (datasets + attrs) from `source_path` and re-writes the raw counts as X
    (float64 CSR, same dtype as today's layers["counts"]), so the worker's
    astype("float64") path is a no-op. No layers/obsp/obsm/uns/varm/varp are
    copied. Writes atomically (tmp + os.replace) and self-verifies the result.
    """
    union_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = union_path.with_name(union_path.name + ".tmp")
    try:
        with h5py.File(source_path, "r") as src, h5py.File(tmp_path, "w") as dst:
            dst.attrs["encoding-type"] = src.attrs.get("encoding-type", "anndata")
            dst.attrs["encoding-version"] = src.attrs.get("encoding-version", "0.1.0")
            for grp in ("obs", "var"):
                src.copy(src[grp], dst, name=grp)
            if "layers/counts" in src:
                counts = src["layers/counts"]
            else:
                print(f"  No layers/counts in {source_path.name}; using its X as "
                      "the counts source (preprocessed views should always "
                      "carry a counts layer).")
                counts = src["X"]
            # Fail closed on non-CSR sources: preprocessed views are CSR by
            # construction (1.1.1_preprocess.py), so anything else means the
            # file was replaced/regenerated elsewhere. Copying a CSC layout
            # under the csr_matrix encoding-type below would silently corrupt
            # the union (the nnz/shape self-check cannot catch it).
            counts_enc = str(counts.attrs.get("encoding-type"))
            if counts_enc != "csr_matrix":
                raise RuntimeError(
                    f"{source_path.name}: counts source encoding-type is "
                    f"'{counts_enc}', expected 'csr_matrix'."
                )
            x = dst.create_group("X")
            x.attrs["encoding-type"] = "csr_matrix"
            x.attrs["encoding-version"] = "0.1.0"
            x.attrs["shape"] = counts.attrs["shape"]
            for ds_name in ("data", "indices", "indptr"):
                src.copy(counts[ds_name], x, name=ds_name)
            expected_nnz = int(counts["indptr"][-1])
        os.replace(tmp_path, union_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    verify_annotation_union(union_path, expected_nnz, sample_col)


def build_union(h5ad_files, union_path, sample_col):
    """
    Build the per-dataset union h5ad (dedup on (sample, barcode)).
    Always writes the minimal annotation union and returns `union_path`.
    """
    if len(h5ad_files) == 1:
        print(f"Single view: writing annotation union from {h5ad_files[0].name} "
              "(streaming copy; the view itself keeps its full layout).")
        write_annotation_union(h5ad_files[0], union_path, sample_col)
        return union_path

    print(f"Building union from {len(h5ad_files)} views: "
          + ", ".join(f.name for f in h5ad_files))

    # 1. Compare (sample, barcode) key sets using backed obs-only reads (no
    #    matrix I/O). If one view already equals the union, stream it.
    key_sets = []
    for f in h5ad_files:
        adata = ad.read_h5ad(f, backed="r")
        keys = set(cell_keys(adata, sample_col))
        key_sets.append(keys)
        adata.file.close()
        print(f"  {f.name}: {len(keys)} cells")

    union_keys = set().union(*key_sets)
    largest_idx = int(np.argmax([len(k) for k in key_sets]))
    print(f"Union: {len(union_keys)} cells")

    if len(key_sets[largest_idx]) == len(union_keys):
        print(f"  Union == view {h5ad_files[largest_idx].name}; streaming copy "
              "(no in-memory concat).")
        write_annotation_union(h5ad_files[largest_idx], union_path, sample_col)
        return union_path

    # 2. Partial-overlap views: concat in memory and dedup. Memory-heavy; the
    #    srun allocation for this dataset must be raised accordingly.
    print("WARNING: views partially overlap; concatenating in memory. "
          "If this OOMs, raise the srun --mem for this dataset in "
          "1_prepare_chunks.sh.")
    adatas = [ad.read_h5ad(f) for f in h5ad_files]
    adata = sc.concat(adatas, join="outer", index_unique=None)
    del adatas
    keys = cell_keys(adata, sample_col)
    # keep first occurrence of each (sample, barcode) key; fancy indexing
    # already allocates fresh matrices, so no extra .copy() is needed
    _, first_idx = np.unique(keys, return_index=True)
    adata = adata[np.sort(first_idx)]

    # Minimal union layout: X = raw counts (CSR), obs/var only — mirrors the
    # streaming writer. anndata writes empty group stubs for emptied mappings,
    # so strip them on disk afterwards.
    counts = adata.layers.get("counts")
    if counts is None:
        counts = adata.X
        if counts is None:
            raise RuntimeError("No counts layer or X to use as the union's X.")
    if not sp.isspmatrix_csr(counts):
        counts = sp.csr_matrix(counts)
    adata.X = counts
    adata.layers = {}
    adata.obsm = {}
    adata.obsp = {}
    adata.uns = {}
    adata.varm = {}
    adata.varp = {}
    expected_nnz = int(counts.nnz)

    union_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(union_path))
    del adata
    with h5py.File(union_path, "r+") as f:
        for grp in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
            if grp in f:
                del f[grp]
    verify_annotation_union(union_path, expected_nnz, sample_col)
    return union_path


def main():
    parser = argparse.ArgumentParser(description="Build per-dataset sample chunk files.")
    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help="Test mode: 1 sample per chunk (production default: 2).",
    )
    args = parser.parse_args()

    project_root = os.environ.get("PROJECT_ROOT")
    if not project_root:
        sys.exit("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

    ds_name = os.environ.get("DS_NAME")
    if not ds_name:
        sys.exit("CRITICAL Error: DS_NAME not set. Ensure it is exported before calling this script.")

    sample_col = os.environ.get("SAMPLE_COLNAME")
    if not sample_col:
        sys.exit("CRITICAL Error: SAMPLE_COLNAME not set. Ensure it is exported before calling this script.")

    hpc_scratch_dir = os.environ.get("HPC_SCRATCH_DIR")
    if not hpc_scratch_dir:
        sys.exit("CRITICAL Error: HPC_SCRATCH_DIR not set. Source slurm_config.sh before calling this script.")

    path_data = Path(hpc_scratch_dir) / ds_name / "output"
    path_union_dir = Path(hpc_scratch_dir) / ds_name / "annotation_union"
    path_output_chunks = path_data / "chunks"

    print(f"Path is: {path_data}")

    # Skip rds->h5ad conversion caches written into the output dir by
    # src.utils.py.preprocess_utils.load_single_input() (they lack the standardized sample
    # column and are not preprocessed views). The union h5ad lives in
    # annotation_union/, outside this glob.
    h5ad_files = sorted(
        f for f in path_data.glob("*.h5ad") if not f.name.endswith("_raw.h5ad")
    )
    print(f"Files found: {', '.join(f.name for f in h5ad_files)}")
    if not h5ad_files:
        sys.exit(f"CRITICAL Error: No preprocessed .h5ad files found in {path_data} "
                 "(run the preprocess array first).")

    # Defensive fail-closed completeness check (mirrors the bash guard in
    # 1_prepare_chunks.sh, mirroring 1.1.1_preprocess.py skip semantics):
    # every view output the preprocess array is expected to produce must
    # already exist. Catches bypasses/drift of the bash check — the dataset
    # lands in FAILED_DATASETS instead of building a partial union.
    datasets_json = os.environ.get("DATASETS_JSON_FILE")
    if not datasets_json:
        print("WARNING: DATASETS_JSON_FILE not set; skipping the expected-view "
              "check (source slurm_config.sh on HPC).")
    else:
        with open(datasets_json) as f:
            ds_entry = json.load(f).get(ds_name, {})
        expected = {
            v.get("output_file_name")
            for v in ds_entry.get("views", {}).values()
            if v.get("input_file_name") is not None and v.get("output_file_name")
        }
        missing = expected - {f.name for f in h5ad_files}
        if missing:
            sys.exit(
                f"CRITICAL Error: preprocessing incomplete for {ds_name}: missing "
                f"expected view file(s) {sorted(missing)} in {path_data} "
                "(run the preprocess array first)."
            )

    # Delete chunk file folder recursively to ensure a perfectly clean start
    shutil.rmtree(path_output_chunks, ignore_errors=True)
    path_output_chunks.mkdir(parents=True, exist_ok=True)

    # Delete any stale union from a previous run (e.g. after --force
    # preprocess re-runs regenerated the views).
    shutil.rmtree(path_union_dir, ignore_errors=True)
    path_union_dir.mkdir(parents=True, exist_ok=True)

    # Number of samples per chunk
    chunk_size = 1 if args.test else 2

    # Build the per-dataset minimal annotation union
    union_file = build_union(h5ad_files, path_union_dir / "union.h5ad", sample_col)

    # Read the union file in backed mode and extract its unique samples
    adata = ad.read_h5ad(union_file, backed="r")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column '{sample_col}' not found in obs of {union_file}")
    file_samples = adata.obs[sample_col].astype(str).unique()
    adata.file.close()
    print(f"Union file {union_file.name} has {len(file_samples)} samples.")

    # Split the union's samples into consecutive groups of chunk_size
    sample_groups = [
        file_samples[i : i + chunk_size] for i in range(0, len(file_samples), chunk_size)
    ]

    # Write each group to a unique chunk file (global counter within the dataset)
    global_chunk_counter = 1
    for group in sample_groups:
        # We write the source file as the VERY FIRST line, followed by the sample IDs
        chunk_path = path_output_chunks / f"chunk_{global_chunk_counter}.txt"
        chunk_path.write_text(
            "\n".join([str(union_file.resolve())] + [str(s) for s in group]) + "\n"
        )

        global_chunk_counter += 1

    print(f"Successfully generated {global_chunk_counter - 1} total chunk files.")

    # Delete stale per-chunk annotation feathers from earlier runs (chunk
    # numbering is per-dataset and changes with chunk size / sample set, so a
    # rerun must not merge leftovers from a previous numbering). Only after
    # all chunks were generated successfully, and never in --test mode
    # (test runs must not destroy production annotation results).
    if not args.test:
        for stale in path_data.glob("annotations_chunk_*.feather"):
            stale.unlink()
            print(f"Removed stale annotations file: {stale.name}")


if __name__ == "__main__":
    main()
