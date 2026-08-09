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
# chunks of 5 (or 1 in --test mode) and chunk_<N>.txt files are written under
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks. Line 1 of each chunk file is the
# absolute union .h5ad path, subsequent lines are sample IDs.
#
# Union construction is memory-lean for the common case:
#   - exactly one view  -> chunk that file directly (no copy)
#   - multiple views, one of which already contains the full union (e.g.
#     Stephenson: benchmark view is a cell subset of the batch-effect view)
#     -> the union is a plain file copy of that view (obs-only reads)
#   - otherwise (partial overlaps between views) -> full in-memory concat +
#     dedup; this path needs a large srun allocation and prints a warning.

import argparse
import os
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
import scanpy as sc


def cell_keys(adata, sample_col):
    """(sample, barcode) composite keys for the obs of `adata`."""
    samples = adata.obs[sample_col].astype(str).values
    return np.char.add(np.char.add(samples, "_"), adata.obs_names.astype(str).values)


def build_union(h5ad_files, union_path, sample_col):
    """
    Build the per-dataset union h5ad (dedup on (sample, barcode)).
    Returns the path of the file to chunk (the union, or a single-view file).
    """
    if len(h5ad_files) == 1:
        print(f"Single view: chunking {h5ad_files[0].name} directly (no union file).")
        return h5ad_files[0]

    print(f"Building union from {len(h5ad_files)} views: "
          + ", ".join(f.name for f in h5ad_files))

    # 1. Compare (sample, barcode) key sets using backed obs-only reads (no
    #    matrix I/O). If one view already equals the union, copy it.
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
        print(f"  Union == view {h5ad_files[largest_idx].name}; hardlinking file "
              "(no in-memory concat, no full copy).")
        union_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(h5ad_files[largest_idx], union_path)
        except OSError:
            # cross-device or unsupported FS: fall back to a plain copy
            print("  Hardlink failed; falling back to file copy.")
            shutil.copyfile(h5ad_files[largest_idx], union_path)
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

    # Force CSR on-disk (backed per-sample subsets in annotation are only
    # selective for CSR), mirroring 1.1.1_preprocess.py::base_preprocessing().
    for mat_name in ("X", "counts"):
        mat = adata.X if mat_name == "X" else adata.layers.get("counts")
        if mat is None:
            continue
        if not sp.issparse(mat):
            mat = sp.csr_matrix(mat)
        else:
            mat = mat.tocsr()
        if mat_name == "X":
            adata.X = mat
        else:
            adata.layers["counts"] = mat

    union_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(union_path))
    return union_path


def main():
    parser = argparse.ArgumentParser(description="Build per-dataset sample chunk files.")
    parser.add_argument(
        "--test",
        action="store_true",
        default=False,
        help="Test mode: 1 sample per chunk (production default: 5).",
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

    # Delete chunk file folder recursively to ensure a perfectly clean start
    shutil.rmtree(path_output_chunks, ignore_errors=True)
    path_output_chunks.mkdir(parents=True, exist_ok=True)

    # Delete any stale union from a previous run (e.g. after --force
    # preprocess re-runs regenerated the views).
    shutil.rmtree(path_union_dir, ignore_errors=True)
    path_union_dir.mkdir(parents=True, exist_ok=True)

    # Number of samples per chunk
    chunk_size = 1 if args.test else 5

    # Build the per-dataset union (or fall back to the single view file)
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
