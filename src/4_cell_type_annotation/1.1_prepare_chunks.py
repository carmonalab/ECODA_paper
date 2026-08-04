#!/usr/bin/env python3
#
# 1.1_prepare_chunks.py — Build dataset chunks (production or test mode)
#
# Native Python replacement for the former reticulate-based 1.1_prepare_chunks.r.
# Reads each preprocessed .h5ad in backed mode, extracts unique sample IDs from
# the SAMPLE_COLNAME column, groups them into consecutive chunks of 5 (or 1 in
# --test mode) and writes chunk_<N>.txt files under
# ${HPC_SCRATCH_DIR}/${DS_NAME}/output/chunks. Line 1 of each chunk file is the
# absolute .h5ad path, subsequent lines are sample IDs.

import argparse
import os
import shutil
import sys
from pathlib import Path

import anndata as ad


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
    path_output_chunks = path_data / "chunks"

    print(f"Path is: {path_data}")

    # Skip rds->h5ad conversion caches written into the output dir by
    # preprocess_utils.load_single_input() (they lack the standardized sample
    # column and are not preprocessed views).
    h5ad_files = sorted(
        f for f in path_data.glob("*.h5ad") if not f.name.endswith("_raw.h5ad")
    )
    print(f"Files found: {', '.join(f.name for f in h5ad_files)}")

    # Delete chunk file folder recursively to ensure a perfectly clean start
    shutil.rmtree(path_output_chunks, ignore_errors=True)
    path_output_chunks.mkdir(parents=True, exist_ok=True)

    # Number of samples per chunk
    chunk_size = 1 if args.test else 5
    global_chunk_counter = 1

    # Loop through each file individually
    for f in h5ad_files:
        print(f"Processing file-specific chunks for: {f.name}")

        # 1. Read the h5ad file in backed mode and extract its unique samples
        adata = ad.read_h5ad(f, backed="r")
        if sample_col not in adata.obs.columns:
            raise KeyError(f"Sample column '{sample_col}' not found in obs of {f}")

        file_samples = adata.obs[sample_col].astype(str).unique()

        # 2. Split ONLY this file's samples into consecutive groups of chunk_size
        sample_groups = [
            file_samples[i : i + chunk_size] for i in range(0, len(file_samples), chunk_size)
        ]

        # 3. Write each group to a unique chunk file
        for group in sample_groups:
            # We write the source file as the VERY FIRST line, followed by the sample IDs
            chunk_path = path_output_chunks / f"chunk_{global_chunk_counter}.txt"
            chunk_path.write_text(
                "\n".join([str(f.resolve())] + [str(s) for s in group]) + "\n"
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
