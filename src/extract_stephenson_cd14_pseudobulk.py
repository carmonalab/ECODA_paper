#!/usr/bin/env python3
"""Extract run-owned, raw-count CD14 low-resolution Stephenson pseudobulk."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

FILENAME = "StephensonE_2021_33879890_preprocessed_benchmark_analysis_ECODAprocessed.h5ad"
INT64_MAX = np.iinfo(np.int64).max
R_INT_MAX = np.iinfo(np.int32).max


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_stat(path):
    stat = path.stat()
    return {"device": stat.st_dev, "inode": stat.st_ino, "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "ctime_ns": stat.st_ctime_ns}


def absolute_path(value):
    path = Path(value)
    if not path.is_absolute():
        raise ValueError(f"Path must be absolute: {value}")
    return path


def extract(h5ad, output_dir, config, helpers_root):
    h5ad, output_dir, config, helpers_root = map(
        absolute_path, (h5ad, output_dir, config, helpers_root)
    )
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if h5ad.name != FILENAME:
        raise ValueError(f"Expected Stephenson benchmark filename {FILENAME}")
    sys.path.insert(0, str(helpers_root / "src" / "utils" / "py"))
    from datasets_io import read_datasets_json
    from h5ad_pseudobulk import (
        _open_obs_and_counts, _read_count_chunk, _validate_obs_columns,
        read_obs_column_values,
    )

    entry = read_datasets_json(str(config), view="benchmark_analysis").get("Stephenson")
    expected = {"sample_col": "Sample", "label_col": "Status",
                "cell_type_low_res": "initial_clustering", "output_file": FILENAME,
                "view_name": "benchmark_analysis"}
    if entry is None or any(entry.get(key) != value for key, value in expected.items()):
        raise ValueError("Config must declare Stephenson benchmark Sample/Status/initial_clustering")
    before = source_stat(h5ad)
    selected_hash = hashlib.sha256()
    with h5py.File(h5ad, "r") as handle:
        obs, counts, n_obs, n_genes, genes, index_name = _open_obs_and_counts(handle, h5ad)
        columns = ["Sample", "Status", "initial_clustering"]
        _validate_obs_columns(obs, columns, n_obs)
        values = {column: read_obs_column_values(obs, column, 0, n_obs) for column in columns}
        for column, array in values.items():
            if any(not isinstance(value, str) or not value.strip() or value != value.strip()
                   for value in array):
                raise ValueError(f"Missing, invalid or padded obs {column} identifier")
        samples = np.asarray(values["Sample"], dtype=str)
        labels = np.asarray(values["Status"], dtype=str)
        annotations = np.asarray(values["initial_clustering"], dtype=str)
        if not set(labels).issubset({"Covid", "Healthy"}):
            raise ValueError("Status must contain only Covid/Healthy")
        sample_ids = sorted(set(samples))
        if "gene" in sample_ids:
            raise ValueError("Sample identifier collides with gene column")
        selected = np.flatnonzero(annotations == "CD14")
        n_cells = []
        status = []
        for sample in sample_ids:
            sample_mask = samples == sample
            sample_labels = set(labels[sample_mask])
            if len(sample_labels) != 1:
                raise ValueError(f"Inconsistent Status across all cells for sample {sample}")
            count = int(np.count_nonzero(sample_mask & (annotations == "CD14")))
            if count < 5:
                raise ValueError(f"Sample {sample} has {count} CD14 cells; minimum is 5")
            n_cells.append(count)
            status.append(next(iter(sample_labels)))
        obs_ids = read_obs_column_values(obs, index_name, 0, n_obs)
        for row in selected:
            selected_hash.update((json.dumps([int(row), str(obs_ids[row]), samples[row],
                                             labels[row], "CD14"], ensure_ascii=False,
                                            separators=(",", ":")) + "\n").encode())
        aggregate = np.zeros((n_genes, len(sample_ids)), dtype=np.int64)
        sample_index = {sample: index for index, sample in enumerate(sample_ids)}
        # Only read contiguous selected runs, bounded to 4096 cells; never X.
        splits = np.flatnonzero(np.diff(selected) != 1) + 1
        for run in np.split(selected, splits):
            for offset in range(0, len(run), 4096):
                rows = run[offset:offset + 4096]
                start, stop = int(rows[0]), int(rows[-1]) + 1
                data, indices, indptr = _read_count_chunk(counts, start, stop, n_genes, n_obs=n_obs)
                for local, row in enumerate(rows):
                    lo, hi = int(indptr[local]), int(indptr[local + 1])
                    idx, val = indices[lo:hi], data[lo:hi]
                    if not len(val):
                        continue
                    order = np.argsort(idx, kind="stable")
                    idx, val = idx[order], val[order]
                    boundaries = np.r_[0, np.flatnonzero(np.diff(idx)) + 1]
                    # Check reduction capacity before NumPy int64 addition (duplicate CSR indices).
                    if int(val.max()) * len(val) > INT64_MAX:
                        ends = np.r_[boundaries[1:], len(val)]
                        sums = [sum(map(int, val[a:b])) for a, b in zip(boundaries, ends)]
                        if any(value > INT64_MAX for value in sums):
                            raise OverflowError("Duplicate count sum exceeds int64 capacity")
                        val = np.asarray(sums, dtype=np.int64)
                    else:
                        val = np.add.reduceat(val, boundaries)
                    idx = idx[boundaries]
                    column = sample_index[samples[row]]
                    current = aggregate[idx, column]
                    if np.any(current > INT64_MAX - val):
                        raise OverflowError("Pseudobulk sum exceeds int64 capacity")
                    aggregate[idx, column] = current + val
        if np.any(aggregate > R_INT_MAX):
            raise OverflowError("Pseudobulk counts exceed maximum R integer")
    if source_stat(h5ad) != before:
        raise RuntimeError("Source H5AD changed during extraction")
    metadata = pd.DataFrame({"Sample": sample_ids, "Status": status,
                             "n_cells": np.asarray(n_cells, dtype=np.int32)})
    total_umis = [sum(map(int, aggregate[:, index])) for index in range(len(sample_ids))]
    sidecars = []
    for suffix in (".md5", ".sha256"):
        path = Path(str(h5ad) + suffix)
        if path.is_file():
            sidecars.append({"path": str(path), "sha256": sha256(path),
                             "content": path.read_text().strip(), "verified": False})
    manifest = {
        "schema_version": 1, "dataset": "Stephenson", "view": "benchmark_analysis",
        "annotation_column": "initial_clustering", "cell_type": "CD14",
        "sample_ids": sample_ids, "n_selected_cells": len(selected), "n_genes": n_genes,
        "total_umis": sum(total_umis), "sample_total_umis": dict(zip(sample_ids, total_umis)),
        "source": {"path": str(h5ad), "stat": before, "full_file_checksum_verified": False,
                   "sidecars": sidecars, "selected_obs_sha256": selected_hash.hexdigest(),
                   "selected_obs_hash_encoding": "UTF-8 JSONL [row,index,Sample,Status,CD14]"},
        "config_sha256": sha256(config), "extraction_code_sha256": sha256(__file__),
    }
    output_dir.mkdir(parents=False, exist_ok=False)
    frame = pd.DataFrame(aggregate.astype(np.int32), columns=sample_ids)
    frame.insert(0, "gene", genes)
    frame.to_feather(output_dir / "counts.feather")
    metadata.to_feather(output_dir / "samples.feather")
    manifest["files"] = {name: sha256(output_dir / name)
                         for name in ("counts.feather", "samples.feather")}
    if source_stat(h5ad) != before:
        raise RuntimeError("Source H5AD changed before publication")
    (output_dir / "aggregation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("h5ad", "output-dir", "config", "helpers-root"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    extract(args.h5ad, args.output_dir, args.config, args.helpers_root)


if __name__ == "__main__":
    main()
