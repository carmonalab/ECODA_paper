#!/usr/bin/env python3
"""Export Adams Authors_HR versus Leiden_res_0.4 counts; read only selected obs.

Run with the repository Pixi Python, --h5ad /mounted/NAS/benchmark.h5ad,
--output-dir /new/run-owned/directory. No expression or embedding values are
read. Existing output directories are never replaced. The manifest records
selected-metadata SHA256, not a freshly verified whole-H5AD checksum.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "utils/py"))
from datasets_io import read_datasets_json
from h5ad_obs_free import load_h5ad_obs_free


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=Path("datasets.json"))
    args = parser.parse_args()
    dataset = read_datasets_json(args.config, view="benchmark_analysis")["Adams"]
    source = args.h5ad.resolve(strict=True)
    if source.name != dataset["output_file"]:
        raise ValueError("Input must match the configured Adams benchmark_analysis H5AD")
    author = dataset["cell_type_high_res"]
    leiden = "leiden_res_0.4_benchmark_analysis_hvg2000"
    before = source.stat()
    obs = load_h5ad_obs_free(source, [author, leiden])
    if obs.empty or obs.index.has_duplicates:
        raise ValueError("Expected nonempty observations with unique cell IDs")
    for column in (author, leiden):
        if obs[column].isna().any() or obs[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing or blank annotations in {column}")
    counts = pd.crosstab(obs[author].astype(str), obs[leiden].astype(str))
    counts.index.name = "Authors_HR"
    counts.columns.name = "Leiden_res_0.4"
    if int(counts.to_numpy().sum()) != len(obs):
        raise ValueError("Contingency counts must include every source cell exactly once")
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Source changed during metadata extraction")
    payload = counts.to_csv().encode()
    digest = hashlib.md5(payload).hexdigest()
    metadata_sha256 = hashlib.sha256(obs.to_csv().encode()).hexdigest()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    target = output / "counts.csv"
    temporary = output / "counts.csv.tmp"
    temporary.write_bytes(payload)
    temporary.rename(target)
    (output / "counts.csv.md5").write_text(
        f"MD5={digest}\nSIZE={len(payload)}\nPATH={target}\n"
    )
    manifest = {
        "dataset": "Adams", "view": "benchmark_analysis",
        "source_h5ad": str(source), "source_size": before.st_size,
        "source_mtime_ns": before.st_mtime_ns,
        "author_column": author, "leiden_column": leiden,
        "n_cells": len(obs), "n_author_types": len(counts),
        "n_leiden_clusters": len(counts.columns),
        "selected_obs_sha256": metadata_sha256,
        "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "counts_md5": digest, "counts_size": len(payload),
        "scope": "selected obs columns and cell IDs only; no expression/embedding values",
    }
    temporary = output / "manifest.json.tmp"
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.rename(output / "manifest.json")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
