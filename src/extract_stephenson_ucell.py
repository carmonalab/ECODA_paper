#!/usr/bin/env python3
"""Export standalone Stephenson benchmark UCell means using obs metadata only.

Usage: python src/extract_stephenson_ucell.py --h5ad /absolute/benchmark.h5ad \
    --output-dir /absolute/new_directory --config /path/to/datasets.json

The H5AD basename must match Stephenson's configured benchmark_analysis view.
Missing/blank annotations are excluded separately for each annotation; missing
scores are excluded from means, not cell counts. Absent groups and groups with
no observed score have NA activity (never invented zero); absent counts are zero.
Every matrix retains every sample, with exact sorted sample and cell-type names.
No expression arrays are opened, no UCell transformation is applied, and labels
are exported separately rather than used for aggregation. Publication exclusively
creates the output directory, moves Feather files, then atomically moves manifest.json
last as the commit marker. Existing paths are never replaced; failed publication
leaves an invalid partial directory without a manifest, preserved as failure evidence.
Existing H5AD checksum sidecars are recorded as unverified evidence, not verified
full-file checksums.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import numbers
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd

SCORES = ("IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell", "cellCycle.G2M_UCell")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _dataset(config: Path) -> tuple[dict, Path]:
    # A copied standalone script can reuse the checkout identified by --config.
    candidates = (config.parent / "src/utils/py/datasets_io.py",
                  Path(__file__).resolve().parent / "utils/py/datasets_io.py")
    helper = next((path for path in candidates if path.is_file()), None)
    if helper is None:
        raise ValueError("Cannot locate src/utils/py/datasets_io.py beside --config or script")
    spec = importlib.util.spec_from_file_location("stephenson_datasets_io", helper)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    dataset = module.read_datasets_json(str(config), view="benchmark_analysis").get("Stephenson")
    if not dataset or dataset["view_name"] != "benchmark_analysis":
        raise ValueError("Config requires Stephenson benchmark_analysis view")
    for key in ("sample_col", "label_col", "cell_type_low_res", "cell_type_high_res"):
        if not isinstance(dataset.get(key), str) or not dataset[key].strip():
            raise ValueError(f"Config requires Stephenson dataset-level {key}")
    return dataset, helper


def _strings(values: pd.Series, column: str, *, required: bool) -> pd.Series:
    missing = values.isna()
    if not values[~missing].map(lambda value: isinstance(value, str)).all():
        raise ValueError(f"obs[{column!r}] must contain string identifiers/labels")
    blank = values.astype(object).map(lambda value: isinstance(value, str) and not value.strip())
    if required and (missing | blank).any():
        raise ValueError(f"obs[{column!r}] contains missing/blank values")
    return values.astype(object).mask(missing | blank, None)


def aggregate(obs: pd.DataFrame, sample_column: str, label_column: str,
              annotations: dict[str, str]) -> tuple[dict[str, pd.DataFrame], dict]:
    """Aggregate already-selected obs; labels never participate in feature means."""
    obs = obs.copy()
    obs[sample_column] = _strings(obs[sample_column], sample_column, required=True)
    obs[label_column] = _strings(obs[label_column], label_column, required=True)
    sample_ids = sorted(obs[sample_column].unique().tolist())
    if not sample_ids:
        raise ValueError("Source has no cells/samples")
    labels = obs.groupby(sample_column, sort=False)[label_column]
    if not labels.nunique().eq(1).all():
        raise ValueError("Each sample must have exactly one unique nonmissing label")
    for score in SCORES:
        present = obs[score].dropna()
        if not present.map(lambda value: isinstance(value, numbers.Real) and
                           not isinstance(value, (bool, np.bool_))).all():
            raise ValueError(f"obs[{score!r}] contains non-numeric scores")
        obs[score] = obs[score].astype(float)
        if np.isinf(obs[score].to_numpy()).any():
            raise ValueError(f"obs[{score!r}] contains infinite scores")
    files = {"samples.feather": pd.DataFrame({
        "Sample": sample_ids, label_column: labels.first().reindex(sample_ids).tolist()
    })}
    excluded, axes = {}, {}
    for annotation, source_column in annotations.items():
        values = _strings(obs[source_column], source_column, required=False)
        valid = values.notna()
        excluded[annotation] = int((~valid).sum())
        cell_types = sorted(values[valid].unique().tolist())
        if "Sample" in cell_types:
            raise ValueError(f"Annotation {source_column!r} has reserved cell-type name 'Sample'")
        axes[annotation] = cell_types
        selected = obs.loc[valid, [sample_column, *SCORES]].copy()
        selected["__cell_type"] = values[valid]
        grouped = selected.groupby([sample_column, "__cell_type"], sort=False, observed=True)
        counts = grouped.size().unstack("__cell_type").reindex(
            index=sample_ids, columns=cell_types).fillna(0).astype("int64")
        counts.insert(0, "Sample", sample_ids)
        files[f"{annotation}__n_cells.feather"] = counts.reset_index(drop=True)
        means = grouped[list(SCORES)].mean()
        for score in SCORES:
            matrix = means[score].unstack("__cell_type").reindex(
                index=sample_ids, columns=cell_types)
            matrix.insert(0, "Sample", sample_ids)
            files[f"{annotation}__{score}.feather"] = matrix.reset_index(drop=True)
    return files, {"sample_ids": sample_ids, "n_cells": len(obs),
                   "excluded_annotation_cells": excluded, "cell_type_columns": axes}


def _publish(staging: Path, output: Path) -> None:
    """Claim a new directory exclusively; publish the manifest last to commit it.

    On failure, preserve the partial output without a manifest; retries must use
    a new output path. Staging and output share a parent so file renames are atomic.
    """
    output.mkdir()
    for feather in sorted(staging.glob("*.feather")):
        feather.rename(output / feather.name)
    (staging / "manifest.json").rename(output / "manifest.json")


def export(h5ad: Path, output_dir: Path, config: Path) -> dict:
    h5ad, output_dir, config = Path(h5ad), Path(output_dir), Path(config)
    if not h5ad.is_absolute() or not output_dir.is_absolute():
        raise ValueError("--h5ad and --output-dir must be absolute paths")
    if os.path.lexists(output_dir):
        raise FileExistsError(f"Output directory must not exist: {output_dir}")
    config = config.resolve(strict=True)
    dataset, helper = _dataset(config)
    if h5ad.name != dataset["output_file"]:
        raise ValueError("H5AD basename does not match configured Stephenson benchmark_analysis output")
    h5ad = h5ad.resolve(strict=True)
    before = h5ad.stat()
    annotations = {"authors_lowres": dataset["cell_type_low_res"],
                   "authors_highres": dataset["cell_type_high_res"],
                   "hitme_layer1": "layer1", "hitme_layer2": "layer2"}
    sample_column, label_column = dataset["sample_col"], dataset["label_col"]
    if label_column == "Sample":
        raise ValueError("Label column cannot use reserved output identifier 'Sample'")
    columns = list(dict.fromkeys([sample_column, label_column, *annotations.values(), *SCORES]))
    with h5py.File(h5ad, "r") as handle:
        if "obs" not in handle:
            raise ValueError("H5AD is missing obs")
        group = handle["obs"]
        missing = [column for column in columns if column not in group]
        if missing:
            raise ValueError(f"H5AD missing required obs fields: {missing}")
        selected = {column: pd.Series(read_elem(group[column])) for column in columns}
        lengths = {len(values) for values in selected.values()}
        index_name = group.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode("utf-8")
        if index_name not in group or len(lengths) != 1 or group[index_name].shape != (next(iter(lengths)),):
            raise ValueError("Selected obs fields must match the one-dimensional obs index")
        obs = pd.DataFrame(selected)
    files, summary = aggregate(obs, sample_column, label_column, annotations)
    # Canonical selected values only: deliberately never hash/read expression.
    digest = hashlib.sha256()
    for column in columns:
        values = obs[column].astype(object).where(obs[column].notna(), None).tolist()
        digest.update(json.dumps([column, values], ensure_ascii=False, allow_nan=False,
                                 separators=(",", ":"), default=lambda value: value.item()).encode("utf-8"))
        digest.update(b"\n")
    selected_digest = digest.hexdigest()
    sidecar = Path(f"{h5ad}.md5")
    source = {"h5ad_path": str(h5ad), "size": before.st_size, "mtime_ns": before.st_mtime_ns,
              "sidecar_text": sidecar.read_text() if sidecar.is_file() else None,
              "full_file_checksum_verified": False, "selected_obs_sha256": selected_digest}
    manifest = {"schema_version": 1, "dataset": "Stephenson", "view": "benchmark_analysis",
                "annotations": annotations, "scores": list(SCORES), "label_column": label_column,
                **summary, "source": source, "config": {"path": str(config), "sha256": sha256(config)},
                "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
                "dataset_accessor_sha256": sha256(helper),
                "runtime": {"executable": sys.executable, "python_version": platform.python_version()},
                "files": {}}
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        for filename, frame in files.items():
            frame.to_feather(staging / filename)
            manifest["files"][filename] = sha256(staging / filename)
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
        after = h5ad.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise ValueError("Source H5AD changed during extraction; output not published")
        _publish(staging, output_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("datasets.json"))
    args = parser.parse_args()
    try:
        manifest = export(args.h5ad, args.output_dir, args.config)
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, f"error: {error}\n")
    print(f"Exported {manifest['n_cells']} cells, {len(manifest['sample_ids'])} samples to {args.output_dir}")
    print("Excluded annotation cells: " + json.dumps(manifest["excluded_annotation_cells"], sort_keys=True))


if __name__ == "__main__":
    main()
