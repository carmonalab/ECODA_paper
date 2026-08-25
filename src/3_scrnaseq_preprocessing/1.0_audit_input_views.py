#!/usr/bin/env python3
"""Audit declared preprocessing views before sample-count filtering.

The audit deliberately stops before ``base_preprocessing``. It applies the
same registry subset and sample-ID standardization as the Stage 3 worker,
then records raw per-sample counts and configured high-resolution annotation
coverage. Biological labels are reported only as metadata; they never enter
selection or validation decisions.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path

import anndata as ad
import h5py
from anndata._io.h5ad import read_dataframe
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.datasets_io import read_datasets_json
from src.utils.py.preprocess_utils import apply_subset_vars, load_input


def read_obs_only(path: Path) -> ad.AnnData:
    """Read only the H5AD observation table, never the count matrix."""
    with h5py.File(path, "r") as handle:
        return ad.AnnData(obs=read_dataframe(handle["obs"]))


def load_audit_input(input_file, input_dir: Path, output_dir: Path):
    """Load only AnnData metadata when a validated raw cache is available."""
    if isinstance(input_file, list):
        metadata = []
        for name in input_file:
            source_path = input_dir / name
            if not source_path.exists() or not str(name).endswith(".h5ad"):
                return load_input(input_file, input_dir, output_dir)
            metadata.append(read_obs_only(source_path).obs)
        combined = pd.concat(metadata, axis=0)
        combined.index = [f"source_{position}" for position in range(len(combined))]
        return ad.AnnData(obs=combined)

    input_path = input_dir / input_file
    if str(input_file).endswith(".rds"):
        raw_cache = output_dir / f"{Path(input_file).stem}_raw.h5ad"
        if raw_cache.exists():
            return read_obs_only(raw_cache)
    elif str(input_file).endswith(".h5ad") and input_path.exists():
        return read_obs_only(input_path)
    return load_input(input_file, input_dir, output_dir)


SENTINEL_ANNOTATIONS = frozenset(
    {
        "",
        "na",
        "n/a",
        "nan",
        "none",
        "null",
        "missing",
        "not available",
        "not_available",
    }
)


def standardize_sample_ids(values: pd.Series, *, dataset: str, view: str) -> list[str]:
    """Validate and standardize sample IDs exactly as Stage 3 does."""
    standardized = []
    for position, value in enumerate(values):
        if pd.isna(value):
            raise ValueError(
                f"{dataset} / {view}: missing sample ID at observation {position}"
            )
        text = str(value)
        if not text.strip():
            raise ValueError(
                f"{dataset} / {view}: empty sample ID at observation {position}"
            )
        standardized.append(f"g{text}" if re.match(r"^\d", text) else text.replace("-", "_"))
    return standardized


def annotation_coverage(series: pd.Series) -> dict[str, int]:
    """Return total, valid, and missing/sentinel annotation counts."""
    total = int(series.shape[0])
    missing = series.isna()
    normalized = series.astype("string").str.strip().str.casefold()
    sentinel = normalized.isin(SENTINEL_ANNOTATIONS)
    invalid = missing | sentinel
    return {
        "total_cells": total,
        "valid_nonmissing": int((~invalid).sum()),
        "invalid_missing_or_sentinel": int(invalid.sum()),
    }


def audit_view(
    dataset: str,
    entry: dict,
    view_name: str,
    view: dict,
    input_root: Path,
    output_root: Path,
) -> dict:
    """Audit one registry view without running preprocessing."""
    input_file = view.get("input_file")
    if not input_file:
        raise ValueError(f"{dataset} / {view_name}: missing input_file_name")

    input_dir = input_root / dataset / "data"
    output_dir = output_root / dataset / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = load_audit_input(input_file, input_dir, output_dir)
    try:
        subset = apply_subset_vars(adata, view.get("subset_vars", {}), copy=False)
        if subset.n_obs == 0:
            raise ValueError(
                f"{dataset} / {view_name}: subset_vars produced an empty view: "
                f"{view.get('subset_vars', {})}"
            )

        sample_col = (view.get("columns") or {}).get("sample") or entry["sample_col"]
        if sample_col not in subset.obs.columns:
            raise KeyError(
                f"{dataset} / {view_name}: sample column {sample_col!r} is missing; "
                f"available columns: {list(subset.obs.columns)}"
            )
        sample_ids = standardize_sample_ids(
            subset.obs[sample_col], dataset=dataset, view=view_name
        )
        sample_counts = Counter(sample_ids)

        high_res_col = (view.get("columns") or {}).get(
            "cell_type_high_res"
        ) or entry["cell_type_high_res"]
        if high_res_col in subset.obs.columns:
            coverage = annotation_coverage(subset.obs[high_res_col])
            coverage["status"] = "present"
        else:
            coverage = {
                "total_cells": int(subset.n_obs),
                "valid_nonmissing": 0,
                "invalid_missing_or_sentinel": int(subset.n_obs),
                "status": "missing_column",
            }

        return {
            "dataset": dataset,
            "view": view_name,
            "input_file": input_file,
            "subset_vars": view.get("subset_vars", {}),
            "sample_column": sample_col,
            "sample_counts": dict(sorted(sample_counts.items())),
            "total_cells": int(subset.n_obs),
            "total_samples": int(len(sample_counts)),
            "high_resolution_column": high_res_col,
            "annotation_coverage": coverage,
        }
    finally:
        backing_file = getattr(adata, "file", None)
        if backing_file is not None:
            backing_file.close()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="datasets.json", type=Path)
    parser.add_argument(
        "--input-root",
        required=True,
        type=Path,
        help="Parent containing <dataset>/data staged inputs.",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        type=Path,
        help="Parent containing/receiving <dataset>/output raw caches.",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--view", default="benchmark_analysis")
    parser.add_argument("--ds-name", action="append", dest="datasets")
    parser.add_argument("--include-underscore", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_root = args.output_root or args.input_root
    config = read_datasets_json(args.config, view=args.view)
    selected = []
    for dataset in sorted(config):
        if args.datasets and dataset not in args.datasets:
            continue
        if not args.include_underscore and dataset.startswith("_"):
            continue
        entry = config[dataset]
        views = entry.get("views") or {}
        if args.view not in views:
            continue
        selected.append((dataset, entry, views[args.view]))
    if not selected:
        raise ValueError(f"No datasets declare view {args.view!r}")

    records = [
        audit_view(dataset, entry, args.view, view, args.input_root, output_root)
        for dataset, entry, view in selected
    ]
    payload = {
        "config": str(args.config.resolve()),
        "view": args.view,
        "threshold": 500,
        "datasets": records,
        "audit_scope": "raw_input_after_declared_subset_before_base_preprocessing",
    }
    write_json_atomic(args.output, payload)
    print(f"Wrote {args.output} ({len(records)} datasets)")


if __name__ == "__main__":
    main()
