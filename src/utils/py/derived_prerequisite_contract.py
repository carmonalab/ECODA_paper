"""Semantic checks for Stage 2 derived preprocessing prerequisites.

These checks are intentionally small and backed/read-only: scheduler submitters
must not skip a hook solely because a checksum-valid artifact is nonempty.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def _read_h5ad(path: str | Path):
    import anndata as ad

    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size == 0:
        raise ValueError(f"h5ad is missing or empty: {artifact}")
    data = ad.read_h5ad(artifact, backed="r")
    return data


def _nonblank(series) -> bool:
    values = series.astype("string")
    return bool(values.notna().any() and (values.str.strip() != "").any())


def _validate_columns(data, path: Path, columns: tuple[str, ...]) -> None:
    if data.n_obs <= 0 or data.n_vars <= 0:
        raise ValueError(f"h5ad is empty: {path}")
    missing = [column for column in columns if column not in data.obs.columns]
    if missing:
        raise ValueError(f"{path}: missing required obs columns {missing}")
    blank = [column for column in columns if not _nonblank(data.obs[column])]
    if blank:
        raise ValueError(f"{path}: required obs columns are empty {blank}")


def validate_myocardial_counts(path: str | Path) -> None:
    artifact = Path(path)
    data = _read_h5ad(artifact)
    try:
        if "counts" not in data.layers:
            raise ValueError(f"{artifact}: layers['counts'] is missing")
        counts = data.layers["counts"]
        if getattr(counts, "shape", None) != (data.n_obs, data.n_vars):
            raise ValueError(f"{artifact}: counts shape does not match h5ad")
        import numpy as np

        values = counts.data if hasattr(counts, "data") else np.asarray(counts)
        if values.size == 0:
            raise ValueError(f"{artifact}: layers['counts'] is empty")
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.round(values)).all():
            raise ValueError(f"{artifact}: layers['counts'] must be nonempty nonnegative integers")
    finally:
        if getattr(data, "isbacked", False):
            data.file.close()

def validate_combinedpbmc(path: str | Path) -> None:
    artifact = Path(path)
    data = _read_h5ad(artifact)
    try:
        _validate_columns(data, artifact, ("Sample", "cond", "batch"))
    finally:
        if getattr(data, "isbacked", False):
            data.file.close()


def validate_joanito_debug(path: str | Path) -> None:
    artifact = Path(path)
    data = _read_h5ad(artifact)
    try:
        sample_column = "Sample" if "Sample" in data.obs.columns else "sample.ID"
        _validate_columns(data, artifact, (sample_column,))
        samples = data.obs[sample_column].astype("string")
        if samples.isna().any():
            raise ValueError(f"{artifact}: sample IDs must not be NA")
        samples = samples.str.strip()
        if (samples == "").any():
            raise ValueError(f"{artifact}: sample IDs must not be blank")
        counts = samples.value_counts(dropna=False)
        if len(counts) != 5 or not (counts == 500).all():
            raise ValueError(f"{artifact}: expected exactly five samples with 500 cells each")
    finally:
        if getattr(data, "isbacked", False):
            data.file.close()
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", required=True)
    parser.add_argument("--kind", required=True, choices=("myocardial", "combinedpbmc", "joanito-debug"))
    args = parser.parse_args()
    if args.kind == "myocardial":
        validate_myocardial_counts(args.path)
    elif args.kind == "combinedpbmc":
        validate_combinedpbmc(args.path)
    else:
        validate_joanito_debug(args.path)
    print(f"derived prerequisite contract OK: {args.path}")


if __name__ == "__main__":
    main()
