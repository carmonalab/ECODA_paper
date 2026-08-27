"""Fail-closed schema, key, and coverage checks for dual annotations."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Iterable

import h5py
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.py.h5ad_source_identity import (  # noqa: E402
    read_obs_column_values,
    read_str_dataset,
)

HITME_REQUIRED_COLUMNS = {"layer1", "layer2", "layer3"}
SCATOMIC_REQUIRED_COLUMNS = {
    "layer_1",
    "layer_2",
    "layer_3",
    "layer_4",
    "layer_5",
    "layer_6",
    "scATOMIC_pred",
    "classification_confidence",
    "S.Score",
    "G2M.Score",
    "Phase",
}
REQUIRED_ANNOTATION_COLUMNS = HITME_REQUIRED_COLUMNS | SCATOMIC_REQUIRED_COLUMNS
OPTIONAL_HITME_COLUMNS = {
    "IFN_UCell",
    "HeatShock_UCell",
    "cellCycle.G1S_UCell",
    "cellCycle.G2M_UCell",
}
KNOWN_ANNOTATION_COLUMNS = REQUIRED_ANNOTATION_COLUMNS | OPTIONAL_HITME_COLUMNS
REQUIRED_KEYS = {"cell_barcode", "Sample"}


def _artifact_path(path: str | Path) -> Path:
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size == 0:
        raise ValueError(f"annotation artifact is missing or empty: {artifact}")
    return artifact


def validate_sidecar(path: str | Path) -> dict:
    """Validate the adjacent MD5/SIZE/PATH sidecar and return its fields."""
    artifact = _artifact_path(path)
    sidecar = Path(f"{artifact}.md5")
    if not sidecar.is_file() or sidecar.stat().st_size == 0:
        raise ValueError(f"annotation artifact checksum sidecar is missing: {sidecar}")
    fields: dict[str, str] = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields[key] = value
    digest = hashlib.md5()
    with artifact.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    if fields.get("MD5", "").lower() != digest.hexdigest():
        raise ValueError(f"annotation artifact MD5 mismatch: {artifact}")
    size = str(artifact.stat().st_size)
    if fields.get("SIZE") != size:
        raise ValueError(f"annotation artifact SIZE mismatch: {artifact}")
    if fields.get("PATH") != str(artifact):
        raise ValueError(f"annotation artifact PATH mismatch: {artifact}")
    return fields


def _is_nonblank(value: object) -> bool:
    if value is None:
        return False
    try:
        missing = pd.isna(value)
        if bool(missing):
            return False
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return bool(text) and text.lower() not in {"nan", "<na>", "na", "none"}


def _as_text(value: object) -> str | None:
    return str(value) if _is_nonblank(value) else None


def _normalize_expected_keys(expected_keys: Iterable[object] | None) -> set[tuple[str, str]] | None:
    if expected_keys is None:
        return None
    normalized: set[tuple[str, str]] = set()
    for item in expected_keys:
        if isinstance(item, (tuple, list)) and len(item) == 2:
            sample, barcode = item
        elif isinstance(item, str):
            if "\t" in item:
                sample, barcode = item.split("\t", 1)
            elif "|" in item:
                sample, barcode = item.split("|", 1)
            else:
                raise ValueError(
                    "expected_keys strings must use '<Sample>\\t<cell_barcode>'"
                )
        else:
            raise ValueError(f"unsupported expected annotation key: {item!r}")
        if not _is_nonblank(sample) or not _is_nonblank(barcode):
            raise ValueError(f"expected annotation keys must be nonblank: {item!r}")
        key = (str(sample), str(barcode))
        if key in normalized:
            raise ValueError(f"expected annotation keys contain duplicates: {key!r}")
        normalized.add(key)
    return normalized


def _coverage_stats(frame: pd.DataFrame, samples: pd.Series) -> tuple[dict, dict]:
    column_stats: dict[str, dict] = {}
    for column in sorted(REQUIRED_ANNOTATION_COLUMNS | (set(frame.columns) & OPTIONAL_HITME_COLUMNS)):
        values = frame[column]
        valid = values.astype(object).map(_is_nonblank).astype(bool)
        column_stats[column] = {
            "n_rows": int(len(values)),
            "n_nonblank": int(valid.sum()),
            "coverage": float(valid.mean()) if len(values) else 0.0,
        }
    sample_stats: dict[str, dict] = {}
    sample_text = samples.map(lambda value: str(value))
    for sample in sorted(sample_text.unique()):
        mask = sample_text == sample
        per_column = {
            column: int(frame.loc[mask, column].astype(object).map(_is_nonblank).astype(bool).sum())
            for column in sorted(REQUIRED_ANNOTATION_COLUMNS | (set(frame.columns) & OPTIONAL_HITME_COLUMNS))
        }
        sample_stats[sample] = {
            "n_rows": int(mask.sum()),
            "columns_nonblank": per_column,
        }
    return column_stats, sample_stats


def _validate_numeric_columns(frame: pd.DataFrame, path: Path) -> None:
    numeric_columns = {
        "classification_confidence",
        "S.Score",
        "G2M.Score",
    } | (OPTIONAL_HITME_COLUMNS & {"IFN_UCell", "HeatShock_UCell",
                                   "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"})
    for column in sorted(numeric_columns & set(frame.columns)):
        values = frame[column]
        converted = pd.to_numeric(values, errors="coerce")
        invalid = values.notna() & converted.isna()
        if bool(invalid.any()):
            raise ValueError(f"annotation column {column!r} contains nonnumeric values: {path}")
        finite = converted.dropna()
        if not finite.empty:
            import numpy as np
            if not np.isfinite(finite.to_numpy(dtype=float)).all():
                raise ValueError(f"annotation column {column!r} contains nonfinite values: {path}")

def _validate_frame(
    frame: pd.DataFrame,
    *,
    path: Path,
    expected_keys: Iterable[object] | None,
    require_dataset_anchor: bool,
) -> dict:
    if frame.empty:
        raise ValueError(f"annotation artifact has no rows: {path}")
    missing_keys = REQUIRED_KEYS - set(frame.columns)
    if missing_keys:
        raise ValueError(f"annotation artifact missing key columns {sorted(missing_keys)}: {path}")
    missing_columns = REQUIRED_ANNOTATION_COLUMNS - set(frame.columns)
    if missing_columns:
        raise ValueError(
            f"annotation artifact missing required dual-method columns "
            f"{sorted(missing_columns)}: {path}"
        )
    _validate_numeric_columns(frame, path)

    samples = frame["Sample"].map(_as_text)
    barcodes = frame["cell_barcode"].map(_as_text)
    if samples.isna().any() or barcodes.isna().any():
        raise ValueError(f"annotation artifact contains blank Sample/cell_barcode keys: {path}")
    keys = list(zip(samples.astype(str), barcodes.astype(str)))
    duplicate_count = len(keys) - len(set(keys))
    if duplicate_count:
        raise ValueError(
            f"annotation artifact has {duplicate_count} duplicate (Sample, cell_barcode) keys: {path}"
        )

    expected = _normalize_expected_keys(expected_keys)
    actual = set(keys)
    if expected is not None:
        missing = expected - actual
        extra = actual - expected
        if missing or extra:
            raise ValueError(
                f"annotation key coverage mismatch for {path}: "
                f"missing={len(missing)}, extra={len(extra)}"
            )

    anchors = {}
    if require_dataset_anchor:
        for column in ("layer1", "scATOMIC_pred"):
            n_nonblank = int(frame[column].astype(object).map(_is_nonblank).astype(bool).sum())
            anchors[column] = n_nonblank
            if n_nonblank == 0:
                raise ValueError(
                    f"annotation dataset anchor {column!r} has no nonblank output: {path}"
                )

    column_stats, sample_stats = _coverage_stats(frame, samples.astype(str))
    return {
        "path": str(path),
        "n_rows": int(len(frame)),
        "n_samples": int(samples.nunique()),
        "required_columns": sorted(REQUIRED_ANNOTATION_COLUMNS),
        "optional_columns_present": sorted(OPTIONAL_HITME_COLUMNS & set(frame.columns)),
        "column_coverage": column_stats,
        "sample_coverage": sample_stats,
        "anchors": anchors,
        "keys": sorted(f"{sample}\t{barcode}" for sample, barcode in actual),
    }


def _persisted_shape(node: object) -> tuple[int, ...]:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "attrs"):
        shape = node.attrs.get("shape")
    if shape is None:
        return ()
    try:
        return tuple(int(value) for value in shape)
    except (TypeError, ValueError):
        return ()


def _read_h5ad_obs_frame(
    path: str | Path,
    columns: Iterable[str] | None = None,
) -> tuple[tuple[int, ...], pd.DataFrame]:
    """Read H5AD shape and selected obs columns without anndata backed-open."""
    artifact = _artifact_path(path)
    with h5py.File(artifact, "r") as handle:
        if "X" not in handle:
            raise ValueError(f"h5ad lacks X: {artifact}")
        shape = _persisted_shape(handle["X"])
        if len(shape) != 2 or any(value <= 0 for value in shape):
            raise ValueError(f"h5ad X is empty or has no persisted shape: {artifact}")
        if "obs" not in handle:
            raise ValueError(f"h5ad lacks obs dataframe: {artifact}")
        obs = handle["obs"]
        encoding = obs.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8")
        if encoding != "dataframe":
            raise ValueError(f"h5ad obs encoding is invalid: {artifact}")
        index_name = obs.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode("utf-8")
        if index_name not in obs:
            raise ValueError(f"h5ad lacks obs index: {artifact}")
        barcodes = read_str_dataset(obs[index_name])
        if len(barcodes) != shape[0]:
            raise ValueError(f"h5ad obs/X row counts disagree: {artifact}")
        requested = set(columns) if columns is not None else {
            str(name) for name in obs.keys() if name != index_name
        }
        values: dict[str, object] = {"cell_barcode": barcodes}
        for column in sorted(requested):
            if column == index_name or column not in obs:
                continue
            values[column] = read_obs_column_values(obs, column)
            if len(values[column]) != shape[0]:
                raise ValueError(f"h5ad obs column length mismatch for {column}: {artifact}")
        return shape, pd.DataFrame(values)


def validate_feather(
    path: str | Path,
    *,
    expected_keys: Iterable[object] | None = None,
    require_dataset_anchor: bool = False,
    require_sidecar: bool = False,
) -> dict:
    """Validate one annotation feather, optionally against an expected key set."""
    artifact = _artifact_path(path)
    if require_sidecar or Path(f"{artifact}.md5").exists():
        validate_sidecar(artifact)
    import pyarrow.feather as feather

    frame = feather.read_table(artifact).to_pandas()
    return _validate_frame(
        frame,
        path=artifact,
        expected_keys=expected_keys,
        require_dataset_anchor=require_dataset_anchor,
    )


def h5ad_key_list(
    path: str | Path, samples: Iterable[object] | None = None
) -> list[tuple[str, str]]:
    """Read ordered (Sample, cell_barcode) keys without anndata backed-open."""
    artifact = _artifact_path(path)
    _shape, frame = _read_h5ad_obs_frame(artifact, {"Sample"})
    if "Sample" not in frame.columns:
        raise ValueError(f"h5ad lacks obs['Sample']: {artifact}")
    selected = None if samples is None else {str(value) for value in samples}
    keys_list: list[tuple[str, str]] = []
    for sample, barcode in zip(frame["Sample"], frame["cell_barcode"]):
        sample_text = str(sample) if _is_nonblank(sample) else ""
        barcode_text = str(barcode) if _is_nonblank(barcode) else ""
        if selected is None or sample_text in selected:
            if not sample_text or not barcode_text:
                raise ValueError(f"h5ad has blank annotation key: {artifact}")
            keys_list.append((sample_text, barcode_text))
    if len(keys_list) != len(set(keys_list)):
        raise ValueError(f"h5ad has duplicate (Sample, cell_barcode) keys: {artifact}")
    return keys_list


def h5ad_keys(path: str | Path, samples: Iterable[object] | None = None) -> set[tuple[str, str]]:
    """Read composite keys from a union/view without requiring annotation columns."""
    return set(h5ad_key_list(path, samples))


def validate_h5ad(
    path: str | Path,
    *,
    expected_keys: Iterable[object] | None = None,
    require_dataset_anchor: bool = True,
    require_sidecar: bool = False,
) -> dict:
    """Validate merged h5ad annotations without anndata backed-open."""
    artifact = _artifact_path(path)
    if require_sidecar or Path(f"{artifact}.md5").exists():
        validate_sidecar(artifact)
    shape, frame = _read_h5ad_obs_frame(
        artifact, REQUIRED_KEYS | KNOWN_ANNOTATION_COLUMNS
    )
    if len(shape) != 2 or shape[0] <= 0 or shape[1] <= 0:
        raise ValueError(f"annotated h5ad is empty: {artifact}")
    return _validate_frame(
        frame,
        path=artifact,
        expected_keys=expected_keys,
        require_dataset_anchor=require_dataset_anchor,
    )
def _keys_from_file(path: str | Path) -> set[tuple[str, str]]:
    key_file = Path(path)
    if not key_file.is_file():
        raise ValueError(f"expected annotation key file is missing: {key_file}")
    return _normalize_expected_keys(
        [line for line in key_file.read_text().splitlines() if line]
    ) or set()

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--path")
    group.add_argument("--h5ad")
    parser.add_argument(
        "--require-sidecar",
        action="store_true",
        help="also require and validate the adjacent MD5/SIZE/PATH sidecar",
    )
    parser.add_argument(
        "--expected-keys-file",
        help="headerless Sample<TAB>cell_barcode keys expected in the artifact",
    )
    parser.add_argument(
        "--expected-union",
        help="minimal union h5ad from which expected annotation keys are read",
    )
    parser.add_argument(
        "--expected-sample",
        action="append",
        default=[],
        help="restrict --expected-union keys to this sample (repeatable)",
    )
    parser.add_argument(
        "--require-dataset-anchor",
        action="store_true",
        help="require nonblank layer1 and scATOMIC_pred anchors",
    )
    args = parser.parse_args()
    artifact = args.h5ad or args.path
    if args.require_sidecar:
        validate_sidecar(artifact)
    expected = None
    if args.expected_keys_file:
        expected = _keys_from_file(args.expected_keys_file)
    if args.expected_union:
        expected = h5ad_keys(args.expected_union, args.expected_sample or None)
    if args.h5ad:
        stats = validate_h5ad(
            args.h5ad,
            expected_keys=expected,
            require_dataset_anchor=True,
        )
        print(f"annotated h5ad contract OK: {args.h5ad} ({stats['n_rows']} rows)")
    else:
        stats = validate_feather(
            args.path,
            expected_keys=expected,
            require_dataset_anchor=args.require_dataset_anchor,
        )
        print(f"annotation feather contract OK: {args.path} ({stats['n_rows']} rows)")


if __name__ == "__main__":
    main()
