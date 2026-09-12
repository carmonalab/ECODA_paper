#!/usr/bin/env python3
"""Export one sample-level Feather table by reading only an H5AD ``obs`` group.

The exporter is deliberately independent of AnnData and the expression matrix:
``X``, ``raw``, and ``layers['counts']`` are never opened.  It retains the first
observation for each configured sample ID in source order, writes the output
and its strict MD5 sidecar atomically, and never publishes an H5AD artifact
record.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
import pandas as pd

try:
    from h5ad_source_identity import read_obs_column_values
except ImportError:  # package import in focused tests
    from .h5ad_source_identity import read_obs_column_values


READ_CHUNK_SIZE = 1024 * 1024


def _decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _obs_column_length(node: Any) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "keys") and "codes" in node:
        shape = node["codes"].shape
    if shape is None or len(shape) != 1:
        raise ValueError(f"H5AD obs column {node.name} is not one-dimensional")
    return int(shape[0])


def _normalise_sample(value: Any) -> str:
    value = _decode(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        raise ValueError("H5AD sample metadata contains a missing sample ID")
    text = str(value).strip()
    if not text or text.casefold() in {"nan", "none", "<na>"}:
        raise ValueError("H5AD sample metadata contains a blank sample ID")
    return text


def _feather_scalar(value: Any) -> Any:
    value = _decode(value)
    if value is None:
        return None
    if pd.isna(value) if not isinstance(value, (list, tuple, dict, np.ndarray)) else False:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _config_entry(config_path: Path, dataset: str, view: str) -> tuple[dict, dict]:
    import json

    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"cannot read datasets config: {config_path}") from exc
    entry = config.get(dataset)
    if not isinstance(entry, dict):
        raise ValueError(f"dataset is missing from config: {dataset}")
    views = entry.get("views")
    if not isinstance(views, dict) or view not in views:
        raise ValueError(f"dataset view is missing from config: {dataset}/{view}")
    return entry, views[view]


def _dataset_spec_columns(config_path: Path, dataset: str) -> list[str]:
    """Return optional candidate columns from the authoritative registry module."""
    module_path = config_path.parent / "notebooks" / "dataset_onboarding" / "dataset_specs.py"
    if not module_path.is_file():
        return []
    spec = importlib.util.spec_from_file_location("ecoda_dataset_specs", module_path)
    if spec is None or spec.loader is None:
        return []
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        dataset_specs = getattr(module, "DATASET_SPECS", {})
        batch_effect_specs = getattr(module, "BATCH_EFFECT_SPECS", {})
        dataset_spec = dataset_specs.get(dataset, {})
        batch_candidates = batch_effect_specs.get(dataset, [])
    except (AttributeError, ImportError, OSError, TypeError, ValueError):
        return []

    columns: list[str] = []

    def add(value: Any) -> None:
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, (list, tuple)):
            values = value
        else:
            return
        columns.extend(
            item for item in values
            if isinstance(item, str) and item.strip()
        )

    # These fields are the complete candidate registry in DATASET_SPECS.  The
    # biological column is included for technical-only/legacy entries whose
    # configured columns may not repeat it.
    if isinstance(dataset_spec, dict):
        add(dataset_spec.get("bio_col"))
        for key in (
            "sample_candidates",
            "sample_stable_cols",
            "batch_cols",
            "cell_type_candidates",
        ):
            add(dataset_spec.get(key, []))
    # Joanito and Stephenson are represented in the technical batch registry
    # even when they have no active DATASET_SPECS entry (notably Site).
    add(batch_candidates)
    return list(dict.fromkeys(columns))
def requested_columns(config_path: Path, dataset: str, entry: dict) -> tuple[str, list[str], list[str]]:
    columns = entry.get("columns")
    if not isinstance(columns, dict):
        raise ValueError(f"dataset columns are malformed: {dataset}")
    raw_sample_column = columns.get("sample")
    label_column = columns.get("label")
    if not isinstance(raw_sample_column, str) or not raw_sample_column.strip():
        raise ValueError(f"dataset sample column is missing: {dataset}")
    if not isinstance(label_column, str) or not label_column.strip():
        raise ValueError(f"dataset primary label column is missing: {dataset}")
    # Stage 3 persists the standardized Sample column. It is the only valid
    # grouping/order key for final Stage 5; the configured raw identity remains
    # an optional source metadata column when it is still present.
    required: list[str] = ["Sample", label_column]
    for key in ("batch", "cell_type_low_res", "cell_type_high_res"):
        value = columns.get(key)
        if isinstance(value, str) and value.strip():
            required.append(value)
        elif isinstance(value, (list, tuple)):
            required.extend(str(item) for item in value if isinstance(item, str) and item.strip())
    required = list(dict.fromkeys(required))
    optional = list(dict.fromkeys([raw_sample_column, *_dataset_spec_columns(config_path, dataset)]))
    optional = [column for column in optional if column not in required]
    all_columns = required + optional
    return "Sample", required, all_columns


def read_obs_metadata(
    input_path: Path,
    sample_column: str,
    required_columns: Iterable[str],
    candidate_columns: Iterable[str],
    chunk_size: int = READ_CHUNK_SIZE,
) -> pd.DataFrame:
    if not input_path.is_file() or input_path.stat().st_size <= 0:
        raise ValueError(f"H5AD input is missing or empty: {input_path}")
    if chunk_size <= 0:
        raise ValueError("chunk size must be positive")
    required = list(dict.fromkeys(str(column) for column in required_columns))
    with h5py.File(input_path, "r") as handle:
        if "obs" not in handle:
            raise ValueError(f"H5AD lacks obs dataframe: {input_path}")
        obs = handle["obs"]
        if obs.attrs.get("encoding-type", "") not in ("dataframe", b"dataframe"):
            raise ValueError(f"H5AD obs is not a dataframe: {input_path}")
        index_name = _decode(obs.attrs.get("_index", "_index"))
        if index_name not in obs:
            raise ValueError(f"H5AD obs index is missing: {input_path}")
        n_obs = _obs_column_length(obs[index_name])
        if n_obs <= 0:
            raise ValueError(f"H5AD obs is empty: {input_path}")
        missing = [column for column in required if column not in obs]
        if missing:
            raise ValueError(f"H5AD is missing configured obs columns: {missing}")
        present_optional = [column for column in candidate_columns if column in obs and column not in required]
        columns = required + present_optional
        values: dict[str, list[Any]] = {column: [] for column in columns}
        seen: dict[str, int] = {}
        sample_ids: list[str] = []
        for start in range(0, n_obs, chunk_size):
            stop = min(n_obs, start + chunk_size)
            chunk = {
                column: read_obs_column_values(obs, column, start, stop)
                for column in columns
            }
            for offset, raw_sample in enumerate(chunk[sample_column]):
                sample = _normalise_sample(raw_sample)
                if sample not in seen:
                    seen[sample] = len(sample_ids)
                    sample_ids.append(sample)
                    for column in columns:
                        value = sample if column == sample_column else chunk[column][offset]
                        values[column].append(_feather_scalar(value))
        if not sample_ids:
            raise ValueError(f"H5AD has no non-empty sample IDs: {input_path}")
    frame = pd.DataFrame(values)
    if sample_column != "Sample":
        frame.insert(0, "Sample", frame[sample_column])
    else:
        frame.insert(0, "Sample", frame.pop("Sample"))
    if frame["Sample"].duplicated().any() or frame["Sample"].astype(str).str.strip().eq("").any():
        raise ValueError("exported sample metadata has blank or duplicate Sample IDs")
    return frame


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_sidecar(path: Path) -> dict[str, str]:
    sidecar = Path(f"{path}.md5")
    lines = sidecar.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3 or [line.split("=", 1)[0] for line in lines] != ["MD5", "SIZE", "PATH"]:
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}")
    fields = {line.split("=", 1)[0]: line.split("=", 1)[1] for line in lines}
    if fields["PATH"] != str(path) or fields["SIZE"] != str(path.stat().st_size):
        raise ValueError(f"metadata checksum sidecar path/size mismatch: {sidecar}")
    if fields["MD5"] != _md5(path):
        raise ValueError(f"metadata checksum sidecar digest mismatch: {sidecar}")
    return fields


def validate_output(path: Path, required_columns: Iterable[str]) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"metadata Feather is missing or empty: {path}")
    _read_sidecar(path)
    try:
        frame = pd.read_feather(path)
    except Exception as exc:  # pyarrow surfaces several exception classes
        raise ValueError(f"metadata Feather is unreadable: {path}") from exc
    required = list(dict.fromkeys(str(column) for column in required_columns))
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"metadata Feather is missing configured columns: {missing}")
    if frame.empty or "Sample" not in frame.columns:
        raise ValueError(f"metadata Feather has no Sample rows: {path}")
    ids = frame["Sample"].astype(str)
    if ids.str.strip().eq("").any() or ids.duplicated().any():
        raise ValueError(f"metadata Feather has blank or duplicate Sample IDs: {path}")


def write_metadata(frame: pd.DataFrame, output: Path) -> None:
    if not output.is_absolute():
        raise ValueError("metadata output path must be absolute")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{output}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    try:
        frame.reset_index(drop=True).to_feather(temporary)
        os.replace(temporary, output)
        digest = _md5(output)
        sidecar_tmp.write_text(
            f"MD5={digest}\nSIZE={output.stat().st_size}\nPATH={output}\n",
            encoding="utf-8",
        )
        os.replace(sidecar_tmp, sidecar)
    finally:
        for path in (temporary, sidecar_tmp):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
    validate_output(output, frame.columns)


def export(args: argparse.Namespace) -> None:
    config_path = args.config.resolve()
    entry, _view = _config_entry(config_path, args.dataset, args.view)
    sample_column, required, candidates = requested_columns(config_path, args.dataset, entry)
    output = args.output.resolve()
    if args.check:
        validate_output(output, required)
        return
    if output.exists() or Path(f"{output}.md5").exists():
        try:
            validate_output(output, required)
            print(f"METADATA_EXPORT=NOOP_VALIDATED PATH={output}")
            return
        except (OSError, ValueError):
            # Rebuild only an invalid output; the source H5AD remains read-only.
            pass
    frame = read_obs_metadata(
        args.input_file.resolve(), sample_column, required, candidates, args.chunk_size
    )
    write_metadata(frame, output)
    print(f"METADATA_EXPORT=OK PATH={output} SAMPLES={len(frame)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--view", required=True)
    parser.add_argument("--input-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=READ_CHUNK_SIZE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        export(args)
    except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
