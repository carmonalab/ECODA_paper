"""Memory-lean H5AD Sample readers and run-owned source identities.

This module intentionally uses h5py only.  Opening an anndata 0.12.x H5AD in
backed mode can eagerly materialize the complete ``layers`` group, so source
sample-ID validation must not use ``anndata.read_h5ad``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np

IDENTITY_SCHEMA = 1
READ_CHUNK_SIZE = 1024 * 1024
_MD5_LENGTH = 32


def _decode_attr(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def read_str_dataset(dataset: Any, start: int | None = None, stop: int | None = None) -> np.ndarray:
    """Read one HDF5 string/numeric slice as a one-dimensional string array."""
    if start is None and stop is None:
        values = dataset[:]
    else:
        values = dataset[start:stop]
    values = np.asarray(values)
    if values.dtype.kind == "S":
        return np.char.decode(values, "utf-8")
    if values.dtype.kind == "U":
        return values.astype(str)
    if values.dtype.kind == "O":
        output = np.empty(values.shape, dtype=object)
        for index, value in enumerate(values):
            if isinstance(value, bytes):
                output[index] = value.decode("utf-8")
            else:
                output[index] = str(value)
        return output.astype(str)
    if values.dtype.kind in "biuf":
        return values.astype(str)
    raise RuntimeError(
        f"Unsupported dtype '{values.dtype}' for HDF5 dataset '{dataset.name}'."
    )


def _obs_column_length(node: Any) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "keys") and "codes" in node:
        shape = node["codes"].shape
    if shape is None or len(shape) != 1:
        raise RuntimeError(f"HDF5 obs column '{node.name}' is not one-dimensional.")
    return int(shape[0])


def read_obs_column(
    obs_group: Any,
    column: str,
    start: int | None = None,
    stop: int | None = None,
) -> np.ndarray:
    """Read one anndata 0.12 obs column without opening an AnnData object."""
    if column not in obs_group:
        raise KeyError(column)
    node = obs_group[column]
    encoding = _decode_attr(node.attrs.get("encoding-type"))
    if encoding in (None, "", "array", "string-array"):
        return read_str_dataset(node, start, stop)
    if encoding == "categorical":
        categories = read_str_dataset(node["categories"])
        codes = np.asarray(
            node["codes"][:] if start is None and stop is None else node["codes"][start:stop]
        ).astype(np.int64)
        if categories.ndim != 1:
            raise RuntimeError(f"HDF5 categorical column '{node.name}' has invalid categories.")
        # -1 is anndata/pandas' missing-category code.  Keep the historical
        # astype(str) spelling here; the public Sample reader rejects it.
        return np.where(codes >= 0, categories[np.clip(codes, 0, None)], "nan")
    raise RuntimeError(
        f"Unsupported encoding-type '{encoding}' for obs column '{column}'."
    )


def _read_dataset_values(node: Any, start: int | None = None, stop: int | None = None) -> np.ndarray:
    values = node[:] if start is None and stop is None else node[start:stop]
    values = np.asarray(values)
    if values.dtype.kind in "SU":
        return values.astype(str)
    if values.dtype.kind == "O":
        output = np.empty(values.shape, dtype=object)
        for index, value in enumerate(values):
            output[index] = value.decode("utf-8") if isinstance(value, bytes) else value
        return output
    return values


def read_obs_column_values(
    obs_group: Any,
    column: str,
    start: int | None = None,
    stop: int | None = None,
) -> np.ndarray:
    """Read one obs column while retaining numeric types and nullable values."""
    if column not in obs_group:
        raise KeyError(column)
    node = obs_group[column]
    encoding = _decode_attr(node.attrs.get("encoding-type"))
    if encoding in (None, "", "array", "string-array"):
        return _read_dataset_values(node, start, stop)
    if encoding == "categorical":
        categories = read_str_dataset(node["categories"])
        codes = np.asarray(
            node["codes"][:] if start is None and stop is None else node["codes"][start:stop]
        ).astype(np.int64)
        if categories.ndim != 1:
            raise RuntimeError(f"HDF5 categorical column '{node.name}' has invalid categories.")
        values = np.empty(codes.shape, dtype=object)
        values[:] = None
        valid = codes >= 0
        values[valid] = categories[np.clip(codes[valid], 0, None)]
        return values
    if encoding in {"nullable-integer", "nullable-boolean", "nullable-string"}:
        values = _read_dataset_values(node["values"], start, stop)
        mask = np.asarray(
            node["mask"][:] if start is None and stop is None else node["mask"][start:stop]
        ).astype(bool)
        if len(values) != len(mask):
            raise RuntimeError(f"HDF5 nullable column '{node.name}' has mismatched values/mask.")
        output = np.asarray(values, dtype=object)
        output[mask] = None
        return output
    raise RuntimeError(
        f"Unsupported encoding-type '{encoding}' for obs column '{column}'."
    )


def _obs_index_name(obs_group: Any) -> str:
    return str(_decode_attr(obs_group.attrs.get("_index", "_index")))


def read_h5ad_sample_ids(path: str | Path, sample_col: str = "Sample") -> list[str]:
    """Return ordered unique Sample IDs using only bounded HDF5 obs reads.

    Only ``obs`` is touched.  In particular, this function never opens the
    H5AD through anndata and never reads ``layers``/``X``.  The complete cell
    column is consumed in bounded slices while only the ordered unique sample
    IDs are retained.
    """
    h5ad_path = Path(path)
    if not h5ad_path.is_file() or h5ad_path.stat().st_size <= 0:
        raise ValueError(f"H5AD path is missing or empty: {h5ad_path}")

    with h5py.File(h5ad_path, "r") as handle:
        if "obs" not in handle:
            raise ValueError(f"H5AD lacks obs dataframe: {h5ad_path}")
        obs = handle["obs"]
        encoding = _decode_attr(obs.attrs.get("encoding-type"))
        if encoding != "dataframe":
            raise ValueError(
                f"H5AD obs encoding-type is {encoding!r}, expected 'dataframe': {h5ad_path}"
            )
        index_name = _obs_index_name(obs)
        if index_name not in obs:
            raise ValueError(f"H5AD lacks obs index {index_name!r}: {h5ad_path}")
        n_obs = _obs_column_length(obs[index_name])
        if n_obs <= 0:
            raise ValueError(f"H5AD obs index is empty: {h5ad_path}")
        if sample_col not in obs:
            raise ValueError(f"H5AD lacks obs[{sample_col!r}]: {h5ad_path}")
        sample_node = obs[sample_col]
        if _obs_column_length(sample_node) != n_obs:
            raise ValueError(f"H5AD Sample column length mismatches obs index: {h5ad_path}")

        ordered: list[str] = []
        seen: set[str] = set()
        for start in range(0, n_obs, READ_CHUNK_SIZE):
            stop = min(n_obs, start + READ_CHUNK_SIZE)
            values = read_obs_column(obs, sample_col, start, stop)
            for value in values:
                sample_id = str(value)
                if not sample_id.strip() or sample_id.casefold() == "nan":
                    raise ValueError(f"H5AD Sample contains a missing/blank value: {h5ad_path}")
                if sample_id not in seen:
                    seen.add(sample_id)
                    ordered.append(sample_id)

    if not ordered:
        raise ValueError(f"H5AD Sample has no nonempty identifiers: {h5ad_path}")
    return ordered


def file_md5(path: str | Path) -> str:
    """Hash a file in bounded chunks."""
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(READ_CHUNK_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def _sidecar_fields(path: Path, *, verify_content: bool = True) -> dict[str, str]:
    sidecar = Path(f"{path}.md5")
    if not sidecar.is_file() or sidecar.stat().st_size <= 0:
        raise ValueError(f"H5AD checksum sidecar is missing: {sidecar}")
    fields: dict[str, str] = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields[key] = value
    digest = fields.get("MD5", "").strip().lower()
    size = fields.get("SIZE", "").strip()
    if os.path.normpath(fields.get("PATH", "")) != os.path.normpath(str(path)):
        raise ValueError(
            f"H5AD checksum PATH mismatch: {sidecar} "
            f"(recorded={fields.get('PATH')!r}, expected={str(path)!r})"
        )
    if len(digest) != _MD5_LENGTH or any(ch not in "0123456789abcdef" for ch in digest):
        raise ValueError(f"H5AD checksum MD5 is malformed: {sidecar}")
    if not size.isdigit() or int(size) != path.stat().st_size:
        raise ValueError(f"H5AD checksum SIZE mismatch: {sidecar}")
    if verify_content and file_md5(path) != digest:
        raise ValueError(f"H5AD checksum MD5 mismatch: {path}")
    return {"MD5": digest, "SIZE": size, "PATH": str(path)}


def resolve_h5ad_path(
    input_root: str | Path,
    config_path: str | Path,
    dataset: str,
    view: str,
) -> Path:
    config = json.loads(Path(config_path).read_text())
    view_spec = config.get(dataset, {}).get("views", {}).get(view, {})
    output = view_spec.get("output_file_name") or view_spec.get("output_file")
    if not isinstance(output, str) or not output:
        raise ValueError(f"missing H5AD output contract for {dataset}/{view}")
    return Path(input_root) / dataset / "output" / output


def _selection_rows(selection: str | Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    seen_keys: set[tuple[str, str]] = set()
    seen_rows: set[tuple[str, ...]] = set()
    for line_number, line in enumerate(Path(selection).read_text().splitlines(), start=1):
        if not line:
            raise ValueError(f"selection contains a blank row at line {line_number}")
        parts = line.split("\t")
        if len(parts) not in (2, 3) or any(not part for part in parts):
            raise ValueError(f"selection row {line_number} must have two or three nonempty columns")
        row = tuple(parts)
        if row in seen_rows:
            raise ValueError(f"selection contains duplicate row at line {line_number}")
        seen_rows.add(row)
        key = (parts[0], parts[1])
        if key not in seen_keys:
            seen_keys.add(key)
            rows.append(key)
    if not rows:
        raise ValueError(f"selection is empty: {selection}")
    return rows


def _identity_payload(entries: Iterable[dict[str, Any]]) -> dict[str, Any]:
    return {"schema": IDENTITY_SCHEMA, "entries": list(entries)}


def build_source_identity(
    selection: str | Path,
    input_root: str | Path,
    config_path: str | Path,
    *,
    validated_sidecars: bool = False,
) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for dataset, view in _selection_rows(selection):
        path = resolve_h5ad_path(input_root, config_path, dataset, view)
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"H5AD path is missing or empty: {path}")
        # A source identity is never built from an unchecked artifact.  The
        # optimized caller may opt out of a second content hash only after it
        # has strictly validated every selected sidecar and source byte pair.
        sidecar = _sidecar_fields(path, verify_content=not validated_sidecars)
        sample_ids = read_h5ad_sample_ids(path)
        entries.append(
            {
                "dataset": dataset,
                "view": view,
                "path": str(path),
                "size": int(sidecar["SIZE"]),
                "md5": sidecar["MD5"],
                "sample_ids": sample_ids,
            }
        )
    return _identity_payload(entries)


def _validate_record(record: Any, key: str) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise ValueError(f"source identity entry is not an object: {key}")
    required = {"dataset", "view", "path", "size", "md5", "sample_ids"}
    if set(record) != required:
        raise ValueError(f"source identity entry has the wrong fields: {key}")
    dataset = record["dataset"]
    view = record["view"]
    path = record["path"]
    size = record["size"]
    digest = record["md5"]
    sample_ids = record["sample_ids"]
    if not all(isinstance(value, str) and value for value in (dataset, view, path)):
        raise ValueError(f"source identity entry has blank names/path: {key}")
    if not isinstance(size, int) or size <= 0:
        raise ValueError(f"source identity entry has invalid size: {key}")
    if not isinstance(digest, str) or len(digest) != _MD5_LENGTH or any(
        char not in "0123456789abcdefABCDEF" for char in digest
    ):
        raise ValueError(f"source identity entry has invalid MD5: {key}")
    if not isinstance(sample_ids, list) or not sample_ids or any(
        not isinstance(value, str) or not value.strip() for value in sample_ids
    ):
        raise ValueError(f"source identity entry has invalid sample IDs: {key}")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"source identity entry has duplicate sample IDs: {key}")
    return {
        "dataset": dataset,
        "view": view,
        "path": path,
        "size": size,
        "md5": digest.lower(),
        "sample_ids": sample_ids,
    }


def load_source_identity(path: str | Path, require_sidecar: bool = True) -> dict[tuple[str, str], dict[str, Any]]:
    identity_path = Path(path)
    if not identity_path.is_file() or identity_path.stat().st_size <= 0:
        raise ValueError(f"source identity manifest is missing or empty: {identity_path}")
    if require_sidecar:
        _sidecar_fields(identity_path)
    payload = json.loads(identity_path.read_text())
    if not isinstance(payload, dict) or payload.get("schema") != IDENTITY_SCHEMA:
        raise ValueError(f"source identity schema is invalid: {identity_path}")
    raw_entries = payload.get("entries")
    if not isinstance(raw_entries, list) or not raw_entries:
        raise ValueError(f"source identity entries are missing: {identity_path}")
    records: dict[tuple[str, str], dict[str, Any]] = {}
    for raw in raw_entries:
        record = _validate_record(raw, str(raw.get("dataset", "?")) if isinstance(raw, dict) else "?")
        key = (record["dataset"], record["view"])
        if key in records:
            raise ValueError(f"source identity contains duplicate row: {key[0]}/{key[1]}")
        records[key] = record
    return records


def verify_source_identity(
    identity_path: str | Path,
    selection: str | Path,
    input_root: str | Path,
    config_path: str | Path,
    *,
    validated_sidecars: bool = False,
) -> None:
    records = load_source_identity(identity_path, require_sidecar=True)
    selected = _selection_rows(selection)
    if list(records) != selected:
        raise ValueError("source identity rows do not exactly match the current selection")
    for dataset, view in selected:
        path = resolve_h5ad_path(input_root, config_path, dataset, view)
        record = records[(dataset, view)]
        if record["path"] != str(path):
            raise ValueError(f"source identity PATH mismatch for {dataset}/{view}")
        sidecar = _sidecar_fields(path, verify_content=not validated_sidecars)
        if (
            int(sidecar["SIZE"]) != record["size"]
            or sidecar["MD5"] != record["md5"]
        ):
            raise ValueError(f"source identity checksum mismatch for {dataset}/{view}")
        if read_h5ad_sample_ids(path) != record["sample_ids"]:
            raise ValueError(f"source identity Sample order mismatch for {dataset}/{view}")


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.build.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--identity", type=Path)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--validated-source-sidecars",
        action="store_true",
        help="trust caller-completed strict source sidecar/content validation",
    )
    args = parser.parse_args()
    if args.output is not None:
        payload = build_source_identity(
            args.selection,
            args.input_root,
            args.config,
            validated_sidecars=args.validated_source_sidecars,
        )
        _write_json_atomic(args.output, payload)
        print(f"source identity written: {args.output}")
    else:
        verify_source_identity(
            args.identity,
            args.selection,
            args.input_root,
            args.config,
            validated_sidecars=args.validated_source_sidecars,
        )
        print(f"source identity verified: {args.identity}")


if __name__ == "__main__":
    main()
