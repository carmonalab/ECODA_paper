#!/usr/bin/env python3
"""Create and validate the Alzheimer donor-by-assay H5AD derivative.

The SEA-AD source is a large H5AD.  This worker deliberately reads only HDF5
shape metadata and the ``obs`` vectors needed for the donor/assay contract;
values in ``X``, ``raw``, and ``layers`` are never loaded.  The derivative is
made by copying the source file and adding one obs column, so every source
matrix and observation is preserved byte-for-byte except for the new column.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import h5py
import numpy as np

# The worker is executed from an immutable source snapshot, not necessarily
# from the repository checkout containing the launcher.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.h5ad_source_identity import read_obs_column_values  # noqa: E402

DERIVED_COLUMN = "donor_id_assay"
DONOR_COLUMN = "donor_id"
ASSAY_COLUMN = "assay"
ASSAY_MAP: Mapping[str, str] = {
    "10x 3' v3": "10x3v3",
    "10x multiome": "10xmultiome",
}
SEX_LABELS = frozenset({"female", "male"})
# These fields are sample-level in the SEA-AD source.  They are optional for
# small synthetic fixtures, but if present they may not vary within a derived
# donor-by-assay sample.
PAIR_METADATA_COLUMNS = (
    "sex",
    "Sex",
    "Cognitive status",
    "disease",
    "tissue",
    "tissue_type",
    "PMI",
)
EXAMPLE_IDS = {
    ("H20.33.001", "10x 3' v3"): "H20.33.001_10x3v3",
    ("H20.33.001", "10x multiome"): "H20.33.001_10xmultiome",
}


class ContractError(ValueError):
    """Raised when a source or derivative violates the strict contract."""


def _decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.tobytes().decode("utf-8")
    return value


def _is_missing(value: Any) -> bool:
    value = _decode(value)
    if value is None:
        return True
    if isinstance(value, (float, np.floating)):
        return not math.isfinite(float(value))
    if isinstance(value, (complex, np.complexfloating)):
        return not math.isfinite(float(value.real)) or not math.isfinite(float(value.imag))
    # pd.NA and similar extension scalars must not be coerced into a boolean
    # expression.  Their canonical text spellings are unambiguously missing.
    try:
        if bool(np.isscalar(value) and np.isnat(value)):
            return True
    except (TypeError, ValueError):
        pass
    text = str(value)
    return text.strip() == "" or text.strip().lower() in {"nan", "nat", "none", "null", "<na>"}


def _required_text(value: Any, column: str, row: int) -> str:
    value = _decode(value)
    if _is_missing(value):
        raise ContractError(f"{column} is missing or blank at obs row {row}")
    text = str(value)
    if text != text.strip():
        raise ContractError(f"{column} has surrounding whitespace at obs row {row}")
    if any(char in text for char in "\r\n\t"):
        raise ContractError(f"{column} contains a record delimiter at obs row {row}")
    return text


def _assay_token(value: Any, row: int) -> tuple[str, str]:
    value = _decode(value)
    if _is_missing(value):
        raise ContractError(f"assay is missing or blank at obs row {row}")
    assay = str(value)
    # Do not trim or sanitize assay values: only the two exact source labels
    # are accepted by the scientific contract.
    token = ASSAY_MAP.get(assay)
    if token is None:
        raise ContractError(
            f"unsupported assay at obs row {row}: {assay!r}; "
            f"expected exactly {sorted(ASSAY_MAP)}"
        )
    return assay, token


def _text_key(value: Any) -> tuple[str, str]:
    value = _decode(value)
    if _is_missing(value):
        return ("missing", "")
    return ("value", str(value))


def _values_equal(left: Any, right: Any) -> bool:
    left = _decode(left)
    right = _decode(right)
    if _is_missing(left) or _is_missing(right):
        return _is_missing(left) and _is_missing(right)
    return str(left) == str(right)


def _sex_label(value: Any) -> str:
    """Return the canonical comparison key used for sample-level sex counts."""
    value = _decode(value)
    if _is_missing(value):
        return "<missing>"
    return str(value).casefold()


def _node_shape(node: Any, label: str, path: Path) -> tuple[int, ...]:
    shape = getattr(node, "shape", None)
    if shape is None:
        shape = node.attrs.get("shape") if hasattr(node, "attrs") else None
    try:
        result = tuple(int(value) for value in shape)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{path}: {label} has no valid persisted shape") from exc
    if len(result) != 2 or any(value <= 0 for value in result):
        raise ContractError(f"{path}: {label} has an empty or invalid shape {result!r}")
    return result


def _node_length(node: Any, label: str, path: Path) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "keys"):
        if "codes" in node:
            shape = node["codes"].shape
        elif "values" in node:
            shape = node["values"].shape
    try:
        result = tuple(int(value) for value in shape)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{path}: {label} has no valid one-dimensional shape") from exc
    if len(result) != 1:
        raise ContractError(f"{path}: {label} is not one-dimensional")
    return result[0]


def _index_name(obs: h5py.Group) -> str:
    value = _decode(obs.attrs.get("_index", "_index"))
    if not isinstance(value, str) or not value.strip() or value not in obs:
        raise ContractError("H5AD obs has no valid index column")
    return value


def _matrix_shapes(handle: h5py.File, path: Path) -> dict[str, tuple[int, ...]]:
    shapes: dict[str, tuple[int, ...]] = {}
    for name in ("X", "raw/X", "layers/counts"):
        if name in handle:
            shapes[name] = _node_shape(handle[name], name, path)
    if not shapes:
        raise ContractError(f"{path}: H5AD has no persisted expression/count matrix")
    return shapes


def _read_obs_column(obs: h5py.Group, column: str, n_obs: int, path: Path) -> np.ndarray:
    if column not in obs:
        raise ContractError(f"{path}: H5AD obs is missing required column {column!r}")
    node = obs[column]
    if _node_length(node, f"obs column {column!r}", path) != n_obs:
        raise ContractError(f"{path}: obs column {column!r} row count mismatch")
    try:
        values = np.asarray(read_obs_column_values(obs, column), dtype=object)
    except (KeyError, RuntimeError, ValueError, OSError) as exc:
        raise ContractError(f"{path}: could not read obs column {column!r}: {exc}") from exc
    if values.ndim != 1 or len(values) != n_obs:
        raise ContractError(f"{path}: obs column {column!r} is not a vector of n_obs rows")
    return values


def _validate_h5ad_layout(
    handle: h5py.File, path: Path, n_obs: int
) -> dict[str, tuple[int, ...]]:
    if "obs" not in handle or not isinstance(handle["obs"], h5py.Group):
        raise ContractError(f"{path}: H5AD obs group is missing")
    if "var" not in handle or not isinstance(handle["var"], h5py.Group):
        raise ContractError(f"{path}: H5AD var group is missing")
    shapes = _matrix_shapes(handle, path)
    wrong_rows = {
        name: shape
        for name, shape in shapes.items()
        if shape[0] != n_obs
    }
    if wrong_rows:
        raise ContractError(
            f"{path}: matrix row count does not match obs ({n_obs}): {wrong_rows!r}"
        )
    n_vars = next(iter(shapes.values()))[1]
    var = handle["var"]
    var_index = _decode(var.attrs.get("_index", "_index"))
    if not isinstance(var_index, str) or not var_index.strip() or var_index not in var:
        raise ContractError(f"{path}: H5AD var index is missing")
    if _node_length(var[var_index], "var index", path) != n_vars:
        raise ContractError(f"{path}: var index length does not match matrix columns")
    top_shape = handle.attrs.get("shape")
    if top_shape is not None:
        try:
            top_shape_tuple = tuple(int(value) for value in top_shape)
        except (TypeError, ValueError) as exc:
            raise ContractError(f"{path}: invalid H5AD top-level shape") from exc
        if top_shape_tuple != (n_obs, n_vars):
            raise ContractError(f"{path}: H5AD top-level shape does not match obs/matrix")
    return shapes


def _obs_columns(obs: h5py.Group) -> list[str]:
    return [str(name) for name in obs.keys()]


def _read_source_contract(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise ContractError(f"source H5AD is missing, empty, or a symlink: {path}")
    with h5py.File(path, "r") as handle:
        obs = handle.get("obs")
        if obs is None or not isinstance(obs, h5py.Group):
            raise ContractError(f"{path}: H5AD obs is missing or not a group")
        encoding = _decode(obs.attrs.get("encoding-type"))
        if encoding != "dataframe":
            raise ContractError(f"{path}: H5AD obs is not a dataframe")
        index_name = _index_name(obs)
        n_obs = _node_length(obs[index_name], "obs index", path)
        if n_obs <= 0:
            raise ContractError(f"{path}: H5AD obs is empty")
        matrix_shapes = _validate_h5ad_layout(handle, path, n_obs)
        index_values = _read_obs_column(obs, index_name, n_obs, path)
        index_text = [_required_text(value, "obs index", row) for row, value in enumerate(index_values)]
        if len(set(index_text)) != n_obs:
            raise ContractError(f"{path}: duplicate observation IDs are not permitted")
        if DERIVED_COLUMN in obs:
            raise ContractError(f"{path}: raw source already contains {DERIVED_COLUMN!r}")
        donor_values = _read_obs_column(obs, DONOR_COLUMN, n_obs, path)
        assay_values = _read_obs_column(obs, ASSAY_COLUMN, n_obs, path)
        donors: list[str] = []
        assays: list[str] = []
        tokens: list[str] = []
        derived: list[str] = []
        pair_by_derived: dict[str, tuple[str, str]] = {}
        for row, (donor_value, assay_value) in enumerate(zip(donor_values, assay_values)):
            donor = _required_text(donor_value, DONOR_COLUMN, row)
            assay, token = _assay_token(assay_value, row)
            identifier = f"{donor}_{token}"
            pair = (donor, assay)
            previous = pair_by_derived.get(identifier)
            if previous is not None and previous != pair:
                raise ContractError(
                    f"derived ID collision for {identifier!r}: {previous!r} versus {pair!r}"
                )
            pair_by_derived[identifier] = pair
            donors.append(donor)
            assays.append(assay)
            tokens.append(token)
            derived.append(identifier)

        pair_metadata: dict[str, np.ndarray] = {}
        for column in PAIR_METADATA_COLUMNS:
            if column in obs:
                pair_metadata[column] = _read_obs_column(obs, column, n_obs, path)
        for column, values in pair_metadata.items():
            first_seen: dict[str, tuple[str, str]] = {}
            for row, (identifier, value) in enumerate(zip(derived, values)):
                key = _text_key(value)
                previous = first_seen.get(identifier)
                if previous is None:
                    first_seen[identifier] = key
                elif previous != key:
                    raise ContractError(
                        f"mixed {column!r} metadata within derived sample {identifier!r} "
                        f"(rows include {row})"
                    )

        sex_column = next(
            (column for column in ("sex", "Sex") if column in pair_metadata),
            None,
        )
        sex_counts: dict[str, int] = {}
        if sex_column is not None:
            sex_by_derived: dict[str, str] = {}
            for identifier, value in zip(derived, pair_metadata[sex_column]):
                if identifier not in sex_by_derived:
                    sex_by_derived[identifier] = _sex_label(value)
            for label in sex_by_derived.values():
                sex_counts[label] = sex_counts.get(label, 0) + 1
            sex_counts = dict(sorted(sex_counts.items()))

        assay_counts = {
            token: len(
                {
                    identifier
                    for identifier, identifier_token in zip(derived, tokens)
                    if identifier_token == token
                }
            )
            for token in sorted(set(tokens))
        }
        return {
            "n_obs": n_obs,
            "index_name": index_name,
            "index_values": index_values,
            "obs_columns": _obs_columns(obs),
            "donor_values": donor_values,
            "assay_values": assay_values,
            "donors": donors,
            "assays": assays,
            "tokens": tokens,
            "derived": derived,
            "pair_metadata": pair_metadata,
            "n_donors": len(set(donors)),
            "n_samples": len(set(derived)),
            "assay_counts": assay_counts,
            "sex_column": sex_column,
            "sex_counts": sex_counts,
            "matrix_shapes": matrix_shapes,
        }


def _validate_expectations(
    contract: Mapping[str, Any],
    *,
    expected_samples: int | None,
    expected_donors: int | None,
    expected_assay_counts: Mapping[str, int] | None,
    expected_sex_counts: Mapping[str, int] | None = None,
    require_example_ids: bool,
) -> None:
    if expected_samples is not None and contract["n_samples"] != expected_samples:
        raise ContractError(
            f"expected {expected_samples} unique donor-by-assay samples, "
            f"observed {contract['n_samples']}"
        )
    if expected_donors is not None and contract["n_donors"] != expected_donors:
        raise ContractError(
            f"expected {expected_donors} donors, observed {contract['n_donors']}"
        )
    if expected_assay_counts is not None:
        observed = dict(contract["assay_counts"])
        if observed != dict(expected_assay_counts):
            raise ContractError(
                f"assay-token sample counts mismatch: expected "
                f"{dict(expected_assay_counts)!r}, observed {observed!r}"
            )
    if expected_sex_counts is not None:
        expected = _normalise_expected_sex_mapping(expected_sex_counts)
        observed = dict(contract.get("sex_counts", {}))
        if observed != expected:
            raise ContractError(
                f"sex sample counts mismatch: expected {expected!r}, observed {observed!r}"
            )
    if require_example_ids:
        pairs = set(zip(contract["donors"], contract["assays"]))
        missing = [pair for pair in EXAMPLE_IDS if pair not in pairs]
        if missing:
            raise ContractError(f"required deterministic example pair(s) missing: {missing!r}")
        for pair, expected in EXAMPLE_IDS.items():
            actual = f"{pair[0]}_{ASSAY_MAP[pair[1]]}"
            if actual != expected:
                raise ContractError(f"deterministic example ID mismatch for {pair!r}")


def _compare_source_obs(input_path: Path, output_path: Path, contract: Mapping[str, Any]) -> None:
    with h5py.File(input_path, "r") as source, h5py.File(output_path, "r") as output:
        source_obs = source["obs"]
        output_obs = output.get("obs")
        if output_obs is None or not isinstance(output_obs, h5py.Group):
            raise ContractError(f"{output_path}: output obs is missing or not a group")
        source_columns = set(contract["obs_columns"])
        output_columns = set(_obs_columns(output_obs))
        expected_columns = source_columns | {DERIVED_COLUMN}
        if output_columns != expected_columns:
            missing = sorted(expected_columns - output_columns)
            extra = sorted(output_columns - expected_columns)
            raise ContractError(
                f"{output_path}: obs columns changed; missing={missing!r}, extra={extra!r}"
            )
        n_obs = int(contract["n_obs"])
        for column in contract["obs_columns"]:
            source_values = _read_obs_column(source_obs, column, n_obs, input_path)
            output_values = _read_obs_column(output_obs, column, n_obs, output_path)
            if len(source_values) != len(output_values) or any(
                not _values_equal(left, right)
                for left, right in zip(source_values, output_values)
            ):
                raise ContractError(f"{output_path}: source obs column changed: {column!r}")
        output_index_name = _index_name(output_obs)
        if output_index_name != contract["index_name"]:
            raise ContractError(f"{output_path}: obs index name changed")
        derived_values = _read_obs_column(output_obs, DERIVED_COLUMN, n_obs, output_path)
        expected = np.asarray(contract["derived"], dtype=object)
        if any(not _values_equal(left, right) for left, right in zip(derived_values, expected)):
            raise ContractError(f"{output_path}: {DERIVED_COLUMN} values do not match source metadata")
        if len(set(str(value) for value in derived_values)) != contract["n_samples"]:
            raise ContractError(f"{output_path}: derived sample IDs are unexpectedly duplicated or changed")

        source_shapes = _validate_h5ad_layout(source, input_path, n_obs)
        output_shapes = _validate_h5ad_layout(output, output_path, n_obs)
        if source_shapes != dict(contract["matrix_shapes"]):
            raise ContractError(f"{input_path}: source matrix layout changed during validation")
        if output_shapes != source_shapes:
            raise ContractError(
                f"{output_path}: matrix shapes changed; expected {source_shapes!r}, observed {output_shapes!r}"
            )
        source_shape = source.attrs.get("shape")
        output_shape = output.attrs.get("shape")
        if source_shape is not None:
            try:
                source_shape_tuple = tuple(int(value) for value in source_shape)
                output_shape_tuple = tuple(int(value) for value in output_shape)
            except (TypeError, ValueError) as exc:
                raise ContractError(f"{output_path}: invalid H5AD top-level shape") from exc
            if source_shape_tuple != output_shape_tuple or not source_shape_tuple or source_shape_tuple[0] != n_obs:
                raise ContractError(f"{output_path}: H5AD top-level shape changed or has wrong row count")


def _validate_output_contract(
    input_path: Path,
    output_path: Path,
    source_contract: Mapping[str, Any],
    *,
    expected_samples: int | None,
    expected_donors: int | None,
    expected_assay_counts: Mapping[str, int] | None,
    expected_sex_counts: Mapping[str, int] | None = None,
    require_example_ids: bool,
) -> None:
    if not output_path.is_file() or output_path.is_symlink() or output_path.stat().st_size <= 0:
        raise ContractError(f"output H5AD is missing, empty, or a symlink: {output_path}")
    with h5py.File(output_path, "r") as handle:
        obs = handle.get("obs")
        if obs is None or not isinstance(obs, h5py.Group):
            raise ContractError(f"{output_path}: output obs is missing")
        encoding = _decode(obs.attrs.get("encoding-type"))
        if encoding != "dataframe":
            raise ContractError(f"{output_path}: output obs is not a dataframe")
        n_obs = int(source_contract["n_obs"])
        _validate_h5ad_layout(handle, output_path, n_obs)
        index_name = _index_name(obs)
        if _node_length(obs[index_name], "output obs index", output_path) != n_obs:
            raise ContractError(f"{output_path}: output row count changed")
        _read_obs_column(obs, DERIVED_COLUMN, n_obs, output_path)
    _compare_source_obs(input_path, output_path, source_contract)
    _validate_expectations(
        source_contract,
        expected_samples=expected_samples,
        expected_donors=expected_donors,
        expected_assay_counts=expected_assay_counts,
        expected_sex_counts=expected_sex_counts,
        require_example_ids=require_example_ids,
    )


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_checksum(path: Path) -> None:
    sidecar = Path(f"{path}.md5")
    if sidecar.is_symlink() or (sidecar.exists() and not sidecar.is_file()):
        raise ContractError(f"checksum sidecar is not a regular file: {sidecar}")
    digest = _md5(path)
    temporary = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(
            f"MD5={digest}\nSIZE={path.stat().st_size}\nPATH={path}\n",
            encoding="utf-8",
        )
        os.replace(temporary, sidecar)
    finally:
        if temporary.exists():
            temporary.unlink()


def _validate_checksum(path: Path) -> None:
    sidecar = Path(f"{path}.md5")
    if sidecar.is_symlink() or not sidecar.is_file():
        raise ContractError(f"checksum sidecar is missing or unsafe: {sidecar}")
    lines = sidecar.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3 or not lines[0].startswith("MD5=") or not lines[1].startswith("SIZE=") or not lines[2].startswith("PATH="):
        raise ContractError(f"checksum sidecar has invalid format: {sidecar}")
    digest = lines[0][4:]
    size_text = lines[1][5:]
    recorded_path = lines[2][5:]
    if len(digest) != 32 or any(char not in "0123456789abcdef" for char in digest):
        raise ContractError(f"checksum sidecar has invalid MD5: {sidecar}")
    try:
        size = int(size_text)
    except ValueError as exc:
        raise ContractError(f"checksum sidecar has invalid SIZE: {sidecar}") from exc
    if size <= 0 or size != path.stat().st_size:
        raise ContractError(f"checksum sidecar SIZE mismatch: {sidecar}")
    if os.path.normpath(recorded_path) != os.path.normpath(str(path)):
        raise ContractError(f"checksum sidecar PATH mismatch: {sidecar}")
    if _md5(path) != digest:
        raise ContractError(f"checksum sidecar MD5 mismatch: {path}")


def _normalise_expected_assay_counts(value: str | None) -> dict[str, int] | None:
    if value is None or value == "":
        return None
    result: dict[str, int] = {}
    for item in value.split(","):
        if "=" not in item:
            raise ContractError(f"invalid --expected-assay-counts item: {item!r}")
        token, count_text = item.split("=", 1)
        token = token.strip()
        if token not in set(ASSAY_MAP.values()) or not count_text.isdigit():
            raise ContractError(f"invalid --expected-assay-counts item: {item!r}")
        result[token] = int(count_text)
    if set(result) != set(ASSAY_MAP.values()):
        raise ContractError("--expected-assay-counts must name both accepted assay tokens")
    return result


def _normalise_expected_sex_mapping(value: Mapping[str, int]) -> dict[str, int]:
    result: dict[str, int] = {}
    for label, count in value.items():
        canonical = str(_decode(label)).strip().casefold()
        if canonical not in SEX_LABELS:
            raise ContractError(
                f"invalid expected sex label {label!r}; expected female and male"
            )
        if canonical in result:
            raise ContractError(f"duplicate expected sex label: {label!r}")
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 0:
            raise ContractError(f"invalid expected sex count for {label!r}: {count!r}")
        result[canonical] = int(count)
    if set(result) != SEX_LABELS:
        raise ContractError("--expected-sex-counts must name both female and male")
    return dict(sorted(result.items()))


def _normalise_expected_sex_counts(value: str | None) -> dict[str, int] | None:
    if value is None or value == "":
        return None
    result: dict[str, int] = {}
    for item in value.split(","):
        if "=" not in item:
            raise ContractError(f"invalid --expected-sex-counts item: {item!r}")
        label, count_text = item.split("=", 1)
        label = label.strip().casefold()
        if (
            label not in SEX_LABELS
            or not count_text.isdigit()
            or label in result
        ):
            raise ContractError(f"invalid --expected-sex-counts item: {item!r}")
        result[label] = int(count_text)
    return _normalise_expected_sex_mapping(result)

def validate_derivative(
    input_path: str | Path,
    output_path: str | Path,
    *,
    expected_samples: int | None = None,
    expected_donors: int | None = None,
    expected_assay_counts: Mapping[str, int] | None = None,
    expected_sex_counts: Mapping[str, int] | None = None,
    require_example_ids: bool = False,
    require_checksum: bool = True,
) -> dict[str, Any]:
    """Validate source/output metadata, schema, sample contract, and checksum."""
    source = Path(input_path)
    output = Path(output_path)
    if source.resolve() == output.resolve():
        raise ContractError("source and derivative paths must differ")
    contract = _read_source_contract(source)
    _validate_output_contract(
        source,
        output,
        contract,
        expected_samples=expected_samples,
        expected_donors=expected_donors,
        expected_assay_counts=expected_assay_counts,
        expected_sex_counts=expected_sex_counts,
        require_example_ids=require_example_ids,
    )
    if require_checksum:
        _validate_checksum(output)
    return contract


def create_derivative(
    input_path: str | Path,
    output_path: str | Path,
    *,
    expected_samples: int | None = None,
    expected_donors: int | None = None,
    expected_assay_counts: Mapping[str, int] | None = None,
    expected_sex_counts: Mapping[str, int] | None = None,
    require_example_ids: bool = False,
    force: bool = False,
) -> bool:
    """Create the derivative atomically; return ``False`` for validated NOOP."""
    source = Path(input_path)
    output = Path(output_path)
    if source.resolve() == output.resolve():
        raise ContractError("source and derivative paths must differ")
    source_contract = _read_source_contract(source)
    _validate_expectations(
        source_contract,
        expected_samples=expected_samples,
        expected_donors=expected_donors,
        expected_assay_counts=expected_assay_counts,
        expected_sex_counts=expected_sex_counts,
        require_example_ids=require_example_ids,
    )
    if output.exists() or output.is_symlink():
        if not force:
            validate_derivative(
                source,
                output,
                expected_samples=expected_samples,
                expected_donors=expected_donors,
                expected_assay_counts=expected_assay_counts,
                expected_sex_counts=expected_sex_counts,
                require_example_ids=require_example_ids,
                require_checksum=True,
            )
            return False
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_symlink():
        raise ContractError(f"output path is a symlink: {output}")
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise ContractError(f"temporary output path already exists: {temporary}")
    try:
        # copyfile streams bytes; it never constructs an AnnData object or
        # reads X/raw/layers values into Python memory.
        shutil.copyfile(source, temporary)
        with h5py.File(temporary, "r+") as handle:
            obs = handle.get("obs")
            if obs is None or not isinstance(obs, h5py.Group):
                raise ContractError(f"{temporary}: copied H5AD obs is missing")
            if DERIVED_COLUMN in obs:
                raise ContractError(f"{temporary}: copied H5AD already has {DERIVED_COLUMN!r}")
            string_dtype = h5py.string_dtype(encoding="utf-8")
            node = obs.create_dataset(
                DERIVED_COLUMN,
                shape=(source_contract["n_obs"],),
                dtype=string_dtype,
            )
            node[:] = np.asarray(source_contract["derived"], dtype=object)
            node.attrs["encoding-type"] = "string-array"
            node.attrs["encoding-version"] = "0.2.0"
            if "column-order" in obs.attrs:
                column_order = obs.attrs["column-order"]
                if isinstance(column_order, (str, bytes, np.str_, np.bytes_)):
                    values = [_decode(column_order)]
                else:
                    values = [_decode(value) for value in np.asarray(column_order).tolist()]
                if DERIVED_COLUMN not in values:
                    values.append(DERIVED_COLUMN)
                obs.attrs["column-order"] = np.asarray(
                    values, dtype=h5py.string_dtype(encoding="utf-8")
                )
        _validate_output_contract(
            source,
            temporary,
            source_contract,
            expected_samples=expected_samples,
            expected_donors=expected_donors,
            expected_assay_counts=expected_assay_counts,
            expected_sex_counts=expected_sex_counts,
            require_example_ids=require_example_ids,
        )
        os.replace(temporary, output)
        _write_checksum(output)
        validate_derivative(
            source,
            output,
            expected_samples=expected_samples,
            expected_donors=expected_donors,
            expected_assay_counts=expected_assay_counts,
            expected_sex_counts=expected_sex_counts,
            require_example_ids=require_example_ids,
            require_checksum=True,
        )
        return True
    finally:
        if temporary.exists():
            temporary.unlink()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-file", "--input_file", "--input", dest="input_file", required=True)
    parser.add_argument("--output-file", "--output_file", "--output", dest="output_file", required=True)
    parser.add_argument("--expected-samples", type=int, default=None)
    parser.add_argument("--expected-donors", type=int, default=None)
    parser.add_argument(
        "--expected-assay-counts",
        default=None,
        help="comma-separated assay-token counts, e.g. 10x3v3=83,10xmultiome=21",
    )

    parser.add_argument(
        "--expected-sex-counts",
        default=None,
        help="comma-separated sex sample counts, e.g. female=59,male=45",
    )
    parser.add_argument("--require-example-ids", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        expected_counts = _normalise_expected_assay_counts(args.expected_assay_counts)
        expected_sex_counts = _normalise_expected_sex_counts(args.expected_sex_counts)
        if args.validate_only:
            contract = validate_derivative(
                args.input_file,
                args.output_file,
                expected_samples=args.expected_samples,
                expected_donors=args.expected_donors,
                expected_assay_counts=expected_counts,
                expected_sex_counts=expected_sex_counts,
                require_example_ids=args.require_example_ids,
                require_checksum=True,
            )
            print(
                f"ALZHEIMER_DONOR_ASSAY_VALIDATED=1 cells={contract['n_obs']} "
                f"samples={contract['n_samples']}",
                flush=True,
            )
            return 0
        created = create_derivative(
            args.input_file,
            args.output_file,
            expected_samples=args.expected_samples,
            expected_donors=args.expected_donors,
            expected_assay_counts=expected_counts,
            expected_sex_counts=expected_sex_counts,
            require_example_ids=args.require_example_ids,
            force=args.force,
        )
        if created:
            print(
                f"ALZHEIMER_DONOR_ASSAY_CREATED=1 output={args.output_file}",
                flush=True,
            )
        else:
            print(
                f"NOOP_VALIDATED=1 output={args.output_file}",
                flush=True,
            )
        return 0
    except (ContractError, OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
