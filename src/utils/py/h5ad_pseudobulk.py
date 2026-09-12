"""Bounded H5AD/CSR reducers for sample and cell-type pseudobulk.

The public sample reducer keeps the historical genes-by-samples ``int64``
contract.  The cell-type reducer has a separate, run-owned HDF5/CSR store
boundary: metadata is discovered first, one raw count pass appends sparse cell
contributions, and a sparse on-disk aggregate is maintained with checked
additions.  No path in this module materializes a dense per-chunk grouped
``U x n_genes`` array.
"""

from __future__ import annotations

import errno
import json
import os
import socket
import time
import uuid
from pathlib import Path
from typing import Any, Iterable, Mapping

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:  # import_from_path exposes this directory as a top-level module path
    from h5ad_source_identity import read_obs_column_values, read_str_dataset
except ImportError:  # package imports used by focused tests
    from .h5ad_source_identity import read_obs_column_values, read_str_dataset

try:  # import_from_path exposes this directory as a top-level module path
    from batch_contract import (
        BatchContractError,
        RESERVED_OBS_NAME,
        normalize_batch_keys,
        serialize_batch_metadata,
        validate_batch_metadata,
    )
except ImportError:  # package imports used by focused tests
    from .batch_contract import (
        BatchContractError,
        RESERVED_OBS_NAME,
        normalize_batch_keys,
        serialize_batch_metadata,
        validate_batch_metadata,
    )


DEFAULT_CHUNK_SIZE = 4096
_INT64_MAX = np.iinfo(np.int64).max
_INT64_EXCLUSIVE = 1 << 63
_INT_MAX = np.iinfo(np.int32).max
_STORE_SCHEMA = 1
_STORE_STAGE = "pseudobulk_ct"
_STORE_STATE_WRITING = "writing"
_STORE_STATE_READY = "ready"
__all__ = [
    "DEFAULT_CHUNK_SIZE",
    "aggregate_h5ad_counts_by_sample",
    "audit_h5ad_ct_group_store",
    "prepare_h5ad_ct_group_store",
    "read_h5ad_ct_group_store",
    "read_h5ad_sample_metadata",
    "validate_h5ad_corrected_batch_metadata",
]



def _decode(value: Any) -> Any:
    return value.decode("utf-8") if isinstance(value, bytes) else value


def _json_safe(value: Any) -> Any:
    """Convert metadata/manifest values to JSON-safe scalar containers."""
    value = _decode(value)
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if value is None:
        return None
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _persisted_shape(node) -> tuple[int, ...]:
    shape = getattr(node, "shape", None)
    if shape is None:
        shape = node.attrs.get("shape")
    if shape is None:
        return ()
    try:
        return tuple(int(value) for value in shape)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"HDF5 node {node.name} has an invalid shape") from exc


def _node_length(node) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and "codes" in node:
        shape = node["codes"].shape
    if shape is None or len(shape) != 1:
        raise ValueError(f"HDF5 node {node.name} is not a one-dimensional vector")
    return int(shape[0])


def _obs_column_length(node) -> int:
    """Return an obs column's vector length, including encoded columns."""
    shape = getattr(node, "shape", None)
    if shape is None and "codes" in node:
        shape = node["codes"].shape
    if shape is None and "values" in node:
        shape = node["values"].shape
    if shape is None or len(shape) != 1:
        raise ValueError(f"HDF5 obs column {node.name} is not one-dimensional")
    return int(shape[0])


def _metadata_scalar(value: Any) -> Any:
    value = _decode(value)
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return None
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        return None
    return value


def _normalise_identifier(value: Any) -> str | None:
    value = _metadata_scalar(value)
    if value is None:
        return None
    text = str(value)
    if not text.strip() or text.casefold() == "nan":
        return None
    return text


def _requested_columns(sample_col: str, metadata_columns: Iterable[str] | None) -> list[str]:
    if isinstance(metadata_columns, str):
        metadata_columns = [metadata_columns]
    columns = [sample_col, *(str(value) for value in (metadata_columns or ()))]
    columns = [value for value in columns if value]
    return list(dict.fromkeys(columns))


def _validate_chunk_size(chunk_size: int) -> int:
    try:
        value = int(chunk_size)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"chunk_size must be a positive integer, got {chunk_size!r}") from exc
    if value <= 0:
        raise ValueError(f"chunk_size must be a positive integer, got {chunk_size!r}")
    return value


def _validate_max_value(max_value: int | None) -> int | None:
    if max_value is None:
        return None
    if isinstance(max_value, (bool, np.bool_)):
        raise ValueError("max_value must be a nonnegative integer")
    if isinstance(max_value, (float, np.floating)):
        if not np.isfinite(max_value) or float(max_value) != float(int(max_value)):
            raise ValueError(f"max_value must be a nonnegative integer, got {max_value!r}")
    try:
        value = int(max_value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"max_value must be a nonnegative integer, got {max_value!r}") from exc
    if value < 0 or value > _INT64_MAX:
        raise ValueError(f"max_value must be between 0 and {_INT64_MAX}")
    return value


def _open_obs_and_counts(handle, artifact: Path):
    x = handle.get("X")
    x_shape = _persisted_shape(x) if x is not None else ()
    if len(x_shape) != 2 or any(value <= 0 for value in x_shape):
        raise ValueError(f"H5AD X is missing or has an invalid shape: {artifact}")
    n_obs, n_vars = x_shape

    obs = handle.get("obs")
    if obs is None or _decode(obs.attrs.get("encoding-type")) != "dataframe":
        raise ValueError(f"H5AD obs is not a dataframe: {artifact}")
    index_name = str(_decode(obs.attrs.get("_index", "_index")))
    if index_name not in obs or _node_length(obs[index_name]) != n_obs:
        raise ValueError(f"H5AD obs index is missing or has the wrong length: {artifact}")

    layers = handle.get("layers")
    counts = layers.get("counts") if layers is not None else None
    counts_shape = _persisted_shape(counts) if counts is not None else ()
    if counts is None or counts_shape != (n_obs, n_vars):
        raise ValueError(
            f"H5AD layers['counts'] is missing or has shape {counts_shape}; "
            f"expected {(n_obs, n_vars)}: {artifact}"
        )
    if _decode(counts.attrs.get("encoding-type")) != "csr_matrix":
        raise ValueError(f"H5AD counts is not CSR: {artifact}")
    if not all(name in counts for name in ("data", "indices", "indptr")):
        raise ValueError(f"H5AD counts is incomplete: {artifact}")
    indptr_length = _node_length(counts["indptr"])
    data_length = _node_length(counts["data"])
    indices_length = _node_length(counts["indices"])
    if indptr_length != n_obs + 1:
        raise ValueError(f"H5AD counts indptr has the wrong length: {artifact}")
    first_indptr = _checked_integer_vector(
        np.asarray(counts["indptr"][:1]), "indptr"
    )
    final_indptr = _checked_integer_vector(
        np.asarray(counts["indptr"][-1:]), "indptr"
    )
    if int(first_indptr[0]) != 0:
        raise ValueError(f"H5AD counts indptr does not start at zero: {artifact}")
    if int(final_indptr[0]) != data_length or int(final_indptr[0]) != indices_length:
        raise ValueError(
            f"H5AD counts indptr does not match CSR data/index lengths: {artifact}"
        )

    var = handle.get("var")
    if var is None:
        raise ValueError(f"H5AD var is missing: {artifact}")
    var_index_name = str(_decode(var.attrs.get("_index", "_index")))
    if var_index_name not in var:
        raise ValueError(f"H5AD var index is missing: {artifact}")
    gene_names = read_str_dataset(var[var_index_name])
    if len(gene_names) != n_vars or any(not str(value).strip() for value in gene_names):
        raise ValueError(f"H5AD var index is invalid: {artifact}")
    if len(set(str(value) for value in gene_names)) != n_vars:
        raise ValueError(f"H5AD var index contains duplicate genes: {artifact}")
    return obs, counts, n_obs, n_vars, np.asarray(gene_names, dtype=str), index_name




def _validate_obs_columns(obs, columns: list[str], n_obs: int) -> None:
    available = set(str(name) for name in obs.keys())
    missing = sorted(set(columns) - available)
    if missing:
        raise ValueError(f"H5AD is missing requested obs columns {missing}")
    for column in columns:
        if _obs_column_length(obs[column]) != n_obs:
            raise ValueError(f"H5AD obs column {column!r} has the wrong length")


def _read_obs_chunks(obs, columns: list[str], n_obs: int, chunk_size: int):
    _validate_obs_columns(obs, columns, n_obs)
    for start in range(0, n_obs, chunk_size):
        stop = min(n_obs, start + chunk_size)
        yield (
            start,
            stop,
            {
                column: read_obs_column_values(obs, column, start, stop)
                for column in columns
            },
        )




def _collect_sample_metadata(
    obs, sample_col: str, metadata_columns: list[str], n_obs: int, chunk_size: int
):
    sample_ids: list[str] = []
    sample_to_index: dict[str, int] = {}
    values = {column: [] for column in metadata_columns}
    for _start, _stop, chunk in _read_obs_chunks(obs, metadata_columns, n_obs, chunk_size):
        sample_values = chunk[sample_col]
        for offset, raw_sample in enumerate(sample_values):
            sample = _normalise_identifier(raw_sample)
            if sample is None:
                raise ValueError(f"H5AD obs[{sample_col!r}] contains a missing/blank sample ID")
            if sample not in sample_to_index:
                sample_to_index[sample] = len(sample_ids)
                sample_ids.append(sample)
                for column in metadata_columns:
                    value = _metadata_scalar(chunk[column][offset])
                    if column == sample_col:
                        value = sample
                    values[column].append(value)
    if not sample_ids:
        raise ValueError("H5AD contains no non-empty sample IDs")
    metadata = pd.DataFrame(
        values,
        index=pd.Index(sample_ids, name=sample_col),
    )
    return sample_ids, sample_to_index, metadata


def validate_h5ad_corrected_batch_metadata(
    path: str | Path,
    batch_keys: str | Iterable[str],
    sample_col: str = "Sample",
    biological_column: str | None = None,
    method_id: str = "Pseudobulk",
    model_id: str = "pseudobulk_composite_v1",
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    near_unique_fraction: float = 0.50,
) -> dict[str, Any]:
    """Validate corrected batch metadata across every selected cell.

    This is the metadata-only boundary for corrected callers.  It reads the
    complete ``obs`` vectors needed by the batch contract in bounded chunks,
    before any sample or cell-type reducer can retain first-observation
    values.  The H5AD expression/count nodes are never opened or materialized.
    """
    keys = normalize_batch_keys(
        batch_keys,
        sample_column=sample_col,
        biological_column=biological_column,
    )
    chunk_size = _validate_chunk_size(chunk_size)
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")

    columns = list(dict.fromkeys([sample_col, *keys]))
    if biological_column is not None:
        columns.append(biological_column)
        columns = list(dict.fromkeys(columns))
    metadata_values: dict[str, list[Any]] = {column: [] for column in columns}

    with h5py.File(artifact, "r") as handle:
        obs = handle.get("obs")
        if obs is None or _decode(obs.attrs.get("encoding-type")) != "dataframe":
            raise ValueError(f"H5AD obs is not a dataframe: {artifact}")
        if RESERVED_OBS_NAME in obs:
            raise BatchContractError(
                f"metadata already contains reserved temporary column {RESERVED_OBS_NAME!r}"
            )
        index_name = str(_decode(obs.attrs.get("_index", "_index")))
        if index_name not in obs:
            raise ValueError(f"H5AD obs index is missing: {artifact}")
        n_obs = _obs_column_length(obs[index_name])
        if n_obs <= 0:
            raise ValueError(f"H5AD obs index is empty: {artifact}")

        for _start, _stop, chunk in _read_obs_chunks(obs, columns, n_obs, chunk_size):
            for column in columns:
                metadata_values[column].extend(list(chunk[column]))

    validation = validate_batch_metadata(
        metadata_values,
        batch_keys,
        sample_column=sample_col,
        biological_column=biological_column,
        near_unique_fraction=near_unique_fraction,
    )
    serialized = serialize_batch_metadata(
        validation,
        method_id=method_id,
        model_id=model_id,
    )
    return _json_safe(serialized)


def read_h5ad_sample_metadata(
    path: str | Path,
    sample_col: str = "Sample",
    metadata_columns: Iterable[str] | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> pd.DataFrame:
    """Return first-observation metadata for each sample without reading counts."""
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    sample_col = str(sample_col)
    if not sample_col:
        raise ValueError("sample_col must be non-empty")
    chunk_size = _validate_chunk_size(chunk_size)
    columns = _requested_columns(sample_col, metadata_columns)
    with h5py.File(artifact, "r") as handle:
        x = handle.get("X")
        x_shape = _persisted_shape(x) if x is not None else ()
        if len(x_shape) != 2 or any(value <= 0 for value in x_shape):
            raise ValueError(f"H5AD X is missing or has an invalid shape: {artifact}")
        obs = handle.get("obs")
        if obs is None or _decode(obs.attrs.get("encoding-type")) != "dataframe":
            raise ValueError(f"H5AD obs is not a dataframe: {artifact}")
        index_name = str(_decode(obs.attrs.get("_index", "_index")))
        if index_name not in obs or _node_length(obs[index_name]) != x_shape[0]:
            raise ValueError(f"H5AD obs index is missing or has the wrong length: {artifact}")
        _sample_ids, _sample_to_index, metadata = _collect_sample_metadata(
            obs, sample_col, columns, x_shape[0], chunk_size
        )
    return metadata


def _checked_integer_vector(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1 or array.dtype.kind not in "iu":
        raise ValueError(f"H5AD counts {name} must be a one-dimensional integer vector")
    if array.size:
        if array.dtype.kind == "i" and int(np.min(array)) < 0:
            raise ValueError(f"H5AD counts {name} contains a negative value")
        if int(np.max(array)) > _INT64_MAX:
            raise ValueError(f"H5AD counts {name} exceeds int64 capacity")
    return np.asarray(array, dtype=np.int64)


def _checked_count_values(values: Any, max_value: int | None = None) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1:
        raise ValueError("H5AD counts data must be a one-dimensional vector")
    kind = array.dtype.kind
    if kind == "f":
        if not np.isfinite(array).all():
            raise ValueError("H5AD counts data must be finite, nonnegative, integer-valued")
        if np.any(array < 0) or np.any(array != np.floor(array)):
            raise ValueError("H5AD counts data must be finite, nonnegative, integer-valued")
        # Comparing against 2**63 rather than float(INT64_MAX) avoids the
        # representational rounding of INT64_MAX to 2**63.
        if np.any(array >= float(_INT64_EXCLUSIVE)):
            raise ValueError("H5AD counts exceed int64 aggregation capacity")
        array = np.rint(array).astype(np.int64)
        if max_value is not None and np.any(array > max_value):
            raise ValueError(f"H5AD counts exceed max_value={max_value}")
    elif kind in "iu":
        if kind == "i" and array.size and int(np.min(array)) < 0:
            raise ValueError("H5AD counts data contains a negative signed-integer value")
        if array.size and int(np.max(array)) > _INT64_MAX:
            raise ValueError("H5AD counts exceed int64 aggregation capacity")
        if max_value is not None and array.size and int(np.max(array)) > max_value:
            raise ValueError(f"H5AD counts exceed max_value={max_value}")
    else:
        raise ValueError(f"H5AD counts data has unsupported dtype {array.dtype}")
    return np.asarray(array, dtype=np.int64)


def _store_integer_vector(node: Any, name: str) -> np.ndarray:
    """Read one persisted CSR integer vector without coercing malformed data."""
    try:
        values = np.asarray(node[:])
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"H5AD CT store {name} is not a readable dataset") from exc
    if values.ndim != 1 or values.dtype.kind not in "iu":
        raise ValueError(
            f"H5AD CT store {name} must be a one-dimensional integer vector"
        )
    return _checked_integer_vector(values, f"store {name}")


def _store_count_vector(
    node: Any, name: str, *, max_value: int | None
) -> np.ndarray:
    values = _store_integer_vector(node, name)
    return _checked_count_values(values, max_value=max_value)


def _store_string_vector(node: Any, name: str) -> list[str]:
    """Read a persisted one-dimensional UTF-8 identifier/JSON vector."""
    try:
        values = np.asarray(node[:])
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(f"H5AD CT store {name} is not a readable dataset") from exc
    if values.ndim != 1 or values.dtype.kind not in "OSU":
        raise ValueError(f"H5AD CT store {name} must be a one-dimensional string vector")
    output: list[str] = []
    for raw_value in values:
        value = _decode(raw_value)
        if isinstance(value, np.str_):
            value = str(value)
        if not isinstance(value, str):
            raise ValueError(f"H5AD CT store {name} contains a non-string value")
        output.append(value)
    return output


def _store_vector_length(node: Any, name: str, kinds: str) -> int:
    """Validate persisted vector shape/dtype without materializing its values."""
    shape = getattr(node, "shape", None)
    dtype = getattr(node, "dtype", None)
    if shape is None or len(shape) != 1:
        raise ValueError(f"H5AD CT store {name} must be a one-dimensional vector")
    if dtype is None or dtype.kind not in kinds:
        raise ValueError(f"H5AD CT store {name} has an invalid dtype")
    return int(shape[0])


def _strict_integer_scalar(
    value: Any, name: str, *, allow_negative: bool = False
) -> int:
    array = np.asarray(value)
    if array.ndim != 0 or array.dtype.kind not in "iu":
        raise ValueError(f"H5AD CT store {name} must be an integer scalar")
    if array.dtype.kind == "i" and int(array) < 0 and not allow_negative:
        raise ValueError(f"H5AD CT store {name} must be nonnegative")
    if int(array) > _INT64_MAX:
        raise ValueError(f"H5AD CT store {name} exceeds int64 capacity")
    return int(array)


def _store_path_without_symlink(store_path: str | Path) -> Path:
    """Reject a final store symlink before any canonicalisation follows it."""
    candidate = Path(store_path).expanduser()
    try:
        is_symlink = candidate.is_symlink()
    except OSError as exc:
        raise ValueError(f"H5AD CT store path cannot be inspected: {candidate}") from exc
    if is_symlink:
        raise ValueError(f"H5AD CT store path must not be a symlink: {candidate}")
    return candidate


def _validate_store_layout(store_path: Path, *, require_ready: bool) -> dict[str, Any]:
    if Path(store_path).is_symlink():
        raise ValueError(f"H5AD CT store path must not be a symlink: {store_path}")
    with h5py.File(store_path, "r") as handle:
        schema = _strict_integer_scalar(
            handle.attrs.get("store_schema", -1), "store_schema"
        )
        stage = str(_decode(handle.attrs.get("stage", "")))
        state = str(_decode(handle.attrs.get("manifest_state", "")))
        if schema != _STORE_SCHEMA or stage != _STORE_STAGE:
            raise ValueError("H5AD CT store schema/stage is invalid")
        if require_ready and state != _STORE_STATE_READY:
            raise ValueError(f"H5AD CT store is not ready (state={state!r})")
        n_vars = _strict_integer_scalar(handle.attrs.get("n_vars", -1), "n_vars")
        if n_vars <= 0:
            raise ValueError("H5AD CT store has an invalid n_vars")
        if "max_value" not in handle.attrs:
            raise ValueError("H5AD CT store max_value attribute is missing")
        max_attr = _strict_integer_scalar(
            handle.attrs["max_value"], "max_value", allow_negative=True
        )
        if max_attr < -1:
            raise ValueError("H5AD CT store max_value attribute is invalid")
        max_value = None if max_attr == -1 else _validate_max_value(max_attr)

        required = {
            "group_ids",
            "group_index",
            "indptr",
            "indices",
            "data",
            "group_catalog_ids",
            "all_sample_ids",
            "sample_ids",
            "cell_type_ids",
            "group_cell_counts",
            "group_metadata_json",
            "gene_names",
            "group_aggregates",
        }
        missing = sorted(required - set(handle.keys()))
        if missing:
            raise ValueError(f"H5AD CT store is missing datasets/groups {missing}")

        contribution_ids = _store_string_vector(handle["group_ids"], "group_ids")
        n_contributions = len(contribution_ids)
        if n_contributions <= 0:
            raise ValueError("H5AD CT store has no contribution rows")
        contribution_group_indices = _store_integer_vector(
            handle["group_index"], "group_index"
        )
        if len(contribution_group_indices) != n_contributions:
            raise ValueError("H5AD CT store contribution group-index length is invalid")

        catalog_ids = _store_string_vector(
            handle["group_catalog_ids"], "group_catalog_ids"
        )
        n_groups = len(catalog_ids)
        if n_groups <= 0:
            raise ValueError("H5AD CT store has no group catalog rows")
        if len(set(catalog_ids)) != n_groups or any(not value for value in catalog_ids):
            raise ValueError("H5AD CT store group catalog IDs are invalid")
        if len(contribution_group_indices) and (
            int(contribution_group_indices.min()) < 0
            or int(contribution_group_indices.max()) >= n_groups
        ):
            raise ValueError("H5AD CT store contribution group-index values are invalid")
        if any(
            value != catalog_ids[int(index)]
            for value, index in zip(contribution_ids, contribution_group_indices)
        ):
            raise ValueError("H5AD CT store contribution IDs disagree with group-index values")

        indptr = _store_integer_vector(handle["indptr"], "indptr")
        if len(indptr) != n_contributions + 1:
            raise ValueError("H5AD CT store contribution indptr length is invalid")
        if not len(indptr) or int(indptr[0]) != 0:
            raise ValueError("H5AD CT store contribution indptr must start at zero")
        if len(indptr) and np.any(np.diff(indptr) < 0):
            raise ValueError("H5AD CT store contribution indptr is not monotonic")
        indices = _store_integer_vector(handle["indices"], "indices")
        data = _store_count_vector(handle["data"], "data", max_value=max_value)
        if int(indptr[-1]) != len(indices) or int(indptr[-1]) != len(data):
            raise ValueError("H5AD CT store contribution CSR lengths are inconsistent")
        if len(indices) and (
            int(indices.min()) < 0 or int(indices.max()) >= n_vars
        ):
            raise ValueError("H5AD CT store contribution indices are out of bounds")

        sample_ids = _store_string_vector(handle["sample_ids"], "sample_ids")
        cell_type_ids = _store_string_vector(
            handle["cell_type_ids"], "cell_type_ids"
        )
        group_counts = _store_integer_vector(
            handle["group_cell_counts"], "group_cell_counts"
        )
        metadata_json_values = _store_string_vector(
            handle["group_metadata_json"], "group_metadata_json"
        )
        if not (
            len(sample_ids)
            == len(cell_type_ids)
            == len(group_counts)
            == len(metadata_json_values)
            == n_groups
        ):
            raise ValueError("H5AD CT store group metadata/catalog lengths are invalid")
        if np.any(group_counts <= 0):
            raise ValueError("H5AD CT store group cell counts are invalid")
        if any(
            catalog_id != _composite_group_id(sample_id, cell_type_id)
            for catalog_id, sample_id, cell_type_id in zip(
                catalog_ids, sample_ids, cell_type_ids
            )
        ):
            raise ValueError("H5AD CT store group catalog metadata disagrees")

        observed_counts = [0] * n_groups
        for raw_group_index in contribution_group_indices:
            group_index = int(raw_group_index)
            observed_counts[group_index] = _checked_add(
                observed_counts[group_index], 1, _INT64_MAX
            )
        if observed_counts != [int(value) for value in group_counts]:
            raise ValueError("H5AD CT store group cell counts disagree with contributions")

        all_sample_ids = _store_string_vector(
            handle["all_sample_ids"], "all_sample_ids"
        )
        if (
            not all_sample_ids
            or len(set(all_sample_ids)) != len(all_sample_ids)
            or any(not value for value in all_sample_ids)
            or any(sample_id not in all_sample_ids for sample_id in sample_ids)
        ):
            raise ValueError("H5AD CT store all_sample_ids are invalid")
        gene_names = _store_string_vector(handle["gene_names"], "gene_names")
        if (
            len(gene_names) != n_vars
            or any(not value for value in gene_names)
            or len(set(gene_names)) != n_vars
        ):
            raise ValueError("H5AD CT store gene_names are invalid")

        aggregate_group = handle["group_aggregates"]
        if not isinstance(aggregate_group, h5py.Group):
            raise ValueError("H5AD CT store group_aggregates is not a group")
        expected_aggregate_keys = {str(index) for index in range(n_groups)}
        if set(aggregate_group.keys()) != expected_aggregate_keys:
            raise ValueError("H5AD CT store aggregate row count is invalid")
        for index in range(n_groups):
            subgroup = aggregate_group[str(index)]
            if not isinstance(subgroup, h5py.Group) or set(subgroup.keys()) != {
                "indices",
                "data",
            }:
                raise ValueError("H5AD CT store aggregate row is incomplete")
            aggregate_id = _decode(subgroup.attrs.get("group_id", ""))
            if not isinstance(aggregate_id, str) or aggregate_id != catalog_ids[index]:
                raise ValueError("H5AD CT store aggregate row ID is invalid")
            aggregate_indices = _store_integer_vector(
                subgroup["indices"], f"aggregate row {index} indices"
            )
            aggregate_data = _store_count_vector(
                subgroup["data"], f"aggregate row {index} data", max_value=max_value
            )
            if len(aggregate_indices) != len(aggregate_data):
                raise ValueError("H5AD CT store aggregate CSR arrays are inconsistent")
            if len(aggregate_indices) and (
                int(aggregate_indices.min()) < 0
                or int(aggregate_indices.max()) >= n_vars
            ):
                raise ValueError("H5AD CT store aggregate indices are out of bounds")
            if len(aggregate_indices) > 1 and np.any(np.diff(aggregate_indices) <= 0):
                raise ValueError(
                    "H5AD CT store aggregate indices are not sorted and unique"
                )

        columns_json = _decode(handle.attrs.get("group_metadata_columns_json", "[]"))
        try:
            parsed_columns = json.loads(columns_json)
            if (
                not isinstance(parsed_columns, list)
                or not parsed_columns
                or any(not isinstance(value, str) or not value for value in parsed_columns)
                or len(set(parsed_columns)) != len(parsed_columns)
            ):
                raise ValueError
            metadata_columns = list(parsed_columns)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("H5AD CT store metadata columns are malformed") from exc
        try:
            for value in metadata_json_values:
                record = json.loads(value)
                if (
                    not isinstance(record, dict)
                    or set(record) != set(metadata_columns)
                ):
                    raise ValueError
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("H5AD CT store group metadata is malformed") from exc

    return {
        "n_vars": n_vars,
        "n_contributions": n_contributions,
        "n_groups": n_groups,
        "max_value": max_value,
        "metadata_columns": metadata_columns,
        "state": state,
    }




def _validate_store_catalog_layout(
    store_path: Path, *, require_ready: bool
) -> dict[str, Any]:
    """Validate store metadata and CSR shapes without reading every row."""
    if Path(store_path).is_symlink():
        raise ValueError(f"H5AD CT store path must not be a symlink: {store_path}")
    with h5py.File(store_path, "r") as handle:
        schema = _strict_integer_scalar(
            handle.attrs.get("store_schema", -1), "store_schema"
        )
        stage = str(_decode(handle.attrs.get("stage", "")))
        state = str(_decode(handle.attrs.get("manifest_state", "")))
        if schema != _STORE_SCHEMA or stage != _STORE_STAGE:
            raise ValueError("H5AD CT store schema/stage is invalid")
        if require_ready and state != _STORE_STATE_READY:
            raise ValueError(f"H5AD CT store is not ready (state={state!r})")
        n_vars = _strict_integer_scalar(handle.attrs.get("n_vars", -1), "n_vars")
        if n_vars <= 0:
            raise ValueError("H5AD CT store has an invalid n_vars")
        if "max_value" not in handle.attrs:
            raise ValueError("H5AD CT store max_value attribute is missing")
        max_attr = _strict_integer_scalar(
            handle.attrs["max_value"], "max_value", allow_negative=True
        )
        if max_attr < -1:
            raise ValueError("H5AD CT store max_value attribute is invalid")
        max_value = None if max_attr == -1 else _validate_max_value(max_attr)

        required = {
            "group_ids",
            "group_index",
            "indptr",
            "indices",
            "data",
            "group_catalog_ids",
            "all_sample_ids",
            "sample_ids",
            "cell_type_ids",
            "group_cell_counts",
            "group_metadata_json",
            "gene_names",
            "group_aggregates",
        }
        missing = sorted(required - set(handle.keys()))
        if missing:
            raise ValueError(f"H5AD CT store is missing datasets/groups {missing}")

        n_contributions = _store_vector_length(
            handle["group_ids"], "group_ids", "OSU"
        )
        if n_contributions <= 0:
            raise ValueError("H5AD CT store has no contribution rows")
        if (
            _store_vector_length(handle["group_index"], "group_index", "iu")
            != n_contributions
        ):
            raise ValueError("H5AD CT store contribution group-index length is invalid")
        if (
            _store_vector_length(handle["indptr"], "indptr", "iu")
            != n_contributions + 1
        ):
            raise ValueError("H5AD CT store contribution indptr length is invalid")
        indptr_node = handle["indptr"]
        indptr_start = _strict_integer_scalar(indptr_node[0], "indptr first value")
        indptr_stop = _strict_integer_scalar(indptr_node[-1], "indptr final value")
        if indptr_start != 0:
            raise ValueError("H5AD CT store contribution indptr must start at zero")
        indices_length = _store_vector_length(handle["indices"], "indices", "iu")
        data_length = _store_vector_length(handle["data"], "data", "iu")
        if indptr_stop != indices_length or indptr_stop != data_length:
            raise ValueError("H5AD CT store contribution CSR lengths are inconsistent")

        catalog_ids = _store_string_vector(
            handle["group_catalog_ids"], "group_catalog_ids"
        )
        n_groups = len(catalog_ids)
        if n_groups <= 0:
            raise ValueError("H5AD CT store has no group catalog rows")
        if len(set(catalog_ids)) != n_groups or any(not value for value in catalog_ids):
            raise ValueError("H5AD CT store group catalog IDs are invalid")
        for name, kinds in (
            ("sample_ids", "OSU"),
            ("cell_type_ids", "OSU"),
            ("group_cell_counts", "iu"),
            ("group_metadata_json", "OSU"),
        ):
            if _store_vector_length(handle[name], name, kinds) != n_groups:
                raise ValueError(f"H5AD CT store {name} length is invalid")
        all_sample_ids = _store_string_vector(
            handle["all_sample_ids"], "all_sample_ids"
        )
        if (
            not all_sample_ids
            or len(set(all_sample_ids)) != len(all_sample_ids)
            or any(not value for value in all_sample_ids)
        ):
            raise ValueError("H5AD CT store all_sample_ids are invalid")
        if _store_vector_length(handle["gene_names"], "gene_names", "OSU") != n_vars:
            raise ValueError("H5AD CT store gene_names length is invalid")

        aggregate_group = handle["group_aggregates"]
        if not isinstance(aggregate_group, h5py.Group):
            raise ValueError("H5AD CT store group_aggregates is not a group")
        expected_aggregate_keys = {str(index) for index in range(n_groups)}
        if set(aggregate_group.keys()) != expected_aggregate_keys:
            raise ValueError("H5AD CT store aggregate row count is invalid")
        for index in range(n_groups):
            subgroup = aggregate_group[str(index)]
            if not isinstance(subgroup, h5py.Group) or set(subgroup.keys()) != {
                "indices",
                "data",
            }:
                raise ValueError("H5AD CT store aggregate row is incomplete")
            aggregate_indices_length = _store_vector_length(
                subgroup["indices"], f"aggregate row {index} indices", "iu"
            )
            aggregate_data_length = _store_vector_length(
                subgroup["data"], f"aggregate row {index} data", "iu"
            )
            if aggregate_indices_length != aggregate_data_length:
                raise ValueError("H5AD CT store aggregate CSR arrays are inconsistent")

        columns_json = _decode(handle.attrs.get("group_metadata_columns_json", "[]"))
        try:
            parsed_columns = json.loads(columns_json)
            if (
                not isinstance(parsed_columns, list)
                or not parsed_columns
                or any(not isinstance(value, str) or not value for value in parsed_columns)
                or len(set(parsed_columns)) != len(parsed_columns)
            ):
                raise ValueError
            metadata_columns = list(parsed_columns)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("H5AD CT store metadata columns are malformed") from exc

    return {
        "n_vars": n_vars,
        "n_contributions": n_contributions,
        "n_groups": n_groups,
        "max_value": max_value,
        "metadata_columns": metadata_columns,
        "state": state,
    }
def _read_count_chunk(
    counts,
    start: int,
    stop: int,
    n_vars: int,
    max_value: int | None = None,
    *,
    n_obs: int | None = None,
):
    raw_indptr = np.asarray(counts["indptr"][start : stop + 1])
    indptr = _checked_integer_vector(raw_indptr, "indptr")
    if len(indptr) != stop - start + 1 or np.any(np.diff(indptr) < 0):
        raise ValueError("H5AD counts indptr is invalid")
    value_start, value_stop = int(indptr[0]), int(indptr[-1])
    data_length = _node_length(counts["data"])
    indices_length = _node_length(counts["indices"])
    if value_stop > data_length or value_stop > indices_length:
        raise ValueError("H5AD counts indptr exceeds CSR data/index lengths")
    if start == 0 and value_start != 0:
        raise ValueError("H5AD counts indptr does not start at zero")
    if n_obs is not None and stop == n_obs:
        if value_stop != data_length or value_stop != indices_length:
            raise ValueError("H5AD counts final indptr does not match CSR lengths")
    data = _checked_count_values(
        np.asarray(counts["data"][value_start:value_stop]), max_value=max_value
    )
    indices = _checked_integer_vector(
        np.asarray(counts["indices"][value_start:value_stop]), "indices"
    )
    if len(data) != value_stop - value_start or len(indices) != len(data):
        raise ValueError("H5AD counts CSR arrays have inconsistent lengths")
    if len(indices) and (int(indices.min()) < 0 or int(indices.max()) >= n_vars):
        raise ValueError("H5AD counts indices are outside the gene dimension")
    indptr -= value_start
    return data, indices, indptr


def _checked_add(current: int, incoming: int, limit: int) -> int:
    current = int(current)
    incoming = int(incoming)
    if current < 0 or incoming < 0 or current > limit or incoming > limit:
        raise ValueError("checked CSR accumulation received an invalid nonnegative value")
    remaining = limit - current
    if incoming > remaining:
        if limit < _INT64_MAX:
            raise ValueError(f"aggregate value exceeds max_value={limit}")
        raise OverflowError("checked CSR accumulation would overflow int64")
    return current + incoming


def _accumulate_sparse_row(
    target: np.ndarray, indices: np.ndarray, data: np.ndarray, limit: int
) -> None:
    """Add one CSR row into a dense final accumulator with a pre-add guard."""
    for raw_index, raw_value in zip(indices, data):
        index = int(raw_index)
        target[index] = _checked_add(int(target[index]), int(raw_value), limit)


def aggregate_h5ad_counts_by_sample(
    path: str | Path,
    sample_col: str = "Sample",
    metadata_columns: Iterable[str] | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    max_value: int | None = None,
):
    """Aggregate raw CSR counts by sample using bounded cell chunks.

    Returns a dict with ``counts`` shaped genes-by-samples (authoritative
    ``int64``), first-seen ``sample_ids``, gene names, and first-observation
    metadata for the requested columns.  ``max_value`` is an optional checked
    aggregate ceiling; values above it are rejected before this function
    returns, which is required before reticulate/DESeq2 transport.  Biological
    labels are never inferred or inserted into model covariates here.
    """
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    sample_col = str(sample_col)
    if not sample_col:
        raise ValueError("sample_col must be non-empty")
    chunk_size = _validate_chunk_size(chunk_size)
    max_value = _validate_max_value(max_value)
    limit = _INT64_MAX if max_value is None else max_value
    columns = _requested_columns(sample_col, metadata_columns)

    with h5py.File(artifact, "r") as handle:
        obs, counts, n_obs, n_vars, gene_names, _index_name = _open_obs_and_counts(
            handle, artifact
        )
        _validate_obs_columns(obs, columns, n_obs)

        sample_ids: list[str] = []
        sample_to_index: dict[str, int] = {}
        metadata_values = {column: [] for column in columns}
        aggregate_rows: list[np.ndarray] = []

        for start in range(0, n_obs, chunk_size):
            stop = min(n_obs, start + chunk_size)
            chunk_obs = {
                column: read_obs_column_values(obs, column, start, stop)
                for column in columns
            }
            group_ids = np.empty(stop - start, dtype=np.int64)
            for offset, raw_sample in enumerate(chunk_obs[sample_col]):
                sample = _normalise_identifier(raw_sample)
                if sample is None:
                    raise ValueError(
                        f"H5AD obs[{sample_col!r}] contains a missing/blank sample ID"
                    )
                group_id = sample_to_index.get(sample)
                if group_id is None:
                    group_id = len(sample_ids)
                    sample_to_index[sample] = group_id
                    sample_ids.append(sample)
                    aggregate_rows.append(np.zeros(n_vars, dtype=np.int64))
                    for column in columns:
                        value = _metadata_scalar(chunk_obs[column][offset])
                        if column == sample_col:
                            value = sample
                        metadata_values[column].append(value)
                group_ids[offset] = group_id

            data, indices, indptr = _read_count_chunk(
                counts, start, stop, n_vars, max_value=max_value, n_obs=n_obs
            )
            for offset, group_id in enumerate(group_ids):
                row_start, row_stop = int(indptr[offset]), int(indptr[offset + 1])
                _accumulate_sparse_row(
                    aggregate_rows[int(group_id)],
                    indices[row_start:row_stop],
                    data[row_start:row_stop],
                    limit,
                )

    if not sample_ids:
        raise ValueError("H5AD contains no non-empty sample IDs")
    metadata = pd.DataFrame(
        metadata_values,
        index=pd.Index(sample_ids, name=sample_col),
    )
    return {
        "counts": np.asarray(aggregate_rows, dtype=np.int64).T,
        "sample_ids": sample_ids,
        "gene_names": gene_names,
        "metadata": metadata,
    }


def _composite_group_id(sample_id: str, cell_type_id: str) -> str:
    """Use the established CT profile spelling while retaining tuple mapping."""
    return f"Sample={sample_id};cell_type={cell_type_id}"


def _collect_ct_metadata(
    obs,
    sample_col: str,
    cell_type_col: str,
    metadata_columns: Iterable[str] | None,
    n_obs: int,
    chunk_size: int,
):
    requested = _requested_columns(sample_col, metadata_columns)
    columns = list(dict.fromkeys([sample_col, cell_type_col, *requested]))
    available = set(str(name) for name in obs.keys())
    missing = sorted(set(columns) - available)
    if missing:
        raise ValueError(f"H5AD is missing requested obs columns {missing}")

    group_ids: list[str] = []
    sample_ids: list[str] = []
    cell_type_ids: list[str] = []
    group_cell_counts: list[int] = []
    all_sample_ids: list[str] = []
    all_sample_seen: set[str] = set()
    group_to_index: dict[tuple[str, str], int] = {}
    group_string_to_key: dict[str, tuple[str, str]] = {}
    metadata_values = {column: [] for column in columns}
    for _start, _stop, chunk in _read_obs_chunks(obs, columns, n_obs, chunk_size):

        for offset, raw_sample in enumerate(chunk[sample_col]):
            sample = _normalise_identifier(raw_sample)
            if sample is None:
                raise ValueError(
                    f"H5AD obs[{sample_col!r}] contains a missing/blank sample ID"
                )
            if sample not in all_sample_seen:
                all_sample_seen.add(sample)
                all_sample_ids.append(sample)
            cell_type = _normalise_identifier(chunk[cell_type_col][offset])
            if cell_type is None:
                # Missing annotations are excluded from the CT contract, but
                # all count values in their rows are still validated later.
                continue
            key = (sample, cell_type)
            group_index = group_to_index.get(key)
            if group_index is None:
                group_index = len(group_ids)
                group_to_index[key] = group_index
                group_id = _composite_group_id(sample, cell_type)
                previous_key = group_string_to_key.get(group_id)
                if previous_key is not None and previous_key != key:
                    raise ValueError(
                        "composite Sample x cell-type identifiers are ambiguous: "
                        f"{group_id!r}"
                    )
                group_string_to_key[group_id] = key
                group_ids.append(group_id)
                sample_ids.append(sample)
                cell_type_ids.append(cell_type)
                group_cell_counts.append(0)
                for column in columns:
                    if column == sample_col:
                        value = sample
                    elif column == cell_type_col:
                        value = cell_type
                    else:
                        value = _metadata_scalar(chunk[column][offset])
                    metadata_values[column].append(value)
            group_cell_counts[group_index] = _checked_add(
                group_cell_counts[group_index], 1, _INT64_MAX
            )

    if not group_ids:
        raise ValueError("H5AD contains no present Sample x cell-type groups")
    group_metadata = pd.DataFrame(
        metadata_values,
        index=pd.Index(group_ids, name="group_id"),
    )
    return {
        "columns": columns,
        "group_ids": group_ids,
        "sample_ids": sample_ids,
        "all_sample_ids": all_sample_ids,
        "cell_type_ids": cell_type_ids,
        "group_cell_counts": group_cell_counts,
        "group_to_index": group_to_index,
        "metadata": group_metadata,
    }


def _store_manifest_path(store_path: Path) -> Path:
    return Path(f"{store_path}.manifest.json")


def _store_lock_path(store_path: Path) -> Path:
    return Path(f"{store_path}.lock")


def _scheduler_identity() -> dict[str, str]:
    keys = (
        "SLURM_JOB_ID",
        "SLURM_STEP_ID",
        "SLURM_ARRAY_JOB_ID",
        "SLURM_ARRAY_TASK_ID",
        "PBS_JOBID",
        "JOB_ID",
        "ECODA_JOB_ID",
    )
    return {key: os.environ[key] for key in keys if os.environ.get(key)}


def _source_identity(artifact: Path, source_identity: Any) -> Any:
    if source_identity is not None:
        return _json_safe(source_identity)
    stat = artifact.stat()
    return {
        "path": str(artifact.resolve()),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _new_manifest(
    store_path: Path,
    artifact: Path,
    run_id: str | None,
    source_identity: Any,
) -> dict[str, Any]:
    resolved_store = str(store_path.resolve())
    resolved_lock = str(_store_lock_path(store_path).resolve())
    source = _source_identity(artifact, source_identity)
    source_checksum = None
    if isinstance(source, Mapping):
        for key in ("checksum", "md5", "MD5", "sha256", "SHA256"):
            if source.get(key) is not None:
                source_checksum = str(source[key])
                break
    return {
        "schema": _STORE_SCHEMA,
        "stage": _STORE_STAGE,
        "state": _STORE_STATE_WRITING,
        "run_id": str(run_id or os.environ.get("ECODA_RUN_ID") or uuid.uuid4().hex),
        "pid": int(os.getpid()),
        "owner_host": socket.gethostname(),
        "scheduler_identity": _scheduler_identity(),
        "source_identity": source,
        "source_checksum": source_checksum,
        "store_path": resolved_store,
        "lock_path": resolved_lock,
        "created_at": float(time.time()),
    }


def _write_manifest(path: Path, manifest: Mapping[str, Any], *, exclusive: bool = False) -> None:
    payload = json.dumps(_json_safe(dict(manifest)), sort_keys=True, separators=(",", ":"))
    if exclusive:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(payload)
    else:
        path.write_text(payload, encoding="utf-8")


def _utf8_dataset(handle, name: str, values: Iterable[str], *, appendable: bool = False):
    dtype = h5py.string_dtype(encoding="utf-8")
    values = [str(value) for value in values]
    if appendable:
        dataset = handle.create_dataset(
            name, shape=(0,), maxshape=(None,), dtype=dtype, chunks=True
        )
    else:
        dataset = handle.create_dataset(name, data=np.asarray(values, dtype=object), dtype=dtype)
    return dataset


def _append_1d(dataset, values: Any, dtype: Any | None = None) -> None:
    if dtype is None:
        array = np.asarray(values)
    elif dtype is object:
        array = np.asarray(values, dtype=object)
    else:
        array = np.asarray(values, dtype=dtype)
    if array.ndim != 1:
        raise ValueError(f"store append for {dataset.name} is not one-dimensional")
    old_length = int(dataset.shape[0])
    new_length = _checked_add(old_length, int(array.shape[0]), _INT64_MAX)
    dataset.resize((new_length,))
    if array.shape[0]:
        dataset[old_length:new_length] = array


def _append_contributions(
    group_ids_dataset,
    group_index_dataset,
    indptr_dataset,
    indices_dataset,
    data_dataset,
    contribution_group_ids: list[str],
    contribution_group_indices: list[int],
    local_indptr: list[int],
    local_indices: list[int],
    local_data: list[int],
) -> None:
    if len(local_indptr) != len(contribution_group_ids) + 1:
        raise ValueError("store contribution CSR row count is inconsistent")
    base_nnz = int(indptr_dataset[-1])
    absolute_indptr = [
        _checked_add(base_nnz, int(offset), _INT64_MAX)
        for offset in local_indptr[1:]
    ]
    _append_1d(group_ids_dataset, contribution_group_ids, dtype=object)
    _append_1d(group_index_dataset, contribution_group_indices, dtype=np.int64)
    _append_1d(indices_dataset, local_indices, dtype=np.int64)
    _append_1d(data_dataset, local_data, dtype=np.int64)
    _append_1d(indptr_dataset, absolute_indptr, dtype=np.int64)


def _update_sparse_store_row(subgroup, additions: dict[int, int], limit: int) -> None:
    old_indices = _checked_integer_vector(np.asarray(subgroup["indices"][:]), "store indices")
    old_data = _checked_count_values(np.asarray(subgroup["data"][:]), max_value=limit)
    if len(old_indices) != len(old_data):
        raise ValueError("store aggregate CSR arrays have inconsistent lengths")
    values: dict[int, int] = {}
    for index, value in zip(old_indices, old_data):
        column = int(index)
        values[column] = _checked_add(values.get(column, 0), int(value), limit)
    for column, value in additions.items():
        values[column] = _checked_add(values.get(column, 0), int(value), limit)
    ordered = sorted(values.items())
    new_indices = np.asarray([column for column, _value in ordered], dtype=np.int64)
    new_data = np.asarray([value for _column, value in ordered], dtype=np.int64)
    subgroup["indices"].resize((len(new_indices),))
    subgroup["data"].resize((len(new_data),))
    if len(new_indices):
        subgroup["indices"][:] = new_indices
        subgroup["data"][:] = new_data




def _read_manifest(store_path: Path) -> dict[str, Any]:
    manifest_path = _store_manifest_path(store_path)
    if not manifest_path.is_file() or manifest_path.stat().st_size <= 0:
        raise ValueError(f"H5AD CT store ownership manifest is missing: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"H5AD CT store ownership manifest is malformed: {manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise ValueError("H5AD CT store ownership manifest is not an object")
    required = {
        "schema",
        "stage",
        "state",
        "run_id",
        "pid",
        "owner_host",
        "scheduler_identity",
        "source_identity",
        "source_checksum",
        "store_path",
        "lock_path",
        "created_at",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"H5AD CT store ownership manifest is missing {missing}")
    if int(manifest["schema"]) != _STORE_SCHEMA or str(manifest["stage"]) != _STORE_STAGE:
        raise ValueError("H5AD CT store ownership manifest schema/stage is invalid")
    if str(manifest["state"]) not in {_STORE_STATE_WRITING, _STORE_STATE_READY}:
        raise ValueError("H5AD CT store ownership manifest state is invalid")
    if not str(manifest["run_id"]):
        raise ValueError("H5AD CT store ownership manifest run_id is empty")
    try:
        int(manifest["pid"])
        created_at = float(manifest["created_at"])
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("H5AD CT store ownership manifest owner/time is invalid") from exc
    if not np.isfinite(created_at) or created_at < 0:
        raise ValueError("H5AD CT store ownership manifest created_at is invalid")
    if not str(manifest.get("owner_host", "")):
        raise ValueError("H5AD CT store ownership owner_host is invalid")
    if not isinstance(manifest["scheduler_identity"], dict):
        raise ValueError("H5AD CT store ownership scheduler_identity is invalid")
    if Path(str(manifest["store_path"])).resolve() != store_path.resolve():
        raise ValueError("H5AD CT store ownership manifest store_path mismatches target")
    if Path(str(manifest["lock_path"])).resolve() != _store_lock_path(store_path).resolve():
        raise ValueError("H5AD CT store ownership manifest lock_path mismatches target")
    return manifest


def _owner_status(manifest: Mapping[str, Any], lock_path: Path) -> str:
    if lock_path.exists():
        return "locked"
    if str(manifest.get("owner_host", "")) != socket.gethostname():
        return "ambiguous"
    try:
        pid = int(manifest["pid"])
    except (TypeError, ValueError, OverflowError):
        return "ambiguous"
    if pid <= 0:
        return "ambiguous"
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        if manifest.get("scheduler_identity"):
            return "ambiguous"
        return "dead"
    except PermissionError:
        return "ambiguous"
    except OSError as exc:
        if exc.errno == errno.ESRCH:
            if manifest.get("scheduler_identity"):
                return "ambiguous"
            return "dead"
        return "ambiguous"
    return "active"


def audit_h5ad_ct_group_store(
    store_path: str | Path,
    expected_run_id: str | None = None,
    max_age_seconds: float | None = None,
    cleanup: bool = False,
):
    """Audit a CT store without computing counts, optionally cleaning stale data.

    A store is eligible for cleanup only when its manifest and CSR layout are
    valid, its run ID matches ``expected_run_id`` when supplied, its age is at
    least ``max_age_seconds``, its lock is absent, and its recorded local PID
    is demonstrably dead with no scheduler identity.  Missing, active,
    ambiguous, malformed, or future-dated stores fail closed and are never
    deleted.  The returned dict contains ``valid``, ``owner_status``,
    ``expired``, and ``cleanup_performed`` diagnostics.
    """
    requested_store = Path(store_path).expanduser()
    result: dict[str, Any] = {
        "store_path": str(requested_store.absolute()),
        "manifest_path": str(_store_manifest_path(requested_store.absolute())),
        "valid": False,
        "owner_status": "missing",
        "expired": False,
        "cleanup_performed": False,
        "reason": None,
    }
    try:
        # Inspect the caller spelling before resolve(); resolving a final
        # symlink could redirect audit and stale cleanup to an external store.
        requested_store = _store_path_without_symlink(requested_store)
        store = requested_store.resolve()
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        result["reason"] = str(exc)
        return result
    result["store_path"] = str(store)
    result["manifest_path"] = str(_store_manifest_path(store))
    if max_age_seconds is not None:
        try:
            max_age_seconds = float(max_age_seconds)
        except (TypeError, ValueError, OverflowError):
            result["reason"] = "invalid max_age_seconds"
            return result
        if not np.isfinite(max_age_seconds) or max_age_seconds < 0:
            result["reason"] = "invalid max_age_seconds"
            return result
    if not store.is_file():
        result["reason"] = "store is missing"
        return result
    try:
        manifest = _read_manifest(store)
        if expected_run_id is not None and str(manifest["run_id"]) != str(expected_run_id):
            result["reason"] = "run_id mismatch"
            result["owner_status"] = _owner_status(manifest, _store_lock_path(store))
            return result
        layout = _validate_store_layout(store, require_ready=False)
        if str(manifest["state"]) != str(layout["state"]):
            raise ValueError("H5AD CT store manifest/layout state mismatch")
    except (OSError, ValueError, KeyError, TypeError, OverflowError) as exc:
        result["reason"] = str(exc)
        return result

    result["valid"] = True
    result["state"] = str(manifest["state"])
    result["schema"] = int(manifest["schema"])
    result["stage"] = str(manifest["stage"])
    result["owner_status"] = _owner_status(manifest, _store_lock_path(store))
    result["layout"] = layout
    now = float(time.time())
    age_raw = now - float(manifest["created_at"])
    age_seconds = max(0.0, age_raw)
    result["age_seconds"] = age_seconds
    result["expired"] = bool(
        max_age_seconds is not None
        and age_raw >= 0.0
        and age_raw >= max_age_seconds
    )

    eligible = (
        cleanup
        and result["valid"]
        and result["expired"]
        and result["owner_status"] == "dead"
    )
    if eligible:
        try:
            # Recheck the unresolved caller path immediately before unlinking.
            # Never unlink the canonical target of a newly introduced symlink.
            if requested_store.is_symlink():
                result["valid"] = False
                result["reason"] = (
                    f"H5AD CT store path must not be a symlink: {requested_store}"
                )
                return result
            requested_store.unlink()
            _store_manifest_path(requested_store).unlink(missing_ok=True)
            result["cleanup_performed"] = True
        except OSError as exc:
            result["reason"] = f"cleanup failed: {exc}"
    return result


def prepare_h5ad_ct_group_store(
    path: str | Path,
    sample_col: str = "Sample",
    *,
    cell_type_col: str,
    metadata_columns: Iterable[str] | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    max_value: int | None = None,
    store_path: str | Path,
    run_id: str | None = None,
    source_identity: Any = None,
):
    """Write one run-owned sparse Sample x cell-type HDF5 contribution store.

    The first pass reads only ``obs`` and assigns stable first-seen composite
    groups plus first-observation metadata.  Exactly one subsequent pass reads
    ``layers['counts']`` CSR rows.  Contributions are append-only CSR rows,
    while one sparse aggregate row per group is updated on disk with a checked
    int64 (or ``max_value``) guard. Missing cell-type annotations are omitted;
    missing sample IDs, invalid counts, negative signed integers, nonfinite or
    noninteger floats, and checked overflow fail closed before return.

    The returned dict has ``store_path``, first-seen ``group_ids``, aligned
    ``sample_ids``/``cell_type_ids``/``group_cell_counts``, ``group_metadata``,
    ``gene_names``, and ``n_vars``.  No dense all-group gene matrix is retained.
    """
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    sample_col = str(sample_col)
    cell_type_col = str(cell_type_col)
    if not sample_col or not cell_type_col:
        raise ValueError("sample_col and cell_type_col must be non-empty")
    chunk_size = _validate_chunk_size(chunk_size)
    max_value = _validate_max_value(max_value)
    limit = _INT64_MAX if max_value is None else max_value
    store_candidate = _store_path_without_symlink(store_path)
    store = store_candidate.resolve()
    store.parent.mkdir(parents=True, exist_ok=True)
    lock_path = _store_lock_path(store)
    manifest_path = _store_manifest_path(store)
    lock_created = False
    manifest: dict[str, Any] | None = None
    result: dict[str, Any] | None = None

    try:
        with lock_path.open("x", encoding="utf-8") as lock_handle:
            lock_handle.write(str(os.getpid()))
        lock_created = True
        manifest = _new_manifest(store, artifact, run_id, source_identity)
        _write_manifest(manifest_path, manifest, exclusive=True)

        with h5py.File(store, "x") as output:
            with h5py.File(artifact, "r") as handle:
                obs, counts, n_obs, n_vars, gene_names, _index_name = _open_obs_and_counts(
                    handle, artifact
                )
                metadata_pass = _collect_ct_metadata(
                    obs,
                    sample_col,
                    cell_type_col,
                    metadata_columns,
                    n_obs,
                    chunk_size,
                )

                output.attrs["store_schema"] = _STORE_SCHEMA
                output.attrs["stage"] = _STORE_STAGE
                output.attrs["manifest_state"] = _STORE_STATE_WRITING
                output.attrs["n_vars"] = int(n_vars)
                output.attrs["n_obs"] = int(n_obs)
                output.attrs["max_value"] = -1 if max_value is None else int(max_value)
                output.attrs["group_metadata_columns_json"] = json.dumps(
                    metadata_pass["columns"], separators=(",", ":")
                )
                output.attrs["run_id"] = manifest["run_id"]
                output.attrs["owner_pid"] = int(manifest["pid"])
                output.attrs["owner_host"] = manifest["owner_host"]
                output.attrs["scheduler_identity_json"] = json.dumps(
                    manifest["scheduler_identity"], separators=(",", ":")
                )
                output.attrs["source_checksum"] = (
                    "" if manifest["source_checksum"] is None else manifest["source_checksum"]
                )
                output.attrs["created_at"] = float(manifest["created_at"])
                output.attrs["source_identity_json"] = json.dumps(
                    _json_safe(manifest["source_identity"]), separators=(",", ":")
                )

                group_ids = metadata_pass["group_ids"]
                _utf8_dataset(output, "group_catalog_ids", group_ids)
                _utf8_dataset(output, "all_sample_ids", metadata_pass["all_sample_ids"])
                _utf8_dataset(output, "sample_ids", metadata_pass["sample_ids"])
                _utf8_dataset(output, "cell_type_ids", metadata_pass["cell_type_ids"])
                group_counts = np.asarray(metadata_pass["group_cell_counts"], dtype=np.int64)
                output.create_dataset("group_cell_counts", data=group_counts, dtype=np.int64)
                metadata_records = []
                for index in range(len(group_ids)):
                    record = {
                        column: _json_safe(metadata_pass["metadata"].iloc[index][column])
                        for column in metadata_pass["columns"]
                    }
                    metadata_records.append(json.dumps(record, separators=(",", ":")))
                _utf8_dataset(output, "group_metadata_json", metadata_records)
                _utf8_dataset(output, "gene_names", [str(value) for value in gene_names])

                _utf8_dataset(output, "group_ids", [], appendable=True)
                output.create_dataset(
                    "group_index", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True
                )
                output.create_dataset(
                    "indptr", data=np.asarray([0], dtype=np.int64), maxshape=(None,), chunks=True
                )
                output.create_dataset(
                    "indices", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True
                )
                output.create_dataset(
                    "data", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True
                )
                aggregate_group = output.create_group("group_aggregates")
                for group_index, group_id in enumerate(group_ids):
                    subgroup = aggregate_group.create_group(str(group_index))
                    subgroup.attrs["group_id"] = group_id
                    subgroup.create_dataset(
                        "indices", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True
                    )
                    subgroup.create_dataset(
                        "data", shape=(0,), maxshape=(None,), dtype=np.int64, chunks=True
                    )

                seen_counts = [0 for _group_id in group_ids]
                # Re-read only obs metadata for the single raw CSR pass.  The
                # count layer itself is opened/read exactly once per chunk.
                raw_columns = [sample_col, cell_type_col]
                for start, stop, chunk_obs in _read_obs_chunks(
                    obs, raw_columns, n_obs, chunk_size
                ):
                    chunk_group_indices = np.full(stop - start, -1, dtype=np.int64)
                    for offset, raw_sample in enumerate(chunk_obs[sample_col]):
                        sample = _normalise_identifier(raw_sample)
                        if sample is None:
                            raise ValueError(
                                f"H5AD obs[{sample_col!r}] contains a missing/blank sample ID"
                            )
                        cell_type = _normalise_identifier(chunk_obs[cell_type_col][offset])
                        if cell_type is None:
                            continue
                        group_index = metadata_pass["group_to_index"].get((sample, cell_type))
                        if group_index is None:
                            raise ValueError(
                                "H5AD metadata changed between composite passes; "
                                f"missing group {(sample, cell_type)!r}"
                            )
                        chunk_group_indices[offset] = group_index

                    data, indices, indptr = _read_count_chunk(
                        counts, start, stop, n_vars, max_value=max_value, n_obs=n_obs
                    )
                    contribution_group_ids: list[str] = []
                    contribution_group_indices: list[int] = []
                    local_indptr = [0]
                    local_indices: list[int] = []
                    local_data: list[int] = []
                    local_accumulators: dict[int, dict[int, int]] = {}
                    for offset, raw_group_index in enumerate(chunk_group_indices):
                        group_index = int(raw_group_index)
                        row_start, row_stop = int(indptr[offset]), int(indptr[offset + 1])
                        if group_index < 0:
                            continue
                        contribution_group_ids.append(group_ids[group_index])
                        contribution_group_indices.append(group_index)
                        row_indices = indices[row_start:row_stop]
                        row_data = data[row_start:row_stop]
                        for raw_index, raw_value in zip(row_indices, row_data):
                            local_indices.append(int(raw_index))
                            local_data.append(int(raw_value))
                            additions = local_accumulators.setdefault(group_index, {})
                            column = int(raw_index)
                            additions[column] = _checked_add(
                                additions.get(column, 0), int(raw_value), limit
                            )
                        local_indptr.append(
                            _checked_add(
                                local_indptr[-1],
                                row_stop - row_start,
                                _INT64_MAX,
                            )
                        )
                        seen_counts[group_index] = _checked_add(
                            seen_counts[group_index], 1, _INT64_MAX
                        )

                    _append_contributions(
                        output["group_ids"],
                        output["group_index"],
                        output["indptr"],
                        output["indices"],
                        output["data"],
                        contribution_group_ids,
                        contribution_group_indices,
                        local_indptr,
                        local_indices,
                        local_data,
                    )
                    for group_index, additions in local_accumulators.items():
                        _update_sparse_store_row(
                            aggregate_group[str(group_index)], additions, limit
                        )

                if seen_counts != metadata_pass["group_cell_counts"]:
                    raise ValueError("H5AD CT metadata/count passes disagreed on group occupancy")
            output.attrs["manifest_state"] = _STORE_STATE_READY
            output.flush()

        _validate_store_layout(store, require_ready=True)
        manifest["state"] = _STORE_STATE_READY
        manifest["completed_at"] = float(time.time())
        _write_manifest(manifest_path, manifest)
        result = {
            "store_path": str(store),
            "group_ids": list(metadata_pass["group_ids"]),
            "sample_ids": list(metadata_pass["sample_ids"]),
            "all_sample_ids": list(metadata_pass["all_sample_ids"]),
            "cell_type_ids": list(metadata_pass["cell_type_ids"]),
            "group_cell_counts": list(metadata_pass["group_cell_counts"]),
            "group_metadata": metadata_pass["metadata"],
            "gene_names": np.asarray(gene_names, dtype=str),
            "n_vars": int(n_vars),
        }
        return result
    finally:
        if lock_created:
            try:
                lock_path.unlink(missing_ok=True)
            except OSError:
                pass

def _requested_group_ids(group_ids: Iterable[str]) -> list[str]:
    if isinstance(group_ids, str):
        values = [group_ids]
    else:
        values = [str(_decode(value)) for value in group_ids]
    if not values:
        raise ValueError("group_ids must contain at least one contribution ID")
    return values


def read_h5ad_ct_group_store(store_path: str | Path, group_ids: Iterable[str]):
    """Read selected composite groups as sparse CSR rows.

    ``group_ids`` are returned in the requested order.  The ``counts``/``rows``
    value is a sparse CSR matrix shaped selected-groups-by-genes; aggregate
    rows are read from the run-owned sparse CSR store without densifying all
    groups.  ``group_metadata``/``metadata`` is first-observation metadata,
    and aligned sample IDs, cell-type IDs, and cell counts are included.
    """
    store_candidate = _store_path_without_symlink(store_path)
    store = store_candidate.resolve()
    if not store.is_file() or store.stat().st_size <= 0:
        raise ValueError(f"H5AD CT store is missing or empty: {store}")
    manifest = _read_manifest(store)
    if str(manifest["state"]) != _STORE_STATE_READY:
        raise ValueError("H5AD CT store ownership manifest is not ready")
    layout = _validate_store_catalog_layout(store, require_ready=True)
    requested = _requested_group_ids(group_ids)

    with h5py.File(store, "r") as handle:
        all_sample_ids = _store_string_vector(
            handle["all_sample_ids"], "all_sample_ids"
        )
        catalog_ids = _store_string_vector(
            handle["group_catalog_ids"], "group_catalog_ids"
        )
        catalog_to_index = {value: index for index, value in enumerate(catalog_ids)}
        missing = [value for value in requested if value not in catalog_to_index]
        if missing:
            raise ValueError(f"H5AD CT store has no requested group IDs {missing}")
        selected_indices = [catalog_to_index[value] for value in requested]

        row_indptr = [0]
        row_indices: list[int] = []
        row_data: list[int] = []
        for group_index in selected_indices:
            subgroup = handle["group_aggregates"][str(group_index)]
            aggregate_id = _decode(subgroup.attrs.get("group_id", ""))
            if not isinstance(aggregate_id, str) or aggregate_id != catalog_ids[group_index]:
                raise ValueError("H5AD CT store aggregate row ID is invalid")
            indices = _store_integer_vector(
                subgroup["indices"], f"aggregate row {group_index} indices"
            )
            data = _store_count_vector(
                subgroup["data"],
                f"aggregate row {group_index} data",
                max_value=layout["max_value"],
            )
            if len(indices) != len(data):
                raise ValueError("H5AD CT store aggregate CSR arrays are inconsistent")
            if len(indices) and (
                int(indices.min()) < 0 or int(indices.max()) >= layout["n_vars"]
            ):
                raise ValueError("H5AD CT store aggregate indices are out of bounds")
            if len(indices) > 1 and np.any(np.diff(indices) <= 0):
                raise ValueError(
                    "H5AD CT store aggregate indices are not sorted and unique"
                )
            for column, value in zip(indices, data):
                row_indices.append(int(column))
                row_data.append(int(value))
            row_indptr.append(_checked_add(row_indptr[-1], len(indices), _INT64_MAX))
        counts = sparse.csr_matrix(
            (
                np.asarray(row_data, dtype=np.int64),
                np.asarray(row_indices, dtype=np.int64),
                np.asarray(row_indptr, dtype=np.int64),
            ),
            shape=(len(requested), layout["n_vars"]),
            dtype=np.int64,
        )
        columns = layout["metadata_columns"]
        records = []
        for index in selected_indices:
            raw_record = _decode(handle["group_metadata_json"][index])
            if not isinstance(raw_record, str):
                raise ValueError("H5AD CT store selected metadata is not a string")
            try:
                record = json.loads(raw_record)
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ValueError("H5AD CT store selected metadata is malformed") from exc
            if not isinstance(record, dict) or set(record) != set(columns):
                raise ValueError("H5AD CT store selected metadata is malformed")
            records.append(record)
        group_metadata = pd.DataFrame(records, index=pd.Index(requested, name="group_id"))
        group_metadata = group_metadata.reindex(columns=columns)
        sample_ids = [
            str(_decode(handle["sample_ids"][index])) for index in selected_indices
        ]
        cell_type_ids = [
            str(_decode(handle["cell_type_ids"][index])) for index in selected_indices
        ]
        group_cell_counts = [
            int(handle["group_cell_counts"][index]) for index in selected_indices
        ]
        if any(
            sample_id not in all_sample_ids
            or not sample_id.strip()
            or not cell_type_id.strip()
            or count <= 0
            or _composite_group_id(sample_id, cell_type_id)
            != catalog_ids[group_index]
            for sample_id, cell_type_id, count, group_index in zip(
                sample_ids, cell_type_ids, group_cell_counts, selected_indices
            )
        ):
            raise ValueError("H5AD CT store selected group metadata is invalid")
        gene_names = _store_string_vector(handle["gene_names"], "gene_names")
        if (
            len(gene_names) != layout["n_vars"]
            or any(not value for value in gene_names)
            or len(set(gene_names)) != layout["n_vars"]
        ):
            raise ValueError("H5AD CT store gene_names are invalid")
        gene_names = np.asarray(gene_names, dtype=str)

    return {
        "counts": counts,
        "rows": counts,
        "group_ids": requested,
        "all_sample_ids": all_sample_ids,
        "group_metadata": group_metadata,
        "metadata": group_metadata,
        "sample_ids": sample_ids,
        "cell_type_ids": cell_type_ids,
        "group_cell_counts": group_cell_counts,
        "gene_names": gene_names,
        "n_vars": int(layout["n_vars"]),
    }
