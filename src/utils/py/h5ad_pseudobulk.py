"""Stream H5AD count aggregation for sample-level pseudobulk preparation.

This module reads the persisted CSR counts layer in bounded cell chunks. It
never constructs the cell-by-gene matrix in memory; only the sample-by-gene
aggregate and selected sample metadata are retained. Count-dependent
pseudobulk still consumes raw counts, while all other metadata/embedding-only
methods use the counts-free loader.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:  # import_from_path exposes this directory as a top-level module path
    from h5ad_source_identity import read_obs_column_values, read_str_dataset
except ImportError:  # package imports used by focused tests
    from .h5ad_source_identity import read_obs_column_values, read_str_dataset


DEFAULT_CHUNK_SIZE = 4096
_INT64_MAX = np.iinfo(np.int64).max


def _decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


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


def _metadata_scalar(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return None
    if isinstance(value, float) and np.isnan(value):
        return None
    return value


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
        raise ValueError(f"H5AD layers['counts'] is not CSR: {artifact}")
    if not all(name in counts for name in ("data", "indices", "indptr")):
        raise ValueError(f"H5AD layers['counts'] is incomplete: {artifact}")
    if _node_length(counts["indptr"]) != n_obs + 1:
        raise ValueError(f"H5AD counts indptr has the wrong length: {artifact}")

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


def _read_obs_chunks(obs, columns: list[str], n_obs: int, chunk_size: int):
    available = set(str(name) for name in obs.keys())
    missing = sorted(set(columns) - available)
    if missing:
        raise ValueError(f"H5AD is missing requested obs columns {missing}")
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


def _collect_sample_metadata(obs, sample_col: str, metadata_columns: list[str], n_obs: int, chunk_size: int):
    sample_ids: list[str] = []
    sample_to_index: dict[str, int] = {}
    values = {column: [] for column in metadata_columns}
    for _start, _stop, chunk in _read_obs_chunks(
        obs, metadata_columns, n_obs, chunk_size
    ):
        sample_values = chunk[sample_col]
        for offset, raw_sample in enumerate(sample_values):
            sample = _metadata_scalar(raw_sample)
            sample = "" if sample is None else str(sample)
            if not sample.strip() or sample.casefold() == "nan":
                raise ValueError(f"H5AD obs[{sample_col!r}] contains a missing/blank sample ID")
            if sample not in sample_to_index:
                sample_to_index[sample] = len(sample_ids)
                sample_ids.append(sample)
                for column in metadata_columns:
                    values[column].append(_metadata_scalar(chunk[column][offset]))
    if not sample_ids:
        raise ValueError("H5AD contains no non-empty sample IDs")
    metadata = pd.DataFrame(
        values,
        index=pd.Index(sample_ids, name=sample_col),
    )
    return sample_ids, sample_to_index, metadata


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


def _read_count_chunk(counts, start: int, stop: int, n_vars: int):
    indptr = np.asarray(counts["indptr"][start : stop + 1], dtype=np.int64)
    if len(indptr) != stop - start + 1 or indptr[0] < 0 or np.any(np.diff(indptr) < 0):
        raise ValueError("H5AD counts indptr is invalid")
    value_start, value_stop = int(indptr[0]), int(indptr[-1])
    data = np.asarray(counts["data"][value_start:value_stop])
    indices = np.asarray(counts["indices"][value_start:value_stop], dtype=np.int64)
    if len(data) != value_stop - value_start or len(indices) != len(data):
        raise ValueError("H5AD counts CSR arrays have inconsistent lengths")
    if len(indices) and (indices.min() < 0 or indices.max() >= n_vars):
        raise ValueError("H5AD counts indices are outside the gene dimension")
    if data.dtype.kind == "f":
        if not np.isfinite(data).all() or np.any(data < 0) or np.any(data != np.floor(data)):
            raise ValueError("H5AD counts data must be finite, nonnegative, integer-valued")
        data = np.rint(data)
    elif data.dtype.kind not in "iu":
        raise ValueError(f"H5AD counts data has unsupported dtype {data.dtype}")
    if len(data) and int(np.max(data)) > _INT64_MAX:
        raise ValueError("H5AD counts exceed int64 aggregation capacity")
    data = np.asarray(data, dtype=np.int64)
    indptr -= value_start
    return data, indices, indptr


def aggregate_h5ad_counts_by_sample(
    path: str | Path,
    sample_col: str = "Sample",
    metadata_columns: Iterable[str] | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
):
    """Aggregate CSR counts by sample using bounded cell chunks.

    Returns a dict with ``counts`` shaped genes-by-samples (int64), ordered
    ``sample_ids``, gene names, and first-observation metadata for the requested
    columns. The metadata is for sample identity/batch bookkeeping; callers
    must keep biological labels out of normalization/model covariates.
    """
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    sample_col = str(sample_col)
    if not sample_col:
        raise ValueError("sample_col must be non-empty")
    chunk_size = _validate_chunk_size(chunk_size)
    columns = _requested_columns(sample_col, metadata_columns)

    with h5py.File(artifact, "r") as handle:
        obs, counts, n_obs, n_vars, gene_names, _index_name = _open_obs_and_counts(
            handle, artifact
        )
        available = set(str(name) for name in obs.keys())
        missing = sorted(set(columns) - available)
        if missing:
            raise ValueError(f"H5AD is missing requested obs columns {missing}")

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
                sample = _metadata_scalar(raw_sample)
                sample = "" if sample is None else str(sample)
                if not sample.strip() or sample.casefold() == "nan":
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
                        metadata_values[column].append(
                            _metadata_scalar(chunk_obs[column][offset])
                        )
                group_ids[offset] = group_id

            data, indices, indptr = _read_count_chunk(counts, start, stop, n_vars)
            chunk_matrix = sparse.csr_matrix(
                (data, indices, indptr), shape=(stop - start, n_vars), dtype=np.int64
            )
            unique_groups, inverse = np.unique(group_ids, return_inverse=True)
            membership = sparse.csr_matrix(
                (
                    np.ones(stop - start, dtype=np.int64),
                    (inverse, np.arange(stop - start, dtype=np.int64)),
                ),
                shape=(len(unique_groups), stop - start),
            )
            grouped = (membership @ chunk_matrix).toarray()
            for local_index, group_id in enumerate(unique_groups):
                aggregate_rows[int(group_id)] += grouped[local_index]

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
