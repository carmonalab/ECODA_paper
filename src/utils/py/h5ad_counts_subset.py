"""Load only selected raw-count genes from a persisted H5AD.

MrVI and scPoli require raw counts, but they operate on the stored HVG
ranking. Reading the complete cell-by-gene counts matrix before selecting
HVGs can exceed hundreds of gigabytes. This module filters the persisted CSR
counts row-by-row in bounded chunks and returns a minimal AnnData containing
only the requested genes and observation columns.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:  # import_from_path and package imports use the same module directory
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


def _validate_chunk_size(chunk_size: int) -> int:
    try:
        value = int(chunk_size)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"chunk_size must be a positive integer, got {chunk_size!r}") from exc
    if value <= 0:
        raise ValueError(f"chunk_size must be a positive integer, got {chunk_size!r}")
    return value


def _requested_columns(obs_columns: Iterable[str] | None) -> list[str]:
    if isinstance(obs_columns, str):
        obs_columns = [obs_columns]
    columns = [str(value) for value in (obs_columns or ()) if str(value)]
    if "Sample" not in columns:
        columns.insert(0, "Sample")
    return list(dict.fromkeys(columns))


def _open_layout(handle, artifact: Path):
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
    if var_index_name not in var or "hvg_rank" not in var:
        raise ValueError(f"H5AD var index or hvg_rank is missing: {artifact}")
    gene_names = np.asarray(read_str_dataset(var[var_index_name]), dtype=str)
    if len(gene_names) != n_vars or any(not value.strip() for value in gene_names):
        raise ValueError(f"H5AD var index is invalid: {artifact}")
    if len(set(gene_names)) != n_vars:
        raise ValueError(f"H5AD var index contains duplicate genes: {artifact}")
    ranks = np.asarray(var["hvg_rank"][:], dtype=float).reshape(-1)
    if len(ranks) != n_vars:
        raise ValueError(f"H5AD hvg_rank has the wrong length: {artifact}")
    return obs, counts, n_obs, n_vars, gene_names, ranks, index_name


def read_h5ad_hvg_genes(
    path: str | Path,
    n_hvg: int,
) -> list[str]:
    """Return the top stored HVG names without opening X or layers."""
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    try:
        requested = int(n_hvg)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"n_hvg must be a positive integer, got {n_hvg!r}") from exc
    if requested <= 0:
        raise ValueError(f"n_hvg must be a positive integer, got {n_hvg!r}")
    with h5py.File(artifact, "r") as handle:
        var = handle.get("var")
        if var is None:
            raise ValueError(f"H5AD var is missing: {artifact}")
        index_name = str(_decode(var.attrs.get("_index", "_index")))
        if index_name not in var or "hvg_rank" not in var:
            raise ValueError(f"H5AD var index or hvg_rank is missing: {artifact}")
        gene_names = np.asarray(read_str_dataset(var[index_name]), dtype=str)
        ranks = np.asarray(var["hvg_rank"][:], dtype=float).reshape(-1)
    valid = np.isfinite(ranks)
    if int(valid.sum()) < requested:
        raise ValueError(
            f"Only {int(valid.sum())} genes have a valid hvg_rank, "
            f"but {requested} were requested: {artifact}"
        )
    order = np.argsort(ranks[valid], kind="stable")[:requested]
    selected = gene_names[valid][order].tolist()
    if len(set(selected)) != len(selected):
        raise ValueError(f"Top hvg_rank genes are not unique: {artifact}")
    return selected


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
    indptr -= value_start
    return np.asarray(data, dtype=np.int64), indices, indptr


def _read_obs_frame(obs, columns: list[str], n_obs: int, chunk_size: int, index_name: str):
    available = set(str(name) for name in obs.keys())
    missing = sorted(set(columns) - available)
    if missing:
        raise ValueError(f"H5AD is missing requested obs columns {missing}")
    chunks = {column: [] for column in columns}
    index_chunks = []
    index_node = obs[index_name]
    for start in range(0, n_obs, chunk_size):
        stop = min(n_obs, start + chunk_size)
        index_chunks.append(np.asarray(read_str_dataset(index_node, start, stop), dtype=str))
        for column in columns:
            chunks[column].append(
                np.asarray(read_obs_column_values(obs, column, start, stop), dtype=object)
            )
    values = {column: np.concatenate(parts) for column, parts in chunks.items()}
    index = np.concatenate(index_chunks)
    if len(index) != n_obs or any(not value.strip() for value in index):
        raise ValueError("H5AD obs index contains missing or blank identifiers")
    frame = pd.DataFrame(values, index=pd.Index(index, name=index_name))
    sample = frame["Sample"].astype("string")
    if sample.isna().any() or (sample.str.strip() == "").any():
        raise ValueError("H5AD obs['Sample'] contains missing or blank values")
    return frame


def load_h5ad_counts_subset(
    path: str | Path,
    gene_names: Iterable[str],
    obs_columns: Iterable[str] | None = None,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
):
    """Return counts-only AnnData for selected genes and obs columns."""
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    chunk_size = _validate_chunk_size(chunk_size)
    selected_genes = [str(value) for value in gene_names if str(value)]
    selected_genes = list(dict.fromkeys(selected_genes))
    if not selected_genes:
        raise ValueError("at least one selected gene is required")
    columns = _requested_columns(obs_columns)

    with h5py.File(artifact, "r") as handle:
        obs, counts, n_obs, n_vars, all_genes, ranks, index_name = _open_layout(
            handle, artifact
        )
        positions = {gene: index for index, gene in enumerate(all_genes)}
        missing = sorted(set(selected_genes) - set(positions))
        if missing:
            raise ValueError(f"H5AD is missing selected genes {missing}: {artifact}")
        selected_indices = np.asarray([positions[gene] for gene in selected_genes], dtype=np.int64)
        selected_lookup = np.full(n_vars, -1, dtype=np.int64)
        selected_lookup[selected_indices] = np.arange(len(selected_genes), dtype=np.int64)
        obs_frame = _read_obs_frame(obs, columns, n_obs, chunk_size, index_name)

        data_parts: list[np.ndarray] = []
        index_parts: list[np.ndarray] = []
        row_count_parts: list[np.ndarray] = []
        for start in range(0, n_obs, chunk_size):
            stop = min(n_obs, start + chunk_size)
            data, indices, indptr = _read_count_chunk(counts, start, stop, n_vars)
            selected_positions = selected_lookup[indices]
            keep = selected_positions >= 0
            if np.any(keep):
                row_ids = np.repeat(np.arange(stop - start, dtype=np.int64), np.diff(indptr))
                kept_rows = row_ids[keep]
                row_counts = np.bincount(kept_rows, minlength=stop - start).astype(np.int64)
                data_parts.append(data[keep])
                index_parts.append(selected_positions[keep])
            else:
                row_counts = np.zeros(stop - start, dtype=np.int64)
            row_count_parts.append(row_counts)

    if data_parts:
        data = np.concatenate(data_parts)
        indices = np.concatenate(index_parts)
    else:
        data = np.asarray([], dtype=np.int64)
        indices = np.asarray([], dtype=np.int64)
    row_counts = np.concatenate(row_count_parts)
    indptr = np.concatenate(([0], np.cumsum(row_counts, dtype=np.int64)))
    selected_counts = sparse.csr_matrix(
        (data, indices, indptr), shape=(n_obs, len(selected_genes)), dtype=np.int64
    )
    selected_ranks = ranks[selected_indices]
    minimal = ad.AnnData(
        X=selected_counts,
        obs=obs_frame,
        var=pd.DataFrame(
            {"hvg_rank": selected_ranks},
            index=pd.Index(selected_genes, name="_index"),
        ),
    )
    minimal.layers["counts"] = selected_counts
    minimal.uns["_ecoda_source_shape"] = np.asarray([n_obs, n_vars], dtype=np.int64)
    return minimal
