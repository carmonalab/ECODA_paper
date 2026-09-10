"""Load requested H5AD observations without materializing matrices or embeddings.

The cell-subsetting analysis consumes cell-level metadata only.  This module
checks the persisted AnnData layout with h5py, reads the requested ``obs``
columns and index, and never touches values in ``X``, ``layers['counts']``, or
``obsm``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd

try:  # import_from_path exposes this directory as a top-level module path
    from h5ad_source_identity import read_obs_column_values, read_str_dataset
except ImportError:  # package imports used by focused tests
    from .h5ad_source_identity import read_obs_column_values, read_str_dataset


def _decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else value


def _persisted_shape(node, artifact: Path, label: str) -> tuple[int, ...]:
    shape = getattr(node, "shape", None)
    if shape is None:
        shape = node.attrs.get("shape")
    try:
        values = tuple(int(value) for value in shape)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"H5AD {label} has no valid persisted shape: {artifact}") from exc
    if len(values) != 2 or any(value <= 0 for value in values):
        raise ValueError(f"H5AD {label} is empty or has an invalid shape: {artifact}")
    return values


def _node_length(node, artifact: Path, label: str) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "keys"):
        if "codes" in node:
            shape = node["codes"].shape
        elif "values" in node:
            shape = node["values"].shape
    if shape is None or len(shape) != 1:
        raise ValueError(f"H5AD {label} is not a one-dimensional vector: {artifact}")
    return int(shape[0])


def _requested_columns(obs_columns: Iterable[str] | None) -> list[str]:
    if obs_columns is None:
        return []
    if isinstance(obs_columns, str):
        obs_columns = [obs_columns]
    columns = [str(value) for value in obs_columns if str(value)]
    return list(dict.fromkeys(columns))


def load_h5ad_obs_free(
    path: str | Path,
    obs_columns: Iterable[str] | None,
) -> pd.DataFrame:
    """Return only requested ``obs`` columns and the original obs index.

    H5AD ``X`` and ``layers['counts']`` are checked for matching, non-empty
    persisted two-dimensional shapes, but their values are never read.  The
    ``var`` index and ``hvg_rank`` node are checked for the required schema;
    ``obsm`` is not inspected.
    """
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    requested = _requested_columns(obs_columns)

    with h5py.File(artifact, "r") as handle:
        x = handle.get("X")
        if x is None:
            raise ValueError(f"H5AD lacks X: {artifact}")
        shape = _persisted_shape(x, artifact, "X")
        n_obs, n_vars = shape

        layers = handle.get("layers")
        counts = layers.get("counts") if layers is not None else None
        if counts is None:
            raise ValueError(f"H5AD lacks layers['counts']: {artifact}")
        counts_shape = _persisted_shape(counts, artifact, "layers['counts']")
        if counts_shape != shape:
            raise ValueError(
                f"H5AD counts shape {counts_shape} does not match X shape {shape}: {artifact}"
            )

        obs_group = handle.get("obs")
        if obs_group is None or _decode(obs_group.attrs.get("encoding-type")) != "dataframe":
            raise ValueError(f"H5AD obs is not a dataframe: {artifact}")
        index_name = str(_decode(obs_group.attrs.get("_index", "_index")))
        if index_name not in obs_group:
            raise ValueError(f"H5AD lacks obs index {index_name!r}: {artifact}")
        index_node = obs_group[index_name]
        if _node_length(index_node, artifact, "obs index") != n_obs:
            raise ValueError(f"H5AD obs index length mismatch: {artifact}")
        obs_index = np.asarray(read_str_dataset(index_node))
        if (
            obs_index.ndim != 1
            or len(obs_index) != n_obs
            or any(not str(value).strip() for value in obs_index)
        ):
            raise ValueError(f"H5AD obs index is invalid: {artifact}")

        available = {str(name) for name in obs_group.keys() if str(name) != index_name}
        missing = sorted(set(requested) - available)
        if missing:
            raise ValueError(f"H5AD is missing requested obs columns {missing}: {artifact}")
        obs_values = {}
        for column in requested:
            node = obs_group[column]
            if _node_length(node, artifact, f"obs column {column!r}") != n_obs:
                raise ValueError(f"H5AD obs column length mismatch for {column}: {artifact}")
            values = read_obs_column_values(obs_group, column)
            if len(values) != n_obs:
                raise ValueError(f"H5AD obs column length mismatch for {column}: {artifact}")
            obs_values[column] = values

        var_group = handle.get("var")
        if var_group is None:
            raise ValueError(f"H5AD lacks var: {artifact}")
        var_index_name = str(_decode(var_group.attrs.get("_index", "_index")))
        if var_index_name not in var_group:
            raise ValueError(f"H5AD lacks var index {var_index_name!r}: {artifact}")
        if "hvg_rank" not in var_group:
            raise ValueError(f"H5AD lacks var['hvg_rank']: {artifact}")
        var_index_node = var_group[var_index_name]
        if _node_length(var_index_node, artifact, "var index") != n_vars:
            raise ValueError(f"H5AD var index length mismatch: {artifact}")
        var_index = np.asarray(read_str_dataset(var_index_node))
        if (
            var_index.ndim != 1
            or len(var_index) != n_vars
            or any(not str(value).strip() for value in var_index)
        ):
            raise ValueError(f"H5AD var index is invalid: {artifact}")

    return pd.DataFrame(obs_values, index=pd.Index(obs_index, name=index_name))
