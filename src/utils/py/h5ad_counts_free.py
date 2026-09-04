"""Load benchmark metadata/embeddings without materializing H5AD counts.

The R and Python benchmark methods that consume only observations and stored
PCA embeddings must not open a pinned AnnData file through ``read_h5ad``:
that path can eagerly materialize ``layers['counts']``.  This module reads
only the HDF5 metadata/embedding nodes, then constructs a minimal in-memory
AnnData object with an empty sparse X matrix.  Count-dependent methods keep
their existing full-count path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:  # import_from_path exposes this directory as a top-level module path
    from h5ad_source_identity import read_obs_column_values, read_str_dataset
except ImportError:  # package imports used by focused tests
    from .h5ad_source_identity import read_obs_column_values, read_str_dataset


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


def _read_var_values(node) -> np.ndarray:
    encoding = _decode(node.attrs.get("encoding-type"))
    if encoding == "categorical":
        categories = read_str_dataset(node["categories"])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        if categories.ndim != 1 or len(codes) == 0:
            raise ValueError(f"invalid categorical HDF5 node: {node.name}")
        values = np.full(codes.shape, np.nan, dtype=float)
        valid = codes >= 0
        try:
            values[valid] = np.asarray(categories[codes[valid]], dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"categorical HDF5 node is not numeric: {node.name}") from exc
        return values
    values = np.asarray(node[:]) if hasattr(node, "shape") else np.asarray(node[()])
    return values.reshape(-1)


def _read_obsm_array(node, n_obs: int, key: str) -> np.ndarray:
    if hasattr(node, "shape"):
        values = np.asarray(node[:])
    elif "values" in node and hasattr(node["values"], "shape"):
        values = np.asarray(node["values"][:])
    elif "data" in node and hasattr(node["data"], "shape"):
        values = np.asarray(node["data"][:])
    else:
        raise ValueError(f"unsupported HDF5 obsm encoding for {key!r}")
    if values.ndim != 2 or values.shape[0] != n_obs or values.shape[1] <= 0:
        raise ValueError(f"obsm[{key!r}] has invalid shape {values.shape}")
    if not np.isfinite(values).all():
        raise ValueError(f"obsm[{key!r}] contains nonfinite values")
    return values


def validate_h5ad_counts_free_input(path: str | Path, view: str, method: str) -> None:
    """Validate the persisted H5AD contract without an AnnData backed-open."""
    try:
        from benchmark_h5ad_contract import validate_benchmark_h5ad_path
    except ImportError:
        from .benchmark_h5ad_contract import validate_benchmark_h5ad_path

    validate_benchmark_h5ad_path(path, view, method)

def load_h5ad_counts_free(
    path: str | Path,
    obs_columns: Iterable[str] | None,
    embedding_keys: Iterable[str],
    obs_prefixes: Iterable[str] | None = None,
):
    """Return a minimal AnnData object containing no persisted count values.

    ``obs_columns`` names the exact cell-level metadata columns required by a
    method.  ``obs_prefixes`` adds all HDF5 obs columns with those prefixes;
    this is used for the generated Leiden resolution columns in composition
    methods.  The returned object retains the original cell and gene shape,
    Sample order, requested metadata, ``var['hvg_rank']``, and requested obsm
    embeddings.  X is an empty sparse matrix and no counts layer is attached.
    """
    artifact = Path(path)
    if not artifact.is_file() or artifact.stat().st_size <= 0:
        raise ValueError(f"H5AD is missing or empty: {artifact}")
    requested_obs = [str(value) for value in (obs_columns or ()) if str(value)]
    prefixes = [str(value) for value in (obs_prefixes or ()) if str(value)]
    requested_embeddings = [str(value) for value in embedding_keys if str(value)]
    if not requested_embeddings:
        raise ValueError("counts-free H5AD loading requires at least one embedding key")
    requested_obs = list(dict.fromkeys(requested_obs))
    requested_embeddings = list(dict.fromkeys(requested_embeddings))

    with h5py.File(artifact, "r") as handle:
        if "X" not in handle:
            raise ValueError(f"H5AD lacks X: {artifact}")
        shape = _persisted_shape(handle["X"])
        if len(shape) != 2 or any(value <= 0 for value in shape):
            raise ValueError(f"H5AD X is empty or has no persisted shape: {artifact}")
        n_obs, n_vars = shape

        layers = handle.get("layers")
        if layers is None or "counts" not in layers:
            raise ValueError(f"H5AD lacks layers['counts']: {artifact}")
        counts = layers["counts"]
        counts_shape = _persisted_shape(counts)
        if counts_shape != shape:
            raise ValueError(
                f"H5AD counts shape {counts_shape} does not match X shape {shape}: {artifact}"
            )
        if _decode(counts.attrs.get("encoding-type")) != "csr_matrix":
            raise ValueError(f"H5AD counts layer is not CSR: {artifact}")
        if not all(name in counts for name in ("data", "indices", "indptr")):
            raise ValueError(f"H5AD counts layer is incomplete: {artifact}")

        obs_group = handle.get("obs")
        if obs_group is None or _decode(obs_group.attrs.get("encoding-type")) != "dataframe":
            raise ValueError(f"H5AD obs is not a dataframe: {artifact}")
        index_name = str(_decode(obs_group.attrs.get("_index", "_index")))
        if index_name not in obs_group:
            raise ValueError(f"H5AD lacks obs index {index_name!r}: {artifact}")
        barcodes = read_str_dataset(obs_group[index_name])
        if len(barcodes) != n_obs or any(not str(value).strip() for value in barcodes):
            raise ValueError(f"H5AD obs index is invalid: {artifact}")

        available_obs = [str(name) for name in obs_group.keys() if str(name) != index_name]
        selected_obs = set(requested_obs)
        for prefix in prefixes:
            selected_obs.update(name for name in available_obs if name.startswith(prefix))
        missing_obs = sorted(selected_obs - set(available_obs))
        if missing_obs:
            raise ValueError(f"H5AD is missing requested obs columns {missing_obs}: {artifact}")
        obs_values = {name: read_obs_column_values(obs_group, name) for name in sorted(selected_obs)}
        for name, values in obs_values.items():
            if len(values) != n_obs:
                raise ValueError(f"H5AD obs column length mismatch for {name}: {artifact}")
        obs = pd.DataFrame(obs_values, index=pd.Index(barcodes, name=index_name))

        var_group = handle.get("var")
        if var_group is None:
            raise ValueError(f"H5AD lacks var: {artifact}")
        var_index_name = str(_decode(var_group.attrs.get("_index", "_index")))
        if var_index_name not in var_group or "hvg_rank" not in var_group:
            raise ValueError(f"H5AD lacks var index or hvg_rank: {artifact}")
        gene_names = read_str_dataset(var_group[var_index_name])
        if len(gene_names) != n_vars or any(not str(value).strip() for value in gene_names):
            raise ValueError(f"H5AD var index is invalid: {artifact}")
        hvg_rank = _read_var_values(var_group["hvg_rank"])
        if len(hvg_rank) != n_vars:
            raise ValueError(f"H5AD hvg_rank length mismatch: {artifact}")
        try:
            hvg_rank = np.asarray(hvg_rank, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"H5AD hvg_rank is not numeric: {artifact}") from exc
        var = pd.DataFrame({"hvg_rank": hvg_rank}, index=pd.Index(gene_names, name=var_index_name))

        obsm_group = handle.get("obsm")
        if obsm_group is None:
            raise ValueError(f"H5AD lacks obsm: {artifact}")
        embeddings = {
            key: _read_obsm_array(obsm_group[key], n_obs, key)
            for key in requested_embeddings
            if key in obsm_group
        }
        missing_embeddings = sorted(set(requested_embeddings) - set(embeddings))
        if missing_embeddings:
            raise ValueError(f"H5AD is missing requested obsm keys {missing_embeddings}: {artifact}")

    empty_x = sparse.csr_matrix((n_obs, n_vars), dtype=np.float32)
    minimal = ad.AnnData(X=empty_x, obs=obs, var=var)
    for key, values in embeddings.items():
        minimal.obsm[key] = values
    return minimal
