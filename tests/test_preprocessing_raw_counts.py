#!/usr/bin/env python3
"""Focused assertions for raw-count adoption and shape alignment."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "ecoda_preprocess", ROOT / "src/3_scrnaseq_preprocessing/1.1.1_preprocess.py"
)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def make_adata():
    n_obs, n_vars = 4, 100
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    gene_names = [f"gene_{i}" for i in range(n_vars)]
    var = pd.DataFrame(
        {"normalized_marker": np.arange(n_vars, dtype=np.int32)},
        index=gene_names,
    )
    raw_var = pd.DataFrame(
        {"raw_marker": np.arange(n_vars, dtype=np.int32)},
        index=gene_names,
    )
    raw_values = np.arange(1, n_obs * n_vars + 1, dtype=np.int32).reshape(
        n_obs, n_vars
    )
    normalized_values = np.full((n_obs, n_vars), 0.25, dtype=np.float32)
    adata = ad.AnnData(
        X=sp.csr_matrix(normalized_values), obs=obs, var=var
    )
    adata.raw = ad.AnnData(
        X=sp.csr_matrix(raw_values), obs=obs.copy(), var=raw_var
    )
    return adata, raw_values, raw_var


def dense(matrix):
    return matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)



def test_sparse_scale_matches_dense_scanpy_after_centering():
    values = np.zeros((120, 4), dtype=np.float64)
    values[0, 0] = 100
    values[1, 0] = 2
    values[0, 1] = 1
    values[1, 1] = 2
    values[2, 1] = 3
    values[0, 3] = -100
    values[1, 3] = -2

    sparse_scaled = module._scale_sparse_for_pca(sp.csr_matrix(values))
    assert sp.isspmatrix_csr(sparse_scaled)

    dense_adata = ad.AnnData(X=values.copy())
    module.sc.pp.scale(dense_adata, max_value=10)

    sparse_centered = dense(sparse_scaled)
    sparse_centered -= sparse_centered.mean(axis=0)
    dense_values = np.asarray(dense_adata.X)
    dense_centered = dense_values - dense_values.mean(axis=0)
    np.testing.assert_allclose(sparse_centered, dense_centered, rtol=1e-10, atol=1e-10)
    assert dense_values[0, 0] == 10
    assert dense_values[0, 3] == -10


def test_compute_pca_uses_lightweight_sparse_subset():
    n_obs, n_vars = 6, 4
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    values = np.arange(n_obs * n_vars, dtype=np.float64).reshape(n_obs, n_vars)
    adata = ad.AnnData(X=sp.csr_matrix(values), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(values)

    original_pca = module.sc.pp.pca
    captured = {}

    def fake_pca(sub, n_comps, svd_solver):
        captured["sub"] = sub
        sub.obsm["X_pca"] = np.zeros((sub.n_obs, n_comps), dtype=np.float64)

    module.sc.pp.pca = fake_pca
    try:
        sub = module.compute_pca_and_store(
            adata, ["gene_0", "gene_2"], "sparse_subset", n_pcs=2
        )
    finally:
        module.sc.pp.pca = original_pca

    assert sub is captured["sub"]
    assert sp.isspmatrix_csr(sub.X)
    assert len(sub.layers) == 0
    assert list(sub.var_names) == ["gene_0", "gene_2"]
    assert "X_pca_sparse_subset" in adata.obsm

def test_integer_raw_counts_are_adopted_and_vaulted():
    module.standardize_gene_symbols = lambda _: None
    adata, raw_values, raw_var = make_adata()

    processed = module.base_preprocessing(adata)

    assert processed is adata
    np.testing.assert_array_equal(dense(processed.layers["counts"]), raw_values)
    assert np.isfinite(dense(processed.X)).all()
    assert not np.array_equal(dense(processed.X), raw_values)
    assert processed.raw is None
    np.testing.assert_array_equal(
        processed.var["raw_marker"].to_numpy(), raw_var["raw_marker"].to_numpy()
    )
    assert "normalized_marker" not in processed.var


def test_integer_raw_shape_mismatch_fails_closed():
    module.standardize_gene_symbols = lambda _: None
    adata, raw_values, raw_var = make_adata()
    adata.raw = ad.AnnData(
        X=sp.csr_matrix(raw_values[:, :-1]),
        obs=adata.obs.copy(),
        var=raw_var.iloc[:-1].copy(),
    )

    try:
        module.base_preprocessing(adata)
    except ValueError as exc:
        assert "shape" in str(exc).lower()
    else:
        raise AssertionError("incompatible raw matrix shape was not rejected")

def main():
    test_sparse_scale_matches_dense_scanpy_after_centering()
    test_compute_pca_uses_lightweight_sparse_subset()
    test_integer_raw_counts_are_adopted_and_vaulted()
    test_integer_raw_shape_mismatch_fails_closed()
    print("raw-count preprocessing contracts OK")


if __name__ == "__main__":
    main()
