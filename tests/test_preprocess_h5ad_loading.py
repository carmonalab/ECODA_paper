#!/usr/bin/env python3
"""Focused contracts for backed H5AD raw-count loading."""
from __future__ import annotations

import tempfile
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.py import preprocess_utils as module


def _write_raw_input(path: Path):
    obs = pd.DataFrame(
        {"Sample": ["sample_a", "sample_a", "sample_b"]},
        index=pd.Index(["cell_1", "cell_2", "cell_3"], name="barcode"),
    )
    normalized_var = pd.DataFrame(
        {"normalized_marker": [0.1, 0.2]},
        index=["gene_1", "gene_2"],
    )
    normalized = sp.csr_matrix(
        np.asarray(
            [
                [0.25, 0.5],
                [0.5, 0.75],
                [0.75, 1.0],
            ],
            dtype=np.float32,
        )
    )
    raw_n_vars = 1205
    raw_var = pd.DataFrame(
        {"raw_marker": np.arange(raw_n_vars, dtype=np.int32)},
        index=[f"raw_gene_{i}" for i in range(raw_n_vars)],
    )
    raw_values = np.zeros((3, raw_n_vars), dtype=np.int32)
    raw_values[0, 0] = 1
    raw_values[0, 1200] = 2
    raw_values[1, 4] = 4
    raw_values[1, 1001] = 5
    raw_values[2, 6] = 6
    raw_values[2, 1204] = 8
    source = ad.AnnData(X=normalized, obs=obs, var=normalized_var)
    source.raw = ad.AnnData(
        X=sp.csr_matrix(raw_values),
        obs=obs.copy(),
        var=raw_var,
    )
    embedding = np.arange(6, dtype=np.float32).reshape(3, 2)
    source.obsm["X_existing"] = embedding
    source.write_h5ad(path)
    return obs, raw_var, raw_values, embedding


def _write_no_raw_input(path: Path):
    obs = pd.DataFrame(
        {"Sample": ["sample_a", "sample_b"]},
        index=["cell_1", "cell_2"],
    )
    values = np.asarray([[0.25, 0.5], [0.75, 1.0]], dtype=np.float32)
    var = pd.DataFrame(index=["gene_1", "gene_2"])
    ad.AnnData(
        X=sp.csr_matrix(values),
        obs=obs,
        var=var,
    ).write_h5ad(path)
    return obs, values


def test_integer_raw_h5ad_is_loaded_without_raw_container():
    with tempfile.TemporaryDirectory(prefix="ecoda-h5ad-loader-") as raw:
        root = Path(raw)
        path = root / "expanded_raw.h5ad"
        obs, raw_var, raw_values, embedding = _write_raw_input(path)

        loaded = module.load_single_input(path.name, root, root)

        assert loaded.shape == raw_values.shape
        assert not loaded.isbacked
        assert loaded.raw is None
        assert sp.issparse(loaded.X)
        np.testing.assert_array_equal(loaded.X.toarray(), raw_values)
        pd.testing.assert_frame_equal(loaded.obs, obs)
        pd.testing.assert_frame_equal(loaded.var, raw_var)
        np.testing.assert_array_equal(loaded.obsm["X_existing"], embedding)


def test_h5ad_without_raw_uses_normal_in_memory_fallback():
    with tempfile.TemporaryDirectory(prefix="ecoda-h5ad-loader-") as raw:
        root = Path(raw)
        path = root / "normalized_only.h5ad"
        obs, values = _write_no_raw_input(path)

        loaded = module.load_single_input(path.name, root, root)

        assert loaded.shape == values.shape
        assert not loaded.isbacked
        assert loaded.raw is None
        np.testing.assert_allclose(loaded.X.toarray(), values)
        pd.testing.assert_frame_equal(loaded.obs, obs)


def test_integer_raw_h5ad_with_misaligned_raw_observations_fails_closed():
    with tempfile.TemporaryDirectory(prefix="ecoda-h5ad-loader-") as raw:
        root = Path(raw)
        path = root / "misaligned_raw.h5ad"
        _write_raw_input(path)

        original_read = module.sc.read_h5ad
        backed = original_read(path, backed="r")

        class RawWithMisalignedObservations:
            X = backed.raw.X
            var = backed.raw.var
            var_names = backed.raw.var_names
            obs_names = pd.Index(["cell_2", "cell_1", "cell_3"])

        class BackedWithMisalignedRaw:
            layers = backed.layers
            obs = backed.obs
            obs_names = backed.obs_names
            obsm = backed.obsm
            raw = RawWithMisalignedObservations()
            file = backed.file

        module.sc.read_h5ad = lambda candidate, **kwargs: BackedWithMisalignedRaw()
        try:
            try:
                module.load_single_input(path.name, root, root)
            except ValueError as exc:
                assert "Raw observations do not align" in str(exc)
            else:
                raise AssertionError("misaligned raw observations were accepted")
        finally:
            module.sc.read_h5ad = original_read


if __name__ == "__main__":
    test_integer_raw_h5ad_is_loaded_without_raw_container()
    test_h5ad_without_raw_uses_normal_in_memory_fallback()
    test_integer_raw_h5ad_with_misaligned_raw_observations_fails_closed()
    print("H5AD loader contracts OK")
