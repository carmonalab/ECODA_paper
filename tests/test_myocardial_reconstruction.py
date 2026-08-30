#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile

import anndata as ad
import numpy as np
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src/2_dataset_specific_preprocessing/1.5.1_reconstruct_myocardial_counts.py"


def load_module():
    spec = importlib.util.spec_from_file_location("myocardial_reconstruction", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    module = load_module()
    raw = sp.csr_matrix(np.array([[1, 2, 0], [0, 1, 2]], dtype=np.float32))

    with tempfile.TemporaryDirectory(prefix="ecoda-myocardial-test-") as directory:
        root = Path(directory)

        reconstructible = ad.AnnData(X=np.log1p(raw.toarray()))
        reconstructible_path = root / "reconstructible.h5ad"
        reconstructible.write_h5ad(reconstructible_path)
        assert module.process_myocardial_file(
            reconstructible_path, reconstructible_path, force=True
        )
        reconstructed = ad.read_h5ad(reconstructible_path)
        np.testing.assert_array_equal(reconstructed.layers["counts"].toarray(), raw.toarray())
        np.testing.assert_array_equal(reconstructed.X.toarray(), raw.toarray())

        valid_layer = ad.AnnData(
            X=np.array([[1000.0, 1001.0, 0.0], [0.0, 1000.0, 1001.0]], dtype=np.float32)
        )
        valid_layer.layers["counts"] = raw.copy()
        valid_layer_path = root / "valid_layer.h5ad"
        valid_layer.write_h5ad(valid_layer_path)
        assert module.process_myocardial_file(valid_layer_path, valid_layer_path, force=True)
        rewritten = ad.read_h5ad(valid_layer_path)
        np.testing.assert_array_equal(rewritten.layers["counts"].toarray(), raw.toarray())
        np.testing.assert_array_equal(rewritten.X.toarray(), raw.toarray())

        overflowing = ad.AnnData(X=np.array([[1000.0, 1001.0]], dtype=np.float32))
        overflowing_path = root / "overflowing.h5ad"
        overflowing.write_h5ad(overflowing_path)
        try:
            module.process_myocardial_file(overflowing_path, overflowing_path, force=True)
        except ValueError as error:
            assert "sanity check failed" in str(error)
        else:
            raise AssertionError("overflowing reconstruction unexpectedly succeeded")

    print("myocardial reconstruction: OK")


if __name__ == "__main__":
    main()
