#!/usr/bin/env python3
"""Focused checks for H5AD observation-index metadata validation."""
from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "ecoda_artifact_contract", ROOT / "src/utils/py/artifact_contract.py"
)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def expect_failure(path: Path, expected: str) -> None:
    try:
        module.validate_artifact(path, "h5ad")
    except ValueError as exc:
        assert expected in str(exc), str(exc)
    else:
        raise AssertionError(f"expected failure containing {expected!r}")


def make_custom_index_h5ad(path: Path) -> None:
    obs = pd.DataFrame(
        {"Sample": ["s1", "s1"]},
        index=pd.Index(["cell-1", "cell-2"], name="barcodes"),
    )
    var = pd.DataFrame(index=["gene-1", "gene-2"])
    ad.AnnData(X=np.ones((2, 2), dtype="float32"), obs=obs, var=var).write_h5ad(path)


def make_literal_index_h5ad(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=np.ones((2, 2), dtype="float32"))
        obs = handle.create_group("obs")
        obs.create_dataset("_index", data=np.asarray(["cell-1", "cell-2"], dtype="S6"))
        handle.create_group("var")
        handle.attrs["shape"] = (2, 2)


def make_counts_only_h5ad(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        layers = handle.create_group("layers")
        counts = layers.create_group("counts")
        counts.attrs["encoding-type"] = "csr_matrix"
        counts.attrs["shape"] = (2, 2)
        counts.create_dataset("data", data=np.asarray([1.0], dtype="float32"))
        counts.create_dataset("indices", data=np.asarray([0], dtype="int32"))
        counts.create_dataset("indptr", data=np.asarray([0, 1, 1], dtype="int32"))

        obs = handle.create_group("obs")
        obs.attrs["_index"] = b"barcodes"
        obs.create_dataset("barcodes", data=np.asarray(["cell-1", "cell-2"], dtype="S6"))
        var = handle.create_group("var")
        var.create_dataset("_index", data=np.asarray(["gene-1", "gene-2"], dtype="S6"))
        handle.attrs["shape"] = (2, 2)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-artifact-contract-") as raw:
        root = Path(raw)
        custom = root / "custom-index.h5ad"
        make_custom_index_h5ad(custom)
        module.validate_artifact(custom, "h5ad")

        with h5py.File(custom, "r+") as handle:
            handle["obs"].attrs["_index"] = "missing"
        expect_failure(custom, "h5ad obs index is missing")

        with h5py.File(custom, "r+") as handle:
            handle["obs"].attrs["_index"] = "   "
        expect_failure(custom, "h5ad obs index is missing")

        with h5py.File(custom, "r+") as handle:
            handle["obs"].attrs["_index"] = 7
        expect_failure(custom, "h5ad obs index is missing")

        with h5py.File(custom, "r+") as handle:
            del handle["obs"].attrs["_index"]
        expect_failure(custom, "h5ad obs index is missing")

        literal = root / "literal-index.h5ad"
        make_literal_index_h5ad(literal)
        module.validate_artifact(literal, "h5ad")

        counts_only = root / "counts-only.h5ad"
        make_counts_only_h5ad(counts_only)
        module.validate_artifact(counts_only, "h5ad")

        missing_matrix = root / "missing-matrix.h5ad"
        make_literal_index_h5ad(missing_matrix)
        with h5py.File(missing_matrix, "r+") as handle:
            del handle["X"]
        expect_failure(missing_matrix, "h5ad missing matrix storage")
    print("artifact contract index checks OK")


if __name__ == "__main__":
    main()
