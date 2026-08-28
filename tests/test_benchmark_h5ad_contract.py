#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from anndata import AnnData

from src.utils.py.benchmark_h5ad_contract import (
    REQUIRED_OBSM,
    validate_benchmark_h5ad_contract,
)


def make_valid_adata(view="benchmark_analysis"):
    n_cells = 4
    n_genes = 3000
    obs = pd.DataFrame(
        {"Sample": ["s1", "s1", "s2", "s2"]},
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        {"hvg_rank": np.arange(n_genes, dtype=float)},
        index=[f"g{i}" for i in range(n_genes)],
    )
    adata = AnnData(
        X=np.ones((n_cells, n_genes), dtype=float),
        obs=obs,
        var=var,
    )
    adata.layers["counts"] = np.ones((n_cells, n_genes), dtype=float)
    for key in REQUIRED_OBSM[view]:
        adata.obsm[key] = np.ones((n_cells, 2), dtype=float)
    return adata


def assert_contract_error(adata, expected):
    try:
        validate_benchmark_h5ad_contract(adata, "benchmark_analysis", "test")
    except ValueError as exc:
        assert expected in str(exc), str(exc)
    else:
        raise AssertionError(f"Expected contract failure mentioning {expected!r}")


def main():
    for view in REQUIRED_OBSM:
        validate_benchmark_h5ad_contract(
            make_valid_adata(view), view, "test"
        )

    missing_counts = make_valid_adata()
    del missing_counts.layers["counts"]
    assert_contract_error(missing_counts, "layers['counts']")

    missing_embedding = make_valid_adata()
    del missing_embedding.obsm["X_pca_benchmark_analysis_hvg3000"]
    assert_contract_error(missing_embedding, "X_pca_benchmark_analysis_hvg3000")
    assert "batch_effect_analysis" not in REQUIRED_OBSM
    try:
        validate_benchmark_h5ad_contract(
            make_valid_adata("batch_effect_uncorrected"),
            "batch_effect_analysis",
            "test",
        )
    except ValueError:
        pass
    else:
        raise AssertionError("legacy batch-effect view was accepted")

    print("benchmark h5ad contract checks OK")


if __name__ == "__main__":
    main()
