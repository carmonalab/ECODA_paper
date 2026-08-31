#!/usr/bin/env python3
"""Focused assertions for the pre-preprocessing sample-count invariant."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.py.preprocess_utils import remove_low_cellcount_samples


def make_adata(sample_ids):
    obs = pd.DataFrame({"Sample": sample_ids}, index=[f"cell_{i}" for i in range(len(sample_ids))])
    return ad.AnnData(X=np.ones((len(obs), 1), dtype=np.float32), obs=obs)


def expect_error(exc_type, callback, message):
    try:
        callback()
    except exc_type:
        return
    raise AssertionError(message)


def main():
    sample_ids = ["A"] * 500 + ["B"] * 499 + ["C"] * 501
    adata = make_adata(sample_ids)
    original_obs = adata.obs.copy(deep=True)
    original_X = adata.X.copy()
    filtered, removed = remove_low_cellcount_samples(adata)

    assert filtered.n_obs == 1001
    assert removed == {"B": 499}
    assert list(filtered.obs["Sample"].unique()) == ["A", "C"]
    assert filtered.obs_names[0] == "cell_0"
    assert filtered.obs_names[-1] == "cell_1499"
    assert adata.n_obs == 1500
    pd.testing.assert_frame_equal(adata.obs, original_obs)
    np.testing.assert_array_equal(adata.X, original_X)
    assert filtered is not adata
    filtered.X[0, 0] = 7
    assert adata.X[0, 0] == original_X[0, 0]
    filtered.obs.iloc[0, 0] = "detached"
    assert adata.obs.iloc[0, 0] == original_obs.iloc[0, 0]

    complete_input = make_adata(["A"] * 500)
    complete, removed_complete = remove_low_cellcount_samples(
        complete_input, min_cells_per_sample=500
    )
    assert complete.n_obs == 500
    assert removed_complete == {}
    assert complete is complete_input

    expect_error(
        KeyError,
        lambda: remove_low_cellcount_samples(make_adata(["A"]), sample_col="missing"),
        "missing sample column did not fail",
    )
    missing_id = make_adata(["A", None])
    expect_error(
        ValueError,
        lambda: remove_low_cellcount_samples(missing_id),
        "missing sample ID did not fail",
    )
    empty = make_adata(["A"] * 500)
    expect_error(
        ValueError,
        lambda: remove_low_cellcount_samples(empty, min_cells_per_sample=501),
        "empty filtered result did not fail",
    )

    print("sample-count filter contracts OK")


if __name__ == "__main__":
    main()
