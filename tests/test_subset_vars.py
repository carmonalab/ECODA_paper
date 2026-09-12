#!/usr/bin/env python3
"""Focused regression checks for subset predicates and sample consistency."""

from __future__ import annotations

from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.py.preprocess_utils import (  # noqa: E402
    apply_subset_vars,
    assert_subset_sample_consistency,
    evaluate_subset_mask,
)


def expect_error(exc_type, callback, message):
    try:
        callback()
    except exc_type:
        return
    raise AssertionError(message)


def make_adata(obs: pd.DataFrame) -> ad.AnnData:
    return ad.AnnData(
        X=np.ones((len(obs), 1), dtype=np.float32),
        obs=obs,
    )


def make_covid_fixture() -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "Sample": ["A", "A", "B", "C", "D", "E"],
            "sampleID": ["sA", "sA", "sB", "sC", "sD", "sE"],
            "PatientID": ["pA", "pA", "pB", "pC", "pD", "pE"],
            "sampling_day": pd.Categorical(
                ["29", "30", "30.5", "control", "unknown", "malformed"]
            ),
        },
        index=pd.Index(["c0", "c1", "c2", "c3", "c4", "c5"], name="cell"),
    )
    return make_adata(obs)


def main() -> None:
    covid = make_covid_fixture()
    covid_rule = {
        "sampling_day": {
            "values": 30,
            "op": "<=",
            "include_values": ["control"],
        }
    }

    mask = evaluate_subset_mask(covid, covid_rule)
    assert isinstance(mask, pd.Series)
    pd.testing.assert_index_equal(mask.index, covid.obs_names)
    assert mask.tolist() == [True, True, False, True, False, False]
    assert list(covid.obs_names[mask.to_numpy()]) == ["c0", "c1", "c3"]
    assert list(covid.obs_names[~mask.to_numpy()]) == ["c2", "c4", "c5"]

    # The two copy modes retain their existing AnnData semantics.
    copied = apply_subset_vars(covid, covid_rule, copy=True)
    assert copied is not covid
    assert not copied.is_view
    copied.X[0, 0] = 99
    assert covid.X[0, 0] == 1

    view = apply_subset_vars(covid, covid_rule, copy=False)
    assert view is not covid
    assert view.is_view
    assert list(view.obs_names) == ["c0", "c1", "c3"]

    split_obs = pd.DataFrame(
        {
            "Sample": ["S", "S"],
            "sampleID": ["sS", "sS"],
            "PatientID": ["pS", "pS"],
            "sampling_day": pd.Categorical(["30", "30.5"]),
        },
        index=["split0", "split1"],
    )
    split = make_adata(split_obs)
    split_mask = pd.Series([True, False], index=split.obs_names, dtype=bool)
    expect_error(
        ValueError,
        lambda: assert_subset_sample_consistency(
            split, split_mask, "sampleID", "Covid19_PBMC / split fixture"
        ),
        "a sample split across retained and dropped cells did not fail",
    )

    repeated_patient_obs = pd.DataFrame(
        {
            "Sample": ["S1", "S2"],
            "sampleID": ["s1", "s2"],
            "PatientID": ["p-repeat", "p-repeat"],
            "sampling_day": pd.Categorical(["30", "31"]),
        },
        index=["repeat0", "repeat1"],
    )
    repeated_patient = make_adata(repeated_patient_obs)
    repeated_mask = pd.Series(
        [True, False], index=repeated_patient.obs_names, dtype=bool
    )
    # Consistency is checked on configured sampleID, not repeated PatientID.
    assert_subset_sample_consistency(
        repeated_patient,
        repeated_mask,
        "sampleID",
        "Covid19_PBMC / repeated PatientID fixture",
    )

    membership = make_adata(
        pd.DataFrame(
            {"label": ["A", "B", "C"]}, index=["m0", "m1", "m2"]
        )
    )
    scalar_in = evaluate_subset_mask(
        membership, {"label": {"values": "A", "op": "in"}}
    )
    scalar_notin = evaluate_subset_mask(
        membership, {"label": {"values": "A", "op": "notin"}}
    )
    assert scalar_in.tolist() == [True, False, False]
    assert scalar_notin.tolist() == [False, True, True]

    comparison = make_adata(
        pd.DataFrame(
            {
                "value": pd.Series(
                    ["", "  ", "not-a-number", "inf", "-inf", None, "2"],
                    index=[f"v{i}" for i in range(7)],
                    dtype="string",
                )
            },
            index=[f"v{i}" for i in range(7)],
        )
    )
    finite_mask = evaluate_subset_mask(
        comparison, {"value": {"values": 3, "op": "<="}}
    )
    assert finite_mask.tolist() == [False, False, False, False, False, False, True]

    expect_error(
        ValueError,
        lambda: evaluate_subset_mask(
            comparison, {"value": {"values": "inf", "op": "<="}}
        ),
        "non-finite comparison threshold did not fail",
    )
    expect_error(
        ValueError,
        lambda: evaluate_subset_mask(
            comparison, {"value": {"values": 3, "op": "approximately"}}
        ),
        "unknown subset operator did not fail",
    )
    expect_error(
        (KeyError, ValueError),
        lambda: evaluate_subset_mask(comparison, {"value": {"values": 3}}),
        "missing subset operator did not fail",
    )
    expect_error(
        KeyError,
        lambda: evaluate_subset_mask(
            comparison, {"missing": {"values": 3, "op": "in"}}
        ),
        "missing subset column did not fail",
    )
    expect_error(
        (KeyError, ValueError),
        lambda: evaluate_subset_mask(comparison, {"value": {"op": "in"}}),
        "malformed subset rule did not fail",
    )

    print("subset predicate and sample-consistency contracts OK")


if __name__ == "__main__":
    main()
