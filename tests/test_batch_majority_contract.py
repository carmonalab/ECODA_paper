#!/usr/bin/env python3
"""Focused regression checks for obs-only H5AD sample metadata export.

The fixtures exercise the exporter at its source boundary: only observation
metadata is reduced to one row per Sample, while the H5AD expression and
count storage remains unrelated to the metadata contract.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import json
from pathlib import Path
import sys
import tempfile

import anndata as ad
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.utils.py.export_h5ad_sample_metadata import (  # noqa: E402
    export,
    read_obs_metadata,
)
from src.utils.py.batch_contract import (  # noqa: E402
    BatchContractError,
    composite_token,
)


ASSAY_POLICY = {
    "sample_aggregation": "majority_v1",
    "majority_keys": ["assay"],
    "accepted_sentinel_values": {},
}

BREAST_POLICY = {
    "sample_aggregation": "majority_v1",
    "majority_keys": ["suspension_dissociation_time"],
    "accepted_sentinel_values": {
        "suspension_dissociation_time": ["unknown"],
    },
}



def _write_h5ad(path: Path, obs: pd.DataFrame) -> dict[str, np.ndarray]:
    """Write a small real AnnData fixture with auxiliary storage present."""

    n_obs = len(obs)
    x = np.arange(n_obs * 2, dtype=np.float32).reshape(n_obs, 2)
    counts = (x + 10).astype(np.int32)
    normalized = (x + 0.5).astype(np.float32)
    raw_x = (x + 100).astype(np.float32)
    var = pd.DataFrame(index=["gene_0", "gene_1"])

    # Keep nullable StringDtype values representable in this small fixture;
    # missing values become Python None while blank and sentinel strings stay
    # observable to the obs-only reader.
    obs = obs.copy()
    for column in obs.columns:
        if isinstance(obs[column].dtype, pd.StringDtype):
            series = obs[column]
            obs[column] = series.astype(object).where(series.notna(), None)

    data = ad.AnnData(X=x.copy(), obs=obs, var=var.copy())
    data.layers["counts"] = counts.copy()
    data.layers["normalized"] = normalized.copy()
    data.raw = ad.AnnData(X=raw_x.copy(), obs=obs.copy(), var=var.copy())
    data.write_h5ad(path)
    return {
        "X": x,
        "counts": counts,
        "normalized": normalized,
        "raw": raw_x,
    }
def _write_obs_only_h5ad(
    path: Path, columns: dict[str, list[object]], index: list[str]
) -> None:
    """Write a valid tiny obs-only H5AD, intentionally omitting X/raw/layers."""

    path.parent.mkdir(parents=True, exist_ok=True)
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "index"
        index_node = obs.create_dataset(
            "index",
            data=np.asarray(index, dtype=object),
            dtype=string_dtype,
        )
        index_node.attrs["encoding-type"] = "string-array"
        index_node.attrs["encoding-version"] = "0.2.0"
        for name, values in columns.items():
            array = np.asarray(values)
            if array.dtype.kind in {"U", "O"}:
                node = obs.create_dataset(
                    name,
                    data=np.asarray(values, dtype=object),
                    dtype=string_dtype,
                )
                node.attrs["encoding-type"] = "string-array"
                node.attrs["encoding-version"] = "0.2.0"
            else:
                node = obs.create_dataset(name, data=array)
                node.attrs["encoding-type"] = "array"
                node.attrs["encoding-version"] = "0.2.0"





def _read(
    path: Path,
    required_columns: list[str],
    *,
    batch_keys: list[str],
    policy: dict[str, object] | None,
    biological_column: str,
) -> pd.DataFrame:
    return read_obs_metadata(
        path,
        "Sample",
        required_columns,
        [],
        chunk_size=2,
        batch_keys=batch_keys,
        batch_metadata_policy=policy,
        biological_column=biological_column,
    )



def _expect_value_error(callback, context: str) -> None:
    try:
        callback()
    except ValueError as exc:
        assert "assay" in str(exc), f"{context} failed for the wrong metadata key: {exc}"
    else:
        raise AssertionError(f"{context} unexpectedly passed")



def _majority_observations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Sample": ["S2"] * 3 + ["S1"] * 5,
            # S1's winner occurs only 2/5 times.  It is unique, but is not a
            # conventional >50% majority and therefore exercises the absence
            # of a fraction threshold.
            "assay": [
                "assay-b-first",
                "assay-b-first",
                "assay-b-other",
                "assay-first",
                "assay-winner",
                "assay-other",
                "assay-winner",
                "assay-last",
            ],
            # ``unknown`` is deliberately on a non-majority technical key.  It
            # must not be interpreted as a missing value for ``assay`` policy.
            "site": [
                "site-s2-first",
                "site-s2-second",
                "site-s2-last",
                "site-s1-first",
                "unknown",
                "site-s1-third",
                "site-s1-fourth",
                "site-s1-last",
            ],
            "cell_type": [
                "label-s2-first",
                "label-s2-second",
                "label-s2-last",
                "label-s1-first",
                "label-s1-second",
                "label-s1-third",
                "label-s1-fourth",
                "label-s1-last",
            ],
        },
        index=[f"cell-{index}" for index in range(8)],
    )



def _check_assay_majority_fixture(root: Path) -> None:
    observations = _majority_observations()
    path = root / "assay-majority.h5ad"
    arrays = _write_h5ad(path, observations)
    before = path.read_bytes()

    required = ["Sample", "assay", "site", "cell_type"]
    baseline = _read(
        path,
        required,
        batch_keys=["assay", "site"],
        policy=None,
        biological_column="cell_type",
    )
    exported = _read(
        path,
        required,
        batch_keys=["assay", "site"],
        policy=ASSAY_POLICY,
        biological_column="cell_type",
    )

    assert list(exported.columns) == required
    assert baseline["Sample"].tolist() == ["S2", "S1"]
    assert exported["Sample"].tolist() == ["S2", "S1"]
    assert baseline["assay"].tolist() == ["assay-b-first", "assay-first"]
    assert exported["assay"].tolist() == ["assay-b-first", "assay-winner"]
    assert exported["site"].tolist() == ["site-s2-first", "site-s1-first"]
    assert exported["cell_type"].tolist() == ["label-s2-first", "label-s1-first"]

    changed_columns = [
        column
        for column in exported.columns
        if not exported[column].equals(baseline[column])
    ]
    assert changed_columns == ["assay"]
    assert "unknown" in observations["site"].tolist()

    # The exporter consumes only obs metadata: expression, raw, and layers are
    # present in the fixture, but neither the file nor their values change.
    assert path.read_bytes() == before
    restored = ad.read_h5ad(path)
    np.testing.assert_array_equal(restored.X, arrays["X"])
    np.testing.assert_array_equal(restored.layers["counts"], arrays["counts"])
    np.testing.assert_array_equal(restored.layers["normalized"], arrays["normalized"])
    assert restored.raw is not None
    np.testing.assert_array_equal(restored.raw.X, arrays["raw"])



def _check_tie_fixture(root: Path) -> None:
    observations = pd.DataFrame(
        {
            "Sample": ["tie"] * 4,
            "assay": ["assay-a", "assay-a", "assay-b", "assay-b"],
            "site": ["site"] * 4,
            "cell_type": ["label"] * 4,
        },
        index=[f"tie-cell-{index}" for index in range(4)],
    )
    path = root / "assay-tie.h5ad"
    _write_h5ad(path, observations)
    _expect_value_error(
        lambda: _read(
            path,
            ["Sample", "assay", "site", "cell_type"],
            batch_keys=["assay", "site"],
            policy=ASSAY_POLICY,
            biological_column="cell_type",
        ),
        "exact assay tie",
    )



def _check_invalid_assay_values(root: Path) -> None:
    invalid_cases = {
        "missing": pd.Series(["assay-valid", pd.NA, "assay-other"], dtype="string"),
        "blank": pd.Series(["assay-valid", "   ", "assay-other"], dtype="string"),
        "unknown": pd.Series(["assay-valid", "unknown", "assay-other"], dtype="string"),
        "nan": pd.Series([1.0, np.nan, 2.0], dtype="float64"),
        "positive-infinity": pd.Series([1.0, np.inf, 2.0], dtype="float64"),
        "negative-infinity": pd.Series([1.0, -np.inf, 2.0], dtype="float64"),
    }
    for name, assay in invalid_cases.items():
        observations = pd.DataFrame(
            {
                "Sample": ["invalid"] * len(assay),
                "site": ["site"] * len(assay),
                "cell_type": ["label"] * len(assay),
            }
        )
        observations["assay"] = assay
        observations.index = [f"{name}-cell-{index}" for index in range(len(assay))]
        path = root / f"invalid-assay-{name}.h5ad"
        _write_h5ad(path, observations)
        _expect_value_error(
            lambda path=path: _read(
                path,
                ["Sample", "assay", "site", "cell_type"],
                batch_keys=["assay", "site"],
                policy=ASSAY_POLICY,
                biological_column="cell_type",
            ),
            f"invalid assay value ({name})",
        )



def _check_breast_fixture(root: Path) -> None:
    observations = pd.DataFrame(
        {
            "Sample": ["breast-unknown"] * 4 + ["breast-recorded"] * 3,
            "suspension_dissociation_time": [
                "unknown",
                "12h",
                "unknown",
                "24h",
                "6h",
                "12h",
                "12h",
            ],
            "disease": [
                "disease-first",
                "disease-second",
                "disease-third",
                "disease-fourth",
                "disease-recorded-first",
                "disease-recorded-second",
                "disease-recorded-third",
            ],
        },
        index=[f"breast-cell-{index}" for index in range(7)],
    )
    path = root / "breast-majority.h5ad"
    _write_h5ad(path, observations)
    exported = _read(
        path,
        ["Sample", "suspension_dissociation_time", "disease"],
        batch_keys=["suspension_dissociation_time"],
        policy=BREAST_POLICY,
        biological_column="disease",
    )
    assert exported["Sample"].tolist() == ["breast-unknown", "breast-recorded"]
    assert exported["suspension_dissociation_time"].tolist() == ["unknown", "12h"]
    assert exported["disease"].tolist() == ["disease-first", "disease-recorded-first"]


def _expect_composite_failure(callback, context: str) -> None:
    try:
        callback()
    except (BatchContractError, ValueError) as exc:
        message = str(exc)
        assert "sentinel" in message or "batch value" in message, (
            f"{context} failed for an unrelated reason: {exc}"
        )
    else:
        raise AssertionError(f"{context} unexpectedly passed")


def _check_breast_three_key_composite() -> None:
    """Allow only Breast's configured literal sentinel on its own key."""

    keys = (
        "assay",
        "sequencing_platform",
        "suspension_dissociation_time",
    )
    accepted = BREAST_POLICY["accepted_sentinel_values"]
    recorded = ("assay-a", "platform-a", "12h")
    with_unknown = ("assay-a", "platform-a", "unknown")

    # The optional mapping changes only the approved literal category; all
    # ordinary values retain the exact default token.
    assert composite_token(keys, recorded) == composite_token(
        keys,
        recorded,
        accepted_sentinel_values=accepted,
    )
    token = composite_token(
        keys,
        with_unknown,
        accepted_sentinel_values=accepted,
    )
    assert ",9:733a756e6b6e6f776e" in token  # normal ``s:unknown`` category

    _expect_composite_failure(
        lambda: composite_token(
            keys,
            ("assay-a", "platform-a", "UNKNOWN"),
            accepted_sentinel_values=accepted,
        ),
        "nonliteral case variant of the configured sentinel",
    )

    # No mapping remains strict, as does a mapping attached to another key.
    _expect_composite_failure(
        lambda: composite_token(keys, with_unknown),
        "Breast unknown without an approved policy",
    )
    _expect_composite_failure(
        lambda: composite_token(
            keys,
            ("unknown", "platform-a", "12h"),
            accepted_sentinel_values=accepted,
        ),
        "unknown on an unapproved Breast key",
    )
    _expect_composite_failure(
        lambda: composite_token(
            keys,
            with_unknown,
            accepted_sentinel_values={"assay": ["unknown"]},
        ),
        "Breast policy attached to an unapproved key",
    )

    # Mapping shape and values are fail-closed rather than broadening the
    # exception to arbitrary strings or datasets.
    _expect_composite_failure(
        lambda: composite_token(
            keys,
            with_unknown,
            accepted_sentinel_values={"not_a_batch_key": ["unknown"]},
        ),
        "accepted sentinel mapping with an unknown key",
    )
    _expect_composite_failure(
        lambda: composite_token(
            keys,
            with_unknown,
            accepted_sentinel_values={
                "suspension_dissociation_time": "unknown"
            },
        ),
        "accepted sentinel mapping with a scalar value",
    )
    _expect_composite_failure(
        lambda: composite_token(
            keys,
            with_unknown,
            accepted_sentinel_values={
                "suspension_dissociation_time": ["not-a-sentinel"]
            },
        ),
        "accepted sentinel mapping with an ordinary category",
    )

def _check_export_boundary(root: Path) -> None:
    """Exercise export() against an H5AD that contains only its obs group."""

    input_path = (root / "BreastCncr_processed.h5ad").resolve()
    output_path = (root / "metadata" / "Breast_cancer_sample_metadata.feather").resolve()
    config_path = (root / "breast-datasets.json").resolve()
    columns = {
        "Sample": [
            "breast-unknown",
            "breast-unknown",
            "breast-unknown",
            "breast-unknown",
            "breast-recorded",
            "breast-recorded",
            "breast-recorded",
        ],
        "disease": [
            "disease-first",
            "disease-second",
            "disease-third",
            "disease-fourth",
            "disease-recorded-first",
            "disease-recorded-second",
            "disease-recorded-third",
        ],
        "assay": [
            "assay-first",
            "assay-second",
            "assay-third",
            "assay-fourth",
            "assay-recorded-first",
            "assay-recorded-second",
            "assay-recorded-third",
        ],
        "sequencing_platform": [
            "platform-first",
            "platform-second",
            "platform-third",
            "platform-fourth",
            "platform-recorded-first",
            "platform-recorded-second",
            "platform-recorded-third",
        ],
        "suspension_dissociation_time": [
            "unknown",
            "12h",
            "unknown",
            "24h",
            "6h",
            "12h",
            "12h",
        ],
        "broad_cell_type": [
            "broad-first",
            "broad-second",
            "broad-third",
            "broad-fourth",
            "broad-recorded-first",
            "broad-recorded-second",
            "broad-recorded-third",
        ],
        "author_cell_type": [
            "author-first",
            "author-second",
            "author-third",
            "author-fourth",
            "author-recorded-first",
            "author-recorded-second",
            "author-recorded-third",
        ],
    }
    index = [f"obs-cell-{number}" for number in range(len(columns["Sample"]))]
    _write_obs_only_h5ad(input_path, columns, index)
    with h5py.File(input_path, "r") as handle:
        assert list(handle.keys()) == ["obs"]

    config = {
        "Breast_cancer": {
            "columns": {
                "sample": "sample_id",
                "label": "disease",
                "batch": [
                    "assay",
                    "sequencing_platform",
                    "suspension_dissociation_time",
                ],
                "cell_type_low_res": "broad_cell_type",
                "cell_type_high_res": "author_cell_type",
            },
            "batch_metadata_policy": BREAST_POLICY,
            "views": {
                "batch_effect_uncorrected": {
                    "input_file_name": input_path.name,
                    "output_file_name": "BreastCncr_processed_batch.h5ad",
                }
            },
        }
    }
    config_path.write_text(json.dumps(config), encoding="utf-8")
    export(
        argparse.Namespace(
            config=config_path,
            dataset="Breast_cancer",
            view="batch_effect_uncorrected",
            input_file=input_path,
            output=output_path,
            chunk_size=1,
            check=False,
        )
    )

    metadata = pd.read_feather(output_path)
    expected_columns = [
        "Sample",
        "disease",
        "assay",
        "sequencing_platform",
        "suspension_dissociation_time",
        "broad_cell_type",
        "author_cell_type",
    ]
    assert list(metadata.columns) == expected_columns
    assert metadata["Sample"].tolist() == ["breast-unknown", "breast-recorded"]
    assert metadata["disease"].tolist() == [
        "disease-first",
        "disease-recorded-first",
    ]
    assert metadata["assay"].tolist() == ["assay-first", "assay-recorded-first"]
    assert metadata["sequencing_platform"].tolist() == [
        "platform-first",
        "platform-recorded-first",
    ]
    assert metadata["suspension_dissociation_time"].tolist() == ["unknown", "12h"]
    assert metadata["broad_cell_type"].tolist() == [
        "broad-first",
        "broad-recorded-first",
    ]
    assert metadata["author_cell_type"].tolist() == [
        "author-first",
        "author-recorded-first",
    ]
    assert not {"X", "raw", "layers", "counts"} & set(metadata.columns)

    sidecar = Path(f"{output_path}.md5")
    digest = hashlib.md5(output_path.read_bytes()).hexdigest()
    assert sidecar.read_text(encoding="utf-8").splitlines() == [
        f"MD5={digest}",
        f"SIZE={output_path.stat().st_size}",
        f"PATH={output_path}",
    ]

def _check_corrected_lupus_highres_contract(root: Path) -> None:
    """Corrected-final export requires configured high-res louvain only."""

    input_path = (root / "Lupus_Perez2022.h5ad").resolve()
    corrected_root = (root / "batch_effect" / "corrected_final").resolve()
    output_path = (
        corrected_root
        / "metadata"
        / "Lupus_PBMC_sample_metadata.feather"
    )
    config_path = (root / "lupus-datasets.json").resolve()
    columns = {
        "Sample": ["sample-a", "sample-a", "sample-b", "sample-b"],
        "sampleID": ["raw-a", "raw-a", "raw-b", "raw-b"],
        "Status": ["case", "case", "control", "control"],
        "batch_cov": ["batch-1", "batch-1", "batch-2", "batch-2"],
        "louvain": ["B", "B", "A", "A"],
    }
    _write_obs_only_h5ad(
        input_path,
        columns,
        [f"lupus-cell-{number}" for number in range(4)],
    )
    config = {
        "Lupus_PBMC": {
            "columns": {
                "sample": "sampleID",
                "label": "Status",
                "batch": "batch_cov",
                "cell_type_low_res": "layer1",
                "cell_type_high_res": "louvain",
            },
            "batch_metadata_policy": {
                "sample_aggregation": "majority_v1",
                "majority_keys": ["batch_cov"],
                "accepted_sentinel_values": {},
            },
            "views": {
                "batch_effect_corrected": {
                    "input_file_name": input_path.name,
                    "output_file_name": "Lupus_corrected.h5ad",
                }
            },
        }
    }
    config_path.write_text(json.dumps(config), encoding="utf-8")
    previous_pass = os.environ.get("ANALYSIS_PASS")
    previous_root = os.environ.get("ANALYSIS_ROOT")
    os.environ["ANALYSIS_PASS"] = "corrected"
    os.environ["ANALYSIS_ROOT"] = str(corrected_root)
    try:
        export(
            argparse.Namespace(
                config=config_path,
                dataset="Lupus_PBMC",
                view="batch_effect_corrected",
                analysis_variant="corrected_final",
                input_file=input_path,
                output=output_path,
                chunk_size=2,
                check=False,
            )
        )
    finally:
        if previous_pass is None:
            os.environ.pop("ANALYSIS_PASS", None)
        else:
            os.environ["ANALYSIS_PASS"] = previous_pass
        if previous_root is None:
            os.environ.pop("ANALYSIS_ROOT", None)
        else:
            os.environ["ANALYSIS_ROOT"] = previous_root
    metadata = pd.read_feather(output_path)
    assert metadata["Sample"].tolist() == ["sample-a", "sample-b"]
    assert metadata["louvain"].tolist() == ["B", "A"]
    assert "layer1" not in metadata.columns
    assert not {"X", "raw", "layers", "counts"} & set(metadata.columns)





def _check_configured_policies() -> None:
    config = json.loads((ROOT / "datasets.json").read_text(encoding="utf-8"))
    expected = {
        "Alzheimer": {
            "sample_aggregation": "majority_v1",
            "majority_keys": ["assay"],
            "accepted_sentinel_values": {},
        },
        "Breast_cancer": {
            "sample_aggregation": "majority_v1",
            "majority_keys": ["suspension_dissociation_time"],
            "accepted_sentinel_values": {
                "suspension_dissociation_time": ["unknown"],
            },
        },
        "Lupus_PBMC": {
            "sample_aggregation": "majority_v1",
            "majority_keys": ["batch_cov"],
            "accepted_sentinel_values": {},
        },
    }
    actual = {
        dataset: entry["batch_metadata_policy"]
        for dataset, entry in config.items()
        if "batch_metadata_policy" in entry
    }
    assert actual == expected
    assert {dataset: policy["majority_keys"] for dataset, policy in actual.items()} == {
        "Alzheimer": ["assay"],
        "Breast_cancer": ["suspension_dissociation_time"],
        "Lupus_PBMC": ["batch_cov"],
    }
    assert set(actual) == set(expected)
    for dataset, entry in config.items():
        if dataset not in expected:
            assert "batch_metadata_policy" not in entry, (
                f"unexpected majority policy configured for {dataset}"
            )
        for view_name, view in entry.get("views", {}).items():
            assert "batch_metadata_policy" not in view, (
                f"majority policy must remain dataset-scoped: {dataset}/{view_name}"
            )



def main() -> None:
    _check_configured_policies()
    with tempfile.TemporaryDirectory(prefix="ecoda-batch-majority-") as temporary:
        root = Path(temporary)
        _check_assay_majority_fixture(root)
        _check_tie_fixture(root)
        _check_invalid_assay_values(root)
        _check_breast_fixture(root)
        _check_breast_three_key_composite()
        _check_export_boundary(root)
        _check_corrected_lupus_highres_contract(root)
    print("H5AD metadata majority contracts: OK")


if __name__ == "__main__":
    main()
