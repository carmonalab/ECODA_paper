#!/usr/bin/env python3

import json
import numpy as np
import pandas as pd
import sys
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from anndata import AnnData

from src.utils.py import benchmark_h5ad_contract as validator
from src.utils.py.batch_contract import (
    build_batch_contract_identity,
    build_historical_batch_contract_identity,
)
from src.utils.py.benchmark_h5ad_contract import (
    REQUIRED_OBSM,
    validate_benchmark_h5ad_contract,
)


def make_valid_adata(view="benchmark_analysis"):
    n_cells = 4
    n_genes = 3000
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s1", "s2", "s2"],
            "batch": ["b1", "b1", "b2", "b2"],
        },
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


def _invoke_validator_argument(argument):
    """Exercise the CLI argument boundary while observing parsed kwargs."""
    captured = []
    original_validator = validator.validate_benchmark_h5ad_path

    def capture(path, view, method, **kwargs):
        captured.append(kwargs)

    validator.validate_benchmark_h5ad_path = capture
    original_argv = sys.argv
    sys.argv = [
        "benchmark_h5ad_contract.py",
        "--path",
        "unused.h5ad",
        "--view",
        "batch_effect_uncorrected",
        "--method",
        "benchmark",
        "--expected-batch-contract",
        argument,
    ]
    try:
        validator.main()
    finally:
        validator.validate_benchmark_h5ad_path = original_validator
        sys.argv = original_argv

    assert len(captured) == 1
    return captured[0]["expected_batch_contract"]


def assert_batch_contract_argument_boundary():
    identity = {"padding": "x" * 1024, "nested": {"value": 17}}
    inline = json.dumps(identity)
    assert len(inline) > 255
    assert _invoke_validator_argument(inline) == identity

    with tempfile.TemporaryDirectory() as directory:
        identity_path = Path(directory) / "identity.json"
        identity_path.write_text(inline, encoding="utf-8")
        assert _invoke_validator_argument(str(identity_path)) == identity

    try:
        _invoke_validator_argument('{"malformed":')
    except ValueError as exc:
        assert "JSON object" in str(exc)
    else:
        raise AssertionError("malformed batch contract JSON was accepted")

    try:
        _invoke_validator_argument("[]")
    except ValueError as exc:
        assert "must be an object" in str(exc)
    else:
        raise AssertionError("non-object batch contract JSON was accepted")


def assert_corrected_final_consumer_policy():
    identity = build_batch_contract_identity(
        ["batch"],
        sample_column="Sample",
        method_id="preprocess",
        model_id="hvg_composite_v1",
    )
    adata = make_valid_adata("batch_effect_corrected")
    adata.uns["batch_contract"] = identity
    validate_benchmark_h5ad_contract(
        adata,
        "batch_effect_corrected",
        "gloscope",
        expected_batch_contract=identity,
        batch_contract=identity,
        corrected_final_consumer=True,
    )
    try:
        validate_benchmark_h5ad_contract(
            adata,
            "batch_effect_corrected",
            "gloscope",
            expected_batch_contract=identity,
            batch_contract=identity,
            require_corrected_summary=False,
        )
    except ValueError as exc:
        assert "missing validation_summary" in str(exc), str(exc)
    else:
        raise AssertionError(
            "ordinary corrected GloScope consumer accepted a summary-free H5AD"
        )


def assert_active_corrected_model_policy():
    """Accept active limma identities and reject retired model IDs."""
    for method_id, model_id in (
        ("ECODA_authors_HR", "limma_fixed_effects_v1"),
        ("Pseudobulk", "pseudobulk_limma_fixed_effects_v1"),
    ):
        identity = build_batch_contract_identity(
            ["batch"],
            sample_column="Sample",
            method_id=method_id,
            model_id=model_id,
        )
        validator.validate_batch_contract_identity(
            batch_contract=identity,
            label=f"{method_id} active identity",
        )

    historical_identities = (
        ("ECODA_authors_HR", "ecoda_additive_random_intercepts_v1"),
        ("Pseudobulk", "pseudobulk_composite_v1"),
    )
    for method_id, model_id in historical_identities:
        identity = build_historical_batch_contract_identity(
            ["batch"],
            sample_column="Sample",
            method_id=method_id,
            model_id=model_id,
        )
        try:
            validator.validate_batch_contract_identity(
                batch_contract=identity,
                label=f"{method_id} historical identity",
            )
        except ValueError as exc:
            assert "unsupported model policy" in str(exc), str(exc)
        else:
            raise AssertionError(
                f"historical model ID {model_id!r} was accepted by the "
                "active H5AD identity validator"
            )


def main():
    assert_batch_contract_argument_boundary()
    assert_corrected_final_consumer_policy()
    assert_active_corrected_model_policy()

    for view in REQUIRED_OBSM:
        kwargs = {}
        method = "test"
        if view == "batch_effect_corrected":
            corrected_identity = build_batch_contract_identity(
                ["batch"],
                sample_column="Sample",
                method_id="preprocess",
                model_id="hvg_composite_v1",
            )
            kwargs = {
                "expected_batch_contract": corrected_identity,
                "batch_contract": corrected_identity,
                "allow_missing_corrected_summary": True,
            }
            method = "preprocessing"
        validate_benchmark_h5ad_contract(
            make_valid_adata(view), view, method, **kwargs
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
