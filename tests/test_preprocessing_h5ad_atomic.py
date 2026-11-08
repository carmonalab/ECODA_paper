#!/usr/bin/env python3
"""Preserve prior artifacts and roll back only a failed new publication."""
from __future__ import annotations

import hashlib
import importlib.util
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "ecoda_preprocess", ROOT / "src/3_scrnaseq_preprocessing/1.1.1_preprocess.py"
)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def make_adata():
    data = ad.AnnData(
        X=np.ones((2, 3000), dtype="float32"),
        obs=pd.DataFrame({"Sample": ["s1", "s2"]}, index=["c1", "c2"]),
        var=pd.DataFrame(
            {"hvg_rank": np.arange(3000, dtype=float)},
            index=[f"g{i}" for i in range(3000)],
        ),
    )
    data.layers["counts"] = np.ones((2, 3000), dtype="float32")
    for key in module.validate_benchmark_h5ad_contract.__globals__["REQUIRED_OBSM"][
        "benchmark_analysis"
    ]:
        data.obsm[key] = np.ones((2, 2), dtype="float32")
    return data


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-preprocess-atomic-") as raw:
        path = Path(raw) / "artifact.h5ad"
        module._write_h5ad_atomic(make_adata(), path, "benchmark_analysis")
        old_bytes = path.read_bytes()
        old_sidecar = Path(f"{path}.md5").read_bytes()
        try:
            module._write_h5ad_atomic(make_adata(), path, "benchmark_analysis")
        except FileExistsError:
            pass
        else:
            raise AssertionError("writer replaced an existing scientific artifact")
        assert path.read_bytes() == old_bytes
        assert Path(f"{path}.md5").read_bytes() == old_sidecar
        validator = module.validate_benchmark_h5ad_path
        calls = {"count": 0}

        def fail_installed(candidate, view, method):
            calls["count"] += 1
            if calls["count"] == 2:
                raise ValueError("installed validation failure")
            return validator(candidate, view, method)

        module.validate_benchmark_h5ad_path = fail_installed
        failed_path = path.with_name("failed-new-publication.h5ad")
        try:
            module._write_h5ad_atomic(make_adata(), failed_path, "benchmark_analysis")
        except ValueError as exc:
            assert str(exc) == "installed validation failure"
        else:
            raise AssertionError("installed-path validation failure was ignored")
        assert not failed_path.exists()
        assert not Path(f"{failed_path}.md5").exists()
        assert path.read_bytes() == old_bytes
        assert Path(f"{path}.md5").read_bytes() == old_sidecar
        assert hashlib.md5(path.read_bytes()).hexdigest() == module._checksum(path)
    print("atomic preprocessing h5ad: OK")


if __name__ == "__main__":
    main()
