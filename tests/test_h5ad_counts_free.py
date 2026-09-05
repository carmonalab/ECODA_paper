#!/usr/bin/env python3
"""Exercise the metadata/embedding-only H5AD benchmark loader."""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from src.utils.py.h5ad_counts_free import load_h5ad_counts_free


def load_python_worker():
    worker_path = (
        ROOT
        / "src/5_run_benchmark_methods/run_python_sample_embedding_methods/"
        "1.1.1_benchmark_methods_py.py"
    )
    spec = importlib.util.spec_from_file_location("ecoda_python_worker", worker_path)
    assert spec and spec.loader
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    return worker


def run_r_loader_check(path: Path) -> None:
    r_code = r"""
args <- commandArgs(trailingOnly = TRUE)
project <- args[[1]]
input <- args[[2]]
Sys.setenv(PROJECT_ROOT = project)
source(file.path(project, "src", "5_run_benchmark_methods", "benchmark_hpc_utils.R"))
adata <- load_h5ad_counts_free(
  input,
  c("Sample", "cell_type"),
  "X_pca_batch_effect_uncorrected_hvg2000",
  "leiden_res_",
  view = "batch_effect_uncorrected",
  method = "gloscope"
)
shape <- reticulate::py_to_r(adata$shape)
stopifnot(as.integer(shape[[1]]) == 3L)
stopifnot(as.integer(shape[[2]]) == 2000L)
stopifnot(length(reticulate::py_to_r(adata$layers$keys())) == 0L)
message("R counts-free H5AD loader: OK")
"""
    subprocess.run(
        [
            "pixi",
            "run",
            "Rscript",
            "--vanilla",
            "-e",
            r_code,
            str(ROOT),
            str(path),
        ],
        check=True,
        env={**os.environ, "RETICULATE_PYTHON": sys.executable},
    )


def run_python_worker_counts_free_check(raw: Path) -> None:
    worker = load_python_worker()
    path = raw / "worker_input.h5ad"
    n_genes = 2000
    obs = pd.DataFrame(
        {"Sample": ["s1", "s1", "s2"], "cell_type": ["T", "B", "T"]},
        index=["cell1", "cell2", "cell3"],
    )
    counts = sparse.csr_matrix((3, n_genes), dtype=np.int32)
    counts[0, 0] = 1
    counts[1, 1] = 2
    counts[2, 2] = 3
    data = ad.AnnData(
        X=counts.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(
            {"hvg_rank": np.arange(1, n_genes + 1, dtype=float)},
            index=[f"g{i}" for i in range(n_genes)],
        ),
    )
    data.layers["counts"] = counts
    data.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.ones((3, 2))
    data.write_h5ad(path)

    captured = {}

    def fake_run_pilot(adata, ct_col, view, n_hvg, output_path):
        captured["adata"] = adata
        frame = pd.DataFrame(
            [[0.0, 1.0], [1.0, 0.0]],
            index=["s1", "s2"],
            columns=["s1", "s2"],
        )
        worker.atomic_to_feather(frame, output_path)

    worker.run_pilot = fake_run_pilot
    worker.log_execution_time = lambda *args, **kwargs: None
    worker.peak_rss_gb = lambda: 0.0
    args = SimpleNamespace(
        view="batch_effect_uncorrected",
        analysis_pass="uncorrected",
        combo=None,
        high_resolution_only=True,
        output_dir=str(raw / "worker_output"),
        input_dir=str(raw),
        method="pilot",
        hvg=[2000],
        force=False,
        device="cpu",
        log_file=str(raw / "worker.log"),
    )
    entry = {
        "views": {"batch_effect_uncorrected": {"output_file": path.name}},
        "cell_type_low_res": "cell_type",
        "cell_type_high_res": "cell_type",
        "batch_col": None,
    }
    worker.process_dataset(args, "Synthetic", entry)
    loaded = captured["adata"]
    assert loaded.X.nnz == 0
    assert "counts" not in loaded.layers
    assert list(loaded.obs.columns) == ["Sample", "cell_type"]


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-counts-free-") as raw:
        path = Path(raw) / "input.h5ad"
        obs = pd.DataFrame(
            {
                "Sample": ["s1", "s1", "s2"],
                "cell_type": ["T", "B", "T"],
                "label": ["case", "case", "control"],
                "leiden_res_0.1_batch_effect_uncorrected_hvg2000": [
                    "0",
                    "1",
                    "0",
                ],
            },
            index=["cell1", "cell2", "cell3"],
        )
        n_genes = 2000
        counts = sparse.csr_matrix((3, n_genes), dtype=np.int32)
        counts[0, 0] = 1
        counts[0, 2] = 2
        counts[1, 1] = 3
        counts[2, 0] = 4
        counts[2, 2] = 5
        data = ad.AnnData(
            X=counts.astype(np.float32),
            obs=obs,
            var=pd.DataFrame(
                {"hvg_rank": np.arange(1, n_genes + 1, dtype=float)},
                index=[f"g{i}" for i in range(1, n_genes + 1)],
            ),
        )
        data.layers["counts"] = counts
        data.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.asarray(
            [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float32
        )
        data.write_h5ad(path)

        minimal = load_h5ad_counts_free(
            path,
            ["Sample", "cell_type"],
            ["X_pca_batch_effect_uncorrected_hvg2000"],
            ["leiden_res_"],
        )
        assert minimal.shape == data.shape
        assert minimal.X.nnz == 0
        assert "counts" not in minimal.layers
        assert "label" not in minimal.obs.columns
        assert "cell_type" in minimal.obs.columns
        assert "leiden_res_0.1_batch_effect_uncorrected_hvg2000" in minimal.obs.columns
        np.testing.assert_allclose(
            minimal.obsm["X_pca_batch_effect_uncorrected_hvg2000"],
            data.obsm["X_pca_batch_effect_uncorrected_hvg2000"],
        )
        assert list(minimal.var.index[:3]) == ["g1", "g2", "g3"]
        np.testing.assert_allclose(
            minimal.var["hvg_rank"].to_numpy()[:3], [1.0, 2.0, 3.0]
        )

        try:
            load_h5ad_counts_free(
                path,
                ["missing"],
                ["X_pca_batch_effect_uncorrected_hvg2000"],
            )
        except ValueError as exc:
            assert "requested obs columns" in str(exc)
        else:
            raise AssertionError("missing requested obs column was accepted")

        run_r_loader_check(path)
        run_python_worker_counts_free_check(Path(raw))

    print("counts-free H5AD loader: OK")


if __name__ == "__main__":
    main()
