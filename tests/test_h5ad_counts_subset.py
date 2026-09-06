#!/usr/bin/env python3
"""Exercise selected-HVG raw-count loading for count-dependent methods."""
from __future__ import annotations

import importlib.util
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

from src.utils.py.h5ad_counts_subset import (
    load_h5ad_counts_subset,
    read_h5ad_hvg_genes,
)


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


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-count-subset-") as raw:
        path = Path(raw) / "input.h5ad"
        counts = sparse.csr_matrix(
            np.asarray(
                [
                    [1, 0, 2, 0, 3, 0],
                    [0, 4, 0, 5, 0, 6],
                    [7, 0, 8, 0, 9, 0],
                ],
                dtype=np.int32,
            )
        )
        obs = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s1"],
                "cell_type": ["T", "B", "T"],
            },
            index=["cell1", "cell2", "cell3"],
        )
        data = ad.AnnData(
            X=counts.astype(np.float32),
            obs=obs,
            var=pd.DataFrame(
                {"hvg_rank": [6.0, 1.0, 5.0, 2.0, 4.0, 3.0]},
                index=["g1", "g2", "g3", "g4", "g5", "g6"],
            ),
        )
        data.layers["counts"] = counts
        data.write_h5ad(path)

        selected = read_h5ad_hvg_genes(path, 3)
        assert selected == ["g2", "g4", "g6"]
        loaded = load_h5ad_counts_subset(
            path,
            selected,
            obs_columns=["Sample", "cell_type"],
            chunk_size=2,
        )
        assert loaded.shape == (3, 3)
        assert loaded.X.nnz == 3
        np.testing.assert_array_equal(
            loaded.layers["counts"].toarray(),
            np.asarray([[0, 0, 0], [4, 5, 6], [0, 0, 0]], dtype=np.int64),
        )
        assert list(loaded.obs.columns) == ["Sample", "cell_type"]
        assert list(loaded.var.index) == selected
        assert tuple(loaded.uns["_ecoda_source_shape"]) == (3, 6)

        worker_path = Path(raw) / "worker_input.h5ad"
        worker_counts = sparse.csr_matrix((3, 2000), dtype=np.int32)
        worker_counts[0, 1] = 2
        worker_counts[1, 3] = 4
        worker_counts[2, 1] = 5
        worker_obs = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s1"],
                "cell_type": ["T", "B", "T"],
            },
            index=["cell1", "cell2", "cell3"],
        )
        worker_data = ad.AnnData(
            X=worker_counts.astype(np.float32),
            obs=worker_obs,
            var=pd.DataFrame(
                {"hvg_rank": np.arange(1, 2001, dtype=float)},
                index=[f"g{i}" for i in range(1, 2001)],
            ),
        )
        worker_data.layers["counts"] = worker_counts
        worker_data.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.ones(
            (3, 2)
        )
        worker_data.write_h5ad(worker_path)

        worker = load_python_worker()
        captured = {}

        def fake_run_mrvi(adata, device, output_path, batch_key=None):
            captured["adata"] = adata
            samples = list(dict.fromkeys(adata.obs["Sample"].astype(str)))
            frame = pd.DataFrame(
                np.eye(len(samples)), index=samples, columns=samples
            )
            worker.atomic_to_feather(frame, output_path)

        worker.run_mrvi = fake_run_mrvi
        args = SimpleNamespace(
            view="batch_effect_uncorrected",
            analysis_pass="uncorrected",
            combo=None,
            high_resolution_only=True,
            output_dir=str(Path(raw) / "worker_output"),
            input_dir=str(Path(raw)),
            method="mrvi",
            hvg=[3],
            force=False,
            device="cpu",
            log_file=str(Path(raw) / "worker.log"),
        )
        entry = {
            "views": {
                "batch_effect_uncorrected": {"output_file": worker_path.name}
            },
            "cell_type_low_res": "cell_type",
            "cell_type_high_res": "cell_type",
            "batch_col": None,
        }
        worker.process_dataset(args, "Synthetic", entry)
        assert captured["adata"].shape == (3, 3)
        assert "counts" in captured["adata"].layers
        assert captured["adata"].X.nnz == 2

    print("selected-HVG H5AD counts: OK")


if __name__ == "__main__":
    main()
