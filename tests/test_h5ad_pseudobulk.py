#!/usr/bin/env python3
"""Exercise bounded CSR sample aggregation for pseudobulk preparation."""
from __future__ import annotations

import os
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from src.utils.py.h5ad_pseudobulk import (
    aggregate_h5ad_counts_by_sample,
    read_h5ad_sample_metadata,
)


def run_r_helper_check(path: Path) -> None:
    r_code = r"""
args <- commandArgs(trailingOnly = TRUE)
project <- args[[1]]
input <- args[[2]]
Sys.setenv(PROJECT_ROOT = project)
source(file.path(project, "src", "5_run_benchmark_methods", "benchmark_hpc_utils.R"))
metadata <- load_h5ad_sample_metadata(
  input,
  sample_col = "Sample",
  metadata_columns = c("batch", "label"),
  chunk_size = 2L
)
stopifnot(identical(rownames(metadata), c("s1", "s2", "s3")))
stopifnot(identical(as.character(metadata$label), c("case", "control", "control")))
seurat <- load_h5ad_pseudobulk_seurat(
  input,
  sample_col = "Sample",
  batch_col = "batch",
  chunk_size = 2L
)
stopifnot(nrow(seurat) == 4L)
stopifnot(ncol(seurat) == 3L)
stopifnot(identical(colnames(seurat), c("s1", "s2", "s3")))
stopifnot(identical(as.character(seurat$batch), c("A", "B", "A")))
message("R streaming pseudobulk helper: OK")
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


def run_prepare_worker_check(raw: Path) -> None:
    path = raw / "worker_input.h5ad"
    counts = sparse.csr_matrix(
        np.asarray(
            [
                [5, 1, 2, 0, 1] + [0] * 1995,
                [2, 3, 1, 4, 0] + [0] * 1995,
                [7, 0, 4, 1, 2] + [0] * 1995,
                [1, 2, 3, 5, 1] + [0] * 1995,
                [4, 1, 0, 2, 6] + [0] * 1995,
            ],
            dtype=np.int32,
        )
    )
    obs = pd.DataFrame(
        {
            "Sample": ["s1", "s2", "s1", "s3", "s2"],
            "label": ["case", "control", "case", "control", "control"],
        },
        index=[f"cell{i}" for i in range(1, 6)],
    )
    data = ad.AnnData(
        X=counts.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(
            {"hvg_rank": np.arange(1, 2001, dtype=float)},
            index=[f"g{i}" for i in range(1, 2001)],
        ),
    )
    data.layers["counts"] = counts
    data.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.ones((5, 2))
    data.write_h5ad(path)

    config_path = raw / "synthetic_datasets.json"
    config_path.write_text(
        json.dumps(
            {
                "Synthetic": {
                    "columns": {"sample": "Sample", "label": "label"},
                    "views": {
                        "batch_effect_uncorrected": {
                            "output_file_name": path.name
                        }
                    },
                }
            }
        )
    )
    pseudobulk_dir = raw / "pseudobulks"
    log_file = raw / "prepare.log"
    worker = (
        ROOT
        / "src/5_run_benchmark_methods/run_r_sample_embedding_methods/"
        "1.1.1_prepare_pseudobulk.R"
    )
    subprocess.run(
        [
            "pixi",
            "run",
            "Rscript",
            "--vanilla",
            str(worker),
            "--config_path",
            str(config_path),
            "--ds_name",
            "Synthetic",
            "--view",
            "batch_effect_uncorrected",
            "--input_dir",
            str(raw),
            "--pseudobulk_dir",
            str(pseudobulk_dir),
            "--log_file",
            str(log_file),
            "--analysis_pass",
            "uncorrected",
        ],
        check=True,
        cwd=ROOT,
        env={**os.environ, "PROJECT_ROOT": str(ROOT), "RETICULATE_PYTHON": sys.executable},
    )
    output = pseudobulk_dir / "Synthetic_batch_effect_uncorrected_pseudobulk_hvg2000.rds"
    assert output.is_file() and output.stat().st_size > 0
    assert Path(f"{output}.md5").is_file()


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-pseudobulk-") as raw:
        path = Path(raw) / "input.h5ad"
        counts = sparse.csr_matrix(
            np.asarray(
                [
                    [1, 0, 2, 0],
                    [0, 3, 0, 1],
                    [4, 0, 5, 0],
                    [2, 2, 0, 0],
                    [1, 1, 1, 1],
                ],
                dtype=np.int32,
            )
        )
        obs = pd.DataFrame(
            {
                "Sample": ["s1", "s2", "s1", "s3", "s2"],
                "batch": ["A", "B", "A", "A", "B"],
                "label": ["case", "control", "case", "control", "control"],
            },
            index=[f"cell{i}" for i in range(1, 6)],
        )
        data = ad.AnnData(
            X=counts.astype(np.float32),
            obs=obs,
            var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
        )
        data.layers["counts"] = counts
        data.write_h5ad(path)

        result = aggregate_h5ad_counts_by_sample(
            path,
            sample_col="Sample",
            metadata_columns=["batch"],
            chunk_size=2,
        )
        assert result["sample_ids"] == ["s1", "s2", "s3"]
        np.testing.assert_array_equal(
            result["counts"],
            np.asarray(
                [
                    [5, 1, 2],
                    [0, 4, 2],
                    [7, 1, 0],
                    [0, 2, 0],
                ],
                dtype=np.int64,
            ),
        )
        assert list(result["gene_names"]) == ["g1", "g2", "g3", "g4"]
        assert list(result["metadata"]["batch"]) == ["A", "B", "A"]

        metadata = read_h5ad_sample_metadata(
            path,
            sample_col="Sample",
            metadata_columns=["batch", "label"],
            chunk_size=2,
        )
        assert list(metadata.index) == ["s1", "s2", "s3"]
        assert list(metadata["label"]) == ["case", "control", "control"]

        try:
            aggregate_h5ad_counts_by_sample(
                path,
                sample_col="Sample",
                metadata_columns=["missing"],
                chunk_size=2,
            )
        except ValueError as exc:
            assert "missing requested obs columns" in str(exc)
        else:
            raise AssertionError("missing metadata column was accepted")

        run_r_helper_check(path)
        run_prepare_worker_check(Path(raw))

    print("streaming H5AD pseudobulk: OK")


if __name__ == "__main__":
    main()
