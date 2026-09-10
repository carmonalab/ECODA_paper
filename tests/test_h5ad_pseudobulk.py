#!/usr/bin/env python3
"""Exercise bounded CSR sample and composite pseudobulk preparation."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from src.utils.py.h5ad_pseudobulk import (
    aggregate_h5ad_counts_by_sample,
    audit_h5ad_ct_group_store,
    prepare_h5ad_ct_group_store,
    read_h5ad_ct_group_store,
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


def _write_counts_fixture(path: Path, counts: sparse.csr_matrix, sample_ids: list[str]) -> None:
    n_obs, n_vars = counts.shape
    obs = pd.DataFrame(
        {"Sample": sample_ids},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    data = ad.AnnData(
        X=sparse.csr_matrix((n_obs, n_vars), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_vars)]),
    )
    data.layers["counts"] = counts
    data.write_h5ad(path)


def _write_short_obs_fixture(path: Path) -> None:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", shape=(2, 1), dtype=np.float32)
        obs = handle.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["_index"] = "index"
        obs.create_dataset(
            "index", data=np.asarray(["cell0", "cell1"], dtype=object), dtype=string_dtype
        )
        obs.create_dataset("Sample", data=np.asarray(["s1"], dtype=object), dtype=string_dtype)
        layers = handle.create_group("layers")
        counts = layers.create_group("counts")
        counts.attrs["encoding-type"] = "csr_matrix"
        counts.attrs["shape"] = (2, 1)
        counts.create_dataset("data", data=np.asarray([1, 2], dtype=np.int32))
        counts.create_dataset("indices", data=np.asarray([0, 0], dtype=np.int32))
        counts.create_dataset("indptr", data=np.asarray([0, 1, 2], dtype=np.int32))
        var = handle.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.asarray(["g1"], dtype=object), dtype=string_dtype)


def _clone_store_fixture(source: Path, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    manifest_path = Path(f"{target}.manifest.json")
    manifest = json.loads(Path(f"{source}.manifest.json").read_text())
    manifest["store_path"] = str(target.resolve())
    manifest["lock_path"] = str(Path(f"{target}.lock").resolve())
    manifest_path.write_text(json.dumps(manifest))
    return target


def _assert_corrupt_store_rejected(
    source: Path, target: Path, mutate, requested_group: str | None
) -> None:
    _clone_store_fixture(source, target)
    with h5py.File(target, "r+") as handle:
        mutate(handle)
    manifest_path = Path(f"{target}.manifest.json")
    manifest = json.loads(manifest_path.read_text())
    manifest["pid"] = 2_147_483_647
    manifest["scheduler_identity"] = {}
    manifest["created_at"] = time.time() - 3600
    manifest_path.write_text(json.dumps(manifest))

    audit = audit_h5ad_ct_group_store(
        target,
        expected_run_id="focused-run",
        max_age_seconds=1,
        cleanup=True,
    )
    assert not audit["valid"]
    assert not audit["cleanup_performed"]
    assert target.exists() and manifest_path.exists()
    if requested_group is None:
        return
    try:
        read_h5ad_ct_group_store(target, [requested_group])
    except (OSError, TypeError, ValueError):
        pass
    else:
        raise AssertionError("malformed CT store was readable")


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
                "cell_type": ["T1", "T2", "T1", "T1", "T2"],
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
        assert result["counts"].dtype == np.int64
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

        # Composite IDs, occupancy, first-observation metadata, and chunk
        # boundaries are all exercised by the same five-cell fixture.
        store = Path(raw) / "ct" / "groups.h5"
        composite = prepare_h5ad_ct_group_store(
            path,
            sample_col="Sample",
            cell_type_col="cell_type",
            metadata_columns=["batch"],
            chunk_size=2,
            max_value=2_147_483_647,
            store_path=store,
            run_id="focused-run",
            source_identity={"checksum": "fixture"},
        )
        assert composite["group_ids"] == [
            "Sample=s1;cell_type=T1",
            "Sample=s2;cell_type=T2",
            "Sample=s3;cell_type=T1",
        ]
        assert composite["sample_ids"] == ["s1", "s2", "s3"]
        assert composite["all_sample_ids"] == ["s1", "s2", "s3"]
        assert composite["cell_type_ids"] == ["T1", "T2", "T1"]
        assert composite["group_cell_counts"] == [2, 2, 1]
        assert list(composite["group_metadata"]["batch"]) == ["A", "B", "A"]
        selected = read_h5ad_ct_group_store(
            store,
            [composite["group_ids"][0], composite["group_ids"][1]],
        )
        assert sparse.isspmatrix_csr(selected["counts"])
        np.testing.assert_array_equal(
            selected["counts"].toarray(),
            np.asarray([[5, 0, 7, 0], [1, 4, 1, 2]], dtype=np.int64),
        )
        assert selected["group_ids"] == composite["group_ids"][:2]
        audit = audit_h5ad_ct_group_store(
            store, expected_run_id="focused-run", max_age_seconds=3600
        )
        assert audit["layout"]["max_value"] == np.iinfo(np.int32).max
        assert audit["valid"] and not audit["cleanup_performed"]
        assert "toarray()" not in (
            ROOT / "src/utils/py/h5ad_pseudobulk.py"
        ).read_text()

        # Every persisted CSR vector and aggregate row is audited before a
        # store can be read or considered eligible for stale cleanup.
        corruption_cases = [
            (
                "bad-contribution-indptr",
                lambda handle: handle["indptr"].__setitem__(0, 1),
            ),
            (
                "bad-contribution-index",
                lambda handle: handle["indices"].__setitem__(
                    0, int(handle.attrs["n_vars"])
                ),
            ),
            (
                "out-of-bounds-aggregate-index",
                lambda handle: handle["group_aggregates"]["0"]["indices"].__setitem__(
                    0, int(handle.attrs["n_vars"])
                ),
            ),
            (
                "duplicate-aggregate-index",
                lambda handle: handle["group_aggregates"]["0"]["indices"].__setitem__(
                    1, handle["group_aggregates"]["0"]["indices"][0]
                ),
            ),
            (
                "negative-aggregate-data",
                lambda handle: handle["group_aggregates"]["0"]["data"].__setitem__(0, -1),
            ),
            (
                "mismatched-group-count",
                lambda handle: handle["group_cell_counts"].__setitem__(
                    0, handle["group_cell_counts"][0] + 1
                ),
            ),
        ]
        for name, mutate in corruption_cases:
            _assert_corrupt_store_rejected(
                store,
                Path(raw) / "ct-corrupt" / f"{name}.h5",
                mutate,
                None
                if name.startswith("bad-contribution") or name == "mismatched-group-count"
                else composite["group_ids"][0],
            )

        unselected_store = Path(raw) / "ct-corrupt" / "unselected-row.h5"
        _clone_store_fixture(store, unselected_store)
        with h5py.File(unselected_store, "r+") as handle:
            handle["group_aggregates"]["1"]["data"][0] = -1
        unselected_audit = audit_h5ad_ct_group_store(
            unselected_store,
            expected_run_id="focused-run",
            max_age_seconds=1,
            cleanup=True,
        )
        assert not unselected_audit["valid"]
        assert not unselected_audit["cleanup_performed"]
        selected_from_unselected = read_h5ad_ct_group_store(
            unselected_store, [composite["group_ids"][0]]
        )
        assert sparse.isspmatrix_csr(selected_from_unselected["counts"])
        assert selected_from_unselected["group_ids"] == [composite["group_ids"][0]]

        # A final-component symlink is rejected before canonicalisation by
        # audit, reader, and writer; the external store remains untouched.
        symlink_store = Path(raw) / "ct-link" / "groups.h5"
        symlink_store.parent.mkdir(parents=True, exist_ok=True)
        symlink_store.symlink_to(store)
        symlink_audit = audit_h5ad_ct_group_store(
            symlink_store,
            expected_run_id="focused-run",
            max_age_seconds=0,
            cleanup=True,
        )
        assert not symlink_audit["valid"]
        assert not symlink_audit["cleanup_performed"]
        assert store.exists() and symlink_store.is_symlink()
        try:
            read_h5ad_ct_group_store(symlink_store, [composite["group_ids"][0]])
        except ValueError:
            pass
        else:
            raise AssertionError("symlinked CT store was readable")
        try:
            prepare_h5ad_ct_group_store(
                path,
                sample_col="Sample",
                cell_type_col="cell_type",
                metadata_columns=["batch"],
                chunk_size=2,
                max_value=2_147_483_647,
                store_path=symlink_store,
                run_id="focused-run",
                source_identity={"checksum": "fixture"},
            )
        except ValueError:
            pass
        else:
            raise AssertionError("symlinked CT store was writable")
        assert store.exists() and symlink_store.is_symlink()
        symlink_store.unlink()

        # A stale but structurally valid store is removable only after the
        # owner is demonstrably dead and the age policy has expired.
        manifest_path = Path(f"{store}.manifest.json")
        manifest = json.loads(manifest_path.read_text())
        manifest["pid"] = 2_147_483_647
        manifest["scheduler_identity"] = {}
        manifest["created_at"] = time.time() - 3600
        manifest_path.write_text(json.dumps(manifest))
        stale = audit_h5ad_ct_group_store(
            store,
            expected_run_id="focused-run",
            max_age_seconds=1,
            cleanup=True,
        )
        assert stale["cleanup_performed"]
        assert not store.exists() and not manifest_path.exists()

        # Boundary regressions: signed negatives, cumulative int64 overflow,
        # and the DESeq2 INT_MAX transport ceiling all fail before return.
        negative_path = Path(raw) / "negative.h5ad"
        _write_counts_fixture(
            negative_path,
            sparse.csr_matrix(np.asarray([[-1]], dtype=np.int64)),
            ["s1"],
        )
        try:
            aggregate_h5ad_counts_by_sample(negative_path)
        except ValueError as exc:
            assert "negative" in str(exc).lower()
        else:
            raise AssertionError("negative signed counts were accepted")

        overflow_path = Path(raw) / "overflow.h5ad"
        _write_counts_fixture(
            overflow_path,
            sparse.csr_matrix(
                np.asarray([[np.iinfo(np.int64).max], [1]], dtype=np.int64)
            ),
            ["s1", "s1"],
        )
        try:
            aggregate_h5ad_counts_by_sample(overflow_path)
        except OverflowError as exc:
            assert "int64" in str(exc)
        else:
            raise AssertionError("cumulative int64 overflow was accepted")

        int_max_path = Path(raw) / "int-max.h5ad"
        _write_counts_fixture(
            int_max_path,
            sparse.csr_matrix(
                np.asarray([[np.iinfo(np.int32).max + 1]], dtype=np.int64)
            ),
            ["s1"],
        )
        try:
            aggregate_h5ad_counts_by_sample(
                int_max_path, max_value=np.iinfo(np.int32).max
            )
        except ValueError as exc:
            assert "max_value" in str(exc)
        else:
            raise AssertionError("INT_MAX aggregate ceiling was accepted")

        float_ceiling_path = Path(raw) / "float-ceiling.h5ad"
        float_ceiling = (1 << 53) + 3
        _write_counts_fixture(
            float_ceiling_path,
            sparse.csr_matrix(np.asarray([[float((1 << 53) + 4)]], dtype=np.float64)),
            ["s1"],
        )
        try:
            aggregate_h5ad_counts_by_sample(
                float_ceiling_path, max_value=float_ceiling
            )
        except ValueError as exc:
            assert "max_value" in str(exc)
        else:
            raise AssertionError("rounded float max_value ceiling was accepted")

        short_metadata_path = Path(raw) / "short-metadata.h5ad"
        _write_short_obs_fixture(short_metadata_path)
        try:
            aggregate_h5ad_counts_by_sample(short_metadata_path)
        except ValueError as exc:
            assert "wrong length" in str(exc)
        else:
            raise AssertionError("short obs metadata was accepted")

        try:
            aggregate_h5ad_counts_by_sample(path, metadata_columns=["missing"], chunk_size=2)
        except ValueError as exc:
            assert "missing requested obs columns" in str(exc)
        else:
            raise AssertionError("missing metadata column was accepted")

        run_r_helper_check(path)
        run_prepare_worker_check(Path(raw))

    print("streaming H5AD pseudobulk: OK")


if __name__ == "__main__":
    main()
