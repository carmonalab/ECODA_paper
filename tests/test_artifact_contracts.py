#!/usr/bin/env python3
"""Standalone checks for fail-closed annotation and h5ad contracts."""
from __future__ import annotations

import hashlib
import importlib.util
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.feather as feather
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils.py.annotation_contract import (  # noqa: E402
    REQUIRED_ANNOTATION_COLUMNS,
    validate_feather,
    validate_h5ad,
)

merge_spec = importlib.util.spec_from_file_location(
    "ecoda_merge_annotations", ROOT / "src/4_cell_type_annotation/3.1_merge_annotations.py"
)
assert merge_spec.loader is not None
merge_module = importlib.util.module_from_spec(merge_spec)
merge_spec.loader.exec_module(merge_module)


def write_sidecar(path: Path) -> None:
    write = f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n"
    path.with_name(f"{path.name}.md5").write_text(write)


def annotation_frame(samples: list[str], barcodes: list[str], *, all_na_layer1=False) -> pd.DataFrame:
    labels = {
        "layer1": [None if all_na_layer1 else "T"] * len(samples),
        "layer2": ["T"] * len(samples),
        "layer3": ["T"] * len(samples),
        "layer_1": ["1"] * len(samples),
        "layer_2": ["2"] * len(samples),
        "layer_3": ["3"] * len(samples),
        "layer_4": ["4"] * len(samples),
        "layer_5": ["5"] * len(samples),
        "layer_6": ["6"] * len(samples),
        "scATOMIC_pred": ["T"] * len(samples),
        "Phase": ["G1"] * len(samples),
    }
    numeric = {
        "classification_confidence": [0.9] * len(samples),
        "S.Score": [0.1] * len(samples),
        "G2M.Score": [0.2] * len(samples),
    }
    return pd.DataFrame({"Sample": samples, "cell_barcode": barcodes, **labels, **numeric})


def make_target(path: Path) -> None:
    obs = pd.DataFrame({"Sample": ["s1", "s1"]}, index=["c1", "c2"])
    var = pd.DataFrame(index=["g1", "g2"])
    matrix = sparse.csr_matrix(np.ones((2, 2), dtype="float32"))
    target = ad.AnnData(X=matrix, obs=obs, var=var)
    target.layers["counts"] = matrix.copy()
    target.write_h5ad(path)


def expect_failure(fn, text: str) -> None:
    try:
        fn()
    except ValueError as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError(f"expected failure containing {text!r}")


def main() -> None:
    assert len(REQUIRED_ANNOTATION_COLUMNS) == 14
    with tempfile.TemporaryDirectory(prefix="ecoda-contract-") as raw:
        root = Path(raw)
        feather_path = root / "annotations_chunk_1.feather"
        frame = annotation_frame(["s1", "s1"], ["c1", "c2"])
        feather.write_feather(frame, feather_path)
        write_sidecar(feather_path)
        stats = validate_feather(feather_path, expected_keys={("s1", "c1"), ("s1", "c2")})
        validated_stats = validate_feather(
            feather_path,
            expected_keys={("s1", "c1"), ("s1", "c2")},
            sidecar_validated=True,
        )
        assert stats["n_rows"] == 2
        assert validated_stats["n_rows"] == 2
        assert stats["column_coverage"]["layer1"]["n_nonblank"] == 2
        categorical = frame.copy()
        categorical["classification_confidence"] = ["confident", "low_confidence"]
        categorical_path = root / "categorical_confidence.feather"
        feather.write_feather(categorical, categorical_path)
        categorical_stats = validate_feather(categorical_path)
        assert categorical_stats["n_rows"] == 2

        duplicate = root / "duplicate.feather"
        feather.write_feather(annotation_frame(["s1", "s1"], ["c1", "c1"]), duplicate)
        expect_failure(lambda: validate_feather(duplicate), "duplicate")

        incomplete = root / "incomplete.feather"
        incomplete_frame = frame.drop(columns=["scATOMIC_pred"])
        feather.write_feather(incomplete_frame, incomplete)
        expect_failure(lambda: validate_feather(incomplete), "required dual-method")

        target = root / "target.h5ad"
        make_target(target)
        write_sidecar(target)
        annot_dir = root / "annotations"
        annot_dir.mkdir()
        source_feather = annot_dir / "annotations_chunk_1.feather"
        feather.write_feather(frame, source_feather)
        write_sidecar(source_feather)
        merged = root / "merged.h5ad"
        merge_module.merge_annotations(str(target), str(annot_dir), str(merged))
        merged_stats = validate_h5ad(merged, expected_keys={("s1", "c1"), ("s1", "c2")})
        assert merged_stats["anchors"]["layer1"] == 2
        assert merged.with_name(f"{merged.name}.md5").is_file()
        merged_data = ad.read_h5ad(merged)
        assert merged_data.obs["classification_confidence"].dtype.kind == "f"

        isolated_frame = frame.copy()
        isolated_frame["layer1"] = pd.Series(["T", None], dtype="string")
        isolated_frame["scATOMIC_pred"] = pd.Series([None, "T"], dtype="string")
        isolated_dir = root / "isolated_annotations"
        isolated_dir.mkdir()
        isolated_file = isolated_dir / "annotations_chunk_1.feather"
        feather.write_feather(isolated_frame, isolated_file)
        write_sidecar(isolated_file)
        isolated_merged = root / "isolated_merged.h5ad"
        merge_module.merge_annotations(str(target), str(isolated_dir), str(isolated_merged))
        isolated_stats = validate_h5ad(isolated_merged)
        assert isolated_stats["anchors"]["layer1"] == 1
        assert isolated_stats["anchors"]["scATOMIC_pred"] == 1

        all_na = root / "all_na.h5ad"
        all_na_dir = root / "all_na_annotations"
        all_na_dir.mkdir()
        all_na_frame = annotation_frame(["s1", "s1"], ["c1", "c2"], all_na_layer1=True)
        all_na_file = all_na_dir / "annotations_chunk_1.feather"
        feather.write_feather(all_na_frame, all_na_file)
        write_sidecar(all_na_file)
        all_na_merged = root / "all_na_merged.h5ad"
        expect_failure(
            lambda: merge_module.merge_annotations(str(target), str(all_na_dir), str(all_na_merged)),
            "dataset anchor",
        )

        extra_dir = root / "extra_annotations"
        extra_dir.mkdir()
        extra_frame = annotation_frame(["s1", "s1", "s2"], ["c1", "c2", "cx"])
        extra_file = extra_dir / "annotations_chunk_1.feather"
        feather.write_feather(extra_frame, extra_file)
        write_sidecar(extra_file)
        expect_failure(
            lambda: merge_module.merge_annotations(str(target), str(extra_dir), str(root / "extra.h5ad")),
            "coverage mismatch",
        )
    print("artifact contracts: OK")


if __name__ == "__main__":
    main()
