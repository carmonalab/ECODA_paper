#!/usr/bin/env python3
"""Standalone deterministic regression: python tests/test_extract_stephenson_cd14_pseudobulk.py."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
FILENAME = "StephensonE_2021_33879890_preprocessed_benchmark_analysis_ECODAprocessed.h5ad"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(path, variant="valid"):
    # Interleave samples; only the first five cells per sample are low-res CD14.
    samples = [sample for _ in range(6) for sample in ("s2", "s1")]
    annotations = ["CD14"] * 10 + ["other", "CD14_mono"]
    labels = ["Covid" if sample == "s2" else "Healthy" for sample in samples]
    if variant == "absent":
        annotations = ["other" if sample == "s1" else ct
                       for sample, ct in zip(samples, annotations)]
    if variant == "below_min":
        annotations[0] = "other"
    if variant == "outside_status_conflict":
        labels[10] = "Healthy"
    if variant == "missing_sample":
        samples[10] = ""
    if variant == "wrong_label":
        labels[10] = "unknown"
    if variant == "wrong_annotation":
        annotations = ["CD14_mono"] * 12
    rows = [[(0, 2), (0, 3), (2, row + 1)] for row in range(12)]
    # Deliberately invalid nonselected counts must neither enter sums nor trigger validation.
    rows[10] = [(1, -0.5)]
    rows[11] = [(1, 9000)]
    if variant == "negative":
        rows[0] = [(0, -1)]
    if variant == "fractional":
        rows[0] = [(0, 1.5)]
    if variant == "r_overflow":
        rows[0] = [(0, 2**31)]
    if variant == "duplicate_overflow":
        rows[0] = [(0, 2**63 - 1), (0, 1)]
    if variant == "sum_overflow":
        rows[0] = [(0, 2**63 - 1)]
        rows[2] = [(0, 1)]
    integer_counts = variant in {"r_overflow", "duplicate_overflow", "sum_overflow"}
    if integer_counts:
        rows[10] = [(1, 0)]
    values, indices, indptr = [], [], [0]
    for row in rows:
        indices.extend(index for index, _ in row)
        values.extend(value for _, value in row)
        indptr.append(len(values))
    values = np.asarray(values, dtype=np.int64 if integer_counts else np.float64)
    indices, indptr = np.asarray(indices), np.asarray(indptr)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=np.full((12, 3), 0.125))
        obs = handle.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["_index"] = "_index"
        for name, column in {"_index": [f"cell{i}" for i in range(12)], "Sample": samples,
                             "Status": labels, "initial_clustering": annotations,
                             "full_clustering": ["CD14_mono"] * 12}.items():
            node = obs.create_dataset(name, data=np.asarray(column, dtype=object),
                                      dtype=h5py.string_dtype())
            node.attrs["encoding-type"] = "string-array"
        if variant == "wrong_column":
            obs.move("initial_clustering", "wrong_column")
        var = handle.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.asarray(["g1", "g2", "g3"], dtype=object),
                           dtype=h5py.string_dtype())
        counts = handle.create_group("layers").create_group("counts")
        counts.attrs["encoding-type"] = "csr_matrix"
        counts.attrs["shape"] = (12, 3)
        for name, value in (("data", values), ("indices", indices), ("indptr", indptr)):
            counts.create_dataset(name, data=value)
        if variant == "no_counts":
            del handle["layers/counts"]
    reference = sparse.csr_matrix((values, indices, indptr), shape=(12, 3))
    return reference, samples, annotations


def test_cd14_extraction():
    with tempfile.TemporaryDirectory() as temporary:
        base = Path(temporary)
        isolated = base / "isolated"
        isolated.mkdir()
        script = isolated / "extract.py"
        shutil.copyfile(ROOT / "src/extract_stephenson_cd14_pseudobulk.py", script)
        frozen = base / "frozen"
        shutil.copytree(ROOT / "src/utils/py", frozen / "src/utils/py",
                        ignore=shutil.ignore_patterns("__pycache__"))
        config = base / "datasets.json"
        shutil.copyfile(ROOT / "datasets.json", config)

        def run(source, output, success):
            before, stat = digest(source), source.stat()
            result = subprocess.run(
                [sys.executable, str(script), "--h5ad", str(source), "--output-dir", str(output),
                 "--config", str(config), "--helpers-root", str(frozen)],
                cwd=isolated, capture_output=True, text=True,
            )
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            assert digest(source) == before
            after = source.stat()
            assert (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns) == (
                after.st_size, after.st_mtime_ns, after.st_ctime_ns)

        source = base / FILENAME
        matrix, samples, annotations = fixture(source)
        output = base / "aggregated"
        run(source, output, True)
        raw = pd.read_feather(output / "counts.feather")
        metadata = pd.read_feather(output / "samples.feather")
        assert list(raw.columns) == ["gene", "s1", "s2"]
        assert list(raw["gene"]) == ["g1", "g2", "g3"]
        for sample in ("s1", "s2"):
            selected = [i for i in range(12) if samples[i] == sample and annotations[i] == "CD14"]
            np.testing.assert_array_equal(raw[sample], np.asarray(matrix[selected].sum(axis=0)).ravel())
            assert raw[sample].dtype == np.dtype("int32")
        assert metadata.to_dict("list") == {
            "Sample": ["s1", "s2"], "Status": ["Healthy", "Covid"], "n_cells": [5, 5]}
        assert metadata.n_cells.dtype == np.dtype("int32")
        manifest = json.loads((output / "aggregation.json").read_text())
        assert manifest["n_selected_cells"] == 10
        assert manifest["cell_type"] == "CD14"
        assert manifest["annotation_column"] == "initial_clustering"
        assert manifest["source"]["full_file_checksum_verified"] is False
        assert manifest["total_umis"] == 105
        for name, expected in manifest["files"].items():
            assert digest(output / name) == expected
        saved = {path.name: digest(path) for path in output.iterdir()}
        run(source, output, False)
        assert saved == {path.name: digest(path) for path in output.iterdir()}
        wrong_name = base / "wrong.h5ad"
        shutil.copyfile(source, wrong_name)
        run(wrong_name, base / "wrong_name_output", False)
        for variant in ("absent", "below_min", "outside_status_conflict", "missing_sample",
                        "wrong_label", "wrong_annotation", "wrong_column", "no_counts",
                        "negative", "fractional", "r_overflow", "duplicate_overflow", "sum_overflow"):
            fixture(source, variant)
            failed_output = base / variant
            run(source, failed_output, False)
            assert not (failed_output / "aggregation.json").exists()
        fixture(source)
        bad_config = json.loads(config.read_text())
        bad_config["Stephenson"]["columns"]["cell_type_low_res"] = "full_clustering"
        config.write_text(json.dumps(bad_config))
        run(source, base / "bad_config", False)


if __name__ == "__main__":
    test_cd14_extraction()
    print("CD14 extraction regression: OK")
