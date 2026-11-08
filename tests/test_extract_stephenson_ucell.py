#!/usr/bin/env python3
"""Run with: pixi run -e default python tests/test_extract_stephenson_ucell.py"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from anndata.io import write_elem
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import extract_stephenson_ucell as extractor

CONFIG = ROOT / "datasets.json"
BASENAME = json.loads(CONFIG.read_text())["Stephenson"]["views"]["benchmark_analysis"]["output_file_name"]


def fixture() -> pd.DataFrame:
    frame = pd.DataFrame({
        "Sample": pd.Categorical(["S2", "S1", "S1", "S2", "S3", "S3", "S4", "S4"],
                                 categories=["S4", "S3", "S2", "S1", "unused"]),
        "Status": pd.Categorical(["Covid", "Healthy", "Healthy", "Covid",
                                  "Healthy", "Healthy", "Covid", "Covid"]),
    }, index=[f"cell_{index}" for index in range(8)])
    for column in ("initial_clustering", "full_clustering", "layer1"):
        frame[column] = pd.Categorical(["A", "A", "A", "B", None, "", "A", "A"],
                                       categories=["B", "A", "", "unused"])
    frame["layer2"] = pd.Categorical(["A", "A", "A", "B", "B", "B", "A", "A"])
    for score in extractor.SCORES:
        frame[score] = [0.0, 2.0, 4.0, np.nan, 8.0, 10.0, np.nan, np.nan]
    frame.loc["cell_2", "HeatShock_UCell"] = np.nan
    return frame


def write_source(path: Path, obs: pd.DataFrame) -> None:
    # Intentionally not a complete AnnData: there is no X, var, raw, layers or obsm.
    with h5py.File(path, "w") as handle:
        write_elem(handle, "obs", obs)


class ExtractorTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name).resolve()
        self.source = self.directory / BASENAME
        self.output = self.directory / "export"
        write_source(self.source, fixture())

    def test_metadata_only_artifacts_and_cli(self):
        source_bytes = self.source.read_bytes()
        Path(f"{self.source}.md5").write_text("MD5=unverified-test-evidence\n")
        real_file, real_read = h5py.File, extractor.read_elem

        def readonly_file(path, mode="r", *args, **kwargs):
            self.assertEqual(mode, "r")
            return real_file(path, mode, *args, **kwargs)

        def obs_child_only(node):
            self.assertTrue(node.name.startswith("/obs/"), node.name)
            return real_read(node)

        with patch.object(extractor, "h5py", SimpleNamespace(File=readonly_file)), \
             patch.object(extractor, "read_elem", side_effect=obs_child_only), \
             patch("anndata.read_h5ad", side_effect=AssertionError("No full/backed reads")):
            manifest = extractor.export(self.source, self.output, CONFIG)
        self.assertEqual(source_bytes, self.source.read_bytes())
        self.assertEqual(manifest["schema_version"], 1)
        self.assertEqual(manifest["dataset"], "Stephenson")
        self.assertEqual(manifest["view"], "benchmark_analysis")
        self.assertEqual(manifest["sample_ids"], ["S1", "S2", "S3", "S4"])
        self.assertEqual(manifest["n_cells"], 8)
        self.assertEqual(len(manifest["files"]), 21)
        self.assertEqual(len(list(self.output.iterdir())), 22)
        self.assertEqual(json.loads((self.output / "manifest.json").read_text()), manifest)
        self.assertFalse(manifest["source"]["full_file_checksum_verified"])
        self.assertEqual(manifest["source"]["sidecar_text"], "MD5=unverified-test-evidence\n")
        self.assertEqual(len(manifest["source"]["selected_obs_sha256"]), 64)
        self.assertEqual(manifest["config"]["sha256"], extractor.sha256(CONFIG))
        self.assertEqual(manifest["script"]["sha256"], extractor.sha256(Path(extractor.__file__)))
        for filename, expected_hash in manifest["files"].items():
            path = self.output / filename
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), expected_hash)
            frame = pd.read_feather(path)
            self.assertEqual(frame.columns[0], "Sample")
            self.assertEqual(frame.Sample.tolist(), manifest["sample_ids"])
            if filename != "samples.feather":
                annotation = filename.split("__")[0]
                self.assertEqual(frame.columns.tolist()[1:], manifest["cell_type_columns"][annotation])
        samples = pd.read_feather(self.output / "samples.feather")
        self.assertEqual(samples.columns.tolist(), ["Sample", "Status"])
        self.assertEqual(samples.Status.tolist(), ["Healthy", "Covid", "Healthy", "Covid"])
        for annotation in manifest["annotations"]:
            self.assertEqual(manifest["excluded_annotation_cells"][annotation],
                             0 if annotation == "hitme_layer2" else 2)
            counts = pd.read_feather(self.output / f"{annotation}__n_cells.feather").set_index("Sample")
            self.assertEqual(counts.loc["S1", "A"], 2)
            self.assertEqual(counts.loc["S2", "A"], 1)
            self.assertEqual(counts.loc["S2", "B"], 1)
            self.assertEqual(counts.loc["S1", "B"], 0)
            self.assertEqual(counts.loc["S4", "A"], 2)
            for score in extractor.SCORES:
                means = pd.read_feather(self.output / f"{annotation}__{score}.feather").set_index("Sample")
                self.assertEqual(means.loc["S1", "A"], 2.0 if score == "HeatShock_UCell" else 3.0)
                self.assertEqual(means.loc["S2", "A"], 0.0)  # Actual score zero is not NA.
                self.assertTrue(np.isnan(means.loc["S1", "B"]))  # Absent group.
                self.assertTrue(np.isnan(means.loc["S2", "B"]))  # Present, missing score.
                self.assertTrue(np.isnan(means.loc["S4", "A"]))  # All scores missing.
                if annotation == "hitme_layer2":
                    self.assertEqual(means.loc["S3", "B"], 9.0)
                else:
                    self.assertTrue(means.loc["S3"].isna().all())
                    self.assertEqual(counts.loc["S3"].sum(), 0)
        # A copied run-owned CLI reuses only the checkout identified by --config.
        second = self.directory / "cli-export"
        copied_script = self.directory / "extract_stephenson_ucell.py"
        copied_script.write_bytes(Path(extractor.__file__).read_bytes())
        subprocess.run([sys.executable, str(copied_script),
                        "--h5ad", str(self.source), "--output-dir", str(second),
                        "--config", str(CONFIG)], check=True, capture_output=True, text=True)
        self.assertEqual(json.loads((second / "manifest.json").read_text())["files"], manifest["files"])
        with self.assertRaises(FileExistsError):
            extractor.export(self.source, self.output, CONFIG)
        empty = self.directory / "already-existing-empty"
        empty.mkdir()
        with self.assertRaises(FileExistsError):
            extractor.export(self.source, empty, CONFIG)

    def test_publish_manifest_last(self):
        real_rename = Path.rename
        published = []

        def recording_rename(source, destination):
            self.assertEqual(destination.parent, self.output)
            self.assertFalse((self.output / "manifest.json").exists())
            if source.name == "manifest.json":
                manifest = json.loads(source.read_text())
                self.assertEqual(set(published), set(manifest["files"]))
                for filename, expected_hash in manifest["files"].items():
                    self.assertEqual(extractor.sha256(self.output / filename), expected_hash)
            else:
                self.assertEqual(source.suffix, ".feather")
            result = real_rename(source, destination)
            published.append(source.name)
            return result

        with patch.object(Path, "rename", recording_rename):
            manifest = extractor.export(self.source, self.output, CONFIG)
        self.assertEqual(published[-1], "manifest.json")
        self.assertEqual(len(published), len(manifest["files"]) + 1)
        self.assertEqual(json.loads((self.output / "manifest.json").read_text()), manifest)

    def test_failed_publication_preserves_partial_and_rejects_retry(self):
        real_rename = Path.rename
        for failure in ("second feather", "manifest"):
            with self.subTest(failure=failure):
                output = self.directory / failure.replace(" ", "-")
                published = []

                def failing_rename(source, destination):
                    self.assertFalse((output / "manifest.json").exists())
                    if (failure == "second feather" and len(published) == 1) or (
                            failure == "manifest" and source.name == "manifest.json"):
                        raise OSError("Injected publication failure")
                    result = real_rename(source, destination)
                    published.append(source.name)
                    return result

                with patch.object(Path, "rename", failing_rename):
                    with self.assertRaisesRegex(OSError, "Injected publication failure"):
                        extractor.export(self.source, output, CONFIG)
                self.assertTrue(published)
                self.assertEqual(len(published), 1 if failure == "second feather" else 21)
                self.assertFalse((output / "manifest.json").exists())
                remaining = {path.name: path.read_bytes() for path in output.iterdir()}
                self.assertEqual(set(remaining), set(published))
                self.assertEqual(list(self.directory.glob(f".{output.name}.*")), [])
                with self.assertRaises(FileExistsError):
                    extractor.export(self.source, output, CONFIG)
                self.assertEqual({path.name: path.read_bytes() for path in output.iterdir()}, remaining)

    def test_publish_never_clobbers_existing_paths(self):
        staging = self.directory / "publish-staging"
        staging.mkdir()
        (staging / "samples.feather").write_bytes(b"staged data")
        (staging / "manifest.json").write_text("{}")
        empty = self.directory / "empty"
        empty.mkdir()
        populated = self.directory / "populated"
        populated.mkdir()
        (populated / "evidence").write_bytes(b"keep")
        existing_file = self.directory / "file"
        existing_file.write_bytes(b"keep")
        symlink = self.directory / "symlink"
        symlink.symlink_to(populated, target_is_directory=True)
        dangling = self.directory / "dangling"
        dangling.symlink_to(self.directory / "missing")
        for target in (empty, populated, existing_file, symlink, dangling):
            with self.subTest(target=target.name):
                with self.assertRaises(FileExistsError):
                    extractor._publish(staging, target)
                with self.assertRaises(FileExistsError):
                    extractor.export(self.source, target, CONFIG)
                self.assertEqual((staging / "samples.feather").read_bytes(), b"staged data")
                self.assertEqual((staging / "manifest.json").read_text(), "{}")
        self.assertEqual(list(empty.iterdir()), [])
        self.assertEqual((populated / "evidence").read_bytes(), b"keep")
        self.assertEqual(existing_file.read_bytes(), b"keep")
        self.assertTrue(symlink.is_symlink())
        self.assertTrue(dangling.is_symlink())

    def test_target_created_after_preflight_is_not_overwritten(self):
        real_publish = extractor._publish

        def concurrent_target(staging, output):
            output.mkdir()
            return real_publish(staging, output)

        with patch.object(extractor, "_publish", side_effect=concurrent_target):
            with self.assertRaises(FileExistsError):
                extractor.export(self.source, self.output, CONFIG)
        self.assertEqual(list(self.output.iterdir()), [])
        self.assertEqual(list(self.directory.glob(".export.*")), [])

    def test_invalid_inputs_do_not_publish(self):
        cases = []
        inconsistent = fixture()
        inconsistent.loc["cell_2", "Status"] = "Covid"
        cases.append(("inconsistent labels", inconsistent, "unique nonmissing label"))
        missing_label = fixture()
        missing_label.loc["cell_0", "Status"] = None
        cases.append(("missing label", missing_label, "missing/blank"))
        missing_sample = fixture()
        missing_sample.loc["cell_0", "Sample"] = None
        cases.append(("missing sample", missing_sample, "missing/blank"))
        infinite = fixture()
        infinite.loc["cell_0", extractor.SCORES[0]] = np.inf
        cases.append(("infinity", infinite, "infinite"))
        nonnumeric = fixture()
        nonnumeric[extractor.SCORES[0]] = ["not-a-number"] * len(nonnumeric)
        cases.append(("nonnumeric", nonnumeric, "non-numeric"))
        for column in fixture().columns:
            cases.append((f"missing {column}", fixture().drop(columns=column), "missing required obs"))
        for name, frame, message in cases:
            with self.subTest(name=name):
                write_source(self.source, frame)
                with self.assertRaisesRegex(ValueError, message):
                    extractor.export(self.source, self.output, CONFIG)
                self.assertFalse(self.output.exists())
        wrong_name = self.directory / "not-the-benchmark.h5ad"
        write_source(wrong_name, fixture())
        with self.assertRaisesRegex(ValueError, "basename"):
            extractor.export(wrong_name, self.output, CONFIG)

    def test_concurrent_source_change_cleans_staging(self):
        real_feather = pd.DataFrame.to_feather
        changed = False

        def changing_source(frame, path, *args, **kwargs):
            nonlocal changed
            result = real_feather(frame, path, *args, **kwargs)
            if not changed:
                changed = True
                with h5py.File(self.source, "r+") as handle:
                    handle.attrs["concurrent_change"] = True
            return result

        with patch.object(pd.DataFrame, "to_feather", changing_source):
            with self.assertRaisesRegex(ValueError, "changed during extraction"):
                extractor.export(self.source, self.output, CONFIG)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.directory.glob(".export.*")), [])


if __name__ == "__main__":
    unittest.main()
