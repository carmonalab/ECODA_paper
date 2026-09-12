#!/usr/bin/env python3
"""Standalone regression for corrected H5AD obs-only source auditing.

The fixture is deliberately self-contained: it creates a tiny H5AD, an
immutable snapshot identity, runtime provenance, and a run-owned report tree
under ``TemporaryDirectory``.  The auditor is invoked with the interpreter
running this test (the Pixi default environment in the supported command).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
AUDITOR = ROOT / "src" / "utils" / "py" / "audit_corrected_h5ad_source.py"
DATASET = "FixtureH5AD"
VIEW = "batch_effect_corrected"
SNAPSHOT_COMMIT = "0123456789abcdef0123456789abcdef01234567"
RUN_ID = "corrected-h5ad-fixture"


def _digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _write_manifest(path: Path, fields: dict[str, str]) -> None:
    path.write_text(
        "".join(f"{key}={value}\n" for key, value in fields.items()),
        encoding="utf-8",
    )


def _make_snapshot(base: Path, config: dict[str, object]) -> dict[str, Path]:
    """Create the canonical snapshot layout required by the auditor."""
    snapshot = base / "snapshot" / SNAPSHOT_COMMIT
    source_root = snapshot / "tree"
    identity_dir = snapshot / "identity"
    runtime_dir = snapshot / "runtime"
    aux_root = snapshot / "aux"
    for directory in (source_root, identity_dir, runtime_dir, aux_root):
        directory.mkdir(parents=True, exist_ok=True)

    config_path = source_root / "datasets.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (source_root / "config_helper.R").write_text(
        "# corrected H5AD metadata fixture helper\n", encoding="utf-8"
    )
    (source_root / "pixi.toml").write_text(
        "[workspace]\nname = 'corrected-h5ad-fixture'\n", encoding="utf-8"
    )
    (source_root / "pixi.lock").write_text(
        "# corrected H5AD metadata fixture lock\n", encoding="utf-8"
    )

    source_archive = identity_dir / "source.tar"
    with tarfile.open(source_archive, mode="w") as archive:
        archive.add(source_root, arcname="tree")

    source_manifest = identity_dir / "source.manifest"
    _write_manifest(
        source_manifest,
        {
            "FORMAT": "1",
            "SOURCE_ROOT": str(source_root),
            "SOURCE_COMMIT": SNAPSHOT_COMMIT,
            "SOURCE_ARCHIVE_PATH": str(source_archive),
            "SOURCE_ARCHIVE_SHA256": _digest(source_archive, "sha256"),
            "CONFIG_HELPER_SHA256": _digest(source_root / "config_helper.R", "sha256"),
            "DATASETS_SHA256": _digest(config_path, "sha256"),
            "PIXI_TOML_SHA256": _digest(source_root / "pixi.toml", "sha256"),
            "PIXI_LOCK_SHA256": _digest(source_root / "pixi.lock", "sha256"),
            "AUX_ROOT": str(aux_root),
            "SCGATE_DB_BRANCH": "fixture",
        },
    )
    (snapshot / "COMPLETE").write_text("COMPLETE\n", encoding="utf-8")

    runtime_image = runtime_dir / "fixture-runtime-image"
    runtime_image.write_bytes(b"corrected H5AD fixture runtime image\n")
    runtime_manifest = runtime_dir / "fixture-runtime.manifest"
    _write_manifest(runtime_manifest, {"FORMAT": "1", "PROFILE": "fixture"})

    run_root = base / "run-root" / RUN_ID
    run_manifests = run_root / "manifests"
    (run_root / "reports").mkdir(parents=True, exist_ok=True)
    run_manifests.mkdir(parents=True, exist_ok=True)
    run_source_manifest = run_manifests / "source.manifest"
    shutil.copyfile(source_manifest, run_source_manifest)
    runtime_identity = run_manifests / "runtime.identity"
    _write_manifest(
        runtime_identity,
        {
            "RUNTIME_IMAGE": str(runtime_image),
            "RUNTIME_MANIFEST": str(runtime_manifest),
            "RUNTIME_IMAGE_SHA256": _digest(runtime_image, "sha256"),
            "RUNTIME_MANIFEST_SHA256": _digest(runtime_manifest, "sha256"),
            "RUNTIME_IMAGE_SIZE": str(runtime_image.stat().st_size),
            "RUNTIME_MANIFEST_SIZE": str(runtime_manifest.stat().st_size),
            "RUNTIME_PROFILE": "fixture",
        },
    )

    # Mark the snapshot payload immutable after every identity-bound file has
    # been written.  The caller restores permissions before TemporaryDirectory
    # cleanup so this remains portable across POSIX filesystems.
    _make_read_only(snapshot)
    return {
        "source_root": source_root,
        "source_manifest": source_manifest,
        "runtime_identity": runtime_identity,
        "runtime_image": runtime_image,
        "runtime_manifest": runtime_manifest,
        "run_root": run_root,
        "config": config_path,
    }


def _make_read_only(root: Path) -> None:
    for path in root.rglob("*"):
        if path.is_dir():
            path.chmod(0o555)
        else:
            path.chmod(0o444)
    root.chmod(0o555)


def _make_writable(root: Path) -> None:
    if not root.exists():
        return
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir():
            path.chmod(0o755)
        else:
            path.chmod(0o644)
    root.chmod(0o755)


def _write_h5ad(path: Path, *, split_sample: bool = False) -> None:
    samples = np.repeat(["S1", "S2", "S3", "S4", "S5"], 2)
    labels = np.repeat(["case", "case", "control", "control", "case"], 2)
    batches = np.repeat(["batch_a", "batch_a", "batch_b", "batch_b", "batch_b"], 2)
    subset = np.array(["keep"] * 8 + ["drop"] * 2, dtype=object)
    if split_sample:
        # One configured sample has one retained and one dropped row.  The
        # shared sample-consistency helper must reject this before batch audit.
        subset[1] = "drop"
    obs = pd.DataFrame(
        {
            "sample_id": samples,
            "label": labels,
            "batch": batches,
            "subset_group": subset,
        },
        index=[f"cell_{index}" for index in range(len(samples))],
    )
    adata = ad.AnnData(
        X=np.zeros((len(obs), 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["gene_0"]),
    )
    adata.write_h5ad(path)


def _invoke(
    fixture: dict[str, Path], input_path: Path, output_path: Path
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(AUDITOR),
        "--config",
        str(fixture["config"]),
        "--input-file",
        str(input_path),
        "--output",
        str(output_path),
        "--dataset",
        DATASET,
        "--view",
        VIEW,
        "--source-root",
        str(fixture["source_root"]),
        "--source-manifest",
        str(fixture["source_manifest"]),
        "--runtime-identity",
        str(fixture["runtime_identity"]),
        "--run-root",
        str(fixture["run_root"]),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        item for item in (str(ROOT), environment.get("PYTHONPATH", "")) if item
    )
    return subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )


def _assert_completed(result: subprocess.CompletedProcess[str], context: str) -> None:
    assert result.returncode == 0, (
        f"{context} failed with status {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )


def _exercise_fixture(
    base: Path, fixture: dict[str, Path], config: dict[str, object]
) -> None:
    input_dir = base / "input"
    negative_dir = base / "negative-input"
    input_dir.mkdir()
    negative_dir.mkdir()
    input_path = input_dir / f"{DATASET}.h5ad"
    negative_input = negative_dir / f"{DATASET}.h5ad"
    _write_h5ad(input_path)
    _write_h5ad(negative_input, split_sample=True)

    report_path = fixture["run_root"] / "reports" / "Fixture_corrected_source.json"
    positive = _invoke(fixture, input_path, report_path)
    _assert_completed(positive, "positive corrected H5AD audit")
    assert report_path.is_file(), "successful audit did not install its report"
    report = json.loads(report_path.read_text(encoding="utf-8"))

    expected_size = input_path.stat().st_size
    expected_md5 = _digest(input_path, "md5")
    expected_sha256 = _digest(input_path, "sha256")
    assert report["status"] == "SOURCE_METADATA_VALIDATED_H5AD"
    assert report["dataset"] == DATASET
    assert report["view"] == VIEW
    assert report["source_type"] == "h5ad"
    assert report["obs_only"] is True
    assert report["read_scope"] == ["obs"]
    assert report["source_path"] == str(input_path)
    assert report["source_identity"] == {
        "path": str(input_path),
        "size": expected_size,
        "md5": expected_md5,
        "sha256": expected_sha256,
    }

    subset_audit = report["subset_audit"]
    assert subset_audit["total_cells"] == 10
    assert subset_audit["retained_cells"] == 8
    assert subset_audit["dropped_cells"] == 2
    assert subset_audit["total_samples"] == 5
    assert subset_audit["retained_samples"] == 4
    assert subset_audit["dropped_samples"] == 1
    assert subset_audit["retained_sample_ids"] == ["S1", "S2", "S3", "S4"]
    assert subset_audit["dropped_sample_ids"] == ["S5"]
    assert subset_audit["split_sample_count"] == 0
    assert subset_audit["split_sample_ids"] == []
    assert subset_audit["configured"] == config[DATASET]["views"][VIEW]["subset_vars"]

    report_config = report["config"]
    assert report_config["sample_column"] == "sample_id"
    assert report_config["label_column"] == "label"
    assert report_config["batch_keys"] == ["batch"]
    assert report_config["subset_vars"] == config[DATASET]["views"][VIEW]["subset_vars"]

    validation = report["validation_summary"]
    assert validation["valid"] is True
    assert validation["estimable"] is True
    assert validation["sample_column"] == "sample_id"
    assert validation["biological_column"] == "label"
    assert validation["batch_keys"] == ["batch"]
    assert validation["n_cells"] == 8
    assert validation["n_samples"] == 4
    assert validation["key_level_counts"] == {"batch": 2}
    assert validation["design_rank"] == 2
    assert validation["design_columns"] == 2
    assert validation["fingerprint"]

    provenance = report["provenance"]
    assert provenance["source_root"] == str(fixture["source_root"])
    assert provenance["source_manifest"]["path"] == str(fixture["source_manifest"])
    assert provenance["source_manifest"]["SOURCE_ROOT"] == str(fixture["source_root"])
    assert provenance["source_manifest_run"]["path"] == str(
        fixture["run_root"] / "manifests" / "source.manifest"
    )
    assert provenance["runtime_identity"]["path"] == str(fixture["runtime_identity"])
    assert provenance["runtime_manifest"]["path"] == str(fixture["runtime_manifest"])
    assert provenance["runtime_image"]["path"] == str(fixture["runtime_image"])

    negative_report = fixture["run_root"] / "reports" / "Fixture_split_sample.json"
    negative = _invoke(fixture, negative_input, negative_report)
    assert negative.returncode != 0, (
        "split-sample H5AD unexpectedly passed audit\n"
        f"stdout:\n{negative.stdout}\nstderr:\n{negative.stderr}"
    )
    assert not negative_report.exists(), "failed audit installed a report"


def main() -> None:
    config: dict[str, object] = {
        DATASET: {
            "display_name": "Corrected H5AD metadata fixture",
            "use_for_batch_effect": True,
            "columns": {
                "sample": "sample_id",
                "label": "label",
                "batch": "batch",
            },
            "views": {
                VIEW: {
                    "input_file_name": f"{DATASET}.h5ad",
                    "output_file_name": f"{DATASET}_{VIEW}_ECODAprocessed.h5ad",
                    "subset_vars": {
                        "subset_group": {"values": ["keep"], "op": "in"}
                    },
                }
            },
        }
    }

    with tempfile.TemporaryDirectory(prefix="corrected-h5ad-source-") as temporary:
        base = Path(temporary).resolve()
        fixture = _make_snapshot(base, config)
        try:
            _exercise_fixture(base, fixture, config)
        finally:
            # Restore the deliberately read-only snapshot before
            # TemporaryDirectory removes it, including when an assertion
            # fails and the fixture exits early.
            _make_writable(base / "snapshot")

    print("corrected H5AD source metadata regression: OK")


if __name__ == "__main__":
    main()
