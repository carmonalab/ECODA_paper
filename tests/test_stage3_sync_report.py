#!/usr/bin/env python3
"""Focused contract test for validator-only corrected Stage 3 sync reports."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from anndata import AnnData

from src.utils.py.batch_contract import build_batch_contract_identity
from src.utils.py.benchmark_h5ad_contract import REQUIRED_OBSM

REPORT_PATH = ROOT / "src" / "3_scrnaseq_preprocessing" / "validate_stage3_sync_report.py"
REPORT_SPEC = importlib.util.spec_from_file_location("stage3_sync_report", REPORT_PATH)
assert REPORT_SPEC is not None and REPORT_SPEC.loader is not None
report = importlib.util.module_from_spec(REPORT_SPEC)
REPORT_SPEC.loader.exec_module(report)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sidecar(path: Path) -> None:
    path.with_name(path.name + ".md5").write_text(
        f"MD5={md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n",
        encoding="utf-8",
    )


def snapshot(root: Path, commit: str) -> Path:
    tree = root / "tree"
    identity = root / "identity"
    tree.mkdir(parents=True)
    identity.mkdir()
    for name in ("datasets.json", "config_helper.R", "pixi.toml", "pixi.lock"):
        shutil.copy2(ROOT / name, tree / name)
    archive = identity / "source.tar"
    with tarfile.open(archive, "w") as handle:
        handle.add(tree, arcname=".")
    manifest = identity / "source.manifest"
    manifest.write_text(
        "\n".join(
            [
                "FORMAT=1",
                f"SOURCE_ROOT={tree}",
                f"SOURCE_COMMIT={commit}",
                f"SOURCE_ARCHIVE_PATH={archive}",
                f"SOURCE_ARCHIVE_SHA256={sha256(archive)}",
                f"CONFIG_HELPER_SHA256={sha256(tree / 'config_helper.R')}",
                f"DATASETS_SHA256={sha256(tree / 'datasets.json')}",
                f"PIXI_TOML_SHA256={sha256(tree / 'pixi.toml')}",
                f"PIXI_LOCK_SHA256={sha256(tree / 'pixi.lock')}",
                f"AUX_ROOT={tree / 'aux'}",
                "SCGATE_DB_BRANCH=test",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (root / "COMPLETE").write_text("complete\n", encoding="utf-8")
    return manifest


def make_h5ad(path: Path) -> None:
    n_obs = 4
    n_vars = 3000
    obs = pd.DataFrame(
        {"Sample": ["s1", "s1", "s2", "s2"], "assay": ["A", "A", "B", "B"], "sex": ["F", "F", "M", "M"]},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"hvg_rank": np.arange(n_vars, dtype=float)},
        index=[f"gene{i}" for i in range(n_vars)],
    )
    adata = AnnData(X=np.ones((n_obs, n_vars), dtype=float), obs=obs, var=var)
    adata.layers["counts"] = np.ones((n_obs, n_vars), dtype=float)
    for key in REQUIRED_OBSM["batch_effect_corrected"]:
        adata.obsm[key] = np.ones((n_obs, 2), dtype=float)
    adata.uns["batch_contract"] = build_batch_contract_identity(
        ["assay", "sex"], sample_column="Sample", method_id="preprocess",
        model_id="hvg_composite_v1"
    )
    adata.write_h5ad(path)
    sidecar(path)


def invoke(args: list[str]) -> None:
    original = sys.argv
    sys.argv = ["validate_stage3_sync_report.py", *args]
    try:
        report.main()
    finally:
        sys.argv = original


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-stage3-report-") as temporary:
        root = Path(temporary)
        run_root = root / "run"
        run_manifests = run_root / "manifests"
        run_manifests.mkdir(parents=True)
        prior_terminal = run_root / "status" / "terminal"
        prior_terminal.parent.mkdir(parents=True)
        prior_terminal.write_text(
            "STATE=FAIL\nRUN_ID=fixture\nREASON=prior failure\n", encoding="utf-8"
        )
        run_snapshot = root / "run-snapshot"
        validator_snapshot = root / "validator-snapshot"
        run_source = snapshot(run_snapshot, "run-source")
        validator_source = snapshot(validator_snapshot, "validator-source")
        run_source_copy = run_manifests / "source.manifest"
        shutil.copy2(run_source, run_source_copy)
        runtime_identity = run_manifests / "runtime.identity"
        runtime_identity.write_text("RUNTIME=fixture\n", encoding="utf-8")
        selection = run_manifests / "selection.tsv"
        selection.write_text("Alzheimer\tbatch_effect_corrected\n", encoding="utf-8")
        input_root = root / "input"
        h5ad = input_root / "Alzheimer" / "output" / "SEAAD_Alzheimer_batch_effect_analysis_corrected_ECODAprocessed.h5ad"
        h5ad.parent.mkdir(parents=True)
        make_h5ad(h5ad)
        output = run_manifests / "validated_sync_report.json"
        invoke(
            [
                "--run-root",
                str(run_root),
                "--run-source-manifest",
                str(run_source_copy),
                "--validator-source-manifest",
                str(validator_source),
                "--runtime-identity",
                str(runtime_identity),
                "--selection",
                str(selection),
                "--input-root",
                str(input_root),
                "--config",
                str(run_snapshot / "tree" / "datasets.json"),
                "--output",
                str(output),
            ]
        )
        payload = json.loads(output.read_text(encoding="utf-8"))
        assert payload["stage"] == "stage3"
        assert payload["run_id"] == run_root.name
        assert payload["run_source_manifest"]["source_commit"] == "run-source"
        assert payload["validator_source_manifest"]["source_commit"] == "validator-source"
        assert payload["selection"]["rows"] == 1
        assert payload["rows"][0]["contract"] == "batch_effect_corrected_h5ad_v1"
        assert payload["rows"][0]["md5"] == md5(h5ad)
        assert payload["rows"][0]["size"] == str(h5ad.stat().st_size)
        assert output.with_name(output.name + ".md5").is_file()


    print("stage3 sync report contract: OK")


if __name__ == "__main__":
    main()
