#!/usr/bin/env python3
"""Focused checks for exact Stage 5 Feather selection validation."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location(
    "matrix_artifact_validator", ROOT / "src/5_run_benchmark_methods/matrix_artifact_validator.py"
)
assert spec and spec.loader
matrix_artifact_validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matrix_artifact_validator)
from src.utils.py.h5ad_source_identity import (
    build_source_identity,
    load_source_identity,
    verify_source_identity,
)


def write_feather(path: Path, frame: pd.DataFrame) -> None:
    frame.to_feather(path)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    path.with_name(f"{path.name}.md5").write_text(
        f"MD5={digest}\nSIZE={path.stat().st_size}\nPATH={path}\n"
    )


def write_h5ad_obs(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.create_dataset("_index", data=np.asarray(["c1", "c2", "c3"], dtype="S2"))
        obs.create_dataset("Sample", data=np.asarray(["s1", "s1", "s2"], dtype="S2"))


def write_file_sidecar(path: Path) -> None:
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    path.with_name(f"{path.name}.md5").write_text(
        f"MD5={digest}\nSIZE={path.stat().st_size}\nPATH={path}\n"
    )


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-matrix-validator-") as raw:
        root = Path(raw)
        (root / "embeddings").mkdir()
        selection = root / "selection.tsv"
        selection.write_text("Adams\tbenchmark_analysis\tmrvi\n")
        frame = pd.DataFrame(
            {"s1": [1.0, 0.0], "s2": [0.0, 1.0]}, index=["s1", "s2"]
        )
        for n in (1000, 2000, 3000):
            write_feather(root / "embeddings" / f"Adams_hvg{n}_mrvi_dists.feather", frame)
        matrix_artifact_validator.validate(
            root, selection, ["mrvi"], batch=False, exact=True
        )
        input_root = root / "input"
        source_h5ad = input_root / "Adams" / "output" / "source.h5ad"
        write_h5ad_obs(source_h5ad)
        write_file_sidecar(source_h5ad)
        config = root / "datasets.json"
        config.write_text(json.dumps({
            "Adams": {"views": {
                "benchmark_analysis": {"output_file_name": "source.h5ad"}
            }}
        }))
        identity = root / "source_identity.json"
        identity.write_text(json.dumps(
            build_source_identity(selection, input_root, config),
            indent=2,
            sort_keys=True,
        ) + "\n")
        write_file_sidecar(identity)
        assert load_source_identity(identity)[("Adams", "benchmark_analysis")]["sample_ids"] == [
            "s1", "s2"
        ]
        verify_source_identity(identity, selection, input_root, config)
        verify_source_identity(
            identity,
            selection,
            input_root,
            config,
            validated_sidecars=True,
        )
        matrix_artifact_validator.validate(
            root,
            selection,
            ["mrvi"],
            batch=False,
            exact=True,
            input_root=input_root,
            config_path=config,
            source_identity=identity,
            source_identity_verified=True,
        )
        source_h5ad.write_bytes(source_h5ad.read_bytes() + b"changed")
        try:
            verify_source_identity(identity, selection, input_root, config)
        except ValueError:
            pass
        else:
            raise AssertionError("changed source identity was accepted")

        wrong_scope = root / "batch-wrong-scope.tsv"
        wrong_scope.write_text("Adams\tbatch_effect_uncorrected\twrong_scope\n")
        try:
            matrix_artifact_validator.validate(
                root,
                wrong_scope,
                ["mrvi"],
                batch=True,
                batch_pass="uncorrected",
            )
        except ValueError:
            pass
        else:
            raise AssertionError("batch scope mismatch was accepted")

        broken = root / "embeddings" / "Adams_hvg2000_mrvi_dists.feather.md5"
        broken.write_text("MD5=00000000000000000000000000000000\n")
        try:
            matrix_artifact_validator.validate(
                root, selection, ["mrvi"], batch=False, exact=True
            )
        except ValueError:
            pass
        else:
            raise AssertionError("invalid checksum was accepted")

        try:
            matrix_artifact_validator.validate(
                root, selection, ["gloscope"], batch=False, exact=True
            )
        except ValueError:
            pass
        else:
            raise AssertionError("unselected exact row scope was accepted")
    print("benchmark matrix validator: OK")


if __name__ == "__main__":
    main()
