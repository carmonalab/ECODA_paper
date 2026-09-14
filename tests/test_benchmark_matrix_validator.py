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
from src.utils.py.gene_utils import standardize_gene_symbols
from src.utils.py.batch_contract import build_batch_contract_identity


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
    class GeneFixture:
        var_names = ["ENSG00000278232", "ENSG00000278232.1", "CRF-R"]

    genes = GeneFixture()
    standardize_gene_symbols(genes)
    assert genes.var_names == ["CRHR1", "CRHR1", "CRHR1"]
    with tempfile.TemporaryDirectory(prefix="ecoda-matrix-validator-") as raw:
        root = Path(raw)
        (root / "embeddings").mkdir()
        selection = root / "selection.tsv"
        selection.write_text("Adams\tbenchmark_analysis\tmrvi\n")
        write_file_sidecar(selection)
        frame = pd.DataFrame(
            {"s1": [1.0, 0.0], "s2": [0.0, 1.0]}, index=["s1", "s2"]
        )
        for n in (1000, 2000, 3000):
            write_feather(root / "embeddings" / f"Adams_hvg{n}_mrvi_dists.feather", frame)
        unrelated_partial = root / "unrelated" / "nested" / "stale.tmp.123"
        unrelated_partial.parent.mkdir(parents=True)
        unrelated_partial.write_text("stale")
        matrix_artifact_validator.validate(
            root, selection, ["mrvi"], batch=False, exact=True
        )
        selected_partial = (
            root / "embeddings" / "Adams_hvg1000_mrvi_dists.feather.tmp.123"
        )
        selected_partial.write_text("stale")
        try:
            matrix_artifact_validator.validate(
                root, selection, ["mrvi"], batch=False, exact=True
            )
        except ValueError as exc:
            assert "partial benchmark artifacts remain" in str(exc)
        else:
            raise AssertionError("selected adjacent partial was accepted")
        finally:
            selected_partial.unlink()
        matrix_artifact_validator.validate(
            root, selection, ["mrvi"], batch=False, exact=True
        )
        pilotgm_selection = root / "pilotgm-selection.tsv"
        pilotgm_selection.write_text("Adams\tbenchmark_analysis\tpilotgm\n")
        write_file_sidecar(pilotgm_selection)
        write_feather(
            root / "embeddings" / "Adams_hvg2000_highres_pilotgm_dists.feather",
            frame,
        )
        matrix_artifact_validator.validate(
            root, pilotgm_selection, ["pilotgm"], batch=False, exact=True
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

        corrected_root = root / "batch_effect" / "corrected_final"
        corrected_root.joinpath("embeddings").mkdir(parents=True)
        corrected_config = root / "corrected-datasets.json"
        corrected_config.write_text(json.dumps({
            "Adams": {
                "use_for_batch_effect": True,
                "columns": {"batch": "batch"},
                "views": {
                    "batch_effect_corrected": {
                        "output_file_name": "source.h5ad"
                    }
                }
            }
        }))
        corrected_identity = build_batch_contract_identity(
            ["batch"],
            sample_column="Sample",
            method_id="PILOT",
            model_id="embedding_consumer_harmony_v1",
        )
        corrected_artifact = (
            corrected_root
            / "embeddings"
            / "Adams_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather"
        )
        write_feather(corrected_artifact, frame)
        corrected_runtime = Path(f"{corrected_artifact}.runtime.json")
        corrected_runtime.write_text(json.dumps({
            "schema_version": 1,
            "artifact_path": str(corrected_artifact),
            "artifact_md5": hashlib.md5(corrected_artifact.read_bytes()).hexdigest(),
            "dataset": "Adams",
            "method": "PILOT_hvg2000",
            "time_secs": 1.0,
            "mem_GB": 1.0,
            "batch_contract": corrected_identity,
        }))
        write_file_sidecar(corrected_runtime)
        corrected_selection = root / "corrected-selection.tsv"
        corrected_selection.write_text(
            "Adams\tbatch_effect_corrected\tbatch_effect_corrected\n"
        )
        write_file_sidecar(corrected_selection)
        matrix_artifact_validator.validate(
            corrected_root,
            corrected_selection,
            ["pilot"],
            batch=True,
            batch_pass="corrected",
            analysis_variant="corrected_final",
            config_path=corrected_config,
        )
        try:
            matrix_artifact_validator.require_nonempty(
                [corrected_artifact],
                "ordinary corrected Feather consumer",
                expected_batch_contract=corrected_identity,
                require_runtime_batch_contract=True,
                require_corrected_summary=True,
            )
        except ValueError as exc:
            assert exc.__cause__ is not None
            assert "missing validation_summary" in str(exc.__cause__), str(exc)
        else:
            raise AssertionError(
                "ordinary corrected Feather consumer accepted a summary-free runtime"
            )

        wrong_scope = root / "batch-wrong-scope.tsv"
        wrong_scope.write_text("Adams\tbatch_effect_uncorrected\twrong_scope\n")
        write_file_sidecar(wrong_scope)
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
        batch_pilotgm = root / "batch-pilotgm.tsv"
        batch_pilotgm.write_text(
            "Adams\tbatch_effect_uncorrected\tbatch_effect_uncorrected\n"
        )
        write_file_sidecar(batch_pilotgm)
        try:
            matrix_artifact_validator.validate(
                root,
                batch_pilotgm,
                ["pilotgm"],
                batch=True,
                batch_pass="uncorrected",
            )
        except ValueError as exc:
            assert "not scheduled" in str(exc)
        else:
            raise AssertionError("batch PILOT-GM-VAE was accepted")

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
