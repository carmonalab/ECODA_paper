#!/usr/bin/env python3
"""Focused contracts for the two-pass batch-effect registry and workers."""

import importlib.util
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "datasets.json"
PY_WORKER = (
    ROOT
    / "src/5_run_benchmark_methods/run_python_sample_embedding_methods/"
    / "1.1.1_benchmark_methods_py.py"
)


def load_worker():
    spec = importlib.util.spec_from_file_location("batch_effect_worker", PY_WORKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    with DATASETS.open() as handle:
        datasets = json.load(handle)

    expected = {
        "Alzheimer": ("donor_id", "Cognitive status", "Subclass", "Supertype"),
        "Breast_cancer": ("sample_id", "disease", "broad_cell_type", "author_cell_type"),
        "Covid19_PBMC": ("sampleID", "CoVID-19 severity", "majorType", "celltype"),
        "Diabetes": ("donor_id", "disease", "cell_type", "cell_type_reannotatedIntegrated"),
        "Kidney_KPMP": ("specimen", "condition.l1", "subclass.l1", "subclass.l3"),
        "Lung": ("sample", "origin", "ann_coarse", "ann_fine"),
        "Lupus_PBMC": ("sampleID", "Status", "layer1", "layer2"),
        "Myocardial_infarction": ("orig_ident", "patient_group", "cell_type", "cell_subtype"),
        "Parkinson": (
            "donor_id",
            "disease",
            "cell_type",
            "leiden_res_5_batch_effect_uncorrected_hvg2000",
        ),
    }
    for name, roles in expected.items():
        entry = datasets[name]
        assert entry["use_for_benchmark"] is False
        assert entry["use_for_batch_effect"] is True
        assert set(entry["views"]) == {"batch_effect_uncorrected", "batch_effect_corrected"}
        assert entry["columns"]["batch"] is None
        cols = entry["columns"]
        assert tuple(cols[key] for key in ("sample", "label", "cell_type_low_res", "cell_type_high_res")) == roles
        for view_name, view in entry["views"].items():
            assert "benchmark" not in view["output_file_name"]
            assert view["output_file_name"].endswith(
                "_batch_effect_analysis_"
                f"{view_name.removeprefix('batch_effect_')}_ECODAprocessed.h5ad"
            )

    assert datasets["Parkinson"]["views"]["batch_effect_corrected"]["columns"][
        "cell_type_high_res"
    ] == "leiden_res_5_batch_effect_corrected_hvg2000_harmony"
    assert datasets["Joanito"]["columns"]["batch"] == "seqtec"
    assert datasets["Stephenson"]["columns"]["batch"] == "Site"

    worker = load_worker()
    fake = SimpleNamespace(
        obsm={
            "X_pca_batch_effect_uncorrected_hvg2000": object(),
            "X_pca_harmony_batch_effect_corrected_hvg2000": object(),
        }
    )
    assert worker.resolve_pass_embedding_key(fake, "batch_effect_uncorrected", 2000) == (
        "X_pca_batch_effect_uncorrected_hvg2000"
    )
    assert worker.resolve_pass_embedding_key(fake, "batch_effect_corrected", 2000) == (
        "X_pca_harmony_batch_effect_corrected_hvg2000"
    )
    try:
        worker.resolve_pass_embedding_key(fake, "batch_effect_corrected", 1000)
    except KeyError:
        pass
    else:
        raise AssertionError("missing exact corrected embedding did not fail")

    # Existing output skips execution, but the pass-qualified name must still
    # be the name the worker recognizes. This catches accidental fallback to
    # ordinary benchmark filenames without loading a cohort-sized h5ad.
    entry = worker.read_datasets_json(str(DATASETS), view="batch_effect_uncorrected")["Alzheimer"]
    with tempfile.TemporaryDirectory() as tmp:
        input_dir = Path(tmp) / "input"
        output_dir = Path(tmp) / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        (input_dir / entry["output_file"]).touch()
        expected_output = output_dir / (
            "Alzheimer_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
        )
        expected_output.touch()
        args = SimpleNamespace(
            view="batch_effect_uncorrected",
            analysis_pass="uncorrected",
            high_resolution_only=True,
            output_dir=str(output_dir),
            input_dir=str(input_dir),
            method="pilot",
            hvg=[2000],
            force=False,
            device="cpu",
        )
        worker.process_dataset(args, "Alzheimer", entry)

    corrected_args = SimpleNamespace(
        view="batch_effect_corrected",
        analysis_pass="corrected",
        high_resolution_only=True,
        output_dir="/tmp/unused-batch-output",
        input_dir="/tmp/unused-batch-input",
        method="pilot",
        hvg=[2000],
        force=False,
        device="cpu",
    )
    corrected_entry = worker.read_datasets_json(str(DATASETS), view="batch_effect_corrected")["Alzheimer"]
    try:
        worker.process_dataset(corrected_args, "Alzheimer", corrected_entry)
    except ValueError as exc:
        assert str(exc) == "corrected batch-effect view requires a confirmed columns.batch"
    else:
        raise AssertionError("corrected null-batch guard did not fail")

    print("batch-effect registry and Python routing OK")


if __name__ == "__main__":
    main()
