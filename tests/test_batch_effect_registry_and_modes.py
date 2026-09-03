#!/usr/bin/env python3
"""Focused contracts for the two-pass batch-effect registry and workers."""

import anndata as ad
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
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
SUBSET_SCRIPT = ROOT / "notebooks/dataset_onboarding/create_subsets_hpc.py"


def load_worker():
    spec = importlib.util.spec_from_file_location("batch_effect_worker", PY_WORKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_subset_worker():
    spec = importlib.util.spec_from_file_location("onboarding_subset_worker", SUBSET_SCRIPT)
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
        "Lung": ("sample", "disease", "ann_coarse", "ann_fine"),
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

    assert datasets["Joanito"]["columns"]["cell_type_low_res"] == "cell.type"
    assert datasets["Joanito"]["columns"]["cell_type_high_res"] == "cell.type_new"
    combined = datasets["CombinedPBMC"]
    assert combined["file_names"] == "combined_pbmc.h5ad"
    assert combined["columns"]["cell_type_low_res"] == "layer1"
    assert combined["columns"]["cell_type_high_res"] == "layer2"
    assert set(combined["views"]) == {"batch_effect_uncorrected", "batch_effect_corrected"}
    assert combined["views"]["batch_effect_uncorrected"] == {
        "input_file_name": "combined_pbmc.h5ad",
        "output_file_name": "combined_pbmc_batch_effect_uncorrected_ECODAprocessed.h5ad",
        "subset_vars": {},
    }
    assert combined["views"]["batch_effect_corrected"] == {
        "input_file_name": "combined_pbmc.h5ad",
        "output_file_name": "combined_pbmc_batch_effect_corrected_ECODAprocessed.h5ad",
        "subset_vars": {},
    }
    spec = importlib.util.spec_from_file_location(
        "dataset_specs", ROOT / "notebooks/dataset_onboarding/dataset_specs.py"
    )
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert len(module.BATCH_EFFECT_DATASET_ORDER) == 12
    assert module.BATCH_EFFECT_DATASET_ORDER[-3:] == ("Joanito", "Stephenson", "CombinedPBMC")
    assert module.BATCH_EFFECT_SPECS["Joanito"] == ["seqtec", "Site"]
    assert module.BATCH_EFFECT_SPECS["Stephenson"] == ["Site"]
    assert module.BATCH_EFFECT_SPECS["CombinedPBMC"] == ["batch"]
    assert module.DEBUG_SPEC["file_name"].endswith(
        "_debug_5samples_batch_effect_uncorrected_ECODAprocessed.h5ad"
    )
    assert module.DEBUG_SPEC["expected_source"]["cells"] == 2500
    assert module.DEBUG_SPEC["expected_source"]["independent_units"] == 5
    stephenson = datasets["Stephenson"]
    assert set(stephenson["views"]) == {
        "benchmark_analysis",
        "batch_effect_uncorrected",
        "batch_effect_corrected",
    }
    assert stephenson["views"]["benchmark_analysis"] == {
        "input_file_name": "StephensonE_2021_33879890_preprocessed.rds",
        "output_file_name": (
            "StephensonE_2021_33879890_preprocessed_"
            "benchmark_analysis_ECODAprocessed.h5ad"
        ),
        "subset_vars": {
            "Site": {"values": ["Ncl"], "op": "in"},
            "Status": {"values": ["Healthy", "Covid"], "op": "in"},
            "Sample": {
                "values": ["BGCV10_CV0198", "MH8919230"],
                "op": "notin",
            },
        },
    }
    expected_stephenson_batch_subset = {
        "Status": {"values": ["Healthy", "Covid"], "op": "in"},
        "Sample": {
            "values": ["BGCV10_CV0198", "MH8919230"],
            "op": "notin",
        },
    }
    for view_name in ("batch_effect_uncorrected", "batch_effect_corrected"):
        view = stephenson["views"][view_name]
        assert view["input_file_name"] == (
            "StephensonE_2021_33879890_preprocessed.rds"
        )
        assert view["output_file_name"] == (
            "StephensonE_2021_33879890_preprocessed_"
            f"batch_effect_analysis_{view_name.removeprefix('batch_effect_')}_"
            "ECODAprocessed.h5ad"
        )
        assert view["subset_vars"] == expected_stephenson_batch_subset
    subset_worker = load_subset_worker()
    safe = subset_worker._json_safe({"nan": float("nan"), "finite": 2.5})
    assert safe == {"nan": None, "finite": 2.5}

    worker = load_worker()
    try:
        worker.validate_gpu_execution("mrvi", "auto")
    except RuntimeError as exc:
        assert "requires --device cuda" in str(exc)
    else:
        raise AssertionError("GPU-backed method accepted implicit auto device")
    original_cuda_available = worker.torch.cuda.is_available
    worker.torch.cuda.is_available = lambda: False
    try:
        try:
            worker.validate_gpu_execution("scpoli", "cuda")
        except RuntimeError as exc:
            assert "torch.cuda.is_available() is False" in str(exc)
        else:
            raise AssertionError("GPU-backed method accepted unavailable CUDA")
    finally:
        worker.torch.cuda.is_available = original_cuda_available
    worker.validate_gpu_execution("pilotgm", "auto")
    worker.validate_gpu_execution("mrvi", "cpu", "hvg1000")
    try:
        worker.validate_gpu_execution("mrvi", "cpu", "hvg2000")
    except RuntimeError as exc:
        assert "default hvg2000 run is H200-only" in str(exc)
    else:
        raise AssertionError("default MrVI accepted CPU execution")

    with tempfile.TemporaryDirectory() as writer_tmp:
        writer_path = Path(writer_tmp) / "embedding.feather"
        writer_frame = pd.DataFrame({"Dim_1": [0.0]}, index=["s1"])
        worker.atomic_to_feather(writer_frame, writer_path)
        assert worker.recorded_feather_valid(writer_path)
        old_bytes = writer_path.read_bytes()
        try:
            worker.atomic_to_feather(pd.DataFrame({"Dim_1": [1.0]}), writer_path)
        except ValueError:
            pass
        else:
            raise AssertionError("Feather writer accepted missing sample identifiers")
        assert writer_path.read_bytes() == old_bytes
        try:
            worker.atomic_to_feather(
                pd.DataFrame({"Dim_1": [float("nan")]}, index=["s1"]),
                writer_path,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Feather writer accepted nonfinite features")
        assert writer_path.read_bytes() == old_bytes
    ordered_adata = SimpleNamespace(
        obs=pd.DataFrame(
            {"Sample": pd.Categorical(["sample_b", "sample_a", "sample_b"])},
            index=["c1", "c2", "c3"],
        )
    )
    assert worker._ordered_sample_ids(ordered_adata) == ["sample_b", "sample_a"]
    square = pd.DataFrame(
        [[0.0, 2.0], [2.0, 0.0]],
        index=["sample_a", "sample_b"],
        columns=["sample_a", "sample_b"],
    )
    aligned_square = worker._align_square_frame(
        square, ["sample_b", "sample_a"], "square.feather"
    )
    assert list(aligned_square.index) == ["sample_b", "sample_a"]
    assert list(aligned_square.columns) == ["sample_b", "sample_a"]
    assert aligned_square.loc["sample_b", "sample_a"] == 2.0
    covariance_adata = SimpleNamespace(
        uns={
            "GMVAE_Representation": {
                "sample_a": {
                    "means": np.zeros((2, 2)),
                    "weights": np.array([0.5, 0.5]),
                    "covariances": np.array(
                        [
                            [[np.nan, 0.0], [0.0, np.nan]],
                            [[1.0, 0.0], [0.0, 1.0]],
                        ]
                    ),
                }
            }
        }
    )
    worker._stabilize_pilotgm_covariances(covariance_adata)
    repaired = covariance_adata.uns["GMVAE_Representation"]["sample_a"][
        "covariances"
    ]
    assert np.isfinite(repaired).all()
    assert (np.linalg.eigvalsh(repaired) >= -1e-10).all()
    benchmark_entries = worker.read_datasets_json(
        str(DATASETS), view="benchmark_analysis"
    )
    benchmark_entry = benchmark_entries["Stephenson"]
    assert benchmark_entry["view_name"] == "benchmark_analysis"
    assert benchmark_entry["input_file"] == (
        "StephensonE_2021_33879890_preprocessed.rds"
    )
    assert benchmark_entry["output_file"] == (
        "StephensonE_2021_33879890_preprocessed_"
        "benchmark_analysis_ECODAprocessed.h5ad"
    )
    assert benchmark_entry["subset_vars"] == stephenson["views"][
        "benchmark_analysis"
    ]["subset_vars"]

    partial = ad.AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"ct": pd.Categorical(["T", None, "B"])},
            index=["c1", "c2", "c3"],
        ),
    )
    worker.fill_unknown_ct(partial, "ct", "test")
    assert list(partial.obs["ct"].astype(str)) == ["T", "Unknown", "B"]

    complete = ad.AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"ct": ["T", "B"]}, index=["c1", "c2"]),
    )
    before = complete.obs["ct"].copy()
    worker.fill_unknown_ct(complete, "ct", "test")
    pd.testing.assert_series_equal(complete.obs["ct"], before)
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
    try:
        worker.resolve_pass_embedding_key(fake, "batch_effect_analysis", 2000)
    except (KeyError, ValueError):
        pass
    else:
        raise AssertionError("legacy analysis embedding fallback was accepted")

    # Existing output skips execution, but the pass-qualified name must still
    # be the name the worker recognizes. This catches accidental fallback to
    # ordinary benchmark filenames without loading a cohort-sized h5ad.
    entry = worker.read_datasets_json(str(DATASETS), view="batch_effect_uncorrected")["Alzheimer"]
    with tempfile.TemporaryDirectory() as tmp:
        input_dir = Path(tmp) / "input"
        output_dir = Path(tmp) / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        input_h5ad = input_dir / entry["output_file"]
        input_adata = ad.AnnData(
            X=np.ones((1, 2000), dtype=np.float32),
            obs=pd.DataFrame({"Sample": ["s1"]}, index=["c1"]),
            var=pd.DataFrame(
                {"hvg_rank": np.arange(2000, dtype=float)},
                index=[f"g{i}" for i in range(2000)],
            ),
        )
        input_adata.layers["counts"] = np.ones((1, 2000), dtype=np.float32)
        input_adata.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.ones(
            (1, 2), dtype=np.float32
        )
        input_adata.write_h5ad(input_h5ad)
        input_digest = hashlib.md5(input_h5ad.read_bytes()).hexdigest()
        input_h5ad.with_name(f"{input_h5ad.name}.md5").write_text(
            f"MD5={input_digest}\nSIZE={input_h5ad.stat().st_size}\nPATH={input_h5ad}\n"
        )
        expected_output = output_dir / (
            "Alzheimer_batch_effect_uncorrected_hvg2000_highres_pilot_dists.feather"
        )
        pd.DataFrame({"Dim_1": [0.0]}, index=["s1"]).to_feather(expected_output)
        output_digest = hashlib.md5(expected_output.read_bytes()).hexdigest()
        expected_output.with_name(f"{expected_output.name}.md5").write_text(
            f"MD5={output_digest}\nSIZE={expected_output.stat().st_size}\nPATH={expected_output}\n"
        )
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
