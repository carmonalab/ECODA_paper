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
FINAL_BATCH_DATASET_ORDER = (
    "Alzheimer",
    "Breast_cancer",
    "Covid19_PBMC",
    "Kidney_KPMP_full",
    "Diabetes",
    "Lupus_PBMC",
    "Lung",
    "Joanito",
    "Stephenson",
)
FINAL_BATCH_METHOD_KEYS = (
    "ECODA_authors_HR",
    "ECODA_seuratres_2",
    "Pseudobulk_hvg2000",
    "GloScope_hvg2000_pcadims30",
    "MrVI_hvg2000",
    "PILOT_hvg2000",
    "QOT_hvg2000",
    "ECODA_authors_HR_NULL",
)


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
    for dataset_name, entry in datasets.items():
        for view_name, view in entry.get("views", {}).items():
            assert "columns" not in view, (
                f"{dataset_name}.{view_name} declares view-level columns"
            )


    expected_final_roles = {
        "Alzheimer": ("donor_id_assay", "Cognitive status", "Subclass", "Supertype"),
        "Breast_cancer": ("sample_id", "disease", "broad_cell_type", "author_cell_type"),
        "Covid19_PBMC": ("sampleID", "CoVID-19 severity", "majorType", "celltype"),
        "Kidney_KPMP_full": ("specimen", "condition.l1", "subclass.l1", "subclass.l3"),
        "Diabetes": ("donor_id", "disease", "cell_type", "cell_type_reannotatedIntegrated"),
        "Lupus_PBMC": ("sampleID", "Status", "layer1", "louvain"),
        "Lung": ("sample", "disease", "ann_coarse", "ann_fine"),
        "Joanito": ("sample.ID", "sample.origin", "cell.type", "cell.type_new"),
        "Stephenson": ("Sample", "Status", "initial_clustering", "full_clustering"),
    }
    assert tuple(expected_final_roles) == FINAL_BATCH_DATASET_ORDER
    assert len(FINAL_BATCH_DATASET_ORDER) == 9
    assert len(set(FINAL_BATCH_DATASET_ORDER)) == 9
    assert all(name in datasets for name in FINAL_BATCH_DATASET_ORDER)
    assert tuple(
        name for name in FINAL_BATCH_DATASET_ORDER if datasets[name]["use_for_batch_effect"]
    ) == FINAL_BATCH_DATASET_ORDER
    assert FINAL_BATCH_METHOD_KEYS == (
        "ECODA_authors_HR",
        "ECODA_seuratres_2",
        "Pseudobulk_hvg2000",
        "GloScope_hvg2000_pcadims30",
        "MrVI_hvg2000",
        "PILOT_hvg2000",
        "QOT_hvg2000",
        "ECODA_authors_HR_NULL",
    )
    final_target_outputs = {
        "Covid19_PBMC": (
            "Covid19_Ren2021_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad",
            "Covid19_Ren2021_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad",
        ),
        "Diabetes": (
            "diabetes_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad",
            "diabetes_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad",
        ),
        "Joanito": (
            "JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad",
            "JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad",
        ),
        "Lung": (
            "lungatlas_batch_effect_analysis_uncorrected_final_ECODAprocessed.h5ad",
            "lungatlas_batch_effect_analysis_corrected_final_ECODAprocessed.h5ad",
        ),
    }
    covid_subset = {
        "Sampling day (Days after symptom onset)": {
            "values": 30,
            "op": "<=",
            "include_values": ["control"],
        }
    }
    assert datasets["Covid19_PBMC"]["views"]["batch_effect_uncorrected"]["subset_vars"] == covid_subset
    assert datasets["Covid19_PBMC"]["views"]["batch_effect_corrected"]["subset_vars"] == covid_subset
    assert datasets["Breast_cancer"]["columns"]["batch"] == [
        "assay",
        "suspension_dissociation_time",
    ]


    legacy_batch_output_names = {
        "Alzheimer": (
            "SEAAD_Alzheimer_donor_assay_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "SEAAD_Alzheimer_donor_assay_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Breast_cancer": (
            "BreastCncr_processed_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "BreastCncr_processed_batch_effect_analysis_corrected_assay_dissociation_ECODAprocessed.h5ad",
        ),
        "Kidney_KPMP": (
            "Kidney_KPMP_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Kidney_KPMP_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Kidney_KPMP_full": (
            "Kidney_KPMP_full_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Kidney_KPMP_full_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Lupus_PBMC": (
            "Lupus_Perez2022_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Lupus_Perez2022_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Myocardial_infarction": (
            "Myocardial_Infarc_2_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Myocardial_Infarc_2_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Parkinson": (
            "Parkinson_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Parkinson_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Stephenson": (
            "StephensonE_2021_33879890_preprocessed_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "StephensonE_2021_33879890_preprocessed_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "CombinedPBMC": (
            "combined_pbmc_batch_effect_uncorrected_ECODAprocessed.h5ad",
            "combined_pbmc_batch_effect_corrected_ECODAprocessed.h5ad",
        ),
    }
    for name, roles in expected_final_roles.items():
        entry = datasets[name]
        assert entry["use_for_batch_effect"] is True
        assert {"batch_effect_uncorrected", "batch_effect_corrected"} <= set(entry["views"])
        cols = entry["columns"]
        assert tuple(cols[key] for key in ("sample", "label", "cell_type_low_res", "cell_type_high_res")) == roles
        output_names = tuple(
            entry["views"][view_name]["output_file_name"]
            for view_name in ("batch_effect_uncorrected", "batch_effect_corrected")
        )
        if name in final_target_outputs:
            assert output_names == final_target_outputs[name]
        else:
            assert output_names == legacy_batch_output_names[name]

    for name in ("CombinedPBMC", "Kidney_KPMP", "Myocardial_infarction", "Parkinson"):
        assert datasets[name]["use_for_benchmark"] is False
        assert datasets[name]["use_for_batch_effect"] is False
    legacy_target_output_names = {
        "Covid19_PBMC": (
            "Covid19_Ren2021_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "Covid19_Ren2021_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Diabetes": (
            "diabetes_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "diabetes_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Joanito": (
            "JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "JoaI_2022_35773407_Nofilt_whole_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
        "Lung": (
            "lungatlas_batch_effect_analysis_uncorrected_ECODAprocessed.h5ad",
            "lungatlas_batch_effect_analysis_corrected_ECODAprocessed.h5ad",
        ),
    }
    final_suffix = "_final_ECODAprocessed.h5ad"
    for name, entry in datasets.items():
        if not {"batch_effect_uncorrected", "batch_effect_corrected"} <= set(entry.get("views", {})):
            continue
        output_names = tuple(
            entry["views"][view_name]["output_file_name"]
            for view_name in ("batch_effect_uncorrected", "batch_effect_corrected")
        )
        assert all((final_suffix in output) == (name in final_target_outputs) for output in output_names)
        if name in legacy_target_output_names:
            assert output_names != legacy_target_output_names[name]
    legacy = datasets["Kidney_KPMP"]
    assert legacy["use_for_benchmark"] is False
    assert legacy["use_for_batch_effect"] is False
    full = datasets["Kidney_KPMP_full"]
    assert full["use_for_benchmark"] is False
    assert full["use_for_batch_effect"] is True
    assert full["file_names"] == "Kidney_KPMP_full.h5ad"
    assert tuple(
        full["columns"][key]
        for key in ("sample", "label", "cell_type_low_res", "cell_type_high_res")
    ) == expected_final_roles["Kidney_KPMP_full"]
    for view in full["views"].values():
        assert view["input_file_name"] == "Kidney_KPMP_full.h5ad"


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
    assert module.BATCH_EFFECT_DATASET_ORDER == (
        "Alzheimer",
        "Breast_cancer",
        "Covid19_PBMC",
        "Kidney_KPMP_full",
        "Myocardial_infarction",
        "Diabetes",
        "Lupus_PBMC",
        "Lung",
        "Parkinson",
        "Joanito",
        "Stephenson",
        "CombinedPBMC",
    )
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
        worker.publish_runtime_metadata(
            expected_output,
            "Alzheimer",
            "PILOT_hvg2000_highres",
            0.0,
            None,
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
    corrected_entry["batch_col"] = None
    try:
        worker.process_dataset(corrected_args, "Alzheimer", corrected_entry)
    except ValueError as exc:
        assert str(exc) == "corrected batch-effect view requires a confirmed columns.batch"
    else:
        raise AssertionError("corrected null-batch guard did not fail")

    print("batch-effect registry and Python routing OK")


if __name__ == "__main__":
    main()
