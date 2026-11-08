#!/usr/bin/env python3
"""Whole-sample membership, approved registry selections, and safe Stage 3 reuse."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.cohort_selection import (
    cohort_configuration, cohort_record, input_file_identities,
    resolve_cohort_selection, validate_cohort_h5ad_path,
)
from src.utils.py.datasets_io import read_datasets_json


def expect_error(callback, fragment):
    try:
        callback()
    except (ValueError, KeyError, FileExistsError) as exc:
        assert fragment in str(exc), str(exc)
    else:
        raise AssertionError(f"expected failure containing {fragment!r}")


def configuration(dataset, entry, view, minimum=1):
    return cohort_configuration(
        dataset=dataset, view=view, sample_column=entry["sample_col"],
        label_column=entry["label_col"], batch_keys=entry["batch_col"],
        subset_vars=entry["views"][view]["subset_vars"],
        batch_metadata_policy=entry.get("batch_metadata_policy"),
        evaluation_labels=entry.get("evaluation_labels"), min_cells_per_sample=minimum,
    )


def check_registered_selections():
    datasets = read_datasets_json(ROOT / "datasets.json")
    active = {key for key, entry in datasets.items()
              if entry["use_for_batch_effect"] and not key.startswith("_")}
    assert active == {"Alzheimer", "Breast_cancer", "Covid19_PBMC", "Joanito",
                      "Kidney_KPMP_full", "Lupus_PBMC", "Lung", "Stephenson"}
    assert sum(bool(entry["use_for_benchmark"]) for key, entry in datasets.items()
               if not key.startswith("_")) == 11
    for dataset in active:
        views = datasets[dataset]["views"]
        assert views["batch_effect_uncorrected"]["subset_vars"] == views["batch_effect_corrected"]["subset_vars"]
    assert datasets["Kidney_KPMP_full"]["evaluation_labels"] == {"condition": "condition.l1", "tissue": "tissue"}
    assert datasets["Adams"]["evaluation_labels"] == {"primary": "Disease_Identity"}
    assert datasets["Alzheimer"]["batch_col"] == "assay"
    assert datasets["Kidney_KPMP_full"]["batch_col"] == "suspension_type"

    def select(dataset, rows, expected, view="batch_effect_uncorrected"):
        entry = datasets[dataset]
        obs = pd.DataFrame([{entry["sample_col"]: f"raw-{index}", **row}
                            for index, row in enumerate(rows)],
                           index=[f"cell_{index}" for index in range(len(rows))])
        mask, _, _ = resolve_cohort_selection(obs, configuration(dataset, entry, view))
        assert mask.tolist() == expected, (dataset, mask.tolist())

    select("Alzheimer", [{"assay": "10x 3' v3", "Cognitive status": label}
                          for label in ["Dementia", "No dementia", "Reference"]], [True, True, False])
    select("Joanito", [{"seqtec": "3' seq", "sample.origin": label}
                        for label in ["Normal", "Tumor", "Tumor-2", "LymphNode"]], [True, True, False, False])
    day = "Sampling day (Days after symptom onset)"
    select("Covid19_PBMC", [
        {"datasets": batch, "CoVID-19 severity": "control" if value == "control" else "mild/moderate", day: value}
        for batch, value in [("d02", "control"), ("d09", "30"), ("d13", "31"),
                             ("d07", "control"), ("d13", "unknown"), ("d13", "23")]
    ], [True, True, False, False, False, True])
    select("Kidney_KPMP_full", [
        {"condition.l1": "Ref", "suspension_type": "nucleus", "region.l1": region, "tissue": tissue}
        for region, tissue in [("Cortex", "cortex of kidney"), ("Medulla", "renal medulla"),
                               ("Papilla", "renal papilla"), ("NA", "kidney")]
    ], [True, True, False, True])
    good_lupus = ["dmx_AbFlare-4", "dmx_YS-JY-20_pool3", "dmx_YS-JY-20_pool4",
                  "dmx_YS-JY-21_pool1", "dmx_YS-JY-21_pool2", "dmx_YS-JY-22_pool5", "dmx_YS-JY-22_pool6"]
    select("Lupus_PBMC", [{"Status": "Case", "batch_cov": batch}
                           for batch in good_lupus + ["dmx_AbFlare-3", "dmx_YS-JY-23_pool1"]], [True] * 7 + [False, False])
    lung = ["Goveia_Carmeliet_2020", "Lambrechts_Thienpont_2018_6149v2", "Lambrechts_Thienpont_2018_6653",
            "Leader_Merad_2021_10x_3p_v2_beads_cite", "Leader_Merad_2021_10x_3p_v2_digest-deadcell_cite",
            "Leader_Merad_2021_10x_5p_v1_CD2"]
    rows = [{"disease": "lung adenocarcinoma", "platform": "10x", "tissue": "lung", "dataset": batch}
            for batch in lung + ["Adams_Kaminski_2020"]]
    rows += [{**rows[0], "platform": "Smart-seq"}, {**rows[0], "tissue": "blood"},
             {**rows[0], "disease": "chronic obstructive pulmonary disease"}]
    select("Lung", rows, [True] * 6 + [False] * 4)
    # The quality exclusions hold in every view; only batch views exclude Sanger.
    rows = [{"Sample": sample, "Status": "Healthy", "Site": site}
            for sample, site in [("good_ncl", "Ncl"), ("good_cambridge", "Cambridge"),
                                 ("good_sanger", "Sanger"), ("MH8919226", "Ncl"),
                                 ("MH8919227", "Ncl"), ("BGCV10_CV0198", "Ncl"), ("MH8919230", "Ncl")]]
    select("Stephenson", rows, [True, True, False, False, False, False, False])
    select("Stephenson", rows, [True, False, False, False, False, False, False], "benchmark_analysis")


def check_majority_and_identity():
    entry = read_datasets_json(ROOT / "datasets.json")["Breast_cancer"]
    config = configuration("Breast_cancer", entry, "batch_effect_uncorrected", minimum=3)
    obs = pd.DataFrame({
        "sample_id": ["S-A"] * 4 + ["S-B"] * 4 + ["S-C"] * 4,
        "disease": ["breast cancer"] * 8 + ["normal"] * 4,
        "assay": ["10x 3' v3"] * 12,
        "suspension_dissociation_time": ["12 hour"] * 3 + ["2 hour"] + ["unknown"] * 3 + ["12 hour"] + ["15 minute"] * 3 + ["12 hour"],
    }, index=[f"cell_{index}" for index in range(12)])
    original = obs.copy(deep=True)
    mask, sample_frame, audit = resolve_cohort_selection(obs, config)
    assert mask.tolist() == [True] * 4 + [False] * 4 + [True] * 4
    assert sample_frame["suspension_dissociation_time"].tolist() == ["12 hour", "unknown", "15 minute"]
    assert audit["retained_sample_ids"] == ["S-A", "S-C"]
    assert audit["cohort_membership"]["retained_after_cellcount_ids"] == ["S_A", "S_C"]
    pd.testing.assert_frame_equal(obs, original)
    reporting_obs = obs.copy()
    reporting_obs["anatomy"] = ["A"] * 8 + ["B"] * 4
    reporting_config = {**config, "evaluation_labels": {"condition": "disease", "anatomy": "anatomy"}}
    _, _, reporting_audit = resolve_cohort_selection(reporting_obs, reporting_config)
    assert reporting_audit["cohort_membership_sha256"] == audit["cohort_membership_sha256"]
    tie = obs.copy()
    tie.loc[tie.index[0], "suspension_dissociation_time"] = "2 hour"
    expect_error(lambda: resolve_cohort_selection(tie, config), "unique majority")
    mixed_label = obs.copy()
    mixed_label.loc[mixed_label.index[0], "disease"] = "normal"
    expect_error(lambda: resolve_cohort_selection(mixed_label, config), "disagree within raw sample")
    missing_id = obs.copy()
    missing_id.loc[missing_id.index[0], "sample_id"] = None
    expect_error(lambda: resolve_cohort_selection(missing_id, config), "missing or blank")
    collision = obs.copy()
    collision.loc[collision.sample_id.eq("S-C"), "sample_id"] = "S_A"
    expect_error(lambda: resolve_cohort_selection(collision, config), "collide")

    with tempfile.TemporaryDirectory(prefix="ecoda-cohort-") as temporary:
        root = Path(temporary)
        audit_spec = importlib.util.spec_from_file_location(
            "cohort_input_audit", ROOT / "src/3_scrnaseq_preprocessing/1.0_audit_input_views.py"
        )
        input_audit = importlib.util.module_from_spec(audit_spec)
        audit_spec.loader.exec_module(input_audit)
        input_audit.load_audit_input = lambda *args: ad.AnnData(obs=obs.copy())
        report = input_audit.audit_view(
            dataset="Breast_cancer", entry=entry, view_name="batch_effect_uncorrected",
            view=entry["views"]["batch_effect_uncorrected"], input_root=root, output_root=root,
        )
        assert report["retained_cells"] == 8
        assert report["subset_audit"]["retained_sample_ids"] == ["S-A", "S-C"]
        assert report["subset_audit"]["cohort_membership_sha256"]
        source = root / "raw.h5ad"
        source.write_bytes(b"immutable raw-source fixture")
        sources = input_file_identities(source.name, root)
        output = obs.loc[mask].copy()
        output["Sample"] = output.sample_id.str.replace("-", "_", regex=False)
        data = ad.AnnData(obs=output)
        data.uns["cohort_selection"] = cohort_record(audit, sources, output)
        path = root / "derived.h5ad"
        data.write_h5ad(path)
        assert validate_cohort_h5ad_path(path, config, sources)
        assert validate_cohort_h5ad_path(path, reporting_config, sources)
        before = path.read_bytes()
        changed = {**config, "subset_vars": {"assay": {"op": "in", "values": ["10x 3' v2"]}}}
        expect_error(lambda: validate_cohort_h5ad_path(path, changed, sources), "selection/policy differs")
        source.write_bytes(b"changed raw-source fixture")
        expect_error(lambda: validate_cohort_h5ad_path(path, config, input_file_identities(source.name, root)), "source identity differs")
        assert path.read_bytes() == before
        data.obs.loc[data.obs.index[0], "Sample"] = "S_C"
        data.write_h5ad(root / "tampered.h5ad")
        expect_error(lambda: validate_cohort_h5ad_path(root / "tampered.h5ad", config, sources), "obs disagrees")


def check_worker_reuse():
    spec = importlib.util.spec_from_file_location("cohort_preprocess", ROOT / "src/3_scrnaseq_preprocessing/1.1.1_preprocess.py")
    worker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(worker)
    view = "batch_effect_uncorrected"
    dataset = "Fixture"
    with tempfile.TemporaryDirectory(prefix="ecoda-cohort-worker-") as temporary:
        root = Path(temporary)
        raw_config = {dataset: {
            "use_for_benchmark": False, "use_for_batch_effect": True,
            "columns": {"sample": "raw_sample", "label": "biology", "batch": "batch"},
            "views": {view: {"input_file_name": "raw.h5ad", "output_file_name": "derived.h5ad", "subset_vars": {}}},
        }}
        config_path = root / "datasets.json"
        config_path.write_text(json.dumps(raw_config))
        obs = pd.DataFrame({"raw_sample": np.repeat(["S-A", "S-B"], 500),
                            "biology": np.repeat(["case", "control"], 500), "batch": "B1"},
                           index=[f"cell_{i}" for i in range(1000)])
        ad.AnnData(obs=obs).write_h5ad(root / "raw.h5ad")
        entry = read_datasets_json(config_path)[dataset]
        config = configuration(dataset, entry, view, minimum=500)
        mask, _, audit = resolve_cohort_selection(obs, config)
        output_obs = obs.loc[mask].copy()
        output_obs["Sample"] = output_obs.raw_sample.str.replace("-", "_", regex=False)
        counts = sp.csr_matrix((np.ones(1000, dtype=np.float32), (np.arange(1000), np.arange(1000))), shape=(1000, 2000))
        data = ad.AnnData(X=counts.copy(), obs=output_obs,
                          var=pd.DataFrame({"hvg_rank": np.arange(2000, dtype=float)}, index=[f"gene_{i}" for i in range(2000)]))
        data.layers["counts"] = counts
        data.obsm["X_pca_batch_effect_uncorrected_hvg2000"] = np.zeros((1000, 2), dtype=np.float32)
        data.uns["cohort_selection"] = cohort_record(audit, input_file_identities("raw.h5ad", root), output_obs)
        path = root / "derived.h5ad"
        worker._write_h5ad_atomic(data, path, view)
        completed = subprocess.run([
            sys.executable, str(ROOT / "src/utils/py/cohort_selection.py"),
            "--config", str(config_path), "--dataset", dataset, "--view", view,
            "--input-dir", str(root), "--path", str(path),
        ], check=True, capture_output=True, text=True)
        summary = json.loads(completed.stdout)
        assert summary["status"] == "COHORT_VALIDATED" and summary["n_samples"] == 2
        assert summary["membership_sha256"] == audit["cohort_membership_sha256"]
        before = path.read_bytes()
        sidecar = Path(f"{path}.md5").read_bytes()
        worker.load_input = lambda *args: (_ for _ in ()).throw(AssertionError("validated reuse loaded expression data"))
        worker.main(config_path, root, root, ds_name=dataset, view=view)
        raw_config[dataset]["views"][view]["subset_vars"] = {"biology": {"values": ["case"], "op": "in"}}
        config_path.write_text(json.dumps(raw_config))
        expect_error(lambda: worker.main(config_path, root, root, ds_name=dataset, view=view), "Existing artifact preserved")
        expect_error(lambda: worker.main(config_path, root, root, ds_name=dataset, view=view, force=True), "--force cannot replace")
        assert path.read_bytes() == before
        assert Path(f"{path}.md5").read_bytes() == sidecar

        # A fresh RDS attempt must never reuse/repair a canonical raw cache.
        legacy_cache = root / "raw_raw.h5ad"
        legacy_cache.write_bytes(b"preserve old raw cache")
        (root / "raw.rds").write_bytes(b"run-owned RDS source fixture")
        raw_config[dataset]["views"][view].update(
            input_file_name="raw.rds", output_file_name="fresh.h5ad", subset_vars={}
        )
        config_path.write_text(json.dumps(raw_config))

        def isolated_rds_loader(names, input_dir, output_dir):
            cache = Path(output_dir)
            assert cache.parent == root and cache != root
            assert cache.name.startswith(f".{dataset}_{view}_source_")
            assert not (cache / "raw_raw.h5ad").exists()
            (cache / "raw_raw.h5ad").write_bytes(b"attempt-owned conversion")
            return ad.AnnData(X=sp.csr_matrix(np.ones((len(obs), 1))), obs=obs.copy())

        worker.load_input = isolated_rds_loader
        worker.process_view = lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("deliberate processing stop"))
        expect_error(lambda: worker.main(config_path, root, root, ds_name=dataset, view=view), "deliberate processing stop")
        assert legacy_cache.read_bytes() == b"preserve old raw cache"
        assert not list(root.glob(f".{dataset}_{view}_source_*"))
        assert not (root / "fresh.h5ad").exists()
        assert path.read_bytes() == before
        assert Path(f"{path}.md5").read_bytes() == sidecar


if __name__ == "__main__":
    check_registered_selections()
    check_majority_and_identity()
    check_worker_reuse()
    print("PASS: approved registry, whole-sample winners, label consistency, source/membership identity, safe reuse")
