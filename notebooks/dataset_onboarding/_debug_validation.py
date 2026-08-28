"""T3.0 -- validate the onboarding helpers + LISI metrics on the _debug dataset.

Runs the exact notebook pipeline (count sanity, metadata exploration, UMAP
panels, cell-level LISI separation) against the Joanito 5-sample `_debug`
batch-effect view on the NAS:

  /Volumes/Shared/Projects/ECODA_paper/_debug/output/
    JoaI_2022_35773407_debug_5samples_batch_effect_uncorrected_ECODAprocessed.h5ad

Known ground truth: bio label `sample.origin` (LymphNode/Normal/Tumor), batch
candidates `seqtec` (constant -> single-label) and `Site` (partially
confounded with sample.origin). The preprocessed view carries only namespaced
PCA obsm keys (`X_pca_batch_effect_uncorrected_hvg2000`; the Harmony key is
only valid for the corrected pass), so this exercises the uncorrected
path used by the metrics step, plus the computed-UMAP path (no X_umap exists).

Usage (from repo root, pixi default env):
  pixi run -e default python notebooks/dataset_onboarding/_debug_validation.py
  or with the onboarding dir as cwd:
  .pixi/envs/default/bin/python _debug_validation.py
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import scanpy as sc  # noqa: E402
import pandas as pd  # noqa: E402

import onboarding_utils as ou  # noqa: E402
from dataset_specs import DATASET_SPECS  # noqa: E402


def _sample_fixture(mixed_condition: bool = False) -> pd.DataFrame:
    rows = []
    for donor_index, donor in enumerate(["D1", "D2"]):
        for sample_index in range(2):
            sample = f"{donor}-S{sample_index + 1}"
            condition = (
                ("case" if sample_index == 0 else "control")
                if mixed_condition and donor == "D1"
                else ("case" if donor == "D1" else "control")
            )
            for _ in range(20):
                rows.append(
                    {
                        "donor": donor,
                        "sample": sample,
                        "condition": condition,
                        "sex": "F" if donor_index == 0 else "M",
                    }
                )
    return pd.DataFrame(rows)


def _hierarchy_fixture(
    two_parents: bool = False,
    missing_high: bool = False,
    n_high: int = 30,
    n_low: int | None = None,
) -> pd.DataFrame:
    rows = []
    for sample_index in range(10):
        for high_index in range(n_high):
            for repeat in range(10):
                high = f"H{high_index}"
                if n_low is None:
                    low = f"L{high_index // 3}"
                else:
                    low = f"L{high_index % n_low}"
                if two_parents and high_index == 0 and repeat < 2:
                    low = "L1"
                if missing_high and sample_index == 0 and high_index < 6:
                    high = None
                rows.append(
                    {
                        "sample": f"S{sample_index}",
                        "status": "case" if sample_index < 5 else "control",
                        "batch": "A" if sample_index < 5 else "B",
                        "low": low,
                        "high": high,
                    }
                )
    return pd.DataFrame(rows)


def _batch_restricted_fixture() -> pd.DataFrame:
    rows = []
    for sample_index in range(10):
        for cell_index in range(300):
            if sample_index < 5 and cell_index < 180:
                label = "restricted"
            else:
                label = f"common{cell_index % 9}"
            rows.append(
                {
                    "sample": f"S{sample_index}",
                    "status": "case" if sample_index < 5 else "control",
                    "batch": "A" if sample_index < 5 else "B",
                    "low": label,
                }
            )
    return pd.DataFrame(rows)


def _role_spec(
    sample: str,
    label: str,
    low: str,
    high: str,
    low_source: str = "author",
    high_source: str = "author",
) -> dict:
    return {
        "key": "fixture",
        "expected_source": {},
        "sample_stable_cols": [label, sample],
        "batch_cols": [],
        "registry_roles": {
            "sample": sample,
            "label": label,
            "cell_type_low_res": low,
            "cell_type_high_res": high,
            "annotation_source": {"low": low_source, "high": high_source},
        },
        "decision_notes": ["fixture decision"],
        "initial_registry_mode": "two_pass_batch_effect",
    }


def run_metadata_fixtures() -> int:
    """Exercise raw, declared-role, and derived-annotation gate contracts."""
    valid_sample = _sample_fixture()
    sample_audit, sample_selection = ou.audit_sample_candidates(
        valid_sample,
        ["donor", "sample"],
        ["condition", "sex", "donor"],
        expected_units=2,
    )
    assert sample_selection["status"] == "PASS", sample_selection
    assert sample_selection["selected"] == "donor", sample_selection
    print("PASS valid donor/sample nesting")

    mixed_sample = _sample_fixture(mixed_condition=True)
    _, mixed_selection = ou.audit_sample_candidates(
        mixed_sample,
        ["donor", "sample"],
        ["condition", "sex", "donor"],
        expected_units=2,
    )
    assert mixed_selection["status"] == "PASS", mixed_selection
    assert mixed_selection["selected"] == "sample", mixed_selection
    print("PASS donor invalidated by mixed condition; specimen/sample selected")

    nested = valid_sample.copy()
    nested["low"] = nested["donor"].map({"D1": "L0", "D2": "L1"})
    nested["high"] = nested["donor"].map({"D1": "H0", "D2": "H1"})
    nested_sample_audit, nested_sample_selection = ou.audit_sample_candidates(
        nested,
        ["donor", "sample"],
        ["condition", "sex"],
        expected_units=4,
    )
    nested_columns, _, nested_ct_selection = ou.audit_cell_type_candidates(
        nested,
        ["low", "high"],
        "sample",
        "condition",
        [],
    )
    nested_gate = ou.build_registry_gate(
        _role_spec("sample", "condition", "low", "high"),
        nested,
        nested_sample_audit,
        nested_sample_selection,
        nested_columns,
        pd.DataFrame(),
        nested_ct_selection,
    )
    assert nested_sample_selection["selected"] == "donor", nested_sample_selection
    assert nested_gate["status"] == "PASS", nested_gate
    assert nested_gate["sample"] == "sample", nested_gate
    assert any(w["role"] == "sample" for w in nested_gate["heuristic_warnings"]), nested_gate
    print("PASS declared nested sample overrides highest-unit heuristic")

    collision = nested.copy()
    collision["sample"] = ["A-1", "A_1"] + ["A-2"] * (len(collision) - 2)
    collision_sample_audit, collision_sample_selection = ou.audit_sample_candidates(
        collision,
        ["sample"],
        ["condition", "sex"],
        expected_units=3,
    )
    collision_gate = ou.build_registry_gate(
        _role_spec("sample", "condition", "low", "high"),
        collision,
        collision_sample_audit,
        collision_sample_selection,
        nested_columns,
        pd.DataFrame(),
        nested_ct_selection,
    )
    assert collision_sample_selection["status"] == "FAIL", collision_sample_selection
    assert collision_gate["status"] == "FAIL", collision_gate
    assert any("collision" in reason for reason in collision_gate["reasons"]), collision_gate
    print("FAIL as expected standardized sample-ID collision survives role override")

    hierarchy = _hierarchy_fixture()
    _, _, hierarchy_selection = ou.audit_cell_type_candidates(
        hierarchy,
        ["low", "high"],
        "sample",
        "status",
        ["batch"],
    )
    assert hierarchy_selection["status"] == "PASS", hierarchy_selection
    assert hierarchy_selection["selected"] == {"low": "low", "high": "high"}, hierarchy_selection
    print("PASS valid 10/30-label hierarchy")

    declared_24 = _hierarchy_fixture(n_high=131, n_low=24)
    tier_sample_audit, tier_sample_selection = ou.audit_sample_candidates(
        declared_24,
        ["sample"],
        ["status", "sample"],
        expected_units=10,
    )
    tier_columns, _, tier_selection = ou.audit_cell_type_candidates(
        declared_24,
        ["low", "high"],
        "sample",
        "status",
        ["batch"],
    )
    assert tier_selection["status"] == "FAIL", tier_selection
    tier_gate = ou.build_registry_gate(
        _role_spec("sample", "status", "low", "high"),
        declared_24,
        tier_sample_audit,
        tier_sample_selection,
        tier_columns,
        pd.DataFrame(),
        tier_selection,
    )
    assert tier_gate["status"] == "PASS", tier_gate
    assert tier_gate["cell_type_low_res"] == "low", tier_gate
    assert tier_gate["cell_type_high_res"] == "high", tier_gate
    print("PASS declared 24-label low tier overrides 3--20 heuristic")

    broken_hierarchy = _hierarchy_fixture(two_parents=True)
    _, _, broken_selection = ou.audit_cell_type_candidates(
        broken_hierarchy,
        ["low", "high"],
        "sample",
        "status",
        ["batch"],
    )
    assert broken_selection["status"] == "FAIL", broken_selection
    assert any(
        row.get("status") == "FAIL"
        for row in broken_selection["hierarchy"]
        if row.get("high") == "high" and row.get("low") == "low"
    ), broken_selection
    print("FAIL as expected high label mapped to two parents")

    restricted = _batch_restricted_fixture()
    restricted_columns, _, restricted_selection = ou.audit_cell_type_candidates(
        restricted,
        ["low"],
        "sample",
        "status",
        ["batch"],
    )
    restricted_row = restricted_columns.iloc[0]
    assert restricted_row["gate_status"] == "FAIL", restricted_row.to_dict()
    assert any("restricted" in reason for reason in restricted_row["gate_reasons"]), restricted_row.to_dict()
    assert restricted_selection["status"] == "FAIL", restricted_selection
    print("FAIL as expected abundant batch-restricted label")

    missing_high = _hierarchy_fixture(missing_high=True)
    missing_columns, _, missing_selection = ou.audit_cell_type_candidates(
        missing_high,
        ["low", "high"],
        "sample",
        "status",
        ["batch"],
    )
    missing_row = missing_columns[missing_columns["candidate"] == "high"].iloc[0]
    assert missing_row["missing_fraction"] > 0.01, missing_row.to_dict()
    assert missing_row["gate_status"] == "FAIL", missing_row.to_dict()
    assert missing_selection["status"] == "FAIL", missing_selection
    print("FAIL as expected high tier with >1% missing values")

    pending_source = _hierarchy_fixture(n_high=12, n_low=4).drop(columns=["low", "high"])
    pending_sample_audit, pending_sample_selection = ou.audit_sample_candidates(
        pending_source,
        ["sample"],
        ["status", "sample"],
        expected_units=10,
    )
    pending_columns, _, pending_selection = ou.audit_cell_type_candidates(
        pending_source,
        ["layer1", "layer2"],
        "sample",
        "status",
        [],
    )
    pending_gate = ou.build_registry_gate(
        _role_spec("sample", "status", "layer1", "layer2", "hitme", "hitme"),
        pending_source,
        pending_sample_audit,
        pending_sample_selection,
        pending_columns,
        pd.DataFrame(),
        pending_selection,
    )
    assert pending_gate["status"] == "PASS_PENDING_DERIVED_ANNOTATION", pending_gate
    processed = pending_source.copy()
    processed["layer1"] = processed["status"].map({"case": "L0", "control": "L1"})
    processed["layer2"] = processed["status"].map({"case": "H0", "control": "H1"})
    finalized = ou.finalize_derived_registry_gate(pending_gate, processed)
    assert finalized["status"] == "PASS", finalized
    assert finalized["derived_pending"] == [], finalized
    print("PASS HiTME/Leiden-style roles require post-processing evidence")
    return 0


def _processed_h5ad_candidates(processed_root: Path, name: str, entry: dict) -> list[Path]:
    candidates = []
    for view in (entry.get("views") or {}).values():
        output_name = view.get("output_file_name")
        if not output_name:
            continue
        candidates.extend(
            [
                processed_root / output_name,
                processed_root / name / output_name,
                processed_root / "output" / output_name,
                processed_root / name / "output" / output_name,
            ]
        )
    return list(dict.fromkeys(candidates))


def _load_or_write_postprocess_gate(
    name: str,
    gate: dict,
    entry: dict,
    processed_root: Path,
) -> dict:
    gate_path = processed_root / f"{name}_postprocess_gate.json"
    if gate_path.exists():
        with gate_path.open() as handle:
            return json.load(handle)
    h5ad_candidates = _processed_h5ad_candidates(processed_root, name, entry)
    h5ad_path = next((path for path in h5ad_candidates if path.exists()), None)
    if h5ad_path is None:
        raise FileNotFoundError(
            f"{name}: no processed h5ad found under {processed_root}; "
            f"looked for {[str(path) for path in h5ad_candidates]}"
        )
    adata = sc.read_h5ad(h5ad_path, backed="r")
    try:
        processed_gate = ou.finalize_derived_registry_gate(gate, adata.obs.copy())
    finally:
        if getattr(adata, "isbacked", False) and hasattr(adata, "file"):
            adata.file.close()
    processed_gate["processed_artifact"] = str(h5ad_path)
    temporary = gate_path.with_suffix(".tmp")
    with temporary.open("w") as handle:
        json.dump(processed_gate, handle, indent=2, allow_nan=False)
    os.replace(temporary, gate_path)
    return processed_gate


def run_registry_audit(
    audit_dir: Path,
    config_path: Path,
    processed_artifact_dir: Path | None = None,
) -> int:
    """Compare registry entries with raw and post-processing gate evidence."""
    with config_path.open() as handle:
        config = json.load(handle)
    errors = []
    for name, spec in DATASET_SPECS.items():
        meta_path = audit_dir / f"{name}_meta.json"
        if not meta_path.exists():
            errors.append(f"{name}: missing audit {meta_path}")
            continue
        with meta_path.open() as handle:
            metadata = json.load(handle)
        gate = metadata.get("registry_gate") or {}
        entry = config.get(name)
        if gate.get("status") == "PASS_PENDING_DERIVED_ANNOTATION":
            if processed_artifact_dir is None:
                print(
                    f"DEFERRED {name}: {gate.get('status')} — "
                    f"{'; '.join(gate.get('pending_reasons', []))}"
                )
                continue
            if entry is None:
                errors.append(f"{name}: pending gate has no datasets.json entry")
                continue
            try:
                gate = _load_or_write_postprocess_gate(
                    name, gate, entry, processed_artifact_dir
                )
            except (FileNotFoundError, OSError, ValueError) as exc:
                errors.append(str(exc))
                continue
        if entry is None:
            if gate.get("status") == "PASS":
                errors.append(f"{name}: PASS gate has no datasets.json entry")
            else:
                print(f"DEFERRED {name}: {gate.get('status')} — {'; '.join(gate.get('reasons', []))}")
            continue
        if gate.get("status") == "FAIL":
            errors.append(f"{name}: registry gate failed: {gate.get('reasons', [])}")
            continue
        if gate.get("status") != "PASS":
            print(f"DEFERRED {name}: {gate.get('status')}")
            continue

        selected = gate.get("selected") or {}
        columns = entry.get("columns") or {}
        expected_columns = {
            "sample": selected.get("sample"),
            "label": selected.get("label"),
            "cell_type_low_res": selected.get("cell_type_low_res"),
            "cell_type_high_res": selected.get("cell_type_high_res"),
        }
        for column, expected in expected_columns.items():
            if columns.get(column) != expected:
                errors.append(f"{name}: columns.{column} != registry_gate.selected.{column}")
        if "batch" not in columns or columns["batch"] is not None:
            errors.append(f"{name}: columns.batch must be explicit null")
        if entry.get("folder_name") != "JooM_2025_41097818":
            errors.append(f"{name}: folder_name is not canonical")
        if entry.get("file_names") != spec["file_name"]:
            errors.append(f"{name}: file_names != canonical input filename")
        if entry.get("use_for_benchmark") is not False:
            errors.append(f"{name}: use_for_benchmark must be false")
        if entry.get("use_for_batch_effect") is not True:
            errors.append(f"{name}: use_for_batch_effect must be true")
        print(f"PASS {name}: registry entry matches {meta_path.name}")

    if errors:
        for error in errors:
            print(f"ERROR {error}")
        return 1
    print("Registry audit comparison PASSED")
    return 0

H5AD = (
    "/Volumes/Shared/Projects/ECODA_paper/_debug/output/"
    "JoaI_2022_35773407_debug_5samples_batch_effect_uncorrected_ECODAprocessed.h5ad"
)
OUT = ROOT / "data" / "new_dataset_checks" / "_debug"
OUT.mkdir(parents=True, exist_ok=True)

CT_COL = "cell.type"
BIO_COL = "sample.origin"
BATCH_COLS = ["seqtec", "Site"]
SAMPLE_COL = "Sample"
EXPECTED = {"cells": 2500, "samples": 5, "ct_types": 10}


def run(cmd: list, **kw):
    print(f"\n$ {' '.join(str(c) for c in cmd)}")
    return subprocess.run(cmd, check=True, **kw)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate onboarding metadata gates and the _debug pipeline")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Run deterministic in-memory sample and annotation gate fixtures",
    )
    parser.add_argument(
        "--registry-audit-dir",
        type=Path,
        default=None,
        help="Compare datasets.json entries with *_meta.json registry gates",
    )
    parser.add_argument(
        "--processed-artifact-dir",
        type=Path,
        default=None,
        help="Processed output root used to finalize pending derived registry gates",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "datasets.json",
        help="Registry JSON used with --registry-audit-dir",
    )
    args = parser.parse_args()
    if args.metadata_only:
        status = run_metadata_fixtures()
        if args.registry_audit_dir is not None:
            status = max(
                status,
                run_registry_audit(
                    args.registry_audit_dir,
                    args.config,
                    args.processed_artifact_dir,
                ),
            )
        return status
    if args.registry_audit_dir is not None:
        return run_registry_audit(
            args.registry_audit_dir,
            args.config,
            args.processed_artifact_dir,
        )

    if not Path(H5AD).exists():
        print(f"ERROR: _debug h5ad not found:\n  {H5AD}\n(NAS not mounted?)")
        return 1

    adata = sc.read_h5ad(H5AD, backed="r")
    print("=== 1. File structure ===")
    print(adata)
    print("obsm:", list(adata.obsm.keys()))
    assert int(adata.n_obs) == EXPECTED["cells"], (
        f"expected {EXPECTED['cells']} cells, got {adata.n_obs}"
    )
    assert SAMPLE_COL in adata.obs.columns, f"missing {SAMPLE_COL} column"
    sample_values = adata.obs[SAMPLE_COL].astype("string")
    assert not sample_values.isna().any(), "Sample contains NA values"
    assert not (sample_values.str.strip() == "").any(), "Sample contains blank values"
    sample_counts = sample_values.value_counts(sort=False)
    assert len(sample_counts) == EXPECTED["samples"], sample_counts.to_dict()
    assert set(sample_counts.astype(int).tolist()) == {500}, sample_counts.to_dict()
    assert "counts" in adata.layers, 'missing layers["counts"]'
    required_obs = {"cell.type_new", "seqtec", "Site"}
    missing_obs = sorted(required_obs - set(adata.obs.columns))
    assert not missing_obs, f"missing required derived metadata: {missing_obs}"
    pca_key = "X_pca_batch_effect_uncorrected_hvg2000"
    assert pca_key in adata.obsm, f"missing {pca_key}"
    print("OK: canonical _debug artifact shape, samples, counts, PCA, and metadata")
    print("layers:", dict(adata.layers))

    print("\n=== 2. Count sanity check ===")
    cs = ou.count_sanity_check(adata)
    assert cs["verdict"] == "PASS", cs
    assert cs["slot"] == "layers/counts", cs
    assert cs["integer_values"] is True, cs
    print("OK: PASS, integer values in layers/counts")

    print("\n=== 3. Metadata exploration ===")
    print(ou.obs_summary(adata).to_string())
    print("\ncandidate columns:")
    print(ou.candidate_col_detection(adata.obs))
    print("\ncells per sample:")
    print(pd.Series(ou.cells_per_sample_stats(adata.obs, SAMPLE_COL)["cells_per_sample"]))
    print("\npaper table compare:")
    print(ou.paper_table_compare(adata.obs, SAMPLE_COL, CT_COL, EXPECTED).to_string())
    print("\nconfounding crosstab:")
    print(ou.confounding_crosstab(adata.obs, BIO_COL, BATCH_COLS).to_string())

    print("\n=== 4. UMAP panels (computed fresh from raw counts) ===")
    um = ou.embed_and_umap_workflow(
        adata,
        label_cols=[BIO_COL] + BATCH_COLS + [CT_COL],
        out_dir=OUT,
        name="debug",
        sample_col=SAMPLE_COL,
    )
    assert um["computed"] is True, "expected fresh computed-UMAP on _debug"
    print("used pca_key:", um["pca_key"])

    print("\n=== 4.5 Sample-first subsetting (T3.1 path validation) ===")
    import resource as _res
    import time as _time

    _t0 = _time.perf_counter()
    sub, ssumm = ou.subset_by_samples(
        adata,
        sample_col=SAMPLE_COL,
        bio_col=BIO_COL,
        batch_cols=BATCH_COLS,
        ct_col=CT_COL,
        config={"MAX_SAMPLES": 4, "MAX_CELLS_PER_SAMPLE": 2000, "CELLS_TARGET": 5000},
    )
    _t1 = _time.perf_counter()
    _rss = _res.getrusage(_res.RUSAGE_SELF).ru_maxrss
    print(f"subset wall time: {_t1 - _t0:.2f}s; peak RSS: {_rss} ({_rss / 1e6:.2f} GB if bytes, {_rss / 1e3:.2f} MB if KB)")
    print("summary:", ssumm)
    # structural checks on the reduced config
    assert ssumm["n_samples_selected"] <= 4, ssumm
    assert ssumm["n_cells_total"] <= 5000, ssumm
    assert max(ssumm["cells_per_sample"].values()) <= 2000, ssumm
    assert ssumm["n_samples_selected"] == len(ssumm["cells_per_sample"]), ssumm
    assert set(ssumm["cells_per_sample"]) == set(ssumm["selected_samples"]), ssumm
    assert sub.obs["Site"].nunique(dropna=True) >= 2, "Site should keep >= 2 levels on the subset"
    assert sub.obs["sample.origin"].nunique(dropna=True) >= 2, "bio groups should all be covered"
    print("OK: subset ran within budget; Site/bio still multi-level")
    print("BATCH-SIGNAL CHECK on the subset: running the LISI metrics on the subset")

    info_sub = ou.write_metrics_input(
        sub,
        out_path=OUT / "metrics_input_subset.feather",
        ct_col=CT_COL,
        bio_col=BIO_COL,
        batch_cols=BATCH_COLS,
        sample_col=SAMPLE_COL,
        max_cells=50_000,
    )
    print("pca source (subset):", info_sub["pc_source"])
    assert "unintegrated" in info_sub["pc_source"] or "fresh" in info_sub["pc_source"]

    run([
        "pixi", "run", "-e", "default", "Rscript", "--vanilla",
        str(HERE / "onboarding_metrics.R"),
        "--input", str(OUT / "metrics_input_subset.feather"),
        "--ct-col", CT_COL,
        "--bio-col", BIO_COL,
        "--batch-cols", ",".join(BATCH_COLS),
        "--out-csv", str(OUT / "metrics_separation_subset.csv"),
        "--out-json", str(OUT / "metrics_separation_subset.json"),
    ])
    tbl_sub = pd.read_csv(OUT / "metrics_separation_subset.csv")
    print(tbl_sub.to_string())
    sc_col = [c for c in tbl_sub.columns if "Site_separation" in c][0]
    st_col = [c for c in tbl_sub.columns if "Site_status" in c][0]
    ok_mask = (tbl_sub[st_col] == "ok") & (tbl_sub["cell_type"] != "<ALL>")
    ok_cts = tbl_sub.loc[ok_mask, :]
    if len(ok_cts) and (ok_cts[sc_col] > 0.5).any():
        print("OK: Site batch separation > 0.5 still detected within CT(s) on the subset")
    else:
        print("WARNING: no Site separation > 0.5 on the reduced 4-sample subset")
        print("         acceptable at this small size; full-size subsets are the notebook default")
    assert len(ok_cts) > 0, "expected at least one CT with an 'ok' Site score on the subset (metrics pipeline ran end-to-end)"

    print("\n=== 5. Metrics input + LISI separation (fresh raw-count PCA) ===")
    info = ou.write_metrics_input(
        adata,
        out_path=OUT / "metrics_input.feather",
        ct_col=CT_COL,
        bio_col=BIO_COL,
        batch_cols=BATCH_COLS,
        sample_col=SAMPLE_COL,
    )
    print("pca source:", info["pc_source"])
    assert "unintegrated" in info["pc_source"] or "fresh" in info["pc_source"]

    rscript = [
        "pixi", "run", "-e", "default", "Rscript", "--vanilla",
        str(HERE / "onboarding_metrics.R"),
        "--input", str(OUT / "metrics_input.feather"),
        "--ct-col", CT_COL,
        "--bio-col", BIO_COL,
        "--batch-cols", ",".join(BATCH_COLS),
        "--out-csv", str(OUT / "metrics_separation.csv"),
        "--out-json", str(OUT / "metrics_separation.json"),
    ]
    run(rscript)

    print("\n=== 6. Separation table ===")
    tbl = pd.read_csv(OUT / "metrics_separation.csv")
    print(tbl.to_string())
    import json

    with open(OUT / "metrics_separation.json") as fh:
        jobj = json.load(fh)
    print("\nverdict:")
    for line in jobj["verdict"]:
        print("  " + line)

    print("\n=== 7. Sanity expectations ===")
    allrow = tbl[tbl["cell_type"] == "<ALL>"].iloc[0]
    assert allrow["bio_status"] == "ok", allrow
    print(f"overall bio separation: {allrow['bio_separation']:.3f}")
    # seqtec is multi-level (both 5' and 3' seq) -> scored
    st_col = [c for c in tbl.columns if "seqtec_status" in c][0]
    sc_col = [c for c in tbl.columns if "seqtec_separation" in c][0]
    scored = tbl[st_col] == "ok"
    assert scored.any(), tbl[[st_col, sc_col]]
    print("seqtec correctly scored across cell types (multi-level batch signal detected)")

    adata.file.close()
    print("\nT3.0 validation PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())