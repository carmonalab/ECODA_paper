"""T3.0 -- validate the onboarding helpers + LISI metrics on the _debug dataset.

Runs the exact notebook pipeline (count sanity, metadata exploration, UMAP
panels, cell-level LISI separation) against the Joanito 5-sample `_debug`
batch-effect view on the NAS:

  /Volumes/Shared/Projects/ECODA_paper/_debug/output/
    JoaI_2022_35773407_debug_5samples_batch_effect_analysis_ECODAprocessed.h5ad

Known ground truth: bio label `sample.origin` (LymphNode/Normal/Tumor), batch
candidates `seqtec` (constant -> single-label) and `Site` (partially
confounded with sample.origin). The preprocessed view carries only namespaced
PCA obsm keys (X_pca_batch_effect_analysis_hvg2000 -- unintegrated; the
X_pca_harmony_* key must never be used), so this exercises the precomputed-X_pca
path used by the metrics step, plus the computed-UMAP path (no X_umap exists).

Usage (from repo root, pixi default env):
  pixi run -e default python notebooks/dataset_onboarding/_debug_validation.py
  or with the onboarding dir as cwd:
  .pixi/envs/default/bin/python _debug_validation.py
"""

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

H5AD = (
    "/Volumes/Shared/Projects/ECODA_paper/_debug/output/"
    "JoaI_2022_35773407_debug_5samples_batch_effect_analysis_ECODAprocessed.h5ad"
)
OUT = ROOT / "data" / "new_dataset_checks" / "_debug"
OUT.mkdir(parents=True, exist_ok=True)

CT_COL = "cell.type"
BIO_COL = "sample.origin"
BATCH_COLS = ["seqtec", "Site"]
SAMPLE_COL = "sample.ID"
EXPECTED = {"cells": 2500, "samples": 5, "ct_types": 10}


def run(cmd: list, **kw):
    print(f"\n$ {' '.join(str(c) for c in cmd)}")
    return subprocess.run(cmd, check=True, **kw)


def main() -> int:
    if not Path(H5AD).exists():
        print(f"ERROR: _debug h5ad not found:\n  {H5AD}\n(NAS not mounted?)")
        return 1

    adata = sc.read_h5ad(H5AD, backed="r")
    print("=== 1. File structure ===")
    print(adata)
    print("obsm:", list(adata.obsm.keys()))
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
    # seqtec is constant -> confounded (single level) or skipped, never scored
    st_col = [c for c in tbl.columns if "seqtec_status" in c][0]
    sc_col = [c for c in tbl.columns if "seqtec_separation" in c][0]
    scored = tbl[st_col] == "ok"
    assert not scored.any(), tbl[[st_col, sc_col]]
    assert tbl.loc[scored, sc_col].isna().all()
    print("seqtec correctly never scored (constant column -> confounded/skipped)")

    adata.file.close()
    print("\nT3.0 validation PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())