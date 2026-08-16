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

    print("\n=== 4. UMAP panels (computed path: no X_umap obsm) ===")
    um = ou.embed_and_umap_workflow(
        adata,
        label_cols=[BIO_COL] + BATCH_COLS + [CT_COL],
        out_dir=OUT,
        name="debug",
        sample_col=SAMPLE_COL,
    )
    assert um["computed"] is True, "expected computed-UMAP path on _debug"
    print("precomputed PCA candidates:", um["precomputed_pca_keys"])

    print("\n=== 5. Metrics input + LISI separation (precomputed X_pca path) ===")
    info = ou.write_metrics_input(
        adata,
        out_path=OUT / "metrics_input.feather",
        ct_col=CT_COL,
        bio_col=BIO_COL,
        batch_cols=BATCH_COLS,
        sample_col=SAMPLE_COL,
    )
    print("pca source:", info["pc_source"])
    assert "obsm" in info["pc_source"], "expected the precomputed-X_pca path"

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