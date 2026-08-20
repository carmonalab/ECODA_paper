#!/usr/bin/env python3
"""create_subsets_hpc.py -- create lightweight diagnostic subsets on the HPC.

Runs on the HPC (login node or debug-cpu/shared-cpu node) using the local
BeeGFS scratch files (${HPC_SCRATCH_DIR}/_downloads/<file>) where I/O is GB/s.
For each dataset, it:
  1. Opens the full .h5ad in backed mode (backed='r').
  2. Runs full-file count_sanity_check() (sampling 200k values on fast storage).
  3. Extracts candidate columns, full obs_summary, cells_per_sample_stats,
     paper_table_compare, and confounding_crosstab on full .obs.
  4. Runs sample-first, RAM-bounded subset_by_samples() (~10k cells, 10-20 samples).
  5. Writes <Name>_subset.h5ad and <Name>_meta.json into the output directory.

Usage (repo root on HPC):
  pixi run -e default python notebooks/dataset_onboarding/create_subsets_hpc.py
  pixi run -e default python notebooks/dataset_onboarding/create_subsets_hpc.py --only breast
  pixi run -e default python notebooks/dataset_onboarding/create_subsets_hpc.py --out-dir /path/to/subsets
"""

from __future__ import annotations

import argparse
import functools
import json
import os
import sys
import time
from pathlib import Path

# Unbuffered stdout for SLURM live log monitoring
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
print = functools.partial(print, flush=True)

# Ensure onboarding_utils can be imported
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import onboarding_utils as ou

DATASET_CATALOG = [
    {
        "name": "Alzheimer",
        "key": "alzheimer",
        "file": "SEAAD_Alzheimer.h5ad",
        "expected": {"cells": 1395601, "samples": 83, "ct_types": 18},
        "sample_col_hint": "donor_id",
        "bio_col_hint": "Cognitive status",
        "batch_cols_hint": ["assay", "tissue_type", "PMI"],
        "ct_col_hint": "cell_type",
    },
    {
        "name": "Breast_cancer",
        "key": "breast",
        "file": "BreastCncr_processed.h5ad",
        "expected": {"cells": 714331, "samples": 126, "ct_types": 10},
        "sample_col_hint": "donor_id",
        "bio_col_hint": "disease",
        "batch_cols_hint": ["assay", "sequencing_platform", "sample_preservation_method", "suspension_dissociation_time"],
        "ct_col_hint": "broad_cell_type",
    },
    {
        "name": "Covid19_PBMC",
        "key": "covid19",
        "file": "Covid19_Ren2021.h5ad",
        "expected": {"cells": 993171, "samples": 151, "ct_types": 10},
        "sample_col_hint": "sampleID",
        "bio_col_hint": "CoVID-19 severity",
        "batch_cols_hint": ["Single cell sequencing platform", "City", "datasets", "Sample type"],
        "ct_col_hint": "celltype",
    },
    {
        "name": "Diabetes",
        "key": "diabetes",
        "file": "diabetes.h5ad",
        "expected": {"cells": 264235, "samples": 52, "ct_types": 13},
        "sample_col_hint": "donor_id",
        "bio_col_hint": "disease",
        "batch_cols_hint": ["dataset", "design", "assay"],
        "ct_col_hint": "cell_type",
    },
    {
        "name": "Kidney_KPMP",
        "key": "kidney",
        "file": "Kidney_KPMP.h5ad",
        "expected": {"cells": 104314, "samples": 45, "ct_types": 14},
        "sample_col_hint": "donor_id",
        "bio_col_hint": "condition.l1",
        "batch_cols_hint": ["assay", "tissue_type", "region.l1", "library"],
        "ct_col_hint": "subclass.l1",
    },
    {
        "name": "Lupus_PBMC",
        "key": "lupus",
        "file": "Lupus_Perez2022.h5ad",
        "expected": {"cells": 1263676, "samples": 261, "ct_types": 11},
        "sample_col_hint": "sampleID",
        "bio_col_hint": "Status",
        "batch_cols_hint": ["batch_cov", "Processing_Cohort"],
        "ct_col_hint": "cell_types",
    },
    {
        "name": "Lung",
        "key": "lung",
        "file": "lungatlas.h5ad",
        "expected": {"cells": 941504, "samples": 165, "ct_types": 12},
        "sample_col_hint": "sample",
        "bio_col_hint": "disease",
        "batch_cols_hint": ["dataset", "study", "platform", "assay"],
        "ct_col_hint": "cell_type",
    },
    {
        "name": "Myocardial_infarction",
        "key": "myocardial",
        "file": "Myocardial_Infarc_2.h5ad",
        "expected": {"cells": 132888, "samples": 23, "ct_types": 11},
        "sample_col_hint": "patient",
        "bio_col_hint": "patient_group",
        "batch_cols_hint": ["batch", "sampleType"],
        "ct_col_hint": "cell_type",
    },
    {
        "name": "Parkinson",
        "key": "parkinson",
        "file": "Parkinson.h5ad",
        "expected": {"cells": 2096155, "samples": 97, "ct_types": 11},
        "sample_col_hint": "donor_id",
        "bio_col_hint": "disease",
        "batch_cols_hint": ["Brain_bank", "assay", "tissue_type"],
        "ct_col_hint": "cell_type",
    },
]


def resolve_default_paths():
    hpc_scratch = os.environ.get("HPC_SCRATCH_DIR")
    if hpc_scratch:
        in_dir = Path(hpc_scratch) / "_downloads"
        out_dir = in_dir / "subsets"
    else:
        # Local / NAS fallback
        in_dir = Path(
            "/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/JooM_2025_41097818/output"
        )
        out_dir = HERE.parents[1] / "data" / "new_dataset_checks" / "subsets"
    return in_dir, out_dir


def pick_effective(guess: str | None, candidates: list[str], obs_cols: list[str]) -> str | None:
    if guess is not None and guess in obs_cols:
        return guess
    for c in candidates:
        if c in obs_cols:
            return c
    return None


def process_dataset(ds_info: dict, in_dir: Path, out_dir: Path, cfg: dict, verbose: bool = True, skip_existing: bool = False) -> bool:
    name = ds_info["name"]
    fname = ds_info["file"]
    fpath = in_dir / fname

    out_h5ad = out_dir / f"{name}_subset.h5ad"
    meta_json_path = out_dir / f"{name}_meta.json"
    if skip_existing and out_h5ad.exists() and meta_json_path.exists() and out_h5ad.stat().st_size > 1000 and meta_json_path.stat().st_size > 100:
        print(f"[{name}] Output subset and metadata JSON already exist ({out_h5ad.name}, {meta_json_path.name}), skipping.")
        return True

    if not fpath.exists():
        print(f"[{name}] ERROR: input file not found: {fpath}")
        return False

    print(f"\n==============================================================================")
    print(f"[{name}] Processing {fpath.name} (size: {fpath.stat().st_size / 1e9:.2f} GB)...")
    print(f"==============================================================================")
    t0 = time.perf_counter()

    # 1. Backed read
    adata = sc.read_h5ad(fpath, backed="r")
    t_read = time.perf_counter()
    if verbose:
        print(f"[{name}] Opened backed h5ad in {t_read - t0:.2f}s: shape {adata.shape}")

    # 2. File structure info
    structure_info = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "layers": list(adata.layers.keys()),
        "obsm": list(adata.obsm.keys()),
        "obsp": list(adata.obsp.keys()),
        "has_raw": adata.raw is not None,
    }

    # 3. Full-file count sanity check
    print(f"[{name}] Running count_sanity_check()...")
    cs = ou.count_sanity_check(adata, n_sample_cells=200_000, verbose=verbose)
    t_cs = time.perf_counter()
    if verbose:
        print(f"[{name}] Count sanity verdict: {cs.get('verdict')} (took {t_cs - t_read:.2f}s)")

    # 4. Column detection & metadata
    obs = adata.obs
    cand = ou.candidate_col_detection(obs)
    sample_col = pick_effective(ds_info.get("sample_col_hint"), cand["sample"], obs.columns.tolist())
    bio_col = pick_effective(ds_info.get("bio_col_hint"), cand["label"], obs.columns.tolist())
    ct_col = pick_effective(ds_info.get("ct_col_hint"), cand["cell_type"], obs.columns.tolist())

    batch_cols = [c for c in ds_info.get("batch_cols_hint", []) if c in obs.columns]
    if not batch_cols:
        batch_cols = [c for c in cand["batch"] if c not in (sample_col, bio_col, ct_col)][:4]

    print(f"[{name}] Effective columns: sample={sample_col!r}, bio={bio_col!r}, ct={ct_col!r}, batch={batch_cols!r}")

    cps = ou.cells_per_sample_stats(obs, sample_col) if sample_col else {}
    obs_sum_df = ou.obs_summary(adata)
    paper_cmp_df = ou.paper_table_compare(obs, sample_col, ct_col, ds_info.get("expected", {})) if (sample_col and ct_col) else pd.DataFrame()
    crosstab_df = ou.confounding_crosstab(obs, bio_col, batch_cols) if (bio_col and batch_cols) else pd.DataFrame()

    # 5. Sample-first subsetting
    print(f"[{name}] Running subset_by_samples()...")
    sub, subset_summary = ou.subset_by_samples(
        adata,
        sample_col=sample_col,
        bio_col=bio_col,
        batch_cols=batch_cols,
        ct_col=ct_col,
        config=cfg,
        verbose=verbose,
    )
    t_sub = time.perf_counter()
    if verbose:
        print(f"[{name}] Extracted subset: {sub.n_obs} cells from {subset_summary['n_samples_selected']} samples (took {t_sub - t_cs:.2f}s)")

    # 6. Save subset .h5ad
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = out_dir / f"{name}_subset.h5ad"
    sub.write_h5ad(out_h5ad)
    print(f"[{name}] Wrote subset h5ad ({out_h5ad.stat().st_size / 1e6:.2f} MB): {out_h5ad}")

    # 7. Save metadata JSON
    meta_json_path = out_dir / f"{name}_meta.json"
    meta_data = {
        "dataset_name": name,
        "canonical_file": fname,
        "full_file_structure": structure_info,
        "count_sanity": cs,
        "candidate_columns": cand,
        "effective_columns": {
            "sample_col": sample_col,
            "bio_col": bio_col,
            "ct_col": ct_col,
            "batch_cols": batch_cols,
        },
        "cells_per_sample_stats": {k: v for k, v in cps.items() if k != "cells_per_sample"},
        "cells_per_sample": cps.get("cells_per_sample", {}),
        "obs_summary": obs_sum_df.to_dict(orient="records"),
        "paper_table_compare": paper_cmp_df.to_dict(orient="records") if not paper_cmp_df.empty else [],
        "confounding_crosstab": crosstab_df.to_dict(orient="records") if not crosstab_df.empty else [],
        "subset_summary": subset_summary,
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    with open(meta_json_path, "w") as f:
        json.dump(meta_data, f, indent=2, default=str)
    print(f"[{name}] Wrote metadata summary: {meta_json_path}")

    # Clean up backed file handle
    try:
        if getattr(adata, "isbacked", False) and hasattr(adata, "file"):
            adata.file.close()
    except Exception:
        pass

    t_end = time.perf_counter()
    print(f"[{name}] COMPLETE in {t_end - t0:.2f}s!")
    return True


def main():
    default_in, default_out = resolve_default_paths()

    parser = argparse.ArgumentParser(description="Create diagnostic subsets on HPC scratch")
    parser.add_argument("--in-dir", type=Path, default=default_in, help=f"Input downloads directory (default: {default_in})")
    parser.add_argument("--out-dir", type=Path, default=default_out, help=f"Output subsets directory (default: {default_out})")
    parser.add_argument("--only", type=str, default=None, help="Process only one dataset key (e.g. alzheimer, breast)")
    parser.add_argument("--max-samples", type=int, default=ou.SUBSET_CONFIG["MAX_SAMPLES"])
    parser.add_argument("--max-cells-per-sample", type=int, default=ou.SUBSET_CONFIG["MAX_CELLS_PER_SAMPLE"])
    parser.add_argument("--cells-target", type=int, default=ou.SUBSET_CONFIG["CELLS_TARGET"])
    parser.add_argument("--cells-max", type=int, default=ou.SUBSET_CONFIG["CELLS_MAX"])
    parser.add_argument("--seed", type=int, default=ou.SUBSET_CONFIG["SEED"])
    parser.add_argument("--skip-existing", action="store_true", default=False, help="Skip datasets with existing subset files")

    args = parser.parse_args()

    cfg = {
        "MAX_SAMPLES": args.max_samples,
        "MAX_CELLS_PER_SAMPLE": args.max_cells_per_sample,
        "CELLS_TARGET": args.cells_target,
        "CELLS_MAX": args.cells_max,
        "SEED": args.seed,
    }

    targets = DATASET_CATALOG
    if args.only:
        only_norm = args.only.lower().replace("-", "_").split("_")[0]
        targets = [
            d for d in DATASET_CATALOG
            if d["key"] == only_norm or d["name"].lower().startswith(only_norm) or d["key"].startswith(only_norm)
        ]
        if not targets:
            print(f"ERROR: Unknown dataset key {args.only!r}. Available: {[d['key'] for d in DATASET_CATALOG]}")
            sys.exit(1)

    print(f"Starting subset generation for {len(targets)} datasets...")
    print(f"Input dir:  {args.in_dir}")
    print(f"Output dir: {args.out_dir}")
    print(f"Config:     {cfg}")

    successes = []
    failures = []
    t_start = time.perf_counter()

    for ds_info in targets:
        ok = process_dataset(ds_info, args.in_dir, args.out_dir, cfg, skip_existing=args.skip_existing)
        if ok:
            successes.append(ds_info["name"])
        else:
            failures.append(ds_info["name"])

    t_total = time.perf_counter() - t_start
    print(f"\n==============================================================================")
    print(f"ALL DONE in {t_total:.2f}s! ({len(successes)} succeeded, {len(failures)} failed)")
    if failures:
        print(f"Failed datasets: {failures}")
    print(f"Output directory: {args.out_dir}")
    print(f"==============================================================================")


if __name__ == "__main__":
    main()
