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
from dataset_specs import DATASET_SPECS, get_dataset_spec


DATASET_NAMES = tuple(DATASET_SPECS)


def _resolve_dataset_name(raw_name: str) -> str | None:
    """Resolve only an exact canonical name or its unambiguous slug."""
    normalized = raw_name.lower().replace("-", "_")
    matches = []
    for name, spec in DATASET_SPECS.items():
        aliases = {
            name.lower(),
            name.lower().replace("_", ""),
            str(spec["key"]).lower(),
            str(spec["key"]).lower().replace("_", ""),
        }
        if normalized in aliases or name.lower().startswith(normalized):
            matches.append(name)
    if len(matches) == 1:
        return matches[0]
    return None



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




def _df_records(df: pd.DataFrame) -> list[dict]:
    if df is None or df.empty:
        return []
    return json.loads(df.to_json(orient="records"))
def _apply_subset_vars(adata, subset_vars: dict):
    """Apply exact categorical inclusion/exclusion before registry audits."""
    if not subset_vars:
        return adata
    mask = pd.Series(True, index=adata.obs_names)
    for column, rule in subset_vars.items():
        if column not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {column}")
        values = rule.get("values", [])
        column_mask = adata.obs[column].isin(values)
        mask &= column_mask if rule.get("op", "in") == "in" else ~column_mask
    return adata[mask]



def process_dataset(
    spec: dict,
    in_dir: Path,
    out_dir: Path,
    cfg: dict,
    verbose: bool = True,
    skip_existing: bool = False,
) -> bool:
    name = spec["key"]
    fname = spec["file_name"]
    fpath = in_dir / fname

    out_h5ad = out_dir / f"{name}_subset.h5ad"
    meta_json_path = out_dir / f"{name}_meta.json"
    if (
        skip_existing
        and out_h5ad.exists()
        and meta_json_path.exists()
        and out_h5ad.stat().st_size > 1000
        and meta_json_path.stat().st_size > 100
    ):
        print(
            f"[{name}] Existing subset and registry audit found "
            f"({out_h5ad.name}, {meta_json_path.name}); skipping by request."
        )
        return True

    if not fpath.exists():
        print(f"[{name}] ERROR: input file not found: {fpath}")
        return False

    print("\n==============================================================================")
    print(f"[{name}] Processing {fpath.name} (size: {fpath.stat().st_size / 1e9:.2f} GB)...")
    print("==============================================================================")
    t0 = time.perf_counter()

    adata = sc.read_h5ad(fpath, backed="r")
    t_read = time.perf_counter()
    if verbose:
        print(f"[{name}] Opened backed h5ad in {t_read - t0:.2f}s: shape {adata.shape}")

    subset_vars = dict(spec.get("subset_vars") or {})
    declared_sample_col = (spec.get("registry_roles") or {}).get("sample")
    pre_filter_obs = adata.obs
    pre_filter_sample_col = (
        declared_sample_col
        if declared_sample_col in pre_filter_obs.columns
        else spec["sample_candidates"][0]
    )
    pre_filter_counts = {
        "n_cells": int(adata.n_obs),
        "n_samples": (
            int(pre_filter_obs[pre_filter_sample_col].nunique(dropna=True))
            if pre_filter_sample_col in pre_filter_obs.columns
            else None
        ),
    }
    if subset_vars:
        adata = _apply_subset_vars(adata, subset_vars)
    post_filter_obs = adata.obs
    post_filter_counts = {
        "n_cells": int(adata.n_obs),
        "n_samples": (
            int(post_filter_obs[pre_filter_sample_col].nunique(dropna=True))
            if pre_filter_sample_col in post_filter_obs.columns
            else None
        ),
    }
    mask_comparison = None
    if {"platform", "assay"}.issubset(pre_filter_obs.columns):
        platform_mask = pre_filter_obs["platform"].astype("string").eq("10x")
        assay_mask = pre_filter_obs["assay"].astype("string").str.contains(
            "10x", regex=False, na=False
        )
        differing = platform_mask != assay_mask
        mask_comparison = {
            "platform_expression": 'platform == "10x"',
            "assay_expression": 'assay contains "10x"',
            "platform_cells": int(platform_mask.sum()),
            "assay_cells": int(assay_mask.sum()),
            "differing_rows": int(differing.sum()),
            "first_differing_obs_names": [
                str(value) for value in pre_filter_obs.index[differing][:20]
            ],
            "platform_samples": int(
                pre_filter_obs.loc[platform_mask, pre_filter_sample_col].nunique(dropna=True)
            ) if pre_filter_sample_col in pre_filter_obs.columns else None,
            "assay_samples": int(
                pre_filter_obs.loc[assay_mask, pre_filter_sample_col].nunique(dropna=True)
            ) if pre_filter_sample_col in pre_filter_obs.columns else None,
            "equivalent": bool(not differing.any()),
        }
    filter_summary = {
        "expression": subset_vars,
        "pre_filter": pre_filter_counts,
        "post_filter": post_filter_counts,
        "mask_comparison": mask_comparison,
    }
    if verbose and subset_vars:
        print(
            f"[{name}] Applied registry subset {subset_vars}: "
            f"{pre_filter_counts} -> {post_filter_counts}"
        )

    structure_info = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "layers": list(adata.layers.keys()),
        "obsm": list(adata.obsm.keys()),
        "obsp": list(adata.obsp.keys()),
        "has_raw": adata.raw is not None,
        "filter": filter_summary,
    }

    print(f"[{name}] Running count_sanity_check()...")
    cs = ou.count_sanity_check(adata, n_sample_cells=200_000, verbose=verbose)
    t_cs = time.perf_counter()
    if verbose:
        print(f"[{name}] Count sanity verdict: {cs.get('verdict')} (took {t_cs - t_read:.2f}s)")

    obs = adata.obs
    expected_source = spec["expected_source"]
    sample_audit, sample_selection = ou.audit_sample_candidates(
        obs,
        spec["sample_candidates"],
        spec["sample_stable_cols"],
        expected_source.get("independent_units"),
    )
    ct_sample = (
        declared_sample_col
        if declared_sample_col in obs.columns
        else sample_selection.get("selected") or spec["sample_candidates"][0]
    )
    ct_column_audit, ct_label_audit, ct_selection = ou.audit_cell_type_candidates(
        obs,
        spec["cell_type_candidates"],
        ct_sample,
        spec["bio_col"],
        spec["batch_cols"],
    )
    registry_gate = ou.build_registry_gate(
        spec,
        obs,
        sample_audit,
        sample_selection,
        ct_column_audit,
        ct_label_audit,
        ct_selection,
    )
    sample_col = registry_gate.get("sample")
    bio_col = registry_gate.get("label")
    low_col = registry_gate.get("cell_type_low_res")
    high_col = registry_gate.get("cell_type_high_res")
    batch_cols = list(spec["batch_cols"])
    existing_batch_cols = [column for column in batch_cols if column in obs.columns]
    print(
        f"[{name}] Registry gate: {registry_gate['status']} "
        f"(sample={sample_col!r}, label={bio_col!r}, low={low_col!r}, high={high_col!r})"
    )

    cps = ou.cells_per_sample_stats(obs, sample_col) if sample_col else {}
    obs_sum_df = ou.obs_summary(adata)
    expected_compare = {
        "cells": expected_source.get("cells"),
        "samples": expected_source.get("independent_units"),
    }
    paper_cmp_df = (
        ou.paper_table_compare(obs, sample_col, high_col or low_col, expected_compare)
        if sample_col and (high_col or low_col)
        else pd.DataFrame()
    )
    crosstab_df = (
        ou.confounding_crosstab(obs, bio_col, existing_batch_cols)
        if bio_col in obs.columns and existing_batch_cols
        else pd.DataFrame()
    )

    subset_summary = {
        "status": "SKIPPED",
        "reason": (
            "declared registry sample role did not pass; full-file metadata gate is authoritative"
            if registry_gate.get("status") == "FAIL"
            else "subset not requested"
        ),
    }
    sub = None
    if registry_gate.get("status") in {"PASS", "PASS_PENDING_DERIVED_ANNOTATION"} and sample_col:
        print(f"[{name}] Running diagnostic subset_by_samples() after full-file audits...")
        sub, subset_summary = ou.subset_by_samples(
            adata,
            sample_col=sample_col,
            bio_col=bio_col if bio_col in obs.columns else None,
            batch_cols=existing_batch_cols,
            ct_col=low_col or high_col,
            config=cfg,
            verbose=verbose,
        )
        if verbose:
            print(
                f"[{name}] Extracted subset: {sub.n_obs} cells from "
                f"{subset_summary['n_samples_selected']} samples "
                f"(took {time.perf_counter() - t_cs:.2f}s)"
            )

    out_dir.mkdir(parents=True, exist_ok=True)
    if sub is not None:
        sub.write_h5ad(out_h5ad)
        print(f"[{name}] Wrote subset h5ad ({out_h5ad.stat().st_size / 1e6:.2f} MB): {out_h5ad}")

    meta_data = {
        "dataset_name": name,
        "canonical_file": fname,
        "full_file_structure": structure_info,
        "filter_summary": filter_summary,
        "count_sanity": cs,
        "candidate_columns": ou.candidate_col_detection(obs),
        "effective_columns": {
            "sample_col": sample_col,
            "bio_col": bio_col,
            "ct_col": high_col or low_col,
            "cell_type_low_res": low_col,
            "cell_type_high_res": high_col,
            "batch_cols": existing_batch_cols,
        },
        "sample_audit": {
            "candidates": _df_records(sample_audit),
            "pairwise": sample_selection.get("pairwise", []),
            "selection": sample_selection,
        },
        "cell_type_column_audit": _df_records(ct_column_audit),
        "cell_type_label_audit": _df_records(ct_label_audit),
        "cell_type_selection": ct_selection,
        "registry_gate": registry_gate,
        "cells_per_sample_stats": {k: v for k, v in cps.items() if k != "cells_per_sample"},
        "cells_per_sample": cps.get("cells_per_sample", {}),
        "obs_summary": _df_records(obs_sum_df),
        "paper_table_compare": _df_records(paper_cmp_df),
        "confounding_crosstab": _df_records(crosstab_df),
        "subset_summary": subset_summary,
        "audit_scope": "full_file_obs; subset is diagnostic only",
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(meta_json_path, "w") as handle:
        json.dump(meta_data, handle, indent=2, default=str, allow_nan=False)
    print(f"[{name}] Wrote metadata and registry gate: {meta_json_path}")

    try:
        if getattr(adata, "isbacked", False) and hasattr(adata, "file"):
            adata.file.close()
    except Exception:
        pass

    print(f"[{name}] COMPLETE in {time.perf_counter() - t0:.2f}s (registry={registry_gate['status']})!")
    return True


def main():
    default_in, default_out = resolve_default_paths()

    parser = argparse.ArgumentParser(
        description="Inspect full metadata, write registry gates, and create diagnostic subsets on HPC scratch"
    )
    parser.add_argument(
        "--in-dir",
        type=Path,
        default=default_in,
        help=f"Input downloads directory (default: {default_in})",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=default_out,
        help=f"Output subsets directory (default: {default_out})",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Process one canonical dataset name or unambiguous slug (e.g. alzheimer, breast)",
    )
    parser.add_argument("--max-samples", type=int, default=ou.SUBSET_CONFIG["MAX_SAMPLES"])
    parser.add_argument(
        "--max-cells-per-sample",
        type=int,
        default=ou.SUBSET_CONFIG["MAX_CELLS_PER_SAMPLE"],
    )
    parser.add_argument("--cells-target", type=int, default=ou.SUBSET_CONFIG["CELLS_TARGET"])
    parser.add_argument("--cells-max", type=int, default=ou.SUBSET_CONFIG["CELLS_MAX"])
    parser.add_argument("--seed", type=int, default=ou.SUBSET_CONFIG["SEED"])
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        default=False,
        help="Skip datasets with existing subset and registry metadata",
    )

    args = parser.parse_args()
    cfg = {
        "MAX_SAMPLES": args.max_samples,
        "MAX_CELLS_PER_SAMPLE": args.max_cells_per_sample,
        "CELLS_TARGET": args.cells_target,
        "CELLS_MAX": args.cells_max,
        "SEED": args.seed,
    }

    targets = list(DATASET_SPECS.values())
    if args.only:
        name = _resolve_dataset_name(args.only)
        if name is None:
            print(
                f"ERROR: Unknown or ambiguous dataset key {args.only!r}. "
                f"Available: {list(DATASET_NAMES)}"
            )
            sys.exit(1)
        targets = [get_dataset_spec(name)]

    print(f"Starting full-file registry audit for {len(targets)} datasets...")
    print(f"Input dir:  {args.in_dir}")
    print(f"Output dir: {args.out_dir}")
    print(f"Config:     {cfg}")

    successes = []
    failures = []
    t_start = time.perf_counter()
    for spec in targets:
        ok = process_dataset(
            spec,
            args.in_dir,
            args.out_dir,
            cfg,
            skip_existing=args.skip_existing,
        )
        if ok:
            successes.append(spec["key"])
        else:
            failures.append(spec["key"])

    t_total = time.perf_counter() - t_start
    print("\n==============================================================================")
    print(
        f"ALL DONE in {t_total:.2f}s! "
        f"({len(successes)} metadata audits succeeded, {len(failures)} failed)"
    )
    if failures:
        print(f"Failed datasets: {failures}")
    print(f"Output directory: {args.out_dir}")
    print("Registry gates are authoritative; PASS is required before datasets.json registration.")
    print("==============================================================================")


if __name__ == "__main__":
    main()
