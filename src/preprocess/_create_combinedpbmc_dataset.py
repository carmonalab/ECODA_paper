"""
Create the Combined PBMC dataset (Stephenson + GongSharma + Zhu) for batch effect analysis.

Usage:
    python _create_combinedpbmc_dataset.py

Output: data/combined_pbmc_batch_effect_analysis.h5ad

Loads each source sequentially, keeps memory in check by trimming obs columns
early and releasing source objects before concat.
"""

import sys
import gc
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.gene_utils import standardize_gene_symbols
from src.datasets_io import read_datasets_json
from src.preprocess._preprocess_utils import load_input, apply_subset_vars


def load_and_prepare_source(ds_name, entry, base_path, output_dir, view_name=None):
    input_names = entry.get("file_names")
    if not input_names:
        raise ValueError(f"No file_names for {ds_name}")

    adata = load_input(input_names, base_path, output_dir)

    subset_vars = {}
    if view_name:
        subset_vars = entry["views"].get(view_name, {}).get("subset_vars", {})
    if subset_vars:
        adata = apply_subset_vars(adata, subset_vars)

    standardize_gene_symbols(adata)
    adata.var_names_make_unique()

    return adata


def keep_only_cols(adata, cols):
    existing = [c for c in cols if c in adata.obs.columns]
    adata.obs = adata.obs[existing]
    return adata


def main():
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "datasets.json"
    base_path = project_root / "data"
    output_dir = project_root / "data"

    config = read_datasets_json(config_path)

    source_names = ["Stephenson", "Gongsharma_cmv_young_males", "Zhu"]
    missing = [n for n in source_names if n not in config]
    if missing:
        raise KeyError(
            f"Missing sources in datasets.json (required for CombinedPBMC): {missing}"
        )

    keep_base = ["Sample", "batch", "cond"]

    # ---- Stephenson ----
    # Subset to the benchmark_analysis view (Site = Ncl only), consistent
    # with previous code.
    print("Loading Stephenson...")
    entry_s = config["Stephenson"]
    adata_s = load_and_prepare_source(
        "Stephenson", entry_s, base_path, output_dir, view_name="benchmark_analysis"
    )
    adata_s.obs["batch"] = "Stephenson"
    adata_s.obs["cond"] = adata_s.obs[entry_s["label_col"]].values
    sample_col_s = entry_s["sample_col"]
    if sample_col_s != "Sample" and sample_col_s in adata_s.obs.columns:
        adata_s.obs["Sample"] = adata_s.obs[sample_col_s].values
    keep_only_cols(adata_s, keep_base)
    print(f"  -> {adata_s.n_obs} cells")

    # ---- GongSharma ----
    print("Loading GongSharma...")
    entry_g = config["Gongsharma_cmv_young_males"]
    adata_g = load_and_prepare_source("Gongsharma_cmv_young_males", entry_g, base_path, output_dir)
    sample_col_g = entry_g["sample_col"]
    unique_samples = adata_g.obs[sample_col_g].unique().tolist()
    rng = np.random.default_rng(123)
    chosen = rng.choice(unique_samples, size=min(15, len(unique_samples)), replace=False)
    print(f"  GongSharma: {len(unique_samples)} total samples, picking {len(chosen)}")
    adata_g = adata_g[adata_g.obs[sample_col_g].isin(chosen)].copy()
    adata_g.obs["batch"] = "GongSharma"
    adata_g.obs["cond"] = "Healthy"
    if sample_col_g != "Sample" and sample_col_g in adata_g.obs.columns:
        adata_g.obs["Sample"] = adata_g.obs[sample_col_g].values
    keep_only_cols(adata_g, keep_base)
    print(f"  -> {adata_g.n_obs} cells")

    # ---- Zhu ----
    print("Loading Zhu...")
    entry_z = config["Zhu"]
    adata_z = load_and_prepare_source("Zhu", entry_z, base_path, output_dir)
    adata_z.obs["batch"] = "Zhu"
    adata_z.obs["cond"] = "Healthy"
    sample_col_z = entry_z["sample_col"]
    if sample_col_z != "Sample" and sample_col_z in adata_z.obs.columns:
        adata_z.obs["Sample"] = adata_z.obs[sample_col_z].values
    keep_only_cols(adata_z, keep_base)
    print(f"  -> {adata_z.n_obs} cells")

    # ---- Compute common genes from standardized var_names ----
    gene_sets = [set(adata.var_names) for adata in [adata_s, adata_g, adata_z]]
    common_genes = sorted(set.intersection(*gene_sets))
    print(f"Common genes across all three sources: {len(common_genes)}")

    if len(common_genes) < 5000:
        print("Warning: fewer than 5000 common genes. Using union instead.")
        common_genes = sorted(set.union(*gene_sets))
        print(f"Union gene count: {len(common_genes)}")

    # Subset each to common genes in-place (no .copy() to avoid doubling memory)
    adata_s = adata_s[:, adata_s.var_names.isin(common_genes)]
    adata_g = adata_g[:, adata_g.var_names.isin(common_genes)]
    adata_z = adata_z[:, adata_z.var_names.isin(common_genes)]

    # ---- Prefix sample IDs to avoid collisions ----
    adata_s.obs["Sample"] = "Stephenson_" + adata_s.obs["Sample"].astype(str)
    adata_g.obs["Sample"] = "GongSharma_" + adata_g.obs["Sample"].astype(str)
    adata_z.obs["Sample"] = "Zhu_" + adata_z.obs["Sample"].astype(str)

    # ---- Concat and release sources ----
    print("Concatenating...")
    adata_combined = sc.concat(
        [adata_s, adata_g, adata_z],
        index_unique=None,
        join="outer",
    )
    print(f"Combined shape: {adata_combined.shape}")

    del adata_s, adata_g, adata_z
    gc.collect()

    output_path = project_root / "data" / "combined_pbmc_batch_effect_analysis.h5ad"
    adata_combined.write_h5ad(str(output_path))
    print(f"Saved combined PBMC dataset to: {output_path}")


if __name__ == "__main__":
    main()
