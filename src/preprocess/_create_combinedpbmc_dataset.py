"""
Create the Combined PBMC dataset (Stephenson + GongSharma + Zhu) for batch effect analysis.

Usage:
    python _create_combinedpbmc_dataset.py

Output: data/combined_pbmc_batch_effect_analysis.h5ad

The script loads each source, applies per-source filters (subset_vars from datasets.json),
standardizes gene symbols, picks a random 15-sample subset of GongSharma (seed 123),
and merges on common genes. No cell type columns are included — those come from the
annotation pipeline running on the combined preprocessed .h5ad later.
"""

import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
import rpy2.robjects as ro
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.gene_utils import standardize_gene_symbols


# ---------------------------------------------------------------------------
# R interop (same as preprocess.py — load .rds files)
# ---------------------------------------------------------------------------
ro.r('source("src/utils/load_all_functions.R")')
ro.r('''
convert_rds_to_raw_h5ad <- function(input_path, output_path) {
  seurat <- readRDS(input_path)
  seurat <- create_clean_seuratv5_object(seurat)
  if (!file.exists(output_path)) {
    seurat@assays$RNA@data <- seurat@assays$RNA@counts
    write_h5ad(seurat, output_path)
    seurat@assays$RNA@data <- NULL
  }
}
''')
convert_rds_to_raw_h5ad_r = ro.globalenv["convert_rds_to_raw_h5ad"]


def load_single(input_name, base_path, output_dir, project_root):
    """Load a single input file (converting .rds to cached .h5ad if needed)."""
    input_path = base_path / input_name
    if str(input_name).endswith(".rds"):
        stem = Path(input_name).stem
        raw_h5ad_path = output_dir / f"{stem}_raw.h5ad"
        convert_rds_to_raw_h5ad_r(str(input_path), str(raw_h5ad_path))
        return sc.read_h5ad(raw_h5ad_path)
    elif str(input_name).endswith(".h5ad"):
        return sc.read_h5ad(input_path)
    else:
        raise ValueError(f"Unsupported format: {input_name}")


def load_input(input_names, base_path, output_dir, project_root):
    """Load one or more input files, concatenating if list."""
    if isinstance(input_names, list):
        adatas = [load_single(n, base_path, output_dir, project_root) for n in input_names]
        return sc.concat(adatas, index_unique="_") if len(adatas) > 1 else adatas[0]
    return load_single(input_names, base_path, output_dir, project_root)


def apply_subset_vars(adata, subset_vars):
    """Apply {"col": {"values": [...], "op": "in"/"notin"}} filters (same as preprocess.py)."""
    if not subset_vars:
        return adata
    mask = pd.Series(True, index=adata.obs_names)
    for col, spec in subset_vars.items():
        if col not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {col}")
        col_mask = adata.obs[col].isin(spec["values"])
        mask &= col_mask if spec.get("op", "in") == "in" else ~col_mask
    return adata[mask].copy()


def load_and_prepare_source(ds_name, view_info, base_path, output_dir, project_root):
    """Load a source dataset and apply per-source filters + gene standardization."""

    input_names = view_info.get("input_file_name")
    if not input_names:
        raise ValueError(f"No input_file_name for {ds_name}")

    adata = load_input(input_names, base_path, output_dir, project_root)

    subset_vars = view_info.get("subset_vars", {})
    if subset_vars:
        adata = apply_subset_vars(adata, subset_vars)

    standardize_gene_symbols(adata)
    adata.var_names_make_unique()

    return adata


def main():
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "datasets.json"
    base_path = project_root / "data"
    output_dir = project_root / "data"

    with open(config_path) as f:
        config = json.load(f)

    source_configs = {
        "Stephenson": {},
        "Gongsharma_cmv_young_males": {},
        "Zhu": {},
    }

    for ds_name, cfg in source_configs.items():
        ds_info = config[ds_name]
        view_info = ds_info["views"]["batch_effect_analysis"]
        cfg["view_info"] = view_info
        cfg["sample_col"] = ds_info.get("columns", {}).get("sample", "Sample")
        cfg["label_col"] = ds_info.get("columns", {}).get("label", "Sample")

    # ---- Stephenson ----
    print("Loading Stephenson...")
    adata_s = load_and_prepare_source(
        "Stephenson", source_configs["Stephenson"]["view_info"],
        base_path, output_dir, project_root
    )
    adata_s.obs["batch"] = "Stephenson"
    adata_s.obs["cond"] = adata_s.obs[source_configs["Stephenson"]["label_col"]].values
    sample_col_s = source_configs["Stephenson"]["sample_col"]
    if sample_col_s != "Sample" and sample_col_s in adata_s.obs.columns:
        adata_s.obs["Sample"] = adata_s.obs[sample_col_s].values
    keep_s = ["Sample", "batch", "cond"]
    adata_s.obs = adata_s.obs[[c for c in keep_s if c in adata_s.obs.columns]]
    print(f"  -> {adata_s.n_obs} cells")

    # ---- GongSharma ----
    print("Loading GongSharma...")
    view_g = source_configs["Gongsharma_cmv_young_males"]["view_info"]
    adata_g = load_and_prepare_source(
        "Gongsharma_cmv_young_males", view_g,
        base_path, output_dir, project_root
    )
    # Subset to 15 random samples (seed 123)
    sample_col_g = source_configs["Gongsharma_cmv_young_males"]["sample_col"]
    unique_samples = adata_g.obs[sample_col_g].unique().tolist()
    rng = np.random.default_rng(123)
    chosen = rng.choice(unique_samples, size=min(15, len(unique_samples)), replace=False)
    print(f"  GongSharma: {len(unique_samples)} total samples, picking {len(chosen)}")
    adata_g = adata_g[adata_g.obs[sample_col_g].isin(chosen)].copy()
    adata_g.obs["batch"] = "GongSharma"
    adata_g.obs["cond"] = "Healthy"
    if sample_col_g != "Sample" and sample_col_g in adata_g.obs.columns:
        adata_g.obs["Sample"] = adata_g.obs[sample_col_g].values
    keep_g = ["Sample", "batch", "cond"]
    adata_g.obs = adata_g.obs[[c for c in keep_g if c in adata_g.obs.columns]]
    print(f"  -> {adata_g.n_obs} cells")

    # ---- Zhu ----
    print("Loading Zhu...")
    view_z = source_configs["Zhu"]["view_info"]
    adata_z = load_and_prepare_source(
        "Zhu", view_z,
        base_path, output_dir, project_root
    )
    adata_z.obs["batch"] = "Zhu"
    adata_z.obs["cond"] = "Healthy"
    sample_col_z = source_configs["Zhu"]["sample_col"]
    if sample_col_z != "Sample" and sample_col_z in adata_z.obs.columns:
        adata_z.obs["Sample"] = adata_z.obs[sample_col_z].values
    keep_z = ["Sample", "batch", "cond"]
    adata_z.obs = adata_z.obs[[c for c in keep_z if c in adata_z.obs.columns]]
    print(f"  -> {adata_z.n_obs} cells")

    # ---- Find common genes ----
    gene_sets = [set(adata.var_names) for adata in [adata_s, adata_g, adata_z]]
    common_genes = sorted(set.intersection(*gene_sets))
    print(f"Common genes across all three sources: {len(common_genes)}")

    if len(common_genes) < 5000:
        print("Warning: fewer than 5000 common genes. Using union instead.")
        common_genes = sorted(set.union(*gene_sets))
        print(f"Union gene count: {len(common_genes)}")

    # ---- Subset each to common genes ----
    adata_s = adata_s[:, adata_s.var_names.isin(common_genes)].copy()
    adata_g = adata_g[:, adata_g.var_names.isin(common_genes)].copy()
    adata_z = adata_z[:, adata_z.var_names.isin(common_genes)].copy()

    # ---- Prefix sample IDs to avoid collisions ----
    adata_s.obs["Sample"] = "Stephenson_" + adata_s.obs["Sample"].astype(str)
    adata_g.obs["Sample"] = "GongSharma_" + adata_g.obs["Sample"].astype(str)
    adata_z.obs["Sample"] = "Zhu_" + adata_z.obs["Sample"].astype(str)

    # ---- Concat ----
    print("Concatenating...")
    adata_combined = sc.concat(
        [adata_s, adata_g, adata_z],
        index_unique=None,
        join="outer",
    )
    print(f"Combined shape: {adata_combined.shape}")

    output_path = project_root / "data" / "combined_pbmc_batch_effect_analysis.h5ad"
    adata_combined.write_h5ad(str(output_path))
    print(f"Saved combined PBMC dataset to: {output_path}")


if __name__ == "__main__":
    main()
