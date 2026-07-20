import os
import json
import scanpy as sc
import rpy2.robjects as ro
import scipy.sparse as sp
import bionty as bt
from pathlib import Path
import pandas as pd
import re


# ---------------------------------------------------------------------------
# R interop
# ---------------------------------------------------------------------------
ro.r('source("src/R/load_all_functions.R")')
convert_rds_to_raw_h5ad_r = ro.globalenv["convert_rds_to_raw_h5ad"]
 
 
# ---------------------------------------------------------------------------
# Cell (row) subsetting for views
# ---------------------------------------------------------------------------
def apply_subset_vars(adata, subset_vars):
    """Apply {"col": {"values": [...], "op": "in"/"notin"}} filters."""
    if not subset_vars:
        return adata
    mask = pd.Series(True, index=adata.obs_names)
    for col, spec in subset_vars.items():
        if col not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {col}")
        col_mask = adata.obs[col].isin(spec["values"])
        mask &= col_mask if spec.get("op", "in") == "in" else ~col_mask
    return adata[mask].copy()
 
 
# ---------------------------------------------------------------------------
# HVG selection
# ---------------------------------------------------------------------------
def compute_hvgs(adata, n_top_genes, batch_key=None, flavor="seurat_v3"):
    """
    One-shot HVG selection: returns the gene names directly.
    Use this when only a single, fixed HVG size is ever needed (e.g. the
    batch view), since there's no need to keep ranks around for later slicing.
    """
    hvg_df = sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key,
        inplace=False,
        check_values=True,
    )
    return adata.var_names[hvg_df["highly_variable"].values]
 
 
def select_hvgs_ranked(adata, n_top_genes, flavor="seurat_v3"):
    """
    Runs HVG selection once and stores per-gene ranks, so multiple
    n_top_genes subsets can later be sliced out correctly via
    top_n_hvg_genes(). Use this only when several HVG sizes are needed
    from the same regime (n_top_genes should be the largest size planned).
    """
    hvg_df = sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=n_top_genes,
        flavor=flavor,
        inplace=False,
        check_values=True,
    )
    adata.var["hvg_rank"] = hvg_df["highly_variable_rank"].values
    return adata
 
 
def top_n_hvg_genes(adata, n):
    """Recover the true top-n genes by sorting on the stored rank."""
    ranks = adata.var["hvg_rank"].dropna().sort_values()
    if len(ranks) < n:
        raise ValueError(
            f"Only {len(ranks)} genes have a stored rank, but {n} were "
            f"requested. Re-run select_hvgs_ranked with a larger n_top_genes."
        )
    return ranks.index[:n]
 
 
# ---------------------------------------------------------------------------
# PCA -> (optional Harmony) -> neighbors -> Leiden, for one gene set
# ---------------------------------------------------------------------------
def run_downstream_for_gene_set(
    adata, genes, key_suffix, resolutions, batch_key=None, use_harmony=False, n_pcs=50
):
    sub = adata[:, adata.var_names.isin(genes)].copy()
 
    sc.pp.scale(sub, max_value=10)
    n_comps = min(n_pcs, sub.n_vars - 1, sub.n_obs - 1)
    sc.pp.pca(sub, n_comps=n_comps, svd_solver="arpack")
 
    pca_key = f"X_pca_{key_suffix}"
    adata.obsm[pca_key] = sub.obsm["X_pca"]
 
    rep_key = pca_key
    if use_harmony:
        if batch_key is None:
            raise ValueError("use_harmony=True requires a batch_key.")
        sc.external.pp.harmony_integrate(
            sub, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony"
        )
        rep_key = f"X_pca_harmony_{key_suffix}"
        adata.obsm[rep_key] = sub.obsm["X_pca_harmony"]
 
    neighbors_key = f"neighbors_{key_suffix}"
    sc.pp.neighbors(
        adata, n_pcs=adata.obsm[rep_key].shape[1], use_rep=rep_key, key_added=neighbors_key
    )
 
    for r in resolutions:
        sc.tl.leiden(
            adata, resolution=r, key_added=f"leiden_res_{r}_{key_suffix}",
            neighbors_key=neighbors_key,
        )
 
    return adata
 
 
# ---------------------------------------------------------------------------
# Shared setup: gene standardization + counts vaulting + normalize/log
# ---------------------------------------------------------------------------
def base_preprocessing(adata):
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var_names = bt.Gene.standardize(adata.var_names, organism="human")
    adata.var_names_make_unique()
 
    if "counts" not in adata.layers:
        adata.layers["counts"] = (
            adata.raw.X.copy() if adata.raw is not None else adata.X.copy()
        )
 
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    if not sp.issparse(adata.layers["counts"]):
        adata.layers["counts"] = sp.csr_matrix(adata.layers["counts"])
 
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata
 
 
# ---------------------------------------------------------------------------
# Orchestration for a single view
# ---------------------------------------------------------------------------
BATCH_VIEW_N_HVG = 2000
BATCH_VIEW_RES = []

def process_view(
    adata, view_name, batch_key, n_hvg_sizes, resolutions, use_harmony, flavor="seurat_v3"
):
    adata = base_preprocessing(adata)
 
    if batch_key is not None:
        # Batch view: only ever needs one fixed HVG size, so skip the
        # rank-storage machinery entirely and go straight to downstream.
        genes = compute_hvgs(adata, n_top_genes=BATCH_VIEW_N_HVG, batch_key=batch_key, flavor=flavor)
        key_suffix = f"{view_name}_batch_hvg{BATCH_VIEW_N_HVG}"
        adata = run_downstream_for_gene_set(
            adata, genes=genes, key_suffix=key_suffix, resolutions=BATCH_VIEW_RES,
            batch_key=batch_key, use_harmony=use_harmony,
        )
        return adata
 
    # Standard (non-batch) view: compute ranks once, slice multiple HVG sizes.
    adata = select_hvgs_ranked(adata, n_top_genes=max(n_hvg_sizes), flavor=flavor)
    for n in n_hvg_sizes:
        genes = top_n_hvg_genes(adata, n=n)
        key_suffix = f"{view_name}_nobatch_hvg{n}"
        adata = run_downstream_for_gene_set(
            adata, genes=genes, key_suffix=key_suffix, resolutions=resolutions,
            batch_key=None, use_harmony=False,
        )
 
    return adata
 
 
# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(config_path="datasets.json", base_path="data", output_dir="data",
         n_hvg_sizes=(3000, 2000, 1000), use_harmony=True):
    
    # Get the project root
    project_root = Path(__file__).resolve().parents[2]

    if config_path is not None:
        # If config_path is already absolute, the '/' operator ignores project_root
        config_path = project_root / config_path
    
    base_path = project_root / base_path
    output_dir = project_root / output_dir
 
    os.makedirs(output_dir, exist_ok=True)
 
    with open(config_path) as f:
        config = json.load(f)["datasets"]
 
    for ds_name, ds_info in config.items():
        input_filename = ds_info.get("file_name")
        if not input_filename:
            print(f"Skipping {ds_name}: No primary file_name specified.")
            continue
 
        input_file_path = os.path.join(base_path, input_filename)
        file_name_no_ext = os.path.splitext(input_filename)[0]
        raw_h5ad_path = os.path.join(output_dir, f"{file_name_no_ext}_raw.h5ad")
        sample_col = ds_info.get("columns", {}).get("sample", "Sample")
        batch_col = ds_info.get("columns", {}).get("batch", sample_col)
 
        # --- load / convert once, reused across views -----------------
        if input_filename.endswith(".rds"):
            convert_rds_to_raw_h5ad_r(input_file_path, raw_h5ad_path)
            adata_full = sc.read_h5ad(raw_h5ad_path)
        elif input_filename.endswith(".h5ad"):
            adata_full = sc.read_h5ad(input_file_path)
        else:
            print(f"  -> Unsupported file format for {input_filename}")
            continue

        if sample_col in adata_full.obs.columns:
            adata_full.obs["Sample"] = [
                f"g{s}" if re.match(r"^\d", str(s)) else str(s).replace("-", "_")
                for s in adata_full.obs[sample_col]
            ]
        else:
            raise ValueError(f"Cannot find {sample_col} in {input_file_path}")
 
        views = ds_info.get("views") or {"default": {"subset_vars": {}}}
 
        for view_name, view_info in views.items():
            processed_file_path = os.path.join(
                output_dir, f"{file_name_no_ext}_{view_name}_ECODAprocessed.h5ad"
            )
            if os.path.exists(processed_file_path):
                print(f"Already processed: {ds_name} / {view_name}")
                continue
 
            is_batch_view = view_name == "batch_effect_analysis" and ds_info.get(
                "use_for_batch_effect", False
            )
            batch_key = batch_col if is_batch_view else None
 
            print(f"Processing {ds_name} / {view_name} (batch_key={batch_key})...")
 
            adata_view = apply_subset_vars(adata_full, view_info.get("subset_vars", {}))
            adata_view = process_view(
                adata_view,
                view_name=view_name,
                batch_key=batch_key,
                n_hvg_sizes=n_hvg_sizes,
                resolutions=(0.1, 0.4, 2, 5, 20),
                use_harmony=use_harmony and batch_key is not None,
            )
 
            adata_view.write_h5ad(processed_file_path)
            print(f"  -> Saved: {processed_file_path}\n")


if __name__ == "__main__":
    main()