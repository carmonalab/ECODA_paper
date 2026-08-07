import os
import sys
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path
import pandas as pd
import re

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.gene_utils import standardize_gene_symbols
from src.datasets_io import read_datasets_json
from src.utils.preprocess_utils import load_input, apply_subset_vars


# ---------------------------------------------------------------------------
# HVG selection
# ---------------------------------------------------------------------------
def select_hvgs_ranked(adata, n_top_genes, batch_key=None, flavor="seurat_v3_paper"):
    """
    Runs HVG selection once (optionally batch-aware) and stores per-gene ranks,
    so multiple n_top_genes subsets can later be sliced out correctly via
    top_n_hvg_genes(). Use this only when several HVG sizes are needed
    from the same regime (n_top_genes should be the largest size planned).

    batch_key: batch column for batch-aware HVG selection (benchmark views:
    the standardized sample column; batch-effect views: the dataset's batch_col).
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
# PCA -> Harmony -> neighbors -> Leiden, for one gene set
# ---------------------------------------------------------------------------
def compute_pca_and_store(adata, genes, key_suffix, n_pcs=50):
    """
    Subset to the gene set, scale, run PCA and store the embedding as
    adata.obsm[f"X_pca_{key_suffix}"]. Runs for every HVG size.
    Returns the subsetted object (used by compute_harmony_and_store).
    """
    sub = adata[:, adata.var_names.isin(genes)].copy()

    sc.pp.scale(sub, max_value=10)
    n_comps = min(n_pcs, sub.n_vars - 1, sub.n_obs - 1)
    sc.pp.pca(sub, n_comps=n_comps, svd_solver="arpack")

    adata.obsm[f"X_pca_{key_suffix}"] = sub.obsm["X_pca"]
    return sub


def compute_harmony_and_store(adata, sub, batch_key, key_suffix):
    """
    Harmony-integrate the PCA embedding of `sub` (in place) and store the
    corrected embedding as adata.obsm[f"X_pca_harmony_{key_suffix}"].
    """
    if batch_key is None:
        raise ValueError("Harmony integration requires a batch_key.")
    import harmonypy

    x = sub.obsm["X_pca"].astype(np.float64)
    harmony_out = harmonypy.run_harmony(x, sub.obs, batch_key)
    z_corr = np.asarray(harmony_out.Z_corr)
    # harmonypy <2 returns Z_corr as (d, N) [PCs x cells]; harmonypy >=2 returns
    # (N, d) [cells x PCs]. scanpy's wrapper assumes the former and applies .T,
    # which breaks with 2.x — call harmonypy directly and normalize here.
    if z_corr.shape[0] != sub.n_obs:
        z_corr = z_corr.T
    sub.obsm["X_pca_harmony"] = z_corr
    adata.obsm[f"X_pca_harmony_{key_suffix}"] = z_corr


def run_clustering(adata, rep_key, key_suffix, resolutions):
    """
    Neighbors graph + Leiden clustering at the given resolutions, stored with
    keys neighbors_{key_suffix} and leiden_res_{r}_{key_suffix}.
    """
    neighbors_key = f"neighbors_{key_suffix}"
    sc.pp.neighbors(
        adata, n_pcs=adata.obsm[rep_key].shape[1], use_rep=rep_key, key_added=neighbors_key
    )

    for r in resolutions:
        sc.tl.leiden(
            adata, resolution=r, key_added=f"leiden_res_{r}_{key_suffix}",
            neighbors_key=neighbors_key,
        )


# ---------------------------------------------------------------------------
# Shared setup: gene standardization + counts vaulting + normalize/log
# ---------------------------------------------------------------------------
def base_preprocessing(adata):
    # anndataR-produced inputs (e.g. the _debug h5ad) can have X=None with the
    # raw counts in layers["counts"]; promote the counts layer to X before any
    # X-dependent step (filter_cells etc.).
    if adata.X is None:
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
        elif adata.raw is not None:
            adata.X = adata.raw.X.copy()
        else:
            raise ValueError(
                "Input has neither X, nor a counts layer, nor raw counts; cannot preprocess."
            )

    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    standardize_gene_symbols(adata)
    adata.var_names_make_unique()
 
    if "counts" not in adata.layers:
        adata.layers["counts"] = (
            adata.raw.X.copy() if adata.raw is not None else adata.X.copy()
        )

    # Force CSR unconditionally (not only for dense inputs): the on-disk sparse
    # format is preserved at write time, and backed-mode per-sample subsets in
    # cell type annotation are only selective for CSR (CSC falls back to a full
    # in-memory read per subset -> OOM). tocsr() on an already-CSR matrix is a
    # no-op (no copy). scanpy ops after this preserve CSR.
    adata.X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    adata.layers["counts"] = (
        adata.layers["counts"].tocsr()
        if sp.issparse(adata.layers["counts"])
        else sp.csr_matrix(adata.layers["counts"])
    )

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata


# ---------------------------------------------------------------------------
# Orchestration for a single view
# ---------------------------------------------------------------------------
# Both view types get the same treatment; only the HVG batch_key differs
# (benchmark: standardized sample column; batch-effect: the dataset's batch_col).
# X_pca is stored for every HVG size; Harmony integration and unsupervised
# clustering (neighbors + Leiden) run ONLY at the CLUSTER_N_HVG pass, on both
# the uncorrected and the harmony-corrected embeddings.
CLUSTER_N_HVG = 2000
BATCH_VIEW_N_HVG = 2000          # single HVG size used for batch-effect views
BENCHMARK_VIEW_N_HVG_SIZES = (3000, 2000, 1000)
RESOLUTIONS = (0.1, 0.4, 2, 5, 20, 50)


def process_view(
    adata, view_name, batch_key, n_hvg_sizes, resolutions, flavor="seurat_v3_paper"
):
    """
    Unified per-view pipeline. Resulting keys (example, benchmark view):
      - X_pca_{view_name}_hvg{n}                     for every HVG size n
      - at the CLUSTER_N_HVG (2000) pass additionally:
          X_pca_harmony_{view_name}_hvg2000
          neighbors_{view_name}_hvg2000              (+ _harmony variant)
          leiden_res_{r}_{view_name}_hvg2000         (+ _harmony variant)
    """
    adata = base_preprocessing(adata)

    adata = select_hvgs_ranked(
        adata, n_top_genes=max(n_hvg_sizes), batch_key=batch_key, flavor=flavor
    )

    for n in n_hvg_sizes:
        genes = top_n_hvg_genes(adata, n=n)
        key_suffix = f"{view_name}_hvg{n}"

        sub = compute_pca_and_store(adata, genes, key_suffix)

        if n == CLUSTER_N_HVG:
            run_clustering(adata, f"X_pca_{key_suffix}", key_suffix, resolutions)
            compute_harmony_and_store(adata, sub, batch_key, key_suffix)
            run_clustering(
                adata, f"X_pca_harmony_{key_suffix}", f"{key_suffix}_harmony", resolutions
            )

    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(config_path, input_dir, output_dir, ds_name=None, force=False):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    config = read_datasets_json(config_path)

    for current_ds, entry in config.items():
        if ds_name is not None and current_ds != ds_name:
            continue

        sample_col = entry["sample_col"]
        batch_col = entry["batch_col"] or entry["sample_col"]
        use_for_batch_effect = entry["use_for_batch_effect"]

        views = entry["views"]
        if not views:
            print(f"Skipping {current_ds}: No views defined.")
            continue

        for view_name, view_info in views.items():
            input_file_name = view_info.get("input_file")
            if not input_file_name:
                print(f"Skipping {current_ds} / {view_name}: No input_file_name.")
                continue

            output_file_name = view_info.get("output_file")
            if not output_file_name:
                print(f"Skipping {current_ds} / {view_name}: No output_file_name.")
                continue

            processed_file_path = output_dir / output_file_name
            if processed_file_path.exists() and not force:
                print(f"Already processed: {current_ds} / {view_name}")
                continue

            print(f"Loading {current_ds} / {view_name} ...")
            adata_full = load_input(input_file_name, input_dir, output_dir)

            # Subset on ORIGINAL sample/label values, before the sample column
            # is standardized (standardization can alter '-' or leading digits
            # that subset values may contain).
            adata_view = apply_subset_vars(adata_full, view_info.get("subset_vars", {}))
            if adata_view.n_obs == 0:
                raise ValueError(
                    f"Subset for {current_ds} / {view_name} is empty after "
                    f"apply_subset_vars. Check subset_vars: {view_info.get('subset_vars', {})}"
                )

            if sample_col in adata_view.obs.columns:
                sample_col_out = os.environ.get("SAMPLE_COLNAME", "Sample")
                adata_view.obs[sample_col_out] = [
                    f"g{s}" if re.match(r"^\d", str(s)) else str(s).replace("-", "_")
                    for s in adata_view.obs[sample_col]
                ]
            else:
                raise ValueError(f"Cannot find {sample_col} in obs for {current_ds} / {view_name}")

            is_batch_view = view_name == "batch_effect_analysis" and use_for_batch_effect
            if is_batch_view:
                batch_key = batch_col
                n_hvg_sizes = (BATCH_VIEW_N_HVG,)
            else:
                batch_key = os.environ.get("SAMPLE_COLNAME", "Sample")
                n_hvg_sizes = BENCHMARK_VIEW_N_HVG_SIZES

            if is_batch_view and batch_col not in adata_view.obs.columns:
                raise ValueError(
                    f"batch_col '{batch_col}' not found in obs for {current_ds} / {view_name}. "
                    f"Available columns: {list(adata_view.obs.columns)}"
                )

            print(f"Processing {current_ds} / {view_name} (batch_key={batch_key})...")

            adata_view = process_view(
                adata_view,
                view_name=view_name,
                batch_key=batch_key,
                n_hvg_sizes=n_hvg_sizes,
                resolutions=RESOLUTIONS,
            )

            adata_view.write_h5ad(str(processed_file_path))
            print(f"  -> Saved: {processed_file_path}\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Standardized scRNA-seq preprocessing")
    parser.add_argument("--config_path", required=True,
                        help="Path to datasets.json")
    parser.add_argument("--input_dir", required=True,
                        help="Directory with raw input files")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for processed .h5ad files")
    parser.add_argument("--ds_name", default=None,
                        help="Only process this dataset (skip all others; default: all datasets)")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Recompute views whose output .h5ad already exists "
                             "(bypasses the 'Already processed' skip; needed for "
                             "debug re-runs or after code changes)")
    args = parser.parse_args()
    main(config_path=args.config_path, input_dir=args.input_dir,
         output_dir=args.output_dir, ds_name=args.ds_name, force=args.force)
