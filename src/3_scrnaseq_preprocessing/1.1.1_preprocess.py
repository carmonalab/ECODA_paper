import hashlib
import os
import sys
import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path
import pandas as pd
import re

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.gene_utils import standardize_gene_symbols
from src.utils.py.datasets_io import read_datasets_json
from src.utils.py.benchmark_h5ad_contract import (
    validate_benchmark_h5ad_contract,
    validate_benchmark_h5ad_path,
)
from src.utils.py.preprocess_utils import (
    load_input,
    apply_subset_vars,
    remove_low_cellcount_samples,
)


def _checksum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_checksum(path):
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    tmp.write_text(
        f"MD5={_checksum(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n"
    )
    os.replace(tmp, sidecar)


def _validate_recorded_checksum(path):
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or path.stat().st_size == 0 or not sidecar.is_file():
        return False
    records = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            records[key] = value
    return (
        records.get("PATH") == str(path)
        and records.get("MD5") == _checksum(path)
        and records.get("SIZE") == str(path.stat().st_size)
    )


def _write_h5ad_atomic(adata, path, view_name):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{path}.md5")
    backup = path.with_name(f".{path.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    had_path = path.is_file()
    had_sidecar = sidecar.is_file()

    validate_benchmark_h5ad_contract(adata, view_name, "preprocessing")
    try:
        adata.write_h5ad(str(tmp))
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"preprocessing produced an empty h5ad: {tmp}")
        validate_benchmark_h5ad_path(str(tmp), view_name, "preprocessing")

        if had_path:
            os.link(path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(tmp, path)
        validate_benchmark_h5ad_path(str(path), view_name, "preprocessing")
        _write_checksum(path)
    except Exception:
        if backup.exists():
            os.replace(backup, path)
        elif not had_path and path.exists():
            path.unlink()
        if sidecar_backup.exists():
            os.replace(sidecar_backup, sidecar)
        elif not had_sidecar and sidecar.exists():
            sidecar.unlink()
        raise
    finally:
        for temporary in (tmp, backup, sidecar_backup):
            if temporary.exists():
                temporary.unlink()


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
    try:
        hvg_df = sc.pp.highly_variable_genes(
            adata,
            layer="counts",
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            inplace=False,
            check_values=True,
        )
    except ValueError as e:
        # seurat_v3_paper fits a loess of log10(var) on log10(mean) per batch;
        # integer counts in small batches quantize log-means into tied
        # x-values -> singular design matrix -> skmisc loess ValueError.
        # Retry once on a deterministic jittered copy of the counts layer.
        print(
            f"seurat_v3 loess fit failed ({e}); retried with deterministic "
            "jittered counts (seed 42)"
        )
        rng = np.random.RandomState(42)
        cnt = adata.layers["counts"].copy()
        cnt.data = (
            cnt.data + rng.random_sample(cnt.data.shape).astype(np.float32) * 1e-6
        )
        adata.layers["counts_jittered"] = cnt
        try:
            hvg_df = sc.pp.highly_variable_genes(
                adata,
                layer="counts_jittered",
                n_top_genes=n_top_genes,
                flavor=flavor,
                batch_key=batch_key,
                inplace=False,
                check_values=False,
            )
        except ValueError:
            del adata.layers["counts_jittered"]
            raise e
        del adata.layers["counts_jittered"]
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
            flavor="igraph", n_iterations=2, directed=False,
        )


# ---------------------------------------------------------------------------
# Shared setup: gene standardization + counts vaulting + normalize/log
# ---------------------------------------------------------------------------
def base_preprocessing(adata):
    """
    Extracts raw integer counts into adata.X upfront, standardizes gene symbols,
    filters cells and genes on raw counts, vaults raw counts into layers['counts'],
    and computes log-normalized expression in adata.X.
    """
    counts_layer_present = "counts" in adata.layers
    if counts_layer_present:
        # Use the existing raw layer as the temporary filter matrix. The copy
        # needed for normalization is deferred until after cell/gene filters,
        # avoiding two full-cohort sparse duplicates for large inputs.
        adata.X = adata.layers["counts"]
    elif adata.raw is not None and adata.raw.X is not None:
        raw_mat = adata.raw.X
        sample_raw = raw_mat.data[:1000] if sp.issparse(raw_mat) else raw_mat.ravel()[:1000]
        sample_raw_pos = sample_raw[sample_raw > 1e-6]
        is_raw_int = bool(sample_raw_pos.size > 0 and np.all(np.abs(sample_raw_pos - np.round(sample_raw_pos)) < 1e-3))
        if is_raw_int:
            raw_var = (
                adata.raw.var.copy()
                if hasattr(adata.raw, "var") and len(adata.raw.var) == raw_mat.shape[1]
                else adata.var.copy()
            )
            # Construct fresh AnnData using raw.X and raw.var to keep gene dimensions aligned
            adata = ad.AnnData(
                X=adata.raw.X.copy(),
                obs=adata.obs.copy(),
                var=raw_var,
                obsm=adata.obsm.copy(),
            )
        elif adata.X is None:
            raise ValueError(
                "Input has raw matrix with non-integer values and X is None; cannot preprocess."
            )
    elif adata.X is None:
        raise ValueError(
            "Input has neither X, nor a counts layer, nor raw counts; cannot preprocess."
        )

    # 2. Ensure CSR format on counts
    adata.X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

    # 3. Filter cells and genes on raw counts
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    # 4. Standardize gene symbols on filtered counts
    standardize_gene_symbols(adata)
    adata.var_names_make_unique()

    if counts_layer_present:
        # Keep layers["counts"] raw and detach X only after filtering.
        adata.X = adata.X.copy()
    else:
        # Raw counts came from X/raw.X; vault one filtered copy before
        # normalizing X.
        adata.layers["counts"] = adata.X.copy()
    # 6. Normalize and log-transform X for PCA/Harmony
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
    adata,
    view_name,
    batch_key,
    n_hvg_sizes,
    resolutions,
    flavor="seurat_v3_paper",
    compute_harmony=True,
):
    """Run one explicit analysis view with pass-qualified output keys."""
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
            if compute_harmony:
                compute_harmony_and_store(adata, sub, batch_key, key_suffix)
                run_clustering(
                    adata,
                    f"X_pca_harmony_{key_suffix}",
                    f"{key_suffix}_harmony",
                    resolutions,
                )

    return adata


def main(config_path, input_dir, output_dir, ds_name=None, force=False, view=None):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    allowed_views = {
        "benchmark_analysis",
        "batch_effect_uncorrected",
        "batch_effect_corrected",
    }
    if view is not None and view not in allowed_views:
        raise ValueError(f"Unknown preprocessing view: {view}")
    config = read_datasets_json(config_path, view=view)
    if view is not None and not any(entry.get("views") for entry in config.values()):
        raise ValueError(f"No dataset declares preprocessing view: {view}")

    for current_ds, entry in config.items():
        if ds_name is not None and current_ds != ds_name:
            continue

        sample_col = entry["sample_col"]
        batch_col = entry.get("batch_col")
        use_for_batch_effect = bool(entry.get("use_for_batch_effect"))
        views = entry.get("views") or {}
        if not views:
            print(f"Skipping {current_ds}: No views defined.")
            continue

        for view_name, view_info in views.items():
            if view_name not in allowed_views:
                raise ValueError(
                    f"Unknown preprocessing view {view_name!r} for dataset {current_ds}"
                )
            is_uncorrected = view_name == "batch_effect_uncorrected"
            is_corrected = view_name == "batch_effect_corrected"
            is_batch_view = is_uncorrected or is_corrected
            if is_batch_view and not use_for_batch_effect:
                raise ValueError(
                    f"Dataset {current_ds} declares {view_name} but "
                    "use_for_batch_effect is false"
                )
            if is_corrected and batch_col is None:
                raise ValueError(
                    "corrected batch-effect view requires a confirmed columns.batch"
                )

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
                try:
                    validate_benchmark_h5ad_path(
                        str(processed_file_path), view_name, "preprocessing"
                    )
                    if _validate_recorded_checksum(processed_file_path):
                        print(f"Already processed and validated: {current_ds} / {view_name}")
                        continue
                except (OSError, ValueError, KeyError) as exc:
                    print(
                        f"Existing artifact invalid for {current_ds} / {view_name}; "
                        f"recomputing ({exc})"
                    )

            print(f"Loading {current_ds} / {view_name} ...")
            adata_full = load_input(input_file_name, input_dir, output_dir)

            # Subset on ORIGINAL sample/label values, before the sample column
            # is standardized (standardization can alter '-' or leading digits).
            adata_view = apply_subset_vars(adata_full, view_info.get("subset_vars", {}))
            if adata_view.n_obs == 0:
                raise ValueError(
                    f"Subset for {current_ds} / {view_name} is empty after "
                    f"apply_subset_vars. Check subset_vars: {view_info.get('subset_vars', {})}"
                )

            if sample_col in adata_view.obs.columns:
                sample_col_out = "Sample" if is_batch_view else os.environ.get("SAMPLE_COLNAME", "Sample")
                adata_view.obs[sample_col_out] = [
                    f"g{s}" if re.match(r"^\d", str(s)) else str(s).replace("-", "_")
                    for s in adata_view.obs[sample_col]
                ]
            else:
                raise ValueError(f"Cannot find {sample_col} in obs for {current_ds} / {view_name}")

            adata_view, removed_samples = remove_low_cellcount_samples(
                adata_view, sample_col=sample_col_out, min_cells_per_sample=500
            )
            removed_summary = (
                ", ".join(
                    f"{sample}={count}" for sample, count in removed_samples.items()
                )
                if removed_samples
                else "none"
            )
            print(
                f"Sample-count filter ({current_ds} / {view_name}): "
                f"threshold=500, removed={len(removed_samples)} [{removed_summary}]"
            )

            if is_uncorrected:
                batch_key = "Sample"
                n_hvg_sizes = (BATCH_VIEW_N_HVG,)
                compute_harmony = False
            elif is_corrected:
                batch_key = batch_col
                n_hvg_sizes = (BATCH_VIEW_N_HVG,)
                compute_harmony = True
            else:
                batch_key = os.environ.get("SAMPLE_COLNAME", "Sample")
                n_hvg_sizes = BENCHMARK_VIEW_N_HVG_SIZES
                compute_harmony = True

            if is_corrected and batch_col not in adata_view.obs.columns:
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
                compute_harmony=compute_harmony,
            )
            _write_h5ad_atomic(adata_view, processed_file_path, view_name)
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
    parser.add_argument("--view", default=None,
                        help="Process exactly one declared view (default: all declared views)")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Recompute views whose output .h5ad already exists "
                             "(bypasses the 'Already processed' skip; needed for "
                             "debug re-runs or after code changes)")
    args = parser.parse_args()
    main(config_path=args.config_path, input_dir=args.input_dir,
         output_dir=args.output_dir, ds_name=args.ds_name, force=args.force,
         view=args.view)
