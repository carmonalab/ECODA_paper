"""Python benchmark methods (MrVI, scPoli, PILOT) as a CLI script.

Replaces the logic of the archived notebook
`1.2_benchmark_methods_py.qmd` (kept as reference; do NOT delete). Consumes
the preprocessed benchmark view h5ad produced by
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`:

- PILOT consumes the stored obsm embedding `X_pca_{view}_hvg{n}` (the qmd
  recomputed PCA from scratch);
- MrVI/scPoli subset genes via the stored `var["hvg_rank"]` (computed by
  `select_hvgs_ranked`, batch-aware) instead of re-running HVG selection —
  subset to HVGs FIRST, then point X at the raw counts layer
  (`layers["counts"]`; X is log-normalized): MrVI keeps the sparse counts,
  scPoli densifies the small HVG subset to float32
  (`layers["counts"].toarray().astype("float32", copy=False)`);
- cell type annotation columns come from datasets.json
  (`cell_type_low_res` / `cell_type_high_res`).

Feather naming, method-string format and data layout are preserved exactly
from the qmd (the R ingest functions `process_mrvi_fig` /
`process_scpoli_fig` / `process_pilot_fig`, `constants.R` label map and the
notebook recodes depend on them): plain `DataFrame.to_feather()` with the
pandas index (sample names) kept — the index is written as the last feather
column, matching R's `column_to_rownames(ncol)`.

Execution time (float seconds, method body only — excluding h5ad loading,
like the qmd and R `exec_time()`) and peak RSS (`mem_GB`) are appended to a
per-task log feather. One process per task writes the file, so no concurrency
issues. Combos run defaults-first (MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
PILOT_hvg2000_highres) so the main-method rows are measured before any
in-process memory bloat (peak RSS is monotonic within a process).
"""

import argparse
import gc
import os
import resource
import sys
import time
from pathlib import Path

# This script sits one level deeper than 1.1.1_preprocess.py, so the repo
# root is parents[3] (preprocess uses parents[2] from src/3_scrnaseq_preprocessing/).
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
from scvi.external import MRVI
import pilotpy as pl
import torch

from src.utils.py.datasets_io import read_datasets_json


# scPoli is imported lazily (get_scpoli) so that MrVI/PILOT runs never touch
# scarches: scarches 0.6.1 (its final release) does `from anndata import
# AnnData, read`, but `anndata.read` was removed in anndata >= 0.12 (the
# pinned 0.12.19). The shim below restores it as the documented alias
# `read_h5ad` — scarches only calls `read()` on .h5ad files, so this is a
# faithful drop-in. See https://github.com/theislab/scarches.
def get_scpoli():
    import anndata as ad

    if not hasattr(ad, "read"):
        ad.read = ad.read_h5ad  # scarches compat (anndata >= 0.12 removed `read`)
    from scarches.models.scpoli import scPoli

    return scPoli


# ---------------------------------------------------------------------------
# Execution time / memory logging
# ---------------------------------------------------------------------------
def peak_rss_gb():
    """Peak resident set size of this process in GB.

    getrusage().ru_maxrss units: KB on Linux, bytes on macOS.
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / 1024.0 / 1024.0 / 1024.0
    return rss / 1024.0 / 1024.0


def log_execution_time(dataset_name, method_str, time_secs, log_file):
    """Append/overwrite one (dataset, method) row in the per-task exec log.

    Read-modify-write on the feather (single process per task). Overwrites
    the row if the (dataset, method) combo already exists (matches the qmd's
    rerun semantics).
    """
    new_row = pd.DataFrame(
        {
            "dataset": [dataset_name],
            "method": [method_str],
            "time_secs": [float(time_secs)],
            "mem_GB": [peak_rss_gb()],
        }
    )
    if os.path.exists(log_file):
        df_existing = pd.read_feather(log_file)
        mask = (df_existing["dataset"] == dataset_name) & (
            df_existing["method"] == method_str
        )
        if mask.any():
            df_final = df_existing[~mask]
        else:
            df_final = df_existing
        df_final = pd.concat([df_final, new_row], ignore_index=True)
    else:
        df_final = new_row
    df_final.reset_index(drop=True).to_feather(log_file)


# ---------------------------------------------------------------------------
# Combo resolution (legacy qmd rules, datasets.json-driven)
# ---------------------------------------------------------------------------
def scpoli_dims_for(n_hvg, res_label):
    """scPoli embedding dims for an (n_hvg, resolution) combo (qmd rules)."""
    if res_label == "_highres":
        return [2, 3, 5, 10, 15] if n_hvg == 2000 else [15]
    if res_label == "_lowres" and n_hvg == 2000:
        return [15]
    return []


def run_pilot_for(n_hvg, res_label):
    """Whether PILOT runs for an (n_hvg, resolution) combo (qmd rules)."""
    if res_label == "_highres":
        return True
    return res_label == "_lowres" and n_hvg == 2000


# Default (main-method) combos — constants.R method_label_map_main and the
# notebook's exec-time figure: MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
# PILOT_hvg2000_highres.
DEFAULT_HVG = 2000
DEFAULT_SCPOLI_DIM = 15
DEFAULT_RES_LABEL = "_highres"


def is_default_combo(method, combo):
    n, res_label, _, payload, _, _ = combo
    if method == "mrvi":
        return n == DEFAULT_HVG
    if method == "scpoli":
        return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL and payload == DEFAULT_SCPOLI_DIM
    if method == "pilot":
        return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL
    return False


def top_n_hvg_genes(adata, n):
    """Top-n genes from the stored hvg_rank (set by 1.1.1_preprocess.py)."""
    ranks = adata.var["hvg_rank"].dropna().sort_values()
    if len(ranks) < n:
        raise ValueError(
            f"Only {len(ranks)} genes have a stored hvg_rank, but {n} were "
            f"requested. Re-run preprocessing with a larger HVG size."
        )
    return list(ranks.index[:n])


def use_counts_layer(sub, method, ds_name):
    """Point X at the raw counts layer for MrVI/scPoli.

    The preprocessed h5ad has log-normalized X; both models need the raw
    counts (vaulted by base_preprocessing as layers["counts"]). scPoli
    additionally needs dense float32 — call this AFTER the HVG subset so
    only the small n_obs x n_hvg matrix is densified. Falls back to X with
    a warning if the counts layer is missing.
    """
    counts = sub.layers.get("counts")
    if counts is None:
        print(f"WARNING: {ds_name}: no 'counts' layer; using log-normalized X "
              f"for {method}.")
        counts = sub.X
    if method == "scpoli":
        if sp.issparse(counts):
            sub.X = counts.toarray().astype("float32", copy=False)
        else:
            sub.X = np.asarray(counts, dtype="float32")
    else:
        sub.X = counts
    return sub


# ---------------------------------------------------------------------------
# Method bodies (qmd semantics preserved)
# ---------------------------------------------------------------------------
def run_mrvi(adata, device, output_path):
    """MrVI local sample distances (lowres only; ct column unused)."""
    adata.obs["dummy_col"] = np.zeros(adata.n_obs)
    MRVI.setup_anndata(adata, sample_key="Sample")
    model = MRVI(adata)
    model.train(max_epochs=50, accelerator=device)
    dists = model.get_local_sample_distances(
        keep_cell=False, groupby="dummy_col", batch_size=32
    )
    df_dists = dists["dummy_col"].isel(dummy_col_name=0).to_pandas()
    df_dists.to_feather(output_path)


def run_scpoli(adata, ct_col, dim, output_path):
    """scPoli conditional sample embeddings for one embedding dim."""
    # scPoli requires a cell-type label for EVERY cell, but datasets whose
    # declared ct columns are the pipeline annotation columns (Lee/Zhang:
    # layer1/layer2) carry NaN for cells the annotators left unclassified
    # (HiTME covers immune subsets only, ~66-88% of cells). scarches'
    # label_encoder runs np.unique() on the mixed str/NaN object column,
    # which crashes under numpy 2.x ("'<' not supported between instances of
    # 'str' and 'float'"). Fill with an explicit "Unknown" class: every cell
    # stays in the output (per-cell embedding rows stay aligned) and
    # unannotated cells get their own scPoli prototype. No-op on complete
    # columns (all other datasets use fully-covered author annotations).
    n_na = int(adata.obs[ct_col].isna().sum())
    if n_na:
        print(f"Filling {n_na}/{adata.n_obs} missing values in '{ct_col}' "
              "with 'Unknown' (scPoli requires a label per cell).")
        col = adata.obs[ct_col]
        # anndata stores string obs columns as pandas Categorical by default;
        # fillna() cannot introduce an undefined category, so register
        # 'Unknown' first.
        if isinstance(col.dtype, pd.CategoricalDtype):
            col = col.cat.add_categories("Unknown")
        adata.obs[ct_col] = col.fillna("Unknown")
    scPoli = get_scpoli()
    scpoli_model = scPoli(
        adata=adata,
        condition_keys="Sample",
        cell_type_keys=ct_col,
        embedding_dims=dim,
        recon_loss="nb",
    )
    scpoli_model.train(
        n_epochs=50,
        pretraining_epochs=40,
        early_stopping_kwargs={
            "early_stopping_metric": "val_prototype_loss",
            "mode": "min",
            "threshold": 0,
            "patience": 20,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        },
        eta=5,
    )
    adata_emb = scpoli_model.get_conditional_embeddings()
    df_embs = pd.DataFrame(adata_emb.X, index=adata_emb.obs_names)
    df_embs.columns = [f"Dim_{i + 1}" for i in range(df_embs.shape[1])]
    df_embs.to_feather(output_path)


def run_pilot(adata, ct_col, view, n_hvg, output_path):
    """PILOT Wasserstein sample distances on the preprocessed obsm PCA."""
    emb_key = f"X_pca_{view}_hvg{n_hvg}"
    emb = adata.obsm[emb_key]
    # PILOT (pilotpy>=2.0.x) requires a named-columns pandas DataFrame in
    # obsm (Trajectory.extract_data_anno_scRNA_from_h5ad accesses .columns);
    # the preprocess step stores scanpy's plain ndarray instead.
    if not hasattr(emb, "columns"):
        emb = pd.DataFrame(
            emb,
            index=adata.obs_names,
            columns=[f"PCA_{i + 1}" for i in range(emb.shape[1])],
        )
    adata.obsm[emb_key] = emb
    pl.tl.wasserstein_distance(
        adata,
        emb_matrix=emb_key,
        clusters_col=ct_col,
        sample_col="Sample",
        status="Sample",
    )
    adata.uns["EMD_df"].to_feather(output_path)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def process_dataset(args, ds_name, entry):
    """Run all combos of the requested method for one dataset.

    Loads the h5ad once per task (not per combo); skips combos whose output
    feather already exists unless --force.
    """
    view_name = args.view
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if view_name not in entry["views"]:
        raise ValueError(
            f"Dataset '{ds_name}' has no '{view_name}' view in datasets.json."
        )
    view_output = entry["views"][view_name]["output_file"]

    input_path = Path(args.input_dir) / view_output
    if not input_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_path}")

    lowres_col = entry.get("cell_type_low_res")
    highres_col = entry.get("cell_type_high_res")

    # combo: (n_hvg, res_label, ct_col, payload, run_fn, out_name)
    combos = []

    if args.method == "mrvi":
        if lowres_col is None:
            print(f"WARNING: {ds_name}: cell_type_low_res is None; skipping MrVI.")
            return
        for n in args.hvg:
            out_name = f"{ds_name}_hvg{n}_mrvi_dists.feather"
            combos.append((n, "_lowres", None, None, run_mrvi, out_name))

    elif args.method == "scpoli":
        for res_label, ct_col in (("_lowres", lowres_col), ("_highres", highres_col)):
            if ct_col is None:
                continue
            for n in args.hvg:
                for dim in scpoli_dims_for(n, res_label):
                    out_name = (
                        f"{ds_name}_hvg{n}{res_label}_scpoli_dims{dim}_embs.feather"
                    )
                    combos.append((n, res_label, ct_col, dim, run_scpoli, out_name))

    elif args.method == "pilot":
        for res_label, ct_col in (("_lowres", lowres_col), ("_highres", highres_col)):
            if ct_col is None:
                continue
            for n in args.hvg:
                if run_pilot_for(n, res_label):
                    out_name = f"{ds_name}_hvg{n}{res_label}_pilot_dists.feather"
                    combos.append((n, res_label, ct_col, None, run_pilot, out_name))

    # Defaults-first ordering: ru_maxrss peak RSS is monotonic within a
    # process, so combos run earlier report the least bloated mem_GB (memory
    # leaks / allocator retention from earlier combos would otherwise inflate
    # the defaults' rows). Stable sort: non-default combos keep their order.
    # If --hvg excludes the default size, no default combo exists and the
    # sort is a no-op — no behavior change.
    combos.sort(key=lambda c: 0 if is_default_combo(args.method, c) else 1)

    pending = []
    for n, res_label, ct_col, payload, run_fn, out_name in combos:
        out_path = output_dir / out_name
        if out_path.exists() and not args.force:
            print(f"Already processed: {out_name}")
            continue
        pending.append((n, res_label, ct_col, payload, run_fn, out_path))

    if not pending:
        return

    print(f"Loading {input_path} ...")
    adata = sc.read_h5ad(str(input_path))
    adata = adata.to_memory()

    if "Sample" not in adata.obs.columns:
        raise ValueError(
            f"Cannot find standardized sample column 'Sample' in obs of {input_path}."
        )

    for n, res_label, ct_col, payload, run_fn, out_path in pending:
        # HVG gene subset per combo, done FIRST (before any dense conversion)
        # so only the small n_obs x n_hvg matrix is materialized. PILOT uses
        # the stored obsm directly, no subset needed.
        if args.method in ("mrvi", "scpoli"):
            genes = top_n_hvg_genes(adata, n)
            sub = adata[:, genes].copy()
            sub = use_counts_layer(sub, args.method, ds_name)
        else:
            sub = adata

        if ct_col is not None and ct_col not in sub.obs.columns:
            raise ValueError(
                f"Cell type column '{ct_col}' not found in obs of {ds_name} "
                f"(available: {list(sub.obs.columns)})."
            )

        # Exact legacy method strings (constants.R + notebook recodes depend
        # on them): MrVI_hvg{n}, scPoli_hvg{n}_dims{d}{res}, PILOT_hvg{n}{res}.
        if args.method == "mrvi":
            method_str = f"MrVI_hvg{n}"
        elif args.method == "scpoli":
            method_str = f"scPoli_hvg{n}_dims{payload}{res_label}"
        else:
            method_str = f"PILOT_hvg{n}{res_label}"

        print(f"Processing {method_str} ...")
        start_time = time.time()
        if args.method == "mrvi":
            run_mrvi(sub, args.device, out_path)
        elif args.method == "scpoli":
            run_scpoli(sub, ct_col, payload, out_path)
        else:
            run_pilot(sub, ct_col, args.view, n, out_path)
        exec_time = time.time() - start_time

        log_execution_time(ds_name, method_str, exec_time, args.log_file)
        print(f"  -> Saved: {out_path} ({exec_time:.2f}s, "
              f"{peak_rss_gb():.2f} GB peak RSS)")
        gc.collect()


def main():
    parser = argparse.ArgumentParser(
        description="Run Python benchmark methods (MrVI/scPoli/PILOT) on a "
                    "preprocessed benchmark view h5ad."
    )
    parser.add_argument("--config_path", required=True,
                        help="Path to datasets.json")
    parser.add_argument("--ds_name", required=True,
                        help="Dataset key in datasets.json")
    parser.add_argument("--view", default="benchmark_analysis",
                        help="View name (default: benchmark_analysis)")
    parser.add_argument("--method", required=True,
                        choices=["mrvi", "scpoli", "pilot"],
                        help="Benchmark method to run")
    parser.add_argument("--input_dir", required=True,
                        help="Directory holding the preprocessed view h5ad")
    parser.add_argument("--output_dir", required=True,
                        help="Feather output dir (created if missing)")
    parser.add_argument("--log_file", default=None,
                        help="Per-task execution-time log feather "
                             "(default: <output_dir>/execution_times.feather "
                             "for local runs)")
    parser.add_argument("--hvg", nargs="+", type=int, default=[1000, 2000, 3000],
                        help="HVG sizes to run (default: 1000 2000 3000)")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Recompute combos whose output feather already exists")
    parser.add_argument("--device", default="auto",
                        choices=["auto", "cpu", "cuda"],
                        help="Train accelerator (default: auto; uses the "
                             "allocated GPU on shared-gpu nodes)")
    args = parser.parse_args()

    args.hvg = sorted(set(args.hvg))
    if args.log_file is None:
        args.log_file = os.path.join(args.output_dir, "execution_times.feather")

    scvi.settings.seed = 0
    print("scvi-tools version:", scvi.__version__)
    print("torch.cuda.is_available():", torch.cuda.is_available())

    config = read_datasets_json(args.config_path, view=args.view)
    if args.ds_name not in config:
        raise ValueError(f"'{args.ds_name}' is not a dataset in {args.config_path} "
                         f"with a '{args.view}' view.")
    process_dataset(args, args.ds_name, config[args.ds_name])
    print("Processing complete!")


if __name__ == "__main__":
    main()
