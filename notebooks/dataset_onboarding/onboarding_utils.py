"""Shared helpers for the dataset-onboarding check notebooks.

Used by the Python ``dataset_check_*.qmd`` notebooks (pixi **default** env,
scanpy 1.12.2 / anndata 0.12.x -- NOT the py-cuda13 env) and by the
standalone ``onboarding_metrics.R`` helper (cell-level LISI separation via
``src/utils/scoring_metrics.R``, run with ``pixi run -e default Rscript``).

Scope: one-off review notebooks for the PILOT-GM-VAE study datasets
downloaded to the NAS folder ``JooM_2025_41097818/output/``. These helpers
are NOT part of the pipeline and must NOT source the R notebook loader
(``load_all_functions.R`` / ``imports.R`` -- those serve only
``benchmark_analysis.rmd`` / ``batch_effect_analysis.rmd``).

Conventions
-----------
* Always read big h5ads with ``sc.read_h5ad(path, backed="r")`` and subsample
  before in-memory work.
* Counts slot chain: ``layers["counts"]`` -> ``X`` -> ``raw.X``. We accept
  float-encoded CSR counts (integer *values*, not dtype).
* Unintegrated space only for metrics/UMAP: never ``X_pca_harmony*`` /
  integrated embeddings.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import scanpy as sc
    import anndata as ad
except ImportError as e:  # pragma: no cover
    raise SystemExit(
        "onboarding_utils.py requires the pixi default env python "
        "(scanpy/anndata available there)"
    ) from e


def locate_counts(adata, verbose: bool = True) -> tuple:
    """Return (counts_matrix, slot_name).

    Slot chain: ``layers["counts"]`` -> ``X`` -> ``raw.X``. Never returns
    ``None``: a missing counts slot raises a clear error so T2 can flag the
    dataset as FAIL and suggest the fallback source.
    """
    for layer_name in ("counts", "counts_cite", "counts_atac"):
        if layer_name in adata.layers:
            if verbose:
                print(f"counts slot: adata.layers[{layer_name!r}]")
            return adata.layers[layer_name], f"layers/{layer_name}"
    if adata.X is not None:
        if verbose:
            print("counts slot: adata.X")
        return adata.X, "X"
    if adata.raw is not None and adata.raw.X is not None:
        if verbose:
            print("counts slot: adata.raw.X")
        return adata.raw.X, "raw.X"
    raise ValueError("No counts found: no 'counts' layer, no X, no raw.X.")


def _sample_values(matrix, n_sample_cells: int = 200_000, seed: int = 0):
    """Draw a flat sample of matrix values (numpy 1-D) for sanity checks."""
    rng = np.random.default_rng(seed)
    n_cells, n_genes = matrix.shape
    total = n_cells * n_genes
    n_pick = min(n_sample_cells, total)
    if total == 0:
        return np.array([], dtype=float)
    if _is_sparse(matrix):
        # Sparse: draw random (row, col) pairs without building a dense copy.
        vals = matrix.data
        if len(vals) >= n_pick:
            return np.asarray(rng.choice(vals, size=n_pick, replace=False), dtype=float)
        return np.asarray(vals, dtype=float)
    flat = np.asarray(matrix).ravel()
    return np.asarray(rng.choice(flat, size=n_pick, replace=False), dtype=float)


def _is_sparse(matrix) -> bool:
    return hasattr(matrix, "toarray") or hasattr(matrix, "tocsr")


def count_sanity_check(
    adata,
    n_sample_cells: int = 200_000,
    tol: float = 1e-6,
    verbose: bool = True,
) -> dict:
    """T2 count sanity check.

    Verdict PASS requires: a counts slot exists, values are non-negative,
    finite, and integer-valued *up to floating-point tolerance* (float-encoded
    CSR like 1.0/2.0 is valid raw counts -- dtype is not checked). Also flags
    X that looks log-normalized while counts live elsewhere.
    """
    try:
        counts, slot = locate_counts(adata, verbose=verbose)
    except ValueError as e:
        res = {"verdict": "FAIL", "reason": str(e), "slot": None}
        if verbose:
            print(json.dumps(res, indent=2))
        return res

    n_cells, n_genes = counts.shape
    sparse = _is_sparse(counts)
    is_csr = hasattr(counts, "format") and counts.format == "csr"
    n_nonzero = counts.nnz if sparse else int(np.count_nonzero(counts))
    sparsity = 1.0 - (n_nonzero / (n_cells * n_genes)) if n_cells * n_genes else float("nan")

    sample = _sample_values(counts, n_sample_cells=n_sample_cells)
    non_negative = bool(np.all(sample >= 0)) if sample.size else True
    finite = bool(np.all(np.isfinite(sample))) if sample.size else True
    # Integer-VALUE check with relative epsilon tolerance: counts encoded as
    # float32 can drift by a few ULP from the exact integer.
    rounded = np.round(sample)
    integer_values = bool(
        sample.size == 0 or np.all(np.abs(sample - rounded) <= tol * np.maximum(1.0, np.abs(sample)))
    )
    max_val = float(sample.max()) if sample.size else float("nan")
    frac_non_integer = float(
        np.mean(np.abs(sample - np.round(sample)) > tol) if sample.size else 0.0
    ) if sample.size else 0.0

    log_norm_flag = None
    if slot != "X" and adata.X is not None and _is_sparse(adata.X):
        x_sample = _sample_values(adata.X, n_sample_cells=n_sample_cells)
        x_max = float(x_sample.max()) if x_sample.size else float("nan")
        x_integer = bool(
            x_sample.size == 0
            or np.all(np.abs(x_sample - np.round(x_sample)) <= tol * np.maximum(1.0, np.abs(x_sample)))
        )
        if not x_integer and x_max < 30:
            log_norm_flag = (
                f"X looks log-normalized (max={x_max:.2f}, fractional values) "
                f"while counts live in {slot}"
            )
            if verbose:
                print("WARNING: " + log_norm_flag)

    checks = {
        "slot": slot,
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "sparse": bool(sparse),
        "is_csr": is_csr,
        "dtype": str(getattr(counts, "dtype", "?")),
        "nnz": int(n_nonzero),
        "sparsity": float(sparsity),
        "non_negative": non_negative,
        "finite": finite,
        "integer_values": integer_values,
        "max_sampled_value": max_val,
        "frac_non_integer_sampled": frac_non_integer,
        "log_normalized_X_flag": log_norm_flag,
    }
    ok = non_negative and finite and integer_values and n_genes > 0 and n_cells > 0
    checks["verdict"] = "PASS" if ok else "FAIL"
    if verbose:
        print(json.dumps(checks, indent=2))
    return checks


def obs_summary(adata, top_n: int = 20) -> pd.DataFrame:
    """Per-column dtype + top value counts (bounded) for obs."""
    rows = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        nunique = int(s.nunique(dropna=True))
        if adata.obs[col].dtype.name in ("category", "object") or nunique <= 30:
            vc = (
                s.astype(str)
                .value_counts(dropna=True)
                .head(top_n)
                .to_dict()
            )
            has_na = bool(s.isna().any())
            rows.append(
                {
                    "column": col,
                    "dtype": str(s.dtype),
                    "n_unique": nunique,
                    "na_count": int(s.isna().sum()),
                    "top_values": json.dumps(vc, default=str)[:400],
                }
            )
        else:
            rows.append(
                {
                    "column": col,
                    "dtype": str(s.dtype),
                    "n_unique": nunique,
                    "na_count": int(s.isna().sum()),
                    "top_values": "",
                }
            )
    return pd.DataFrame(rows)


def candidate_col_detection(
    obs: pd.DataFrame, sample_terms=None, label_terms=None, batch_terms=None, ct_terms=None
) -> dict:
    """Heuristic column-role detection for the study-summary card.

    Purely advisory -- the notebook always lets the curated sets (from the
    paper/feasibility table) override the heuristics.
    """
    sample_terms = sample_terms or ["sample", "specimen", "donor", "patient", "participant", "subject", "individual"]
    label_terms = label_terms or ["status", "disease", "condition", "class", "group", "label", "diagnosis", "origin", "tissue"]
    batch_terms = batch_terms or ["batch", "site", "center", "study", "dataset", "seqtec", "platform", "assay", "chemistry", "library"]
    ct_terms = ct_terms or ["cell", "cluster", "annotation", "subtype", "subclass", "class", "type", "identity", "level", "label"]

    cols = [str(c) for c in obs.columns]
    low = {c: c.lower() for c in cols}

    # Cluster/embedding/QC-like columns are never sample/bio/batch candidates
    # (e.g. preprocessed-leiden columns, UMAP/PCA keys, UCell scores, nCount/
    # nFeature). Cell-type candidates may legitimately contain 'cluster'.
    noise_terms = (
        "leiden", "kmeans", "seurat_clusters", "scvi", "scanorama",
        "harmony", "umap", "pca", "phase", "ucell", "_score", "score_",
        "ncount", "nfeature", "n_genes", "cell_barcode", "barcode",
        "classification_confidence",
    )

    def is_noise(c):
        lc = low[c]
        return any(n in lc for n in noise_terms)

    out = {"sample": [], "label": [], "batch": [], "cell_type": []}
    for c in cols:
        lc = low[c]
        if is_noise(c):
            continue
        if any(t in lc for t in sample_terms) and obs[c].nunique(dropna=True) >= 2:
            out["sample"].append(c)
        if any(t in lc for t in label_terms) and 2 <= obs[c].nunique(dropna=True) <= 50:
            out["label"].append(c)
        if any(t in lc for t in batch_terms) and 2 <= obs[c].nunique(dropna=True) <= 200:
            out["batch"].append(c)
        if any(t in lc for t in ct_terms) and 2 <= obs[c].nunique(dropna=True) <= 200:
            out["cell_type"].append(c)
    # drop exact duplicates across roles (e.g. a plain 'Sample' column matching
    # both sample and batch terms)
    return out


def cells_per_sample_stats(obs: pd.DataFrame, sample_col: str) -> dict:
    vc = obs[sample_col].value_counts()
    stats = {
        "n_samples": int(vc.size),
        "total_cells": int(vc.sum()),
        "min": int(vc.min()) if vc.size else None,
        "median": float(vc.median()) if vc.size else None,
        "mean": float(vc.mean()) if vc.size else None,
        "std": float(vc.std(ddof=0)) if vc.size else None,
        "max": int(vc.max()) if vc.size else None,
        "n_samples_lt_50_cells": int((vc < 50).sum()),
        "cells_per_sample": vc.to_dict(),
    }
    return stats


def paper_table_compare(obs: pd.DataFrame, sample_col: str, ct_col: str | None, expected: dict) -> pd.DataFrame:
    n_cells = len(obs)
    n_samples = int(obs[sample_col].nunique(dropna=True)) if sample_col in obs else None
    n_cts = int(obs[ct_col].nunique(dropna=True)) if ct_col and ct_col in obs else None
    exp_cells = expected.get("cells")
    exp_samples = expected.get("samples")
    exp_cts = expected.get("ct_types")
    row = {
        "metric": ["n cells", "n samples", "n cell types"],
        "observed": [n_cells, n_samples, n_cts],
        "expected (paper)": [exp_cells, exp_samples, exp_cts],
        "delta": [None if exp_cells is None else n_cells - exp_cells,
                  None if exp_samples is None else n_samples - exp_samples,
                  None if exp_cts is None else n_cts - exp_cts],
    }
    return pd.DataFrame(row)


def confounding_crosstab(obs: pd.DataFrame, bio_col: str, batch_cols, max_uniq: int = 60) -> pd.DataFrame:
    """bio x batch contingency tables + collinearity warnings.

    Returns a DataFrame with one row per batch candidate: chi2/p-value of
    independence plus a simple collinearity rule (a bio group or batch level
    occupying a single cell while another level is absent elsewhere). A
    (near-)perfectly collinear design makes bio and batch statistically
    indistinguishable -- loud warning, metrics/UMAP interpretation must note it.
    """
    import scipy.stats as stats

    rows = []
    for bcol in batch_cols:
        if bcol not in obs.columns:
            rows.append({"batch_candidate": bcol, "note": "column not found"})
            continue
        bio = obs[bio_col].astype(str).replace({"nan": "<NA>"})
        bat = obs[bcol].astype(str).replace({"nan": "<NA>"})
        if bio.nunique() > max_uniq or bat.nunique() > max_uniq:
            rows.append(
                {
                    "batch_candidate": bcol,
                    "n_bio_uniq": int(bio.nunique()),
                    "n_batch_uniq": int(bat.nunique()),
                    "note": f"too many unique levels (>{max_uniq}); crosstab skipped",
                }
            )
            continue
        ct = pd.crosstab(bio, bat)
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            rows.append(
                {
                    "batch_candidate": bcol,
                    "n_bio_uniq": int(bio.nunique()),
                    "n_batch_uniq": int(bat.nunique()),
                    "chi2_pvalue": None,
                    "perfectly_collinear": True,
                    "warning": "batch is (near-)perfectly collinear with the bio label -- bio and batch are statistically indistinguishable",
                    "crosstab": ct.to_dict(),
                }
            )
            continue
        chi2, p, dof, _ = stats.chi2_contingency(ct)
        # rule: at least one batch level is present in exactly one bio group
        # (or vice versa) -> strongly confounded
        collinear = bool(((ct > 0).sum(axis=0) == 1).any() or ((ct > 0).sum(axis=1) == 1).any())
        rows.append(
            {
                "batch_candidate": bcol,
                "n_bio_uniq": int(bio.nunique()),
                "n_batch_uniq": int(bat.nunique()),
                "chi2_pvalue": float(p),
                "perfectly_collinear": False,
                "strongly_confounded": collinear,
                "warning": (
                    "bio x batch strongly confounded (a bio group or batch level exists in only one level of the other)"
                    if collinear
                    else ""
                ),
                "crosstab": ct.to_dict(),
            }
        )
    return pd.DataFrame(rows)


def detect_unintegrated_obsm(adata) -> dict:
    """Find precomputed non-harmony PCA and UMAP keys in obsm.

    Excludes any key containing 'harmony' (integrated space). For the paper's
    own preprocessed views the keys are namespaced
    (``X_pca_batch_effect_analysis_hvg2000``) -- still unintegrated and valid.
    Preference for metrics: exact ``X_pca``, then a key containing ``hvg2000``,
    then the lexicographically-smallest remaining candidate.
    """
    keys = list(adata.obsm.keys())
    pca_keys = [k for k in keys if (k == "X_pca" or k.startswith("X_pca_")) and "harmony" not in k]
    umap_keys = [k for k in keys if (k == "X_umap" or k.startswith("X_umap_")) and "harmony" not in k]

    def pick(cands, preferred_substr=None):
        if not cands:
            return None
        for c in cands:
            if c == "X_pca" or c == "X_umap":
                return c
        if preferred_substr:
            for c in cands:
                if preferred_substr in c:
                    return c
        return sorted(cands)[0]

    return {
        "pca_key": pick(pca_keys, "hvg2000"),
        "umap_key": pick(umap_keys),
        "pca_candidates": pca_keys,
        "umap_candidates": umap_keys,
        "all_keys": keys,
    }


def get_counts_matrix_in_memory(adata, subset_idx, verbose: bool = False):
    """Return a dense/sparse counts matrix for a cell subset, from raw counts."""
    if subset_idx is not None:
        sub = adata[subset_idx].to_memory()
    else:
        sub = adata.to_memory()
    counts, slot = locate_counts(sub, verbose=verbose)
    if slot.startswith("raw."):
        # align raw var index to var_names if needed
        counts = sub.raw[:, sub.var_names].X
    return counts, slot


def plot_umap_panels(
    umap_coords: np.ndarray,
    obs_subset: pd.DataFrame,
    label_cols,
    out_dir: Path,
    name: str,
    dpi: int = 150,
    max_cells: int = 200_000,
    seed: int = 0,
    point_size: float | None = None,
) -> list:
    """Scatter panels of the UMAP embedding colored by each label column.

    ``umap_coords`` and ``obs_subset`` must be cell-aligned and already
    subsampled by the caller. Saves one PNG per label column.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n_cells = umap_coords.shape[0]
    if point_size is None:
        point_size = 4.0 if n_cells <= 50_000 else 1.0

    saved = []
    for col in label_cols:
        if col not in obs_subset.columns:
            warnings.warn(f"plot_umap_panels: label column {col!r} not found; skipped")
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        cats = obs_subset[col].astype(str).fillna("NA")
        for cat in pd.unique(cats):
            m = cats.values == cat
            ax.scatter(
                umap_coords[m, 0],
                umap_coords[m, 1],
                s=point_size,
                alpha=0.6,
                linewidths=0,
                label=cat if len(pd.unique(cats)) <= 40 else None,
            )
        if len(pd.unique(cats)) <= 40:
            ax.legend(
                bbox_to_anchor=(1.02, 1), loc="upper left",
                frameon=False, fontsize=8, markerscale=2,
            )
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in col)
        fig.tight_layout()
        path = out_dir / f"{name}_umap_{safe}.png"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        saved.append(str(path))
        print(f"saved {path}")
    return saved


def compute_embed_umap(
    adata,
    compute_subset_max: int = 100_000,
    pca_n_comps: int = 50,
    hvg_n_top: int = 2000,
    seed: int = 0,
    verbose: bool = True,
):
    """Compute a fresh PCA+UMAP on a subsample of RAW counts.

    Returns ``(pca_key, umap_key, n_cells_used, subset_idx)`` where the keys
    point into the returned subset AnnData (``subset.obsm[pca_key]`` etc.).
    The subset object is returned as third element.
    """
    rng = np.random.default_rng(seed)
    n_cells = adata.n_obs
    compute_max = min(n_cells, compute_subset_max)
    subset_idx = rng.choice(n_cells, size=compute_max, replace=False)
    subset_idx = np.sort(subset_idx)
    sub = adata[subset_idx].to_memory()

    counts, slot = locate_counts(sub, verbose=False)
    if slot.startswith("raw."):
        counts = sub.raw[:, sub.var_names].X
    sub.X = counts
    sub.layers["counts"] = counts

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(sub, n_top_genes=hvg_n_top, flavor="seurat", subset=False)
    hvg_mask = sub.var["highly_variable"].values
    if not hvg_mask.any():
        hvg_mask = np.ones(len(sub.var), dtype=bool)
    sub_hvg = sub[:, hvg_mask]
    sc.pp.pca(sub_hvg, n_comps=pca_n_comps, svd_solver="arpack")
    sub.obsm["X_pca_onboard"] = sub_hvg.obsm["X_pca"]
    sub.obsm["X_pca_onboard_hvg"] = sub_hvg.obsm["X_pca"]
    sc.pp.neighbors(sub, use_rep="X_pca_onboard", n_neighbors=15, n_pcs=pca_n_comps)
    sc.tl.umap(sub, min_dist=0.3)
    if verbose:
        print(f"computed PCA+UMAP on {compute_max} cells (raw counts, {slot})")
    return "X_pca_onboard", "X_umap", sub


def embed_and_umap_workflow(
    adata,
    label_cols,
    out_dir: Path,
    name: str,
    sample_col: str | None = None,
    plot_max: int = 200_000,
    compute_subset_max: int = 100_000,
    seed: int = 0,
) -> dict:
    """Driver: plot UMAP panels (precomputed if present, else computed fresh).

    Returns a dict with the resolved pca/umap keys (potentially inside a
    computed subset object), used by the metrics step.
    """
    obsm_info = detect_unintegrated_obsm(adata)
    umap_key = obsm_info["umap_key"]
    pca_key = obsm_info["pca_key"]

    if umap_key is not None:
        rng = np.random.default_rng(seed)
        plot_max = min(adata.n_obs, plot_max)
        idx = np.sort(rng.choice(adata.n_obs, size=plot_max, replace=False))
        coords = np.asarray(adata.obsm[umap_key][idx])
        obs_subset = adata.obs.iloc[idx]
        saved = plot_umap_panels(
            coords, obs_subset, label_cols, out_dir, name, seed=seed
        )
        print(f"used precomputed UMAP obsm key {umap_key!r} ({coords.shape[0]} cells)")
        return {
            "pca_key": pca_key,
            "umap_key": umap_key,
            "computed": False,
            "subset": None,
            "saved_pngs": saved,
            "precomputed_pca_keys": obsm_info["pca_candidates"],
        }

    # compute fresh on raw-count subsample
    pca_key, umap_key, sub = compute_embed_umap(
        adata, compute_subset_max=compute_subset_max, seed=seed
    )
    obs_subset = sub.obs.reset_index(drop=True)
    coords = np.asarray(sub.obsm[umap_key])
    saved = plot_umap_panels(coords, obs_subset, label_cols, out_dir, name, seed=seed)
    return {
        "pca_key": pca_key,
        "umap_key": umap_key,
        "computed": True,
        "subset": sub,
        "saved_pngs": saved,
        "precomputed_pca_keys": obsm_info["pca_candidates"],
    }


def write_metrics_input(
    adata,
    out_path,
    ct_col,
    bio_col,
    batch_cols,
    pca_key: str | None = None,
    sample_col: str | None = None,
    max_cells: int = 300_000,
    seed: int = 0,
    n_pcs: int = 50,
    verbose: bool = True,
) -> dict:
    """Write the feather consumed by ``onboarding_metrics.R``.

    Global random subsample (fixed seed), unintegrated PCA embedding (top
    ``n_pcs`` PCs) + CT/bio/batch labels per cell. If no usable precomputed PCA
    exists, computes PCA on a raw-count subsample. Never uses Harmony keys.
    """
    if pca_key is None:
        pca_key = detect_unintegrated_obsm(adata)["pca_key"]

    if pca_key is not None and pca_key in adata.obsm:
        rng = np.random.default_rng(seed)
        max_cells = min(adata.n_obs, max_cells)
        idx = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        emb = np.asarray(adata.obsm[pca_key][idx])[:, :n_pcs]
        obs_used = adata.obs.iloc[idx]
        pc_source = f"obsm[{pca_key!r}]"
        sub = None
    else:
        pca_key, _umap_key, sub = compute_embed_umap(adata, seed=seed)
        emb = np.asarray(sub.obsm[pca_key])[:, :n_pcs]
        obs_used = sub.obs.reset_index(drop=True)
        pc_source = "computed X_pca (raw-count subsample)"

    df = pd.DataFrame(emb, columns=[f"PC_{i + 1}" for i in range(emb.shape[1])])
    df["cell_index"] = np.arange(len(df))
    for col in [ct_col, bio_col] + list(batch_cols or []):
        if col in obs_used.columns:
            df[col] = obs_used[col].astype(str).fillna("<NA>").values
        else:
            df[col] = "<NA>"
    if sample_col and sample_col in obs_used.columns:
        df[sample_col] = obs_used[sample_col].astype(str).values

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_feather(out_path)
    if verbose:
        print(
            f"wrote metrics input to {out_path}: {df.shape[0]} cells x "
            f"{df.shape[1]} cols (PCA from {pc_source})"
        )
    return {"path": str(out_path), "n_cells": int(len(df)), "pc_source": pc_source}


def save_png(fig, path, dpi: int = 150):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt_close(fig)
    print(f"saved {path}")


def plt_close(fig):
    import matplotlib.pyplot as plt

    plt.close(fig)