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


# ---------------------------------------------------------------------------
# T3.1 -- sample-first, RAM-bounded subsetting (onboarding diagnostics only)
# ---------------------------------------------------------------------------

SUBSET_CONFIG = {
    # Number of samples to select (range 10-20; files with fewer keep all).
    "MAX_SAMPLES": 15,
    # Samples per biological group (range 3-5), before the global cap applies.
    "N_PER_BIO": 4,
    # Random per-sample cap (range 500-5000); bounds every per-sample read so
    # peak memory is independent of the file size.
    "MAX_CELLS_PER_SAMPLE": 2000,
    # Optional per-cell-type cap applied AFTER the concat, stratified by the
    # batch candidates (e.g. 100) to avoid artificial within-CT imbalances.
    "MAX_CELLS_PER_CT": None,
    # Overall soft target: the subset is trimmed (stratified by sample) when it
    # exceeds this many cells.
    "CELLS_TARGET": 10_000,
    # Absolute ceiling the subset never exceeds (used when CELLS_TARGET is 0).
    "CELLS_MAX": 50_000,
    "SEED": 0,
}


def _strata_codes(obs: pd.DataFrame, strat_cols) -> np.ndarray:
    """Combine several categorical/string columns into one dense int array."""
    if not strat_cols:
        return np.zeros(len(obs), dtype=int)
    out = np.zeros(len(obs), dtype=np.int64)
    basis = 1
    for col in strat_cols:
        if col not in obs.columns:
            continue
        codes, _ = pd.factorize(obs[col].astype(str))
        codes = np.where(codes < 0, 0, codes)
        out += codes * basis
        basis *= int(codes.max()) + 1
    return out


def _cap_group_values(obs: pd.DataFrame, group_col, cap, strat_cols, seed) -> np.ndarray:
    """Downsample obs rows so each value of ``group_col`` holds at most ``cap``
    rows, distributing the cap across the ``strat_cols`` strata proportionally
    to their within-group counts (fractional remainder filled round-robin by
    largest fraction). Returns a boolean keep-mask over ``obs``.
    """
    n = len(obs)
    if cap is None or cap <= 0:
        return np.ones(n, dtype=bool)
    rng = np.random.default_rng(seed)
    gcodes, guniq = pd.factorize(obs[group_col].astype(str))
    strat = _strata_codes(obs, strat_cols) if strat_cols else np.zeros(n, dtype=int)
    keep = np.zeros(n, dtype=bool)
    for g in range(len(guniq)):
        rows = np.flatnonzero(gcodes == g)
        total = len(rows)
        if total <= cap:
            keep[rows] = True
            continue
        strata, counts = np.unique(strat[rows], return_counts=True)
        floors = np.floor(cap * counts / total).astype(int)
        leftover = int(cap - floors.sum())
        kept = np.zeros(total, dtype=bool)
        for i in range(len(strata)):
            sel = np.flatnonzero(strat[rows] == strata[i])
            take = int(floors[i])
            if take >= len(sel):
                kept[sel] = True
            else:
                p = rng.choice(len(sel), size=take, replace=False)
                kept[sel[np.sort(p)]] = True
        if leftover > 0:
            taken = floors.astype(int).copy()
            frac = cap * counts / total - floors
            order = np.argsort(-frac, kind="stable")
            for i in order:
                if leftover <= 0:
                    break
                sel = np.flatnonzero(strat[rows] == strata[i])
                avail = sel[~kept[sel]]
                if len(avail) == 0 or taken[i] >= len(sel):
                    continue
                kept[avail[rng.integers(0, len(avail))]] = True
                taken[i] += 1
                leftover -= 1
        keep[rows] = kept
    return keep


def _cap_total_stratified(obs: pd.DataFrame, strat_col, cap, seed) -> np.ndarray:
    """Cap the TOTAL number of rows to ``cap``, distributed across the
    ``strat_col`` values (strata) proportionally to their counts (fractional
    remainder filled round-robin by largest fraction) so every stratum stays
    represented. Without a valid strat column, a plain seeded random sample.
    Returns a boolean keep-mask over ``obs``.
    """
    n = len(obs)
    if cap is None or cap <= 0 or n <= cap:
        return np.ones(n, dtype=bool)
    rng = np.random.default_rng(seed)
    if strat_col is None or strat_col not in obs.columns:
        idx = rng.choice(n, size=int(cap), replace=False)
        keep = np.zeros(n, dtype=bool)
        keep[np.sort(idx)] = True
        return keep
    codes, _ = pd.factorize(obs[strat_col].astype(str))
    counts = np.bincount(codes, minlength=int(codes.max()) + 1)
    floors = np.floor(cap * counts / n).astype(int)
    leftover = int(cap - floors.sum())
    keep = np.zeros(n, dtype=bool)
    for s in range(len(counts)):
        sel = np.flatnonzero(codes == s)
        take = int(floors[s])
        if take >= len(sel):
            keep[sel] = True
        elif take > 0:
            p = rng.choice(len(sel), size=take, replace=False)
            keep[sel[np.sort(p)]] = True
    if leftover > 0:
        taken = floors.astype(int).copy()
        frac = cap * counts / n - floors
        order = np.argsort(-frac, kind="stable")
        for s in order:
            if leftover <= 0:
                break
            sel = np.flatnonzero(codes == s)
            avail = sel[~keep[sel]]
            if len(avail) == 0 or taken[s] >= len(sel):
                continue
            keep[avail[rng.integers(0, len(avail))]] = True
            taken[s] += 1
            leftover -= 1
    return keep


def subset_by_samples(
    adata,
    sample_col: str | None = None,
    bio_col: str | None = None,
    batch_cols: list | None = None,
    ct_col: str | None = None,
    config: dict | None = None,
    verbose: bool = True,
) -> tuple:
    """Sample-first, RAM-bounded cell subsetting for onboarding diagnostics.

    Flow: ``obs``-only read -> select ``MAX_SAMPLES`` samples stratified by the
    bio condition and (round-robin over the batch-candidate values) the batch
    candidates -> read the selected samples' cells PER SAMPLE via
    ``.to_memory()`` (peak memory bounded by ``MAX_CELLS_PER_SAMPLE`` regardless
    of file size) -> optional per-CT cap (stratified by batch) -> overall
    target cap (stratified by sample, every selected sample stays represented).
    Small files (< 10 samples) keep all samples.

    The returned AnnData is in memory (<= ``CELLS_MAX`` cells) and carries the
    parent's precomputed UNINTEGRATED ``obsm`` arrays (``X_pca*``/``X_umap*``,
    never harmony) re-sliced to the subset rows, so the downstream UMAP and
    LISI metrics stay in the original global embedding space. DIAGNOSTIC ONLY
    -- the HPC pipeline always uses the full data.

    Returns ``(subset, summary_dict)`` with samples-per-bio-group,
    cells-per-sample, selected samples and the config actually used.
    """
    cfg = dict(SUBSET_CONFIG)
    if config:
        cfg.update(config)
    if cfg.get("CELLS_TARGET") and cfg.get("CELLS_MAX"):
        if cfg["CELLS_TARGET"] > cfg["CELLS_MAX"]:
            raise ValueError("SUBSET_CONFIG: CELLS_TARGET must be <= CELLS_MAX")
    seed = int(cfg["SEED"])
    rng = np.random.default_rng(seed)

    if sample_col is None or sample_col not in adata.obs.columns:
        cand = candidate_col_detection(adata.obs)
        raise ValueError(
            f"subset_by_samples: sample_col {sample_col!r} not in obs columns; "
            f"candidate sample columns: {cand['sample']}"
        )

    obs = adata.obs

    if bio_col and bio_col in obs.columns and obs[bio_col].nunique(dropna=True) >= 2:
        bio = bio_col
    else:
        if bio_col:
            warnings.warn(
                f"subset_by_samples: bio_col {bio_col!r} missing or single-valued; "
                f"treating all samples as one group"
            )
        bio = None

    present_batch = [c for c in (batch_cols or []) if c in obs.columns]
    if batch_cols and len(present_batch) != len(batch_cols):
        warnings.warn(
            f"subset_by_samples: batch columns not found in obs: "
            f"{sorted(set(batch_cols) - set(present_batch))}"
        )

    # --- per-sample bookkeeping (obs only; cheap even for 2M cells) ----------
    vc = obs[sample_col].value_counts(dropna=True)
    n_all_samples = int(vc.size)
    codes, uniq_samples = pd.factorize(obs[sample_col].astype(str))
    sample_code = {str(u): i for i, u in enumerate(uniq_samples)}
    sample_list = [str(s) for s in vc.index]
    n_na = int((codes < 0).sum())
    if n_na:
        warnings.warn(
            f"subset_by_samples: excluding {n_na} cells with missing {sample_col!r}"
        )

    valid = codes >= 0
    smp = pd.DataFrame({"sample": obs[sample_col].astype(str).values})[valid]
    if bio is not None:
        smp["bio"] = obs[bio].astype(str).values[valid]
    if present_batch:
        for bcol in present_batch:
            smp[bcol] = obs[bcol].astype(str).values[valid]

    bio_of = (
        smp.groupby("sample")["bio"].agg(lambda s: s.mode().iloc[0] if len(s.mode()) else "<NA>").to_dict()
        if bio is not None
        else {}
    )
    batch_tuple_of = {}
    if present_batch:
        modes = smp.groupby("sample")[present_batch].agg(
            lambda s: s.mode().iloc[0] if len(s.mode()) else "<NA>"
        )
        batch_tuple_of = {s: tuple(row) for s, (_, row) in zip(modes.index, modes.iterrows())}

    # --- sample selection: bio stratification + round-robin over batches -----
    if bio is not None:
        groups = sorted({bio_of.get(s, "<NA>") for s in sample_list})
        n_group = {g: sum(1 for s in sample_list if bio_of.get(s, "<NA>") == g) for g in groups}
        alloc = {g: min(1, n_group[g]) for g in groups}
        remaining = int(cfg["MAX_SAMPLES"]) - sum(alloc.values())
        order_groups = sorted(groups, key=lambda g: -n_group[g])
        while remaining > 0:
            advanced = False
            for g in order_groups:
                if alloc[g] >= min(int(cfg["N_PER_BIO"]), n_group[g]):
                    continue
                alloc[g] += 1
                remaining -= 1
                advanced = True
                if remaining == 0:
                    break
            if not advanced:
                break
    else:
        groups = [None]
        alloc = {None: min(int(cfg["MAX_SAMPLES"]), len(sample_list))}
    if sum(alloc.values()) < 1:
        raise ValueError("subset_by_samples: no samples available to select")

    def ordered_samples(grp_samples):
        """Round-robin over the batch-value tuples so distinct batch values are
        covered early; deterministic shuffle fallback when no batch columns."""
        if not batch_tuple_of:
            out = list(grp_samples)
            rng.shuffle(out)
            return out
        buckets: dict = {}
        for s in grp_samples:
            buckets.setdefault(batch_tuple_of.get(s), []).append(s)
        for lst in buckets.values():
            rng.shuffle(lst)
        keys = list(buckets)
        out = []
        i = 0
        guard = 0
        while any(buckets.values()):
            k = keys[i % len(keys)]
            if buckets[k]:
                out.append(buckets[k].pop())
            i += 1
            guard += 1
            if guard > (len(keys) + len(grp_samples)) * 2:
                break
        return out

    selected = []
    for g in groups:
        if bio is None:
            grp_samples = list(sample_list)
        else:
            grp_samples = [s for s in sample_list if bio_of.get(s, "<NA>") == g]
        selected.extend(ordered_samples(grp_samples)[: alloc[g]])

    # --- per-sample reads (bounded by MAX_CELLS_PER_SAMPLE) ------------------
    cap_sample = int(cfg["MAX_CELLS_PER_SAMPLE"])
    parts = []
    parent_positions = []
    for s in selected:
        c = sample_code.get(s)
        if c is None:
            continue
        idx = np.flatnonzero(codes == c)
        if len(idx) > cap_sample:
            pick = rng.choice(len(idx), size=cap_sample, replace=False)
            idx = np.sort(idx[pick])
        else:
            idx = np.sort(idx)
        parent_positions.append(idx)
        part = adata[idx].to_memory()
        part.obsm = {}
        part.obsp = {}
        parts.append(part)
    if not parts:
        raise ValueError("subset_by_samples: no cells selected")
    parent_pos = np.concatenate(parent_positions)

    sub = ad.concat(parts, join="outer", merge="same")
    del parts

    # Re-slice the parent's unintegrated obsm arrays to the subset rows so the
    # downstream UMAP/metrics stay in the original (global) embedding space.
    # Backed h5ad obsm arrays live in the file as h5py Datasets but
    # `adata.obsm[k]` materializes the full array on access -- read the
    # dataset directly and fancy-index it in the file (memory bounded by the
    # subset size, independent of file size). Fall back to the plain access
    # for in-memory parents / non-dataset (e.g. sparse) entries.
    try:
        import h5py
    except ImportError:  # pragma: no cover
        h5py = None
    h5_obsm = None
    if getattr(adata, "isbacked", False):
        try:
            h5_obsm = adata.file["obsm"]
        except (KeyError, TypeError):
            h5_obsm = None
    obsm_info = detect_unintegrated_obsm(adata)
    for k in obsm_info["pca_candidates"] + obsm_info["umap_candidates"]:
        if h5_obsm is not None and h5py is not None:
            ds = h5_obsm.get(k)
            if isinstance(ds, h5py.Dataset) and ds.ndim == 2 and ds.shape[0] == adata.n_obs:
                sub.obsm[k] = ds[parent_pos]
                continue
        arr = np.asarray(adata.obsm[k])
        if arr.ndim == 2 and arr.shape[0] == adata.n_obs:
            sub.obsm[k] = arr[parent_pos]

    # --- post-concat caps -----------------------------------------------------
    if cfg.get("MAX_CELLS_PER_CT") and ct_col and ct_col in sub.obs.columns:
        keep = _cap_group_values(
            sub.obs, ct_col, int(cfg["MAX_CELLS_PER_CT"]), present_batch, seed
        )
        sub = sub[keep]
    cap_max = cfg.get("CELLS_MAX")
    cap_target = cfg.get("CELLS_TARGET") or cap_max
    if cap_max and sub.n_obs > cap_max:
        keep = _cap_total_stratified(sub.obs, sample_col, int(cap_max), seed)
        sub = sub[keep]
    if cap_target and sub.n_obs > cap_target:
        keep = _cap_total_stratified(sub.obs, sample_col, int(cap_target), seed)
        sub = sub[keep]

    # --- summary --------------------------------------------------------------
    # Keep the keys str-normalized like every earlier stage (sample_list,
    # bio_of, sample_code) so samples_per_bio/cells_per_sample stay consistent
    # even for int-typed sample columns.
    cells_per_sample = (
        sub.obs[sample_col].astype(str).value_counts().astype(int).to_dict()
    )
    samples_selected = list(cells_per_sample)
    samples_per_bio = {}
    for g in groups:
        if bio is None:
            samples_per_bio["all"] = len(samples_selected)
        else:
            samples_per_bio[g] = sum(
                1 for s in samples_selected if bio_of.get(s, "<NA>") == g
            )
    summary = {
        "n_samples_total": n_all_samples,
        "n_samples_selected": len(samples_selected),
        "samples_per_bio": samples_per_bio,
        "selected_samples": samples_selected,
        "cells_per_sample": cells_per_sample,
        "n_cells_total": int(sub.n_obs),
        "n_cells_full": int(obs.shape[0]),
        "config": dict(cfg),
    }
    if verbose:
        print(
            f"subset_by_samples: {summary['n_cells_total']} cells from "
            f"{summary['n_samples_selected']}/{summary['n_samples_total']} samples "
            f"(full file: {summary['n_cells_full']} cells) -- diagnostic subset "
            f"only; the HPC pipeline uses the full data."
        )
    return sub, summary


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


def load_onboarding_dataset_or_subset(name: str, nas_file: str | Path, here: Path) -> dict:
    """Find and load either a precomputed diagnostic subset or the full NAS file.

    Search priority for precomputed subset:
      1. data/new_dataset_checks/subsets/<name>_subset.h5ad
      2. data/new_dataset_checks/<name>/<name>_subset.h5ad
      3. <nas_dir>/subsets/<name>_subset.h5ad

    Returns a dict:
      mode: 'subset' or 'nas'
      sub: in-memory subset AnnData (if mode=='subset')
      adata: backed full AnnData (if mode=='nas')
      meta: dict from <name>_meta.json (if mode=='subset')
    """
    root = here.parent.parent
    nas_path = Path(nas_file)
    cands = [
        root / "data" / "new_dataset_checks" / "subsets" / f"{name}_subset.h5ad",
        root / "data" / "new_dataset_checks" / name / f"{name}_subset.h5ad",
        nas_path.parent.parent / "subsets" / f"{name}_subset.h5ad",
    ]

    for subset_p in cands:
        meta_p = subset_p.with_name(f"{name}_meta.json")
        if subset_p.exists() and meta_p.exists():
            print(f"[{name}] Loading precomputed subset from: {subset_p}")
            sub = sc.read_h5ad(subset_p)
            with open(meta_p) as f:
                meta = json.load(f)
            print(f"[{name}] Precomputed metadata loaded (created {meta.get('created_at', 'unknown')})")
            return {
                "mode": "subset",
                "sub": sub,
                "adata": None,
                "meta": meta,
                "subset_path": str(subset_p),
                "meta_path": str(meta_p),
            }

    if nas_path.exists():
        print(f"[{name}] Opening full file from NAS in backed mode: {nas_path}")
        adata = sc.read_h5ad(nas_path, backed="r")
        return {
            "mode": "nas",
            "sub": None,
            "adata": adata,
            "meta": None,
            "nas_path": str(nas_path),
        }

    err_msg = (
        f"[{name}] ERROR: Neither precomputed subset nor NAS file found.\n"
        f"  Looked for subset at: {[str(p) for p in cands]}\n"
        f"  Looked for NAS file at: {nas_path}\n\n"
        f"To generate subsets on the HPC and pull them locally:\n"
        f"  1. Run on HPC: ./notebooks/dataset_onboarding/run_subset_hpc.sh --only {name.lower().split('_')[0]}\n"
        f"  2. Pull to Mac: rsync -avP bamboo:scratch/ECODA_paper/_downloads/subsets/ data/new_dataset_checks/subsets/"
    )
    raise FileNotFoundError(err_msg)