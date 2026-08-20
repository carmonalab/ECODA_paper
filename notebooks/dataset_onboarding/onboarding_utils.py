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
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure full visibility without clipping in notebook tables
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", 1000)

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

    Slot chain: ``layers["counts"]`` -> ``raw.X`` (if raw has integer counts while X is normalized) -> ``X`` -> ``raw.X``.
    """
    for layer_name in ("counts", "counts_cite", "counts_atac"):
        if layer_name in adata.layers:
            if verbose:
                print(f"counts slot: adata.layers[{layer_name!r}]")
            return adata.layers[layer_name], f"layers/{layer_name}"

    # If raw.X exists, check if X is scaled or log-normalized while raw.X holds integer counts
    if adata.raw is not None and adata.raw.X is not None:
        sample_x = _sample_values(adata.X, n_sample_cells=1000) if adata.X is not None else np.array([])
        sample_raw = _sample_values(adata.raw.X, n_sample_cells=1000)
        x_is_int = bool(sample_x.size > 0 and np.all(np.abs(sample_x - np.round(sample_x)) < 1e-4))
        raw_is_int = bool(sample_raw.size > 0 and np.all(np.abs(sample_raw - np.round(sample_raw)) < 1e-4))
        if len(sample_x) > 0 and (np.any(sample_x < 0) or (not x_is_int and raw_is_int)):
            if verbose:
                print("counts slot: adata.raw.X (adata.raw.X contains raw integer counts while adata.X is normalized/scaled)")
            return adata.raw.X, "raw.X"

    if adata.X is not None:
        if verbose:
            print("counts slot: adata.X")
        return adata.X, "X"
    if adata.raw is not None and adata.raw.X is not None:
        if verbose:
            print("counts slot: adata.raw.X")
        return adata.raw.X, "raw.X"
    raise ValueError("No counts found: no 'counts' layer, no X, no raw.X.")


def _extract_sparse_data_array(matrix):
    """Extract the 1-D array / dataset of non-zero values from any sparse matrix."""
    if hasattr(matrix, "format") and hasattr(matrix, "data"):
        # In scipy.sparse (csr_matrix, csc_matrix, etc.): matrix.data is a 1D ndarray
        return matrix.data
    if hasattr(matrix, "_data"):
        # In anndata backed sparse dataset (_CSRDataset, _CSCDataset): matrix._data is an h5py.Dataset
        return matrix._data
    if hasattr(matrix, "group") and hasattr(matrix.group, "__getitem__") and "data" in matrix.group:
        return matrix.group["data"]
    return None


def _is_sparse(matrix) -> bool:
    if hasattr(matrix, "format") or hasattr(matrix, "tocsr") or hasattr(matrix, "toarray") or hasattr(matrix, "_data"):
        return True
    if hasattr(matrix, "group") and hasattr(matrix.group, "__getitem__") and "data" in matrix.group:
        return True
    return False


def _sample_values(matrix, n_sample_cells: int = 200_000, seed: int = 0):
    """Draw a flat sample of matrix values (numpy 1-D) for sanity checks."""
    rng = np.random.default_rng(seed)
    n_cells, n_genes = matrix.shape
    total = n_cells * n_genes
    n_pick = min(n_sample_cells, total)
    if total == 0:
        return np.array([], dtype=float)

    sp_data = _extract_sparse_data_array(matrix)
    if sp_data is not None:
        n_data = len(sp_data)
        if n_data == 0:
            return np.array([], dtype=float)
        take = min(n_pick, n_data)
        if isinstance(sp_data, np.ndarray):
            idx = rng.choice(n_data, size=take, replace=False)
            return np.asarray(sp_data[np.sort(idx)], dtype=float)
        # H5PY / on-disk dataset: sample random contiguous blocks to avoid
        # 200,000 non-contiguous random seeks on compressed chunks
        n_blocks = 10
        block_size = max(1000, take // n_blocks)
        max_start = max(0, n_data - block_size)
        if max_start == 0:
            return np.asarray(sp_data[:take], dtype=float)
        starts = np.sort(rng.integers(0, max_start, size=n_blocks))
        samples = [sp_data[s : s + block_size] for s in starts]
        res = np.concatenate(samples)[:take]
        return np.asarray(res, dtype=float)

    # In-memory dense numpy array
    if isinstance(matrix, np.ndarray):
        flat = matrix.ravel()
        take = min(n_pick, len(flat))
        idx = np.sort(rng.choice(len(flat), size=take, replace=False))
        return np.asarray(flat[idx], dtype=float)

    # Dense backed h5py dataset (sample rows to avoid pulling full matrix)
    if hasattr(matrix, "__getitem__"):
        n_rows = min(1000, n_cells)
        row_idx = np.sort(rng.choice(n_cells, size=n_rows, replace=False))
        sub_chunk = np.asarray(matrix[row_idx, :]).ravel()
        take = min(n_pick, len(sub_chunk))
        idx = np.sort(rng.choice(len(sub_chunk), size=take, replace=False))
        return np.asarray(sub_chunk[idx], dtype=float)

    return np.array([], dtype=float)


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

    if hasattr(counts, "nnz"):
        n_nonzero = counts.nnz
    elif hasattr(counts, "_data"):
        n_nonzero = len(counts._data)
    elif hasattr(counts, "group") and hasattr(counts.group, "__getitem__") and "data" in counts.group:
        n_nonzero = len(counts.group["data"])
    elif sparse:
        sp_data = _extract_sparse_data_array(counts)
        n_nonzero = len(sp_data) if sp_data is not None else 0
    else:
        n_nonzero = int(np.count_nonzero(counts))

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


def _fast_nunique(s: pd.Series) -> int:
    """Fast nunique check on Series avoiding full scans on categoricals."""
    if isinstance(s.dtype, pd.CategoricalDtype):
        return int(len(s.cat.categories))
    return int(s.nunique(dropna=True))


def obs_summary(adata, top_n: int = 20) -> pd.DataFrame:
    """Per-column dtype + top value counts (bounded) for obs."""
    rows = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        nunique = _fast_nunique(s)
        if isinstance(s.dtype, pd.CategoricalDtype):
            vc = s.value_counts(dropna=True).head(top_n).to_dict()
            rows.append({
                "column": col,
                "dtype": str(s.dtype),
                "n_unique": nunique,
                "na_count": int(s.isna().sum()),
                "top_values": json.dumps(vc, default=str)[:400],
            })
        elif s.dtype.name == "object" or nunique <= 30:
            if len(s) > 200_000 and nunique > 500:
                # High-cardinality string/identifier column (e.g. cell barcodes)
                vc = {"<unique_ids>": nunique}
            else:
                vc = s.value_counts(dropna=True).head(top_n).to_dict()
            rows.append({
                "column": col,
                "dtype": str(s.dtype),
                "n_unique": nunique,
                "na_count": int(s.isna().sum()),
                "top_values": json.dumps(vc, default=str)[:400],
            })
        else:
            rows.append({
                "column": col,
                "dtype": str(s.dtype),
                "n_unique": nunique,
                "na_count": int(s.isna().sum()),
                "top_values": "",
            })
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
        nu = _fast_nunique(obs[c])
        if any(t in lc for t in sample_terms) and nu >= 2:
            out["sample"].append(c)
        if any(t in lc for t in label_terms) and 2 <= nu <= 50:
            out["label"].append(c)
        if any(t in lc for t in batch_terms) and 2 <= nu <= 200:
            out["batch"].append(c)
        if any(t in lc for t in ct_terms) and 2 <= nu <= 200:
            out["cell_type"].append(c)
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
            rows.append({
                "batch_candidate": bcol,
                "n_bio_uniq": 0,
                "n_batch_uniq": 0,
                "chi2_pvalue": None,
                "perfectly_collinear": False,
                "strongly_confounded": False,
                "warning": "column not found",
                "note": "column not found",
                "crosstab": {},
            })
            continue
        bio = obs[bio_col].astype(str).replace({"nan": "<NA>"})
        bat = obs[bcol].astype(str).replace({"nan": "<NA>"})
        if bio.nunique() > max_uniq or bat.nunique() > max_uniq:
            rows.append(
                {
                    "batch_candidate": bcol,
                    "n_bio_uniq": int(bio.nunique()),
                    "n_batch_uniq": int(bat.nunique()),
                    "chi2_pvalue": None,
                    "perfectly_collinear": False,
                    "strongly_confounded": False,
                    "warning": f"too many unique levels (>{max_uniq}); crosstab skipped",
                    "note": f"too many unique levels (>{max_uniq}); crosstab skipped",
                    "crosstab": {},
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
                    "strongly_confounded": True,
                    "warning": "batch is (near-)perfectly collinear with the bio label -- bio and batch are statistically indistinguishable",
                    "note": "",
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
                "note": "",
                "crosstab": ct.to_dict(),
            }
        )
    return pd.DataFrame(rows)


def cell_type_harmonization_check(
    obs: pd.DataFrame,
    ct_cols: list[str],
    batch_col: str | None = None,
    sample_col: str | None = None,
) -> pd.DataFrame:
    """Evaluate cross-batch / cross-study cell type harmonization and label sharing.

    For each candidate cell type annotation column, assesses whether annotations
    are atlas-wide harmonized (present across batches) or study-specific.
    """
    rows = []
    has_batch = batch_col is not None and batch_col in obs.columns
    n_batches = int(obs[batch_col].nunique(dropna=True)) if has_batch else 1

    for col in ct_cols:
        if col not in obs.columns:
            continue
        nu = int(obs[col].nunique(dropna=True))
        if nu == 0:
            continue

        if not has_batch or n_batches <= 1:
            rows.append({
                "Candidate Column": col,
                "N Cell Types": nu,
                "Shared (≥80% Batches)": f"{nu}/{nu} (100.0%)",
                "Single-Batch Only": f"0/{nu} (0.0%)",
                "Sharing Status": "Harmonized (Single Batch / Study)",
                "Recommendation / Notes": "Single cohort context; all labels shared.",
            })
            continue

        ct = pd.crosstab(obs[col], obs[batch_col])
        b_per_ct = (ct > 0).sum(axis=1)
        thr_80 = max(2, int(np.ceil(0.8 * n_batches)))
        n_shared_80 = int((b_per_ct >= thr_80).sum())
        n_single = int((b_per_ct == 1).sum())
        pct_shared = (n_shared_80 / nu * 100) if nu else 0.0
        pct_single = (n_single / nu * 100) if nu else 0.0

        if pct_shared >= 75.0 and pct_single <= 20.0:
            status = "Harmonized (Atlas-Wide)"
            notes = "Atlas-wide harmonized labels; optimal for ECODA patient stratification."
        elif pct_shared >= 35.0 or n_shared_80 >= 4:
            status = "Partially Harmonized"
            notes = "Core lineages shared across studies; fine subtypes partially batch-restricted."
        else:
            status = "Study-Specific / Unharmonized"
            notes = "Disjoint label sets per study; unsuited for cross-study ECODA stratification."

        rows.append({
            "Candidate Column": col,
            "N Cell Types": nu,
            "Shared (≥80% Batches)": f"{n_shared_80}/{nu} ({pct_shared:.1f}%)",
            "Single-Batch Only": f"{n_single}/{nu} ({pct_single:.1f}%)",
            "Sharing Status": status,
            "Recommendation / Notes": notes,
        })
    return pd.DataFrame(rows)


def categorize_obs_columns(obs: pd.DataFrame, config: dict | None = None) -> pd.DataFrame:
    """Classify all obs columns into 7 canonical roles with metadata summaries.

    Roles:
      1. Main Biological Condition (PILOT-GM-VAE / Primary Contrast)
      2. Secondary / Demographic Biological Covariates
      3. Batch Effect Candidates (Technical Covariates)
      4. Cell Type Annotations
      5. Sample & Donor Identifiers
      6. Technical QC Metrics & Single-Cell Artifacts
      7. Uninformative / Constant Columns
    """
    cfg = config or {}
    role_map = {}

    for col in cfg.get("bio_primary", []):
        if col in obs.columns:
            role_map[col] = "1. Main Biological Condition"
    for col in cfg.get("bio_secondary", []):
        if col in obs.columns:
            role_map[col] = "2. Secondary / Demographic Biological Covariates"
    for col in cfg.get("batch_candidates", []):
        if col in obs.columns:
            role_map[col] = "3. Batch Effect Candidates (Technical Covariates)"
    for col in cfg.get("cell_types", []):
        if col in obs.columns:
            role_map[col] = "4. Cell Type Annotations"
    for col in cfg.get("sample_ids", []):
        if col in obs.columns:
            role_map[col] = "5. Sample & Donor Identifiers"
    for col in cfg.get("qc_metrics", []):
        if col in obs.columns:
            role_map[col] = "6. Technical QC Metrics & Single-Cell Artifacts"
    for col in cfg.get("uninformative", []):
        if col in obs.columns:
            role_map[col] = "7. Uninformative / Constant Columns"

    uninformative_exact = {
        "is_primary_data", "organism", "observation_joinid", "mapped_reference_assembly",
        "mapped_reference_annotation", "alignment_software", "organism_ontology_term_id",
        "suspension_uuid", "library_uuid", "sample_uuid",
    }
    qc_terms = [
        "ncount", "nfeature", "n_genes", "n_counts", "total_counts", "percent_mito", "percent_rb",
        "percent.mt", "pct_counts_mito", "doublet_score", "doublet", "emptydrops", "matrisome",
        "collagen", "s.score", "g2m.score", "rin", "low_q", "degen.score", "percent.er", "glycoprotein",
        "proteoglycan", "cyc.score", "aepi.score", "astr.score", "genes detected", "number of umis",
        "fraction mitochrondrial", "log10_n_counts", "mt_frac", "nfeaturess_rna", "dissociation_s1",
    ]
    sample_terms = [
        "donor_id", "sample_id", "sampleid", "sample", "specimen", "patientid", "patient",
        "file_id", "orig_ident", "ind_cov", "accsample", "patient_region_id",
    ]
    ct_terms = [
        "cell_type", "celltype", "cell_types", "majortype", "broad_cell_type", "author_cell_type",
        "subclass", "class", "ann_fine", "ann_coarse", "cell_subtype", "cg_cov", "opt_clust", "supertype",
        "cell_state", "louvain", "leiden", "accannot", "acclabel", "structure", "state", "major_cell_types",
    ]
    batch_terms = [
        "assay", "batch", "batch_cov", "sequencing_platform", "platform", "single cell sequencing platform",
        "city", "brain_bank", "hospital", "study", "dataset", "tissue_type", "tissue_location",
        "region", "dissociation", "preservation", "sample_source", "library", "sample_prep", "site",
        "experiment", "processing_cohort", "sampletype", "tissue_sampling_method", "donor_status",
    ]
    bio_terms = [
        "disease", "condition", "status", "severity", "patient_group", "origin", "cognitive status",
        "age", "sex", "bmi", "braak", "adnc", "cerad", "outcome", "comorbidities", "diabetes",
        "hypertension", "egfr", "smoker", "tumor_stage", "mutation", "pop_cov", "sle_status", "pathology",
        "uicc_stage", "ever_smoker", "thal phase", "apoe4", "late-nc", "lewy body", "microinfarct",
        "neurotypical", "breast_density", "sars-cov-2", "sampling day", "strain", "diabetes_model",
    ]

    rows = []
    for col in obs.columns:
        s = obs[col]
        lc = col.lower()
        nu = int(s.nunique(dropna=True))
        na_pct = float(s.isna().mean() * 100)
        dtype_str = str(s.dtype)

        if nu == 0:
            val_summary = "<Empty / All NA>"
        elif pd.api.types.is_numeric_dtype(s) and nu > 10:
            val_summary = f"[{s.min():.2g}, {s.max():.2g}] (mean: {s.mean():.2g})"
        else:
            top_vals = s.dropna().astype(str).unique()
            if len(top_vals) <= 4:
                val_summary = ", ".join(top_vals)
            else:
                val_summary = ", ".join(top_vals[:3]) + f", ... ({nu} total)"

        if col in role_map:
            role = role_map[col]
        elif nu <= 1 or col in uninformative_exact or "_ontology_term_id" in lc:
            role = "7. Uninformative / Constant Columns"
        elif any(t in lc for t in qc_terms):
            role = "6. Technical QC Metrics & Single-Cell Artifacts"
        elif any(t in lc for t in sample_terms) or nu > len(obs) * 0.8:
            role = "5. Sample & Donor Identifiers"
        elif any(t in lc for t in ct_terms):
            role = "4. Cell Type Annotations"
        elif any(t in lc for t in batch_terms):
            role = "3. Batch Effect Candidates (Technical Covariates)"
        elif any(t in lc for t in bio_terms):
            role = "2. Secondary / Demographic Biological Covariates"
        else:
            role = "2. Secondary / Demographic Biological Covariates"

        rows.append({
            "Column": col,
            "Assigned Role": role,
            "Data Type": dtype_str,
            "Unique Levels": nu,
            "Missing %": f"{na_pct:.1f}%",
            "Sample Values / Categories": val_summary,
        })

    df = pd.DataFrame(rows)
    return df.sort_values(by=["Assigned Role", "Column"]).reset_index(drop=True)


def summarize_obs_categories(cat_df: pd.DataFrame) -> pd.DataFrame:
    """Concise executive overview of metadata categories and their assigned columns."""
    rows = []
    for role, group in cat_df.groupby("Assigned Role", sort=True):
        cols = group["Column"].tolist()
        rows.append({
            "Assigned Role": role,
            "Column Count": len(cols),
            "Assigned Columns": ", ".join(cols),
        })
    return pd.DataFrame(rows)


def compute_clr_composition(
    obs: pd.DataFrame,
    sample_col: str,
    ct_col: str,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """Compute centered log-ratio (CLR) abundance per sample and cell type.

    Counts are imputed with ``+pseudocount`` (matching ECODA counts_all default),
    normalized to proportions, and CLR-transformed.
    """
    counts = pd.crosstab(obs[sample_col].astype(str), obs[ct_col].astype(str)).astype(float)
    counts_plus = counts + float(pseudocount)
    props = counts_plus.div(counts_plus.sum(axis=1), axis=0)
    log_props = np.log(props)
    clr = log_props.sub(log_props.mean(axis=1), axis=0)
    return clr


def plot_clr_abundance_boxplots(
    clr_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_col: str,
    group_col: str,
    max_cts: int = 20,
    title: str | None = None,
):
    """Plot grouped boxplots of sample-level CLR abundances with Wilcoxon/Kruskal-Wallis tests.

    Annotates Benjamini-Hochberg FDR-corrected significance stars above cell types.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy.stats as stats
    import statsmodels.stats.multitest as smt

    meta_dedup = meta_df[[sample_col, group_col]].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    common_samples = clr_df.index.intersection(meta_dedup.index)
    clr_sub = clr_df.loc[common_samples]
    meta_sub = meta_dedup.loc[common_samples]

    valid_mask = (
        meta_sub[group_col].notna()
        & (meta_sub[group_col].astype(str) != "<NA>")
        & (meta_sub[group_col].astype(str) != "nan")
    )
    clr_sub = clr_sub[valid_mask]
    meta_sub = meta_sub[valid_mask]

    if len(meta_sub) == 0 or clr_sub.shape[1] == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, f"No valid samples for {group_col}", ha="center", va="center")
        return fig

    s_grp = meta_sub[group_col].astype(str)
    uniq_levels = np.sort(s_grp.unique())

    if clr_sub.shape[1] > max_cts:
        variances = clr_sub.var(axis=0)
        selected_cts = variances.nlargest(max_cts).index.tolist()
    else:
        selected_cts = clr_sub.columns.tolist()

    raw_pvals = []
    for ct in selected_cts:
        vals = clr_sub[ct]
        if len(uniq_levels) < 2:
            raw_pvals.append(np.nan)
        elif len(uniq_levels) == 2:
            g1 = vals[s_grp == uniq_levels[0]]
            g2 = vals[s_grp == uniq_levels[1]]
            if len(g1) > 0 and len(g2) > 0:
                try:
                    _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                    raw_pvals.append(p)
                except Exception:
                    raw_pvals.append(np.nan)
            else:
                raw_pvals.append(np.nan)
        else:
            grouped = [vals[s_grp == lvl] for lvl in uniq_levels if (s_grp == lvl).sum() > 0]
            if len(grouped) >= 2:
                try:
                    _, p = stats.kruskal(*grouped)
                    raw_pvals.append(p)
                except Exception:
                    raw_pvals.append(np.nan)
            else:
                raw_pvals.append(np.nan)

    raw_pvals = np.array(raw_pvals, dtype=float)
    valid_p = np.isfinite(raw_pvals)
    adj_p = np.full_like(raw_pvals, np.nan)
    if valid_p.sum() > 0:
        _, corrected, _, _ = smt.multipletests(raw_pvals[valid_p], method="fdr_bh")
        adj_p[valid_p] = corrected

    def star_p(p):
        if np.isnan(p):
            return ""
        if p < 0.0001:
            return "****"
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    stars = [star_p(p) for p in adj_p]

    plot_df = clr_sub[selected_cts].copy()
    plot_df[group_col] = s_grp.values
    plot_df[sample_col] = plot_df.index
    melted = pd.melt(
        plot_df, id_vars=[sample_col, group_col], value_vars=selected_cts,
        var_name="cell_type", value_name="clr_abundance"
    )

    fig_w = max(8.0, len(selected_cts) * 0.75 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))

    palette = sns.color_palette("tab10", n_colors=max(1, len(uniq_levels)))
    sns.boxplot(
        data=melted, x="cell_type", y="clr_abundance", hue=group_col,
        ax=ax, showmeans=False, fliersize=0, boxprops=dict(alpha=0.6), palette=palette
    )
    sns.stripplot(
        data=melted, x="cell_type", y="clr_abundance", hue=group_col,
        ax=ax, dodge=True, jitter=0.2, alpha=0.8, size=4.5, edgecolor="black", linewidth=0.5,
        palette=palette, legend=False
    )

    y_max = float(melted["clr_abundance"].max()) if len(melted) else 1.0
    y_min = float(melted["clr_abundance"].min()) if len(melted) else -1.0
    y_range = y_max - y_min if y_max > y_min else 1.0
    star_y = y_max + 0.08 * y_range
    ax.set_ylim(y_min - 0.05 * y_range, y_max + 0.22 * y_range)

    for i, (ct, st, p) in enumerate(zip(selected_cts, stars, adj_p)):
        if st:
            color = "#D62728" if st != "ns" else "#7F7F7F"
            weight = "bold" if st != "ns" else "normal"
            ax.text(i, star_y, st, ha="center", va="bottom", fontsize=10, color=color, fontweight=weight)

    ax.set_title(title or f"Cell-Type Compositional CLR Abundance by {group_col}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Cell Type", fontweight="bold", fontsize=10)
    ax.set_ylabel("CLR Abundance", fontweight="bold", fontsize=10)
    plt.xticks(rotation=45, ha="right", fontsize=9)

    handles, labels = ax.get_legend_handles_labels()
    if len(uniq_levels) > 0:
        ax.legend(
            handles[: len(uniq_levels)], labels[: len(uniq_levels)],
            title=group_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True
        )
    fig.tight_layout()
    return fig


def compute_compositional_significance_matrix(
    clr_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_col: str,
    group_cols: list[str],
):
    """Compute FDR-corrected p-values for all cell types across multiple clinical/technical covariates.

    Returns (sig_matrix_df, heatmap_fig).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy.stats as stats
    import statsmodels.stats.multitest as smt

    dedup_group_cols = list(dict.fromkeys([c for c in group_cols if c in meta_df.columns]))
    meta_dedup = meta_df[[sample_col] + dedup_group_cols].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    common_samples = clr_df.index.intersection(meta_dedup.index)
    clr_sub = clr_df.loc[common_samples]
    meta_sub = meta_dedup.loc[common_samples]

    cts = clr_sub.columns.tolist()
    pvals_dict = {g: [] for g in dedup_group_cols if g in meta_sub.columns}

    for g in pvals_dict:
        s_grp = meta_sub[g]
        if isinstance(s_grp, pd.DataFrame):
            s_grp = s_grp.iloc[:, 0]
        s_grp = s_grp.astype(str)
        valid_mask = (s_grp != "<NA>") & (s_grp != "nan") & (s_grp != "None")
        valid_grp = s_grp[valid_mask]
        uniq_levels = valid_grp.unique()

        for ct in cts:
            vals = clr_sub[ct][valid_mask]
            if len(uniq_levels) < 2:
                pvals_dict[g].append(np.nan)
            elif len(uniq_levels) == 2:
                g1 = vals[valid_grp == uniq_levels[0]]
                g2 = vals[valid_grp == uniq_levels[1]]
                if len(g1) > 0 and len(g2) > 0:
                    try:
                        _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                        pvals_dict[g].append(p)
                    except Exception:
                        pvals_dict[g].append(np.nan)
                else:
                    pvals_dict[g].append(np.nan)
            else:
                grouped_vals = [vals[valid_grp == lvl] for lvl in uniq_levels if (valid_grp == lvl).sum() > 0]
                if len(grouped_vals) >= 2:
                    try:
                        _, p = stats.kruskal(*grouped_vals)
                        pvals_dict[g].append(p)
                    except Exception:
                        pvals_dict[g].append(np.nan)
                else:
                    pvals_dict[g].append(np.nan)

    sig_mat = pd.DataFrame(index=cts)
    for g in pvals_dict:
        raw_p = np.array(pvals_dict[g], dtype=float)
        valid_p_mask = np.isfinite(raw_p)
        adj_p = np.full_like(raw_p, np.nan)
        if valid_p_mask.sum() > 0:
            _, corrected, _, _ = smt.multipletests(raw_p[valid_p_mask], method="fdr_bh")
            adj_p[valid_p_mask] = corrected
        sig_mat[g] = -np.log10(np.clip(adj_p, 1e-10, 1.0))

    n_rows = len(sig_mat)
    n_cols = len(sig_mat.columns)
    max_label_len = max([len(str(x)) for x in sig_mat.index] + [10])
    label_margin = max(4.0, max_label_len * 0.14)
    fig_w = max(6.0, label_margin + n_cols * 1.8 + 1.5)
    fig_h = max(4.5, n_rows * 0.35 + 1.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    max_val = max(3.0, float(np.nanmax(sig_mat.values))) if not np.all(np.isnan(sig_mat.values)) else 3.0
    sns.heatmap(
        sig_mat, annot=True, fmt=".2f", cmap="YlOrRd", vmin=0, vmax=max_val,
        cbar_kws={"label": "-log10(FDR p-value)", "shrink": 0.6}, ax=ax
    )
    ax.set_title("Cell Type Compositional Shift Significance (-log10 FDR p)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Cell Type", fontsize=10, fontweight="bold")
    ax.set_xlabel("Covariate", fontsize=10, fontweight="bold")
    fig.tight_layout()

    return sig_mat, fig


def compute_variance_partition(
    clr_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_col: str,
    bio_cols: list[str],
    tech_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate the fraction of total variance in cell-type abundance explained by each covariate.

    Follows Sikkema et al. 2023 (Nature Medicine, Fig 4a) variance partitioning.
    Returns (var_tech_df, var_bio_df).
    """
    import statsmodels.formula.api as smf

    dedup_bio = list(dict.fromkeys([c for c in bio_cols if c in meta_df.columns]))
    dedup_tech = list(dict.fromkeys([c for c in tech_cols if c in meta_df.columns]))
    all_cols = list(dict.fromkeys(dedup_bio + dedup_tech))

    meta_dedup = meta_df[[sample_col] + all_cols].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    common_samples = clr_df.index.intersection(meta_dedup.index)
    clr_sub = clr_df.loc[common_samples]
    meta_sub = meta_dedup.loc[common_samples].copy()

    cts = ["Whole atlas"] + clr_sub.columns.tolist()
    var_bio = pd.DataFrame(index=cts, columns=dedup_bio, dtype=float)
    var_tech = pd.DataFrame(index=cts, columns=dedup_tech, dtype=float)

    for cdf, cols in [(var_tech, dedup_tech), (var_bio, dedup_bio)]:
        for col in cols:
            s_grp = meta_sub[col]
            if isinstance(s_grp, pd.DataFrame):
                s_grp = s_grp.iloc[:, 0]
            s_clean = s_grp.dropna().astype(str)
            s_clean = s_clean[(s_clean != "<NA>") & (s_clean != "nan") & (s_clean != "None")]
            if s_clean.nunique() < 2 or len(s_clean) < 3:
                continue

            df_reg = pd.DataFrame({"x": s_grp.loc[s_clean.index].astype("category")})
            r2_list = []
            for ct in clr_sub.columns:
                df_reg["y"] = clr_sub[ct].loc[s_clean.index].values
                try:
                    fit = smf.ols("y ~ C(x)", data=df_reg).fit()
                    r2 = max(0.0, float(fit.rsquared))
                except Exception:
                    r2 = 0.0
                cdf.loc[ct, col] = r2
                r2_list.append(r2)
            cdf.loc["Whole atlas", col] = np.mean(r2_list) if r2_list else np.nan

    return var_tech, var_bio


def plot_variance_partition_heatmap(
    var_tech: pd.DataFrame,
    var_bio: pd.DataFrame,
    title: str | None = None,
    vmax: float = 0.40,
) -> plt.Figure:
    """Plot dual-panel variance partition heatmap matching Sikkema et al. 2023 Fig 4a."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    n_rows = len(var_tech)
    n_tech = max(1, var_tech.shape[1])
    n_bio = max(1, var_bio.shape[1])

    max_label_len = max([len(str(x)) for x in var_tech.index] + [10])
    label_margin = max(2.5, max_label_len * 0.12)
    fig_w = max(7.5, label_margin + (n_tech + n_bio) * 0.65 + 1.2)
    fig_h = max(4.5, n_rows * 0.38 + 2.0)

    fig, (ax_tech, ax_bio) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [n_tech, n_bio], "wspace": 0.08},
        sharey=True,
    )

    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color="#7f7f7f")

    # Left: Technical covariates
    if var_tech.shape[1] > 0:
        sns.heatmap(
            var_tech, ax=ax_tech, cmap=cmap, vmin=0.0, vmax=vmax,
            cbar=False, annot=False, linewidths=0.5, linecolor="white"
        )
        ax_tech.set_xticklabels(ax_tech.get_xticklabels(), rotation=90, ha="center", fontsize=9)
    ax_tech.set_title("Covariate (technical)", fontsize=11, fontweight="bold", pad=10)
    ax_tech.set_xlabel("")
    ax_tech.set_ylabel("Cell type", fontsize=11, fontweight="bold")
    ax_tech.tick_params(left=True, bottom=True)

    # Right: Biological covariates
    if var_bio.shape[1] > 0:
        sns.heatmap(
            var_bio, ax=ax_bio, cmap=cmap, vmin=0.0, vmax=vmax,
            cbar=True, annot=False, linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Fraction of total variance", "shrink": 0.75, "extend": "max"}
        )
        ax_bio.set_xticklabels(ax_bio.get_xticklabels(), rotation=90, ha="center", fontsize=9)
    ax_bio.set_title("Covariate (biological)", fontsize=11, fontweight="bold", pad=10)
    ax_bio.set_xlabel("")
    ax_bio.set_ylabel("")
    ax_bio.tick_params(left=False, bottom=True)

    if title:
        fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)

    fig.tight_layout()
    return fig


def compute_compositional_joint_lm(
    clr_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    sample_col: str,
    group_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Fit a multi-variable OLS linear model (CLR ~ cov1 + cov2 + ...) across all cell types simultaneously.

    Returns (covariate_anova_summary_df, per_celltype_fits_df, top_model_summary_text),
    providing joint Type-II ANOVA covariate attribution, full model fit statistics,
    and the comprehensive statsmodels OLS summary for the top affected cell type.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    import statsmodels.stats.multitest as smt

    dedup_group_cols = list(dict.fromkeys([c for c in group_cols if c in meta_df.columns]))
    meta_dedup = meta_df[[sample_col] + dedup_group_cols].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    common_samples = clr_df.index.intersection(meta_dedup.index)
    clr_sub = clr_df.loc[common_samples]
    meta_sub = meta_dedup.loc[common_samples].copy()

    # Filter to covariates with >= 2 distinct non-NA categories
    valid_cols = []
    rename_map = {}
    for i, col in enumerate(dedup_group_cols):
        if col in meta_sub.columns:
            s_grp = meta_sub[col]
            if isinstance(s_grp, pd.DataFrame):
                s_grp = s_grp.iloc[:, 0]
            s_clean = s_grp.dropna().astype(str)
            s_clean = s_clean[(s_clean != "<NA>") & (s_clean != "nan") & (s_clean != "None")]
            if s_clean.nunique() >= 2:
                safe_name = f"cov_{i}"
                rename_map[col] = safe_name
                valid_cols.append(col)

    if not valid_cols or clr_sub.shape[1] == 0:
        return pd.DataFrame(), pd.DataFrame(), ""

    df_reg = pd.DataFrame(index=common_samples)
    for col in valid_cols:
        s_grp = meta_sub[col]
        if isinstance(s_grp, pd.DataFrame):
            s_grp = s_grp.iloc[:, 0]
        df_reg[rename_map[col]] = s_grp.astype("category")

    formula = "clr_val ~ " + " + ".join([f"C({rename_map[c]})" for c in valid_cols])

    fit_summaries = []
    anova_rows = []
    models_dict = {}

    for ct in clr_sub.columns:
        df_reg["clr_val"] = clr_sub[ct].values
        try:
            model = smf.ols(formula, data=df_reg).fit()
            r2 = float(model.rsquared) if np.isfinite(model.rsquared) else 0.0
            adj_r2 = float(model.rsquared_adj) if np.isfinite(model.rsquared_adj) else 0.0
            f_stat = float(model.fvalue) if np.isfinite(model.fvalue) else 0.0
            p_val = float(model.f_pvalue) if np.isfinite(model.f_pvalue) else 1.0

            fit_summaries.append({
                "Cell Type": ct,
                "R-squared": round(r2, 4),
                "Adj R-squared": round(adj_r2, 4),
                "F-statistic": round(f_stat, 2),
                "Prob (F-stat)": p_val,
                "AIC": round(float(model.aic), 1),
                "BIC": round(float(model.bic), 1),
            })
            models_dict[ct] = model

            try:
                a_tbl = anova_lm(model, typ=2)
                for col in valid_cols:
                    safe = f"C({rename_map[col]})"
                    if safe in a_tbl.index:
                        f_cov = float(a_tbl.loc[safe, "F"]) if np.isfinite(a_tbl.loc[safe, "F"]) else 0.0
                        p_cov = float(a_tbl.loc[safe, "PR(>F)"]) if np.isfinite(a_tbl.loc[safe, "PR(>F)"]) else 1.0
                        ss_cov = float(a_tbl.loc[safe, "sum_sq"])
                        df_cov = int(a_tbl.loc[safe, "df"])
                        anova_rows.append({
                            "Cell Type": ct,
                            "Covariate": col,
                            "Sum of Squares": round(ss_cov, 3),
                            "df": df_cov,
                            "F-statistic": round(f_cov, 2),
                            "p-value": p_cov,
                        })
            except Exception:
                pass
        except Exception:
            pass

    df_fits = pd.DataFrame(fit_summaries)
    df_anova = pd.DataFrame(anova_rows)

    if df_anova.empty:
        return pd.DataFrame(), df_fits, ""

    # Multi-testing correction on ANOVA p-values
    for col in valid_cols:
        sub_c = df_anova[df_anova["Covariate"] == col]
        raw_p = sub_c["p-value"].values
        valid_p_mask = np.isfinite(raw_p)
        adj_p = np.full_like(raw_p, 1.0)
        if valid_p_mask.sum() > 0:
            _, corr, _, _ = smt.multipletests(raw_p[valid_p_mask], method="fdr_bh")
            adj_p[valid_p_mask] = corr
        df_anova.loc[df_anova["Covariate"] == col, "FDR p-value"] = adj_p
        df_anova.loc[df_anova["Covariate"] == col, "Significance"] = [
            "****" if p < 0.0001 else ("***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns")))
            for p in adj_p
        ]

    cov_rows = []
    for col in valid_cols:
        sub_c = df_anova[df_anova["Covariate"] == col]
        raw_p = sub_c["p-value"].values
        adj_p = sub_c["FDR p-value"].values
        n_sig_raw = int((raw_p < 0.05).sum())
        n_sig_fdr = int((adj_p < 0.05).sum())
        min_fdr = float(np.nanmin(adj_p)) if len(adj_p) else 1.0
        mean_f = float(sub_c["F-statistic"].mean()) if len(sub_c) else 0.0
        best_ct = sub_c.loc[sub_c["p-value"].idxmin(), "Cell Type"] if len(sub_c) else "-"

        signif_code = (
            "****" if min_fdr < 0.0001 else ("***" if min_fdr < 0.001 else ("**" if min_fdr < 0.01 else ("*" if min_fdr < 0.05 else "ns")))
        )
        cov_rows.append({
            "Covariate": col,
            "Joint Model Term": f"C({col})",
            "Mean F-stat": round(mean_f, 2),
            "Sig CTs (p < 0.05)": f"{n_sig_raw}/{len(clr_sub.columns)}",
            "Sig CTs (FDR < 0.05)": f"{n_sig_fdr}/{len(clr_sub.columns)}",
            "Min FDR p-value": f"{min_fdr:.3e}" if min_fdr < 0.01 else f"{min_fdr:.3f}",
            "Top Affected Cell Type": best_ct,
            "Joint Model Significance": signif_code,
        })

    df_cov_summary = pd.DataFrame(cov_rows)

    # Pick top affected cell type for detailed statsmodels OLS summary text
    top_ct = df_fits.sort_values(by="F-statistic", ascending=False)["Cell Type"].iloc[0] if not df_fits.empty else None
    top_summary_text = models_dict[top_ct].summary().as_text() if top_ct and top_ct in models_dict else ""

    return df_cov_summary, df_fits, top_summary_text


# Backward-compatible alias
compute_compositional_lm_summary = compute_compositional_joint_lm


def select_balanced_samples(
    obs: pd.DataFrame,
    sample_col: str,
    bio_col: str | None = None,
    batch_cols: list | None = None,
    target_samples: int = 25,
    min_cells: int = 200,
    seed: int = 0,
) -> list[str]:
    """Joint balanced sample allocation: min(N_samples, target_samples).

    Filters candidate samples with >= min_cells, then stratifies across bio_col
    and round-robins across batch_cols combinations.
    """
    rng = np.random.default_rng(seed)
    vc = obs[sample_col].value_counts(dropna=True)
    eligible = vc[vc >= min_cells].index.astype(str).tolist()
    if not eligible:
        eligible = vc.index.astype(str).tolist()
    if len(eligible) <= target_samples:
        return eligible

    sample_obs = obs[obs[sample_col].astype(str).isin(eligible)].copy()
    sample_obs[sample_col] = sample_obs[sample_col].astype(str)
    present_batch = [c for c in (batch_cols or []) if c in sample_obs.columns]

    if bio_col and bio_col in sample_obs.columns:
        bio_mode = sample_obs.groupby(sample_col)[bio_col].agg(
            lambda s: s.mode().iloc[0] if len(s.mode()) else "<NA>"
        ).to_dict()
    else:
        bio_mode = {s: "all" for s in eligible}

    if present_batch:
        modes = sample_obs.groupby(sample_col)[present_batch].agg(
            lambda s: s.mode().iloc[0] if len(s.mode()) else "<NA>"
        )
        batch_dict = {s: tuple(row) for s, (_, row) in zip(modes.index, modes.iterrows())}
    else:
        batch_dict = {s: () for s in eligible}

    groups = sorted(set(bio_mode.values()))
    base_per_group = target_samples // len(groups)
    rem = target_samples % len(groups)
    alloc = {g: base_per_group + (1 if i < rem else 0) for i, g in enumerate(groups)}

    selected = []
    for g in groups:
        grp_samples = [s for s in eligible if bio_mode.get(s) == g]
        if not grp_samples:
            continue
        buckets = {}
        for s in grp_samples:
            buckets.setdefault(batch_dict.get(s), []).append(s)
        for b in buckets.values():
            rng.shuffle(b)
        keys = list(buckets.keys())
        grp_selected = []
        i = 0
        while len(grp_selected) < alloc[g] and any(buckets.values()):
            k = keys[i % len(keys)]
            if buckets[k]:
                grp_selected.append(buckets[k].pop())
            i += 1
        selected.extend(grp_selected)

    remaining = [s for s in eligible if s not in selected]
    if len(selected) < target_samples and remaining:
        rng.shuffle(remaining)
        selected.extend(remaining[: target_samples - len(selected)])

    return selected[:target_samples]


# ---------------------------------------------------------------------------
# T3.1 -- sample-first, RAM-bounded subsetting (onboarding diagnostics only)
# ---------------------------------------------------------------------------

SUBSET_CONFIG = {
    # Number of samples to select (range 15-30; files with fewer keep all).
    "MAX_SAMPLES": 25,
    # Samples per biological group (range 3-6), before the global cap applies.
    "N_PER_BIO": 5,
    # Random per-sample cap (1000 cells/sample for new dataset onboarding; 500 for _debug).
    "MAX_CELLS_PER_SAMPLE": 1000,
    # Optional per-cell-type cap applied AFTER the concat.
    "MAX_CELLS_PER_CT": None,
    # Overall soft target: the subset is trimmed (stratified by sample) when it
    # exceeds this many cells.
    "CELLS_TARGET": 25_000,
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


def _fast_backed_csr_slice(h5_group, parent_pos: np.ndarray, n_vars: int):
    """Fast bulk slice of on-disk CSR matrix from an HDF5 group or dataset.

    Merges nearby row ranges to replace thousands of single-row HDF5 seeks
    with a small number of bulk contiguous reads.
    """
    import scipy.sparse as sp

    indptr_ds = h5_group["indptr"]
    data_ds = h5_group["data"]
    indices_ds = h5_group["indices"]

    indptr = np.asarray(indptr_ds)
    row_starts = indptr[parent_pos]
    row_ends = indptr[parent_pos + 1]
    row_lengths = row_ends - row_starts

    sub_indptr = np.zeros(len(parent_pos) + 1, dtype=np.int64)
    sub_indptr[1:] = np.cumsum(row_lengths)
    total_nnz = sub_indptr[-1]

    sub_data = np.empty(total_nnz, dtype=data_ds.dtype)
    sub_indices = np.empty(total_nnz, dtype=indices_ds.dtype)

    if total_nnz == 0:
        return sp.csr_matrix((sub_data, sub_indices, sub_indptr), shape=(len(parent_pos), n_vars))

    merged_intervals = []
    cur_start = row_starts[0]
    cur_end = row_ends[0]
    cur_rows = []

    for i in range(len(parent_pos)):
        rs = row_starts[i]
        re = row_ends[i]
        dst_s = sub_indptr[i]
        dst_e = sub_indptr[i + 1]
        if rs == re:
            continue
        if rs <= cur_end + 2_000_000:
            cur_end = max(cur_end, re)
            cur_rows.append((rs, re, dst_s, dst_e))
        else:
            merged_intervals.append((cur_start, cur_end, cur_rows))
            cur_start = rs
            cur_end = re
            cur_rows = [(rs, re, dst_s, dst_e)]
    if cur_rows:
        merged_intervals.append((cur_start, cur_end, cur_rows))

    for blk_start, blk_end, rows in merged_intervals:
        blk_data = data_ds[blk_start:blk_end]
        blk_indices = indices_ds[blk_start:blk_end]
        for rs, re, dst_s, dst_e in rows:
            offset_s = rs - blk_start
            offset_e = re - blk_start
            sub_data[dst_s:dst_e] = blk_data[offset_s:offset_e]
            sub_indices[dst_s:dst_e] = blk_indices[offset_s:offset_e]

    return sp.csr_matrix((sub_data, sub_indices, sub_indptr), shape=(len(parent_pos), n_vars))


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
    if not parent_positions:
        raise ValueError("subset_by_samples: no cells selected")
    parent_pos = np.sort(np.concatenate(parent_positions))

    # Fast slice on backed AnnData (orders of magnitude faster than scanpy default row seeks)
    extracted_fast = False
    if getattr(adata, "isbacked", False):
        try:
            import h5py
            import anndata as ad
            fpath = getattr(adata, "filename", None)
            if fpath and os.path.exists(fpath):
                with h5py.File(fpath, "r") as hf:
                    if "X" in hf and isinstance(hf["X"], h5py.Group) and "indptr" in hf["X"]:
                        sub_X = _fast_backed_csr_slice(hf["X"], parent_pos, adata.n_vars)
                        sub_obs = adata.obs.iloc[parent_pos].copy()
                        sub_var = adata.var.copy()
                        sub = ad.AnnData(X=sub_X, obs=sub_obs, var=sub_var)
                        extracted_fast = True
        except Exception:
            extracted_fast = False

    if not extracted_fast:
        sub = adata[parent_pos].to_memory()
    sub.obsm = {}
    sub.obsp = {}
    sub.varm = {}
    sub.varp = {}
    sub.uns = {}

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
    label_cols: list,
    out_dir: Path | None = None,
    name: str = "dataset",
    dpi: int = 150,
    point_size: float | None = None,
    seed: int = 0,
) -> list:
    """Scatter panels of the UMAP embedding colored by each label column.

    Returns a list of generated matplotlib Figure objects (and saves PNGs if out_dir is provided).
    """
    import matplotlib.pyplot as plt

    n_cells = umap_coords.shape[0]
    if point_size is None:
        point_size = 6.0 if n_cells <= 10_000 else (3.0 if n_cells <= 50_000 else 1.0)

    saved_figs = []
    for col in label_cols:
        if col not in obs_subset.columns:
            warnings.warn(f"plot_umap_panels: label column {col!r} not found; skipped")
            continue
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        cats = obs_subset[col].astype(str).fillna("NA")
        uniq_cats = pd.unique(cats)
        n_uniq = len(uniq_cats)
        for cat in uniq_cats:
            m = cats.values == cat
            ax.scatter(
                umap_coords[m, 0],
                umap_coords[m, 1],
                s=point_size,
                alpha=0.7,
                linewidths=0,
                label=cat if n_uniq <= 30 else None,
            )
        if n_uniq <= 30:
            ncol = 2 if n_uniq > 15 else 1
            fsize = 7 if n_uniq > 15 else 8
            ax.legend(
                bbox_to_anchor=(1.02, 1), loc="upper left",
                frameon=False, fontsize=fsize, markerscale=2,
                ncol=ncol, title=col,
            )
        ax.set_title(f"UMAP colored by {col}", fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        fig.tight_layout()

        if out_dir is not None:
            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)
            safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in col)
            p = out_path / f"{name}_umap_{safe}.png"
            fig.savefig(p, dpi=dpi, bbox_inches="tight")
        saved_figs.append(fig)
    return saved_figs


def compute_embed_umap(
    adata,
    compute_subset_max: int = 100_000,
    pca_n_comps: int = 50,
    hvg_n_top: int = 2000,
    seed: int = 0,
    verbose: bool = True,
):
    """Compute a fresh unintegrated PCA+UMAP on a subsample of RAW counts.

    Mirrors the standard preprocessing pipeline in src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:
    1. Extracts raw counts (adata.raw, layers['counts'], or X).
    2. Filter cells (min_genes=100) & filter genes (min_cells=3).
    3. Standard total-count normalization (sc.pp.normalize_total(target_sum=1e4)) + log1p.
    4. HVG selection (n_top_genes=2000, flavor='seurat_v3_paper' with jitter fallback).
    5. Scaling (max_value=10) and PCA (n_comps=50, svd_solver='arpack').
    6. Neighborhood graph (n_pcs=50, n_neighbors=15).
    7. UMAP embedding (Scanpy defaults: min_dist=0.5, spread=1.0).
    """
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    n_cells = adata.n_obs
    compute_max = min(n_cells, compute_subset_max)
    subset_idx = rng.choice(n_cells, size=compute_max, replace=False)
    subset_idx = np.sort(subset_idx)
    sub = adata[subset_idx].to_memory()
    sub.obsm = {}
    sub.obsp = {}
    sub.varm = {}
    sub.varp = {}
    sub.uns = {}

    counts, slot = locate_counts(sub, verbose=False)
    if slot.startswith("raw."):
        common_vars = sub.var_names.intersection(sub.raw.var_names)
        if len(common_vars) > 0:
            counts = sub.raw[:, common_vars].X
            sub = sub[:, common_vars].copy()
        else:
            counts = sub.raw.X
    sub.X = counts.copy() if hasattr(counts, "copy") else counts
    sub.layers["counts"] = sub.X.copy()

    if not sp.issparse(sub.X):
        sub.X = sp.csr_matrix(sub.X)
    if not sp.issparse(sub.layers["counts"]):
        sub.layers["counts"] = sp.csr_matrix(sub.layers["counts"])

    if sub.n_obs > 20 and sub.n_vars > 100:
        sc.pp.filter_cells(sub, min_genes=100)
        sc.pp.filter_genes(sub, min_cells=3)

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            sc.pp.highly_variable_genes(
                sub, layer="counts", n_top_genes=hvg_n_top, flavor="seurat_v3_paper", check_values=True
            )
        except Exception:
            try:
                sc.pp.highly_variable_genes(sub, n_top_genes=hvg_n_top, flavor="seurat", subset=False)
            except Exception:
                sub.var["highly_variable"] = True

    hvg_mask = sub.var["highly_variable"].values
    if not hvg_mask.any():
        hvg_mask = np.ones(len(sub.var), dtype=bool)

    sub_hvg = sub[:, hvg_mask].copy()
    sc.pp.scale(sub_hvg, max_value=10)
    if sp.issparse(sub_hvg.X):
        sub_hvg.X.data = np.nan_to_num(sub_hvg.X.data, nan=0.0, posinf=10.0, neginf=-10.0)
    elif isinstance(sub_hvg.X, np.ndarray):
        sub_hvg.X = np.nan_to_num(sub_hvg.X, nan=0.0, posinf=10.0, neginf=-10.0)
    n_comps = min(pca_n_comps, sub_hvg.n_vars - 1, sub_hvg.n_obs - 1)
    try:
        sc.pp.pca(sub_hvg, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    except Exception:
        sc.pp.pca(sub_hvg, n_comps=n_comps, svd_solver="randomized", random_state=seed)
    sub.obsm["X_pca"] = sub_hvg.obsm["X_pca"]

    n_pcs_use = min(n_comps, sub.obsm["X_pca"].shape[1])
    sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=15, n_pcs=n_pcs_use)
    sc.tl.umap(sub, min_dist=0.5, spread=1.0, random_state=seed)

    if verbose:
        print(f"computed unintegrated PCA+UMAP on {sub.n_obs} cells (raw counts, {slot})")
    return "X_pca", "X_umap", sub


def embed_and_umap_workflow(
    adata,
    label_cols: list,
    out_dir: Path | None = None,
    name: str = "dataset",
    sample_col: str | None = None,
    compute_subset_max: int = 100_000,
    seed: int = 0,
    verbose: bool = True,
) -> dict:
    """Driver: always compute fresh unintegrated PCA+UMAP on raw counts and generate panels."""
    pca_key, umap_key, sub = compute_embed_umap(
        adata, compute_subset_max=compute_subset_max, seed=seed, verbose=verbose
    )
    obs_subset = sub.obs.reset_index(drop=True)
    coords = np.asarray(sub.obsm[umap_key])

    effective_cols = list(label_cols)
    if sample_col and sample_col in obs_subset.columns and sample_col not in effective_cols:
        effective_cols.insert(0, sample_col)

    figs = plot_umap_panels(coords, obs_subset, effective_cols, out_dir=out_dir, name=name, seed=seed)
    return {
        "pca_key": pca_key,
        "umap_key": umap_key,
        "subset": sub,
        "figs": figs,
        "computed": True,
    }


def write_metrics_input(
    adata,
    out_path,
    ct_col,
    bio_col,
    batch_cols,
    sample_col: str | None = None,
    max_cells: int = 300_000,
    seed: int = 0,
    n_pcs: int = 50,
    verbose: bool = True,
) -> dict:
    """Write the feather consumed by onboarding_metrics.R on unintegrated PCA."""
    if "X_pca" in adata.obsm:
        sub = adata
        pca_key = "X_pca"
    elif "X_pca_onboard" in adata.obsm:
        sub = adata
        pca_key = "X_pca_onboard"
    else:
        pca_key, _umap_key, sub = compute_embed_umap(adata, seed=seed, verbose=verbose)

    rng = np.random.default_rng(seed)
    max_cells = min(sub.n_obs, max_cells)
    idx = np.sort(rng.choice(sub.n_obs, size=max_cells, replace=False))
    emb = np.asarray(sub.obsm[pca_key][idx])[:, :n_pcs]
    obs_used = sub.obs.iloc[idx]

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
        print(f"wrote metrics input to {out_path}: {df.shape[0]} cells x {df.shape[1]} cols")
    return {"path": str(out_path), "n_cells": int(len(df)), "pc_source": "fresh unintegrated X_pca"}


def run_onboarding_metrics(
    sub,
    ct_col: str,
    bio_col: str,
    batch_cols: list,
    sample_col: str | None = None,
    temp_dir: Path | None = None,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Runs cell-level LISI separation metrics and returns (sep_df, json_summary)."""
    import json
    import shutil
    import subprocess
    import tempfile

    if temp_dir is None:
        t_dir = Path(tempfile.mkdtemp(prefix="onboard_metrics_"))
        cleanup = True
    else:
        t_dir = Path(temp_dir)
        t_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False

    feather_p = t_dir / "metrics_input.feather"
    csv_p = t_dir / "metrics_separation.csv"
    json_p = t_dir / "metrics_separation.json"

    write_metrics_input(
        sub, out_path=feather_p, ct_col=ct_col, bio_col=bio_col,
        batch_cols=batch_cols, sample_col=sample_col, seed=seed, verbose=False
    )

    here = Path(__file__).resolve().parent
    script_path = here / "onboarding_metrics.R"
    pixi = shutil.which("pixi") or f"{Path.home()}/.pixi/bin/pixi"
    cmd = [
        pixi, "run", "-e", "default", "Rscript", "--vanilla",
        str(script_path),
        "--input", str(feather_p),
        "--ct-col", str(ct_col),
        "--bio-col", str(bio_col),
        "--batch-cols", ",".join(batch_cols or []),
        "--out-csv", str(csv_p),
        "--out-json", str(json_p),
        "--seed", str(seed),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    sep_df = pd.read_csv(csv_p)
    with open(json_p) as f:
        meta_json = json.load(f)

    if cleanup:
        shutil.rmtree(t_dir, ignore_errors=True)

    return sep_df, meta_json


def plot_separation_heatmap(sep_df: pd.DataFrame, name: str = "Dataset", out_path: Path | None = None):
    """Generate and return a seaborn heatmap Figure for cell-type x label separation."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    import re

    mat = sep_df.set_index("cell_type")
    score_cols = [c for c in mat.columns if c.endswith("_separation")]
    lab = {c: re.sub("_separation$", "", c).replace("batch_", "") for c in score_cols}
    hm = mat[score_cols].rename(columns=lab)

    n_rows = len(hm)
    n_cols = len(score_cols)
    max_label_len = max([len(str(x)) for x in hm.index] + [10])
    label_margin = max(4.0, max_label_len * 0.14)
    col_width = 1.8
    fig_w = label_margin + n_cols * col_width + 1.8
    fig_h = max(5.0, n_rows * 0.38 + 1.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        hm.astype(float),
        annot=True,
        fmt=".2f",
        annot_kws={"size": 9, "weight": "bold"},
        cmap="RdYlGn_r",
        vmin=0,
        vmax=1,
        mask=hm.isna(),
        linewidths=0.5,
        cbar_kws={"shrink": 0.5, "aspect": 20, "label": "LISI Separation (1=Separated, 0=Mixed)"},
        ax=ax,
    )
    ax.set_title(f"{name}: Per-Cell-Type LISI Separation (Unintegrated PCA)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Cell Type", fontsize=10, fontweight="bold")
    ax.set_xlabel("Label / Covariate", fontsize=10, fontweight="bold")
    fig.tight_layout()

    if out_path is not None:
        out_p = Path(out_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_p, dpi=150, bbox_inches="tight")
    return fig


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