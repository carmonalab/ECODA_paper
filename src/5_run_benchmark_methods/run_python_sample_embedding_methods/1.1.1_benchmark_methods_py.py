"""Python benchmark methods (MrVI, scPoli, PILOT, QOT, PILOT-GM-VAE) as a CLI script.

Replaces the logic of the archived notebook
`1.2_benchmark_methods_py.qmd` (kept as reference; do NOT delete). Consumes
the preprocessed benchmark view h5ad produced by
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`:

- PILOT/QOT/PILOT-GM-VAE consume the stored obsm embedding `X_pca_{view}_hvg{n}`
  (the qmd recomputed PCA from scratch);
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
`process_scpoli_fig` / `process_pilot_fig` / `process_qot_fig` /
`process_pilotgm_fig`, `constants.R` label map and the notebook recodes
depend on them): plain `DataFrame.to_feather()` with the pandas index (sample
names) kept — the index is written as the last feather column, matching R's
`column_to_rownames(ncol)`.

QOT and PILOT-GM-VAE are extended methods (see the implementation plan
`.kilo/plans/1786651957910-pilotgm-qot-benchmark-implementation.md`): QOT runs
the vendored `qot_utils_re.py` (PennShenLab/QOT @ 28cd529880c1, one bug fix
in `Gaussian_Mixture_Representation`), PILOT-GM-VAE runs the `pilotgm` PyPI
package (CostaLab/PILOT-GM-VAE, BIB 2025). Both receive a distinct temp obs
column (`_bench_prog` / `_bench_status`) instead of `"Sample"` for the
status/progession argument: their `rename()` dicts collapse when the sample
and status column are the same (duplicate dict key -> both columns renamed
to 'status', no 'sampleID' survives -> KeyError in the GMM groupby). The
temp column also keeps any bio label out of the distance path (no-leakage).

Execution time (float seconds, method body only — excluding h5ad loading,
like the qmd and R `exec_time()`) and peak RSS (`mem_GB`) are appended to a
per-task log feather. One process per task writes the file, so no concurrency
issues. Combos run defaults-first (MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
PILOT_hvg2000_highres, QOT_hvg2000_highres, PILOT-GM-VAE_hvg2000_highres) so
the main-method rows are measured before any in-process memory bloat (peak
RSS is monotonic within a process).
"""

import argparse
import gc
import hashlib
import os
import resource
import sys
import tempfile
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
from src.utils.py.benchmark_h5ad_contract import validate_benchmark_h5ad_contract


def _file_md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_feather_frame(frame, path):
    if frame.empty:
        raise ValueError(f"Feather artifact is empty: {path}")
    if isinstance(frame.index, pd.RangeIndex):
        raise ValueError(f"Feather artifact lacks sample identifiers: {path}")
    if frame.index.hasnans or not frame.index.is_unique:
        raise ValueError(f"Feather artifact has invalid sample identifiers: {path}")
    if any(str(value).strip() == "" for value in frame.index):
        raise ValueError(f"Feather artifact has blank sample identifiers: {path}")
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.empty or not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"Feather artifact has no finite numeric features: {path}")


def _ordered_sample_ids(adata):
    """Return unique Sample IDs in their first-appearance obs order."""
    if "Sample" not in adata.obs.columns:
        raise ValueError("benchmark AnnData is missing the 'Sample' column")
    values = adata.obs["Sample"]
    if values.isna().any():
        raise ValueError("benchmark AnnData Sample contains missing values")
    sample_ids = [str(value) for value in values]
    if any(not value.strip() for value in sample_ids):
        raise ValueError("benchmark AnnData Sample contains blank values")
    return list(dict.fromkeys(sample_ids))


def _align_square_frame(frame, sample_ids, path):
    """Align a sample-by-sample frame to the canonical obs sample order."""
    expected = [str(value) for value in sample_ids]
    row_ids = [str(value) for value in frame.index]
    column_ids = [str(value) for value in frame.columns]
    if (
        len(row_ids) != len(expected)
        or len(column_ids) != len(expected)
        or len(set(row_ids)) != len(row_ids)
        or len(set(column_ids)) != len(column_ids)
        or set(row_ids) != set(expected)
        or set(column_ids) != set(expected)
    ):
        raise ValueError(
            f"sample-by-sample Feather output identifiers do not match "
            f"canonical samples: {path}"
        )
    aligned = frame.copy()
    aligned.index = row_ids
    aligned.columns = column_ids
    return aligned.loc[expected, expected]


def recorded_feather_valid(path):
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or path.stat().st_size == 0 or not sidecar.is_file():
        return False
    records = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            records[key] = value
    if (
        records.get("PATH") != str(path)
        or records.get("MD5") != _file_md5(path)
        or records.get("SIZE") != str(path.stat().st_size)
    ):
        return False
    try:
        frame = pd.read_feather(path)
        _validate_feather_frame(frame, path)
        return True
    except Exception:
        return False


def atomic_to_feather(frame, path):
    """Write a complete Feather file and checksum before publication."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{path}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    backup = path.with_name(f".{path.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    had_path = path.is_file()
    had_sidecar = sidecar.is_file()
    try:
        _validate_feather_frame(frame, path)
        frame.to_feather(tmp)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"empty Feather output: {tmp}")
        if had_path:
            os.link(path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(tmp, path)
        sidecar_tmp.write_text(
            f"MD5={_file_md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n"
        )
        os.replace(sidecar_tmp, sidecar)
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
        for temporary in (tmp, sidecar_tmp, backup, sidecar_backup):
            if temporary.exists():
                temporary.unlink()

def _validate_execution_measurements(frame, path):
    for column, allow_missing in (("time_secs", False), ("mem_GB", True)):
        raw = frame[column]
        values = pd.to_numeric(raw, errors="coerce")
        missing = raw.isna()
        if (not allow_missing and missing.any()) or values[~missing].isna().any():
            raise ValueError(f"execution log has invalid numeric values: {path}")
        if not np.isfinite(values[~missing].to_numpy(dtype=float)).all():
            raise ValueError(f"execution log has invalid numeric values: {path}")


def _validate_execution_log_frame(frame, path):
    required = {"dataset", "method", "time_secs", "mem_GB"}
    if set(frame.columns) != required or frame.empty:
        raise ValueError(f"execution log has an invalid schema: {path}")
    if frame[["dataset", "method"]].isna().any().any() or \
       frame[["dataset", "method"]].astype(str).apply(lambda column: column.str.strip() == "").any().any():
        raise ValueError(f"execution log has blank identifiers: {path}")
    if frame[["dataset", "method"]].duplicated().any():
        raise ValueError(f"execution log has duplicate identifiers: {path}")
    _validate_execution_measurements(frame, path)

def execution_log_atomic_to_feather(frame, path):
    """Atomically write the documented indexless execution-log schema."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{path}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    backup = path.with_name(f".{path.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    had_path = path.is_file()
    had_sidecar = sidecar.is_file()
    try:
        _validate_execution_log_frame(frame, path)
        frame.reset_index(drop=True).to_feather(tmp)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"empty execution log: {tmp}")
        if had_path:
            os.link(path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(tmp, path)
        sidecar_tmp.write_text(
            f"MD5={_file_md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n"
        )
        os.replace(sidecar_tmp, sidecar)
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
        for temporary in (tmp, sidecar_tmp, backup, sidecar_backup):
            if temporary.exists():
                temporary.unlink()


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


def get_pilotgm():
    """Lazy import of pilotgm with a shim for its non-relative internal imports.

    Upstream packaging bug (pilotgm 0.1.1): `pilotgm/model/GMVAE.py` does
    `from networks.Networks import *` (also `losses`/`metrics`), but the
    modules ship inside the package (`pilotgm/networks/...`), so a plain
    `import pilotgm` fails with ModuleNotFoundError. The shim inserts the
    pilotgm package directory into sys.path so the top-level names resolve.

    The entry is left on sys.path ON PURPOSE: gmmvae_wasserstein_distance
    dispatches `compute_emd` to loky workers via joblib, and unpickling the
    task re-imports pilotgm in the worker. loky spawns copy the parent's
    sys.path, so removing the entry after the import would break the worker
    import ("failed to un-serialize" / ModuleNotFoundError: networks).
    Collision check (py-cpu env, 2026-08-14): no other installed package
    imports `core`, `model`, `networks`, `losses` or `metrics` top-level, so
    the extra path entry is inert.
    """
    import importlib.metadata as md

    for f in md.files("pilotgm"):
        if str(f).endswith("pilotgm/__init__.py"):
            pkg_dir = f.locate().parent
            break
    if pkg_dir is None:
        raise ImportError("pilotgm package files not found (importlib.metadata)")
    if str(pkg_dir) not in sys.path:
        sys.path.insert(0, str(pkg_dir))
    import pilotgm

    return pilotgm


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

def report_gpu_memory(method_str):
    """Print peak CUDA memory and device capacity for resource profiling."""
    if not torch.cuda.is_available():
        return
    try:
        torch.cuda.synchronize()
        props = torch.cuda.get_device_properties(torch.cuda.current_device())
        gib = float(1024**3)
        print(
            f"GPU_MEMORY_PROFILE method={method_str} "
            f"device={props.name} total_GiB={props.total_memory / gib:.2f} "
            f"peak_allocated_GiB={torch.cuda.max_memory_allocated() / gib:.2f} "
            f"peak_reserved_GiB={torch.cuda.max_memory_reserved() / gib:.2f}",
            flush=True,
        )
    except Exception as exc:
        print(f"GPU_MEMORY_PROFILE_ERROR method={method_str}: {exc}", flush=True)


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
    execution_log_atomic_to_feather(df_final, log_file)


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


def run_wass_combo_for(n_hvg, res_label):
    """Whether a Wasserstein-distance method (PILOT/QOT/PILOT-GM-VAE) runs
    for an (n_hvg, resolution) combo (qmd rules)."""
    if res_label == "_highres":
        return True
    return res_label == "_lowres" and n_hvg == 2000


# Default (main-method) combos — constants.R method_label_map_main and the
# notebook's exec-time figure: MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
# PILOT_hvg2000_highres, QOT_hvg2000_highres, PILOT-GM-VAE_hvg2000_highres.
DEFAULT_HVG = 2000
DEFAULT_SCPOLI_DIM = 15
DEFAULT_RES_LABEL = "_highres"


def is_default_combo(method, combo):
    n, res_label, _, payload, _, _ = combo
    if method == "mrvi":
        return n == DEFAULT_HVG
    if method == "scpoli":
        return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL and payload == DEFAULT_SCPOLI_DIM
    if method in ("pilot", "qot", "pilotgm"):
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
    """Point X at the validated raw counts layer for MrVI/scPoli.

    The preprocessed h5ad has log-normalized X; both models need the raw
    counts vaulted by base_preprocessing in layers["counts"]. The contract is
    checked before this function runs, so missing counts are an error rather
    than a normalized-expression fallback.
    """
    counts = sub.layers["counts"]
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
def run_mrvi(adata, device, output_path, batch_key=None):
    """MrVI local sample distances with an optional technical batch key."""
    adata.obs["dummy_col"] = np.zeros(adata.n_obs)
    setup_kwargs = {"sample_key": "Sample"}
    if batch_key is not None:
        setup_kwargs["batch_key"] = batch_key
    MRVI.setup_anndata(adata, **setup_kwargs)
    model = MRVI(adata)
    model.train(max_epochs=50, accelerator=device)
    dists = model.get_local_sample_distances(
        keep_cell=False, groupby="dummy_col", batch_size=32
    )
    df_dists = dists["dummy_col"].isel(dummy_col_name=0).to_pandas()
    df_dists = _align_square_frame(
        df_dists, _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path)


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
    atomic_to_feather(df_embs, output_path)


def resolve_pass_embedding_key(adata, view, n_hvg):
    """Resolve one explicitly declared embedding without compatibility fallback."""
    if view == "benchmark_analysis":
        key = f"X_pca_benchmark_analysis_hvg{n_hvg}"
    elif view == "batch_effect_uncorrected":
        key = f"X_pca_batch_effect_uncorrected_hvg{n_hvg}"
    elif view == "batch_effect_corrected":
        key = f"X_pca_harmony_batch_effect_corrected_hvg{n_hvg}"
    else:
        raise ValueError(f"Unknown preprocessing view: {view}")
    if key not in adata.obsm:
        raise KeyError(f"Required embedding {key!r} is missing from adata.obsm")
    return key


def run_pilot(adata, ct_col, view, n_hvg, output_path):
    """PILOT Wasserstein sample distances on the exact view embedding."""
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
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
    # Same NaN cell-type guard as QOT/PILOT-GM-VAE/scPoli (see run_scpoli):
    # PILOT's cost_matrix() treats NaN labels as a pseudo-cell-type whose
    # centroid is the median of zero cells -> NaN entries in the cost matrix
    # -> ot.emd2 collapses (all-zero EMD distances written on HPC; segfault
    # locally). Lee/Zhang (HiTME immune-only annotation, 12-33% NaN) hit this.
    fill_unknown_ct(adata, ct_col, "PILOT")
    # pilotpy writes hard-coded Results_PILOT/plots relative to cwd.  A
    # read-only project bind cannot host those scratch plots, and they are not
    # benchmark artifacts; isolate only this third-party side effect.
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="pilot_") as temp_dir:
        try:
            os.chdir(temp_dir)
            pl.tl.wasserstein_distance(
                adata,
                emb_matrix=emb_key,
                clusters_col=ct_col,
                sample_col="Sample",
                status="Sample",
            )
        finally:
            os.chdir(cwd)
    df_dists = _align_square_frame(
        adata.uns["EMD_df"], _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path)


def fill_unknown_ct(adata, ct_col, method):
    """Fill NaN cell-type labels with an explicit "Unknown" class.

    Same rationale as scPoli (see run_scpoli): the QOT script filters
    `Cell_type != 'Unknown'` and NaN would otherwise leak into the GMM
    groupby keys. No-op on complete columns.
    """
    n_na = int(adata.obs[ct_col].isna().sum())
    if n_na:
        print(f"Filling {n_na}/{adata.n_obs} missing values in '{ct_col}' "
              f"with 'Unknown' ({method}).")
        col = adata.obs[ct_col]
        if isinstance(col.dtype, pd.CategoricalDtype):
            col = col.cat.add_categories("Unknown")
        adata.obs[ct_col] = col.fillna("Unknown")


def run_qot(adata, ct_col, view, n_hvg, output_path):
    """QOT Wasserstein sample distances on the preprocessed obsm PCA.

    Runs the vendored qot_utils_re.py (PennShenLab/QOT @ 28cd529880c1, two
    hotfixes in Gaussian_Mixture_Representation — see the file header and
    docs/qot_hotfixes.md). Lazy import so the phate dependency is only
    touched for QOT runs.
    """
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
    import qot_utils_re

    fill_unknown_ct(adata, ct_col, "QOT")
    # Distinct temp column: Extract_Info renames {type_cell, id, progession}
    # in ONE dict — passing "Sample" for both id and progession collapses the
    # duplicate dict key, BOTH columns get renamed to 'status', no 'sampleID'
    # survives and the GMM groupby raises KeyError. The temp column also
    # keeps bio labels out of the distance path.
    adata.obs["_bench_prog"] = adata.obs["Sample"]
    qot_utils_re.Run_QOT(
        adata,
        gene_matrix=emb_key,
        type_cell=ct_col,
        id_col="Sample",
        progession="_bench_prog",
        dataset_type="rna",
        num_components_list=[1],
        random_state=2,
        min_samples_for_gmm=0,
        qot_method="cosine",
        normalized_set=False,
    )
    samples = _ordered_sample_ids(adata)
    # Plain object strings: anndata obs columns are categorical by default and
    # a categorical DataFrame index would be written to the feather as
    # categorical (pyarrow 24/25 pandas-compat cannot read that back:
    # "data type 'categorical' not understood"). PILOT's EMD_df index is
    # plain object — keep the identical layout.
    samples = np.asarray(samples, dtype=object)
    df_dists = pd.DataFrame(
        adata.uns["QOT_Distance"], index=samples, columns=samples
    )
    df_dists = _align_square_frame(
        df_dists, samples, output_path
    )
    atomic_to_feather(df_dists, output_path)


def _stabilize_pilotgm_covariances(adata):
    """Regularize undefined empirical covariances before PILOT-GM distances.

    The upstream representation helper uses ``np.cov`` for each assigned
    component. A component with one cell has an undefined covariance and
    yields NaNs; interpreting that component as zero empirical covariance,
    followed by the package's diagonal regularization, is the finite
    single-observation limit. Means and weights remain untouched and must
    already be finite.
    """
    representations = adata.uns.get("GMVAE_Representation")
    if not isinstance(representations, dict) or not representations:
        raise ValueError("PILOT-GM-VAE did not produce sample representations")
    for sample_id, representation in representations.items():
        means = np.asarray(representation["means"], dtype=float)
        weights = np.asarray(representation["weights"], dtype=float)
        if (
            not np.isfinite(means).all()
            or not np.isfinite(weights).all()
            or weights.ndim != 1
            or weights.size == 0
            or float(weights.sum()) <= 0
        ):
            raise ValueError(
                f"PILOT-GM-VAE produced invalid means/weights for sample {sample_id!r}"
            )
        covariances = np.asarray(representation["covariances"], dtype=float)
        if (
            covariances.ndim != 3
            or covariances.shape[1] != covariances.shape[2]
            or covariances.shape[0] != means.shape[0]
        ):
            raise ValueError(
                f"PILOT-GM-VAE produced invalid covariances for sample {sample_id!r}"
            )
        stable = []
        for covariance in covariances:
            if np.isinf(covariance).any():
                raise ValueError(
                    f"PILOT-GM-VAE produced infinite covariance for sample {sample_id!r}"
                )
            covariance = np.nan_to_num(covariance, nan=0.0)
            covariance = (covariance + covariance.T) / 2.0
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.clip(eigenvalues, 0.0, None)
            stable.append((eigenvectors * eigenvalues) @ eigenvectors.T)
        representation["covariances"] = np.asarray(stable, dtype=float)


def _run_pilotgm_distance_with_stable_covariances(pilotgm, adata, *, emb_key):
    """Run PILOT-GM after finite PSD covariance stabilization."""
    import importlib

    pilotgm_core = importlib.import_module("pilotgm.core")
    original_representation = pilotgm_core.gaussian_mixture_vae_representation

    def stabilized_representation(*args, **kwargs):
        result = original_representation(*args, **kwargs)
        _stabilize_pilotgm_covariances(result)
        return result

    pilotgm_core.gaussian_mixture_vae_representation = stabilized_representation
    try:
        pilotgm.gmmvae_wasserstein_distance(
            adata,
            emb_matrix=emb_key,
            clusters_col="component_assignment",
            sample_col="Sample",
            status="_bench_status",
            wass_dis=True,
            covariance_type="full",
            epsilon=1e-3,
        )
    finally:
        pilotgm_core.gaussian_mixture_vae_representation = original_representation


def run_pilotgm(adata, ct_col, view, n_hvg, output_path, ds_name, device):
    """PILOT-GM-VAE Wasserstein sample distances on the preprocessed obsm PCA.

    Runs the `pilotgm` PyPI package (CostaLab/PILOT-GM-VAE, BIB 2025):
    `train_gmvae` (50 epochs, num_classes = n unique cell types) then
    `gmmvae_wasserstein_distance`. The whole pilotgm block runs inside a
    node-local tempdir: `train_gmvae` hardcodes `./trained_models/<ds>/` and
    saves weights — running it from the repo root (worker cwd on HPC) would
    pollute the repo, and running it from the scratch output dir would
    pollute the NAS sync (the submit tail rsyncs benchmark/ wholesale).
    Weights are ephemeral by design (load_weights=False; retries re-train
    from scratch).
    """
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
    emb = adata.obsm[emb_key]
    # train_gmvae needs torch.tensor(obsm[key]) (fails on a pandas
    # DataFrame: "could not determine the shape of object type 'DataFrame'"
    # with torch >= 2.x), while extract_data_anno_scRNA_from_h5ad (in the
    # distance step) needs `.columns`. Store the plain ndarray for training,
    # swap in the named-columns DataFrame only for
    # gmmvae_wasserstein_distance — a DataFrame at that point is also
    # joblib-picklable (a __main__-defined ndarray subclass broke loky's
    # task serialization).
    if hasattr(emb, "columns"):
        emb = np.asarray(emb)
    adata.obsm[emb_key] = emb

    # Unlike PILOT/QOT, PILOT-GM's distance API uses the model-generated
    # `component_assignment` column, not the biological cell-type labels. The
    # configured annotation only preserves the historical component-count
    # choice; do not turn missing labels into an extra rare component.
    ct_values = adata.obs[ct_col].astype("string")
    valid_ct = ct_values.notna() & (ct_values.str.strip() != "")
    num_classes = max(2, int(ct_values[valid_ct].nunique()))
    # Distinct temp column for `status`: same duplicate-key rename bug as QOT
    # (gmmvae_wasserstein_distance renames the last three columns via a dict;
    # sample_col == status == "Sample" would collapse the keys).
    adata.obs["_bench_status"] = adata.obs["Sample"]

    pilotgm = get_pilotgm()
    # Plain `device == "cuda"` would silently run CPU on GPU nodes under the
    # default --device auto.
    use_cuda = device == "cuda" or (device == "auto" and torch.cuda.is_available())

    cwd = os.getcwd()
    tmp_dir = tempfile.mkdtemp(prefix="pilotgm_")
    try:
        os.chdir(tmp_dir)
        pilotgm.train_gmvae(
            adata,
            dataset_name=ds_name,
            pca_key=emb_key,
            labels_column=None,
            epochs=50,
            num_classes=num_classes,
            cuda=use_cuda,
            gpuID=0,
            load_weights=False,
            save_model=True,
            seed=1,
        )
        # Swap in the named-columns DataFrame for the distance step (see
        # the comment at the top of this function).
        adata.obsm[emb_key] = pd.DataFrame(
            emb,
            index=adata.obs_names,
            columns=[f"PCA_{i + 1}" for i in range(emb.shape[1])],
        )
        # The configured cell-type column determines the requested number of
        # mixture components, but the package's distance step uses the
        # model-generated component_assignment column, not those annotations.
        _run_pilotgm_distance_with_stable_covariances(
            pilotgm,
            adata,
            emb_key=emb_key,
        )
    finally:
        os.chdir(cwd)
    df_dists = _align_square_frame(
        adata.uns["EMD_df"], _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def process_dataset(args, ds_name, entry):
    """Run all combos of the requested method for one dataset.

    Loads the h5ad once per task (not per combo); skips combos whose output
    feather already exists unless --force.
    """
    view_name = args.view
    analysis_pass = getattr(args, "analysis_pass", None)
    high_resolution_only = bool(getattr(args, "high_resolution_only", False)) or analysis_pass is not None
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if view_name not in entry["views"]:
        raise ValueError(
            f"Dataset '{ds_name}' has no '{view_name}' view in datasets.json."
        )
    expected_view_for_pass = {
        "uncorrected": {"batch_effect_uncorrected"},
        "corrected": {"batch_effect_corrected"},
    }
    if analysis_pass is not None:
        if analysis_pass not in expected_view_for_pass:
            raise ValueError(f"Unknown analysis pass: {analysis_pass}")
        if view_name not in expected_view_for_pass[analysis_pass]:
            raise ValueError(
                f"analysis pass {analysis_pass!r} requires one of "
                f"{sorted(expected_view_for_pass[analysis_pass])!r}"
            )
        if analysis_pass == "corrected" and entry.get("batch_col") is None:
            raise ValueError(
                "corrected batch-effect view requires a confirmed columns.batch"
            )
    view_output = entry["views"][view_name]["output_file"]
    input_path = Path(args.input_dir) / view_output
    if not input_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_path}")

    lowres_col = entry.get("cell_type_low_res")
    highres_col = entry.get("cell_type_high_res")
    technical_batch = entry.get("batch_col") if analysis_pass == "corrected" else None

    def output_name(suffix, n, res_label=None, extension="dists"):
        if analysis_pass is not None:
            return (
                f"{ds_name}_batch_effect_{analysis_pass}_hvg{n}_highres_"
                f"{suffix}_{extension}.feather"
            )
        return f"{ds_name}_hvg{n}{res_label or ''}_{suffix}_{extension}.feather"

    tiers = [("_highres", highres_col)] if high_resolution_only else [
        ("_lowres", lowres_col),
        ("_highres", highres_col),
    ]
    combos = []
    if args.method == "mrvi":
        if lowres_col is None:
            print(f"WARNING: {ds_name}: cell_type_low_res is None; skipping MrVI.")
            return
        for n in args.hvg:
            out_name = (
                output_name("mrvi", n, extension="dists")
                if analysis_pass is not None
                else f"{ds_name}_hvg{n}_mrvi_dists.feather"
            )
            combos.append((n, "_highres" if analysis_pass else "_lowres",
                           None, None, run_mrvi, out_name))

    elif args.method == "scpoli":
        for res_label, ct_col in tiers:
            if ct_col is None:
                continue
            for n in args.hvg:
                for dim in scpoli_dims_for(n, res_label):
                    out_name = (
                        f"{ds_name}_hvg{n}{res_label}_scpoli_dims{dim}_embs.feather"
                    )
                    combos.append((n, res_label, ct_col, dim, run_scpoli, out_name))

    elif args.method in ("pilot", "qot", "pilotgm"):
        suffix = {"pilot": "pilot", "qot": "qot", "pilotgm": "pilotgm"}[args.method]
        for res_label, ct_col in tiers:
            if ct_col is None:
                continue
            for n in args.hvg:
                if run_wass_combo_for(n, res_label):
                    out_name = (
                        output_name(suffix, n, res_label=res_label)
                        if analysis_pass is not None
                        else f"{ds_name}_hvg{n}{res_label}_{suffix}_dists.feather"
                    )
                    combos.append((n, res_label, ct_col, None, None, out_name))

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
        if recorded_feather_valid(out_path) and not args.force:
            print(f"Already processed and validated: {out_name}")
            continue
        pending.append((n, res_label, ct_col, payload, run_fn, out_path))

    if not pending:
        return

    print(f"Loading {input_path} ...")
    adata = sc.read_h5ad(str(input_path), backed="r")
    validate_benchmark_h5ad_contract(adata, args.view, args.method)
    adata = adata.to_memory()
    print(
        f"INPUT_PROFILE bytes={input_path.stat().st_size} "
        f"cells={adata.n_obs} genes={adata.n_vars}",
        flush=True,
    )

    if "Sample" not in adata.obs.columns:
        raise ValueError(
            f"Cannot find standardized sample column 'Sample' in obs of {input_path}."
        )
    if technical_batch is not None and technical_batch not in adata.obs.columns:
        raise ValueError(
            f"Confirmed batch column '{technical_batch}' not found in obs of {input_path}"
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
        # on them): MrVI_hvg{n}, scPoli_hvg{n}_dims{d}{res},
        # PILOT_hvg{n}{res}, QOT_hvg{n}{res}, PILOT-GM-VAE_hvg{n}{res}.
        if args.method == "mrvi":
            method_str = f"MrVI_hvg{n}"
        elif args.method == "scpoli":
            method_str = f"scPoli_hvg{n}_dims{payload}{res_label}"
        elif args.method == "qot":
            method_str = f"QOT_hvg{n}{res_label}"
        elif args.method == "pilotgm":
            method_str = f"PILOT-GM-VAE_hvg{n}{res_label}"
        else:
            method_str = f"PILOT_hvg{n}{res_label}"

        print(f"Processing {method_str} ...")
        start_time = time.time()
        profile_gpu = args.method in ("mrvi", "scpoli") and torch.cuda.is_available()
        if profile_gpu:
            torch.cuda.reset_peak_memory_stats()
        try:
            if args.method == "mrvi":
                run_mrvi(sub, args.device, out_path, batch_key=technical_batch)
            elif args.method == "scpoli":
                run_scpoli(sub, ct_col, payload, out_path)
            elif args.method == "qot":
                run_qot(sub, ct_col, args.view, n, out_path)
            elif args.method == "pilotgm":
                run_pilotgm(sub, ct_col, args.view, n, out_path, ds_name, args.device)
            else:
                run_pilot(sub, ct_col, args.view, n, out_path)
        finally:
            if profile_gpu:
                report_gpu_memory(method_str)
        exec_time = time.time() - start_time

        log_execution_time(ds_name, method_str, exec_time, args.log_file)
        print(f"  -> Saved: {out_path} ({exec_time:.2f}s, "
              f"{peak_rss_gb():.2f} GB peak RSS)")
        gc.collect()


def main():
    parser = argparse.ArgumentParser(
        description="Run Python benchmark methods (MrVI/scPoli/PILOT/QOT/"
                    "PILOT-GM-VAE) on a preprocessed benchmark view h5ad."
    )
    parser.add_argument("--config_path", required=True,
                        help="Path to datasets.json")
    parser.add_argument("--ds_name", required=True,
                        help="Dataset key in datasets.json")
    parser.add_argument("--view", default="benchmark_analysis",
                        help="View name (default: benchmark_analysis)")
    parser.add_argument("--analysis_pass", default=None,
                        choices=["uncorrected", "corrected"],
                        help="Batch-effect pass; requires the matching explicit view")
    parser.add_argument("--high_resolution_only", action="store_true",
                        help="Run only the configured high-resolution tier")
    parser.add_argument("--method", required=True,
                        choices=["mrvi", "scpoli", "pilot", "qot", "pilotgm"],
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
    if args.analysis_pass is not None:
        args.hvg = [2000]

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
