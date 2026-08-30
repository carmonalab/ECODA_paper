"""
1.5.1_reconstruct_myocardial_counts.py — Reconstruct raw integer count matrix for Myocardial Infarction.

The staged Myocardial Infarction dataset (Kuppe et al. 2022 Nature, PMID 35948637;
Myocardial_Infarc_2.h5ad) only provides log1p-normalized expression in .X and lacks
an explicit raw count layer or raw.X.

This script reconstructs the exact raw UMI count matrix via cell-wise minimum positive
step inversion:
    s_i = min_{j: X_{ij} > 0}(expm1(X_{ij}))
    c_{ij} = round(expm1(X_{ij}) / s_i)

The reconstructed raw integer counts are vaulted to adata.layers["counts"] and set as
adata.X (CSR format) in place (atomic write via temp file + os.replace).

Usage (HPC, via 1.5_submit_myocardial.sh):
    ${PYTHON_BIN} 1.5.1_reconstruct_myocardial_counts.py [--config_path ...] [--data_dir ...]
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.datasets_io import read_datasets_json

DS_NAME = "Myocardial_infarction"
DEFAULT_FILE_NAME = "Myocardial_Infarc_2.h5ad"


def reconstruct_raw_counts_matrix(X: sp.spmatrix | np.ndarray) -> sp.csr_matrix:
    """Reconstruct exact raw count matrix from log1p-normalized expression."""
    if not sp.isspmatrix_csr(X):
        X = sp.csr_matrix(X)

    exp_data = np.expm1(X.data)
    indptr = X.indptr
    n_obs = X.shape[0]

    recovered_data = np.empty_like(exp_data, dtype=np.float32)

    for i in range(n_obs):
        start, end = indptr[i], indptr[i + 1]
        if start == end:
            continue
        cell_vals = exp_data[start:end]
        pos_vals = cell_vals[cell_vals > 1e-7]
        if len(pos_vals) == 0:
            recovered_data[start:end] = 0.0
            continue
        step = np.min(pos_vals)
        recovered_data[start:end] = np.round(cell_vals / step)

    return sp.csr_matrix((recovered_data, X.indices, X.indptr), shape=X.shape)

def _valid_counts_layer(counts, shape: tuple[int, int]) -> bool:
    """Return whether an existing raw-count layer satisfies the full contract."""
    if getattr(counts, "shape", None) != shape:
        return False
    values = counts.data if sp.issparse(counts) else np.asarray(counts).reshape(-1)
    try:
        return bool(
            values.size > 0
            and np.isfinite(values).all()
            and (values >= 0).all()
            and np.equal(values, np.round(values)).all()
        )
    except (TypeError, ValueError):
        return False


def process_myocardial_file(
    input_path: Path, output_path: Path, force: bool = False
) -> bool:
    print(f"=== Processing Myocardial Infarction Counts: {input_path.name} ===", flush=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    t0 = time.perf_counter()
    adata = ad.read_h5ad(str(input_path))
    print(f"Loaded {input_path.name}: shape {adata.shape}, layers: {list(adata.layers.keys())}", flush=True)

    # A valid counts layer is already the authoritative raw-count artifact.
    # --force invalidates the sidecar and forces revalidation, but must not
    # replace valid counts with an overflow-prone inversion of X.
    if "counts" in adata.layers and _valid_counts_layer(
        adata.layers["counts"], adata.shape
    ):
        print(
            "Layer 'counts' already satisfies the finite nonnegative integer "
            "contract; retaining it after full validation.",
            flush=True,
        )
        return True

    if adata.X is None:
        raise ValueError("Cannot reconstruct counts: adata.X is None.")

    print("Reconstructing raw count matrix via cell-wise minimum step inversion...", flush=True)
    counts_csr = reconstruct_raw_counts_matrix(adata.X)

    # Verification of reconstructed counts
    rec_data = counts_csr.data
    is_non_neg = bool(np.all(rec_data >= 0))
    is_int = bool(np.all(np.abs(rec_data - np.round(rec_data)) < 1e-3))
    print(
        f"Reconstruction sanity: non_negative={is_non_neg}, integer_values={is_int}, "
        f"min={rec_data.min():.1f}, max={rec_data.max():.1f}, nnz={counts_csr.nnz}",
        flush=True,
    )

    if not is_non_neg or not is_int:
        raise ValueError("Count reconstruction sanity check failed!")

    adata.layers["counts"] = counts_csr
    adata.X = counts_csr.copy()

    # Atomic write to temp file and replace
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(f"{output_path.name}.tmp_{os.getpid()}")
    try:
        print(f"Writing updated h5ad to temporary file: {tmp_path.name}...", flush=True)
        adata.write_h5ad(str(tmp_path))
        os.replace(str(tmp_path), str(output_path))
        print(f"Atomic update complete -> {output_path} (took {time.perf_counter() - t0:.2f}s)", flush=True)
    finally:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct raw integer counts for Myocardial Infarction"
    )
    parser.add_argument(
        "--config_path",
        default=os.environ.get("DATASETS_JSON_FILE"),
        help="Path to datasets.json (defaults to $DATASETS_JSON_FILE)",
    )
    parser.add_argument(
        "--data_dir",
        default=None,
        help="Path to dataset data directory (defaults to ${HPC_SCRATCH_DIR}/Myocardial_infarction/data)",
    )
    parser.add_argument(
        "--input_file",
        default=None,
        help="Explicit input file path (overrides data_dir resolution)",
    )
    parser.add_argument(
        "--output_file",
        default=None,
        help="Explicit output file path (defaults to overwriting input_file in place)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Revalidate existing counts; recompute only when invalid or missing",
    )
    args = parser.parse_args()

    if args.input_file:
        input_path = Path(args.input_file)
        output_path = Path(args.output_file) if args.output_file else input_path
    else:
        scratch_dir = os.environ.get("HPC_SCRATCH_DIR")
        if args.data_dir:
            data_dir = Path(args.data_dir)
        elif scratch_dir:
            data_dir = Path(scratch_dir) / DS_NAME / "data"
        else:
            raise ValueError(
                "Neither --input_file, --data_dir, nor HPC_SCRATCH_DIR was provided."
            )

        fname = DEFAULT_FILE_NAME
        if args.config_path and Path(args.config_path).exists():
            cfg = read_datasets_json(args.config_path)
            if DS_NAME in cfg:
                fnames = cfg[DS_NAME].get("file_names")
                if isinstance(fnames, list) and len(fnames) > 0:
                    fname = fnames[0]
                elif isinstance(fnames, str):
                    fname = fnames

        input_path = data_dir / fname
        output_path = input_path

    process_myocardial_file(input_path, output_path, force=args.force)


if __name__ == "__main__":
    main()
