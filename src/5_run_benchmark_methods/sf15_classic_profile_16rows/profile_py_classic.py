#!/usr/bin/env python3
"""Scoped 16-row classic full-H5AD profiling worker for Python benchmark methods.

Reproduces the pre-streaming (fdf137d^) full-H5AD-in-memory loading path
(`sc.read_h5ad(..., backed='r')` + `adata.to_memory()`) and records method
wall-clock runtime (`time_secs`) and peak RSS (`ru_maxrss / 1024^2` -> `mem_GB`).
"""

import argparse
import csv
import importlib.util
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

import scanpy as sc
import scvi
import torch

from src.utils.py.benchmark_h5ad_contract import validate_benchmark_h5ad_contract
from src.utils.py.datasets_io import read_datasets_json


def _load_py_benchmark_module():
    mod_path = (
        REPO_ROOT
        / "src"
        / "5_run_benchmark_methods"
        / "run_python_sample_embedding_methods"
        / "1.1.1_benchmark_methods_py.py"
    )
    if str(mod_path.parent) not in sys.path:
        sys.path.insert(0, str(mod_path.parent))
    spec = importlib.util.spec_from_file_location("ecoda_bench_py", mod_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load benchmark module from {mod_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ds_name", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    args = parser.parse_args()

    method_map = {
        "PILOT_hvg2000": ("pilot", 2000, None),
        "QOT_hvg2000": ("qot", 2000, None),
        "PILOT-GM-VAE_hvg2000": ("pilotgm", 2000, None),
        "MrVI_hvg2000": ("mrvi", 2000, None),
        "scPoli_hvg2000_dims15_highres": ("scpoli", 2000, 15),
    }
    if args.method not in method_map:
        raise ValueError(f"Unsupported Python profiling method: {args.method}")
    py_method, n_hvg, scpoli_dim = method_map[args.method]

    bench_py = _load_py_benchmark_module()
    if py_method in bench_py.GPU_BACKED_METHODS:
        bench_py.validate_gpu_execution(
            py_method, args.device, combo=f"hvg{n_hvg}", view="benchmark_analysis"
        )
    if py_method == "pilotgm" and args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("PILOT-GM-VAE requested CUDA but torch.cuda.is_available() is False")

    scvi.settings.seed = 0

    scratch_dir = Path(os.environ["HPC_SCRATCH_DIR"])
    config = read_datasets_json(str(REPO_ROOT / "datasets.json"), view="benchmark_analysis")
    if args.ds_name not in config:
        raise KeyError(f"Dataset not found in datasets.json: {args.ds_name}")
    entry = config[args.ds_name]
    input_path = scratch_dir / args.ds_name / "output" / entry["views"]["benchmark_analysis"]["output_file"]
    if not input_path.is_file():
        raise FileNotFoundError(f"Input H5AD not found: {input_path}")

    highres_col = entry.get("cell_type_high_res")
    if py_method != "mrvi" and not highres_col:
        raise ValueError(f"{args.ds_name} is missing cell_type_high_res")

    # Classic full-H5AD-in-memory loading (identical to pre-Nov-9 fdf137d^)
    t_load0 = time.time()
    adata = sc.read_h5ad(str(input_path), backed="r")
    validate_benchmark_h5ad_contract(adata, "benchmark_analysis", py_method)
    adata = adata.to_memory()
    load_secs = time.time() - t_load0

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_feather = work_dir / f"{args.ds_name}_{args.method}.feather"

    if py_method in ("mrvi", "scpoli"):
        genes = bench_py.top_n_hvg_genes(adata, n_hvg)
        sub = adata[:, genes].copy()
        sub = bench_py.use_counts_layer(sub, py_method, args.ds_name)
    else:
        sub = adata

    start_time = time.time()
    if args.device == "cuda" and torch.cuda.is_available():
        torch.cuda.reset_peak_memory_stats()

    if py_method == "mrvi":
        bench_py.run_mrvi(sub, args.device, out_feather)
    elif py_method == "scpoli":
        bench_py.run_scpoli(sub, highres_col, scpoli_dim, out_feather)
    elif py_method == "pilot":
        bench_py.run_pilot(sub, highres_col, "benchmark_analysis", n_hvg, out_feather)
    elif py_method == "qot":
        bench_py.run_qot(sub, highres_col, "benchmark_analysis", n_hvg, out_feather)
    elif py_method == "pilotgm":
        bench_py.run_pilotgm(
            sub, highres_col, "benchmark_analysis", n_hvg, out_feather, args.ds_name, args.device
        )

    if args.device == "cuda" and torch.cuda.is_available():
        bench_py.report_gpu_memory(args.method)

    time_secs = time.time() - start_time
    mem_gb = bench_py.peak_rss_gb()
    total_wall_secs = load_secs + time_secs

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp_csv = out_csv.with_name(f".{out_csv.name}.tmp.{os.getpid()}")
    with open(tmp_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "dataset",
                "method",
                "lang",
                "resource_mode",
                "time_secs",
                "load_secs",
                "total_wall_secs",
                "mem_GB",
                "slurm_job_id",
                "slurm_task_id",
                "node",
                "cpu_model",
                "gpu_model",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "dataset": args.ds_name,
                "method": args.method,
                "lang": "py",
                "resource_mode": "gpu" if args.device == "cuda" else "cpu",
                "time_secs": f"{time_secs:.6f}",
                "load_secs": f"{load_secs:.6f}",
                "total_wall_secs": f"{total_wall_secs:.6f}",
                "mem_GB": f"{mem_gb:.6f}",
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "slurm_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
                "node": os.environ.get("ECODA_PROF_NODE", os.uname().nodename),
                "cpu_model": os.environ.get("ECODA_PROF_CPU_MODEL", ""),
                "gpu_model": os.environ.get("ECODA_PROF_GPU_MODEL", "none"),
            }
        )
    os.replace(tmp_csv, out_csv)
    print(
        f"PROFILE_PY_DONE dataset={args.ds_name} method={args.method} "
        f"time_secs={time_secs:.3f} load_secs={load_secs:.3f} mem_GB={mem_gb:.4f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
