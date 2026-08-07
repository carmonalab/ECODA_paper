"""
Create the Combined PBMC dataset (Stephenson + GongSharma + Zhu) for batch effect analysis.

Usage:
    Local (default):      python 1.1.1_create_combinedpbmc_dataset.py
    HPC (per-dataset):    python 1.1.1_create_combinedpbmc_dataset.py --layout per-dataset
                          [--base-path ${HPC_SCRATCH_DIR} --output-dir ${HPC_SCRATCH_DIR}/CombinedPBMC/data]

Output:
    Local (flat):         PROJECT_ROOT/data/combined_pbmc_batch_effect_analysis.h5ad
    HPC (per-dataset):    ${HPC_SCRATCH_DIR}/CombinedPBMC/data/combined_pbmc_batch_effect_analysis.h5ad

Defaults are driven by the environment: if HPC_SCRATCH_DIR is set, base-path
defaults to ${HPC_SCRATCH_DIR} with layout=per-dataset (per-source raw inputs at
${HPC_SCRATCH_DIR}/<ds>/data, e.g. Stephenson/data, Gongsharma_cmv_young_males/data,
Zhu/data); otherwise the local flat layout (PROJECT_ROOT/data) is used.
output_dir doubles as the rds->h5ad conversion cache ({stem}_raw.h5ad).

Parallel workers:
- The three sources are loaded concurrently in 3 fork worker processes
  (ProcessPoolExecutor, max_workers=3); wall time is ~max(source) instead of the
  sum. Each worker writes a trimmed per-source intermediate to
  output_dir/_intermediates/ (<source>_subset.h5ad, overwritten on rerun); the
  main process then loads the small intermediates for the shared final steps
  (common-gene intersection with <5000 union fallback, Sample prefixing,
  concat, write).
- Stephenson and Zhu load via the R interop (rds -> cached {stem}_raw.h5ad).
  rpy2 is imported lazily INSIDE the worker functions only, so the parent
  process never initializes R before forking.
- GongSharma is read in backed mode (sc.read_h5ad(..., backed="r")): non-counts
  layers are dropped on disk and only the HDF5 chunks covering the ~15 picked
  samples are materialized via to_memory() (identical rng(123) pick logic as
  before); on any failure it falls back to the full in-memory load. The raw
  GongSharma sample column (specimen.specimenGuid, datasets.json columns.sample)
  is standardized to "Sample" by finalize_and_write_intermediate, consistent
  with src/3_scrnaseq_preprocessing/1.1.1_preprocess.py.

HPC notes:
- Must run AFTER src/1_stage_data/1_stage_data.sh (which stages raw inputs per
  dataset to ${HPC_SCRATCH_DIR}/<ds>/data) and BEFORE
  src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh (staging skips datasets
  with folder_name: null, and the preprocess array task reads the combined file
  from ${HPC_SCRATCH_DIR}/CombinedPBMC/data).
- Requires `module load GCCcore/12.2.0` for the R interop (rds->h5ad conversion).
  CWD-independent: preprocess_utils.py pins the embedded R working directory to
  ${PROJECT_ROOT} at import time.
- Heavy loads (GongSharma is huge) may warrant running via a single sbatch job
  instead of interactively on the login node if OOM occurs.
"""

import os
import sys
import gc
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import scanpy as sc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.gene_utils import standardize_gene_symbols
from src.datasets_io import read_datasets_json

KEEP_BASE = ["Sample", "batch", "cond"]
GONGSHARMA_N_SAMPLES = 15
SAMPLE_RNG_SEED = 123


def _resolve_source_paths(ds_name, entry, base_path, layout):
    if layout == "per-dataset":
        base_path = base_path / ds_name / "data"

    input_names = entry.get("file_names")
    if not input_names:
        raise ValueError(f"No file_names for {ds_name}")

    if layout == "per-dataset":
        input_names_list = input_names if isinstance(input_names, list) else [input_names]
        missing = [n for n in input_names_list if not (base_path / n).exists()]
        if missing:
            raise FileNotFoundError(
                f"Raw input(s) not staged for {ds_name}: {missing}. "
                f"Run src/1_stage_data/1_stage_data.sh first (expected in {base_path})."
            )

    return input_names, base_path


def load_and_prepare_source(ds_name, entry, base_path, output_dir, view_name=None, layout="flat"):
    # Lazy import: preprocess_utils.py initializes rpy2 (embedded R session) at
    # module import time. Only worker processes (post-fork) may trigger it.
    from src.utils.preprocess_utils import load_input, apply_subset_vars

    input_names, src_path = _resolve_source_paths(ds_name, entry, base_path, layout)
    adata = load_input(input_names, src_path, output_dir)

    subset_vars = {}
    if view_name:
        subset_vars = entry["views"].get(view_name, {}).get("subset_vars", {})
    if subset_vars:
        adata = apply_subset_vars(adata, subset_vars)

    standardize_gene_symbols(adata)
    adata.var_names_make_unique()

    return adata


def keep_only_cols(adata, cols):
    existing = [c for c in cols if c in adata.obs.columns]
    adata.obs = adata.obs[existing]
    return adata


def finalize_and_write_intermediate(adata, sample_col, batch_label, cond_value, intermediate_path):
    adata.obs["batch"] = batch_label
    adata.obs["cond"] = cond_value
    if sample_col != "Sample" and sample_col in adata.obs.columns:
        adata.obs["Sample"] = adata.obs[sample_col].values
    keep_only_cols(adata, KEEP_BASE)
    print(f"  -> {adata.n_obs} cells", flush=True)
    adata.write_h5ad(str(intermediate_path))
    return adata.n_obs, str(intermediate_path)


def worker_stephenson(entry, base_path, output_dir, view_name, layout, intermediate_path):
    # Subset to the benchmark_analysis view (Site = Ncl only), consistent
    # with previous code.
    print("Loading Stephenson...", flush=True)
    adata = load_and_prepare_source(
        "Stephenson", entry, base_path, output_dir,
        view_name=view_name, layout=layout,
    )
    cond = adata.obs[entry["label_col"]].values
    return finalize_and_write_intermediate(
        adata, entry["sample_col"], "Stephenson", cond, intermediate_path
    )


def worker_gongsharma(entry, base_path, output_dir, layout, intermediate_path):
    print("Loading GongSharma...", flush=True)
    sample_col = entry["sample_col"]
    input_names, src_path = _resolve_source_paths(
        "Gongsharma_cmv_young_males", entry, base_path, layout
    )
    names = input_names if isinstance(input_names, list) else [input_names]
    try:
        backed = []
        unique_samples = []
        for name in names:
            adata = sc.read_h5ad(str(src_path / name), backed="r")
            backed.append(adata)
            for s in adata.obs[sample_col].unique().tolist():
                if s not in unique_samples:
                    unique_samples.append(s)
        rng = np.random.default_rng(SAMPLE_RNG_SEED)
        chosen = rng.choice(unique_samples, size=min(GONGSHARMA_N_SAMPLES, len(unique_samples)), replace=False)
        print(f"  GongSharma: {len(unique_samples)} total samples, picking {len(chosen)}", flush=True)
        subsets = []
        for adata in backed:
            # Backed layer deletion drops the HDF5 dataset without loading X;
            # keeps heavy scaled layers from being materialized.
            for layer_name in list(adata.layers.keys()):
                if layer_name != "counts":
                    del adata.layers[layer_name]
            mask = adata.obs[sample_col].isin(chosen)
            subsets.append(adata[mask].to_memory())
            try:
                adata.file.close()
            except Exception:
                pass
        adata = sc.concat(subsets, index_unique="_") if len(subsets) > 1 else subsets[0]
    except Exception as exc:
        print(
            f"  Warning: backed-mode GongSharma read failed ({exc!r}); "
            "falling back to full in-memory load.",
            flush=True,
        )
        adata = load_and_prepare_source(
            "Gongsharma_cmv_young_males", entry, base_path, output_dir, layout=layout
        )
        unique_samples = adata.obs[sample_col].unique().tolist()
        rng = np.random.default_rng(SAMPLE_RNG_SEED)
        chosen = rng.choice(unique_samples, size=min(GONGSHARMA_N_SAMPLES, len(unique_samples)), replace=False)
        print(f"  GongSharma: {len(unique_samples)} total samples, picking {len(chosen)}", flush=True)
        adata = adata[adata.obs[sample_col].isin(chosen)].copy()

    standardize_gene_symbols(adata)
    adata.var_names_make_unique()
    return finalize_and_write_intermediate(adata, sample_col, "GongSharma", "Healthy", intermediate_path)


def worker_zhu(entry, base_path, output_dir, layout, intermediate_path):
    print("Loading Zhu...", flush=True)
    adata = load_and_prepare_source("Zhu", entry, base_path, output_dir, layout=layout)
    return finalize_and_write_intermediate(adata, entry["sample_col"], "Zhu", "Healthy", intermediate_path)


def main():
    project_root = Path(os.environ.get("PROJECT_ROOT") or Path(__file__).resolve().parents[2])
    hpc_scratch_dir = os.environ.get("HPC_SCRATCH_DIR", "")

    parser = argparse.ArgumentParser(description="Create the Combined PBMC dataset")
    parser.add_argument("--config-path", default="datasets.json",
                        help="Path to datasets.json (relative to PROJECT_ROOT unless absolute)")
    parser.add_argument("--base-path", default=None,
                        help="Root of raw source data. HPC default: ${HPC_SCRATCH_DIR}; "
                             "local default: PROJECT_ROOT/data.")
    parser.add_argument("--output-dir", default=None,
                        help="Output dir for combined .h5ad + rds->h5ad cache. "
                             "HPC default: ${HPC_SCRATCH_DIR}/CombinedPBMC/data; "
                             "local default: PROJECT_ROOT/data.")
    parser.add_argument("--layout", choices=["per-dataset", "flat"], default=None,
                        help="per-dataset: sources at <base-path>/<ds>/data. "
                             "Default: per-dataset if HPC_SCRATCH_DIR set, else flat.")
    args = parser.parse_args()

    if args.base_path:
        base_path = Path(args.base_path)
    elif hpc_scratch_dir:
        base_path = Path(hpc_scratch_dir)
    else:
        base_path = project_root / "data"

    if args.output_dir:
        output_dir = Path(args.output_dir)
    elif hpc_scratch_dir:
        output_dir = Path(hpc_scratch_dir) / "CombinedPBMC" / "data"
    else:
        output_dir = project_root / "data"

    layout = args.layout or ("per-dataset" if hpc_scratch_dir else "flat")

    config_path = Path(args.config_path)
    if not config_path.is_absolute():
        config_path = project_root / config_path

    config = read_datasets_json(config_path)

    source_names = ["Stephenson", "Gongsharma_cmv_young_males", "Zhu"]
    missing = [n for n in source_names if n not in config]
    if missing:
        raise KeyError(
            f"Missing sources in datasets.json (required for CombinedPBMC): {missing}"
        )

    # ---- Load the three sources concurrently (3 fork workers) ----
    # Each worker writes its own small intermediate into output_dir/_intermediates/.
    # The parent never imports src.utils.preprocess_utils, so rpy2/R init happens
    # only inside the R workers, after forking.
    intermediate_dir = output_dir / "_intermediates"
    os.makedirs(intermediate_dir, exist_ok=True)
    mp_context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=4, mp_context=mp_context) as executor:
        futures = {
            "Stephenson": executor.submit(
                worker_stephenson, config["Stephenson"], base_path, output_dir,
                view_name="benchmark_analysis", layout=layout,
                intermediate_path=intermediate_dir / "Stephenson_subset.h5ad",
            ),
            "GongSharma": executor.submit(
                worker_gongsharma, config["Gongsharma_cmv_young_males"], base_path,
                output_dir, layout=layout,
                intermediate_path=intermediate_dir / "GongSharma_subset.h5ad",
            ),
            "Zhu": executor.submit(
                worker_zhu, config["Zhu"], base_path, output_dir, layout=layout,
                intermediate_path=intermediate_dir / "Zhu_subset.h5ad",
            ),
        }
        results = {name: future.result() for name, future in futures.items()}

    path_s = results["Stephenson"][1]
    path_g = results["GongSharma"][1]
    path_z = results["Zhu"][1]

    # ---- Reload the small per-source intermediates ----
    print("Loading per-source intermediates...")
    adata_s = sc.read_h5ad(path_s)
    adata_g = sc.read_h5ad(path_g)
    adata_z = sc.read_h5ad(path_z)

    # ---- Compute common genes from standardized var_names ----
    gene_sets = [set(adata.var_names) for adata in [adata_s, adata_g, adata_z]]
    common_genes = sorted(set.intersection(*gene_sets))
    print(f"Common genes across all three sources: {len(common_genes)}")

    if len(common_genes) < 5000:
        print("Warning: fewer than 5000 common genes. Using union instead.")
        common_genes = sorted(set.union(*gene_sets))
        print(f"Union gene count: {len(common_genes)}")

    # Subset each to common genes in-place (no .copy() to avoid doubling memory)
    adata_s = adata_s[:, adata_s.var_names.isin(common_genes)]
    adata_g = adata_g[:, adata_g.var_names.isin(common_genes)]
    adata_z = adata_z[:, adata_z.var_names.isin(common_genes)]

    # ---- Prefix sample IDs to avoid collisions ----
    adata_s.obs["Sample"] = "Stephenson_" + adata_s.obs["Sample"].astype(str)
    adata_g.obs["Sample"] = "GongSharma_" + adata_g.obs["Sample"].astype(str)
    adata_z.obs["Sample"] = "Zhu_" + adata_z.obs["Sample"].astype(str)

    # ---- Concat and release sources ----
    print("Concatenating...")
    adata_combined = sc.concat(
        [adata_s, adata_g, adata_z],
        index_unique=None,
        join="outer",
    )
    print(f"Combined shape: {adata_combined.shape}")

    del adata_s, adata_g, adata_z
    gc.collect()

    output_path = output_dir / "combined_pbmc_batch_effect_analysis.h5ad"
    os.makedirs(output_dir, exist_ok=True)
    adata_combined.write_h5ad(str(output_path))
    print(f"Saved combined PBMC dataset to: {output_path}")


if __name__ == "__main__":
    main()
