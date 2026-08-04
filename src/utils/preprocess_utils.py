import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import rpy2.robjects as ro
from pathlib import Path


# ---------------------------------------------------------------------------
# R interop — convert .rds to cached raw .h5ad
# ---------------------------------------------------------------------------
ro.r('source("src/utils/load_all_functions.R")')
ro.r('''
convert_rds_to_raw_h5ad <- function(input_path, output_path) {
  seurat <- readRDS(input_path)
  seurat <- create_clean_seuratv5_object(seurat)

  if (!file.exists(output_path)) {
    seurat@assays$RNA@data <- seurat@assays$RNA@counts
    write_h5ad(seurat, output_path)
    seurat@assays$RNA@data <- NULL
  }
}
''')
convert_rds_to_raw_h5ad_r = ro.globalenv["convert_rds_to_raw_h5ad"]


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------
def load_single_input(input_name, input_dir, output_dir):
    input_path = input_dir / input_name
    if str(input_name).endswith(".rds"):
        stem = Path(input_name).stem
        raw_h5ad_path = output_dir / f"{stem}_raw.h5ad"
        convert_rds_to_raw_h5ad_r(str(input_path), str(raw_h5ad_path))
        return sc.read_h5ad(raw_h5ad_path)
    elif str(input_name).endswith(".h5ad"):
        return sc.read_h5ad(input_path)
    else:
        raise ValueError(f"Unsupported file format: {input_name}")


def load_input(input_names, input_dir, output_dir):
    if isinstance(input_names, list):
        adatas = [load_single_input(n, input_dir, output_dir) for n in input_names]
        return sc.concat(adatas, index_unique="_") if len(adatas) > 1 else adatas[0]
    return load_single_input(input_names, input_dir, output_dir)


# ---------------------------------------------------------------------------
# Cell (row) subsetting for views
# ---------------------------------------------------------------------------
def apply_subset_vars(adata, subset_vars):
    if not subset_vars:
        return adata
    mask = pd.Series(True, index=adata.obs_names)
    for col, spec in subset_vars.items():
        if col not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {col}")
        col_mask = adata.obs[col].isin(spec["values"])
        mask &= col_mask if spec.get("op", "in") == "in" else ~col_mask
    return adata[mask].copy()
