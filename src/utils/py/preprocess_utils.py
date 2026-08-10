import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import rpy2.robjects as ro
from pathlib import Path


# ---------------------------------------------------------------------------
# R interop — convert .rds to cached raw .h5ad
# ---------------------------------------------------------------------------
# The R import is pinned to PROJECT_ROOT at import time (setwd) so the
# repo-relative source() below works regardless of caller CWD. All R interop
# calls use absolute paths; this makes callers CWD-independent.
PROJECT_ROOT = Path(__file__).resolve().parents[3]
ro.r('setwd')(str(PROJECT_ROOT))
ro.r('''
# Lighter than load_all_functions.R: preprocessing only converts .rds -> h5ad
# (readRDS -> create_clean_seuratv5_object -> write_h5ad); the benchmark-only
# packages (HiTME, MOFA2, scITD, GloScope, scECODA, ...) must not be pulled in
# (same pattern as the annotation worker 2.1.1_process_chunk.R).
source("src/utils/seurat_utils.R")
library(Seurat)
library(anndataR)
''')
ro.r('''
convert_rds_to_raw_h5ad <- function(input_path, output_path) {
  seurat <- readRDS(input_path)
  seurat <- create_clean_seuratv5_object(seurat)

  if (!file.exists(output_path)) {
    # anndataR records string columns with encoding='ascii' regardless of
    # content; any non-ASCII byte then breaks sc.read_h5ad (UnicodeDecodeError).
    # Lossy replacement is deliberate: clinical/sample metadata is ASCII-safe.
    # Factors carry their values in levels (also character data), so both
    # types are sanitized; factor-ness (incl. ordered) is preserved.
    md <- seurat@meta.data
    md[] <- lapply(md, function(x) {
      if (is.factor(x)) {
        x_conv <- iconv(as.character(x), from = "latin1", to = "ASCII", sub = " ")
        lvls <- iconv(levels(x), from = "latin1", to = "ASCII", sub = " ")
        if (is.ordered(x)) ordered(x_conv, levels = lvls) else factor(x_conv, levels = lvls)
      } else if (is.character(x)) {
        iconv(x, from = "latin1", to = "ASCII", sub = " ")
      } else {
        x
      }
    })
    seurat@meta.data <- md

    # Defensive layer alignment: anndataR validates that layer colnames are
    # identical to the assay features (order- and content-sensitive); some RDSes
    # carry layers misaligned in the middle of the gene list (head/tail
    # identical), which would fail the anndataR strict-identity check. Fresh
    # Assay5 layers have no rownames (features live in @features and are
    # resolved by anndataR), so only layers that actually carry rownames are
    # aligned.
    rna <- seurat[["RNA"]]
    if (inherits(rna, "Assay5") && length(rna@layers) > 0) {
      features <- rownames(rna)
      for (lyr in names(rna@layers)) {
        m <- rna@layers[[lyr]]
        if (!is.null(rownames(m)) && !identical(rownames(m), features)) {
          if (!setequal(rownames(m), features)) {
            stop(sprintf(
              "Layer '%s' rownames are not the same gene set as the assay features (%d missing, %d extra); cannot align.",
              lyr,
              sum(!features %in% rownames(m)),
              sum(!rownames(m) %in% features)
            ))
          }
          n_diff <- sum(rownames(m) != features)
          first_diff <- rownames(m)[which(rownames(m) != features)[seq_len(min(20, n_diff))]]
          message(sprintf(
            "convert_rds_to_raw_h5ad: reindexing layer '%s' rownames to assay features (%d/%d positions differ; first diffs: %s)",
            lyr, n_diff, length(features), paste(first_diff, collapse = ", ")
          ))
          rna@layers[[lyr]] <- m[features, , drop = FALSE]
        }
      }
      seurat[["RNA"]] <- rna
    }

    write_h5ad(seurat, output_path)
  }
}
''')
convert_rds_to_raw_h5ad_r = ro.globalenv["convert_rds_to_raw_h5ad"]


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------
def load_single_input(input_name, input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    input_path = input_dir / input_name
    if str(input_name).endswith(".rds"):
        stem = Path(input_name).stem
        raw_h5ad_path = output_dir / f"{stem}_raw.h5ad"
        convert_rds_to_raw_h5ad_r(str(input_path), str(raw_h5ad_path))
        try:
            adata = sc.read_h5ad(raw_h5ad_path)
        except UnicodeDecodeError:
            # Legacy caches written before meta.data ASCII sanitization carry
            # non-ASCII bytes tagged encoding='ascii' (anndataR); delete the
            # broken cache and regenerate it once with the sanitized converter.
            print(
                f"UnicodeDecodeError reading cached {raw_h5ad_path.name}; "
                "deleting and re-converting."
            )
            raw_h5ad_path.unlink(missing_ok=True)
            convert_rds_to_raw_h5ad_r(str(input_path), str(raw_h5ad_path))
            adata = sc.read_h5ad(raw_h5ad_path)
        # anndataR writes X=None with the raw counts in layers["counts"];
        # promote the counts layer to X (mirror of base_preprocessing).
        if adata.X is None and "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
        return adata
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
