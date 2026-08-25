import scipy.sparse as sp
import pandas as pd
import scanpy as sc
import rpy2.robjects as ro
from numbers import Integral
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

    # Defensive layer alignment (generalized): anndataR validates that layer
    # dimnames are identical to the assay features (order- and content-
    # sensitive); some RDSes carry layers misaligned in the middle of the gene
    # list (head/tail identical), which would fail the anndataR strict-identity
    # check. Fresh Assay5 layers have no rownames (features live in @features
    # and are resolved by anndataR). Covers ALL assays (not just "RNA") and
    # cell-major layers whose gene names sit in the COLNAMES (transposed
    # convention, e.g. Wu): those are transposed back to the Seurat
    # (features x cells) convention first. Layers with a genuinely different
    # gene set still fail closed.
    align_assay_layers <- function(a, assay_name) {
      if (!inherits(a, "Assay5") || length(a@layers) == 0) return(a)
      # Canonical feature vector, stripped of any stray 'names' attribute
      # (Wu quirk: named Assay5 features propagate names into layer dimnames,
      # breaking anndataR's identical()-based layer validation). unname() is a
      # no-op for clean objects.
      features <- unname(rownames(a))
      for (lyr in names(a@layers)) {
        m <- a@layers[[lyr]]
        if (!is.null(rownames(m))) {
          # Genes in rownames (Seurat convention)
          if (identical(rownames(m), features)) next
          if (setequal(rownames(m), features)) {
            n_diff <- sum(rownames(m) != features)
            first_diff <- rownames(m)[which(rownames(m) != features)[seq_len(min(20, n_diff))]]
            message(sprintf(
              "convert_rds_to_raw_h5ad: reindexing assay '%s' layer '%s' rownames to assay features (%d/%d positions differ; first diffs: %s)",
              assay_name, lyr, n_diff, length(features), paste(first_diff, collapse = ", ")
            ))
            m <- m[features, , drop = FALSE]
          } else if (!is.null(colnames(m)) && setequal(colnames(m), features)) {
            # Cell-major layer with both dimnames: genes are in the colnames
            message(sprintf(
              "convert_rds_to_raw_h5ad: transposing cell-major assay '%s' layer '%s' (genes in colnames) to Seurat convention",
              assay_name, lyr
            ))
            m <- Matrix::t(m)[features, , drop = FALSE]
          } else {
            stop(sprintf(
              "Layer '%s' of assay '%s' has neither rownames nor colnames equal to the assay features (%d missing, %d extra in rownames); cannot align.",
              lyr, assay_name,
              sum(!features %in% rownames(m)),
              sum(!rownames(m) %in% features)
            ))
          }
        } else if (!is.null(colnames(m)) && setequal(colnames(m), features)) {
          # Cell-major layer without rownames (e.g. Wu): transposed convention
          message(sprintf(
            "convert_rds_to_raw_h5ad: transposing cell-major assay '%s' layer '%s' (genes in colnames) to Seurat convention",
            assay_name, lyr
          ))
          m <- Matrix::t(m)[features, , drop = FALSE]
        } else {
          # Fresh layer without gene dimnames (or barcodes only): anndataR
          # resolves features from @features
          next
        }
        a@layers[[lyr]] <- m
      }
      return(a)
    }
    for (assay_name in names(seurat@assays)) {
      seurat[[assay_name]] <- align_assay_layers(seurat[[assay_name]], assay_name)
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
        # Cache auto-repair (2): a *raw.h5ad cache left behind by a write_h5ad
        # crash in an earlier run reads fine but has neither X nor a counts
        # layer (anndataR always writes counts into layers["counts"], X=None).
        # The R converter only writes when the cache file is missing, so such a
        # stale partial cache blocks regeneration and surfaces as
        # base_preprocessing's "Input has neither X, nor a counts layer, nor
        # raw counts". Delete and re-convert once (same pattern as the
        # UnicodeDecodeError repair above).
        if adata.X is None and "counts" not in adata.layers and adata.raw is None:
            print(
                f"Cached {raw_h5ad_path.name} has no X/counts/raw (partial write "
                "from an earlier run); deleting and re-converting."
            )
            raw_h5ad_path.unlink(missing_ok=True)
            convert_rds_to_raw_h5ad_r(str(input_path), str(raw_h5ad_path))
            adata = sc.read_h5ad(raw_h5ad_path)
            if adata.X is None and "counts" not in adata.layers and adata.raw is None:
                raise ValueError(
                    f"Re-converted {raw_h5ad_path.name} still has no X, counts "
                    "layer, or raw; the source RDS appears to carry no count data."
                )
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
        adata = sc.concat(adatas, index_unique="_") if len(adatas) > 1 else adatas[0]
    else:
        adata = load_single_input(input_names, input_dir, output_dir)

    # anndata's write_h5ad rejects obs.index.name colliding with an obs column
    # whose values differ. sc.concat(index_unique="_") de-duplicates the index
    # with _N suffixes when inputs share barcode values, which breaks an
    # index==column equality that held per-file (e.g. the GongSharma SoundLife
    # files' redundant 'barcodes' column). The index holds the canonical cell
    # barcodes; drop the redundant duplicate column if present.
    if adata.obs.index.name is not None and adata.obs.index.name in adata.obs.columns:
        del adata.obs[adata.obs.index.name]
    return adata


# ---------------------------------------------------------------------------
# Cell (row) subsetting for views
# ---------------------------------------------------------------------------
def apply_subset_vars(adata, subset_vars, copy=True):
    if not subset_vars:
        return adata
    mask = pd.Series(True, index=adata.obs_names)
    for col, spec in subset_vars.items():
        if col not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {col}")
        col_mask = adata.obs[col].isin(spec["values"])
        mask &= col_mask if spec.get("op", "in") == "in" else ~col_mask
    subset = adata[mask]
    return subset.copy() if copy else subset


def remove_low_cellcount_samples(
    adata, sample_col="Sample", min_cells_per_sample=500
):
    """Remove samples with fewer than ``min_cells_per_sample`` observations.

    The returned metadata maps each removed sample ID to its observation
    count. Counts are computed before any per-cell or per-gene preprocessing;
    missing and blank sample IDs are invalid input rather than a removable
    sample.
    """
    if sample_col not in adata.obs.columns:
        raise KeyError(
            f"sample column {sample_col!r} is missing; "
            f"available columns: {list(adata.obs.columns)}"
        )
    if isinstance(min_cells_per_sample, bool) or not isinstance(
        min_cells_per_sample, Integral
    ):
        raise TypeError("min_cells_per_sample must be a positive integer")
    if min_cells_per_sample <= 0:
        raise ValueError("min_cells_per_sample must be a positive integer")

    sample_values = adata.obs[sample_col]
    sample_ids = sample_values.astype("string")
    invalid = sample_values.isna() | sample_ids.str.strip().eq("")
    if bool(invalid.any()):
        invalid_rows = list(adata.obs_names[invalid][:5])
        raise ValueError(
            f"sample column {sample_col!r} contains missing or blank IDs; "
            f"first invalid observations: {invalid_rows}"
        )

    counts = sample_ids.value_counts(sort=False)
    removed = counts[counts < min_cells_per_sample].sort_index()
    keep_ids = counts[counts >= min_cells_per_sample].index
    keep_mask = sample_ids.isin(keep_ids).to_numpy()
    if not bool(keep_mask.any()):
        raise ValueError(
            f"sample-count filtering removed all {adata.n_obs} observations; "
            f"threshold={min_cells_per_sample}"
        )

    filtered = adata[keep_mask].copy()
    removed_counts = {str(sample): int(count) for sample, count in removed.items()}
    return filtered, removed_counts
