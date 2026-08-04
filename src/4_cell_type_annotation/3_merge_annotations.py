import os
import pandas as pd
import anndata as ad
import glob
import sys
from pathlib import Path


def merge_annotations(h5ad_path: str, annot_dir: str, output_path: str | None = None):
    """
    Merge all chunk annotation feather files into the input .h5ad's obs.

    Parameters
    ----------
    h5ad_path
        Path to the input .h5ad file (same as used by cell type annotation).
    annot_dir
        Directory containing annotations_chunk_*.feather files.
    output_path
        Output .h5ad path. If None, overwrites h5ad_path.
    """
    if output_path is None:
        output_path = h5ad_path

    annot_files = sorted(glob.glob(f"{annot_dir}/annotations_chunk_*.feather"))
    if not annot_files:
        print(f"No annotation feather files found in {annot_dir}")
        return

    print(f"Found {len(annot_files)} annotation chunk files")
    annotations = pd.concat(
        [pd.read_feather(f) for f in annot_files], ignore_index=True
    )
    # Barcodes repeat across samples, and multi-view datasets produce one
    # feather set per view in the same dir; join on (sample, barcode) instead
    # of cell_barcode alone to avoid duplicate-index row explosion.
    sample_col = os.environ.get("SAMPLE_COLNAME", "Sample")
    if "sample" not in annotations.columns:
        print(f"'sample' column missing in annotation feathers of {annot_dir}")
        return
    annotations["_key"] = (
        annotations["sample"].astype(str) + "_" + annotations["cell_barcode"].astype(str)
    )
    n_dup = annotations["_key"].duplicated().sum()
    if n_dup:
        print(f"WARNING: {n_dup} duplicate (sample, barcode) keys across chunk feathers; keeping first.")
        annotations = annotations.drop_duplicates("_key", keep="first")
    annotations = annotations.set_index("_key").drop(columns=["cell_barcode", "sample"])
    print(f"Total annotation entries: {len(annotations)}")

    adata = ad.read_h5ad(h5ad_path)
    print(f"h5ad obs entries before merge: {adata.n_obs}")

    obs = adata.obs.copy()
    if sample_col not in obs.columns:
        raise ValueError(
            f"Sample column '{sample_col}' not found in obs of {h5ad_path}. "
            f"Available columns: {list(obs.columns)}"
        )
    obs["_key"] = obs[sample_col].astype(str) + "_" + adata.obs_names.astype(str)

    original_n_obs = adata.n_obs
    adata.obs = obs.join(annotations, how="left", on="_key").drop(columns="_key")
    merged_count = adata.obs["layer1"].notna().sum()
    print(f"Rows with annotations after merge: {merged_count} / {original_n_obs}")

    # Subset obs to only whitelisted annotation columns
    hitme_cols_keep = ["IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                       "cellCycle.G2M_UCell", "layer1", "layer2", "layer3"]
    scatomic_cols = ["layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                     "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence"]
    annot_cols = set(hitme_cols_keep + scatomic_cols)
    orig_cols = [c for c in adata.obs.columns if c not in annot_cols]
    existing_annot = [c for c in annot_cols if c in adata.obs.columns]
    adata.obs = adata.obs[orig_cols + existing_annot]

    adata.write_h5ad(output_path)
    print(f"Saved annotated .h5ad to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python 3_merge_annotations.py <h5ad_path> [annot_dir] [output_path]")
        sys.exit(1)

    h5ad_path = sys.argv[1]
    annot_dir = sys.argv[2] if len(sys.argv) > 2 else str(Path(h5ad_path).parent)
    output_path = sys.argv[3] if len(sys.argv) > 3 else None
    merge_annotations(h5ad_path, annot_dir, output_path)
