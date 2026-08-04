import os
import argparse
import pandas as pd
import anndata as ad
import glob
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

    # Whitelisted annotation columns (must mirror 2.1.1_process_chunk.R)
    hitme_cols_keep = ["IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                       "cellCycle.G2M_UCell", "layer1", "layer2", "layer3"]
    scatomic_cols = ["layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                     "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence"]
    annot_cols = set(hitme_cols_keep + scatomic_cols)

    annot_files = sorted(glob.glob(f"{annot_dir}/annotations_chunk_*.feather"))
    if not annot_files:
        print(f"No annotation feather files found in {annot_dir}")
        return

    print(f"Found {len(annot_files)} annotation chunk files")
    annotations = pd.concat(
        [pd.read_feather(f) for f in annot_files], ignore_index=True
    )
    # Barcodes repeat across samples, and (pre-union) multi-view datasets could
    # produce one feather set per view in the same dir; join on (sample,
    # barcode) instead of cell_barcode alone to avoid duplicate-index row
    # explosion.
    sample_col = os.environ.get("SAMPLE_COLNAME", "Sample")
    if sample_col not in annotations.columns:
        print(f"'{sample_col}' column missing in annotation feathers of {annot_dir}")
        return
    annotations["_key"] = (
        annotations[sample_col].astype(str) + "_" + annotations["cell_barcode"].astype(str)
    )
    n_dup = annotations["_key"].duplicated().sum()
    if n_dup:
        print(f"WARNING: {n_dup} duplicate (sample, barcode) keys across chunk feathers; keeping first.")
        annotations = annotations.drop_duplicates("_key", keep="first")
    annotations = annotations.set_index("_key").drop(columns=["cell_barcode", sample_col])
    print(f"Total annotation entries: {len(annotations)}")

    adata = ad.read_h5ad(h5ad_path)
    print(f"h5ad obs entries before merge: {adata.n_obs}")

    obs = adata.obs.copy()
    if sample_col not in obs.columns:
        raise ValueError(
            f"Sample column '{sample_col}' not found in obs of {h5ad_path}. "
            f"Available columns: {list(obs.columns)}"
        )

    # Idempotency: if this view was already merged (e.g. a previous merge run
    # merged other views first and failed), drop the existing whitelisted
    # annotation columns before the join — pandas join would otherwise raise
    # "columns overlap but no suffix specified".
    existing = [c for c in annot_cols if c in obs.columns]
    if existing:
        print(f"Dropping {len(existing)} existing annotation columns before re-merge: {existing}")
        obs = obs.drop(columns=existing)

    obs["_key"] = obs[sample_col].astype(str) + "_" + adata.obs_names.astype(str)

    original_n_obs = adata.n_obs
    adata.obs = obs.join(annotations, how="left", on="_key").drop(columns="_key")

    # Subset obs to only whitelisted annotation columns
    orig_cols = [c for c in adata.obs.columns if c not in annot_cols]
    existing_annot = [c for c in annot_cols if c in adata.obs.columns]
    adata.obs = adata.obs[orig_cols + existing_annot]

    merged_count = int(adata.obs[existing_annot].notna().any(axis=1).sum()) if existing_annot else 0
    print(f"Rows with annotations after merge: {merged_count} / {original_n_obs}")

    # Atomic write: write to a temp file in the same directory, then rename —
    # an interrupted run must not truncate the only scratch copy of the view.
    tmp_path = str(Path(output_path).with_name(Path(output_path).name + ".tmp"))
    adata.write_h5ad(tmp_path)
    os.replace(tmp_path, output_path)
    print(f"Saved annotated .h5ad to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge per-chunk annotation feathers into a view .h5ad's obs."
    )
    parser.add_argument("--h5ad-path", required=True,
                        help="Path to the view .h5ad file to annotate.")
    parser.add_argument("--annot-dir", default=None,
                        help="Directory containing annotations_chunk_*.feather "
                             "(default: the .h5ad's parent directory).")
    parser.add_argument("--output-path", default=None,
                        help="Output .h5ad path (default: overwrite --h5ad-path).")
    args = parser.parse_args()

    if args.annot_dir is None:
        args.annot_dir = str(Path(args.h5ad_path).parent)
    merge_annotations(args.h5ad_path, args.annot_dir, args.output_path)
