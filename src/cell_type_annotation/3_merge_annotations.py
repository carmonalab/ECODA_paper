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
    annotations = annotations.set_index("cell_barcode")
    print(f"Total annotation entries: {len(annotations)}")

    adata = ad.read_h5ad(h5ad_path)
    print(f"h5ad obs entries before merge: {adata.n_obs}")

    original_n_obs = adata.n_obs
    adata.obs = adata.obs.join(annotations, how="left")
    merged_count = adata.obs["layer1"].notna().sum()
    print(f"Rows with annotations after merge: {merged_count} / {original_n_obs}")

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
