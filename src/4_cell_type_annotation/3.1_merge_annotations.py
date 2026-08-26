import os
import re
import sys
import argparse
import pandas as pd
import anndata as ad
import glob
from pathlib import Path


def _safe_annotation_scalar(value):
    """Convert mixed annotation objects to h5ad-safe nullable strings."""
    if value is None:
        return pd.NA
    try:
        if bool(pd.isna(value)):
            return pd.NA
    except (TypeError, ValueError):
        pass
    return str(value)


def _coerce_annotation_columns(obs, columns):
    """Keep numeric scores numeric; normalize label columns for HDF5 writes."""
    for column in columns:
        series = obs[column]
        if not pd.api.types.is_numeric_dtype(series):
            obs[column] = series.map(_safe_annotation_scalar).astype("string")
    return obs

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

    # Fresh annotation output columns: the ONLY annotation columns this
    # pipeline emits, taken fresh from the annotation feathers. NOT a
    # keep-whitelist — pre-existing versions of these names are always
    # dropped (see drop_cols below); obs may only carry what this annotation
    # run freshly produced (or study metadata). Must mirror
    # 2.1.1_process_chunk.R.
    HITME_OUTPUT_COLS = ["IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                         "cellCycle.G2M_UCell", "layer1", "layer2", "layer3"]
    SCATOMIC_OUTPUT_COLS = ["layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                            "scATOMIC_pred", "S.Score", "G2M.Score", "Phase",
                            "classification_confidence"]
    ANNOT_OUTPUT_COLS = set(HITME_OUTPUT_COLS + SCATOMIC_OUTPUT_COLS)

    # Legacy annotation column classification: pre-existing columns that must
    # never survive into the final obs (the annotation worker wipes them at
    # source; this is the merge-time backstop). Tier 1: unconditional
    # pattern matches. Tier 2: exact-name matches, dropped with a loud
    # warning — exact names only, so author columns like a plain "cellCycle"
    # or "Phase_variant" never collide. Must mirror 2.1.1_process_chunk.R.
    LEGACY_ANNOT_TIER1 = [r"^scGate", r"^functional\.cluster", r"_UCell$",
                          r"^scATOMIC", r"^layer_?\d"]
    LEGACY_ANNOT_TIER2 = ["S.Score", "G2M.Score", "Phase", "classification_confidence",
                          "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"]

    annot_files = sorted(glob.glob(f"{annot_dir}/annotations_chunk_*.feather"))
    if not annot_files:
        print(f"No annotation feather files found in {annot_dir}")
        sys.exit(1)

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
        sys.exit(1)
    annotations["__annot_key"] = (
        annotations[sample_col].astype(str) + "_" + annotations["cell_barcode"].astype(str)
    )
    n_dup = annotations["__annot_key"].duplicated().sum()
    if n_dup:
        print(f"WARNING: {n_dup} duplicate (sample, barcode) keys across chunk feathers; keeping first.")
        annotations = annotations.drop_duplicates("__annot_key", keep="first")
    annotations = annotations.set_index("__annot_key").drop(columns=["cell_barcode", sample_col])
    print(f"Total annotation entries: {len(annotations)}")

    adata = ad.read_h5ad(h5ad_path)
    print(f"h5ad obs entries before merge: {adata.n_obs}")

    obs = adata.obs.copy()
    # Pre-existing annotation columns are dropped unconditionally: the fresh
    # output names (ANNOT_OUTPUT_COLS — old whitelist semantics never apply to
    # pre-existing columns), Tier-1 pattern matches (any scGate_*,
    # functional.cluster*, *_UCell, scATOMIC*, layer<digit>*/layer_<digit>*
    # column), and Tier-2 exact names present in obs.
    tier1_matches = {c for c in obs.columns if any(re.search(p, c) for p in LEGACY_ANNOT_TIER1)}
    tier2_present = set(LEGACY_ANNOT_TIER2) & set(obs.columns)
    drop_cols = set(ANNOT_OUTPUT_COLS) | tier1_matches | tier2_present
    if sample_col not in obs.columns:
        raise ValueError(
            f"Sample column '{sample_col}' not found in obs of {h5ad_path}. "
            f"Available columns: {list(obs.columns)}"
        )

    # Idempotency: if this view was already merged (e.g. a previous merge run
    # merged other views first and failed), drop the existing annotation
    # columns (fresh-output names + legacy tiers) before the join — pandas
    # join would otherwise raise "columns overlap but no suffix specified".
    existing = [c for c in drop_cols if c in obs.columns]
    if existing:
        tier2_flagged = sorted(set(LEGACY_ANNOT_TIER2) & set(existing))
        msg = (f"Dropping {len(existing)} pre-existing annotation column(s) "
               f"before re-merge: {existing}")
        if tier2_flagged:
            msg += (f" — {tier2_flagged} matched Tier-2 exact names "
                    f"(possibly author metadata, not pipeline output)")
        print(msg)
        obs = obs.drop(columns=existing)

    obs["__annot_key"] = obs[sample_col].astype(str) + "_" + adata.obs_names.astype(str)

    original_n_obs = adata.n_obs
    adata.obs = obs.join(annotations, how="left", on="__annot_key").drop(columns="__annot_key")

    # Subset obs to the original (non-annotation) metadata plus only the
    # fresh annotation columns produced by THIS run (all pre-existing
    # annotation columns were dropped above; the join cannot bring back any).
    orig_cols = [c for c in adata.obs.columns if c not in drop_cols]
    existing_annot = [c for c in ANNOT_OUTPUT_COLS if c in adata.obs.columns]
    adata.obs = adata.obs[orig_cols + existing_annot]
    adata.obs = _coerce_annotation_columns(adata.obs, existing_annot)

    if not existing_annot:
        print(
            f"ERROR: no annotation output columns in the merged obs of "
            f"{h5ad_path} — the annotation run produced no results (all "
            f"scATOMIC/HiTME attempts failed in 2.1.1_process_chunk.R, so the "
            f"feathers only hold cell barcodes). Not saving an unannotated h5ad."
        )
        sys.exit(1)

    merged_count = int(adata.obs[existing_annot].notna().any(axis=1).sum())
    print(f"Rows with annotations after merge: {merged_count} / {original_n_obs}")
    if merged_count == 0:
        print(
            f"ERROR: merge matched 0/{original_n_obs} cells. This indicates "
            f"(sample, barcode) key drift between the view h5ad and the "
            f"annotation feathers — check that SAMPLE_COLNAME "
            f"('{os.environ.get('SAMPLE_COLNAME', 'Sample')}') is consistent "
            f"between slurm_config.sh and 2.1.1_process_chunk.R."
        )
        sys.exit(1)
    print(f"Match rate: {merged_count / original_n_obs:.2%} of cells matched")

    # Post-merge invariant (fail loudly): every final obs column matching the
    # legacy tier patterns must have come fresh from THIS run's feathers —
    # the worker wipes pre-existing columns at source, so any tier-matching
    # column that is not in the feathers means stale data leaked in (e.g. a
    # drop pattern was missed).
    feather_cols = set(annotations.columns)
    tier_matching = [c for c in adata.obs.columns
                     if any(re.search(p, c) for p in LEGACY_ANNOT_TIER1)
                     or c in LEGACY_ANNOT_TIER2]
    not_from_feather = [c for c in tier_matching if c not in feather_cols]
    if not_from_feather:
        print("ERROR: post-merge invariant violated — obs column(s) matching "
              "the legacy annotation patterns were NOT produced by this "
              f"annotation run: {not_from_feather}")
        sys.exit(1)

    # Provenance marker: "annotated here" — every merged view records the
    # annotation version that wrote its annotation columns.
    adata.uns["ecoda_annotation_version"] = "1"

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
