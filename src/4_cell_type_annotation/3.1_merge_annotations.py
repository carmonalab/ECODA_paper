"""Merge run-owned dual-method annotation feathers into a preprocessed h5ad."""

from __future__ import annotations

import argparse
import glob
import hashlib
import os
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from src.utils.py.annotation_contract import (  # noqa: E402
    KNOWN_ANNOTATION_COLUMNS,
    OPTIONAL_HITME_COLUMNS,
    REQUIRED_ANNOTATION_COLUMNS,
    _normalize_expected_keys,
    _validate_frame,
    h5ad_keys,
    validate_feather,
    validate_h5ad,
    validate_sidecar,
)

HITME_OUTPUT_COLS = ["IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                     "cellCycle.G2M_UCell", "layer1", "layer2", "layer3"]
SCATOMIC_OUTPUT_COLS = ["layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                        "scATOMIC_pred", "S.Score", "G2M.Score", "Phase",
                        "classification_confidence"]
ANNOT_OUTPUT_COLS = HITME_OUTPUT_COLS + SCATOMIC_OUTPUT_COLS
NUMERIC_ANNOTATION_COLS = {
    "IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell", "cellCycle.G2M_UCell",
    "S.Score", "G2M.Score",
}
LEGACY_ANNOT_TIER1 = [r"^scGate", r"^functional\.cluster", r"_UCell$",
                      r"^scATOMIC", r"^layer_?\d"]
LEGACY_ANNOT_TIER2 = ["S.Score", "G2M.Score", "Phase", "classification_confidence",
                      "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"]
CONFIDENCE_LABELS = {"confident", "low_confidence"}


def _safe_annotation_scalar(value):
    """Convert mixed annotation objects to h5ad-safe nullable strings."""
    if value is None:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def _coerce_annotation_columns(obs: pd.DataFrame, columns) -> pd.DataFrame:
    """Keep scores numeric; preserve scATOMIC confidence labels or scores."""
    for column in columns:
        if column not in obs:
            continue
        if column == "classification_confidence":
            values = obs[column]
            converted = pd.to_numeric(values, errors="coerce")
            labels = values.astype(object).map(
                lambda value: str(value).strip() if _safe_annotation_scalar(value) is not None else ""
            )
            invalid = values.notna() & converted.isna() & ~labels.isin(CONFIDENCE_LABELS)
            if bool(invalid.any()):
                raise ValueError(
                    "classification_confidence contains values outside the "
                    "numeric/confidence-label contract"
                )
            if bool(labels.isin(CONFIDENCE_LABELS).any()):
                obs[column] = values.map(_safe_annotation_scalar).astype("string")
            else:
                obs[column] = converted.astype("float64")
        elif column in NUMERIC_ANNOTATION_COLS:
            obs[column] = pd.to_numeric(obs[column], errors="coerce").astype("float64")
        else:
            obs[column] = obs[column].map(_safe_annotation_scalar).astype("string")
    return obs


def _keys_from_frame(frame: pd.DataFrame) -> list[tuple[str, str]]:
    if "Sample" not in frame or "cell_barcode" not in frame:
        raise ValueError("annotation frame lacks Sample/cell_barcode keys")
    samples = frame["Sample"].map(_safe_annotation_scalar)
    barcodes = frame["cell_barcode"].map(_safe_annotation_scalar)
    if samples.isna().any() or barcodes.isna().any():
        raise ValueError("annotation frame contains blank keys")
    return list(zip(samples.astype(str), barcodes.astype(str)))


def _write_sidecar(path: Path) -> None:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    sidecar = Path(f"{path}.md5")
    temporary = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    temporary.write_text(
        f"MD5={digest.hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n"
    )
    os.replace(temporary, sidecar)


def _validate_annotation_feather(path: Path) -> dict:
    return validate_feather(path, require_sidecar=True)


def _atomic_write_h5ad(adata, output_path: Path, expected_keys) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(f".{output_path.name}.tmp.{os.getpid()}")
    backup = output_path.with_name(f".{output_path.name}.previous.{os.getpid()}")
    sidecar = Path(f"{output_path}.md5")
    sidecar_backup = output_path.with_name(f".{output_path.name}.md5.previous.{os.getpid()}")
    had_output = output_path.exists()
    had_sidecar = sidecar.exists()
    installed = False
    sidecar_installed = False
    try:
        adata.write_h5ad(temporary)
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise ValueError(f"atomic annotation merge produced an empty file: {temporary}")
        validate_h5ad(temporary, expected_keys=expected_keys)
        if had_output:
            os.link(output_path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(temporary, output_path)
        installed = True
        _write_sidecar(output_path)
        sidecar_installed = True
        validate_h5ad(output_path, expected_keys=expected_keys)
        validate_sidecar(output_path)
    except Exception:
        if installed and output_path.exists():
            output_path.unlink()
        if had_output and backup.exists():
            os.replace(backup, output_path)
        if sidecar_installed and sidecar.exists():
            sidecar.unlink()
        if had_sidecar and sidecar_backup.exists():
            os.replace(sidecar_backup, sidecar)
        elif not had_sidecar and sidecar.exists():
            sidecar.unlink()
        raise
    finally:
        for path in (temporary, backup, sidecar_backup):
            if path.exists():
                path.unlink()

def merge_annotations(
    h5ad_path: str,
    annot_dir: str,
    output_path: str | None = None,
    *,
    expected_keys=None,
    union_path: str | None = None,
):
    """Validate one run-owned annotation union, then project it to a view."""
    source_path = Path(h5ad_path)
    destination = Path(output_path or h5ad_path)
    validate_sidecar(source_path)
    annot_files = sorted(glob.glob(f"{annot_dir}/annotations_chunk_*.feather"))
    if not annot_files:
        raise ValueError(f"No annotation feather files found in {annot_dir}")

    frames = []
    for filename in annot_files:
        path = Path(filename)
        _validate_annotation_feather(path)
        frames.append(pd.read_feather(path))
    annotations = pd.concat(frames, ignore_index=True)
    if annotations.empty:
        raise ValueError(f"annotation feathers contain no rows: {annot_dir}")

    annotation_keys = _keys_from_frame(annotations)
    if len(annotation_keys) != len(set(annotation_keys)):
        raise ValueError("duplicate (Sample, cell_barcode) keys across annotation chunks")
    target_expected = (
        h5ad_keys(source_path)
        if expected_keys is None
        else _normalize_expected_keys(expected_keys)
    )
    if target_expected is None:
        raise ValueError("target annotation keys could not be derived")
    if union_path is None:
        union_expected = target_expected
    else:
        union = Path(union_path)
        validate_sidecar(union)
        union_expected = h5ad_keys(union)
        if not target_expected.issubset(union_expected):
            raise ValueError("target h5ad keys are not contained in the annotation union")
    _validate_frame(
        annotations,
        path=Path(annot_files[0]),
        expected_keys=union_expected,
        require_dataset_anchor=True,
    )

    adata = ad.read_h5ad(source_path)
    sample_col = "Sample"
    if sample_col not in adata.obs.columns:
        raise ValueError(f"h5ad lacks obs['Sample']: {source_path}")
    current_keys = h5ad_keys(source_path)
    if current_keys != target_expected:
        raise ValueError("expected annotation keys do not match target h5ad")

    obs = adata.obs.copy()
    stale = set(KNOWN_ANNOTATION_COLUMNS)
    stale.update(c for c in obs.columns if any(re.search(p, c) for p in LEGACY_ANNOT_TIER1))
    stale.update(set(LEGACY_ANNOT_TIER2) & set(obs.columns))
    obs = obs.drop(columns=[c for c in stale if c in obs.columns])
    obs_keys = [
        (str(sample), str(barcode))
        for sample, barcode in zip(obs[sample_col], adata.obs_names)
    ]
    if len(obs_keys) != len(set(obs_keys)) or set(obs_keys) != target_expected:
        raise ValueError("target h5ad has duplicate or mismatched annotation keys")

    annotations["__annot_key"] = [f"{sample}\x1f{barcode}" for sample, barcode in annotation_keys]
    obs["__annot_key"] = [f"{sample}\x1f{barcode}" for sample, barcode in obs_keys]
    fresh_columns = [column for column in ANNOT_OUTPUT_COLS if column in annotations.columns]
    fresh = annotations.set_index("__annot_key")[fresh_columns]
    obs = obs.join(fresh, on="__annot_key", how="left", validate="one_to_one")
    obs = obs.drop(columns=["__annot_key"])
    if not REQUIRED_ANNOTATION_COLUMNS.issubset(obs.columns):
        raise ValueError("merged h5ad is missing required dual-method columns")
    adata.obs = _coerce_annotation_columns(obs, fresh_columns)
    adata.uns["ecoda_annotation_version"] = "1"
    _atomic_write_h5ad(adata, destination, target_expected)
    return destination
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge per-chunk dual-method annotations into an h5ad."
    )
    parser.add_argument("--h5ad-path", required=True)
    parser.add_argument("--annot-dir", required=True)
    parser.add_argument("--output-path", default=None)
    parser.add_argument("--union-path", default=None)
    args = parser.parse_args()
    merge_annotations(args.h5ad_path, args.annot_dir, args.output_path, union_path=args.union_path)
