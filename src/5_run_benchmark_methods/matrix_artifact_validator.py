"""Validate the exact selected Pipeline 5 artifacts before synchronization."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.py.h5ad_source_identity import (  # noqa: E402
    load_source_identity,
    read_h5ad_sample_ids,
    resolve_h5ad_path,
    verify_source_identity,
)

PB_VARIANTS = ("schvg2000", "hvg2000", "hvg500", "hvg2000_bl", "hvg1000", "hvg3000")
BATCH_DATASET_ORDER = (
    "Alzheimer",
    "Breast_cancer",
    "Covid19_PBMC",
    "Kidney_KPMP",
    "Myocardial_infarction",
    "Diabetes",
    "Lupus_PBMC",
    "Lung",
    "Parkinson",
    "Joanito",
    "Stephenson",
    "CombinedPBMC",
)
PYTHON_METHODS = {"mrvi", "scpoli", "pilot", "qot", "pilotgm"}
R_METHODS = {"gloscope", "mofa", "pseudobulk", "composition", "scitd"}
CONSUMES_PSEUDOBULK = {"mofa", "pseudobulk", "composition"}


def checksum_ok(path: Path) -> bool:
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or path.stat().st_size == 0 or not sidecar.is_file():
        return False
    records = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            records[key] = value
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return (
        records.get("PATH") == str(path)
        and records.get("MD5") == digest.hexdigest()
        and records.get("SIZE") == str(path.stat().st_size)
    )

def require_nonempty(
    paths: list[Path],
    description: str,
    expected_samples: list[str] | None = None,
) -> None:
    if not paths:
        raise ValueError(f"missing/invalid {description}: []")
    expected = None if expected_samples is None else list(expected_samples)
    for path in paths:
        if not checksum_ok(path):
            raise ValueError(f"missing/invalid {description}: {path}")
        if path.suffix != ".feather":
            continue
        frame = pd.read_feather(path)
        if frame.empty:
            raise ValueError(f"empty Feather output for {description}: {path}")
        id_columns = [
            column
            for column in ("__index_level_0__", "Sample", "sample")
            if column in frame.columns
        ]
        id_column = id_columns[0] if id_columns else None
        if id_column is None:
            if isinstance(frame.index, pd.RangeIndex):
                raise ValueError(
                    f"Feather output has no sample identifier for {description}: {path}"
                )
            ids = [str(value) if value is not None else "" for value in frame.index]
        else:
            ids = [str(value) if value is not None else "" for value in frame[id_column]]
        if any(not value.strip() for value in ids):
            raise ValueError(f"Feather output has blank sample identifiers: {path}")
        if len(ids) != len(set(ids)):
            raise ValueError(f"Feather output has duplicate sample identifiers: {path}")
        if expected is not None and ids != expected:
            raise ValueError(
                f"Feather sample identifiers do not match ordered h5ad samples for {path}"
            )
        feature_columns = [column for column in frame.columns if column != id_column]
        if not feature_columns:
            raise ValueError(f"Feather output has no feature columns: {path}")
        values = frame[feature_columns].apply(pd.to_numeric, errors="coerce")
        if values.isna().all(axis=None):
            raise ValueError(f"Feather output has no numeric finite features: {path}")
        if not np.isfinite(values.to_numpy(dtype=float, na_value=np.nan)).all():
            raise ValueError(f"Feather output has nonfinite features: {path}")
        if "_dists.feather" in path.name:
            if len(feature_columns) != len(ids) or feature_columns != ids:
                raise ValueError(f"distance Feather is not square with ordered IDs: {path}")


def expected_artifacts(root: Path, ds: str, label: str, batch: bool, batch_pass: str | None) -> list[Path]:
    if label == "prepare_pseudobulk":
        if batch:
            return [
                root
                / "pseudobulks"
                / f"{ds}_batch_effect_{batch_pass or 'uncorrected'}_pseudobulk_hvg2000.rds"
            ]
        return [
            root / "pseudobulks" / f"{ds}_pseudobulk_{variant}.rds"
            for variant in PB_VARIANTS
        ]
    if label in {"trans", "zeroimp"}:
        return [root / "results" / f"{ds}_{label}.rds"]
    if label in R_METHODS:
        stem = f"{ds}_batch_effect_{batch_pass or 'uncorrected'}" if batch else ds
        return [root / "results" / f"{stem}_{label}.rds"]
    if label == "mrvi":
        if batch:
            return [
                root
                / "embeddings"
                / f"{ds}_batch_effect_{batch_pass or 'uncorrected'}_hvg2000_highres_mrvi_dists.feather"
            ]
        return [root / "embeddings" / f"{ds}_hvg{n}_mrvi_dists.feather" for n in (1000, 2000, 3000)]
    if label == "scpoli":
        if batch:
            raise ValueError("scPoli is not supported in batch-effect mode")
        paths = [root / "embeddings" / f"{ds}_hvg2000_lowres_scpoli_dims15_embs.feather"]
        paths.extend(
            root / "embeddings" / f"{ds}_hvg{n}_highres_scpoli_dims15_embs.feather"
            for n in (1000, 3000)
        )
        paths.extend(
            root / "embeddings" / f"{ds}_hvg2000_highres_scpoli_dims{dim}_embs.feather"
            for dim in (2, 3, 5, 10, 15)
        )
        return paths
    if label in {"pilot", "qot", "pilotgm"}:
        if batch:
            return [
                root
                / "embeddings"
                / f"{ds}_batch_effect_{batch_pass or 'uncorrected'}_hvg2000_highres_{label}_dists.feather"
            ]
        paths = [root / "embeddings" / f"{ds}_hvg2000_lowres_{label}_dists.feather"]
        paths.extend(
            root / "embeddings" / f"{ds}_hvg{n}_highres_{label}_dists.feather"
            for n in (1000, 2000, 3000)
        )
        return paths
    raise ValueError(f"unsupported benchmark output label: {label}")


def read_selection(selection: Path) -> list[tuple[str, str, str]]:
    rows = []
    seen = set()
    for line_number, line in enumerate(selection.read_text().splitlines(), start=1):
        if not line:
            raise ValueError(f"selection contains a blank row at line {line_number}")
        parts = line.split("\t")
        if len(parts) != 3 or any(not part for part in parts):
            raise ValueError(f"selection row {line_number} must have three non-empty columns")
        row = tuple(parts)
        if row in seen:
            raise ValueError(f"selection contains duplicate row: {'/'.join(row)}")
        seen.add(row)
        rows.append(row)
    if not rows:
        raise ValueError(f"selection is empty: {selection}")
    return rows


def expected_sample_ids(
    input_root: Path | None,
    config_path: Path | None,
    ds: str,
    view: str,
    source_identity_records: dict[tuple[str, str], dict] | None = None,
) -> list[str] | None:
    if input_root is None:
        return None
    if config_path is None or not config_path.is_file():
        raise ValueError("--input-root requires an existing --config")
    h5ad_path = resolve_h5ad_path(input_root, config_path, ds, view)
    if source_identity_records is not None:
        record = source_identity_records.get((ds, view))
        if record is None or record["path"] != str(h5ad_path):
            raise ValueError(f"source identity is missing or mismatched for {ds}/{view}")
        return list(record["sample_ids"])
    if not checksum_ok(h5ad_path):
        raise ValueError(f"missing or checksum-invalid input h5ad: {h5ad_path}")
    return read_h5ad_sample_ids(h5ad_path)


def validate(
    root: Path,
    selection: Path,
    labels: list[str],
    batch: bool,
    batch_pass: str | None = None,
    exact: bool = False,
    input_root: Path | None = None,
    config_path: Path | None = None,
    source_identity: Path | None = None,
    source_identity_verified: bool = False,
) -> None:
    rows = read_selection(selection)
    allowed = list(dict.fromkeys(labels))
    source_identity_records = (
        load_source_identity(source_identity) if source_identity is not None else None
    )
    if source_identity is not None and not source_identity_verified:
        if input_root is None or config_path is None:
            raise ValueError("source identity verification requires --input-root and --config")
        verify_source_identity(source_identity, selection, input_root, config_path)
    if not allowed:
        raise ValueError("no selected benchmark labels")
    if batch and batch_pass not in {"uncorrected", "corrected"}:
        raise ValueError("batch validation requires --batch-pass")
    if batch and exact:
        expected_rows = [(ds, "batch_effect_uncorrected", "batch_effect_uncorrected")
                        for ds in BATCH_DATASET_ORDER]
        if rows != expected_rows or batch_pass != "uncorrected":
            raise ValueError("batch exact selection is not the literal twelve-row uncorrected matrix")
    for ds, view, scope in rows:
        if batch:
            expected_view = f"batch_effect_{batch_pass}"
            if view != expected_view:
                raise ValueError(f"batch selection view mismatch for {ds}: {view}")
            if scope != expected_view:
                raise ValueError(f"batch selection scope mismatch for {ds}: {scope}")
            selected_labels = allowed
        else:
            selected_labels = [scope] if exact else allowed
            if exact and scope not in allowed:
                raise ValueError(f"selection scope {scope!r} is not in --labels")
        expected = expected_sample_ids(
            input_root, config_path, ds, view, source_identity_records
        )
        for label in selected_labels:
            paths = expected_artifacts(root, ds, label, batch, batch_pass)
            require_nonempty(paths, f"{ds}/{view}/{label}", expected)
    partials = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and any(token in path.name for token in (".tmp.", ".build.", ".partial"))
    ]
    if partials:
        raise ValueError(f"partial benchmark artifacts remain: {partials}")



def validate_single(path: Path) -> None:
    require_nonempty([path], "benchmark artifact")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--artifact", type=Path)
    group.add_argument("--root", type=Path)
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--batch-pass", default=None)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--source-identity", type=Path, default=None)
    parser.add_argument("--source-identity-verified", action="store_true")
    args = parser.parse_args()
    if args.artifact is not None:
        validate_single(args.artifact)
    else:
        if args.selection is None or not args.labels:
            parser.error("--root requires --selection and --labels")
        validate(
            args.root,
            args.selection,
            args.labels,
            args.batch,
            args.batch_pass,
            args.exact,
            args.input_root,
            args.config,
            args.source_identity,
            args.source_identity_verified,
        )
    print("matrix artifact contract OK")




if __name__ == "__main__":
    main()
