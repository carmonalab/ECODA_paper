"""Validate the exact selected Pipeline 5 artifacts before synchronization."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
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
    "Kidney_KPMP_full",
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

_CHECKSUM_FIELDS = ("MD5", "SIZE", "PATH")
_ARTIFACT_RECORD_FIELDS = ("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
_MD5_RE = re.compile(r"^[0-9a-f]{32}$")
_RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
_PARTIAL_MARKERS = (".tmp", ".build", ".partial")


def _partial_name_patterns(path: Path) -> tuple[str, ...]:
    """Return only atomic-temp names adjacent to one selected path."""
    bases = (Path(path), Path(f"{path}.md5"))
    patterns: list[str] = []
    for base in bases:
        for name in (base.name, f".{base.name}"):
            for marker in _PARTIAL_MARKERS:
                patterns.extend((f"{name}{marker}", f"{name}{marker}.*"))
    return tuple(patterns)


def _selected_partial_paths(
    paths: list[Path],
    producer_run_id: str | None = None,
) -> list[Path]:
    """Find partial names for selected outputs without discovering a root."""
    selected = [Path(path) for path in paths]
    for path in tuple(selected):
        record = _record_candidate(path, producer_run_id)
        if record is not None:
            selected.append(record)
    partials: set[Path] = set()
    for path in selected:
        try:
            for pattern in _partial_name_patterns(path):
                partials.update(path.parent.glob(pattern))
        except (OSError, RuntimeError) as exc:
            raise ValueError(
                f"unable to inspect adjacent partial artifacts for {path}"
            ) from exc
    return sorted(
        (path for path in partials if path.exists() or path.is_symlink()),
        key=str,
    )


def _reject_selected_partials(
    paths: list[Path],
    producer_run_id: str | None = None,
) -> None:
    partials = _selected_partial_paths(paths, producer_run_id)
    if partials:
        raise ValueError(f"partial benchmark artifacts remain: {partials}")




def _read_checksum_sidecar(path: Path) -> dict[str, str]:
    """Read one exact MD5/SIZE/PATH sidecar without hashing ``path``."""
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or path.stat().st_size <= 0 or not sidecar.is_file():
        raise ValueError(f"checksum sidecar is missing: {sidecar}")
    try:
        lines = sidecar.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"checksum sidecar is unreadable: {sidecar}") from exc
    if len(lines) != len(_CHECKSUM_FIELDS):
        raise ValueError(f"checksum sidecar has an invalid schema: {sidecar}")
    records: dict[str, str] = {}
    for key, line in zip(_CHECKSUM_FIELDS, lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in records:
            raise ValueError(f"checksum sidecar has an invalid schema: {sidecar}")
        records[key] = line[len(prefix):]
    if not _MD5_RE.fullmatch(records["MD5"]):
        raise ValueError(f"checksum sidecar has an invalid MD5: {sidecar}")
    try:
        size = int(records["SIZE"])
    except ValueError as exc:
        raise ValueError(f"checksum sidecar has an invalid SIZE: {sidecar}") from exc
    if size <= 0 or str(size) != records["SIZE"]:
        raise ValueError(f"checksum sidecar has an invalid SIZE: {sidecar}")
    if records["PATH"] != str(path):
        raise ValueError(f"checksum sidecar has the wrong PATH: {sidecar}")
    if size != path.stat().st_size:
        raise ValueError(f"checksum sidecar has the wrong SIZE: {sidecar}")
    return records


def _full_checksum(path: Path) -> dict[str, str]:
    """Strictly verify a sidecar and the bytes it describes."""
    records = _read_checksum_sidecar(path)
    digest = hashlib.md5()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ValueError(f"cannot read artifact: {path}") from exc
    if digest.hexdigest() != records["MD5"]:
        raise ValueError(f"checksum sidecar does not match artifact: {path}")
    if path.stat().st_size != int(records["SIZE"]):
        raise ValueError(f"artifact changed during checksum validation: {path}")
    return records


def checksum_ok(path: Path) -> bool:
    """Return whether ``path`` has an exact, content-matching MD5 sidecar."""
    try:
        _full_checksum(Path(path))
    except (OSError, ValueError):
        return False
    return True


def _record_context(run_id: str | None = None) -> tuple[Path, str] | None:
    """Resolve explicit/current run metadata without scanning other runs."""
    selected_run_id = (
        run_id
        or os.environ.get("ECODA_PRODUCER_RUN_ID")
        or os.environ.get("ECODA_RUN_ID")
    )
    runs_root = os.environ.get("ECODA_RUNS_ROOT", "")
    if not runs_root:
        scratch_root = os.environ.get("HPC_SCRATCH_DIR") or os.environ.get(
            "ECODA_SCRATCH_ROOT", ""
        )
        if scratch_root:
            runs_root = str(Path(scratch_root) / "_ecoda_runs")
    if not selected_run_id or not runs_root:
        return None
    if not _RUN_ID_RE.fullmatch(selected_run_id):
        raise ValueError(f"invalid producer run ID: {selected_run_id!r}")
    return Path(runs_root), selected_run_id


def artifact_record_path(path: Path, run_id: str) -> Path:
    """Return the bounded record path for one canonical artifact path."""
    context = _record_context(run_id)
    if context is None:
        raise ValueError("artifact records require a run ID and run root")
    runs_root, selected_run_id = context
    canonical = Path(path).resolve()
    path_digest = hashlib.sha256(str(canonical).encode("utf-8")).hexdigest()[:32]
    return (
        runs_root
        / selected_run_id
        / "manifests"
        / "artifacts"
        / f"{path_digest}.record"
    )


def _record_candidate(path: Path, run_id: str | None = None) -> Path | None:
    context = _record_context(run_id)
    if context is None:
        return None
    return artifact_record_path(path, context[1])


def _read_artifact_record(
    path: Path,
    producer: str | None = None,
    run_id: str | None = None,
    *,
    require: bool = False,
) -> dict[str, str] | None:
    """Validate a run-owned record and its sidecar fields, without hashing bytes."""
    candidate = _record_candidate(path, run_id)
    if candidate is None:
        if require:
            raise ValueError(f"artifact record context is missing: {path}")
        return None
    if not candidate.exists() and not candidate.is_symlink():
        if require:
            raise ValueError(f"artifact record is missing: {candidate}")
        return None
    if not candidate.is_file():
        raise ValueError(f"artifact record is not a regular file: {candidate}")
    try:
        lines = candidate.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"artifact record is unreadable: {candidate}") from exc
    if len(lines) != len(_ARTIFACT_RECORD_FIELDS):
        raise ValueError(f"artifact record has an invalid schema: {candidate}")
    records: dict[str, str] = {}
    for key, line in zip(_ARTIFACT_RECORD_FIELDS, lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in records:
            raise ValueError(f"artifact record has an invalid schema: {candidate}")
        records[key] = line[len(prefix):]
    canonical = str(path.resolve())
    expected_run_id = (_record_context(run_id) or (None, ""))[1]
    if records["PATH"] != canonical:
        raise ValueError(f"artifact record has the wrong PATH: {candidate}")
    if not _MD5_RE.fullmatch(records["MD5"]):
        raise ValueError(f"artifact record has an invalid MD5: {candidate}")
    try:
        size = int(records["SIZE"])
    except ValueError as exc:
        raise ValueError(f"artifact record has an invalid SIZE: {candidate}") from exc
    if size <= 0 or str(size) != records["SIZE"]:
        raise ValueError(f"artifact record has an invalid SIZE: {candidate}")
    if records["RUN_ID"] != expected_run_id:
        raise ValueError(f"artifact record has the wrong RUN_ID: {candidate}")
    if not records["PRODUCER"] or any(
        char in records["PRODUCER"] for char in "\t\r\n"
    ):
        raise ValueError(f"artifact record has an invalid PRODUCER: {candidate}")
    if producer is not None and records["PRODUCER"] != producer:
        raise ValueError(f"artifact record has the wrong PRODUCER: {candidate}")
    if records["STATE"] != "PUBLISHED":
        raise ValueError(f"artifact record is not published: {candidate}")
    if not path.is_file() or path.stat().st_size != size:
        raise ValueError(f"artifact record SIZE does not match artifact: {path}")
    sidecar = _read_checksum_sidecar(path)
    if (
        sidecar["PATH"] != str(path)
        or sidecar["SIZE"] != records["SIZE"]
        or sidecar["MD5"] != records["MD5"]
    ):
        raise ValueError(f"artifact record does not match checksum sidecar: {path}")
    return records


def validate_artifact_record(
    path: Path,
    producer: str,
    run_id: str,
) -> dict[str, str]:
    """Require and validate one exact run-owned artifact record."""
    record = _read_artifact_record(path, producer, run_id, require=True)
    assert record is not None
    return record

def require_nonempty(
    paths: list[Path],
    description: str,
    expected_samples: list[str] | None = None,
    producer: str | None = None,
    producer_run_id: str | None = None,
) -> None:
    if not paths:
        raise ValueError(f"missing/invalid {description}: []")
    expected = None if expected_samples is None else list(expected_samples)
    for raw_path in paths:
        path = Path(raw_path)
        candidate = _record_candidate(path, producer_run_id)
        record_present = candidate is not None and (
            candidate.exists() or candidate.is_symlink()
        )
        # Feather is deserialized below, so its sidecar/content pair must be
        # proved immediately before the read.  A run-owned record can
        # additionally bind the digest/producer without replacing this hash.
        # RDS and other non-deserialized artifacts may use the run-owned
        # record, but only when that record is actually present.
        try:
            if path.suffix.lower() == ".feather":
                _full_checksum(path)
                if record_present:
                    _read_artifact_record(
                        path, producer, producer_run_id, require=True
                    )
            elif record_present:
                _read_artifact_record(
                    path, producer, producer_run_id, require=True
                )
            else:
                _full_checksum(path)
        except (OSError, ValueError) as exc:
            raise ValueError(f"missing/invalid {description}: {path}") from exc
        if path.suffix.lower() != ".feather":
            continue
        try:
            frame = pd.read_feather(path)
        except Exception as exc:
            raise ValueError(f"invalid Feather output for {description}: {path}") from exc
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
            if label == "pilotgm":
                raise ValueError("PILOT-GM-VAE is not scheduled in batch-effect mode")
            return [
                root
                / "embeddings"
                / f"{ds}_batch_effect_{batch_pass or 'uncorrected'}_hvg2000_highres_{label}_dists.feather"
            ]
        if label == "pilotgm":
            return [
                root / "embeddings" / f"{ds}_hvg2000_highres_pilotgm_dists.feather"
            ]
        paths = [root / "embeddings" / f"{ds}_hvg2000_lowres_{label}_dists.feather"]
        paths.extend(
            root / "embeddings" / f"{ds}_hvg{n}_highres_{label}_dists.feather"
            for n in (1000, 2000, 3000)
        )
        return paths
    raise ValueError(f"unsupported benchmark output label: {label}")


def read_selection(selection: Path) -> list[tuple[str, str, str]]:
    selection = Path(selection)
    if not checksum_ok(selection):
        raise ValueError(f"selection checksum is missing or invalid: {selection}")
    try:
        lines = selection.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"selection is unreadable: {selection}") from exc
    rows = []
    seen = set()
    for line_number, line in enumerate(lines, start=1):
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
    producer_run_id: str | None = None,
    producer: str | None = None,
) -> None:
    rows = read_selection(selection)
    selected_paths = [Path(selection)]
    if source_identity is not None:
        selected_paths.append(Path(source_identity))
    allowed = list(dict.fromkeys(labels))
    if source_identity is not None:
        try:
            _full_checksum(Path(source_identity))
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"source identity checksum is missing or invalid: {source_identity}"
            ) from exc
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
        expected_rows = [
            (ds, "batch_effect_uncorrected", "batch_effect_uncorrected")
            for ds in BATCH_DATASET_ORDER
        ]
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
            selected_paths.extend(paths)
            require_nonempty(
                paths,
                f"{ds}/{view}/{label}",
                expected,
                producer=producer or label,
                producer_run_id=producer_run_id,
            )
    _reject_selected_partials(selected_paths, producer_run_id)



def validate_single(
    path: Path,
    producer: str | None = None,
    producer_run_id: str | None = None,
) -> None:
    require_nonempty(
        [path],
        "benchmark artifact",
        producer=producer,
        producer_run_id=producer_run_id,
    )
    _reject_selected_partials([Path(path)], producer_run_id)
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
    parser.add_argument("--producer", default=None)
    parser.add_argument("--producer-run-id", default=None)
    args = parser.parse_args()
    if args.artifact is not None:
        validate_single(
            args.artifact,
            producer=args.producer,
            producer_run_id=args.producer_run_id,
        )
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
            producer_run_id=args.producer_run_id,
            producer=args.producer,
        )
    print("matrix artifact contract OK")




if __name__ == "__main__":
    main()
