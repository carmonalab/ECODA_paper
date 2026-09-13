#!/usr/bin/env python3
"""Audit declared preprocessing views before sample-count filtering.

The audit deliberately stops before ``base_preprocessing``. It applies the
same registry subset and sample-ID standardization as the Stage 3 worker,
then records raw per-sample counts and configured high-resolution annotation
coverage. Biological labels are reported only as metadata; they never enter
selection or validation decisions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path

import anndata as ad
import h5py
from anndata._io.h5ad import read_dataframe
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.utils.py.datasets_io import read_datasets_json
from src.utils.py.preprocess_utils import (
    assert_subset_sample_consistency,
    evaluate_subset_mask,
    load_input,
)


def read_obs_only(path: Path) -> ad.AnnData:
    """Read only the H5AD observation table, never the count matrix."""
    with h5py.File(path, "r") as handle:
        return ad.AnnData(obs=read_dataframe(handle["obs"]))


def _file_digest(path: Path, algorithm="sha256") -> str | None:
    """Return a file digest, or ``None`` when an optional identity is absent."""
    if not path.is_file():
        return None
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def _file_identity(path_value: str | None) -> dict:
    if not path_value:
        return {"path": None, "size": None, "md5": None, "sha256": None}
    path = Path(path_value).expanduser()
    identity = {
        "path": str(path.resolve(strict=False)),
        "size": None,
        "md5": None,
        "sha256": None,
    }
    if path.is_file():
        identity["size"] = int(path.stat().st_size)
        identity["md5"] = _file_digest(path, "md5")
        identity["sha256"] = _file_digest(path, "sha256")
    return identity


def _runtime_source_provenance() -> dict:
    source_root = os.environ.get("ECODA_SOURCE_ROOT")
    source_manifest = os.environ.get("ECODA_SOURCE_MANIFEST")
    source_manifest_run = os.environ.get("ECODA_SOURCE_MANIFEST_RUN")
    runtime_identity = os.environ.get("ECODA_RUNTIME_IDENTITY")
    runtime_image = os.environ.get("ECODA_RUNTIME_IMAGE")
    runtime_manifest = os.environ.get("ECODA_RUNTIME_MANIFEST")
    return {
        "run_id": os.environ.get("ECODA_RUN_ID"),
        "run_root": os.environ.get("ECODA_RUN_ROOT"),
        "source_root": source_root,
        "source_manifest": _file_identity(source_manifest),
        "source_manifest_run": _file_identity(source_manifest_run),
        "runtime_identity": _file_identity(runtime_identity),
        "runtime_image": _file_identity(runtime_image),
        "runtime_manifest": _file_identity(runtime_manifest),
    }


def _json_scalar(value):
    if value is None or value is pd.NA:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _sampling_day_audit(series: pd.Series, include_values) -> dict:
    """Describe raw sampling-day values without changing their representation."""
    trimmed = series.astype("string").str.strip()
    numeric = pd.to_numeric(trimmed, errors="coerce")
    numeric_values = np.asarray(numeric.fillna(np.nan), dtype=float)
    missing = series.isna()
    blank = (~missing) & trimmed.eq("")
    non_finite = (
        (~missing) & (~blank) & numeric.notna() & ~np.isfinite(numeric_values)
    )
    non_numeric = (~missing) & (~blank) & numeric.isna()
    if isinstance(include_values, str):
        include_values = [include_values]
    include_text = {
        str(value).strip().casefold()
        for value in (include_values or [])
        if value is not None and not pd.isna(value)
    }
    unknown = non_numeric & ~trimmed.str.casefold().isin(include_text)
    malformed = unknown & ~trimmed.str.casefold().eq("unknown")
    raw_unique = [_json_scalar(value) for value in pd.unique(series.astype(object))]
    return {
        "column": str(series.name),
        "dtype": str(series.dtype),
        "raw_unique_values": raw_unique,
        "missing_count": int(missing.sum()),
        "blank_count": int(blank.sum()),
        "non_finite_count": int(non_finite.sum()),
        "non_numeric_count": int(non_numeric.sum()),
        "unknown_count": int(
            (unknown & trimmed.str.casefold().eq("unknown")).sum()
        ),
        "malformed_count": int(malformed.sum()),
        "include_values_count": int(
            trimmed.str.casefold().isin(include_text).sum()
        ),
        "finite_numeric_count": int(
            ((~missing) & (~blank) & np.isfinite(numeric_values)).sum()
        ),
    }


def load_audit_input(input_file, input_dir: Path, output_dir: Path):
    """Load only AnnData metadata when a validated raw cache is available."""
    if isinstance(input_file, list):
        metadata = []
        for name in input_file:
            source_path = input_dir / name
            if not source_path.exists() or not str(name).endswith(".h5ad"):
                return load_input(input_file, input_dir, output_dir)
            metadata.append(read_obs_only(source_path).obs)
        combined = pd.concat(metadata, axis=0)
        combined.index = [f"source_{position}" for position in range(len(combined))]
        return ad.AnnData(obs=combined)

    input_path = input_dir / input_file
    if str(input_file).endswith(".rds"):
        raw_cache = output_dir / f"{Path(input_file).stem}_raw.h5ad"
        if raw_cache.exists():
            return read_obs_only(raw_cache)
    elif str(input_file).endswith(".h5ad") and input_path.exists():
        return read_obs_only(input_path)
    return load_input(input_file, input_dir, output_dir)


SENTINEL_ANNOTATIONS = frozenset(
    {
        "",
        "na",
        "n/a",
        "nan",
        "none",
        "null",
        "missing",
        "not available",
        "not_available",
    }
)


def standardize_sample_ids(values: pd.Series, *, dataset: str, view: str) -> list[str]:
    """Validate and standardize sample IDs exactly as Stage 3 does."""
    standardized = []
    for position, value in enumerate(values):
        if pd.isna(value):
            raise ValueError(
                f"{dataset} / {view}: missing sample ID at observation {position}"
            )
        text = str(value)
        if not text.strip():
            raise ValueError(
                f"{dataset} / {view}: empty sample ID at observation {position}"
            )
        standardized.append(f"g{text}" if re.match(r"^\d", text) else text.replace("-", "_"))
    return standardized


def annotation_coverage(series: pd.Series) -> dict[str, int]:
    """Return total, valid, and missing/sentinel annotation counts."""
    total = int(series.shape[0])
    missing = series.isna()
    normalized = series.astype("string").str.strip().str.casefold()
    sentinel = normalized.isin(SENTINEL_ANNOTATIONS)
    invalid = missing | sentinel
    return {
        "total_cells": total,
        "valid_nonmissing": int((~invalid).sum()),
        "invalid_missing_or_sentinel": int(invalid.sum()),
    }


def audit_view(
    dataset: str,
    entry: dict,
    view_name: str,
    view: dict,
    input_root: Path,
    output_root: Path,
    direct_input_path: Path | None = None,
) -> dict:
    """Audit one registry view without running preprocessing."""
    input_file = view.get("input_file")
    if not input_file:
        raise ValueError(f"{dataset} / {view_name}: missing input_file_name")

    if direct_input_path is None:
        input_dir = input_root / dataset / "data"
        output_dir = output_root / dataset / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        adata = load_audit_input(input_file, input_dir, output_dir)
        input_identity = None
    else:
        input_path = Path(direct_input_path).resolve()
        input_identity = _file_identity(str(input_path))
        adata = read_obs_only(input_path)

    try:
        sample_col = (view.get("columns") or {}).get("sample") or entry["sample_col"]
        subset_vars = view.get("subset_vars", {})
        subset_mask = evaluate_subset_mask(adata, subset_vars)
        subset_audit = assert_subset_sample_consistency(
            adata,
            subset_mask,
            sample_col,
            context=f"{dataset} / {view_name}",
        )
        subset = adata[subset_mask]
        if subset.n_obs == 0:
            raise ValueError(
                f"{dataset} / {view_name}: subset_vars produced an empty view: "
                f"{subset_vars}"
            )

        sample_ids = standardize_sample_ids(
            subset.obs[sample_col], dataset=dataset, view=view_name
        )
        sample_counts = Counter(sample_ids)

        high_res_col = (view.get("columns") or {}).get(
            "cell_type_high_res"
        ) or entry["cell_type_high_res"]
        if high_res_col in subset.obs.columns:
            coverage = annotation_coverage(subset.obs[high_res_col])
            coverage["status"] = "present"
        else:
            coverage = {
                "total_cells": int(subset.n_obs),
                "valid_nonmissing": 0,
                "invalid_missing_or_sentinel": int(subset.n_obs),
                "status": "missing_column",
            }

        record = {
            "dataset": dataset,
            "view": view_name,
            "input_file": input_file,
            "input_identity": input_identity,
            "subset_vars": subset_vars,
            "sample_column": sample_col,
            "sample_counts": dict(sorted(sample_counts.items())),
            # ``total_cells``/``total_samples`` retain their historical
            # post-subset meaning; explicit fields describe the complete
            # source and predicate result.
            "total_cells": int(subset.n_obs),
            "total_samples": int(len(sample_counts)),
            "input_total_cells": subset_audit["total_cells"],
            "input_total_samples": subset_audit["total_samples"],
            "retained_cells": subset_audit["retained_cells"],
            "dropped_cells": subset_audit["dropped_cells"],
            "retained_samples": subset_audit["retained_samples"],
            "dropped_samples": subset_audit["dropped_samples"],
            "split_sample_count": subset_audit["split_sample_count"],
            "subset_audit": subset_audit,
            "high_resolution_column": high_res_col,
            "annotation_coverage": coverage,
        }
        if direct_input_path is not None:
            sampling_column = "Sampling day (Days after symptom onset)"
            rule = subset_vars.get(sampling_column)
            if not isinstance(rule, dict):
                raise ValueError(
                    f"{dataset} / {view_name}: missing configured sampling-day rule"
                )
            include_values = rule.get("include_values", [])
            record["sampling_day_column"] = sampling_column
            record["sampling_day_audit"] = _sampling_day_audit(
                adata.obs[sampling_column], include_values
            )
            configured_cardinalities = {}
            for column in ("sampleID", "PatientID"):
                if column not in adata.obs.columns:
                    raise KeyError(
                        f"{dataset} / {view_name}: required Covid obs column "
                        f"{column!r} is missing"
                    )
                configured_cardinalities[column] = int(
                    adata.obs[column].nunique(dropna=False)
                )
            record["configured_cardinalities"] = configured_cardinalities
            record["sampleID_cardinality"] = configured_cardinalities["sampleID"]
            record["PatientID_cardinality"] = configured_cardinalities["PatientID"]
        return record
    finally:
        backing_file = getattr(adata, "file", None)
        if backing_file is not None:
            backing_file.close()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def write_checksum_atomic(path: Path) -> None:
    """Write the run-owned checksum sidecar for an audit report."""
    digest = _file_digest(path, "md5")
    if digest is None:
        raise ValueError(f"cannot checksum missing audit report: {path}")
    sidecar = Path(f"{path}.md5")
    temporary = sidecar.with_name(f".{sidecar.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            f"MD5={digest}\nSIZE={path.stat().st_size}\nPATH={path}\n",
            encoding="utf-8",
        )
        os.replace(temporary, sidecar)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _validate_path_ancestry(
    path: Path,
    label: str,
    *,
    boundary: Path | None = None,
) -> None:
    """Reject traversal and unsafe ancestry before any mkdir.

    When a trusted run-root boundary is available, symlink aliases above that
    canonical boundary are allowed (for example macOS ``/var``), while every
    existing component below the boundary is checked strictly.
    """
    if not path.is_absolute():
        raise ValueError(f"{label} must be an absolute path")
    if ".." in path.parts:
        raise ValueError(f"{label} must not contain '..' path traversal")

    boundary_real = None
    if boundary is not None:
        if not boundary.is_absolute() or ".." in boundary.parts:
            raise ValueError("ECODA_RUN_ROOT must be an absolute, canonical path")
        try:
            boundary_real = boundary.resolve(strict=False)
        except (OSError, RuntimeError) as exc:
            raise ValueError("ECODA_RUN_ROOT cannot be resolved safely") from exc

    current = Path(path.anchor)
    parts = path.parts[1:] if path.anchor else path.parts
    boundary_reached = boundary_real is None
    for index, part in enumerate(parts):
        current /= part
        try:
            current.lstat()
        except FileNotFoundError:
            # A missing ancestor means all remaining descendants are also
            # absent; they will be created only after the full validation.
            break
        except OSError as exc:
            raise ValueError(f"{label} ancestry cannot be inspected") from exc

        if boundary_real is not None and not boundary_reached:
            try:
                current_real = current.resolve(strict=False)
            except (OSError, RuntimeError) as exc:
                raise ValueError(f"{label} ancestry cannot be resolved") from exc
            if current_real in boundary_real.parents:
                # Existing ancestors above the bound are outside this
                # output's security domain and may use platform aliases.
                continue
            if current_real == boundary_real or boundary_real in current_real.parents:
                # The lexical path is now at or below the trusted boundary;
                # every component from here onward must be non-symlinked.
                boundary_reached = True

        if current.is_symlink():
            raise ValueError(f"{label} ancestry contains a symlink: {current}")
        if index < len(parts) - 1 and not current.is_dir():
            raise ValueError(f"{label} ancestry is not directory-backed: {current}")


def _validate_direct_output_scope(
    output_root: Path,
    output: Path,
    run_root: Path,
) -> tuple[Path, Path]:
    """Validate an obs-only report path without creating any directories."""
    output_root = output_root.expanduser()
    output = output.expanduser()
    run_root = run_root.expanduser()
    if not output_root.is_absolute() or not output.is_absolute():
        raise ValueError("obs-only output-root and output must be absolute paths")
    if not run_root.is_absolute():
        raise ValueError("ECODA_RUN_ROOT must be an absolute path")
    if ".." in run_root.parts:
        raise ValueError("ECODA_RUN_ROOT must not contain '..' path traversal")
    if not run_root.is_dir() or run_root.is_symlink():
        raise ValueError("ECODA_RUN_ROOT must be a real directory")

    try:
        run_real = run_root.resolve(strict=False)
        checksum = Path(f"{output}.md5")
        _validate_path_ancestry(
            output_root,
            "obs-only output-root",
            boundary=run_real,
        )
        _validate_path_ancestry(
            output,
            "obs-only report",
            boundary=run_real,
        )
        _validate_path_ancestry(
            checksum,
            "obs-only report checksum",
            boundary=run_real,
        )
        root_real = output_root.resolve(strict=False)
        output_real = output.resolve(strict=False)
        checksum_real = checksum.resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        raise ValueError("obs-only output path cannot be resolved safely") from exc

    if output_root.exists() and not output_root.is_dir():
        raise ValueError("obs-only output-root must be a directory")
    if output.exists() and not output.is_file():
        raise ValueError("obs-only report path must be a regular file")
    if output.exists() and output.is_symlink():
        raise ValueError("obs-only report must not replace a symlink")
    if checksum.exists() and checksum.is_symlink():
        raise ValueError("obs-only report checksum must not replace a symlink")

    if output_real == root_real or root_real not in output_real.parents:
        raise ValueError("obs-only report path escapes output-root")
    if checksum_real == root_real or root_real not in checksum_real.parents:
        raise ValueError("obs-only report checksum escapes output-root")
    if root_real != run_real and run_real not in root_real.parents:
        raise ValueError("obs-only report root escapes ECODA_RUN_ROOT")
    for candidate in (output_real, checksum_real):
        if candidate == run_real or run_real not in candidate.parents:
            raise ValueError(
                "obs-only report and checksum must remain below ECODA_RUN_ROOT"
            )

    # All checks above intentionally precede directory creation.  The caller
    # may now let the atomic writers create missing run-owned parents.
    return root_real, output


def _required_identity_path(name: str, *, directory=False) -> Path:
    value = os.environ.get(name, "")
    if not value or not Path(value).is_absolute():
        raise ValueError(f"obs-only audit requires absolute {name}")
    path = Path(value)
    if directory:
        valid = path.is_dir() and not path.is_symlink()
    else:
        valid = path.is_file() and not path.is_symlink()
    if not valid:
        raise ValueError(f"obs-only audit identity path is missing or unsafe: {name}")
    return path.resolve()


def _validate_direct_identity(args: argparse.Namespace, output_root: Path) -> dict:
    source_root = _required_identity_path("ECODA_SOURCE_ROOT", directory=True)
    source_manifest = _required_identity_path("ECODA_SOURCE_MANIFEST")
    runtime_identity = _required_identity_path("ECODA_RUNTIME_IDENTITY")
    runtime_image = _required_identity_path("ECODA_RUNTIME_IMAGE")
    runtime_manifest = _required_identity_path("ECODA_RUNTIME_MANIFEST")
    run_root = _required_identity_path("ECODA_RUN_ROOT", directory=True)
    run_id = os.environ.get("ECODA_RUN_ID", "")
    if not run_id or run_root.name != run_id:
        raise ValueError("obs-only audit run identity is missing or inconsistent")
    if source_root.name != "tree":
        raise ValueError("obs-only audit source root is not a snapshot tree")
    expected_manifest = source_root.parent / "identity" / "source.manifest"
    if source_manifest != expected_manifest.resolve():
        raise ValueError("obs-only audit source manifest is not bound to source tree")
    source_manifest_run = os.environ.get("ECODA_SOURCE_MANIFEST_RUN")
    if source_manifest_run:
        run_manifest = (run_root / "manifests" / "source.manifest").resolve()
        if Path(source_manifest_run).resolve() != run_manifest:
            raise ValueError("obs-only audit run source manifest binding mismatch")
        if not run_manifest.is_file() or run_manifest.is_symlink():
            raise ValueError("obs-only audit run source manifest is missing or unsafe")
        if run_manifest.read_bytes() != source_manifest.read_bytes():
            raise ValueError("obs-only audit run source manifest differs from snapshot")
    runtime_run_identity = (run_root / "manifests" / "runtime.identity").resolve()
    if runtime_identity != runtime_run_identity:
        raise ValueError("obs-only audit runtime identity is not run-bound")
    _validate_direct_output_scope(
        output_root,
        args.output.expanduser(),
        run_root,
    )
    config_real = args.config.expanduser().resolve()
    if config_real != (source_root / "datasets.json").resolve():
        raise ValueError("obs-only audit config is not the immutable snapshot datasets.json")
    return _runtime_source_provenance()


def _validate_covid_direct_rule(view: dict) -> None:
    column = "Sampling day (Days after symptom onset)"
    subset_vars = view.get("subset_vars")
    rule = subset_vars.get(column) if isinstance(subset_vars, dict) else None
    if not isinstance(rule, dict) or rule.get("op") != "<=":
        raise ValueError("Covid obs-only audit requires the configured sampling-day <= rule")
    values = rule.get("values")
    if isinstance(values, (list, tuple)):
        if len(values) != 1:
            raise ValueError("Covid sampling-day threshold must have one value")
        values = values[0]
    if values is None or isinstance(values, bool):
        raise ValueError("Covid sampling-day threshold must be finite numeric 30")
    threshold = pd.to_numeric(
        pd.Series([str(values).strip()]), errors="coerce"
    ).iloc[0]
    if pd.isna(threshold) or not np.isfinite(float(threshold)) or float(threshold) != 30:
        raise ValueError("Covid sampling-day threshold must be exactly 30")
    include_values = rule.get("include_values")
    if isinstance(include_values, str):
        include_values = [include_values]
    if not isinstance(include_values, (list, tuple)) or [
        str(value).strip().casefold() for value in include_values
    ] != ["control"]:
        raise ValueError("Covid sampling-day rule must include exactly 'control'")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="datasets.json", type=Path)
    parser.add_argument(
        "--input-root",
        required=False,
        type=Path,
        help="Parent containing <dataset>/data staged inputs.",
    )
    parser.add_argument(
        "--input-file",
        default=None,
        type=Path,
        help="Direct H5AD path for strict obs-only auditing.",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        type=Path,
        help="Parent containing/receiving <dataset>/output raw caches.",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--view", default="benchmark_analysis")
    parser.add_argument("--ds-name", action="append", dest="datasets")
    parser.add_argument("--include-underscore", action="store_true")
    parser.add_argument(
        "--obs-only",
        action="store_true",
        help="Audit only H5AD obs through h5py; requires direct Covid input.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.obs_only:
        if args.input_root is not None:
            raise ValueError(
                "--obs-only requires --input-file and rejects --input-root loading"
            )
        if args.input_file is None:
            raise ValueError("--obs-only requires --input-file")
        if args.output_root is None:
            raise ValueError("--obs-only requires --output-root")
        if args.datasets != ["Covid19_PBMC"]:
            raise ValueError("--obs-only requires exactly --ds-name Covid19_PBMC")
        if args.view not in {
            "batch_effect_uncorrected",
            "batch_effect_corrected",
        }:
            raise ValueError(
                "--obs-only is restricted to a Covid batch-effect view"
            )
        input_argument = args.input_file.expanduser()
        if not input_argument.is_absolute():
            raise ValueError("--obs-only --input-file must be an absolute path")
        if input_argument.is_symlink() or not input_argument.is_file():
            raise ValueError(
                f"obs-only input H5AD is missing or unsafe: {input_argument}"
            )
        if input_argument.suffix.lower() != ".h5ad":
            raise ValueError("--obs-only --input-file must name an .h5ad")
        direct_input = input_argument.resolve()
        config = read_datasets_json(args.config, view=args.view)
        try:
            entry = config["Covid19_PBMC"]
            view = entry["views"][args.view]
        except KeyError as exc:
            raise ValueError(
                f"Covid19_PBMC does not declare {args.view!r}"
            ) from exc
        configured_input = view.get("input_file")
        if not isinstance(configured_input, str) or not configured_input:
            raise ValueError("Covid obs-only view has no configured input_file_name")
        configured_name = Path(configured_input).name
        if direct_input.name != configured_name:
            raise ValueError(
                "direct Covid H5AD filename does not match configured input_file_name"
            )
        scratch_root = os.environ.get("HPC_SCRATCH_DIR")
        if scratch_root:
            expected_path = (
                Path(scratch_root)
                / "Covid19_PBMC"
                / "data"
                / configured_name
            ).resolve(strict=False)
            if direct_input != expected_path:
                raise ValueError(
                    "direct Covid H5AD path does not match HPC scratch source path"
                )
        _validate_covid_direct_rule(view)
        output_root = args.output_root.expanduser()
        _validate_direct_identity(args, output_root)
        record = audit_view(
            "Covid19_PBMC",
            entry,
            args.view,
            view,
            output_root,
            output_root,
            direct_input_path=direct_input,
        )
        record["input_path"] = record["input_identity"]["path"]
        provenance = _runtime_source_provenance()
        payload = {
            "config": str(args.config.expanduser().resolve()),
            "view": args.view,
            "threshold": 500,
            "obs_only": True,
            "audit_scope": "direct_h5ad_obs_only_after_declared_subset_before_base_preprocessing",
            "provenance": provenance,
            "immutable_source_runtime": provenance,
            "datasets": [record],
        }
        output = args.output.expanduser()
        if not output.is_absolute():
            raise ValueError("obs-only --output must be an absolute path")
        write_json_atomic(output, payload)
        write_checksum_atomic(output)
        print(f"Wrote {output} (1 dataset)")
        return

    if args.input_file is not None:
        raise ValueError("--input-file is only valid with --obs-only")
    if args.input_root is None:
        raise ValueError("--input-root is required unless --obs-only is used")
    output_root = args.output_root or args.input_root
    config = read_datasets_json(args.config, view=args.view)
    selected = []
    for dataset in sorted(config):
        if args.datasets and dataset not in args.datasets:
            continue
        if not args.include_underscore and dataset.startswith("_"):
            continue
        entry = config[dataset]
        views = entry.get("views") or {}
        if args.view not in views:
            continue
        selected.append((dataset, entry, views[args.view]))
    if not selected:
        raise ValueError(f"No datasets declare view {args.view!r}")

    records = [
        audit_view(dataset, entry, args.view, view, args.input_root, output_root)
        for dataset, entry, view in selected
    ]
    payload = {
        "config": str(args.config.resolve()),
        "view": args.view,
        "threshold": 500,
        "datasets": records,
        "audit_scope": "raw_input_after_declared_subset_before_base_preprocessing",
    }
    write_json_atomic(args.output, payload)
    print(f"Wrote {args.output} ({len(records)} datasets)")




if __name__ == "__main__":
    main()
