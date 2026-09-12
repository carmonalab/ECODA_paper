#!/usr/bin/env python3
"""Audit corrected Stage 3 H5AD sources without reading expression data.

This source-bound boundary opens an H5AD with :mod:`h5py` and reads only the
``obs`` dataframe columns required by the immutable configuration.  It never
accesses ``X``, ``raw``, ``layers``, or any other HDF5 group.  Subset semantics
come from the shared preprocessing evaluator, and corrected batch metadata is
validated by the shared Python batch contract before the run-owned report is
installed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np
import pandas as pd

# The evaluator import is deliberately from the immutable snapshot package.
# Keeping this boundary shared prevents a second implementation of subset
# operators while this script remains usable as a direct source-tree script.
try:
    from src.utils.py.batch_contract import (
        RESERVED_OBS_NAME,
        build_batch_contract_identity,
        normalize_batch_keys,
        validate_batch_metadata,
    )
    from src.utils.py.h5ad_source_identity import (
        read_obs_column_values,
        read_str_dataset,
    )
    from src.utils.py.preprocess_utils import (
        assert_subset_sample_consistency,
        evaluate_subset_mask,
    )
except ModuleNotFoundError:  # direct execution with src/utils/py on sys.path
    from batch_contract import (  # type: ignore[no-redef]
        RESERVED_OBS_NAME,
        build_batch_contract_identity,
        normalize_batch_keys,
        validate_batch_metadata,
    )
    from h5ad_source_identity import (  # type: ignore[no-redef]
        read_obs_column_values,
        read_str_dataset,
    )
    from preprocess_utils import (  # type: ignore[no-redef]
        assert_subset_sample_consistency,
        evaluate_subset_mask,
    )


REPORT_SCHEMA_VERSION = 1
CORRECTED_VIEW = "batch_effect_corrected"
CONTRACT_METHOD = "preprocess"
CONTRACT_MODEL = "hvg_composite_v1"
_SOURCE_DIGEST_SIZE = 1024 * 1024
_HEX_SHA256 = re.compile(r"^[0-9A-Fa-f]{64}$")
_MANIFEST_KEY = re.compile(r"^[A-Za-z0-9_]+$")
_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


class _ObsOnly:
    """Small AnnData-compatible view consumed by shared subset helpers."""

    __slots__ = ("obs", "obs_names")

    def __init__(self, obs: pd.DataFrame):
        self.obs = obs
        self.obs_names = obs.index


def _safe_text(value: object, label: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be one string")
    if any(character in value for character in "\r\n\t"):
        raise ValueError(f"{label} contains record-delimiter characters")
    if not allow_empty and (not value or value != value.strip()):
        raise ValueError(f"{label} must be nonblank and free of surrounding whitespace")
    return value


def _safe_absolute(value: object, label: str) -> Path:
    if isinstance(value, Path):
        text = str(value)
    else:
        text = _safe_text(value, label)
    path = Path(text)
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    if ".." in path.parts:
        raise ValueError(f"{label} must not contain parent traversal")
    return path


def _assert_no_symlink_ancestors(path: Path, label: str) -> None:
    current = path
    while True:
        if current.is_symlink():
            raise ValueError(f"{label} must not be a symlink: {current}")
        if current.parent == current:
            break
        current = current.parent


def _existing_file(value: object, label: str) -> Path:
    path = _safe_absolute(value, label)
    _assert_no_symlink_ancestors(path, label)
    if not path.is_file():
        raise ValueError(f"{label} is not an existing regular file: {path}")
    if path.stat().st_size <= 0:
        raise ValueError(f"{label} is empty: {path}")
    resolved = path.resolve(strict=True)
    if resolved != path:
        raise ValueError(f"{label} is not canonical: {path}")
    return path


def _existing_directory(value: object, label: str) -> Path:
    path = _safe_absolute(value, label)
    _assert_no_symlink_ancestors(path, label)
    if not path.is_dir():
        raise ValueError(f"{label} is not an existing directory: {path}")
    resolved = path.resolve(strict=True)
    if resolved != path:
        raise ValueError(f"{label} is not canonical: {path}")
    return path


def _run_owned_output(value: object, run_root: Path) -> Path:
    path = _safe_absolute(value, "output path")
    parent = path.parent
    _assert_no_symlink_ancestors(path, "output path")
    if not parent.is_dir() or parent.is_symlink():
        raise ValueError(f"output parent directory is missing or unsafe: {parent}")
    parent_resolved = parent.resolve(strict=True)
    if not _path_below(parent_resolved, run_root):
        raise ValueError(f"output path escapes the run root: {path}")
    if path.exists() or path.is_symlink():
        if path.is_symlink():
            raise ValueError(f"output path is a symlink: {path}")
        if not path.is_file():
            raise ValueError(f"output path is not a regular file: {path}")
    return path


def _path_below(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_SOURCE_DIGEST_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_identity(path: Path, label: str) -> dict[str, Any]:
    path = _existing_file(path, label)
    return {
        "path": str(path),
        "size": int(path.stat().st_size),
        "md5": _digest_file(path, "md5"),
        "sha256": _digest_file(path, "sha256"),
    }


def _read_manifest(path: Path, label: str) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"could not read {label}: {path}") from exc
    if not lines:
        raise ValueError(f"{label} is empty: {path}")
    values: dict[str, str] = {}
    for line_number, line in enumerate(lines, start=1):
        if "=" not in line:
            raise ValueError(f"{label} has malformed line {line_number}")
        key, value = line.split("=", 1)
        _safe_text(key, f"{label} key {line_number}")
        if not _MANIFEST_KEY.fullmatch(key):
            raise ValueError(f"{label} has an invalid key: {key!r}")
        if key in values:
            raise ValueError(f"{label} duplicates key {key}")
        values[key] = _safe_text(
            value, f"{label} value {key}", allow_empty=True
        )
    return values


def _require_manifest_field(
    manifest: Mapping[str, str], key: str, label: str
) -> str:
    if key not in manifest:
        raise ValueError(f"{label} is missing {key}")
    value = _safe_text(manifest[key], f"{label} {key}")
    if not value:
        raise ValueError(f"{label} has an empty {key}")
    return value


def _validate_digest(value: str, label: str, pattern: re.Pattern[str]) -> str:
    if not pattern.fullmatch(value):
        raise ValueError(f"{label} is not a valid digest")
    return value.lower()


def _runtime_file_identity(path: Path, label: str) -> dict[str, Any]:
    identity = _file_identity(path, label)
    return {
        "path": identity["path"],
        "sha256": identity["sha256"],
        "size": identity["size"],
    }


def _validate_identity_manifests(
    source_root: Path,
    source_manifest_path: Path,
    runtime_identity_path: Path,
    run_root: Path,
    config_path: Path,
) -> dict[str, Any]:
    """Validate the immutable source/runtime binding and return provenance."""
    if source_root.name != "tree":
        raise ValueError("source root must end in /tree")
    if not re.fullmatch(r"[0-9A-Fa-f]{40}", source_root.parent.name):
        raise ValueError("source root is not commit keyed")

    expected_source_manifest = source_root.parent / "identity" / "source.manifest"
    if source_manifest_path != expected_source_manifest:
        raise ValueError("source manifest is not the immutable snapshot manifest")
    expected_runtime_identity = run_root / "manifests" / "runtime.identity"
    if runtime_identity_path != expected_runtime_identity:
        raise ValueError("runtime identity is not run-bound")

    run_source_manifest = run_root / "manifests" / "source.manifest"
    run_source_manifest = _existing_file(
        run_source_manifest, "run-bound source manifest"
    )
    if run_source_manifest.read_bytes() != source_manifest_path.read_bytes():
        raise ValueError("run-bound source manifest differs from immutable source manifest")

    source = _read_manifest(source_manifest_path, "source manifest")
    runtime = _read_manifest(runtime_identity_path, "runtime identity")
    expected_source_keys = {
        "FORMAT",
        "SOURCE_ROOT",
        "SOURCE_COMMIT",
        "SOURCE_ARCHIVE_PATH",
        "SOURCE_ARCHIVE_SHA256",
        "CONFIG_HELPER_SHA256",
        "DATASETS_SHA256",
        "PIXI_TOML_SHA256",
        "PIXI_LOCK_SHA256",
        "AUX_ROOT",
        "SCGATE_DB_BRANCH",
    }
    if set(source) != expected_source_keys:
        raise ValueError("source manifest must contain exactly its immutable identity fields")
    if _require_manifest_field(source, "FORMAT", "source manifest") != "1":
        raise ValueError("source manifest FORMAT must be 1")
    if _require_manifest_field(source, "SOURCE_ROOT", "source manifest") != str(source_root):
        raise ValueError("source manifest SOURCE_ROOT does not match source root")
    source_commit = _require_manifest_field(source, "SOURCE_COMMIT", "source manifest")
    if not re.fullmatch(r"[0-9A-Fa-f]{40}", source_commit):
        raise ValueError("source manifest SOURCE_COMMIT is not a full Git commit")
    if source_commit.casefold() != source_root.parent.name.casefold():
        raise ValueError("source manifest SOURCE_COMMIT does not match source root")

    expected_snapshot_root = source_root.parent
    expected_identity_dir = expected_snapshot_root / "identity"
    expected_source_archive = expected_identity_dir / "source.tar"
    source_archive_path = _safe_absolute(
        source["SOURCE_ARCHIVE_PATH"], "source archive path"
    )
    if source_archive_path != expected_source_archive:
        raise ValueError("source archive is not bound to the immutable snapshot identity")
    aux_root = _safe_absolute(source["AUX_ROOT"], "source aux root")
    if aux_root != expected_snapshot_root / "aux":
        raise ValueError("source aux root is not bound to the immutable snapshot")
    complete_marker = expected_snapshot_root / "COMPLETE"
    if not complete_marker.is_file() or complete_marker.is_symlink():
        raise ValueError("immutable snapshot completion marker is missing or unsafe")
    try:
        complete_value = complete_marker.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError) as exc:
        raise ValueError("immutable snapshot completion marker is unreadable") from exc
    if complete_value != "COMPLETE":
        raise ValueError("immutable snapshot completion marker is invalid")
    if (source_root / ".git").exists():
        raise ValueError("immutable source tree must not contain a Git directory")

    for key in (
        "SOURCE_ARCHIVE_PATH",
        "SOURCE_ARCHIVE_SHA256",
        "CONFIG_HELPER_SHA256",
        "DATASETS_SHA256",
        "PIXI_TOML_SHA256",
        "PIXI_LOCK_SHA256",
        "AUX_ROOT",
        "SCGATE_DB_BRANCH",
    ):
        _require_manifest_field(source, key, "source manifest")
    source_archive_sha = _validate_digest(
        source["SOURCE_ARCHIVE_SHA256"],
        "source manifest SOURCE_ARCHIVE_SHA256",
        _HEX_SHA256,
    )
    for key in (
        "CONFIG_HELPER_SHA256",
        "DATASETS_SHA256",
        "PIXI_TOML_SHA256",
        "PIXI_LOCK_SHA256",
    ):
        _validate_digest(source[key], f"source manifest {key}", _HEX_SHA256)

    source_archive = _existing_file(source_archive_path, "source archive")
    if _digest_file(source_archive, "sha256") != source_archive_sha:
        raise ValueError("source archive SHA-256 does not match source manifest")
    expected_files = {
        "CONFIG_HELPER_SHA256": source_root / "config_helper.R",
        "DATASETS_SHA256": source_root / "datasets.json",
        "PIXI_TOML_SHA256": source_root / "pixi.toml",
        "PIXI_LOCK_SHA256": source_root / "pixi.lock",
    }
    expected_config = _existing_file(expected_files["DATASETS_SHA256"], "immutable datasets.json")
    if config_path != expected_config:
        raise ValueError("configuration is not the immutable snapshot datasets.json")
    for field, bound_path in expected_files.items():
        bound_file = _existing_file(bound_path, f"immutable source {field}")
        if _digest_file(bound_file, "sha256") != source[field].lower():
            raise ValueError(f"immutable source {field} does not match source manifest")

    runtime_image = _existing_file(
        _require_manifest_field(runtime, "RUNTIME_IMAGE", "runtime identity"),
        "runtime image",
    )
    runtime_manifest = _existing_file(
        _require_manifest_field(runtime, "RUNTIME_MANIFEST", "runtime identity"),
        "runtime manifest",
    )
    runtime_image_sha = _validate_digest(
        _require_manifest_field(runtime, "RUNTIME_IMAGE_SHA256", "runtime identity"),
        "runtime identity RUNTIME_IMAGE_SHA256",
        _HEX_SHA256,
    )
    runtime_manifest_sha = _validate_digest(
        _require_manifest_field(runtime, "RUNTIME_MANIFEST_SHA256", "runtime identity"),
        "runtime identity RUNTIME_MANIFEST_SHA256",
        _HEX_SHA256,
    )
    runtime_image_size = _require_manifest_field(
        runtime, "RUNTIME_IMAGE_SIZE", "runtime identity"
    )
    runtime_manifest_size = _require_manifest_field(
        runtime, "RUNTIME_MANIFEST_SIZE", "runtime identity"
    )
    if not re.fullmatch(r"[1-9][0-9]*", runtime_image_size):
        raise ValueError("runtime identity RUNTIME_IMAGE_SIZE is invalid")
    if not re.fullmatch(r"[1-9][0-9]*", runtime_manifest_size):
        raise ValueError("runtime identity RUNTIME_MANIFEST_SIZE is invalid")
    runtime_image_info = _runtime_file_identity(runtime_image, "runtime image")
    runtime_manifest_info = _runtime_file_identity(runtime_manifest, "runtime manifest")
    if runtime_image_info["sha256"] != runtime_image_sha:
        raise ValueError("runtime image SHA-256 does not match runtime identity")
    if runtime_manifest_info["sha256"] != runtime_manifest_sha:
        raise ValueError("runtime manifest SHA-256 does not match runtime identity")
    if runtime_image_info["size"] != int(runtime_image_size):
        raise ValueError("runtime image size does not match runtime identity")
    if runtime_manifest_info["size"] != int(runtime_manifest_size):
        raise ValueError("runtime manifest size does not match runtime identity")
    runtime_manifest_fields = _read_manifest(runtime_manifest, "runtime manifest")
    runtime_format = _require_manifest_field(
        runtime_manifest_fields, "FORMAT", "runtime manifest"
    )
    if runtime_format not in {"1", "2"}:
        raise ValueError("runtime manifest FORMAT is unsupported")

    return {
        "source_root": str(source_root),
        "source_manifest": {"path": str(source_manifest_path), **source},
        "source_manifest_run": {"path": str(run_source_manifest)},
        "runtime_identity": {"path": str(runtime_identity_path), **runtime},
        "runtime_manifest": {
            "path": str(runtime_manifest),
            "sha256": runtime_manifest_info["sha256"],
            "size": runtime_manifest_info["size"],
            **runtime_manifest_fields,
        },
        "runtime_image": runtime_image_info,
    }


def _merge_columns(entry: Mapping[str, Any], view: Mapping[str, Any], dataset: str) -> dict[str, Any]:
    base = entry.get("columns")
    override = view.get("columns")
    if base is None:
        base = {}
    if override is None:
        override = {}
    if not isinstance(base, Mapping) or not isinstance(override, Mapping):
        raise ValueError(f"{dataset} columns must be objects")
    columns = dict(base)
    columns.update(override)
    return columns


def _resolve_config(
    config_path: Path,
    dataset: str,
    view_name: str,
    input_path: Path,
    expected_sample: str | None,
    expected_label: str | None,
    expected_batch_json: str | None,
    expected_subset_json: str | None,
) -> dict[str, Any]:
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read immutable configuration: {config_path}") from exc
    if not isinstance(config, Mapping):
        raise ValueError("datasets.json must be an object")
    _safe_text(dataset, "dataset")
    if not dataset or dataset.startswith("_"):
        raise ValueError(f"dataset is not a production dataset: {dataset}")
    entry = config.get(dataset)
    if not isinstance(entry, Mapping):
        raise ValueError(f"dataset is not configured: {dataset}")
    if entry.get("use_for_batch_effect") is not True:
        raise ValueError(f"dataset is not enabled for batch-effect processing: {dataset}")
    views = entry.get("views")
    if not isinstance(views, Mapping) or not isinstance(views.get(view_name), Mapping):
        raise ValueError(f"corrected view is not configured: {dataset}/{view_name}")
    if view_name != CORRECTED_VIEW:
        raise ValueError(f"audit requires {CORRECTED_VIEW} view")
    view = views[view_name]
    input_name = view.get("input_file_name")
    output_name = view.get("output_file_name")
    if input_name is None:
        input_name = view.get("input_file")
    if output_name is None:
        output_name = view.get("output_file")
    input_name = _safe_text(input_name, "configured input_file_name")
    output_name = _safe_text(output_name, "configured output_file_name")
    if input_path.name != Path(input_name).name:
        raise ValueError("source path basename does not match configured input_file_name")
    if "/" in output_name or output_name.startswith("."):
        raise ValueError(f"configured output_file_name is unsafe: {output_name}")

    columns = _merge_columns(entry, view, dataset)
    sample_column = _safe_text(columns.get("sample"), "configured sample column")
    label_column = _safe_text(columns.get("label"), "configured label column")
    if sample_column == label_column:
        raise ValueError("configured sample and label columns must differ")
    try:
        batch_keys = list(
            normalize_batch_keys(
                columns.get("batch"),
                sample_column=sample_column,
                biological_column=label_column,
            )
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{dataset} configured batch keys are invalid") from exc
    subset_vars = view.get("subset_vars")
    if subset_vars is None:
        subset_vars = {}
    if not isinstance(subset_vars, Mapping):
        raise ValueError("subset_vars must be an object")
    subset_vars = dict(subset_vars)

    if expected_sample is not None and expected_sample != sample_column:
        raise ValueError("configured sample column argument disagrees with immutable config")
    if expected_label is not None and expected_label != label_column:
        raise ValueError("configured label column argument disagrees with immutable config")
    if expected_batch_json is not None:
        try:
            expected_batch = json.loads(expected_batch_json)
        except json.JSONDecodeError as exc:
            raise ValueError("configured batch keys argument is not JSON") from exc
        try:
            expected_batch_keys = list(
                normalize_batch_keys(
                    expected_batch,
                    sample_column=sample_column,
                    biological_column=label_column,
                )
            )
        except (TypeError, ValueError) as exc:
            raise ValueError("configured batch keys argument is invalid") from exc
        if expected_batch_keys != batch_keys:
            raise ValueError("configured batch keys argument disagrees with immutable config")
    if expected_subset_json is not None:
        try:
            expected_subset = json.loads(expected_subset_json)
        except json.JSONDecodeError as exc:
            raise ValueError("configured subset argument is not JSON") from exc
        if expected_subset != subset_vars:
            raise ValueError("configured subset argument disagrees with immutable config")

    return {
        "config": config,
        "entry": entry,
        "view": view,
        "input_file_name": input_name,
        "output_file_name": output_name,
        "columns": columns,
        "sample_column": sample_column,
        "label_column": label_column,
        "batch_keys": batch_keys,
        "subset_vars": subset_vars,
    }


def _decode_index(values: Any, path: Path) -> pd.Index:
    index = np.asarray(values)
    if index.ndim != 1 or index.size == 0:
        raise ValueError(f"H5AD obs index is empty or malformed: {path}")
    index = index.astype(str)
    if any(not value.strip() for value in index):
        raise ValueError(f"H5AD obs index contains blank values: {path}")
    if len(set(index.tolist())) != len(index):
        raise ValueError(f"H5AD obs index contains duplicate values: {path}")
    return pd.Index(index)


def _read_obs_metadata(path: Path, columns: list[str]) -> pd.DataFrame:
    """Read exactly the requested obs vectors; no non-obs HDF5 group is touched."""
    with h5py.File(path, "r") as handle:
        obs = handle.get("obs")
        if obs is None or obs.attrs.get("encoding-type") not in ("dataframe", b"dataframe"):
            raise ValueError(f"H5AD obs is not a dataframe: {path}")
        if RESERVED_OBS_NAME in obs:
            raise ValueError(
                f"H5AD obs contains reserved temporary column {RESERVED_OBS_NAME!r}: {path}"
            )
        index_name_value = obs.attrs.get("_index", "_index")
        if isinstance(index_name_value, bytes):
            index_name_value = index_name_value.decode("utf-8")
        index_name = str(index_name_value)
        if index_name not in obs:
            raise ValueError(f"H5AD obs index is missing: {path}")
        obs_index = _decode_index(read_str_dataset(obs[index_name]), path)
        n_obs = len(obs_index)
        values: dict[str, Any] = {}
        for column in columns:
            if column == index_name:
                continue
            if column not in obs:
                raise ValueError(f"H5AD is missing requested obs column {column!r}: {path}")
            column_values = np.asarray(read_obs_column_values(obs, column), dtype=object)
            if column_values.ndim != 1 or len(column_values) != n_obs:
                raise ValueError(f"H5AD obs column {column!r} has the wrong length: {path}")
            values[column] = column_values
    frame = pd.DataFrame(values, index=obs_index)
    frame.index.name = index_name
    return frame


def _compact_validation(validation: Any, identity: Mapping[str, Any]) -> dict[str, Any]:
    keys = list(validation.keys)
    key_level_counts = {key: int(len(validation.levels[key])) for key in keys}
    key_near_unique_fraction = {
        key: float(len(validation.levels[key]) / validation.n_samples) for key in keys
    }
    return {
        "valid": bool(validation is not None),
        "sample_column": validation.sample_column,
        "biological_column": validation.biological_column,
        "batch_keys": keys,
        "n_cells": int(validation.n_obs),
        "n_samples": int(validation.n_samples),
        "key_level_counts": key_level_counts,
        "key_near_unique_fraction": key_near_unique_fraction,
        "near_unique_fraction": float(validation.near_unique_fraction),
        "estimable": bool(validation.estimable),
        "design_rank": int(validation.design_rank),
        "design_columns": int(validation.design_columns),
        "composite_design_rank": int(validation.composite_design_rank),
        "composite_design_columns": int(validation.composite_design_columns),
        "fingerprint": identity["fingerprint"],
        "method_id": CONTRACT_METHOD,
        "model_id": CONTRACT_MODEL,
    }


def _audit(
    *,
    config_path: Path,
    input_path: Path,
    output_path: Path,
    source_root: Path,
    source_manifest_path: Path,
    runtime_identity_path: Path,
    run_root: Path,
    dataset: str,
    view_name: str,
    expected_sample: str | None,
    expected_label: str | None,
    expected_batch_json: str | None,
    expected_subset_json: str | None,
) -> dict[str, Any]:
    provenance = _validate_identity_manifests(
        source_root,
        source_manifest_path,
        runtime_identity_path,
        run_root,
        config_path,
    )
    resolved = _resolve_config(
        config_path,
        dataset,
        view_name,
        input_path,
        expected_sample,
        expected_label,
        expected_batch_json,
        expected_subset_json,
    )
    sample_column = resolved["sample_column"]
    label_column = resolved["label_column"]
    batch_keys = resolved["batch_keys"]
    subset_vars = resolved["subset_vars"]
    required_columns = list(
        dict.fromkeys([sample_column, label_column, *batch_keys, *subset_vars.keys()])
    )
    metadata = _read_obs_metadata(input_path, required_columns)
    adata = _ObsOnly(metadata)

    # These calls are intentionally the shared production implementations.  In
    # particular, no comparison operator is reproduced in this source audit.
    subset_mask = evaluate_subset_mask(adata, subset_vars)
    subset_audit = assert_subset_sample_consistency(
        adata,
        subset_mask,
        sample_column,
        context=f"{dataset}/{view_name} corrected source subset",
    )
    if subset_audit["retained_cells"] < 1:
        raise ValueError("corrected source subset retained no cells")
    retained_metadata = metadata.loc[subset_mask]
    validation_columns = list(dict.fromkeys([sample_column, *batch_keys, label_column]))
    validation = validate_batch_metadata(
        retained_metadata[validation_columns],
        batch_keys,
        sample_column=sample_column,
        biological_column=label_column,
    )
    contract_identity = build_batch_contract_identity(
        batch_keys,
        sample_column=sample_column,
        method_id=CONTRACT_METHOD,
        model_id=CONTRACT_MODEL,
    )
    validation_summary = _compact_validation(validation, contract_identity)
    if not validation_summary["valid"] or not validation_summary["estimable"]:
        raise ValueError("corrected source metadata validation did not produce a valid estimable contract")

    # Keep the configured rule in the subset section as well as config: this
    # makes the actual audited predicate explicit to shell-side evidence checks.
    subset_audit = {
        **subset_audit,
        "configured": subset_vars,
        "validation_scope": "full_source_metadata",
    }
    source_identity = _file_identity(input_path, "H5AD source")
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "status": "SOURCE_METADATA_VALIDATED_H5AD",
        "dataset": dataset,
        "view": view_name,
        "source_type": "h5ad",
        "obs_only": True,
        "read_scope": ["obs"],
        "source_path": source_identity["path"],
        "source_identity": source_identity,
        "config": {
            "dataset": dataset,
            "view": view_name,
            "input_file_name": resolved["input_file_name"],
            "output_file_name": resolved["output_file_name"],
            "sample_column": sample_column,
            "label_column": label_column,
            "batch_keys": batch_keys,
            "subset_vars": subset_vars,
        },
        "subset_audit": subset_audit,
        "validation_summary": validation_summary,
        "provenance": provenance,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def _write_json_atomic(payload: Mapping[str, Any], output_path: Path) -> None:
    temporary = None
    try:
        fd, temporary = tempfile.mkstemp(
            prefix=f".{output_path.name}.build.",
            dir=str(output_path.parent),
            text=True,
        )
        os.close(fd)
        temporary_path = Path(temporary)
        with temporary_path.open("w", encoding="utf-8") as handle:
            json.dump(_json_safe(payload), handle, sort_keys=True, indent=2, allow_nan=False)
            handle.write("\n")
        os.chmod(temporary_path, 0o600)
        os.replace(temporary_path, output_path)
        temporary = None
    finally:
        if temporary is not None:
            try:
                Path(temporary).unlink()
            except FileNotFoundError:
                pass


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", "--config-path", dest="config", required=True)
    parser.add_argument("--input-file", "--input", dest="input_file", required=True)
    parser.add_argument("--output", "--output-file", dest="output", required=True)
    parser.add_argument("--dataset", "--ds-name", dest="dataset", required=True)
    parser.add_argument("--view", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--source-manifest", required=True)
    parser.add_argument("--runtime-identity", required=True)
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--sample-column")
    parser.add_argument("--label-column")
    parser.add_argument("--batch-keys-json")
    parser.add_argument("--subset-vars-json")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        source_root = _existing_directory(args.source_root, "source root")
        run_root = _existing_directory(args.run_root, "run root")
        config_path = _existing_file(args.config, "configuration")
        input_path = _existing_file(args.input_file, "H5AD source")
        output_path = _run_owned_output(args.output, run_root)
        source_manifest_path = _existing_file(args.source_manifest, "source manifest")
        runtime_identity_path = _existing_file(args.runtime_identity, "runtime identity")
        if input_path.suffix.casefold() != ".h5ad":
            raise ValueError("H5AD source path must end in .h5ad")
        if not _RUN_ID.fullmatch(run_root.name):
            raise ValueError("run root basename is not a valid run ID")
        report = _audit(
            config_path=config_path,
            input_path=input_path,
            output_path=output_path,
            source_root=source_root,
            source_manifest_path=source_manifest_path,
            runtime_identity_path=runtime_identity_path,
            run_root=run_root,
            dataset=args.dataset,
            view_name=args.view,
            expected_sample=args.sample_column,
            expected_label=args.label_column,
            expected_batch_json=args.batch_keys_json,
            expected_subset_json=args.subset_vars_json,
        )
        _write_json_atomic(report, output_path)
    except (OSError, KeyError, TypeError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=os.sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
