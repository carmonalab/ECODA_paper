#!/usr/bin/env python3
"""Export one sample-level Feather table by reading only an H5AD ``obs`` group.

``X``, ``raw``, and ``layers['counts']`` are never opened. It retains the first
observation for each configured sample ID in source order by default. An
explicit ``majority_v1`` policy replaces only its configured technical winners
with unique sample-level winners; every other metadata field remains first-row.
Batch-effect exports require the configured high-resolution cell-type column
when a cell-type column is needed; low-resolution annotations are not inputs to
batch-effect calculations and are carried only when present.
It writes the output and its strict MD5 sidecar atomically, and never publishes
an H5AD artifact record.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import os
from pathlib import Path
from typing import Any, Hashable, Iterable, Mapping

import h5py
import numpy as np
import pandas as pd

try:
    from h5ad_source_identity import read_obs_column_values
except ImportError:  # package import in focused tests
    from .h5ad_source_identity import read_obs_column_values


READ_CHUNK_SIZE = 1024 * 1024


def _decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _obs_column_length(node: Any) -> int:
    shape = getattr(node, "shape", None)
    if shape is None and hasattr(node, "keys") and "codes" in node:
        shape = node["codes"].shape
    if shape is None or len(shape) != 1:
        raise ValueError(f"H5AD obs column {node.name} is not one-dimensional")
    return int(shape[0])


def _normalise_sample(value: Any) -> str:
    value = _decode(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        raise ValueError("H5AD sample metadata contains a missing sample ID")
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        raise ValueError("H5AD sample metadata contains a non-finite sample ID")
    text = str(value).strip()
    if not text or text.casefold() in {"nan", "none", "<na>"}:
        raise ValueError("H5AD sample metadata contains a blank sample ID")
    return text


def _feather_scalar(value: Any) -> Any:
    value = _decode(value)
    if value is None:
        return None
    if pd.isna(value) if not isinstance(value, (list, tuple, dict, np.ndarray)) else False:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _config_entry(config_path: Path, dataset: str, view: str) -> tuple[dict, dict]:
    import json

    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"cannot read datasets config: {config_path}") from exc
    entry = config.get(dataset)
    if not isinstance(entry, dict):
        raise ValueError(f"dataset is missing from config: {dataset}")
    views = entry.get("views")
    if not isinstance(views, dict) or view not in views:
        raise ValueError(f"dataset view is missing from config: {dataset}/{view}")
    view_config = views[view]
    if not isinstance(view_config, dict):
        raise ValueError(f"dataset view is malformed: {dataset}/{view}")
    return entry, view_config


def _dataset_spec_columns(config_path: Path, dataset: str) -> list[str]:
    """Return optional candidate columns from the authoritative registry module."""
    module_path = config_path.parent / "notebooks" / "dataset_onboarding" / "dataset_specs.py"
    if not module_path.is_file():
        return []
    spec = importlib.util.spec_from_file_location("ecoda_dataset_specs", module_path)
    if spec is None or spec.loader is None:
        return []
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        dataset_specs = getattr(module, "DATASET_SPECS", {})
        batch_effect_specs = getattr(module, "BATCH_EFFECT_SPECS", {})
        dataset_spec = dataset_specs.get(dataset, {})
        batch_candidates = batch_effect_specs.get(dataset, [])
    except (AttributeError, ImportError, OSError, TypeError, ValueError):
        return []

    columns: list[str] = []

    def add(value: Any) -> None:
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, (list, tuple)):
            values = value
        else:
            return
        columns.extend(
            item for item in values
            if isinstance(item, str) and item.strip()
        )

    # These fields are the complete candidate registry in DATASET_SPECS.  The
    # biological column is included for technical-only/legacy entries whose
    # configured columns may not repeat it.
    if isinstance(dataset_spec, dict):
        add(dataset_spec.get("bio_col"))
        for key in (
            "sample_candidates",
            "sample_stable_cols",
            "batch_cols",
            "cell_type_candidates",
        ):
            add(dataset_spec.get(key, []))
    # Joanito and Stephenson are represented in the technical batch registry
    # even when they have no active DATASET_SPECS entry (notably Site).
    add(batch_candidates)
    return list(dict.fromkeys(columns))

def _dataset_columns(entry: dict) -> dict:
    """Return the authoritative dataset-level metadata column contract."""
    columns = entry.get("columns")
    if not isinstance(columns, dict):
        raise ValueError("dataset columns are malformed")
    return columns


def requested_columns(
    config_path: Path,
    dataset: str,
    entry: dict,
    view: dict | None = None,
    *,
    batch_view: bool = False,
) -> tuple[str, list[str], list[str]]:
    columns = _dataset_columns(entry)
    raw_sample_column = columns.get("sample")
    label_column = columns.get("label")
    if not isinstance(raw_sample_column, str) or not raw_sample_column.strip():
        raise ValueError(f"dataset sample column is missing: {dataset}")
    if not isinstance(label_column, str) or not label_column.strip():
        raise ValueError(f"dataset primary label column is missing: {dataset}")
    # Stage 3 persists the standardized Sample column. It is the only valid
    # grouping/order key for final Stage 5; the configured raw identity remains
    # an optional source metadata column when it is still present.
    required: list[str] = ["Sample", label_column]
    required_keys = ["batch", "cell_type_low_res", "cell_type_high_res"]
    if batch_view:
        # Batch-effect methods use only the configured high-resolution
        # annotation. The low-resolution column remains optional carry-through.
        required_keys.remove("cell_type_low_res")
    for key in required_keys:
        value = columns.get(key)
        if isinstance(value, str) and value.strip():
            required.append(value)
        elif isinstance(value, (list, tuple)):
            required.extend(
                str(item) for item in value if isinstance(item, str) and item.strip()
            )
    required = list(dict.fromkeys(required))
    optional_configured = []
    if batch_view:
        value = columns.get("cell_type_low_res")
        if isinstance(value, str) and value.strip():
            optional_configured.append(value)
        elif isinstance(value, (list, tuple)):
            optional_configured.extend(
                str(item) for item in value if isinstance(item, str) and item.strip()
            )
    optional = list(
        dict.fromkeys(
            [
                *optional_configured,
                raw_sample_column,
                *_dataset_spec_columns(config_path, dataset),
            ]
        )
    )
    optional = [column for column in optional if column not in required]
    all_columns = required + optional
    return "Sample", required, all_columns

def _ordered_metadata_columns(
    config_path: Path,
    dataset: str,
    entry: dict,
    view: dict,
    frame_columns: Iterable[str],
) -> list[str]:
    """Preserve the configured metadata order while omitting absent optional fields."""
    columns = _dataset_columns(entry)
    available = [str(column) for column in frame_columns]
    ordered: list[str] = []
    seen: set[str] = set()

    def add(value: Any) -> None:
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, (list, tuple)):
            values = value
        else:
            return
        for column in values:
            if isinstance(column, str) and column in available and column not in seen:
                seen.add(column)
                ordered.append(column)

    add("Sample")
    for key in ("label", "batch", "cell_type_low_res", "cell_type_high_res", "sample"):
        add(columns.get(key))
    add(_dataset_spec_columns(config_path, dataset))
    for column in available:
        add(column)
    return ordered

_APPROVED_MAJORITY_DATASETS = frozenset(
    {"Alzheimer", "Breast_cancer", "Lupus_PBMC"}
)
_APPROVED_MAJORITY_KEYS = {
    "Alzheimer": ("assay",),
    "Breast_cancer": ("suspension_dissociation_time",),
    "Lupus_PBMC": ("batch_cov",),
}
_BREAST_UNKNOWN_BATCH_KEY = "suspension_dissociation_time"
_MISSING_BATCH_TEXTS = frozenset(
    {
        "na",
        "nan",
        "none",
        "<na>",
        "n/a",
        "null",
        "unknown",
        "inf",
        "+inf",
        "-inf",
        "infinity",
        "+infinity",
        "-infinity",
    }
)


def _normalise_batch_keys(
    raw_batch_keys: Any,
    *,
    sample_column: str = "Sample",
    biological_column: str | None = None,
) -> list[str]:
    if raw_batch_keys is None:
        return []
    if isinstance(raw_batch_keys, str):
        values = [raw_batch_keys]
    elif isinstance(raw_batch_keys, (bytes, bytearray, Mapping, set, frozenset)):
        raise ValueError("configured batch keys are invalid")
    else:
        try:
            values = list(raw_batch_keys)
        except TypeError as exc:
            raise ValueError("configured batch keys are invalid") from exc
    if not isinstance(sample_column, str) or not sample_column.strip():
        raise ValueError("sample column is invalid")
    if biological_column is not None and (
        not isinstance(biological_column, str) or not biological_column.strip()
    ):
        raise ValueError("biological label column is invalid")
    result: list[str] = []
    seen: set[str] = set()
    for key in values:
        if not isinstance(key, str) or not key.strip():
            raise ValueError("configured batch keys must be nonblank strings")
        if key in seen:
            raise ValueError(f"configured batch keys contain duplicate {key!r}")
        if key == sample_column:
            raise ValueError("configured batch keys must not contain the sample column")
        if biological_column is not None and key == biological_column:
            raise ValueError("configured batch keys must not contain the biological label")
        seen.add(key)
        result.append(key)
    return result


def _normalise_batch_metadata_policy(
    raw_policy: Any,
    *,
    dataset: str | None = None,
) -> dict[str, Any] | None:
    if raw_policy is None:
        return None
    if not isinstance(raw_policy, Mapping):
        raise ValueError("batch metadata policy is invalid")
    allowed_fields = {
        "sample_aggregation",
        "majority_keys",
        "accepted_sentinel_values",
    }
    unknown_fields = set(raw_policy) - allowed_fields
    if unknown_fields:
        raise ValueError(
            f"batch metadata policy has unsupported fields: {sorted(unknown_fields)!r}"
        )
    if raw_policy.get("sample_aggregation") != "majority_v1":
        raise ValueError("only sample_aggregation=majority_v1 is supported")
    if dataset is not None and dataset not in _APPROVED_MAJORITY_DATASETS:
        raise ValueError(
            f"majority_v1 batch metadata policy is not approved for {dataset!r}"
        )

    raw_majority_keys = raw_policy.get("majority_keys")
    if (
        isinstance(raw_majority_keys, (str, bytes, bytearray))
        or not isinstance(raw_majority_keys, (list, tuple))
        or not raw_majority_keys
    ):
        raise ValueError("majority_keys must be a non-empty list")
    majority_keys = _normalise_batch_keys(raw_majority_keys)
    if dataset is not None and tuple(majority_keys) != _APPROVED_MAJORITY_KEYS[dataset]:
        raise ValueError(
            f"majority_keys are not approved for {dataset!r}: {majority_keys!r}"
        )

    raw_accepted = raw_policy.get("accepted_sentinel_values", {})
    if not isinstance(raw_accepted, Mapping):
        raise ValueError("accepted_sentinel_values must be a mapping")
    accepted: dict[str, tuple[str, ...]] = {}
    for key, values in raw_accepted.items():
        if not isinstance(key, str) or not key.strip():
            raise ValueError("accepted sentinel batch keys must be nonblank strings")
        if isinstance(values, str) or not isinstance(values, (list, tuple)):
            raise ValueError(f"accepted sentinel values for {key!r} are invalid")
        values_list = list(values)
        if key != _BREAST_UNKNOWN_BATCH_KEY or values_list != ["unknown"]:
            raise ValueError(
                "only literal 'unknown' for "
                f"{_BREAST_UNKNOWN_BATCH_KEY!r} may be accepted"
            )
        if key not in majority_keys:
            raise ValueError(
                f"accepted sentinel key {key!r} must be a majority key"
            )
        accepted[key] = ("unknown",)
    if dataset in {"Alzheimer", "Lupus_PBMC"} and accepted:
        raise ValueError(f"accepted sentinel values are not approved for {dataset!r}")
    return {
        "sample_aggregation": "majority_v1",
        "majority_keys": tuple(majority_keys),
        "accepted_sentinel_values": accepted,
    }


def _resolved_batch_configuration(
    entry: dict,
    view: dict,
    label_column: str,
    dataset: str | None = None,
) -> tuple[list[str], dict[str, Any] | None]:
    """Resolve configured technical keys and the authoritative dataset policy."""
    columns = _dataset_columns(entry)
    batch_keys = _normalise_batch_keys(
        columns.get("batch"),
        sample_column="Sample",
        biological_column=label_column,
    )
    # Batch aggregation is a dataset-level contract.  Do not allow a
    # view-local object to replace the authoritative datasets.json policy.
    policy = _normalise_batch_metadata_policy(
        entry.get("batch_metadata_policy"),
        dataset=dataset,
    )
    if policy is not None:
        majority_keys = _normalise_batch_keys(
            policy["majority_keys"],
            sample_column="Sample",
            biological_column=label_column,
        )
        if any(key not in batch_keys for key in majority_keys):
            raise ValueError("majority_keys must be configured batch keys")
        policy = {**policy, "majority_keys": tuple(majority_keys)}
    return batch_keys, policy


def _validate_variant_output_binding(
    output: Path,
    dataset: str,
    view: str,
    variant: str | None = None,
) -> None:
    """Reject an explicit Stage 5 output that escapes its variant lane."""

    environment_variant = os.environ.get("ANALYSIS_VARIANT", "")
    if (
        variant is not None
        and environment_variant
        and str(variant) != environment_variant
    ):
        raise ValueError(
            "metadata export variant disagrees with ANALYSIS_VARIANT"
        )
    selected_variant = (
        environment_variant if variant is None else variant
    )
    selected_variant = str(selected_variant or "")
    if selected_variant not in ("", "final", "corrected_final"):
        raise ValueError(f"unsupported Stage 5 analysis variant: {selected_variant}")
    if not selected_variant:
        return
    expected_pass = (
        "uncorrected" if selected_variant == "final" else "corrected"
    )
    configured_pass = os.environ.get("ANALYSIS_PASS", "")
    if configured_pass and configured_pass != expected_pass:
        raise ValueError(
            f"{selected_variant} metadata export requires the {expected_pass} pass"
        )
    expected_view = f"batch_effect_{expected_pass}"
    if view != expected_view:
        raise ValueError(
            f"{selected_variant} metadata export requires view {expected_view}"
        )
    analysis_root = os.environ.get("ANALYSIS_ROOT", "")
    if not analysis_root or not Path(analysis_root).is_absolute():
        raise ValueError(
            f"{selected_variant} metadata export requires an absolute ANALYSIS_ROOT"
        )
    root_text = analysis_root.rstrip("/")
    expected_suffix = f"/batch_effect/{expected_pass}_final"
    if not root_text.endswith(expected_suffix):
        raise ValueError(
            f"{selected_variant} metadata export is not bound to "
            "the variant-qualified output root"
        )
    expected = Path(root_text) / "metadata" / f"{dataset}_sample_metadata.feather"
    if output != expected:
        raise ValueError(
            "metadata output is not the exact variant-qualified path: "
            f"{output} (expected {expected})"
        )


def _canonical_batch_value(
    value: Any,
    *,
    key: str,
    accepted_sentinel_values: Mapping[str, tuple[str, ...]],
) -> tuple[Hashable, Any]:
    """Validate one technical class and return a hashable counting key."""
    value = _decode(value)
    if value is None or value is pd.NA or value is pd.NaT:
        raise ValueError(f"batch metadata column {key!r} contains a missing value")
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, (bool, np.bool_)) and bool(missing):
        raise ValueError(f"batch metadata column {key!r} contains a missing value")

    if isinstance(value, str):
        if value == "unknown" and key in accepted_sentinel_values:
            return ("str", value), value
        if not value.strip() or value.strip().casefold() in _MISSING_BATCH_TEXTS:
            raise ValueError(
                f"batch metadata column {key!r} contains a missing/blank value"
            )
        return ("str", value), value
    if isinstance(value, (bool, np.bool_)):
        return ("bool", bool(value)), bool(value)
    if isinstance(value, (int, np.integer)):
        return ("int", int(value)), int(value)
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(
                f"batch metadata column {key!r} contains a non-finite value"
            )
        return ("float", numeric), value
    try:
        hash(value)
    except TypeError as exc:
        raise ValueError(
            f"batch metadata column {key!r} contains an unsupported value"
        ) from exc
    return (
        type(value).__module__,
        type(value).__qualname__,
        value,
    ), value


def _majority_batch_winners(
    sample_ids: Iterable[str],
    batch_keys: Iterable[str],
    batch_counts: Mapping[str, Mapping[str, Mapping[Hashable, int]]],
    batch_first_values: Mapping[
        str, Mapping[str, Mapping[Hashable, Any]]
    ],
) -> dict[str, list[Any]]:
    """Reduce streamed per-Sample technical counts to unique class winners."""
    keys = list(batch_keys)
    winners: dict[str, list[Any]] = {key: [] for key in keys}
    for sample in sample_ids:
        sample_counts = batch_counts.get(sample)
        sample_first_values = batch_first_values.get(sample)
        if sample_counts is None or sample_first_values is None:
            raise ValueError(f"batch metadata counts are missing for Sample {sample!r}")
        for key in keys:
            counts = sample_counts.get(key)
            first_values = sample_first_values.get(key)
            if not counts or first_values is None:
                raise ValueError(
                    f"batch metadata counts are missing for {key!r} in Sample {sample!r}"
                )
            top_count = max(counts.values())
            top_classes = [
                category for category, count in counts.items() if count == top_count
            ]
            if len(top_classes) != 1:
                raise ValueError(
                    f"batch key {key!r} has no unique majority within Sample {sample!r}"
                )
            category = top_classes[0]
            if category not in first_values:
                raise ValueError(
                    f"batch metadata first value is missing for {key!r} "
                    f"in Sample {sample!r}"
                )
            winners[key].append(first_values[category])
    return winners


def read_obs_metadata(
    input_path: Path,
    sample_column: str,
    required_columns: Iterable[str],
    candidate_columns: Iterable[str],
    chunk_size: int = READ_CHUNK_SIZE,
    *,
    batch_keys: Iterable[str] | None = None,
    batch_metadata_policy: Mapping[str, Any] | None = None,
    biological_column: str | None = None,
) -> pd.DataFrame:
    if not input_path.is_file() or input_path.stat().st_size <= 0:
        raise ValueError(f"H5AD input is missing or empty: {input_path}")
    if chunk_size <= 0:
        raise ValueError("chunk size must be positive")
    if not isinstance(sample_column, str) or not sample_column.strip():
        raise ValueError("sample column is invalid")
    required = list(dict.fromkeys(str(column) for column in required_columns))
    configured_batch_keys = _normalise_batch_keys(
        batch_keys,
        sample_column=sample_column,
        biological_column=biological_column,
    )
    policy = _normalise_batch_metadata_policy(batch_metadata_policy)
    majority_keys: list[str] = []
    if policy is not None:
        majority_keys = _normalise_batch_keys(
            policy["majority_keys"],
            sample_column=sample_column,
            biological_column=biological_column,
        )
        if any(key not in configured_batch_keys for key in majority_keys):
            raise ValueError("majority_keys must be configured batch keys")
        accepted = policy["accepted_sentinel_values"]
        if any(key not in majority_keys for key in accepted):
            raise ValueError("accepted sentinel keys must be majority keys")
    use_majority = bool(majority_keys)
    with h5py.File(input_path, "r") as handle:
        if "obs" not in handle:
            raise ValueError(f"H5AD lacks obs dataframe: {input_path}")
        obs = handle["obs"]
        if obs.attrs.get("encoding-type", "") not in ("dataframe", b"dataframe"):
            raise ValueError(f"H5AD obs is not a dataframe: {input_path}")
        index_name = _decode(obs.attrs.get("_index", "_index"))
        if index_name not in obs:
            raise ValueError(f"H5AD obs index is missing: {input_path}")
        n_obs = _obs_column_length(obs[index_name])
        if n_obs <= 0:
            raise ValueError(f"H5AD obs is empty: {input_path}")
        needed = list(dict.fromkeys([sample_column, *required, *configured_batch_keys]))
        missing = [column for column in needed if column not in obs]
        if missing:
            raise ValueError(f"H5AD is missing configured obs columns: {missing}")
        if isinstance(candidate_columns, str):
            candidate_values: Iterable[str] = (candidate_columns,)
        else:
            candidate_values = candidate_columns
        present_optional = [
            column
            for column in candidate_values
            if column in obs and column not in required and column not in configured_batch_keys
        ]
        columns = list(
            dict.fromkeys([sample_column, *required, *configured_batch_keys, *present_optional])
        )
        values: dict[str, list[Any]] = {column: [] for column in columns}
        seen: set[str] = set()
        sample_ids: list[str] = []
        batch_counts: dict[str, dict[str, dict[Hashable, int]]] = {}
        batch_first_values: dict[
            str, dict[str, dict[Hashable, Any]]
        ] = {}
        accepted = (
            policy["accepted_sentinel_values"] if policy is not None else {}
        )
        for start in range(0, n_obs, chunk_size):
            stop = min(n_obs, start + chunk_size)
            chunk = {
                column: read_obs_column_values(obs, column, start, stop)
                for column in columns
            }
            for offset, raw_sample in enumerate(chunk[sample_column]):
                sample = _normalise_sample(raw_sample)
                if sample not in seen:
                    seen.add(sample)
                    sample_ids.append(sample)
                    for column in columns:
                        value = sample if column == sample_column else chunk[column][offset]
                        values[column].append(_feather_scalar(value))
                    if use_majority:
                        batch_counts[sample] = {
                            key: {} for key in majority_keys
                        }
                        batch_first_values[sample] = {
                            key: {} for key in majority_keys
                        }
                if use_majority:
                    for key in majority_keys:
                        category, raw_value = _canonical_batch_value(
                            chunk[key][offset],
                            key=key,
                            accepted_sentinel_values=accepted,
                        )
                        if key not in batch_counts[sample]:
                            continue
                        counts = batch_counts[sample][key]
                        first_values = batch_first_values[sample][key]
                        counts[category] = counts.get(category, 0) + 1
                        first_values.setdefault(category, raw_value)
        if not sample_ids:
            raise ValueError(f"H5AD has no non-empty sample IDs: {input_path}")
    frame = pd.DataFrame(values)
    if sample_column != "Sample":
        frame.insert(0, "Sample", frame[sample_column])
    else:
        frame.insert(0, "Sample", frame.pop("Sample"))
    if frame["Sample"].duplicated().any() or frame["Sample"].astype(str).str.strip().eq("").any():
        raise ValueError("exported sample metadata has blank or duplicate Sample IDs")
    if use_majority:
        winners = _majority_batch_winners(
            sample_ids,
            majority_keys,
            batch_counts,
            batch_first_values,
        )
        for key in majority_keys:
            frame[key] = [_feather_scalar(value) for value in winners[key]]
    return frame


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_sidecar(path: Path) -> dict[str, str]:
    sidecar = Path(f"{path}.md5")
    if path.is_symlink() or sidecar.is_symlink():
        raise ValueError(f"metadata checksum path is a symlink: {path}")
    try:
        lines = sidecar.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}") from exc
    if len(lines) != 3:
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}")
    fields: dict[str, str] = {}
    for key, line in zip(("MD5", "SIZE", "PATH"), lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in fields:
            raise ValueError(f"invalid metadata checksum sidecar: {sidecar}")
        fields[key] = line[len(prefix):]
    digest = fields["MD5"]
    if (
        len(digest) != 32
        or digest != digest.lower()
        or any(char not in "0123456789abcdef" for char in digest)
    ):
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}")
    try:
        size = int(fields["SIZE"])
    except ValueError as exc:
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}") from exc
    if size <= 0 or str(size) != fields["SIZE"]:
        raise ValueError(f"invalid metadata checksum sidecar: {sidecar}")
    if fields["PATH"] != str(path) or size != path.stat().st_size:
        raise ValueError(f"metadata checksum sidecar path/size mismatch: {sidecar}")
    if digest != _md5(path):
        raise ValueError(f"metadata checksum sidecar digest mismatch: {sidecar}")
    return fields


def validate_output(path: Path, required_columns: Iterable[str]) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"metadata Feather is missing or empty: {path}")
    _read_sidecar(path)
    try:
        frame = pd.read_feather(path)
    except Exception as exc:  # pyarrow surfaces several exception classes
        raise ValueError(f"metadata Feather is unreadable: {path}") from exc
    required = list(dict.fromkeys(str(column) for column in required_columns))
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"metadata Feather is missing configured columns: {missing}")
    if frame.empty or "Sample" not in frame.columns:
        raise ValueError(f"metadata Feather has no Sample rows: {path}")
    ids = frame["Sample"].astype(str)
    if ids.str.strip().eq("").any() or ids.duplicated().any():
        raise ValueError(f"metadata Feather has blank or duplicate Sample IDs: {path}")


def write_metadata(frame: pd.DataFrame, output: Path) -> None:
    if not output.is_absolute():
        raise ValueError("metadata output path must be absolute")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{output}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    backup = output.with_name(f".{output.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    if output.is_symlink() or sidecar.is_symlink():
        raise ValueError(f"metadata output path is a symlink: {output}")
    if output.exists() and not output.is_file():
        raise ValueError(f"metadata output is not a regular file: {output}")
    if sidecar.exists() and not sidecar.is_file():
        raise ValueError(f"metadata checksum is not a regular file: {sidecar}")
    had_output = output.is_file()
    had_sidecar = sidecar.is_file()
    try:
        frame.reset_index(drop=True).to_feather(temporary)
        if not temporary.is_file() or temporary.stat().st_size <= 0:
            raise ValueError(f"metadata Feather is empty: {temporary}")
        if had_output:
            os.link(output, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(temporary, output)
        digest = _md5(output)
        sidecar_tmp.write_text(
            f"MD5={digest}\nSIZE={output.stat().st_size}\nPATH={output}\n",
            encoding="utf-8",
        )
        os.replace(sidecar_tmp, sidecar)
        validate_output(output, frame.columns)
    except Exception:
        if backup.exists():
            os.replace(backup, output)
        elif not had_output and output.exists():
            output.unlink()
        if sidecar_backup.exists():
            os.replace(sidecar_backup, sidecar)
        elif not had_sidecar and sidecar.exists():
            sidecar.unlink()
        raise
    finally:
        for path in (temporary, sidecar_tmp, backup, sidecar_backup):
            try:
                path.unlink()
            except FileNotFoundError:
                pass


def export(args: argparse.Namespace) -> None:
    config_path = args.config.resolve()
    entry, view = _config_entry(config_path, args.dataset, args.view)
    selected_variant = (
        getattr(args, "analysis_variant", None)
        or os.environ.get("ANALYSIS_VARIANT", "")
    )
    sample_column, required, candidates = requested_columns(
        config_path,
        args.dataset,
        entry,
        view,
        batch_view=args.view
        in {"batch_effect_uncorrected", "batch_effect_corrected"},
    )
    columns = _dataset_columns(entry)
    label_column = columns["label"]
    batch_keys, policy = _resolved_batch_configuration(
        entry,
        view,
        label_column,
        dataset=args.dataset,
    )
    output = args.output.resolve()
    _validate_variant_output_binding(
        output,
        args.dataset,
        args.view,
        getattr(args, "analysis_variant", None),
    )
    if args.view == "batch_effect_uncorrected" and selected_variant == "final":
        # The explicit uncorrected-final lane retains the first source row for
        # every field; corrected-final is the only final lane using the
        # approved majority_v1 policy.
        policy = None
    if args.check:
        validate_output(output, required)
        return
    if output.exists() or Path(f"{output}.md5").exists():
        try:
            validate_output(output, required)
            print(f"METADATA_EXPORT=NOOP_VALIDATED PATH={output}")
            return
        except (OSError, ValueError):
            # Rebuild only an invalid output; the source H5AD remains read-only.
            pass
    frame = read_obs_metadata(
        args.input_file.resolve(),
        sample_column,
        required,
        candidates,
        args.chunk_size,
        batch_keys=batch_keys,
        batch_metadata_policy=policy,
        biological_column=label_column,
    )
    frame = frame.loc[
        :,
        _ordered_metadata_columns(
            config_path,
            args.dataset,
            entry,
            view,
            frame.columns,
        ),
    ]
    write_metadata(frame, output)
    print(f"METADATA_EXPORT=OK PATH={output} SAMPLES={len(frame)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--view", required=True)
    parser.add_argument(
        "--analysis-variant",
        default=None,
        choices=["final", "corrected_final"],
    )
    parser.add_argument("--input-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=READ_CHUNK_SIZE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        export(args)
    except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
