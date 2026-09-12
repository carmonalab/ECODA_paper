"""Deterministic corrected-mode batch metadata and composite contract.

The corrected pipeline has several consumers that cannot all accept the same
batch representation.  This module is the small, dependency-shared boundary:
configuration is normalized here, cell metadata is validated before any
sample-level reduction, and scalar composite values use one byte-exact
cross-language encoding.

Nothing in this module writes an AnnData object or mutates the caller's
metadata.  :func:`build_batch_composite` returns a copied, in-memory frame for
callers that need a temporary scalar column; the reserved column is never
persisted by this module.
"""

from __future__ import annotations

import datetime as _datetime
import json
import hashlib
import math
import struct
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


TOKEN_VERSION = "ecoda_batch_composite_v1"
FINGERPRINT_VERSION = "ecoda_batch_contract_v1"
RESERVED_OBS_NAME = "__ecoda_batch_combined_v1"
DIRECT_SCALARIZATION = "direct_v1"
COMPOSITE_SCALARIZATION = "composite_v1"

# Compact validated summaries are deliberately versioned independently from
# the source/configuration identity.  The identity fingerprint remains
# configuration-only; this field records only the bounded result of validating
# the full cell metadata.
VALIDATION_SUMMARY_SCHEMA_VERSION = 1
PREPROCESS_CORRECTION_MODE = "preprocess_hvg_harmony"
NATIVE_HARMONY_CORRECTION_MODE = "native_harmony_embedding"
MRVI_CORRECTION_MODE = "mrvi_composite_batch"
PSEUDOBULK_CORRECTION_MODE = "batch_only_pseudobulk"
ECODA_CORRECTION_MODE = "additive_random_intercepts"
_VALIDATION_SUMMARY_FIELDS = frozenset(
    {
        "schema_version",
        "validated_before_reduction",
        "sample_constancy",
        "per_key_levels",
        "key_level_counts",
        "composite_levels",
        "composite_level_count",
        "n_cells",
        "n_samples",
        "correction_mode",
        "correction_formula",
    }
)

# These are compared after stripping surrounding whitespace and case-folding.
# The original non-sentinel text is otherwise preserved byte-for-byte after
# UTF-8 encoding.
MISSING_SENTINELS = frozenset(
    {"na", "nan", "none", "<na>", "n/a", "null", "unknown"}
)

# Corrected method IDs are part of the cache/run identity contract.  The
# preprocessing boundary uses ``preprocess``; corrected Stage 5 uses exactly
# the eight method labels below.
METHOD_IDS = frozenset(
    {
        "preprocess",
        "ECODA_authors_HR",
        "ECODA_seuratres_2",
        "ECODA_authors_HR_NULL",
        "Pseudobulk",
        "GloScope",
        "PILOT",
        "MrVI",
        "QOT",
    }
)

# Model IDs are likewise fixed policy tokens; changing one changes cache
# identity and therefore requires an explicit contract change.
MODEL_IDS = frozenset(
    {
        "hvg_composite_v1",
        "harmony_native_list_v1",
        "ecoda_additive_random_intercepts_v1",
        "pseudobulk_composite_v1",
        "mrvi_composite_v1",
        "embedding_consumer_harmony_v1",
    }
)

__all__ = [
    "TOKEN_VERSION",
    "FINGERPRINT_VERSION",
    "RESERVED_OBS_NAME",
    "DIRECT_SCALARIZATION",
    "COMPOSITE_SCALARIZATION",
    "VALIDATION_SUMMARY_SCHEMA_VERSION",
    "PREPROCESS_CORRECTION_MODE",
    "NATIVE_HARMONY_CORRECTION_MODE",
    "MRVI_CORRECTION_MODE",
    "PSEUDOBULK_CORRECTION_MODE",
    "ECODA_CORRECTION_MODE",
    "MISSING_SENTINELS",
    "METHOD_IDS",
    "MODEL_IDS",
    "BatchContractError",
    "BatchValidation",
    "BatchComposite",
    "normalize_batch_keys",
    "canonicalize_batch_value",
    "canonicalize_batch_values",
    "composite_token",
    "validate_batch_metadata",
    "build_batch_composite",
    "batch_contract_fingerprint",
    "build_batch_contract_identity",
    "build_batch_validation_summary",
    "validate_batch_validation_summary",
    "read_h5ad_validation_summary",
    "batch_correction_spec_for_keys",
    "batch_correction_spec",
    "augment_batch_contract",
    "serialize_batch_metadata",
]


class BatchContractError(ValueError):
    """Raised when corrected batch configuration or metadata is invalid."""


@dataclass(frozen=True)
class BatchValidation:
    """Validated full-cell batch metadata.

    ``canonical_values`` and ``row_tokens`` are ordered by input observation;
    ``sample_ids``, ``sample_group_ids``, and ``sample_composite_values`` are
    ordered by first observation of each validated sample.  ``levels`` and
    ``composite_levels`` are sorted by UTF-8 bytes.  The maps and tuples are
    newly allocated and never alias the input frame.
    """

    keys: tuple[str, ...]
    sample_column: str
    biological_column: str | None
    scalarization: str
    token_version: str
    n_obs: int
    sample_ids: tuple[str, ...]
    levels: dict[str, tuple[str, ...]]
    canonical_values: dict[str, tuple[str, ...]]
    row_tokens: tuple[str, ...]
    composite_levels: tuple[str, ...]
    design_rank: int
    design_columns: int
    composite_design_rank: int
    composite_design_columns: int
    sample_constancy: bool = True
    estimable: bool = True
    near_unique_fraction: float = 0.50
    # Appended defaults preserve the legacy positional/manual constructor
    # shape; validated instances always populate both tuples explicitly.
    sample_group_ids: tuple[str, ...] = ()
    sample_composite_values: tuple[str, ...] = ()

    @property
    def key_count(self) -> int:
        """Return the number of configured batch columns."""

        return len(self.keys)

    @property
    def n_samples(self) -> int:
        """Return the number of distinct standardized samples."""

        return len(self.sample_group_ids or self.sample_ids)

    @property
    def composite_level_count(self) -> int:
        """Return the number of composite levels, or zero for direct mode."""

        return len(self.composite_levels)

    @property
    def scalarized_values(self) -> tuple[str, ...]:
        """Return direct canonical values or full composite tokens by cell."""

        if self.scalarization == DIRECT_SCALARIZATION:
            return self.canonical_values[self.keys[0]]
        return self.row_tokens

    @property
    def sample_scalarized_values(self) -> tuple[str, ...]:
        """Return direct canonical values or composite tokens by sample."""

        return self.sample_composite_values


@dataclass(frozen=True)
class BatchComposite:
    """An ephemeral in-memory composite column and its validation result.

    ``frame`` is a deep copy when the input is a pandas frame.  Deleting or
    writing this returned copy cannot alter the source metadata.  The module
    deliberately provides no file-writing operation, so the reserved column
    cannot become a persisted H5AD field through this API.
    """

    frame: Any
    validation: BatchValidation
    values: tuple[str, ...]
    column_name: str = RESERVED_OBS_NAME

    @property
    def keys(self) -> tuple[str, ...]:
        """Return the ordered source keys."""

        return self.validation.keys

    @property
    def levels(self) -> tuple[str, ...]:
        """Return sorted full-token composite levels."""

        return self.validation.composite_levels

    @property
    def tokens(self) -> tuple[str, ...]:
        """Return composite tokens in input-observation order."""

        return self.values

    def metadata(
        self,
        *,
        method_id: str,
        model_id: str,
        include_tokens: bool = False,
    ) -> dict[str, Any]:
        """Serialize the run-owned contract metadata for this composite."""

        return serialize_batch_metadata(
            self.validation,
            method_id=method_id,
            model_id=model_id,
            include_tokens=include_tokens,
        )


def _validate_column_name(value: object, role: str) -> str:
    if not isinstance(value, str) or not value or not value.strip():
        raise BatchContractError(f"{role} must be a nonempty string")
    return value


def _validate_key_vector(
    batch_keys: object,
    *,
    sample_column: str | None = None,
    biological_column: str | None = None,
) -> tuple[str, ...]:
    """Validate an ordered key vector without converting its elements."""

    if isinstance(batch_keys, str):
        raw_keys: list[object] = [batch_keys]
    elif batch_keys is None:
        raise BatchContractError("corrected columns.batch must not be null")
    elif isinstance(batch_keys, (bytes, bytearray, Mapping, set, frozenset)):
        raise BatchContractError(
            "corrected columns.batch must be a string or an ordered list of strings"
        )
    elif isinstance(batch_keys, Sequence):
        raw_keys = list(batch_keys)
    else:
        raise BatchContractError(
            "corrected columns.batch must be a string or an ordered list of strings"
        )

    if not raw_keys:
        raise BatchContractError("corrected columns.batch must contain at least one key")

    keys: list[str] = []
    seen: set[str] = set()
    for position, key in enumerate(raw_keys):
        if not isinstance(key, str):
            raise BatchContractError(
                f"corrected batch key at position {position} must be a string"
            )
        if not key or not key.strip():
            raise BatchContractError(
                f"corrected batch key at position {position} must be nonblank"
            )
        if key in seen:
            raise BatchContractError(f"corrected batch keys contain duplicate {key!r}")
        seen.add(key)
        keys.append(key)

    if sample_column is not None and any(key == sample_column for key in keys):
        raise BatchContractError(
            f"corrected batch key must not equal standardized sample column {sample_column!r}"
        )
    if biological_column is not None and any(key == biological_column for key in keys):
        raise BatchContractError(
            f"corrected batch key must not equal biological label column {biological_column!r}"
        )
    return tuple(keys)


def normalize_batch_keys(
    batch_keys: str | Sequence[str] | None,
    *,
    sample_column: str = "Sample",
    biological_column: str | None = None,
) -> tuple[str, ...]:
    """Return corrected ``columns.batch`` as an ordered, validated key tuple.

    A scalar string becomes a one-element tuple.  ``None`` and an empty list
    are rejected in corrected mode; no raw configuration object is changed.
    """

    sample_column = _validate_column_name(sample_column, "sample_column")
    if biological_column is not None:
        biological_column = _validate_column_name(
            biological_column, "biological_column"
        )
    return _validate_key_vector(
        batch_keys,
        sample_column=sample_column,
        biological_column=biological_column,
    )


def _is_date_like(value: object) -> bool:
    if isinstance(value, (_datetime.date, _datetime.datetime, _datetime.time)):
        return True
    if isinstance(value, (np.datetime64, np.timedelta64)):
        return True
    if isinstance(value, (pd.Timestamp, pd.Timedelta)):
        return True
    return False


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if _is_date_like(value):
        # NaT is missing, while ordinary date/time values are rejected by the
        # date-like branch in canonicalize_batch_value.
        try:
            return bool(pd.isna(value))
        except (TypeError, ValueError):
            return False
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    if isinstance(missing, (bool, np.bool_)):
        return bool(missing)
    return False


def _reject_text_placeholder(text: str, *, value: object) -> None:
    if not text.strip():
        raise BatchContractError(f"batch value is blank: {value!r}")
    if text.strip().casefold() in MISSING_SENTINELS:
        raise BatchContractError(f"batch value is a missing sentinel: {value!r}")


def _is_signed_integer(value: object) -> bool:
    if isinstance(value, bool) or isinstance(value, np.bool_):
        return False
    if isinstance(value, int):
        return True
    return isinstance(value, np.integer) and getattr(value.dtype, "kind", "") == "i"


def _is_float(value: object) -> bool:
    if isinstance(value, bool) or isinstance(value, np.bool_):
        return False
    if isinstance(value, float):
        return True
    if isinstance(value, np.floating):
        dtype = np.dtype(value.dtype)
        return dtype.kind == "f" and dtype.itemsize == np.dtype(np.float64).itemsize
    return False


def _is_supported_factor_scalar(value: object) -> bool:
    if isinstance(value, (str, np.str_)):
        return True
    if isinstance(value, (bool, np.bool_)):
        return True
    if _is_signed_integer(value) or _is_float(value):
        return True
    return False


def canonicalize_batch_value(value: object, *, factor: bool = False) -> str:
    """Canonicalize one nonmissing corrected batch category value.

    Strings and factor labels use ``s:<exact label>``; booleans use
    ``b:true``/``b:false``; signed integers use ``i:<base-10>``; and finite
    IEEE-754 binary64 values use ``f64:<16 lowercase hex bits>``.  NumPy
    floating dtypes other than ``float64`` are rejected rather than widened.
    Lists, dates, arbitrary objects, unsigned integers, missing values, and
    sentinels fail closed rather than being stringified.
    """

    if _is_missing(value):
        raise BatchContractError(f"batch value is missing: {value!r}")
    if _is_date_like(value):
        raise BatchContractError(f"date/time batch values are unsupported: {value!r}")

    if factor:
        if not _is_supported_factor_scalar(value):
            raise BatchContractError(
                f"factor batch value has unsupported type {type(value).__name__}"
            )
        # Factor labels are textual categories even when a Python categorical
        # happens to store numeric category objects.
        text = str(value)
        _reject_text_placeholder(text, value=value)
        return f"s:{text}"

    if isinstance(value, (bool, np.bool_)):
        return "b:true" if bool(value) else "b:false"

    if _is_signed_integer(value):
        return f"i:{int(value):d}"

    if _is_float(value):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise BatchContractError(f"float batch value is not finite: {value!r}")
        bits = struct.pack(">d", numeric).hex()
        return f"f64:{bits}"

    if isinstance(value, (str, np.str_)):
        text = str(value)
        _reject_text_placeholder(text, value=value)
        return f"s:{text}"

    raise BatchContractError(
        f"batch value has unsupported type {type(value).__name__}; "
        "only strings, factors, booleans, signed integers, and finite floats are allowed"
    )


def canonicalize_batch_values(
    values: Sequence[object], *, factor: bool = False
) -> tuple[str, ...]:
    """Canonicalize a sequence in input order without modifying it."""

    if isinstance(values, (str, bytes, bytearray)):
        raise BatchContractError("batch values must be a sequence of scalar values")
    try:
        return tuple(
            canonicalize_batch_value(value, factor=factor) for value in values
        )
    except TypeError as exc:
        raise BatchContractError("batch values must be a sequence") from exc

def _sample_identifier(value: object, *, factor: bool) -> str:
    """Validate and return one standardized Sample label without coercion."""

    if _is_missing(value):
        raise BatchContractError(f"standardized Sample value is missing: {value!r}")
    if not factor and not isinstance(value, (str, np.str_)):
        raise BatchContractError(
            "standardized Sample values must be strings or factor labels"
        )
    if factor and not _is_supported_factor_scalar(value):
        raise BatchContractError(
            f"standardized Sample factor has unsupported value type {type(value).__name__}"
        )
    text = str(value)
    _reject_text_placeholder(text, value=value)
    return text


def _encode_canonical_token(keys: tuple[str, ...], values: tuple[str, ...]) -> str:
    if len(keys) != len(values):
        raise BatchContractError("batch key/value lengths differ")
    pairs: list[str] = []
    for key, value in zip(keys, values):
        key_bytes = key.encode("utf-8")
        value_bytes = value.encode("utf-8")
        pairs.append(
            f"{len(key_bytes)}:{key_bytes.hex()},"
            f"{len(value_bytes)}:{value_bytes.hex()}"
        )
    return f"{TOKEN_VERSION}|{len(keys)}|" + ";".join(pairs)


def _ordered_values(
    keys: tuple[str, ...],
    values: Sequence[object] | Mapping[str, object] | object,
) -> tuple[object, ...]:
    if isinstance(values, Mapping):
        value_keys = tuple(values.keys())
        if set(value_keys) != set(keys) or len(value_keys) != len(keys):
            raise BatchContractError("batch value mapping keys do not match batch keys")
        if any(not isinstance(key, str) for key in value_keys):
            raise BatchContractError("batch value mapping keys must be strings")
        return tuple(values[key] for key in keys)
    if len(keys) == 1 and not isinstance(values, (Sequence, np.ndarray)):
        return (values,)
    if isinstance(values, (str, bytes, bytearray)):
        if len(keys) == 1:
            return (values,)
        raise BatchContractError("batch values must provide one value per key")
    try:
        ordered = tuple(values)  # type: ignore[arg-type]
    except TypeError as exc:
        raise BatchContractError("batch values must provide one value per key") from exc
    if len(ordered) != len(keys):
        raise BatchContractError("batch key/value lengths differ")
    return ordered


def composite_token(
    keys: str | Sequence[str],
    values: Sequence[object] | Mapping[str, object] | object,
    *,
    factors: Sequence[bool] | Mapping[str, bool] | None = None,
) -> str:
    """Encode one ordered composite token using the exact UTF-8 contract.

    Composite construction is defined only for two or more ordered keys.
    Scalar and one-key corrected configurations use direct canonical values.
    """

    normalized = normalize_batch_keys(keys)
    if len(normalized) < 2:
        raise BatchContractError(
            "composite token construction requires at least two batch keys"
        )
    ordered = _ordered_values(normalized, values)

    factor_by_key: dict[str, bool] = {key: False for key in normalized}
    if factors is not None:
        if isinstance(factors, Mapping):
            if set(factors) != set(normalized):
                raise BatchContractError("factor mapping keys do not match batch keys")
            factor_by_key = {}
            for key in normalized:
                marker = factors[key]
                if not isinstance(marker, (bool, np.bool_)):
                    raise BatchContractError("factor markers must be boolean")
                factor_by_key[key] = bool(marker)
        else:
            if isinstance(
                factors, (str, bytes, bytearray, set, frozenset)
            ):
                raise BatchContractError("factor markers must be an ordered sequence")
            markers = tuple(factors)
            if len(markers) != len(normalized):
                raise BatchContractError("batch key/factor lengths differ")
            for key, marker in zip(normalized, markers):
                if not isinstance(marker, (bool, np.bool_)):
                    raise BatchContractError("factor markers must be boolean")
                factor_by_key[key] = bool(marker)

    canonical = tuple(
        canonicalize_batch_value(value, factor=factor_by_key[key])
        for key, value in zip(normalized, ordered)
    )
    return _encode_canonical_token(normalized, canonical)


def _data_columns(data: object) -> tuple[list[object], int]:
    if isinstance(data, pd.DataFrame):
        columns = list(data.columns)
        if len(set(columns)) != len(columns):
            raise BatchContractError("metadata contains duplicate column names")
        return columns, len(data)
    if isinstance(data, Mapping):
        columns = list(data.keys())
        if len(set(columns)) != len(columns):
            raise BatchContractError("metadata contains duplicate column names")
        lengths: set[int] = set()
        for column in columns:
            try:
                lengths.add(len(data[column]))
            except (TypeError, KeyError) as exc:
                raise BatchContractError(
                    f"metadata column {column!r} is not a sized sequence"
                ) from exc
        if len(lengths) > 1:
            raise BatchContractError("metadata columns have different lengths")
        return columns, next(iter(lengths), 0)
    columns = getattr(data, "columns", None)
    if columns is None:
        raise BatchContractError("metadata must be a pandas DataFrame or column mapping")
    columns = list(columns)
    if len(set(columns)) != len(columns):
        raise BatchContractError("metadata contains duplicate column names")
    try:
        n_obs = len(data)
    except TypeError as exc:
        raise BatchContractError("metadata has no observation length") from exc
    return columns, n_obs


def _column_values(data: object, column: str, n_obs: int) -> tuple[tuple[object, ...], bool]:
    try:
        values = data[column]  # type: ignore[index]
    except (KeyError, TypeError) as exc:
        raise BatchContractError(f"metadata is missing column {column!r}") from exc
    factor = isinstance(getattr(values, "dtype", None), pd.CategoricalDtype)
    if isinstance(values, (str, bytes, bytearray)):
        raise BatchContractError(f"metadata column {column!r} is not a sequence")
    try:
        sequence = tuple(values)
    except TypeError as exc:
        raise BatchContractError(f"metadata column {column!r} is not a sequence") from exc
    if len(sequence) != n_obs:
        raise BatchContractError(
            f"metadata column {column!r} has {len(sequence)} rows; expected {n_obs}"
        )
    return sequence, factor


def _sorted_levels(values: Sequence[str]) -> tuple[str, ...]:
    return tuple(sorted(set(values), key=lambda value: value.encode("utf-8")))


def _design_rank(
    rows: Sequence[Mapping[str, str]],
    keys: Sequence[str],
    levels: Mapping[str, Sequence[str]],
) -> tuple[int, int]:
    """Return rank/column count for an intercept plus additive factors."""

    n_rows = len(rows)
    columns: list[np.ndarray] = [np.ones(n_rows, dtype=np.float64)]
    for key in keys:
        # The first byte-sorted level is the deterministic reference level.
        for level in levels[key][1:]:
            columns.append(
                np.fromiter(
                    (1.0 if row[key] == level else 0.0 for row in rows),
                    dtype=np.float64,
                    count=n_rows,
                )
            )
    matrix = np.column_stack(columns)
    try:
        rank = int(np.linalg.matrix_rank(matrix))
    except (TypeError, ValueError, np.linalg.LinAlgError) as exc:
        raise BatchContractError("corrected batch design rank could not be estimated") from exc
    return rank, int(matrix.shape[1])


def _single_factor_rank(values: Sequence[str], levels: Sequence[str]) -> tuple[int, int]:
    rows = [{"__factor__": value} for value in values]
    return _design_rank(rows, ("__factor__",), {"__factor__": levels})


def validate_batch_metadata(
    data: pd.DataFrame | Mapping[str, Sequence[object]],
    batch_keys: str | Sequence[str] | None,
    *,
    sample_column: str = "Sample",
    biological_column: str | None = None,
    near_unique_fraction: float = 0.50,
) -> BatchValidation:
    """Validate every cell's corrected batch metadata before reduction.

    The check rejects missing/sentinel values, within-``Sample`` disagreement,
    constant or near-unique factors, disconnected/rank-deficient additive
    designs, and non-estimable composite levels.  It consumes the entire cell
    table and never chooses a first observation as a fallback.
    """

    sample_column = _validate_column_name(sample_column, "sample_column")
    if biological_column is not None:
        biological_column = _validate_column_name(
            biological_column, "biological_column"
        )
    if (
        isinstance(near_unique_fraction, (bool, np.bool_))
        or not isinstance(near_unique_fraction, (int, float, np.number))
        or not math.isfinite(float(near_unique_fraction))
        or not 0.0 < float(near_unique_fraction) < 1.0
    ):
        raise BatchContractError("near_unique_fraction must be finite and in (0, 1)")
    near_unique_fraction = float(near_unique_fraction)

    keys = normalize_batch_keys(
        batch_keys,
        sample_column=sample_column,
        biological_column=biological_column,
    )
    columns, n_obs = _data_columns(data)
    available = set(columns)
    if RESERVED_OBS_NAME in available:
        raise BatchContractError(
            f"metadata already contains reserved temporary column {RESERVED_OBS_NAME!r}"
        )
    required = (sample_column, *keys)
    missing_columns = [column for column in required if column not in available]
    if missing_columns:
        raise BatchContractError(f"metadata is missing required columns {missing_columns}")
    if n_obs <= 0:
        raise BatchContractError("metadata has no observations")

    sample_values, sample_factor = _column_values(data, sample_column, n_obs)
    raw_sample_labels = tuple(
        _sample_identifier(value, factor=sample_factor) for value in sample_values
    )
    canonical_samples = tuple(
        canonicalize_batch_value(value, factor=sample_factor) for value in sample_values
    )
    sample_ids: list[str] = []
    sample_group_ids: list[str] = []
    sample_seen: set[str] = set()
    for label, sample in zip(raw_sample_labels, canonical_samples):
        if sample not in sample_seen:
            sample_seen.add(sample)
            sample_group_ids.append(sample)
            sample_ids.append(label)
    if len(sample_group_ids) < 2:
        raise BatchContractError("corrected batch metadata requires at least two samples")
    canonical_by_key: dict[str, tuple[str, ...]] = {}
    for key in keys:
        values, factor = _column_values(data, key, n_obs)
        canonical_by_key[key] = tuple(
            canonicalize_batch_value(value, factor=factor) for value in values
        )

    # Full-cell constancy check.  A sample's first value is used only as a
    # comparison anchor after every value has been read and canonicalized; a
    # disagreement is an error, never a silently selected value.
    first_by_sample: dict[str, dict[str, str]] = {}
    first_row_by_sample: dict[str, int] = {}
    for row_number, sample in enumerate(canonical_samples):
        first_row_by_sample.setdefault(sample, row_number)
        sample_first = first_by_sample.setdefault(sample, {})
        for key in keys:
            value = canonical_by_key[key][row_number]
            prior = sample_first.get(key)
            if prior is None:
                sample_first[key] = value
            elif prior != value:
                raise BatchContractError(
                    f"batch key {key!r} disagrees within Sample {sample!r} "
                    f"(rows {first_row_by_sample[sample]} and {row_number})"
                )

    levels = {
        key: _sorted_levels(canonical_by_key[key])
        for key in keys
    }
    for key in keys:
        key_levels = levels[key]
        if len(key_levels) < 2:
            raise BatchContractError(
                f"corrected batch key {key!r} has fewer than two observed levels"
            )
        ratio = len(key_levels) / len(sample_group_ids)
        if ratio > near_unique_fraction:
            raise BatchContractError(
                f"corrected batch key {key!r} is near-unique: "
                f"{len(key_levels)}/{len(sample_group_ids)} levels "
                f"(threshold {near_unique_fraction:g})"
            )
    sample_rows = [
        {key: first_by_sample[sample][key] for key in keys}
        for sample in sample_group_ids
    ]
    design_rank, design_columns = _design_rank(sample_rows, keys, levels)
    if design_rank < design_columns:
        raise BatchContractError(
            "corrected additive batch design is rank-deficient or disconnected"
        )

    scalarization = (
        COMPOSITE_SCALARIZATION if len(keys) >= 2 else DIRECT_SCALARIZATION
    )
    if len(keys) >= 2:
        row_tokens = tuple(
            _encode_canonical_token(
                keys,
                tuple(canonical_by_key[key][row_number] for key in keys),
            )
            for row_number in range(n_obs)
        )
        composite_levels = _sorted_levels(row_tokens)
        if len(composite_levels) < 2:
            raise BatchContractError(
                "corrected composite batch has fewer than two observed levels"
            )
        composite_ratio = len(composite_levels) / len(sample_group_ids)
        if composite_ratio > near_unique_fraction:
            raise BatchContractError(
                "corrected composite batch is near-unique: "
                f"{len(composite_levels)}/{len(sample_group_ids)} levels "
                f"(threshold {near_unique_fraction:g})"
            )
        sample_composite_values = tuple(
            _encode_canonical_token(
                keys, tuple(first_by_sample[sample][key] for key in keys)
            )
            for sample in sample_group_ids
        )
        composite_design_rank, composite_design_columns = _single_factor_rank(
            sample_composite_values, composite_levels
        )
    else:
        sample_composite_values = tuple(
            first_by_sample[sample][keys[0]] for sample in sample_group_ids
        )
        row_tokens = ()
        composite_levels = ()
        composite_design_rank = design_rank
        composite_design_columns = design_columns
    return BatchValidation(
        keys=keys,
        sample_column=sample_column,
        biological_column=biological_column,
        scalarization=scalarization,
        token_version=TOKEN_VERSION,
        n_obs=n_obs,
        sample_ids=tuple(sample_ids),
        sample_group_ids=tuple(sample_group_ids),
        levels=levels,
        canonical_values=canonical_by_key,
        row_tokens=row_tokens,
        composite_levels=composite_levels,
        sample_composite_values=sample_composite_values,
        design_rank=design_rank,
        design_columns=design_columns,
        composite_design_rank=composite_design_rank,
        composite_design_columns=composite_design_columns,
        near_unique_fraction=near_unique_fraction,
    )


def _copy_with_composite(data: object, values: Sequence[str]) -> object:
    if isinstance(data, pd.DataFrame):
        copied = data.copy(deep=True)
        if RESERVED_OBS_NAME in copied.columns:
            raise BatchContractError(
                f"metadata already contains reserved temporary column {RESERVED_OBS_NAME!r}"
            )
        copied.insert(len(copied.columns), RESERVED_OBS_NAME, list(values))
        return copied
    if isinstance(data, Mapping):
        copied: dict[object, object] = {}
        for key, column in data.items():
            try:
                copied[key] = tuple(column)
            except TypeError as exc:
                raise BatchContractError(
                    f"metadata column {key!r} is not a sequence"
                ) from exc
        copied[RESERVED_OBS_NAME] = tuple(values)
        return copied
    raise BatchContractError("metadata must be a pandas DataFrame or column mapping")


def build_batch_composite(
    data: pd.DataFrame | Mapping[str, Sequence[object]],
    batch_keys: str | Sequence[str] | None,
    *,
    sample_column: str = "Sample",
    biological_column: str | None = None,
    near_unique_fraction: float = 0.50,
) -> BatchComposite:
    """Build a validated, copied, ephemeral composite column for 2+ keys.

    Scalar and one-key configurations intentionally remain direct columns and
    therefore fail here instead of accidentally changing representation.
    """

    validation = validate_batch_metadata(
        data,
        batch_keys,
        sample_column=sample_column,
        biological_column=biological_column,
        near_unique_fraction=near_unique_fraction,
    )
    if validation.scalarization != COMPOSITE_SCALARIZATION:
        raise BatchContractError(
            "composite construction requires at least two corrected batch keys"
        )
    copied = _copy_with_composite(data, validation.row_tokens)
    return BatchComposite(
        frame=copied,
        validation=validation,
        values=validation.row_tokens,
    )


def _field_bytes(name: str, value: str | bytes) -> bytes:
    name_bytes = name.encode("utf-8")
    value_bytes = value if isinstance(value, bytes) else value.encode("utf-8")
    return (
        str(len(name_bytes)).encode("ascii")
        + b":"
        + name_bytes.hex().encode("ascii")
        + b","
        + str(len(value_bytes)).encode("ascii")
        + b":"
        + value_bytes.hex().encode("ascii")
        + b";"
    )


def _ordered_key_vector(keys: tuple[str, ...]) -> str:
    parts = [f"{len(keys)}|"]
    for key in keys:
        key_bytes = key.encode("utf-8")
        parts.append(f"{len(key_bytes)}:{key_bytes.hex()};")
    return "".join(parts)


def _validate_method_model(method_id: object, model_id: object) -> tuple[str, str]:
    if not isinstance(method_id, str) or not method_id or not method_id.strip():
        raise BatchContractError("method_id must be a nonempty string")
    if method_id not in METHOD_IDS:
        raise BatchContractError(f"unsupported corrected batch method_id {method_id!r}")
    if not isinstance(model_id, str) or not model_id or not model_id.strip():
        raise BatchContractError("model_id must be a nonempty string")
    if model_id not in MODEL_IDS:
        raise BatchContractError(f"unsupported corrected batch model_id {model_id!r}")
    return method_id, model_id


def _batch_contract_payload(
    keys: tuple[str, ...],
    scalarization: str,
    method_id: str,
    model_id: str,
) -> bytes:
    payload = bytearray(FINGERPRINT_VERSION.encode("utf-8") + b"\0")
    payload.extend(_field_bytes("encoding", TOKEN_VERSION))
    payload.extend(_field_bytes("keys", _ordered_key_vector(keys)))
    payload.extend(_field_bytes("scalarization", scalarization))
    payload.extend(_field_bytes("method", method_id))
    payload.extend(_field_bytes("model", model_id))
    return bytes(payload)


def batch_contract_fingerprint(
    batch_keys: str | Sequence[str],
    scalarization: str | None = None,
    method_id: str | None = None,
    model_id: str | None = None,
) -> str:
    """Return the lowercase SHA-256 key/configuration fingerprint.

    The hash input is the exact byte stream specified by the contract; no JSON,
    locale, map iteration, or platform-native integer representation is used.
    """

    keys = normalize_batch_keys(batch_keys)
    expected_scalarization = (
        COMPOSITE_SCALARIZATION if len(keys) >= 2 else DIRECT_SCALARIZATION
    )
    if scalarization is None:
        scalarization = expected_scalarization
    if scalarization != expected_scalarization:
        raise BatchContractError(
            f"scalarization {scalarization!r} does not match {len(keys)} batch key(s)"
        )
    if method_id is None or model_id is None:
        raise BatchContractError("method_id and model_id are required for fingerprinting")
    method_id, model_id = _validate_method_model(method_id, model_id)
    return hashlib.sha256(
        _batch_contract_payload(keys, scalarization, method_id, model_id)
    ).hexdigest()


def build_batch_contract_identity(
    batch_keys: str | Sequence[str] | None,
    sample_column: str = "Sample",
    method_id: str | None = None,
    model_id: str | None = None,
) -> dict[str, Any]:
    """Build lightweight source/configuration identity for corrected batches.

    This builder deliberately consumes configuration only.  It never reads
    cell metadata and therefore cannot include per-cell or per-sample values.
    The returned mapping is JSON-safe and uses the canonical identity field
    names consumed by corrected artifact validators.
    """

    sample_column = _validate_column_name(sample_column, "sample_column")
    if sample_column != sample_column.strip():
        raise BatchContractError(
            "sample_column must be a nonblank string without surrounding whitespace"
        )
    keys = normalize_batch_keys(batch_keys, sample_column=sample_column)
    if any(key != key.strip() for key in keys):
        raise BatchContractError(
            "corrected batch keys must be nonblank strings without surrounding whitespace"
        )
    if sample_column == RESERVED_OBS_NAME:
        raise BatchContractError(
            f"sample_column uses the reserved temporary name {RESERVED_OBS_NAME!r}"
        )
    if RESERVED_OBS_NAME in keys:
        raise BatchContractError(
            f"corrected batch keys contain the reserved temporary name {RESERVED_OBS_NAME!r}"
        )

    method_id, model_id = _validate_method_model(method_id, model_id)
    scalarization = (
        COMPOSITE_SCALARIZATION if len(keys) >= 2 else DIRECT_SCALARIZATION
    )
    fingerprint = hashlib.sha256(
        _batch_contract_payload(keys, scalarization, method_id, model_id)
    ).hexdigest()
    return {
        "contract_version": FINGERPRINT_VERSION,
        "token_version": TOKEN_VERSION,
        "ordered_source_keys": list(keys),
        "scalarization": scalarization,
        "method_id": method_id,
        "model_id": model_id,
        "required_source_obs_columns": [sample_column, *keys],
        "reserved_obs_name": RESERVED_OBS_NAME,
        "reserved_obs_absent": True,
        "fingerprint": fingerprint,
    }

def _validated_batch_value(
    validation: BatchValidation | BatchComposite,
) -> BatchValidation:
    """Return one validated result for compact-summary construction."""

    if isinstance(validation, BatchComposite):
        validation = validation.validation
    if not isinstance(validation, BatchValidation):
        raise BatchContractError(
            "batch validation summary expects BatchValidation"
        )
    if not validation.keys:
        raise BatchContractError("batch validation summary requires configured keys")
    if validation.n_obs <= 0 or validation.n_samples < 2:
        raise BatchContractError(
            "batch validation summary requires nonempty validated metadata"
        )
    for key in validation.keys:
        try:
            levels = validation.levels[key]
        except (KeyError, TypeError) as exc:
            raise BatchContractError(
                f"batch validation summary is missing levels for {key!r}"
            ) from exc
        if tuple(levels) != _sorted_levels(levels) or len(levels) < 2:
            raise BatchContractError(
                f"batch validation summary has invalid sorted levels for {key!r}"
            )
    if validation.scalarization == DIRECT_SCALARIZATION:
        if validation.composite_levels:
            raise BatchContractError(
                "direct batch validation cannot carry composite levels"
            )
    elif validation.scalarization == COMPOSITE_SCALARIZATION:
        if not validation.composite_levels:
            raise BatchContractError(
                "composite batch validation requires composite levels"
            )
        if tuple(validation.composite_levels) != _sorted_levels(
            validation.composite_levels
        ):
            raise BatchContractError(
                "batch validation summary has invalid sorted composite levels"
            )
    else:
        raise BatchContractError(
            f"unsupported batch validation scalarization {validation.scalarization!r}"
        )
    return validation


def build_batch_validation_summary(
    validation: BatchValidation | BatchComposite,
    correction_mode: str,
    correction_formula: str,
) -> dict[str, Any]:
    """Build the compact, vector-free summary of validated batch metadata.

    The source ``BatchValidation`` contains full cell and sample vectors for
    in-memory consumers.  This function intentionally projects only bounded
    levels/counts and validation facts needed by a persisted corrected
    artifact.  Its output has a fixed key set and fresh ordered maps/lists, so
    callers can safely store it in ``AnnData.uns`` or JSON metadata.
    """

    validation = _validated_batch_value(validation)
    for value, label in (
        (correction_mode, "correction_mode"),
        (correction_formula, "correction_formula"),
    ):
        if (
            not isinstance(value, str)
            or not value
            or value != value.strip()
        ):
            raise BatchContractError(
                f"{label} must be a nonblank string without surrounding whitespace"
            )

    per_key_levels = {
        key: list(validation.levels[key]) for key in validation.keys
    }
    key_level_counts = {
        key: len(per_key_levels[key]) for key in validation.keys
    }
    sample_constancy = {key: True for key in validation.keys}
    return {
        "schema_version": VALIDATION_SUMMARY_SCHEMA_VERSION,
        "validated_before_reduction": True,
        "sample_constancy": sample_constancy,
        "per_key_levels": per_key_levels,
        "key_level_counts": key_level_counts,
        "composite_levels": list(validation.composite_levels),
        "composite_level_count": len(validation.composite_levels),
        "n_cells": int(validation.n_obs),
        "n_samples": int(validation.n_samples),
        "correction_mode": correction_mode,
        "correction_formula": correction_formula,
    }


def batch_correction_spec_for_keys(
    method_id: str,
    batch_keys: str | Sequence[str],
) -> tuple[str, str]:
    """Return one deterministic correction policy from configured key order."""

    keys = normalize_batch_keys(batch_keys)
    if not isinstance(method_id, str) or not method_id:
        raise BatchContractError("method_id must be a nonblank string")
    scalar = keys[0] if len(keys) == 1 else RESERVED_OBS_NAME
    ordered_keys = f"[{','.join(keys)}]"
    if method_id == "preprocess":
        return (
            PREPROCESS_CORRECTION_MODE,
            f"HVG batch_key={scalar}; Harmony vars_use={ordered_keys}",
        )
    if method_id in {"GloScope", "PILOT", "QOT"}:
        return (
            NATIVE_HARMONY_CORRECTION_MODE,
            "embedding=X_pca_harmony_batch_effect_corrected_hvg2000",
        )
    if method_id == "MrVI":
        return (
            MRVI_CORRECTION_MODE,
            f"MRVI.setup_anndata(batch_key={scalar})",
        )
    if method_id == "Pseudobulk":
        return (
            PSEUDOBULK_CORRECTION_MODE,
            f"DESeq2 design=~ 1; limma removeBatchEffect(batch={scalar})",
        )
    if method_id in {
        "ECODA_authors_HR",
        "ECODA_seuratres_2",
        "ECODA_authors_HR_NULL",
    }:
        aliases = (
            ["batch"]
            if len(keys) == 1
            else [f"batch_key_{index}" for index in range(1, len(keys) + 1)]
        )
        formula = "y ~ 1 + " + " + ".join(
            f"(1 | {alias})" for alias in aliases
        )
        return ECODA_CORRECTION_MODE, formula
    raise BatchContractError(
        f"unsupported corrected batch consumer {method_id!r}"
    )


def batch_correction_spec(
    method_id: str,
    validation: BatchValidation | BatchComposite,
) -> tuple[str, str]:
    """Return the explicit corrected-mode summary policy for one consumer."""

    validation = _validated_batch_value(validation)
    return batch_correction_spec_for_keys(method_id, validation.keys)


def augment_batch_contract(
    contract: Mapping[str, Any],
    validation: BatchValidation | BatchComposite | Mapping[str, Any],
    correction_mode: str,
    correction_formula: str,
) -> dict[str, Any]:
    """Attach one compact validated summary without changing identity fields."""

    if not isinstance(contract, Mapping):
        raise BatchContractError("batch contract identity must be a mapping")
    forbidden = {
        "composite_values",
        "sample_composite_values",
        "sample_ids",
        "sample_group_ids",
        "tokens",
        "canonical_values",
        "canonical_cell_values",
        "canonical_sample_metadata",
        "sample_metadata",
        "row_tokens",
        "scalarized_values",
        "sample_scalarized_values",
    }
    if isinstance(validation, Mapping):
        summary = validate_batch_validation_summary(validation)
        for value, label in (
            (correction_mode, "correction_mode"),
            (correction_formula, "correction_formula"),
        ):
            if (
                not isinstance(value, str)
                or not value
                or value != value.strip()
            ):
                raise BatchContractError(
                    f"{label} must be a nonblank string without surrounding whitespace"
                )
        summary["correction_mode"] = correction_mode
        summary["correction_formula"] = correction_formula
    else:
        summary = build_batch_validation_summary(
            validation,
            correction_mode,
            correction_formula,
        )
    identity_keys = contract.get("ordered_source_keys")
    if identity_keys is not None:
        summary = validate_batch_validation_summary(summary, identity_keys)
    augmented = {
        key: value for key, value in contract.items() if key not in forbidden
    }
    augmented["validation_summary"] = summary
    return augmented

def validate_batch_validation_summary(
    summary: Mapping[str, Any],
    batch_keys: str | Sequence[str] | None = None,
) -> dict[str, Any]:
    """Validate and copy one persisted vector-free validation summary."""

    if not isinstance(summary, Mapping):
        raise BatchContractError("validation_summary must be a mapping")
    if set(summary) != _VALIDATION_SUMMARY_FIELDS:
        raise BatchContractError(
            "validation_summary has an invalid field set"
        )

    def _normalize_numpy_scalar(value: Any) -> Any:
        if isinstance(value, np.bool_):
            return bool(value)
        if isinstance(value, np.integer):
            return int(value)
        return value

    def _normalize_numpy_vector(value: Any) -> Any:
        if isinstance(value, np.ndarray) and value.ndim == 1:
            return value.tolist()
        return value

    normalized = dict(summary)
    for field in (
        "schema_version",
        "validated_before_reduction",
        "composite_level_count",
        "n_cells",
        "n_samples",
    ):
        normalized[field] = _normalize_numpy_scalar(summary[field])
    constancy = summary["sample_constancy"]
    if isinstance(constancy, Mapping):
        normalized["sample_constancy"] = {
            key: _normalize_numpy_scalar(value)
            for key, value in constancy.items()
        }
    levels = summary["per_key_levels"]
    if isinstance(levels, Mapping):
        normalized["per_key_levels"] = {
            key: _normalize_numpy_vector(value)
            for key, value in levels.items()
        }
    counts = summary["key_level_counts"]
    if isinstance(counts, Mapping):
        normalized["key_level_counts"] = {
            key: _normalize_numpy_scalar(value)
            for key, value in counts.items()
        }
    normalized["composite_levels"] = _normalize_numpy_vector(
        summary["composite_levels"]
    )
    summary = normalized
    if (
        type(summary["schema_version"]) is not int
        or summary["schema_version"] != VALIDATION_SUMMARY_SCHEMA_VERSION
    ):
        raise BatchContractError("validation_summary has an invalid schema_version")
    if summary["validated_before_reduction"] is not True:
        raise BatchContractError(
            "validation_summary must assert validated_before_reduction"
        )

    constancy = summary["sample_constancy"]
    levels = summary["per_key_levels"]
    counts = summary["key_level_counts"]
    if not all(isinstance(value, Mapping) for value in (constancy, levels, counts)):
        raise BatchContractError(
            "validation_summary key fields must be ordered maps"
        )
    constancy_keys = list(constancy)
    level_keys = list(levels)
    count_keys = list(counts)
    if (
        len(constancy_keys) != len(set(constancy_keys))
        or len(level_keys) != len(set(level_keys))
        or len(count_keys) != len(set(count_keys))
        or set(constancy_keys) != set(level_keys)
        or set(level_keys) != set(count_keys)
    ):
        raise BatchContractError(
            "validation_summary key maps must contain the same configured keys"
        )
    if batch_keys is not None:
        expected_keys = normalize_batch_keys(batch_keys)
        if set(level_keys) != set(expected_keys):
            raise BatchContractError(
                "validation_summary keys do not match configured batch keys"
            )
        # HDF5 group iteration may be lexical even though the source summary
        # was written in configured order.  Return a canonical configured map
        # order after checking that no key was added or dropped.
        keys = list(expected_keys)
    else:
        if level_keys != constancy_keys or level_keys != count_keys:
            raise BatchContractError(
                "validation_summary key maps must use one configured order"
            )
        keys = level_keys
    if not keys or any(
        not isinstance(key, str) or not key or key != key.strip() for key in keys
    ):
        raise BatchContractError("validation_summary contains invalid batch keys")
    for key in keys:
        if constancy[key] is not True:
            raise BatchContractError(
                f"validation_summary sample_constancy is false for {key!r}"
            )
        key_levels = levels[key]
        if (
            isinstance(key_levels, str)
            or not isinstance(key_levels, (list, tuple))
            or any(not isinstance(level, str) for level in key_levels)
            or tuple(key_levels) != _sorted_levels(key_levels)
            or len(key_levels) < 2
        ):
            raise BatchContractError(
                f"validation_summary has invalid levels for {key!r}"
            )
        level_count = counts[key]
        if (
            isinstance(level_count, bool)
            or not isinstance(level_count, int)
            or level_count != len(key_levels)
        ):
            raise BatchContractError(
                f"validation_summary has invalid level count for {key!r}"
            )
    composite_levels = summary["composite_levels"]
    if (
        isinstance(composite_levels, str)
        or not isinstance(composite_levels, (list, tuple))
        or any(not isinstance(level, str) for level in composite_levels)
        or tuple(composite_levels) != _sorted_levels(composite_levels)
    ):
        raise BatchContractError("validation_summary has invalid composite_levels")
    expected_composite_count = len(composite_levels)
    composite_count = summary["composite_level_count"]
    if (
        isinstance(composite_count, bool)
        or not isinstance(composite_count, int)
        or composite_count != expected_composite_count
    ):
        raise BatchContractError(
            "validation_summary has an invalid composite_level_count"
        )
    if len(keys) == 1 and composite_levels:
        raise BatchContractError(
            "direct validation_summary cannot carry composite levels"
        )
    if len(keys) >= 2 and not composite_levels:
        raise BatchContractError(
            "composite validation_summary requires composite levels"
        )

    n_cells = summary["n_cells"]
    n_samples = summary["n_samples"]
    if (
        isinstance(n_cells, bool)
        or not isinstance(n_cells, int)
        or n_cells <= 0
        or isinstance(n_samples, bool)
        or not isinstance(n_samples, int)
        or n_samples < 2
    ):
        raise BatchContractError("validation_summary has invalid observation counts")
    for field in ("correction_mode", "correction_formula"):
        value = summary[field]
        if (
            not isinstance(value, str)
            or not value
            or value != value.strip()
        ):
            raise BatchContractError(
                f"validation_summary has invalid {field}"
            )
    return {
        "schema_version": int(summary["schema_version"]),
        "validated_before_reduction": True,
        "sample_constancy": {key: True for key in keys},
        "per_key_levels": {key: list(levels[key]) for key in keys},
        "key_level_counts": {key: int(counts[key]) for key in keys},
        "composite_levels": list(composite_levels),
        "composite_level_count": int(composite_count),
        "n_cells": int(n_cells),
        "n_samples": int(n_samples),
        "correction_mode": summary["correction_mode"],
        "correction_formula": summary["correction_formula"],
    }

def _h5_summary_value(node: Any) -> Any:
    """Decode an HDF5 summary while preserving configured underscore keys."""

    if hasattr(node, "keys"):
        return {
            str(name): _h5_summary_value(node[name])
            for name in node.keys()
        }
    try:
        value = node[()]
    except (OSError, TypeError, ValueError):
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, np.ndarray):
        values = value.tolist()
        if isinstance(values, list):
            return [
                item.decode("utf-8") if isinstance(item, bytes) else item
                for item in values
            ]
        value = values
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, str) and value[:1] in "[{":
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            pass
    return value


def read_h5ad_validation_summary(
    path: str,
    batch_keys: str | Sequence[str],
) -> dict[str, Any]:
    """Read a corrected H5AD summary without dropping underscore key names."""

    try:
        import h5py

        with h5py.File(path, "r") as handle:
            uns = handle.get("uns")
            contract_node = uns.get("batch_contract") if uns is not None else None
            if contract_node is None:
                raise ValueError("batch_contract is missing")
            if hasattr(contract_node, "keys"):
                summary_node = contract_node.get("validation_summary")
                if summary_node is None:
                    raise ValueError("validation_summary is missing")
                summary = _h5_summary_value(summary_node)
            else:
                contract = _h5_summary_value(contract_node)
                summary = (
                    contract.get("validation_summary")
                    if isinstance(contract, Mapping)
                    else None
                )
    except (ImportError, OSError, TypeError, UnicodeError, ValueError) as exc:
        raise BatchContractError(
            f"could not read corrected H5AD validation_summary: {path}"
        ) from exc
    return validate_batch_validation_summary(summary, batch_keys)


def serialize_batch_metadata(
    validation: BatchValidation | BatchComposite,
    *,
    method_id: str,
    model_id: str,
    include_tokens: bool = False,
) -> dict[str, Any]:
    """Return fresh JSON-safe run-owned metadata for a validated contract."""

    if isinstance(validation, BatchComposite):
        validation = validation.validation
    if not isinstance(validation, BatchValidation):
        raise BatchContractError("serialize_batch_metadata expects BatchValidation")
    method_id, model_id = _validate_method_model(method_id, model_id)
    fingerprint_payload = _batch_contract_payload(
        validation.keys,
        validation.scalarization,
        method_id,
        model_id,
    )
    fingerprint = hashlib.sha256(fingerprint_payload).hexdigest()
    scalarized_values = validation.scalarized_values
    sample_group_ids = validation.sample_group_ids
    sample_composite_values = validation.sample_composite_values
    sample_count = len(validation.sample_ids)
    if (
        not sample_group_ids
        or not sample_composite_values
        or len(sample_group_ids) != sample_count
        or len(sample_composite_values) != sample_count
    ):
        raise BatchContractError(
            "serialize_batch_metadata requires complete validated sample identity fields"
        )
    metadata: dict[str, Any] = {
        "contract_version": FINGERPRINT_VERSION,
        "token_version": validation.token_version,
        "encoding": validation.token_version,
        "encoding_version": validation.token_version,
        "ordered_keys": list(validation.keys),
        "keys": list(validation.keys),
        "key_count": validation.key_count,
        "scalarization": validation.scalarization,
        "per_key_levels": {
            key: list(validation.levels[key]) for key in validation.keys
        },
        "levels": {
            key: list(validation.levels[key]) for key in validation.keys
        },
        "composite_levels": list(validation.composite_levels),
        "composite_level_count": validation.composite_level_count,
        "composite_values": list(scalarized_values),
        "sample_composite_values": list(sample_composite_values),
        "method": method_id,
        "method_id": method_id,
        "model": model_id,
        "model_id": model_id,
        "fingerprint_payload": fingerprint_payload.decode("utf-8"),
        "fingerprint_payload_hex": fingerprint_payload.hex(),
        "fingerprint": fingerprint,
        "sample_column": validation.sample_column,
        "biological_column": validation.biological_column,
        "required_source_obs_columns": [validation.sample_column, *validation.keys],
        "reserved_obs_name": RESERVED_OBS_NAME,
        "reserved_obs_absent": True,
        "n_obs": validation.n_obs,
        "n_samples": validation.n_samples,
        "sample_ids": list(validation.sample_ids),
        "sample_group_ids": list(sample_group_ids),
        "sample_constancy": validation.sample_constancy,
        "estimable": validation.estimable,
        "near_unique_fraction": validation.near_unique_fraction,
        "design_rank": validation.design_rank,
        "design_columns": validation.design_columns,
        "composite_design_rank": validation.composite_design_rank,
        "composite_design_columns": validation.composite_design_columns,
    }
    if include_tokens:
        metadata["tokens"] = list(scalarized_values)
    return metadata
