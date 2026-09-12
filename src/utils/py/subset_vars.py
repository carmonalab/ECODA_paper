"""Lightweight shared subset predicates and sample-consistency checks.

This module intentionally depends only on the Python standard library,
NumPy, and pandas so metadata-only auditors can reuse preprocessing semantics
without importing Scanpy, rpy2, or R-backed preprocessing code.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd


def _normalize_subset_values(value, *, field):
    """Normalize scalar/list rule values without changing membership types."""
    if isinstance(value, Mapping):
        raise ValueError(f"{field} must be a scalar or sequence, not a mapping")
    if isinstance(value, (str, bytes)) or np.isscalar(value):
        values = [value]
    elif isinstance(value, (list, tuple, set, frozenset, np.ndarray, pd.Index, pd.Series)):
        if isinstance(value, np.ndarray):
            values = np.asarray(value).reshape(-1).tolist()
        elif isinstance(value, (pd.Index, pd.Series)):
            values = value.tolist()
        else:
            values = list(value)
    else:
        raise ValueError(f"{field} must be a scalar or sequence")
    if not values:
        raise ValueError(f"{field} must not be empty")
    return values


def evaluate_subset_mask(adata, subset_vars):
    """Evaluate declared row filters and return a mask indexed by obs names.

    Membership rules are intentionally exact.  Numeric comparison rules parse
    trimmed strings and reject malformed values by leaving those rows out;
    optional ``include_values`` exceptions are local to their comparison
    rule and use case-insensitive trimmed matching.
    """
    if subset_vars is None:
        subset_vars = {}
    if not isinstance(subset_vars, Mapping):
        raise ValueError("subset_vars must be a mapping of obs columns to rules")

    obs_names = pd.Index(adata.obs_names)
    mask = pd.Series(True, index=obs_names, dtype=bool)
    operators = {
        "<=": lambda values, threshold: values <= threshold,
        "<": lambda values, threshold: values < threshold,
        ">=": lambda values, threshold: values >= threshold,
        ">": lambda values, threshold: values > threshold,
    }

    for column, rule in subset_vars.items():
        if not isinstance(column, str) or not column:
            raise ValueError(f"subset_vars column names must be non-empty strings: {column!r}")
        if column not in adata.obs.columns:
            raise KeyError(f"subset_vars references missing obs column: {column}")
        if not isinstance(rule, Mapping):
            raise ValueError(f"subset_vars rule for {column!r} must be a mapping")
        if "op" not in rule:
            raise KeyError(f"subset_vars rule for {column!r} is missing operator 'op'")
        operator = rule["op"]
        if not isinstance(operator, str):
            raise ValueError(f"subset_vars operator for {column!r} must be a string")
        if operator not in {"in", "notin", *operators}:
            raise ValueError(
                f"unknown subset_vars operator {operator!r} for column {column!r}"
            )
        if "values" not in rule:
            raise KeyError(f"subset_vars rule for {column!r} is missing 'values'")
        values = _normalize_subset_values(rule["values"], field=f"{column}.values")
        include_values = None
        if "include_values" in rule:
            include_values = _normalize_subset_values(
                rule["include_values"], field=f"{column}.include_values"
            )
            if operator in {"in", "notin"}:
                raise ValueError(
                    f"include_values is only valid for numeric comparison rules: {column!r}"
                )

        series = adata.obs[column]
        if operator in {"in", "notin"}:
            col_mask = series.isin(values)
            if operator == "notin":
                col_mask = ~col_mask
        else:
            if len(values) != 1:
                raise ValueError(
                    f"comparison rule for {column!r} requires exactly one threshold"
                )
            threshold_value = values[0]
            if isinstance(threshold_value, (bool, np.bool_)):
                raise ValueError(
                    f"comparison threshold for {column!r} must be finite numeric"
                )
            threshold = pd.to_numeric(
                pd.Series([str(threshold_value).strip()]), errors="coerce"
            ).iloc[0]
            if pd.isna(threshold) or not np.isfinite(float(threshold)):
                raise ValueError(
                    f"comparison threshold for {column!r} must be finite numeric"
                )

            trimmed = series.astype("string").str.strip()
            numeric = pd.to_numeric(trimmed, errors="coerce")
            numeric_values = np.asarray(numeric.fillna(np.nan), dtype=float)
            finite = np.isfinite(numeric_values)
            col_mask = pd.Series(False, index=obs_names, dtype=bool)
            col_mask.iloc[:] = finite & operators[operator](numeric_values, float(threshold))
            if include_values is not None:
                include_text = []
                for value in include_values:
                    if value is None or pd.isna(value):
                        raise ValueError(
                            f"include_values for {column!r} must not contain missing values"
                        )
                    text = str(value).strip().casefold()
                    if not text:
                        raise ValueError(
                            f"include_values for {column!r} must not contain blank values"
                        )
                    include_text.append(text)
                col_mask |= trimmed.str.casefold().isin(include_text).fillna(False)

        col_mask = pd.Series(col_mask, index=obs_names, dtype=bool)
        mask &= col_mask

    mask.index = obs_names
    return mask


def assert_subset_sample_consistency(
    adata, mask, sample_col, context="subset"
):
    """Reject filters that split a configured sample across retained rows."""
    obs_names = pd.Index(adata.obs_names)
    if not isinstance(mask, pd.Series):
        mask = pd.Series(mask, index=obs_names)
    if not mask.index.equals(obs_names):
        raise ValueError(
            f"{context}: subset mask index must exactly match adata.obs_names"
        )
    if mask.isna().any():
        raise ValueError(f"{context}: subset mask contains missing values")
    try:
        mask = mask.astype(bool)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context}: subset mask must be boolean") from exc
    if sample_col not in adata.obs.columns:
        raise KeyError(
            f"{context}: sample column {sample_col!r} is missing; "
            f"available columns: {list(adata.obs.columns)}"
        )

    sample_values = adata.obs[sample_col]
    sample_ids = sample_values.astype("string")
    invalid = sample_values.isna() | sample_ids.str.strip().eq("")
    if bool(invalid.any()):
        invalid_rows = list(obs_names[invalid][:5])
        raise ValueError(
            f"{context}: sample column {sample_col!r} contains missing or blank IDs; "
            f"first invalid observations: {invalid_rows}"
        )

    retained = mask.to_numpy(dtype=bool)
    sample_array = sample_ids.to_numpy(dtype=object)
    retained_samples = []
    dropped_samples = []
    split_samples = []
    for sample in pd.unique(sample_array):
        sample_rows = sample_array == sample
        retained_rows = bool(np.any(sample_rows & retained))
        dropped_rows = bool(np.any(sample_rows & ~retained))
        if retained_rows:
            retained_samples.append(str(sample))
        if dropped_rows:
            dropped_samples.append(str(sample))
        if retained_rows and dropped_rows:
            split_samples.append(str(sample))
    if split_samples:
        raise ValueError(
            f"{context}: subset splits {len(split_samples)} sample(s) in "
            f"{sample_col!r}: {split_samples[:5]}"
        )
    return {
        "context": context,
        "sample_column": sample_col,
        "total_cells": int(len(mask)),
        "retained_cells": int(retained.sum()),
        "dropped_cells": int((~retained).sum()),
        "total_samples": int(len(pd.unique(sample_array))),
        "retained_samples": int(len(retained_samples)),
        "dropped_samples": int(len(dropped_samples)),
        "split_sample_count": 0,
        "retained_sample_ids": retained_samples,
        "dropped_sample_ids": dropped_samples,
        "split_sample_ids": [],
    }
