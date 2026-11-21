"""Whole-sample registry filtering and source-bound cohort identity.

Only approved technical ``majority_v1`` fields may be reduced within a
sample. Biological labels and other subset fields must be sample-consistent.
The resulting mask retains every cell of a selected raw sample identifier.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
import re
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd

try:
    from .export_h5ad_sample_metadata import (
        _canonical_batch_value, _feather_scalar, _majority_batch_winners,
        _normalise_batch_keys, _normalise_batch_metadata_policy,
    )
    from .h5ad_source_identity import read_obs_column
    from .subset_vars import assert_subset_sample_consistency, evaluate_subset_mask
except ImportError:  # direct source-tree scripts
    from export_h5ad_sample_metadata import (
        _canonical_batch_value, _feather_scalar, _majority_batch_winners,
        _normalise_batch_keys, _normalise_batch_metadata_policy,
    )
    from h5ad_source_identity import read_obs_column
    from subset_vars import assert_subset_sample_consistency, evaluate_subset_mask


COHORT_IDENTITY_VERSION = 1


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False)


def _sha256(value):
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _selection_configuration(configuration):
    # Adding a reporting outcome must not invalidate preprocessing/features.
    return {key: value for key, value in configuration.items() if key != "evaluation_labels"}


def _selection_columns(configuration):
    return list(dict.fromkeys([
        configuration["sample_column"], configuration["label_column"],
        *configuration["batch_keys"], *configuration["subset_vars"],
    ]))


def cohort_configuration(*, dataset, view, sample_column, label_column, batch_keys=None,
                         subset_vars=None, batch_metadata_policy=None,
                         evaluation_labels=None, min_cells_per_sample=500):
    """Return the exact selection/policy configuration bound to each artifact."""
    if subset_vars is None:
        subset_vars = {}
    if not isinstance(subset_vars, Mapping):
        raise ValueError("subset_vars must be an object")
    if evaluation_labels is None:
        evaluation_labels = {"primary": label_column}
    if (not isinstance(evaluation_labels, Mapping) or not evaluation_labels
            or any(not isinstance(name, str) or not name.strip()
                   or not isinstance(column, str) or not column.strip()
                   for name, column in evaluation_labels.items())):
        raise ValueError("evaluation_labels must map nonblank names to metadata columns")
    keys = _normalise_batch_keys(batch_keys, sample_column=sample_column,
                                 biological_column=label_column)
    policy = _normalise_batch_metadata_policy(batch_metadata_policy, dataset=dataset)
    if policy is not None and any(key not in keys for key in policy["majority_keys"]):
        raise ValueError("majority_keys must be configured technical batch keys")
    if isinstance(min_cells_per_sample, bool) or not isinstance(min_cells_per_sample, int) or min_cells_per_sample < 1:
        raise ValueError("min_cells_per_sample must be a positive integer")
    return {
        "schema_version": COHORT_IDENTITY_VERSION,
        "dataset": dataset, "view": view,
        "sample_column": sample_column, "label_column": label_column,
        "batch_keys": keys, "subset_vars": dict(subset_vars),
        "batch_metadata_policy": policy,
        "evaluation_labels": dict(evaluation_labels),
        "min_cells_per_sample": min_cells_per_sample,
    }


def sample_metadata_for_selection(obs, configuration):
    """Reduce technical winners while rejecting inconsistent selection/outcome fields."""
    sample_column = configuration["sample_column"]
    if len(obs) == 0:
        raise ValueError("cohort source has no observations")
    if not obs.index.is_unique:
        raise ValueError("cohort source observation IDs must be unique")
    fields = list(dict.fromkeys([*_selection_columns(configuration),
                                *configuration["evaluation_labels"].values()]))
    missing = [column for column in fields if column not in obs]
    if missing:
        raise KeyError(f"cohort metadata is missing configured columns: {missing}")
    data = SimpleNamespace(obs=obs, obs_names=obs.index)
    assert_subset_sample_consistency(
        data, pd.Series(True, index=obs.index), sample_column, "cohort source"
    )
    samples = obs[sample_column].astype("string")
    policy = configuration["batch_metadata_policy"]
    majority_keys = list(policy["majority_keys"]) if policy else []
    strict = [column for column in fields if column != sample_column and column not in majority_keys]
    if strict:
        cardinalities = obs[strict].groupby(samples, sort=False, observed=True).nunique(dropna=False)
        conflicts = [(str(sample), column) for sample, row in cardinalities.iterrows()
                     for column in strict if row[column] > 1]
        if conflicts:
            raise ValueError(f"cohort fields disagree within raw sample: {conflicts[:5]}")
    frame = obs.loc[~samples.duplicated(), fields].copy()
    frame[sample_column] = samples.loc[~samples.duplicated()].to_numpy()
    frame.index = pd.Index(frame[sample_column], name=None)
    if majority_keys:
        sample_ids = frame.index.tolist()
        counts = {sample: {key: {} for key in majority_keys} for sample in sample_ids}
        first = {sample: {key: {} for key in majority_keys} for sample in sample_ids}
        for values in obs[[sample_column, *majority_keys]].itertuples(index=False, name=None):
            sample = str(values[0])
            for key, value in zip(majority_keys, values[1:]):
                category, raw = _canonical_batch_value(
                    value, key=key, accepted_sentinel_values=policy["accepted_sentinel_values"]
                )
                counts[sample][key][category] = counts[sample][key].get(category, 0) + 1
                first[sample][key].setdefault(category, raw)
        winners = _majority_batch_winners(sample_ids, majority_keys, counts, first)
        for key in majority_keys:
            frame[key] = winners[key]
    return frame


def canonical_sample_ids(raw_ids):
    values = [f"g{sample}" if re.match(r"^\d", str(sample)) else str(sample).replace("-", "_")
              for sample in raw_ids]
    if len(set(values)) != len(values):
        raise ValueError("raw sample IDs collide after canonical Sample standardization")
    return values


def _observation_digest(rows):
    digest = hashlib.sha256(b"ecoda_obs_sample_v1\n")
    for cell, sample in rows:
        digest.update((_json([str(cell), str(sample)]) + "\n").encode("utf-8"))
    return digest.hexdigest()


def _canonical_category_map(series, *, key, accepted_sentinel_values):
    """Map raw cell values in a series to canonical categories."""
    canon_cache = {}

    def canonicalize(value):
        if isinstance(value, np.generic):
            value = value.item()
        if value is None or value is pd.NA or value is pd.NaT:
            return _canonical_batch_value(value, key=key, accepted_sentinel_values=accepted_sentinel_values)[0]
        if isinstance(value, float) and not np.isfinite(value):
            return _canonical_batch_value(value, key=key, accepted_sentinel_values=accepted_sentinel_values)[0]
        cache_key = (type(value), value)
        if cache_key not in canon_cache:
            canon_cache[cache_key] = _canonical_batch_value(
                value, key=key, accepted_sentinel_values=accepted_sentinel_values
            )[0]
        return canon_cache[cache_key]

    return [canonicalize(value) for value in series.tolist()]


def resolve_cohort_selection(obs, configuration):
    """Return a cell selection mask, resolved sample metadata, and an auditable identity."""
    sample_column = configuration["sample_column"]
    frame = sample_metadata_for_selection(obs, configuration)
    selected = evaluate_subset_mask(
        SimpleNamespace(obs=frame, obs_names=frame.index), configuration["subset_vars"]
    )
    retained = frame.index[selected].tolist()
    sample_series = obs[sample_column].astype("string")
    sample_mask = sample_series.isin(retained)
    audit = assert_subset_sample_consistency(
        SimpleNamespace(obs=obs, obs_names=obs.index), sample_mask, sample_column,
        f"{configuration['dataset']} / {configuration['view']}",
    )
    mask = sample_mask.copy()
    policy = configuration["batch_metadata_policy"]
    majority_cell_filter = policy is not None and policy.get("cell_filter") == "majority_only"
    counts = sample_series.value_counts(sort=False)
    if majority_cell_filter and retained:
        retained_positions = mask.to_numpy(dtype=bool)
        retained_samples = sample_series.iloc[retained_positions].astype(str).tolist()
        for key in policy["majority_keys"]:
            winner_cats = _canonical_category_map(
                frame.loc[retained, key], key=key,
                accepted_sentinel_values=policy["accepted_sentinel_values"],
            )
            winner_by_sample = dict(zip([str(s) for s in retained], winner_cats))
            cell_cats = _canonical_category_map(
                obs.loc[mask, key], key=key,
                accepted_sentinel_values=policy["accepted_sentinel_values"],
            )
            keep_flags = [
                cell_cats[idx] == winner_by_sample[retained_samples[idx]]
                for idx in range(len(retained_samples))
            ]
            new_mask_arr = mask.to_numpy(dtype=bool).copy()
            new_mask_arr[retained_positions] = keep_flags
            mask = pd.Series(new_mask_arr, index=obs.index, dtype=bool)
            retained_positions = new_mask_arr
            retained_samples = sample_series.iloc[retained_positions].astype(str).tolist()
    if majority_cell_filter:
        post_counts = sample_series.loc[mask].value_counts(sort=False)
        retained_after_mask = (
            sample_series.loc[mask].astype(str).drop_duplicates().tolist()
        )
        retained_after_mask_set = set(retained_after_mask)
        retained = retained_after_mask + [
            sample for sample in retained if sample not in retained_after_mask_set
        ]
        after_count_filter = [
            sample for sample in retained
            if int(post_counts.get(sample, 0)) >= configuration["min_cells_per_sample"]
        ]
        audit["retained_sample_ids"] = list(retained)
        audit["cell_filter"] = "majority_only"
        audit["pre_cell_filter_retained_cells"] = int(sample_mask.sum())
        audit["retained_cells"] = int(mask.sum())
        audit["dropped_cells"] = int(len(mask) - mask.sum())
        audit["dropped_non_majority_cells"] = int(sample_mask.sum() - mask.sum())
    else:
        post_counts = None
        after_count_filter = [sample for sample in retained
                              if counts[sample] >= configuration["min_cells_per_sample"]]
    membership = {
        "configuration_sha256": _sha256(_selection_configuration(configuration)),
        "source_observations_sha256": _observation_digest(zip(obs.index, obs[sample_column])),
        "source_sample_counts": {sample: int(counts[sample]) for sample in frame.index},
        "source_sample_metadata": [
            {column: _feather_scalar(value) for column, value in row.items()}
            for row in frame[_selection_columns(configuration)].to_dict("records")
        ],
        "retained_raw_sample_ids": retained,
        "dropped_raw_sample_ids": frame.index[~selected].tolist(),
        "retained_after_cellcount_raw_ids": after_count_filter,
        "retained_after_cellcount_ids": canonical_sample_ids(after_count_filter),
    }
    if majority_cell_filter:
        membership["post_cell_filter_sample_counts"] = {
            sample: int(post_counts.get(sample, 0)) for sample in retained
        }
    audit["cohort_configuration"] = configuration
    audit["cohort_membership"] = membership
    audit["cohort_membership_sha256"] = _sha256(membership)
    return mask, frame, audit


def input_file_identities(input_names, input_dir):
    """Bind raw inputs without opening expression arrays or using mutable caches."""
    names = input_names if isinstance(input_names, list) else [input_names]
    identities = []
    for name in names:
        path = Path(input_dir) / name
        if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"cohort source is missing, empty, or a symlink: {path}")
        digest = hashlib.md5()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        identities.append({"path": str(path.resolve()), "size": path.stat().st_size,
                           "md5": digest.hexdigest()})
    return identities


def cohort_record(audit, source_identities, output_obs, output_sample_column="Sample"):
    membership = audit["cohort_membership"]
    actual_ids = output_obs[output_sample_column].astype("string").drop_duplicates().tolist()
    if actual_ids != membership["retained_after_cellcount_ids"]:
        raise ValueError("processed Sample universe/order differs from the selected cohort")
    configuration = _selection_configuration(audit["cohort_configuration"])
    return {
        "schema_version": COHORT_IDENTITY_VERSION,
        "configuration_json": _json(configuration), "configuration_sha256": _sha256(configuration),
        "source_identity_json": _json(source_identities),
        "membership_json": _json(membership), "membership_sha256": _sha256(membership),
        "output_sample_column": output_sample_column,
        "output_observations_sha256": _observation_digest(zip(output_obs.index, output_obs[output_sample_column])),
    }


def validate_cohort_h5ad_path(path, expected_configuration, expected_sources):
    """Reject missing/stale cohort identity; read only its provenance and obs IDs."""
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"cohort artifact is missing, empty, or a symlink: {path}")
    with h5py.File(path, "r") as handle:
        if "uns/cohort_selection" not in handle:
            raise ValueError("H5AD lacks cohort identity; use a new derivative destination")
        record = handle["uns/cohort_selection"]

        def scalar(name):
            if name not in record:
                raise ValueError(f"H5AD cohort identity lacks {name}")
            value = record[name][()]
            return value.decode("utf-8") if isinstance(value, bytes) else value

        if scalar("schema_version") != COHORT_IDENTITY_VERSION:
            raise ValueError("unsupported H5AD cohort identity version")
        configuration = json.loads(scalar("configuration_json"))
        if configuration != json.loads(_json(_selection_configuration(expected_configuration))) or scalar("configuration_sha256") != _sha256(configuration):
            raise ValueError("H5AD cohort selection/policy differs from the current configuration")
        if json.loads(scalar("source_identity_json")) != expected_sources:
            raise ValueError("H5AD cohort source identity differs from the current raw inputs")
        membership = json.loads(scalar("membership_json"))
        if scalar("membership_sha256") != _sha256(membership) or membership["configuration_sha256"] != _sha256(configuration):
            raise ValueError("H5AD cohort membership identity is inconsistent")
        obs = handle["obs"]
        index_name = obs.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode("utf-8")
        sample_column = scalar("output_sample_column")
        n_obs = obs[index_name].shape[0]
        ordered = []
        seen = set()
        digest = hashlib.sha256(b"ecoda_obs_sample_v1\n")
        for start in range(0, n_obs, 1024 * 1024):
            stop = min(n_obs, start + 1024 * 1024)
            cells = read_obs_column(obs, index_name, start, stop)
            samples = read_obs_column(obs, sample_column, start, stop)
            for cell, sample in zip(cells, samples):
                digest.update((_json([str(cell), str(sample)]) + "\n").encode("utf-8"))
                if sample not in seen:
                    seen.add(sample)
                    ordered.append(str(sample))
        if ordered != membership["retained_after_cellcount_ids"] or digest.hexdigest() != scalar("output_observations_sha256"):
            raise ValueError("H5AD obs disagrees with its cohort membership identity")
    return True


if __name__ == "__main__":
    import argparse
    try:
        from .datasets_io import read_datasets_json
    except ImportError:
        from datasets_io import read_datasets_json

    parser = argparse.ArgumentParser(description="Validate a Stage 3 cohort without reading expression arrays")
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--view", required=True, choices=("benchmark_analysis", "batch_effect_uncorrected", "batch_effect_corrected"))
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--path", required=True)
    args = parser.parse_args()
    try:
        entry = read_datasets_json(args.config, view=args.view)[args.dataset]
        enabled = entry["use_for_benchmark"] if args.view == "benchmark_analysis" else entry["use_for_batch_effect"]
        if not enabled:
            raise ValueError(f"dataset is not enabled for {args.view}: {args.dataset}")
        view = entry["views"][args.view]
        configuration = cohort_configuration(
            dataset=args.dataset, view=args.view, sample_column=entry["sample_col"],
            label_column=entry["label_col"], batch_keys=entry["batch_col"],
            subset_vars=view["subset_vars"], batch_metadata_policy=entry["batch_metadata_policy"],
            evaluation_labels=entry["evaluation_labels"],
        )
        sources = input_file_identities(view["input_file"], args.input_dir)
        validate_cohort_h5ad_path(args.path, configuration, sources)
        with h5py.File(args.path, "r") as handle:
            record = handle["uns/cohort_selection"]
            membership = json.loads(record["membership_json"].asstr()[()])
            print(_json({"status": "COHORT_VALIDATED", "configuration_sha256": _sha256(_selection_configuration(configuration)),
                         "membership_sha256": _sha256(membership),
                         "n_samples": len(membership["retained_after_cellcount_ids"])}))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
