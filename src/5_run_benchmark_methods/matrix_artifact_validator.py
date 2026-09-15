"""Validate the exact selected Pipeline 5 artifacts before synchronization."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.py.batch_contract import (  # noqa: E402
    build_batch_contract_identity,
    batch_correction_spec_for_keys,
    validate_batch_validation_summary,
)
from src.utils.py.datasets_io import read_datasets_json  # noqa: E402
from src.utils.py.h5ad_source_identity import (  # noqa: E402
    load_source_identity,
    read_h5ad_sample_ids,
    resolve_h5ad_path,
    verify_source_identity,
)
from src.utils.py.benchmark_h5ad_contract import (  # noqa: E402
    validate_batch_contract_identity as _legacy_validate_batch_contract_identity,
)
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
_CORRECTED_METHOD_IDENTITIES = {
    "prepare_pseudobulk": (
        "Pseudobulk",
        "pseudobulk_limma_fixed_effects_v1",
    ),
    "pseudobulk": (
        "Pseudobulk",
        "pseudobulk_limma_fixed_effects_v1",
    ),
    "composition": (
        "ECODA_authors_HR",
        "limma_fixed_effects_v1",
    ),
    "gloscope": ("GloScope", "embedding_consumer_harmony_v1"),
    "mrvi": ("MrVI", "mrvi_composite_v1"),
    "pilot": ("PILOT", "embedding_consumer_harmony_v1"),
    "qot": ("QOT", "embedding_consumer_harmony_v1"),
}

_NEW_LIMMA_MODEL_IDS = frozenset(
    {"limma_fixed_effects_v1", "pseudobulk_limma_fixed_effects_v1"}
)
_HISTORICAL_LIMMA_MODEL_IDS = frozenset(
    {"ecoda_additive_random_intercepts_v1", "pseudobulk_composite_v1"}
)
_FORBIDDEN_LIMMA_CORRECTION_TOKENS = (
    "__ecoda_batch_combined_v1",
    "additive_random_intercepts",
    "pseudobulk_composite",
    "limma::lmFit",
    "remove technical contribution",
    "(1 |",
)

FINAL_BATCH_METHODS = (
    "prepare_pseudobulk",
    "pseudobulk",
    "gloscope",
    "composition",
    "mrvi",
    "pilot",
    "qot",
)

FINAL_UNCORRECTED_DATASETS = (
    "Covid19_PBMC",
    "Diabetes",
    "Joanito",
    "Lung",
    "Kidney_KPMP_full",
)


def _identity_value(identity: dict, fields: tuple[str, ...], label: str):
    present = [field for field in fields if field in identity]
    if not present:
        raise ValueError(f"{label} is missing")
    value = identity[present[0]]
    if any(identity[field] != value for field in present[1:]):
        raise ValueError(f"{label} aliases disagree")
    return value


def _identity_model(identity: object) -> str | None:
    if not isinstance(identity, dict):
        return None
    value = identity.get("model_id", identity.get("model"))
    return value if isinstance(value, str) else None


def _validate_new_limma_identity(
    identity: object,
    label: str,
    *,
    require_metadata: bool,
    require_summary: bool,
) -> dict:
    if not isinstance(identity, dict):
        raise ValueError(f"{label} must be a JSON object")
    keys = tuple(
        _identity_value(
            identity,
            (
                "ordered_source_keys",
                "ordered_keys",
                "source_keys",
                "batch_keys",
                "keys",
            ),
            f"{label} ordered source keys",
        )
    )
    if (
        not keys
        or any(not isinstance(key, str) or not key or key != key.strip() for key in keys)
        or len(set(keys)) != len(keys)
    ):
        raise ValueError(f"{label} has invalid ordered source keys")
    method_id = _identity_value(
        identity, ("method_id", "method", "method_policy"), f"{label} method"
    )
    model_id = _identity_value(
        identity, ("model_id", "model", "model_policy"), f"{label} model"
    )
    expected_model = {
        "Pseudobulk": "pseudobulk_limma_fixed_effects_v1",
        "ECODA_authors_HR": "limma_fixed_effects_v1",
        "ECODA_seuratres_2": "limma_fixed_effects_v1",
        "ECODA_authors_HR_NULL": "limma_fixed_effects_v1",
    }.get(method_id)
    if expected_model is None or model_id != expected_model:
        raise ValueError(
            f"{label} must use the active fixed-effect model for its method"
        )
    expected = build_batch_contract_identity(
        keys,
        sample_column="Sample",
        method_id=method_id,
        model_id=model_id,
    )
    for field in (
        "contract_version",
        "token_version",
        "ordered_source_keys",
        "scalarization",
        "method_id",
        "model_id",
        "required_source_obs_columns",
        "reserved_obs_name",
        "reserved_obs_absent",
        "fingerprint",
    ):
        if identity.get(field) != expected[field]:
            raise ValueError(f"{label} has mismatched {field}")

    metadata_fields = (
        "effective_batch_keys",
        "non_estimable_batch_keys",
        "correction_state",
        "correction_mode",
        "correction_formula",
        "fixed_effect_aliases",
        "correction_design_formula",
        "design_rank",
        "design_columns",
        "design_residual_df",
    )
    has_metadata = any(field in identity for field in metadata_fields)
    if require_metadata and not all(field in identity for field in metadata_fields):
        missing = [field for field in metadata_fields if field not in identity]
        raise ValueError(f"{label} is missing fixed-effect metadata: {missing}")
    if has_metadata:
        effective = identity.get("effective_batch_keys")
        non_estimable = identity.get("non_estimable_batch_keys")
        if (
            not isinstance(effective, (list, tuple))
            or not isinstance(non_estimable, (list, tuple))
            or any(not isinstance(key, str) for key in (*effective, *non_estimable))
            or len(set(effective)) != len(effective)
            or len(set(non_estimable)) != len(non_estimable)
            or list(effective) != [key for key in keys if key in effective]
            or list(non_estimable) != [key for key in keys if key in non_estimable]
            or set(effective).union(non_estimable) != set(keys)
            or set(effective).intersection(non_estimable)
        ):
            raise ValueError(f"{label} has invalid effective/non-estimable keys")
        expected_mode, expected_formula = batch_correction_spec_for_keys(
            method_id,
            keys,
            effective_batch_keys=effective,
            non_estimable_batch_keys=non_estimable,
        )
        expected_state = "BATCH_CORRECTION" if effective else "NO_CORRECTION"
        if (
            identity["correction_state"] != expected_state
            or identity["correction_mode"] != expected_mode
            or identity["correction_formula"] != expected_formula
        ):
            raise ValueError(f"{label} has the wrong fixed-effect correction policy")
        expected_aliases = {
            key: f"batch_key_{index}"
            for index, key in enumerate(keys, start=1)
            if key in effective
        }
        if identity["fixed_effect_aliases"] != expected_aliases:
            raise ValueError(f"{label} has the wrong fixed-effect aliases")
        expected_design = (
            "~1"
            if not effective
            else "~1 + " + " + ".join(expected_aliases.values())
        )
        design = identity["correction_design_formula"]
        if (
            not isinstance(design, str)
            or "".join(design.split()) != "".join(expected_design.split())
        ):
            raise ValueError(f"{label} has the wrong separate-covariate design")
        for field in ("design_rank", "design_columns", "design_residual_df"):
            value = identity[field]
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or value < 0
            ):
                raise ValueError(f"{label} has invalid {field}")
        if identity["design_rank"] != identity["design_columns"]:
            raise ValueError(f"{label} design rank is not full rank")
        if identity["design_residual_df"] <= 0:
            raise ValueError(f"{label} design has no residual degrees of freedom")
        if identity["correction_mode"] in {
            "limma_fixed_effects",
            "limma_fixed_effects_pseudobulk",
        } and any(
            token in identity["correction_formula"]
            for token in _FORBIDDEN_LIMMA_CORRECTION_TOKENS
        ):
            raise ValueError(f"{label} advertises a historical/scalarized correction")

    if "validation_summary" in identity:
        if require_summary and identity["validation_summary"] is None:
            raise ValueError(f"{label} validation_summary is missing")
        if identity["validation_summary"] is not None:
            try:
                summary = validate_batch_validation_summary(
                    identity["validation_summary"], keys
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{label} has an invalid validation_summary") from exc
            if has_metadata and (
                summary["correction_mode"] != identity["correction_mode"]
                or summary["correction_formula"] != identity["correction_formula"]
            ):
                raise ValueError(f"{label} summary correction policy disagrees")
    elif require_summary:
        raise ValueError(f"{label} is missing validation_summary")
    return identity


def validate_batch_contract_identity(
    expected_batch_contract=None,
    batch_contract=None,
    *,
    source_obs_columns=None,
    reserved_absent=None,
    require_recorded=False,
    label="batch contract",
    require_summary=None,
    historical_compatibility=False,
    require_effective_metadata=False,
):
    """Validate active limma identities, with an explicit historical branch."""
    source = (
        expected_batch_contract
        if expected_batch_contract is not None
        else batch_contract
    )
    model = _identity_model(source)
    if model in _HISTORICAL_LIMMA_MODEL_IDS:
        if not historical_compatibility:
            raise ValueError(
                f"{label} uses a historical correction identity; "
                "pass historical_compatibility only for read-only legacy validation"
            )
        # Historical identities may carry the old correction formula.  Strip
        # only that optional summary before the legacy structural reader; the
        # explicit compatibility flag prevents this branch from validating a
        # variant-qualified/new-policy output.
        historical_expected = (
            dict(expected_batch_contract)
            if isinstance(expected_batch_contract, dict)
            else expected_batch_contract
        )
        historical_recorded = (
            dict(batch_contract)
            if isinstance(batch_contract, dict)
            else batch_contract
        )
        for identity in (historical_expected, historical_recorded):
            if isinstance(identity, dict):
                identity.pop("validation_summary", None)
        return _legacy_validate_batch_contract_identity(
            historical_expected,
            historical_recorded,
            source_obs_columns=source_obs_columns,
            reserved_absent=reserved_absent,
            require_recorded=require_recorded,
            label=label,
            require_summary=False,
        )
    if model in _NEW_LIMMA_MODEL_IDS:
        required_summary = bool(require_summary)
        if expected_batch_contract is not None:
            _validate_new_limma_identity(
                expected_batch_contract,
                f"{label} expected identity",
                require_metadata=False,
                require_summary=False,
            )
        if batch_contract is not None:
            metadata_fields = (
                "effective_batch_keys",
                "non_estimable_batch_keys",
                "correction_state",
                "correction_mode",
                "correction_formula",
                "fixed_effect_aliases",
                "correction_design_formula",
                "design_rank",
                "design_columns",
                "design_residual_df",
            )
            recorded_metadata = any(
                field in batch_contract for field in metadata_fields
            )
            _validate_new_limma_identity(
                batch_contract,
                f"{label} recorded identity",
                require_metadata=(
                    require_effective_metadata or recorded_metadata
                ),
                require_summary=required_summary,
            )
        elif require_recorded:
            raise ValueError(f"{label} is missing recorded corrected batch identity")
        if expected_batch_contract is not None and batch_contract is not None:
            for field in (
                "ordered_source_keys",
                "scalarization",
                "method_id",
                "model_id",
                "fingerprint",
            ):
                if expected_batch_contract.get(field) != batch_contract.get(field):
                    raise ValueError(f"{label} source/config identity does not match")
        if source_obs_columns is not None:
            required_columns = expected_batch_contract or batch_contract
            required = required_columns["required_source_obs_columns"]
            missing = [column for column in required if column not in source_obs_columns]
            if missing:
                raise ValueError(f"{label} source obs columns are missing: {missing}")
        if reserved_absent is not None and reserved_absent is not True:
            raise ValueError(f"{label} source obs contains a reserved temporary column")
        return batch_contract or expected_batch_contract
    return _legacy_validate_batch_contract_identity(
        expected_batch_contract,
        batch_contract,
        source_obs_columns=source_obs_columns,
        reserved_absent=reserved_absent,
        require_recorded=require_recorded,
        label=label,
        require_summary=require_summary,
    )

def _normalized_absolute_path(value: str, label: str) -> Path:
    if not isinstance(value, str) or not value or "\n" in value or "\t" in value:
        raise ValueError(f"{label} must be a non-empty absolute path")
    candidate = Path(os.path.normpath(value))
    if not candidate.is_absolute():
        raise ValueError(f"{label} must be an absolute path: {value}")
    return candidate


def _variant_root_identity(root: Path) -> str | None:
    normalized = Path(os.path.normpath(str(root)))
    if not normalized.is_absolute():
        return None
    if normalized.name == "uncorrected_final" and normalized.parent.name == "batch_effect":
        return "uncorrected_final"
    if normalized.name == "corrected_final" and normalized.parent.name == "batch_effect":
        return "corrected_final"
    if (
        normalized.name == "recovery_35row"
        and normalized.parent.name == "corrected_final"
        and normalized.parent.parent.name == "batch_effect"
    ):
        return "corrected_final/recovery_35row"
    return None


def _validate_variant_root_binding(
    root: Path, analysis_variant: str
) -> None:
    """Bind a variant to one physical scratch/NAS root identity.

    The historical direct ``corrected_final`` root remains valid when no
    replacement-root metadata is present.  The disjoint recovery root is
    intentionally fail-closed: a path below it is accepted only when the
    scheduler-bound root version/identity says exactly
    ``corrected_final/recovery_35row``.  This keeps direct historical
    artifacts immutable while preventing a matrix run from mixing roots.
    """

    if analysis_variant not in {"final", "corrected_final"}:
        raise ValueError(f"unknown analysis variant: {analysis_variant}")
    root = _normalized_absolute_path(str(root), "variant validation root")
    actual_identity = _variant_root_identity(root)
    if actual_identity is None:
        raise ValueError(
            f"{analysis_variant} artifact validation requires a supported "
            "batch_effect variant root"
        )

    version_values: list[tuple[str, str]] = []
    for field in (
        "ANALYSIS_ROOT_VERSION",
        "ECODA_STAGE5_CORRECTED_FINAL_ROOT_VERSION",
    ):
        value = os.environ.get(field, "")
        if value:
            if "\n" in value or "\t" in value:
                raise ValueError(f"{field} contains a record delimiter")
            version_values.append((field, value))
    distinct_versions = {value for _, value in version_values}
    if distinct_versions - {"recovery_35row"}:
        raise ValueError(
            "corrected-final root version must be recovery_35row: "
            + ", ".join(sorted(distinct_versions - {"recovery_35row"}))
        )
    if len(distinct_versions) > 1:
        raise ValueError("corrected-final root version declarations disagree")

    declared_identity = os.environ.get("ANALYSIS_ROOT_IDENTITY", "")
    if declared_identity and ("\n" in declared_identity or "\t" in declared_identity):
        raise ValueError("ANALYSIS_ROOT_IDENTITY contains a record delimiter")
    if analysis_variant == "final":
        if distinct_versions:
            raise ValueError(
                "ANALYSIS_ROOT_VERSION is only valid for corrected_final"
            )
        if declared_identity and declared_identity != "uncorrected_final":
            raise ValueError(
                "final analysis cannot use a corrected-final root identity"
            )
        expected_identity = "uncorrected_final"
    else:
        if declared_identity not in {
            "",
            "corrected_final",
            "corrected_final/recovery_35row",
        }:
            raise ValueError(
                f"invalid ANALYSIS_ROOT_IDENTITY: {declared_identity}"
            )
        if distinct_versions:
            expected_identity = "corrected_final/recovery_35row"
            if (
                declared_identity
                and declared_identity != expected_identity
            ):
                raise ValueError(
                    "corrected-final root version and identity disagree"
                )
        elif declared_identity:
            expected_identity = declared_identity
        else:
            expected_identity = "corrected_final"

    env_root_paths: list[tuple[str, Path]] = []
    for field in ("ANALYSIS_ROOT", "ANALYSIS_NAS_ROOT"):
        value = os.environ.get(field, "")
        if value:
            env_root_paths.append(
                (field, _normalized_absolute_path(value, field))
            )
    env_identities: list[tuple[str, str]] = []
    for field, env_path in env_root_paths:
        env_identity = _variant_root_identity(env_path)
        if env_identity is None:
            raise ValueError(f"{field} is not a supported variant root: {env_path}")
        env_identities.append((field, env_identity))
    if (
        analysis_variant == "corrected_final"
        and not distinct_versions
        and not declared_identity
        and any(identity == "corrected_final/recovery_35row"
                for _, identity in env_identities)
    ):
        raise ValueError(
            "replacement corrected_final root requires a bound root identity"
        )
    if env_identities and len({identity for _, identity in env_identities}) > 1:
        raise ValueError("ANALYSIS_ROOT and ANALYSIS_NAS_ROOT root identities disagree")
    if env_identities and env_identities[0][1] != expected_identity:
        raise ValueError(
            "variant root metadata disagrees with ANALYSIS_ROOT/ANALYSIS_NAS_ROOT"
        )
    for field, env_path in env_root_paths:
        if field == "ANALYSIS_ROOT" and env_path != root:
            raise ValueError(
                "variant artifact root mixing is not allowed: "
                f"expected {env_path} but found {root}"
            )
    if actual_identity != expected_identity:
        if (
            actual_identity == "corrected_final/recovery_35row"
            and analysis_variant == "corrected_final"
        ):
            raise ValueError(
                "replacement corrected_final root requires a bound root identity"
            )
        raise ValueError(
            f"{analysis_variant} artifact validation requires the "
            f"{expected_identity} root"
        )


def _validate_variant_root(root: Path, analysis_variant: str) -> None:
    _validate_variant_root_binding(root, analysis_variant)

def batch_artifact_stem(
    ds: str,
    batch_pass: str | None,
    analysis_variant: str | None = None,
) -> str:
    """Return the single variant-qualified stem used by every batch artifact."""
    if not batch_pass:
        if analysis_variant:
            raise ValueError("analysis variant requires a batch-effect pass")
        return ds
    if batch_pass not in {"uncorrected", "corrected"}:
        raise ValueError(f"invalid batch pass: {batch_pass}")
    variant = analysis_variant or ""
    if variant not in {"", "final", "corrected_final"}:
        raise ValueError(f"unknown analysis variant: {variant}")
    if variant == "final":
        if batch_pass != "uncorrected":
            raise ValueError("final analysis variant requires uncorrected batch pass")
        return f"{ds}_batch_effect_uncorrected_final"
    if variant == "corrected_final":
        if batch_pass != "corrected":
            raise ValueError(
                "corrected_final analysis variant requires corrected batch pass"
            )
        return f"{ds}_batch_effect_corrected_final"
    return f"{ds}_batch_effect_{batch_pass}"


def _corrected_method_identity(label: str) -> tuple[str, str]:
    try:
        return _CORRECTED_METHOD_IDENTITIES[label]
    except KeyError as exc:
        raise ValueError(
            f"unsupported corrected batch method for identity: {label}"
        ) from exc


def _build_corrected_batch_contract(
    config_entries: dict,
    ds: str,
    view: str,
    label: str,
) -> dict:
    entry = config_entries.get(ds)
    if not isinstance(entry, dict):
        raise ValueError(f"dataset {ds!r} is missing from the selected config")
    raw_keys = entry.get("batch_col")
    if raw_keys is None:
        raise ValueError(
            f"corrected batch config is missing columns.batch for {ds}/{view}"
        )
    method_id, model_id = _corrected_method_identity(label)
    try:
        return build_batch_contract_identity(
            raw_keys,
            sample_column="Sample",
            method_id=method_id,
            model_id=model_id,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"invalid corrected batch config for {ds}/{view}/{label}: {exc}"
        ) from exc


_RUNTIME_METADATA_FIELDS = frozenset(
    {
        "schema_version",
        "artifact_path",
        "artifact_md5",
        "dataset",
        "method",
        "time_secs",
        "mem_GB",
    }
)
_RUNTIME_METADATA_IDENTITY_FIELDS = frozenset({"batch_contract"})
_RUNTIME_CHECKSUM_FIELDS = ("MD5", "SIZE", "PATH")


def _reject_nonfinite_json_constant(value):
    raise ValueError(f"runtime metadata contains non-finite JSON value: {value}")


def _reject_duplicate_json_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"runtime metadata has duplicate key: {key}")
        result[key] = value
    return result


def _runtime_number(value, field, *, allow_none=False):
    if value is None and allow_none:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"runtime metadata has invalid {field}")
    if not np.isfinite(float(value)) or value < 0:
        raise ValueError(f"runtime metadata has invalid {field}")


def _read_feather_batch_contract(
    path: Path,
    *,
    require_summary: bool = True,
    historical_compatibility: bool = False,
) -> dict | None:
    """Read and verify the corrected identity beside one Feather artifact."""
    metadata_path = Path(f"{path}.runtime.json")
    if not metadata_path.is_file():
        return None
    try:
        _full_checksum(metadata_path)
    except (OSError, ValueError) as exc:
        raise ValueError(
            f"Feather runtime metadata checksum is invalid: {metadata_path}"
        ) from exc
    try:
        payload = json.loads(
            metadata_path.read_text(encoding="utf-8"),
            parse_constant=_reject_nonfinite_json_constant,
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(
            f"Feather runtime metadata is malformed: {metadata_path}"
        ) from exc
    required_fields = _RUNTIME_METADATA_FIELDS | _RUNTIME_METADATA_IDENTITY_FIELDS
    if not isinstance(payload, dict) or set(payload) != required_fields:
        raise ValueError(
            f"Feather runtime metadata has an invalid schema: {metadata_path}"
        )
    if type(payload["schema_version"]) is not int or payload["schema_version"] != 1:
        raise ValueError(
            f"Feather runtime metadata has an invalid schema_version: {metadata_path}"
        )
    if payload["artifact_path"] != str(path):
        raise ValueError(
            f"Feather runtime metadata has the wrong artifact_path: {metadata_path}"
        )
    artifact_checksum = _read_checksum_sidecar(path)
    if (
        not isinstance(payload["artifact_md5"], str)
        or payload["artifact_md5"] != artifact_checksum["MD5"]
    ):
        raise ValueError(
            f"Feather runtime metadata has the wrong artifact_md5: {metadata_path}"
        )
    if (
        not isinstance(payload["dataset"], str)
        or not payload["dataset"]
        or not isinstance(payload["method"], str)
        or not payload["method"]
    ):
        raise ValueError(
            f"Feather runtime metadata has invalid dataset/method: {metadata_path}"
        )
    _runtime_number(payload["time_secs"], "time_secs")
    _runtime_number(payload["mem_GB"], "mem_GB", allow_none=True)
    identity = payload["batch_contract"]
    if not isinstance(identity, dict):
        raise ValueError(
            f"Feather runtime metadata batch_contract is not an object: "
            f"{metadata_path}"
        )
    validate_batch_contract_identity(
        identity,
        identity,
        require_recorded=True,
        require_summary=require_summary,
        historical_compatibility=historical_compatibility,
        require_effective_metadata=True,
        label=f"Feather runtime metadata {metadata_path}",
    )
    return identity

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
    *,
    expected_batch_contract=None,
    batch_contract=None,
    require_runtime_batch_contract: bool = False,
    require_corrected_summary: bool = True,
    historical_compatibility: bool = False,
) -> None:
    if not paths:
        raise ValueError(f"missing/invalid {description}: []")
    if expected_batch_contract is not None or batch_contract is not None:
        validate_batch_contract_identity(
            expected_batch_contract,
            batch_contract,
            require_recorded=(
                expected_batch_contract is not None and batch_contract is not None
            ),
            require_summary=False,
            historical_compatibility=historical_compatibility,
            label=description,
        )
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
                if require_runtime_batch_contract:
                    runtime_batch_contract = _read_feather_batch_contract(
                        path,
                        require_summary=require_corrected_summary,
                        historical_compatibility=historical_compatibility,
                    )
                    if runtime_batch_contract is None:
                        raise ValueError(
                            "missing recorded corrected batch contract identity"
                        )
                    if batch_contract is not None:
                        validate_batch_contract_identity(
                            batch_contract,
                            runtime_batch_contract,
                            require_recorded=True,
                            require_summary=require_corrected_summary,
                            historical_compatibility=historical_compatibility,
                            label=f"{description} embedded identity",
                        )
                    validate_batch_contract_identity(
                        expected_batch_contract,
                        runtime_batch_contract,
                        require_recorded=expected_batch_contract is not None,
                        require_summary=require_corrected_summary,
                        historical_compatibility=historical_compatibility,
                        label=description,
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

def expected_artifacts(
    root: Path,
    ds: str,
    label: str,
    batch: bool,
    batch_pass: str | None,
    analysis_variant: str | None = None,
) -> list[Path]:
    stem = batch_artifact_stem(
        ds, batch_pass if batch else None, analysis_variant if batch else None
    )
    if label == "prepare_pseudobulk":
        if batch:
            return [root / "pseudobulks" / f"{stem}_pseudobulk_hvg2000.rds"]
        return [
            root / "pseudobulks" / f"{ds}_pseudobulk_{variant}.rds"
            for variant in PB_VARIANTS
        ]
    if label in {"trans", "zeroimp"}:
        return [root / "results" / f"{ds}_{label}.rds"]
    if label in R_METHODS:
        return [root / "results" / f"{stem}_{label}.rds"]
    if label == "mrvi":
        if batch:
            return [
                root / "embeddings" / f"{stem}_hvg2000_highres_mrvi_dists.feather"
            ]
        return [
            root / "embeddings" / f"{ds}_hvg{n}_mrvi_dists.feather"
            for n in (1000, 2000, 3000)
        ]
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
                / f"{stem}_hvg2000_highres_{label}_dists.feather"
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
    *,
    analysis_variant: str | None = None,
    expected_batch_contract=None,
    batch_contract=None,
    historical_compatibility: bool = False,
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
    if analysis_variant not in (None, "", "final", "corrected_final"):
        raise ValueError(f"unknown analysis variant: {analysis_variant}")
    if historical_compatibility and analysis_variant:
        raise ValueError(
            "historical compatibility is read-only and cannot validate a "
            "variant-qualified corrected output"
        )
    variant_datasets: tuple[str, ...] | None = None
    if analysis_variant in {"final", "corrected_final"}:
        expected_pass = (
            "uncorrected" if analysis_variant == "final" else "corrected"
        )
        if not batch or batch_pass != expected_pass:
            raise ValueError(
                f"{analysis_variant} analysis variant requires the "
                f"{expected_pass} batch-effect pass"
            )
        forbidden = sorted(set(allowed) - set(FINAL_BATCH_METHODS))
        if forbidden:
            raise ValueError(
                f"{analysis_variant} analysis variant has unsupported methods: "
                + ", ".join(forbidden)
            )
        _validate_variant_root(root, analysis_variant)
        if analysis_variant == "final":
            variant_datasets = FINAL_UNCORRECTED_DATASETS
        else:
            if config_path is None or not config_path.is_file():
                raise ValueError(
                    "corrected_final artifact validation requires an existing --config"
                )
            variant_datasets = _configured_batch_effect_datasets(
                config_path, "batch_effect_corrected"
            )
            if not variant_datasets:
                raise ValueError(
                    "corrected_final config has no selected corrected datasets"
                )
    if batch and batch_pass not in {"uncorrected", "corrected"}:
        raise ValueError("batch validation requires --batch-pass")
    corrected = batch and batch_pass == "corrected"
    config_by_view: dict[str, dict] = {}
    if corrected:
        if config_path is None or not config_path.is_file():
            raise ValueError(
                "corrected batch artifact validation requires an existing --config"
            )
    if (
        variant_datasets is not None
        and len(rows) > 1
        and tuple(row[0] for row in rows) != variant_datasets
    ):
        raise ValueError(
            f"{analysis_variant} selection must use its exact configured dataset order"
        )
    if batch and exact:
        expected_rows = [
            (ds, "batch_effect_uncorrected", "batch_effect_uncorrected")
            for ds in BATCH_DATASET_ORDER
        ]
        if rows != expected_rows or batch_pass != "uncorrected":
            raise ValueError("batch exact selection is not the literal twelve-row uncorrected matrix")
    for ds, view, scope in rows:
        if variant_datasets is not None and ds not in variant_datasets:
            raise ValueError(
                f"{analysis_variant} selection contains an unapproved dataset: {ds}"
            )
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
        if corrected and view not in config_by_view:
            config_by_view[view] = read_datasets_json(config_path, view=view)
        for label in selected_labels:
            row_expected_batch_contract = expected_batch_contract
            if corrected:
                derived_batch_contract = _build_corrected_batch_contract(
                    config_by_view[view], ds, view, label
                )
                row_expected_batch_contract = derived_batch_contract
                if expected_batch_contract is not None:
                    validate_batch_contract_identity(
                        derived_batch_contract,
                        expected_batch_contract,
                        require_recorded=True,
                        require_summary=False,
                        historical_compatibility=historical_compatibility,
                        label=f"{ds}/{view}/{label} supplied identity",
                    )
                    row_expected_batch_contract = expected_batch_contract
            paths = expected_artifacts(
                root, ds, label, batch, batch_pass, analysis_variant
            )
            selected_paths.extend(paths)
            require_nonempty(
                paths,
                f"{ds}/{view}/{label}",
                expected,
                producer=producer or label,
                producer_run_id=producer_run_id,
                expected_batch_contract=row_expected_batch_contract,
                batch_contract=batch_contract,
                require_runtime_batch_contract=corrected,
                require_corrected_summary=analysis_variant != "corrected_final",
                historical_compatibility=historical_compatibility,
            )
    _reject_selected_partials(selected_paths, producer_run_id)



def validate_single(
    path: Path,
    producer: str | None = None,
    producer_run_id: str | None = None,
    *,
    corrected: bool = False,
    analysis_variant: str | None = None,
    expected_batch_contract=None,
    batch_contract=None,
    historical_compatibility: bool = False,
) -> None:
    if analysis_variant not in (None, "", "final", "corrected_final"):
        raise ValueError(f"unknown analysis variant: {analysis_variant}")
    if historical_compatibility and analysis_variant:
        raise ValueError(
            "historical compatibility is read-only and cannot validate a "
            "variant-qualified corrected output"
        )
    if analysis_variant == "final" and corrected:
        raise ValueError("final analysis variant cannot validate corrected Stage 5 artifacts")
    if analysis_variant == "corrected_final" and not corrected:
        raise ValueError(
            "corrected_final analysis variant requires corrected Stage 5 artifacts"
        )
    if analysis_variant:
        # Every selected artifact must stay under the same bound variant root.
        # ``expected_artifacts`` places outputs one directory below that root
        # (results/, embeddings/, or pseudobulks/).
        _validate_variant_root(Path(path).parent.parent, analysis_variant)
        marker = (
            "_batch_effect_uncorrected_final_"
            if analysis_variant == "final"
            else "_batch_effect_corrected_final_"
        )
        if marker not in path.name:
            raise ValueError(
                f"{analysis_variant} artifact path is not variant-qualified: {path}"
            )
    if (
        corrected
        and expected_batch_contract is None
        and batch_contract is None
    ):
        raise ValueError(
            "corrected batch artifact validation requires explicit contract identity"
        )
    require_nonempty(
        [path],
        "benchmark artifact",
        producer=producer,
        producer_run_id=producer_run_id,
        expected_batch_contract=expected_batch_contract,
        batch_contract=batch_contract,
        require_runtime_batch_contract=corrected,
        require_corrected_summary=analysis_variant != "corrected_final",
        historical_compatibility=historical_compatibility,
    )
    _reject_selected_partials([Path(path)], producer_run_id)
def _load_batch_contract_argument(value):
    if value is None:
        return None
    candidate = Path(value)
    try:
        text = (
            candidate.read_text(encoding="utf-8")
            if candidate.is_file()
            else value
        )
        identity = json.loads(text)
    except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as exc:
        raise ValueError(
            "batch contract identity must be a JSON object or a readable JSON path"
        ) from exc
    if not isinstance(identity, dict):
        raise ValueError("batch contract identity JSON must be an object")
    return identity
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--artifact", type=Path)
    group.add_argument("--root", type=Path)
    parser.add_argument("--analysis-variant", default=None,
                        choices=["final", "corrected_final"],
                        help="variant-qualified batch artifact paths")
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--batch-pass", default=None)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--source-identity", type=Path, default=None)
    parser.add_argument("--source-identity-verified", action="store_true")
    parser.add_argument("--expected-batch-contract", default=None)
    parser.add_argument("--batch-contract", default=None)
    parser.add_argument("--producer", default=None)
    parser.add_argument("--producer-run-id", default=None)
    parser.add_argument(
        "--historical-compatibility",
        action="store_true",
        help="read-only compatibility for immutable historical corrected artifacts",
    )
    args = parser.parse_args()
    expected_batch_contract = _load_batch_contract_argument(
        args.expected_batch_contract
    )
    batch_contract = _load_batch_contract_argument(args.batch_contract)
    if args.artifact is not None:
        validate_single(
            args.artifact,
            producer=args.producer,
            producer_run_id=args.producer_run_id,
            corrected=args.batch_pass == "corrected",
            historical_compatibility=args.historical_compatibility,
            analysis_variant=args.analysis_variant,
            expected_batch_contract=expected_batch_contract,
            batch_contract=batch_contract,
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
            analysis_variant=args.analysis_variant,
            expected_batch_contract=expected_batch_contract,
            batch_contract=batch_contract,
            historical_compatibility=args.historical_compatibility,
        )
    print("matrix artifact contract OK")




if __name__ == "__main__":
    main()
