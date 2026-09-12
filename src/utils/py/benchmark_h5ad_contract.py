"""Content contract for benchmark-analysis AnnData artifacts."""

import os
import re
import numpy as np

REQUIRED_OBSM = {
    "benchmark_analysis": {
        "X_pca_benchmark_analysis_hvg1000",
        "X_pca_benchmark_analysis_hvg2000",
        "X_pca_benchmark_analysis_hvg3000",
        "X_pca_harmony_benchmark_analysis_hvg2000",
    },
    "batch_effect_uncorrected": {
        "X_pca_batch_effect_uncorrected_hvg2000",
    },
    "batch_effect_corrected": {
        "X_pca_batch_effect_corrected_hvg2000",
        "X_pca_harmony_batch_effect_corrected_hvg2000",
    },
}
from collections.abc import Mapping
import hashlib
import json
from pathlib import Path

try:
    from src.utils.py.batch_contract import (  # noqa: E402
        batch_correction_spec_for_keys as _batch_correction_spec_for_keys,
        validate_batch_validation_summary as _validate_batch_validation_summary,
    )
except ModuleNotFoundError:  # pragma: no cover - direct script execution
    from batch_contract import (  # type: ignore[no-redef]
        batch_correction_spec_for_keys as _batch_correction_spec_for_keys,
        validate_batch_validation_summary as _validate_batch_validation_summary,
    )


# These values are duplicated here deliberately rather than inferred from an
# audit/checksum sidecar.  The identity is an explicit corrected-mode caller
# contract; ordinary and uncorrected validation never enters this code.
_BATCH_TOKEN_VERSION = "ecoda_batch_composite_v1"
_BATCH_CONTRACT_VERSION = "ecoda_batch_contract_v1"
_BATCH_RESERVED_OBS_NAME = "__ecoda_batch_combined_v1"
_BATCH_DIRECT_SCALARIZATION = "direct_v1"
_BATCH_COMPOSITE_SCALARIZATION = "composite_v1"
_BATCH_H5AD_METHOD_ID = "preprocess"
_BATCH_H5AD_MODEL_ID = "hvg_composite_v1"
_BATCH_METHOD_IDS = frozenset(
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
_BATCH_MODEL_IDS = frozenset(
    {
        "hvg_composite_v1",
        "harmony_native_list_v1",
        "ecoda_additive_random_intercepts_v1",
        "pseudobulk_composite_v1",
        "mrvi_composite_v1",
        "embedding_consumer_harmony_v1",
    }
)
_BATCH_KEY_FIELDS = (
    "ordered_source_keys",
    "ordered_keys",
    "source_keys",
    "batch_keys",
    "keys",
)
_BATCH_TOKEN_FIELDS = ("token_version", "encoding_version", "encoding")
_BATCH_METHOD_FIELDS = ("method_id", "method", "method_policy")
_BATCH_MODEL_FIELDS = ("model_id", "model", "model_policy")
_BATCH_FINGERPRINT_FIELDS = (
    "fingerprint",
    "batch_contract_fingerprint",
    "key_set_fingerprint",
)
_BATCH_OBS_FIELDS = (
    "required_source_obs_columns",
    "required_obs_columns",
    "source_obs_columns",
    "obs_columns",
)
_BATCH_IDENTITY_VECTOR_FIELDS = frozenset(
    {
        "composite_values",
        "sample_composite_values",
        "sample_ids",
        "sample_group_ids",
        "canonical_values",
        "canonical_cell_values",
        "canonical_sample_metadata",
        "sample_metadata",
        "row_tokens",
        "tokens",
        "scalarized_values",
        "sample_scalarized_values",
    }
)
_BATCH_RESERVED_ABSENCE_FIELDS = (
    "reserved_obs_absent",
    "reserved_column_absent",
    "reserved_absent",
)


def _identity_alias(identity, fields, label):
    """Return one consistent value from equivalent identity field aliases."""
    present = [field for field in fields if field in identity]
    if not present:
        raise ValueError(
            f"corrected batch contract identity is missing {label}"
        )
    value = identity[present[0]]
    for field in present[1:]:
        if identity[field] != value:
            raise ValueError(
                f"corrected batch contract identity has mismatched {label} aliases"
            )
    return value


def _identity_string_list(value, label, *, nonempty=True):
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise ValueError(
            f"corrected batch contract identity {label} must be an ordered list of strings"
        )
    values = list(value)
    if nonempty and not values:
        raise ValueError(
            f"corrected batch contract identity {label} must be nonempty"
        )
    if any(
        not isinstance(item, str)
        or not item
        or item != item.strip()
        for item in values
    ):
        raise ValueError(
            f"corrected batch contract identity {label} contains a blank or "
            "whitespace-padded value"
        )
    if len(values) != len(set(values)):
        raise ValueError(
            f"corrected batch contract identity {label} contains duplicates"
        )
    return tuple(values)

def validate_batch_validation_summary(summary, batch_keys=None):
    """Validate one persisted vector-free corrected batch summary."""
    try:
        return _validate_batch_validation_summary(summary, batch_keys)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"validation_summary is invalid: {exc}") from exc

def _validate_identity_summary(
    identity,
    keys,
    label,
    *,
    method_id=None,
    required=False,
    validate_optional=True,
):
    """Validate an identity summary and its method-specific correction policy."""
    if not isinstance(identity, Mapping):
        raise ValueError(
            f"corrected batch contract identity for {label} must be a mapping"
        )
    has_summary = "validation_summary" in identity
    vector_fields = sorted(_BATCH_IDENTITY_VECTOR_FIELDS.intersection(identity))
    if (required or (validate_optional and has_summary)) and vector_fields:
        raise ValueError(
            f"{label} contains forbidden per-cell/sample vector fields: "
            f"{', '.join(vector_fields)}"
        )
    if not has_summary:
        if required:
            raise ValueError(f"{label} is missing validation_summary")
        return None
    if not required and not validate_optional:
        return None
    try:
        normalized = validate_batch_validation_summary(
            identity["validation_summary"],
            keys,
        )
        if method_id is not None:
            expected_mode, expected_formula = _batch_correction_spec_for_keys(
                method_id,
                keys,
            )
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} has an invalid validation_summary") from exc
    if method_id is not None and (
        normalized["correction_mode"] != expected_mode
        or normalized["correction_formula"] != expected_formula
    ):
        raise ValueError(
            f"{label} has the wrong correction policy for method {method_id!r}"
        )
    return normalized


def _identity_fingerprint(keys, scalarization, method_id, model_id):
    """Compute the byte-exact fingerprint without serializing an object."""
    def field(name, value):
        name_bytes = name.encode("utf-8")
        value_bytes = value.encode("utf-8")
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

    key_parts = [f"{len(keys)}|".encode("ascii")]
    for key in keys:
        key_bytes = key.encode("utf-8")
        key_parts.extend(
            (
                str(len(key_bytes)).encode("ascii"),
                b":",
                key_bytes.hex().encode("ascii"),
                b";",
            )
        )
    payload = bytearray(_BATCH_CONTRACT_VERSION.encode("utf-8") + b"\0")
    payload.extend(field("encoding", _BATCH_TOKEN_VERSION))
    payload.extend(field("keys", b"".join(key_parts).decode("ascii")))
    payload.extend(field("scalarization", scalarization))
    payload.extend(field("method", method_id))
    payload.extend(field("model", model_id))
    return hashlib.sha256(bytes(payload)).hexdigest()


def _normalize_batch_contract_identity(identity, label="batch contract"):
    """Normalize and validate one explicit corrected identity mapping."""
    if not isinstance(identity, Mapping):
        raise ValueError(
            f"corrected batch contract identity for {label} must be a mapping"
        )

    keys = _identity_string_list(
        _identity_alias(identity, _BATCH_KEY_FIELDS, "ordered source keys"),
        "ordered source keys",
    )
    if _BATCH_RESERVED_OBS_NAME in keys:
        raise ValueError(
            "corrected batch contract identity ordered source keys contain "
            f"reserved column {_BATCH_RESERVED_OBS_NAME!r}"
        )

    token_values = _identity_alias(
        identity, _BATCH_TOKEN_FIELDS, "token/encoding version"
    )
    if (
        not isinstance(token_values, str)
        or token_values != _BATCH_TOKEN_VERSION
    ):
        raise ValueError(
            "corrected batch contract identity token/encoding version must be "
            f"{_BATCH_TOKEN_VERSION!r}"
        )

    scalarization = identity.get("scalarization")
    if not isinstance(scalarization, str):
        raise ValueError(
            "corrected batch contract identity is missing scalarization"
        )
    expected_scalarization = (
        _BATCH_COMPOSITE_SCALARIZATION
        if len(keys) >= 2
        else _BATCH_DIRECT_SCALARIZATION
    )
    if scalarization != expected_scalarization:
        raise ValueError(
            "corrected batch contract identity scalarization does not match "
            "the ordered source key count"
        )

    method_id = _identity_alias(identity, _BATCH_METHOD_FIELDS, "method policy")
    model_id = _identity_alias(identity, _BATCH_MODEL_FIELDS, "model policy")
    if not isinstance(method_id, str) or method_id not in _BATCH_METHOD_IDS:
        raise ValueError(
            f"corrected batch contract identity has unsupported method policy: {method_id!r}"
        )
    if not isinstance(model_id, str) or model_id not in _BATCH_MODEL_IDS:
        raise ValueError(
            f"corrected batch contract identity has unsupported model policy: {model_id!r}"
        )

    contract_version = identity.get("contract_version", _BATCH_CONTRACT_VERSION)
    if contract_version != _BATCH_CONTRACT_VERSION:
        raise ValueError(
            "corrected batch contract identity has an unsupported contract version"
        )

    required_obs_columns = _identity_string_list(
        _identity_alias(
            identity, _BATCH_OBS_FIELDS, "required source obs columns"
        ),
        "required source obs columns",
    )
    if _BATCH_RESERVED_OBS_NAME in required_obs_columns:
        raise ValueError(
            "corrected batch contract identity required source obs columns "
            f"contain reserved column {_BATCH_RESERVED_OBS_NAME!r}"
        )
    missing_keys = [key for key in keys if key not in required_obs_columns]
    if missing_keys:
        raise ValueError(
            "corrected batch contract identity required source obs columns "
            f"omit ordered source keys: {', '.join(missing_keys)}"
        )

    reserved_absent = _identity_alias(
        identity,
        _BATCH_RESERVED_ABSENCE_FIELDS,
        "reserved-column absence",
    )
    if type(reserved_absent) is not bool or not reserved_absent:
        raise ValueError(
            "corrected batch contract identity must assert absence of "
            f"{_BATCH_RESERVED_OBS_NAME!r}"
        )
    reserved_name_fields = ("reserved_obs_name", "reserved_name")
    present_reserved_name_fields = [
        field for field in reserved_name_fields if field in identity
    ]
    if present_reserved_name_fields:
        reserved_name = _identity_alias(
            identity, reserved_name_fields, "reserved-column name"
        )
        if reserved_name != _BATCH_RESERVED_OBS_NAME:
            raise ValueError(
                "corrected batch contract identity has an unsupported "
                "reserved-column name"
            )

    fingerprint = _identity_alias(identity, _BATCH_FINGERPRINT_FIELDS, "fingerprint")
    expected_fingerprint = _identity_fingerprint(
        keys, scalarization, method_id, model_id
    )
    if (
        not isinstance(fingerprint, str)
        or not re.fullmatch(r"[0-9a-f]{64}", fingerprint)
        or fingerprint != expected_fingerprint
    ):
        raise ValueError(
            "corrected batch contract identity fingerprint is missing, malformed, "
            "or does not match its source/configuration fields"
        )
    normalized = {
        "ordered_source_keys": keys,
        "token_version": token_values,
        "scalarization": scalarization,
        "method_id": method_id,
        "model_id": model_id,
        "fingerprint": fingerprint,
        "required_source_obs_columns": required_obs_columns,
        "reserved_obs_absent": True,
        "contract_version": contract_version,
    }
    if "validation_summary" in identity:
        _validate_identity_summary(
            identity,
            keys,
            label,
            method_id=method_id,
        )
    return normalized


def _normalize_h5ad_expected_identity(identity, label):
    expected = validate_batch_contract_identity(
        identity,
        require_recorded=False,
        label=label,
    )
    if (
        expected["method_id"] != _BATCH_H5AD_METHOD_ID
        or expected["model_id"] != _BATCH_H5AD_MODEL_ID
    ):
        raise ValueError(
            f"{label} must use the corrected H5AD preprocessing identity "
            f"{_BATCH_H5AD_METHOD_ID}/{_BATCH_H5AD_MODEL_ID}"
        )
    summary = _validate_identity_summary(
        identity,
        expected["ordered_source_keys"],
        f"{label} expected identity",
        method_id=expected["method_id"],
    )
    if summary is not None:
        expected["validation_summary"] = summary
    return expected

def validate_batch_contract_identity(
    expected_batch_contract=None,
    batch_contract=None,
    *,
    source_obs_columns=None,
    reserved_absent=None,
    require_recorded=False,
    label="batch contract",
    require_summary=None,
):
    """Validate corrected identity fields and any persisted validation summary.

    Configuration-only expected identities may omit a summary.  A summary is
    validated, including its method-specific correction policy, whenever it
    is present; explicit summaries are compared only when both identities
    provide one.
    """
    if expected_batch_contract is None and batch_contract is None:
        return None
    expected_source = (
        expected_batch_contract
        if expected_batch_contract is not None
        else batch_contract
    )
    expected = _normalize_batch_contract_identity(
        expected_source, f"{label} (expected)"
    )
    summary_required = (
        require_recorded if require_summary is None else require_summary
    )
    expected_summary = _validate_identity_summary(
        expected_source,
        expected["ordered_source_keys"],
        f"{label} expected identity",
        method_id=expected["method_id"],
        validate_optional=True,
    )
    if batch_contract is None:
        if require_recorded:
            raise ValueError(
                f"{label} is missing recorded corrected batch contract identity"
            )
        recorded = None
        recorded_summary = None
    else:
        recorded = _normalize_batch_contract_identity(
            batch_contract, f"{label} (recorded)"
        )
        recorded_summary = _validate_identity_summary(
            batch_contract,
            recorded["ordered_source_keys"],
            f"{label} recorded identity",
            method_id=recorded["method_id"],
            required=bool(summary_required),
            validate_optional=True,
        )
        if recorded != expected:
            raise ValueError(
                f"{label} source/config identity does not match the expected "
                "corrected batch contract"
            )
        if (
            expected_summary is not None
            and recorded_summary is not None
            and recorded_summary != expected_summary
        ):
            raise ValueError(
                f"{label} validation_summary does not match the expected "
                "corrected metadata"
            )

    if source_obs_columns is not None:
        try:
            observed_columns = list(source_obs_columns)
        except TypeError:
            raise ValueError(
                f"{label} source obs columns are not iterable"
            ) from None
        if any(not isinstance(column, str) for column in observed_columns):
            raise ValueError(f"{label} source obs columns are not strings")
        missing = [
            column
            for column in expected["required_source_obs_columns"]
            if column not in observed_columns
        ]
        if missing:
            raise ValueError(
                f"{label} source obs columns are missing: {', '.join(missing)}"
            )
        if _BATCH_RESERVED_OBS_NAME in observed_columns:
            raise ValueError(
                f"{label} source obs contains reserved temporary column "
                f"{_BATCH_RESERVED_OBS_NAME!r}"
            )
    if reserved_absent is not None and (
        type(reserved_absent) is not bool or not reserved_absent
    ):
        raise ValueError(
            f"{label} source obs contains reserved temporary column "
            f"{_BATCH_RESERVED_OBS_NAME!r}"
        )
    return expected if recorded is None else recorded
def _embedded_batch_contract(mapping):
    """Return an explicitly recorded contract from an AnnData ``uns`` map."""
    if not isinstance(mapping, Mapping):
        return None
    for field in (
        "batch_contract",
        "batch_contract_identity",
        "ecoda_batch_contract",
        "_ecoda_batch_contract",
    ):
        value = mapping.get(field)
        if value is not None:
            return value
    return None


def _reject_duplicate_json_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"embedded identity has duplicate key: {key}")
        result[key] = value
    return result


def _h5_identity_value(node):
    """Decode the small JSON-like values used for ``uns`` identity metadata."""
    if hasattr(node, "keys"):
        return {
            str(name): _h5_identity_value(node[name])
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
            decoded = json.loads(
                value,
                object_pairs_hook=_reject_duplicate_json_keys,
            )
        except json.JSONDecodeError:
            pass
        else:
            return decoded
    return value


def _embedded_h5_batch_contract(handle):
    """Read a recorded contract from ``uns`` without treating checksums as proof."""
    uns = handle.get("uns")
    if uns is None:
        return None
    for field in (
        "batch_contract",
        "batch_contract_identity",
        "ecoda_batch_contract",
        "_ecoda_batch_contract",
    ):
        if field in uns:
            return _h5_identity_value(uns[field])
    return None


def _h5_shape(node):
    """Return a node's persisted shape, including group-backed values."""
    if node is None:
        return ()
    shape = getattr(node, "shape", None)
    if shape is None:
        shape = node.attrs.get("shape")
    if shape is None:
        return ()
    try:
        return tuple(int(value) for value in shape)
    except (TypeError, ValueError):
        return ()


def _h5_nonempty_node(node):
    shape = _h5_shape(node)
    return bool(shape) and all(value > 0 for value in shape)


def _h5_vector_non_na(node):
    """Return the number of non-NA values in a scalar/vector HDF5 field."""
    value_node = node
    if hasattr(node, "keys"):
        for name in ("values", "data", "_values", "codes"):
            if name in node:
                value_node = node[name]
                break
    try:
        values = np.asarray(value_node[()]).reshape(-1)
    except (TypeError, ValueError, OSError):
        return 0
    if "mask" in getattr(node, "keys", lambda: ())():
        try:
            mask = np.asarray(node["mask"][()]).reshape(-1).astype(bool)
            return int((~mask).sum())
        except (TypeError, ValueError, OSError):
            return 0
    if getattr(value_node, "name", "").endswith("/codes"):
        try:
            return int((values >= 0).sum())
        except TypeError:
            return 0
    if np.issubdtype(values.dtype, np.number):
        return int((~np.isnan(values)).sum())
    return int(values.size)


def _h5_index_name(group):
    index_name = group.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    return str(index_name)


def _validate_count_values(values, label):
    values = np.asarray(values)
    if values.size == 0:
        raise ValueError(f"{label} is empty")
    try:
        finite = np.isfinite(values).all()
        nonnegative = (values >= 0).all()
        integral = np.equal(values, np.floor(values)).all()
    except (TypeError, ValueError):
        raise ValueError(f"{label} is not numeric") from None
    if not finite:
        raise ValueError(f"{label} contains nonfinite values")
    if not nonnegative or not integral:
        raise ValueError(f"{label} must be finite, nonnegative, integer-valued counts")


def _validate_counts_layer(layer, shape, label="layers['counts']"):
    import scipy.sparse as sparse

    layer_shape = _h5_shape(layer)
    if layer_shape and tuple(shape) != layer_shape:
        raise ValueError(f"{label} shape {layer_shape} does not match {tuple(shape)}")
    if sparse.issparse(layer):
        values = layer.data
    else:
        values = np.asarray(layer)
    _validate_count_values(values, label)
def _validate_h5_count_node(node, label, max_values=1024 * 1024):
    shape = _h5_shape(node)
    if not shape:
        raise ValueError(f"{label} has no persisted shape")
    if hasattr(node, "shape"):
        if len(shape) == 1:
            for start in range(0, shape[0], max_values):
                _validate_count_values(
                    node[start:start + max_values], label
                )
        else:
            rows_per = max(1, max_values // max(1, int(np.prod(shape[1:])))
                           )
            for start in range(0, shape[0], rows_per):
                _validate_count_values(node[start:start + rows_per], label)
    else:
        raise ValueError(f"{label} is not a readable HDF5 dataset")


def _validate_h5_counts_layer(layer, expected_shape, label="layers['counts']"):
    shape = _h5_shape(layer)
    if tuple(shape) != tuple(expected_shape):
        raise ValueError(
            f"{label} shape {shape} does not match {tuple(expected_shape)}"
        )
    if hasattr(layer, "keys"):
        encoding = layer.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise ValueError(f"{label} has unsupported encoding {encoding!r}")
        if "data" not in layer:
            raise ValueError(f"{label} has no data values")
        _validate_h5_count_node(layer["data"], label)
    else:
        _validate_h5_count_node(layer, label)


def validate_benchmark_h5ad_contract(
    adata,
    view,
    method,
    *,
    expected_batch_contract=None,
    batch_contract=None,
):
    """Reject incomplete AnnData artifacts before benchmark computation."""
    if view not in REQUIRED_OBSM:
        raise ValueError(f"Unknown preprocessing view for h5ad contract: {view}")
    if (
        view == "batch_effect_corrected"
        and expected_batch_contract is None
    ):
        raise ValueError(
            "corrected batch h5ad validation requires explicit contract identity "
            "(expected config identity)"
        )
    if view == "batch_effect_corrected":
        expected_batch_contract = _normalize_h5ad_expected_identity(
            expected_batch_contract,
            f"h5ad {method} ({view}) expected identity",
        )
    if expected_batch_contract is not None or batch_contract is not None:
        obs_columns = list(getattr(adata.obs, "columns", ()))
        embedded_batch_contract = _embedded_batch_contract(
            getattr(adata, "uns", {})
        )
        if batch_contract is not None and embedded_batch_contract is not None:
            validate_batch_contract_identity(
                batch_contract,
                embedded_batch_contract,
                require_recorded=True,
                require_summary=view == "batch_effect_corrected",
                label=f"h5ad {method} ({view}) embedded identity",
            )
        recorded_batch_contract = (
            batch_contract
            if batch_contract is not None
            else embedded_batch_contract
        )
        validate_batch_contract_identity(
            expected_batch_contract,
            recorded_batch_contract,
            source_obs_columns=obs_columns,
            reserved_absent=(
                _BATCH_RESERVED_OBS_NAME not in obs_columns
            ),
            require_recorded=expected_batch_contract is not None,
            require_summary=view == "batch_effect_corrected",
            label=f"h5ad {method} ({view})",
        )
    missing = []
    try:
        x_present = adata.X is not None
    except Exception:
        x_present = False
    if not x_present:
        missing.append("X")
    if "counts" not in adata.layers:
        missing.append("layers['counts']")
    if "counts" in adata.layers:
        _validate_counts_layer(
            adata.layers["counts"], (adata.n_obs, adata.n_vars)
        )

    obsm = getattr(adata, "obsm", {})
    missing_obsm = sorted(REQUIRED_OBSM[view] - set(obsm.keys()))
    missing.extend(f"obsm['{key}']" for key in missing_obsm)

    if "hvg_rank" not in adata.var.columns:
        missing.append("var['hvg_rank']")
    else:
        required_hvg = 3000 if view == "benchmark_analysis" else 2000
        n_ranked = int(adata.var["hvg_rank"].notna().sum())
        if n_ranked < required_hvg:
            missing.append(
                f"var['hvg_rank'] with at least {required_hvg} non-NA ranks "
                f"(found {n_ranked})"
            )

    if "Sample" not in adata.obs.columns:
        missing.append("obs['Sample']")
    if getattr(adata, "n_obs", 0) <= 0 or getattr(adata, "n_vars", 0) <= 0:
        missing.append("non-empty obs/var")

    if missing:
        raise ValueError(
            f"h5ad content contract failed for {method} ({view}): "
            f"missing or invalid {', '.join(missing)}. "
            "Re-run 1.1.1_preprocess.py with --force and use the "
            "authoritative processed h5ad."
        )


def validate_benchmark_h5ad_path(
    path,
    view,
    method,
    *,
    expected_batch_contract=None,
    batch_contract=None,
):
    """Validate persisted h5ad structure without materializing AnnData."""
    import h5py

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError(f"h5ad path is missing or empty: {path}")
    if view not in REQUIRED_OBSM:
        raise ValueError(f"Unknown preprocessing view for h5ad contract: {view}")
    if (
        view == "batch_effect_corrected"
        and expected_batch_contract is None
    ):
        raise ValueError(
            "corrected batch h5ad validation requires explicit contract identity "
            "(expected config identity)"
        )
    if view == "batch_effect_corrected":
        expected_batch_contract = _normalize_h5ad_expected_identity(
            expected_batch_contract,
            f"h5ad {method} ({view}) expected identity",
        )

    identity_active = (
        expected_batch_contract is not None or batch_contract is not None
    )
    missing = []
    with h5py.File(path, "r") as handle:
        obs_columns = (
            list(handle["obs"].keys()) if "obs" in handle else None
        )
        if identity_active:
            embedded_batch_contract = _embedded_h5_batch_contract(handle)
            if (
                batch_contract is not None
                and embedded_batch_contract is not None
            ):
                validate_batch_contract_identity(
                    batch_contract,
                    embedded_batch_contract,
                    require_recorded=True,
                    require_summary=view == "batch_effect_corrected",
                    label=f"h5ad {method} ({view}) embedded identity",
                )
            recorded_batch_contract = (
                batch_contract
                if batch_contract is not None
                else embedded_batch_contract
            )
            validate_batch_contract_identity(
                expected_batch_contract,
                recorded_batch_contract,
                source_obs_columns=obs_columns,
                reserved_absent=(
                    obs_columns is not None
                    and _BATCH_RESERVED_OBS_NAME not in obs_columns
                ),
                require_recorded=expected_batch_contract is not None,
                require_summary=view == "batch_effect_corrected",
                label=f"h5ad {method} ({view})",
            )

        x_shape = _h5_shape(handle["X"]) if "X" in handle else ()
        if "X" not in handle:
            missing.append("X")
        elif len(x_shape) != 2 or any(value <= 0 for value in x_shape):
            missing.append("X with a non-empty persisted shape")
        if "layers" not in handle or "counts" not in handle["layers"]:
            missing.append("layers['counts']")
        else:
            counts = handle["layers"]["counts"]
            counts_shape = _h5_shape(counts)
            expected_shape = x_shape if len(x_shape) == 2 else counts_shape
            try:
                _validate_h5_counts_layer(counts, expected_shape)
            except (OSError, TypeError, ValueError) as exc:
                missing.append(str(exc))
        if "obsm" not in handle:
            missing.extend(f"obsm['{key}']" for key in sorted(REQUIRED_OBSM[view]))
        else:
            for key in sorted(REQUIRED_OBSM[view]):
                if key not in handle["obsm"]:
                    missing.append(f"obsm['{key}']")

        required_hvg = 3000 if view == "benchmark_analysis" else 2000
        if "var" not in handle or "hvg_rank" not in handle["var"]:
            missing.append("var['hvg_rank']")
        else:
            hvg_rank = handle["var"]["hvg_rank"]
            if _h5_vector_non_na(hvg_rank) < required_hvg:
                missing.append(
                    f"var['hvg_rank'] with at least {required_hvg} non-NA ranks"
                )

        if "obs" not in handle or "Sample" not in handle["obs"]:
            missing.append("obs['Sample']")
        if "obs" not in handle:
            missing.append("non-empty obs index")
        else:
            obs = handle["obs"]
            index_name = _h5_index_name(obs)
            if index_name not in obs or not _h5_nonempty_node(obs[index_name]):
                missing.append("non-empty obs index")
        if "var" not in handle:
            missing.append("non-empty var index")
        else:
            var = handle["var"]
            index_name = _h5_index_name(var)
            if index_name not in var or not _h5_nonempty_node(var[index_name]):
                missing.append("non-empty var index")

    if missing:
        raise ValueError(
            f"h5ad content contract failed for {method} ({view}): "
            f"missing or invalid {', '.join(dict.fromkeys(missing))}. "
            "Re-run 1.1.1_preprocess.py with --force and use the "
            "authoritative processed h5ad."
        )



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
    if not isinstance(identity, Mapping):
        raise ValueError("batch contract identity JSON must be an object")
    return identity
def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", required=True)
    parser.add_argument("--view", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--expected-batch-contract", default=None)
    parser.add_argument("--batch-contract", default=None)
    args = parser.parse_args()
    validate_benchmark_h5ad_path(
        args.path,
        args.view,
        args.method,
        expected_batch_contract=_load_batch_contract_argument(
            args.expected_batch_contract
        ),
        batch_contract=_load_batch_contract_argument(args.batch_contract),
    )
    print(f"h5ad contract OK: {args.path}")


if __name__ == "__main__":
    main()
