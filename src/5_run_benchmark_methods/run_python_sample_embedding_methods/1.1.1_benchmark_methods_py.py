"""Python benchmark methods (MrVI, scPoli, PILOT, QOT, PILOT-GM-VAE) as a CLI script.

Replaces the logic of the archived notebook
`1.2_benchmark_methods_py.qmd` (kept as reference; do NOT delete). Consumes
the preprocessed benchmark view h5ad produced by
`src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`:

- PILOT/QOT/PILOT-GM-VAE consume only selected obs columns and the stored
  obsm embedding `X_pca_{view}_hvg{n}`. They use the h5py/minimal-AnnData
  loader and never materialize `X` or `layers["counts"]`;
- MrVI/scPoli use the stored `var["hvg_rank"]` (computed by
  `select_hvgs_ranked`) to stream only the requested HVG columns from the raw
  CSR counts layer into a minimal AnnData; MrVI keeps those counts sparse and
  scPoli densifies only its selected HVG subset to float32;
- cell type annotation columns come from datasets.json
  (`cell_type_low_res` / `cell_type_high_res`).

Feather naming, method-string format and data layout are preserved exactly
from the qmd (the R ingest functions `process_mrvi_fig` /
`process_scpoli_fig` / `process_pilot_fig` / `process_qot_fig` /
`process_pilotgm_fig`, `constants.R` label map and the notebook recodes
depend on them): plain `DataFrame.to_feather()` with the pandas index (sample
names) kept — the index is written as the last feather column, matching R's
`column_to_rownames(ncol)`.

QOT and PILOT-GM-VAE are extended methods (see the implementation plan
`.kilo/plans/1786651957910-pilotgm-qot-benchmark-implementation.md`): QOT runs
the vendored `qot_utils_re.py` (PennShenLab/QOT @ 28cd529880c1, one bug fix
in `Gaussian_Mixture_Representation`), PILOT-GM-VAE runs the `pilotgm` PyPI
package (CostaLab/PILOT-GM-VAE, BIB 2025). Both receive a distinct temp obs
column (`_bench_prog` / `_bench_status`) instead of `"Sample"` for the
status/progession argument: their `rename()` dicts collapse when the sample
and status column are the same (duplicate dict key -> both columns renamed
to 'status', no 'sampleID' survives -> KeyError in the GMM groupby). The
temp column also keeps any bio label out of the distance path (no-leakage).

Execution time (float seconds, method body only — excluding h5ad loading,
like the qmd and R `exec_time()`) and peak RSS (`mem_GB`) are appended to a
per-task log feather. One process per task writes the file, so no concurrency
issues. Combos run defaults-first (MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
PILOT_hvg2000_highres, QOT_hvg2000_highres, PILOT-GM-VAE_hvg2000_highres) so
the main-method rows are measured before any in-process memory bloat (peak
RSS is monotonic within a process).
"""

import argparse
import gc
import hashlib
import json
import os
import re
import resource
import sys
import tempfile
import time
from pathlib import Path

# This script sits one level deeper than 1.1.1_preprocess.py, so the repo
# root is parents[3] (preprocess uses parents[2] from src/3_scrnaseq_preprocessing/).
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
from scvi.external import MRVI
import pilotpy as pl
import torch

from src.utils.py.datasets_io import read_datasets_json
from src.utils.py.batch_contract import (
    augment_batch_contract,
    batch_correction_spec_for_keys,
    build_batch_composite,
    build_batch_contract_identity,
    build_batch_validation_summary,
    normalize_batch_keys,
    read_h5ad_validation_summary,
    validate_batch_metadata,
    validate_batch_validation_summary,
)
from src.utils.py.h5ad_counts_free import load_h5ad_counts_free
from src.utils.py.h5ad_counts_subset import (
    load_h5ad_counts_subset,
    read_h5ad_hvg_genes,
)
from src.utils.py.benchmark_h5ad_contract import (
    validate_batch_contract_identity,
    validate_benchmark_h5ad_contract,
    validate_benchmark_h5ad_path,
)
_CORRECTED_H5AD_METHOD_ID = "preprocess"
_CORRECTED_H5AD_MODEL_ID = "hvg_composite_v1"
_CORRECTED_STAGE5_METHOD_MODELS = {
    "mrvi": ("MrVI", "mrvi_composite_v1"),
    "pilot": ("PILOT", "embedding_consumer_harmony_v1"),
    "qot": ("QOT", "embedding_consumer_harmony_v1"),
}


def _file_md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

_CHECKSUM_FIELDS = ("MD5", "SIZE", "PATH")
_ARTIFACT_RECORD_FIELDS = ("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
_MD5_RE = re.compile(r"^[0-9a-f]{32}$")
_RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
def _read_corrected_validation_summary(path, batch_keys):
    """Read and validate the compact summary embedded in a corrected H5AD."""

    try:
        normalized = read_h5ad_validation_summary(str(path), batch_keys)
    except ValueError as exc:
        raise ValueError(
            f"corrected H5AD validation_summary is invalid: {path}"
        ) from exc
    expected_mode, expected_formula = batch_correction_spec_for_keys(
        "preprocess",
        batch_keys,
    )
    if (
        normalized["correction_mode"] != expected_mode
        or normalized["correction_formula"] != expected_formula
    ):
        raise ValueError(
            f"corrected H5AD validation_summary has the wrong correction policy: {path}"
        )
    return normalized


STAGE5_ARTIFACT_PRODUCERS = {
    "mrvi": "stage5_mrvi",
    "scpoli": "stage5_scpoli",
    "pilot": "stage5_pilot",
    "qot": "stage5_qot",
    "pilotgm": "stage5_pilotgm",
}


def stage5_artifact_producer(method):
    """Return the canonical run-bound producer for one Stage 5 method."""
    try:
        return STAGE5_ARTIFACT_PRODUCERS[method]
    except (KeyError, TypeError) as exc:
        raise ValueError(
            f"unsupported Stage 5 benchmark method: {method!r}"
        ) from exc


def _corrected_h5ad_identity(batch_keys):
    """Build the source identity required for corrected H5AD validation."""
    return build_batch_contract_identity(
        batch_keys,
        sample_column="Sample",
        method_id=_CORRECTED_H5AD_METHOD_ID,
        model_id=_CORRECTED_H5AD_MODEL_ID,
    )


def _corrected_stage5_identity(method, batch_keys, validation_summary=None):
    """Build Stage 5 identity and attach the source validation summary."""

    method_model = _CORRECTED_STAGE5_METHOD_MODELS.get(method)
    if method_model is None:
        return None
    method_id, model_id = method_model
    identity = build_batch_contract_identity(
        batch_keys,
        sample_column="Sample",
        method_id=method_id,
        model_id=model_id,
    )
    if validation_summary is None:
        return identity
    correction_mode, correction_formula = batch_correction_spec_for_keys(
        method_id,
        batch_keys,
    )
    return augment_batch_contract(
        identity,
        validation_summary,
        correction_mode,
        correction_formula,
    )


def _validate_h5ad_path(path, view, method, expected_batch_contract):
    """Validate an H5AD path, preserving legacy call shape when unbound."""
    if expected_batch_contract is None:
        return validate_benchmark_h5ad_path(path, view, method)
    return validate_benchmark_h5ad_path(
        path,
        view,
        method,
        expected_batch_contract=expected_batch_contract,
    )


def _validate_corrected_batch_contract(expected, recorded, label):
    """Validate identity and exact compact summary for corrected caches."""

    validate_batch_contract_identity(
        expected,
        recorded,
        require_recorded=True,
        label=label,
    )
    if not isinstance(expected, dict) or not isinstance(recorded, dict):
        raise ValueError(f"{label} must be a mapping")
    expected_summary = expected.get("validation_summary")
    recorded_summary = recorded.get("validation_summary")
    if expected_summary is None or recorded_summary is None:
        raise ValueError(f"{label} is missing validation_summary")
    forbidden = {
        "composite_values",
        "sample_composite_values",
        "sample_ids",
        "sample_group_ids",
        "tokens",
    }
    for mapping_name, mapping in (("expected", expected), ("recorded", recorded)):
        vectors = sorted(forbidden.intersection(mapping))
        if vectors:
            raise ValueError(
                f"{label} {mapping_name} contains per-cell/sample vectors: "
                f"{', '.join(vectors)}"
            )
    try:
        expected_normalized = validate_batch_validation_summary(
            expected_summary,
            expected["ordered_source_keys"],
        )
        recorded_normalized = validate_batch_validation_summary(
            recorded_summary,
            expected["ordered_source_keys"],
        )
    except ValueError as exc:
        raise ValueError(f"{label} has an invalid validation_summary") from exc
    if recorded_normalized != expected_normalized:
        raise ValueError(
            f"{label} validation_summary does not match expected corrected metadata"
        )


def _validate_h5ad_content(adata, view, method, expected_batch_contract):
    """Validate loaded H5AD content with corrected identity when required."""
    if expected_batch_contract is None:
        return validate_benchmark_h5ad_contract(adata, view, method)
    return validate_benchmark_h5ad_contract(
        adata,
        view,
        method,
        expected_batch_contract=expected_batch_contract,
    )

def _read_checksum_sidecar(path):
    """Read an exact MD5/SIZE/PATH sidecar without hashing artifact bytes."""
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or path.stat().st_size <= 0 or not sidecar.is_file():
        raise ValueError(f"checksum sidecar is missing: {sidecar}")
    try:
        lines = sidecar.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"checksum sidecar is unreadable: {sidecar}") from exc
    if len(lines) != len(_CHECKSUM_FIELDS):
        raise ValueError(f"checksum sidecar has an invalid schema: {sidecar}")
    records = {}
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


def _full_checksum(path):
    """Strictly verify an artifact sidecar and the bytes it describes."""
    path = Path(path)
    records = _read_checksum_sidecar(path)
    digest = _file_md5(path)
    if digest != records["MD5"]:
        raise ValueError(f"checksum sidecar does not match artifact: {path}")
    if path.stat().st_size != int(records["SIZE"]):
        raise ValueError(f"artifact changed during checksum validation: {path}")
    return records


def _record_context(run_id=None):
    """Resolve one explicit/current run without scanning other run roots."""
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


def artifact_record_path(path, run_id):
    """Return the bounded record path for a canonical artifact path."""
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


def _record_candidate(path, run_id=None):
    context = _record_context(run_id)
    if context is None:
        return None
    return artifact_record_path(path, context[1])


def _read_artifact_record(path, producer=None, run_id=None, *, require=False):
    """Validate run-owned record metadata and sidecar fields without rehashing."""
    path = Path(path)
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
    records = {}
    for key, line in zip(_ARTIFACT_RECORD_FIELDS, lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in records:
            raise ValueError(f"artifact record has an invalid schema: {candidate}")
        records[key] = line[len(prefix):]
    context = _record_context(run_id)
    expected_run_id = "" if context is None else context[1]
    if records["PATH"] != str(path.resolve()):
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
        character in records["PRODUCER"] for character in "\t\r\n"
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


def validate_artifact_record(path, producer, run_id):
    """Require and validate one exact run-owned artifact record."""
    record = _read_artifact_record(path, producer, run_id, require=True)
    if record is None:
        raise ValueError(f"artifact record is missing: {path}")
    return record


def _write_artifact_record(path, producer, run_id=None, checksum=None):
    """Publish a run-owned record after a strict artifact publication."""
    context = _record_context(run_id)
    if context is None:
        return None
    if not isinstance(producer, str) or not producer or any(
        character in producer for character in "\t\r\n"
    ):
        raise ValueError(f"invalid artifact producer: {producer!r}")
    path = Path(path)
    sidecar = _full_checksum(path) if checksum is None else checksum
    destination = artifact_record_path(path, context[1])
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
    serialized = (
        f"PATH={path.resolve()}\n"
        f"SIZE={sidecar['SIZE']}\n"
        f"MD5={sidecar['MD5']}\n"
        f"RUN_ID={context[1]}\n"
        f"PRODUCER={producer}\n"
        "STATE=PUBLISHED\n"
    )
    try:
        temporary.write_text(serialized, encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def _validate_feather_frame(frame, path):
    if frame.empty:
        raise ValueError(f"Feather artifact is empty: {path}")
    if isinstance(frame.index, pd.RangeIndex):
        raise ValueError(f"Feather artifact lacks sample identifiers: {path}")
    if frame.index.hasnans or not frame.index.is_unique:
        raise ValueError(f"Feather artifact has invalid sample identifiers: {path}")
    if any(str(value).strip() == "" for value in frame.index):
        raise ValueError(f"Feather artifact has blank sample identifiers: {path}")
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.empty or not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"Feather artifact has no finite numeric features: {path}")


def _ordered_sample_ids(adata):
    """Return unique Sample IDs in their first-appearance obs order."""
    if "Sample" not in adata.obs.columns:
        raise ValueError("benchmark AnnData is missing the 'Sample' column")
    values = adata.obs["Sample"]
    if values.isna().any():
        raise ValueError("benchmark AnnData Sample contains missing values")
    sample_ids = [str(value) for value in values]
    if any(not value.strip() for value in sample_ids):
        raise ValueError("benchmark AnnData Sample contains blank values")
    return list(dict.fromkeys(sample_ids))


def _align_square_frame(frame, sample_ids, path):
    """Align a sample-by-sample frame to the canonical obs sample order."""
    expected = [str(value) for value in sample_ids]
    row_ids = [str(value) for value in frame.index]
    column_ids = [str(value) for value in frame.columns]
    if (
        len(row_ids) != len(expected)
        or len(column_ids) != len(expected)
        or len(set(row_ids)) != len(row_ids)
        or len(set(column_ids)) != len(column_ids)
        or set(row_ids) != set(expected)
        or set(column_ids) != set(expected)
    ):
        raise ValueError(
            f"sample-by-sample Feather output identifiers do not match "
            f"canonical samples: {path}"
        )
    aligned = frame.copy()
    aligned.index = row_ids
    aligned.columns = column_ids
    return aligned.loc[expected, expected]


def recorded_feather_valid(
    path, producer=None, producer_run_id=None, expected_batch_contract=None
):
    """Validate a cache Feather and its semantic frame before reading callers use it."""
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    present = path.exists() or sidecar.exists()
    if not present:
        return False
    if not path.is_file() or not sidecar.is_file():
        raise ValueError(f"Feather cache is incomplete: {path}")
    if expected_batch_contract is None:
        # This full check is mandatory immediately before Feather deserialization.
        _full_checksum(path)
    candidate = _record_candidate(path, producer_run_id)
    if candidate is not None and (candidate.exists() or candidate.is_symlink()):
        _read_artifact_record(path, producer, producer_run_id, require=True)
    if expected_batch_contract is not None:
        recorded_batch_contract = _read_runtime_batch_contract(path)
        _validate_corrected_batch_contract(
            expected_batch_contract,
            recorded_batch_contract,
            "runtime metadata",
        )
        # This full check is mandatory immediately before Feather deserialization.
        _full_checksum(path)
    try:
        frame = pd.read_feather(path)
        _validate_feather_frame(frame, path)
    except Exception as exc:
        if isinstance(exc, ValueError):
            raise
        raise ValueError(f"Feather cache is malformed: {path}") from exc
    return True


def atomic_to_feather(frame, path, producer=None):
    """Write a complete Feather file, checksum, and optional run record."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{path}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    backup = path.with_name(f".{path.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    had_path = path.is_file()
    had_sidecar = sidecar.is_file()
    try:
        _validate_feather_frame(frame, path)
        frame.to_feather(tmp)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"empty Feather output: {tmp}")
        if had_path:
            os.link(path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(tmp, path)
        checksum = {
            "MD5": _file_md5(path),
            "SIZE": str(path.stat().st_size),
            "PATH": str(path),
        }
        sidecar_tmp.write_text(
            f"MD5={checksum['MD5']}\n"
            f"SIZE={checksum['SIZE']}\n"
            f"PATH={checksum['PATH']}\n",
            encoding="utf-8",
        )
        os.replace(sidecar_tmp, sidecar)
        effective_producer = producer or os.environ.get(
            "ECODA_ARTIFACT_PRODUCER"
        )
        if _record_context() is not None:
            if not effective_producer:
                raise ValueError(
                    "run-bound artifact producer is required for Feather output"
                )
            _write_artifact_record(path, effective_producer, checksum=checksum)
    except Exception:
        if backup.exists():
            os.replace(backup, path)
        elif not had_path and path.exists():
            path.unlink()
        if sidecar_backup.exists():
            os.replace(sidecar_backup, sidecar)
        elif not had_sidecar and sidecar.exists():
            sidecar.unlink()
        raise
    finally:
        for temporary in (tmp, sidecar_tmp, backup, sidecar_backup):
            if temporary.exists():
                temporary.unlink()


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
_RUNTIME_CHECKSUM_KEYS = ("MD5", "SIZE", "PATH")


def runtime_metadata_path(output_path):
    """Return the runtime metadata path associated with a Feather output."""
    return Path(f"{Path(output_path)}.runtime.json")


def runtime_metadata_checksum_path(output_path):
    """Return the checksum sidecar path for a runtime metadata JSON file."""
    return Path(f"{runtime_metadata_path(output_path)}.md5")

def _recorded_feather_md5(path):
    """Return the verified MD5 for a recorded Feather artifact.

    Runtime metadata is only meaningful when the output and its checksum
    sidecar were published together.  Verify the exact sidecar schema and
    digest here without parsing the Feather payload; cache validation performs
    the more expensive frame-level checks separately.
    """
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    try:
        output_size = path.stat().st_size
    except OSError as exc:
        raise ValueError(f"cannot publish runtime metadata for missing output: {path}") from exc
    if not path.is_file() or output_size == 0:
        raise ValueError(f"cannot publish runtime metadata for missing output: {path}")
    if not sidecar.is_file():
        raise ValueError(f"cannot publish runtime metadata without checksum sidecar: {sidecar}")
    try:
        lines = sidecar.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read Feather checksum sidecar: {sidecar}") from exc
    if len(lines) != len(_RUNTIME_CHECKSUM_KEYS):
        raise ValueError(f"Feather checksum sidecar has an invalid schema: {sidecar}")
    records = {}
    for key, line in zip(_RUNTIME_CHECKSUM_KEYS, lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in records:
            raise ValueError(f"Feather checksum sidecar has an invalid schema: {sidecar}")
        records[key] = line[len(prefix):]
    artifact_md5 = _file_md5(path)
    if records["PATH"] != str(path):
        raise ValueError(f"Feather checksum sidecar has the wrong PATH: {sidecar}")
    if records["SIZE"] != str(output_size):
        raise ValueError(f"Feather checksum sidecar has the wrong SIZE: {sidecar}")
    if records["MD5"] != artifact_md5:
        raise ValueError(f"Feather checksum sidecar does not match: {sidecar}")
    if path.stat().st_size != output_size:
        raise ValueError(f"Feather artifact changed during checksum validation: {path}")
    return artifact_md5

def _runtime_number(value, field, allow_none=False):
    if value is None and allow_none:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"runtime metadata has invalid {field}: {value!r}")
    number = float(value)
    if not np.isfinite(number) or number < 0:
        raise ValueError(f"runtime metadata has invalid {field}: {value!r}")
    return number

def _runtime_metadata_payload(
    output_path,
    dataset_name,
    method_str,
    time_secs,
    mem_gb,
    batch_contract=None,
):
    output_path = Path(output_path)
    artifact_md5 = _recorded_feather_md5(output_path)
    if not isinstance(dataset_name, str) or not dataset_name:
        raise ValueError(f"runtime metadata has invalid dataset: {dataset_name!r}")
    if not isinstance(method_str, str) or not method_str:
        raise ValueError(f"runtime metadata has invalid method: {method_str!r}")
    time_value = _runtime_number(time_secs, "time_secs")
    memory_value = _runtime_number(mem_gb, "mem_GB", allow_none=True)
    payload = {
        "schema_version": 1,
        "artifact_path": str(output_path),
        "artifact_md5": artifact_md5,
        "dataset": dataset_name,
        "method": method_str,
        "time_secs": time_value,
        "mem_GB": memory_value,
    }
    if batch_contract is not None:
        # Corrected runtime metadata must carry the compact summary; ordinary
        # and uncorrected calls leave the legacy payload untouched.
        if not isinstance(batch_contract, dict):
            raise ValueError("runtime metadata batch_contract must be a mapping")
        _validate_corrected_batch_contract(
            batch_contract,
            batch_contract,
            "runtime metadata",
        )
        forbidden = {
            "composite_values",
            "sample_composite_values",
            "sample_ids",
            "sample_group_ids",
            "tokens",
        }
        vectors = sorted(forbidden.intersection(batch_contract))
        if vectors:
            raise ValueError(
                "runtime metadata batch_contract contains per-cell/sample "
                f"vectors: {', '.join(vectors)}"
            )
        payload["batch_contract"] = dict(batch_contract)
    return payload



def _reject_nonfinite_json_constant(value):
    raise ValueError(f"runtime metadata contains non-finite JSON value: {value}")
def _reject_duplicate_json_keys(pairs):
    """Reject duplicate object keys instead of silently keeping the last."""
    payload = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"runtime metadata has duplicate key: {key}")
        payload[key] = value
    return payload


def _read_runtime_checksum(metadata_path):
    checksum_path = Path(f"{metadata_path}.md5")
    if not checksum_path.is_file():
        raise ValueError(f"runtime metadata checksum is missing: {checksum_path}")
    try:
        lines = checksum_path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read runtime metadata checksum: {checksum_path}") from exc
    if len(lines) != len(_RUNTIME_CHECKSUM_KEYS):
        raise ValueError(f"runtime metadata checksum has an invalid schema: {checksum_path}")
    records = {}
    for key, line in zip(_RUNTIME_CHECKSUM_KEYS, lines):
        prefix = f"{key}="
        if not line.startswith(prefix) or key in records:
            raise ValueError(f"runtime metadata checksum has an invalid schema: {checksum_path}")
        records[key] = line[len(prefix):]
    expected_md5 = records["MD5"]
    if (
        len(expected_md5) != 32
        or expected_md5 != expected_md5.lower()
        or any(char not in "0123456789abcdef" for char in expected_md5)
    ):
        raise ValueError(f"runtime metadata checksum has an invalid MD5: {checksum_path}")
    try:
        expected_size = int(records["SIZE"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"runtime metadata checksum has an invalid SIZE: {checksum_path}") from exc
    if expected_size < 0 or str(expected_size) != records["SIZE"]:
        raise ValueError(f"runtime metadata checksum has an invalid SIZE: {checksum_path}")
    if records["PATH"] != str(metadata_path):
        raise ValueError(f"runtime metadata checksum has the wrong PATH: {checksum_path}")
    if expected_md5 != _file_md5(metadata_path):
        raise ValueError(f"runtime metadata checksum does not match: {checksum_path}")
    if expected_size != metadata_path.stat().st_size:
        raise ValueError(f"runtime metadata checksum has the wrong SIZE: {checksum_path}")


def _read_runtime_batch_contract(output_path):
    """Read and verify the corrected identity beside one Feather artifact."""
    output_path = Path(output_path)
    metadata_path = runtime_metadata_path(output_path)
    if not metadata_path.is_file() or metadata_path.stat().st_size == 0:
        raise ValueError(f"runtime metadata is missing: {metadata_path}")
    _read_runtime_checksum(metadata_path)
    try:
        payload = json.loads(
            metadata_path.read_text(encoding="utf-8"),
            parse_constant=_reject_nonfinite_json_constant,
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"runtime metadata is malformed: {metadata_path}") from exc
    required_fields = _RUNTIME_METADATA_FIELDS | _RUNTIME_METADATA_IDENTITY_FIELDS
    if not isinstance(payload, dict) or set(payload) != required_fields:
        raise ValueError(f"runtime metadata has an invalid schema: {metadata_path}")
    if type(payload["schema_version"]) is not int or payload["schema_version"] != 1:
        raise ValueError(f"runtime metadata has an invalid schema_version: {metadata_path}")
    if payload["artifact_path"] != str(output_path):
        raise ValueError(f"runtime metadata has the wrong artifact_path: {metadata_path}")
    artifact_checksum = _read_checksum_sidecar(output_path)
    if (
        not isinstance(payload["artifact_md5"], str)
        or payload["artifact_md5"] != artifact_checksum["MD5"]
    ):
        raise ValueError(f"runtime metadata has the wrong artifact_md5: {metadata_path}")
    if (
        not isinstance(payload["dataset"], str)
        or not payload["dataset"]
        or not isinstance(payload["method"], str)
        or not payload["method"]
    ):
        raise ValueError(f"runtime metadata has invalid dataset/method: {metadata_path}")
    _runtime_number(payload["time_secs"], "time_secs")
    _runtime_number(payload["mem_GB"], "mem_GB", allow_none=True)
    identity = payload["batch_contract"]
    if not isinstance(identity, dict):
        raise ValueError(
            f"runtime metadata batch_contract is not an object: {metadata_path}"
        )
    return identity

def publish_runtime_metadata(
    output_path,
    dataset_name,
    method_str,
    time_secs,
    mem_gb,
    producer=None,
    producer_run_id=None,
    batch_contract=None,
):
    """Atomically publish runtime metadata after its output is complete."""
    output_path = Path(output_path)
    effective_producer = producer or os.environ.get("ECODA_ARTIFACT_PRODUCER")
    record_context = _record_context(producer_run_id)
    if record_context is not None:
        if not effective_producer:
            raise ValueError(
                "run-bound artifact producer is required for runtime metadata"
            )
        _read_artifact_record(
            output_path,
            effective_producer,
            producer_run_id,
            require=True,
        )
    metadata_path = runtime_metadata_path(output_path)
    checksum_path = runtime_metadata_checksum_path(output_path)
    metadata_tmp = metadata_path.with_name(
        f".{metadata_path.name}.tmp.{os.getpid()}"
    )
    checksum_tmp = checksum_path.with_name(
        f".{checksum_path.name}.tmp.{os.getpid()}"
    )
    metadata_backup = metadata_path.with_name(
        f".{metadata_path.name}.previous.{os.getpid()}"
    )
    checksum_backup = checksum_path.with_name(
        f".{checksum_path.name}.previous.{os.getpid()}"
    )
    payload = _runtime_metadata_payload(
        output_path,
        dataset_name,
        method_str,
        time_secs,
        mem_gb,
        batch_contract=batch_contract,
    )
    serialized = (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    had_metadata = metadata_path.is_file()
    had_checksum = checksum_path.is_file()
    try:
        metadata_tmp.write_text(serialized, encoding="utf-8")
        checksum_tmp.write_text(
            f"MD5={_file_md5(metadata_tmp)}\n"
            f"SIZE={metadata_tmp.stat().st_size}\n"
            f"PATH={metadata_path}\n",
            encoding="utf-8",
        )
        if had_metadata:
            os.link(metadata_path, metadata_backup)
        if had_checksum:
            os.link(checksum_path, checksum_backup)
        os.replace(metadata_tmp, metadata_path)
        os.replace(checksum_tmp, checksum_path)
        if record_context is not None:
            metadata_checksum = {
                "MD5": _file_md5(metadata_path),
                "SIZE": str(metadata_path.stat().st_size),
                "PATH": str(metadata_path),
            }
            _write_artifact_record(
                metadata_path,
                effective_producer,
                run_id=producer_run_id,
                checksum=metadata_checksum,
            )
    except Exception:
        if metadata_backup.exists():
            os.replace(metadata_backup, metadata_path)
        elif not had_metadata and metadata_path.exists():
            metadata_path.unlink()
        if checksum_backup.exists():
            os.replace(checksum_backup, checksum_path)
        elif not had_checksum and checksum_path.exists():
            checksum_path.unlink()
        raise
    finally:
        for temporary in (
            metadata_tmp,
            checksum_tmp,
            metadata_backup,
            checksum_backup,
        ):
            if temporary.exists():
                temporary.unlink()


def read_runtime_metadata(
    output_path,
    dataset_name,
    method_str,
    producer=None,
    producer_run_id=None,
    expected_batch_contract=None,
):
    """Read and strictly validate metadata for a valid Feather cache hit."""
    output_path = Path(output_path)
    effective_producer = producer or os.environ.get("ECODA_ARTIFACT_PRODUCER")
    record_context = _record_context(producer_run_id)
    if record_context is not None and not effective_producer:
        raise ValueError(
            "run-bound artifact producer is required for runtime metadata"
        )
    if expected_batch_contract is None:
        feather_valid = recorded_feather_valid(
            output_path,
            producer=effective_producer,
            producer_run_id=producer_run_id,
        )
    else:
        feather_valid = recorded_feather_valid(
            output_path,
            producer=effective_producer,
            producer_run_id=producer_run_id,
            expected_batch_contract=expected_batch_contract,
        )
    if not feather_valid:
        raise ValueError(f"output Feather is not a valid recorded artifact: {output_path}")
    metadata_path = runtime_metadata_path(output_path)
    if not metadata_path.is_file() or metadata_path.stat().st_size == 0:
        raise ValueError(f"runtime metadata is missing: {metadata_path}")
    _read_runtime_checksum(metadata_path)
    metadata_record = _record_candidate(metadata_path, producer_run_id)
    metadata_record_present = metadata_record is not None and (
        metadata_record.exists() or metadata_record.is_symlink()
    )
    if metadata_record_present:
        _read_artifact_record(
            metadata_path,
            effective_producer,
            producer_run_id,
            require=True,
        )
    elif producer_run_id is not None:
        raise ValueError(f"runtime metadata artifact record is missing: {metadata_path}")
    try:
        payload = json.loads(
            metadata_path.read_text(encoding="utf-8"),
            parse_constant=_reject_nonfinite_json_constant,
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"runtime metadata is malformed: {metadata_path}") from exc
    expected_runtime_fields = _RUNTIME_METADATA_FIELDS
    if expected_batch_contract is not None:
        expected_runtime_fields |= _RUNTIME_METADATA_IDENTITY_FIELDS
    if not isinstance(payload, dict) or set(payload) != expected_runtime_fields:
        raise ValueError(f"runtime metadata has an invalid schema: {metadata_path}")
    if payload["schema_version"] != 1 or type(payload["schema_version"]) is not int:
        raise ValueError(f"runtime metadata has an invalid schema_version: {metadata_path}")
    if payload["artifact_path"] != str(output_path):
        raise ValueError(f"runtime metadata has the wrong artifact_path: {metadata_path}")
    artifact_md5 = payload["artifact_md5"]
    if (
        not isinstance(artifact_md5, str)
        or len(artifact_md5) != 32
        or artifact_md5 != artifact_md5.lower()
        or any(char not in "0123456789abcdef" for char in artifact_md5)
        or artifact_md5 != _file_md5(output_path)
    ):
        raise ValueError(f"runtime metadata has the wrong artifact_md5: {metadata_path}")
    if payload["dataset"] != dataset_name or not isinstance(payload["dataset"], str):
        raise ValueError(f"runtime metadata has the wrong dataset: {metadata_path}")
    if payload["method"] != method_str or not isinstance(payload["method"], str):
        raise ValueError(f"runtime metadata has the wrong method: {metadata_path}")
    if expected_batch_contract is not None:
        _validate_corrected_batch_contract(
            expected_batch_contract,
            payload["batch_contract"],
            f"runtime metadata {method_str}",
        )
    _runtime_number(payload["time_secs"], "time_secs")
    _runtime_number(payload["mem_GB"], "mem_GB", allow_none=True)
    return payload
def _validate_execution_measurements(frame, path):
    for column, allow_missing in (("time_secs", False), ("mem_GB", True)):
        raw = frame[column]
        values = pd.to_numeric(raw, errors="coerce")
        missing = raw.isna()
        if (not allow_missing and missing.any()) or values[~missing].isna().any():
            raise ValueError(f"execution log has invalid numeric values: {path}")
        numeric = values[~missing].to_numpy(dtype=float)
        if not np.isfinite(numeric).all() or (numeric < 0).any():
            raise ValueError(f"execution log has invalid numeric values: {path}")


def _validate_execution_log_frame(frame, path):
    required = {"dataset", "method", "time_secs", "mem_GB"}
    if set(frame.columns) != required or frame.empty:
        raise ValueError(f"execution log has an invalid schema: {path}")
    if frame[["dataset", "method"]].isna().any().any() or \
       frame[["dataset", "method"]].astype(str).apply(lambda column: column.str.strip() == "").any().any():
        raise ValueError(f"execution log has blank identifiers: {path}")
    if frame[["dataset", "method"]].duplicated().any():
        raise ValueError(f"execution log has duplicate identifiers: {path}")
    _validate_execution_measurements(frame, path)

def _read_recorded_execution_log(
    path, producer=None, producer_run_id=None, *, require_record=False
):
    """Validate an execution log fully before deserializing it."""
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    present = path.exists() or sidecar.exists()
    if not present:
        return None
    if not path.is_file() or not sidecar.is_file():
        raise ValueError(f"execution log is incomplete: {path}")
    _full_checksum(path)
    candidate = _record_candidate(path, producer_run_id)
    if candidate is not None and (candidate.exists() or candidate.is_symlink()):
        _read_artifact_record(path, producer, producer_run_id, require=True)
    elif require_record:
        raise ValueError(f"execution log artifact record is missing: {path}")
    try:
        frame = pd.read_feather(path)
    except Exception as exc:
        raise ValueError(f"execution log is malformed: {path}") from exc
    _validate_execution_log_frame(frame, path)
    return frame

def execution_log_atomic_to_feather(frame, path, producer=None):
    """Atomically write the documented indexless execution-log schema."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    sidecar = Path(f"{path}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    backup = path.with_name(f".{path.name}.previous.{os.getpid()}")
    sidecar_backup = sidecar.with_name(f".{sidecar.name}.previous.{os.getpid()}")
    had_path = path.is_file()
    had_sidecar = sidecar.is_file()
    try:
        _validate_execution_log_frame(frame, path)
        frame.reset_index(drop=True).to_feather(tmp)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"empty execution log: {tmp}")
        if had_path:
            os.link(path, backup)
        if had_sidecar:
            os.link(sidecar, sidecar_backup)
        os.replace(tmp, path)
        checksum = {
            "MD5": _file_md5(path),
            "SIZE": str(path.stat().st_size),
            "PATH": str(path),
        }
        sidecar_tmp.write_text(
            f"MD5={checksum['MD5']}\n"
            f"SIZE={checksum['SIZE']}\n"
            f"PATH={checksum['PATH']}\n",
            encoding="utf-8",
        )
        os.replace(sidecar_tmp, sidecar)
        effective_producer = producer or os.environ.get(
            "ECODA_EXECUTION_LOG_PRODUCER", "stage5_execution_log"
        )
        if _record_context() is not None:
            _write_artifact_record(path, effective_producer, checksum=checksum)
    except Exception:
        if backup.exists():
            os.replace(backup, path)
        elif not had_path and path.exists():
            path.unlink()
        if sidecar_backup.exists():
            os.replace(sidecar_backup, sidecar)
        elif not had_sidecar and sidecar.exists():
            sidecar.unlink()
        raise
    finally:
        for temporary in (tmp, sidecar_tmp, backup, sidecar_backup):
            if temporary.exists():
                temporary.unlink()


# scPoli is imported lazily (get_scpoli) so that MrVI/PILOT runs never touch
# scarches: scarches 0.6.1 (its final release) does `from anndata import
# AnnData, read`, but `anndata.read` was removed in anndata >= 0.12 (the
# pinned 0.12.19). The shim below restores it as the documented alias
# `read_h5ad` — scarches only calls `read()` on .h5ad files, so this is a
# faithful drop-in. See https://github.com/theislab/scarches.
def get_scpoli():
    import anndata as ad

    if not hasattr(ad, "read"):
        ad.read = ad.read_h5ad  # scarches compat (anndata >= 0.12 removed `read`)
    from scarches.models.scpoli import scPoli

    return scPoli


def get_pilotgm():
    """Lazy import of pilotgm with a shim for its non-relative internal imports.

    Upstream packaging bug (pilotgm 0.1.1): `pilotgm/model/GMVAE.py` does
    `from networks.Networks import *` (also `losses`/`metrics`), but the
    modules ship inside the package (`pilotgm/networks/...`), so a plain
    `import pilotgm` fails with ModuleNotFoundError. The shim inserts the
    pilotgm package directory into sys.path so the top-level names resolve.

    The entry is left on sys.path ON PURPOSE: gmmvae_wasserstein_distance
    dispatches `compute_emd` to loky workers via joblib, and unpickling the
    task re-imports pilotgm in the worker. loky spawns copy the parent's
    sys.path, so removing the entry after the import would break the worker
    import ("failed to un-serialize" / ModuleNotFoundError: networks).
    Collision check (py-cpu env, 2026-08-14): no other installed package
    imports `core`, `model`, `networks`, `losses` or `metrics` top-level, so
    the extra path entry is inert.
    """
    import importlib.metadata as md

    for f in md.files("pilotgm"):
        if str(f).endswith("pilotgm/__init__.py"):
            pkg_dir = f.locate().parent
            break
    if pkg_dir is None:
        raise ImportError("pilotgm package files not found (importlib.metadata)")
    if str(pkg_dir) not in sys.path:
        sys.path.insert(0, str(pkg_dir))
    import pilotgm

    return pilotgm


# ---------------------------------------------------------------------------
# Execution time / memory logging
# ---------------------------------------------------------------------------
def peak_rss_gb():
    """Peak resident set size of this process in GB.

    getrusage().ru_maxrss units: KB on Linux, bytes on macOS.
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / 1024.0 / 1024.0 / 1024.0
    return rss / 1024.0 / 1024.0

def report_gpu_memory(method_str):
    """Print peak CUDA memory and device capacity for resource profiling."""
    if not torch.cuda.is_available():
        return
    try:
        torch.cuda.synchronize()
        props = torch.cuda.get_device_properties(torch.cuda.current_device())
        gib = float(1024**3)
        print(
            f"GPU_MEMORY_PROFILE method={method_str} "
            f"device={props.name} total_GiB={props.total_memory / gib:.2f} "
            f"peak_allocated_GiB={torch.cuda.max_memory_allocated() / gib:.2f} "
            f"peak_reserved_GiB={torch.cuda.max_memory_reserved() / gib:.2f}",
            flush=True,
        )
    except Exception as exc:
        print(f"GPU_MEMORY_PROFILE_ERROR method={method_str}: {exc}", flush=True)



GPU_BACKED_METHODS = frozenset(("mrvi", "scpoli"))


def validate_gpu_execution(method, device, combo=None):
    """Validate the resource class selected by the Stage 5 submitter."""
    if method == "mrvi" and device == "cpu":
        if combo is None or combo == "hvg2000":
            raise RuntimeError(
                "MrVI CPU execution requires an explicit non-default "
                "--combo; the default hvg2000 run is H200-only"
            )
        return
    if method not in GPU_BACKED_METHODS:
        return
    if device != "cuda":
        raise RuntimeError(
            f"{method} is GPU-backed and requires --device cuda; "
            "refusing implicit CPU/auto execution"
        )
    if not torch.cuda.is_available():
        raise RuntimeError(
            f"{method} requested CUDA but torch.cuda.is_available() is False; "
            "refusing CPU fallback"
        )

_USE_PEAK_RSS = object()


def log_execution_time(
    dataset_name,
    method_str,
    time_secs,
    log_file,
    mem_gb=_USE_PEAK_RSS,
    producer=None,
    producer_run_id=None,
):
    """Append/overwrite one (dataset, method) row in the per-task exec log."""
    if mem_gb is _USE_PEAK_RSS:
        mem_gb = peak_rss_gb()
    new_row = pd.DataFrame(
        {
            "dataset": [dataset_name],
            "method": [method_str],
            "time_secs": [float(time_secs)],
            "mem_GB": [mem_gb],
        }
    )
    effective_producer = producer or os.environ.get(
        "ECODA_EXECUTION_LOG_PRODUCER", "stage5_execution_log"
    )
    df_existing = _read_recorded_execution_log(
        log_file,
        producer=effective_producer,
        producer_run_id=producer_run_id,
    )
    if df_existing is None:
        df_final = new_row
    else:
        mask = (df_existing["dataset"] == dataset_name) & (
            df_existing["method"] == method_str
        )
        if mask.any():
            df_existing = df_existing[~mask]
        df_final = pd.concat([df_existing, new_row], ignore_index=True)
    execution_log_atomic_to_feather(
        df_final, log_file, producer=effective_producer
    )


def replay_runtime_metadata(
    output_path,
    dataset_name,
    method_str,
    log_file,
    producer=None,
    producer_run_id=None,
    expected_batch_contract=None,
):
    payload = read_runtime_metadata(
        output_path,
        dataset_name,
        method_str,
        producer=producer,
        producer_run_id=producer_run_id,
        expected_batch_contract=expected_batch_contract,
    )
    log_execution_time(
        dataset_name,
        method_str,
        payload["time_secs"],
        log_file,
        mem_gb=payload["mem_GB"],
        producer=(
            os.environ.get("ECODA_EXECUTION_LOG_PRODUCER")
            or "stage5_execution_log"
        ),
        producer_run_id=producer_run_id,
    )
    return payload


# ---------------------------------------------------------------------------
# Combo resolution (legacy qmd rules, datasets.json-driven)
# ---------------------------------------------------------------------------
def scpoli_dims_for(n_hvg, res_label):
    """scPoli embedding dims for an (n_hvg, resolution) combo (qmd rules)."""
    if res_label == "_highres":
        return [2, 3, 5, 10, 15] if n_hvg == 2000 else [15]
    if res_label == "_lowres" and n_hvg == 2000:
        return [15]
    return []


def run_wass_combo_for(n_hvg, res_label):
    """Whether a Wasserstein-distance method (PILOT/QOT/PILOT-GM-VAE) runs
    for an (n_hvg, resolution) combo (qmd rules)."""
    if res_label == "_highres":
        return True
    return res_label == "_lowres" and n_hvg == 2000


# Default (main-method) combos — constants.R method_label_map_main and the
# notebook's exec-time figure: MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
# PILOT_hvg2000_highres, QOT_hvg2000_highres, PILOT-GM-VAE_hvg2000_highres.
DEFAULT_HVG = 2000
DEFAULT_SCPOLI_DIM = 15
DEFAULT_RES_LABEL = "_highres"


def legacy_method_label(method, n_hvg, res_label, payload):
    """Return the exact execution-log label used by the legacy benchmark."""
    if method == "mrvi":
        return f"MrVI_hvg{n_hvg}"
    if method == "scpoli":
        return f"scPoli_hvg{n_hvg}_dims{payload}{res_label}"
    if method == "qot":
        return f"QOT_hvg{n_hvg}{res_label}"
    if method == "pilotgm":
        return f"PILOT-GM-VAE_hvg{n_hvg}{res_label}"
    return f"PILOT_hvg{n_hvg}{res_label}"


def is_default_combo(method, combo):
    n, res_label, _, payload, _, _ = combo
    if method == "mrvi":
        return n == DEFAULT_HVG
    if method == "scpoli":
        return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL and payload == DEFAULT_SCPOLI_DIM
    if method in ("pilot", "qot", "pilotgm"):
        return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL
    return False


def canonical_combo_token(method, combo):
    """Return the strict matrix-manifest token for one generated combo."""
    n_hvg, res_label, _, payload, _, _ = combo
    try:
        resolution = {"_lowres": "lowres", "_highres": "highres"}[res_label]
    except KeyError as exc:
        raise ValueError(
            f"Cannot derive a canonical combo token from resolution {res_label!r}"
        ) from exc

    if method == "mrvi":
        return f"hvg{n_hvg}"
    if method == "scpoli":
        if payload is None:
            raise ValueError("scPoli combo is missing its embedding dimension")
        return f"hvg{n_hvg}_{resolution}_dims{payload}"
    if method in ("pilot", "qot", "pilotgm"):
        return f"hvg{n_hvg}_{resolution}"
    raise ValueError(f"Cannot derive a canonical combo token for method {method!r}")


def select_requested_combo(method, combos, requested_combo):
    """Select exactly one generated combo, failing closed on bad tokens."""
    matches = [
        combo
        for combo in combos
        if canonical_combo_token(method, combo) == requested_combo
    ]
    if not matches:
        available = sorted(
            {canonical_combo_token(method, combo) for combo in combos}
        )
        raise ValueError(
            f"Unknown --combo {requested_combo!r} for method {method!r}; "
            f"available canonical combos: {available}"
        )
    if len(matches) != 1:
        raise ValueError(
            f"Ambiguous --combo {requested_combo!r} for method {method!r}: "
            f"{len(matches)} generated combos match"
        )
    return matches


def top_n_hvg_genes(adata, n):
    """Top-n genes from the stored hvg_rank (set by 1.1.1_preprocess.py)."""
    ranks = adata.var["hvg_rank"].dropna().sort_values()
    if len(ranks) < n:
        raise ValueError(
            f"Only {len(ranks)} genes have a stored hvg_rank, but {n} were "
            f"requested. Re-run preprocessing with a larger HVG size."
        )
    return list(ranks.index[:n])




def use_counts_layer(sub, method, ds_name):
    """Point X at the validated raw counts layer for MrVI/scPoli.

    The preprocessed h5ad has log-normalized X; both models need the raw
    counts vaulted by base_preprocessing in layers["counts"]. The contract is
    checked before this function runs, so missing counts are an error rather
    than a normalized-expression fallback.
    """
    counts = sub.layers["counts"]
    if method == "scpoli":
        if sp.issparse(counts):
            sub.X = counts.toarray().astype("float32", copy=False)
        else:
            sub.X = np.asarray(counts, dtype="float32")
    else:
        sub.X = counts
    return sub


# ---------------------------------------------------------------------------
# Method bodies (qmd semantics preserved)
# ---------------------------------------------------------------------------
def run_mrvi(adata, device, output_path, batch_key=None, producer=None):
    """MrVI local sample distances with an optional technical batch key."""
    adata.obs["dummy_col"] = np.zeros(adata.n_obs)
    setup_kwargs = {"sample_key": "Sample"}
    if batch_key is not None:
        setup_kwargs["batch_key"] = batch_key
    MRVI.setup_anndata(adata, **setup_kwargs)
    model = MRVI(adata)
    model.train(max_epochs=50, accelerator=device)
    dists = model.get_local_sample_distances(
        keep_cell=False, groupby="dummy_col", batch_size=32
    )
    df_dists = dists["dummy_col"].isel(dummy_col_name=0).to_pandas()
    df_dists = _align_square_frame(
        df_dists, _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path, producer=producer)


def run_scpoli(adata, ct_col, dim, output_path, producer=None):
    """scPoli conditional sample embeddings for one embedding dim."""
    # scPoli requires a cell-type label for EVERY cell, but datasets whose
    # declared ct columns are the pipeline annotation columns (Lee/Zhang:
    # layer1/layer2) carry NaN for cells the annotators left unclassified
    # (HiTME covers immune subsets only, ~66-88% of cells). scarches'
    # label_encoder runs np.unique() on the mixed str/NaN object column,
    # which crashes under numpy 2.x ("'<' not supported between instances of
    # 'str' and 'float'"). Fill with an explicit "Unknown" class: every cell
    # stays in the output (per-cell embedding rows stay aligned) and
    # unannotated cells get their own scPoli prototype. No-op on complete
    # columns (all other datasets use fully-covered author annotations).
    n_na = int(adata.obs[ct_col].isna().sum())
    if n_na:
        print(f"Filling {n_na}/{adata.n_obs} missing values in '{ct_col}' "
              "with 'Unknown' (scPoli requires a label per cell).")
        col = adata.obs[ct_col]
        # anndata stores string obs columns as pandas Categorical by default;
        # fillna() cannot introduce an undefined category, so register
        # 'Unknown' first.
        if isinstance(col.dtype, pd.CategoricalDtype):
            col = col.cat.add_categories("Unknown")
        adata.obs[ct_col] = col.fillna("Unknown")
    scPoli = get_scpoli()
    scpoli_model = scPoli(
        adata=adata,
        condition_keys="Sample",
        cell_type_keys=ct_col,
        embedding_dims=dim,
        recon_loss="nb",
    )
    scpoli_model.train(
        n_epochs=50,
        pretraining_epochs=40,
        early_stopping_kwargs={
            "early_stopping_metric": "val_prototype_loss",
            "mode": "min",
            "threshold": 0,
            "patience": 20,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        },
        eta=5,
    )
    adata_emb = scpoli_model.get_conditional_embeddings()
    df_embs = pd.DataFrame(adata_emb.X, index=adata_emb.obs_names)
    df_embs.columns = [f"Dim_{i + 1}" for i in range(df_embs.shape[1])]
    atomic_to_feather(df_embs, output_path, producer=producer)


def resolve_pass_embedding_key(adata, view, n_hvg):
    """Resolve one explicitly declared embedding without compatibility fallback."""
    if view == "benchmark_analysis":
        key = f"X_pca_benchmark_analysis_hvg{n_hvg}"
    elif view == "batch_effect_uncorrected":
        key = f"X_pca_batch_effect_uncorrected_hvg{n_hvg}"
    elif view == "batch_effect_corrected":
        key = f"X_pca_harmony_batch_effect_corrected_hvg{n_hvg}"
    else:
        raise ValueError(f"Unknown preprocessing view: {view}")
    if key not in adata.obsm:
        raise KeyError(f"Required embedding {key!r} is missing from adata.obsm")
    return key


def run_pilot(adata, ct_col, view, n_hvg, output_path, producer=None):
    """PILOT Wasserstein sample distances on the exact view embedding."""
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
    emb = adata.obsm[emb_key]
    # PILOT (pilotpy>=2.0.x) requires a named-columns pandas DataFrame in
    # obsm (Trajectory.extract_data_anno_scRNA_from_h5ad accesses .columns);
    # the preprocess step stores scanpy's plain ndarray instead.
    if not hasattr(emb, "columns"):
        emb = pd.DataFrame(
            emb,
            index=adata.obs_names,
            columns=[f"PCA_{i + 1}" for i in range(emb.shape[1])],
        )
    adata.obsm[emb_key] = emb
    # Same NaN cell-type guard as QOT/PILOT-GM-VAE/scPoli (see run_scpoli):
    # PILOT's cost_matrix() treats NaN labels as a pseudo-cell-type whose
    # centroid is the median of zero cells -> NaN entries in the cost matrix
    # -> ot.emd2 collapses (all-zero EMD distances written on HPC; segfault
    # locally). Lee/Zhang (HiTME immune-only annotation, 12-33% NaN) hit this.
    fill_unknown_ct(adata, ct_col, "PILOT")
    # pilotpy writes hard-coded Results_PILOT/plots relative to cwd.  A
    # read-only project bind cannot host those scratch plots, and they are not
    # benchmark artifacts; isolate only this third-party side effect.
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="pilot_") as temp_dir:
        try:
            os.chdir(temp_dir)
            pl.tl.wasserstein_distance(
                adata,
                emb_matrix=emb_key,
                clusters_col=ct_col,
                sample_col="Sample",
                status="Sample",
            )
        finally:
            os.chdir(cwd)
    df_dists = _align_square_frame(
        adata.uns["EMD_df"], _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path, producer=producer)


def fill_unknown_ct(adata, ct_col, method):
    """Fill NaN cell-type labels with an explicit "Unknown" class.

    Same rationale as scPoli (see run_scpoli): the QOT script filters
    `Cell_type != 'Unknown'` and NaN would otherwise leak into the GMM
    groupby keys. No-op on complete columns.
    """
    n_na = int(adata.obs[ct_col].isna().sum())
    if n_na:
        print(f"Filling {n_na}/{adata.n_obs} missing values in '{ct_col}' "
              f"with 'Unknown' ({method}).")
        col = adata.obs[ct_col]
        if isinstance(col.dtype, pd.CategoricalDtype):
            col = col.cat.add_categories("Unknown")
        adata.obs[ct_col] = col.fillna("Unknown")


def run_qot(adata, ct_col, view, n_hvg, output_path, producer=None):
    """QOT Wasserstein sample distances on the preprocessed obsm PCA.

    Runs the vendored qot_utils_re.py (PennShenLab/QOT @ 28cd529880c1, two
    hotfixes in Gaussian_Mixture_Representation — see the file header and
    docs/qot_hotfixes.md). Lazy import so the phate dependency is only
    touched for QOT runs.
    """
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
    import qot_utils_re

    fill_unknown_ct(adata, ct_col, "QOT")
    # Distinct temp column: Extract_Info renames {type_cell, id, progession}
    # in ONE dict — passing "Sample" for both id and progession collapses the
    # duplicate dict key, BOTH columns get renamed to 'status', no 'sampleID'
    # survives and the GMM groupby raises KeyError. The temp column also
    # keeps bio labels out of the distance path.
    adata.obs["_bench_prog"] = adata.obs["Sample"]
    qot_utils_re.Run_QOT(
        adata,
        gene_matrix=emb_key,
        type_cell=ct_col,
        id_col="Sample",
        progession="_bench_prog",
        dataset_type="rna",
        num_components_list=[1],
        random_state=2,
        min_samples_for_gmm=0,
        qot_method="cosine",
        normalized_set=False,
    )
    samples = _ordered_sample_ids(adata)
    # Plain object strings: anndata obs columns are categorical by default and
    # a categorical DataFrame index would be written to the feather as
    # categorical (pyarrow 24/25 pandas-compat cannot read that back:
    # "data type 'categorical' not understood"). PILOT's EMD_df index is
    # plain object — keep the identical layout.
    samples = np.asarray(samples, dtype=object)
    df_dists = pd.DataFrame(
        adata.uns["QOT_Distance"], index=samples, columns=samples
    )
    df_dists = _align_square_frame(
        df_dists, samples, output_path
    )
    atomic_to_feather(df_dists, output_path, producer=producer)


def _stabilize_pilotgm_covariances(adata):
    """Regularize undefined empirical covariances before PILOT-GM distances.

    The upstream representation helper uses ``np.cov`` for each assigned
    component. A component with one cell has an undefined covariance and
    yields NaNs; interpreting that component as zero empirical covariance,
    followed by the package's diagonal regularization, is the finite
    single-observation limit. Means and weights remain untouched and must
    already be finite.
    """
    representations = adata.uns.get("GMVAE_Representation")
    if not isinstance(representations, dict) or not representations:
        raise ValueError("PILOT-GM-VAE did not produce sample representations")
    for sample_id, representation in representations.items():
        means = np.asarray(representation["means"], dtype=float)
        weights = np.asarray(representation["weights"], dtype=float)
        if (
            not np.isfinite(means).all()
            or not np.isfinite(weights).all()
            or weights.ndim != 1
            or weights.size == 0
            or float(weights.sum()) <= 0
        ):
            raise ValueError(
                f"PILOT-GM-VAE produced invalid means/weights for sample {sample_id!r}"
            )
        covariances = np.asarray(representation["covariances"], dtype=float)
        if (
            covariances.ndim != 3
            or covariances.shape[1] != covariances.shape[2]
            or covariances.shape[0] != means.shape[0]
        ):
            raise ValueError(
                f"PILOT-GM-VAE produced invalid covariances for sample {sample_id!r}"
            )
        stable = []
        for covariance in covariances:
            if np.isinf(covariance).any():
                raise ValueError(
                    f"PILOT-GM-VAE produced infinite covariance for sample {sample_id!r}"
                )
            covariance = np.nan_to_num(covariance, nan=0.0)
            covariance = (covariance + covariance.T) / 2.0
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            eigenvalues = np.clip(eigenvalues, 0.0, None)
            stable.append((eigenvectors * eigenvalues) @ eigenvectors.T)
        representation["covariances"] = np.asarray(stable, dtype=float)


def _run_pilotgm_distance_with_stable_covariances(pilotgm, adata, *, emb_key):
    """Run PILOT-GM after finite PSD covariance stabilization."""
    import importlib

    pilotgm_core = importlib.import_module("pilotgm.core")
    original_representation = pilotgm_core.gaussian_mixture_vae_representation

    def stabilized_representation(*args, **kwargs):
        result = original_representation(*args, **kwargs)
        _stabilize_pilotgm_covariances(result)
        return result

    pilotgm_core.gaussian_mixture_vae_representation = stabilized_representation
    try:
        pilotgm.gmmvae_wasserstein_distance(
            adata,
            emb_matrix=emb_key,
            clusters_col="component_assignment",
            sample_col="Sample",
            status="_bench_status",
            wass_dis=True,
            covariance_type="full",
            epsilon=1e-3,
        )
    finally:
        pilotgm_core.gaussian_mixture_vae_representation = original_representation


def run_pilotgm(adata, ct_col, view, n_hvg, output_path, ds_name, device, producer=None):
    """PILOT-GM-VAE Wasserstein sample distances on the preprocessed obsm PCA.

    Runs the `pilotgm` PyPI package (CostaLab/PILOT-GM-VAE, BIB 2025):
    `train_gmvae` (50 epochs, num_classes = n unique cell types) then
    `gmmvae_wasserstein_distance`. The whole pilotgm block runs inside a
    node-local tempdir: `train_gmvae` hardcodes `./trained_models/<ds>/` and
    saves weights — running it from the repo root (worker cwd on HPC) would
    pollute the repo, and running it from the scratch output dir would
    pollute the NAS sync (the submit tail rsyncs benchmark/ wholesale).
    Weights are ephemeral by design (load_weights=False; retries re-train
    from scratch).
    """
    emb_key = resolve_pass_embedding_key(adata, view, n_hvg)
    emb = adata.obsm[emb_key]
    # train_gmvae needs torch.tensor(obsm[key]) (fails on a pandas
    # DataFrame: "could not determine the shape of object type 'DataFrame'"
    # with torch >= 2.x), while extract_data_anno_scRNA_from_h5ad (in the
    # distance step) needs `.columns`. Store the plain ndarray for training,
    # swap in the named-columns DataFrame only for
    # gmmvae_wasserstein_distance — a DataFrame at that point is also
    # joblib-picklable (a __main__-defined ndarray subclass broke loky's
    # task serialization).
    if hasattr(emb, "columns"):
        emb = np.asarray(emb)
    adata.obsm[emb_key] = emb

    # Unlike PILOT/QOT, PILOT-GM's distance API uses the model-generated
    # `component_assignment` column, not the biological cell-type labels. The
    # configured annotation only preserves the historical component-count
    # choice; do not turn missing labels into an extra rare component.
    ct_values = adata.obs[ct_col].astype("string")
    valid_ct = ct_values.notna() & (ct_values.str.strip() != "")
    num_classes = max(2, int(ct_values[valid_ct].nunique()))
    # Distinct temp column for `status`: same duplicate-key rename bug as QOT
    # (gmmvae_wasserstein_distance renames the last three columns via a dict;
    # sample_col == status == "Sample" would collapse the keys).
    adata.obs["_bench_status"] = adata.obs["Sample"]

    pilotgm = get_pilotgm()
    # Plain `device == "cuda"` would silently run CPU on GPU nodes under the
    # default --device auto.
    use_cuda = device == "cuda" or (device == "auto" and torch.cuda.is_available())

    cwd = os.getcwd()
    tmp_dir = tempfile.mkdtemp(prefix="pilotgm_")
    try:
        os.chdir(tmp_dir)
        pilotgm.train_gmvae(
            adata,
            dataset_name=ds_name,
            pca_key=emb_key,
            labels_column=None,
            epochs=50,
            num_classes=num_classes,
            cuda=use_cuda,
            gpuID=0,
            load_weights=False,
            save_model=True,
            seed=1,
        )
        # Swap in the named-columns DataFrame for the distance step (see
        # the comment at the top of this function).
        adata.obsm[emb_key] = pd.DataFrame(
            emb,
            index=adata.obs_names,
            columns=[f"PCA_{i + 1}" for i in range(emb.shape[1])],
        )
        # The configured cell-type column determines the requested number of
        # mixture components, but the package's distance step uses the
        # model-generated component_assignment column, not those annotations.
        _run_pilotgm_distance_with_stable_covariances(
            pilotgm,
            adata,
            emb_key=emb_key,
        )
    finally:
        os.chdir(cwd)
    df_dists = _align_square_frame(
        adata.uns["EMD_df"], _ordered_sample_ids(adata), output_path
    )
    atomic_to_feather(df_dists, output_path, producer=producer)


# ---------------------------------------------------------------------------
# Orchestration
def process_dataset(args, ds_name, entry):
    """Run all combos of the requested method for one dataset.

    Loads the h5ad once per task (not per combo); skips only a valid,
    recorded Feather cache unless --force. Existing incomplete or invalid
    artifacts fail closed instead of being silently recomputed.
    """
    artifact_producer = stage5_artifact_producer(args.method)
    # The method supplied by the Stage 5 submitter is authoritative for all
    # benchmark artifacts, including calls into helpers that omit producer.
    os.environ["ECODA_ARTIFACT_PRODUCER"] = artifact_producer
    view_name = args.view
    analysis_pass = getattr(args, "analysis_pass", None)
    requested_combo = getattr(args, "combo", None)
    if requested_combo is not None and analysis_pass is not None:
        raise ValueError("--combo is only supported for ordinary benchmark runs")

    high_resolution_only = bool(getattr(args, "high_resolution_only", False)) or analysis_pass is not None
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = getattr(args, "log_file", None)
    if log_file is None:
        log_file = output_dir / "execution_times.feather"

    if view_name not in entry["views"]:
        raise ValueError(
            f"Dataset '{ds_name}' has no '{view_name}' view in datasets.json."
        )
    expected_view_for_pass = {
        "uncorrected": {"batch_effect_uncorrected"},
        "corrected": {"batch_effect_corrected"},
    }
    if analysis_pass is not None:
        if analysis_pass not in expected_view_for_pass:
            raise ValueError(f"Unknown analysis pass: {analysis_pass}")
        if view_name not in expected_view_for_pass[analysis_pass]:
            raise ValueError(
                f"analysis pass {analysis_pass!r} requires one of "
                f"{sorted(expected_view_for_pass[analysis_pass])!r}"
            )
        if analysis_pass == "corrected" and entry.get("batch_col") is None:
            raise ValueError(
                "corrected batch-effect view requires a confirmed columns.batch"
            )
    view_output = entry["views"][view_name]["output_file"]
    input_path = Path(args.input_dir) / view_output
    if not input_path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_path}")

    corrected_batch_keys = ()
    corrected_validation_summary = None
    expected_h5ad_batch_contract = None
    expected_stage5_batch_contract = None
    if analysis_pass == "corrected":
        corrected_batch_keys = normalize_batch_keys(entry.get("batch_col"))
        corrected_validation_summary = _read_corrected_validation_summary(
            input_path,
            corrected_batch_keys,
        )
        expected_h5ad_batch_contract = _corrected_h5ad_identity(
            corrected_batch_keys
        )
        expected_stage5_batch_contract = _corrected_stage5_identity(
            args.method,
            corrected_batch_keys,
            corrected_validation_summary,
        )
        if args.method == "pilotgm":
            raise ValueError(
                "PILOT-GM-VAE is not scheduled for corrected batch-effect runs"
            )

    lowres_col = entry.get("cell_type_low_res")
    highres_col = entry.get("cell_type_high_res")
    technical_batch = (
        corrected_batch_keys[0] if len(corrected_batch_keys) == 1 else None
    )

    def output_name(suffix, n, res_label=None, extension="dists"):
        if analysis_pass is not None:
            return (
                f"{ds_name}_batch_effect_{analysis_pass}_hvg{n}_highres_"
                f"{suffix}_{extension}.feather"
            )
        return f"{ds_name}_hvg{n}{res_label or ''}_{suffix}_{extension}.feather"

    tiers = [("_highres", highres_col)] if high_resolution_only else [
        ("_lowres", lowres_col),
        ("_highres", highres_col),
    ]
    combos = []
    if args.method == "mrvi":
        if lowres_col is None:
            if requested_combo is not None:
                raise ValueError(
                    f"Unknown --combo {requested_combo!r} for method 'mrvi': "
                    "MrVI has no runnable combo without cell_type_low_res"
                )

            print(f"WARNING: {ds_name}: cell_type_low_res is None; skipping MrVI.")
            return
        for n in args.hvg:
            out_name = (
                output_name("mrvi", n, extension="dists")
                if analysis_pass is not None
                else f"{ds_name}_hvg{n}_mrvi_dists.feather"
            )
            combos.append((n, "_highres" if analysis_pass else "_lowres",
                           None, None, run_mrvi, out_name))

    elif args.method == "scpoli":
        for res_label, ct_col in tiers:
            if ct_col is None:
                continue
            for n in args.hvg:
                for dim in scpoli_dims_for(n, res_label):
                    out_name = (
                        f"{ds_name}_hvg{n}{res_label}_scpoli_dims{dim}_embs.feather"
                    )
                    combos.append((n, res_label, ct_col, dim, run_scpoli, out_name))

    elif args.method in ("pilot", "qot", "pilotgm"):
        suffix = {"pilot": "pilot", "qot": "qot", "pilotgm": "pilotgm"}[args.method]
        for res_label, ct_col in tiers:
            if ct_col is None:
                continue
            for n in args.hvg:
                if run_wass_combo_for(n, res_label):
                    out_name = (
                        output_name(suffix, n, res_label=res_label)
                        if analysis_pass is not None
                        else f"{ds_name}_hvg{n}{res_label}_{suffix}_dists.feather"
                    )
                    combos.append((n, res_label, ct_col, None, None, out_name))

    if requested_combo is not None:
        combos = select_requested_combo(
            args.method, combos, requested_combo
        )

    # Defaults-first ordering: ru_maxrss peak RSS is monotonic within a
    # process, so combos run earlier report the least bloated mem_GB (memory
    # leaks / allocator retention from earlier combos would otherwise inflate
    # the defaults' rows). Stable sort: non-default combos keep their order.
    # If --hvg excludes the default size, no default combo exists and the
    # sort is a no-op — no behavior change.
    combos.sort(key=lambda c: 0 if is_default_combo(args.method, c) else 1)

    pending = []
    method_labels = {}
    for n, res_label, ct_col, payload, run_fn, out_name in combos:
        out_path = output_dir / out_name
        method_str = legacy_method_label(args.method, n, res_label, payload)
        method_labels[out_path] = method_str
        if expected_stage5_batch_contract is None:
            cache_valid = (
                not args.force
                and recorded_feather_valid(out_path, producer=artifact_producer)
            )
        else:
            cache_valid = (
                not args.force
                and recorded_feather_valid(
                    out_path,
                    producer=artifact_producer,
                    expected_batch_contract=expected_stage5_batch_contract,
                )
            )
        if cache_valid:
            replay_runtime_metadata(
                out_path,
                ds_name,
                method_str,
                log_file,
                producer=artifact_producer,
                expected_batch_contract=expected_stage5_batch_contract,
            )
            print(f"Already processed and validated: {out_name}")
            continue
        pending.append((n, res_label, ct_col, payload, run_fn, out_path))

    if not pending:
        return

    print(f"Loading {input_path} ...")
    source_shape = None
    if args.method in ("pilot", "qot", "pilotgm"):
        _validate_h5ad_path(
            input_path,
            args.view,
            args.method,
            expected_h5ad_batch_contract,
        )
        obs_columns = {"Sample"}
        obs_columns.update(
            str(ct_col) for _, _, ct_col, _, _, _ in pending if ct_col is not None
        )
        embedding_keys = []
        for n, _, _, _, _, _ in pending:
            if args.view == "batch_effect_corrected":
                key = f"X_pca_harmony_batch_effect_corrected_hvg{n}"
            elif args.view == "batch_effect_uncorrected":
                key = f"X_pca_batch_effect_uncorrected_hvg{n}"
            else:
                key = f"X_pca_benchmark_analysis_hvg{n}"
            if key not in embedding_keys:
                embedding_keys.append(key)
        adata = load_h5ad_counts_free(
            input_path,
            sorted(obs_columns),
            embedding_keys,
        )
        print("COUNTS_ACCESS=none; loaded selected obs/obsm into minimal AnnData")
    elif args.method in ("mrvi", "scpoli"):
        _validate_h5ad_path(
            input_path,
            args.view,
            args.method,
            expected_h5ad_batch_contract,
        )
        obs_columns = {"Sample"}
        if args.method == "scpoli":
            obs_columns.update(
                str(ct_col)
                for _, _, ct_col, _, _, _ in pending
                if ct_col is not None
            )
        if args.method == "mrvi" and corrected_batch_keys:
            obs_columns.update(corrected_batch_keys)
        elif technical_batch is not None:
            obs_columns.add(str(technical_batch))
        max_hvg = max(n for n, _, _, _, _, _ in pending)
        selected_genes = read_h5ad_hvg_genes(input_path, max_hvg)
        adata = load_h5ad_counts_subset(
            input_path,
            selected_genes,
            sorted(obs_columns),
        )
        source_shape = tuple(
            int(value) for value in adata.uns.pop("_ecoda_source_shape")
        )
        print(
            "COUNTS_ACCESS=selected; loaded stored HVG counts into minimal AnnData"
        )
    else:
        adata = sc.read_h5ad(str(input_path), backed="r")
        _validate_h5ad_content(
            adata,
            args.view,
            args.method,
            expected_h5ad_batch_contract,
        )
        adata = adata.to_memory()
    profile_shape = source_shape or (adata.n_obs, adata.n_vars)
    selected_suffix = (
        f" selected_genes={adata.n_vars}" if source_shape is not None else ""
    )
    print(
        f"INPUT_PROFILE bytes={input_path.stat().st_size} "
        f"cells={profile_shape[0]} genes={profile_shape[1]}{selected_suffix}",
        flush=True,
    )

    if "Sample" not in adata.obs.columns:
        raise ValueError(
            f"Cannot find standardized sample column 'Sample' in obs of {input_path}."
        )
    if technical_batch is not None and technical_batch not in adata.obs.columns:
        raise ValueError(
            f"Confirmed batch column '{technical_batch}' not found in obs of {input_path}"
        )
    if args.method == "mrvi" and analysis_pass == "corrected":
        validation = validate_batch_metadata(adata.obs, corrected_batch_keys)
        loaded_summary = build_batch_validation_summary(
            validation,
            corrected_validation_summary["correction_mode"],
            corrected_validation_summary["correction_formula"],
        )
        if loaded_summary != corrected_validation_summary:
            raise ValueError(
                "loaded corrected batch metadata differs from H5AD validation_summary"
            )

    for n, res_label, ct_col, payload, run_fn, out_path in pending:
        # HVG gene subset per combo, done FIRST (before any dense conversion)
        # so only the small n_obs x n_hvg matrix is materialized. PILOT uses
        # the stored obsm directly, no subset needed.
        if args.method in ("mrvi", "scpoli"):
            genes = top_n_hvg_genes(adata, n)
            sub = adata[:, genes].copy()
            sub = use_counts_layer(sub, args.method, ds_name)
        else:
            sub = adata

        if ct_col is not None and ct_col not in sub.obs.columns:
            raise ValueError(
                f"Cell type column '{ct_col}' not found in obs of {ds_name} "
                f"(available: {list(sub.obs.columns)})."
            )
        run_batch_key = technical_batch
        temporary_batch_key = None

        # Exact legacy method strings (constants.R + notebook recodes depend
        # on them); these were derived before the cache scan as well.
        method_str = method_labels[out_path]

        print(f"Processing {method_str} ...")
        start_time = time.time()
        profile_gpu = args.method in GPU_BACKED_METHODS and args.device == "cuda"
        if profile_gpu:
            torch.cuda.reset_peak_memory_stats()
        try:
            if (
                args.method == "mrvi"
                and analysis_pass == "corrected"
                and len(corrected_batch_keys) >= 2
            ):
                batch_composite = build_batch_composite(
                    sub.obs,
                    corrected_batch_keys,
                )
                sub.obs = batch_composite.frame
                temporary_batch_key = batch_composite.column_name
                run_batch_key = temporary_batch_key
            if args.method == "mrvi":
                run_mrvi(
                    sub,
                    args.device,
                    out_path,
                    batch_key=run_batch_key,
                )
            elif args.method == "scpoli":
                run_scpoli(sub, ct_col, payload, out_path)
            elif args.method == "qot":
                run_qot(sub, ct_col, args.view, n, out_path)
            elif args.method == "pilotgm":
                run_pilotgm(
                    sub,
                    ct_col,
                    args.view,
                    n,
                    out_path,
                    ds_name,
                    args.device,
                )
            else:
                run_pilot(sub, ct_col, args.view, n, out_path)
        finally:
            if (
                temporary_batch_key is not None
                and temporary_batch_key in sub.obs.columns
            ):
                del sub.obs[temporary_batch_key]
            if profile_gpu:
                report_gpu_memory(method_str)
        exec_time = time.time() - start_time
        mem_gb = peak_rss_gb()
        publish_runtime_metadata(
            out_path,
            ds_name,
            method_str,
            exec_time,
            mem_gb,
            producer=artifact_producer,
            batch_contract=expected_stage5_batch_contract,
        )
        log_execution_time(
            ds_name,
            method_str,
            exec_time,
            log_file,
            mem_gb=mem_gb,
            producer=(
                os.environ.get("ECODA_EXECUTION_LOG_PRODUCER")
                or "stage5_execution_log"
            ),
        )
        print(f"  -> Saved: {out_path} ({exec_time:.2f}s, "
              f"{mem_gb:.2f} GB peak RSS)")
        gc.collect()
def main():
    parser = argparse.ArgumentParser(
        description="Run Python benchmark methods (MrVI/scPoli/PILOT/QOT/"
                    "PILOT-GM-VAE) on a preprocessed benchmark view h5ad."
    )
    parser.add_argument("--config_path", required=True,
                        help="Path to datasets.json")
    parser.add_argument("--ds_name", required=True,
                        help="Dataset key in datasets.json")
    parser.add_argument("--view", default="benchmark_analysis",
                        help="View name (default: benchmark_analysis)")
    parser.add_argument("--analysis_pass", default=None,
                        choices=["uncorrected", "corrected"],
                        help="Batch-effect pass; requires the matching explicit view")
    parser.add_argument("--high_resolution_only", action="store_true",
                        help="Run only the configured high-resolution tier")
    parser.add_argument("--method", required=True,
                        choices=["mrvi", "scpoli", "pilot", "qot", "pilotgm"],
                        help="Benchmark method to run")
    parser.add_argument("--input_dir", required=True,
                        help="Directory holding the preprocessed view h5ad")
    parser.add_argument("--output_dir", required=True,
                        help="Feather output dir (created if missing)")
    parser.add_argument("--log_file", default=None,
                        help="Per-task execution-time log feather "
                             "(default: <output_dir>/execution_times.feather "
                             "for local runs)")
    parser.add_argument("--hvg", nargs="+", type=int, default=[1000, 2000, 3000],
                        help="HVG sizes to run (default: 1000 2000 3000)")
    parser.add_argument(
        "--combo",
        default=None,
        help=(
            "Run exactly one generated canonical parameter combo "
            "(for example, hvg2000_highres_dims15)"
        ),
    )
    parser.add_argument("--force", action="store_true", default=False,
                        help="Recompute combos whose output feather already exists")
    parser.add_argument("--device", default="auto",
                        choices=["auto", "cpu", "cuda"],
                        help="Training device; default MrVI/scPoli combos "
                             "require CUDA, non-default MrVI may use CPU")
    args = parser.parse_args()

    args.hvg = sorted(set(args.hvg))
    if args.analysis_pass is not None:
        args.hvg = [2000]
    validate_gpu_execution(args.method, args.device, args.combo)

    scvi.settings.seed = 0
    print("scvi-tools version:", scvi.__version__)
    print("torch.cuda.is_available():", torch.cuda.is_available())

    config = read_datasets_json(args.config_path, view=args.view)
    if args.ds_name not in config:
        raise ValueError(f"'{args.ds_name}' is not a dataset in {args.config_path} "
                         f"with a '{args.view}' view.")
    process_dataset(args, args.ds_name, config[args.ds_name])
    print("Processing complete!")


if __name__ == "__main__":
    main()
