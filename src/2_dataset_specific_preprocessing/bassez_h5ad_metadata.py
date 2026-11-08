#!/usr/bin/env python3
"""Create a source-bound Bassez HR-only H5AD derivative on HPC scratch.

Expression arrays are copied as file bytes and compared in bounded HDF5
slices, never materialized as AnnData or whole matrices. Only cellSubType and
new, explicitly separate repair provenance may differ from the parent H5AD.
The HPC wrapper must call ecoda_require_input_ownership before and after this
worker, and ecoda_validate_input_artifact for a supplied producer record.
Legacy checksum-proof mode never invents a producer or PUBLISHED state.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.utils.py.h5ad_source_identity import read_obs_column_values

HR = "cellSubType"
LR = "cellType"
REPAIR = "uns/author_annotation_repair"
POLICY = "bassez_cellSubType_fallback_v1"
MISSING = ("", "na", "nan", "none", "unknown")
BUFFER_BYTES = 8 * 1024 * 1024


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def _digest(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(BUFFER_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def _absolute(path, *, exists=True):
    path = Path(path)
    if (not path.is_absolute() or any(char in str(path) for char in "\t\r\n")
            or os.path.normpath(path) != str(path) or path.resolve(strict=exists) != path):
        raise ValueError(f"path must be canonical and absolute: {path}")
    return path


def _record(path, keys):
    path = _absolute(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"provenance record is not a regular file: {path}")
    raw = path.read_bytes()
    if b"\r" in raw:
        raise ValueError(f"invalid record line endings: {path}")
    text = raw.decode("utf-8")
    lines = text.splitlines()
    if not text.endswith("\n") or len(lines) != len(keys):
        raise ValueError(f"invalid provenance record: {path}")
    fields = {}
    for line, key in zip(lines, keys):
        if not line.startswith(f"{key}=") or not line[len(key) + 1:]:
            raise ValueError(f"invalid {key} in provenance record: {path}")
        fields[key] = line[len(key) + 1:]
    return fields


def _source_identity(source, expected_md5, parent_record, parent_run_id, parent_producer,
                     scratch_root, legacy_input_proof=None):
    source = _absolute(source)
    scratch_root = _absolute(scratch_root)
    if not source.is_relative_to(scratch_root) or not source.is_file() or source.is_symlink():
        raise ValueError("Bassez input must be a regular H5AD on the explicit HPC scratch root")
    if not re.fullmatch(r"[0-9a-f]{32}", expected_md5):
        raise ValueError("expected input MD5 must contain 32 lowercase hexadecimal digits")
    if (parent_record is None) == (legacy_input_proof is None):
        raise ValueError("select exactly one parent producer record or explicit legacy input proof")
    if parent_record is not None:
        if not isinstance(parent_run_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", parent_run_id):
            raise ValueError("invalid parent producer run ID")
        key = hashlib.sha256(str(source).encode()).hexdigest()[:32]
        canonical_record = scratch_root / "_ecoda_runs" / parent_run_id / "manifests/artifacts" / f"{key}.record"
        if _absolute(parent_record) != canonical_record:
            raise ValueError("parent producer record is not the canonical run-owned artifact record")
        producer = parent_producer or "stage3"
        fields = _record(parent_record, ("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE"))
        expected = {"PATH": str(source), "SIZE": str(source.stat().st_size), "MD5": expected_md5,
                    "RUN_ID": parent_run_id, "PRODUCER": producer, "STATE": "PUBLISHED"}
        provenance = {"provenance_kind": "published_producer", "producer_run_id": parent_run_id,
                      "producer": producer, "record": str(_absolute(parent_record)),
                      "record_sha256": _digest(parent_record),
                      "ownership_validation": "wrapper_requires_validate_input_artifact_before_and_after"}
    else:
        if parent_run_id is not None or parent_producer is not None:
            raise ValueError("legacy inputs must not claim a producer run or PUBLISHED record")
        fields = _record(legacy_input_proof, ("MD5", "SIZE", "PATH"))
        expected = {"MD5": expected_md5, "SIZE": str(source.stat().st_size), "PATH": str(source)}
        provenance = {"provenance_kind": "legacy_validated_input",
                      "checksum_proof": str(_absolute(legacy_input_proof)),
                      "checksum_proof_sha256": _digest(legacy_input_proof),
                      "ownership_validation": "wrapper_requires_no_active_writer_before_and_after"}
    if fields != expected or source.stat().st_size <= 0:
        raise ValueError("parent H5AD provenance differs from the selected input")
    if _digest(source, "md5") != expected_md5:
        raise ValueError("parent H5AD checksum differs from the selected input")
    return {"path": str(source), "size": source.stat().st_size, "md5": expected_md5, **provenance}


def _annotation(value):
    if value is None or (isinstance(value, (float, np.floating)) and math.isnan(float(value))):
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if not isinstance(value, str):
        raise ValueError("Bassez author annotations must contain strings or missing values")
    return value


def _missing(value):
    return value is None or value.strip(" \t\r\n").lower() in MISSING


def annotations(path, sample_column="Sample"):
    """Match bassez_cellsubtype_utils.R's sentinel rules; keep valid text verbatim."""
    with h5py.File(path, "r") as handle:
        if "obs" not in handle or "uns" not in handle:
            raise ValueError("processed Bassez H5AD requires obs and uns groups")
        obs = handle["obs"]
        index = obs.attrs.get("_index", "_index")
        if isinstance(index, bytes):
            index = index.decode("utf-8")
        cells = read_obs_column_values(obs, index)
        samples = read_obs_column_values(obs, sample_column)
        high = [_annotation(value) for value in read_obs_column_values(obs, HR)]
        low = [_annotation(value) for value in read_obs_column_values(obs, LR)]
    n = len(cells)
    if not n or len(set(map(str, cells))) != n or any(len(values) != n for values in (samples, high, low)):
        raise ValueError("Bassez cell IDs or metadata lengths are invalid")
    if any(not isinstance(value, str) or not value.strip() or value.strip().lower() == "nan"
           for value in [*cells, *samples]):
        raise ValueError("Bassez cell/sample IDs are missing, blank, or non-text")
    replace = [_missing(value) for value in high]
    invalid = [str(cells[index]) for index in range(n) if replace[index] and _missing(low[index])]
    if invalid:
        raise ValueError(f"invalid cellType fallback for missing cellSubType at cells: {invalid[:10]}")
    filled = [low[index] if replace[index] else high[index] for index in range(n)]
    coverage = {"n_cells": n, "n_samples": len(set(map(str, samples))),
                "n_replaced": sum(replace), "before_valid_categories": len(set(value for value in high if not _missing(value))),
                "after_categories": len(set(filled)), "sample_column": sample_column,
                "fallback_counts": dict(sorted(Counter(low[index] for index in range(n) if replace[index]).items()))}
    return filled, coverage


def _dtype(dtype):
    if h5py.check_ref_dtype(dtype) is not None:
        raise ValueError("HDF5 object references cannot be proven unchanged by this repair")
    string = h5py.check_string_dtype(dtype)
    variable = h5py.check_vlen_dtype(dtype)
    if variable is not None and variable not in (str, bytes):
        raise ValueError("variable-length numeric datasets cannot be streamed with a bounded buffer")
    if dtype.names:
        for name in dtype.names:
            _dtype(dtype.fields[name][0])
    return {"dtype": dtype.str, "descr": repr(dtype.descr),
            "string": (string.encoding, string.length) if string else None,
            "vlen": repr(h5py.check_vlen_dtype(dtype)), "enum": h5py.check_enum_dtype(dtype)}


def _feed(digest, value):
    """Hash typed values, including variable-length strings, without pointer bytes."""
    if value is None:
        digest.update(b"N")
    elif isinstance(value, (str, bytes)):
        raw = value.encode("utf-8") if isinstance(value, str) else value
        digest.update(("S" if isinstance(value, str) else "B").encode() + str(len(raw)).encode() + b":" + raw)
    elif isinstance(value, dict):
        digest.update(b"D")
        _feed(digest, sorted(value))
        for key in sorted(value):
            _feed(digest, value[key])
    elif isinstance(value, (list, tuple)):
        digest.update(f"L{len(value)}:".encode())
        for item in value:
            _feed(digest, item)
    elif isinstance(value, h5py.Empty):
        _feed(digest, {"null_dtype": _dtype(value.dtype)})
    else:
        array = np.asarray(value)
        digest.update(b"A")
        _feed(digest, _dtype(array.dtype))
        _feed(digest, list(array.shape))
        if array.dtype.names:
            # Structured Scanpy metadata can contain vlen string fields. Hash
            # field values, never the object-pointer bytes or padding.
            for name in array.dtype.names:
                _feed(digest, name)
                _feed(digest, array[name])
        elif array.dtype.kind in "OU":
            for item in array.flat:
                _feed(digest, item)
        else:
            digest.update(array.tobytes(order="C"))


def _slices(shape, itemsize, buffer_bytes=BUFFER_BYTES):
    if not all(shape):
        return
    block = list(shape)
    target = max(1, buffer_bytes // max(1, itemsize))
    for axis in range(len(block)):
        rest = math.prod(block[:axis] + block[axis + 1:])
        block[axis] = min(block[axis], max(1, target // max(1, rest)))
    for start in product(*(range(0, size, step) for size, step in zip(shape, block))):
        yield tuple(slice(offset, min(offset + step, size)) for offset, step, size in zip(start, block, shape))


def protected_tree(path, buffer_bytes=BUFFER_BYTES):
    """Prove every other HDF5 dataset, storage property, group and attribute unchanged."""
    digest = hashlib.sha256()
    counts = {"objects": 0, "datasets": 0}

    def visit(node, name):
        excluded = any(name == root or name.startswith(root + "/") for root in (f"obs/{HR}", REPAIR))
        if not excluded:
            counts["objects"] += 1
            _feed(digest, {"path": name, "kind": "dataset" if isinstance(node, h5py.Dataset) else "group"})
            for key in sorted(node.attrs):
                _feed(digest, key)
                _feed(digest, node.attrs[key])
        if isinstance(node, h5py.Group):
            for key in sorted(node):
                if not isinstance(node.get(key, getlink=True), h5py.HardLink):
                    raise ValueError("external and soft HDF5 links are not supported for metadata repair")
                child = node[key]
                if h5py.h5o.get_info(child.id).rc != 1:
                    raise ValueError("aliased HDF5 objects are not supported for metadata repair")
                visit(child, f"{name}/{key}" if name else key)
        elif not excluded:
            if node.is_virtual:
                raise ValueError("virtual HDF5 datasets are not supported for metadata repair")
            counts["datasets"] += 1
            _feed(digest, {"shape": node.shape, "dtype": _dtype(node.dtype), "chunks": node.chunks,
                           "compression": node.compression, "compression_opts": node.compression_opts,
                           "shuffle": node.shuffle, "fletcher32": node.fletcher32,
                           "scaleoffset": node.scaleoffset, "fillvalue": node.fillvalue, "maxshape": node.maxshape})
            if node.shape is None:
                return
            if not node.shape:
                _feed(digest, node[()])
            elif all(node.shape):
                # Vlen text is metadata; cap the number of objects per read as
                # pointer itemsize does not bound the decoded string lengths.
                budget = min(buffer_bytes, 1024 * max(1, node.dtype.itemsize)) if node.dtype.hasobject else buffer_bytes
                for selection in _slices(node.shape, node.dtype.itemsize, budget):
                    _feed(digest, node[selection])

    with h5py.File(path, "r") as handle:
        visit(handle, "")
    return {"schema": "h5ad_protected_tree_v1", "sha256": digest.hexdigest(),
            "buffer_bytes": buffer_bytes, **counts}


def _destinations(output, run_root, scratch_root):
    run_root = _absolute(run_root)
    scratch_root = _absolute(scratch_root)
    if run_root.parent != scratch_root / "_ecoda_runs" or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_root.name):
        raise ValueError("repair run root must be an exact run below scratch/_ecoda_runs")
    output = _absolute(output, exists=False)
    expected_parent = run_root / "artifacts/datasets/Bassez/output"
    if output.parent != expected_parent or output.suffix != ".h5ad":
        raise ValueError("repair output must be a run-owned Bassez H5AD derivative")
    identities = {}
    for name in ("source.manifest", "runtime.identity"):
        path = _absolute(run_root / "manifests" / name)
        if not path.is_file() or path.stat().st_size <= 0:
            raise ValueError(f"missing repair run identity: {path}")
        identities[str(path)] = _digest(path)
    return output, identities


def _write_repaired_column(handle, filled):
    categories = sorted(set(filled))
    codes = {value: index for index, value in enumerate(categories)}
    del handle["obs"][HR]  # Only the new temporary derivative is writable.
    group = handle["obs"].create_group(HR)
    group.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
    values = group.create_dataset("categories", data=np.asarray(categories, dtype=object), dtype=h5py.string_dtype("utf-8"))
    values.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
    values = group.create_dataset("codes", data=np.asarray([codes[value] for value in filled], dtype=np.int32))
    values.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})


def create_derivative(*, input_file, output_file, expected_md5, scratch_root, run_root,
                      parent_record=None, parent_run_id=None, parent_producer=None,
                      legacy_input_proof=None, sample_column="Sample"):
    output, run_identities = _destinations(output_file, run_root, scratch_root)
    sidecar, report_path = Path(f"{output}.md5"), Path(f"{output}.repair.json")
    if any(path.exists() or path.is_symlink() for path in (output, sidecar, report_path)):
        raise FileExistsError("existing repair outputs are preserved; validate them or select a new run")
    parent = _source_identity(input_file, expected_md5, parent_record, parent_run_id, parent_producer,
                              scratch_root, legacy_input_proof)
    filled, coverage = annotations(input_file, sample_column)
    if coverage["n_replaced"] == 0:
        raise ValueError("Bassez HR annotations are already complete; reuse the validated parent")
    with h5py.File(input_file, "r") as handle:
        if REPAIR in handle:
            raise ValueError("parent already contains author-annotation repair provenance; validate/reuse it")
    unchanged = protected_tree(input_file)
    identity = {"schema": "ecoda_bassez_metadata_repair_v1", "dataset": "Bassez", "view": "benchmark_analysis",
                "run_id": Path(run_root).name, "parent": parent, "run_identities": run_identities,
                "worker_sha256": _digest(__file__), "policy": POLICY, "missing_tokens": list(MISSING),
                "changed_obs_column": HR, "coverage": coverage, "unchanged_protected_tree": unchanged}
    identity_sha256 = hashlib.sha256(_json(identity).encode()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    owned = []

    def temporary(payload=None):
        descriptor, name = tempfile.mkstemp(prefix=f".{output.name}.repair-", dir=output.parent)
        path = Path(name)
        owned.append(path)
        with os.fdopen(descriptor, "wb") as handle:
            if payload is not None:
                handle.write(payload)
        return path

    try:
        candidate = temporary()
        shutil.copyfile(input_file, candidate)
        with h5py.File(candidate, "r+") as handle:
            _write_repaired_column(handle, filled)
            group = handle["uns"].create_group("author_annotation_repair")
            group.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
            for key, value in (("identity_json", _json(identity)), ("identity_sha256", identity_sha256)):
                node = group.create_dataset(key, data=value, dtype=h5py.string_dtype("utf-8"))
                node.attrs.update({"encoding-type": "string", "encoding-version": "0.2.0"})
        actual, after = annotations(candidate, sample_column)
        if actual != filled or after["n_replaced"] != 0 or protected_tree(candidate) != unchanged:
            raise ValueError("Bassez repair changed protected data or failed HR completeness")
        if _source_identity(input_file, expected_md5, parent_record, parent_run_id, parent_producer,
                            scratch_root, legacy_input_proof) != parent:
            raise ValueError("parent provenance changed during repair")
        if any(_digest(path) != digest for path, digest in run_identities.items()):
            raise ValueError("repair run identity changed during repair")
        output_md5 = _digest(candidate, "md5")
        report = {"identity": identity, "identity_sha256": identity_sha256,
                  "output": {"path": str(output), "size": candidate.stat().st_size, "md5": output_md5}}
        checksum = temporary(f"MD5={output_md5}\nSIZE={candidate.stat().st_size}\nPATH={output}\n".encode())
        proof = temporary((_json(report) + "\n").encode())
        # Publish the H5AD last. A reader never sees it before both proofs are
        # complete. Failed partial publication remains immutable evidence.
        for source, destination in ((proof, report_path), (checksum, sidecar), (candidate, output)):
            os.link(source, destination)
        return report
    finally:
        for path in owned:
            path.unlink(missing_ok=True)


def validate_derivative(**kwargs):
    output, run_identities = _destinations(kwargs["output_file"], kwargs["run_root"], kwargs["scratch_root"])
    parent = _source_identity(kwargs["input_file"], kwargs["expected_md5"], kwargs.get("parent_record"),
                              kwargs.get("parent_run_id"), kwargs.get("parent_producer"), kwargs["scratch_root"],
                              kwargs.get("legacy_input_proof"))
    report_path = _absolute(f"{output}.repair.json")
    report = json.loads(report_path.read_text())
    identity = report["identity"]
    expected_hash = hashlib.sha256(_json(identity).encode()).hexdigest()
    if (report["identity_sha256"] != expected_hash or identity["parent"] != parent
            or identity["schema"] != "ecoda_bassez_metadata_repair_v1"
            or identity["dataset"] != "Bassez" or identity["view"] != "benchmark_analysis"
            or identity["run_id"] != Path(kwargs["run_root"]).name or identity["run_identities"] != run_identities
            or identity["policy"] != POLICY or identity["changed_obs_column"] != HR
            or identity["missing_tokens"] != list(MISSING)):
        raise ValueError("repair provenance does not match the selected parent/run/policy")
    size, md5 = output.stat().st_size, _digest(output, "md5")
    if report["output"] != {"path": str(output), "size": size, "md5": md5}:
        raise ValueError("repair output identity differs from its proof")
    if _record(f"{output}.md5", ("MD5", "SIZE", "PATH")) != {"MD5": md5, "SIZE": str(size), "PATH": str(output)}:
        raise ValueError("repair checksum sidecar differs")
    with h5py.File(output, "r") as handle:
        if (handle[f"{REPAIR}/identity_json"].asstr()[()] != _json(identity)
                or handle[f"{REPAIR}/identity_sha256"].asstr()[()] != expected_hash):
            raise ValueError("H5AD repair provenance differs from its proof")
    expected, coverage = annotations(kwargs["input_file"], kwargs.get("sample_column", "Sample"))
    actual, after = annotations(output, kwargs.get("sample_column", "Sample"))
    if expected != actual or after["n_replaced"] or identity["coverage"] != coverage:
        raise ValueError("repaired HR metadata differs from the approved fallback")
    if protected_tree(kwargs["input_file"]) != identity["unchanged_protected_tree"] or protected_tree(output) != identity["unchanged_protected_tree"]:
        raise ValueError("protected H5AD datasets/attributes differ from the parent")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("audit", "create", "validate"))
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--expected-md5", required=True)
    parent = parser.add_mutually_exclusive_group(required=True)
    parent.add_argument("--parent-record")
    parent.add_argument("--legacy-input-proof")
    parser.add_argument("--parent-run-id")
    parser.add_argument("--parent-producer", help="Strict producer-record mode only; defaults to stage3")
    parser.add_argument("--scratch-root", required=True)
    parser.add_argument("--sample-column", default="Sample")
    parser.add_argument("--output-file")
    parser.add_argument("--run-root")
    parser.add_argument("--require-complete", action="store_true")
    args = vars(parser.parse_args())
    command, require_complete = args.pop("command"), args.pop("require_complete")
    if command == "audit":
        _source_identity(args["input_file"], args["expected_md5"], args["parent_record"], args["parent_run_id"],
                         args["parent_producer"], args["scratch_root"], args["legacy_input_proof"])
        _, coverage = annotations(args["input_file"], args["sample_column"])
        if require_complete and coverage["n_replaced"]:
            raise ValueError(f"Bassez producer still has {coverage['n_replaced']} missing HR labels")
        print(_json(coverage))
    else:
        if args["output_file"] is None or args["run_root"] is None:
            parser.error("create/validate require --output-file and --run-root")
        result = create_derivative(**args) if command == "create" else validate_derivative(**args)
        print(_json(result))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, TypeError) as exc:
        raise SystemExit(f"ERROR: {exc}")
