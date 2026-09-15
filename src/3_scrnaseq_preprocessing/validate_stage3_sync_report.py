#!/usr/bin/env python3
"""Validate corrected Stage 3 H5ADs and write a run-bound sync report."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

SOURCE_ROOT = Path(__file__).resolve().parents[2]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))


from src.utils.py.batch_contract import build_batch_contract_identity
from src.utils.py.benchmark_h5ad_contract import validate_benchmark_h5ad_path


_REPORT_CONTRACT = "batch_effect_corrected_h5ad_v1"
_FIELD_RE = re.compile(r"^[A-Z][A-Z0-9_]*=[^\r\n]*$")
_COMPONENT_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _absolute(path: str, label: str) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        raise ValueError(f"{label} must be absolute: {path}")
    return candidate


def _regular_file(path: Path, label: str, *, nonempty: bool = True) -> Path:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be a regular non-symlink file: {path}")
    if nonempty and path.stat().st_size <= 0:
        raise ValueError(f"{label} must be nonempty: {path}")
    return path


def _regular_dir(path: Path, label: str) -> Path:
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"{label} must be a regular non-symlink directory: {path}")
    return path


def _within(path: Path, root: Path, label: str) -> None:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
    except ValueError as exc:
        raise ValueError(f"{label} escapes its run root: {path}") from exc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_identity(path: Path, label: str) -> dict[str, Any]:
    _regular_file(path, label)
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "size": path.stat().st_size,
    }

def _prior_terminal_identity(path: Path, run_id: str) -> dict[str, Any]:
    identity = _file_identity(path, "prior Stage 3 terminal status")
    lines = path.read_text(encoding="utf-8").splitlines()
    state_lines = [line for line in lines if line.startswith("STATE=")]
    run_lines = [line for line in lines if line.startswith("RUN_ID=")]
    if state_lines != ["STATE=FAIL"] or run_lines != [f"RUN_ID={run_id}"]:
        raise ValueError(
            "prior Stage 3 terminal status must be STATE=FAIL for the target run"
        )
    return identity

def _manifest(path: Path, label: str) -> dict[str, str]:
    _regular_file(path, label)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"{label} is empty: {path}")
    result: dict[str, str] = {}
    for line in lines:
        if not _FIELD_RE.fullmatch(line):
            raise ValueError(f"{label} has malformed identity row: {path}")
        key, value = line.split("=", 1)
        if not value or key in result:
            raise ValueError(f"{label} has duplicate or empty identity fields: {path}")
        result[key] = value
    return result


def _snapshot_identity(
    path: Path, label: str, *, run_copy: bool = False
) -> tuple[dict[str, str], Path]:
    fields = _manifest(path, label)
    if fields.get("FORMAT") != "1":
        raise ValueError(f"{label} must use FORMAT=1: {path}")
    source_root = _absolute(fields.get("SOURCE_ROOT", ""), f"{label} SOURCE_ROOT")
    if source_root.name != "tree":
        raise ValueError(f"{label} SOURCE_ROOT must end in /tree: {path}")
    snapshot_root = source_root.parent
    expected_manifest = snapshot_root / "identity" / "source.manifest"
    if not run_copy and path.resolve(strict=False) != expected_manifest.resolve(strict=False):
        raise ValueError(f"{label} is not the snapshot's immutable manifest: {path}")
    _regular_file(snapshot_root / "COMPLETE", f"{label} COMPLETE")
    archive = _absolute(fields.get("SOURCE_ARCHIVE_PATH", ""), f"{label} SOURCE_ARCHIVE_PATH")
    _regular_file(archive, f"{label} archive")
    if fields.get("SOURCE_ARCHIVE_SHA256") != _sha256(archive):
        raise ValueError(f"{label} archive checksum does not match its manifest: {path}")
    if path.read_bytes() != expected_manifest.read_bytes():
        raise ValueError(f"{label} differs from its immutable snapshot copy: {path}")
    return fields, snapshot_root


def _strict_sidecar(path: Path, label: str) -> dict[str, str]:
    sidecar = _regular_file(Path(f"{path}.md5"), f"{label} checksum sidecar")
    lines = sidecar.read_text(encoding="utf-8").splitlines()
    if len(lines) != 3 or [line.split("=", 1)[0] for line in lines] != ["MD5", "SIZE", "PATH"]:
        raise ValueError(f"{label} checksum sidecar schema is invalid: {sidecar}")
    values = {line.split("=", 1)[0]: line.split("=", 1)[1] for line in lines}
    if not re.fullmatch(r"[0-9a-f]{32}", values["MD5"]):
        raise ValueError(f"{label} checksum sidecar MD5 is invalid: {sidecar}")
    if not re.fullmatch(r"[1-9][0-9]*", values["SIZE"]):
        raise ValueError(f"{label} checksum sidecar SIZE is invalid: {sidecar}")
    if values["PATH"] != str(path):
        raise ValueError(f"{label} checksum sidecar PATH is not canonical: {sidecar}")
    actual_size = path.stat().st_size
    actual_md5 = _md5(path)
    if values["SIZE"] != str(actual_size) or values["MD5"] != actual_md5:
        raise ValueError(f"{label} checksum sidecar does not match its file: {path}")
    return values


def _selection(path: Path, run_root: Path, config: dict[str, Any], input_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    _within(path, run_root, "selection")
    _regular_file(path, "selection")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"selection is empty: {path}")
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for line in lines:
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"selection row must have exactly two columns: {line}")
        dataset, view = fields
        if not _COMPONENT_RE.fullmatch(dataset) or view != "batch_effect_corrected":
            raise ValueError(f"selection row is outside corrected Stage 3 scope: {line}")
        key = (dataset, view)
        if key in seen:
            raise ValueError(f"selection contains a duplicate row: {line}")
        seen.add(key)
        entry = config.get(dataset)
        if not isinstance(entry, dict):
            raise ValueError(f"selection dataset is not configured: {dataset}")
        views = entry.get("views")
        view_entry = views.get(view) if isinstance(views, dict) else None
        if not isinstance(view_entry, dict):
            raise ValueError(f"corrected view is not configured: {dataset}/{view}")
        output_name = view_entry.get("output_file_name")
        if not isinstance(output_name, str) or not _COMPONENT_RE.fullmatch(output_name):
            raise ValueError(f"corrected output filename is unsafe: {dataset}/{view}")
        path_value = input_root / dataset / "output" / output_name
        _regular_file(path_value, f"{dataset}/{view} H5AD")
        batch_keys = (entry.get("columns") or {}).get("batch")
        expected_identity = build_batch_contract_identity(
            batch_keys,
            sample_column="Sample",
            method_id="preprocess",
            model_id="hvg_composite_v1",
        )
        validate_benchmark_h5ad_path(
            str(path_value),
            view,
            "Stage 3 preprocessing",
            expected_batch_contract=expected_identity,
            allow_missing_corrected_summary=True,
        )
        sidecar = _strict_sidecar(path_value, f"{dataset}/{view} H5AD")
        rows.append(
            {
                "dataset": dataset,
                "view": view,
                "path": str(path_value),
                "md5": sidecar["MD5"],
                "size": sidecar["SIZE"],
                "contract": _REPORT_CONTRACT,
            }
        )
    return rows, {"path": str(path), "sha256": _sha256(path), "size": path.stat().st_size, "rows": len(rows)}


def _write_atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    os.replace(temporary, path)
    sidecar = Path(f"{path}.md5")
    sidecar_tmp = sidecar.with_name(f".{sidecar.name}.tmp.{os.getpid()}")
    sidecar_tmp.write_text(
        f"MD5={_md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n",
        encoding="utf-8",
    )
    os.replace(sidecar_tmp, sidecar)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--run-source-manifest", required=True)
    parser.add_argument("--validator-source-manifest", required=True)
    parser.add_argument("--runtime-identity", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    run_root = _regular_dir(_absolute(args.run_root, "run root"), "run root")
    run_source_manifest = _absolute(args.run_source_manifest, "run source manifest")
    validator_source_manifest = _absolute(args.validator_source_manifest, "validator source manifest")
    runtime_identity = _absolute(args.runtime_identity, "runtime identity")
    selection = _absolute(args.selection, "selection")
    input_root = _regular_dir(_absolute(args.input_root, "input root"), "input root")
    config_path = _regular_file(_absolute(args.config, "config"), "config")
    output = _absolute(args.output, "output")
    _within(run_source_manifest, run_root, "run source manifest")
    _within(runtime_identity, run_root, "runtime identity")
    _within(selection, run_root, "selection")
    _within(output, run_root, "output")

    run_source_fields, _ = _snapshot_identity(
        run_source_manifest, "run source manifest", run_copy=True
    )
    validator_source_fields, validator_snapshot = _snapshot_identity(
        validator_source_manifest, "validator source manifest"
    )
    runtime_identity_info = _file_identity(runtime_identity, "runtime identity")
    config_info = _file_identity(config_path, "config")
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"config is not valid JSON: {config_path}") from exc
    if not isinstance(config, dict):
        raise ValueError("config must contain a JSON object")
    rows, selection_info = _selection(selection, run_root, config, input_root)
    prior_terminal = _prior_terminal_identity(
        run_root / "status" / "terminal", run_root.name
    )
    payload = {
        "format": 1,
        "stage": "stage3",
        "run_id": run_root.name,
        "run_source_manifest": {
            **_file_identity(run_source_manifest, "run source manifest"),
            "source_commit": run_source_fields["SOURCE_COMMIT"],
        },
        "validator_source_manifest": {
            **_file_identity(validator_source_manifest, "validator source manifest"),
            "source_commit": validator_source_fields["SOURCE_COMMIT"],
            "snapshot_root": str(validator_snapshot),
        },
        "runtime_identity": runtime_identity_info,
        "prior_terminal": prior_terminal,
        "selection": selection_info,
        "config": config_info,
        "rows": rows,
    }
    _write_atomic_json(output, payload)
    print(f"validated Stage 3 sync report: {output} ({len(rows)} rows)")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, TypeError, KeyError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
