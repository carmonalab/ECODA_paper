"""Merge per-task benchmark execution-time logs into one feather.

Concatenates the per-task logs matching the given (label x dataset) cross
product (labels are the benchmark method names or 'trans'/'zeroimp'
analyses; each log file is `execution_times_<label>_<ds>.feather`) from the
benchmark embeddings output dir into `execution_times.feather`, deduplicating
on (dataset, method) with the last occurrence kept (matches the qmd's
overwrite-on-rerun semantics). Runs on the login node after the benchmark
arrays complete.

Scoping to the run's label x dataset cross product keeps stale logs from
previous failed runs out of the merge; `--existing-log` preserves the NAS log
across partial runs (e.g. `--ds_name _debug`), so a subset run extends the
full log instead of overwriting it. Per-task log deletion is `--cleanup`
(default on); the submit script passes `--no-cleanup` and deletes the logs
itself only after the NAS rsync has succeeded.

Usage:
    python 1.1.2_merge_execution_times.py [--output_dir <dir>]
        [--labels <name>... --datasets <ds>...] [--existing-log <path>]
        [--no-cleanup]
"""

import argparse
import glob
import hashlib
import os
import re
import sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


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
    records = _read_checksum_sidecar(path)
    digest = _file_md5(path)
    if digest != records["MD5"]:
        raise ValueError(f"checksum sidecar does not match artifact: {path}")
    if path.stat().st_size != int(records["SIZE"]):
        raise ValueError(f"artifact changed during checksum validation: {path}")
    return records


def _record_context(run_id=None):
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
    record = _read_artifact_record(path, producer, run_id, require=True)
    if record is None:
        raise ValueError(f"artifact record is missing: {path}")
    return record


def _write_artifact_record(path, producer, run_id=None, checksum=None):
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

def _current_run_root(run_id=None):
    """Return the canonical producer root for the current run, if bound."""
    configured_root = os.environ.get("ECODA_RUN_ROOT", "")
    if configured_root:
        try:
            return Path(configured_root).resolve()
        except (OSError, RuntimeError):
            return None
    context = _record_context(run_id)
    if context is None:
        return None
    return (context[0] / context[1]).resolve()


def _path_is_current_run_owned(path, run_id=None):
    root = _current_run_root(run_id)
    if root is None:
        return False
    try:
        Path(path).resolve().relative_to(root)
    except (OSError, RuntimeError, ValueError):
        return False
    return True


def _recorded_checksum_ok(
    path, producer=None, run_id=None, *, require_record=False
):
    """Require a strict checksum and any record required by path ownership."""
    path = Path(path)
    record_required = require_record or _path_is_current_run_owned(path, run_id)
    try:
        _full_checksum(path)
        if record_required:
            _read_artifact_record(path, producer, run_id, require=True)
    except (OSError, ValueError):
        return False
    return True


def _validate_log_frame(frame, path):
    required = {"dataset", "method", "time_secs", "mem_GB"}
    if set(frame.columns) != required or frame.empty:
        raise ValueError(
            f"execution log has an invalid schema or no rows: {path}"
        )
    if frame[["dataset", "method"]].isna().any().any() or \
       frame[["dataset", "method"]].astype(str).apply(lambda column: column.str.strip() == "").any().any():
        raise ValueError(f"execution log has blank identifiers: {path}")
    if frame[["dataset", "method"]].duplicated().any():
        raise ValueError(f"execution log has duplicate identifiers: {path}")
    for column, allow_missing in (("time_secs", False), ("mem_GB", True)):
        raw = frame[column]
        values = pd.to_numeric(raw, errors="coerce")
        missing = raw.isna()
        if (not allow_missing and missing.any()) or values[~missing].isna().any():
            raise ValueError(f"execution log has invalid numeric values: {path}")
        if (not np.isfinite(values[~missing].to_numpy(dtype=float)).all() or (values[~missing] < 0).any()):
            raise ValueError(f"execution log has invalid numeric values: {path}")

def _read_recorded_execution_log(
    path, producer=None, run_id=None, *, require_record=False
):
    """Validate a log's checksum, record (when required), and schema."""
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    present = path.exists() or sidecar.exists()
    if not present:
        return None
    if not path.is_file() or not sidecar.is_file():
        raise ValueError(f"execution log is incomplete: {path}")
    record_required = require_record or _path_is_current_run_owned(path, run_id)
    if not _recorded_checksum_ok(
        path,
        producer=producer,
        run_id=run_id,
        require_record=record_required,
    ):
        if record_required:
            candidate = _record_candidate(path, run_id)
            if candidate is None or (
                not candidate.exists() and not candidate.is_symlink()
            ):
                raise ValueError(f"execution log artifact record is missing: {path}")
        raise ValueError(f"execution log checksum or record failed: {path}")
    try:
        frame = pd.read_feather(path)
    except Exception as exc:
        raise ValueError(f"execution log is malformed: {path}") from exc
    _validate_log_frame(frame, path)
    return frame


def atomic_feather(frame, path, producer=None, *, write_record=True):
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
        _validate_log_frame(frame, path)
        frame.reset_index(drop=True).to_feather(tmp)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise RuntimeError(f"empty merged execution log: {tmp}")
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
        if write_record and _record_context() is not None:
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

def _atomic_text(lines, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        tmp.write_text("".join(lines))
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()

def main():
    parser = argparse.ArgumentParser(
        description="Merge per-task benchmark exec-time logs into "
                    "execution_times.feather."
    )
    parser.add_argument(
        "--output_dir",
        default="benchmark/embeddings",
        help="Directory with per-task logs and output feather "
             "(default: benchmark/embeddings)",
    )
    parser.add_argument(
        "--log-dir",
        default=None,
        help="Run-owned directory containing per-task logs "
             "(default: --output_dir)",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help="Method/analysis labels whose per-task logs "
             "(execution_times_<label>_<ds>.feather) to merge "
             "(default: all task logs)",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        help="Dataset names to merge logs for, crossed with --labels "
             "(default: all datasets)",
    )
    parser.add_argument(
        "--existing-log",
        default=None,
        help="Existing execution_times.feather (e.g. the NAS copy) to merge "
             "with: its rows are kept unless overridden by a row of this run",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        default=True,
        help="Delete merged per-task logs (default: on)",
    )
    parser.add_argument(
        "--no-cleanup",
        action="store_false",
        dest="cleanup",
        help="Keep per-task logs after merging",
    )
    parser.add_argument(
        "--filename_prefix",
        default="execution_times_",
        help="Prefix before <label>_<dataset> in per-task log filenames "
             "(default: execution_times_)",
    )
    parser.add_argument(
        "--migrate-existing-log",
        default=None,
        help="Validate and sidecar-migrate an existing indexless execution log in place",
    )
    parser.add_argument(
        "--cleanup-manifest",
        default=None,
        help="Write the exact task-log paths eligible for cleanup",
    )
    args = parser.parse_args()
    log_producer = os.environ.get(
        "ECODA_EXECUTION_LOG_PRODUCER", "stage5_execution_log"
    )
    if args.migrate_existing_log:
        migrated = Path(args.migrate_existing_log)
        frame = _read_recorded_execution_log(
            migrated, producer=log_producer
        )
        if frame is None:
            parser.error(f"existing execution log is missing: {migrated}")
        atomic_feather(
            frame,
            migrated,
            producer=log_producer,
            write_record=_path_is_current_run_owned(migrated),
        )
        print(f"Execution log sidecar migrated: {migrated}")
        return

    output_dir = Path(args.output_dir)
    log_dir = Path(args.log_dir or args.output_dir)
    task_logs = []
    if args.datasets is not None and args.labels is None:
        parser.error("--datasets requires --labels")
    if args.labels is not None:
        for label in args.labels:
            if args.datasets is not None:
                for ds in args.datasets:
                    base = log_dir / f"{args.filename_prefix}{label}_{ds}.feather"
                    shard_pattern = (
                        log_dir / f"{args.filename_prefix}{label}_{ds}_*.feather"
                    )
                    task_logs.extend(
                        sorted(
                            set(
                                glob.glob(str(base))
                                + glob.glob(str(shard_pattern))
                            )
                        )
                    )
            else:
                task_logs.extend(
                    sorted(
                        glob.glob(
                            str(log_dir / f"{args.filename_prefix}{label}_*.feather")
                        )
                    )
                )
        task_logs = sorted(set(task_logs))
    else:
        # The merged execution_times.feather has no task suffix and is never
        # matched by a prefix ending in an underscore.
        task_logs = sorted(
            glob.glob(str(log_dir / f"{args.filename_prefix}*.feather"))
        )

    out_path = output_dir / "execution_times.feather"

    if not task_logs:
        print(
            f"WARNING: No per-task execution logs found in {log_dir} "
            f"for the requested labels/datasets; nothing new to merge."
        )
        if args.cleanup_manifest:
            _atomic_text([], args.cleanup_manifest)
        if not args.existing_log:
            raise ValueError(
                "no per-task execution logs or validated existing log were found"
            )
        existing = _read_recorded_execution_log(
            args.existing_log, producer=log_producer
        )
        if existing is None:
            raise ValueError(f"existing execution log is missing: {args.existing_log}")
        atomic_feather(existing, out_path, producer=log_producer)
        print(
            f"No new rows; wrote existing log unchanged -> {out_path} "
            f"({len(existing)} rows)"
        )
        return

    frames = []
    for task_log in task_logs:
        frame = _read_recorded_execution_log(
            task_log,
            producer=log_producer,
        )
        if frame is None:
            raise ValueError(f"execution log is missing: {task_log}")
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    # Merge with the existing log (NAS continuity): this run's rows win.
    if args.existing_log:
        existing = _read_recorded_execution_log(
            args.existing_log, producer=log_producer
        )
        if existing is not None:
            merged = pd.concat([existing, merged], ignore_index=True)
            merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    merged = merged.reset_index(drop=True)
    _validate_log_frame(merged, out_path)
    atomic_feather(merged, out_path, producer=log_producer)
    if args.cleanup_manifest:
        _atomic_text(
            [f"{task_log}\n" for task_log in task_logs],
            args.cleanup_manifest,
        )
    print(f"Merged {len(task_logs)} task logs -> {out_path} "
          f"({len(merged)} rows)")

    if args.cleanup:
        for task_log in task_logs:
            os.remove(task_log)
            sidecar = Path(f"{task_log}.md5")
            if sidecar.exists():
                sidecar.unlink()
        print(f"Deleted {len(task_logs)} per-task logs.")


if __name__ == "__main__":
    main()
