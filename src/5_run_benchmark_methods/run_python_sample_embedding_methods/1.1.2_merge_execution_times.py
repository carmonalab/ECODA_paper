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

def _recorded_checksum_ok(path):
    path = Path(path)
    sidecar = Path(f"{path}.md5")
    if not path.is_file() or not sidecar.is_file():
        return False
    records = {}
    for line in sidecar.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            records[key] = value
    return (
        records.get("PATH") == str(path)
        and records.get("MD5") == _file_md5(path)
        and records.get("SIZE") == str(path.stat().st_size)
    )


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
        if not np.isfinite(values[~missing].to_numpy(dtype=float)).all():
            raise ValueError(f"execution log has invalid numeric values: {path}")

def atomic_feather(frame, path):
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
        sidecar_tmp.write_text(
            f"MD5={_file_md5(path)}\nSIZE={path.stat().st_size}\nPATH={path}\n"
        )
        os.replace(sidecar_tmp, sidecar)
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
    if args.migrate_existing_log:
        migrated = Path(args.migrate_existing_log)
        if not migrated.is_file():
            parser.error(f"existing execution log is missing: {migrated}")
        frame = pd.read_feather(migrated)
        _validate_log_frame(frame, migrated)
        atomic_feather(frame, migrated)
        print(f"Execution log sidecar migrated: {migrated}")
        return

    output_dir = Path(args.output_dir)
    log_dir = Path(args.log_dir or args.output_dir)
    if args.datasets is not None and args.labels is None:
        parser.error("--datasets requires --labels")
    if args.labels is not None:
        task_logs = []
        for label in args.labels:
            if args.datasets is not None:
                for ds in args.datasets:
                    task_logs.extend(
                        sorted(
                            glob.glob(
                                str(
                                    log_dir
                                    / f"{args.filename_prefix}{label}_{ds}.feather"
                                )
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
        if args.existing_log and os.path.exists(args.existing_log):
            existing = pd.read_feather(args.existing_log)
            _validate_log_frame(existing, args.existing_log)
            atomic_feather(existing, out_path)
            print(
                f"No new rows; wrote existing log unchanged -> {out_path} "
                f"({len(existing)} rows)"
            )
        return

    frames = []
    for task_log in task_logs:
        if not _recorded_checksum_ok(task_log):
            raise ValueError(f"execution log checksum failed: {task_log}")
        frame = pd.read_feather(task_log)
        _validate_log_frame(frame, task_log)
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    # Merge with the existing log (NAS continuity): this run's rows win.
    if args.existing_log and os.path.exists(args.existing_log):
        existing = pd.read_feather(args.existing_log)
        _validate_log_frame(existing, args.existing_log)
        merged = pd.concat([existing, merged], ignore_index=True)
        merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    merged = merged.reset_index(drop=True)
    _validate_log_frame(merged, out_path)
    atomic_feather(merged, out_path)
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
