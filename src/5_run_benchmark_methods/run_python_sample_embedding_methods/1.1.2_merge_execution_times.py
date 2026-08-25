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
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


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
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
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
                                    output_dir
                                    / f"{args.filename_prefix}{label}_{ds}.feather"
                                )
                            )
                        )
                    )
            else:
                task_logs.extend(
                    sorted(
                        glob.glob(
                            str(output_dir / f"{args.filename_prefix}{label}_*.feather")
                        )
                    )
                )
        task_logs = sorted(set(task_logs))
    else:
        # The merged execution_times.feather has no task suffix and is never
        # matched by a prefix ending in an underscore.
        task_logs = sorted(
            glob.glob(str(output_dir / f"{args.filename_prefix}*.feather"))
        )

    out_path = output_dir / "execution_times.feather"

    if not task_logs:
        print(f"WARNING: No per-task execution logs found in {output_dir} "
              f"for the requested labels/datasets; nothing new to merge.")
        if args.existing_log and os.path.exists(args.existing_log):
            existing = pd.read_feather(args.existing_log)
            existing.reset_index(drop=True).to_feather(out_path)
            print(f"No new rows; wrote existing log unchanged -> {out_path} "
                  f"({len(existing)} rows)")
        return

    frames = [pd.read_feather(f) for f in task_logs]
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    # Merge with the existing log (NAS continuity): this run's rows win
    # (keep="last" after concat), untouched rows from previous runs remain.
    if args.existing_log and os.path.exists(args.existing_log):
        existing = pd.read_feather(args.existing_log)
        merged = pd.concat([existing, merged], ignore_index=True)
        merged = merged.drop_duplicates(subset=["dataset", "method"], keep="last")

    merged = merged.reset_index(drop=True)
    merged.to_feather(out_path)
    print(f"Merged {len(task_logs)} task logs -> {out_path} "
          f"({len(merged)} rows)")

    if args.cleanup:
        for f in task_logs:
            os.remove(f)
        print(f"Deleted {len(task_logs)} per-task logs.")


if __name__ == "__main__":
    main()
