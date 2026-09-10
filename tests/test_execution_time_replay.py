#!/usr/bin/env python3
"""Focused runtime-sidecar replay and execution-log validation checks."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PYTHON_WORKER = (
    ROOT
    / "src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py"
)
MERGE_WORKER = (
    ROOT
    / "src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py"
)
R_WORKER = ROOT / "src/5_run_benchmark_methods/benchmark_hpc_utils.R"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_checksum(path: Path) -> None:
    path.with_name(f"{path.name}.md5").write_text(
        f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\n"
        f"SIZE={path.stat().st_size}\n"
        f"PATH={path}\n"
    )


def write_feather(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_feather(path)
    write_checksum(path)


def output_frame(value: float = 1.0) -> pd.DataFrame:
    return pd.DataFrame(
        {"s1": [value, 0.0], "s2": [0.0, value]},
        index=["s1", "s2"],
    )


def expect_value_error(callable_) -> None:
    try:
        callable_()
    except ValueError:
        return
    raise AssertionError("expected a fail-closed ValueError")


def python_args(root: Path, log_file: Path) -> SimpleNamespace:
    input_dir = root / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    # process_dataset checks that the configured input path exists before it
    # evaluates cache state, but a valid cache must not load this placeholder.
    (input_dir / "view.h5ad").write_bytes(b"not an h5ad; cache path must skip it")
    return SimpleNamespace(
        view="benchmark_analysis",
        analysis_pass=None,
        combo=None,
        high_resolution_only=True,
        output_dir=str(root / "embeddings"),
        input_dir=str(input_dir),
        hvg=[2000],
        method="pilotgm",
        force=False,
        log_file=str(log_file),
        device="cpu",
    )


def python_entry() -> dict:
    return {
        "views": {"benchmark_analysis": {"output_file": "view.h5ad"}},
        "cell_type_high_res": "cell_type_high",
    }


def _check_python_cache_hit_replays_without_work(worker, root: Path) -> None:
    output = root / "embeddings" / "Adams_hvg2000_highres_pilotgm_dists.feather"
    log_file = root / "logs" / "execution_times.feather"
    write_feather(output, output_frame())
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=12.345678,
        mem_gb=2.75,
    )

    args = python_args(root, log_file)
    calls: list[str] = []

    def forbidden(*_args, **_kwargs):
        calls.append("work")
        raise AssertionError("a valid cache hit invoked method or input work")

    with (
        patch.object(worker, "run_pilotgm", forbidden),
        patch.object(worker, "load_h5ad_counts_free", forbidden),
        patch.object(worker, "validate_benchmark_h5ad_path", forbidden),
    ):
        worker.process_dataset(args, "Adams", python_entry())

    assert calls == []
    replayed = pd.read_feather(log_file)
    assert list(replayed.columns) == ["dataset", "method", "time_secs", "mem_GB"]
    assert len(replayed) == 1
    row = replayed.iloc[0]
    assert row["dataset"] == "Adams"
    assert row["method"] == "PILOT-GM-VAE_hvg2000_highres"
    assert row["time_secs"] == 12.345678
    assert row["mem_GB"] == 2.75

def _check_python_cache_hit_uses_default_log_file(worker, root: Path) -> None:
    output = root / "embeddings" / "Adams_hvg2000_highres_pilotgm_dists.feather"
    write_feather(output, output_frame())
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=21.987654,
        mem_gb=4.125,
    )

    args = python_args(root, root / "unused" / "execution_times.feather")
    args.log_file = None
    calls: list[str] = []

    def forbidden(*_args, **_kwargs):
        calls.append("work")
        raise AssertionError("a valid cache hit invoked method or input work")

    with (
        patch.object(worker, "run_pilotgm", forbidden),
        patch.object(worker, "load_h5ad_counts_free", forbidden),
        patch.object(worker, "validate_benchmark_h5ad_path", forbidden),
    ):
        worker.process_dataset(args, "Adams", python_entry())

    assert calls == []
    default_log = Path(args.output_dir) / "execution_times.feather"
    assert default_log.is_file()
    replayed = pd.read_feather(default_log)
    assert len(replayed) == 1
    row = replayed.iloc[0]
    assert row["dataset"] == "Adams"
    assert row["method"] == "PILOT-GM-VAE_hvg2000_highres"
    assert row["time_secs"] == 21.987654
    assert row["mem_GB"] == 4.125


def _check_python_nullable_memory_replays_as_missing(worker, root: Path) -> None:
    output = root / "embeddings" / "Adams_hvg2000_highres_pilotgm_dists.feather"
    log_file = root / "logs" / "execution_times.feather"
    write_feather(output, output_frame())
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=8.5,
        mem_gb=None,
    )
    metadata = json.loads(Path(f"{output}.runtime.json").read_text())
    assert metadata["mem_GB"] is None
    assert worker.read_runtime_metadata(
        output, "Adams", "PILOT-GM-VAE_hvg2000_highres"
    )["mem_GB"] is None

    args = python_args(root, log_file)
    calls: list[str] = []

    def forbidden(*_args, **_kwargs):
        calls.append("work")
        raise AssertionError("a valid cache hit invoked method or input work")

    with (
        patch.object(worker, "run_pilotgm", forbidden),
        patch.object(worker, "load_h5ad_counts_free", forbidden),
        patch.object(worker, "validate_benchmark_h5ad_path", forbidden),
    ):
        worker.process_dataset(args, "Adams", python_entry())

    assert calls == []
    replayed = pd.read_feather(log_file)
    assert len(replayed) == 1
    row = replayed.iloc[0]
    assert row["time_secs"] == 8.5
    assert pd.isna(row["mem_GB"])



def _check_python_metadata_rejection(worker, root: Path) -> None:
    output = root / "embeddings" / "Adams_hvg2000_highres_pilotgm_dists.feather"
    write_feather(output, output_frame())
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=3.0,
        mem_gb=1.5,
    )
    metadata = Path(f"{output}.runtime.json")
    metadata_checksum = Path(f"{metadata}.md5")

    # Missing metadata is an error, not a cache miss that silently drops the
    # execution-time row.
    metadata.unlink()
    metadata_checksum.unlink()
    expect_value_error(
        lambda: worker.read_runtime_metadata(
            output, "Adams", "PILOT-GM-VAE_hvg2000_highres"
        )
    )

    # A metadata checksum must cover the metadata JSON itself.
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=3.0,
        mem_gb=1.5,
    )
    payload = json.loads(metadata.read_text())
    payload["time_secs"] = 4.0
    metadata.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    expect_value_error(
        lambda: worker.read_runtime_metadata(
            output, "Adams", "PILOT-GM-VAE_hvg2000_highres"
        )
    )
    # A semantically wrong label must fail even when the metadata JSON has a
    # fresh, valid checksum.
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=3.0,
        mem_gb=1.5,
    )
    payload = json.loads(metadata.read_text())
    payload["method"] = "PILOT-GM-VAE_hvg1000_highres"
    metadata.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    write_checksum(metadata)
    expect_value_error(
        lambda: worker.read_runtime_metadata(
            output, "Adams", "PILOT-GM-VAE_hvg2000_highres"
        )
    )


    # Keep the output checksum valid but change its bytes.  The metadata still
    # names the old artifact MD5 and must therefore be rejected.
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=3.0,
        mem_gb=1.5,
    )
    write_feather(output, output_frame(2.0))
    expect_value_error(
        lambda: worker.read_runtime_metadata(
            output, "Adams", "PILOT-GM-VAE_hvg2000_highres"
        )
    )

    # The process-level cache branch also fails before an input loader or
    # method can run when metadata is absent.
    worker.publish_runtime_metadata(
        output,
        "Adams",
        "PILOT-GM-VAE_hvg2000_highres",
        time_secs=3.0,
        mem_gb=1.5,
    )
    metadata.unlink()
    metadata_checksum.unlink()
    args = python_args(root, root / "logs" / "execution_times.feather")

    def forbidden(*_args, **_kwargs):
        raise AssertionError("invalid cache metadata fell through to computation")

    with (
        patch.object(worker, "run_pilotgm", forbidden),
        patch.object(worker, "load_h5ad_counts_free", forbidden),
        patch.object(worker, "validate_benchmark_h5ad_path", forbidden),
    ):
        expect_value_error(lambda: worker.process_dataset(args, "Adams", python_entry()))


def _check_merge_rejects_duplicate_and_invalid_runtime_rows(merge_worker) -> None:
    valid = pd.DataFrame(
        {
            "dataset": ["Adams"],
            "method": ["PILOT-GM-VAE_hvg2000_highres"],
            "time_secs": [12.0],
            "mem_GB": [2.0],
        }
    )
    shared = pd.DataFrame(
        {
            "dataset": ["Adams"],
            "method": ["prepare_pseudobulk_shared"],
            "time_secs": [4.0],
            "mem_GB": [np.nan],
        }
    )
    # Legacy rows and the one shared preparation row coexist in the canonical
    # four-column log; a second shared row is still a duplicate identifier.
    merge_worker._validate_log_frame(pd.concat([valid, shared]), "mixed")
    duplicate = pd.concat([valid, valid], ignore_index=True)
    expect_value_error(lambda: merge_worker._validate_log_frame(duplicate, "duplicate"))
    duplicate_shared = pd.concat([shared, shared], ignore_index=True)
    expect_value_error(
        lambda: merge_worker._validate_log_frame(duplicate_shared, "duplicate-shared")
    )

    for column, value in (
        ("time_secs", np.nan),
        ("time_secs", np.inf),
        ("time_secs", -1.0),
        ("mem_GB", np.inf),
        ("mem_GB", -1.0),
        ("time_secs", "not-a-number"),
    ):
        invalid = valid.copy()
        invalid.loc[0, column] = value
        expect_value_error(lambda invalid=invalid: merge_worker._validate_log_frame(invalid, "invalid"))
    dynamic_shared = pd.DataFrame(
        {
            "dataset": ["Adams", "Adams", "Adams"],
            "method": [
                "prepare_pseudobulk_ct_shared_LR",
                "prepare_pseudobulk_ct_shared_HR",
                "prepare_pseudobulk_ct_shared_LR",
            ],
            "time_secs": [4.0, 5.0, 4.0],
            "mem_GB": [np.nan, np.nan, np.nan],
            "shared_time_secs": [4.0, 5.0, 4.0],
            "timing_id": ["ct-lr", "ct-hr", "ct-lr"],
            "timing_schema": [2, 2, 2],
        }
    )
    # Per-task shards validate independently, then deduplication keeps both
    # dynamic CT shared methods and removes only the repeated LR identity.
    merge_worker._validate_log_frame(dynamic_shared.iloc[:2], "dynamic-shared")
    merge_worker._validate_log_frame(
        dynamic_shared.iloc[[2]], "dynamic-shared-repeat"
    )
    dynamic_dedup = merge_worker._deduplicate_log_frame(
        pd.concat(
            [dynamic_shared.iloc[:2], dynamic_shared.iloc[[2]]],
            ignore_index=True,
        )
    )
    assert len(dynamic_dedup) == 2
    assert set(dynamic_dedup["method"]) == {
        "prepare_pseudobulk_ct_shared_LR",
        "prepare_pseudobulk_ct_shared_HR",
    }
    assert set(dynamic_dedup["time_secs"]) == {4.0, 5.0}

    malformed_schema = valid.copy()
    malformed_schema["timing_schema"] = ["not-a-schema"]
    expect_value_error(
        lambda: merge_worker._validate_log_frame(
            malformed_schema, "malformed-schema"
        )
    )

    local_missing_variant = dynamic_shared.iloc[[0]].copy()
    local_missing_variant.loc[:, "method"] = "Pseudobulk_hvg500"
    expect_value_error(
        lambda: merge_worker._validate_log_frame(
            local_missing_variant, "missing-variant"
        )
    )

    shared_base_mismatch = dynamic_shared.iloc[[0]].copy()
    shared_base_mismatch.loc[:, "time_secs"] = 3.0
    expect_value_error(
        lambda: merge_worker._validate_log_frame(
            shared_base_mismatch, "shared-base-mismatch"
        )
    )

    inconsistent_timing = pd.DataFrame(
        {
            "dataset": ["Adams", "Adams"],
            "method": [
                "prepare_pseudobulk_ct_shared_LR",
                "Pseudobulk_hvg500",
            ],
            "time_secs": [4.0, 0.5],
            "mem_GB": [np.nan, np.nan],
            "aggregate_time_secs": [1.25, 1.5],
            "shared_fit_time_secs": [2.75, 2.5],
            "shared_time_secs": [4.0, 4.0],
            "variant_time_secs": [np.nan, 0.5],
            "shared_mem_GB": [np.nan, np.nan],
            "timing_id": ["ct-id", "ct-id"],
            "timing_schema": [2, 2],
        }
    )
    expect_value_error(
        lambda: merge_worker._validate_log_frame(
            inconsistent_timing, "inconsistent-timing"
        )
    )


def run_r_replay_fixture(root: Path) -> pd.DataFrame:
    script = root / "r_runtime_replay_fixture.R"
    script.write_text(
        r'''args <- commandArgs(trailingOnly = TRUE)
helper <- normalizePath(args[[1]], mustWork = TRUE)
root <- normalizePath(args[[2]], mustWork = FALSE)
source(helper)
helpers <- normalizePath(
  file.path(dirname(helper), "..", "utils", "helpers.R"),
  mustWork = TRUE
)
source(helpers)
results <- file.path(root, "results")
dir.create(results, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(root, "execution_times.feather")

specs <- list(
  list(suffix = "trans", method = "trans_analysis", time = 17.25, mem = 2.75),
  list(suffix = "zeroimp", method = "zeroimp_analysis", time = 23.5, mem = 3.125)
)
for (spec in specs) {
  output <- file.path(results, paste0("Adams_", spec$suffix, ".rds"))
  save_rds_atomic(list(result = 1L), output)
  write_runtime_metadata(
    output, "Adams", spec$method, spec$time, spec$mem
  )
  validated <- validate_runtime_metadata(output, "Adams", spec$method)
  stopifnot(is.list(validated))
  stopifnot(isTRUE(all.equal(as.numeric(validated$time_secs), spec$time)))
  stopifnot(isTRUE(all.equal(as.numeric(validated$mem_GB), spec$mem)))
  replay_runtime_metadata(output, "Adams", spec$method, log_file)
}

rows <- arrow::read_feather(log_file)
stopifnot(nrow(rows) == 2L)
for (spec in specs) {
  row <- rows[rows$method == spec$method, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  stopifnot(identical(as.character(row$dataset[[1L]]), "Adams"))
  stopifnot(identical(as.character(row$method[[1L]]), spec$method))
  stopifnot(isTRUE(all.equal(as.numeric(row$time_secs[[1L]]), spec$time)))
  stopifnot(isTRUE(all.equal(as.numeric(row$mem_GB[[1L]]), spec$mem)))
}

# Schema-2 pseudobulk cache timing is replayed as one shared row plus one
# variant-local row per cache.  Keep this in a separate log so the legacy
# runtime fixture above remains a four-column two-row log.
schema2_variants <- list(
  hvg500 = list(
    pb = matrix(1, nrow = 1L, ncol = 1L,
                dimnames = list("s1", "g1")),
    time_secs = 0.25,
    mem_GB = NA_real_,
    aggregate_time_secs = 1.25,
    shared_fit_time_secs = 2.75,
    shared_time_secs = 4,
    variant_time_secs = 0.25,
    shared_mem_GB = NA_real_,
    timing_id = "run-1:Adams:benchmark_analysis:none",
    timing_schema = 2L
  ),
  hvg2000 = list(
    pb = matrix(2, nrow = 1L, ncol = 1L,
                dimnames = list("s1", "g1")),
    time_secs = 0.75,
    mem_GB = NA_real_,
    aggregate_time_secs = 1.25,
    shared_fit_time_secs = 2.75,
    shared_time_secs = 4,
    variant_time_secs = 0.75,
    shared_mem_GB = NA_real_,
    timing_id = "run-1:Adams:benchmark_analysis:none",
    timing_schema = 2L
  )
)
schema2_fields <- c(
  "pb", "time_secs", "mem_GB", "aggregate_time_secs",
  "shared_fit_time_secs", "shared_time_secs", "variant_time_secs",
  "shared_mem_GB", "timing_id", "timing_schema"
)
stopifnot(identical(sort(names(schema2_variants[[1L]])), sort(schema2_fields)))
shared_log_file <- file.path(root, "pseudobulk-execution-times.feather")
emit_pseudobulk_timing_rows(
  schema2_variants,
  ds = "Adams",
  log_file = shared_log_file
)
shared_rows <- arrow::read_feather(shared_log_file)
stopifnot(
  nrow(shared_rows) == 3L,
  sum(shared_rows$method == "prepare_pseudobulk_shared") == 1L,
  sum(shared_rows$time_secs) == 4 + 0.25 + 0.75,
  all(c("dataset", "method", "time_secs", "mem_GB") %in%
      colnames(shared_rows))
)
emit_pseudobulk_timing_rows(
  schema2_variants,
  ds = "Adams",
  log_file = shared_log_file
)
shared_rows_replayed <- arrow::read_feather(shared_log_file)
stopifnot(
  nrow(shared_rows_replayed) == 3L,
  sum(shared_rows_replayed$method == "prepare_pseudobulk_shared") == 1L,
  sum(shared_rows_replayed$time_secs) == 4 + 0.25 + 0.75
)
bad_timing_ids <- schema2_variants
bad_timing_ids$hvg2000$timing_id <- "run-2:Adams:benchmark_analysis:none"
bad <- try(
  emit_pseudobulk_timing_rows(
    bad_timing_ids,
    ds = "Adams",
    log_file = file.path(root, "invalid-timing-ids.feather")
  ),
  silent = TRUE
)
stopifnot(inherits(bad, "try-error"))
bad_timing_ids$hvg2000$timing_id <- ""
bad <- try(
  emit_pseudobulk_timing_rows(
    bad_timing_ids,
    ds = "Adams",
    log_file = file.path(root, "blank-timing-id.feather")
  ),
  silent = TRUE
)
stopifnot(inherits(bad, "try-error"))

# A legacy cache payload remains inclusive and does not fabricate shared time.
legacy_timing_log <- file.path(root, "legacy-execution-times.feather")
log_exec_row("Adams", "prepare_pseudobulk_hvg500", 5.5, legacy_timing_log)
legacy_rows <- arrow::read_feather(legacy_timing_log)
stopifnot(nrow(legacy_rows) == 1L, legacy_rows$time_secs[[1L]] == 5.5)

# Helper timing validation is fail-closed for vectors and malformed schema
# markers rather than silently coercing them to legacy rows.
expect_helper_error <- function(fn) {
  value <- try(fn(), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
}
expect_helper_error(function() .execution_scalar(c(1, 2), "vector"))
malformed_schema <- data.frame(
  dataset = c("Adams", "Adams"),
  method = c("legacy_a", "legacy_b"),
  time_secs = c(1, 2),
  mem_GB = c(NA_real_, NA_real_),
  timing_schema = c("2", "not-a-schema"),
  stringsAsFactors = FALSE
)
expect_helper_error(function() normalize_exec_times(malformed_schema))

# A timing_id without timing_schema remains one legacy-inclusive row.
legacy_bundle <- execution_time_rows_from_bundle(
  "Adams",
  "legacy_method",
  list(exec_time = 6, variant_time_secs = 2, timing_id = "id-without-schema")
)
stopifnot(
  nrow(legacy_bundle) == 1L,
  legacy_bundle$time_secs[[1L]] == 6,
  !"timing_schema" %in% names(legacy_bundle)
)

ct_bundle_rows <- execution_time_rows_from_bundle(
  "Adams",
  "Pseudobulk_CT_HR_hvg2000",
  schema2_variants$hvg500
)
stopifnot(
  sum(ct_bundle_rows$method == "prepare_pseudobulk_ct_shared") == 1L,
  sum(ct_bundle_rows$method == "Pseudobulk_CT_HR_hvg2000") == 1L
)

ct_shared_log <- data.frame(
  dataset = rep("Adams", 3L),
  method = c(
    "prepare_pseudobulk_ct_shared_LR",
    "prepare_pseudobulk_ct_shared_HR",
    "prepare_pseudobulk_ct_shared_LR"
  ),
  time_secs = c(4, 5, 4),
  mem_GB = c(NA_real_, NA_real_, NA_real_),
  shared_time_secs = c(4, 5, 4),
  shared_mem_GB = c(NA_real_, NA_real_, NA_real_),
  timing_id = c("ct-lr", "ct-hr", "ct-lr"),
  timing_schema = c(2L, 2L, 2L),
  stringsAsFactors = FALSE
)
ct_shared_dedup <- deduplicate_exec_times(ct_shared_log)
stopifnot(
  nrow(ct_shared_dedup) == 2L,
  setequal(
    ct_shared_dedup$method,
    c(
      "prepare_pseudobulk_ct_shared_LR",
      "prepare_pseudobulk_ct_shared_HR"
    )
  ),
  setequal(ct_shared_dedup$time_secs, c(4, 5))
)
dynamic_ct_bundle <- schema2_variants$hvg500
dynamic_ct_bundle$shared_timing_method <-
  "prepare_pseudobulk_ct_shared_LR"
dynamic_ct_rows <- execution_time_rows_from_bundle(
  "Adams",
  "Pseudobulk_CT_LR_hvg2000",
  dynamic_ct_bundle
)
stopifnot(
  sum(
    dynamic_ct_rows$method ==
      "prepare_pseudobulk_ct_shared_LR"
  ) == 1L
)

# The persisted execution log remains the canonical four-column format while
# retaining separate LR/HR CT-column shared rows.
ct_four_column_log <- file.path(root, "ct-shared-execution-times.feather")
log_exec_row(
  "Adams", "prepare_pseudobulk_ct_shared_LR", 4, ct_four_column_log
)
log_exec_row(
  "Adams", "prepare_pseudobulk_ct_shared_HR", 5, ct_four_column_log
)
ct_four_column_rows <- arrow::read_feather(ct_four_column_log)
stopifnot(
  identical(
    colnames(ct_four_column_rows),
    c("dataset", "method", "time_secs", "mem_GB")
  ),
  nrow(ct_four_column_rows) == 2L,
  setequal(
    ct_four_column_rows$method,
    c(
      "prepare_pseudobulk_ct_shared_LR",
      "prepare_pseudobulk_ct_shared_HR"
    )
  )
)

# Legacy preparation rows retain their inclusive time, and a dataset mixing
# schema-2 and legacy rows classifies each row by its own marker.
legacy_summary <- summarize_exec_times(data.frame(
  dataset = c("Adams", "Adams"),
  method = c("prepare_pseudobulk_hvg500", "trans_analysis"),
  time_secs = c(5.5, 2),
  mem_GB = c(NA_real_, 1),
  stringsAsFactors = FALSE
))
legacy_row <- legacy_summary[legacy_summary$dataset == "Adams", , drop = FALSE]
stopifnot(
  nrow(legacy_row) == 1L,
  legacy_row$variant_time_secs[[1L]] == 0,
  legacy_row$shared_time_secs[[1L]] == 0,
  legacy_row$legacy_inclusive_time_secs[[1L]] == 7.5,
  legacy_row$total_time_secs[[1L]] == 7.5
)

mixed_rows <- data.frame(
  dataset = rep("Adams", 4L),
  method = c(
    "legacy_method", "Pseudobulk_hvg500",
    "prepare_pseudobulk_hvg500", "prepare_pseudobulk_shared"
  ),
  time_secs = c(3, 0.75, 2, 4),
  mem_GB = c(NA_real_, NA_real_, NA_real_, NA_real_),
  aggregate_time_secs = c(NA_real_, 1.25, NA_real_, 1.25),
  shared_fit_time_secs = c(NA_real_, 2.75, NA_real_, 2.75),
  shared_time_secs = c(NA_real_, 4, NA_real_, 4),
  variant_time_secs = c(NA_real_, 0.75, NA_real_, NA_real_),
  shared_mem_GB = c(NA_real_, NA_real_, NA_real_, NA_real_),
  timing_id = c(NA_character_, "id-a", NA_character_, "id-a"),
  timing_schema = c(NA_integer_, 2L, NA_integer_, 2L),
  stringsAsFactors = FALSE
)
mixed_summary <- summarize_exec_times(mixed_rows)
mixed_row <- mixed_summary[mixed_summary$dataset == "Adams", , drop = FALSE]
stopifnot(
  mixed_row$variant_time_secs[[1L]] == 0.75,
  mixed_row$shared_time_secs[[1L]] == 4,
  mixed_row$legacy_inclusive_time_secs[[1L]] == 5,
  mixed_row$total_time_secs[[1L]] == 9.75
)

# Missing shared rows are derived independently for distinct timing IDs.
derived_rows <- data.frame(
  dataset = c("Adams", "Adams"),
  method = c("Pseudobulk_hvg500", "Pseudobulk_hvg2000"),
  time_secs = c(0.5, 0.75),
  mem_GB = c(NA_real_, NA_real_),
  aggregate_time_secs = c(1.25, 2.25),
  shared_fit_time_secs = c(2.75, 3.75),
  shared_time_secs = c(4, 6),
  variant_time_secs = c(0.5, 0.75),
  shared_mem_GB = c(NA_real_, NA_real_),
  timing_id = c("id-a", "id-b"),
  timing_schema = c(2L, 2L),
  stringsAsFactors = FALSE
)
derived_summary <- summarize_exec_times(derived_rows)
derived_row <- derived_summary[derived_summary$dataset == "Adams", , drop = FALSE]
stopifnot(
  derived_row$variant_time_secs[[1L]] == 1.25,
  derived_row$shared_time_secs[[1L]] == 10,
  derived_row$total_time_secs[[1L]] == 11.25
)
bad_decomposition <- derived_rows
bad_decomposition$shared_time_secs[[1L]] <- 5
expect_helper_error(function() summarize_exec_times(bad_decomposition))
bad_shared_row <- mixed_rows
bad_shared_row$variant_time_secs[[4L]] <- 0
expect_helper_error(function() summarize_exec_times(bad_shared_row))

# Keep the canonical runtime log assertions independent from the schema-2
# pseudobulk timing log above.

# A metadata artifact-MD5 mismatch is rejected while the RDS checksum remains
# valid.  Re-write the metadata sidecar so this is not merely a stale JSON
# checksum case.
trans_output <- file.path(results, "Adams_trans.rds")
trans_metadata <- paste0(trans_output, ".runtime.json")
payload <- jsonlite::fromJSON(trans_metadata)
payload$artifact_md5 <- strrep("0", 32)
jsonlite::write_json(payload, trans_metadata, auto_unbox = TRUE, pretty = FALSE)
writeLines(c(
  paste0("MD5=", unname(tools::md5sum(trans_metadata))),
  paste0("SIZE=", file.info(trans_metadata)$size),
  paste0("PATH=", trans_metadata)
), paste0(trans_metadata, ".md5"))
bad <- try(
  replay_runtime_metadata(trans_output, "Adams", "trans_analysis", log_file),
  silent = TRUE
)
stopifnot(inherits(bad, "try-error"))

zeroimp_output <- file.path(results, "Adams_zeroimp.rds")
zeroimp_metadata <- paste0(zeroimp_output, ".runtime.json")
unlink(c(zeroimp_metadata, paste0(zeroimp_metadata, ".md5")))
bad <- try(
  replay_runtime_metadata(zeroimp_output, "Adams", "zeroimp_analysis", log_file),
  silent = TRUE
)
stopifnot(inherits(bad, "try-error"))
stopifnot(nrow(arrow::read_feather(log_file)) == 2L)
'''
    )
    subprocess.run(
        ["pixi", "run", "Rscript", "--vanilla", str(script), str(R_WORKER), str(root)],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    return pd.read_feather(root / "execution_times.feather")


def _check_pipeline_b_replay_with_rds_and_log() -> None:
    with tempfile.TemporaryDirectory(prefix="ecoda-runtime-replay-") as raw:
        rows = run_r_replay_fixture(Path(raw))
        assert list(rows.columns) == ["dataset", "method", "time_secs", "mem_GB"]
        assert set(rows["method"]) == {"trans_analysis", "zeroimp_analysis"}
        expected = {
            "trans_analysis": (17.25, 2.75),
            "zeroimp_analysis": (23.5, 3.125),
        }
        for method, (time_secs, mem_gb) in expected.items():
            row = rows.loc[rows["method"] == method].iloc[0]
            assert row["dataset"] == "Adams"
            assert row["time_secs"] == time_secs
            assert row["mem_GB"] == mem_gb


def main() -> None:
    worker = load_module("ecoda_python_benchmark_worker", PYTHON_WORKER)
    merge_worker = load_module("ecoda_merge_execution_times", MERGE_WORKER)
    with tempfile.TemporaryDirectory(prefix="ecoda-python-runtime-") as raw:
        root = Path(raw)
        _check_python_cache_hit_replays_without_work(worker, root)
    with tempfile.TemporaryDirectory(prefix="ecoda-python-runtime-default-log-") as raw:
        root = Path(raw)
        _check_python_cache_hit_uses_default_log_file(worker, root)
    with tempfile.TemporaryDirectory(prefix="ecoda-python-runtime-null-memory-") as raw:
        root = Path(raw)
        _check_python_nullable_memory_replays_as_missing(worker, root)
    with tempfile.TemporaryDirectory(prefix="ecoda-python-runtime-invalid-") as raw:
        root = Path(raw)
        _check_python_metadata_rejection(worker, root)
    _check_merge_rejects_duplicate_and_invalid_runtime_rows(merge_worker)
    _check_pipeline_b_replay_with_rds_and_log()
    print("execution-time replay contracts: OK")


if __name__ == "__main__":
    main()
