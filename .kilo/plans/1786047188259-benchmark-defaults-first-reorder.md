# Reorder Python benchmark combos so default methods run first

## Context

`1.1.1_benchmark_methods_py.py` runs all combos of one method sequentially in a
single process (one SLURM task per dataset per method). `mem_GB` is logged as
the **process-wide peak RSS** (`resource.getrusage().ru_maxrss`, monotonic
within a process — `log_execution_time()` at 1.1.1_benchmark_methods_py.py:71
calls `peak_rss_gb()`). Memory not returned to the OS between combos (torch/
scvi allocator retention, densified scPoli matrices, cuDNN workspaces) bloats
every later combo's `mem_GB`, and can slow later combos via GC/memory pressure.

The paper's **default (main-method) combos** — the ones in `constants.R`
`method_label_map_main` and the notebook's exec-time figure
(`benchmark_analysis.rmd` ~line 1833, plus the `PILOT_hvg2000_highres` →
`PILOT_hvg2000` recode at line 1793) — currently run **after** other combos:

| Method | Default combo (method string) | Current position in task |
|---|---|---|
| mrvi | `MrVI_hvg2000` (hvg 2000, lowres) | 2nd of 3 (hvg sorted 1000, 2000, 3000) |
| scpoli | `scPoli_hvg2000_dims15_highres` (hvg 2000, highres, dims 15) | last of 6 |
| pilot | `PILOT_hvg2000_highres` (hvg 2000, highres) | 3rd of 4 |

Their timing/RAM is therefore bloated by earlier combos.

## Decision (user-confirmed)

**Reorder only.** No new CLI flags, no changes to `1.1_run_worker.sh` /
`1_submit_hpc_array.sh` / `1.1.2_merge_execution_times.py`. The default combo
runs **first** within each task, so its `mem_GB` = peak of (h5ad load + the
default combo itself) and its `time_secs` is measured without prior-combo
memory pressure. Non-default combos keep their current relative order (stable
sort) and remain subject to the existing in-process bloat (accepted; they are
supplementary parameter-screening combos).

Explicitly out of scope: a `--defaults` run mode / separate defaults pass, and
per-combo process isolation (would change task granularity).

## Task 1 — Reorder combos in `1.1.1_benchmark_methods_py.py`

In `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`:

1. Add constants + a helper near the existing combo helpers (`scpoli_dims_for`
   / `run_pilot_for`, ~line 104–118):

   ```python
   # Default (main-method) combos — constants.R method_label_map_main and the
   # notebook's exec-time figure: MrVI_hvg2000, scPoli_hvg2000_dims15_highres,
   # PILOT_hvg2000_highres.
   DEFAULT_HVG = 2000
   DEFAULT_SCPOLI_DIM = 15
   DEFAULT_RES_LABEL = "_highres"


   def is_default_combo(method, combo):
       n, res_label, _, payload, _, _ = combo
       if method == "mrvi":
           return n == DEFAULT_HVG
       if method == "scpoli":
           return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL and payload == DEFAULT_SCPOLI_DIM
       if method == "pilot":
           return n == DEFAULT_HVG and res_label == DEFAULT_RES_LABEL
       return False
   ```

2. In `process_dataset()`, after the `combos` list is built (after the
   `elif args.method == "pilot":` block, before `pending = []`), insert:

   ```python
   # Defaults-first ordering: ru_maxrss peak RSS is monotonic within a
   # process, so combos run earlier report the least bloated mem_GB (memory
   # leaks / allocator retention from earlier combos would otherwise inflate
   # the defaults' rows). Stable sort: non-default combos keep their order.
   combos.sort(key=lambda c: 0 if is_default_combo(args.method, c) else 1)
   ```

   Note: if `--hvg` excludes 2000 (custom HVG list), no default combo exists
   and the sort is a no-op — no behavior change.

3. Update the module docstring (lines 26–29, the execution-time/memory
   paragraph) with one sentence: combos run defaults-first so the main-method
   rows (`MrVI_hvg2000`, `scPoli_hvg2000_dims15_highres`,
   `PILOT_hvg2000_highres`) are measured before any in-process memory bloat.

Resulting per-task order:
- mrvi: hvg2000, hvg1000, hvg3000
- scpoli: dims15_highres (2000), dims15_lowres (2000), dims2/3/5/10_highres (2000)
- pilot: highres (2000), lowres (2000), highres (1000), highres (3000)

## Task 2 — Update docs

- `docs/ARCHITECTURE.md`, row `1.1.1_benchmark_methods_py.py` (~line 231),
  after the combos description ("...h5ad loaded ONCE per task."): add that
  combos run **defaults-first** (the three main-method combos above) because
  `ru_maxrss` peak RSS is monotonic within a process — earlier combos report
  the least bloated `mem_GB`.
- `TODO.md`: add a short note to the Phase 3.1 bullet (line 42–50) that combos
  run defaults-first for accurate `mem_GB`/`time_secs` of the main methods,
  and add a changelog entry under the Changelog section.

No changes to README.md (does not document combo ordering), the worker, the
submitter, the merge script, or the archived qmd.

## Validation

- `python -m py_compile src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`
- No bash files touched → no `bash -n` needed.
- HPC smoke test (user-run, per AGENTS.md; NOT part of implementation):
  `./1_submit_hpc_array.sh --ds_name _debug --methods mrvi` (then scpoli,
  pilot) — check the merged `execution_times.feather`: the default rows'
  `mem_GB` must be ≤ the same-method non-default rows (i.e. they were measured
  first, unbloated).
