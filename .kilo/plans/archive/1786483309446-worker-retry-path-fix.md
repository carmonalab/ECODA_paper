# Fix worker_retry.sh source path in benchmark workers (off-by-one `..`)

## Context / root cause

Commit `5387784` ("HPC worker hardening: self-healing retry, R-lib staging, per-sample checkpoints", added 2026-08-11) made all five array workers source `src/utils/bash/worker_retry.sh`. The three benchmark workers (which sit **two** levels below `src/`, at `src/5_run_benchmark_methods/<dir>/`) were given a **three-level** relative path:

```bash
source "${SCRIPT_DIR}/../../../utils/bash/worker_retry.sh"   # resolves to <PROJECT_ROOT>/utils/bash/... — does not exist
```

Correct path (two levels up → `src/utils/bash/worker_retry.sh`):

```bash
source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"
```

Consequence: every benchmark array task dies at the `source` line with exit 1 (`No such file or directory`), **before** `worker_retry.sh` is loaded, so the transient-failure self-requeue never engages. Annotation/preprocess workers (one level below `src/`) use `../utils/bash/worker_retry.sh` and are correct — do not touch them.

## Changes (3 files, 1 line each)

Replace `source "${SCRIPT_DIR}/../../../utils/bash/worker_retry.sh"` with `source "${SCRIPT_DIR}/../../utils/bash/worker_retry.sh"` in:

1. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh` (line 65)
2. `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh` (line 71)
3. `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh` (line 75)

No other `../../../` occurrences exist anywhere under `src/5_run_benchmark_methods/` (verified by grep) — scope is exactly these three lines.

## Validation (static only — do NOT run pipeline scripts, per AGENTS.md)

1. `bash -n` on all three edited worker scripts.
2. `grep -rn 'worker_retry.sh' src/5_run_benchmark_methods/` → exactly three `../../utils/bash/worker_retry.sh` matches (plus comments); no `../../../` remaining.
3. Path resolution check: from each worker's directory, `../../utils/bash/worker_retry.sh` must resolve to the existing file `src/utils/bash/worker_retry.sh` (e.g. `test -f src/utils/bash/worker_retry.sh` from repo root after mentally applying the relative path).
4. Sanity: `src/4_cell_type_annotation/2.1_run_worker.sh` and `src/3_scrnaseq_preprocessing/1.1_run_worker.sh` keep their correct `../utils/bash/worker_retry.sh` — no diff.

## Commit & push (per AGENTS.md task-completion workflow)

1. Move this plan to `.kilo/plans/archive/`.
2. `git add .`
3. Commit (style: concise, e.g. `Fix worker_retry.sh source path in benchmark workers (off-by-one ..)`).
4. Push.

## User steps after the fix is pushed (context, not implementer work)

1. Ctrl-C the two stuck monitor terminals (they are blocked on the down slurmdbd anyway and would fail closed without syncing).
2. On HPC: `git pull` (if the HPC clone has local uncommitted edits to these worker files, the pull will conflict — resolve or stash first).
3. If the temporary symlink `~/ECODA_paper/utils/bash/worker_retry.sh` was created, remove it (`rm` the link; it is untracked and would otherwise persist).
4. Resubmit both submitters.

## Parallel launch decision (user question: run python + R submitters simultaneously?)

**Yes — launch both in parallel.** Rationale:
- Disjoint resources: python methods wait on H200 GPUs (gpu006); R methods use the EPYC-7742 CPU pool. No interference.
- PILOT (part of the python submission) already shares the CPU class with the R methods regardless of launch order, so serializing buys nothing.
- Only interaction: both submitters end with the shared merge/sync tail (`benchmark/embeddings/execution_times.feather` read-modify-write + `checksums.md5` + rsync). A true simultaneous merge race is unlikely (R run starts immediately and finishes before the H200 wait typically ends) and is recoverable via `--sync-only <job-id>`; the checksum/rsync steps are last-writer-wins. Nothing is lost on failure (fail-closed design, per-task logs deleted only after successful rsync).
- If the H200 node is free, submitting python now lets it start immediately; R uses otherwise-idle CPU capacity.

## Risks / notes

- **slurmdbd was down** (`sacct: Connection refused` on `lunihpcslurm1.admin:6819`) during the incident — check it is back (`sacct -j <id>`) before relying on the final sync gate of a resubmitted run; the arrays themselves run fine while it is down.
- The failure mode tonight (`No such file or directory` at the source line, exit 1, no requeue) is exactly the documented stale-BeeGFS/startup-failure class — after this fix the worker self-healing engages for such cases.
- Out of scope: no other worker or config changes; docs (ARCHITECTURE.md / AGENTS.md) already describe the two-level convention correctly.

## Validation of the actual fix (by user, on HPC)

After `git pull` + resubmission: tasks must pass the `source worker_retry.sh` line (no `No such file or directory` in `~/ECODA_paper/logs/5_benchmark_*.err`), and the R array (CPU, starts immediately) should reach "running {method} on {dataset}" in its `.log`.
