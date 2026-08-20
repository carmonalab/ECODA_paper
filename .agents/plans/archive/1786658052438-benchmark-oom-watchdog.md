# Compute-node watchdog for benchmark OOM auto-escalation

## Context

- The OOM auto-escalation loop (`benchmark_wait_oom_retry`) runs in the **login-node submitter tail**. An SSH drop kills the tail (SIGHUP), so an `OUT_OF_MEMORY` task is never retried and nothing is synced — observed 2026-08-14: `4313942_3` (scitd, Gongsharma) OOM'd at 128G, no retry array was ever submitted, no sync email arrived.
- Workers and their transient self-requeue (`worker_retry.sh`) already survive SSH loss (they run on compute nodes). OOM kills the whole task cgroup, so the worker cannot self-requeue — the escalation must be **external** to the task.
- The final NAS sync can never leave the login node (worker nodes have no NAS mount). The goal is narrower: **the escalation must survive a login-tail death**; the sync is recovered as today (`--sync-only` or an idempotent re-run).
- Login-node policy (AGENTS.md): no long-running login processes, no tmux — so the watchdog must run as a **Slurm batch job on a compute node**.

## Goal

One small "watchdog" batch job per method array. The watchdog owns the terminal-wait + per-task gate + OOM escalation (re-submitting only OOM'd datasets with doubled `--mem`, clamped to `BENCHMARK_MEM_MAX`), and writes a **status file** the login tail reads. The login tail only waits for the watchdog job and then does the merge/sync — so even if the tail dies, the escalation continues unattended and the morning recovery is a re-run or `--sync-only <watchdog_id>`.

## Design decisions

- **Watchdog owns gate + escalation; login tail waits on the watchdog id.** Prevents double-escalation (tail and watchdog both detecting OOM would double-resubmit the same dataset).
- **One watchdog per method array**, submitted by the submitter right after each array (matches the current per-array monitor structure; parallel-friendly).
- **prepare_pseudobulk gets a watchdog too**, with a `soft-gate` mode mirroring today's `benchmark_wait_prep_array` (artifact presence pass first; strict OOM-aware gate fallback). Keeps the gate uniform; the login side no longer holds any gate logic.
- **Status-file protocol**: watchdog writes `${HPC_SCRATCH_DIR}/_benchmark_watchdog/<watchdog_job_id>.status` before exiting (lines: `STATE=OK|FAIL`, `LABEL=`, `JOB_REPORT=<label>|<id>|<wall>` per gated array, `FAIL_REASON=`, `REPORT=<multiline>`). Tail polls for the file (≤2 min grace) after the watchdog job is terminal.
- **Reuse `benchmark_wait_oom_retry` as the watchdog's engine**: add an optional 5th arg `STATUS_FILE`; when set, each terminal path (pass / ceiling-OOM / non-OOM fail / empty sacct / unparseable mem / cap) writes the status file, and the `notify_sync_status` calls are skipped (the login tail emails instead; compute nodes have no mail CLI).
- **Watchdog resubmit closure**: `watchdog_main.sh` receives the method's partition/throttle/log-prefix/worker-path/flags as positional args (after `--`) and defines its own `watchdog_resubmit <label> <ds_csv> <mem> <manifest>` (mirrors `submit_method_array_retry`). `METHOD`/`BENCHMARK_MANIFEST`/`FORCE_BENCHMARK` propagate via sbatch env inheritance (submitter → watchdog → retry workers).
- **Watchdog job specs**: 1 task, 1 cpu, `--mem=2G`, `--time=${WATCHDOG_TIME_LIMIT:-12:00:00}` — **12h is the `shared-*` partition MaxTime** (the workers' `#SBATCH --time=12:00:00` is documented as the partition max); a higher value is rejected by sbatch at submit time (fail-fast, visible in the submitter output), so `WATCHDOG_TIME_LIMIT` must never be set above the target partition's MaxTime. Default partition (not the pinned benchmark class), logs to `${LOGS_DIR}/5_benchmark_watchdog_<label>_%A.log/.err`. Pure bash + slurm CLI (`squeue`/`sacct`/`sbatch`/`scontrol`) — no pixi/R/Python needed.
- **`--sync-only`**: accepts watchdog job ids (wait terminal → read status → OK: proceed to merge/sync; FAIL: fail closed with the report). Array ids keep today's strict `benchmark_wait_for_array` gate (unchanged).
- **`benchmark_wait_for_array` and `benchmark_merge_sync_cleanup` unchanged**; trans/zeroimp submitter untouched.

## Implementation tasks

### 1. `src/5_run_benchmark_methods/benchmark_submit_common.sh`
- Add `WATCHDOG_STATUS_DIR="${HPC_SCRATCH_DIR}/_benchmark_watchdog"` (derived, like `BENCHMARK_MERGE_SCRIPT`); never rsync'd (outside `benchmark/`).
- Add `benchmark_pb_variant_names <benchmark_hpc_utils_path>` (extract the existing sed/grep parse from the R submitter's `benchmark_wait_prep_array`).
- Extend `benchmark_wait_oom_retry <job_id> <label> <resubmit_fn> <manifest> [status_file]`:
  - When `status_file` is set: write the status file at every terminal path (OK with one `JOB_REPORT=` line per gated array incl. the final retry id; FAIL with `FAIL_REASON` + the per-task report) and skip `notify_sync_status` calls. Behavior without the arg unchanged.
- Add `benchmark_submit_watchdog <array_id> <label> <manifest> <mode:strict|soft-gate> <partition> <throttle> <log_prefix> <worker_script> <flags...>`:
  - `mkdir -p` the status dir, `sbatch` `watchdog_main.sh` with the args (1 cpu/2G/`WATCHDOG_TIME_LIMIT`, default 12h = shared-* MaxTime), echo ONLY the watchdog job id on stdout (progress to stderr). No status-file arg: the watchdog names its own file `${WATCHDOG_STATUS_DIR}/<watchdog_job_id>.status` from `SLURM_JOB_ID` at runtime (the id is unknowable at submit time).
- Add `benchmark_wait_watchdog <watchdog_id> <label>`:
  - `benchmark_wait_array_terminal` on the watchdog id → poll for `${WATCHDOG_STATUS_DIR}/<watchdog_id>.status` (≤2 min grace) → parse `STATE`:
    - OK → append its `JOB_REPORT=` lines to `JOB_REPORTS`, return 0.
    - FAIL → print report, `notify_sync_status` ("NOT synced — watchdog failed", with the report), exit 1.
    - Watchdog job non-COMPLETED (FAILED/TIMEOUT/PREEMPTED/CANCELLED — e.g. node panic or preemption before the status file was written) → fail closed reporting its sacct `State,ExitCode` + a pointer to `${LOGS_DIR}/5_benchmark_watchdog_<label>_<id>.*`.
    - Watchdog COMPLETED but status file missing/empty after the grace → fail closed ("watchdog exited without a status file; check its logs").

### 2. New `src/5_run_benchmark_methods/watchdog_main.sh`
- Standard worker boilerplate: `set -euo pipefail`, `SCRIPT_DIR` recovery via `scontrol show job` `Command=` (sbatch spool) with `BASH_SOURCE` fallback, `source slurm_config.sh`, `cd "${PROJECT_ROOT}"`, `source benchmark_submit_common.sh`.
- Args: `<array_id> <label> <manifest> <mode> -- <partition> <throttle> <log_prefix> <worker_script> <flags...>`.
- Defines its own status path `${WATCHDOG_STATUS_DIR}/${SLURM_JOB_ID}.status` (same derivation as the login tail's `benchmark_wait_watchdog`).
- Define `watchdog_resubmit <label> <ds_csv> <mem> <manifest>`: writes the reduced manifest (`benchmark_manifest_<label>_retry_$$.txt`), exports `METHOD`/`BENCHMARK_MANIFEST`, `sbatch`es the retry array with the given partition/flags/throttle/`--mem`/log patterns/worker script + `--mail-user`, echoes the job id only.
- `mode=strict`: call `benchmark_wait_oom_retry <array_id> <label> watchdog_resubmit <manifest> <status_path>`.
- `mode=soft-gate` (prepare_pseudobulk): `benchmark_wait_array_terminal` → `benchmark_pb_variant_names "${SCRIPT_DIR}/benchmark_hpc_utils.R"` → all variants present per dataset → write OK status (JOB_REPORT for the prep array) → exit 0; else fall through to `benchmark_wait_oom_retry` with the status file.
- No emailing, no NAS access, no exec-log merging.

### 3. `run_r_sample_embedding_methods/1_submit_hpc_array.sh`
- After each method array submission: `WATCHDOG_ID="$(benchmark_submit_watchdog "${ARRAY_JOB_ID}" "${METHOD}" "${MANIFEST}" strict "${PARTITION}" "${MAX_NUM_CHUNKS_PARALLEL}" "5_benchmark_r_${METHOD}" "${SCRIPT_DIR}/1.1_run_worker.sh" "${EXTRA_FLAGS[@]}")"`; track per-method watchdog ids (parallel to `ARRAY_JOB_IDS`/`ARRAY_JOB_METHODS`).
- Prep flow: submit prep array → submit prep watchdog (`soft-gate` mode) → `benchmark_wait_watchdog <prep_watchdog_id> prepare_pseudobulk` → only then submit method arrays.
- Monitor tail: replace `benchmark_wait_oom_retry ...` per job with `benchmark_wait_watchdog <watchdog_id> <method>`.
- Delete `benchmark_wait_prep_array` (its soft-gate + strict fallback logic moves to `watchdog_main.sh` soft-gate mode). `submit_method_array_retry` stays for the normal submission path (the watchdog uses its own closure).
- `--sync-only`: if the id has a status file → watchdog path (wait + read); else strict `benchmark_wait_for_array` (unchanged).

### 4. `run_python_sample_embedding_methods/1_submit_hpc_array.sh`
- Same pattern as (3), with per-method `PARTITION`/`EXTRA_FLAGS`/`THROTTLE` (GPU flags preserved) passed to `benchmark_submit_watchdog`; monitor tail switches to `benchmark_wait_watchdog` per method.

### 5. Docs
- `AGENTS.md` (worker self-healing paragraph): benchmark OOM escalation runs as a **compute-node watchdog job** (`watchdog_main.sh`, one per method array) that survives login-tail/SSH loss; the login tail only waits for the watchdog and syncs; `--sync-only <watchdog_id>` resumes after a tail death. Note `WATCHDOG_TIME_LIMIT` env.
- `docs/ARCHITECTURE.md`: R/python workflow blocks (submit array → submit watchdog → gate on watchdog id → merge/sync), `benchmark_submit_common.sh` row (`benchmark_submit_watchdog`/`benchmark_wait_watchdog`/status-file arg), new `watchdog_main.sh` row, `_benchmark_watchdog/` dir entry in the HPC layout.

## Validation

1. `bash -n` on all edited/new scripts.
2. Mocked unit tests (extend the existing `test_oom_retry.sh` harness; stubs for `sacct`/`squeue`/`sbatch`/`notify_sync_status`):
   - `benchmark_wait_oom_retry` with `status_file`: OK path writes `STATE=OK` + `JOB_REPORT=` incl. the FINAL retry id; FAIL paths (ceiling OOM with MaxRSS, non-OOM, empty sacct, unparseable mem) write `STATE=FAIL` + reason + report; no `notify_sync_status` calls when a status file is set.
   - `watchdog_main.sh` (extracted): strict mode → `watchdog_resubmit` builds the correct sbatch (stubbed `sbatch`; assert `--mem` doubling + clamp, reduced manifest contents, flags preserved); soft-gate mode → variants present = OK status without gate; variants missing = strict gate path.
   - `benchmark_wait_watchdog`: OK status → passes + `JOB_REPORTS` merged from file; FAIL → email + exit 1; watchdog non-COMPLETED → fail closed with state + log pointer; watchdog COMPLETED but no status file → fail closed.
   - `benchmark_submit_watchdog`: stubbed sbatch → echoes only the watchdog id; args correctly forwarded.
3. **Critical gate — sbatch from a compute node (user, HPC; run in parallel with items 1–2)**: submit a trivial debug-cpu job that performs the full client round-trip the watchdog needs — `sbatch --wrap 'echo ok'`, then `squeue`/`sacct` the inner job to COMPLETED, plus `scontrol show job <id>` — and reports the inner job id. If compute nodes cannot reach slurmctld as a client, **stop and reassess**: the design is not viable (fallback: keep the escalation login-side and document the nohup/tmux exception for unattended runs).
4. Smoke test (user, HPC): `--ds_name _debug --methods pilot` (python) and `--methods gloscope --ds_name _debug` (R) — full cycle with watchdogs; verify status files, no sbatch rejection of the watchdog's `--time`, and the sync email.
5. **Acceptance test** (user, HPC): run `--methods gloscope,scitd`, then Ctrl-C the login tail right after submission (simulating SSH drop). The watchdog must still: gate gloscope (all COMPLETED), detect Gongsharma scitd OOM at 128G, re-submit at 256G, complete, write `STATE=OK`. Then `--sync-only <watchdog_id>` syncs to NAS and emails.

## Risks / notes

- **sbatch-from-compute-node is a hard prerequisite** (validation item 3). If the cluster forbids it, fall back to documenting the nohup/tmux exception for unattended runs.
- Retry arrays submitted by the watchdog are ordinary method arrays — workers unchanged (skip-if-exists idempotency applies).
- Watchdog job specs must be modest (1 cpu/2G) so it never competes with method tasks or distorts the pinned-class picture.
- If the watchdog job itself fails (node issue), the login tail fails closed on the watchdog id with a pointer to its `.err` — recover with `--sync-only <array ids>`.
- If the tail dies and the user re-runs the same command while the old watchdog is still escalating, the new run's arrays skip existing results and its new watchdogs gate cleanly; at worst a benign double-compute of the OOM'd dataset (atomic writes, same exec-log row).
- **`WATCHDOG_TIME_LIMIT` (default 12h = shared-* MaxTime) bounds the watchdog.** An escalation chain can outlive it: an array running near its full 12h before an OOM could leave a retry (up to 12h more) beyond the watchdog's own limit. The retry array is an independent Slurm job and keeps running — only the gating/status is lost; the login tail fails closed on the watchdog id (state-aware message, see task 1) and recovery is `--sync-only` with the retry array id (from `sacct`) or a fresh re-run (idempotent). Observed OOMs died early (4 min / 1h45 in), so this is a tail-risk, not the norm. Setting `WATCHDOG_TIME_LIMIT` above the partition MaxTime is rejected by sbatch at submit time (fail-fast, visible in the submitter output).
- Stale status files in `_benchmark_watchdog/` are harmless (job-id named) and never synced.

## Out of scope

- NAS sync on compute nodes (impossible — login-only mount); the sync tail remains login-bound with `--sync-only`/re-run recovery.
- trans/zeroimp submitter (fixed 32G, obs-only reads).
- Watchdogs for stages 1–4 arrays (preprocess/annotation have their own mechanisms).
- Emails from the watchdog (login tail only).
