# Speed up SLURM array-wait in submit scripts (scontrol wait + poll-until-terminal)

## Goal

Replace the `sleep 60` squeue polling loop + fixed `sleep 30` sacct buffer in the 4 login-node submit scripts with:
1. `scontrol wait <job_id>` — event-driven blocking (Slurm >= 23.02; cluster confirmed at 26.05.1), zero polling load on the scheduler, exact job-id match.
2. A bounded poll-until-terminal loop over `sacct` replacing the blind `sleep 30` — returns as soon as accounting has settled (typically 1-5s vs. a fixed 30s) and waits *longer* than 30s when needed, removing today's false-failure risk on large arrays.

The fail-closed sacct gate (every state row must be COMPLETED before any NAS sync) is **unchanged** and remains the single authority.

## Files to change

| File | Change |
|---|---|
| `src/5_run_benchmark_methods/benchmark_submit_common.sh` | Rewrite `benchmark_wait_for_array` (lines 66-92); update header doc (lines 17-20) |
| `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` | Replace lines 108-114 (squeue loop + `sleep 30`); keep gate (115-127) |
| `src/4_cell_type_annotation/2_submit_hpc_array.sh` | Replace lines 122-129 (squeue loop + `sleep 30`); keep gate (131-142) |
| `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` | Replace lines 48-63 (multi-job wait loop + `sleep 30`); keep summary + per-job `-X` COMPLETED check (65-81) |
| `docs/ARCHITECTURE.md` | Update the two explicit "squeue poll" descriptions (lines 239, 367) |

The 3 benchmark callers (`run_python_sample_embedding_methods/`, `run_r_sample_embedding_methods/`, `run_transformation_zeroimp_analysis/` submitters) need **no code changes** — they call the shared helper. Optionally adjust their header comments that say "polled to completion" (e.g. `run_r_sample_embedding_methods/1_submit_hpc_array.sh:25-28`).

## Design

### 1. Wait block (replaces the squeue loop)

```bash
# Event-driven block until the job leaves the scheduler (no polling).
# Exit code deliberately ignored: the fail-closed sacct gate below is the
# authoritative check (covers cancellation, failure, purged controller records).
scontrol wait "${JOB_ID}" > /dev/null 2>&1 || true
```

- For `1_submit_hpc.sh` (multiple independent jobs): `for job_id in "${JOB_IDS[@]}"; do scontrol wait "${job_id}" > /dev/null 2>&1 || true; done` (sequential blocking; already-finished jobs return immediately).
- `scontrol wait` uses exact job ids — also removes the old substring-match risk of `grep -q "${JOB_ID}"` against full `squeue` output.

### 2. Poll-until-terminal (replaces `sleep 30`)

```bash
# sacct may lag a few seconds behind the job leaving the scheduler; poll
# (bounded) until every state row is terminal instead of a blind fixed sleep.
local TAIL_ITER=0
local STATES
while (( TAIL_ITER < 60 )); do  # max 5 min at 5s
  STATES="$(sacct -j "${JOB_ID}" --format=State -n 2>/dev/null || true)"
  if [[ -n "${STATES//[[:space:]]/}" ]] \
     && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
    break
  fi
  sleep 5
  TAIL_ITER=$((TAIL_ITER + 1))
done
```

For `1_submit_hpc.sh`: same loop but per-job, reusing the existing per-job state fetch:
`sacct -j "${job_id}" -X -n -o "State" 2>/dev/null | tr -d ' \n'`; job is "settled" when state is non-empty and not in the non-terminal set. Break when all jobs settled.

### 3. Gate — unchanged

Keep the existing fail-closed checks verbatim in each file: empty-states → `exit 1` without syncing; any row not exactly `COMPLETED` → print `sacct --format=JobID,JobName,State,ExitCode` + `exit 1` without syncing.

## Safety analysis (can this lead to incomplete file writes?)

**No.** Key points:

- **Write completeness is guaranteed by the gate, not by timing.** A task only reaches state `COMPLETED` after its batch step process exited — process exit flushes buffers and closes file descriptors, so all output files (RDS bundles, feathers, exec logs) are fully written before the state is recorded. The gate is unchanged; we only change *when* we notice the job finished, never *what* we verify.
- **The old `sleep 30` was never a write-completion buffer** — it was a sacct-propagation buffer. The poll-until-terminal is strictly safer: on large arrays where accounting takes >30s to finalize all rows, today's code **false-fails and skips the NAS sync**; the poll waits as long as needed and only then runs the gate.
- **Failed/cancelled jobs behave identically**: `scontrol wait` returns non-zero, which we ignore; the gate then sees FAILED/CANCELLED/TIMEOUT and aborts without syncing — same fail-closed outcome as today, just faster.
- **`scontrol wait` errors are harmless**: if the controller has no record (purged) or restarts mid-wait, the `|| true` swallow + accounting-DB-based gate still yields the correct verdict.
- **No change to the sync path**: `benchmark_merge_sync_cleanup` and the rsync steps are untouched; sync still happens only after the gate passes.
- **Login-node etiquette improves**: zero scheduler polling during long runs (AGENTS.md: login node must not execute heavy work; squeue every 60s was light but unnecessary).

## Validation

1. **Syntax**: `bash -n` on all 4 modified scripts (local, no HPC needed).
2. **Login-node smoke test (user runs, needs HPC access)**:
   - `sbatch --parsable --partition=debug-cpu --wrap="sleep 5"` → time `scontrol wait <id>`; check exit code; immediately `sacct -j <id> --format=State -n` to measure the actual propagation lag (validates the 5s poll interval / 5-min bound).
   - Array variant: `sbatch --parsable --partition=debug-cpu --array=1-2 --wrap="sleep 3"` → confirm `scontrol wait` accepts the array master id and blocks until all tasks done.
   - Negative test: `sbatch --parsable --partition=debug-cpu --wrap="exit 1"` → confirm `scontrol wait` returns promptly (non-zero) and sacct shows FAILED so the gate aborts.
3. **End-to-end**: existing `_debug` benchmark smoke run (`--ds_name _debug`) covering submit → wait → gate → NAS sync, per the user's usual workflow.
4. **Fallback if the smoke test shows `scontrol wait` misbehaves** (unexpected exit semantics, unsupported on login node): keep the squeue loop at `sleep 15` + the poll-until-terminal tail. Still removes the fixed 30s and most tail latency; no new primitive needed.

## Out of scope

- Changing the fail-closed gate logic or the merge/sync/cleanup tail.
- `scontrol wait` exit-code semantics as a *fast gate* — the sacct gate stays the single authority.
- Adding timeouts to `scontrol wait` (arrays already carry `--time` limits; old loop also blocked indefinitely).
