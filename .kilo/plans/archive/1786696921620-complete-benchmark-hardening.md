# Complete benchmark hardening: TRANSIENT_REQEX + squeue-poll (post-watchdog review, v2)

## Review verdict (updated 2026-08-14, after HPC validation)

The OOM-watchdog plan was implemented (`12fbd13`), reviewed (bash -n clean, 103/103 mocked
assertions), and **validated on the HPC**. Two post-implementation bugs were found and fixed
during validation:

1. **`watchdog_main.sh:42` relative-depth bug** — `source "${SCRIPT_DIR}/../../slurm_config.sh"`
   resolves to the repo root for a script one level below `src/` (the workers sit two levels
   below and `../../` is correct for them). Symptom: watchdog jobs died at startup (`ExitCode 1`,
   `Run time 00:00:00`, empty `.log`, ~163-byte `.err`), so the tails fail-closed ("watchdog
   lost"). **Fixed in `450d1cd`** (`../slurm_config.sh`). The mocked suite could not catch it
   (it awk-extracts only the watchdog *functions*, never executing line 42; `bash -n` is
   syntax-only) — the HPC smoke test was the right detector.
2. **`--sync-only` unbound-variable crash in both submitters** — in `--sync-only` mode the
   gating loop inside the branch runs, then the shared tail loop after the `if/else` iterates
   `ARRAY_JOB_IDS` (non-empty) and dereferences `WATCHDOG_JOB_IDS[$i]`/`ARRAY_JOB_METHODS[$i]`,
   which are declared only in the submission branch → `unbound variable` under `set -u` AFTER a
   successful gate, BEFORE merge/sync. Symptom: `--sync-only 4314217` gated `STATE=OK` then
   crashed at `line 295: WATCHDOG_JOB_IDS[$i]: unbound variable`; nothing synced. **Fixed in
   `dfdf56b`** (tail gating loop wrapped in `if [[ -z "${SYNC_ONLY_IDS}" ]]`; python:293,
   R:334) and re-validated on the HPC (gate → merge → sync → email → exit 0).

**HPC validation status:**

| Item | Status |
|---|---|
| sbatch-from-compute-node gate (probe 4314199/4314200, inner job from cpu002, squeue/sacct/scontrol round-trip) | ✅ PASSED |
| Smoke: `--ds_name _debug --methods pilot` (4314216/4314217) | ✅ PASSED (STATE=OK, synced, emailed) |
| Smoke: `--ds_name _debug --methods gloscope` (4314218/4314219) | ✅ PASSED |
| Smoke: `--ds_name _debug --methods prepare_pseudobulk` soft-gate watchdog (4314235/4314236) | ✅ PASSED (soft-gate mode, STATE=OK) |
| `--sync-only <watchdog_id>` resume (4314217) | ✅ PASSED (after `dfdf56b`) |
| **Acceptance: `--methods gloscope,scitd` + Ctrl-C the tail + OOM escalation + `--sync-only` recovery** | ⏳ **PENDING — the only remaining validation** |

Remaining tasks of this plan:

1. `TRANSIENT_REQEX` extension — `src/utils/bash/worker_retry.sh:51` still lacks the arrow
   `.onLoad` signatures (observed 2026-08-14: `Error: .onLoad failed in loadNamespace() for
   'arrow' ... attempt to apply non-function` — stale BeeGFS node view; not matched today).
2. Squeue-poll hardening — `benchmark_wait_array_terminal`
   (`benchmark_submit_common.sh`) still treats a transient empty/error `squeue` response as
   "job left the scheduler". The watchdog conversion already made the downstream bounded sacct
   poll (20 min) the safety net, but a false early exit still misreports, and squeue+sacct both
   misbehaving > 20 min would read a RUNNING task as `BAD_TASK`. The 2026-08-14 prep-gate
   incident was the same root cause.
3. **NEW (fold-in): `test_oom_retry.sh` squeue stub + poll assertions** — required for task 2,
   see below.

## Tasks

### 1. TRANSIENT_REQEX: add the arrow .onLoad stale-view signatures

`src/utils/bash/worker_retry.sh:51` — append two alternatives to the existing pattern:

```bash
TRANSIENT_REQEX='No such file or directory|cannot open file|cannot open shared object|package or namespace load failed|there is no package called|missing from the pixi environment|cannot open connection|failed to load|No module named|\.onLoad failed in loadNamespace|attempt to apply non-function'
```

Also:
- Update the header comment of worker_retry.sh (the "TRANSIENT_REQEX" description) to mention
  the arrow `.onLoad` stale-view signature class.
- Update the AGENTS.md worker self-healing paragraph (it enumerates example signatures) to
  include `.onLoad failed in loadNamespace` / `attempt to apply non-function`.

Deliberately NOT added: `error reading from connection` / `unknown input format`
(corrupt-RDS signatures) — a genuinely corrupt cache file is not healed by requeue; the retry
cap (WORKER_MAX_RETRIES=6) and the visible `.err` keep failures diagnosable.

### 2. Harden the squeue poll in benchmark_wait_array_terminal

`src/5_run_benchmark_methods/benchmark_submit_common.sh` — replace the raw
`while squeue ... | grep -qx` loop with a consecutive-clean-miss loop:

- squeue non-zero exit (or non-empty stderr) → treat as **still present** (warning to stderr,
  reset the miss counter, keep polling);
- clean response containing the id → reset the miss counter, keep polling;
- clean response without the id → increment the counter; **exit the loop only after 2
  consecutive clean misses** (adds ≤ 60 s to the normal completion path).

Sketch (fits the existing function; `set -euo pipefail` safe — the command substitution sits in
an `if` condition, and `(( MISS >= 2 )) && break` is exempt as the first command of an `&&`
list):

```bash
local MISS=0 OUT=""
while :; do
  if ! OUT="$(squeue -u "$USER" -h -o "%A" 2>/dev/null)"; then
    echo "WARNING: squeue query failed while waiting for ${JOB_ID}; keeping polling." >&2
    MISS=0
  elif grep -qx "${JOB_ID}" <<< "${OUT}"; then
    MISS=0
  else
    MISS=$((MISS + 1))
  fi
  (( MISS >= 2 )) && break
  sleep 60
done
```

Keep the existing "left the scheduler" echo and the unchanged bounded sacct terminal-state
poll below it. This single fix point covers the login tail AND the watchdog (both call
`benchmark_wait_array_terminal` via `benchmark_wait_for_array` / `benchmark_wait_oom_retry` /
`benchmark_wait_watchdog`).

### 3. test_oom_retry.sh: update the squeue stub + add 2-miss poll assertions

**Required — without the stub change the suite HANGS under task 2:** the harness stubs
`squeue() { return 1; }` ("immediately terminal"), but the new poll treats a failing squeue as
"still present" with `MISS=0` reset every iteration and `sleep` stubbed → infinite loop.

- Change the stub to a clean empty response: `squeue() { return 0; }` (comment: "clean empty
  response -> 2 consecutive clean misses -> terminal"). The existing 103 assertions must still
  pass unchanged (the poll now exits after 2 instant misses).
- Add assertions for the miss logic by redefining `squeue()` inside the test block (functions
  resolve at call time) and counting invocations via a capture file (`${CAPTURE_DIR}/squeue_calls.txt`):
  1. clean-empty from the start → exits after exactly **2** squeue calls.
  2. id present once, then clean-empty → exits after **3** calls (2 misses after the id
     disappears).
  3. one squeue failure (return 1) then clean-empty → failure resets the counter → exits after
     **3** calls.
  Assert exit 0 and the "left the scheduler." line in each case.

Out of scope (note only): the identical raw pattern in the stage-2/3/4 submitters
(`src/2_dataset_specific_preprocessing/1_submit_hpc.sh:137`,
`src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh:191`,
`src/4_cell_type_annotation/2_submit_hpc_array.sh:274`) — same latent early-exit risk but
separate downstream gates; not touched in this round.

## Validation

1. `bash -n src/utils/bash/worker_retry.sh src/5_run_benchmark_methods/benchmark_submit_common.sh src/5_run_benchmark_methods/test_oom_retry.sh`
2. Signature regression check:
   `printf 'Error: .onLoad failed in loadNamespace() for %s '\''arrow'\'', details:\n  error: attempt to apply non-function\n' '' | grep -Eiq "$(grep -o "TRANSIENT_REQEX='[^']*'" src/utils/bash/worker_retry.sh | cut -d"'" -f2)"` → exit 0.
3. `bash src/5_run_benchmark_methods/test_oom_retry.sh` → ALL assertions pass (103 + the 3 new
   poll assertions), no hang.
4. User (HPC) — acceptance test (the only remaining validation; run AFTER this round is
   committed + pulled):
   - `squeue -u $USER` — if in-flight old-code arrays are still listed, wait (concurrent sync
     tails could race on the exec-log merge).
   - `source src/slurm_config.sh && cd "${PROJECT_ROOT}" && ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods gloscope,scitd`
   - Ctrl-C the login tail after BOTH `gloscope watchdog job ID:` and `scitd watchdog job ID:`
     print (if too early, re-run — idempotent).
   - The watchdogs must: gate gloscope (all COMPLETED → `STATE=OK`), detect the Gongsharma
     scitd OOM at 128G, re-submit at 256G, complete, write `STATE=OK` with TWO `JOB_REPORT=`
     lines (original + retry id). Escalation evidence:
     `grep "OOM escalation" "${HOME}"/ECODA_paper/logs/5_benchmark_watchdog_scitd_*.log`
   - Recover: `--sync-only <gloscope_wd>,<scitd_wd> --methods gloscope,scitd` → merges, rsyncs
     to NAS, sync email. Verify the Gongsharma scitd RDS bundle on the NAS.
   - Fallback if Gongsharma does not OOM at 128G this time: the run still validates the
     watchdog happy path; to force escalation re-run with `BENCHMARK_MEM=96G ... --methods scitd --force`.
5. Worker retry sanity (optional): if a task fails with the arrow `.onLoad` signature again,
   confirm the retry counter bumps (self-requeue) instead of an immediate FAILED.

## Rollout / git hygiene

- HPC clone is currently clean (last pull was a fast-forward of `450d1cd..dfdf56b`); no local
  modifications are expected. Before pulling this round: `git status --short` (if anything
  shows as modified, `git stash` or `git checkout -- <file>` first).
- Ordering: implement + commit + push this round FIRST, then pull on the HPC, then run the
  acceptance test (final full-stack validation of the committed state).

## Commit / archive (AGENTS.md task-completion workflow)

- Archive this plan file to `.kilo/plans/archive/`.
- `git add` only: `src/utils/bash/worker_retry.sh`,
  `src/5_run_benchmark_methods/benchmark_submit_common.sh`,
  `src/5_run_benchmark_methods/test_oom_retry.sh`, `AGENTS.md`, this plan (archived).
  Do NOT touch the submitters (sync-only fix already committed in `dfdf56b`) nor the unrelated
  untracked plan `.kilo/plans/1786651957910-pilotgm-qot-benchmark-implementation.md`.
- Commit (repo style, e.g. "Retry arrow .onLoad stale-view signatures + harden squeue poll
  against transient empty responses") and push.

## Out of scope / decisions

- **Mixed OOM + non-OOM failures** in `benchmark_wait_oom_retry` (the first `BAD_TASK` still
  short-circuits before the OOM escalation): conscious non-goal this round — recovery is an
  idempotent rerun. Revisit if it recurs.
- Corrupt-RDS retry signatures: not added (see Task 1 rationale).
- Stage-2/3/4 submitter squeue polls: not touched (own gates; see Task 2 note).
- Exec-log bookkeeping quirk observed during validation (a `--sync-only` merge may report "No
  per-task execution logs found" and sync a stale leftover per-task feather from an earlier
  failed-closed run): benign — the merged `execution_times.feather` is the authoritative log
  and was written unchanged; no action this round.
- The watchdog itself: no further code changes expected; remaining work is the acceptance test.
