# Auto-escalate BENCHMARK_MEM on OOM (benchmark submitters)

## Context

- `4313847_3` (scitd, Gongsharma_cmv_young_males) was OOM-killed at the 128G task cap: sacct shows `State=OUT_OF_MEMORY`, `ExitCode 0:125`, `MaxRSS 127,505,612 KB` (~122 GiB). All other 10 scitd tasks COMPLETED.
- Root cause: `process_scitd_fig()` (`benchmark_methods_r.R:252`) passes the **full counts matrix** into scITD's `make_new_container` (which densifies) plus a 14-core parallel Tucker/ICA, on top of the already-loaded full Seurat object. Large datasets (Gongsharma; future additions) exceed 128G. gloscope/mofa fit because they don't add scITD's working set.
- Current behavior: an OOM'd task fails the submitter's fail-closed gate → the whole submission aborts without syncing → the user must notice, manually re-run with a higher `BENCHMARK_MEM` (which today is also not env-overridable — `slurm_config.sh:118` assigns unconditionally).
- Goal: within ONE submission, detect OOM tasks via sacct, automatically re-submit **only those tasks' datasets** with doubled `--mem`, up to a ceiling; fail closed only when the ceiling is reached or a non-OOM failure occurs.

## Decisions

- **Escalation: double each retry** (user-confirmed): 128G → 256G → 512G. Ceiling `BENCHMARK_MEM_MAX` (env-overridable, default `512G`, validated against actual node RAM — see Validation).
- **Scope**: `run_r_sample_embedding_methods/` (gloscope/mofa/pseudobulk/scitd + prepare_pseudobulk) and `run_python_sample_embedding_methods/` (mrvi/scpoli/pilot) — both use `BENCHMARK_MEM`. `run_transformation_zeroimp_analysis/` (hardcoded 32G, obs-only reads) is out of scope.
- **Retry granularity**: per-method reduced manifest containing only the OOM'd datasets; array task ids map 1:1 to it (existing manifest mechanism). Workers are idempotent (skip-if-exists), so non-OOM'd datasets are never re-run.
- **Non-OOM failures** (FAILED/CANCELLED/TIMEOUT/…): fail closed exactly as today (per-task report + sync-status email). OOM is the only auto-retried state.
- **`--sync-only` mode**: unchanged (no submissions happen there; provided ids are gated strictly as today).
- **No cross-run persistence** of observed memory needs: each OOM costs minutes (tasks die early, e.g. 4:28 for Gongsharma scitd), so per-run escalation suffices.
- **No partition switch** (e.g. `shared-bigmem`) for retries: it would break the pinned-hardware runtime-comparability premise (EPYC-7742 constraint).

## Implementation tasks

### 1. `src/slurm_config.sh`
- Line 118: `BENCHMARK_MEM="128G"` → `BENCHMARK_MEM="${BENCHMARK_MEM:-128G}"` (honors the existing "All env-overridable" comment; currently a prefix override is clobbered on source).
- Add next to it: `BENCHMARK_MEM_MAX="${BENCHMARK_MEM_MAX:-512G}"` with a comment (doubling ceiling for OOM auto-escalation; per-command override, e.g. `BENCHMARK_MEM_MAX=256G ./1_submit_hpc_array.sh`).

### 2. `src/5_run_benchmark_methods/benchmark_submit_common.sh`
- Add helpers:
  - `benchmark_bump_mem <mem>` — parse `<N>G`/`<N>T`, return `2N` with the same suffix; unparseable → non-zero exit (caller fails closed).
  - `benchmark_mem_ge <a> <b>` — compare two mem strings.
- Extract the shared "wait until terminal" poll (squeue exact-id poll at 60 s + bounded sacct poll-until-terminal) into a helper used by `benchmark_wait_for_array` (behavior unchanged) and the new function.
- Add `benchmark_wait_oom_retry <job_id> <label> <resubmit_fn> <manifest>`:
  1. Wait until the array leaves the scheduler (shared poll).
  2. Read per-task states: `sacct -j <id> --format=JobID,State -n`, keeping only task rows (`<jobid>_<n>`, excluding `.batch`/`.extern` and the master row).
  3. If all COMPLETED → record in `JOB_REPORTS` (as `benchmark_wait_for_array` does), return 0.
  4. If any non-COMPLETED, non-OOM state → fail closed with the existing `benchmark_task_report` + `notify_sync_status` path (exit 1, no sync).
  5. If OOM tasks exist:
     - If `benchmark_mem_ge` current mem `BENCHMARK_MEM_MAX` (or bump fails) → fail closed with an OOM report including per-task `MaxRSS` (`sacct --format=JobID,State,MaxRSS,Elapsed`), exit 1.
     - Else: bump mem (×2), map OOM task ids → datasets via `sed -n "${t}p" "${manifest}"`, echo an escalation line (label, datasets, old→new mem), call `${resubmit_fn} <label> <comma-separated-datasets> <new_mem> <new_manifest_path>`, which echoes the new array job id; loop back to step 1 with the new id/manifest. Cap total attempts (e.g. 4) as belt-and-braces.
- Keep `benchmark_wait_for_array` untouched (still used by `--sync-only` paths and trans/zeroimp).

### 3. `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh`
- Factor the sbatch block of `submit_method_array()` into a shared `submit_method_array_retry <method> <ds_csv> <mem> <manifest_path>` (same partition/flags/constraint/throttle as the method's normal submission; `--mem` from the argument; manifest = only the given datasets; echoes the job id). The normal `submit_method_array()` delegates to it with the full dataset list and `${BENCHMARK_MEM}`.
- Method arrays (gloscope/mofa/pseudobulk/scitd): replace `benchmark_wait_for_array` in the monitor tail with `benchmark_wait_oom_retry <id> <method> submit_method_array_retry <manifest>`.
- `benchmark_wait_prep_array()`: keep the artifact soft-pass first (all variants present → pass, as implemented). When variants are missing, replace the strict `benchmark_wait_for_array` fallback with the OOM loop (escalate OOM'd prep tasks; fail closed on non-OOM bad states; all-COMPLETED passes).

### 4. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`
- Same pattern: a retry-capable submit function per method (preserving the per-method partition/`--gpus`/`--constraint`/throttle flags), reduced manifest, and `benchmark_wait_oom_retry` in the monitor tail instead of `benchmark_wait_for_array`.

### 5. Documentation
- `AGENTS.md`: extend the worker-self-healing paragraph with one sentence: OOM kills cannot self-requeue (task is dead), so the benchmark submitters auto-detect `OUT_OF_MEMORY` tasks via sacct and re-submit those datasets with doubled `--mem` (up to `BENCHMARK_MEM_MAX`, default 512G) before failing closed.
- `docs/ARCHITECTURE.md`: R/python submitter workflow blocks ("poll + sacct gate") → note the OOM auto-escalation; `benchmark_submit_common.sh` row → document `benchmark_wait_oom_retry` + mem helpers.

## Validation

1. `bash -n` on all edited scripts.
2. Mocked unit tests (pattern of the prep-gate test used during implementation):
   - bump/compare helpers (128G→256G, 512G ceiling, unparseable input fails);
   - OOM loop with stubbed `sacct`/`squeue`/resubmit callback: all-COMPLETED → pass; OOM → escalate + resubmit with only the OOM'd datasets → pass on retry; OOM at ceiling → fail closed with MaxRSS report; mixed OOM+FAILED → fail closed immediately; task-row filtering excludes `.batch`/`.extern`/master rows.
3. `sinfo -p shared-cpu -o "%n %m %e"` and `sinfo -p shared-gpu -o "%n %m %e"` on the HPC → confirm the 512G default ceiling fits node RAM (adjust `BENCHMARK_MEM_MAX` default if nodes are smaller; a 256G request must be schedulable).
4. Live validation (after commit/push/pull on the HPC): re-run `./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods gloscope,scitd` — gloscope + 10 scitd datasets skip (results exist), Gongsharma scitd must auto-escalate 128G→256G and complete, then the merge/sync tail runs normally. Also covers the pending Gongsharma scitd gap.
5. If the in-flight mofa/pseudobulk run OOMs on Gongsharma pseudobulk (same code path loads the full Seurat): its rerun (`--methods mofa,pseudobulk`) auto-escalates the same way.

## Risks / notes

- Escalated `--mem` requests occupy more of a node → the retried single task may sit PENDING a while; acceptable and visible in the escalation log lines.
- The OOM'd attempt's partial exec-log rows are overwritten by the retry worker (merge dedup keep=last) — safe.
- The retried array's logs use its own job id → no collision with the original array's logs.
- Do not switch retries to `shared-bigmem` (breaks pinned-hardware comparability).
- Keep `benchmark_wait_for_array` and the prep artifact soft-gate semantics intact; the OOM loop composes on top.
- Git hygiene: the user commits/pushes locally and pulls on the HPC (repo convention); the current running submissions are unaffected by the code change.

## Out of scope

- trans/zeroimp submitter (fixed 32G, obs-only reads).
- Cross-run persistence of per-dataset/method memory observations.
- `--sync-only` OOM handling (no submissions happen there).
- TIMEOUT retries (would need `--time` escalation, not memory).
