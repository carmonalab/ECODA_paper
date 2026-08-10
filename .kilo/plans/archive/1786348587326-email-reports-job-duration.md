# Email reports with job durations + per-dataset/per-method breakdowns

## Goal

1. Every completion email must report **how long each job took** (wall-clock, from
   `sacct` `Elapsed`; array = per-task rows + total array wall time).
2. Add the detailed report format everywhere: **finished? synced to NAS? per
   dataset (or per method) completed/failed**, with durations.
3. Add emails to the pipelines where they are missing entirely
   (`1_stage_data.sh`, `2_dataset_specific_preprocessing/1_submit_hpc.sh`).

## Current state (verified in code)

- `notify_sync_status` (`src/utils/bash/sync_status_email.sh`) is sourced by:
  preprocess array submitter, annotation array submitter, merge script,
  `benchmark_submit_common.sh` (all 3 benchmark submitters). Emails exist but
  contain only dataset/method name lists + sacct excerpt on gate failure — no
  durations, no per-item status (except `3_submit_merge.sh` default-all
  summary: OK/FAILED dataset lists).
- Missing entirely: `src/1_stage_data/1_stage_data.sh` (login-node staging,
  no jobs) and `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
  (dispatcher; header comment says "no status email is sent" — stale after
  this change).
- `sacct Elapsed` is available per task and for the array master (row with the
  bare job id). Task→name mappings:
  - Preprocess: task `i` → i-th line of `jq -r 'keys[]'` (full sorted list;
    valid in default-all AND single-dataset modes, non-underscore keys are a
    prefix).
  - Annotation: task `i` → i-th line of `${HPC_SCRATCH_DIR}/chunks_manifest.txt`
    (field 1 = DS_NAME); manifest persists for `--sync-only` resume.
  - Benchmark (all 3 submitters): task `i` → `DATASET_NAMES[i-1]`
    (manifest order).
  - Dispatcher: plain jobs, no mapping needed (JobName + State + Elapsed).

## Design decisions

- **Duration = wall-clock** (`sacct Elapsed`, HH:MM:SS as-is). Per-task rows
  filtered to `^<job>_[0-9]+[[:space:]]` (excludes `.batch`/`.extern`/`.0` step
  rows); array total from the bare-`<job>` master row. `n/a` when sacct is
  empty/purged (existing fail-closed paths unchanged).
- **Annotation granularity**: per-dataset aggregation (chunk count
  COMPLETED/FAILED + total task elapsed), because chunk arrays can be
  thousands of tasks (spam). Failed chunk task ids listed, bounded to 20 +
  "+ N more".
- **Email format** (success and failure):

  ```
  Pipeline: <label> (datasets: A, B, C)
  Job: <id>
  Status: all tasks COMPLETED — synced to NAS / NOT synced — <reason>
  Per-task report (task → dataset, state, elapsed, exit code):
    1 → Joanito       COMPLETED    00:12:33
    2 → Smillie       FAILED (1)   01:04:11
  Array wall time: 03:12:45
  ```

- Failure emails keep the current "NOT synced — reason" subject line style and
  replace the raw `sacct --format=JobID,JobName,State,ExitCode -n` dump with
  the mapped report (ExitCode kept in the per-task line).
- Success emails get the same per-item breakdown + array wall time appended to
  the existing body (merge script: per-view/per-dataset durations measured
  with local `date +%s` around each `srun`, since it is login-node driven).

## Changes (ordered task list)

### 1. `src/utils/bash/sync_status_email.sh` — add report builders

Two new functions (no new file; every consumer already sources this):

- `array_task_report <job_id>` → prints per-task lines
  `TASKID<TAB>STATE<TAB>ELAPSED<TAB>EXITCODE` for rows matching
  `^<job_id>_[0-9]+[[:space:]]` (via `sacct -j <id> -n --format=JobID,State,Elapsed,ExitCode`),
  numeric-sorted by task number; empty if `sacct` unavailable (guard with
  `command -v sacct`) or returns nothing.
- `array_wall_time <job_id>` → prints the master row Elapsed
  (`sacct -j <id> -X -n --format=Elapsed`), else `n/a`.

Both print to stdout for callers to embed in email bodies. Header doc updated
(one paragraph: duration helpers, per-task filtering rule).

### 2. `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`

- Build `DS_BY_TASK=( $(jq -r 'keys[]' "${DATASETS_JSON_FILE}") )` once; helper
  `task_dataset <n>` = `DS_BY_TASK[n-1]`.
- Function `preprocess_report <job_id>` → renders the "Per-task report" +
  "Array wall time" block (task → dataset, state, elapsed, exit code).
- Replace the raw sacct dump in the gate-failure email with
  `preprocess_report`; append it to the success email; keep subject lines.
- `--sync-only` works unchanged (DATASET_NAMES/jq list rebuilt before the
  branch; sacct + mapping available).

### 3. `src/4_cell_type_annotation/2_submit_hpc_array.sh`

- Function `annotation_report <job_id>`: reads `array_task_report`, maps each
  task → dataset from `chunks_manifest.txt` line <n> (field 1; guard: manifest
  missing → fall back to task-id-only lines), aggregates per dataset
  (COMPLETED/FAILED counts + summed Elapsed), collects failed chunk task ids
  (bounded 20), appends array wall time.
- Use in gate-failure + success emails (replaces raw sacct dump on failure).
- `--sync-only` works (manifest persists on scratch).

### 4. `src/5_run_benchmark_methods/benchmark_submit_common.sh`

- In `benchmark_wait_for_array`: capture per-job report at the gate
  (task → `DATASET_NAMES[i-1]`); failure email body = mapped per-task report
  + array wall time (replaces raw sacct dump); after the success gate, store
  `JOB_REPORTS[<job_id>]="<label>|array wall"` for the final email.
- In `benchmark_merge_sync_cleanup`: append a "Job durations" block to the
  success email (per submitted array: label, job id, array wall time) and to
  the NAS-unreachable email (durations from the gated jobs).
- Signature unchanged; the 3 benchmark submitters need no edits.

### 5. `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — add emails

- Source `sync_status_email.sh`.
- Add `--mail-user="${USER_EMAIL}"` to both sbatch calls (cap step + loop
  steps) for consistency with the other submitters (step scripts only carry
  `#SBATCH --mail-type=END,FAIL`; recipient would otherwise be the cluster
  default user).
- After the per-job sacct gate (and in `--sync-only` mode): collect
  `JobID,JobName,State,Elapsed,ExitCode` for all jobs; send ONE email:
  - success: "ECODA: stage-2 steps COMPLETED (N steps) — no NAS sync by
    design" + per-step lines (name, state, elapsed);
  - failure: "ECODA: stage-2 steps FAILED (k/N)" + per-step lines + failed
    step names.
- Update the stale header comment ("no status email is sent" → sends
  completion report; no NAS sync by design).

### 6. `src/1_stage_data/1_stage_data.sh` — add completion email

- Source `sync_status_email.sh`; track `SECONDS`, count staged files and
  missing-source warnings (replace `echo "  -> [WARNING] ..."` with
  accumulation into `WARNINGS[]`).
- Final email: "ECODA: data staging completed (N files staged, M warnings,
  duration HH:MM:SS)" listing any missing sources; also a `--ds_name`-scoped
  subject when single-dataset mode. Login-node script → same best-effort mail
  CLI (skipped silently if absent).

### 7. `src/4_cell_type_annotation/3_submit_merge.sh` — add durations

- Time each per-view `srun` merge with `date +%s`; append "Merge duration:
  HH:MM:SS" to the per-event success email; time each `merge_one_ds` call in
  default-all mode and append per-dataset durations to the summary email
  (OK and FAILED lists unchanged).
- No sacct needed (srun runs are tracked by wall clock here).

### 8. Documentation

- `docs/ARCHITECTURE.md`: update the submit-script rows for preprocess,
  annotation array, merge, `benchmark_submit_common.sh`, stage-2 dispatcher,
  and the NAS↔scratch flow paragraph — one clause each: emails now include
  per-task/per-dataset (or per-method) status + durations (array wall time).
  Add a row note for `1_stage_data.sh` (completion email).
- `AGENTS.md`: extend the sync-status email bullet (line ~120) with
  "emails include per-dataset/method completion + job durations".
- README.md: no change (kept short by convention).

## Edge cases

- `--sync-only` with purged/unknown id → existing empty-sacct fail-closed
  path; durations report `n/a`; no crash.
- Annotation manifest missing on `--sync-only` after scratch wipe → fall back
  to task-id-only lines (dataset column blank).
- Arrays with `.batch`/`.extern`/`.0` step rows → filtered out by the
  `^<job>_[0-9]+[[:space:]]` regex (also excludes the bare master row, which
  is used for array wall time only).
- `sacct` unavailable (local/dev) → helpers print empty/`n/a`; scripts keep
  their current behavior.
- No mail binary → existing silent skip (notify_sync_status returns 0).

## Validation

1. `bash -n` on all changed bash scripts.
2. Local helper test: stub `sacct` on PATH (prints fixed rows) + a tiny driver
   sourcing `sync_status_email.sh` → `array_task_report`/`array_wall_time`
   output verified; `NOTIFY_MAIL_BIN=cat notify_sync_status ...` prints the
   mail with the report block (no sending).
3. HPC smoke (user, debug-cpu): `--ds_name _debug` run of preprocess and
   annotation; `--sync-only` re-run of the same ids → emails contain the
   per-dataset report + durations; a bogus `--sync-only` id → failure email
   with `n/a` durations.

## Out of scope

- Changing Slurm's own `--mail-type=END,FAIL` emails (kept as-is).
- datasets.json; benchmark submitters' code (helpers only); adding email to
  worker scripts (they are covered by the submitter reports).

## Final step (repo convention)

After implementation: move this plan to `.kilo/plans/archive/`, `git add .`,
commit, push (AGENTS.md Task Completion Workflow).
