# Disconnect-proof submit scripts: --sync-only resume + best-effort sync-status email

## Context and decisions

- **NAS is login-node-only** (user-verified): compute nodes cannot reach `/srv/smednas515.unige.ch`, so the fail-closed sacct gate + rsync tails stay on the login node; the dependent-sync-job (fire-and-forget) design is infeasible.
- **User constraints**: minimal complexity; no tmux; no new tool dependencies, email accounts, or credentials.
- **Decisions**:
  1. Add `--sync-only <job-id>` resume mode to the 6 blocking submitters. Default (blocking) behavior unchanged — when the SSH session stays open, nothing is different.
  2. Add a best-effort sync-status email sent by the tails via the cluster's existing MTA (`mailx`/`mail`/`sendmail`, guarded, silently skipped if absent). Slurm cannot send this email (the sync runs on the login node after the job, not as a Slurm job), so the scripts send it directly — the same MTA that already delivers the existing `--mail-user="${USER_EMAIL}"` job emails. No new dependencies.
  3. Drop afterok and tmux from the previous plan.
- Slurm already emails job END/FAIL (`--mail-user` in every submitter, `--mail-type=END,FAIL` in workers) — the new script email complements it with the "synced to NAS?" answer.

## Out of scope

- NAS-on-compute / DTN routes; `--dependency=afterok` ordering; tmux/nohup; annotation-array → `3_submit_merge.sh` chaining; `--sync-only` for `3_submit_merge.sh` (re-running it is already safe and idempotent); `--parsable` job-id consistency.

## Changes (ordered task list)

### 1. New helper `src/utils/bash/sync_status_email.sh`
- Defines `notify_sync_status <subject> <body>`:
  - Binary discovery: `NOTIFY_MAIL_BIN` env override (test hook), else first available of `mailx`, `mail`, `sendmail` via `command -v`; if none → print "no mail binary; skipping email" and return 0.
  - `sendmail`: pipe `To: ${USER_EMAIL}` + Subject header + body to `sendmail -t`. `mailx`/`mail`: `printf '%s\n' "${BODY}" | ${BIN} -s "${SUBJECT}" "${USER_EMAIL}"`.
  - Sourced (after `slurm_config.sh`) only by the login-node scripts that need it; uses `USER_EMAIL`.

### 2. `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`
- Parse `--sync-only <JOB_ID>` (error on missing id; error on combination with `--force`; update header usage comment).
- When set: skip submission block; set `ARRAY_JOB_ID` from the flag; run the existing tail unchanged.
- Tail: on sync success → `notify_sync_status` "preprocess synced N datasets to NAS"; on gate failure / empty sacct / NAS unreachable → notify "NOT synced — <reason>" (with sacct excerpt on gate failure) before `exit 1`.

### 3. `src/4_cell_type_annotation/2_submit_hpc_array.sh`
- Same as #2; in `--sync-only` mode also skip the scGate DB cache `srun` block and the chunk-manifest build.

### 4. `src/5_run_benchmark_methods/benchmark_submit_common.sh`
- `benchmark_wait_for_array`: on gate failure → notify "NOT synced — task states" before `exit 1`.
- `benchmark_merge_sync_cleanup`: on NAS-unreachable → notify before `exit 1`; on rsync success → notify "synced benchmark/ to NAS".
- Header doc: note the `--sync-only` usage (helpers reused as-is).

### 5. Benchmark submitters (3 files) — add `--sync-only <id1,id2,...>`
- `run_python_sample_embedding_methods/1_submit_hpc_array.sh`: skip submission loop; gate each provided id via `benchmark_wait_for_array`; then `benchmark_merge_sync_cleanup "${METHODS[@]}"`. Header update.
- `run_r_sample_embedding_methods/1_submit_hpc_array.sh`: same (no afterok change).
- `run_transformation_zeroimp_analysis/1_submit_hpc_array.sh`: same with `--analysis` labels.
- Flag rules identical to #2 (comma-split ids; reject `--force` combination; repeat original `--ds_name`/`--methods`/`--analysis` flags so merge labels and sync scoping match).

### 6. `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`
- Add `--sync-only <id1,id2,...>`: skip the submission loop; run the existing wait + `-X` sacct gate + summary for the provided ids. No sync here → no email. Header update.

### 7. `src/4_cell_type_annotation/3_submit_merge.sh`
- Add `notify_sync_status` on sync success and before the NAS-unreachable/coverage-gate failure exits (cheap consistency; script stays interactive and re-runnable). No `--sync-only`.

### 8. Documentation
- **README.md** — in `### Workflow`: short subsection "SSH disconnects are safe": tails are login-bound by design (NAS login-node-only); if the connection drops the pipeline keeps running, and on reconnect you re-run the same command with `--sync-only <job-id>` to re-check + sync without re-running compute. You receive Slurm job emails plus a script email "synced to NAS" / "NOT synced — reason" (best-effort).
- **docs/ARCHITECTURE.md** — update the submit-script rows (preprocess, annotation, 3 benchmark submitters, stage-2 dispatcher, merge) with `--sync-only` and the status email; one sentence in the NAS↔scratch flow paragraphs.
- **AGENTS.md** — one bullet under HPC general information: blocking tails are login-bound; disconnect recovery = re-run with `--sync-only <job-id>` (same flags), never a bare re-submit; sync-status emails are best-effort via the login node's mail CLI.

## Edge cases

- `--sync-only` with unknown/purged job id → empty sacct → fail-closed `exit 1`, no sync, failure email sent.
- `--sync-only` used while the array is still running → `scontrol wait` blocks until completion, then gates (safe).
- No mail binary on the login node → email skipped with a stdout note; `--sync-only` output remains the source of truth.
- Email body includes: script + dataset(s)/method(s) + job id(s) + outcome (+ sacct excerpt on gate failure).
- Overlapping submissions unaffected (benchmark manifests already carry the submit `$$` PID).
- Existing conventions: source `slurm_config.sh`, `cd "${PROJECT_ROOT}"`; flags on the sbatch command line (no `#SBATCH` expansion).

## Validation

1. `bash -n` on all changed bash scripts (local, no HPC needed).
2. Local helper test: `NOTIFY_MAIL_BIN=cat` → `notify_sync_status "test" "body"` prints the mail to stdout (verifies invocation without sending).
3. HPC smoke tests (user runs, debug-cpu):
   - `--sync-only` negative: bogus job id → fail-closed `exit 1`, no sync, failure email.
   - `--sync-only` positive: completed `_debug` run (e.g. preprocess `_debug` or benchmark `--ds_name _debug`) → gate passes, rsync runs, success email.
4. End-to-end `_debug` pipeline runs per the usual workflow.
