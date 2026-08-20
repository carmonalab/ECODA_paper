# Docs update: watchdog 4325218 false-fail + Adams sample-name standardization

Status: Completed / Archived. Sync for array 4325217 completed in previous session; technical recovery notes and sample-name standardization facts documented.

## Context (session ses_feaf38333ffe1VGa7wsuSBdZmI, 2026-08-18)

`./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods composition`
submitted:
- prepare_pseudobulk array 4324976 + soft-gate watchdog 4324977 → **OK**
- composition array 4325217 (11 datasets) + strict watchdog 4325218 → **watchdog FAILED (exit 1)**

The work **actually completed** (Slurm email: `4325217_* COMPLETED, ExitCode [0-0]`, all 11
tasks), but the watchdog 4325218 **failed closed on a stale sacct state**: its per-task
report shows `1 -> Adams RUNNING (0:0)` — a RUNNING state captured while the last task was
still draining / SlurmDBD accounting had not flushed the terminal state. The watchdog's
bounded sacct poll (`benchmark_wait_array_terminal`) timed out at its 20-min cap (15 min +
5 min grace, `benchmark_submit_common.sh:294-318`), then `benchmark_wait_oom_retry`
(`benchmark_submit_common.sh:463-477`) treated the leftover RUNNING row as a
non-COMPLETED, non-OOM task and wrote `STATE=FAIL` → sync to NAS skipped.

Not a real task failure; **not** related to the Adams `-`/`_` sample-name question
(verified: standardization happens once in `1.1.1_preprocess.py`; R workers do not
re-standardize; the composition worker is obs-only and would never see unstandardized
names; and the array summary shows COMPLETED anyway).

## Verified facts (from code, 2026-08-18)

1. **Adams `-`/`_` standardization site**: `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:285-290`
   — before writing each view h5ad, `obs["Sample"]` gets `f"g{s}"` for names starting
   with a digit and `str(s).replace("-", "_")` otherwise. Subsetting (from `subset_vars`)
   runs on ORIGINAL values BEFORE standardization (`1.1.1_preprocess.py:275-278`).
2. **R workers must NOT re-standardize** (would diverge hyphen→underscore from obs names
   for h5ads predating the python change, e.g. Adams → broken bundle label match):
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R:59-62` (comment present)
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:128-131, 181-183` (comments present)
   - Legacy Seurat-only path (notebook): `src/5_run_benchmark_methods/benchmark_pipeline.R:438-441` — only place that still calls `standardize_sample_names()`
   - GloScope worker still standardizes its own rownames: `benchmark_methods_r.R:359` (intentional, GloScope returns its own rownames)
3. **Docs coverage**: ARCHITECTURE.md documents the rationale in the prep-worker row
   (`docs/ARCHITECTURE.md:439`) and the R-method-worker row (`:438`). AGENTS.md does
   NOT mention it.
4. **Recovery for a FAIL watchdog status file**: `benchmark_wait_watchdog`
   (`benchmark_submit_common.sh:661-721`) re-reads the status file; a `STATE=FAIL` file
   always fails closed (`:688-697`) — so `--sync-only 4325218` FAILS again.
   The working recovery is `--sync-only 4325217` (the **array** id): no status file →
   strict all-rows gate (`benchmark_wait_for_array`, `benchmark_submit_common.sh:574-599`)
   against current sacct (all COMPLETED) → sync runs.
   Caveat: `benchmark_wait_for_array` includes the master + batch rows in the
   all-COMPLETED gate (in contrast to the watchdog's per-task gate). The master row for
   a completed array is COMPLETED; if sacct is still laggy it fails closed again and the
   command can simply be re-run (idempotent).

## Scope & decisions

- **Docs-only + one comment**: update AGENTS.md, ARCHITECTURE.md, TODO.md; add one
  clarifying comment in `benchmark_wait_array_terminal` about the stale-RUNNING race
  (root cause of 4325218).
- No pipeline logic changes (no watchdog fix, no resubmission code). The user gets the
  recovery command for HPC.

## Ordered tasks

### 1. AGENTS.md — Adams sample-name standardization note
Add to the "Worker environment invariants" bullet list (or the R Modules section), e.g.:
"Sample-name standardization (`g`-prefix for digit-leading names, `-`→`_` in `Sample`)
happens ONCE in `1.1.1_preprocess.py` (before view h5ad write; subsetting runs on
original values first). R benchmark/prep workers must NOT re-apply
`standardize_sample_names()` — it would diverge (hyphen→underscore) from the obs names
for h5ads predating the python change (e.g. Adams) and break the bundle label match.
Remaining call sites: legacy Seurat path `benchmark_pipeline.R` and GloScope rownames
(`benchmark_methods_r.R`)."

### 2. ARCHITECTURE.md — watchdog false-fail + standardization
- In the `1_submit_hpc_array.sh` row (`docs/ARCHITECTURE.md:436`) or the watchdog row
  (`:441`), add: watchdog gates can fail closed on a stale `RUNNING` sacct row (task
  left scheduler while SlurmDBD still flushes; the bounded 20-min poll can time out).
  Observed 2026-08-18: composition watchdog 4325218 FAILED while array 4325217 was
  fully COMPLETED; sync recovered via `--sync-only <array-id>` (not the watchdog id,
  whose FAIL status file re-fails in `benchmark_wait_watchdog`).
- In `1.1.1_preprocess.py` row (`:98`): mention the sample-name standardization
  (`g`-prefix + `-`→`_`), applied per-view before write, subsetting on original values.

### 3. TODO.md — optional cleanup note
Under "Ideas for later" (or a new note): watchdog stale-RUNNING false-fail — consider a
short post-gate grace + re-read of sacct before writing FAIL (no code change now).

### 4. Code comment (docs-quality, 1 location)
`benchmark_submit_common.sh:294` (the sacct-lag comment block): add one sentence noting
that when the bounded poll times out with a stale RUNNING, the strict gate below can
false-fail a COMPLETED array; recover via `--sync-only <array-id>` (the FAIL watchdog
status file re-fails).

### 5. Verify
- No HPC/pipeline run required (docs + comment only).
- `git grep` to confirm the docs text is consistent.
- Stage + commit + push per the Task Completion Workflow (plan file goes to
  `.kilo/plans/archive/`).

## User commands (run on HPC, from `~/ECODA_paper`, `$HOME/ECODA_paper` = clone)

Verify true terminal states (optional but recommended):
```bash
sacct -j 4325217 --format=JobID,JobName,State,ExitCode,Elapsed -n
```
Recover the sync (repeat original flags; uses the array id, NOT the watchdog id):
```bash
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods composition --sync-only 4325217
```
If the gate still sees laggy sacct rows, re-run the same command (idempotent).

## Risks / open questions
- None material: docs-only; the only judgment call is wording placement (AGENTS.md
  section), which the implementer resolves by following the existing style.
