# Centralize jq module loading in slurm_config.sh

## Context & root cause

`./src/4_cell_type_annotation/1_prepare_chunks.sh test _debug` failed on the Bamboo login node:

```
Lmod has detected the following error: These module(s) or extension(s) exist but cannot
be loaded as requested: "jq/1.6"
```

Lmod on UNIGE clusters uses a **hierarchical module tree**: `jq/1.6` exists but is hidden
behind its toolchain prerequisite. Verified on bamboo (2026-08-07): `module spider jq` shows

```
You will need to load all module(s) on any one of the lines below before the "jq/1.6" module is available to load.
  GCCcore/12.2.0
```

The toolchain match is **exact** — loading `GCCcore/12.3.0` does NOT expose `jq/1.6`.
Fix verified interactively: `module load GCCcore/12.2.0` then `module load jq/1.6` works.

Repo state today (all committed, working tree clean):
- 3 scripts already load `GCCcore/12.2.0` + `jq/1.6` (bare, no guard): `src/1_stage_data/1_stage_data.sh`, `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`
- 2 scripts load bare `jq/1.6` only (the ones that failed): `src/4_cell_type_annotation/1_prepare_chunks.sh`, `src/4_cell_type_annotation/2_submit_hpc_array.sh`
- 2 scripts load guarded `jq/1.6` only: `src/4_cell_type_annotation/2.1_run_worker.sh`, `src/4_cell_type_annotation/3_submit_merge.sh`
- 1 script loads bare `jq/1.6` only: `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`

All 8 jq-using scripts source `src/slurm_config.sh` (login-node scripts and sbatch
workers alike — module loads do not propagate through `sbatch`, which is exactly why the
central config is the right home for the load). This also fixes the workers
(`1.1_run_worker.sh`, `2.1_run_worker.sh`), which had the same latent failure on nodes
with a bare environment.

## Design decision

**Centralize guarded module loads in `src/slurm_config.sh`; delete all per-script
`module load` lines (jq-related); keep the existing `command -v jq` guards as the
fail-closed backstop.**

- Rationale: single source of truth for the toolchain env (the file already owns
  `PYTHON_BIN`/`PIXI_RSCRIPT`/PATH prepend); removes the current 3-pattern inconsistency
  (bare pair / bare jq / guarded jq).
- `|| true` + `>/dev/null 2>&1` is mandatory: a failing module load must never abort a
  `set -euo pipefail` script nor print Lmod noise; scripts then fail closed with their
  own clear "jq not available" message.
- macOS-safe: `module` does not exist on macOS, the guarded load silently no-ops
  (`command not found` → 127 → `|| true`), so `slurm_config.sh` remains sourceable
  anywhere.
- Rejected alternative: finish the per-script `GCCcore/12.2.0` + `jq/1.6` pattern in the
  5 remaining files — 8 duplicated pairs to maintain, no single place to bump when the
  module tree changes.
- Deferred (documented, not implemented): adding `jq` to the py-cuda13 pixi env
  (`pixi.toml`) would remove the module dependency entirely (env bin is already on PATH
  everywhere); module lines would become harmless fallbacks. Needs pixi.toml change +
  HPC `pixi install --environment py-cuda13`.

## File changes

### 1. `src/slurm_config.sh` — add jq section

Insert after the PATH-prepend block (after line 39, before the reticulate block), as a
new section:

```bash
# --- jq (JSON parsing) ---
# Lmod's hierarchical tree on UNIGE clusters hides jq/1.6 behind its toolchain
# prerequisite GCCcore/12.2.0 (verified via `module spider jq` on bamboo, 2026-08-07).
# Both loads are guarded: a failing module load must never abort a `set -e` script
# nor print noise — the `command -v jq` guards in each consumer script are the
# fail-closed backstop. If the module tree updates (jq 1.6 no longer builds on newer
# GCCcore; EasyBuild pairs jq/1.7.1-1.8.1 with GCCcore/13.x and jq/1.8.1 with
# GCCcore/14.x), re-run `module spider jq` and bump both lines — or add jq to the
# py-cuda13 pixi env to drop the module dependency entirely.
module load GCCcore/12.2.0 >/dev/null 2>&1 || true
module load jq/1.6 >/dev/null 2>&1 || true
```

Note: `set -u` compatibility — the added lines reference no unset vars. `slurm_config.sh`
is sourced after `set -euo pipefail` in every consumer; the `|| true` guards make this safe.

### 2. Remove per-script module loads (7 files)

For each file below, delete the exact `module load` line(s) listed; keep/insert the
`command -v jq` guard where noted. No other logic changes.

| File | Lines to delete | Guard status |
|---|---|---|
| `src/1_stage_data/1_stage_data.sh` | 22-23 (`GCCcore/12.2.0` + `jq/1.6`) | **none today — add** guard after `cd "${PROJECT_ROOT}"` (line 20) |
| `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` | 7-8 (`GCCcore/12.2.0` + `jq/1.6`) | **none today — add** guard after `cd "${PROJECT_ROOT}"` (line 5) |
| `src/3_scrnaseq_preprocessing/1.1_run_worker.sh` | 24-25 (`GCCcore/12.2.0` + `jq/1.6`) | keep existing guard (lines 26-29) |
| `src/4_cell_type_annotation/1_prepare_chunks.sh` | 35 (`jq/1.6`) | keep existing guard (lines 36-39) |
| `src/4_cell_type_annotation/2_submit_hpc_array.sh` | 15 (`jq/1.6`) | keep existing guard (lines 16-19) |
| `src/4_cell_type_annotation/2.1_run_worker.sh` | 22-23 (comment + guarded `jq/1.6`) | keep existing guard (lines 24-27) |
| `src/4_cell_type_annotation/3_submit_merge.sh` | 34 (guarded `jq/1.6`) | **none today — add** guard before the `if ! jq -e` check (line 35) |
| `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh` | 34 (`jq/1.6`) | **none today — add** guard after `cd "${PROJECT_ROOT}"` (line 32) |

Guard snippet (match existing style in `1_prepare_chunks.sh:36-39`):

```bash
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq not available; cannot read ${DATASETS_JSON_FILE}."
  exit 1
fi
```

For `2.1_run_worker.sh`, also delete the now-obsolete comment on line 22 ("jq module
loads in the submit script do not propagate to workers; load here.") — the load now
comes from the centrally sourced `slurm_config.sh`.

Adding guards to the 4 files that lack them is part of this change: the design's
backstop only works where it exists. Without a guard, a future module-tree update would
surface as a misleading "'<ds>' is not a dataset in datasets.json" error instead of a
clear jq-missing message.

### 3. Documentation

- `AGENTS.md` line 102: replace "(via jq, `module load jq/1.6` on the worker)" with a
  reference to the centralized guarded jq load in `slurm_config.sh` (GCCcore/12.2.0 +
  jq/1.6, `|| true`; consumers' `command -v jq` guards fail closed). Keep the sentence
  structure intact.
- `docs/ARCHITECTURE.md`:
  - Line 143 (`2.1_run_worker.sh` row): replace "loads jq (module loads do not propagate
    from the submit script)" with something like "gets jq from `slurm_config.sh`
    (sourced in-worker; guarded `module load GCCcore/12.2.0` + `jq/1.6` — module loads
    do not propagate through sbatch, hence the central location)".
  - Line 154 (Environment propagation bullet): update "(via jq, loaded on the worker —
    module loads do not propagate either)" to reference the centralized load, and add a
    short clause in the bullet noting `slurm_config.sh` now also performs guarded module
    loads for jq (toolchain + version pinned, Lmod hierarchical prerequisite).
- `README.md`, `TODO.md`: no jq/module-load mentions — no changes.

## Explicitly out of scope (do not touch)

- `src/2_dataset_specific_preprocessing/1.1_submit_combinedpbmc.sh` line 22
  (`module load GCCcore/12.2.0`): this load is for rpy2/R interop (rds→h5ad conversion,
  see docstring in `1.1.1_create_combinedpbmc_dataset.py:41`), NOT jq. Keep it — it must
  stay unguarded-explicit since the worker needs the R toolchain to actually work, and
  prior plan `1786048138782` decided to keep it. (Post-change it becomes idempotent with
  the slurm_config.sh load; harmless.)
- `src/2_dataset_specific_preprocessing/1.2_submit_joanito.sh`,
  `1.3_submit_kfoury_lowres_ct.sh`: no jq usage, no module loads — untouched.
- Adding `jq` to `pixi.toml` (py-cuda13 env): future hardening, documented in the
  slurm_config.sh comment, not part of this change.

## Validation

Per AGENTS.md, do not run pipeline scripts as part of implementation validation unless
the user asks; HPC acceptance runs are performed by the user. Implementer-verifiable
steps:

1. `bash -n` on all 9 edited files (8 scripts + `slurm_config.sh`).
2. macOS smoke (proves the centralized load cannot break local sourcing):
   `bash -c 'source src/slurm_config.sh && echo OK'` — must print `OK` with no error
   (no `module` command on macOS → silently no-ops via `|| true`).
3. Confirm no remaining `module load` lines in the 8 jq scripts:
   `git grep -n "module load" -- src` should return only
   `src/2_dataset_specific_preprocessing/1.1_submit_combinedpbmc.sh:22`.

User acceptance on HPC (fresh SSH session, NO manual `module load` — this is the
regression scenario):

1. `./src/4_cell_type_annotation/1_prepare_chunks.sh test _debug` → must print
   "Datasets to process: _debug" with no Lmod error (chunks re-created in test mode,
   harmless).
2. `./src/1_stage_data/1_stage_data.sh --ds_name _debug` → clean no-op smoke.
3. Optional: `module purge && bash -n <any edited script>` on the login node.
4. Optionally verify a worker path: submit the annotation array for `_debug` with
   `2_submit_hpc_array.sh _debug` (or defer to the regular pipeline run).

## Risks

- **Module-tree update** (jq/1.6 or GCCcore/12.2.0 removed): guarded loads no-op; scripts
  fail closed with the explicit "jq not available" message. Mitigation: comment in
  slurm_config.sh documents the `module spider jq` re-check procedure; pixi-jq is the
  permanent fix (deferred).
- **GCCcore auto-swap**: if a user/session already has a different GCCcore (e.g. 12.3.0)
  loaded, Lmod auto-swaps to 12.2.0 ("reloaded with a version change", observed
  behavior). Benign — no script depends on a specific GCCcore version.
- **Silent jq absence**: guarded load hides failures by design; the added
  `command -v jq` guards in all 8 consumer scripts ensure the failure is loud and
  unambiguous.
