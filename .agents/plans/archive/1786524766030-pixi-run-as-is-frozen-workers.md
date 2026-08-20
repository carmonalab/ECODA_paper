# Harden PIXI_RSCRIPT: `pixi run --as-is` so workers never mutate the env at runtime

## Context / problem

`PIXI_RSCRIPT` (src/slurm_config.sh:34) is the R entry point for every compute-node
R invocation in the pipeline: annotation workers, R benchmark workers,
transformation/zeroimp workers, the srun scGate-DB session, and the dataset-specific
prep sbatch jobs (`1.3_submit_joanito.sh`, `1.4_submit_kfoury_lowres_ct.sh`).

Today it expands to `pixi run -e py-cuda13 Rscript --vanilla`. A plain `pixi run`
"will also update the lock file and install the environment if it is required"
(pixi docs) — i.e. any worker whose installed prefix does not match the lockfile
silently triggers a solve + full env install at runtime. The 2026-08-11 incident
("~12:35 pixi prefix update via `pixi run`", immediately preceding the
stale-view/ENOENT array failures) is exactly this mutation class. The existing
guards (env_refresh.lock, no-active-jobs check) cover the *deliberate* refresh
paths but not implicit worker-side updates.

Goal: a worker can never trigger a lockfile update or an environment (prefix)
install. Env refreshes happen only via the two guarded entry points
(`setup_env_sbatch.sh` / `refresh_env.sh`), which use `${PIXI_BIN}` directly and
are unaffected by this change.

## Verified flag semantics (pixi 0.70.2 help + source, main branch)

Source: `crates/pixi_cli/src/cli_config.rs` (`LockAndInstallConfig::allow_installs()`
= `!as_is && !no_install`), `crates/pixi_cli/src/run.rs` (env-install block gated on
`allow_installs()`, calls `lock_file.prefix(..., UpdateMode::QuickValidate)`),
`crates/pixi_core/src/lock_file/update.rs` (QuickValidate → full `update_prefix`
when the prefix hash mismatches the lockfile).

| Flags | Lockfile update | Out-of-date check | Env (prefix) install |
|---|---|---|---|
| (default) | yes | yes | yes |
| `--frozen` | **no** | no | **yes** |
| `--locked` | no | yes (abort on manifest drift) | **yes** |
| `--no-install` | yes | yes | no |
| `--as-is` (= `--no-install --frozen`) | **no** | **no** | **no** |

Key finding: `--frozen` alone does NOT prevent env installs — it only freezes the
lockfile. The correct flag for "never mutate at runtime" is **`--as-is`**
(user-confirmed decision; validated locally with pixi 0.70.2:
`pixi run --as-is -e py-cuda13 Rscript --vanilla -e 'cat("ok")'` → exit 0).

## Implementation tasks

1. **`src/slurm_config.sh`** — change line 34 to:
   ```bash
   export PIXI_RSCRIPT="${HOME}/.pixi/bin/pixi run --as-is -e py-cuda13 Rscript --vanilla"
   ```
   Update the comment block (lines 30–32) to document `--as-is`: pixi `run` must
   never update the lockfile or install/repair the env from a worker — that is the
   "prefix update via `pixi run`" mutation class behind the 2026-08-11 corrupted
   library/ENOENT failures; env changes only via the guarded
   `setup_env_sbatch.sh`/`refresh_env.sh`. Note explicitly that `--frozen` is NOT
   sufficient (it still installs the prefix when it mismatches the lockfile) and
   that `--as-is` = `--no-install --frozen`.

2. **Comment updates where the command string is quoted** (cosmetic, keep accurate):
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh` line 74
   - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh` line 78
   - `src/4_cell_type_annotation/2.1.1_process_chunk.R` line 4
   - `src/4_cell_type_annotation/2.0_create_scgate_db.R` line 4
   (`pixi run ...` → `pixi run --as-is ...`)

3. **`docs/ARCHITECTURE.md`** — line ~165 ("PIXI_RSCRIPT expands to `pixi run
   -e py-cuda13 Rscript --vanilla`"): add `--as-is` and one clause: workers never
   update the lockfile or install the env at runtime; refreshes only via the
   guarded entry points.

4. **`AGENTS.md`** — in the "Worker environment invariants" bullet (line ~143):
   add that `PIXI_RSCRIPT` runs with `--as-is` (no runtime lockfile/env mutation;
   `--frozen` alone insufficient), and in the env-refresh paragraph (line ~141)
   note that drift now fails as stale-env errors instead of auto-installing.

## Behavior changes / risks

- **Env missing on a worker** (fresh clone, env build skipped): previously pixi
  auto-installed; now the task fails fast at activation with a pixi error. This is
  desired — the workflow mandates building the env first via the guarded script.
- **Env stale vs a new lockfile** (e.g. git pull with a regenerated `pixi.lock`,
  no HPC reinstall): workers run the OLD env; missing new packages surface as the
  familiar `package or namespace load failed` / `No module named` signatures
  (transient retry requeues up to 3×, then fails). Env stays untouched — the
  corruption vector is gone; the "refresh before arrays" workflow becomes
  mechanically enforced.
- **Rattler NFS-cache warning**: with `--as-is` the update check is skipped, so the
  repodata-cache warning should stop appearing in worker logs (cosmetic; not a goal).
- **No impact**: on macOS (slurm_config.sh is HPC-only), on
  `setup_env_sbatch.sh`/`refresh_env.sh` (use `${PIXI_BIN}` directly), on Python
  workers (bare `${PYTHON_BIN}`, no pixi), on rpy2 workers (env `Rscript` via PATH).
  No pixi.toml/pixi.lock changes, no env rebuild needed.

## Validation

1. Local smoke (re-run by implementer): `pixi run --as-is -e py-cuda13 true` and
   `pixi run --as-is -e py-cuda13 Rscript --vanilla -e 'cat("ok")'` → exit 0.
2. HPC login node (non-compute): `~/.pixi/bin/pixi --version` — confirm `--as-is`
   support. If the HPC pixi lacks it (unlikely; present in 0.70.2), fall back to
   `pixi run --no-install --frozen ...` (equivalent, older flags) and keep the
   slurm_config.sh comment in sync.
3. HPC sbatch smoke on `debug-cpu`: a tiny job that sources `slurm_config.sh`,
   echoes `${PIXI_RSCRIPT}`, runs it with `-e 'print(R.version.string)'`, and
   checks the env was untouched — record mtime/hash of
   `.pixi/envs/py-cuda13/conda-meta/history` and `pixi.lock` before/after and
   compare (or `git status --short pixi.lock`).
4. Real-world: next array run (any stage) exercises the new string end-to-end.

## Workflow (per AGENTS.md)

Implement, then move this plan to `.kilo/plans/archive/`, `git add .`, commit
(message e.g. "Harden PIXI_RSCRIPT with pixi run --as-is (no runtime env
mutation)"), and push — unless the user wants to review first.
