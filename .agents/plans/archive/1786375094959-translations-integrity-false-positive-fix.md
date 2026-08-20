# Plan: fix `translations` false positive in R-library integrity check + failure diagnostics in env-install scripts

## Context / problem (user-reported, root cause confirmed locally)

Failed HPC login-node install (`rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup`, old approach, in tmux) ended with:

```
Error: R library integrity check FAILED: 1 package(s) missing Meta/package.rds (corrupt/partial install): translations. ...
Execution halted
```

**Root cause (verified on this Mac, read-only checks):** `translations` is NOT a failed package. It is R's own message-catalog component shipped inside `lib/R/library/` by conda r-base (`DESCRIPTION` says `Package: translations`, plus locale dirs `ar/ de/ fr/ ...`; no `Meta/package.rds`). It is copied at R build time (`make install`), never installed via `R CMD INSTALL`, so it never has a `Meta/` dir. Evidence:

- `.pixi/envs/{default,py-cpu,py-tools,r-tools}/lib/R/library/translations/` all contain `DESCRIPTION` and no `Meta/package.rds`.
- Simulating the exact `check_lib_integrity()` from pixi.toml over `.pixi/envs/default/lib/R/library` (422 dirs): exactly `translations` is flagged.

**Why now / why the env is fine:** the integrity check is the LAST step of `pixi run setup`. The user's run was the first full clean rebuild since the check was added (commit ec85cb2); every remotes install had completed — only the final verification false-positived. **The HPC env is almost certainly healthy: NO `rm -rf` is needed for recovery.** A re-run of the (fixed) setup is fast because remotes skips SHA-pinned already-installed packages.

**Answer to "can we check before re-running":** yes — done locally (above). Additionally, this plan makes the install scripts self-diagnosing so any future real failure is identifiable on the worker node without a second guess.

## Decisions (locked)

- **Fix in pixi.toml `[tasks.setup]`**: skip the `translations` directory by name in `check_lib_integrity()` (with a comment). Name-based is precise; do NOT switch to a `Built:`-field heuristic — an interrupted `R CMD INSTALL` can leave DESCRIPTION without `Built:`, which would mask real partial installs.
- **Enrich the failure output**: flagged dirs are reported with their FULL PATH + the `Package:` field from their DESCRIPTION, so a future non-package dir is recognizable and a future corruption is actionable. Keep the existing wipe-and-reinstall repair text.
- **Automatic failure diagnostics** (no opt-in flag; no-op on success) in BOTH entry-point scripts: on `pixi run -e py-cuda13 setup` failure, run a diagnostic Rscript that prints R version, `.libPaths()`, per-library counts, every dir with DESCRIPTION missing `Meta/package.rds` (path + `Package:` field), and the "translations is the only legitimate case" hint. Same wrap added to `refresh_env.sh` (login-node path) for consistency.
- **No pixi.toml dependency change** → no `pixi install` side effects, pixi.lock untouched.
- **Recovery path**: HPC `git pull` then re-run `sbatch src/utils/bash/setup_env_sbatch.sh` (it runs `pixi install` — idempotent/fast — then setup, which skips installed packages). Do NOT `rm -rf` again.

## Tasks (ordered)

### 1. pixi.toml — fix `check_lib_integrity` + enriched failure message
Replace lines ~190–204 (current block) with:
- In `check_lib_integrity()`: after the DESCRIPTION check, `if (basename(p) == 'translations') next` with a comment (R message-catalog component, copied at R build time, never installed via R CMD INSTALL, no Meta/package.rds by design).
- Collect FULL dir paths (`bad <- c(bad, p)`, not `basename(p)`).
- On failure: before `stop()`, per-dir detail line via `cat()`:
  `pkg_line <- grep('^Package:', readLines(file.path(p, 'DESCRIPTION'), n = 20), value = TRUE)`
  `cat('Missing Meta/package.rds:', p, paste(pkg_line, collapse = '; '))`
  (`cat()` appends a newline per call — no `'\n'` escapes needed.)
- `stop(sprintf(...))` keeps the repair text; mention that the only legitimate case is `.../translations` (R message catalogs).
- **Hard constraints (block lives inside TOML `"""` AND inside a shell `Rscript -e "..."` wrapper):** single quotes only, no `"`, no `$`, no backticks, no `\n` inside R string literals. (Validated every time via the tomllib/R-parse check below.)

### 2. src/utils/bash/setup_env_sbatch.sh — failure diagnostics
- Replace the bare step-3 call with:
  ```bash
  if ! "${PIXI_BIN}" run -e py-cuda13 setup; then
    echo "ERROR: pixi run setup failed — running env diagnostics..." >&2
    run_env_diagnostics
    exit 1
  fi
  ```
  (keep the success `echo ... elapsed` after the `if`).
- Add `run_env_diagnostics()`: a `"${PIXI_BIN}" run -e py-cuda13 Rscript --vanilla -e '...'` single-quoted block (R uses DOUBLE quotes there — opposite of the TOML constraint) that prints:
  - `R.version.string` and each `.libPaths()` entry
  - per library: total DESCRIPTION-dir count and every dir missing `Meta/package.rds` (full path + `Package:` field via `grep('^Package:', ...)`)
  - hint: only legitimate case is `.../translations` (R message catalogs, false positive → update pixi.toml); otherwise the env is genuinely corrupt → wipe-and-reinstall repair.

### 3. src/utils/bash/refresh_env.sh — same failure-diagnostics wrap
- Same `if ! "${PIXI_BIN}" run -e py-cuda13 setup; then run_env_diagnostics; exit 1; fi` pattern around its step `[2/3]`, reusing the identical diagnostic function (keeps login-node and worker-node failure output consistent).

### 4. AGENTS.md — one-line doc touch
- In the env-refresh bullet, update the integrity-check description to "(Meta/package.rds presence, skipping R's own `translations` component, + critical packages load)".

### 5. Validation (no pipeline runs — per AGENTS.md)
- `bash -n` on both scripts.
- tomllib parse of pixi.toml; extract the setup `cmd`, strip the `Rscript -e "` wrapper; assert no `"`/`$`/backtick and no backslash-n in the raw block; `Rscript --vanilla -e 'parse(...)'` the extracted R code.
- Simulate the FIXED `check_lib_integrity()` on `.pixi/envs/default/lib/R/library` → expect **zero** flagged dirs; re-run the UNFIXED version → expect exactly `translations` (documents the false positive).
- R-parse the diagnostic R snippet (extract from the script or replicate in a temp file).
- HPC run (user-driven, when no other jobs are active): `git pull` → `sbatch src/utils/bash/setup_env_sbatch.sh` — NO `rm -rf`; expect fast completion (pixi install idempotent, setup skips installed SHA-pinned packages, "R library integrity check OK" then "All critical packages load OK", exit COMPLETED). Alternative login-node path: `tmux new -s env-refresh` → `src/utils/bash/refresh_env.sh`.
- Follow the repo's Task Completion Workflow: archive this plan to `.kilo/plans/archive/`, `git add .`, commit (repo message style), push; then HPC `git pull` and the sbatch run.

## Risks / notes

- Name-based `translations` skip is safe: it is the only non-package dir in a conda r-base library (verified across 4 local envs, 422 dirs). The enriched error output (full path + `Package:` field) makes any future exception identifiable.
- TOML-block constraints are the main footgun (no `"`, `$`, backticks, `'\n'` escapes); the sbatch/refresh diagnostics have the OPPOSITE constraint (single-quoted shell block → double quotes in R). Validated per Task 5.
- The failed previous run left the env fully installed — re-running setup WITHOUT `rm -rf` is both safe and fast; mention this in the commit message / README note so the user does not re-wipe.

## Out of scope

- HPC-side inspection before the re-run (local evidence is conclusive; the re-run self-verifies and the diagnostics cover any real failure).
- Changing the check strategy beyond the `translations` skip; `pixi.lock` regeneration; ccache.
