# Plan: pixi py-cuda13 env switch, script exec bits, HPC setup docs, commit & push

## Context

- HPC (Bamboo) had no pixi binary (no module); fixed by user-space install (`~/.pixi/bin/pixi`).
- HPC clone had "unstaged changes" blocking `git pull`; root cause: all `src/*.sh` are committed as `100644` while AGENTS.md instructs `chmod +x` on HPC → mode-change noise. Fixed on HPC with `git config core.fileMode false`; `git pull` succeeded.
- User then ran on the Mac: `find src -name "*.sh" -exec git update-index --chmod=+x {} \;` (staged, NOT committed).
- User wants: (a) everything works now and in the future, (b) docs updated, (c) commit & push, (d) point all pixi interpreter pointers to the `py-cuda13` env, (e) validate.

## Verified facts (current state)

- Mac: `git status` shows `MM` on all 12 `src/*.sh` — index staged `100644 → 100755` (mode only, 0 content lines), but **worktree files are still 644** → they must be `chmod +x`-ed on disk, otherwise the Mac repo is dirty again right after commit.
- All interpreter usage routes through `slurm_config.sh` — no other pipeline file hardcodes env paths:
  - `PYTHON_BIN` used (quoted) in: `1.1_submit_combinedpbmc.sh`, `3.1_submit_merge.sh`, `1_prepare_chunks.sh`, `1.1_run_worker.sh`.
  - `PIXI_RSCRIPT` used (unquoted, relies on word-splitting) in: `1.2_submit_joanito.sh`, `1.4_submit_kfoury_lowres_ct.sh`, `2.1_run_worker.sh`, `2_submit_hpc_array.sh`.
  - `RETICULATE_PYTHON` consumed by R workers (reticulate).
- `.Rprofile:11,14` references `.pixi/envs/default` — Mac/IDE-only fallback. `py-cuda13` is linux-64-target-scoped (pixi.toml: `[feature.cuda13.target.linux-64...]`), so it does NOT exist on macOS — **do not change `.Rprofile`**.
- Doc touchpoints: `README.md` lines 61–76 (Installation), 100–110 (HPC execution); `docs/ARCHITECTURE.md` line 152 (Environment propagation) + HPC folder layout `.pixi/` note; `AGENTS.md` line 78 (worker invariants) + HPC general-info section.
- `pixi.lock` (v7, in sync) records envs `default`, `py-cpu`, `py-cuda13` with full linux-64 package sets → installs are lock-replayed (no re-solve).

## Decisions

- **D1 (user):** All HPC interpreter pointers switch to `py-cuda13`; only `src/slurm_config.sh` needs the change (verified all call sites use the env vars).
- **D2:** `.Rprofile` unchanged — Mac fallback; `py-cuda13` does not exist on osx-arm64.
- **D3:** On HPC only the `py-cuda13` env is required after the switch. A stale `default` env from an earlier install is harmless; leave it.
- **D4:** Keep the committed exec bits (`100755`) for all `src/*.sh`; also `chmod +x` the Mac worktree files so the tree stays clean. HPC clone keeps `core.fileMode false` (still needed for any pre-existing chmod'ed files; harmless).
- **D5:** Two logical commits, one push.
- **D6:** HPC R-package setup must run `pixi run -e py-cuda13 setup` so R libs land in the env `PIXI_RSCRIPT` uses.
- **D7:** Validation scope = env/interpreter checks via srun on `debug-cpu`; full pipeline validation with the `_debug` dataset stays deferred (AGENTS.md convention).

## Tasks

### 1. Local repo fixes (Mac)

1. Make worktree modes match the staged index:
   ```bash
   find src -name "*.sh" -exec chmod +x {} \;
   git status --porcelain   # expect clean listing of only the 12 mode changes, no MM
   ```

### 2. Switch interpreter env in `src/slurm_config.sh`

2. Edit lines 30–31 and 39:
   - `export PYTHON_BIN="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"`
   - `export PIXI_RSCRIPT="${HOME}/.pixi/bin/pixi run -e py-cuda13 Rscript --vanilla"`
   - `export RETICULATE_PYTHON="${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python"`
   - Update the adjacent comments (lines 26–29, 33–38) to say the HPC env is `py-cuda13` and that `pixi run` uses `-e py-cuda13` (default env otherwise).

### 3. Documentation updates

3. `README.md`:
   - Installation: correct "creates the default `py-cpu` environment" phrasing (pixi creates the `default` env); for HPC: install the pixi binary in user space (`curl -fsSL https://pixi.sh/install.sh | bash` → `~/.pixi/bin/pixi`, no module exists), then `pixi install --environment py-cuda13` and `pixi run -e py-cuda13 setup`; note HPC uses only `py-cuda13` (all interpreter pointers in `src/slurm_config.sh`).
   - HPC execution block: add a "first-time setup" step (pixi binary + env + setup; heavy installs run via sbatch on `shared-cpu`, NOT the login node — login-node policy).
4. `docs/ARCHITECTURE.md`:
   - HPC folder layout: `.pixi/` line → note envs live at `.pixi/envs/py-cuda13/` on HPC.
   - "Environment propagation" paragraph (line 152): name `py-cuda13` as the env behind `PYTHON_BIN` / `PIXI_RSCRIPT` / `RETICULATE_PYTHON` on HPC.
5. `AGENTS.md`:
   - Worker environment invariants: state HPC interpreters come from the `py-cuda13` pixi env (via `slurm_config.sh`); any `pixi run`/`pixi run setup` on HPC must use `-e py-cuda13`.
   - HPC general info: pixi is a user-space binary at `~/.pixi/bin/pixi` (not a module); HPC clones must set `git config core.fileMode false` (scripts committed `100755`; avoids mode-change noise blocking `git pull`); first-time env setup via the sbatch job in README.

### 4. Commit & push

6. Commit 1 (mode changes only):
   ```bash
   git add src/slurm_config.sh src/*/1_*.sh src/3_scrnaseq_preprocessing/1.1_run_worker.sh src/4_cell_type_annotation/*.sh
   git commit -m "Make pipeline shell scripts executable"
   ```
   (or `git add -u` after verifying no other files are staged — `git diff --cached --summary` must show only the 12 mode changes).
7. Commit 2: `src/slurm_config.sh` + `README.md` + `docs/ARCHITECTURE.md` + `AGENTS.md`
   → message like `Point HPC interpreter env to py-cuda13; document HPC pixi setup`.
8. `git push`.

### 5. HPC verification (user runs on cluster)

9. `cd ~/ECODA_paper && git pull` — must succeed with `core.fileMode false` set; `git status --porcelain` must be empty.
10. Check pixi binary: `~/.pixi/bin/pixi --version` (if not installed yet, do the user-space install first).
11. Env setup — conditional on whether it already ran:
    - If `~/ECODA_paper/.pixi/envs/py-cuda13/bin/python` does not exist, submit the setup sbatch job (README version), containing:
      ```bash
      export PATH="${HOME}/.pixi/bin:${PATH}"
      cd "${HOME}/ECODA_paper"
      pixi install --environment py-cuda13
      pixi run -e py-cuda13 setup
      ```
      (`--partition=shared-cpu`, `--mem=32G`, long `--time`; workers share `$HOME`).
12. Validate via srun on `debug-cpu`:
    ```bash
    srun --partition=debug-cpu --mem=8G --time=30:00 bash -lc '
      ~/ECODA_paper/.pixi/envs/py-cuda13/bin/python -c "import scanpy, anndata, harmonypy; print(scanpy.__version__)"
      ~/.pixi/bin/pixi run -e py-cuda13 Rscript --vanilla -e "cat(R.version.string); library(Seurat); library(anndataR)"
      ~/ECODA_paper/.pixi/envs/py-cuda13/bin/python -c "import torch; print(torch.__version__)"
    '
    ```
13. Smoke-test one pipeline entry point using the env vars (login node submit only, e.g. `bash -lc 'source src/slurm_config.sh; echo "$PYTHON_BIN"; echo "$PIXI_RSCRIPT"'` via srun) to confirm paths resolve. Full pipeline runs remain deferred to the `_debug` validation (out of scope).

### 6. Future-proofing

14. Grep guard: no `envs/default` left in `src/` or docs (`.Rprofile` Mac fallback is the only allowed exception — documented).
15. Conventions now recorded in AGENTS.md: use `PYTHON_BIN`/`PIXI_RSCRIPT` in any new script; `pixi run -e py-cuda13` on HPC; `core.fileMode false` on HPC clones; setup order = pull → env install → run.

## Risks / notes

- Missing Mac worktree `chmod +x` (Task 1) would leave the Mac dirty after commit — do it first.
- If the HPC setup job previously ran without `-e`, R packages are in the now-orphaned `default` env; re-running `pixi run -e py-cuda13 setup` reinstalls them into `py-cuda13` (per-package idempotent). Leave the old env in place.
- `py-cuda13` pulls several GB (pytorch-cuda, jax cuda13); check home quota on HPC (`quota` / `lquota`). CPU-only pipeline stages are unaffected (torch/cuda libs only load when used).
- `pixi run -e py-cuda13 ...` fails with a clear error until the env is installed — pipeline scripts must run only after Task 11.
- Per AGENTS.md, no pipeline scripts are executed for validation here; the interpreter/env checks (Task 12) are the validation boundary.

## Out of scope

- Full pipeline validation on the `_debug` dataset (deferred by convention).
- Changing `.Rprofile` (Mac-only; `py-cuda13` is linux-only).
- Deleting the stale `default` env on HPC.
