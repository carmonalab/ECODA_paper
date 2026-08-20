# Fix USER_EMAIL, switch scanpy leiden to igraph, harden sacct monitor

Three independent fixes, confirmed with the user:

1. **USER_EMAIL**: `slurm_config.sh` currently hardcodes `USER_EMAIL="${USER}@unige.ch"` which silently produces a non-deliverable address (`halterc@unige.ch`). The user's real address must NOT appear anywhere in the repo — solve via an env var set in the user's HPC shell profile.
2. **igraph**: `sc.tl.leiden` emits a FutureWarning because scanpy's igraph flavor is not used and the Python `igraph` package is not installed (`pixi.toml:35` `r-igraph` is the R package only; `pixi.lock` contains no Python igraph, only `leidenalg`). Switch to `flavor="igraph"` — no mass re-run of already-preprocessed datasets (`leiden_res_*` has no consumers anywhere; verified by grep + archived plan note).
3. **sacct monitor**: job 4294752's NAS sync was skipped because the bounded sacct poll (60×5s = 5 min) in `1_submit_hpc_array.sh` expired while SlurmDBD still reported RUNNING (accounting lag). Harden the poll in all 4 submit scripts.

## Tasks (ordered)

### 1. `src/slurm_config.sh` — USER_EMAIL env-respecting default (lines 116-117)

Replace the unconditional assignment with a fallback that respects a pre-set env var (sbatch propagates the login shell's environment to the submit script; workers inherit it too):

```bash
# --- User Info ---
# USER_EMAIL is the recipient for Slurm --mail-user and sync-status emails
# (notify_sync_status). It must be set by the user in their HPC shell profile
# (e.g. ~/.bashrc) — personal addresses must NOT be hardcoded in this repo.
# The default below is a best-effort guess and may be non-deliverable.
if [[ -z "${USER_EMAIL:-}" ]]; then
  export USER_EMAIL="${USER}@unige.ch"
  echo "WARNING: USER_EMAIL unset — falling back to ${USER_EMAIL}. Set USER_EMAIL in your HPC ~/.bashrc to receive Slurm + sync-status emails." >&2
else
  export USER_EMAIL
fi
```

(`${USER_EMAIL:-}` is safe under callers' `set -u`.)

### 2. `README.md` — one-line note

In the HPC setup/usage section: state that `USER_EMAIL` must be exported in the HPC shell profile (`~/.bashrc`) to receive `--mail-user` and sync-status emails. Do NOT put any real address in the repo.

### 3. `pixi.toml` — add Python igraph

In the CORE PYTHON STACK block (after `scanpy = "==1.12.2"`, line 57):

```toml
# Python igraph backend for sc.tl.leiden(flavor="igraph"); r-igraph above is
# the R package (unrelated)
python-igraph = ">=0.11,<1"
```

Then run `pixi install` locally (updates `pixi.lock` — commit it; installs python-igraph into the local `py-cpu` env).

### 4. `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` — leiden flavor (line 104)

In `run_clustering`, add `flavor="igraph"`, `n_iterations=2`, `directed=False` to the `sc.tl.leiden` call (params per the FutureWarning; neighbors graph is undirected):

```python
sc.tl.leiden(
    adata, resolution=r, key_added=f"leiden_res_{r}_{key_suffix}",
    neighbors_key=neighbors_key,
    flavor="igraph", n_iterations=2, directed=False,
)
```

### 5. sacct monitor hardening — 4 submit scripts

Same pattern in all 4 files: extend the first bounded poll cap from 60 to 180 iterations (15 min at 5s); if it exhausts with non-terminal states, print an explicit warning and run a 60-iteration (5 min) grace window before falling through to the existing fail-closed gate (the gate itself — `grep -qvE '^ *COMPLETED *$'` — is UNCHANGED).

- `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` (poll loop, lines 146-155)
- `src/4_cell_type_annotation/2_submit_hpc_array.sh` (poll loop, lines 165-174)
- `src/5_run_benchmark_methods/benchmark_submit_common.sh` (`benchmark_wait_for_array`, ~lines 92-105)
- `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` (per-job poll loop, lines 97-111; same pattern around its per-job `state` variable)

Snippet for the array submitters (adapt var names for 2_dataset_specific's per-job `state` loop):

```bash
TAIL_ITER=0
while (( TAIL_ITER < 180 )); do  # max 15 min at 5s
    STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
    if [[ -n "${STATES//[[:space:]]/}" ]] \
       && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
        break
    fi
    sleep 5
    TAIL_ITER=$((TAIL_ITER + 1))
done
if (( TAIL_ITER >= 180 )); then
    echo "WARNING: sacct still reporting non-terminal states after 15 min; extending wait by a 5 min grace window..." >&2
    TAIL_ITER=0
    while (( TAIL_ITER < 60 )); do
        STATES="$(sacct -j "${ARRAY_JOB_ID}" --format=State -n 2>/dev/null || true)"
        if [[ -n "${STATES//[[:space:]]/}" ]] \
           && ! grep -qE 'PENDING|RUNNING|REQUEUED|CONFIGURING|SUSPENDED|RESIZING' <<< "${STATES}"; then
            break
        fi
        sleep 5
        TAIL_ITER=$((TAIL_ITER + 1))
    done
fi
```

## User actions on HPC (after commit + pull)

1. Add `export USER_EMAIL="[REDACTED_EMAIL]"` to `~/.bashrc` on the HPC login node (login shells source it; sbatch propagates it).
2. `git pull` in `[REDACTED_PATH]/ECODA_paper`, then `~/.pixi/bin/pixi install -e py-cuda13` (installs allowed on the login node; adds python-igraph to worker env).
3. Resume the Kfoury sync (if not already done): once `sacct -j 4294752 --format=JobID,JobName,State,ExitCode` shows COMPLETED, run
   `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name Kfoury --sync-only 4294752`

## Risks / constraints

- igraph's leiden optimizer may yield slightly different clusters than leidenalg; `leiden_res_*` has no downstream consumers (grep + archived plan confirm), so the paper is unaffected. Already-synced h5ads keep leidenalg labels until re-run (`--force`); acceptable per user decision.
- Do NOT add `--mail-user` to `2_dataset_specific_preprocessing/1_submit_hpc.sh`: its current Slurm-default address (`${USER}@<submit-host>`) is forwarded to the real address by the cluster's alias (proven by the combine_pbmc .eml) and works even when USER_EMAIL falls back. Adding `--mail-user` there would break delivery in the fallback case.
- Extending the poll only affects login-node wait time in pathological DBD-lag cases; the fail-closed gate and all sync semantics are untouched.
- Do not run pipeline scripts for validation (AGENTS.md rule).

## Validation

- `bash -n` on the 5 modified bash scripts.
- `python -m py_compile src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`.
- Locally: `pixi install` then `pixi run -e py-cpu python -c "import igraph; print(igraph.__version__)"`.
- Confirm `git diff pixi.lock` contains only the python-igraph additions (no unrelated churn).

## Completion workflow (AGENTS.md)

After implementation: move this plan to `.kilo/plans/archive/`, `git add .`, commit summarizing the three fixes, push.
