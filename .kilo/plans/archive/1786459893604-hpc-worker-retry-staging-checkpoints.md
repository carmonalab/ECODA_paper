# HPC worker hardening: self-healing retry, R-lib staging, per-sample checkpoints

## Context

Annotation array `4308096` (re-run of `2_submit_hpc_array.sh`) failed/slow workers:

- **Task 128** (Gongsharma, chunk 28): `cannot open file '.../S4Arrays/R/S4Arrays.rdb': No such file or directory` — BeeGFS stale client-cache view on that worker node (file exists on login node). ABI warnings in other tasks' logs (`SeuratObject` built under R 4.5.1, current R 4.5.2; Matrix 1.7.4 vs 1.7.5) confirm the env was mutated recently, which opens the staleness window (AGENTS.md documents up to ~1h).
- **Task 378** (Pelka, chunk 48): still RUNNING at 47+ min against a **1h wall limit** (`#SBATCH --time=01:00:00` in `2.1_run_worker.sh`) — will very likely die as Slurm TIMEOUT. `.err` shows `scATOMIC attempt 1 with 643 s timeout` then `Error in namedrop(args) : reached CPU time limit` (`namedrop` confirmed from `bbmle`, scATOMIC's mle2 path).

### Test evidence (reticulate/timeout escape) — verified locally 2026-08-11

Reproduced the annotation pattern (`tryCatch` + `TimeoutException`/`error` handlers around `R.utils::withTimeout`) with reticulate python calls in the local pixi env:

- **Timeout fires at an R checkpoint** (e.g. inside `optim`/`mle2` loops): R.utils catches it correctly → `TimeoutException` handled → retry path returns NULL, time limits reset cleanly. The common scATOMIC failure mode is therefore *catchable*.
- **R blocked inside a python call (reticulate) when the budget expires**: the timeout does **not** fire during the python call and, if no R loop/interrupt checkpoint occurs afterwards, may never fire at all (silent overrun past the budget). `setTimeLimit` is only checked at R evaluation points — a python call can overrun the whole budget.
- Conclusion: `withTimeout` per-attempt budgets are **best-effort only** around scATOMIC. The plan therefore enforces the real budget at R level (wall-clock checks between attempts/samples vs `SLURM_TIME_LIMIT`), catches escapes at the **sample** level, and makes any residual overrun cost at most one in-flight sample via per-sample checkpoints. Slurm's wall limit remains the backstop.

## Decisions (user-confirmed)

1. **Rollout scope: all array workers** — annotation (`4_cell_type_annotation`), preprocess (`3_scrnaseq_preprocessing`), and the three benchmark workers.
2. **Wall-time policy (annotation)**: `--time` 1h → **2h**, per-attempt timeout capped at **1800 s**, attempts stop when remaining wall time < timeout + 300 s margin, BLAS/OMP threads pinned to 1 so CPU time ≈ wall time. Thread pinning is **annotation-only** (benchmark workers are hardware-pinned for runtime comparability — pinning would invalidate cross-method runtime comparisons; preprocess is legitimately multi-threaded).
3. Staging: R library **only** (`lib/R/library`), per-task copy to node-local `/scratch`, with a size guard (skip staging + warn if too large). Python site-packages are NOT staged (torch etc. too large; retry covers flakes). Node-shared staging (`/srv/share/users/...`) is a documented follow-up, not part of this plan.
4. Retry: transient-failure signature grep on the Slurm `.err` file + `scontrol requeue <array_job>_<task>` from the RUNNING task, capped by a counter file in scratch (default max 3 retries).

## Implementation tasks

### A. New shared helper: `src/utils/bash/worker_retry.sh`

Sourced by all 5 array workers (after `slurm_config.sh`). Follows repo conventions (bash, `set -euo pipefail`-compatible, no comments unless needed). Provides:

- `TRANSIENT_REQEX` (grep -Ei pattern):
  `No such file or directory|cannot open file|cannot open shared object|package or namespace load failed|there is no package called|cannot open connection|failed to load|No module named`
- `worker_retry_count_file` → `${HPC_SCRATCH_DIR}/_worker_retries/<jobid>_<taskid>.count` where jobid/taskid = `${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}` and `${SLURM_ARRAY_TASK_ID:-0}`.
- `worker_bump_retry_count <max>` → read count (0 if missing), return 1 if count >= max, else write count+1 atomically (tmp+`mv`) and echo the new count.
- `worker_clear_retry_count` → `rm -f` the counter file.
- `worker_requeue_if_transient <err_file> [max]`:
  - greps `<err_file>` (the Slurm `--error` file for this task; must exist and be non-empty) with `TRANSIENT_REQEX`.
  - If match and `worker_bump_retry_count` succeeds → echo "Transient failure detected; requeueing (attempt N)" → `scontrol requeue "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"` (or `${SLURM_JOB_ID}` for non-array) → `sleep 2` → return 0 (caller exits 0; the task restarts).
  - Else return 1 (caller propagates the real exit code).
- `stage_env_rlib [prefix]`:
  - Source: `${PROJECT_ROOT}/.pixi/envs/py-cuda13/lib/R/library`. Guard: if `WORKER_STAGE_R_LIB` != `1` (default 1) or size > `WORKER_R_LIB_MAX_MB` (default 10240 MB, measured via `du -sm`) → warn and skip (fall back to BeeGFS; retry still covers flakes).
  - Target: `${STAGE_ROOT}/${prefix}_R_library`, `STAGE_ROOT="${WORKER_STAGE_ROOT:-/scratch}/${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}"` (fall back to `/tmp` if `/scratch` absent).
  - `cp -a` (log start/end + seconds + size); on failure return non-zero (the caller's transient check covers stale-view `No such file` errors).
  - On success: `export R_LIBS="${STAGE_ROOT}/${prefix}_R_library${R_LIBS:+:${R_LIBS}}"` (default env library remains as fallback in `.libPaths()`).
- `export_worker_thread_env` → `export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`. Annotation worker only.

### B. Annotation worker: `src/4_cell_type_annotation/2.1_run_worker.sh`

1. `#SBATCH --time=02:00:00` (was 01:00:00). Do NOT add `--requeue` yet — validation step V2 determines whether manual requeue of a RUNNING task works on this cluster without it; if not, add it.
2. After the feather-skip check (lines 61-68, unchanged — skip stays fast and before staging):
   - source `"${SCRIPT_DIR}/../utils/bash/worker_retry.sh"`
   - `export_worker_thread_env`
   - run R with staging and unified retry handling:
     ```bash
     set +e
     stage_env_rlib "annotation" && ${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.1.1_process_chunk.R" "${CHUNK_FILE}"
     RC=$?
     set -e
     if [[ ${RC} -eq 0 ]]; then
       worker_clear_retry_count
       echo "Task ${SLURM_ARRAY_TASK_ID}: chunk processing complete."
       exit 0
     fi
     ERR_FILE="${LOGS_DIR}/4_cell_type_annotation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
     if worker_requeue_if_transient "${ERR_FILE}" "${WORKER_MAX_RETRIES:-3}"; then
       exit 0   # requeued; script will restart on (likely) another node
     fi
     exit ${RC}
     ```
   - Both staging and R stderr land in the Slurm `.err` (worker stderr + R stderr), so one grep covers both.
   - Echo the staged path + `R_LIBS` for log visibility.

### C. Annotation R script: `src/4_cell_type_annotation/2.1.1_process_chunk.R`

1. **Wall-clock budget plumbing** (near the top, after args parsing):
   ```r
   wall_limit_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_TIME_LIMIT", "7200")))
   if (is.na(wall_limit_s) || wall_limit_s <= 0) wall_limit_s <- 7200
   job_start_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_START_TIME", "")))
   wall_elapsed <- function() if (!is.na(job_start_s)) as.numeric(Sys.time()) - job_start_s else proc.time()[3]
   wall_left <- function() wall_limit_s - wall_elapsed()
   ```
2. **Budget-aware attempts** (scATOMIC block lines 185-211 and HiTME block lines 214-245):
   - `attempt_timeout <- min(timeout, 1800)` (timeout = existing ncol formula, line 182).
   - Inside each attempt loop, before `withTimeout`: `eff <- min(attempt_timeout, wall_left() - 300); if (eff < 60) { message("...insufficient wall time remaining; skipping attempt..."); break }`, then `withTimeout(..., timeout = eff)`.
   - Keep the existing `TimeoutException`/`error` tryCatch handlers unchanged (best-effort fast-fail; verified working at R checkpoints).
3. **Per-sample checkpoints** (replace the tail of the sample loop, lines 247-262):
   - `tmp_dir <- file.path(paths$path_output, "annotation_tmp", sub("\\.txt$", "", basename(args$chunk_file)))` (e.g. `output/annotation_tmp/chunk_48`). Must NOT match the merge/coverage glob `annotations_chunk_*.feather` (it is a subdirectory — safe).
   - Per sample `i` (position in `samples_to_process`): intermediate `file.path(tmp_dir, sprintf("sample_%02d.feather", i))`.
     - If the intermediate already exists → `message("resume: sample ... already annotated")` → `next`.
     - Else wrap the existing per-sample body (scATOMIC + HiTME + `meta`/`annot` extraction) in a sample-level `tryCatch(error = function(e) { message("SAMPLE FAILED (", target_sample, "): ", conditionMessage(e), " - continuing"); FALSE })`; on success write the intermediate atomically (`write_feather` to `*.tmp` + `file.rename`) and return TRUE.
     - Record failed samples in a vector; keep `rm(seurat_obj); gc()` per sample.
   - After the loop:
     - If no failures: read all intermediates in index order, `do.call(rbind, ...)` (same columns as today: annot cols + `cell_barcode` + sample colname), write the final feather atomically (`annot_file` unchanged), `unlink(tmp_dir, recursive = TRUE)`, exit normally.
     - If any failures: `message("chunk INCOMPLETE - failed samples: ...")`, `quit(status = 1)` WITHOUT writing the final feather. Intermediates persist → a later re-run resumes only the failed sample(s) (worker's feather-skip check requires the final feather, so the chunk is correctly re-submitted).
   - This makes any timeout escape / late-fire cost at most the in-flight sample, and turns `scATOMIC`-deterministic failures into cheap resume cycles instead of full-chunk re-runs.

### D. Chunk rebuild cleanup: `src/4_cell_type_annotation/1.1_prepare_chunks.py`

In the stale-feather cleanup block (lines 424-431, production mode only), also delete `annotation_tmp/` under `path_data` (mirror logic/condition of the feather deletion): chunk numbering changes on rebuild, so old intermediates must not be resumed.

### E. NAS rsync excludes (annotation_tmp must never reach NAS)

- `src/4_cell_type_annotation/2_submit_hpc_array.sh` line 341: add `--exclude='annotation_tmp/'`.
- `src/4_cell_type_annotation/3_submit_merge.sh` line ~232: add `--exclude='annotation_tmp/'`.

### F. Preprocess worker: `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`

- Source `worker_retry.sh`; wrap the `${PYTHON_BIN} 1.1.1_preprocess.py ...` invocation (lines 49-54) with the same `set +e`/RC/requeue-if-transient pattern, err file `${LOGS_DIR}/3_scrnaseq_preprocessing_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err`.
- No staging, no thread pinning (python env too large; preprocess is multi-threaded by design).
- **Pre-check before finalizing**: verify `1.1.1_preprocess.py` skip/overwrite semantics — if a re-run after a mid-write crash could exit 0 with partial outputs, requeue is NOT safe (would mask incompleteness). If unsafe, restrict retry to failures that occur before any output writing is possible (e.g., only when the error is a load/import failure), or skip retry for preprocess and note it. Validation step V3 covers this.

### G. Benchmark workers (retry; R workers also stage)

- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`: source helper; wrap the main python invocation; err file `${LOGS_DIR}/5_benchmark_${METHOD}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err`.
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh` and `run_transformation_zeroimp_analysis/1.1_run_worker.sh`: same + `stage_env_rlib "benchmark"` before the R invocation (env `R_LIBS` export), err files `${LOGS_DIR}/5_benchmark_r_${METHOD}_...err` and `${LOGS_DIR}/5_transzeroimp_${ANALYSIS}_...err`.
- Do NOT pin threads, do NOT change `--time` (pinned resource classes).
- Pre-check: benchmark method outputs are overwritten on re-run (idempotent) — verify per method before enabling retry; if a method's partial output could be mistaken for complete, keep retry but ensure the worker re-runs produce fresh outputs (most write per-dataset output files; confirm).

### H. Docs

- `AGENTS.md`: extend the BeeGFS section — workers now self-requeue on transient signatures (`src/utils/bash/worker_retry.sh`, counter dir `${HPC_SCRATCH_DIR}/_worker_retries/`, `WORKER_MAX_RETRIES` default 3); R workers stage `lib/R/library` to node-local `/scratch` (`R_LIBS`, `WORKER_STAGE_R_LIB`, `WORKER_R_LIB_MAX_MB`); annotation tasks: 2h wall, per-attempt cap 1800 s, budget-aware attempts, per-sample checkpoints (`output/annotation_tmp/`); thread pinning is annotation-only (benchmark comparability).
- `docs/ARCHITECTURE.md`: annotation worker flow (staging → R run → retry), per-sample checkpoint contract (intermediates + atomic final feather; merge coverage gate unchanged), `annotation_tmp/` + `_worker_retries/` in the HPC folder layout, rsync excludes.
- No README/TODO changes.

## Validation

- **V1 (done)**: withTimeout/reticulate escape test (documented above) — evidence for the design, not a code change.
- **V2 (HPC, before the real run)**: `scontrol show job <any-running-id> -o | grep Requeue` and a tiny 2-task debug array that fails once with an echo of `cannot open file ... No such file or directory` + `exit 1`, then succeeds on the second attempt via the counter → verify REQUEUED → COMPLETED, counter cleaned, and that manual `scontrol requeue <job>_<task>` of a RUNNING task works without `--requeue` (if not, add `#SBATCH --requeue` to the annotation worker).
- **V3 (HPC)**: `du -sm ~/ECODA_paper/.pixi/envs/py-cuda13/lib/R/library` (size; drives the staging guard) and one small sbatch test measuring `cp -a` time to `/scratch` on a worker node + `df -h /scratch` (node disk headroom). Verify preprocess/benchmark re-run idempotency (tasks F/G pre-checks).
- **V4 (HPC, end-to-end)**: `./src/4_cell_type_annotation/2_submit_hpc_array.sh _debug` — check: staging echo + `R_LIBS` in logs, intermediates appear under `output/annotation_tmp/chunk_*/` during the run, final feather written, `annotation_tmp/` removed on success, merge via `3_submit_merge.sh _debug` passes the coverage gate, no `annotation_tmp/` on NAS.
- **V5 (HPC, resume)**: `scancel` one `_debug` task mid-run; resubmit; verify only the missing sample(s) re-ran (intermediates skipped).
- **V6**: real run, one dataset first (`./src/4_cell_type_annotation/2_submit_hpc_array.sh Gongsharma_cmv_young_males`), then full. Watch a full cycle for TIMEOUT-free completion and zero stale-view failures.

## Risks / notes

- Requeued tasks may land on the same node (no node exclusion) — staging makes env-staleness moot; scratch-data staleness gets another chance via the retry cap; `SBATCH_EXCLUDE` remains available to the user.
- Staging copy load on BeeGFS metadata servers (up to 1000 concurrent copies of the same tree): the size guard + measured copy time (V3) decide whether the node-shared `/srv/share/users/...` optimization becomes a follow-up.
- Stale retry counter files (task killed after success) are harmless: filenames embed job+task ids, so they never collide across submissions.
- Pre-existing env inconsistency (ABI warnings: SeuratObject/Matrix built vs loaded versions): the documented definitive repair (`rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup`) is OUT OF SCOPE here — it re-opens the staleness window and must be done as a separate, planned step with the ~1h wait afterwards.
- The currently-running array `4308096` is unaffected; task 378 will likely TIMEOUT at 1h — after this implementation, its dataset (`Pelka`) re-run is cheap thanks to per-sample resume.

# IMPORTANT NOTE (by the user):
Do not commit and push anything after implementing this plan.
Notify the user to double-check.
Only commit and push all the changes (including the plan) once the user checked and confirmed the changes.