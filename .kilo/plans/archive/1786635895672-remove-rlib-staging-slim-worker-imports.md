# Plan: Remove R-library staging + disentangle worker imports

## Context & diagnosis

Symptom (session ses_005c, 2026-08-13): transzeroimp arrays (4313320/4313321) hit a
recurring BeeGFS metadata storm — ~22 tasks simultaneously run `du -sm` + `cp -a` on
`.pixi/envs/py-cuda13/lib/R/library` (~1.3 GB, ~500k files), most `cp` reads hit
stale-cache ENOENT, tasks self-requeue (retry budget 6), and the storms repeat.

Two amplifiers identified and confirmed in code:

1. **Per-task full-library staging** (`stage_env_rlib` in `src/utils/bash/worker_retry.sh:121-179`):
   copies the entire R library to node-local `/scratch` per task. The storm source.
2. **Monolithic imports**: every R benchmark/transzeroimp worker sources
   `src/utils/load_all_functions.R` → `src/utils/imports.R`, attaching all 42 packages
   (Seurat, MOFA2, scITD, GloScope, HiTME, tidyverse, EPIC, …) even for a 30 s
   trans/zeroimp task. Runtime BeeGFS footprint far larger than needed.

`WORKER_STAGE_R_LIB=0` (env-var kill switch) already skips both the `du` and the `cp`
— this plan makes that path the only path, and slims the imports so direct-env
execution is cheap. The transient-signature self-requeue mechanism
(`worker_requeue_if_transient`, `WORKER_MAX_RETRIES=6`) stays as the safety net.

Decisions (user, 2026-08-13):
- Staging removal is **complete**: annotation worker included (its package set is
  already slim and per-sample checkpoints make retries cheap).
- **No array throttling** in this pass (deferred; see "Deferred follow-ups").
- **One implementation pass** (staging removal + import disentanglement), then
  `_debug` validation on HPC, then resubmit the pending production runs.

## Constraints

- `src/utils/imports.R` stays **untouched**: it is the env-verification list for the
  guarded env-refresh smoke checks (`src/utils/bash/setup_env_sbatch.sh:150-159`,
  `src/utils/bash/refresh_env.sh:102-111`) and the notebook loader. New worker import
  files are *subsets* of it.
- `src/utils/load_all_functions.R` stays **untouched** (notebooks:
  `benchmark_analysis.rmd`, `batch_effect_analysis.rmd` keep full behavior).
- Utils are already mostly namespace-qualified (`DESeq2::`, `limma::`,
  `zCompositions::`, `vegan::`, `mclust::`, `cluster::`, `igraph::`, `Matrix::`,
  `GloScope::`, `arrow::`, `parallelly::`, `SummarizedExperiment::`,
  `MatrixGenerics::`, `jsonlite::`, `tibble::`, `tools::`, `proxy::`). Bare-attach
  packages actually required by the executed worker code paths:
  - benchmark/prepare_pseudobulk workers: **Seurat** (bare `CreateSeuratObject`,
    `FindVariableFeatures`, `ScaleData`, `RunPCA`, `CreateDimReducObject` in
    `src/utils/seurat_utils.R`), **dplyr** (bare verbs + `%>%`), **reticulate**
    (bare `import()`/`py_to_r()`).
  - per method: **MOFA2** (bare `create_mofa` in
    `benchmark_methods_r.R:230`, mofa only), **scITD** (bare `initialize_params`,
    `make_new_container` in `process_scitd_fig`, scitd only). GloScope already
    namespaced (`GloScope::gloscope`) — attach not needed.
  - trans/zeroimp workers: **dplyr**, **reticulate**, **foreach**, **doParallel**.
    **No Seurat** (obs-only backed reads; `get_ct_comp_df` is base + pipe only).
    Note: `datrans` is a **local function** (`benchmark_pipeline.R:5`), NOT an
    scECODA call — scECODA is never referenced in worker paths (only
    `imports.R:35` + a commented-out notebook line), so it is NOT attached.
    However, local `datrans` bare-calls `makeCluster` (`parallel`),
    `registerDoParallel` (`doParallel`) and `foreach` — attaching `foreach` +
    `doParallel` pulls `parallel`/`iterators` via Depends, exactly as the current
    full `imports.R` does.
- No env refresh needed: no packages are added/removed (only load lists change).

## Task 1 — Remove the staging machinery (complete)

1. `src/utils/bash/worker_retry.sh`:
   - Delete `stage_env_rlib()` (lines 121–179).
   - Delete its header-comment block (lines 37–51) and the "R library staging"
     mention in the file header (lines 3–4).
   - Keep `TRANSIENT_REQEX`, retry counters, `worker_requeue_if_transient`,
     `export_worker_thread_env` unchanged.
2. Call sites (drop the `stage_env_rlib "..." &&` prefix, update the surrounding
   "Staging + unified retry handling" comments to "unified retry handling"):
   - `src/4_cell_type_annotation/2.1_run_worker.sh:79` →
     plain `${PIXI_RSCRIPT} "${SCRIPT_DIR}/2.1.1_process_chunk.R" "${CHUNK_FILE}"`.
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh:78`.
   - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh:82`.
   - (Preprocess + python-benchmark workers already never call it — no change.)
3. `src/slurm_config.sh`: remove `WORKER_STAGE_R_LIB` and `WORKER_R_LIB_MAX_MB`
   (lines 149–159 block); keep `WORKER_MAX_RETRIES=6`, reword the block comment to
   "worker self-healing (retry only)".

## Task 2 — Disentangle worker imports

New files in `src/utils/` (mirror the `imports.R` load-and-stop pattern;
`suppressPackageStartupMessages` + missing-package `stop()`):

1. `imports_worker_core.R` — attach `Seurat`, `dplyr`, `reticulate` (this order
   matches the existing relative order in `imports.R`). Header comment: subset of
   `imports.R` for R benchmark workers; `imports.R` remains the canonical
   env-verification list; namespace-qualified packages load on demand and must not
   be attached here.
2. `imports_worker_transzeroimp.R` — attach `doParallel`, `foreach`, `reticulate`,
   `dplyr` (NOT scECODA: `datrans` is local at `benchmark_pipeline.R:5`; this
   order mirrors `imports.R`'s relative order, lines 6/8/27/46).
   Header comment: obs-only backed reads; **no Seurat**.
3. `load_worker_functions.R` — sources the util files in `load_all_functions.R`
   dependency order **minus** `imports.R` and `plotting.R`:
   `datasets_io.R`, `constants.R`, `helpers.R`, `math_utils.R`, `scoring_metrics.R`,
   `pseudobulk.R`, `hvcs.R`, `seurat_utils.R`,
   `src/5_run_benchmark_methods/benchmark_methods_r.R`,
   `src/5_run_benchmark_methods/benchmark_pipeline.R`. Header comment: worker-only
   loader; notebooks use `load_all_functions.R`. (Implementer: verify no top-level
   executable statements beyond function definitions exist in these files — only
   `source()` order should matter.)
4. Update the 4 entry R scripts (replace the `load_all_functions.R` source + drop
   the now-redundant `library(reticulate)`; keep sourcing `benchmark_hpc_utils.R`):
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:31-34`:
     source `imports_worker_core.R` + `load_worker_functions.R` +
     `benchmark_hpc_utils.R`; after the method is parsed (line ~51), add conditional
     attaches: `if (method == "mofa") library(MOFA2)`, `if (method == "scitd")
     library(scITD)` (gloscope needs only the installed namespace).
   - `.../1.1.1_prepare_pseudobulk.R:20-23`: source `imports_worker_core.R` +
     `load_worker_functions.R` + `benchmark_hpc_utils.R`.
   - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1.1_run_transformation_analysis.R:21-24`
     and `.../1.1.1_run_zeroimp_analysis.R:21-24`: source
     `imports_worker_transzeroimp.R` + `load_worker_functions.R` +
     `benchmark_hpc_utils.R`.
   - `src/4_cell_type_annotation/2.1.1_process_chunk.R`: **no change** (already the
     slim-import model).
5. Extend the env smoke checks so env drift in the subsets fails loudly at refresh
   time: in `src/utils/bash/setup_env_sbatch.sh` (after line 158) and
   `src/utils/bash/refresh_env.sh` (after line 110), also
   `source("src/utils/imports_worker_core.R")` and
   `source("src/utils/imports_worker_transzeroimp.R")` in the same Rscript smoke
   check (MOFA2/scITD remain covered by the full `imports.R` list).

## Task 3 — Documentation

- `AGENTS.md` (~line 166, "Worker self-healing" paragraph): remove the
  `stage_env_rlib` / `WORKER_STAGE_R_LIB` / `WORKER_R_LIB_MAX_MB` description;
  state that R workers run directly against the env library with slim per-class
  imports (`imports_worker_core.R` / `imports_worker_transzeroimp.R` /
  `load_worker_functions.R`; `imports.R` = env-verification list only), and that
  the retry mechanism covers residual stale-view flakes.
- `docs/ARCHITECTURE.md`: update the staging mentions in the `2.1_run_worker.sh`
  row (~line 153), the worker self-healing section (~line 163), the R benchmark
  worker row (~line 389) and the transzeroimp worker row (~line 506); add the new
  import files to the R-module file table.
- `TODO.md`: add deferred follow-up — CPU benchmark array throttling
  (`BENCHMARK_CPU_ARRAY_THROTTLE`, e.g. 4, for the R/transzeroimp/PILOT submitters;
  `MAX_NUM_CHUNKS_PARALLEL` untouched) as a separate design if direct-env imports
  still show startup storms; node-shared `/srv/share/users/...` staging as the
  documented long-term alternative.

## Task 4 — Static validation (local, no pipeline runs)

- `bash -n` on the three edited worker scripts.
- Parse-only R checks (no package loads):
  `Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f))' <files>` for all
  edited/created R files.
- Confirm `imports.R`, `load_all_functions.R`, and `2.1.1_process_chunk.R` are
  byte-identical to HEAD (git diff).

## Task 5 — HPC validation + resubmit (user-run)

Commit + push (AGENTS.md task-completion workflow: archive the plan, `git add .`,
commit, push), then on the HPC login node:

```bash
cd ~/ECODA_paper && source src/slurm_config.sh && cd "${PROJECT_ROOT}"
git pull                     # no env refresh needed (no new packages)

# 1. (optional, cheap) Load-time before/after numbers for the docs — run on the
#    debug-cpu partition (never the login node); expect the worker subsets to be
#    several× faster than the 42-package imports.R:
sbatch --partition=debug-cpu --mem=16G -t 00:20:00 -o "${HOME}/import_timing_%j.out" <<'EOF'
#!/bin/bash
source "${SLURM_SUBMIT_DIR}/src/slurm_config.sh" && cd "${PROJECT_ROOT}"
for f in src/utils/imports.R src/utils/imports_worker_core.R src/utils/imports_worker_transzeroimp.R; do
  /usr/bin/time -f "$f: %e s" "${HOME}/.pixi/bin/pixi" run --as-is -e py-cuda13 Rscript --vanilla \
    -e "t <- system.time(source('$f')); cat('$f:', t[['elapsed']], 's\n')" 2>&1 | tail -1
done
EOF
#   (Runs twice per file — bash `time` and R `system.time` — take the R number.)

# 2. _debug smoke validation (starts fast; validates slim imports + retry path):
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --ds_name _debug --analysis trans,zeroimp --force
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name _debug --methods prepare_pseudobulk,pseudobulk,gloscope,mofa,scitd
#   (mofa/scitd combos are skipped on 5 samples, but MOFA2/scITD attach at startup
#    → the conditional library() loads are exercised. Check the .err logs for the
#    tasks that ran.)

# 3. Only after _debug passes, the pending production work:
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --analysis trans,zeroimp --force          # 130bde6 breaking-change re-run (required)
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh
#   (re-runs prepare_pseudobulk fast (cached), then gloscope/mofa/pseudobulk/scitd;
#    tail syncs everything incl. execution_times.feather + checksums.md5 to NAS)
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Lee   --methods scpoli --partition "${SLURM_PARTITION_PRIVATE}"
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Zhang --methods scpoli --partition "${SLURM_PARTITION_PRIVATE}"
#   (gpu004: evening 18:00+ window; daytime blocked by the TRES=cpu-64 reservation —
#    IT email about adding gres/gpu=1 remains a separate, user-owned follow-up)
```

Verification: Step-1/Step-2 status tables from session ses_005c (one line per
dataset); `_worker_retries/` counters empty for successful tasks; NAS
`execution_times.feather` + `checksums.md5` present after the R submitter tail;
scpoli feathers for Lee/Zhang on NAS.

## Risks & rollback

- Missing package in a slim set → loud `stop()` at startup (mirrors `imports.R`
  pattern), caught by the `_debug` run; add the package to the relevant worker
  import file (no env refresh needed). The highest-risk runtime gap is
  trans/zeroimp's local `datrans` (`foreach`/`doParallel`/`parallel` bare calls) —
  the `_debug` **trans** run specifically exercises it (zeroimp never calls
  `datrans`), so do not skip the trans half of the smoke run.
- Load-order regression → mitigated by mirroring the existing relative order
  (Seurat → … → MOFA2 → scITD; doParallel/foreach → reticulate → dplyr).
- Residual stale-view flakes → `worker_requeue_if_transient` (budget 6) unchanged;
  a failed task still fails closed and its `.err` tail shows the real error.
- Rollback: revert the commit; `WORKER_STAGE_R_LIB` path is gone by design — the
  pre-refactor code is restorable from git, no env changes involved.

## Deferred follow-ups (explicitly out of scope)

- CPU array throttling (`BENCHMARK_CPU_ARRAY_THROTTLE=4`) — only if storms persist
  (see TODO.md entry).
- Node-shared `/srv/share/users/...` staging — long-term alternative, documented,
  not implemented.
- gpu004 reservation TRES (GPU) — IT ticket, user-owned.
