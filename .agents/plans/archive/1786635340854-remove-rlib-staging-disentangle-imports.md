# Plan: Remove R-library staging + disentangle worker imports

## Goal

Eliminate the BeeGFS metadata-storm failure mode observed 2026-08-13 (11–22 array tasks concurrently `cp -a`-ing the ~1.3 GB / ~500k-file pixi R library caused rotating `cp: cannot open ... No such file or directory` ENOENTs on worker nodes, transient cluster incident). Two changes:

1. **Remove the staging machinery entirely** (`stage_env_rlib`, `WORKER_STAGE_R_LIB`, `WORKER_R_LIB_MAX_MB`) — workers read packages directly from BeeGFS (exactly the already-proven `WORKER_STAGE_R_LIB=0` behavior; Python workers already do this with zero staging and succeeded for all datasets).
2. **Disentangle `src/utils/imports.R`** — today every benchmark/trans/zeroimp R task attaches all 42 packages via `load_all_functions.R`, pulling a huge fraction of the library tree on every task start. Replace with slim per-class import files so a 30-second trans task no longer imports Seurat/MOFA2/scITD/funkyheatmap/EPIC.

The transient-retry mechanism (`worker_retry.sh` TRANSIENT_REQEX grep + counter-capped self-requeue, `WORKER_MAX_RETRIES=6`) stays — it covers residual flakes of direct-BeeGFS reads.

## Context (verified in code — corrects the second opinion)

- `stage_env_rlib()` lives at `src/utils/bash/worker_retry.sh:121-179` (header doc lines 37-51), called from exactly **3** workers:
  - `src/4_cell_type_annotation/2.1_run_worker.sh:79` (`stage_env_rlib "annotation" && ${PIXI_RSCRIPT} ...`)
  - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh:78` (`stage_env_rlib "benchmark" && ...`)
  - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh:82` (`stage_env_rlib "benchmark" && ...`)
- Config vars: `src/slurm_config.sh:149-159` (`WORKER_STAGE_R_LIB=1`, `WORKER_R_LIB_MAX_MB=10240`, `WORKER_MAX_RETRIES=6`). **Not** referenced in `benchmark_hpc_utils.R` (second opinion was wrong there) or anywhere else.
- `R_LIBS` is only ever set by `stage_env_rlib`; after removal, R resolves the env library via `.libPaths()` by default (env `lib/R/library` is the env's own library) — nothing else needed.
- `imports.R` (`src/utils/imports.R:5-48`, 42 packages via `require()`) is sourced only through `src/utils/load_all_functions.R:7`, which the **4** worker R scripts source:
  - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R:31`
  - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R:20`
  - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1.1_run_transformation_analysis.R:21`
  - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1.1_run_zeroimp_analysis.R:21`
  - All 4 additionally source `benchmark_hpc_utils.R` explicitly (keep that line).
- **`imports.R` must stay untouched**: it is the env-verification canonical list used by `src/utils/bash/setup_env_sbatch.sh:150-158` and `src/utils/bash/refresh_env.sh:102-110` smoke checks.
- **Annotation worker is already the lightweight model**: `2.1.1_process_chunk.R:11-25` sources only `seurat_utils.R` + explicit `library()` calls; it does not use `load_all_functions.R`. Its `stage_env_rlib` call is removed in Phase 1 but its imports stay as-is (scATOMIC/HiTME/scGate genuinely need those packages).
- Package usage in utils/benchmark files is **mostly namespaced** (`DESeq2::` ×9, `limma::` ×4, `zCompositions::` ×5, `Hotelling::` ×3, `vegan::`, `Matrix::`, `arrow::`, `plotly::`, ...). Attached-only (bare-call) heavy packages in the benchmark pipeline: **MOFA2, scITD, edgeR, scECODA, funkyheatmap, EPIC, anndataR** (zero `pkg::` occurrences). `%>%` is used unnamespaced in utils (`seurat_utils.R:153`, etc.) → **dplyr must be attached** in every worker class.
- `run_transformation_analysis`/`run_zeroimp_analysis` live in `benchmark_pipeline.R:1184/1212` (not in a separate file); `run_ct_comps_analysis_worker` in `benchmark_hpc_utils.R:233`; trans/zeroimp use `datrans` (scECODA), `get_ct_comp_df` (seurat_utils.R:166, data.frame variant — no Seurat object needed), `calc_perc_df`, `impute_zeros`, `clr`, `calc_sep_score` (math_utils/scoring_metrics).
- `plotting.R` is notebook-only (no plotting/funkyheatmap/ggplot references in benchmark_methods_r.R or benchmark_pipeline.R).
- CPU arrays are unthrottled (`--array=1-N%${MAX_NUM_CHUNKS_PARALLEL}`, =1000) in all 4 submitters; GPU arrays use `BENCHMARK_GPU_ARRAY_THROTTLE` (=4).

## Decisions (user-confirmed)

- **Array throttling: DEFERRED** — follow-up decision after direct-BeeGFS imports are proven under concurrency (~2 weeks of runs). Not in this plan.
- **Import granularity: per-class files** (not per-method, not inline-in-script):
  - Mirrors the existing `imports.R` + `load_all_functions.R` repo convention.
  - Avoids duplicating package blocks across the 4 worker scripts.
  - Per-method splits are unnecessary: the 4 heavy methods run 12h tasks where a few extra attached packages are negligible; the tiny tasks (trans/zeroimp) are the ones that need slimming.
  - `imports.R` stays as the env-verification list; `load_all_functions.R` stays as the notebook/interactive everything-loader.

## Phase 1 — Remove the staging machinery (small, independent, safe)

1. **`src/utils/bash/worker_retry.sh`**: delete `stage_env_rlib()` (lines 121-179) and its header doc block (lines 37-51, the `#   stage_env_rlib [prefix]` bullet). Keep `TRANSIENT_REQEX`, retry counter helpers, `worker_requeue_if_transient`, `export_worker_thread_env`. Update the file header (line 3-4: drop "+ R library staging" wording).
2. **Call sites** — drop the `stage_env_rlib "benchmark" &&` / `stage_env_rlib "annotation" &&` prefix (keep `set +e`, `RC=$?`, `set -e`, retry block as-is):
   - `src/4_cell_type_annotation/2.1_run_worker.sh:79`
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh:78`
   - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh:82`
   - Also update the "Staging + unified retry handling" comments above each call (now "Direct-env R library reads + unified retry handling"; the `.err` single-grep rationale still holds — R stderr still lands there).
3. **`src/slurm_config.sh:149-159`**: delete the `WORKER_STAGE_R_LIB` / `WORKER_R_LIB_MAX_MB` lines + their comments; keep `WORKER_MAX_RETRIES=6`.
4. **Docs**: replace staging language in `AGENTS.md` (~line 166) and `docs/ARCHITECTURE.md` (rows ~153, 163, 389, 506) with "R workers read the pixi env library directly from BeeGFS; transient flakes covered by self-requeue". Keep the `/srv/share/users` node-shared-cache line as a documented follow-up (unchanged).

## Phase 2 — Per-class slim import files

### 2a. Dependency map (implementation-time verification task)

For each worker class, read the sourced files and list every package: attach (bare calls / `%>%`) vs `requireNamespace` (namespaced-only). Expected outcome (verify, don't trust blindly):

| Class | Attach | Namespaced-only checks (requireNamespace) |
|---|---|---|
| **benchmark R** (prepare_pseudobulk + gloscope/mofa/pseudobulk/scitd) | dplyr, Seurat, anndataR, reticulate, MOFA2, scITD, edgeR, scECODA | DESeq2, limma, Matrix, arrow, zCompositions, Hotelling, GloScope (1 namespaced use), vegan, igraph, cluster, proxy, factoextra, mclust, parallelly, BiocParallel, SummarizedExperiment, MatrixGenerics, tibble, jsonlite, stats |
| **transzeroimp** | dplyr, reticulate, scECODA | zCompositions, DESeq2 (via pseudobulk.R sourcing? verify), math/scoring namespaced deps, jsonlite, arrow |

Notes: `benchmark_pipeline.R` is sourced but only *defines* functions (no top-level package calls) — safe to source without attaching its packages. `plotting.R` is NOT sourced by workers (notebook-only). `pseudobulk.R`/`hvcs.R`/`scoring_metrics.R` use namespaced DESeq2/vegan etc. — included via the namespaced check list. Annotation worker (`2.1.1_process_chunk.R`) is left as-is.

### 2b. New files

1. **`src/utils/imports_core.R`** — minimal shared loader:
   - Attach: `dplyr` (provides `%>%`), `jsonlite` (if utils need it bare — verify), `arrow` (benchmark_hpc_utils uses `arrow::read_feather` namespaced — check-only).
   - `requireNamespace`-check the namespaced-only list with the same fail-loudly `stop()` style as `imports.R:60-66` ("missing from the pixi environment" message).
   - Reuse the `my_packages` + `load_my_packages` pattern (copy, renamed e.g. `core_packages`/`check_packages`).
2. **`src/utils/imports_benchmark.R`** — `source("src/utils/imports_core.R")` (repo-relative, matching `load_all_functions.R` style) then attach benchmark class packages (Seurat, anndataR, reticulate, MOFA2, scITD, edgeR, scECODA) + their checks.
3. **`src/utils/imports_transzeroimp.R`** — `source("src/utils/imports_core.R")` + attach (reticulate, scECODA) + checks.
4. **`src/utils/load_benchmark_functions.R`** — class loader mirroring `load_all_functions.R` but with slim imports and no `plotting.R`:
   `imports_benchmark.R` → `datasets_io.R` → `constants.R` → `helpers.R` → `math_utils.R` → `scoring_metrics.R` → `pseudobulk.R` → `hvcs.R` → `seurat_utils.R` → `benchmark_methods_r.R` → `benchmark_pipeline.R` (verify this order against `load_all_functions.R:7-18` minus `imports.R`/`plotting.R`).
5. **`src/utils/load_transzeroimp_functions.R`** — `imports_transzeroimp.R` → the utils subset trans/zeroimp actually need (expected: `datasets_io.R`, `constants.R`, `helpers.R`, `math_utils.R`, `scoring_metrics.R`, `seurat_utils.R`, `benchmark_pipeline.R` for `run_transformation_analysis`/`run_zeroimp_analysis`; verify `pseudobulk.R`/`hvcs.R` need).

### 2c. Switch the 4 worker R scripts

Replace `source(file.path(project_root, "src/utils/load_all_functions.R"))` with:
- `1.1.1_run_benchmark_methods_r.R:31` and `1.1.1_prepare_pseudobulk.R:20` → `source(file.path(project_root, "src/utils/load_benchmark_functions.R"))`
- `1.1.1_run_transformation_analysis.R:21` and `1.1.1_run_zeroimp_analysis.R:21` → `source(file.path(project_root, "src/utils/load_transzeroimp_functions.R"))`

Keep the existing `source(... benchmark_hpc_utils.R)` and `library(reticulate)` lines. Do NOT touch `imports.R` or `load_all_functions.R` (notebook + env-verification roles).

## Phase 3 — Validation

1. **Local sanity** (no package loads on macOS — py-cuda13 is linux-only): `Rscript -e 'parse("src/utils/imports_core.R")'` etc. for syntax; `bash -n` the edited shell scripts.
2. **HPC smoke test** (debug-cpu partition, per AGENTS.md convention — never the login node):
   ```
   sbatch --partition=debug-cpu --time=00:10:00 --mem=8G -o ${HOME}/import_smoke_%j.out \
     -e ${HOME}/import_smoke_%j.err <<'EOF'
   #!/bin/bash
   cd "${HOME}/ECODA_paper" && source src/slurm_config.sh && cd "${PROJECT_ROOT}"
   ${PIXI_RSCRIPT} -e 't <- system.time(source("src/utils/load_benchmark_functions.R")); cat("benchmark OK", t[["elapsed"]], "s\n")'
   ${PIXI_RSCRIPT} -e 't <- system.time(source("src/utils/load_transzeroimp_functions.R")); cat("transzeroimp OK", t[["elapsed"]], "s\n")'
   EOF
   ```
   Expect: both loaders succeed; record elapsed (transzeroimp should be a small fraction of the previous 42-package load).
3. **`_debug` pipeline smoke runs** (the repo's established validation target, per AGENTS.md):
   ```
   ./src/4_cell_type_annotation/2_submit_hpc_array.sh --ds_name _debug
   ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug
   ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name _debug --force
   ```
   Verify outputs: annotation feathers complete; `_debug` R bundles (3: gloscope/mofa/pseudobulk, no scITD for 5 samples) with `exec_time` set; trans/zeroimp bundles present. Compare structure against pre-change bundles.
4. **Regression check**: `_debug` + Kfoury R-method bundles byte-identical structure to previous runs (or recomputed identically with `--force`).

## Phase 4 — Docs

- Update `AGENTS.md` (~line 166 worker-self-healing paragraph) and `docs/ARCHITECTURE.md` (rows 153, 163, 389, 506 + any `worker_retry.sh` references): staging removed; workers read the env library directly from BeeGFS; slim per-class import files (`imports_core/benchmark/transzeroimp`, `load_benchmark_functions.R`, `load_transzeroimp_functions.R`); `imports.R` remains the env-verification list.
- Add a TODO.md note: re-evaluate array throttling (`BENCHMARK_CPU_ARRAY_THROTTLE`) after ~2 weeks of direct-BeeGFS imports if flakiness reappears.

## Risks & rollback

- **Phase 1**: pure removal, behaviorally identical to the proven `WORKER_STAGE_R_LIB=0` path (used for the successful trans/zeroimp re-runs). Rollback = revert commit.
- **Phase 2**: a missing package in a slim class file → loud `stop()` at task start (same style as today's `imports.R`), caught by the Phase 3 smoke + `_debug` runs. Wrong load order → mirror `load_all_functions.R:7-18` order per class. `%>%` availability guaranteed by attaching dplyr in `imports_core.R`.
- Residual direct-BeeGFS flakes → existing `worker_requeue_if_transient` + `WORKER_MAX_RETRIES=6` self-healing unchanged.
- Do NOT delete `imports.R` or `load_all_functions.R` (setup_env_sbatch.sh/refresh_env.sh smoke checks + notebooks depend on them).

## Out of scope / follow-ups

- Array throttling (deferred by user decision; TODO.md note in Phase 4).
- `tidyverse` → explicit dplyr/tidyr/forcats/stringr hygiene in `imports.R` itself (not needed for the storm; imports.R is verification-only now).
- Node-shared `/srv/share/users/...` R-library cache (already a documented follow-up).
- Per-method import splits (rejected — per-class chosen).
