# R imports refactor: per-file explicit imports instead of centralized imports.R

## Goal

Replace the centralized hard-fail package loader `src/utils/imports.R` (41 packages, `stop()` if any missing) with per-file `load_my_packages()` calls, so each R file documents its own requirements and failures are localized. Motivation: session `ses_022e` — a gloscope-only HPC job died because `HiTME`/`MOFA2` (needed by nothing in that job) were missing; and `preprocess_utils.py` (Python preprocessing) currently attaches all 41 packages, including benchmark methods, via `load_all_functions.R`.

## Core design decisions (confirmed with user)

1. **Per-file imports in `src/utils/*.R`**: each module calls `invisible(load_my_packages(c(...)))` with exactly the packages it needs for its *unqualified* calls. Calls already written as `pkg::fn` stay namespaced (no attach needed — they are their own documentation). This keeps the aggregate missing-package error (all missing packages reported at once, `require()` loop) while making lists per-file.
2. **Per-method imports live in the new pipeline worker** `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R` (the per-method dispatcher), **not** at the top of the shared `benchmark_methods_r.R` / `benchmark_pipeline.R` files. Rationale: in R, attaching only matters at *call time* — sourcing a file executes nothing in it. A gloscope job sources all benchmark files but only needs `GloScope` attached (worker switches on `--method`). The shared benchmark files get header comments listing per-method package requirements instead of attaches for MOFA2/scITD/GloScope.
3. **`load_all_functions.R`** keeps sourcing all modules in the current order (entry-point convenience for notebooks/workers), minus the `imports.R` line. Attach order follows source order; keep dplyr/tidyverse attached *after* Seurat to preserve current masking behavior.
4. **`preprocess_utils.py` slims to `src/utils/seurat_utils.R` only** (it only uses `create_clean_seuratv5_object` + `write_h5ad`). Decouples the Python preprocessing pipeline from benchmark packages entirely.
5. **Notebooks restore explicit `load_my_packages()` calls** in a setup chunk (they were stripped when imports.R was centralized) — this also makes the notebooks self-verifying when rendered. Both `notebooks/benchmark_analysis.rmd` and `notebooks/batch_effect_analysis.rmd` keep `source("src/utils/load_all_functions.R")` for functions.
6. **Delete `src/utils/imports.R`**; move `%notin%` definition to `src/utils/helpers.R` (`notebooks/batch_effect_analysis.rmd:39` already defines its own locally — harmless, can stay).
7. **`load_my_packages` is kept, moved to a tiny new `src/utils/package_utils.R`** (not helpers.R — see ordering gotcha in Risks), with a **corrected error message**: the current "add to pixi.toml as 'r-packagename'" advice is wrong for GitHub/Bioc/CRAN-pinned packages (that is exactly the ses_022e misdiagnosis). New message explains both install paths: conda-available packages → `r-*` dep in pixi.toml `[dependencies]` + `pixi install`; GitHub/Bioc/CRAN-pinned (HiTME, MOFA2, scITD, GloScope, Seurat, anndataR, ...) → `pixi run -e py-cuda13 setup`.
8. **Missing-package prevention = manual `check-r-deps` pixi task** (user chose manual, NOT auto-run in submit scripts): a fresh-session Rscript that sources `load_all_functions.R` (every module's per-file check executes) plus a documented extras vector for packages attached outside modules (GloScope, MOFA2, scITD, HiTME, scATOMIC, cutoff.scATOMIC, ProjecTILs, STACAS, SignatuR, scGate, R.utils, reticulate). Run on the login node before submitting; it would have caught the ses_022e failure pre-submission.
9. **Install mechanism stays the single organized `setup` task in pixi.toml** (user chose single task, not per-scope tasks). Conda deps are impossible for most of these packages (conda can't pin GitHub commits — Seurat `e4cc89238`, anndataR `07612e4f`, HiTME, STACAS, ProjecTILs, EPIC, scECODA — nor CRAN/Bioc versions absent from conda-forge/bioconda: MOFA2 1.20.2, GloScope 2.0.1, scITD, limma 3.68.4). Reorganize the `[tasks.setup]` Rscript into scoped comment sections (base / annotation / benchmark) with per-package install-source comments; no new install tasks.
10. `datasets.json` is **untouched** (ground truth).

Key fact to communicate: repeated `load_my_packages`/`library()` calls across files cost nothing — after the first attach in a session they are no-ops. Per-file imports have zero redundancy cost.

## Task list

### 1. New `src/utils/package_utils.R`

Move `load_my_packages()` from imports.R here, unchanged logic (`require()` loop, `suppressPackageStartupMessages`, aggregate `missing_pkgs` stop), but rewrite the `stop()` message per decision #7 (name both install paths; mention `pixi run -e py-cuda13 setup`).

### 2. `src/utils/helpers.R`

Add `%notin% <- Negate(%in%)` (currently at imports.R:75).

### 3. `src/utils/load_all_functions.R`

- Remove `source("src/utils/imports.R")` (line 7).
- Add `source("src/utils/package_utils.R")` as the FIRST source (before helpers.R), so every module's top-level `load_my_packages()` call resolves.
- Keep the rest, same order. Update header comment: each module now declares its own imports.

### 4. Per-file `load_my_packages()` calls in `src/utils/` modules (top of file)

Starter audit (implementer must verify with `grep -nE "\w+\(" <file>` for unqualified calls and `%>%`):

- `datasets_io.R`: none needed (`jsonlite::fromJSON` namespaced) — verify.
- `constants.R`: none.
- `math_utils.R`: `load_my_packages(c("magrittr"))` (uses `%>%`).
- `scoring_metrics.R`: magrittr if `%>%` present (vegan/cluster/igraph/Matrix/mclust/thisutils are `::`-prefixed).
- `pseudobulk.R`: magrittr if `%>%` present (DESeq2/limma/SummarizedExperiment/MatrixGenerics `::`-prefixed).
- `hvcs.R`: magrittr if `%>%` present (dplyr/tidyr/ggrepel `::`-prefixed).
- `seurat_utils.R`: `load_my_packages(c("Seurat", "anndataR", "arrow", "dplyr", "magrittr"))` — unqualified `Idents/subset/VariableFeatures/write_h5ad/read_h5ad/set_cpu_count/mutate/if_else/%>%`.
- `plotting.R`: `load_my_packages(c("ggplot2", "dplyr", "magrittr"))` (unqualified geoms after `ggplot2::ggplot(...)`, `bind_rows`, `%>%`).

Rule: a file attaches only packages used *unqualified*; `::` calls are already explicit documentation.

### 5. Shared benchmark files — no top-level heavy attaches

- `src/5_run_benchmark_methods/benchmark_methods_r.R` and `benchmark_pipeline.R`: do **not** attach MOFA2/scITD/GloScope. Add header comments listing per-method package requirements (MOFA2: `create_mofa`, `run_mofa`, `prepare_mofa`, `get_default_*`; scITD: `initialize_params`, `make_new_container`, `form_tensor`, `run_tucker_ica`; GloScope: `GloScope::gloscope` + unqualified `gloscopeProp`; plus dplyr/magrittr/tibble/arrow/parallelly/BiocParallel/proxy across both files).

### 6. HPC worker `run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R` (the "new pipeline")

After `method` is validated, add a method-scoped attach block using the helper (aggregate error if a method package is missing):

```r
load_my_packages(c("reticulate", "dplyr", "magrittr"))   # audit + extend with the file's unqualified needs
if (method == "gloscope") load_my_packages(c("GloScope"))
if (method == "mofa")     load_my_packages(c("MOFA2"))
if (method == "scitd")    load_my_packages(c("scITD"))
# pseudobulk: DESeq2/limma are ::-namespaced in pseudobulk.R — no attach needed
```

This is the actual ses_022e fix: a gloscope job now requires only `GloScope`, not `HiTME`/`MOFA2`/`scITD`.

### 7. Other HPC workers (audit + attach what they use unqualified)

- `run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R` (already `library(reticulate)`; convert to `load_my_packages(c("reticulate", ...))`).
- `run_transformation_zeroimp_analysis/1.1.1_run_zeroimp_analysis.R` and `1.1.1_run_transformation_analysis.R` (same treatment).
- All keep `source(...load_all_functions.R)` for functions.

### 8. Annotation worker `src/4_cell_type_annotation/2.1.1_process_chunk.R`

- Add `source(file.path(project_root, "src/utils/package_utils.R"))` BEFORE the existing `source(...seurat_utils.R)` (seurat_utils.R's top-level `load_my_packages` needs the helper defined; the worker sources modules directly, not load_all_functions.R).
- Its own `library(scGate/ProjecTILs/SignatuR/HiTME/...)` calls may stay as-is or convert to `load_my_packages` (implementer's choice; keep behavior identical).

### 9. `src/utils/preprocess_utils.py`

Replace `ro.r('source("src/utils/load_all_functions.R")')` (line 16) with:

```python
ro.r('source("src/utils/package_utils.R")')
ro.r('source("src/utils/seurat_utils.R")')
```

and update the comment block above (it only needs `create_clean_seuratv5_object` + `write_h5ad`; seurat_utils.R now attaches its own deps at source time).

### 10. Notebooks — restore explicit imports in a setup chunk

In `notebooks/benchmark_analysis.rmd` and `notebooks/batch_effect_analysis.rmd`, right after `source("src/utils/load_all_functions.R")`, add `invisible(load_my_packages(c(...)))` derived from a usage audit (grep unqualified calls). Known entries from sampling (implementer must complete): ggplot2, dplyr, tidyr, stringr, forcats, patchwork (`wrap_plots`), funkyheatmap (`funky_heatmap`), pheatmap, plotly, scales, ggrepel, ggtext, mclust (`adjustedRandIndex`), Hotelling, zCompositions, EPIC, Seurat, arrow, jsonlite, writexl (`write_xlsx`), plus whatever the audit finds (rstatix? ggpubr? pROC? zoo? gtools? RColorBrewer? ncdf4? scECODA? reticulate? anndataR?). Prefer explicit tidyverse sub-packages over `library(tidyverse)`. The chunk makes each notebook self-verifying.

### 11. Delete `src/utils/imports.R`

After all consumers are converted. Note: the `%notin%` definition (imports.R:75) moves to helpers.R first (task 2).

### 12. pixi.toml — organized `setup` task + manual `check-r-deps` task

- Reorganize `[tasks.setup]` Rscript into scoped comment sections (base / annotation / benchmark) with per-package install-source comments (GitHub SHA pin vs `install_version` vs conda dep). No new install tasks, no new conda deps.
- Add a manual check task (exact extras vector finalized during implementation; run on login node before submitting arrays — per user decision, NOT wired into submit scripts):

```toml
[tasks.check-r-deps]
cmd = "Rscript -e 'source(\"src/utils/load_all_functions.R\"); load_my_packages(c(\"GloScope\", \"MOFA2\", \"scITD\", \"HiTME\", \"scATOMIC\", \"cutoff.scATOMIC\", \"ProjecTILs\", \"STACAS\", \"SignatuR\", \"scGate\", \"R.utils\", \"reticulate\"))'"
```

### 13. Docs

- `AGENTS.md`: add one line under "R Modules" noting the convention: every `src/utils/*.R` and worker declares its own `load_my_packages()` calls at the top; heavy per-method packages (MOFA2/scITD/GloScope) are attached by the dispatcher worker `1.1.1_run_benchmark_methods_r.R` based on `--method`; `pkg::fn` calls need no attach; `pixi run -e py-cuda13 check-r-deps` verifies the env. (AGENTS.md's "11 utility files loaded by load_all_functions.R" stays true — imports.R is gone, 11 modules remain.)
- Check `docs/ARCHITECTURE.md` for stale imports.R/load_all_functions references (grep found none for imports.R, but verify).

## Dead-weight cleanup (opportunistic, do not touch pixi.toml deps)

Many `my_packages` entries are likely unused anywhere (doParallel, foreach, ggtext, GGally, rstatix, RColorBrewer, forcats, gtools, ncdf4, pROC, Rfast, zoo, HiTME, scECODA...). They simply disappear from all attach lists; do **not** remove them from `pixi.toml`/`setup` (out of scope, and some are installed for future/QC uses). The audit may keep some for notebooks if actually used there.

## Risks / gotchas

- **`load_my_packages` ordering**: modules call it at top-level, so `package_utils.R` must be sourced before any module — in load_all_functions.R (task 3), in the annotation worker (task 8), and in preprocess_utils.py (task 9). Missing this = `could not find function "load_my_packages"` at source time.
- **Unqualified `%>%`** is the most common hidden dep — every file using it needs magrittr in its list (or migrate to `|>`; prefer minimal diff: keep `%>%`).
- **Masking order**: dplyr/tidyverse must attach after Seurat (as today) so `filter/mutate` masks behave unchanged.
- **Annotation worker**: seurat_utils.R's new attaches are idempotent after the worker's own `library()` calls and make seurat_utils.R self-sufficient (fixes a latent fragility: today its unqualified `mutate/if_else` only work if some attached package happens to expose dplyr).
- **`process_scitd_fig`/`process_mofa_fig` etc. are only called on HPC** — notebooks don't need MOFA2/scITD/GloScope attached.
- Do not add attaches to the QC_filtering notebooks (`notebooks/QC_filtering/*.Rmd`) — they are already self-contained; out of scope.

## Validation (user runs; per AGENTS.md agents don't run pipeline scripts)

1. `pixi run -e py-cuda13 check-r-deps` on the login node — must complete without missing-package errors (this is the "how to make sure packages are not missing" step).
2. Local smoke test: fresh `Rscript -e 'source("src/utils/load_all_functions.R")'` and `Rscript -e 'source("src/utils/package_utils.R"); source("src/utils/seurat_utils.R")'` — both succeed without missing-package stops.
3. HPC rerun of the ses_022e repro: `./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods gloscope,mofa,pseudobulk` — all three arrays must COMPLETE and sync to NAS (no HiTME/MOFA2 dependency failure on the gloscope job).
4. Render `notebooks/benchmark_analysis.rmd` and `notebooks/batch_effect_analysis.rmd` (setup chunk must cover all attached-package usage).
5. Python interop: run one preprocessing step using `preprocess_utils.py` (rpy2 path) — must not attempt to load benchmark packages.

## Rollout

All changes land in one commit (single atomic change); no intermediate broken states matter. Suggested commit message: `refactor(R): replace centralized imports.R with per-file package imports`.
