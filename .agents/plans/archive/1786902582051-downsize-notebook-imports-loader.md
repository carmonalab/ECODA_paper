# Downsize the notebook "god import" (imports.R / load_all_functions.R) + fix stale docs

## Goal

The two analysis notebooks (`notebooks/benchmark_analysis.rmd`, `notebooks/batch_effect_analysis.rmd`) still attach all 42 packages through `src/utils/imports.R` (via `src/utils/load_all_functions.R`). Verify every package is actually needed, trim the attach list to ~20 packages, add a full repo-wide env-verification file, and sweep stale documentation that still references the old "42-package" loader.

## Verified findings (from exploration)

### Consumers of the loaders (grep-verified, 2026-08-16)
- `load_all_functions.R` is sourced by **exactly 2 files**: `notebooks/benchmark_analysis.rmd:46` and `notebooks/batch_effect_analysis.rmd:21`.
- `imports.R` is sourced by `load_all_functions.R:7` plus the guarded env-refresh smoke checks `src/utils/bash/setup_env_sbatch.sh:159` and `src/utils/bash/refresh_env.sh:111` (the "canonical env-verification list" role).
- All HPC workers were already disentangled in commit `c477852` (2026-08-11): `src/utils/imports_worker_core.R` (Seurat/reticulate/dplyr), `src/utils/imports_worker_transzeroimp.R` (doParallel/foreach/reticulate/dplyr), `src/utils/load_worker_functions.R` (util files minus `imports.R`/`plotting.R`). **Do not touch the workers.**
- The archived plan `1786190762455-r-per-file-imports-refactor.md` (per-file `load_my_packages()` in every util module, delete `imports.R`) was **abandoned** — its approach would re-god-import inside HPC workers (a top-level attach in `seurat_utils.R` would hit trans/zeroimp workers). Its two useful ideas are incorporated here: a manual `check-r-deps` pixi task and an explicit `requireNamespace` check for packages used namespaced-only.

### Package usage audit (exports-scan of the notebooks' R code chunks + the 11 sourced util files; false positives manually resolved — see appendix for methodology)

| Class | Packages | Action |
|---|---|---|
| **Bare use in notebook-executed code** (keep attached) | ggplot2, ggpubr, ggrepel, ggtext, GGally, scales, Seurat, stringr, tidyr, dplyr, purrr, tibble, forcats, funkyheatmap, RColorBrewer, gtools, scECODA, pheatmap, patchwork, writexl | attach list (20) |
| **Namespaced-only** (no attach; keep in env check) | arrow, jsonlite, zCompositions, limma, DESeq2, mclust, cluster, igraph, vegan, Matrix, proxy, parallelly, BiocParallel, factoextra, Hotelling | env_check.R |
| **Worker/annotation-only** (no attach; keep in env check) | anndataR, reticulate, doParallel, foreach, MOFA2, scITD, EPIC, GloScope, HiTME, scATOMIC, cutoff.scATOMIC, scGate, ProjecTILs, SignatuR, R.utils, robCompositions | env_check.R |
| **Dead in the notebooks** (drop entirely) | plotly, pROC, Rfast, zoo, ncdf4, rstatix, tidyverse (meta; members listed explicitly), mclust (attach only) | remove |

Key resolutions (why borderline hits are NOT real usage): `min_max()`/`clr()` are local functions in `math_utils.R:68`; `highlight()` is local in `helpers.R:26` (helpers.R:26, not plotly); `sym` is re-exported by ggplot2/dplyr; `all_of`/`contains` by dplyr; `%>%` by dplyr; `map` at `benchmark_analysis.rmd:2825` is purrr (mclust dropped, no conflict); `rows`/`var` matched Rfast/pROC exports but no real calls. **writexl is missing from imports.R today** — `benchmark_analysis.rmd:3443` works only because of a notebook-local `library(writexl)`.

### Latent gaps found (fix in this plan)
- DESeq2, vegan, robCompositions + all annotation extras (scATOMIC, scGate, ...) are **not covered by any current env check** — this was the failure mode of session `ses_022e` (missing package discovered at job runtime).
- `pROC::var` currently masks `stats::var` for the bare `var(values, na.rm=TRUE)` call in `hvcs.R:23` (identical output; resolution changes on drop — verify in validation).
- Redundant local `%notin%` definition in `batch_effect_analysis.rmd:39` (loader defines it).

### Stale docs to update (see Task 6)
`AGENTS.md:132,162,166`; `docs/ARCHITECTURE.md:105,439`; `README.md:79`; `src/slurm_config.sh:72` (comment says "imports.R" pins reticulate python — imports.R will no longer attach reticulate); `src/5_run_benchmark_methods/benchmark_hpc_utils.R:6-7` (says "after load_all_functions.R" — workers actually source `load_worker_functions.R`); `src/utils/imports_worker_core.R:8-11` + `src/utils/imports_worker_transzeroimp.R:6-10` ("imports.R remains the canonical env-verification list"); optional: `src/4_cell_type_annotation/2.1.1_process_chunk.R:14-15`, `src/utils/py/preprocess_utils.py:17`. TODO.md: no stale references (checked).

## Design decision (confirmed with user)

**Single downsized notebook loader + new env-check file.** Keep `load_all_functions.R`/`imports.R` names and the loader pattern (notebooks keep `source("src/utils/load_all_functions.R")`); trim `imports.R` to the 20-package attach list; add `src/utils/env_check.R` as the repo-wide `requireNamespace()` verification (attach ∪ namespaced-only ∪ worker-only). No `pkg::` conversion sweep, no per-file `load_my_packages()` in shared utils, no pixi.toml dependency changes (AGENTS.md: never drop defined versions), QC_filtering notebooks out of scope (self-contained, don't source the loader).

## Tasks

### 1. Downsize `src/utils/imports.R`
Replace `my_packages` with the 20-package list **in the current effective attach order** (verified against the present order incl. tidyverse's member attach; dplyr attaches right after tidyr as tidyverse does today; writexl joins at the end):

```
ggplot2, ggpubr, ggrepel, ggtext, GGally, scales, Seurat, stringr, tidyr,
dplyr, purrr, tibble, forcats, funkyheatmap, RColorBrewer, gtools, scECODA,
pheatmap, patchwork, writexl
```

- Keep `load_my_packages()` and `%notin% <- Negate(%in%)` exactly as-is.
- Rewrite the `stop()` message in `load_my_packages` (per archived plan): conda-available packages → `r-*` dep in pixi.toml `[dependencies]` + `pixi install`; GitHub/Bioc/CRAN-pinned (Seurat, anndataR, MOFA2, scITD, HiTME, ...) → `pixi run -e py-cuda13 setup`.
- Update the header comment: attach list for the notebooks only; repo-wide env verification lives in `src/utils/env_check.R`.

### 2. New `src/utils/env_check.R`
- `env_check_packages <- c(` the 20 attach packages + the namespaced-only 15 + the worker/annotation-only 16 listed in the table above (51 total; includes `HiTME` — annotation worker).
- `check_env_packages(pkgs)` using `requireNamespace(pkg, quietly = TRUE)` in a loop (NOT `require` — must not attach), aggregate missing list, fail loudly with the same two-path install message as Task 1.
- Top-level `invisible(check_env_packages(env_check_packages))`.
- Header comment: intended to run on the HPC (login node / setup smoke checks); may report "missing" packages on the macOS env where HPC-only packages are not installed — expected.

### 3. Env-refresh smoke checks (`setup_env_sbatch.sh` ~151-162, `refresh_env.sh` ~102-114)
- Add `source("src/utils/env_check.R")` to the smoke-check Rscript, after the existing `imports.R` + worker-subset sources (keep the existing integrity check).
- Update the success message to something like: `All packages in src/utils/imports.R + worker subsets + env_check.R load OK` (must match Task 6's README update).

### 4. pixi.toml — manual `check-r-deps` task
Add (manual only, NOT wired into submitters — archived-plan decision):
```toml
[tasks.check-r-deps]
cmd = "Rscript -e 'source(\"src/utils/env_check.R\")'"
```
Check how sibling tasks declare their env (e.g. `pixi run -e py-cuda13 ...`) and match it. Do not touch `[dependencies]`/`[tasks.setup]`.

### 5. Notebooks (minimal edits)
- `notebooks/benchmark_analysis.rmd:3443`: remove `library(writexl)` (loader now provides it).
- `notebooks/batch_effect_analysis.rmd:39`: remove the local `%notin%` redefinition (loader provides it).
- Keep both `source("src/utils/load_all_functions.R")` lines unchanged.

### 6. Docs sweep (stale → new wording)
- `AGENTS.md:132` — "the canonical 42-package env-verification list" → "the notebook-loader attach list (~20 pkgs; repo-wide env verification via `src/utils/env_check.R`)".
- `AGENTS.md:162` (env-refresh paragraph) — "full loader package check via `src/utils/imports.R`" → "via `src/utils/env_check.R` (attach ∪ namespaced-only ∪ worker packages)".
- `AGENTS.md:166` (worker self-healing paragraph) — "`imports.R` remains the env-verification list only" → "`imports.R` is the notebook attach list; env verification via `env_check.R`".
- `docs/ARCHITECTURE.md:105` — "the 42-package `load_all_functions.R`" → downsized description ("20-package notebook loader").
- `docs/ARCHITECTURE.md:439` — worker-loader table row: "`imports.R` (canonical 42-package env-verification list; notebooks only)" → update to the new roles of `imports.R` + `env_check.R`.
- `README.md:79` — sync the setup success message with Task 3's new text.
- `src/slurm_config.sh:72` — comment mentions "imports.R" as a reticulate consumer; after Task 1 reticulate is no longer attached there → reword to reference the actual `RETICULATE_PYTHON` consumers (annotation worker `2.1.1_process_chunk.R`, rpy2 python workers).
- `src/5_run_benchmark_methods/benchmark_hpc_utils.R:6-7` — "sourced explicitly after load_all_functions.R" → "after `src/utils/load_worker_functions.R`".
- `src/utils/imports_worker_core.R:8-11` and `src/utils/imports_worker_transzeroimp.R:6-10` — "`src/utils/imports.R` remains the canonical env-verification list" → "`imports.R` is the notebook attach list; env verification via `src/utils/env_check.R`".
- Optional rewording: `src/4_cell_type_annotation/2.1.1_process_chunk.R:14-15`, `src/utils/py/preprocess_utils.py:17`.
- Final grep sweep: `git grep -n "42-package\|42 packages\|imports.R" -- '*.md' '*.R' '*.sh' '*.py'` — no remaining stale references.

### 7. Validation (user runs; per AGENTS.md agents don't run pipeline scripts)
1. Local fresh-session: `pixi run -e default Rscript --vanilla -e 'source("src/utils/load_all_functions.R"); cat("loader OK\n")'` — must succeed (all 20 packages exist in the macOS env).
2. Render both notebooks end-to-end (`notebooks/benchmark_analysis.rmd`, `notebooks/batch_effect_analysis.rmd`) — no missing-function errors; figures unchanged.
3. Verify the pROC→stats `var` resolution change: `pixi run -e default Rscript --vanilla -e 'x <- rnorm(100); stopifnot(identical(unname(stats::var(x, na.rm=TRUE)), unname(stats::var(x, na.rm=TRUE))))'` plus spot-check `get_hvcs()` output on one dataset before/after (optional).
4. HPC login node: `pixi run -e py-cuda13 check-r-deps` (or re-run `setup_env_sbatch.sh` smoke) — must pass, covering the previously unverified packages (DESeq2, vegan, robCompositions, annotation extras).
5. Optional: one `_debug` benchmark array (`--ds_name _debug --methods composition`) to confirm workers are untouched.

## Risks / gotchas
- **Attach-order preservation is critical** for masking: Seurat before dplyr, dplyr right after tidyr (its tidyverse-era position), dplyr's `as_factor`/`all_of` re-exports before forcats/tidyselect — the Task 1 order encodes this. Do not alphabetize.
- `purrr::map` replaces `mclust::map` resolution at `benchmark_analysis.rmd:2825` (previously resolved via tidyverse→purrr) — unchanged behavior; mclust is namespaced-only in `scoring_metrics.R`.
- `pROC::var` → `stats::var` at `hvcs.R:23` — same result; confirm via Task 7.3.
- Do NOT remove any package from `pixi.toml` (installed-but-unused is fine: plotly, zoo, ncdf4, ... stay installed; AGENTS.md rule).
- `env_check.R` intentionally fails on macOS for HPC-only packages — document it as an HPC-side check (Task 2 header).
- The `check-r-deps` task must not run while jobs are active on HPC (it only reads the env — actually safe; no install happens — but keep it as a manual login-node command per archived decision).

## Commit / workflow
Single commit, e.g. `refactor(R): downsize notebook imports.R to 20 packages + add env_check.R; fix stale loader docs`. After implementation: archive this plan to `.kilo/plans/archive/`, stage the changed files, commit, push (per AGENTS.md Task Completion Workflow).

## Appendix — audit methodology (for re-running)
1. Extract R code chunks from both `.rmd` files; strip `#` comments; scan the 11 files sourced by `load_all_functions.R` the same way.
2. For each candidate package: count `pkg::` occurrences; for each genuine export of the package (resolved via `getExportedValue` → `environmentName(environment(fn))`, dropping re-exports like plotly/ggpubr re-exporting dplyr verbs), test word-boundary `name(` call matches.
3. Manually resolve collisions with local functions (`clr`, `min_max`, `highlight` in `helpers.R:26`/`math_utils.R:68`) and re-exports (`sym` via dplyr/ggplot2, `all_of` via dplyr, `%>%` via dplyr).
4. Audit script used for this plan (Rscript, run with `pixi run -e default Rscript --vanilla audit_imports.R` from repo root): embedded below.

```r
pkgs <- c("doParallel","factoextra","foreach","ggplot2","ggpubr","ggrepel",
  "ggtext","GGally","mclust","plotly","scales","anndataR","Seurat","limma",
  "Matrix","HiTME","stringr","tidyr","tidyverse","rstatix","EPIC","reticulate",
  "MOFA2","scITD","GloScope","funkyheatmap","RColorBrewer","forcats","gtools",
  "scECODA","Hotelling","zCompositions","pheatmap","ncdf4","arrow","pROC",
  "Rfast","zoo","parallelly","patchwork","dplyr","jsonlite",
  "purrr","readr","tibble","magrittr","rlang","tidyselect",
  "vegan","DESeq2","cluster","igraph","BiocParallel","proxy","writexl",
  "SummarizedExperiment","robCompositions")
notebook_files <- c("notebooks/benchmark_analysis.rmd","notebooks/batch_effect_analysis.rmd")
util_files <- c("src/utils/datasets_io.R","src/utils/constants.R","src/utils/helpers.R",
  "src/utils/math_utils.R","src/utils/scoring_metrics.R","src/utils/pseudobulk.R",
  "src/utils/hvcs.R","src/utils/seurat_utils.R","src/utils/plotting.R",
  "src/5_run_benchmark_methods/benchmark_methods_r.R",
  "src/5_run_benchmark_methods/benchmark_pipeline.R")
strip_comments <- function(txt) gsub("#.*$", "", txt, perl = TRUE)
extract_r_chunks <- function(f) {
  lines <- readLines(f, warn = FALSE); in_chunk <- FALSE; out <- character(0)
  for (ln in lines) {
    if (grepl("^```\\{r", ln)) { in_chunk <- TRUE; next }
    if (grepl("^```", ln)) { in_chunk <- FALSE; next }
    if (in_chunk) out <- c(out, ln)
  }
  paste(out, collapse = "\n")
}
esc <- function(nm) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", nm)
owning_ns <- function(pkg) {
  ex <- getNamespaceExports(pkg)
  out <- stats::setNames(rep(NA_character_, length(ex)), ex)
  for (nm in ex) {
    v <- tryCatch(getExportedValue(pkg, nm), error = function(e) NULL)
    if (is.function(v)) out[[nm]] <- environmentName(environment(v))
  }
  out
}
scan_text <- function(txt, pkg) {
  txt <- strip_comments(txt)
  ns <- 0L
  m <- gregexpr(paste0("(^|[^.\\w])", pkg, "::"), txt, perl = TRUE)[[1]]
  if (m[1] != -1) ns <- length(m)
  bare <- character(0); load_ok <- requireNamespace(pkg, quietly = TRUE)
  if (load_ok) {
    own <- owning_ns(pkg); keep <- own[own == pkg]
    hits <- vapply(names(keep), function(nm) {
      pat <- paste0("(^|[^A-Za-z0-9._:])", esc(nm), "\\s*\\(")
      grepl(pat, txt, perl = TRUE)
    }, logical(1))
    bare <- names(hits)[hits]
  }
  list(ns = ns, bare = bare, load_ok = load_ok)
}
notebook_txt <- paste(vapply(notebook_files, extract_r_chunks, character(1)), collapse = "\n")
util_txt <- paste(vapply(util_files, function(f) {
  strip_comments(paste(readLines(f, warn = FALSE), collapse = "\n"))
}, character(1)), collapse = "\n")
for (pkg in unique(pkgs)) {
  r_nb <- scan_text(notebook_txt, pkg); r_ut <- scan_text(util_txt, pkg)
  if (!r_nb$load_ok || !r_ut$load_ok) { cat(sprintf("%-18s | NS-LOAD-FAIL\n", pkg)); next }
  if (r_nb$ns == 0 && r_ut$ns == 0 && length(r_nb$bare) == 0 && length(r_ut$bare) == 0) {
    cat(sprintf("%-18s | NONE\n", pkg)); next
  }
  cat(sprintf("%-18s | ns_nb=%d ns_util=%d\n", pkg, r_nb$ns, r_ut$ns))
  if (length(r_nb$bare) > 0) cat(sprintf("    notebook bare: %s\n", paste(r_nb$bare, collapse = ", ")))
  if (length(r_ut$bare) > 0) cat(sprintf("    util     bare: %s\n", paste(r_ut$bare, collapse = ", ")))
}
```
