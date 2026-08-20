# Toggle "old" vs "new" pipeline + select methods/datasets in benchmark_analysis.rmd

## Context & decisions (verified against the current notebook)

- User request: a toggle in `notebooks/benchmark_analysis.rmd` between the **"old" pipeline** (published analysis: current datasets, WITHOUT QOT and PILOT-GM-VAE) and the **"new" pipeline** (all methods and datasets currently in the pipeline), plus the ability to select which methods and datasets the analysis includes.
- **Mechanism (user decision)**: code-level scope parameters in the setup chunk — same pattern as the existing `benchmark_metrics` parameter (notebook line 82). NOT interactive Shiny (keeps the static `rmdformats::downcute` output).
- The notebook currently contains NO QOT / PILOT-GM-VAE anywhere (verified by grep: zero matches) — the R-ingest side (constants.R, `benchmark_methods_r.R`, `benchmark_pipeline.R` `QOT_hvg{i}`/`PILOT-GM-VAE_hvg{i}` blocks) is already implemented and committed. This plan therefore ALSO adds the "new-pipeline" notebook content (methods in all figure lists, recode entries, screening branches, exec-key normalization) — all of it **inert under "old"** because the scope filter removes the rows before any figure sees them.
- **Scope is enforced at the data level**, not per-ggplot: `df_results` is filtered once after flattening, and `result_list[["bmark"]]`/`trans`/`zeroimp` are pruned once after the dataset loop. Every downstream chunk — including those without explicit `filter_methods` lists ("All methods (overview)" line 322, Supp fig 21 loops line 1820/1856, PB-CT comparison line 2435, the Excel export) — is scoped automatically. Per-figure `filter_methods`/`methods` lists are additionally wrapped with scope helpers so *restricted* scopes (user-selected subsets) propagate everywhere.
- **"old" must reproduce today's notebook render exactly** (paper figures are sacred; output PDF filenames unchanged). Defaults (`pipeline_version = "new"`, both scope vectors `NULL`) equal today's full behavior (all non-debug datasets, all methods present in the results).
- Recode maps, `default_methods` regexes and the screening `case_when`/factor levels get the QOT/PILOT-GM-VAE entries **unconditionally**: `dplyr::recode` on an absent key is a no-op, unmatched `default_methods` regexes match nothing, and `facet_grid` drops unused factor levels by default (`drop = TRUE`) — no empty panels in "old" mode.
- `df_results$method` is factorized at line 309-312 (`levs`); the scope block must run BEFORE that so levels are exactly the scoped set.
- `exec_times` comes from the NAS merged feather (all datasets) or the `data/exec_times.rds` cache; it is dataset-scoped once after the cache load. Supp fig 14B (RAM plot, line 2078) currently filters ONLY `!is.na(mem_GB)` — it needs an explicit method-scope filter or excluded methods would leak in (behavior-preserving under default scope).
- Cache semantics: `result_list.rds` (saved per dataset inside `load_hpc_benchmark_results`) keeps excluded datasets' bundles — they are re-pruned in memory every run, never deleted. `data/exec_times.rds` is filtered after load (both fresh and cached paths).
- `analysis_datasets` explicitly allows `_debug` (the `_`-prefix skip applies only when `analysis_datasets` is NULL); `_debug` has no HPC R bundles → `load_hpc_benchmark_results` warns, figures omit them (acceptable for smoke testing).
- Files touched: **`notebooks/benchmark_analysis.rmd` only**. No `src/utils/` changes (helpers live in the notebook; `prepare_benchmark_data`/`build_funky_heatmap` receive already-scoped `filter_methods`).

## Design

### New setup-chunk block (insert at line ~52, replacing the current `_`-filter line)

```r
# ============================================================
# ANALYSIS SCOPE — select which methods + datasets the analysis includes.
#   pipeline_version: "old" = published set (excludes the extended methods
#                     below); "new" = everything currently in the pipeline.
#   analysis_methods:  NULL = all methods found in the loaded results;
#                      or exact result keys, e.g.
#                      c("ECODA_authors_HR", "MrVI_hvg2000", "QOT_hvg2000")
#   analysis_datasets: NULL = all non-debug datasets from datasets.json;
#                      or dataset keys, e.g. c("Adams", "Pelka") (allows "_debug").
# Defaults reproduce the published notebook exactly.
# ============================================================
pipeline_version <- "new"   # "old" | "new"
analysis_methods <- NULL
analysis_datasets <- NULL

# Result-key regex of methods added after the published analysis (excluded
# when pipeline_version == "old"; extend when new methods are added).
extended_method_regex <- "^QOT_|^PILOT-GM-VAE_"
# Dataset keys added after the published analysis (excluded when "old").
extended_datasets <- character(0)

stopifnot(pipeline_version %in% c("old", "new"))
if (!is.null(analysis_datasets)) {
  datasets <- datasets[names(datasets) %in% analysis_datasets]
} else {
  datasets <- datasets[!grepl("^_", names(datasets))]
}
stopifnot(length(datasets) > 0)

# Scope helpers (closed over analysis_methods, resolved after df_results is
# built; only called from figure chunks, which run later).
methods_in_scope <- function(x) intersect(x, analysis_methods)
methods_in_scope_named <- function(m) m[names(m) %in% analysis_methods]
```

`default_methods` (line 88-99): append `"QOT_hvg2000_HR\\b"` and `"PILOT-GM-VAE_hvg2000_HR\\b"` unconditionally (regexes are matched against the recoded screening names; inert when the methods are excluded).

### Scope resolution + dataset pruning

1. **After the dataset loop** (end of the "Process datasets" chunk, after the `gc()` at line 251): prune cached entries of excluded datasets so every later `names(result_list[["bmark"]])` loop is scoped:
   ```r
   result_list[["bmark"]] <- result_list[["bmark"]][names(datasets)]
   result_list[["trans"]] <- result_list[["trans"]][names(datasets)]
   result_list[["zeroimp"]] <- result_list[["zeroimp"]][names(datasets)]
   ```
2. **In the flattening chunk**, after the `mod_score` drop (line 303), BEFORE the `levs`/factor block (line 309):
   ```r
   if (is.null(analysis_methods)) {
     analysis_methods <- unique(as.character(df_results$method))
   }
   if (pipeline_version == "old") {
     analysis_methods <- analysis_methods[
       !grepl(extended_method_regex, analysis_methods)
     ]
   }
   df_results <- df_results[df_results$method %in% analysis_methods, ]
   ```
   (Existing `levs`/factor code then builds levels from the scoped set.)
3. **In the exec-times chunk**, after the cache-load `if/else` (after line 1956), so both the fresh and cached paths are scoped:
   ```r
   # NAS feather covers all benchmark datasets; keep only this run's datasets
   # (n_cells_by_ds keys are display keys, matching exec_times$dataset).
   exec_times <- exec_times[exec_times$dataset %in% names(n_cells_by_ds), ]
   ```

### Per-figure scope wrapping + new-method additions

For every site below: wrap the vector with the scope helper AND add the new-method entries (`"QOT_hvg2000"`, `"PILOT-GM-VAE_hvg2000"` / `= "QOT"`, `= "PILOT-GM-VAE"`). Wrapping is behavior-preserving under the default scope; the additions are inert under "old".

- `filter_methods <- c(...)` → `filter_methods <- methods_in_scope(c(...))`, add the two keys: lines **344** (Supp fig 17, + recode 372-386), **444** (minmax-scaled, + recode 464-476), **1962** (Supp fig 14 exec times, + recode 1987-1995), **2121** (Figure 2A funky, + method_recode 2138-2148), **2205** (rank-based, + recode 2233-2250), **2373** (Supp fig 15, + method_recode 2395-2412), **2470** (presentation, + method_recode 2490-2506), **2530** (Supp fig 1, + method_recode 2547-2557).
- `methods <- list(...)` → after the definition add `methods <- methods_in_scope_named(methods)`, and add the two entries where PILOT is present: **661** (Pelka PCA), **761** (Supp fig 3-13 MDS all datasets), **835** (Pelka features vs dist). Do NOT add entries to 696 (PB+ECODA only), 1150/1189 (annotation-method comparison), 2163 (Figure 2B Bassez — no PILOT) — but still wrap those with `methods_in_scope_named` so restricted scopes apply.
- **Screening chunk** (1247-1340): in the `case_when` (1278-1289) add `grepl("^QOT_", method) ~ "QOT"` and `grepl("^PILOT-GM-VAE_", method) ~ "PILOT-GM-VAE"` (the existing `^PILOT_` branch does NOT match the dashed name — the char after `PILOT` is `-`); extend the `factor(levels = ...)` (1291-1293) with `"QOT"`, `"PILOT-GM-VAE"`; next to the `PILOT_hvg(\d+)` normalization regex (1303) add `stringr::str_replace("QOT_hvg(\\d+)", "QOT_hvg\\1_HR")` and the `PILOT-GM-VAE_hvg(\\d+)` analogue. The unused `filter_methods` at 1247 is dead code (the `filter(method %in% ...)` is commented out) — leave as-is.
- **Supp fig 14B** (2078): add `filter(method %in% analysis_methods)` after the `!is.na(mem_GB)` filter (currently method-unfiltered; behavior-preserving under default scope, excludes QOT/PILOT-GM-VAE under "old").
- **normalize_exec_keys** (1933-1939): add `method = gsub("^(QOT_hvg[0-9]+)_highres$", "\\1", method)` and `method = gsub("^(PILOT-GM-VAE_hvg[0-9]+)_highres$", "\\1", method)` (feather rows carry the `_highres` suffix, bundle keys do not — same as PILOT); update the comment.

Skip the `eval=FALSE` legacy ECODA+PB chunks (1386, 1453) — dead code.

## Task list (ordered)

1. Add the ANALYSIS SCOPE block to the setup chunk (replace the `_`-filter line; keep `datasets`/`dataset_display_map` construction order).
2. Append the two `default_methods` regex entries.
3. Add the bmark/trans/zeroimp pruning after the dataset loop.
4. Add the scope-resolution + df_results filter in the flattening chunk (before `levs`).
5. Add the exec_times dataset filter after the cache-load block.
6. Apply the scope wrapping + new-method entries to all filter_methods/methods/recode sites listed above; add the screening branches/levels/regexes; fix Supp fig 14B; extend normalize_exec_keys.
7. Validation (below).

## Validation

1. **Default render** (`pipeline_version = "new"`, scopes NULL): knit with the local data mirror; QOT + PILOT-GM-VAE appear in df_results and all figures where feathers/bundles exist (missing results → rows absent, figures omit them, no errors).
2. **"old" render** (`pipeline_version = "old"`): no QOT/PILOT-GM-VAE anywhere in df_results, exec tables, or figures; outputs match the current (pre-change) render (same method set, same PDFs).
3. **Restricted scope**: `analysis_datasets <- c("Adams", "Pelka")`, `analysis_methods <- c("ECODA_authors_HR", "MrVI_hvg2000", "QOT_hvg2000")` — every figure (incl. exec-time plots, funky heatmaps, MDS/PCA loops, Excel export) respects the selection; no NA/empty-panel errors.
4. Optional: `analysis_datasets <- "_debug"` smoke check (warnings for missing HPC bundles expected).

## Risks / notes

- Any new method added later must: (a) be covered by the recode maps it appears in, (b) get a `default_methods` regex + screening branch if it is a main method, (c) be added to `extended_method_regex` (or the old/new distinction silently stops covering it) — document this convention in a comment next to `extended_method_regex`.
- Recode of a factor with unmatched keys yields NA with a warning — the scope filter removes rows first, and every scoped method has a recode entry; keep the recode maps as the union of all possible methods.
- `result_list.rds` is not pruned on disk (excluded datasets' bundles remain cached; harmless, re-pruned in memory each run).
- No HPC pipeline runs are needed for validation; notebook knits only. If the local `data/` mirror lacks QOT/PILOT-GM-VAE feathers, "new" renders simply omit them (validation item 1 can be done on a restricted scope).
