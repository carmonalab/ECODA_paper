# Plan: Simplify constants.R `dataset_label_map`

## Summary

Remove the hardcoded `dataset_label_map` from `constants.R` and replace its three usages in `benchmark_analysis.rmd` with a dynamic lookup built from `datasets.json` via the existing `read_datasets_json()` function.

## Context

- `dataset_label_map` is a named character vector in `src/R/constants.R:50-66` mapping dataset short names → display labels (e.g., `"Adams" = "Adams (pulmonary fibrosis)"`)
- Used in `benchmark_analysis.rmd` at 3 locations (lines 1996, 2435, 2759) via `recode(dataset, !!!dataset_label_map)`
- `read_datasets_json()` (in `src/R/datasets_io.R:1`) already reads `display_name` from `datasets.json` for each dataset — this is the same info as `dataset_label_map`
- `batch_effect_analysis.rmd` does NOT use `dataset_label_map`
- One mismatch: `datasets.json` key is `Gongsharma_cmv_young_males`, but `benchmark_analysis.rmd` manually renames it to `GongSharma` on line 105-106 before building the flat dataframe. The `dataset_display_map` must account for this.

## Changes

### 1. `src/R/constants.R` — Remove `dataset_label_map`

Delete lines 50–66 (the entire `dataset_label_map` definition + its comment header).

### 2. `benchmark_analysis.rmd` — Add dynamic display map + replace recodes

#### 2a. Add `dataset_display_map` definition (after line 39)

After `datasets <- read_datasets_json(view = "benchmark_analysis")` on line 39, insert:

```r
dataset_display_map <- sapply(datasets, `[[`, "display_name")
dataset_display_map[["GongSharma"]] <- dataset_display_map[["Gongsharma_cmv_young_males"]]
```

This builds a named vector from JSON (key → display_name), then adds an alias entry for the manually-renamed `GongSharma`.

#### 2b. Replace three `recode` calls

At **line 1996**, replace:
```r
  dataset = recode(dataset, !!!dataset_label_map))
```
with:
```r
  dataset = recode(dataset, !!!dataset_display_map))
```

At **line 2435**, same replacement.

At **line 2759**, same replacement.

## Verification

- Search the codebase for any remaining references to `dataset_label_map` — should be zero after the change.
- Confirm `batch_effect_analysis.rmd` and all other `.R`/`.rmd`/`.qmd` files have no references.
- The three recode sites in `benchmark_analysis.rmd` should produce the same display labels for all datasets that exist in both the old map and `datasets.json`. The only difference is `GongSharma` → `"GongSharma (healthy young males, CMV)"` (from JSON) vs. the old `"GongSharma (healthy, CMV)"` — this is the intentionally updated/improved label in `datasets.json`.
