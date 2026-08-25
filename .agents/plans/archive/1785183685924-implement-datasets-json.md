# Plan: Replace `Processed_dataset_metadata.R` with `datasets.json`

## Overview

`Processed_dataset_metadata.R` is sourced by `benchmark_analysis.rmd` and `batch_effect_analysis.rmd` to obtain dataset metadata. It is superseded by `datasets.json`. We will:
1. Restructure `datasets.json` so each **view** declares its own `input_file_name` and `output_file_name` (instead of a single dataset-level `file_name`)
2. Add `Gongsharma_cmv_young_males` as a new entry
3. Update all consumers to read from the new structure
4. Add a new R module `src/R/datasets_io.R` to read datasets.json into R
5. Update `preprocess.py` and `copy_data_from_nas_to_hpc_scratch.sh` to use per-view file names
6. Update the Rmd files to reference `output_file_name` (`.h5ad`) from the relevant view
7. Update documentation and TODOs (add h5ad loading and multi-file preprocessing as out-of-scope tasks)
8. Delete `Processed_dataset_metadata.R`

---

## Design Decisions (confirmed with user)

- **Per-view file names**: Each view has `input_file_name` (raw input) and `output_file_name` (preprocessed output). No dataset-level `file_name`.
- **Array `input_file_name`**: `input_file_name` can be a string or array of strings (for multi-file inputs like Gongsharma).
- **Stephenson correction**: The correct input filename is `StephensonE_2021_33879890_preprocessed.rds` (already correct in JSON; old R file's `MitchelJ_2023_...` was outdated).
- **No `processed_file_stem`**: Removed from plan.
- **h5ad loading mechanism**: Out of scope — Rmd files will reference `.h5ad` output files but loading will be implemented in a separate task.
- **Multi-file preprocessing (Gongsharma → preprocess.py)**: Out of scope — add to TODO.md.

---

## Step-by-step tasks

### Task 1: Restructure `datasets.json` with per-view file names

**Remove** the `file_name` field from every dataset-level entry. **Add** `input_file_name` and `output_file_name` to every view.

**Adams** (only has `benchmark_analysis`):
```json
"Adams": {
  "display_name": "Adams (pulmonary fibrosis)",
  "folder_name": "AdamsT_2020_32832599",
  "tissue": "Lung",
  "normal_tissue": true,
  "use_for_benchmark": true,
  "use_for_batch_effect": false,
  "columns": { ... unchanged ... },
  "views": {
    "benchmark_analysis": {
      "input_file_name": "AdamsT_2020_32832599.rds",
      "output_file_name": "AdamsT_2020_32832599_benchmark_analysis_ECODAprocessed.h5ad",
      "subset_vars": { ... unchanged ... }
    }
  }
}
```

**Stephenson** (has both `benchmark_analysis` and `batch_effect_analysis`):
```json
"Stephenson": {
  ...
  "views": {
    "benchmark_analysis": {
      "input_file_name": "StephensonE_2021_33879890_preprocessed.rds",
      "output_file_name": "StephensonE_2021_33879890_preprocessed_benchmark_analysis_ECODAprocessed.h5ad",
      "subset_vars": { ... unchanged ... }
    },
    "batch_effect_analysis": {
      "input_file_name": "StephensonE_2021_33879890_preprocessed.rds",
      "output_file_name": "StephensonE_2021_33879890_preprocessed_batch_effect_analysis_ECODAprocessed.h5ad",
      "subset_vars": { ... unchanged ... }
    }
  }
}
```

**Gongsharma** (umbrella entry, not directly used in analysis — set `use_for_benchmark` and `use_for_batch_effect` to `false`):
```json
"Gongsharma": {
  ...
  "use_for_benchmark": false,
  "use_for_batch_effect": false,
  "views": {
    "benchmark_analysis": {
      "input_file_name": [
        "Sound_Life_YoungAdult_Male_CMVneg.h5ad",
        "Sound_Life_YoungAdult_Male_CMVpos.h5ad"
      ],
      "output_file_name": "Gongsharma_benchmark_analysis_ECODAprocessed.h5ad",
      "subset_vars": {}
    },
    "batch_effect_analysis": {
      "input_file_name": [
        "Sound_Life_YoungAdult_Male_CMVneg.h5ad",
        "Sound_Life_YoungAdult_Male_CMVpos.h5ad"
      ],
      "output_file_name": "Gongsharma_batch_effect_analysis_ECODAprocessed.h5ad",
      "subset_vars": {}
    }
  }
}
```

**Gongsharma_cmv_young_males** (new entry, the actual dataset used):
```json
"Gongsharma_cmv_young_males": {
  "display_name": "GongSharma (healthy young males, CMV)",
  "folder_name": "GongSharma_2024_PrePrintTBD",
  "tissue": "Blood",
  "normal_tissue": true,
  "use_for_benchmark": true,
  "use_for_batch_effect": true,
  "columns": {
    "sample": "Sample",
    "label": "subject.cmv",
    "cell_type_low_res": "AIFI_L1",
    "cell_type_high_res": "AIFI_L3"
  },
  "views": {
    "benchmark_analysis": {
      "input_file_name": "Gongsharma_cmv_young_males_ECODAprocessed.rds",
      "output_file_name": "Gongsharma_cmv_young_males_benchmark_analysis_ECODAprocessed.h5ad",
      "subset_vars": {}
    },
    "batch_effect_analysis": {
      "input_file_name": "Gongsharma_cmv_young_males_ECODAprocessed.rds",
      "output_file_name": "Gongsharma_cmv_young_males_batch_effect_analysis_ECODAprocessed.h5ad",
      "subset_vars": {}
    }
  }
}
```

Apply the same pattern to all other datasets (Bassez, Joanito, Kfoury, Kim, Lee, Pelka, Smillie, Wu, Zhang, Zhu). `output_file_name` follows the convention `{input_file_stem}_{view_name}_ECODAprocessed.h5ad` where `input_file_stem` is the `input_file_name` without extension (or the dataset key for multi-file inputs).

### Task 2: Create `src/R/datasets_io.R`

New R module with a single exported function `read_datasets_json()`:

```r
read_datasets_json(path = "datasets.json", view = NULL) -> list
```

Behavior:
- Reads `datasets.json` via `jsonlite::fromJSON(simplifyVector = FALSE)`
- If `view` is specified, filters to datasets where:
  - The dataset has the requested view
  - The view has an `output_file_name`
  - (`use_for_benchmark`/`use_for_batch_effect` flag could also be checked, but simpler to just check view existence)
- For each matching dataset, returns a flat list with:
  - `output_file`: the view's `output_file_name`
  - `label_col`: `columns.label`
  - `low_res_ct_col`: `columns.cell_type_low_res`
  - `hi_res_ct_col`: `columns.cell_type_high_res`
  - `display_name`: `display_name`
  - `tissue`: `tissue`
- Returns a named list (names = dataset keys)

Add `jsonlite` to `my_packages` in `src/R/imports.R`.
Add `source("src/R/datasets_io.R")` to `src/R/load_all_functions.R` (after imports.R, before constants.R).

### Task 3: Update `preprocess.py` to use per-view file names

Currently `preprocess.py`:
```python
for ds_name, ds_info in config.items():
    input_filename = ds_info.get("file_name")
    if not input_filename:
        continue
    # ... convert, then iterate views ...
    for view_name, view_info in views.items():
        processed_file_path = f"{file_name_no_ext}_{view_name}_ECODAprocessed.h5ad"
```

Update to:
```python
for ds_name, ds_info in config.items():
    views = ds_info.get("views", {})
    for view_name, view_info in views.items():
        input_filename = view_info.get("input_file_name")
        output_filename = view_info.get("output_file_name")
        if not input_filename or not output_filename:
            continue
        # Handle array input_file_name (Gongsharma case — for now, skip arrays with a message)
        if isinstance(input_filename, list):
            print(f"Skipping {ds_name}/{view_name}: array input_file_name not yet supported")
            continue
        # Process single input file -> output file
        ...
```

This is a transitional implementation — multi-file inputs (Gongsharma) are skipped until the `preprocess_gongsharma.qmd → preprocess.py` adaptation (Task 9 in TODO.md).

### Task 4: Update `src/bash/copy_data_from_nas_to_hpc_scratch.sh`

Currently the script uses `jq` to extract `folder_name` and `file_name` at the dataset level. With per-view file names, it needs to extract `folder_name` (dataset level) and `input_file_name` (per view) and collect unique files.

Update the jq query to extract all unique `input_file_name` references across all views, paired with their dataset's `folder_name` for the NAS path. Handle both string and array `input_file_name`.

```bash
# Extract all (folder_name, input_file_name) pairs from all views
jq -r '.datasets | to_entries[] | 
  .value as $ds | 
  ($ds.folder_name) as $folder |
  $ds.views | to_entries[] | 
  .value.input_file_name | 
  if type == "array" then .[] else . end |
  "\($folder) \(.)"' "${DATASETS_JSON_FILE}"
```

### Task 5: Update `benchmark_analysis.rmd`

**Line 39**: Replace `source("Processed_dataset_metadata.R")` with:
```r
datasets <- read_datasets_json(view = "benchmark_analysis")
```

**Lines 67 (markdown note)**: Update reference from `"Processed_dataset_metadata.R"` to `"datasets.json"`.

**Lines 70-84 (processing loop)**: Update to use `output_file` from the helper instead of constructing RDS paths:
```r
# Load the preprocessed h5ad file
# h5ad_file <- file.path(path_data, datasets[[ds]]$output_file)
# Loading mechanism TBD — see TODO.md task for h5ad→R integration
seurat <- readRDS(file.path(path_data, paste0(datasets[[ds]][["ds_name"]], seurat_file_ending)))
```

Keep the existing `readRDS` call as a **temporary fallback** but add a comment noting it will be replaced. The `ds_name` field is kept for backward compatibility with the fallback, derived from the view's `output_file_name` by stripping `_{view_name}_ECODAprocessed.h5ad`.

Wait — the helper should still provide `ds_name` for the fallback. Let me adjust the helper to include it.

**Lines 1826-1835 (exec_times loop)**: Remove dead code (commented-out GongSharma keys that don't exist).

**Line 73**: Keep the `Site == "Ncl"` Stephenson subset with a comment that this is now a safety net (handled at preprocessing).

### Task 6: Update `batch_effect_analysis.rmd`

**Line 22**: Replace `source("Processed_dataset_metadata.R")` with:
```r
datasets <- read_datasets_json(view = "batch_effect_analysis")
```

**Lines 430-434 (Stephenson in Combined PBMC)**: Update to use `output_file`:
```r
ds <- "Stephenson"
# h5ad_file <- file.path(path_data, datasets[[ds]]$output_file)
seurat <- readRDS(file.path(path_data, paste0(datasets[[ds]][["ds_name"]], seurat_file_ending)))
```

**Lines 468-472 (GongSharma in Combined PBMC)**: The key `Gongsharma_cmv_young_males` now exists in JSON. Use the view's `output_file`:
```r
ds <- "Gongsharma_cmv_young_males"
# h5ad_file <- file.path(path_data, datasets[[ds]]$output_file)
seurat <- readRDS(file.path(path_data, paste0(datasets[[ds]][["ds_name"]], seurat_file_ending)))
```

### Task 7: Update documentation

**`docs/ARCHITECTURE.md`** (lines 113-124): Replace LAYER 6 section:

```markdown
### LAYER 6: Configuration & Data
```
datasets.json (read by src/R/datasets_io.R::read_datasets_json())
└── per-view entries
    ├── "Adams" → benchmark_analysis view
    │   ├── output_file: AdamsT_2020_32832599_benchmark_analysis_ECODAprocessed.h5ad
    │   ├── label_col, low_res_ct_col, hi_res_ct_col (from dataset-level columns)
    ├── ...
```
Each view declares input_file_name (raw input) and output_file_name (preprocessed output).
Dataset-level metadata (columns, display_name, tissue) is shared across views.
```

**`README.md`**: Line 38 already mentions `datasets.json`. Update if the description is stale (it says "sample/label columns, subsetting rules, batch information" — still accurate). No change needed.

**`AGENTS.md`**: Lines 28-29 say "This acts as ground truth for the datasets evaluated in this study." Still accurate. Lines 53 mentions datasets.json driving preprocessing — still accurate. No change needed.

### Task 8: Update `TODO.md`

Replace the "Implement new `datasets.json`" section (lines 1-4) to mark it done. Add the following new tasks:

```
# Implement h5ad loading in R (for preprocessed .h5ad files produced by preprocess.py)
- Evaluate and implement the best approach to load h5ad files in benchmark_analysis.rmd and batch_effect_analysis.rmd:
  - anndataR (R package for native AnnData loading)
  - reticulate (Python anndata via rpy2)
  - Convert h5ad to Seurat object (since some benchmarked R methods may need Seurat objects)
- Update the loading code in both .rmd files once the approach is chosen

# Integrate multi-file preprocessing (e.g. Gongsharma h5ad files) into preprocess.py
- Move the downsample_by_group() logic from preprocess_gongsharma.qmd into preprocess.py
- Support array input_file_name in datasets.json views
- This may require additional fields in datasets.json (out of scope for the initial per-view migration)
```

### Task 9: Delete `Processed_dataset_metadata.R`

Remove the file.

### Task 10: Verify

- Run `python src/py/preprocess.py --help` (no crash from JSON parsing)
- Run `jq` query on updated datasets.json to confirm valid JSON
- Source `src/R/load_all_functions.R` in R console, call `read_datasets_json()` and inspect result
- Verify `names(datasets)` includes only datasets with the requested view (e.g. `benchmark_analysis` excludes `Gongsharma` umbrella entry)
- Verify `datasets[["Stephenson"]]$output_file` resolves to the correct h5ad path
- Verify `datasets[["Gongsharma_cmv_young_males"]]$output_file` resolves to the correct h5ad path

---

## Edge cases & risks

| Risk | Mitigation |
|---|---|
| `Gongsharma` umbrella entry appears in benchmark loop | Set `use_for_benchmark=false` in its entry; `read_datasets_json(view="benchmark_analysis")` only includes datasets with that view |
| Stephenson RDS has old filename on disk (`MitchelJ_2023_...`) while JSON uses new name (`StephensonE_2021_...`) | The Rmd files keep a temporary `readRDS` fallback with the old naming; h5ad loading task will resolve properly |
| `jsonlite::fromJSON(simplifyVector=FALSE)` required to prevent auto-conversion of nested lists | Use this explicitly in the helper |
| Array `input_file_name` in Gongsharma views — preprocess.py can't handle yet | preprocess.py skips array inputs with a message; marked as out-of-scope TODO |
| Rmd files will fail to load h5ad files until loading mechanism is implemented | Expected — the h5ad loading is a separate task. Rmd files retain the old RDS fallback during transition. |
