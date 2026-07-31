# Plan: Shared `read_datasets_json()` in Python (and enriched R version) + `file_names` in datasets.json

## Goal

Implement a Python `read_datasets_json()` analogous to `src/utils/datasets_io.R:1`, usable by `src/preprocess/1.1.1_preprocess.py` and `src/preprocess/_create_combinedpbmc_dataset.py` to replace their ad-hoc `json.load` config reading. Both the Python and R functions return **all fields defined in datasets.json** (all views, all variables) for future-proof reuse. Add a root-level `file_names` field to every dataset entry and re-add the `Zhu` dataset to datasets.json.

## Additional notes
See for more information on dataset.json file paths and NAS connection in:
AGENTS.md:36-36
```
## datasets.json
```


## Confirmed design decisions (from user)

1. **`file_names` at dataset root level** (string or list of strings): the dataset's raw input file(s), independent of views. This lets the combine script pull raw sources (Stephenson, GongSharma, Zhu) even for datasets that have no analysis views (Zhu) and are never preprocessed standalone.
2. **View-level `input_file_name` stays** in existing views (backward compatible; may diverge per view).
3. **Zhu re-added to datasets.json** (restored from git history `c38ff02`, adapted to current entry format) with `"use_for_benchmark": false`, `"use_for_batch_effect": false`, and **no views** (never preprocessed/analyzed standalone; only a raw source for CombinedPBMC).
4. **GongSharma gets no new view** — the `file_names` design makes it unnecessary.
5. Both functions return all dataset-level fields + all views per dataset. Python returns view-less datasets too (combine script needs Zhu); R keeps skipping view-less datasets (its callers always pass a view filter, so behavior is unchanged).

## Current state (verified)

- `datasets.json`: 13 entries; no `Zhu` (removed in `ec5dfc1`); original Zhu entry in `c38ff02` used old format (`file_name`, no `batch`, view without `input_file_name`).
- `data/ZhuH_2023_37379396whole.rds`
- `_create_combinedpbmc_dataset.py` currently crashes against datasets.json (no Zhu entry; GongSharma has no `batch_effect_analysis` view).
- R callers: `benchmark_analysis.rmd:39` (uses `display_name`, `label_col`, `low_res_ct_col`, `hi_res_ct_col`, legacy `ds_name` path building = known TODO) and `batch_effect_analysis.rmd:22` (sparse usage). Both pass a view filter; enriched entries stay backward compatible.
- `load_all_functions.R:8` sources `datasets_io.R`; `docs/ARCHITECTURE.md:9` documents it.

## Tasks

### 1. datasets.json — add `file_names` to all entries + re-add Zhu

- Add `"file_names"` to every existing entry (mirrors the view `input_file_name` today):
  - Adams `AdamsT_2020_32832599.rds`, Bassez `BassezA_2021_33958794whole.rds`, Gongsharma `["Sound_Life_YoungAdult_Male_CMVneg.h5ad", "Sound_Life_YoungAdult_Male_CMVpos.h5ad"]`, Joanito `JoaI_2022_35773407_Nofilt_whole.rds`, Kfoury `Kfoury_2021_34719426.rds`, Kim `KimN_2020_32385277whole.rds`, Lee `LeeA_2021_34836966p_tumor_seurat.rds`, Pelka `PelkaK_2022_34450029whole.rds`, Smillie `SmillieC_2019_31348891.rds`, Stephenson `StephensonE_2021_33879890_preprocessed.rds`, Wu `WuS_2021_34493872.rds`, Zhang `ZhangY_2022_34653365.rds`, CombinedPBMC `combined_pbmc_batch_effect_analysis.h5ad`.
- Re-add Zhu (after Zhang, alphabetical) in the current entry format:
  ```json
  "Zhu": {
    "display_name": "Zhu",
    "file_names": "ZhuH_2023_37379396whole.rds",
    "folder_name": "ZhuH_2023_37379396",
    "tissue": "Blood",
    "normal_tissue": true,
    "use_for_benchmark": false,
    "use_for_batch_effect": false,
    "columns": {
      "sample": "Sample",
      "label": "Sample",
      "cell_type_low_res": null,
      "cell_type_high_res": null
    },
    "meta_cols_keep": ["Sample"],
    "views": {}
  }
  ```

### 2. New `src/datasets_io.py` — Python `read_datasets_json()`

New lightweight module (stdlib `json` only — deliberately NOT in `_preprocess_utils.py`, which imports rpy2 at module level; benchmark worker scripts must be able to import it without R).

```python
def read_datasets_json(path="datasets.json", view=None):
    """Python analog of src/utils/datasets_io.R::read_datasets_json().

    Returns {ds_name: entry} with ALL datasets.json fields (future-proof).
    Entry summary fields (output_file, view_name, input_file, subset_vars)
    refer to the FIRST view matching `view` (with an output_file_name),
    mirroring the R version. `views` holds ALL matching views. Datasets
    without any matching view (e.g. Zhu, view-less raw source) are included
    with summary fields set to None/{} so raw inputs (file_names) are
    discoverable; the R version skips them (its callers use view filters).
    """
```

Entry keys (both languages, identical names):
- Dataset-level: `display_name`, `folder_name`, `file_names`, `tissue`, `normal_tissue`, `use_for_benchmark`, `use_for_batch_effect`, `sample_col`, `label_col`, `batch_col`, `low_res_ct_col`, `hi_res_ct_col`, `meta_cols_keep`
- First-matching-view summary (R-compat): `view_name`, `output_file`, `input_file`, `subset_vars`
- `views`: `{view_name: {"input_file": …, "output_file": …, "subset_vars": …}}` (all views matching filter; skip views without `output_file_name`, like R)

Semantics: no views → included with `views={}`, summary `None`/`{}`; `view` filter applies only to view matching (summary + `views`).

### 3. `src/utils/datasets_io.R` — enrich the R function

Rewrite `read_datasets_json()`:
- Keep existing inclusion rule (skip datasets with no views) and the `break`-equivalent (first matching view with `output_file_name`) for summary fields — existing callers unchanged.
- Keep existing entry keys (`output_file`, `label_col`, `batch_col`, `low_res_ct_col`, `hi_res_ct_col`, `display_name`, `tissue`) exactly as-is.
- Add the new keys listed in Task 2 (use `ds[["file_names"]]`, `ds[["columns"]][["sample"]]`, `ds[["meta_cols_keep"]]`, etc.) and a `views` list of all matching views (each with `view_name`, `input_file_name`, `output_file_name`, `subset_vars`).

### 4. Refactor `src/preprocess/1.1.1_preprocess.py` main()

- Remove `import json`; replace the `json.load` block with `config = read_datasets_json(config_path)` (import from `src.datasets_io`).
- Dataset loop: `sample_col = entry["sample_col"]`; `batch_col = entry["batch_col"] or entry["sample_col"]`; `use_for_batch_effect = entry["use_for_batch_effect"]`; skip datasets with empty `entry["views"]` (replaces `if not views: continue`; also skips view-less Zhu).
- View loop: iterate `entry["views"]`; skip views without `input_file`/`output_file`; `subset_vars = view_info["subset_vars"]`.
- Keep all other logic identical (`is_batch_view`, batch_key, Harmony, HVG, resolution parameters, output-exists skip, error messages).

### 5. Refactor `src/preprocess/_create_combinedpbmc_dataset.py` main()

- Remove `import json`; `config = read_datasets_json(config_path)` (view=None).
- Replace `source_configs` construction: ordered list `["Stephenson", "Gongsharma_cmv_young_males", "Zhu"]`; look up `config[ds_name]` and raise a clear `KeyError` listing missing sources (user adds them to datasets.json).
- Per source: `file_names = entry["file_names"]`; `sample_col = entry["sample_col"]`; `label_col = entry["label_col"]`; `subset_vars = entry["views"].get("batch_effect_analysis", {}).get("subset_vars", {})` (Stephenson's view subsets Status/samples; GongSharma and Zhu have no such view → `{}`).
- `load_and_prepare_source()` takes the entry instead of view_info and loads `file_names` via `load_input()` (already handles str and list).
- Keep all per-source behavior (GongSharma 15-sample RNG pick, `batch`/`cond` assignment, Sample prefixing, common-gene subsetting with `<5000` union fallback, obs trimming, concat order Stephenson→GongSharma→Zhu, output path `data/combined_pbmc_batch_effect_analysis.h5ad`).

### 6. R callers — verify, no functional changes

- `benchmark_analysis.rmd:39` and `batch_effect_analysis.rmd:22`: entries remain backward compatible; no edits required.
- Do NOT touch the legacy `datasets[[ds]][["ds_name"]]` RDS-path construction (lines 75/435; loading mechanism is an open TODO per TODO.md:174) — out of scope.

### 7. Documentation

- Update `docs/ARCHITECTURE.md:9` to mention both readers (`src/utils/datasets_io.R` and new `src/datasets_io.py`), the enriched entry structure, and `file_names`.
- Note in ARCHITECTURE.md: Zhu is a view-less raw source for CombinedPBMC (not preprocessed standalone).

## Validation

- `python3 -c "from src.datasets_io import read_datasets_json; import json; c=read_datasets_json(); assert 'Zhu' in c and c['Zhu']['file_names']=='ZhuH_2023_37379396whole.rds' and c['Zhu']['views']=={}; assert len(c['Stephenson']['views'])==2; assert read_datasets_json(view='benchmark_analysis')['Stephenson']['output_file'].endswith('benchmark_analysis_ECODAprocessed.h5ad')"`
- `python3 -m py_compile src/datasets_io.py src/preprocess/1.1.1_preprocess.py src/preprocess/_create_combinedpbmc_dataset.py`
- R: `pixi run Rscript -e 'source("src/utils/load_all_functions.R"); d <- read_datasets_json(view="benchmark_analysis"); stopifnot(all(c("output_file","label_col","batch_col","low_res_ct_col","hi_res_ct_col","display_name","tissue","file_names","sample_col","meta_cols_keep","views") %in% names(d[[1]])))'` — plus confirm the result is identical to the old version for previously available keys.
- Dry-run `1.1.1_preprocess.py` on a dataset whose output already exists (prints "Already processed", no crash).
- Full combine run (heavy: loads Stephenson rds + 2 GongSharma h5ads + ZhuH_2023_37379396whole.rds via rpy2): confirm `data/combined_pbmc_batch_effect_analysis.h5ad` produced, 3 batches, GongSharma contributes exactly 15 samples, obs limited to Sample/batch/cond. **Check `ZhuH_2023_37379396whole.rds` has a `Sample` obs column and raw counts**

## Out of scope / notes

- Legacy RDS-path loading in the .rmd files (TODO.md:174 h5ad loading task).
- HPC preprocess array (`1.1_run_worker.sh` `jq -r 'keys[]'`) will include Zhu; its worker exits quickly since Zhu has no views (harmless; can exclude later).
- R/Python divergence: R skips view-less datasets, Python includes them (documented in docstrings + ARCHITECTURE.md).
- `1.1.1_preprocess.py` also gets `read_datasets_json` for CombinedPBMC processing — CombinedPBMC entry already has its batch_effect_analysis view, no change needed there.

## Risks

- `Zhu.rds` content unverified (raw vs processed): mitigated by the validation step; if it is the processed object, swap `file_names` to the correct raw file once the user provides it.
- R behavior regression from enriched entries: mitigated by keeping all existing keys/values identical and validating with the Rscript check above.
