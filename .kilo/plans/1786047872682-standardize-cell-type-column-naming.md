# Standardize cell-type column naming to datasets.json ground truth

## Goal

Standardize the cell-type column identifiers across the codebase to match the
datasets.json ground-truth naming (`columns.cell_type_low_res` /
`columns.cell_type_high_res`), replacing the enriched key names
`low_res_ct_col` / `hi_res_ct_col` everywhere.

Confirmed scope decisions (user):
- Rename **also** the `seurat@misc$` slots (notebook + `benchmark_pipeline.R`).
- Keep as-is: generic `ct_col` function params, `sample_col`/`label_col`/`batch_col`
  enriched keys, `Freq_highres`/`Freq_lowres` method labels, `_lowres|_highres`
  feather file suffixes, real data column values (`cells_lowres`, `layer2`,
  `scATOMIC_pred`), `.kilo/plans/*.md` (historical records), AGENTS.md.

## Naming chain (all sites to touch)

```
datasets.json columns.cell_type_low_res / cell_type_high_res   <- ground truth
├── src/utils/datasets_io.R:47-48         enriched keys (R reader)
├── src/datasets_io.py:58-59              enriched keys (Python reader)
│   ├── 1.1.1_benchmark_methods_py.py:17,256-257,264
│   └── notebooks/benchmark_analysis.rmd:83-84,3694
│       └── seurat@misc$low_res_ct_col / hi_res_ct_col (set at rmd:83-84)
│           └── src/5_run_benchmark_methods/benchmark_pipeline.R (~23 sites)
└── docs/ARCHITECTURE.md:14,231,262
```

## Tasks (ordered)

1. **`src/utils/datasets_io.R:47-48`** — rename enriched keys:
   - `low_res_ct_col = columns[["cell_type_low_res"]]` → `cell_type_low_res = columns[["cell_type_low_res"]]`
   - `hi_res_ct_col = columns[["cell_type_high_res"]]` → `cell_type_high_res = columns[["cell_type_high_res"]]`

2. **`src/datasets_io.py:58-59`** — same rename (keep key order mirroring datasets.json).

3. **`src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`**:
   - Line 17 docstring: ``(`low_res_ct_col` / `hi_res_ct_col`)`` → ``(`cell_type_low_res` / `cell_type_high_res`)``
   - Lines 256-257: `entry.get("low_res_ct_col")` → `entry.get("cell_type_low_res")`; `entry.get("hi_res_ct_col")` → `entry.get("cell_type_high_res")` (local vars `lowres_col`/`highres_col` stay)
   - Line 264 print message: `low_res_ct_col is None` → `cell_type_low_res is None`

4. **`notebooks/benchmark_analysis.rmd`**:
   - Line 83: `seurat@misc$low_res_ct_col <- datasets[[ds]][["low_res_ct_col"]]` → `seurat@misc$cell_type_low_res <- datasets[[ds]][["cell_type_low_res"]]`
   - Line 84: same pattern for high res
   - Line 3694: `datasets[[ds]][["hi_res_ct_col"]]` → `datasets[[ds]][["cell_type_high_res"]]`

5. **`src/5_run_benchmark_methods/benchmark_pipeline.R`** — replaceAll on the two exact
   strings (do NOT touch `ct_col = "layer2"`, `Freq_highres`, feather names):
   - `seurat@misc$low_res_ct_col` → `seurat@misc$cell_type_low_res`
   - `seurat@misc$hi_res_ct_col` → `seurat@misc$cell_type_high_res`

6. **`docs/ARCHITECTURE.md`**:
   - Line 14: `label_col, low_res_ct_col, hi_res_ct_col` → `label_col, cell_type_low_res, cell_type_high_res`
   - Line 231: ``(`low_res_ct_col`/`hi_res_ct_col`)`` → ``(`cell_type_low_res`/`cell_type_high_res`)``
   - Line 262: ``MrVI without `low_res_ct_col` skips`` → ``MrVI without `cell_type_low_res` skips``

## Validation (config/syntax checks only — no pipeline runs, per AGENTS.md)

1. Grep for leftovers: `low_res_ct_col|hi_res_ct_col` must match **only** `.kilo/plans/*.md`.
2. R config check (precedent: prior datasets.json plan):
   `pixi run Rscript -e 'source("src/utils/load_all_functions.R"); d <- read_datasets_json(view="benchmark_analysis"); stopifnot(all(c("cell_type_low_res","cell_type_high_res") %in% names(d[[1]]))); stopifnot(!any(c("low_res_ct_col","hi_res_ct_col") %in% names(d[[1]])))'`
3. Python config check:
   `python3 -c "import sys; sys.path.insert(0, \"src\"); from datasets_io import read_datasets_json; e = next(iter(read_datasets_json(view=\"benchmark_analysis\").values())); assert all(k in e for k in (\"cell_type_low_res\", \"cell_type_high_res\")); assert not any(k in e for k in (\"low_res_ct_col\", \"hi_res_ct_col\"))"`
4. Syntax checks: `python3 -m py_compile src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py` and `pixi run Rscript -e 'parse("src/5_run_benchmark_methods/benchmark_pipeline.R")'`.

## Risks / notes

- Breaking change to the enriched-entry API, but the only consumers are the two
  in-repo files updated in the same change; no external consumers exist.
- `1.2_benchmark_methods_py.qmd` (legacy notebook) does not reference the keys — no change.
- Key order in both readers should keep mirroring datasets.json `columns` order
  (`sample`, `label`, `batch`, `cell_type_low_res`, `cell_type_high_res`).
