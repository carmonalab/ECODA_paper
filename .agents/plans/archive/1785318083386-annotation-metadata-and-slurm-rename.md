# Annotation Metadata Columns & slurm_bash_config.env Rename

## Task A — Implement scATOMIC & Subset Annotation Columns

### Problem

1. `2.2_process_chunk.sh:220` extracts only `c("layer1", "layer2", "layer3")`, dropping HiTME UCell scores and never running scATOMIC.
2. `3_merge_annotations.py` applies no column filtering.
3. The intent (from `TODO.md` line 8) is to **subset to the correct columns** defined in the legacy `Preprocess_datasets.Rmd`, not keep all columns.

### Required Columns (from `Preprocess_datasets.Rmd:502-524`)

```
HiTME_cols_keep: IFN_UCell, HeatShock_UCell, cellCycle.G1S_UCell, cellCycle.G2M_UCell, layer1, layer2, layer3
scATOMIC_cols:   layer_1, layer_2, layer_3, layer_4, layer_5, layer_6, scATOMIC_pred, S.Score, G2M.Score, Phase, classification_confidence
```

### Changes to `src/cell_type_annotation/2.2_process_chunk.sh`

#### 1. Add scATOMIC R library loading

After the existing `library(HiTME)` / `library(Seurat)` block (line ~69), add:

```r
library(scATOMIC)
library(cutoff.scATOMIC)
```

⚠️ Verify scATOMIC is installed in the pixi env (`pixi.toml` only lists `cutoff.scATOMIC`). The package `abelson-lab/scATOMIC` may be auto-pulled as a dependency, or may need to be added to `[tasks.install-scatomic-deps]`.

#### 2. Add NORMAL_TISSUE env var

After the existing env var parsing block (line ~112-121), add:

```r
env_normal_tissue <- Sys.getenv("NORMAL_TISSUE")
```

To `defaults` list (line ~116-121):

```r
normal_tissue = if (env_normal_tissue != "") as.logical(env_normal_tissue) else TRUE
```

#### 3. Define column whitelists

After the defaults parsing block, add:

```r
hitme_cols_keep <- c("IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                     "cellCycle.G2M_UCell", "layer1", "layer2", "layer3")
scatomic_cols <- c("layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                   "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence")
annot_cols <- c(hitme_cols_keep, scatomic_cols)
```

#### 4. Add scATOMIC execution before HiTME in per-sample loop

Before the existing `Run.HiTME()` block (line ~200), add:

```r
### scATOMIC annotation ####
if (is.null(seurat_obj@meta.data[["layer_1"]])) {
  sca_preds <- run_scATOMIC(seurat_obj@assays$RNA$counts)
  sca_results <- create_summary_matrix(
    prediction_list = sca_preds,
    raw_counts = seurat_obj@assays$RNA$counts,
    normal_tissue = args$normal_tissue
  )
  sca_cols <- intersect(scatomic_cols, colnames(sca_results))
  seurat_obj <- AddMetaData(seurat_obj, sca_results[, sca_cols, drop = FALSE])
}
```

Match the retry/timeout pattern from `Preprocess_datasets.Rmd:633-677` with at most 5 retries on the HPC.

#### 5. Update annotation extraction (line 220)

Replace:
```r
annot <- meta[, c("layer1", "layer2", "layer3"), drop = FALSE]
```

With:
```r
keep_cols <- intersect(annot_cols, colnames(meta))
annot <- meta[, keep_cols, drop = FALSE]
```

### Changes to `src/cell_type_annotation/3_merge_annotations.py`

After the join (line 40) and before/after the sentinel check, add column subsetting:

```python
hitme_cols_keep = ["IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
                    "cellCycle.G2M_UCell", "layer1", "layer2", "layer3"]
scatomic_cols = ["layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
                  "scATOMIC_pred", "S.Score", "G2M.Score", "Phase", "classification_confidence"]
annot_cols = set(hitme_cols_keep + scatomic_cols)
# Keep original obs columns plus whitelisted annotation columns
orig_cols = [c for c in adata.obs.columns if c not in annot_cols]
existing_annot = [c for c in annot_cols if c in adata.obs.columns]
adata.obs = adata.obs[orig_cols + existing_annot]
```

### Dependencies / Risks

- **scATOMIC R package**: Must be installed. If missing from pixi env, add to `pixi.toml` under `[tasks.install-scatomic-deps]`:
  ```toml
  remotes::install_github('abelson-lab/scATOMIC')
  ```
- **Rmagic / dlm**: Already installed via `pixi.toml` `install-scatomic-deps`.
- **NORMAL_TISSUE env var**: Must be set per-dataset before running `1_prepare_chunks.sh`. If unset, default to `TRUE` (safe for blood-derived datasets).
- **`normal_tissue` -> `NORMAL_TISSUE`**: Legacy code stored this in `datasets[[ds]]$normal_tissue`. In the new pipeline, it's an env var like `TISSUE_TYPE`. Update any caller scripts (if they exist in the HPC env) to export `NORMAL_TISSUE`. (should be called from datasets.json and passed through the pipeline).

### Validation (SKIP, as creating a debug dataset is not yet implemented. Will be done later. Ignore here.)

1. Run annotation pipeline on a debug dataset with `NORMAL_TISSUE=FALSE` for a tumor dataset.
2. Check feather output columns — should contain only `hitme_cols_keep` + `scatomic_cols` + `cell_barcode` + `sample`.
3. Check final merged `.h5ad` `obs` columns — no extra annotation columns beyond the whitelist.
4. Confirm `obs["layer1"].notna().sum()` still reports correct count.
5. Confirm scATOMIC columns (`layer_1`, `scATOMIC_pred`, etc.) are populated.

---

## Task B — Rename `slurm_bash_config.env` to `slurm_config.sh`

### Rationale
File is a bash script (`#!/bin/bash`, `export`, variable expansion) — `.env` is misleading. Rename to `.sh` per convention.

### Sub-question: Run-once vs. source in all scripts
**Must be sourced in every script.** Each `sbatch --array` task starts a new shell with no inherited environment. Config variables (`PROJECT_ROOT`, `NAS_TARGET_DIR`, `SLURM_PARTITION`, etc.) must be set in every script that uses them. There is no "run once" mechanism for SLURM array jobs.

### Files to modify

#### 1. Rename config file
`git mv src/slurm_bash_config.env src/slurm_config.sh`

Update self-referencing comment in `src/slurm_config.sh` (line 6) to reflect new filename.

#### 2. Update source paths in shell scripts

| File | Line | Change |
|---|---|---|
| `src/cell_type_annotation/1_prepare_chunks.sh` | 13 | `../slurm_bash_config.env` → `../slurm_config.sh` |
| `src/cell_type_annotation/2_submit_hpc_array.sh` | 4 | `../slurm_bash_config.env` → `../slurm_config.sh` |
| `src/cell_type_annotation/2.1_run_worker.sh` | 16 | `../slurm_bash_config.env` → `../slurm_config.sh` |
| `src/cell_type_annotation/2.2_process_chunk.sh` | 13 | `../slurm_bash_config.env` → `../slurm_config.sh` |
| `src/preprocess/1_submit_hpc_array.sh` | 4 | `../slurm_bash_config.env` → `../slurm_config.sh` |
| `src/preprocess/1.1_run_worker.sh` | 16 | `../slurm_bash_config.env` → `../slurm_config.sh` |

#### 3. Update path in 1 R script

| File | Line | Change |
|---|---|---|
| `src/cell_type_annotation/1_prepare_chunks.r` | 4 | `"src", "slurm_bash_config.env"` → `"src", "slurm_config.sh"` |

`readRenviron()` reads `KEY=VALUE` lines — works identically regardless of extension.

#### 4. Update documentation references

| File | Action |
|---|---|
| `TODO.md` | Bulk replace all `slurm_bash_config.env` → `slurm_config.sh` (lines: 6, 9, 23, 42, 69, 84, 189, 197, 302, 342, 358, 409, 445, 451, 496) |
| `AGENTS.md` line 56 | `slurm_bash_config.env` → `slurm_config.sh` |
| `docs/ARCHITECTURE.md` (if references exist) | Check and update |
| `.kilo/plans/1785265434719-pipeline-overhaul-plan.md` | Historical plan record — update or leave as-is |

### Order of operations
1. `git mv src/slurm_bash_config.env src/slurm_config.sh`
2. Update comment in `src/slurm_config.sh` (line 6 path)
3. Update all 6 shell script source paths
4. Update 1 R script path
5. Bulk-replace all documentation references
6. `git grep slurm_bash_config.env` → should return 0 matches

### Validation
1. `git grep slurm_bash_config.env` returns 0 matches.
2. Each updated script parses without syntax error: `bash -n <script>`

---

## TODO.md Update

After implementation:

### Under "Step 3d Follow-up — Completed" (line 1-10)

- [x] Line 8: Mark as done after Task A implementation
- [x] Line 9: Mark as done after Task B implementation
- [x] Line 10: Mark as done — resolved: must be sourced in every script

### All other references

Bulk replace `slurm_bash_config.env` → `slurm_config.sh` across `TODO.md` (15 occurrences).
