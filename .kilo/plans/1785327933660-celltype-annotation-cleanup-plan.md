# Plan: Cleanup & Simplify Cell Type Annotation Pipeline

## Problems Found

### 1. `config_helper.R` does not exist in project root
- Both `1_prepare_chunks.r` and `2.2_process_chunk.sh` `source()` it, but it only exists in `DEPRECATED_LEGACY_CODE/`
- Runtime would fail at this line
- The old version used `readRenviron("config.env")` (deleted), `pipeline_config.json` (gone), and generated paths under `~/DS_NAME/` (old pattern, not `HPC_SCRATCH_DIR`)

### 2. renv is deprecated; pixi handles the environment
- `2.2_process_chunk.sh` line 35: `RENV_CONFIG_EXTERNAL_LIBRARIES`
- `2.2_process_chunk.sh` R-code lines 48-53: `local_renv_lib` path munging + `.libPaths()` manipulation
- `2.2_process_chunk.sh` R-code line 46: `Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "false")`
- `1_prepare_chunks.sh` lines 73-74: `RENV_CONFIG_EXTERNAL_LIBRARIES` + `R_LIBS_SITE`
- `1_prepare_chunks.r` lines 6-9: renv library path munging
- All of these are no-ops with pixi (`pixi run Rscript` already sets the correct `R_LIBS_SITE`)

### 3. 260 lines of R embedded as bash heredoc in `2.2_process_chunk.sh`
- Hard to debug, edit, lint, or maintain
- No syntax highlighting, no `R CMD check`, no IDE support
- The heredoc argument `"chunk_file__${CHUNK_FILE}"` after `<<'RS'` is unusual bash

### 4. `1_prepare_chunks.sh` uses undefined variables
- `HOME_DATA_DIR`, `NAS_DATA_DIR`, `HOME_REF_DIR`, `NAS_REF_DIR`, `GENE_REF_FILE`, `GENE_REF_URL` are all referenced but not defined in `slurm_config.sh`

### 5. `1_prepare_chunks.r` has broken `readRenviron()` call
- Line 4: `readRenviron(file.path(project_root, "src", "slurm_config.sh"))` — `readRenviron` expects `NAME=VALUE` lines, cannot parse `export NAME=VALUE` from bash scripts
- Also sources `config_helper.R` which doesn't exist in project root

### 6. Duplicated path `PIXI_R_LIB` across 3 scripts
- `2.2_process_chunk.sh` line 35, `1_prepare_chunks.sh` line 70, `DEPRECATED_LEGACY_CODE/1_prepare_chunks.sh` line 70 — all compute `PROJECT_ROOT/.pixi/envs/default/lib/R/library`

### 7. `2_submit_hpc_array.sh` references `DS_NAME` but does not set it
- The caller must export `DS_NAME`; this is fragile

---

## Proposed Changes

### Task A: Update `slurm_config.sh` — add missing shared variables

Add the following to `src/slurm_config.sh`:

```bash
# Pixi R library path (used by cell type annotation scripts)
PIXI_R_LIB="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"

# Reference atlas paths (cell type annotation)
NAS_REF_DIR="${NAS_PREFIX}/DataCollections/reference_atlases/sketched_200ct/"
HOME_REF_DIR="${HOME}/reference_atlases/sketched_200ct/"

# Gene reference (cell type annotation) (-> shouldn't this be "${PROJECT_ROOT}/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz" ?)
# Wait, actually, the GENE_REF_FILE is not needed anymore. It was only used for gene standardization with STACAS, which is now implemented in `src/preprocess/1.2_preprocess.py` which is run before cell type annotation. So can GENE_REF_FILE be removed?
# Move script to download the GENE_REF_FILE from GENE_REF_URL README.md?
GENE_REF_FILE="${PROJECT_ROOT}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
GENE_REF_URL="https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"
```

**Rationale**: These variables are used by multiple scripts; centralizing them eliminates duplication and prevents "undefined variable" runtime errors.

---

### Task B: Create `config_helper.R` in project root

Write `config_helper.R` to project root. It must:

- **NOT** call `readRenviron()` — rely on env vars already exported by bash wrappers
- **NOT** use `renv` path manipulation
- Derive paths from `Sys.getenv("HPC_SCRATCH_DIR")` and `Sys.getenv("DS_NAME")` or project root paths (depending on what's appropriate for HPC scratch)
- Keep the same interface (`get_pipeline_config()` returning a list with `path_data`, `path_output`, `path_output_chunks`, `path_ref`, etc.) so the consuming R code doesn't need internal changes

Key path mappings (new pipeline):

| Old path (`~/DS_NAME/...`) | New path |
|---|---|
| `path_data` | `HPC_SCRATCH_DIR/data` |
| `path_output` | `SCRATCH_OUTPUT_DIR` |
| `path_output_chunks` | `SCRATCH_OUTPUT_DIR/chunks` |
| `path_ref` | `HOME_REF_DIR` (from slurm_config.sh) |
| `gene_ref` | `GENE_REF_FILE` (from slurm_config.sh) |

---

### Task C: Extract R code from `2.2_process_chunk.sh` into `2.2_process_chunk.R`

1. Create `src/cell_type_annotation/2.2_process_chunk.R` with:
   - The `get_sample_seurat_obj()` helper
   - The argument parsing block
   - The scGate/ProjecTILs ref map loading
   - The processing loop (annotations_list loop)
   - The feather output writing
   - NO renv bootstrap code
   - NO `Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED)`
   - NO `local_renv_lib` in `.libPaths()`
   - Keep the `library()` calls (scGate, ProjecTILs, HiTME, Seurat, arrow, reticulate, etc.)

2. Rewrite `2.2_process_chunk.sh` as a thin wrapper:
   ```bash
   #!/bin/bash
   set -euo pipefail
   TASK_ID="$1"
   export PROJECT_ROOT="$2"
   CHUNK_FILE="$3"
   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
   source "${SCRIPT_DIR}/../slurm_config.sh"
   cd "${PROJECT_ROOT}"
   # Use pixi directly — no renv vars needed
   "${HOME}/.pixi/bin/pixi" run Rscript --vanilla \
     "${SCRIPT_DIR}/2.2_process_chunk.R" \
     "chunk_file__${CHUNK_FILE}"
   ```

---

### Task D: Remove renv remnants from remaining scripts

1. **`1_prepare_chunks.sh`** (lines 70-74): Replace:
   ```bash
   PIXI_R_LIB="${PROJECT_ROOT}/.pixi/envs/default/lib/R/library"
   export RENV_CONFIG_EXTERNAL_LIBRARIES="${PIXI_R_LIB}"
   export R_LIBS_SITE="${PIXI_R_LIB}:${R_LIBS_SITE:-}"
   ```
   With:
   ```bash
   export R_LIBS_SITE="${PIXI_R_LIB}:${R_LIBS_SITE:-}"
   ```
   (Where `PIXI_R_LIB` now comes from `slurm_config.sh`.)

2. **`1_prepare_chunks.r`** (lines 6-10): Remove renv library path resolution block. The `.libPaths()` manipulation is handled by pixi.

3. **`2.2_process_chunk.R`**: No renv bootstrap block at all.

---

### Task E: Fix `1_prepare_chunks.r`

1. Line 4: Replace `readRenviron(file.path(project_root, "src", "slurm_config.sh"))` with direct `Sys.getenv()` calls for the needed variables. Since the bash wrapper (`1_prepare_chunks.sh`) already sources `slurm_config.sh` and exports the variables before calling `srun`, the R script can just use `Sys.getenv("HPC_SCRATCH_DIR")`, `Sys.getenv("DS_NAME")`, etc.

2. Line 40: Check if `config_helper.R` now exists in project root (Task B), and if the working directory is `PROJECT_ROOT`, then `source("config_helper.R")` should work.

---

### Task F: Fix undefined variables in `1_prepare_chunks.sh`

The variables `HOME_DATA_DIR`, `NAS_DATA_DIR` are used for rsyncing raw data. In the new pipeline, this script should copy from the preprocessed data location (`HPC_SCRATCH_DIR`) instead of from raw NAS paths. However, since the cell type annotation pipeline may need the original h5ad files (not just preprocessed output), this needs clarification.

**Proposed approach**: Update the rsync source to use the new staging approach:
- Source raw h5ad from `NAS_SC_DIR` (already in `slurm_config.sh`) — this is the standardized single-cell datasets
- Destination: `HPC_SCRATCH_DIR` (already in `slurm_config.sh`)

So `HOME_DATA_DIR` → `HPC_SCRATCH_DIR` and `NAS_DATA_DIR` is replaced by constructing the path from `NAS_SC_DIR/${DS_NAME}/data/`.

Similarly for ref maps: `HOME_REF_DIR` → `HOME_REF_DIR` (now defined in `slurm_config.sh`), `NAS_REF_DIR` → `NAS_REF_DIR` (now defined).

---

### Task G: Clean up `2_submit_hpc_array.sh`

No major changes needed aside from:
- Remove the `pipeline_config.json` dependency (line 15-16) — replace with env-var-based config check or remove test mode flag logic
- The `DS_NAME` variable: document that it must be exported before calling, or read from datasets.json like the preprocessing submit script does

---

## Order of Execution

1. **Task A**: Update `slurm_config.sh` — add missing variables (foundational)
2. **Task B**: Write `config_helper.R` to project root
3. **Task C**: Extract `2.2_process_chunk.R` + rewrite `2.2_process_chunk.sh`
4. **Task D**: Remove renv from `1_prepare_chunks.sh` and `1_prepare_chunks.r`
5. **Task E**: Fix `1_prepare_chunks.r` (`readRenviron` → `Sys.getenv`)
6. **Task F**: Fix undefined variables in `1_prepare_chunks.sh`
7. **Task G**: Clean up `2_submit_hpc_array.sh`
8. **Delete** `DEPRECATED_LEGACY_CODE/config_helper.R` if the new one is complete (keep the rest of DEPRECATED_LEGACY_CODE for reference)
9. **Update AGENTS.md** to note the changes

## Validation

For each task:
- **Shell syntax**: `bash -n <file>` for all modified `.sh` files
- **R syntax**: `R -f <file>` for `.R` files (or parse check via `parse(text = readLines(...))`)
- **Consistency**: Verify that all env vars referenced in `.R` files are exported by the calling `.sh` wrapper
- **Integration**: The `2.2_process_chunk.R` script, when called via `pixi run Rscript`, should load all required R packages without `.libPaths()` manipulation

## Risk Register

| Risk | Mitigation |
|---|---|
| `config_helper.R` path mappings differ from what downstream R code expects | Keep the same list keys (`path_data`, `path_output`, `path_ref`, etc.) and update only the values |
| `1_prepare_chunks.sh` still copies raw data from NAS instead of preprocessed data | Confirm the data staging flow: cell type annotation reads from preprocessed h5ads, which are already staged by the preprocessing pipeline |
| Removing `readRenviron()` may break scripts that depend on specific env vars being set | `slurm_config.sh` is sourced by the bash wrapper before calling R, so all exports are already in the environment; explicit `Sys.getenv()` is more reliable |
