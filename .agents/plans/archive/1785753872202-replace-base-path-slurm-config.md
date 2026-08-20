# Plan: Replace `base_path`/`project_root` with slurm_config.sh env vars (per-dataset layout)

Implements TODO.md "## Replace base_path" (lines 3-12). User decisions:
- Per-dataset folder structure is the target convention (confirmed already in place for preprocess).
- Also align the cell type annotation pipeline paths so it becomes functional end-to-end (not just minimal renames).
- `_create_combinedpbmc_dataset.py` should be made HPC-capable (read raw sources from HPC scratch, write combined h5ad into scratch).

## Verified current state

- **Preprocess already per-dataset** (no changes needed):
  - `1_submit_hpc_array.sh:30-36` stages `NAS_SC_DIR/<folder>/output/<file>` → `${HPC_SCRATCH_DIR}/${KEY}/data`
  - `1.1_run_worker.sh:31-38` reads `${HPC_SCRATCH_DIR}/${DS_NAME}/data` (`--base_path`), writes `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` (`--output_dir`), both absolute
  - Post-array rsync `${SCRATCH_OUTPUT_DIR}/` → `${NAS_TARGET_DIR}/output/` (line 85)
- **Annotation pipeline path seams (currently non-functional)**:
  - `config_helper.R:28` `path_data = ${HPC_SCRATCH_DIR}/data` (flat) vs preprocessed outputs in `${SCRATCH_OUTPUT_DIR}/${DS_NAME}/`
  - `config_helper.R:32` `path_output_chunks = ${SCRATCH_OUTPUT_DIR}/chunks` (flat) vs `2_submit_hpc_array.sh:17` `HOME_CHUNKS_DIR=${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks`
  - `slurm_config.sh` exports only `PROJECT_ROOT` + `NAS_PREFIX`; `config_helper.R` reads `HPC_SCRATCH_DIR`/`SCRATCH_OUTPUT_DIR`/`DS_NAME` via `Sys.getenv()`, which never propagates through `srun` (1_prepare_chunks.sh) → script stops immediately
  - `SAMPLE_COLNAME` exported nowhere; R scripts default to `"sample"` but `1.1.1_preprocess.py:210-213` standardizes to `"Sample"`
  - `1.1_prepare_chunks.r:1` and `2.1.1.1_process_chunk.R:8` use `project_root <- getwd()` (works only because bash `cd`s to PROJECT_ROOT; fragile, has TODO comments)
- **`1.1.1_preprocess.py:168-174`**: `project_root = Path(__file__).resolve().parents[2]`, then `project_root / base_path` — only correct for absolute CLI args because pathlib collapses absolute right operands (implicit, fragile).
- **Leftover `base_path`/`project_root`** (checked via grep): `DEPRECATED_LEGACY_CODE/` (dead), `src/gene_utils.py`, `src/datasets_io.py` (repo-relative resources: `aux/`, `datasets.json` — keep), `src/preprocess/_create_combinedpbmc_dataset.py` (handled below), `.kilo/plans/*` (historical).
- **Benchmark**: TODO says decision deferred to user; runs locally or on HPC later → out of scope here (keep TODO note).

## Changes

### 1. `src/slurm_config.sh` — export core env vars

Add `export` to the following (required so `srun`/`sbatch` children and `Sys.getenv()`/`os.environ` see them):
`DATASETS_JSON_FILE`, `NAS_SC_DIR`, `NAS_TARGET_DIR`, `NAS_REF_DIR`, `HOME_REF_DIR`, `GENE_REF_FILE`, `PIXI_R_LIB`, `HPC_SCRATCH_DIR`, `SCRATCH_OUTPUT_DIR`.
Add new export: `SAMPLE_COLNAME="Sample"` (preprocess.py standardizes every dataset's sample column to `"Sample"`).

Do **not** add the flat `SCRATCH_DATA_DIR=${HPC_SCRATCH_DIR}/data` / `SCRATCH_CHUNKS_DIR=${HPC_SCRATCH_DIR}/chunks` from the TODO example: they conflict with the per-dataset convention. Dataset-specific dirs are derived in scripts (`${HPC_SCRATCH_DIR}/${DS_NAME}/data`, `${SCRATCH_OUTPUT_DIR}/${DS_NAME}/chunks`); `HOME_CHUNKS_DIR` stays the per-dataset chunk-dir mechanism exported by `2_submit_hpc_array.sh`. Document this deviation in TODO.md.

### 2. `config_helper.R` — per-dataset paths (make annotation functional)

`get_pipeline_config()` (keep interface, keep DS_NAME/HPC_SCRATCH_DIR guards):
- `path_data = file.path(scratch_output_dir, ds_name)` (preprocessed output = annotation input)
- `path_output = file.path(scratch_output_dir, ds_name)`
- `path_output_samples = file.path(scratch_output_dir, ds_name, "samples")`
- `path_output_chunks = file.path(scratch_output_dir, ds_name, "chunks")` (now matches `2_submit_hpc_array.sh:17` HOME_CHUNKS_DIR)
- `path_output_ecoda = file.path(scratch_output_dir, ds_name, "ecoda")`
- `path_ref`, `gene_ref`: unchanged.

### 3. R scripts — resolve `project_root` from env, drop `getwd()` reliance

- `src/cell_type_annotation/1.1_prepare_chunks.r:1`: `project_root <- Sys.getenv("PROJECT_ROOT"); if (project_root == "") project_root <- getwd()`; line 31: `source(file.path(project_root, "config_helper.R"))`; remove stale TODO comments (lines 1, 27). Add `if (sample_col == "") stop(...)` guard for SAMPLE_COLNAME.
- `src/cell_type_annotation/2.1.1.1_process_chunk.R:8-10`: same `project_root` + explicit `source(file.path(project_root, "config_helper.R"))`.

### 4. `src/preprocess/1.1.1_preprocess.py` — explicit path resolution

- `project_root = Path(os.environ.get("PROJECT_ROOT") or Path(__file__).resolve().parents[2])`
- Make the currently-implicit pathlib behavior explicit: for `config_path`/`base_path`/`output_dir`, only join against `project_root` when the arg is relative; keep absolute args as-is. CLI interface unchanged → worker bash untouched.

### 5. `src/preprocess/_create_combinedpbmc_dataset.py` — HPC-capable

- Add CLI args (argparse): `--base-path` (input root), `--output-dir` (combined output + rds→h5ad cache), `--layout {per-dataset,flat}`.
- Defaults: if `HPC_SCRATCH_DIR` env set → `base-path=${HPC_SCRATCH_DIR}`, `layout=per-dataset` (per-source inputs `${HPC_SCRATCH_DIR}/<ds>/data`, e.g. `Stephenson/data`, `Gongsharma_cmv_young_males/data`, `Zhu/data`), `output-dir=${HPC_SCRATCH_DIR}/CombinedPBMC/data`. Local fallback (env unset) → current behavior (`PROJECT_ROOT/data` flat, output `PROJECT_ROOT/data`).
- `load_and_prepare_source()` gets a per-source `base_path` argument.
- Must run **before** `1_submit_hpc_array.sh` (CombinedPBMC has `folder_name: null`, staging skips it; the preprocess array task reads `${HPC_SCRATCH_DIR}/CombinedPBMC/data/combined_pbmc_batch_effect_analysis.h5ad`). Run from repo root (`cd ${PROJECT_ROOT}`) on the HPC login node because `_preprocess_utils.py` R-interop sources `src/utils/load_all_functions.R` relative to cwd. Note in script docstring: heavy loads (GongSharma) may warrant running via a single `sbatch` job instead of interactively on the login node.

### 6. No changes needed (verify only)

- `src/preprocess/1_submit_hpc_array.sh`, `src/preprocess/1.1_run_worker.sh`: already per-dataset with absolute args.
- `src/cell_type_annotation/2_submit_hpc_array.sh`, `2.1_run_worker.sh`, `2.1.1_process_chunk.sh`: paths already align after change 2 (HOME_CHUNKS_DIR pattern). Optional cleanup (not required): `2.1_run_worker.sh:24` could stop passing `$2` PROJECT_ROOT since `2.1.1_process_chunk.sh` sources slurm_config.sh itself.
- `src/preprocess/_preprocess_utils.py`, `src/datasets_io.py`, `src/utils/datasets_io.R`, `src/gene_utils.py`: path-taking helpers / repo-relative resources — keep.
- Optional (mark clearly, small + functional): `2_submit_hpc_array.sh` could auto-export `TISSUE_TYPE` and `NORMAL_TISSUE` from datasets.json (`jq -r --arg ds "$DS_NAME" '.[$ds].tissue'` / `.normal_tissue`) since `2.1.1.1_process_chunk.R:51,53` reads them and nothing exports them today. `AUTHOR_ANNOT_COLNAMES` remains user-exported (unknown mapping to label columns).

### 7. Docs

- `docs/ARCHITECTURE.md`: update annotation pipeline section (path table: `path_data`/`path_output*` per-dataset under `SCRATCH_OUTPUT_DIR/<DS_NAME>`), note env vars now exported by slurm_config.sh, note combine-script HPC layout + ordering before preprocess array.
- `AGENTS.md`: one line noting the annotation path convention + SAMPLE_COLNAME export.
- `TODO.md`: move "## Replace base_path" to `# Completed` (with summary incl. decision NOT to add flat SCRATCH_DATA_DIR/SCRATCH_CHUNKS_DIR); keep benchmark decision as pending item; note combine script is HPC-capable and must run before the preprocess array.

## Out of scope

- `src/benchmark/` (NAS_SC_DIR vs SCRATCH_OUTPUT_DIR, local vs HPC) — user decision later, keep TODO note.
- `DEPRECATED_LEGACY_CODE/` — dead code, untouched.
- Full "simplify cell type annotation pipeline" refactor (TODO Step 1) — only path alignment done here.
- Running/validating on the HPC cluster (requires user login) — validation steps below are local unless marked HPC.

## Risks

| Risk | Mitigation |
|---|---|
| config_helper.R per-dataset change breaks a consumer expecting old flat dirs | Grep for `path_data`/`path_output`/`path_output_chunks` consumers (only `1.1_prepare_chunks.r`, `2.1.1.1_process_chunk.R`); alignment with `2_submit_hpc_array.sh` verified |
| srun env propagation still fails | slurm_config.sh now exports; keep the existing stop() guards for fast diagnosis; test export+`Sys.getenv` locally via `srun`-less R run |
| Combine script memory (GongSharma) on login node | Note in docstring to run via single sbatch job if OOM |
| Preprocess array picks up CombinedPBMC before file exists | Document ordering (combine → submit array); staging already skips null-folder datasets |

## Validation

1. Grep: remaining `base_path`/`project_root`/`getwd()` only in `DEPRECATED_LEGACY_CODE/`, `gene_utils.py`, `datasets_io.py` (allowed) and `.kilo/plans/`.