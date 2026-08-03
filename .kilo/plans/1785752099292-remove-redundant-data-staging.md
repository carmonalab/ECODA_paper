# Plan: Move raw-data staging out of cell_type_annotation into preprocess

Implements TODO.md item (lines 16-17): *"Staging of raw data should be handled in
`preprocess/1_submit_hpc_array.sh`, as this is the first HPC pipeline script to be run,
and it should not be part of `cell_type_annotation/` (remove it there)"*.

## Context (verified against code)

- `src/preprocess/1_submit_hpc_array.sh` already has a Phase 1 NAS→scratch staging step,
  but it is broken:
  - Uses `jq` (Phase 1, line 13) **before** `module load jq/1.6` (Phase 2, line 47) — fails
    on login nodes where jq is not preloaded.
  - Stages files flat into `${HPC_SCRATCH_DIR}` (line 32), while the worker
    `src/preprocess/1.1_run_worker.sh:31` reads from `DATA_DIR="${HPC_SCRATCH_DIR}/${DS_NAME}/data"`
    (`--base_path`). Files would not be found → worker exits non-zero.
- `src/cell_type_annotation/1_prepare_chunks.sh:34-48` still stages
  `NAS_SC_DIR/${DS_NAME}/data/` → `${HPC_SCRATCH_DIR}/data`, marked with an in-script
  TODO (line 37) to remove it. The NAS `data/` dirs contain no `.h5ad` files
  (e.g. Joanito: `.h5`, zips, csv), so this staging is redundant/wrong: the annotation
  pipeline's real input is the **preprocessed** `.h5ad` written by preprocess.py to
  `${SCRATCH_OUTPUT_DIR}/${DS_NAME}/`.
- `config_helper.R` `path_data = ${HPC_SCRATCH_DIR}/data` and the
  `2_submit_hpc_array.sh` vs `config_helper.R` chunk-dir seam are **deliberately out of
  scope** (user decision): the annotation pipeline will remain non-functional until the
  later "simplify cell type annotation pipeline" step fixes its input paths.
- Doc convention: `AGENTS.md`, `ARCHITECTURE.md`, `README.md` kept up to date.

## Changes

### 1. Fix raw-data staging in `src/preprocess/1_submit_hpc_array.sh`

- Move `module load GCCcore/12.2.0` + `module load jq/1.6` to the top of the script,
  before Phase 1.
- Extend the Phase 1 `jq` to emit the dataset **key** alongside folder and file, and
  stage each file into a per-dataset dir matching the worker:
  ```bash
  jq -r '
    to_entries[] |
    .key as $key |
    .value.folder_name as $folder |
    .value.views | to_entries[] |
    .value.input_file_name |
    if type == "array" then .[] else . end |
    "\($key) \($folder) \(.)"
  ' "${DATASETS_JSON_FILE}" | sort -u | \
  while read -r KEY FOLDER_NAME RAW_FILE_NAME; do
      if [ "$FOLDER_NAME" == "null" ] || [ -z "$FOLDER_NAME" ]; then
          continue
      fi
      NAS_FILE_PATH="${NAS_SC_DIR}/${FOLDER_NAME}/output/${RAW_FILE_NAME}"
      DEST_DIR="${HPC_SCRATCH_DIR}/${KEY}/data"
      echo "Dataset folder: ${FOLDER_NAME}, file: ${RAW_FILE_NAME}"
      if [ -f "$NAS_FILE_PATH" ]; then
          mkdir -p "${DEST_DIR}"
          rsync -ah --progress "$NAS_FILE_PATH" "${DEST_DIR}/"
      else
          echo "  -> [WARNING] Source not found: ${NAS_FILE_PATH}"
      fi
  done
  ```
  The `mkdir -p "${HPC_SCRATCH_DIR}"` before the loop is superseded by the per-dataset
  `mkdir -p "${DEST_DIR}"`; keep or drop (either is fine — drop it to avoid dead code).
- No other changes to the submit script (Phase 2/3 unchanged).

### 2. Remove raw-data staging from `src/cell_type_annotation/1_prepare_chunks.sh`

- Delete lines 34-48 (the `HOME_DATA_DIR`/`NAS_DATA_DIR` variables and the whole
  "STAGE DATA" rsync block, including its "Skipping rsync file transfer" branch).
- Delete the stale TODO comment on line 37.
- Keep the reference-maps staging (`NAS_REF_DIR` → `HOME_REF_DIR`) and the gene-annotation
  download (`GENE_REF_FILE`) blocks unchanged.
- Note: `DS_NAME` is still required by this script (used by chunk generation) — do not
  remove the check.

### 3. Documentation updates

- `docs/ARCHITECTURE.md`:
  - Table row for `1_prepare_chunks.sh` (line 59): change "stages data + ref maps from
    NAS → scratch" to "stages ref maps from NAS → scratch (raw data is staged by
    `src/preprocess/1_submit_hpc_array.sh`)".
  - "NAS ↔ Scratch data flow" note (line 72): clarify that raw-data staging happens only
    in the preprocess submit script; cell type annotation consumes preprocessed output.
- `TODO.md`: mark the item done (move to the `# Completed` section, e.g.
  `- [x] Raw-data staging consolidated in src/preprocess/1_submit_hpc_array.sh; redundant staging removed from src/cell_type_annotation/1_prepare_chunks.sh`).
- `AGENTS.md`: no change needed (already describes `1_prepare_chunks.sh` as reading
  preprocessed `.h5ad` files).

## Out of scope (explicitly deferred)

- Repointing `config_helper.R` `path_data` to `${SCRATCH_OUTPUT_DIR}/${DS_NAME}` and
  aligning annotation output paths per-dataset.
- Exporting `SAMPLE_COLNAME="Sample"` (preprocess.py standardizes the sample column to
  `"Sample"`; the annotation scripts still default to `"sample"`).
- These will be handled by the later "simplify cell type annotation pipeline" step;
  until then the annotation pipeline is known non-functional.
- CombinedPBMC (`folder_name: null`) is not staged — its input is created locally by
  `_create_combinedpbmc_dataset.py` (TODO.md line 18). Pre-existing gap; the preprocess
  array task for CombinedPBMC will fail until the file is placed on scratch/NAS manually.
  Not addressed here.

## Risks

| Risk | Mitigation |
|---|---|
| Staging jq change introduces a syntax/path regression | Test the jq expression and loop locally against `datasets.json` (below) |
| Worker fails if a staged file is missing | Staging warns on missing NAS files; `1.1_run_worker.sh` will still error for genuinely absent inputs (CombinedPBMC, see above) |
| Docs drift | Section 3 updated in the same change |

## Validation

- Run the new jq expression locally and check output is `KEY FOLDER FILE` triples,
  deduplicated, e.g.:
  ```bash
  jq -r 'to_entries[] | .key as $key | .value.folder_name as $folder |
         .value.views | to_entries[] | .value.input_file_name |
         if type == "array" then .[] else . end |
         "\($key) \($folder) \(.)"' datasets.json | sort -u
  ```
- Spot-check that `NAS_SC_DIR/<folder>/output/<file>` exists for all emitted triples
  (NAS currently mounted at `/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets`),
  and that `${HPC_SCRATCH_DIR}/${KEY}/data/<file>` matches what
  `src/preprocess/1.1_run_worker.sh` passes as `--base_path` + `input_file_name`.
- Grep `cell_type_annotation/` for any remaining references to `HOME_DATA_DIR`,
  `NAS_DATA_DIR`, or the removed staging block (should be none outside
  `DEPRECATED_LEGACY_CODE/`).