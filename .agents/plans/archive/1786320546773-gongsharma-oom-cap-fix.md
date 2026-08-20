# GongSharma OOM fix: per-sample 5000-cell cap step (TODO.md:146-157)

## Goal

Fix the preprocess-array OOM for `Gongsharma_cmv_young_males` (task 4, 128G worker) by
implementing a new auto-discovered dataset-specific preprocessing step that caps each
sample at 5000 cells (historical `downsample_by_group` strategy, seed 42), then re-run
the pipeline for that dataset.

**Decision (user-confirmed): option (a) — in-place overwrite of the two staged
SoundLife h5ads in `${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data`. No
datasets.json change (ground truth untouched; NAS remains the pristine source).**

## Context (verified facts)

- Staged inputs: `SoundLife_YoungAdult_Male_CMVneg.h5ad` (1,712,244 cells, 107 samples,
  max 33,900/sample) and `SoundLife_YoungAdult_Male_CMVpos.h5ad` (1,206,761 cells,
  73 samples, max 30,849/sample); 33,538 genes; X = uint16 CSR counts, no layers, no raw.
- OOM cause: the preprocess worker holds both files + `sc.concat` (2.92M cells) and
  densifies `sc.pp.scale` matrices per HVG pass (2.92M × 3000 × float64 ≈ 70 GB) —
  exceeds 128G.
- Cap outcome (seed 42, 5000/sample): 531,291 + 365,000 = **896,291 cells, 180 samples,
  max 5000/sample** — exactly matches the historical NAS artifact
  `Gongsharma_2024_PrePrintTBD/output/Gongsharma_cmv_young_males.h5ad` (896,291 cells /
  180 samples / max 5000). That file is the validation target for the sampling.
- **Race hazard**: `1.1.1_create_combinedpbmc_dataset.py` reads the same staged
  SoundLife files (backed mode) and `1_submit_hpc.sh` submits all steps in parallel.
  An in-place overwrite racing the CombinedPBMC read would nondeterminize the
  CombinedPBMC dataset. The dispatcher must serialize cap → CombinedPBMC.
- Historical reference: `downsample_by_group` at git 3a4711e
  (`src/py/preprocess_gongsharma.qmd`, identical copy kept as draft at
  `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd`).

## Implementation tasks

### 1. New `src/2_dataset_specific_preprocessing/1.4.1_subset_gongsharma.py`

Standalone script (scanpy/anndata/numpy only, no R interop). Run with `${PYTHON_BIN}`.

- CLI: `--config_path` (default `${DATASETS_JSON_FILE}` env) and `--data_dir`
  (default `${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data`).
- Read the `Gongsharma_cmv_young_males` entry via
  `src.utils.py.datasets_io.read_datasets_json`; take `file_names` (must be exactly
  the two SoundLife names — fail loudly on drift) and the sample column
  `columns.sample` (`specimen.specimenGuid`).
- Per file (process sequentially, one at a time — peak memory ≈ one capped matrix):
  1. `adata = sc.read_h5ad(path, backed='r')` (obs/var eager; X lazy).
  2. Build keep-mask exactly mirroring historical `downsample_by_group`:
     - `groups = adata.obs[group_col].unique()` (pandas order of first appearance —
       same call as historical, so deterministic per file).
     - single `rng = np.random.RandomState(42)` created before the loop; per group:
       `pos = np.where(adata.obs[group_col] == group)[0]`;
       if `len(pos) > 5000`: `keep[rng.choice(pos, 5000, replace=False)] = True`;
       else `keep[pos] = True`.
     - Use a boolean mask (not the historical index list) so the output preserves the
       **original row order**; the selected cell SET is identical to historical
       (same seed, same group order, same `RandomState.choice` sequence).
  3. `capped = adata[mask].to_memory()` (anndata backed slice materializes only the
     needed HDF5 chunks — same proven pattern as the CombinedPBMC GongSharma worker).
  4. `close()` the backed handle, then write to a temp file in the same dir
     (`<name>.capped_tmp.h5ad`) and `os.replace()` onto the staged path (atomic —
     never leave a corrupted staged file on failure).
- Print per-file summary (cells before/after, samples, max cells/sample) and the total
  — expected 531,291 / 365,000 / 896,291 (compare against NAS historical artifact).
- Idempotent by construction: re-running on an already-capped file keeps every cell
  (all groups ≤ 5000) → byte-identical re-cap. No skip logic needed.
- Preserve all obs/var/obsm/uns/varp as-is (anndata subsetting handles row alignment;
  obsm X_umap / uns colors are harmless to keep).

### 2. New `src/2_dataset_specific_preprocessing/1.4_submit_gongsharma.sh`

- `#SBATCH`: `--job-name=gongsharma_cap`, `--time=02:00:00`, `--nodes=1`,
  `--ntasks=1`, `--cpus-per-task=4`, `--mem=128G`, `--mail-type=END,FAIL`.
- Mandatory boilerplate (AGENTS.md): SCRIPT_DIR from `scontrol show job`
  (`Command=` field) with `BASH_SOURCE` fallback, `source "${SCRIPT_DIR}/../slurm_config.sh"`,
  `cd "${PROJECT_ROOT}"`.
- Body: `${PYTHON_BIN} "${SCRIPT_DIR}/1.4.1_subset_gongsharma.py"`.
- Executable bit `100755` (repo convention). Header comment: purpose, in-place
  overwrite semantics, re-staging caveat (re-stage restores uncapped files; re-run
  `1_submit_hpc.sh` before the preprocess array), and the CombinedPBMC ordering note.

### 3. Modify `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` (serialization)

The cap step is auto-discovered by the existing `1.*_submit_*.sh` glob — no discovery
change needed. But steps run in parallel and the cap must finish before CombinedPBMC:

- Before the submission loop: if `1.4_submit_gongsharma.sh` exists in `SCRIPT_DIR`,
  submit it first (same sbatch args as the loop), capture its id as `CAP_JOB_ID`,
  and prepend it to `JOB_IDS` (keeps wait/sacct gate/summary semantics unchanged).
- In the loop, when the step basename is `1.1_submit_combinedpbmc.sh` and
  `CAP_JOB_ID` is non-empty, add `--dependency=afterok:${CAP_JOB_ID}`.
  `afterok` = fail-closed: if the cap fails, CombinedPBMC is not submitted, the gate
  reports non-COMPLETED, dispatcher exits 1.
- The `-f` existence guard keeps the dispatcher working even if the step is absent.
- Update the header comment (parallelism note: 1.1 depends on 1.4). `--sync-only`
  mode is unaffected.

### 4. Documentation

- `docs/ARCHITECTURE.md`: add the `1.4` step row (role: in-place per-sample 5000-cell
  cap; inputs/outputs; expected counts; re-staging caveat) and note the
  cap → CombinedPBMC serialization in the `1_submit_hpc.sh` description and the
  CombinedPBMC step blurb (GongSharma component now reads capped files; rebuilt
  automatically on every `1_submit_hpc.sh` run — content changes to ≤5000 cells/sample,
  deterministic seed 42).
- `TODO.md`: mark the GongSharma OOM task done (lines 146-157), noting option (a)
  chosen and the re-run commands.

## Validation (user, HPC — agents do not run pipeline scripts per AGENTS.md)

1. Check no stale output blocks the preprocess skip:
   `ls "${HOME}/scratch/ECODA_paper/Gongsharma_cmv_young_males/output/"` — the failed
   array (gate failed closed) should have left nothing; use `--force` only if a
   partial output exists.
2. `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — verify the
   `gongsharma_cap` job COMPLETED, its log shows 531,291 / 365,000 cells and max
   5000/sample, and CombinedPBMC ran after the cap (job summary + squeue ordering).
3. (Optional, strong check) per-sample barcode sets of the capped staged files equal
   those of the NAS historical `Gongsharma_cmv_young_males.h5ad` (same seed/strategy).
4. `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name Gongsharma_cmv_young_males`
   — must reach the NAS sync gate (sacct all COMPLETED; sync email; verify
   `Gongsharma_cmv_young_males_benchmark_analysis_ECODAprocessed.h5ad` on NAS).

## Risks / fallbacks

- **Capped preprocess still tight at 128G** (~896k cells; estimated peak 90-110 GB:
  X + counts-layer copies ≈ 50 GB, scale-dense ≈ 21 GB/pass, harmony/leiden/neighbors).
  If OOM recurs: bump `--mem` for this dataset only in `1.1_run_worker.sh` (per-DS_NAME
  conditional) — do NOT change the cap without approval (alters dataset content).
- **Re-staging caveat**: `1_stage_data.sh` rsync overwrites capped files with NAS
  originals; the cap step re-caps on the next `1_submit_hpc.sh` (documented in the
  step header). Pipeline convention already runs staging → dataset-specific steps →
  preprocess array in order.
- **CombinedPBMC content change**: its GongSharma component becomes capped (≤5000
  cells/sample, seed-42 random subset) — deterministic, no composition bias introduced
  beyond random per-sample subsampling; batch-effect view rebuilt automatically.
- **pandas `.unique()` order** on the categorical obs column: deterministic per file
  and identical to the historical call; the optional NAS-barcode comparison catches
  any drift.

## Out of scope

- Option (b) (single h5ad + datasets.json change) — rejected by user.
- Changing the cap value / seed / group key.
- Annotation and benchmark stages for GongSharma (follow after preprocess succeeds).
- Other Phase 6 recovery items (Joanito/Smillie/Wu/Zhang/Stephenson re-runs, etc.).
