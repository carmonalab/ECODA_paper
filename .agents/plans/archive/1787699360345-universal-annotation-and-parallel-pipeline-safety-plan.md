# Universal Dataset Annotation & Parallel Multi-Dataset Pipeline Execution with Safety Checks

## Problem & Context
Previously, `src/4_cell_type_annotation/1_submit_onboarding_stage.sh` had a hardcoded assumption that only `Lupus_PBMC` required automated cell-type annotation, under the mistaken notion that author-annotated cohorts could bypass automated annotation.

In reality:
1. **Mandatory Automated Annotation for Benchmark Views & Suitable Cohorts:** Every single-cell benchmark dataset and all suitable batch-effect cohorts in the ECODA study must undergo automated cell-type annotation via **HiTME** (layers 1–3) and **scATOMIC** (layers 1–6, predicted labels, confidence, cell cycle scores). Original author annotations are preserved as baseline ground truth metadata in `obs`, while HiTME and scATOMIC provide standardized, uniform cross-cohort cell-type annotations.
2. **Explicit Exemption via `not_suitable_for_auto_annotation`:** Datasets explicitly marked in [datasets.json](file:///Users/christianhalter/Desktop/ECODA_paper/datasets.json) with `"not_suitable_for_auto_annotation": ["hitme", "scatomic"]` (e.g., `Alzheimer`, `Diabetes`, `Parkinson` due to tissue-incompatible reference maps or non-applicable models) must **not** run through the automated annotation pipeline and must be cleanly skipped with explicit logging.
3. **Multi-Dataset Parallel Execution with Fail-Closed Safety:** All pipeline stages (preprocessing, cell-type annotation, and benchmark method evaluations) must execute across **all eligible datasets in parallel SLURM arrays**, protected by rigorous **safety and idempotency checks** that verify artifact integrity before skipping already-completed runs, with full `--force` recomputation support.
4. **Test Suite Cleanup:** Stale test scripts created as temporary one-off checks during past bug fixes (such as `tests/test_bassez_and_benchmark_regressions.R` and fragile commit-specific submitter stubs) represent unnecessary bloat and will be removed.

---

## Proposed Changes

### 1. Repository Invariants & Documentation

#### [MODIFY] [AGENTS.md](file:///Users/christianhalter/Desktop/ECODA_paper/AGENTS.md)
- **Cell Type Annotation Invariant:**
  - Codify in *Non-negotiable scientific and repository rules* and *Architecture & Data Flow* that **all benchmark datasets and suitable cohorts** must be annotated with both HiTME (layers 1–3) and scATOMIC (layers 1–6, predicted labels, confidence, cell cycle phases).
  - Explicitly document that datasets declaring `"not_suitable_for_auto_annotation"` in [datasets.json](file:///Users/christianhalter/Desktop/ECODA_paper/datasets.json) (e.g., `Alzheimer`, `Diabetes`, `Parkinson`) are exempt and must be skipped by the automated annotation pipeline.
  - Clarify that author annotations serve as ground-truth metadata labels in `obs`, while HiTME and scATOMIC provide standardized, uniform cross-cohort annotations.
- **Multi-Dataset Parallel Execution Standard:**
  - Codify that all pipeline stages must dispatch all target datasets concurrently in parallel SLURM arrays rather than running single datasets sequentially.
- **Fail-Closed Idempotency & Safety Invariant:**
  - Define exact skip-check criteria: check file non-emptiness, valid schema/structure, and checksum integrity.
  - Require atomic writes (`.tmp.<pid>` -> atomic rename) so partial runs never leave misleading artifacts.
  - Require `--force` support across all submission and worker entry points.
- **Testing & QA Section Update:**
  - Remove references to deleted stale test scripts (`tests/test_bassez_and_benchmark_regressions.R`, etc.) from development commands and documentation.

#### [MODIFY] [docs/ARCHITECTURE.md](file:///Users/christianhalter/Desktop/ECODA_paper/docs/ARCHITECTURE.md)
- Update Stage 4 section to state that cell-type annotation is applied across all benchmark and suitable batch-effect cohorts, respecting `not_suitable_for_auto_annotation` flags in `datasets.json`.

---

### 2. Stage 4: Cell-Type Annotation Pipeline (`src/4_cell_type_annotation/`)

#### [MODIFY] [src/4_cell_type_annotation/1_submit_onboarding_stage.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/1_submit_onboarding_stage.sh)
- **Multi-Dataset Default with Suitability Filter:** Update default dataset resolution to include all active datasets defined in [datasets.json](file:///Users/christianhalter/Desktop/ECODA_paper/datasets.json) (skipping `_debug` unless explicitly requested), automatically filtering out datasets flagged with `"not_suitable_for_auto_annotation"` for `hitme` and `scatomic`.
- **Remove Outdated Comments:** Clean up comments that claimed author-annotated cohorts do not need annotation.
- **Unified Multi-Dataset Array Manifest:** Collect chunk files across all eligible target datasets into a single combined manifest, launching them in a single parallel SLURM array with the compute-node watchdog (`1.2_annotation_watchdog.sh`).
- **Idempotency & Clean-Entry Gate:** Skip chunk creation and array execution for datasets whose chunk annotations already exist and pass clean-entry checks (`1.1_prepare_chunks.py --check-clean`), unless `--force` is specified.
- **Multi-Dataset Atomic Merge:** Execute [3_submit_merge.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/3_submit_merge.sh) across all datasets whose annotation chunks are validated by the watchdog.

#### [MODIFY] [src/4_cell_type_annotation/1_prepare_chunks.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/1_prepare_chunks.sh)
- Respect `not_suitable_for_auto_annotation` from [datasets.json](file:///Users/christianhalter/Desktop/ECODA_paper/datasets.json) during dataset iteration: skip chunk creation for flagged unsuitable cohorts (logging a clean skip notice) and prepare chunks for all suitable datasets.

---

### 3. Stage 3: Preprocessing Pipeline (`src/3_scrnaseq_preprocessing/`)

#### [MODIFY] [src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh) & [1_submit_batch_effect_stage.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/3_scrnaseq_preprocessing/1_submit_batch_effect_stage.sh)
- Ensure all cohorts and declared views are submitted in parallel arrays.
- Validate safety checks in [1.1.1_preprocess.py](file:///Users/christianhalter/Desktop/ECODA_paper/src/3_scrnaseq_preprocessing/1.1.1_preprocess.py):
  - Check that existing output `.h5ad` views are non-empty and valid before skipping.
  - Recompute if missing, incomplete, or if `--force` is provided.

---

### 4. Stage 5: Benchmark & Embedding Methods (`src/5_run_benchmark_methods/`)

#### [MODIFY] [src/5_run_benchmark_methods/benchmark_submit_common.sh](file:///Users/christianhalter/Desktop/ECODA_paper/src/5_run_benchmark_methods/benchmark_submit_common.sh)
- Maintain parallel multi-dataset SLURM array dispatch.
- Ensure result-bundle safety checks:
  - Check existence and MD5 checksum integrity of `.rds` bundles and distance matrices in `results/`.
  - Skip completed method evaluations on a per-dataset, per-method, per-parameter level unless `--force` is specified.
  - Retain compute-node OOM watchdog with bounded memory auto-escalation.

---

### 5. Cleanup of Stale Test Scripts (`tests/`)

#### [DELETE] [tests/test_bassez_and_benchmark_regressions.R](file:///Users/christianhalter/Desktop/ECODA_paper/tests/test_bassez_and_benchmark_regressions.R)
- Remove stale one-off regression script originally written for a historical Bassez metadata fix.

#### [DELETE] [tests/test_preprocessing_stage_submitter.sh](file:///Users/christianhalter/Desktop/ECODA_paper/tests/test_preprocessing_stage_submitter.sh)
- Remove fragile mock submitter test hardcoded to historical onboarding row order.

---

## Verification Plan

### Automated Checks
1. **Repository Health & Environment Dependencies:**
   ```bash
   pixi run check-r-deps
   bash src/5_run_benchmark_methods/test_oom_retry.sh
   ```
2. **Shell Syntax Verification:**
   ```bash
   bash -n src/4_cell_type_annotation/1_submit_onboarding_stage.sh
   bash -n src/4_cell_type_annotation/1_prepare_chunks.sh
   bash -n src/4_cell_type_annotation/1.2_annotation_watchdog.sh
   bash -n src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh
   bash -n src/5_run_benchmark_methods/benchmark_submit_common.sh
   ```
3. **Dataset Suitability & Skip Logic Verification:**
   - Verify that `not_suitable_for_auto_annotation` datasets (`Alzheimer`, `Diabetes`, `Parkinson`) are properly skipped by `1_prepare_chunks.sh` and `1_submit_onboarding_stage.sh`.
   - Verify on `_debug` dataset that chunking, clean-entry checking, and parallel array manifest generation work seamlessly.
