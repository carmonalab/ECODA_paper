# ECODA Pipelines 2–5 Structural Review, Simplification & Snakemake Architectural Plan

## 1. Status & Operational Context

**Planning & Documentation Reference Only.**
This plan is a strategic architectural roadmap and follow-up document to [`.agents/plans/20260912-final-batch-effect-subset-plan.md`](file:///Users/christianhalter/Desktop/ECODA_paper/.agents/plans/20260912-final-batch-effect-subset-plan.md). It is **not authorized for implementation during the active batch-effect execution session**.

The current production priority is completing the batch-effect analysis on Bamboo using the existing hardened Bash/durable-gate infrastructure. This document captures structural findings, architectural evaluation of workflow management (specifically Snakemake), and concrete simplification proposals according to the **Ponytail** engineering philosophy, to be executed in a dedicated refactoring phase after paper deliverables are secured.

---

## 2. Diagnosis: Pipeline Bloat & Architectural Anti-Patterns

A comprehensive structural audit across `src/` (Pipelines 2–5 and `src/utils/`) revealed that the codebase has grown to **over 73,500 lines of code**, with **18+ files exceeding 1,000 LOC** (and several exceeding 3,000 LOC).

### The Primary Anti-Pattern: "An In-House Distributed OS in Bash"
The dominant driver of this growth is that the Bash submitters and shared scripts are not merely job launchers; they implement a custom distributed runtime environment in shell:
* **Custom locking & state tracking:** Locking via directory creation, PID probing, trap handlers, and `.tombstone` markers.
* **Manifest parsing in shell:** Slicing multi-column TSVs with nested `while IFS= read` loops, regex checks, and `sed -n "${TASK_ID}p"`.
* **Resource watchdogs:** Background subshells polling `sacct` in loops to detect OOM and trigger resubmissions.
* **Path canonicalization hacks:** Complex symlink-resolution logic working around differences between host filesystems (`$HOME/scratch` vs `/srv/beegfs/...`), Slurm spool copies (`/var/spool/slurmd/...`), and Apptainer container isolation.

### Key Code Duplications & Fragilities
1. **Copy-Pasted Submitter Logic (~1,200 LOC):** `stageX_manifest_value`, `stageX_copy_atomic`, `stageX_record_identity_metadata`, and Slurm `#SBATCH` header synthesis are duplicated across Stage 2, 3, 4, and 5 submitters.
2. **Dual-Language Batch Contract (~2,840 LOC):** [`src/utils/py/batch_contract.py`](file:///Users/christianhalter/Desktop/ECODA_paper/src/utils/py/batch_contract.py) (1,589 LOC) and [`src/utils/batch_contract.R`](file:///Users/christianhalter/Desktop/ECODA_paper/src/utils/batch_contract.R) (1,253 LOC) independently implement identical batch key tokenization, sorting, hashing, and rank checking in two languages.
3. **Quadruple MD5 Checksum Implementations:** MD5 sidecar reading and verification are implemented in 4 different places (Bash `ecoda_run_common.sh`, Python `1.1.1_benchmark_methods_py.py`, R workers `1.1.1_run_benchmark_methods_r.R`, and R loader `batch_effect_analysis.R`). The R implementation strictly requires `c("MD5", "SIZE", "PATH")` line ordering, creating cross-language fragility.
4. **Candidate Dead / Orphaned Scripts:** 
   * Stage 4 stubs ([`1_prepare_chunks.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/1_prepare_chunks.sh), [`2_submit_hpc_array.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/2_submit_hpc_array.sh), [`3_submit_merge.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/3_submit_merge.sh)) exiting unconditionally with code 64.
   * Forwarding shims in Stage 5 ([`run_batch_effect_methods/1_submit_hpc_array.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/5_run_benchmark_methods/run_batch_effect_methods/1_submit_hpc_array.sh), [`run_transformation_zeroimp_analysis/1_submit_hpc_array.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh)).
   * 6 single-dataset Stage 2 wrapper scripts ([`1.1_submit_*.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/2_dataset_specific_preprocessing/) through `1.6_submit_*.sh`) each duplicating 70 lines around `1_submit_hpc.sh`.

---

## 3. The Ponytail Simplicity Audit

Applying the **Ponytail** engineering philosophy ([DietrichGebert/ponytail](https://github.com/DietrichGebert/ponytail), *"The best code is the code never written"*):

### The Simplicity Ladder for ECODA:
* **Rung 1 (YAGNI):** Hand-rolled distributed locking, PID probing, and TSV manifest slicing in Bash do not need to exist.
* **Rung 2 (Reuse):** A single canonical batch key computed once in Stage 3 Python replaces 1,253 lines of R batch contract parsing.
* **Rung 3/4 (Native / Platform):** Slurm and Linux handle job dependencies, process accounting, and file isolation natively.
* **Rung 5 (Installed/Standard Tool):** A declarative workflow engine (Snakemake) solves DAG dependencies, cluster arrays, Apptainer bindings, and resource retries.

### Ponytail Review Log (`ponytail-review`):
```text
src/4_cell_type_annotation/1_prepare_chunks.sh:L1-7: delete: 7-line stub exiting 64. Deleted, use 1_submit_onboarding_stage.sh.
src/4_cell_type_annotation/2_submit_hpc_array.sh:L1-7: delete: 7-line stub exiting 64. Nothing replaces it.
src/4_cell_type_annotation/3_submit_merge.sh:L1-7: delete: 7-line stub exiting 64. Nothing replaces it.
src/5_run_benchmark_methods/run_batch_effect_methods/1_submit_hpc_array.sh:L1-20: delete: forwarding shim. Call parent submitter with explicit args.
src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh:L1-20: delete: forwarding shim. Call parent submitter directly.
src/3_scrnaseq_preprocessing/diagnose_layer_alignment.R:L1-85: delete: orphaned one-off debug script. Nothing replaces it.
src/2_dataset_specific_preprocessing/1.1_submit_azimuth_references.sh:L1-72: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/2_dataset_specific_preprocessing/1.2_submit_covid19_pbmc.sh:L1-68: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/2_dataset_specific_preprocessing/1.3_submit_diabetes.sh:L1-68: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/2_dataset_specific_preprocessing/1.4_submit_joanito.sh:L1-68: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/2_dataset_specific_preprocessing/1.5_submit_lung.sh:L1-68: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/2_dataset_specific_preprocessing/1.6_submit_lupus_pbmc.sh:L1-68: yagni: single-dataset wrapper. Parameterize 1_submit_hpc.sh --datasets.
src/utils/batch_contract.R:L1-1253: yagni: 1,253-line duplicate batch tokenizer. Persist composite batch key in Python Stage 3 H5AD.
src/5_run_benchmark_methods/1.1.1_run_benchmark_methods_r.R:L63-88: shrink: copy-pasted ecoda_local_checksum_fields. Centralize in src/utils/helpers.R.
src/5_run_benchmark_methods/1.1.1_prepare_pseudobulk.R:L35-60: shrink: copy-pasted ecoda_local_checksum_fields. Centralize in src/utils/helpers.R.
src/5_run_benchmark_methods/consolidate_gloscope_results.R:L42-67: shrink: copy-pasted ecoda_local_checksum_fields. Centralize in src/utils/helpers.R.
```
**Estimated Direct Savings:** `-3,650 lines` from shell/R deduplication alone, reaching `>5,000 lines` with batch contract consolidation.

---

## 4. Target Architecture: The 3-Tier Model with Snakemake

Adopting Snakemake eliminates the 15,000 lines of custom Bash orchestration while preserving 100% of the scientific rules and trust-boundary validation.

```text
┌─────────────────────────────────────────────────────────────────────────┐
│  Tier 1: Durable HPC Gate (Outer Boundary)                              │
│  - Commit snapshot & runtime image validation                           │
│  - Launch envelope, session monitoring, and terminal Luna Max audit    │
│  - Synchronization ownership (NAS transfer gate)                       │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │ executes
┌────────────────────────────────────▼────────────────────────────────────┐
│  Tier 2: Declarative Workflow Engine (Snakemake)                        │
│  - Replaces 15,000 LOC of custom Bash submitters & watchdogs            │
│  - Direct DAG resolution from datasets.json                             │
│  - Native Slurm resource management & array scheduling                  │
│  - Native Apptainer container binding (ecoda-py-cuda13.sif)             │
│  - Incomplete-output automatic cleanup and atomic rule execution       │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │ dispatches
┌────────────────────────────────────▼────────────────────────────────────┐
│  Tier 3: Scientific Workers & Trust Boundaries (R & Python)             │
│  - Unchanged science: Scanpy/AnnData, GloScope, DESeq2, MOFA2, scITD    │
│  - Zero label-leakage enforcement                                       │
│  - Explicit preflights: 1.0_audit_input_views.py                        │
│  - Strict artifact validators: matrix_artifact_validator.py             │
│  - Atomic MD5 / SHA-256 sidecar production                              │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 5. Architectural Evaluation: Where Snakemake Fits vs. What It Must Not Replace

### A. What Snakemake Solves Natively
1. **DAG Dependencies:** Coordinates dataflow cleanly across stages (Raw $\to$ Stage 2 $\to$ Stage 3 Python `.h5ad` $\to$ Stage 4 Annotation $\to$ Stage 5 R `.rds` $\to$ Stage 6 Analysis).
2. **Wildcard Resolution from `datasets.json`:** Direct Python parsing of dataset lists, views, and method suites without intermediate TSV manifest scripts.
3. **Slurm Array & Resource Declarations:** Snakemake’s Slurm executor replaces dynamic `#SBATCH` generation, array bound math, and manual log-path routing.
4. **Apptainer / Singularity Execution:**
   `container: "scratch/_ecoda_runtime/ecoda-py-cuda13.sif"`
   Snakemake handles path binding, environment propagation, and worker execution natively, eliminating container bootstrap workarounds.
5. **Atomic Cleanup:** If a task crashes, Snakemake automatically deletes incomplete or truncated output files (`--delete-all-output`), eliminating corrupt artifacts.

### B. What Snakemake Does NOT Solve Out of the Box (Crucial Trust Boundaries)
As identified in the technical review, Snakemake’s default heuristics must not replace ECODA’s scientific contracts:
1. **OOM-Only Retries vs. Generic `--restart-times`:**
   * *The Risk:* Snakemake’s `--restart-times N` retries *any* failure (including code bugs, schema errors, and corrupt data).
   * *The Invariant:* ECODA requires **fail-closed** behavior on scientific errors and code defects; **only** exit code `137` (OOM / SIGKILL) may be retried with scaled memory.
   * *The Solution:* A thin execution adapter or cluster error handler that inspects the exit code before allowing Snakemake to trigger a retry.
2. **File Existence (`mtime`) vs. Cryptographic Integrity:**
   * *The Risk:* Snakemake treats a file as complete if it exists and has an updated timestamp.
   * *The Invariant:* ECODA requires non-empty arrays, schema validation, and cryptographic MD5 sidecars.
   * *The Solution:* Snakemake rules must declare the `.md5` sidecar as an explicit rule output, produced only after the worker passes internal validation.
3. **No Label Leakage:**
   * *The Invariant:* Biological labels (`Status`, `cond`, `Disease_Identity`) must never enter normalization, HVG selection, batch correction, or embeddings.
   * *The Solution:* Parameter separation remains enforced in the Python/R CLI arguments; Snakemake rules only pass unlabelled count matrices.
4. **Stage 5 Shared Sync Ownership:**
   * *The Invariant:* Stage 5 runs targeting the same analysis root must not concurrently synchronize to NAS.
   * *The Solution:* The synchronization step remains governed by the Tier 1 Durable HPC Gate.

---

## 6. Adapter Design: Integrating Existing Validators into Snakemake

Existing validation helpers can be integrated into Snakemake rules without rewrites:

### Example: Stage 3 Preprocessing Rule with Native Memory Escalation & Sidecar Validation
```python
# rules/stage3_preprocessing.smk
rule stage3_preprocess:
    input:
        raw_h5ad = "scratch/raw/{dataset}.h5ad",
        preflight_done = "scratch/audit/{dataset}_{view}.audit.done"
    output:
        proc_h5ad = "scratch/preprocessed/{dataset}_{view}_final.h5ad",
        sidecar   = "scratch/preprocessed/{dataset}_{view}_final.h5ad.md5"
    params:
        view = "{view}",
        config = "datasets.json"
    container:
        "scratch/_ecoda_runtime/ecoda-py-cuda13.sif"
    resources:
        # Dynamic memory escalation for OOM recovery
        mem_mb = lambda wildcards, attempt: 32000 * (2 ** (attempt - 1)),
        time   = "04:00:00"
    retries: 2
    shell:
        """
        # Wrapper enforces fail-closed on logic bugs, allows retry only on exit 137
        src/utils/bash/ecoda_rule_guard.sh \
            python src/3_scrnaseq_preprocessing/1.1.1_preprocess.py \
                --input "{input.raw_h5ad}" \
                --output "{output.proc_h5ad}" \
                --view "{params.view}" \
                --config "{params.config}"
        
        # Cryptographic integrity validation before rule completion
        python src/utils/py/matrix_artifact_validator.py \
            --input "{output.proc_h5ad}" \
            --write-md5 "{output.sidecar}"
        """
```

### Example: Stage 5 R Benchmark Method Rule
```python
# rules/stage5_benchmark.smk
rule stage5_benchmark_r:
    input:
        h5ad = "scratch/preprocessed/{dataset}_{view}_final.h5ad"
    output:
        rds     = "scratch/benchmarks/{view}_final/{dataset}/{dataset}_{method}.rds",
        sidecar = "scratch/benchmarks/{view}_final/{dataset}/{dataset}_{method}.rds.md5"
    params:
        method = "{method}"
    container:
        "scratch/_ecoda_runtime/ecoda-py-cuda13.sif"
    resources:
        mem_mb = lambda wildcards, attempt: 16000 * (2 ** (attempt - 1)),
        time   = "02:00:00"
    retries: 2
    shell:
        """
        src/utils/bash/ecoda_rule_guard.sh \
            Rscript src/5_run_benchmark_methods/1.1.1_run_benchmark_methods_r.R \
                --input "{input.h5ad}" \
                --method "{params.method}" \
                --output "{output.rds}"
        
        Rscript src/5_run_benchmark_methods/validate_benchmark_rds_contract.R \
            --input "{output.rds}" \
            --write-md5 "{output.sidecar}"
        """
```

---

## 7. Phased Migration Sequence (Post-Paper Implementation)

To safely execute this transformation without risking active deliverables:

### Phase 1: Call-Graph & Entry-Point Inventory
1. Catalog all live entry points in Pipelines 2–5.
2. Confirm deprecation of candidate orphan scripts ([`1_prepare_chunks.sh`](file:///Users/christianhalter/Desktop/ECODA_paper/src/4_cell_type_annotation/1_prepare_chunks.sh), forwarding shims, single-dataset wrappers).
3. Confirm that no active scheduler manifest or gate references them.

### Phase 2: Standalone Prototype on `_debug` / Synthetic Fixtures
1. Implement a minimal `Snakefile` covering Stage 3 preprocessing and Stage 5 pseudobulk/gloscope for the 5-sample `_debug` Joanito subset.
2. Run strictly on scratch outside canonical output roots (`scratch/test_snakemake/`).
3. Verify that `ecoda-py-cuda13.sif` executes inside Snakemake with identical R/Python package resolution.

### Phase 3: Contract Parity & OOM Guard Validation
1. Verify that `ecoda_rule_guard.sh` correctly aborts on syntax/schema errors and triggers retry only on exit 137.
2. Verify that output `.h5ad`, `.rds`, and `.md5` sidecars match the exact schemas of the current pipeline.
3. Test compatibility with the downstream analysis loader ([`batch_effect_analysis.R`](file:///Users/christianhalter/Desktop/ECODA_paper/src/utils/batch_effect_analysis.R)).

### Phase 4: Stage-by-Stage Production Cutover
1. Migrate Stage 3 (preprocessing and view generation).
2. Migrate Stage 4 (HiTME & scATOMIC chunk scheduling and merging).
3. Migrate Stage 5 (R and Python benchmark method array).
4. Retire the 15,000 lines of custom Bash submitters (`1_submit_hpc_array.sh`, `benchmark_submit_common.sh`, etc.).

---

## 8. Non-Goals for the Current Batch-Effect Run

During the active batch-effect execution:
* **Do NOT install Snakemake or add it to `pixi.toml`.**
* **Do NOT replace any active submitter script with Snakemake rules.**
* **Do NOT modify `datasets.json`, artifact schemas, or gate structures.**
* **Do NOT refactor the dual R/Python batch contracts.**
* Keep all active blockers focused on the execution plan: fixing the missing import in [`1.0_audit_input_views.py`](file:///Users/christianhalter/Desktop/ECODA_paper/src/3_scrnaseq_preprocessing/1.0_audit_input_views.py), resolving the Joanito RDS inspection memory ceiling, and completing the approved 4-row Stage 3 and Stage 5 runs.
