# Repository Guidelines

## Project Overview

ECODA (Exploratory Compositional Data Analysis) is a reproducible R/Python
workflow for unsupervised patient stratification from single-cell cohorts. It
compares CLR-based cell-type composition, pseudobulk, and sample-embedding
methods, then scores recovery of known biological groups.

### Non-negotiable scientific and repository rules

- **No label leakage.** Biological labels such as `Status`, `sample.origin`,
  `cond`, and `Disease_Identity` are ground truth only. Never pass them to
  preprocessing, HVG selection, normalization, batch correction, embeddings,
  or model covariates.
- `DESeq2.normalize()` benchmark defaults are `blind=TRUE`, `batch_col=NULL`,
  `correct_batch=FALSE` (`~ 1`). Batch-effect mode is batch-only:
  `blind=FALSE`, `batch_col=<batch>`, `correct_batch=TRUE`; never protect
  biological labels in `removeBatchEffect`.
- `datasets.json` is the dataset/view ground truth. **Do not modify it without
  explicit user confirmation.**
- **Universal cell-type annotation with HiTME and scATOMIC.** All
  benchmark-analysis datasets and all suitable cohorts used by benchmark
  workflows must undergo dual automated annotation with HiTME (layers 1–3)
  and scATOMIC (layers 1–6, predicted labels, confidence, and cell-cycle
  scores). Author annotations remain baseline ground-truth metadata in `obs`.
  Datasets marked `not_suitable_for_auto_annotation` in `datasets.json` are
  exempt and must be cleanly skipped by the annotation pipeline.
- **Batch-effect-only annotation exception.** The
  `batch_effect_uncorrected` and `batch_effect_corrected` views do not require
  HiTME/scATOMIC Pipeline 4. Do not schedule Pipeline 4 solely for
  batch-effect analysis; preserve and use configured source/author cell-type
  columns. This view-scoped exception does not weaken the universal
  annotation rule for benchmark-analysis data or configured exemptions.
- **Fail-closed idempotency.** Every stage verifies existing outputs for
  integrity, non-emptiness, schema, ownership, and checksums before reuse.
  `--force` is reserved for explicitly scoped recomputation with a recorded
  dependency or integrity reason; routine checks must not trigger broad
  reruns.
- Files beginning with `Figure` or `Supp_fig` are publication figures: fix
  them, never remove them. `Figure 2A` uses default/main settings, `Supp fig
  15` contains extended methods, and `Supp fig 2` is parameter screening.
  Exclude legacy `ECODA_PB_combo_*` from publication figures.
- Preserve all version constraints in `pixi.toml` and the resolved `pixi.lock`.
- Use focused tests and the `_debug` Joanito five-sample subset for routine
  verification. Do not launch full cohorts for routine checks.

## Architecture & Data Flow

1. **Configuration:** `datasets.json` defines datasets, metadata columns,
   views, and filenames. `src/utils/datasets_io.R` and
   `src/utils/py/datasets_io.py` are the language-specific access layer.
2. **Data staging:** `src/1_stage_data/1_stage_data.sh` stages raw data from
   NAS to scratch; `src/2_dataset_specific_preprocessing/` performs
   cohort-specific conversion and harmonization.
3. **Canonical preprocessing:** `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`
   filters data, preserves raw counts in `layers["counts"]`, normalizes
   `X`, ranks HVGs, computes PCA, and creates Harmony/neighbors/Leiden
   outputs.
4. **Cell-type annotation:** `src/4_cell_type_annotation/` prepares chunks,
   runs the configured annotation workers, checkpoints per-sample Feather
   output, validates completeness, and merges annotations into eligible
   benchmark-analysis `.h5ad` files. Batch-effect views use the exception
   above and retain configured source/author columns.
5. **Benchmarking:** `src/5_run_benchmark_methods/` runs R and Python methods
   through SLURM arrays. Methods converge on sample feature matrices,
   distance matrices, or `create_result_bundle(feat_mat, labels, dist_mat)`
   bundles.
6. **Scoring and persistence:** `src/utils/scoring_metrics.R` computes
   silhouette, modularity, ANOSIM, ARI, and LISI. Results are written
   atomically as `.rds` bundles with checksum sidecars; Feather carries
   cross-language embeddings, distances, and execution logs.
7. **Analysis:** local notebooks consume precomputed results and generate
   publication figures; they do not rerun cohort preprocessing.

Operational concurrency is explicit rather than application-async: R uses
`foreach`/`doParallel`, Python/R workers run in SLURM arrays, and shell
watchdogs gate synchronization on scheduler accounting. Missing status,
checksum mismatch, worker failure, or exhausted OOM retry must fail closed.

### Pipeline ordering and durable HPC execution

- Numbered pipelines remain ordered: Pipeline 1 is the NAS-bound serial
  staging step; Pipeline 2 performs dataset-specific preprocessing, Pipeline
  3 canonical preprocessing, Pipeline 4 annotation, and Pipeline 5
  benchmarking. Independent selected rows within a stage MAY and SHOULD run
  concurrently in one SLURM array; only documented dependencies may serialize
  work.
- Every full-cohort preprocessing, annotation, benchmark, evidence, or
  correction run MUST use the checked-in `durable-hpc-gate-ecoda` profile.
  Direct SSH-launched long-running wrappers are not an acceptable substitute.
- **Canonical host policy:** `bamboo` is the default and canonical host for
  normal ECODA Pipeline 1–5 durable compute and gates. `yggdrasil`, reached as
  `ssh yggdrasil`, is never an implicit compute host. Only for the
  2026-09-15–18 Bamboo maintenance window may it be used as the explicitly
  named temporary backup destination for repository/scratch clone and restore
  checks. This does not change the durable-gate `remote_host=bamboo` policy,
  source/runtime contracts, or NAS synchronization. Bamboo↔Yggdrasil transfer
  requires an explicitly configured SSH key or agent forwarding; passwords
  must never be stored or used. Every backup command states its host,
  source, destination, and scope.
- Arrays that can OOM MUST use a compute-node watchdog with OOM-only
  resubmission of affected manifest rows, bounded memory escalation, and
  fail-closed handling of non-OOM failures or an exhausted ceiling.

#### Stage 5 synchronization and gate contract

- Independent Stage 5 datasets and methods sharing an `ANALYSIS_ROOT` MUST be
  submitted together in one explicit selection-manifest/SLURM wave. A single
  global synchronization owner protects `sync/${ANALYSIS_ROOT}`; separate
  gates using the same root are serialized through the durable profile's
  `ecoda-benchmark` group.
- Do not invent a second serialization group or bypass the shared owner.
  Queue later same-root gates only after the current gate's terminal wait,
  inspect, and reviewer approval. Atomic ownership and checksum checks remain
  mandatory.
- Before each approved launch, record the exact wrapper, selected
  dataset/view/method rows, and expected row count. If emitted rows exceed
  that scope, cancel every scheduler ID and the durable runner, preserve
  evidence, and terminal-inspect the run as `FAILED`.
- Repair, validation, checksum, and reviewer/release audits are no-compute
  operations by default. Inspect existing artifacts, manifests, statuses, and
  checksums first. Valid rows remain outside recomputation selections.
- Targeted recovery is mandatory: rerun only failed, missing, or invalid
  dataset/view/method/parameter rows, and record the dependency or integrity
  reason. Never use broad `--force` selection for repair or recovery.
- After `launch`, arm exactly one unbounded durable `wait`; do not repeatedly
  poll `squeue` or `sacct` from the agent session. Perform one terminal
  `inspect` with every scheduler/watchdog ID emitted by the wrapper, followed
  by the required Luna Max reviewer approval. Use `status` only for
  non-mutating recovery checks, and never rerun a wrapper after an ambiguous
  launch.

### Current plans and immutable run identity
- **Baseline provenance:** The processing plan was prepared from repository
  anchor `5302671ad94556edcf9acccf372d2dc34121d714`; record later
  implementation and documentation commits separately.

- Current dataset/view/method selections, output roots and stems, row counts,
  gate IDs, incidents, maintenance details, and one-off recovery history are
  task-specific. Keep them in the authoritative active plan, e.g.
  (this is just an example, as suggested by the user. use specific plans as needed and discussed with the user)
  `.agents/plans/20260912-final-batch-effect-subset-plan.md`, not in this
  durable guide. Use the checked-in
  `.agents/skills/durable-hpc-gate-ecoda/SKILL.md` and its durable-hpc-gate
  references for operational command details.
- Full-cohort runs require a commit-keyed source snapshot and a versioned
  runtime identity before durable-gate `prepare`. Create the snapshot from a
  clean full-commit checkout and execute from the snapshot, never from the
  mutable canonical checkout or an unversioned runtime.
- The run-owned manifest MUST bind the immutable source manifest, runtime
  image and manifest, auxiliary-root identity, exact run ID, scratch/log
  roots, wrapper, and explicit selected-row scope. Resolve Bamboo's home
  before composing remote paths (for example with
  `ssh bamboo 'printf %s "$HOME"'`).
- The run-scoped `ecoda_run_audit.sh` checks the selected stage, manifest,
  source identity, runtime identity, and artifacts before terminal review; it
  must not scan all run roots or submit repair compute.
- Unpinned direct stage submitters are legacy/validator-only. They may submit
  new work only when invoked by the immutable snapshot executor with the
  required source/runtime/auxiliary identities, exact run ID, and explicit
  row scope. Validator checks may inspect existing artifacts but must not
  create scheduler work.
- A user-explicit temporary noncanonical scoring run may bypass the durable
  profile only when its scope, inputs, and output location are clear; it must
  use named datasets/columns, write outside canonical outputs, and preserve
  validated artifacts. This exception does not authorize unpinned submitter
  commands.

### Local resource boundary

- Full-cohort H5ADs and whole expression/count matrices MUST NOT be copied,
  staged, or retained on the local workstation. Operate against authoritative
  HPC scratch or NAS and read only the metadata needed by the local result.
- Aggregate derived results dataset-by-dataset and release each dataset's
  metadata before processing the next. Local subset mirrors are diagnostic
  only and never replace authoritative full-cohort sources.
- The only permitted local transfer is an explicitly scoped small debug or
  metadata artifact. Existing H5ADs, RDS bundles, pseudobulks, manifests,
  checksums, and gates are immutable; new derived outputs use separate
  run-owned roots.

## Key Directories

- `src/1_stage_data/` — NAS-to-scratch staging.
- `src/2_dataset_specific_preprocessing/` — cohort-specific conversion.
- `src/3_scrnaseq_preprocessing/` — shared Scanpy preprocessing and H5AD
  production.
- `src/4_cell_type_annotation/` — annotation preparation, workers, and merge.
- `src/5_run_benchmark_methods/` — benchmark workers, submitters, watchdogs,
  and result synchronization.
- `src/utils/` — dataset I/O, scoring, imports, environment checks,
  preprocessing utilities, and shell environment setup.
- `notebooks/` — local analysis, publication figures, and onboarding reports.
- `tests/` — focused standalone regressions.
- `data/` — large/gitignored data; never scan recursively or delete
  recursively without explicit confirmation.
- On `bamboo`, `$HOME/scratch/ECODA_paper` is data storage, not a git clone;
  the HPC repository is `$HOME/ECODA_paper`.

## Runtime and tooling

- Pixi is the package/environment manager. Do not introduce Conda, renv,
  pip-only, npm, or a second lockfile.
- Do not assume a system `python` or `python3`. Prefer
  `pixi run -e default python ...` for local Python work. OMP Python eval
  must use this repository's absolute Pixi interpreter.
- Remote binary artifacts (`.feather`, `.parquet`, `.arrow`, `.h5`,
  `.h5ad`, `.rds`, `.sqlite`, and similar) must be inspected with the
  configured Pixi interpreter over SSH or copied with an explicitly scoped
  `scp`/`rsync`; text-only remote reads are insufficient.
- Required runtimes are R `4.5.2` and Python `3.13.*`; HPC workers use the
  `py-cuda13` Pixi environment.
- Worker jobs must source `src/slurm_config.sh`, enter `${PROJECT_ROOT}`, and
  invoke the immutable configured commands:

  ```bash
  source src/slurm_config.sh
  cd "${PROJECT_ROOT}"
  "${PYTHON_BIN}" path/to/worker.py
  ${PIXI_RSCRIPT} path/to/worker.R
  ```

  Never use bare `python`, `Rscript`, or ordinary `pixi run` inside jobs;
  `PYTHON_BIN` and `PIXI_RSCRIPT` resolve to the pinned worker environment.
- Environment setup and refresh serialize on `logs/env_refresh.lock` and must
  not run while arrays are active. Login nodes are for editing, compilation,
  staging, NAS sync, and SLURM submission only.
- Shared benchmark shell helpers and tests support Bash 3.2; avoid newer Bash
  syntax unless a target script explicitly requires it. Submitted jobs may
  start in `/var/spool/slurmd/`; recover the source directory before sourcing
  configuration.
- Never run `rm -rf` against `$HOME/scratch` or `data/` without explicit user
  confirmation. No project-wide formatter, linter, Makefile, or CI workflow
  is configured; match surrounding style.

## Code conventions

- Resolve datasets and columns through `datasets.json` helpers; source
  `src/slurm_config.sh` in every HPC shell entry point.
- Preserve cross-language contracts: H5AD raw counts in `layers["counts"]`,
  Feather tabular sample identity, and RDS feature/label/distance bundles.
  Preserve row names and sample order; reject NA or mismatched identifiers.
- scITD may legitimately drop samples during cell-type filtering or tensor
  factorization. Validators must report dropped IDs and apply this exception
  only to scITD; other methods retain the current sample universe.
- Use `stop()`/`stopifnot()`, parallel `.errorhandling="stop"`, strict shell
  status gates, and required-column validation. Warn-and-skip is only for
  explicitly optional artifacts.
- Write artifacts through temporary files plus atomic rename, checkpoint
  per-sample work, and verify MD5 before `readRDS`. Never weaken ownership,
  provenance, schema, or checksum checks.
- Numbered scripts encode pipeline order; dataset-specific code stays under
  stage 2; shared helpers belong in `src/utils/`; shell variables are
  uppercase. Pass configuration/state explicitly through functions and CLI
  arguments rather than adding a global state manager.
- Retain sparse matrices and subset before densifying. Use namespaced R calls
  where established and lazy imports for optional heavyweight Python methods.
- Document scientific invariants and non-obvious scheduler behavior, not
  line-by-line mechanics. Save plans in `.agents/plans/` with a
  date-prefixed descriptive name; archive completed implementation plans
  under `.agents/plans/archive/`.

## Important files

- `datasets.json` — central dataset, metadata, and view contract.
- `pixi.toml`, `pixi.lock` — pinned environments and dependencies.
- `src/slurm_config.sh` — canonical HPC paths, interpreters, modules,
  resources, and retry ceilings.
- `src/utils/datasets_io.R`, `src/utils/py/datasets_io.py` — shared config
  access.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` — canonical H5AD
  preprocessing entry point.
- `src/5_run_benchmark_methods/benchmark_pipeline.R`,
  `benchmark_methods_r.R`, and `benchmark_submit_common.sh` — benchmark
  transforms, method contracts, submission, retry, sync, and checksums.
- `src/utils/scoring_metrics.R` — benchmark metrics.
- `src/utils/bash/setup_env_sbatch.sh` and `refresh_env.sh` — serialized
  environment mutation and smoke checks.
- `README.md` and `docs/ARCHITECTURE.md` — operator workflow and pipeline map.
- `notebooks/` — execute relevant code chunks as needed; do not knitr-render
  benchmark notebooks merely to produce their PDF outputs.

## Testing and QA

QA is focused and script-based; there is no aggregate test runner, CI gate,
coverage threshold, `pytest`, or `testthat` suite.

- Deterministic shell watchdog tests cover memory escalation, scheduler
  states, status files, notifications, retry ceilings, and fail-closed
  behavior.
- Focused stage submitter, watchdog, artifact-contract, chunk-preparation,
  benchmark-matrix, selection, ownership, and synchronization tests cover
  the corresponding contracts.
- For changes, run the narrow contract test and exercise the `_debug` path
  when appropriate. Add tests only for new observable behavior, boundaries,
  failure handling, or scientific invariants. Keep fixtures temporary and
  deterministic, and stub Slurm/NAS effects instead of requiring live
  infrastructure.