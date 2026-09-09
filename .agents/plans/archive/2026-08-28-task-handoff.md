# ECODA Task Handoff — 2026-08-28

## Objective

Complete the approved benchmark preprocessing consistency plan:

- enforce the raw-view `<500` cell-count rule;
- refresh the four affected ordinary benchmark datasets (Bassez, Lee, Smillie, Zhang);
- regenerate annotations and benchmark artifacts through the supported Stage 3–5 paths;
- validate ordering, finite values, checksums, provenance, and sample-universe contracts;
- release only a reviewed durable Stage 5 gate;
- then audit the remaining ordinary `benchmark_analysis` datasets and rerun only targeted invalid or missing method artifacts.

The approved plan is:
`local://1787689540000-benchmark-preprocessing-consistency-plan.md`.

## Scientific and repository constraints

- Never use biological labels (`Status`, `sample.origin`, `cond`, `Disease_Identity`, or equivalent) as preprocessing, HVG, normalization, batch-correction, embedding, or model covariates.
- Retain samples with at least 500 observations; remove samples with fewer than 500. Exactly 500 is retained.
- Preserve `layers["counts"]`, canonical sample order, row/column/distance alignment, atomic writes, checksum validation, and fail-closed behavior.
- scITD is the only accepted reduced-universe exception. All other methods must retain the complete source sample universe.
- `datasets.json` is ground truth and must not be changed without explicit approval.
- Do not manually create checksum sidecars or relabel failed artifacts as valid.
- Full-cohort work must use the checked-in `durable-hpc-gate-ecoda` profile.
- After a durable launch, use exactly one durable unbounded wait, then one terminal inspect with all emitted scheduler IDs, followed by reviewer approval only after a passing audit.
- Preserve `.gate/`, `.kilo/gates/`, logs, caches, rendered outputs, and temporary evidence outside Git commits.
- Do not stage the user's unrelated WIP.

## Completed implementation

### Sample-count invariant and registry

- Added `remove_low_cellcount_samples()` to `src/utils/py/preprocess_utils.py`.
- Integrated it into `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` before `process_view()`/`base_preprocessing()`.
- Added the raw/view audit tool `src/3_scrnaseq_preprocessing/1.0_audit_input_views.py`.
- Authoritative raw low-cell inventory:
  - Bassez: `BIOKEY_8_Pre=365`, `BIOKEY_25_Pre=296`.
  - Lee: `LB4180T=496`.
  - Smillie: `N58.LPB2=498`, `N19.LPB=485`, `N8.LPB=482`, `N12.LPB=441`, `N14.LPA=432`, `N12.LPA=243`.
  - Zhang: `Pre_P010_t=8`, `Pre_P018_b=437`.
- Restored Stephenson's ordinary `benchmark_analysis` view while retaining both batch-effect views.
- The approved registry change was pushed previously in `0865bb4`; do not modify the registry again without approval.

### Stage 3 and Stage 4

- Durable Stage 3 preprocessing gate completed and reviewed for Bassez, Lee, Smillie, and Zhang.
- Released H5AD summaries:
  - Bassez: 74,948 cells, 26 samples, minimum 569 cells/sample.
  - Lee: 113,865 cells, 27 samples, minimum 630.
  - Smillie: 197,706 cells, 64 samples, minimum 536.
  - Zhang: 201,442 cells, 31 samples, minimum 2,076.
- Stage 4 annotation refresh completed through the supported chunk/worker/merge path.
- Durable recovery/release evidence exists at:
  `.gate/ecoda_annotation_refresh_recovery_20260826T004041Z.*`.

### Benchmark artifact contracts

- R metadata/sample-order fixes:
  - `collapse_sample_metadata()` in `src/utils/seurat_utils.R`;
  - `prepare_mofa_metadata()` in `benchmark_methods_r.R`;
  - canonical pseudobulk ordering;
  - first-appearance factor levels in composition helpers;
  - canonical matrix/distance alignment in `align_result_samples()`;
  - canonical null-control labels and GloScope/GloProp label construction.
- Python fixes:
  - canonical sample-order and square-frame alignment;
  - finite numeric Feather validation;
  - MRVI, PILOT, QOT, and PILOT-GM alignment fixes;
  - PILOT-GM singleton covariance stabilization;
  - missing annotation values excluded from configured component-cardinality counting.
- H5AD source-sidecar repair path added to `src/5_run_benchmark_methods/1_submit_hpc_array.sh`; it validates scratch/NAS identity before recording a missing sidecar and never overwrites an invalid sidecar.

## Environment investigation findings

### Historical evidence

The archived investigation plan `.agents/plans/archive/1786440267018-setup-lazyload-integrity-check.md` records the original pattern:

- random cross-node `ENOENT` or empty reads after environment mutations;
- files healthy from the login node;
- only approximately one random task failing per array;
- failures across different nodes and different packages;
- the prior best explanation was stale NFS client views after mutations, not a package-specific trigger.

The HPC filesystem is currently reported as BeeGFS rather than generic NFS:

- `$HOME` and `.pixi`: `beegfs_home`;
- `$HOME/scratch/ECODA_paper`: `beegfs_scratch`.

There is no per-node “sync” operation. Nodes access the shared BeeGFS namespace through their own clients and metadata/data caches.

The Bamboo shell history also contains historical direct writer commands, including:

- manual removal of `HiTME`, `GloScope`, and `ProjecTILs` package directories;
- `pixi run -e py-cuda13 Rscript -e 'install.packages("abind", ...)'`;
- direct `pixi run -e py-cuda13 setup`;
- direct `pixi install` commands.

These are valid historical writer paths, but the history does not prove that each one overlapped a benchmark array.

### Confirmed Pixi/R activation writer

The conda `r-base` package owns:
`etc/conda/activate.d/activate-r-base.sh`.

That hook runs:

```text
R CMD javareconf
```

`strace` of `pixi run --as-is` showed `sed`/`mv` writes to:

```text
.pixi/envs/py-cuda13/lib/R/etc/Makeconf
.pixi/envs/py-cuda13/lib/R/etc/ldpaths
```

Therefore `--as-is` prevents solving/installing but does not make Pixi activation read-only.

The worker runtime was changed to the direct absolute Rscript binary in `src/slurm_config.sh`, bypassing activation. A direct-R `strace` showed no writes under the environment prefix.

### Confirmed systemic shared-prefix read failure

A direct-R 16-node concurrent load test (no Pixi activation, no installer) failed on five nodes:

```text
cpu027, cpu030, cpu032, cpu033, cpu037
```

The missing files were unrelated R files, including:

```text
methods/Meta/package.rds
utils/Meta/package.rds
S4Vectors/Meta/package.rds
RColorBrewer/R/RColorBrewer
pheatmap/Meta/package.rds
cli/Meta/nsInfo.rds
```

A clean direct-Python 16-node test failed on four nodes:

```text
cpu032, cpu035, cpu036, cpu040
```

Examples included missing:

```text
numpy/version.py
numpy/dtypes.py
numpy/_core/einsumfunc.py
sklearn/externals/array_api_compat/numpy/_info.py
sklearn/metrics/_plot/precision_recall_curve.py
```

The corresponding files exist from the login-node view. The conda-prefix inventory found:

```text
81,711 declared files
0 missing files
711 packages
```

The known `pilotgm` `networks` import problem was excluded from the clean Python result.

### Concurrency correlation

The same previously failing nodes passed with four concurrent direct-R workers and four concurrent direct-Python workers.

Additional evidence:

- CPU009 passed an eight-worker direct package-load stress test.
- A static 16-node `md5sum` test passed for one file on `beegfs_home` and one file on `beegfs_scratch`.
- A sequential direct package retest on CPU027 passed.
- The 16-node direct-R failure happened even after bypassing the Pixi/r-base activation writer.

This rules out one permanently bad node and one problematic R package. The best current explanation is a concurrency-dependent BeeGFS metadata/lookup visibility problem when many processes cold-start against the large shared HOME prefix. A later lookup can refresh a client’s view, making the package appear to “repair itself”; `library()` itself is not reinstalling anything.

The 52-node and 40-node tests were attempted but could not allocate because of partition node limits/reservations; they were cancelled rather than left queued. The 16-node tests are already sufficient to establish the systemic pattern.

### Important limitation

The direct R runtime and atomic mutator lock are **partial mitigations**, not a complete root-cause fix:

- the atomic lock protects competing environment mutators;
- the no-active-job check prevents mutation while the check observes jobs absent;
- a worker can still be submitted after that check and overlap a later mutation;
- BeeGFS is configured with `tuneUseGlobalFileLocks = false`, so do not assume advisory `flock`/fcntl is globally visible across clients.

Do not claim the benchmark environment is fully fixed yet.

## Environment changes completed

Current source-of-truth commits after the earlier `63c331b` baseline:

- `bd0f358` — direct Rscript worker runtime to avoid r-base activation writes.
- `66a8164` — atomic fail-closed environment mutation directory lock.
- `566745c` — Slurm query failures now abort before `pixi install`/setup.
- `25a855d` — extracted guarded `src/utils/setup_r_packages.R`, removed `[tasks.setup]` from `pixi.toml`, migrated both environment wrappers and documentation.

Current local `HEAD` and Bamboo `origin/master` are expected to be:

```text
25a855d8c1add5c37a87c9b2d81307593516f0fd
```

The guarded R setup path is now:

```text
src/utils/bash/refresh_env.sh
src/utils/bash/setup_env_sbatch.sh
  -> pixi install (under atomic lock and successful squeue check)
  -> direct PIXI_RSCRIPT src/utils/setup_r_packages.R
     with ECODA_ENV_MUTATION_GUARD=1
  -> direct-R smoke checks
```

Direct `pixi run setup` is no longer a task and is no longer documented as authoritative. The setup script refuses to run unless `ECODA_ENV_MUTATION_GUARD=1` is supplied by the guarded wrappers.

## Verification completed

Passed local/remote checks include:

```text
bash tests/test_env_mutation_lock.sh
bash tests/test_r_environment_preflight.sh
bash tests/test_ecoda_run_common.sh
bash tests/test_benchmark_matrix_submitter.sh
bash tests/test_benchmark_matrix_watchdog.sh
bash tests/test_h5ad_preflight.sh
bash src/5_run_benchmark_methods/test_oom_retry.sh
```

Also passed:

- `setup_r_packages.R` parse check;
- unguarded setup rejection check;
- direct compute-node R package preflight;
- package-load tracing showing activation writes for Pixi and no prefix writes for direct R.

The lock regression test covers:

- active lock rejection;
- fail-closed behavior for invalid/stale lock state;
- concurrent acquisition with exactly one winner;
- failed `squeue` query rejection;
- active-job rejection and current-job exclusion.

## Durable gate history and current state

### Old retry5 gate

`ecoda_benchmark_methods_repair_retry5_20260828T071924Z`:

- scheduler IDs included `4364884`, `4364888`/`4364889`, `4364890`/`4364893`, `4364894`/`4364895`, `4364896`/`4364897`, `4364898`/`4364899`, `4364900`/`4364901`, `4364902`/`4364903`;
- all emitted IDs and array children were later inspected as terminal;
- composition failed because of missing `DelayedArray.rdb`/`arrow.rdb` on early tasks;
- the run had no durable terminal status when its runner was stopped;
- its 25 orphaned owner records were marked `FAIL` only after all IDs were terminal and no runner/writer remained;
- evidence was preserved at:
  `.gate/ecoda_benchmark_methods_repair_retry5_20260828T071924Z.owner-recovery.tsv` and the remote gate logs.

### Gate that was stopped for missing `--force`

`ecoda_benchmark_methods_repair_direct_rscript_20260828T082857Z` launched the exact selection without `--force` and skipped already-valid rows. It was stopped before release to prevent a mixed old/new artifact scope. Its orphaned owners were recovered after terminal job checks. Do not reuse this gate.

### Latest forced gate

`ecoda_benchmark_methods_repair_forced_direct_rscript_20260828T083643Z` used direct Rscript and `--force`.

Remote run root:

```text
/home/users/h/halterc/scratch/ECODA_paper/_ecoda_runs/stage5_20260828103720_2696790
```

It emitted:

```text
H5AD preflight: 4364986
R preflight:    4364990
Arrays:         4364991, 4364996, 4364998, 4365000, 4365002,
                4365004, 4365006, 4365011, 4365013
Watchdogs:      4364992, 4364997, 4364999, 4365001, 4365003,
                4365005, 4365010, 4365012, 4365014
Aggregate gate visible in scheduler output: 4365015
```

Composition array `4365000` failed Bassez on `cpu009` with:

```text
cannot open .../SparseArray/R/SparseArray.rdb
```

Lee, Smillie, and Zhang composition tasks on CPU010 succeeded.

The gate was cancelled because GPU jobs were blocked for approximately a day and the composition failure made success impossible. Its durable status is `FAILED` with exit code `143`, and the one terminal inspect was completed with all recorded/preflight IDs plus visible aggregate ID. No reviewer approval exists and no artifacts from this gate are releasable.

The latest 34 orphaned owner records were recovered only after the terminal audit and no-active-job/writer checks. The gate’s owner-recovery evidence is:

```text
/home/users/h/halterc/scratch/ECODA_paper/gates/ecoda_benchmark_methods_repair_forced_direct_rscript_20260828T083643Z.owner-recovery.tsv
```

There is currently no benchmark gate, no active Slurm job, no environment lock, and no tmux session on Bamboo.

## Local WIP that must not be staged

The previous status check showed unrelated user work:

```text
modified:   notebooks/benchmark_analysis.rmd
untracked:  .agents/plans/1787783466080-onboarding-uncorrected-batch-plan.md
untracked:  .agents/plans/archive/1787736504341-supp-fig-1-combination-plan.md
untracked:  .agents/plans/review_findings.md
untracked:  .gate/
```

Do not reset, clean, or commit these paths without explicit user direction.

## Exact next steps

1. **Do not launch Stage 5 yet.** The shared HOME prefix remains unsafe under high concurrent cold-start imports.
2. Build a genuinely isolated runtime A/B test. The previous copy attempt was invalid because the packaged R binary hardcodes the HOME prefix; verify `R.home()` and Python `sys.prefix` point to the local/container runtime before loading packages.
3. Preferred A/B design:
   - Arm A: direct R/Python binaries from the shared HOME prefix, using the already demonstrated 16-node concurrent test.
   - Arm B: the exact same package loads from an immutable node-local copy or Apptainer image.
   - Repeat on at least 16 nodes and record host, runtime root, package path, and exact failure. Do not install anything in either arm.
4. If Arm B passes, adopt an immutable node-local/container runtime for workers. Do not rely only on a startup delay; delay is a mitigation, not proof of correctness.
5. If Arm B fails, provide HPC administrators the exact evidence that both direct R and direct Python produce random `ENOENT` on `beegfs_home` under fan-out while the login conda inventory is complete. Ask them to inspect BeeGFS client/metadata-server logs and cache behavior on the failing nodes.
6. Do not add more preflight/error-handling layers unless the A/B test fails to identify a usable runtime fix. The existing compute preflight is diagnostic, not the root fix.
7. Keep the direct Rscript runtime, guarded setup script, atomic lock, and fail-closed `squeue` check. Do not introduce cross-node `flock` without proving global lock semantics because `tuneUseGlobalFileLocks=false`.
8. Once the immutable-prefix/FS issue is resolved, prepare a **new** durable gate ID using the exact 30-row selection file:

   ```bash
   cd "$HOME/ECODA_paper" &&
   source src/slurm_config.sh &&
   ./src/5_run_benchmark_methods/1_submit_hpc_array.sh \
     --selection-file "$HOME/ECODA_paper/.gate/stage5_repair_selection_20260828.tsv" \
     --force
   ```

   The durable profile is:
   `.agents/skills/durable-hpc-gate-ecoda/references/profile.json`

   Profile digest:
   `4daca1336ef9d294e4e4db229916d58a86a800754ede5566d9875022c1d29ae0`

   Remote workdir:
   `/home/users/h/halterc/ECODA_paper`

   Serialization group:
   `ecoda-benchmark`

9. For the next gate, confirm before launch:
   - Bamboo `HEAD == origin/master == 25a855d...`;
   - no active Slurm jobs;
   - no environment lock/writer;
   - current selection file has 30 rows;
   - immutable-prefix A/B has passed;
   - the first scheduler job is the compute-node R preflight.
10. After a passing reviewed four-dataset gate, audit every remaining ordinary `benchmark_analysis` dataset read-only and generate exact per-dataset/per-method rerun selections. Do not use broad `--force` runs.
11. Only after the artifact/gate review is closed, organize remaining WIP into deliberate commits. Keep gate/evidence/log/cache/rendered-output paths outside Git.

## Current completion status

The implementation and released-artifact verification tasks remain blocked, not complete:

- Stage 3/4 work is complete and reviewed.
- Stage 5 has not reached a successful reviewed release.
- The shared runtime failure is reproduced in both R and Python.
- The durable next action is the immutable-prefix/filesystem A/B test, followed by HPC/BeeGFS diagnosis or a proven isolated runtime.
