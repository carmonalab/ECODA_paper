# Bassez Annotation Fix and Benchmark Regeneration

> **Authoritative execution status (2026-08-25): COMPLETED_WITH_AUDIT_NOTE.** Bassez annotation correction and all planned benchmark/result calculations are complete; final artifact audit, local mirror refresh, canonical render, and Bamboo checksum verification passed. This `.kilo` plan remains the sole source of truth.
> **Composition retry checkpoint:** gate `ecoda_composition_retry_20260825T160000Z` launched once with Bamboo source `e4a8d3f541cf744f7e8879ff8e17ef5a1bc7d06e`; prepare array/watchdog `4351259`/`4351269` completed, composition array/watchdog `4351288`/`4351289` failed, and the wrapper exited `1` without NAS synchronization. Wait evidence is `.kilo/gates/ecoda_composition_retry_20260825T160000Z.wait.json`; terminal audit is `.kilo/gates/ecoda_composition_retry_20260825T160000Z.inspect.json`.
> **Composition failure diagnosis:** `get_pb_deseq2()` produced pseudobulk sample rownames such as `BIOKEY-2-Pre`, while the canonical h5ad obs labels used `BIOKEY_2_Pre`; `create_result_bundle()` therefore rejected all labels for the affected datasets. Commit `4f9610b` adds one-to-one exact/standardized sample-ID reconciliation, passes the focused regression test, and is deployed to Bamboo. A fresh composition retry is required; no dependent benchmark gate may start before its terminal inspection and reviewer approval.
> **Composition retry 2 checkpoint:** gate `ecoda_composition_retry2_20260825T143410Z` launched once against Bamboo source `4f9610bdf0cc7397653d945105d458d6948d7dbc`; exact wrapper is unchanged. Local manifest is `.kilo/gates/ecoda_composition_retry2_20260825T143410Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_composition_retry2_20260825T143410Z.wait.json`. No dependent benchmark gate may start before terminal inspection and reviewer approval.
> **Composition retry 2 terminal result:** scheduler arrays/watchdogs `4351389`/`4351390` and `4351457`/`4351458` all completed `0:0`; the wrapper log records 20 task logs merged, checksum sidecar generation, and successful NAS synchronization. The durable inspect correctly rejected release because Bamboo HEAD changed from `4f9610b` at launch to concurrent preprocessing-watchdog commit `4d6c2da` before terminal audit, and it observed a transient checksum-root path failure. The source diff is disjoint from benchmark code; a documented `--sync-only` recovery gate will re-audit the completed outputs under the stabilized checkout before pseudobulk proceeds.
> **Composition sync-only recovery checkpoint:** gate `ecoda_composition_sync_20260825T153000Z` launched once against stabilized Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`, using documented `--sync-only 4351457,4351458` recovery and pre-attached scheduler IDs. Local manifest is `.kilo/gates/ecoda_composition_sync_20260825T153000Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_composition_sync_20260825T153000Z.wait.json`.
> **Composition sync-only recovery 2 checkpoint:** the ECODA profile was corrected to validate the actual Bassez registry contract and prefer the canonical `benchmark/checksums.md5` root. Gate `ecoda_composition_sync2_20260825T153500Z` launched once at stable Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`; local manifest is `.kilo/gates/ecoda_composition_sync2_20260825T153500Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_composition_sync2_20260825T153500Z.wait.json`.
> **Composition sync-only recovery 2 audit:** accounting, repository HEAD, Bassez registry, and terminal roots passed; the profile's NAS checksum command returned exit `1` with no output during the single audit query. An immediate read-only `md5sum -c --quiet` from `$NAS_TARGET_DIR/benchmark` now returns `0`, so a final explicit sync-only gate will re-audit without changing benchmark calculations.
> **Composition sync-only recovery 3 checkpoint:** gate `ecoda_composition_sync3_20260825T153800Z` launched once at stable Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`; it uses the same documented sync-only command and pre-attached successful scheduler IDs. Local manifest is `.kilo/gates/ecoda_composition_sync3_20260825T153800Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_composition_sync3_20260825T153800Z.wait.json`.
> **Composition sync-only recovery 4 checkpoint:** the profile NAS audit now sources `src/slurm_config.sh` before using `NAS_TARGET_DIR`. Gate `ecoda_composition_sync4_20260825T154000Z` launched once at stable Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`; local manifest is `.kilo/gates/ecoda_composition_sync4_20260825T154000Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_composition_sync4_20260825T154000Z.wait.json`.
> **Composition release checkpoint:** gate `ecoda_composition_sync4_20260825T154000Z` terminal accounting for `4351457`/`4351458` passed, configured-root/checksum/registry/HEAD audits passed, and reviewer `Luna Max` approved at `2026-08-25T15:19:07Z`; evidence is `.kilo/gates/ecoda_composition_sync4_20260825T154000Z.review.json`. Composition outputs are released for the next serialized benchmark gate.
> **Bassez pseudobulk checkpoint:** gate `ecoda_bassez_pseudobulk_20260825T152000Z` launched once with the exact plan wrapper at Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`, dependency-bound to the reviewed composition gate. Local manifest is `.kilo/gates/ecoda_bassez_pseudobulk_20260825T152000Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_bassez_pseudobulk_20260825T152000Z.wait.json`.
> **Bassez pseudobulk release checkpoint:** gate `ecoda_bassez_pseudobulk_20260825T152000Z` terminal accounting for `4351531`/`4351532` and `4351543`/`4351544` passed, configured-root/checksum/registry/HEAD audits passed, and reviewer `Luna Max` approved at `2026-08-25T15:34:35Z`; evidence is `.kilo/gates/ecoda_bassez_pseudobulk_20260825T152000Z.review.json`.
> **Transform/zero-imputation checkpoint:** gate `ecoda_transzeroimp_20260825T154500Z` launched once at Bamboo HEAD `4d6c2daece92bc3dc711b21ff672f596448c7d0`, dependency-bound to the reviewed Bassez pseudobulk gate. Exact wrapper is `cd "$HOME/ECODA_paper" && source src/slurm_config.sh && ./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --analysis trans,zeroimp --force`; local manifest is `.kilo/gates/ecoda_transzeroimp_20260825T154500Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_transzeroimp_20260825T154500Z.wait.json`.
> **Transform/zero-imputation failure diagnosis:** trans array `4351566` completed all tasks, but zeroimp array `4351567` failed all ten datasets before NAS synchronization because `zCompositions::cenGeoMean()` called bare `survreg`; `imports_worker_transzeroimp.R` did not attach `survival`. Commit `664dd5f` adds the dependency to the worker attach list, the local import smoke test passed, and the fix is deployed to Bamboo. A fresh serialized retry is required; no partial trans sync was accepted.
> **Zeroimp retry checkpoint:** gate `ecoda_transzeroimp_retry_20260825T160000Z` launched once at Bamboo HEAD `664dd5f44e7399a68efb99c7872ea83cd7f0e4e2`, dependency-bound to the reviewed Bassez pseudobulk gate. It uses the exact transformation/zero-imputation wrapper; local manifest is `.kilo/gates/ecoda_transzeroimp_retry_20260825T160000Z.manifest.json`; one durable waiter is armed at `.kilo/gates/ecoda_transzeroimp_retry_20260825T160000Z.wait.json`.
> **Transform/zero-imputation release checkpoint:** gate `ecoda_transzeroimp_retry_20260825T160000Z` terminal accounting for arrays `4351616`/`4351617` passed, configured-root/checksum/registry/HEAD audits passed, and reviewer `Luna Max` approval is recorded in `.kilo/gates/ecoda_transzeroimp_retry_20260825T160000Z.review.json`. The prior zeroimp failure was isolated to the missing `survival` attach and is superseded; all ten trans and zeroimp results synchronized to NAS.
> **Final completion checkpoint:** Bassez h5ad audit passed shape `(75609, 19366)`, required metadata/layer/HVG/PCA keys, zero sentinel/missing `cellSubType`, and 40 categories; scratch/NAS MD5 remained `2c103d6a2b22e4519562e7d43b9c9822`. All 20 maintained Bassez Python feathers opened with 28 unique samples and expected distance/scPoli shapes; Bassez composition/pseudobulk bundles opened with required HR keys and nonmissing labels; all 12 configured result datasets have readable composition/trans/zeroimp RDS bundles; `checksums.md5` passed; `execution_times.feather` has 1,364 rows and zero duplicate `(dataset, method)` keys. Local benchmark mirror and Bassez h5ad were refreshed; canonical `analysis_scope <- "preprint"` render completed and required publication PDFs are nonempty. Final Bamboo checkout is `a2d1849fb1fd9f9f70612167e9bb8760367eb81c`, and the final NAS checksum audit passed.
> Source hardening is complete in commit `7a2d591`,
> which is pushed to `origin/master` and present in the bamboo clone. The
> canonical preflight passed. RDS patch job `4334702` completed on `shared-cpu`
> with exit `0:0`; it filled 162,753 of 226,454 missing `cellSubType` values,
> changed subtype cardinality from 33 to 40, and atomically installed the
> validated staged RDS. Compute step `4334703` then atomically patched the
> canonical processed `.h5ad`; validation step `4334704` confirmed shape
> `(75609, 19366)`, zero invalid `cellSubType` values, and 40 categories. The
> **Section 3 artifact equality: PASS.** The canonical scratch and NAS files
> each have size `3,324,722,085` bytes and MD5
> `2c103d6a2b22e4519562e7d43b9c9822`. The sync log names job `4334705`, but
> `sacct -j 4334705 -X` returned unrelated `4318816_28`; this provenance
> discrepancy remains an audit fact. The user explicitly accepted that
> discrepancy as safe to proceed and authorized continuation into Section 4.
> The prior “no benchmark wrapper started” statement records the verified
> prelaunch state and is superseded by the single gate below.
> **Audit note (retained):** The exact Bassez Python methods wrapper was launched once
> after reconciling Bamboo HEAD `7a2d5918f74ed80207e98925543282f39c186616`,
> with no duplicate active gate and no later gate started:
> `./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Bassez --methods scpoli,pilot,qot,pilotgm --force`.
> Array IDs are `scpoli=4334711`, `pilot=4334713`, `qot=4334715`,
> `pilotgm=4334717`; watchdog IDs are `scpoli=4334712`, `pilot=4334714`,
> `qot=4334716`, `pilotgm=4334718`. The persistent tmux session is
> `bassez_py_gate_20260824T022421Z`; durable log and exit-status files are
> `/home/users/h/halterc/ECODA_paper/logs/bassez_py_gate_20260824T022421Z.log`
> and `/home/users/h/halterc/ECODA_paper/logs/bassez_py_gate_20260824T022421Z.status`.
> The durable launcher status is `STATE=FAILED`, `EXIT_STATUS=1`,
> `END_UTC=2026-08-24T03:12:35Z`. Exactly one accounting query covered all
> eight recorded IDs; every task row (`4334711_1`, `4334713_1`, `4334715_1`,
> `4334717_1`) and watchdog (`4334712`, `4334714`, `4334716`, `4334718`) is
> `COMPLETED` with `0:0`. The final log records all four watchdogs `STATE=OK`,
> checksum-sidecar finalization (`Wrote benchmark/checksums.md5`), four
> execution logs merged into 1,364 rows, NAS synchronization, and post-sync
> per-task cleanup. On NAS, read-only `md5sum -c --quiet checksums.md5` returned
> `CHECKSUMS_OK`. A read-only feather audit opened all 20 expected Bassez Python
> artifacts: each has 28 rows; distance outputs have 29 columns; scPoli outputs
> have `N+1` columns for dimensions 2, 3, 5, 10, and 15.
> The launcher transcript shows the run was sent to a fresh `bash -l` without
> exported `LOG`, `STATUS`, `SESSION`, or `CHANNEL`: `tee -a ""` failed at
> capture setup and the final empty `${STATUS}` redirect emitted
> `bash: : No such file or directory`. Therefore durable `FAILED/1` and the
> launcher exit `1` are monitoring/status-write evidence, not proof of the
> wrapper result; the actual wrapper return code was captured only in the dead
> shell and is unrecoverable. The unresolved wrapper rc is retained as this
> audit note, while the independently verified scheduler, finalization, sync,
> checksum, and artifact evidence supports completion.
> **Reviewer verdict: PASS as COMPLETED_WITH_AUDIT_NOTE.** The user explicitly
> accepted the verified outputs (“Accept verified outputs”) despite the
> unrecoverable wrapper rc audit note. The next pending gate is **Composition,
> all datasets**, which may proceed serially; no other later gate has started.
> **Transport recovery checkpoint (2026-08-24): STOP — composition gate remains ambiguous; do not relaunch, inspect, query Slurm, or start a later gate.** The sole harness waiter `bg_7` ended after a remote completion-event read failed with `ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Undefined error: 0` (`.kilo/gates/ecoda_composition_20260824T093714Z.wait.json`, discrepancy `d0001` / `completion_transport_error`). One explicit `ssh bamboo 'printf %s "$HOME"'` recovery check reached `/home/users/h/halterc`; read-only durable `reconcile` and `status` then completed without mutation.
> Reconcile evidence is `AMBIGUOUS` with `identity_match=true`, remote manifest identity matching gate `ecoda_composition_20260824T093714Z` and command digest `1b45b492c60d990bfa8532093a4ac7c38c3f55e8237c3fd6323c0b2a84de0e28`, but no remote terminal status (`remote_status=null`). The exact runner process remains PID `3363302`; the named tmux session is present, but its pane commands include `"bash -l"` before the runner and fail the exact `pane_start_command` match (`d0002` / `tmux_runner_mismatch`). Status evidence records local `PRELAUNCH_STOP`, `terminal=true`, `remote_status=null`, and `release_eligible=false`. This is transport/identity ambiguity, not gate evidence of completion.
> The durable wrapper log was read without scheduler access. It emitted array/watchdog IDs: prepare_pseudobulk array `4334869`, prepare_pseudobulk watchdog `4334870`, composition array `4334887`, and composition watchdog `4334888`. IDs are recorded here for handoff only; no `inspect` or accounting query is legal until a later explicit reconciliation resolves the ambiguity.
> Gate ID: `ecoda_composition_20260824T093714Z`; serialization group: `ecoda-benchmark`; profile: `.agents/skills/durable-hpc-gate-ecoda/references/profile.json`; local manifest: `.kilo/gates/ecoda_composition_20260824T093714Z.manifest.json`.
> Exact wrapper: `cd "$HOME/ECODA_paper" && source src/slurm_config.sh && ./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --methods composition --force`.
> Persistent tmux session: `ecoda_composition_20260824T093714Z`; remote manifest: `/home/users/h/halterc/scratch/ECODA_paper/_durable_gates/ecoda_composition_20260824T093714Z/manifest.json`; remote runner: `/home/users/h/halterc/scratch/ECODA_paper/_durable_gates/ecoda_composition_20260824T093714Z/runner.sh`.
> Remote log: `/home/users/h/halterc/scratch/ECODA_paper/_durable_gates/ecoda_composition_20260824T093714Z/runner.log`; remote status: `/home/users/h/halterc/scratch/ECODA_paper/_durable_gates/ecoda_composition_20260824T093714Z/status.json`; completion event: `file:/home/users/h/halterc/scratch/ECODA_paper/_durable_gates/ecoda_composition_20260824T093714Z/done`. Recovery result records: `.kilo/gates/ecoda_composition_20260824T093714Z.reconcile.json` and `.kilo/gates/ecoda_composition_20260824T093714Z.status.json`.

> **Composition failure checkpoint (2026-08-24; evidence-only, no rerun or remote mutation):** The existing durable runner log records prepare_pseudobulk array `4334869` and watchdog `4334870` completing before composition array `4334887` and watchdog `4334888`. The existing watchdog status `/home/users/h/halterc/scratch/ECODA_paper/_benchmark_watchdog/4334888.status` is `STATE=FAIL`, with all 11 tasks (Adams, Bassez, Gongsharma_cmv_young_males, Kfoury, Kim, Lee, Pelka, Smillie, Stephenson, Wu, Zhang) `FAILED (1:0)` and NAS sync skipped. Representative Bamboo worker stderr `5_benchmark_r_composition_4334887_1.err` and task stderr files for tasks 2, 3, 4, and 11 have the same exact error: `load_pb_variants(NULL, ...)` reports all six variants (`schvg2000, hvg2000, hvg500, hvg2000_bl, hvg1000, hvg3000`) missing and no Seurat fallback. The prepare task-1 stderr records all six Adams files saved, and the read-only Bamboo cache listing contains all six files for every failed dataset.
>
> **Diagnosed root cause and local fix:** Bamboo's source still has the old composition branch calling `load_pb_variants(..., force = force)` with `seurat = NULL`; `pb_variants_missing()` intentionally treats `force=TRUE` as every pseudobulk variant missing. Since `--force` was supplied to recompute composition result bundles, every obs-only composition worker converted the prepared cache into an impossible no-Seurat fallback and exited `1`. Local source now routes composition through `load_composition_pb_variants()`, which always reuses the prepared pseudobulk cache (`force=FALSE`) while `run_composition_methods_hpc(..., force=force)` continues to force composition result recomputation. Focused temporary-fixture regression covers both distinctions and passes with `pixi run Rscript tests/test_bassez_and_benchmark_regressions.R`; only expected optional-cell-type warnings were emitted. Local changed files are `src/5_run_benchmark_methods/benchmark_hpc_utils.R`, `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R`, and `tests/test_bassez_and_benchmark_regressions.R`. Bamboo remains on the pre-fix source and must receive these files (or the containing commit) before a future composition-only rerun; no deployment, scheduler query, retry, or NAS sync was performed here.


## Verdict

The proposed scope is broadly correct, but the original execution plan is unsafe as written:

- The four benchmark submitters share `execution_times.feather`, `checksums.md5`, and whole-tree NAS synchronization. Serialize top-level submitters.
- `composition` and `pseudobulk` also both auto-run `prepare_pseudobulk`, writing shared caches/logs.
- Stage 5 reads only `${HPC_SCRATCH_DIR}/Bassez/output/...`; remove the unsupported fallback path.
- Directly overwriting the `.h5ad` is not atomic. Validate a same-directory temporary file before `os.replace`.
- Benchmark synchronization does not persist `${HPC_SCRATCH_DIR}/Bassez/output/`; sync the patched dataset separately.
- `--force` is method-wide, not high-resolution-only.
- Commit `76be859` also selected scPoli dimensions 15 as the default. Audit those artifacts for all datasets.
- The notebook reads the local `data/` mirror. Refresh it before rendering.
- Keep the canonical render at `analysis_scope <- "preprint"`; an all-scope render writes the same sacred figure filenames.

## Decisions

- Use the fast metadata-only path now: patch the staged RDS for future preprocessing and patch the existing processed Bassez `.h5ad` without rerunning Stages 3/4.
- Do not modify `datasets.json`; its Bassez columns and canonical output contract are correct.
- Do not pass biological labels into preprocessing/normalization.
- Do not rerun MrVI, scITD, GloScope, or MOFA for Bassez. They do not consume `cellSubType` (scITD uses `cellType`).
- Rerun Bassez scPoli/PILOT/QOT/PILOT-GM-VAE and pseudobulk. Their wrappers also recompute unaffected variants because no resolution/combo selector exists.
- Rerun `composition`, `trans`, and `zeroimp` for every benchmark dataset for the sentinel-filter rollout.
- QOT/PILOT-GM-VAE remain outside the preprint render, but regenerate their maintained Bassez artifacts.
- The prior-session contract review confirms the canonical Bassez paths and
  `Sample`/`expansion` metadata contract. A full preprocessing rerun would
  require explicit invalidation of the cached `<stem>_raw.h5ad`; the approved
  fast path does not rerun preprocessing and therefore does not touch that
  cache.

## 1. Source Hardening

Before HPC execution, make and validate these minimal changes:

1. `src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh`: use the exact `SLURM_JOB_ID`/`scontrol show job` `Command=` spool-recovery pattern required by `AGENTS.md`.
2. `src/2_dataset_specific_preprocessing/1.6.1_fill_bassez_cellsubtype.R`:
   - Treat actual NA plus stripped `""`, `"NA"`, `"nan"`, `"None"`, and `"Unknown"` as missing.
   - Fail if `cellType` is absent or invalid on any replacement row.
   - Preserve valid subtype values exactly.
   - Write/reopen a same-directory temporary RDS, validate cell count and metadata, then atomically rename it over the staged scratch input.
3. `src/5_run_benchmark_methods/benchmark_pipeline.R`: load `checksums.md5` from `dirname(path_results_nas)`. Finalization writes it at `benchmark/checksums.md5`, while the loader currently looks under `benchmark/results/` and silently disables verification.
4. Add focused tests where existing test conventions permit: mixed factor/string sentinels, preservation of valid labels, rejection of invalid fallback values, and checksum lookup from a `benchmark/results/` path. Do not run a full cohort pipeline.
5. Commit/push these changes and update `~/ECODA_paper` on `bamboo`. Do not stage unrelated onboarding work.

## 2. HPC Preflight

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
git rev-parse HEAD
git merge-base --is-ancestor 76be859cb3cc98cf3c5feb1f472a70a276bce4f8 HEAD
squeue -u "${USER}"
ls -lh "${HPC_SCRATCH_DIR}/Bassez/data/BassezA_2021_33958794whole.rds"
ls -lh "${HPC_SCRATCH_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad"
ls -lh "${NAS_TARGET_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad"
```

Gates:

- `HEAD` contains `76be859` and Section 1's hardening commit.
- No preprocessing, annotation merge, or benchmark job is writing Bassez/shared benchmark files.
- Canonical scratch and NAS files exist. The NAS copy remains rollback until the new scratch file is validated/synced.
- Do not accept `${HPC_SCRATCH_DIR}/BassezA_2021_33958794/...`; workers do not read it.

## 3. Patch Bassez

### Staged source RDS

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
sbatch --partition="${SLURM_PARTITION_BENCHMARK_CPU}" \
  src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh
```

The explicit partition is required: an unqualified submission can inherit the
15-minute `debug-cpu` default even though this wrapper requests one hour. The
completed fast-path submission was job `4334702`; do not resubmit it.

For a future submission, record the job ID. Do not poll in a loop. After
notification, inspect once:

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
sacct -j <RDS_PATCH_JOB_ID> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

Require `COMPLETED`/`0:0`. This RDS establishes future reproducibility but is not consumed by the current fast path.

### Canonical processed `.h5ad`

Run an inline `${PYTHON_BIN}` patch under `srun` on a compute node:

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
srun --partition="${SLURM_PARTITION_BENCHMARK_CPU}" --mem=64G --cpus-per-task=4 \
  "${PYTHON_BIN}" -c '<atomic validated patch described below>'
```

The patch must:

1. Open only `${HPC_SCRATCH_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad`.
2. Require `Sample`, `expansion`, `cellType`, `cellSubType`; `layers["counts"]`; `var["hvg_rank"]`; and `obsm["X_pca_benchmark_analysis_hvg{1000,2000,3000}"]`.
3. Record shape, obs/var indices, metadata columns, layer/obsm keys, and valid original subtype values.
4. Define missing values as actual NA or stripped strings in `{"", "NA", "nan", "None", "Unknown"}`.
5. Fail if replacement rows have invalid `cellType`; replace only missing subtypes; preserve valid subtypes exactly; restore categorical dtype.
6. Write a unique same-directory temporary `.h5ad`, reopen it, verify all structural invariants and zero invalid subtypes, then install it with `os.replace`.
7. Remove only the unique temporary file on failure. Print replacement count and before/after subtype cardinality. Zero replacements are acceptable when final invariants pass.

### Persist dataset output

Benchmark wrappers do not sync dataset outputs. In persistent `tmux` on the login node:

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
mkdir -p "${NAS_TARGET_DIR}/Bassez/output"
rsync -rlptDv \
  "${HPC_SCRATCH_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad" \
  "${NAS_TARGET_DIR}/Bassez/output/"
md5sum \
  "${HPC_SCRATCH_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad" \
  "${NAS_TARGET_DIR}/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad"
```

Require equal hashes before Stage 5.

## 4. Regenerate Benchmarks Serially

Each wrapper internally parallelizes, gates, merges logs, regenerates checksums, and syncs. Do not start the next command until the previous wrapper exits successfully.

### Bassez Python methods

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Bassez \
  --methods scpoli,pilot,qot,pilotgm \
  --force
```

This recomputes all configured resolutions/HVGs/parameter variants for those methods.

### Composition, all datasets

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --methods composition \
  --force
```

This auto-prepends forced `prepare_pseudobulk` and refreshes the complete composition bundle, not only author-HR/GloProp/frequency outputs.

### Bassez pseudobulk

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name Bassez \
  --methods pseudobulk \
  --force
```

This necessarily prepares Bassez pseudobulk again and recomputes plain/LR/HR/HVG/PCA variants. Keep wrapper guarantees rather than bypassing this redundancy.

### Transformations and zero imputation, all datasets

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh \
  --analysis trans,zeroimp \
  --force
```

Both analyses use `cell_type_high_res` through `get_ct_comp_df()`, so all datasets require regeneration.

Recovery rules:

- Record array and watchdog IDs.
- After a dropped login session, use the matching wrapper's `--sync-only <IDs>` with the same selectors and without `--force`.
- Prefer watchdog IDs for Python/R wrappers; transformation/zero-imputation uses array IDs.
- Inspect failures once with `sacct -j <ID> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS`; do not poll.

## 5. Audit scPoli Dimensions 15

Against `${NAS_TARGET_DIR}/benchmark/embeddings`, derive all non-debug benchmark datasets from `datasets.json` and require readable:

- `<DS>_hvg1000_highres_scpoli_dims15_embs.feather`
- `<DS>_hvg2000_highres_scpoli_dims15_embs.feather`
- `<DS>_hvg3000_highres_scpoli_dims15_embs.feather`
- `<DS>_hvg2000_lowres_scpoli_dims15_embs.feather` when a low-resolution column exists

Verify unique sample IDs, finite numeric values, and sample counts matching each metadata bundle. Filename existence is insufficient.

For each non-Bassez dataset with missing scPoli artifacts, run serially without `--force` so only missing combos execute:

```bash
source src/slurm_config.sh
cd "${PROJECT_ROOT}"
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh \
  --ds_name <DATASET_KEY> \
  --methods scpoli
```

## 6. Final Artifact Audit

On `bamboo`, require:

1. All recorded jobs are successful.
2. Bassez `.h5ad` still passes Section 3 invariants.
3. Bassez scPoli/PILOT/QOT/PILOT-GM-VAE feathers open with expected samples/dimensions.
4. `Bassez_composition.rds` and `Bassez_pseudobulk.rds` open and contain expected HR keys.
5. Every benchmark dataset has readable, newly regenerated `<DS>_composition.rds`, `<DS>_trans.rds`, and `<DS>_zeroimp.rds`.
6. `execution_times.feather` has current rerun rows and no duplicate `(dataset, method)` keys.
7. From `${NAS_TARGET_DIR}/benchmark`, `md5sum -c checksums.md5` succeeds for every listed RDS.

Rsync is intentionally non-destructive. Validate required current files explicitly; stale filenames may remain.

## 7. Refresh Local Mirror and Render

On macOS with `/Volumes/Shared/Projects/ECODA_paper` mounted:

```bash
cd /Users/christianhalter/Desktop/ECODA_paper
rsync -rlptDv "/Volumes/Shared/Projects/ECODA_paper/benchmark/" "data/benchmark/"
rsync -rlptDv \
  "/Volumes/Shared/Projects/ECODA_paper/Bassez/output/BassezA_2021_33958794whole_benchmark_analysis_ECODAprocessed.h5ad" \
  "data/Bassez/output/"
```

Keep `analysis_scope <- "preprint"`, then render locally:

```bash
cd /Users/christianhalter/Desktop/ECODA_paper
pixi run --as-is Rscript --vanilla -e "source('src/utils/load_all_functions.R')"
pixi run --as-is Rscript --vanilla -e "rmarkdown::render('notebooks/benchmark_analysis.rmd')"
```

Require no missing-bundle/method, checksum, or unreadable-feather warnings. Inspect without deleting/renaming:

- `plots/Figure_2_A_benchmark_mainmethods.pdf`
- `plots/Figure_2_B_MDS_Bassez.pdf`
- `plots/Supp_fig_2_methodparameter_screening.pdf`
- `plots/Supp_fig_15_benchmark_extended.pdf`
- `plots/Supp_fig_16A_transformations_anosim.pdf`
- `plots/Supp_fig_16B_zero_replacement_anosim.pdf`
- `plots/Supp_fig_16AB_trans_zero_anosim.pdf`

Check that Bassez has no NA/Unknown HR category, Figure 2A retains all datasets and dimensions-15 scPoli default, and Supp fig 16 contains all expected datasets.

If an all-scope diagnostic render is later needed, first parameterize output filenames/directories or render an isolated copy outside `plots/`; never overwrite canonical publication figures with another scope.

## 8. Completion and Failure Boundaries

- If `.h5ad` writing fails before `os.replace`, retain the original and delete only the unique temporary file.
- If post-replacement validation fails before NAS sync, restore scratch from unchanged NAS and stop.
- If dataset sync fails, block benchmarks until hashes match.
- If a wrapper fails, do not start the next wrapper or manually sync partial benchmark outputs.
- If rendering fails, regenerate the bad artifacts; never remove a `Figure*` or `Supp_fig*` output.
- After validation, archive this plan under `.agents/plans/archive/`, inspect status/diff/log, stage only intended source/tests/tracked publication outputs and the archived plan, then commit and push per `AGENTS.md`. Exclude local data mirrors and unrelated worktree changes.
