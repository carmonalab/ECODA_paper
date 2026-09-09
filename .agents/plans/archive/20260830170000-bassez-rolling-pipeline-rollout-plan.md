# Bassez-led rolling execution through pipelines 2–5

Plan slug: `20260830170000-bassez-rolling-pipeline-rollout`
Canonical plan artifact: `local://20260830170000-bassez-rolling-pipeline-rollout-plan.md`

## Context

Run the Bassez benchmark dataset through pipelines 2, 3, 4, and 5 as the control lane. After each Bassez stage reaches a reviewed terminal `OK`, promote the other datasets into that same pipeline stage without submitting Bassez again. Run only the uncorrected batch-effect view in this rollout. Pipeline 1 remains unchanged.

`datasets.json` is the sole eligibility and view authority. The production benchmark set is the eleven `use_for_benchmark=true` datasets excluding `_debug`: `Adams`, `Bassez`, `Gongsharma_cmv_young_males`, `Kfoury`, `Kim`, `Lee`, `Pelka`, `Smillie`, `Stephenson`, `Wu`, and `Zhang`. The production uncorrected batch-effect set is the twelve `use_for_batch_effect=true` datasets excluding `_debug`: `Joanito`, `Stephenson`, `CombinedPBMC`, `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Kidney_KPMP`, `Myocardial_infarction`, `Diabetes`, `Lupus_PBMC`, `Lung`, and `Parkinson`.

`Bassez`, `Lee`, `Smillie`, and `Zhang` are benchmark-only and each declares only `benchmark_analysis`; never request a batch-effect view for them. `Stephenson` belongs to both the benchmark and batch-effect sets and must be represented by two separately declared view rows.

The current remote evidence is not sufficient to treat old production gates as release evidence. Historical status files report `COMPLETED`, but their manifests remain `PREPARED` with `audit_state=NOT_STARTED` and no reviewer approval. Bassez/Lee/Smillie/Zhang benchmark H5ADs and sidecars are present but unverified against the current source. Myocardial's expected uncorrected output is present without a sidecar. CombinedPBMC's declared uncorrected output is missing; a differently named legacy artifact is not an acceptable substitute. Treat these artifacts as stale or unverified and force the relevant rebuilds.

## Approach

### 1. Establish the rollout identity and readiness ledger

1. Use the current Bamboo checkout and reviewed path-preserving runtime only after checking that Bamboo `HEAD` equals the runtime manifest `GIT_REVISION`. The reviewed image is currently `/home/users/h/halterc/scratch/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-path-preserving-c293.sif` with its matching `.manifest`; if either the checkout or manifest differs, stop and rebuild/review the image before launching any pipeline gate.
2. Confirm no active environment mutation lock or conflicting ECODA gate exists. Preserve `datasets.json`, `pixi.toml`, and `pixi.lock` unchanged.
3. Create explicit rollout selection records under the rollout's remote gate/run root. Use these row formats exactly:
   - Stage 3/4: `DATASET<TAB>VIEW`.
   - Stage 5 ordinary benchmark: `DATASET<TAB>benchmark_analysis<TAB>METHOD_OR_ANALYSIS`.
   - Stage 5 batch: `DATASET<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected`.
4. Record Pipeline 2 applicability explicitly. The current Stage 2 dispatcher has hooks only for `Gongsharma_cmv_young_males`, `CombinedPBMC`, `Joanito`, `_debug`, `Kfoury`, `Myocardial_infarction`, and `Bassez`. For this production rollout, Pipeline 2 is required for:
   - Bassez control lane: `Bassez`.
   - Hook-backed rest wave: `Gongsharma_cmv_young_males`, `Kfoury`, `Joanito`, `CombinedPBMC`, `Myocardial_infarction`.
   - `CombinedPBMC` automatically adds `gongsharma_cap` and submits `combinedpbmc` after `afterok:<gongsharma_cap>`.
5. Mark Pipeline 2 as not applicable—not passed—for benchmark datasets `Adams`, `Kim`, `Lee`, `Pelka`, `Smillie`, `Stephenson`, `Wu`, and `Zhang`, and for batch datasets `Stephenson`, `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Kidney_KPMP`, `Diabetes`, `Lupus_PBMC`, `Lung`, and `Parkinson`. Their downstream declared inputs still require current Stage 3 validation; old file presence does not establish current readiness.
6. Because current production gate manifests are unaudited and sidecar/provenance evidence is incomplete, use `--force` for every selected Pipeline 2, 3, 4, and 5 run below. Do not reuse the malformed CombinedPBMC legacy filename.

### 2. Durable gate contract for every rollout wave

Run every listed wave as its own checked-in `durable-hpc-gate-ecoda` gate from the Bamboo repository clone `$HOME/ECODA_paper`. Each gate must:

1. Use a unique gate ID, remote manifest, runner, log, status path, and serialization-group name.
2. Source `src/slurm_config.sh` before invoking the canonical submitter and set the reviewed `ECODA_RUNTIME_MODE=apptainer`, image, and manifest.
3. Run exactly one `prepare`, one `reconcile`, and one `launch`, then arm exactly one unbounded durable `wait`.
4. After terminal status, run exactly one first `inspect` with every scheduler array/watchdog/preflight/aggregate ID emitted by that wrapper, then one separate Luna Max reviewer approval inspect.
5. Treat a wave as passed only when its wrapper state, watchdog/owner/run state, selected sync, artifact contracts, checksums, source identity, terminal accounting, first audit, and reviewer approval all pass. A scheduler `COMPLETED` row without the saved audit and reviewer result is not a promotion barrier.
6. Never use `squeue`/`sacct` polling, manual partial synchronization, or a broad no-argument default selection. A `--sync-only RUN_ID` recovery is allowed only after the run-owned immutable manifests and terminal gates validate.

Sibling gates may run in parallel only after their declared predecessor gates are reviewed `COMPLETED`; give each sibling a distinct serialization group. Never submit a row in two sibling gates or include Bassez in a rest-wave selection.

### 3. Bassez Pipeline 2 control gate

Prepare and launch gate `B2-BASSEZ` with this exact scientific wrapper:

```bash
cd "$HOME/ECODA_paper" && \
export ECODA_RUNTIME_MODE=apptainer \
  ECODA_RUNTIME_IMAGE="$HOME/scratch/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-path-preserving-c293.sif" \
  ECODA_RUNTIME_MANIFEST="$HOME/scratch/ECODA_paper/_ecoda_runtime/ecoda-py-cuda13-path-preserving-c293.sif.manifest" \
  ECODA_RUNTIME_PROFILE=stage2 && \
source src/slurm_config.sh && \
bash src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Bassez --force
```

Require the `bassez_cellsubtype` hook output and its checksum/semantic validation to pass. Do not treat the old Bassez benchmark H5AD as the Pipeline 2 output; Pipeline 2's Bassez contract is the dataset-specific raw/prerequisite artifact defined by `step_outputs()`.

After B2-BASSEZ receives reviewer approval, launch the next two gates in parallel:

- `B3-BASSEZ`, defined in Section 4.
- `B2-HOOKED-REST`, defined in Section 5.

If B2-BASSEZ fails, do not launch either dependent gate. Preserve the failed run and create a new reconciled repair gate; do not manually delete or resynchronize artifacts.

### 4. Bassez Pipeline 3 and Pipeline 4/5 control lane

After B2-BASSEZ is reviewed `COMPLETED`, run `B3-BASSEZ`:

```bash
bash src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --datasets Bassez --views benchmark_analysis --force
```

Require the current `benchmark_analysis` H5AD contract, raw-count layer, required embeddings, declared observation columns, checksums, selected NAS sync, and run-owned terminal state. Do not require a second Bassez view.

Immediately after B3-BASSEZ is reviewed `COMPLETED`, launch `B4-BASSEZ`:

```bash
bash src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  --datasets Bassez --views benchmark_analysis --force
```

Require preparation, annotation, merge, annotation-column, sample-identity, checksum, and sync contracts. After B4-BASSEZ is reviewed `COMPLETED`, launch `B5-BASSEZ`:

```bash
bash src/5_run_benchmark_methods/1_submit_hpc_array.sh \
  --datasets Bassez \
  --methods gloscope,mofa,pseudobulk,composition,scitd,mrvi,scpoli,pilot,qot,pilotgm \
  --analyses trans,zeroimp \
  --force
```

B5-BASSEZ is the final benchmark control barrier. Require every selected method/analysis watchdog, aggregate gate, Feather/RDS result contract, execution log, sample order, checksum, source identity, and NAS synchronization contract to pass before releasing the corresponding rest-wave Pipeline 5 gates.

### 5. Pipeline 2 hook-backed rest wave

After B2-BASSEZ is reviewed, launch `B2-HOOKED-REST` without Bassez:

```bash
bash src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  --datasets Gongsharma_cmv_young_males,Kfoury,Joanito,CombinedPBMC,Myocardial_infarction \
  --force
```

The submitter must retain the built-in `gongsharma_cap -> combinedpbmc` `afterok` dependency. Require:

- Gongsharma cap outputs and checksums.
- Kfoury prerequisite output and checksum.
- Joanito prerequisite output and inline semantic RDS/H5AD checks.
- Myocardial reconstructed counts with `layers['counts']`, matching shape, finite nonnegative integer values, and checksum.
- CombinedPBMC canonical `combined_pbmc.h5ad` with nonblank `Sample`, `cond`, and `batch`; reject `combined_pbmc_batch_effect_analysis_batch_effect_analysis_ECODAprocessed.h5ad` as a substitute.

If this gate fails, block only the hook-backed Stage 3/4/5 rows. The benchmark and batch no-hook Stage 3 waves may proceed after B3-BASSEZ if their declared inputs pass current Stage 3 preflight.

### 6. Pipeline 3 rolling waves after Bassez validation

Once B3-BASSEZ is reviewed, launch these no-hook waves in parallel; each excludes Bassez and uses `--force`:

**Benchmark no-hook Stage 3 wave (`B3-BENCHMARK-NOHOOK`):**

```bash
bash src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --datasets Adams,Kim,Lee,Pelka,Smillie,Stephenson,Wu,Zhang \
  --views benchmark_analysis --force
```

**Batch no-hook Stage 3 wave (`B3-BATCH-NOHOOK`):**

```bash
bash src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --datasets Stephenson,Alzheimer,Breast_cancer,Covid19_PBMC,Kidney_KPMP,Diabetes,Lupus_PBMC,Lung,Parkinson \
  --views batch_effect_uncorrected --force
```

After both B2-HOOKED-REST and B3-BASSEZ are reviewed, launch these hook-backed Stage 3 waves in parallel:

**Benchmark hook Stage 3 wave (`B3-BENCHMARK-HOOK`):**

```bash
bash src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --datasets Gongsharma_cmv_young_males,Kfoury \
  --views benchmark_analysis --force
```

**Batch hook Stage 3 wave (`B3-BATCH-HOOK`):**

```bash
bash src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  --datasets Joanito,CombinedPBMC,Myocardial_infarction \
  --views batch_effect_uncorrected --force
```

For every Stage 3 gate, require the declared input/output view names from `datasets.json`, compute-node H5AD preflight, `benchmark_h5ad_contract.py`, raw counts, required PCA/graph keys, checksums, selected sync, and complete run-owned status. Do not use `--exact-batch-selection` for these split rolling waves; that switch is reserved for the immutable twelve-row batch contract.

### 7. Pipeline 4 rolling waves

`B4-BASSEZ` may run as soon as B3-BASSEZ passes, while the rest of Stage 3 executes. After B4-BASSEZ and all corresponding benchmark Stage 3 gates are reviewed, launch `B4-BENCHMARK-REST`:

```bash
bash src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  --datasets Adams,Gongsharma_cmv_young_males,Kfoury,Kim,Lee,Pelka,Smillie,Stephenson,Wu,Zhang \
  --views benchmark_analysis --force
```

After B4-BASSEZ and all corresponding batch Stage 3 gates are reviewed, launch `B4-BATCH-REST`:

```bash
bash src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  --datasets Joanito,Stephenson,CombinedPBMC,Alzheimer,Breast_cancer,Covid19_PBMC,Kidney_KPMP,Myocardial_infarction,Diabetes,Lupus_PBMC,Lung,Parkinson \
  --views batch_effect_uncorrected --force
```

The batch Stage 4 gate must record the three explicit auto-annotation exemptions. `Alzheimer`, `Diabetes`, and `Parkinson` are unsuitable for both HiTME and scATOMIC and must be cleanly omitted from the runnable annotation selection rather than assigned fabricated annotations. Their batch-effect Stage 3 and Stage 5 rows remain in scope. Require the Stage 4 `runnable_selection`, preparation/chunk/annotation/merge watchdog chain, valid annotation coverage for runnable rows, checksums, per-view sample identity, and selected sync.

### 8. Pipeline 5 rolling waves

After B4-BASSEZ is reviewed, launch B5-BASSEZ. After B5-BASSEZ and `B4-BENCHMARK-REST` are reviewed, launch `B5-BENCHMARK-REST`:

```bash
bash src/5_run_benchmark_methods/1_submit_hpc_array.sh \
  --datasets Adams,Gongsharma_cmv_young_males,Kfoury,Kim,Lee,Pelka,Smillie,Stephenson,Wu,Zhang \
  --methods gloscope,mofa,pseudobulk,composition,scitd,mrvi,scpoli,pilot,qot,pilotgm \
  --analyses trans,zeroimp --force
```

After `B4-BATCH-REST` is reviewed, launch `B5-BATCH-REST` independently of
`B5-BASSEZ`; Bassez is not in the batch dataset universe and the batch suite
uses the separate flexible any-GPU resource class:
```bash
bash src/5_run_benchmark_methods/1_submit_hpc_array.sh \
  --datasets Joanito,Stephenson,CombinedPBMC,Alzheimer,Breast_cancer,Covid19_PBMC,Kidney_KPMP,Myocardial_infarction,Diabetes,Lupus_PBMC,Lung,Parkinson \
  --pass uncorrected \
  --methods prepare_pseudobulk,pseudobulk,gloscope,composition,mrvi,pilot,qot \
  --force
```

The ordinary benchmark gates use `benchmark_analysis` only. The batch gate uses `batch_effect_uncorrected` only and the fixed pass suite. Never pass Bassez, Lee, Smillie, or Zhang to the batch `--pass` command.

Stage 5 resource behavior and parameter screening remain explicit:

- Ordinary GloScope, MrVI, scPoli, and PILOT-GM-VAE use fixed parameter
  shards only for the heavy methods; small methods remain one dataset task to
  avoid repeated H5AD loading.
- The ordinary PILOT-GM-VAE suite is default-only:
  `hvg2000_highres`; no pilotgm parameter screening is scheduled.
- Ordinary default MrVI (`hvg2000`) and default scPoli
  (`hvg2000_highres_dims15`) use the H200-pinned class.
- Ordinary non-default MrVI shards (`hvg1000`, `hvg3000`) use CPU; ordinary
  non-default scPoli shards use the relaxed any-GPU class.
- Batch-effect PILOT-GM-VAE is not scheduled. Batch MrVI uses the flexible
  any-GPU class automatically through `--pass uncorrected`.
- Preserve input-size and CUDA peak-memory telemetry; do not silently move
  standard/default benchmark results off H200.

### 9. Promotion, failure, and completion rules

1. Bassez is the promotion barrier for the ordinary benchmark pipeline. The
   batch Pipeline 5 rest wave does not wait for ordinary B5-BASSEZ because
   Bassez is not a batch dataset; it depends on reviewed B4-BATCH-REST and its
   own batch Stage 3/resource contracts.
2. Stage-specific rest gates depend on their own upstream gate plus the
   applicable Bassez control gate. No-hook datasets do not wait for a
   nonexistent Pipeline 2 hook; hook-backed datasets wait for B2-HOOKED-REST.
3. A failure blocks only its dependent branch unless it is a Bassez control failure. Sibling gates already launched may finish and must retain their evidence; do not rerun successful siblings.
4. On a failed wave, create a new durable repair gate with a new run ID. Do not manually remove stale artifacts, partially synchronize outputs, or reuse a failed gate manifest.
5. A method result is reusable only after current checksum, schema, source identity, and ordered sample-universe validation. For every method except scITD, sample IDs and first-appearance order must exactly match the selected H5AD. Apply the documented scITD dropped-sample exception only where the scITD validator reports it.
6. Complete the rollout only after every Bassez control gate and every selected benchmark/batch rest gate has a passing first audit and Luna Max reviewer approval, with all selected NAS checksums and run-owned terminal statuses present. Record Pipeline 2 N/A rows explicitly and do not label them as passed.

## Critical files & anchors

- `datasets.json` — dataset flags, exact declared views, output names, subset filters, and the three auto-annotation exemptions.
- `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — fixed Pipeline 2 hook mapping, checksum/semantic reuse, force invalidation, and the CombinedPBMC dependency.
- `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh` — explicit two-column selection files, H5AD preflight, force/reuse, watchdog, and selected sync.
- `src/4_cell_type_annotation/1_submit_onboarding_stage.sh` — declared-view selection, auto-annotation exclusions, preparation/annotation/merge dependencies, and sync-only validation.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh` — ordinary versus batch namespaces, fixed uncorrected suite, GPU policy, per-method arrays/watchdogs, aggregate gate, and result validation.

## Verification

### Before the first gate

- Confirm Bamboo checkout and runtime manifest source identity match.
- Confirm the current reviewed image/manifest pair is readable and checksum-valid.
- Confirm no active environment mutation lock, competing durable gate, or unresolved owner exists for the selected rows.
- Verify the rollout lists against `datasets.json`; reject any benchmark-only dataset assigned a batch view.

### Per-gate checks

- **Pipeline 2:** validate every selected hook output with its checksum and `ecoda_validate_stage2_output`; run the myocardial, CombinedPBMC, and Joanito semantic contracts; verify `gongsharma_cap` precedes `combinedpbmc`.
- **Pipeline 3:** validate every selected declared H5AD view with `benchmark_h5ad_contract.py`, required counts/embeddings, nonempty schema, checksum, source identity, run-owned status, and selected NAS sync.
- **Pipeline 4:** verify `runnable_selection`, chunk completeness, annotation key coverage, sample identity, checksum, merge status, and exact exemption of Alzheimer/Diabetes/Parkinson from automated annotation.
- **Pipeline 5:** validate every selected Feather/RDS result, execution log, ordered sample IDs, source identity, matrix/RDS aggregate contracts, remote checksums, and final sync. Confirm the batch output namespace is `batch_effect/uncorrected`, never ordinary `benchmark/`.
- **Rolling promotion:** before each dependent gate launch, verify the predecessor's saved reviewer-approved release evidence and confirm the new selection manifest has no duplicate Bassez row and no overlap with a running sibling gate.
- **Final release:** all wave manifests must be current-source, terminal `COMPLETED`, audit-passed, reviewer-approved, and synchronized. Old `PREPARED`/`audit_state=NOT_STARTED` gate manifests remain historical and cannot satisfy this barrier.

## Assumptions & contingencies

- The rollout intentionally excludes `_debug` from production selections; retain the reviewed `_debug` gate as runtime evidence only.
- **Superseding assumption (2026-09-04):** Corrected batch-effect views are
  deferred until the uncorrected branch is terminal, inspected, reviewed, and
  complete; they are not omitted. The targeted Joanito Pipeline 3 refresh
  already covers both `batch_effect_uncorrected` and
  `batch_effect_corrected`. Corrected-mode Pipeline 5 remains pending after
  the uncorrected recovery and requires the confirmed `columns.batch` contract
  before submission.
- Current evidence justifies forcing all hook-backed Pipeline 2 rows and all selected downstream rows. If a future current audited run proves a row complete before its wave launches, use that run's reviewed `sync-only` path instead of inventing a second run; otherwise retain the forced rebuild rule.
- If Bamboo `HEAD` changes after plan approval or no longer matches the runtime manifest, stop before launching and rebuild/review the immutable image. Do not run the gate with a dirty or mismatched source identity.
- If the H200 queue delays a standard/default benchmark row, keep it H200-pinned. Use flexible GPUs only for batch-effect or explicitly non-default GPU work; do not change the scientific comparability class to improve queue time.

## Execution checkpoint — rolling repairs after source fixes

Checkpoint updated after the interrupted batch Stage 3 gate, the Stage 3
memory/raw-variable/raw-loader repairs, three failed hooked Stage 2 attempts,
the reviewed hooked Stage 2 success, the failed full batch retry, the strict
accounting failure of the raw-loader retry, the lost completion transport of
the 500G retry, its later terminal recovery, the reviewed nine-row 500G batch
success, and the reviewed B4-BASSEZ annotation success. This section is
operational evidence, not a replacement for the promotion rules above.

### Source, gate implementation, and runtime

- Bamboo was fast-forwarded from `e147c248d948449486f72321e5354c7deec36697`
  through `f5f5b3c`, `f876db3`, `10ed7ac`, `3f16171`, `fc52240`, `ad1e98c`,
  `eaa4030`, and `6fcf1fb`. The latest source includes Myocardial repair,
  Stage 3 raw-count ownership, bounded sparse-PCA allocation, named/counts-only
  H5AD contracts, 64G Stage 2 watchdog memory, raw-variable expansion, reduced
  preprocessing copy peaks, and a backed raw-only H5AD loader.
- Commit `46f6e3f` adds the cross-pipeline preflight optimization: one strict
  digest/size read is reused for sidecar and scratch/NAS checks, duplicate
  local sync hashes are removed, repeated Stage 4 union hashes are cached per
  run, and annotation contracts expose an explicit caller-validated-sidecar
  mode that skips only a duplicate sidecar hash after `ecoda_validate_checksum`.
  The Stage 5 source-repair loop now validates each scratch/NAS H5AD once and
  records progress; no invalid-sidecar overwrite or status-only success path
  is allowed.
- The durable gate implementation now treats `dependency_manifests` as
  reviewed predecessor lineage independent of `serialization_group`.
  `serialization_group` remains the explicit deterministic mutex: an active
  `RUNNING` sibling with the same project/profile/profile digest/group is
  rejected as `resource_lock_conflict` and moved to `PRELAUNCH_STOP`; distinct
  groups sharing a reviewed predecessor may launch concurrently. The focused
  regression passed for B4 and benchmark-P3 distinct groups, malformed and
  mismatched predecessors, and same-group active conflict.
- The gate contract documentation was updated in the global durable-gate
  skill, its lifecycle/schema references, and the repository ECODA skill.
  The repository regression is committed as `af4e66f` and pushed. The global
  implementation and global skill references remain outside the repository
  checkout and are the local harness runtime.
- The pre-pull full Bamboo `git status --porcelain` checks before the 500G
  batch and B4 launches reported only preserved untracked `logs/`; a later
  check also found the copied audit helper under remote `.gate`, which was
  relocated to the scratch gate evidence root before B4 launch. `HEAD` was
  `6fcf1fb`, with no source/config changes, Git lock, environment lock, or
  lockfile mutation.
- The source-matched runtime gate
  `ecoda_runtime_build_bassez_rolling_6fcf1fb_20260831T061240Z` completed with
  scheduler build job `4368128`, passed its first audit, and received Luna Max
  approval. Its image SHA-256 is
  `37b173419b7a31bfa1b873ff6268da33eb32fc8d5801c537f44dfa2abdaed9a0`;
  its runtime manifest SHA-256 is
  `b3e1e861be20c8c5b18d99f72ac0cac743cb82100a751bc32c0560c31f22d9b3`.
- The reviewed runtime manifest binds Git revision
  `6fcf1fb1d4e71759bb7fa6e4e826d882384ae011`, path-preserving layout,
  `py-cuda13`, and the immutable `pixi.lock` SHA-256.
- The Stage 4 worker-timeout hardening commit `8471bd5` raises the
  annotation/merge worker submission and retry limit to the explicit
  `ANNOTATION_WORKER_TIME_LIMIT` default of `12:00:00`, changes both worker
  directives to `12:00:00`, and classifies terminal non-OOM annotation states
  (including `TIMEOUT`) as failures. Bamboo was synchronized to this commit.
- The first source-matched `8471bd5` runtime rebuild gate
  `ecoda_runtime_build_bassez_rolling_8471bd5_20260902T054254Z` failed before
  image mutation because `build_ecoda_runtime.sh` correctly rejected active
  Slurm jobs (`4371377`, `4371379`, and the build job `4371743`). Its first
  inspect covered `4371743` exactly once and recorded failed accounting; the
  gate is immutable failed evidence and must not be reused. A fresh runtime
  gate is required after the active jobs drain.
- The optimized runtime gate
  `ecoda_runtime_build_bassez_rolling_46f6e3f_20260902T194508Z` completed
  with build job `4374019`; its first inspect and artifact audit passed, and
  Luna Max approved it at `2026-09-02T19:57:55Z`. It is release-eligible for
  commit `46f6e3f` only. The subsequent Ensembl stable-ID correction is
  intentionally excluded from this image and requires a second runtime build.

### B5 preflight failure and cross-pipeline checksum optimization

- `ecoda_bassez_rolling_b5_benchmark_rest_80f71c8_20260902T081614Z` reached
  terminal `FAILED` at `2026-09-02T18:44:54Z` after source identity creation
  and H5AD compute preflight array `4373935`. Its ten task statuses were
  `STATE=OK` with task IDs 1–10, and the one required first inspect covered
  array parent/tasks exactly once (`COMPLETED|0:0`); the wrapper still failed
  with exit code 1 and no benchmark arrays were submitted. The failed gate is
  preserved and is not release-eligible.
- The local preflight optimization keeps failure closed: a nonzero
  `sbatch --wait` result is logged and rejected; per-task statuses now carry
  the run ID, dataset, view, and task ID and must match the run-owned
  preflight manifest. No status-only success fallback is introduced.
- Shared checksum primitives now expose the exact digest/size from a strict
  validation and provide a no-write sidecar-record confirmation. Pipeline 2
  watchdog output, Pipeline 3/4 NAS comparisons, Pipeline 4 merge markers,
  and Pipeline 5 selected-result manifests reuse those records instead of
  rereading payloads.
- Pipeline 5 source repair validates each scratch/NAS H5AD path once, reuses
  the validated digest/size for scratch/NAS identity, and uses the already
  validated sidecars when building source identity. Existing invalid sidecars
  remain fatal; only genuinely missing sidecars may be created after the H5AD
  semantic contract passes.
- Focused local verification passed for checksum reuse, H5AD run/task binding,
  source identity, Stage 2/3/4/5 submitters/watchdogs, benchmark sync,
  annotation merge safety, Ensembl stable-ID mapping, and the explicit
  sidecar-validation path. Commit `46f6e3f` carries the preflight speedup;
  commit `0145b24` adds the source-matched Ensembl correction and was built
  into a separately reviewed runtime before affected Stage 3/4 work.
- `ecoda_bassez_rolling_b4_batch_rest_80f71c8_20260902T081547Z` reached
  terminal `FAILED` at `2026-09-02T19:33:55Z` after preparation and annotation
  succeeded but merge array/watchdog `4373971`/`4373972` both returned
  `FAILED|1:0`. Its one first inspect covered SCGATE `4372893` from wrapper
  evidence, arrays `4372896`, `4372917`, `4373971`, and watchdogs `4372897`,
  `4372918`, `4373972` exactly once; artifact/immutable audits passed where
  applicable, but accounting failed and no review was performed. Preserve the
  gate as failed evidence; no repair launch has been made.
- The B4 failed run's merge tasks 4, 6, and 9
  (`Breast_cancer`, `Kidney_KPMP`, and `Lung`) failed the canonical annotation
  contract because the merged `layer1` dataset anchor was entirely blank.
  This is a real dual-annotation semantic failure, not a checksum or scheduler
  issue; do not loosen anchor requirements or sidecar validation. A fresh
  reviewed B4 repair gate for the affected batch rows is required before
  releasing B5 batch-rest.
- Inspection of the three failed batch inputs found Ensembl stable IDs in
  `var_names` (for example `ENSG00000278232`) while the annotation signatures
  require gene symbols. The canonical `gene_utils` map previously omitted its
  `Gene stable ID` column, leaving those identifiers unmapped. The repair adds
  stable-ID and version-suffix normalization before existing symbol/alias
  mappings; this is an upstream preprocessing correction, not an annotation
  anchor relaxation. A fresh Stage 3 rebuild for the affected rows is required
  before Stage 4 repair.
- The optimized B5 retry
  `ecoda_bassez_rolling_b5_benchmark_rest_0145b24_20260902T201538Z` reduced
  source-check logging to one pass per selected scratch/NAS H5AD, but reached
  terminal `FAILED` at `2026-09-02T20:26:29Z`: all ten H5AD preflight array
  tasks (`4374121` parent, tasks `4374122`–`4374130`) audited
  `COMPLETED|0:0`, while the submitter checked status files before the shared
  filesystem exposed `Adams__benchmark_analysis.status`. The gate's one first
  inspect passed accounting and generic audits, but the wrapper failure is
  preserved; a bounded local status-settle grace is required and does not
  mask nonzero `sbatch --wait` results.
- The follow-up status-settle patch waits a bounded
  `H5AD_PREFLIGHT_STATUS_GRACE_SECONDS` (default 60 seconds) for run-owned
  status files after a successful `sbatch --wait`. It only waits on local
  filesystem publication, still rejects a nonzero scheduler return, and
  preserves run/task-bound status validation. This prevents shared-filesystem
  publication races without converting status files into scheduler evidence.
- The same local status-settle barrier is now applied to Pipeline 3 before
  consuming compute-node H5AD statuses; the pending source correction is
  intentionally kept separate from the already reviewed `de621c6` runtime and
  requires one final source-matched runtime rebuild before any dependent gate.
- Fresh Stage 3 gene-fix gate
  `ecoda_bassez_rolling_b3_batch_gene_fix_0145b24_20260902T201359Z` rebuilt
  `Breast_cancer`, `Kidney_KPMP`, and `Lung` with stable-ID-to-symbol mapping.
- Array `4374117` and watchdog `4374118` completed, the fixed three-row H5AD
  audit passed with strict checksums and reported residual Ensembl IDs that
  have no `Gene name` in the reference table; downstream annotation anchors
  remain the required success criterion. Luna Max approved the gate at
  `2026-09-02T20:48:23Z`.
- The third runtime
  `ecoda_runtime_build_bassez_rolling_0145b24_20260902T195959Z` completed
  with build job `4374113`, passed its first audit, and received Luna Max
  approval at `2026-09-02T20:12:25Z`; its image is
  `ecoda-py-cuda13-path-preserving-0145b24.sif` with SHA-256
  `96e64d6e94f944c47ca9f9ec17b7fcb6b05ccca9bbd70ea9d35268404fe6da6c`.

### B4/B5 final-gate failures and parameter sharding

- The source-matched final runtime
  `ecoda_runtime_build_bassez_rolling_01c8d93_20260902T211859Z` completed with
  build job `4374197`, passed its first audit, and received Luna Max approval
  at `2026-09-02T21:31:24Z`. Its immutable source is
  `01c8d93928629dddb5ae388b0e8f02c77394d0aa`.
- `ecoda_bassez_rolling_b4_batch_repair_01c8d93_20260902T213629Z` reached
  terminal `FAILED` at `2026-09-03T00:00:41Z`. Preparation, annotation, and
  merge arrays completed, but merge watchdog `4375215` was `OUT_OF_MEMORY`
  (`0:125`) while running `annotation_contract.py --h5ad` with a hard-coded
  2G allocation. Its first inspect covered every recorded array/watchdog ID
  exactly once and found only that OOM; no reviewer approval was performed.
  The checked-in repair raises the merge-watchdog allocation to 32G.
- `ecoda_bassez_rolling_b5_benchmark_rest_01c8d93_20260902T213757Z` reached
  terminal `FAILED` at `2026-09-03T09:52:23Z`. GloScope task 1 (Adams)
  rejected one-sided sample-name normalization; scPoli watchdog `4374284`
  and PILOT-GM-VAE watchdog `4374290` hit the shared 12-hour time limit, and
  scPoli array `4374283` was still `RUNNING` at first inspection. The first
  inspect covered two preflight arrays, thirteen method arrays, and thirteen
  watchdogs exactly once; no reviewer approval was performed.
- The prior matrix watchdogs were submitted without an `afterany` dependency,
  so their 12-hour wall clock included queued-array time. The submitter now
  binds each watchdog to its own array with `--dependency=afterany:<array>`,
  preserving the full watchdog window for terminal task accounting.
- Remote execution-time evidence shows the costly methods are GloScope,
  MrVI, scPoli, and PILOT-GM-VAE; short PILOT, QOT, pseudobulk, and composition
  methods remain one dataset task to avoid repeated H5AD loading.
- The Stage 5 submitter now emits fixed four-column parameter-shard manifests
  for those heavy ordinary methods. Default MrVI (`hvg2000`) and default
  scPoli (`hvg2000_highres_dims15`) use H200; non-default MrVI uses CPU and
  non-default scPoli uses any-GPU. GloScope shards consolidate into the
  canonical method RDS only after their watchdog gate.
- Ordinary PILOT-GM-VAE is default-only (`hvg2000_highres`); it is not
  screened across HVG/resolution combinations. PILOT-GM-VAE is excluded from
  the batch-effect method suite and batch candidate evidence.
- Focused verification passed for parameter-shard submission/dependencies,
  four-column OOM retry preservation, GloScope consolidation, execution-log
  shard merging, default-only PILOT-GM-VAE artifacts, batch exclusion, and
  MrVI CPU/default-H200 policy. A new runtime build and fresh durable B4/B5
  gates are required for release.

### Historical failed batch Stage 3 gates

- The original c9 batch no-hook gate
  `ecoda_bassez_rolling_b3_batch_nohook_c9f9398_20260830T194232Z` was
  reconciled after its runner exited. Its first inspect covered array
  `4367815`, watchdog `4367816`, retry arrays `4367825` and `4367835` exactly
  once; the gate remains immutable `PRELAUNCH_STOP`/failed evidence with no
  reviewer approval.
- The first post-memory-fix nine-row gate
  `ecoda_bassez_rolling_b3_batch_nohook_3f16171_20260831T014424Z` failed at
  task 7 (`Lupus_PBMC`) with array `4368041` and watchdog `4368042`. Its first
  inspect covered both IDs exactly once. The worker rejected a legitimate
  raw-variable expansion (`raw.X` 32,738 genes versus current X 1,999).
- The `fc52240` nine-row gate
  `ecoda_bassez_rolling_b3_batch_nohook_fc52240_20260831T024957Z` reached
  terminal `FAILED` at `2026-08-31T05:55:05Z` after initial array `4368061`,
  watchdog `4368062`, and retry arrays `4368071`/`4368115`; the final retry
  hit the 500G OOM ceiling for Alzheimer. Its first inspect covered all four
  IDs exactly once. Failure evidence and note remain preserved.
- The raw-loader nine-row gate
  `ecoda_bassez_rolling_b3_batch_nohook_6fcf1fb_20260831T062617Z` reached
  wrapper `STATE=OK`, but its first inspect failed strict accounting because
  initial array `4368129` and retry `4368203` were OOM; final retry `4368207`
  and watchdog `4368130` completed. Failure evidence and note remain
  preserved; output was not promoted.
- The 500G initial-ceiling gate
  `ecoda_bassez_rolling_b3_batch_nohook_6fcf1fb_20260831T130325Z` initially
  entered immutable `PRELAUNCH_STOP` after transport loss. Recovery later
  confirmed no runner/tmux session and remote `COMPLETED|0` at
  `2026-08-31T22:02:13Z`; its one first inspect covered array `4368609`,
  all accounting rows and profile audit passed, but state remained
  `PRELAUNCH_STOP`, so it was not release-eligible and was not reused.

### Stage 3 batch no-hook success

- Fresh gate
  `ecoda_bassez_rolling_b3_batch_nohook_6fcf1fb_20260901T095035Z` rebuilt all
  nine rows with the reviewed runtime, `--force`, `--mem 500G`,
  `--max-mem 500G`, and `--throttle 1`. It completed at
  `2026-09-01T19:55:35Z` with array `4369394` and watchdog `4369395`; no OOM
  retry was emitted.
- Its one first inspect covered both IDs exactly once; accounting, terminal
  profile checks, immutable fingerprints, and artifact-contract checks passed.
  The independent read-only Stage 3 batch audit passed all nine declared
  `batch_effect_uncorrected` H5AD outputs, required `Sample`/declared-label
  metadata, strict scratch/NAS MD5+size+PATH sidecars, and byte identity.
  Alzheimer completed with 1,395,601 cells and an identical
  191,350,083,466-byte scratch/NAS artifact.
- Luna Max reviewer approval passed and the gate is release-eligible. This
  reviewed gate is the current batch Stage 3 predecessor; no batch output is
  skipped.

### Stage 3 source repairs

- `base_preprocessing()` adopts validated integer raw matrices by reference
  when dimensions match and replaces the local AnnData container by reference
  when raw variables expand or shrink. It validates raw observation order and
  raw variable metadata, releases obsolete normalized/raw ownership before
  replacement, and preserves the required filtered counts vault plus
  normalized/log X.
- `load_single_input()` opens H5AD inputs backed, validates a bounded raw
  matrix sample and raw obs/var alignment, materializes only the authoritative
  integer raw matrix plus required metadata, closes the backed handle, and
  falls back to the established eager path for counts-layer or non-integer
  inputs. This targets the residual Alzheimer peak without changing the
  counts-layer contract.
- `compute_pca_and_store()` constructs a lightweight selected-X/obs/var
  object without parent layers. Sparse centered/clipped scaling uses bounded
  CSR arithmetic and implicit-zero baselines; dense behavior remains unchanged.
- Generic H5AD artifact validation accepts valid top-level-X and established
  counts-only `layers["counts"]` layouts while failing closed on missing
  matrix storage or malformed index metadata.
- Focused verification passed:
  `tests/test_preprocess_h5ad_loading.py`,
  `tests/test_preprocessing_raw_counts.py`,
  `tests/test_preprocessing_sample_filter.py`,
  `tests/test_preprocessing_h5ad_atomic.py`,
  `tests/test_durable_hpc_gate_parallelism.py`,
  `tests/test_durable_profile_stage_neutral.sh`, Python compilation, and the
  prior artifact/Stage 2 watchdog/submitter contracts. The B4 audit helper
  additionally passed its partial-coverage H5AD smoke and production audit.
  Source commits are `9b3f790`, `f5f5b3c`, `f876db3`, `10ed7ac`, `3f16171`,
  `fc52240`, `ad1e98c`, `eaa4030`, `6fcf1fb`, `af4e66f`, and `8471bd5`.

### Bassez control lane

- `B2-BASSEZ`
  (`ecoda_bassez_rolling_b2_bassez_c9f9398_20260830T194232Z`) remains
  `COMPLETED`, first-audit passed, Bassez RDS checksum/size/path validation
  passed, pinned-R semantic validation passed (`226454` metadata rows), and
  Luna Max approval passed.
- `B3-BASSEZ`
  (`ecoda_bassez_rolling_b3_bassez_c9f9398_r2_20260830T194232Z`) remains
  `COMPLETED`, first-audit passed, the `benchmark_analysis` H5AD contract and
  declared `Sample`/`expansion` columns passed, scratch/NAS checksums matched,
  and Luna Max approval passed.
- Fresh `B4-BASSEZ`
  (`ecoda_bassez_rolling_b4_bassez_6fcf1fb_20260901T205324Z`) selected only
  Bassez/`benchmark_analysis`, ran 13 contiguous annotation chunks, and
  completed at `2026-09-01T21:16:31Z`. Its first inspect covered SCGATE
  `4371168`, preparation array/watchdog `4371169`/`4371170`, annotation
  array/watchdog `4371171`/`4371184`, and merge array/watchdog
  `4371194`/`4371195` exactly once; all accounting and profile checks passed.
- The independent B4 artifact audit passed run-owned manifests/status/owner/
  merge/sync contracts, required dual-method schema, exact keys, checksums,
  and scratch/NAS identity. HiTME layer1/2/3 had 60,336 nonblank rows and
  14,612 intentional NA rows; scATOMIC_pred had 74,948 nonblank rows. This
  matches the canonical partial-coverage contract and is not a failure.
  Luna Max approval passed; B4 is release-eligible.
- The first B4 benchmark-rest gate
  `ecoda_bassez_rolling_b4_benchmark_rest_8723d91_20260902T025830Z` failed
  after preparation and annotation array `4371403`; annotation watchdog
  `4371393` retried OOM rows as `4371494`, then failed validation because
  Adams `annotations_chunk_34.feather` was missing after annotation array
  task 27 reached terminal `TIMEOUT` (`JobIDRaw=4371430`). Its first inspect
  covered wrapper IDs `4371391`, `4371392`, `4371393`, and `4371403` exactly
  once; accounting for those IDs and the profile audit passed, but the
  wrapper state is failed and no review was performed. The failure exposed
  the two-hour Stage 4 worker limit; the repair is committed as `8471bd5`.
- `B5-BASSEZ` is now eligible after reviewed B4. It must use a fresh gate and
  distinct serialization group from the independent Stage 3 benchmark/batch
  waves. It remains the final Bassez control barrier before corresponding
  Pipeline 5 rest waves.

### Pipeline 2 hooked rest

- Earlier c9/f5/f876/10ed hooked Stage 2 attempts are failed historical
  evidence. Their emitted scheduler IDs were first-inspected exactly once,
  and none received reviewer approval.
- The memory-sized gate
  `ecoda_bassez_rolling_b2_hooked_rest_3f16171_20260831T012450Z` ran step jobs
  `4368034`–`4368038` and watchdog `4368039` with
  `STAGE2_WATCHDOG_MEM=64G`. Its watchdog and independent artifact audit
  validated Gongsharma, CombinedPBMC, Joanito, Kfoury, and Myocardial; the
  first audit and Luna Max review passed.
- Gongsharma, Kfoury, Joanito, CombinedPBMC, and Myocardial are released from
  the Pipeline 2 prerequisite. No current-source Stage 2 rerun is needed
  solely because the later Stage 3 code changed; its source-level contract
  remains validated by the reviewed 3f gate.

### Ordinary benchmark rows

- The benchmark no-hook Stage 3 gate
  `ecoda_bassez_rolling_b3_benchmark_nohook_c9f9398_20260830T194232Z` remains
  `COMPLETED` and Luna Max approved for Adams, Kim, Lee, Pelka, Smillie,
  Stephenson, Wu, and Zhang.
- The benchmark hook Stage 3 rows (Gongsharma, Kfoury) are now eligible from
  reviewed B2-HOOKED-REST and B3-BASSEZ. They may launch in a distinct
  serialization group while B5-BASSEZ and the batch-hook branch run, because
  lineage no longer requires group equality.
- Benchmark Pipeline 4/5 rest waves remain pending until their corresponding
  Stage 3 hook/no-hook and Bassez control predecessors pass.

### Current parallel wave

- B4-BASSEZ is reviewed `COMPLETED`; the previously active batch no-hook gate
  is also reviewed `COMPLETED`. The next eligible independent gates are:
  1. B3 benchmark hook Stage 3 for Gongsharma/Kfoury, with dependencies
     B2-HOOKED-REST and B3-BASSEZ, distinct group
     `ecoda-bassez-rolling-b3-benchmark-hook`.
  2. B3 batch hook Stage 3 for Joanito/CombinedPBMC/Myocardial, with the same
     reviewed Pipeline 2/B3-BASSEZ lineage and distinct group
     `ecoda-bassez-rolling-b3-batch-hook`.
  3. B5-BASSEZ after reviewed B4, with distinct group
     `ecoda-bassez-rolling-b5-bassez`.
- These three gates may be prepared and launched in parallel after clean
  Bamboo status checks. They have disjoint output namespaces/resources and
  share only reviewed predecessor lineage. Each still gets its own durable
  waiter, first inspect, gate-specific artifact audit, and Luna Max review.
- B4 benchmark-rest and B4 batch-rest are not yet launchable: each also needs
  its corresponding Stage 3 hook gate. Their Stage 4 scGate database/model/
  ontology cache is a true shared resource; keep both on the same explicit
  Stage 4 mutex group and run them one at a time unless the cache is made
  safely immutable. They remain independent of B5-BASSEZ.

### Ordered next operations

1. Prepare/reconcile/launch B3 benchmark hook, B3 batch hook, and B5-BASSEZ
   as independent fresh gates in parallel, using reviewed lineage paths and
   distinct serialization groups. Before each launch, run full Bamboo
   `git status --porcelain`; only preserved `?? logs/` is allowed.
2. After each terminal wave, inspect every emitted scheduler ID exactly once,
   run its independent artifact/checksum/NAS audit, and obtain Luna Max review.
   A failed branch is preserved and repaired independently; do not block or
   rerun a successful sibling.
3. Once B3 benchmark hook passes, launch B4 benchmark-rest when B4-BASSEZ and
   B3 benchmark no-hook are already reviewed. Once B3 batch hook passes, launch
   B4 batch-rest when B4-BASSEZ and B3 batch no-hook are reviewed. Serialize
   those two Stage 4 gates on the shared scGate mutex only.
4. After reviewed B5-BASSEZ and each corresponding B4 rest gate, launch the
   benchmark-rest and batch-rest Pipeline 5 gates independently; B5-BASSEZ
   may run concurrently with either Pipeline 4 branch.
5. Before every future launch, verify the current source/runtime identity,
   immutable dataset/Pixi fingerprints, no competing same-resource gate, and
   no environment or Git lock. Preserve all old failed, PRELAUNCH_STOP, and
   stale PREPARED manifests as discrepancy evidence.
- **Current source/runtime lineage:** commit `34a315f` introduced selective
  Pipeline 5 parameter sharding, exact GloScope sample IDs, default-only
  ordinary PILOT-GM-VAE, and batch exclusion of PILOT-GM-VAE. Commit `c5da3c2`
  fixed GloScope consolidation when R resolves the configured scratch symlink
  differently from checksum sidecars; its symlink regression test passed.
  Commit `d558b26` added the mandatory targeted-recovery rule to `AGENTS.md`.
  The source-matched `d558b26` path-preserving runtime was built by
  `ecoda_runtime_build_bassez_rolling_d558b26_20260903174000Z`, first-audited,
  and Luna Max approved. `datasets.json`, `pixi.toml`, and `pixi.lock` were
  unchanged.

- **Stage 4 recovery completed:** The old B4 batch repair run
  `stage4_20260902233731_3919808` had successful preparation, annotation, and
  merge arrays; only its 2G merge watchdog OOMed. The reviewed 32G
  merge-only recovery passed for array `4375214` and watchdog `4375917`.
  A guarded owner recovery and `--sync-only` gate then completed and was
  reviewed. It transferred the nine eligible batch H5ADs to NAS at roughly
  18--26 GB each and left the three automated-annotation exemptions untouched.

- **Ordinary Pipeline 5 recovery completed:** The current-source ordinary
  reconciliation skipped every already-valid row except
  `Gongsharma_cmv_young_males/benchmark_analysis/scitd`; its targeted
  no-force recovery also skipped the valid artifact and was reviewed. Adams
  GloScope's five parameter shards were consolidated without recomputation
  after the sidecar-path fix and its RDS contract passed. The initial ordinary
  gate's first inspect correctly preserved the failed OOM attempt
  `4376109`; retry `4376194` completed. No broad ordinary rerun was used.

- **Joanito stale-cache finding:** Pipeline 2's
  `1.3.1_prepare_joanito.R` derives `cell.type_new` and
  `ecoda_validate_stage2_output` requires it. Stage 2 Joanito owner
  `stage2_20260831032828_3836504` is `OK`. Pipeline 3's
  `preprocess_utils.py` creates `JoaI_..._raw.h5ad` only when absent; the
  cached raw H5AD was timestamped `2026-08-10`, before the derivation
  hardening commit `087eee1` (`2026-08-28`). Pipeline 3 therefore propagated a
  stale cache. The current uncorrected H5AD lacks `cell.type_new`, and the
  corrected Joanito H5AD is absent. This is a Pipeline 3 cache-reuse defect,
  not a missing Pipeline 2 derivation.

- **Batch uncorrected attempts and exact failures:** The first forced batch
  gate failed before method submission because the Alzheimer and Parkinson
  preflight status files exceeded the 60-second publication grace; their
  H5ADs are approximately 179 GB and 150 GB. The retry used a 1800-second
  grace and submitted the fixed seven-method suite without PILOT-GM-VAE or
  scPoli. Its run root was terminal `FAIL` at
  `2026-09-04T05:58:33Z`. The final recovery accounting found:
  - PILOT task 10 (`Joanito`) failed non-OOM because `cell.type_new` was
    missing;
  - QOT task 10 failed for the same reason;
  - prepare-pseudobulk task 1 (`Alzheimer`) reached the 500G OOM ceiling;
  - pseudobulk/composition did not obtain terminal matrix rows after that
    dependency failure;
  - MrVI completed successfully;
  - GloScope's remaining Alzheimer retry was pending for node availability
    and was canceled by the agent; completed GloScope artifacts and all
    manifests remain preserved.

- **Next steps before corrected mode:** Do not rerun Pipeline 2. After the
  current failed gate is preserved, perform a targeted Pipeline 3 refresh for
  Joanito's `batch_effect_uncorrected` and `batch_effect_corrected` views by
  invalidating only the stale raw H5AD cache and regenerating from the
  validated Pipeline 2 RDS. Revalidate/remerge existing Stage 4 annotation
  checkpoints only for Joanito if regenerated H5ADs require it. Then run
  targeted Pipeline 5 recovery for missing/invalid Joanito methods and
  Alzheimer GloScope, with Alzheimer launched at the maximum configured
  memory. Preserve all valid rows and use distinct durable serialization
  groups for independent Joanito and Alzheimer repairs. Run corrected-mode
  Pipeline 3--5 only after this uncorrected branch is reviewed and complete.

- **Evidence and safety state:** All failed, PRELAUNCH_STOP, and unlaunched
  PREPARED gate manifests remain preserved under `.gate/`; no historical gate
  is being reused as a promotion barrier. Focused local verification passed
  for matrix sharding/watchdog behavior, GloScope consolidation, benchmark
  synchronization, batch-method exclusion, and the symlink-safe consolidation
  regression. The next repair must be a new durable gate with an explicit
  failure scope; no canonical Stage 3 artifact is to be rewritten without
  preserving its prior evidence and recording the downstream revalidation
  dependency.
- **Joanito Pipeline 3 cache refresh completed:** The targeted gate
  `ecoda_bassez_rolling_b3_joanito_batch_refresh_d558b26_20260904070000Z`
  backed up the pre-refresh raw cache and existing uncorrected H5AD by
  hardlink, invalidated only the stale raw cache, and reran both
  `batch_effect_uncorrected` and `batch_effect_corrected` with the fixed
  source-matched d558b26 runtime. The wrapper emitted preprocessing array
  `4378942` and watchdog `4378943`; both completed. The single terminal
  inspect passed at `2026-09-04T09:40:49Z` with one accounting query and no
  discrepancies, and the Luna Max reviewer approved at `2026-09-04T09:42:11Z`.
  The gate is release-eligible; evidence is in the corresponding `.gate`
  wait, inspect, review, and manifest records.

- **Immediate downstream action:** Revalidate/remerge only Joanito's existing
  Stage 4 dual-annotation checkpoints into the two regenerated H5ADs, using a
  new run-owned Stage 4 merge gate and the prior validated annotation union
  and per-sample Feather checkpoints. Do not rerun annotation workers or
  touch any other dataset. After that gate is inspected and reviewer-approved,
  repair only the missing/invalid uncorrected Pipeline 5 rows: Joanito
  PILOT, QOT, and composition, plus both missing GloScope rows for Alzheimer
  and Parkinson. Preserve valid MrVI, GloScope, pseudobulk, and all other
  dataset/method artifacts. The Alzheimer prepare-pseudobulk 500G OOM remains
  a separate unresolved dependency and must not receive a blind retry.
- **Joanito Stage 4 remerge completed:** Gate
  `ecoda_bassez_rolling_b4_joanito_remerge_d558b26_20260904094500Z` launched
  at `2026-09-04T10:03:26Z` with command digest
  `29d5e0c70ca11bbc596d9dc3e902c083aed3021724aa6b9d537cd6aeabda27f1` and
  completed at `2026-09-04T10:54:26Z`. It created run
  `stage4_20260904094500_joanito_remerge`, built a fresh union with
  `758,928,057` nonzeros and 84 chunks covering 168 samples, reused only the
  84 prior validated annotation Feather checkpoints, merged both regenerated
  views without annotation workers, and guardedly synced both H5ADs to NAS.
  The single terminal inspect passed at `2026-09-04T10:55:24Z` with one
  accounting query covering arrays `4379393`/`4379428` and watchdogs
  `4379394`/`4379429`; all four were `COMPLETED|0:0`, all five artifact
  contracts and immutable fingerprints passed, and no discrepancies were
  recorded. Luna Max approved at `2026-09-04T10:57:05Z`; the gate is
  release-eligible.
- **Joanito uncorrected Pipeline 5 recovery completed, with a reconciled
  strict-validation expansion:** Gate
  `ecoda_bassez_rolling_b5_joanito_uncorrected_recovery_d558b26_20260904110000Z`
  completed at `2026-09-04T11:36:36Z` with command digest
  `eca798a68d8e26b1fb66408f6c7af834ee1c23cedf7f2ee8ac4e41e57b55592d`.
  The exact one-row selection remained Joanito/
  `batch_effect_uncorrected`/`batch_effect_uncorrected`, and no force was
  used. The canonical no-force submitter skipped only
  `prepare_pseudobulk` and `gloscope`; strict current-source validation
  classified `pseudobulk`, `composition`, `mrvi`, `pilot`, and `qot` as
  needing work, so those five rows were submitted and all completed. The
  resulting matrix/RDS contracts passed and the selected outputs were
  guardedly synced to NAS. The extra pseudobulk and MrVI recomputation is
  preserved as evidence, not hidden: post-remerge Joanito H5AD identity
  changed from the prior run's MD5/size
  `5663f60e47ab096caf22a6a39ea016c5`/`18514536708` to
  `65061f83948239642ea394db149eb1aa`/`18514914522`, and the prior MrVI row
  was no longer accepted by the strict current-source selection. No other
  dataset row was selected.

- The gate's single terminal inspect passed at `2026-09-04T11:37:48Z` with
  one accounting query covering preflight `4379518`, method arrays
  `4379521`, `4379523`, `4379525`, `4379527`, `4379529`, watchdogs
  `4379522`, `4379524`, `4379526`, `4379528`, `4379530`, and aggregate
  `4379531`; all matched `COMPLETED|0:0`, all five artifact contracts and
  immutable fingerprints passed, and no discrepancies were recorded. Luna
  Max approved the reconciled gate in the approval-only review phase, so it
  is release-eligible. Alzheimer/Parkinson GloScope remain the only
  uncorrected GloScope rows missing from the canonical result root.
- **Alzheimer/Parkinson GloScope attempt was canceled and audited
  fail-closed:** Gate
  `ecoda_bassez_rolling_b5_gloscope_alzheimer_parkinson_d558b26_20260904121500Z`
  launched with command digest
  `e82af5d86041dc5dbdeb1415bcb7788235fdcee84febb580e9c1f03147a4e61f`
  and reached a durable transport failure at `2026-09-04T14:47:00Z`.
  Preflight array `4379822` completed successfully for both H5AD rows;
  GloScope array `4379874`, watchdog `4379875`, and aggregate `4379876`
  were submitted, but no GloScope task started before cancellation.
  Following the user's explicit authorization, preflight `4379822`, GloScope
  array `4379874`, and watchdog `4379875` were sent to `scancel`; aggregate
  `4379876` was not canceled and subsequently reached `FAILED|1:0`. The named
  tmux session was killed and verified runner PID `2975098` received
  `SIGTERM`. Reconciliation then found no runner process or tmux session, and
  the remote terminal status was `FAILED` with exit
  `143` at `2026-09-04T15:34:08Z`.

- The required single first inspect completed at `2026-09-04T15:35:29Z`
  with one accounting query over all four IDs. It matched preflight
  `4379822` as `COMPLETED|0:0`, array `4379874` and watchdog `4379875` as
  `CANCELLED|0:0`, and aggregate `4379876` as `FAILED|1:0`; therefore
  `audit_passed=false` and `release_eligible=false`, while all individual
  artifact/terminal/fingerprint checks passed. No reviewer approval was
  issued and this failed `PRELAUNCH_STOP` gate is not a predecessor. The
  replacement must be a new explicitly reconciled gate using the reviewed
  Joanito predecessor.

- **New reroute implementation:** The user authorized a reusable targeted
  batch recovery path and a rerun on `shared-bigmem`, rather than waiting for
  the shared-cpu queue. Commit `fdf137d` now preserves the fixed seven-method
  batch default while adding selection-file-scoped
  `--target-methods LIST` recovery with `--pass`, explicit partition, and
  existing owner/preflight/watchdog/validation/sync machinery. It also makes
  GloScope, composition, PILOT, QOT, and PILOT-GM-VAE use a genuinely
  counts-free h5py/minimal-AnnData loader; MrVI, scPoli, and pseudobulk retain
  counts only where their algorithms require them. Focused submitter,
  Python-loader/worker, R-loader, syntax, and H5AD contract checks passed.
  The source is pushed and a source-matched runtime gate is now running
  before any replacement launch.

- **Alzheimer prepare-pseudobulk remains blocked by a real memory floor:** The
  failed batch attempts load the complete 1,395,601-cell by 34,800-gene counts
  matrix into Seurat through `load_benchmark_seurat`, then
  `get_pb_deseq2` calls `AggregateExpression` on the full object before
  selecting the requested hvg2000 result. The affected task OOMed at the
  configured 500G ceiling. A blind retry would repeat the same full-object
  allocation; the next correction must preserve all-gene aggregation and
  DESeq2 semantics through a bounded/streaming path, then use a new
  source-matched runtime and targeted Alzheimer dependency repair.
- **Runtime build retries are fail-closed and preserved:** The first fdf137d
  runtime gate was invalid because it captured Bamboo HEAD `d558b263` before
  a post-launch fast-forward; its runner was terminated and no image from it
  is source-matched. The next gate
  `ecoda_runtime_build_bassez_rolling_fdf137d_final_20260904185000Z` correctly
  captured Bamboo HEAD
  `fdf137d0c486e5edf884c4d1851d62aa6e72e655`, but its guarded build refused
  environment mutation because stale runtime-build job `4380369` was still
  active. Build job `4380402` failed with `1:0`; the verified jobs were
  canceled, and the required single inspect ran at `2026-09-04T16:28:12Z`
  with one accounting query. That gate is terminal `FAILED` and not
  release-eligible. No failed image is adopted.

- **Source-matched runtime build stopped fail-closed after transport recovery:** Gate
  `ecoda_runtime_build_bassez_rolling_fdf137d_clean_20260904194000Z` launched
  at `2026-09-04T16:34:06Z` with command digest
  `214fec748034e5055452bb0623c6dcda62220b5a1d1c05dd6fa4753b90def278`.
  Prepare/reconcile captured Bamboo HEAD
  `fdf137d0c486e5edf884c4d1851d62aa6e72e655`, confirmed no stale runtime-build
  process, and launched the one 128G shared-cpu build with an in-wrapper HEAD
  assertion. The durable waiter hit one SSH completion transport error at
  `2026-09-04T16:39:59Z`, so the local manifest entered terminal
  `PRELAUNCH_STOP` rather than guessing. Recovery status later proved remote
  `COMPLETED` at `2026-09-04T16:42:17Z`; scheduler ID `4380472` was extracted
  from the durable wrapper log.

- The required single terminal inspect ran at `2026-09-04T19:07:51Z` with one
  accounting query for `4380472`, which matched `COMPLETED|0:0`; all five
  artifact contracts, terminal commands, and immutable fingerprints passed.
  Because the local lifecycle is irreversibly `PRELAUNCH_STOP`, reviewer
  approval is prohibited and the gate remains release-ineligible. The produced
  image is retained as evidence but is not adopted as a reviewed runtime.
  A fresh source-matched runtime gate must rebuild/review before the
  shared-bigmem replacement.
- **Recovery runtime gate is now running:** New gate
  `ecoda_runtime_build_bassez_rolling_fdf137d_recovery_20260904190825Z` was
  prepared and launched after the stopped gate's terminal evidence was
  recorded. It asserts Bamboo HEAD
  `fdf137d0c486e5edf884c4d1851d62aa6e72e655`, rebuilds with
  `build_ecoda_runtime.sh --layout path-preserving --force`, and writes
  `ecoda-py-cuda13-path-preserving-fdf137d-reviewed.sif`. Its single durable
  waiter is armed. No benchmark replacement may launch until this gate's
  terminal inspect passes and Luna Max approves it.
- **Recovery runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_fdf137d_recovery_20260904190825Z`
  completed at `2026-09-04T19:17:26Z`. Its single terminal inspect at
  `2026-09-04T19:18:44Z` covered scheduler job `4381441` with one accounting
  query; accounting, all five artifact contracts, terminal commands, and
  immutable fingerprints passed. Luna Max approved at `2026-09-04T19:19:44Z`;
  the gate is release-eligible. The reviewed image is
  `ecoda-py-cuda13-path-preserving-fdf137d-reviewed.sif`.

- **GloScope replacement scope expanded to both missing rows:** The user
  selected one shared-bigmem gate for only GloScope, counts-free, across
  Alzheimer and Parkinson. The earlier prepared Alzheimer-only manifest
  `ecoda_bassez_rolling_b5_alzheimer_gloscope_bigmem_fdf137d_20260904192119Z`
  is superseded and remains unlaunched; it is not a predecessor.
  The replacement selection file
  `/home/users/h/halterc/scratch/ECODA_paper/gates/ecoda_bassez_rolling_b5_gloscope_alzheimer_parkinson_bigmem_fdf137d_20260905155916Z.selection.tsv`
  contains exactly:
  `Alzheimer<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected` and
  `Parkinson<TAB>batch_effect_uncorrected<TAB>batch_effect_uncorrected`.
  Its SHA-256 is
  `877a7437fddb1b9d72f97d68aac0fe1a28d06c846b6dfd4f7d18da2ae02554fe`.
  The replacement gate will invoke the canonical submitter with
  `--pass uncorrected --target-methods gloscope --partition shared-bigmem`,
  `500G` memory, and throttle `1`, using the reviewed runtime and the
  reviewed Joanito Stage 4 and uncorrected Pipeline 5 predecessors. No other
  method or dataset is selected.
- **Dual GloScope bigmem submission failed closed on an invalid CPU request:**
  Gate
  `ecoda_bassez_rolling_b5_gloscope_alzheimer_parkinson_bigmem_fdf137d_20260905155916Z`
  reached durable `FAILED` at `2026-09-05T16:52:24Z`. Both selected H5ADs
  passed source checksum validation; H5AD preflight `4382272` and R
  environment preflight `4382274` were submitted, but the GloScope array
  submission requested the configured 16 CPUs on `shared-bigmem`, whose
  effective nodes expose 14 allocatable CPUs (`CPUSpecList=7,15`). Slurm
  rejected that request with `Requested node configuration is not available`;
  no GloScope array, watchdog, or aggregate gate was submitted.

- The required single terminal inspect ran at `2026-09-05T16:54:08Z`; its one
  accounting query covered the scheduler ID emitted in the durable log
  (`4382274`), which was `COMPLETED|0:0`, and all five artifact contracts,
  terminal commands, and immutable fingerprints passed. The gate remains
  `FAILED` and release-ineligible because the wrapper failed before method
  submission. A non-mutating `sbatch --test-only` confirmed that the same
  shared-bigmem request is schedulable with `--cpus-per-task=14`.
  The next replacement keeps the reviewed runtime and exact dual selection,
  and exports `BENCHMARK_CPU_CPUS_PER_TASK=14` before sourcing the canonical
  config; no source or lockfile change is needed.

- **CPU-adjusted GloScope gate failed on the pre-existing R-loader import
  bug:** Gate
  `ecoda_bassez_rolling_b5_gloscope_alzheimer_parkinson_bigmem_cpu14_fdf137d_20260905165638Z`
  reached durable `FAILED` at `2026-09-05T17:52:35Z`. Both H5ADs passed
  source validation; preflights `4382279`/`4382282`, GloScope array
  `4382283`, watchdog `4382284`, and aggregate `4382285` were recorded.
  The method array failed twice because the fdf137d runtime's R
  `load_h5ad_counts_free()` imported `benchmark_h5ad_contract` through a
  relative-import fallback with no package parent. The aggregate and
  watchdogs consequently failed; no GloScope artifact was accepted.

- The required single terminal inspect ran at `2026-09-05T17:54:01Z` with one
  accounting query covering all five recorded scheduler IDs. Preflights
  `4382279`/`4382282` were `COMPLETED|0:0`; array `4382283`, watchdog
  `4382284`, and aggregate `4382285` were `FAILED|1:0`. The audit therefore
  failed and the gate is not a predecessor. Extra accounting rows
  `4382280` and `4382286` were preserved as scheduler evidence.

- **Loader fix committed and synchronized:** Commit `bc021d9` adds the
  module directory to reticulate's Python `sys.path` before
  `import_from_path`, and extends the R integration regression to call the
  persisted H5AD contract path. `pixi run python tests/test_h5ad_counts_free.py`
  passes locally, and Bamboo is synchronized to
  `bc021d92b6bf5e157c5bb96315d5eba08ddde1f6`. The reviewed fdf137d runtime is
  source-stale for this fix; a new source-matched runtime must be rebuilt and
  reviewed before another GloScope launch.

- **Loader-fix runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_bc021d9_loaderfix_20260905180014Z`
  completed at `2026-09-05T18:10:49Z`; its single inspect at
  `2026-09-05T18:11:56Z` covered scheduler job `4382288`, with accounting,
  all five artifact contracts, terminal commands, and immutable fingerprints
  passing. Luna Max approved at `2026-09-05T18:12:59Z`; the gate is
  release-eligible and its image is
  `ecoda-py-cuda13-path-preserving-bc021d9-reviewed.sif`.

- **Loader-fixed GloScope recovery completed and reviewed:** Gate
  `ecoda_bassez_rolling_b5_gloscope_alzheimer_parkinson_bigmem_cpu14_bc021d9_20260905181320Z`
  completed at `2026-09-05T23:25:05Z`. Its single inspect at
  `2026-09-05T23:26:19Z` covered preflights `4382320`/`4382322`, method array
  `4382323`, watchdog `4382324`, and aggregate gate `4382325`; every required
  accounting row was `COMPLETED|0:0`, and all five artifact contracts,
  terminal audits, immutable fingerprints, and the overall audit passed.
  Luna Max approved at `2026-09-05T23:27:29Z`; the gate is release-eligible.

- The canonical wrapper merged two GloScope task logs into a 16-row execution
  log, validated the Alzheimer and Parkinson RDS bundles, and synchronized
  both outputs plus checksums to NAS:
  `Alzheimer_batch_effect_uncorrected_gloscope.rds` and
  `Parkinson_batch_effect_uncorrected_gloscope.rds`. The selected R worker
  used the bc021d9 counts-free H5AD loader; no counts layer was materialized
  for GloScope. No other method or dataset was selected. Uncorrected
  GloScope is now closed; the remaining uncorrected blocker is Alzheimer
  `prepare_pseudobulk` and its dependent pseudobulk/composition rows.
- **Alzheimer pseudobulk memory-safe correction committed:** Commit
  `d83ac1d` adds an h5py-only, bounded CSR aggregator that retains only the
  sample-by-gene pseudobulk matrix and selected sample metadata. The prepare
  worker now builds a sample-level Seurat object from that aggregate, so it
  preserves the existing `AggregateExpression`/DESeq2 normalization semantics
  without materializing the 1,395,601-cell count matrix. MOFA and batch-mode
  pseudobulk result workers also use sample metadata/count-free paths when
  cached pseudobulks suffice, while ordinary CT pseudobulk and scITD retain
  their required cell counts. Focused Python, R-helper, and actual prepare
  worker regressions passed; Bamboo is synchronized to
  `d83ac1d`.

- **Pseudobulk source-matched runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_d83ac1d_pseudobulk_20260905234243Z`
  completed at `2026-09-05T23:53:44Z`; its single inspect at
  `2026-09-05T23:54:33Z` covered scheduler job `4382417` with one passing
  accounting query, all five artifact contracts, terminal commands, and
  immutable fingerprints. Luna Max approved at `2026-09-05T23:55:36Z`; the
  gate is release-eligible and its image is
  `ecoda-py-cuda13-path-preserving-d83ac1d-reviewed.sif`.

- **Alzheimer prepare-pseudobulk recovery completed and reviewed:** Gate
  `ecoda_bassez_rolling_b5_alzheimer_prepare_pseudobulk_bigmem_d83ac1d_20260905235603Z`
  completed at `2026-09-06T00:46:26Z`. Its single inspect at
  `2026-09-06T00:47:31Z` covered preflights `4382422`/`4382423`, method array
  `4382424`, watchdog `4382425`, and aggregate gate `4382426`; every required
  accounting row was `COMPLETED|0:0`, and all five artifact contracts,
  terminal audits, immutable fingerprints, and the overall audit passed.
  Luna Max approved at `2026-09-06T00:48:23Z`; the gate is release-eligible.
  The wrapper produced and synchronized
  `Alzheimer_batch_effect_uncorrected_pseudobulk_hvg2000.rds` and its
  checksum. The bounded CSR aggregation replaced the prior full-count OOM.

- **Alzheimer dependent pseudobulk/composition recovery completed and
  reviewed:** Gate
  `ecoda_bassez_rolling_b5_alzheimer_pseudobulk_composition_bigmem_d83ac1d_20260906004850Z`
  completed at `2026-09-06T01:36:09Z`; its single inspect at
  `2026-09-06T01:37:13Z` covered preflights `4382433`/`4382441`, arrays
  `4382442`/`4382444`, watchdogs `4382443`/`4382445`, and aggregate
  `4382446`. All required accounting rows were `COMPLETED|0:0`; all five
  artifact contracts, terminal audits, immutable fingerprints, and the
  overall audit passed. Luna Max approved at `2026-09-06T01:38:23Z`; the gate
  is release-eligible. It synchronized the Alzheimer pseudobulk,
  composition, and metadata RDS bundles plus checksums.

- **Full uncorrected batch audit exposed two preserved failures:** Gate
  `ecoda_bassez_rolling_b5_batch_uncorrected_audit_d83ac1d_20260906014011Z`
  reached durable `FAILED` at `2026-09-06T05:47:48Z`. Strict validation
  skipped intact rows and submitted the remaining invalid artifacts. The
  composition array `4382483`/watchdog `4382484` failed because Parkinson's
  configured high-resolution column
  `leiden_res_5_batch_effect_uncorrected_hvg2000` was renamed in place to
  `RNA_snn_res.5` before the caller used the configured name. The MrVI array
  `4382485` hit `OUT_OF_MEMORY|0:125`; its watchdog recorded a completed
  `4382711` retry, but the original array remains non-`COMPLETED` and cannot
  release this gate. Other retry rows are preserved.

- The required single terminal inspect ran at `2026-09-06T05:49:22Z` with one
  accounting query covering all 15 recorded IDs. It matched the composition
  and original MrVI failures above plus aggregate gate `4382491 FAILED|1:0`;
  the audit failed and no reviewer approval was issued. Commit `eafff66`
  now preserves configured Leiden source columns while adding the legacy
  aliases; focused mapping, pseudobulk, and mode regressions passed. Bamboo
  is synchronized to `eafff66`. A source-matched runtime rebuild/review is
  required before a new no-force full uncorrected audit; no unchanged rerun
  is allowed.

- **Composition-fix runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_eafff66_compositionfix_20260906055522Z`
  completed at `2026-09-06T06:05:54Z`; its single inspect at
  `2026-09-06T06:06:48Z` covered scheduler job `4382785`, with accounting,
  all five artifact contracts, terminal commands, and immutable fingerprints
  passing. Luna Max approved at `2026-09-06T06:07:53Z`; the gate is
  release-eligible.

- **Second full uncorrected audit failed only on missing Alzheimer MrVI:**
  Gate
  `ecoda_bassez_rolling_b5_batch_uncorrected_audit_eafff66_retry_20260906060835Z`
  reached durable `FAILED` at `2026-09-06T11:53:31Z`. The corrected
  composition mapping passed; strict validation and reruns produced the
  remaining method artifacts, but the Alzheimer MrVI Feather remained
  missing. The original MrVI array `4382873` had two `OUT_OF_MEMORY|0:125`
  rows, and retry array `4383121` completed while its retry task log
  processed Joanito even though
  `batch_effect_uncorrected__mrvi.retry_1.tsv` contained Alzheimer/Parkinson.
  This retry-manifest/worker mismatch is preserved as a failed recovery, not
  treated as successful Alzheimer completion.

- The required single terminal inspect ran at `2026-09-06T11:55:54Z` with one
  accounting query covering all 18 emitted scheduler IDs. It matched 90
  `COMPLETED|0:0` rows and two preserved `OUT_OF_MEMORY|0:125` rows
  (`4382873` and extra task `4382883`); the gate audit failed and no reviewer
  approval was issued. Commit `707468e` now binds OOM retries explicitly to
  `MATRIX_RETRY_MANIFEST` and streams only selected HVG raw counts for
  MrVI/scPoli instead of materializing all genes. Focused count-loader,
  retry-manifest, pseudobulk, and batch-mode regressions passed; Bamboo is
  synchronized to `707468e`.

- **MrVI-subset runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_707468e_mrvi_subset_20260906120918Z`
  completed at `2026-09-06T12:20:47Z`; its single inspect at
  `2026-09-06T12:22:00Z` covered scheduler job `4383277`, with accounting,
  all five artifact contracts, terminal commands, and immutable fingerprints
  passing. Luna Max approved at `2026-09-06T12:23:25Z`; the gate is
  release-eligible and its image is
  `ecoda-py-cuda13-path-preserving-707468e-reviewed.sif`.

- **Targeted Alzheimer MrVI recovery completed and reviewed:** Gate
  `ecoda_bassez_rolling_b5_alzheimer_mrvi_gpu_707468e_20260906122358Z`
  completed at `2026-09-06T15:43:07Z`. Its single inspect at
  `2026-09-06T15:44:55Z` covered preflight `4383294`, method array `4383311`,
  watchdog `4383312`, and aggregate `4383313`; all required accounting rows
  were `COMPLETED|0:0`, and all artifact contracts, terminal audits,
  immutable fingerprints, and the overall audit passed. Luna Max approved at
  `2026-09-06T15:45:59Z`; the gate is release-eligible. It produced and
  synchronized `Alzheimer_batch_effect_uncorrected_hvg2000_highres_mrvi_dists.feather`
  with its checksum using the selected-HVG loader.

- The corrected retry binding was exercised by the targeted run without an
  OOM retry; the full current-source uncorrected audit can now be repeated
  without force to close the branch. No broad recomputation is allowed.

- **Final uncorrected audit wrapper failed on Stephenson composition:** Gate
  `ecoda_bassez_rolling_b5_batch_uncorrected_final_707468e_20260906154643Z`
  reached remote terminal `FAILED` at `2026-09-06T22:36:11Z` after all method
  scheduler chains completed. RDS validation rejected
  `Stephenson_batch_effect_uncorrected_composition.rds` because its three
  combo keys omitted the expected `ECODA_HiTME_HR_layer2` and
  `ECODA_scATOMIC_HR` entries. The source H5AD contains both annotation
  columns; the counts-free R worker simply had not requested them.

- The completion transport outage was recovered through the durable status
  record. The required single inspect ran at `2026-09-07T10:08:17Z` with one
  accounting query covering preflights `4383592`/`4383906`, arrays
  `4383937`/`4383946`/`4383948`/`4383950`/`4383952`/`4383954`/`4383956`,
  watchdogs `4383938`/`4383947`/`4383949`/`4383951`/`4383953`/`4383955`/
  `4383957`, and aggregate `4383958`; all required rows were
  `COMPLETED|0:0`, and the profile audit passed. The local lifecycle remains
  `PRELAUNCH_STOP` because the wrapper exited nonzero, so no reviewer approval
  was issued and this gate is not release-eligible.

- **Batch contract correction committed:** The earlier `6c05613` runtime
  build completed at `2026-09-07T10:23:21Z` but was not inspected, reviewed,
  or adopted. After clarification, commit `0c25455` narrows batch composition
  generation to `ECODA_authors_HR`, `ECODA_authors_HR_NULL`, and
  `ECODA_seuratres_2`; the validator requires those three and permits only
  the recognized legacy extras `ECODA_HiTME_HR_layer2` and
  `ECODA_scATOMIC_HR`. Ordinary benchmark composition retains its annotation
  outputs, and batch workers no longer request those optional columns.
  Existing extra bundles are ignored for validity beyond that allowlist. The
  focused batch-composition and RDS-contract regressions passed; Bamboo is
  synchronized to `0c25455`.

- **Batch-contract source fix is now complete:** Commit `12373b3` restores
  the root-mode validator's batch `stem`/result-file path construction and
  adds root-mode coverage proving base-only and allowlisted legacy-extra
  composition bundles pass while missing/unknown keys fail. The focused
  regression passes, and Bamboo is synchronized to
  `12373b3d3409c06c2caef48a2db448a6570f6cbb`.

- **Batch-contract runtime reviewed:** Gate
  `ecoda_runtime_build_bassez_rolling_12373b3_batch_contract_20260907104734Z`
  completed at `2026-09-07T10:59:03Z`; its single inspect at
  `2026-09-07T11:00:21Z` covered scheduler job `4385036`, with accounting,
  all five artifact contracts, terminal commands, and immutable fingerprints
  passing. Luna Max approved at `2026-09-07T11:01:45Z`; the gate is
  release-eligible and its image is
  `ecoda-py-cuda13-path-preserving-12373b3-reviewed.sif`.

- **Final no-force uncorrected audit was canceled as overbroad:** Gate
  `ecoda_bassez_rolling_b5_batch_uncorrected_final_12373b3_20260907110212Z`
  launched at `2026-09-07T11:03:01Z` with the reviewed 12373b3 runtime and
  exact twelve-row selection, but its exact wrapper included all seven
  methods. It skipped the 48 validated `prepare_pseudobulk`, `pseudobulk`,
  `gloscope`, and `composition` artifacts, then submitted MRVI/PILOT/QOT
  arrays and watchdogs (`4385850`/`4385852`, `4385853`/`4385865`, and
  `4385866`/`4385867`). The H5AD and R preflights were `4385714` and
  `4385842`. This was an unnecessary compute launch for a contract repair.

- **Overbroad gate contained:** After explicit user approval, `scancel` was
  issued for `4385714`, `4385842`, `4385850`, `4385852`, `4385853`, `4385865`,
  `4385866`, and `4385867`, and the named durable tmux runner was terminated.
  The existing durable waiter recorded terminal `FAILED` with exit code 143
  at `2026-09-07T13:17:13Z`. The single first inspect ran at
  `2026-09-07T13:18:43Z`, issued one accounting query over the six
  array/watchdog IDs, and recorded `audit_state: COMPLETED`,
  `audit.passed: false`, `state: FAILED`, and `release_eligible: false`.
  PILOT array/watchdog accounting rows were completed; MRVI and QOT
  array/watchdog rows were canceled. No reviewer approval or dependent gate
  was started. Run-owned logs, manifests, and partial artifacts were
  preserved; the preflight IDs were canceled but were not included in the
  typed array/watchdog accounting list.

- **No-unnecessary-compute guardrail strengthened:** `AGENTS.md` now makes
  repair/validation no-compute by default, prohibits broad method or
  all-dataset selections and `--force` for repair work, requires explicit
  dataset/view/method row counts and rationale before compute, and mandates
  immediate cancellation plus one failed inspect when emitted scope exceeds
  approval.
- **Batch-unccorrected local snapshot transfer completed:** Per the user's
  explicit request, the complete current scratch tree
  `$HOME/scratch/ECODA_paper/batch_effect/uncorrected/` was copied with
  `rsync` to the local mirror
  `/Users/christianhalter/Desktop/ECODA_paper/data/batch_effect/uncorrected/`.
  The transfer included `results/`, `embeddings/`, `pseudobulks/`,
  `gloscope_dists/`, and the existing `checksums.md5` plus sidecars:
  348 files total, 343 transferred, 101801831 bytes. This records completion
  of the requested local snapshot, not a new HPC computation or a passing
  durable benchmark gate; existing partial/unpromoted artifacts remain
  preserved as-is.
- A NAS overlay was attempted after the initial scratch copy and transferred
  10 changed files; it was not treated as authoritative. The final
  authoritative local sync used `rsync -a --delete` from the HPC scratch
  source to the local destination. It removed only local NAS-only entries
  (`@eaDir/`, its Synology metadata, `.DS_Store`, and the NAS-only Alzheimer
  MRVI sidecar); it did not contact or modify NAS. No file bytes needed
  retransmission, and a final `rsync -ani --delete` reported no differences.
  The local destination now follows the HPC scratch tree; existing partial or
  unpromoted source artifacts remain preserved as-is.
