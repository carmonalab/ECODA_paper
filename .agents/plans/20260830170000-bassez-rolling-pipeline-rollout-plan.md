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
- The rollout intentionally excludes corrected batch-effect views. If a corrected view is requested later, require a separate plan and confirm `columns.batch` before submission.
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
