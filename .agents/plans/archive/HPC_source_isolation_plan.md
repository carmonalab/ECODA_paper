# Parallel HPC Source Isolation and Check Simplification

## Context

The current repository HEAD captured during planning is `5302671ad94556edcf9acccf372d2dc34121d714`. The user-confirmed production benchmark (`benchmark_analysis`) baseline has completed pipelines 1–5. The `batch_effect_uncorrected` baseline is **not terminally complete** after a recent user update to `datasets.json`: `Lupus_PBMC` and `Parkinson` changed `cell_type_high_res`, so methods depending on that input must be rerun, and the batch selection must use the new `Kidney_KPMP_full` dataset key. This status is a user-provided operational fact and must not be inferred from stale gate files. The requested end state is that a running job executes from immutable source, dependency, auxiliary-data, and runtime-image identities, a later commit can be pulled and launched independently, redundant full-file hashing and global stale-gate checks are removed, and future methods remain explicitly pending until selected and run. `datasets.json` is read-only for this change and must not be edited.

## Approach

### 1. Record the completed baseline and pending-method rule

Update the root `AGENTS.md` with a `Current processing baseline` section immediately after the durable-HPC rules. Record exactly:

- Baseline anchor: `5302671ad94556edcf9acccf372d2dc34121d714` (the HEAD captured while preparing this plan); preserve this value even if implementation occurs at a later commit.
- User-confirmed completed production benchmark view (`benchmark_analysis`) datasets, derived from the current `datasets.json` flags: `Adams`, `Bassez`, `Gongsharma_cmv_young_males`, `Kfoury`, `Kim`, `Lee`, `Pelka`, `Smillie`, `Stephenson`, `Wu`, and `Zhang`.
- The historical batch-effect selection was processed, but `batch_effect_uncorrected` is not terminally complete after the user's recent `datasets.json` update. The changed inputs are `cell_type_high_res` for `Lupus_PBMC` and `Parkinson`; only downstream methods that consume that high-resolution field must be explicitly rerun. The required current batch dataset key is `Kidney_KPMP_full`, not `Kidney_KPMP`. Record the configured batch-effect dataset list as `Joanito`, `Stephenson`, `CombinedPBMC`, `Alzheimer`, `Breast_cancer`, `Covid19_PBMC`, `Kidney_KPMP_full`, `Myocardial_infarction`, `Diabetes`, `Lupus_PBMC`, `Lung`, and `Parkinson`, with the two changed datasets and the new Kidney key marked pending targeted validation/rerun.
- `_debug` is a five-sample verification fixture, not a production cohort; record it separately as the routine verification dataset for both configured views.
- `Alzheimer`, `Diabetes`, and `Parkinson` remain covered by the existing `not_suitable_for_auto_annotation` exemption; their historical batch processing status does not imply that automatic HiTME/scATOMIC annotation was required.
- Baseline benchmark methods are `gloscope`, `mofa`, `pseudobulk`, `composition`, `scitd`, `mrvi`, `scpoli`, `pilot`, `qot`, and `pilotgm`; baseline analyses are `trans` and `zeroimp`; the batch-effect suite is `prepare_pseudobulk`, `pseudobulk`, `gloscope`, `composition`, `mrvi`, `pilot`, and `qot`.
- This baseline does not include benchmark methods or scripts added later. The Stage 5 default method list must remain the baseline list above; a newly registered method is not a default. A post-baseline method is runnable only when named explicitly with `--methods` or an explicit selection manifest. Existing valid benchmark rows and unaffected batch rows are skipped individually without `--force`; changed `Lupus_PBMC`/`Parkinson` high-resolution consumers, `Kidney_KPMP_full`, and absent new-method rows are submitted only by explicit targeted selection.
- Gate history is evidence only. Submitters and recovery paths must decide reuse/recompute from the selected artifact contract and explicit user scope, never from the existence of an old `FAILED`, `PRELAUNCH_STOP`, or stale gate manifest.

Replace the current direct full-cohort command examples in `AGENTS.md` with the source-snapshot plus durable-gate wrapper workflow. Label direct submitter commands as legacy/validator-only unless they carry a required immutable source/runtime manifest. Document that a `datasets.json` change can invalidate only dependent downstream rows; it must not trigger a blanket rerun. Do not copy the dataset lists into `datasets.json` or change its flags. If HEAD has advanced before this documentation edit, retain the baseline anchor above and add the actual documentation commit separately; do not replace the baseline anchor.

### 2. Create a self-verifying immutable source snapshot

Add `src/utils/bash/ecoda_source_snapshot.sh`; no existing utility creates a source snapshot. Implement two explicit subcommands:

- `create --source-root ABSOLUTE_REPOSITORY --snapshot-parent ABSOLUTE_PARENT --commit FULL_COMMIT`
- `exec --source-root ABSOLUTE_SNAPSHOT_TREE --source-manifest ABSOLUTE_MANIFEST --host-env-prefix ABSOLUTE_ENV_PREFIX --runtime-image ABSOLUTE_IMAGE --runtime-manifest ABSOLUTE_IMAGE_MANIFEST --run-id RUN_ID --scratch-root ABSOLUTE_SCRATCH_ROOT --logs-root ABSOLUTE_LOG_ROOT --script RELATIVE_SCRIPT -- [SCRIPT_ARGS...]`

Use the full source commit as the immutable snapshot ID so multiple runs on the same source reuse one verified snapshot:

```text
${HPC_SCRATCH_DIR}/_ecoda_source_snapshots/<SOURCE_COMMIT>/
├── COMPLETE
├── tree/                 # complete committed source, including tracked aux/
└── identity/
    ├── source.manifest
    └── source.tar
```

`create` must:

1. Require absolute paths, a safe snapshot parent, an existing Git repository, and a full 40-hex commit resolved by `git -C SOURCE_ROOT rev-parse --verify`. For each required file, require `git cat-file -e <COMMIT>:<PATH>` so an older requested commit cannot produce a manifest from files that exist only in the current HEAD.
2. Require a clean source tree (`git status --porcelain --untracked-files=all` must be empty) so newly added methods/scripts are committed before they become runnable.
3. Require `src`, `datasets.json`, `config_helper.R`, `pixi.toml`, `pixi.lock`, and the three tracked auxiliary files `aux/scGateDB.rds`, `aux/genes.blocklist.rds`, and `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` in the requested commit/tree. Keep `tree/aux` as the sole frozen auxiliary root; do not create a sibling aux copy or a second aux checksum pipeline.
4. Create `source.tar` from the requested commit, persist it in the temporary snapshot directory, hash it with SHA-256, extract it into `tree` (including `tree/aux`), and compare the extracted content/file list with a fresh extraction before applying read-only permissions. Hash `config_helper.R`, `datasets.json`, `pixi.toml`, and `pixi.lock` from the verified extraction, not from the mutable checkout.
5. Install `tree`, `identity/source.tar`, `identity/source.manifest`, and `COMPLETE` with one same-filesystem temporary-parent rename. Do not publish `tree` before its manifest/archive exists. Make the entire published snapshot read-only. An existing commit snapshot must fail closed unless its manifest and archive digest match; never overwrite it.
6. Write a strict key/value manifest with exactly these fields: `FORMAT=1`, `SOURCE_ROOT`, `SOURCE_COMMIT`, `SOURCE_ARCHIVE_PATH`, `SOURCE_ARCHIVE_SHA256`, `CONFIG_HELPER_SHA256`, `DATASETS_SHA256`, `PIXI_TOML_SHA256`, `PIXI_LOCK_SHA256`, `AUX_ROOT`, and `SCGATE_DB_BRANCH`. Runtime-image and host-environment identities are run-bound, not part of the commit-keyed source snapshot cache key.

`exec` must:

- require `SOURCE_ROOT` to be `<snapshot-id>/tree` and the manifest/archive to be under the same `<snapshot-id>/identity`; require `AUX_ROOT` to be exactly `${SOURCE_ROOT}/aux`; reject symlinks or paths escaping that parent;
- require `COMPLETE`, a read-only source tree, no `.git`, a valid relative script under `tree`, and an absolute existing host environment prefix when host-side execution is used;
- rehash the retained `source.tar`, extract it into a temporary directory with the archive’s normal modes, compare every source/aux file and directory against `tree`, and only then check that the published tree is read-only. Never compare raw archive modes with post-`chmod` modes;
- validate the versioned runtime image/manifest paths and required image manifest fields before executing;
- export `ECODA_SOURCE_ROOT`, `ECODA_SOURCE_MANIFEST`, `ECODA_SOURCE_SNAPSHOT_REQUIRED=1`, `ECODA_HOST_ENV_PREFIX`, `ECODA_RUNTIME_IMAGE`, `ECODA_RUNTIME_MANIFEST`, `ECODA_RUN_ID`, `HPC_SCRATCH_DIR=${scratch-root}`, `ECODA_SCRATCH_ROOT=${scratch-root}`, `ECODA_LOGS_DIR=${logs-root}`, and `ECODA_AUX_ROOT=${ECODA_SOURCE_ROOT}/aux`;
- execute `/bin/bash "${ECODA_SOURCE_ROOT}/${script}"` and never source the mutable canonical checkout.

It must not create a scheduler job or create `${HPC_SCRATCH_DIR}/_ecoda_runs/<RUN_ID>`; the stage submitter still creates that run root later. It may create only the explicitly supplied run-specific log directory. A missing/mismatched snapshot, manifest, archive, runtime identity, environment prefix, or log path fails before submission.

Create a commit-keyed snapshot before durable-gate `prepare`; the exact wrapper references the snapshot’s own `ecoda_source_snapshot.sh`, not the mutable checkout’s copy. Keep one immutable snapshot per full source commit. A new runtime image or host environment does not invalidate that source snapshot; its identity is checked separately for each run.

### 3. Separate image-build, bound-source, dependency, host-environment, and run-image identities

Update `src/utils/bash/build_ecoda_runtime.sh`, `src/utils/bash/ecoda_runtime.sh`, `src/slurm_config.sh`, and `src/utils/py/gene_utils.py` together. The current builder copies the realized environment into the Apptainer image but records the build Git revision while later binding live source; a format-2 image therefore records its own build identity separately from the source revision executed by a job. The current `gene_utils.py:11–13` also derives `aux/` from the source tree, so the source/aux layout must remain internally consistent.

Change new runtime manifests to `FORMAT=2` with these fields:

- `IMAGE_BUILD_GIT_REVISION` — provenance only for the checkout used to build the environment image;
- `IMAGE_SHA256`, `IMAGE_PATH`, and the existing runtime/layout/toolchain fields;
- `IMAGE_PIXI_TOML_SHA256` and `IMAGE_PIXI_LOCK_SHA256` — dependency identity of the embedded environment.

During image construction, hash `pixi.toml` and `pixi.lock` before and after `pixi containerize`; abort if either changes. New format-2 images require `--runtime-id RUNTIME_ID` and an output under `${HPC_SCRATCH_DIR}/_ecoda_runtime/<RUNTIME_ID>/`; the builder rejects an existing format-2 output even when `--force` is supplied. A later image build must use a new runtime ID/path. The image and its manifest are published read-only and their parent directory is non-writable for the run. Keep `FORMAT=1` parsing and behavior for existing images/current legacy gates; do not reinterpret old `GIT_REVISION` fields.

At the first submission boundary, write `${ECODA_RUN_ROOT}/manifests/runtime.identity` with exactly these run-bound values: `RUNTIME_IMAGE`, `RUNTIME_MANIFEST`, `RUNTIME_IMAGE_SHA256`, `RUNTIME_MANIFEST_SHA256`, `RUNTIME_IMAGE_SIZE`, and `RUNTIME_MANIFEST_SIZE`, plus the image dependency fields. Do not record or compare `st_dev`/device IDs or prescribe GNU-only `stat` fields: distributed filesystems can report different device values on different compute nodes. Submission performs the sole full hash of the large image and records `RUNTIME_IMAGE_SHA256`; later workers rehash only the small runtime manifest and check the exact image/manifest paths, sizes, read-only permissions, and non-writable parent directory. They do not claim a current image-content SHA comparison. A replacement or manifest edit therefore fails before retry through the versioned non-overwritable path, manifest digest, size, or permissions.

In `ecoda_runtime.sh`:

1. Replace the current live-checkout comparison in `_ecoda_runtime_validate_source_identity` (`:116–141`) with validation of the self-verifying `ECODA_SOURCE_MANIFEST`/`ECODA_SOURCE_ROOT`, source/config/aux digests, host-environment identity, and dependency compatibility. For Apptainer, source `PIXI_TOML_SHA256` and `PIXI_LOCK_SHA256` must equal the image dependency fields. A source-only commit change with identical dependency fields is allowed; a lock/runtime dependency change requires a new matching image/environment.
2. Keep `ecoda_runtime_validate_submission MODE` as the one full submission-time image/environment validation. It validates image bytes, image manifest, versioned image path, source archive, dependency compatibility, and runtime layout, then writes `runtime.identity`. For `FORMAT=2`, it must not compare the mutable canonical checkout’s `git HEAD` with the image build revision.
3. Add `ecoda_runtime_validate_bound_run()` with no positional arguments. It validates the run-bound source/archive and runtime identity, exact containment, read-only permissions, the small runtime-manifest SHA-256, image/manifest path and size recorded at submission, and dependency fields without rehashing the large image. Submission performs the sole full image hash and records `RUNTIME_IMAGE_SHA256`; bound workers must not claim to compare the current image bytes to that digest. A changed/replaced image or image manifest is rejected through the versioned path, non-writable parent, size, manifest digest, and read-only checks; a changed source archive or source file fails through archive comparison.
4. Make `ecoda_runtime_export_csv PROFILE NV` export source, run ID, versioned runtime image/manifest paths, recorded SHA-256 values, sizes, and permissions alongside existing runtime fields. It must not initiate another full image hash; retain comma/newline rejection and exact-value comparison.
5. Make `ecoda_runtime_reexec_worker PROFILE SCRIPT` require `ECODA_SOURCE_SNAPSHOT_REQUIRED=1` and `ecoda_runtime_validate_bound_run()` for new runs. Host mode may execute only a script under the immutable source root and a recorded host environment prefix; it must reject a live-checkout script or mutable/nonmatching environment. Apptainer mode must bind the immutable source tree read-only and execute the snapshot script.
6. Use the relocated runtime layout for format-2 snapshot-backed jobs. Reject format-2 path-preserving execution when the bound source root differs from the build-time project root; retain path-preserving only for format-1 compatibility.

Update `src/slurm_config.sh` so snapshot execution keeps source and operational paths separate:

- Preserve `PROJECT_ROOT` as the directory containing the executing snapshot script.
- Honor absolute `ECODA_HOST_ENV_PREFIX`; it is exactly the `${ENV_ROOT}/.pixi/envs/py-cuda13` directory, not the project root. Validate `bin/python` and `bin/Rscript`, derive `ECODA_HOST_PYTHON_BIN`, `ECODA_HOST_PIXI_RSCRIPT`, `PYTHON_BIN`, `PIXI_RSCRIPT`, `PATH`, `LD_LIBRARY_PATH`, and `RETICULATE_PYTHON` from that same prefix, and compare the recorded binary digests. Snapshot-backed full-cohort workers use the relocated Apptainer environment; host-mode tests/legacy runs require an explicitly pinned, non-mutating prefix.
- Honor absolute `ECODA_LOGS_DIR` for writable run logs instead of placing logs inside read-only source. `ecoda_source_snapshot.sh exec` must export it; it must never pre-create the stage run root.
- Honor absolute `ECODA_AUX_ROOT` and set `SCGATE_DB_PATH` from the frozen snapshot `PROJECT_ROOT/aux` directory. Preserve the existing `SCGATE_DB_BRANCH` and auxiliary-file paths from the snapshot.
- Export `PYTHONDONTWRITEBYTECODE=1` for snapshot workers so Python cannot create `__pycache__` files in the read-only tree.
- Preserve inherited `HPC_SCRATCH_DIR`, NAS paths, and container-side interpreter overrides.

Update `src/utils/py/gene_utils.py:_load_ensembl105_map` so it first uses `Path(os.environ["ECODA_AUX_ROOT"]) / "EnsemblGenes105_Hsa_GRCh38.p13.txt.gz"` when the variable is set, validates that the file is readable, and otherwise retains the current `project_root / "aux"` fallback for legacy execution. Import `os` explicitly. This keeps both host and container code independent of the caller’s source-tree layout.

Update `ecoda_runtime_build_bind_args` so format-2 execution binds the immutable source tree, including its tracked `aux/`, read-only and binds run-specific logs/temp/scratch/reference paths with their existing modes. The image path and manifest path always come from the run identity, never from the mutable canonical default after submission.

Stage 4 must never write the frozen `aux/scGateDB.rds`. In `src/4_cell_type_annotation/1_submit_onboarding_stage.sh:749–776`, snapshot mode validates the pre-existing snapshot file with the already implemented `2.0_create_scgate_db.R --validate-only` path and does not submit the `--wrap` writer or pass `--force`. If validation fails, abort before scheduler submission. Do not change `2.0_create_scgate_db.R`; its existing `--validate-only` interface is the required invariant. Existing prevalidated `scGateDB.rds`, CellOntology cache, Ensembl map, and gene blocklist are immutable inputs for the run.

### 4. Wire snapshots through every submitter, worker, watchdog, retry, merge, and recovery path

Create the commit-keyed source snapshot and versioned runtime image before durable-gate `prepare`. The durable wrapper must invoke the executor from the snapshot itself, not from the mutable canonical checkout. Keep `--remote-workdir` passed to durable-gate `prepare` at `${BAMBOO_HOME}/ECODA_paper` so the existing profile path constraint passes; the exact command executes an absolute snapshot path from that work directory.

Use this exact wrapper shape, changing only fully resolved paths, stage script, and explicitly scoped arguments:

```bash
/bin/bash \
  "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_source_snapshots/<SOURCE_COMMIT>/tree/src/utils/bash/ecoda_source_snapshot.sh" exec \
  --source-root "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_source_snapshots/<SOURCE_COMMIT>/tree" \
  --source-manifest "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_source_snapshots/<SOURCE_COMMIT>/identity/source.manifest" \
  --host-env-prefix "/home/users/h/halterc/ECODA_paper/.pixi/envs/py-cuda13" \
  --runtime-image "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_runtime/<RUNTIME_ID>/ecoda-py-cuda13.sif" \
  --runtime-manifest "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_runtime/<RUNTIME_ID>/ecoda-py-cuda13.sif.manifest" \
  --run-id "<RUN_ID>" \
  --scratch-root "/home/users/h/halterc/scratch/ECODA_paper" \
  --logs-root "/home/users/h/halterc/scratch/ECODA_paper/_ecoda_logs/<RUN_ID>" \
  --script "src/<stage>/<submitter>.sh" -- <explicit-selection-arguments>
```

The implementation must not hard-code `/home/users/h/halterc`; the operator resolves remote `$HOME` before constructing the exact command. The durable manifest stores the fully resolved wrapper, remote work directory, snapshot, source manifest, runtime image, runtime manifest, run ID, selection, and output scope. Because the bootstrap executor is inside the immutable snapshot, a later pull cannot change it between gate preparation and execution.

After each stage opens/creates its run root, copy `source.manifest` and `runtime.identity` into `${ECODA_RUN_ROOT}/manifests/` atomically and record the same source/runtime/dependency fields in run metadata. For `--sync-only` and `--reuse-run`, invoke the exact recorded source snapshot and runtime identity. If a legacy run lacks either manifest, exit with `legacy_source_unpinned` before any worker, retry, or scheduler submission; only a separately invoked validator-only command may inspect its existing artifacts.

The executor may create only the separately named `${HPC_SCRATCH_DIR}/_ecoda_logs/<RUN_ID>` directory. It must not create `${ECODA_RUNS_ROOT}/${RUN_ID}` or `${ECODA_RUNS_ROOT}/${RUN_ID}/logs`. The new-run stage submitter then calls `ecoda_init_run STAGE RUN_ID`, which requires the run root to be absent and creates `manifests/`, `status/`, and `logs/` itself. `--sync-only`/`--reuse-run` calls `ecoda_open_run` on an existing exact run. A pre-existing run root in a new-run path is a hard conflict, not an adoption case.

Keep top-level `ecoda_runtime_validate_submission` only at the first submission boundary for each new run. Replace later full runtime checks with `ecoda_runtime_validate_bound_run()` and reuse the exported identity. Migrate these concrete call sites:

- Stage 2: `src/2_dataset_specific_preprocessing/1_submit_hpc.sh:390–393`, `stage2_watchdog.sh:102–106`, and worker wrappers `1.1_submit_gongsharma.sh:49`, `1.2_submit_combinedpbmc.sh:22`, `1.3_submit_joanito.sh:22`, `1.4_submit_kfoury_lowres_ct.sh:22`, `1.5_submit_myocardial.sh:31`, and `1.6_submit_bassez.sh:20`.
- Stage 3: `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh:99–106`, `1.2_preprocess_watchdog.sh:84–88`, `1.1_run_worker.sh:18`, and `src/utils/bash/h5ad_preflight_worker.sh:32`.
- Stage 4: `src/4_cell_type_annotation/1_submit_onboarding_stage.sh:780–785`, `1.2_annotation_watchdog.sh:37–41`, `1.3_prepare_chunks_watchdog.sh:43–47`, `3.3_merge_watchdog.sh:34–38`, `1.2_prepare_chunks_worker.sh:18`, `2.1_run_worker.sh:18`, and `3.2_merge_worker.sh:18`.
- Stage 5: `src/5_run_benchmark_methods/1_submit_hpc_array.sh:111–118` and its per-method `ecoda_runtime_export_csv` call at `:1289`, `matrix_watchdog.sh:88–94`, `watchdog_main.sh:94–101`, Python/R/trans-zero-imputation worker wrappers at their `1.1_run_worker.sh:18` sites, and `src/utils/bash/r_environment_preflight_worker.sh:32`.

All scheduler script paths must be checked by a shared `ecoda_require_source_script_path CANDIDATE SOURCE_ROOT` helper immediately before every scheduler boundary. Add that check to `src/utils/bash/h5ad_preflight_submit.sh` before its `sbatch` at `:28`, `src/5_run_benchmark_methods/benchmark_submit_common.sh` before watchdog submission and forwarded worker retry paths around `:693–705`, `src/5_run_benchmark_methods/matrix_watchdog.sh` before retry `sbatch` at `:147–149`, `src/5_run_benchmark_methods/watchdog_main.sh` before retry `sbatch` at `:110–120`, and every Stage 4 `sbatch` including the scGate `--wrap` at `src/4_cell_type_annotation/1_submit_onboarding_stage.sh:762–765, 792–810, 866–884, 918–936`. The helper canonicalizes both paths, rejects symlink escapes and non-snapshot paths, and is called again for every retry submission.

Retain Stage 2’s GongSharma-cap `afterok` edge until the CombinedPBMC intermediate is uniquely versioned. It protects a real in-place overwrite dependency and is unrelated to repository commit pinning.

Add `src/utils/bash/ecoda_run_audit.sh` with the exact CLI `ecoda_run_audit.sh --run-root ABSOLUTE_RUN_ROOT --stage STAGE --selection ABSOLUTE_SELECTION --source-manifest ABSOLUTE_SOURCE_MANIFEST --runtime-identity ABSOLUTE_RUNTIME_IDENTITY`. It validates only the specified run: run ID/stage, source/runtime manifest copies and digests, selection checksum and row count, scheduler-ID manifest, terminal status, artifact records, owners, and selected artifact contracts. It must never glob all run roots or submit/repair workers. The durable completion task invokes this validator on the recorded run before the one terminal `inspect`.

For a new snapshot-backed Stage 4 run, `2.0_create_scgate_db.R` is invoked only with `--validate-only`; no scGate writer job is submitted and no file under the frozen snapshot `aux/` is replaced.

### Deployment boundary for legacy running jobs

Before changing the canonical Bamboo checkout, identify active format-1 jobs and keep the exact checkout/runtime source they already use available until terminal audit. Do not pull the new implementation over that checkout while workers or watchdogs can re-enter it. Develop the new implementation in a separate checkout, create a format-2 source/runtime snapshot from that checkout, and launch the first snapshot-backed run from that wrapper. After legacy jobs finish, the canonical checkout may advance normally; old gate manifests remain readable evidence but are never used as rerun triggers.

### 5. Remove current-checkout and stale-global gate coupling

Update `.agents/skills/durable-hpc-gate-ecoda/references/profile.json` for snapshot-backed runs:

- Remove `datasets-pixi-lock-sha256` and `repository-head` from `immutable_fingerprints`; they fingerprint the mutable canonical checkout and cause false terminal staleness.
- Delete `terminal-run-owned-manifests-status` and `terminal-emitted-scheduler-status-records` from `policy.audit_commands`.
- Delete `stage-neutral-run-root-contract` and `stage-neutral-scheduler-status-contract` from `policy.artifact_contracts`.
- Do not replace any of these four entries with another broad `${HPC_SCRATCH_DIR}/_ecoda_runs/*` scan. The exact `ecoda_run_audit.sh` invocation, current run-owned terminal status, source/runtime manifests, exact wrapper identity, and durable gate’s exact scheduler-ID accounting are authoritative.
- Retain canonical-root existence, configured scratch/NAS paths, required wrapper presence, exact command digest, one terminal accounting query, reviewer approval, and fail-closed lifecycle rules.
- Keep same-group mutual exclusion. Distinct serialization groups may launch concurrently only after concrete artifact disjointness is proven.
- After removing both live-checkout fingerprints, persist `"immutable_fingerprints": []`; the durable-gate schema and parser must accept an explicit empty list. The run-bound source/runtime manifests and `ecoda_run_audit.sh` provide the identity checks instead.

Update `.agents/skills/durable-hpc-gate-ecoda/SKILL.md` sections `Enforce Bamboo and repository invariants`, `Coordinate benchmark waves and reviewed lineage`, `Profile severity and audit contracts`, and `Apply fail-closed recovery` so they describe run-bound source/runtime/aux manifests instead of current-checkout HEAD/pixi fingerprints and do not require global run-root scans. Keep exact scheduler-ID accounting, one waiter, reviewer approval, no-compute repair, and fail-closed lifecycle behavior.

The source snapshot and stage runtime manifests are provenance authorities. Old terminal gate manifests remain immutable evidence and are never deleted or used as automatic rerun instructions.

Change no-op behavior explicitly: first perform a validator-only selection preflight. If every requested artifact row is valid, write `NOOP_VALIDATED` to the specified run report and do not call durable-gate `prepare/launch`; the ECODA profile retains `require_scheduler_ids=true` for actual compute gates. If at least one row is missing/invalid, submit only those rows and send that run through the durable gate.

Update submitter/recovery behavior so:

- existing artifact reuse is decided only by the selected artifact’s nonempty/schema/identity/checksum contract and explicit `--force` scope;
- old `FAILED`/`PRELAUNCH_STOP` gates cannot create a new scheduler selection;
- `--sync-only` and annotation reuse require the exact requested run ID and run-owned manifests;
- legacy runs without source/runtime manifests cannot submit compute or retries;
- validated baseline artifacts remain outside every recomputation selection;
- status messages say `baseline artifact valid; no rerun selected` rather than treating a stale gate as evidence that computation is missing.

Do not manually mark old gates complete, delete their evidence, or broad-force pipelines 1–5. Repair/validation remains no-compute by default.

### 6. Reduce redundant full-file checks without weakening artifact safety

Use the existing `ecoda_validate_checksum_record()` in `src/utils/bash/ecoda_run_common.sh:336–358` whenever a prior strict digest/size record is available. It checks sidecar fields and current size without calculating another MD5. Keep `ecoda_validate_checksum()` for first acceptance of an artifact or an untrusted boundary.

Add shared artifact records in `src/utils/bash/ecoda_run_common.sh` with exact functions `ecoda_artifact_record_path PATH RUN_ID`, `ecoda_write_artifact_record PATH PRODUCER RUN_ID`, and `ecoda_validate_artifact_record PATH PRODUCER RUN_ID`. `ecoda_artifact_record_path` must derive the producer run root as `${ECODA_RUNS_ROOT}/${RUN_ID}` from the supplied RUN_ID, not from the caller’s current `ECODA_RUN_ROOT`; publication callers pass their own run ID and downstream consumers pass the upstream producer run ID. It writes under `<producer-run-root>/manifests/artifacts/<first-32-hex-of-SHA256(canonical-path)>.record`; the record filename is bounded and the record `PATH` field retains the full canonical path. Each atomic record contains exactly `PATH`, `SIZE`, `MD5`, `RUN_ID`, `PRODUCER`, and `STATE=PUBLISHED`. `ecoda_write_artifact_record` may run only after a full checksum and semantic contract pass. `ecoda_validate_artifact_record` validates the record/run/producer binding and calls `ecoda_validate_checksum_record`; it must never be used for a newly produced artifact before its first strict checksum. Stage-specific status files may reference this record but must not invent another digest format.

Retain full checksum verification at these boundaries:

- artifact creation/publication before declaring a producer successful;
- H5AD source acceptance and source-identity creation when the source sidecar has not already been strictly validated;
- scratch-to-NAS synchronization and remote comparison;
- deserialization boundaries that must prove the sidecar/content pair before loading.

Collapse only duplicate rehashes after the artifact is run-owned and its artifact record is published:

- Stage 3 preflight/worker/watchdog H5AD checks, using the preflight artifact record downstream;
- Stage 4 source/union/Feather/merge checks, using artifact records plus the existing merge source/union relationship;
- Stage 5 source preflight before source identity, then `h5ad_source_identity.py --validated-source-sidecars`; retain ordered `Sample`-ID validation;
- benchmark cache/watchdog/final validators for non-deserializing checks, using recorded digest/size values.

Stage 5 must reorder its current `stage5_prepare_source_identity`/`stage5_compute_h5ad_preflight` sequence (`1_submit_hpc_array.sh:785–788`): first call `ecoda_require_input_ownership` only to prove that no run, including the current run, owns the input as an ACTIVE writer; then run strict compute-node H5AD preflight and publish source artifact records; then create/verify source identity with `--validated-source-sidecars`. After preflight, call `ecoda_validate_artifact_record` and the semantic H5AD contract before trusting the source identity. A source identity must never become trusted before the first strict preflight, and the preflight must fail on any same-size mutation.

`ecoda_require_input_ownership PATH RUN_ID` must reject an ACTIVE writer belonging to any run, including `RUN_ID`; readers never bypass an active writer. It may allow a missing record before preflight. Add `ecoda_validate_input_artifact PATH PRODUCER RUN_ID` for the later record/schema/identity requirement after the writer is terminal and the source preflight has published its record. This separation prevents a preflight/ownership circular dependency; only `ecoda_validate_output_ownership` has same-run re-entrancy for OOM retries.

Use the one-full-hash-per-boundary rule as follows:

- newly produced H5AD/Feather/RDS artifacts: full hash plus semantic validation at publication, then artifact record;
- Stage 3/4/5 watchdogs and matrix validators: record/size/path checks plus semantic checks, unless they deserialize an artifact;
- `src/5_run_benchmark_methods/matrix_artifact_validator.py`: full sidecar/content verification before any RDS load, record-only checks for non-deserializing artifacts;
- `src/5_run_benchmark_methods/validate_benchmark_rds_contract.R`: preserve strict checksum-before-load behavior for RDS artifacts;
- `src/5_run_benchmark_methods/benchmark_pipeline.R:278–301`: verify a checksum sidecar before `readRDS` whenever the file is listed. The current legacy no-sidecar path must return explicit `legacy_unverified` and skip without calling `readRDS`; it must not be silently accepted;
- `src/5_run_benchmark_methods/benchmark_submit_common.sh:1080` and `:1168`: retain strict remote `md5sum -c` checks after synchronization;
- `src/5_run_benchmark_methods/benchmark_hpc_utils.R:artifact_checksum_ok`, `benchmark_methods_r.R`, `run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R`, `run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R`, `consolidate_gloscope_results.R`, the Python benchmark worker cache, `matrix_artifact_validator.py`, and execution-time merge paths: migrate every cache/read caller to the artifact-record interface, but keep full checksum verification immediately before RDS/Feather deserialization when no earlier verified record is carried into that same read boundary;
- `src/utils/py/h5ad_source_identity.py:202–233`, `:277–304`, and `:365–389`: retain source size/MD5/sample-order identity and use `validated_sidecars` only after strict preflight;
- retain literal reference-map MD5 verification in `src/4_cell_type_annotation/1.0_stage_reference_maps.sh` at acquisition/publication.

Do not replace MD5 with a weaker existence-only test. If an artifact can still be modified by another run, rehash it before reuse; record-only validation is allowed only after output ownership, atomic publication, and the immutable run record make the earlier digest trustworthy.

### 7. Prevent overlapping writers and protect cross-stage readers

Add these shared functions to `src/utils/bash/ecoda_run_common.sh`:

- `ecoda_artifact_owner_key PATH` — canonical absolute final artifact path, never a selection-file digest. Use a readable basename prefix only when the resulting key is at most 100 bytes; otherwise use `artifact_<first-32-hex-of-SHA256(canonical-path)>`. This keeps every owner filename below filesystem component limits and preserves the exact canonical path inside `owner`.
- `ecoda_validate_output_ownership STAGE SELECTION_MANIFEST RUN_ID` — expands every selected dataset/view/method/analysis/parameter row into exact scratch and NAS output paths, rejects duplicate paths in the selection, and rejects any path with an ACTIVE global artifact writer belonging to another run. If the active owner belongs to `RUN_ID`, allow re-validation for that run’s OOM retry without creating a second owner or changing ownership state.
- `ecoda_require_input_ownership PATH RUN_ID` — rejects an ACTIVE writer for the input path regardless of owner run; readers never bypass an active writer and may proceed with a missing record only to let strict preflight create it.
- `ecoda_validate_input_artifact PATH PRODUCER RUN_ID` — after preflight or when consuming an upstream-published artifact, requires the published artifact record plus schema/identity contract.

Store cross-stage writer owners under `${HPC_SCRATCH_DIR}/_ecoda_owners/artifact/<safe-canonical-path-key>/owner`, with `RUN_ID`, stage, path, state, and PID. Preserve existing stage-scoped owner files for compatibility, but acquire/release the global artifact owner for every canonical writer. A Stage 5 reader of a Stage 3/4 H5AD must prove that the writer owner is terminal `OK` before trusting a record-only checksum. Different selection digests must not bypass path intersection.

Call `ecoda_validate_output_ownership` before every stage submits any worker and before every OOM retry array. The same-run ACTIVE-owner exception above is mandatory for retries; another run’s ACTIVE owner remains a hard failure. Call `ecoda_require_input_ownership` before Stage 4 annotation preparation and before Stage 5 source preflight reads a canonical H5AD; it must block any ACTIVE writer, including the current run. Call `ecoda_validate_input_artifact` only after the source preflight publishes its record or when an upstream stage already published one. An overlap fails before scheduler submission, records the concrete path and existing owner, and does not reclaim or overwrite the owner.

Use distinct durable `serialization_group` values only after this exact path set is disjoint. A job writing the same canonical H5AD, Feather, RDS, execution-time log, or NAS path remains in the canonical group and is rejected by the owner check. New methods use unique method-specific output keys. Keep the Stage 2 CombinedPBMC dependency until its write path is versioned.

Leave `src/1_stage_data/1_stage_data.sh` unchanged. It is the intentional serial NAS-to-scratch staging step and has no repository-SHA, MD5 sidecar, or durable-gate bloat.

### 8. Align exact batch selection and post-baseline method selection

Update the hard-coded exact batch dataset identity from `Kidney_KPMP` to the current `datasets.json` key `Kidney_KPMP_full` in `src/utils/bash/ecoda_run_common.sh:512–513`, `src/5_run_benchmark_methods/matrix_artifact_validator.py:28–32`, and `src/5_run_benchmark_methods/validate_benchmark_rds_contract.R:495–499`. Update their exact-batch fixtures in `tests/test_ecoda_run_common.sh`, `tests/test_preprocessing_stage_submitter.sh`, `tests/test_annotation_stage_submitter.sh`, and `tests/test_benchmark_matrix_submitter.sh`; update other batch fixtures only where they represent this same pipeline selection contract. Do not add a compatibility alias for the old dataset key.

In `src/5_run_benchmark_methods/1_submit_hpc_array.sh:789–805`, define a separate `BASELINE_METHODS` list with the ten baseline ordinary methods and use it when `--methods` is omitted. Keep `BATCH_EFFECT_METHODS` and the fixed `EXPECTED_BATCH_METHODS` suite for the current batch pass. A method added to `method_spec` is not added to either default list. Require `--methods` or an explicit selection manifest for that method. Preserve `--force` as an explicit row-scoped recomputation flag; never use it merely because a historical gate is stale.

Document the post-`datasets.json` batch correction state in the baseline section: `Lupus_PBMC` and `Parkinson` have changed `cell_type_high_res`, so only methods whose feature/input path consumes that field are pending targeted rerun; unaffected batch rows remain reusable. The new `Kidney_KPMP_full` rows are pending targeted validation/rerun. Build the explicit Stage 5 selection from the actual method/input dependency paths and record those rows in the run manifest; do not infer the affected method list from a stale gate or rerun every batch method.

### 9. Add focused regression coverage

Extend `tests/test_runtime_container.sh` to cover format-2 identity and replacement behavior:

1. Build an image manifest with `IMAGE_BUILD_GIT_REVISION=A` and a source manifest with `SOURCE_COMMIT=B`, while keeping identical dependency digests; assert submission validation succeeds without live `git HEAD` equality.
2. Change the source lock/config digest; assert validation fails before scheduler execution.
3. Change a source-tree file, retained archive, or source manifest after creation; assert archive/tree or manifest validation fails.
4. Replace the versioned image or image manifest after initial validation; assert `ecoda_runtime_validate_bound_run()` rejects the changed path/size/manifest digest/read-only contract before worker execution. Do not use device IDs in the fixture.
5. Assert relocated Apptainer binds the snapshot source/aux and separate run-log paths and executes the snapshot script.
6. Assert `PYTHONDONTWRITEBYTECODE=1` is propagated and a snapshot worker cannot create `__pycache__`.
7. Assert host-mode required-snapshot execution rejects a live-checkout script and a mutable/nonmatching host environment prefix.
8. Assert old `FORMAT=1` manifests retain their legacy validation path for compatibility.

Add `tests/test_source_snapshot.sh` for a temporary clean Git repository and auxiliary/runtime trees: commit A, create a complete commit-keyed snapshot A, commit B in the source repository, prove snapshot A still contains A, mutate a snapshot file/archive and assert execution fails, test dirty/untracked source rejection, test a commit-snapshot identity conflict, and test that the bootstrap executor runs from snapshot A rather than the mutable checkout. Assert that `tree/aux` contains the three tracked auxiliary files and that `ECODA_AUX_ROOT` resolves them.

Add `tests/test_gene_utils_aux_path.py` with a temporary source tree and auxiliary directory: set `ECODA_AUX_ROOT` to the auxiliary directory, load the Ensembl map through `_load_ensembl105_map`, and assert the environment-selected file is used; unset the variable and assert the legacy `project_root/aux` fallback remains available.

Extend `tests/test_checksum_reuse.sh` to assert that a strict `ecoda_validate_checksum()` followed by repeated `ecoda_validate_checksum_record()`/artifact-record calls invokes the hash command only once, while changed artifact size, sidecar digest, record producer, or run ID fails. Keep the existing remote comparison assertion that hashes the remote boundary once.

Extend `tests/test_stage2_submitter.sh`, `tests/test_stage2_watchdog.sh`, `tests/test_preprocessing_stage_submitter.sh`, `tests/test_annotation_stage_submitter.sh`, `tests/test_benchmark_matrix_submitter.sh`, `tests/test_benchmark_matrix_watchdog.sh`, `tests/test_h5ad_preflight.sh`, and `tests/test_benchmark_sync.sh` so captured scheduler commands contain snapshot source/runtime paths and run manifests. Assert shared preflight/watchdog/retry helpers reject a script outside the source root, and assert an ACTIVE owner belonging to the same run is accepted for an OOM retry while another run is rejected.

Extend `tests/test_annotation_stage_submitter.sh` to assert that snapshot mode uses the existing `2.0_create_scgate_db.R --validate-only` path, never submits the scGate `--wrap` writer or passes `--force` to the frozen aux tree, and rejects a missing/invalid frozen `aux/scGateDB.rds` before scheduler submission. Assert the source snapshot’s Ensembl map is used by the preprocessing fixture.

Update `tests/test_durable_profile_stage_neutral.sh` to assert that both live repository immutable fingerprints and all four broad run-glob contracts are absent, that an empty `immutable_fingerprints` array is accepted, and that canonical-root, exact-command, scheduler/accounting, and reviewer contracts remain. Keep `test_durable_hpc_gate_parallelism.py` for same-group blocking and distinct-group parallelism, and add two different-selection fixtures whose expanded artifact paths overlap; assert the path-level ownership check rejects them before `sbatch`.

Extend the existing Stage 5 submitter test with an existing valid baseline method row plus an absent new method row: assert only the new method is submitted, no `--force` is propagated to the valid row, and a stale/failed gate fixture does not add any row. Add fixtures for the changed `Lupus_PBMC`/`Parkinson` `cell_type_high_res` dependency selection and the new `Kidney_KPMP_full` key; assert unaffected batch rows remain outside the targeted rerun.

Extend `tests/test_benchmark_rds_contract.R` and a direct `benchmark_pipeline.R` loader fixture with a missing-sidecar case that reports `legacy_unverified` and proves `readRDS` is not called; keep valid sidecar-before-load and malformed-sidecar failures.

### 10. Verify the migration on isolated `_debug` roots only

Run focused tests through the repository Pixi environment, not a system Python/R:

- `bash tests/test_runtime_container.sh`
- `bash tests/test_source_snapshot.sh`
- `bash tests/test_checksum_reuse.sh`
- `bash tests/test_stage2_submitter.sh`
- `bash tests/test_stage2_watchdog.sh`
- `bash tests/test_h5ad_preflight.sh`
- `bash tests/test_preprocessing_stage_submitter.sh`
- `bash tests/test_annotation_stage_submitter.sh`
- `bash tests/test_benchmark_matrix_submitter.sh`
- `bash tests/test_benchmark_matrix_watchdog.sh`
- `bash tests/test_benchmark_sync.sh`
- `bash tests/test_durable_profile_stage_neutral.sh`
- `pixi run -e default python tests/test_gene_utils_aux_path.py`
- `pixi run -e default python tests/test_durable_hpc_gate_parallelism.py`
- `pixi run -e default Rscript --vanilla tests/test_benchmark_rds_contract.R`

Run the end-to-end source-isolation smoke test with two temporary scratch roots populated only with the existing Joanito five-sample `_debug` fixture, so the test does not touch canonical baseline artifacts:

1. Use the first committed format-2 implementation revision `IMPLEMENTATION_COMMIT_A` (not the documentation baseline commit `5302671ad94556edcf9acccf372d2dc34121d714`, which cannot contain the new executor) to create source/runtime snapshot A and launch an actual Stage 3/4 `_debug` wrapper against scratch root A.
2. Commit a harmless new method/script fixture as `IMPLEMENTATION_COMMIT_B`, create source/runtime snapshot B against scratch root B, and launch the explicitly scoped Stage 5 new-method wrapper in a distinct serialization group while A remains active.
3. Keep scratch roots, run IDs, global artifact owners, and output namespaces distinct. Do not run both jobs against canonical `${HPC_SCRATCH_DIR}/_debug/output` paths.
4. Confirm A’s delayed worker/watchdog/retry reports source/runtime A, B reports source/runtime B, changing the development checkout does not fail either gate, and no output/input owner is crossed.
5. Replace or attempt to rebuild A’s versioned image path; confirm A rejects the replacement before retry.
6. Confirm a valid baseline artifact is skipped without `--force`, the absent new method row is submitted without touching successful rows, an old stale/failed gate does not create a scheduler row, and missing/malformed checksum, source, runtime, schema, or ownership contracts fail closed.

Run validator-only checks for already completed production benchmark artifacts and unaffected batch-effect rows. Explicitly target only the `Lupus_PBMC`/`Parkinson` downstream methods that consume changed `cell_type_high_res` and the new `Kidney_KPMP_full` batch rows when those reruns are approved. Do not launch full cohorts for this migration test and do not broad-recompute the completed benchmark baseline.

## Critical files & anchors

- `src/utils/bash/ecoda_runtime.sh:116–141, 226–320, 337–454` — separate image/source/dependency validation, run-bound identity, source-aware binds, export, and worker re-exec.
- `src/utils/bash/build_ecoda_runtime.sh:110–236` — versioned format-2 image identity and dependency checks during image construction.
- `src/slurm_config.sh:10–67, 102–177` — exact host environment prefix, source/aux/log paths, and container interpreter handling.
- `src/utils/bash/ecoda_run_common.sh:171–207, 275–408, 410–453` — run metadata, artifact records, ownership, and checksum-record reuse.
- `AGENTS.md` and `.agents/skills/durable-hpc-gate-ecoda/references/profile.json` — completed baseline, pending-method rule, snapshot operator workflow, and removal of live/global gate coupling.

## Verification

The migration is complete only when focused tests pass and isolated `_debug` runs demonstrate simultaneous source/runtime A and source/runtime B execution while the development checkout changes. The observable proof must include self-verifying source archives containing the sole frozen `tree/aux`, run-bound nonreplaceable image identity, source/runtime identities in both run manifests, exact source paths at every scheduler/retry boundary, distinct ownership, no terminal failure from the main checkout advancing, no-op handling without scheduler-ID errors, and fail-closed rejection for changed dependencies, changed/replaced images, mutated snapshots, malformed artifacts, missing required checksums, invalid schema, legacy unverified RDS, and overlapping output ownership. Completed benchmark artifacts and unaffected batch rows must remain untouched and outside recomputation selections; changed high-resolution consumers and new `Kidney_KPMP_full` rows must be explicitly represented as pending targeted work.

## Assumptions & contingencies

- The user-confirmed status is authoritative even if existing gate manifests disagree: the benchmark baseline is complete, while `batch_effect_uncorrected` is not terminally complete after the recent `datasets.json` update. The changed `Lupus_PBMC`/`Parkinson` `cell_type_high_res` consumers and new `Kidney_KPMP_full` rows are pending targeted work. If `datasets.json` flags change again before implementation, preserve the recorded baseline and record the later configuration separately; do not infer completion or trigger a blanket rerun from gates.
- The current commit anchor is `5302671ad94556edcf9acccf372d2dc34121d714`. If the repository advances before documentation, keep this as `baseline anchor` and record the later implementation HEAD separately.
- New snapshot-backed full-cohort jobs use versioned runtime images and Apptainer relocated layout by default. If the relocated image cannot run the isolated `_debug` worker, stop before full-cohort launch and fix the binding contract; do not fall back to the mutable path-preserving checkout.
- Old `FORMAT=1` images/runs remain readable and may finish under their existing contract. They are not converted in place; new runs use `FORMAT=2`, versioned runtime paths, self-verifying source/aux snapshots, and run-bound identity.
- Legacy RDS files without a checksum sidecar are not automatically trusted. Normal loading skips them as `legacy_unverified`; any later checksum registration is an explicit validator-only operation with no scheduler submission.
- Existing canonical output paths remain shared. If an intended new method cannot obtain a disjoint output namespace and concrete artifact owner key, keep it in the existing serialization group rather than bypassing the lock.
- No-op validator-only runs are not durable compute gates and therefore do not need scheduler IDs; every run that submits compute still uses the existing exact-ID terminal accounting and reviewer flow.
- The user selected commit-keyed source snapshots and the tracked `tree/aux` files as the sole frozen auxiliary root. Multiple runs may reuse an unchanged source snapshot; each run still records separate runtime identity, manifests, logs, owners, and outputs. A new commit creates a new snapshot and must not invalidate successful older jobs or artifacts. A requested dependent rerun remains explicit and targeted.
