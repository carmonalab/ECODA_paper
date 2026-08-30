# HPC Runtime Environment Containerization — Implementation Plan

Plan slug: `20260829211352-hpc-environment-containerization`

Repository destination: `.agents/plans/20260829211352-hpc-environment-containerization-plan.md`

## Context

Concurrent Bamboo workers currently import Python and R packages directly from `${PROJECT_ROOT}/.pixi/envs/py-cuda13` on BeeGFS. The handoff records reproducible 16-node failures in both direct-R and direct-Python tests even after Pixi activation was bypassed, while smaller fan-out and sequential reads passed. The working technical document identifies the likely failure mechanism as a BeeGFS metadata/small-file storm and selects an immutable Apptainer/SIF runtime as the primary mitigation.

The end state is one shared scientific script set for pipelines 2–5, with an opt-in `ECODA_RUNTIME_MODE=host|apptainer` boundary. Host mode keeps the current direct Pixi binaries and guarded environment mutation. Apptainer mode enters the image once per scheduled scientific worker, then exposes ordinary in-image `PYTHON_BIN` and `PIXI_RSCRIPT` paths so R-to-Python calls remain direct and never nest Apptainer. Pipeline 1, local macOS analysis, `refresh_env.sh`, `setup_env_sbatch.sh`, `datasets.json`, `pixi.toml`, `pixi.lock`, and duplicate pipeline trees remain unchanged except where a phase explicitly names a shared runtime export callsite.

## Approach

### Phase 0 — Persist the plan and correct the recorded build example

**Parallel work unit 0A — `DocsOwner` (exclusive ownership of the documentation file).**

1. Before any source edit, read this canonical plan and write the identical content to `.agents/plans/20260829211352-hpc-environment-containerization-plan.md`. Do not overwrite or stage existing user WIP reported by `git status` (`notebooks/benchmark_analysis.rmd`, unrelated plan/review files, or `.gate/`).
2. In `docs/hpc_environment_containerization.md`, change only the candidate Pixitainer line under `## 5. Pixitainer assessment` from the equals form to the parser-compatible separate-argument form:

   ```bash
   --base-image rockylinux:9 \
   ```

   Keep `rockylinux:9` explicit; it is the production base-image version pin. Do not rewrite earlier findings.
3. Append a short `### Entry 3 — Pixitainer option syntax` below Entry 2 stating that the feasibility command passes `--base-image` and `rockylinux:9` as separate arguments and retains the explicit Rocky Linux 9 pin.

**Verification:** confirm the repository plan exists at the exact destination; inspect the documentation diff and verify only the command syntax plus the appended discussion entry changed. Confirm unrelated WIP remains untouched.

**Barrier B0:** the exact plan and corrected command are persisted before implementation files are modified.

### Phase 1 — Establish the immutable runtime and configuration contracts

**Parallel work unit 1A — `RuntimeFoundation` (exclusive ownership of the new runtime/build files and central config).**

Targets: `src/utils/bash/ecoda_runtime.sh`, `src/utils/bash/build_ecoda_runtime.sh`, `src/slurm_config.sh`.

1. Extend `src/slurm_config.sh` without changing the existing host values for `PYTHON_BIN`, `PIXI_RSCRIPT`, `PATH`, `LD_LIBRARY_PATH`, or `RETICULATE_PYTHON`:
   - Export `ECODA_RUNTIME_MODE` with default `host`.
   - Export `ECODA_RUNTIME_IMAGE` with default `${HPC_SCRATCH_DIR}/_ecoda_runtime/ecoda-py-cuda13.sif` and `ECODA_RUNTIME_MANIFEST` with default `${ECODA_RUNTIME_IMAGE}.manifest`.
   - Export `ECODA_RUNTIME_IN_CONTAINER` with default `0`, `ECODA_RUNTIME_PROFILE` with default `default`, and `ECODA_APPTAINER_NV` with default `0`.
   - In the `ECODA_RUNTIME_IN_CONTAINER=1` branch, use `ECODA_RUNTIME_PREFIX` supplied by the launcher/manifest and set direct in-image paths: `${prefix}/bin/python`, `${prefix}/bin/Rscript --vanilla`, `${prefix}/bin` first in `PATH`, `${prefix}/lib` first in `LD_LIBRARY_PATH`, `${prefix}/lib/R` as `R_HOME`, and `${prefix}/bin/python` as `RETICULATE_PYTHON`.
   - In that branch set `PYTHONNOUSERSITE=1` and unset `PYTHONHOME`, `PYTHONPATH`, `R_LIBS_USER`, `R_LIBS_SITE`, `R_ENVIRON_USER`, and `R_PROFILE_USER` so host package discovery cannot override the image. Keep the existing module loads and host branch unchanged.
   - The host branch is selected whenever `ECODA_RUNTIME_IN_CONTAINER` is not `1`; `ECODA_RUNTIME_MODE=apptainer` alone must never replace host `PYTHON_BIN` or `PIXI_RSCRIPT`. Login-side submitters/watchdogs validate the image with host tools, then only the worker re-exec boundary sets `ECODA_RUNTIME_IN_CONTAINER=1` and the in-image prefix.
   - Preserve the logical `HPC_SCRATCH_DIR` and `HOME_REF_DIR` values explicitly inside the image even when `--no-home`/containment changes `HOME`; use node-local `TMPDIR` for `NUMBA_DISABLE_CACHING`, `XDG_CACHE_HOME`, and `MPLCONFIGDIR` so read-only SIF layers receive no cache writes.
2. Implement `src/utils/bash/ecoda_runtime.sh`, sourced after `slurm_config.sh`, with Bash 3.2-compatible functions and no associative arrays:
   - `ecoda_runtime_validate_submission` accepts only `host` and `apptainer`. Host mode returns without requiring an image. Apptainer mode requires the configured image and key-value manifest to be readable, nonempty, regular files; validates the manifest SHA-256 against the SIF; invokes `apptainer inspect`; validates `RUNTIME_ENV=py-cuda13`, `BASE_IMAGE=rockylinux:9`, `PIXITAINER_VERSION=0.8.3`, a recognized bind layout, a nonempty `CONTAINER_ENV_PREFIX`, and the recorded lock/source/runtime identity. Missing, unreadable, corrupt, mismatched, or uninspectable images fail closed before any `sbatch` call.
   - `ecoda_runtime_export_csv <profile> <nv>` prints a comma-separated export fragment containing `ECODA_RUNTIME_MODE`, `ECODA_RUNTIME_IMAGE`, `ECODA_RUNTIME_MANIFEST`, `ECODA_RUNTIME_PROFILE`, and `ECODA_APPTAINER_NV`. Reject commas in values so `sbatch --export` cannot be ambiguous. Validate `<nv>` as `0` or `1`.
   - `ecoda_runtime_build_bind_args <profile>` populates the global `ECODA_RUNTIME_BIND_ARGS` array. Every profile binds the resolved project source, resolved real BeeGFS scratch to the original `HPC_SCRATCH_DIR` destination, `LOGS_DIR` read/write, and the actual node-local `${TMPDIR:-/tmp}` read/write. The `stage4` profile additionally binds `HOME_REF_DIR` read-only. `relocated` layout binds the project root read-only; `path-preserving` layout does not bind the project root because that would mask the embedded environment, and instead binds `${PROJECT_ROOT}/src`, `datasets.json`, `config_helper.R`, `aux`, and `LOGS_DIR` at their original absolute destinations. Resolve symlinked source paths with `realpath -e`; fail for missing required sources. Reject any destination that overlaps `CONTAINER_ENV_PREFIX`.
   - `ecoda_runtime_reexec_worker <profile> <absolute-script> [args...]` validates the mode, script, image, manifest, and bind profile; returns without action in host mode or when `ECODA_RUNTIME_IN_CONTAINER=1`; otherwise `exec`s exactly one `apptainer exec` with `--containall`, `--no-home`, `--no-mount home,cwd,hostfs,bind-paths`, the explicit bind array, `--nv` only for `ECODA_APPTAINER_NV=1`, `--env ECODA_RUNTIME_IN_CONTAINER=1`, the manifest-provided `ECODA_RUNTIME_PREFIX`, the image, `/bin/bash`, the source script path, and its arguments. Do not put a multiword Apptainer command in `PYTHON_BIN` or `PIXI_RSCRIPT`, and do not use `--cleanenv`; dynamic worker/Slurm variables are inherited while the config branch sanitizes runtime-sensitive host variables. An Apptainer startup failure replaces the worker process and is terminal, not a transient worker retry.
3. Implement `src/utils/bash/build_ecoda_runtime.sh` as the only checked-in builder:
   - CLI: `build_ecoda_runtime.sh --layout relocated|path-preserving --output ABSOLUTE_SIF [--force]`. Reject relative output, unknown layout, missing `pixi.lock`, missing realized `.pixi/envs/py-cuda13`, active environment locks, failed Slurm queries, and active jobs other than the current build allocation. Require execution on a Bamboo compute allocation, not a login node or macOS.
   - Pin `PIXITAINER_VERSION=0.8.3` and `BASE_IMAGE=rockylinux:9`. Derive the exact host Pixi version from `pixi -V` and pass that value through `--pixi-version`; record it rather than silently selecting a newer image toolchain. Use the separate Pixitainer syntax `--base-image rockylinux:9`, never `--base-image=rockylinux:9`.
   - Run a quiet dry-run first and retain its generated definition as `${IMAGE}.dryrun.def`. Require the definition to show the requested Rocky base, `--manual` direct shell entrypoint, no fresh `py-cuda13` installation, and an explicit `--add-file` of the realized environment. Build with `--manual --no-install --env py-cuda13 --base-image rockylinux:9 --pixi-version <pinned-version> --add-file <source>:<layout-destination> --keep-def --output <image>` using `pixitainer=0.8.3`; do not invoke Pixitainer or Pixi in production workers.
     The pinned build also installs the minimal Rocky utilities required by base R/Python startup (`dnf install -y which jq`) through a fixed post-command; this is not a Pixi environment installation.
   - Set `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` to approved node-local/temporary storage, build atomically, run `apptainer inspect`, compute SHA-256, and write `${IMAGE}.manifest` atomically. The manifest format is strict `KEY=VALUE` with `FORMAT=1`, `IMAGE_PATH`, `IMAGE_SHA256`, `RUNTIME_ENV=py-cuda13`, `RUNTIME_LAYOUT`, `CONTAINER_ENV_PREFIX`, `CONTAINER_PROJECT_ROOT` when path-preserving, `BASE_IMAGE=rockylinux:9`, `PIXITAINER_VERSION=0.8.3`, `PIXI_VERSION`, `APPTAINER_VERSION`, `GIT_REVISION`, and `PIXI_LOCK_SHA256`. Never leave a success manifest for a failed or partial SIF.
     Source-artifact fingerprints record the logical absolute path plus physical MD5/size, so a `$HOME/scratch` symlink namespace is stable between container workers and host-side watchdog validators.
   - `relocated` copies the realized environment to `/opt/ecoda/py-cuda13`. If its relocation smoke fails, build the predetermined `path-preserving` fallback using the same pinned base/tool versions and copy destination `${PROJECT_ROOT}/.pixi/envs/py-cuda13`; its manifest records the exact build-time project root. Do not substitute multiple BeeGFS environment directories, startup delays, or an unvalidated tar extraction path.

**Parallel work unit 1B — `RuntimeContractTests` (exclusive ownership of the new runtime test and host preflight regression).**

Targets: `tests/test_runtime_container.sh`, `tests/test_r_environment_preflight.sh`.

1. Add `tests/test_runtime_container.sh` using temporary roots and fake `apptainer`, `python`, `Rscript`, and scheduler binaries. Cover host mode direct paths; successful relocated and path-preserving image manifests; missing/unreadable/truncated image; SHA-256 mismatch; failed `apptainer inspect`; missing bind source; project/runtime bind collision; symlinked scratch resolution; read-only versus read/write bind arguments; exactly one container boundary; internal `PYTHON_BIN`, `PIXI_RSCRIPT`, `R_HOME`, `RETICULATE_PYTHON`, and sanitized import paths; and `ECODA_APPTAINER_NV=1` versus `0`. Assert all invalid cases return nonzero, print the image/check name, and leave no worker marker or scheduler submission.
2. Update `tests/test_r_environment_preflight.sh` to assert the direct host contract specifically under `ECODA_RUNTIME_MODE=host`; retain the exact host `Rscript --vanilla` path assertion and prove no `pixi run` is used. Add only the configuration-level container branch assertions that do not require a real SIF.

**Verification:** run `bash -n` on both new scripts and changed config; run `bash tests/test_runtime_container.sh` and `bash tests/test_r_environment_preflight.sh`. Expected result: all fake launch cases pass, host paths remain exact, and no fake worker runs after invalid image/bind input.

**Barrier B1:** runtime helper/config tests are green before any pipeline submitter or worker is changed.

### Phase 2 — Bamboo image feasibility and immutable-prefix A/B gate

**Work unit 2A — `RuntimeFeasibility` (no production caller edits).**

Use a Bamboo compute allocation and the approved existing host mutation lifecycle. Do not run the builder on macOS or a login node; do not modify `refresh_env.sh`, `setup_env_sbatch.sh`, `datasets.json`, `pixi.toml`, or `pixi.lock`.

1. Confirm no active Slurm jobs/writers and validate the host `py-cuda13` environment with the existing direct-R preflight. Install or select the exact build-time helper with `pixi global install -c https://prefix.dev/raphaelribes pixitainer=0.8.3`; record `pixi -V`, `apptainer --version`, lock hash, and source revision.
2. Run `build_ecoda_runtime.sh --layout relocated --output <approved-scratch-image>`. Inspect the dry-run and retained `.def`; reject any definition that reconstructs the environment from only the lockfile, uses Pixitainer seamless `pixi run` at runtime, omits the copied realized environment, or does not retain the explicit Rocky Linux 9 pin.
3. Run a CPU image probe that imports representative Python packages (`numpy`, `scanpy`, `anndata`, `torch`) and representative realized/custom R packages (`arrow`, `DESeq2`, `EPIC`, `GloScope`, `MOFA2`, `scITD`, `HiTME`, `scATOMIC`, `ProjecTILs`, `scGate`). Record `sys.prefix`, `R.home()`, `.libPaths()`, `RETICULATE_PYTHON`, `PATH`, and `LD_LIBRARY_PATH`; require all runtime paths to resolve inside the selected image prefix and no command to resolve the host `.pixi` tree. Check Python shebangs, R compiled/default home, shared-library dependencies, escaping symlinks, R lazy-load databases, and the annotation package stack.
4. Exercise the actual Stage 4 R-to-Python child call with a five-sample Joanito/debug chunk. Require the R process to call the ordinary in-image `PYTHON_BIN` path, produce the expected Feather/checksum contract, and show no nested Apptainer invocation.
5. If the relocated image fails any prefix/import/child-process check, run `build_ecoda_runtime.sh --layout path-preserving --output <approved-scratch-preserved-image>` and repeat steps 3–4. If both layouts fail, stop before production routing and preserve the exact logs for HPC administrators; do not implement a tarball, duplicate BeeGFS directories, or a delay-only mitigation.
6. Run a GPU probe on the actual Stage 5 `shared-gpu`/`nvidia_h200_nvl` allocation. With `ECODA_APPTAINER_NV=1`, require NVIDIA device visibility and CUDA/PyTorch compatibility. With CPU methods, require no `--nv`. A missing driver, unavailable `--nv`, or incompatible CUDA stack is a terminal feasibility failure; never silently run an approved GPU method on CPU.
7. Run the host Arm A and immutable-image Arm B direct package-load stress at at least 16 concurrent Bamboo nodes, recording node, runtime root, package path, image SHA, and exact error for every task. Then run the bounded `_debug` pipelines 3–5 path. Arm B must have no intermittent missing-library/`ENOENT` failures and must preserve artifact schemas, checksums, absolute manifest paths, run-owned logs/status, and retry state. Arm A is a diagnostic comparison and is not a release criterion.

The multi-node A/B and `_debug` evidence run must use the checked-in `durable-hpc-gate-ecoda` profile from `.agents/skills/durable-hpc-gate-ecoda/references/profile.json`, remote workdir `$HOME/ECODA_paper` on `bamboo`, canonical scratch `$HOME/scratch/ECODA_paper`, and serialization group `ecoda-benchmark`. Prepare/reconcile/launch the exact wrapper once, arm exactly one asynchronous unbounded `wait`, then perform one terminal `inspect` with every scheduler array/watchdog ID emitted by the wrapper followed by the Luna Max reviewer approval inspect. Never poll `squeue`/`sacct` from the agent session and never treat the completion event as evidence.

**Verification:** the chosen layout has a valid SIF/manifest pair, the CPU/GPU/import/R-to-Python probes pass, the ≥16-node Arm B stress has zero missing-library failures, and the reviewed durable gate is `COMPLETED`. If the gate fails, no dependent production gate starts.

**Barrier B2:** no production submitter is switched to Apptainer until one layout passes all image, bind, GPU, child-process, concurrency, and artifact checks and receives durable-gate reviewer approval.

### Phase 3 — Shared scheduler/runtime propagation primitives

**Parallel work unit 3A — `SchedulerPropagation` (exclusive ownership of shared preflight/watchdog helpers).**

Targets: `src/utils/bash/h5ad_preflight_submit.sh`, `src/utils/bash/h5ad_preflight_worker.sh`, `src/utils/bash/r_environment_preflight_worker.sh`, `src/5_run_benchmark_methods/benchmark_submit_common.sh`, `src/5_run_benchmark_methods/watchdog_main.sh`.
Targets also include `src/5_run_benchmark_methods/matrix_watchdog.sh`, the canonical Stage 5 host-side watchdog/retry path.

1. Extend `ecoda_submit_h5ad_preflight` with a required final `runtime_export` argument and append it to the existing `--export="ALL,..."` string. The preflight worker must call `ecoda_runtime_reexec_worker stage3|stage5 <source-script>` before reading the manifest so its contract checks use the same image as the gated worker.
2. Extend `benchmark_submit_watchdog` with a required `runtime_export` argument after `WORKER_SCRIPT`; pass it explicitly to the host-side watchdog job. Keep the watchdog shell/Slurm-only; it must not start Python/R or Apptainer merely to wait.
3. In `watchdog_main.sh::watchdog_resubmit`, validate the inherited runtime configuration immediately before `sbatch`, export the same `ECODA_RUNTIME_MODE`, image, manifest, profile, and method-specific `ECODA_APPTAINER_NV`, and preserve all existing worker flags. A missing/mismatched image must write the existing fail-closed status and must not submit a retry array.
4. Preserve existing `worker_retry.sh` transient signatures and OOM-only watchdog classification. Container startup failures terminate the worker process before the transient retry block; scientific missing-library failures after a successful container entry retain the existing residual requeue behavior.

**Verification:** syntax-check the shared helpers; use fake `sbatch`/`apptainer` captures to prove preflight and watchdog retry commands carry identical runtime exports and no nested launch. Do not run the Stage 5 OOM suite until its caller signature updates land.

**Barrier B3:** shared function signatures and runtime-export behavior are stable; stage-specific owners may now update callers without editing these shared files.

### Phase 4 — Integrate pipelines 2–5 through the worker boundary

All four work units below are independent after B3. Each owner has exclusive write access to its listed paths. Watchdogs remain host-side shell/Slurm processes; only scientific workers and package preflight workers cross the image boundary.

**Parallel work unit 4A — `Stage2Runtime` (exclusive Stage 2 ownership).**

Targets: `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`, `src/2_dataset_specific_preprocessing/stage2_watchdog.sh`, `src/2_dataset_specific_preprocessing/1.1_submit_gongsharma.sh`, `1.2_submit_combinedpbmc.sh`, `1.3_submit_joanito.sh`, `1.4_submit_kfoury_lowres_ct.sh`, `1.5_submit_myocardial.sh`, `1.6_submit_bassez.sh`.

1. Source `ecoda_runtime.sh` after `slurm_config.sh` in every six scheduled hook script and call `ecoda_runtime_reexec_worker stage2 <literal-source-hook-path>` before any data load or interpreter call. Add the in-container guard to each hook's existing Slurm spool recovery so `scontrol` is not required inside the exact-copy image.
2. In `1_submit_hpc.sh`, validate the runtime once before the first Stage 2 `sbatch`; append `ecoda_runtime_export_csv stage2 0` to each hook export. Add the same export to `stage2_watchdog.sh` OOM retry submissions. Preserve the GongSharma-cap dependency and immediate submission of unrelated hooks.
3. Do not route `stage2_watchdog.sh` itself, Pipeline 1 staging, or host environment mutation through the image. Their host-side validators and scheduler calls remain unchanged.
4. CombinedPBMC's host-only module load must be guarded when `ECODA_RUNTIME_IN_CONTAINER=1`; no host `module` command may be required inside the SIF.

**Verification:** run `bash tests/test_stage2_submitter.sh`; captured hook arrays must include the runtime export, preserve the CombinedPBMC-only dependency, and retain `FORCE_PREPROCESS`. A fake invalid image must cause zero hook submissions.

**Parallel work unit 4B — `Stage3Runtime` (exclusive Stage 3 ownership).**

Targets: `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, `src/3_scrnaseq_preprocessing/1.1_run_worker.sh`, `src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh`.

1. Source the runtime helper and call `ecoda_runtime_reexec_worker stage3 <literal-source-worker-path>` before manifest parsing and Python execution. Skip the existing `scontrol` spool-recovery query only when `ECODA_RUNTIME_IN_CONTAINER=1`; retain it on the host outer invocation.
2. Validate runtime before the first array submission. Pass `ecoda_runtime_export_csv stage3 0` to the H5AD preflight helper and preprocessing array. Add it to every OOM retry export in `1.2_preprocess_watchdog.sh` so reduced manifests retain the same image/profile.
3. Leave watchdog H5AD validation host-side, but ensure it continues to see the same absolute scratch/log/run paths and fails closed on missing status or invalid artifacts.

**Verification:** run `bash tests/test_preprocessing_stage_submitter.sh` and `bash tests/test_h5ad_preflight.sh`; captured arrays and retries must preserve exact manifest paths, force flags, scheduler dependencies, and runtime exports.

**Parallel work unit 4C — `Stage4Runtime` (exclusive Stage 4 ownership).**

Targets: `src/4_cell_type_annotation/1_submit_onboarding_stage.sh`, `src/4_cell_type_annotation/1.2_prepare_chunks_worker.sh`, `2.1_run_worker.sh`, `3.2_merge_worker.sh`, `1.3_prepare_chunks_watchdog.sh`, `1.2_annotation_watchdog.sh`, `3.3_merge_watchdog.sh`.

1. Source and invoke the runtime helper before manifest parsing/valid-artifact checks in all three workers, using profiles `stage4` and their literal source script paths. Preserve all run-owned absolute paths, Feather/H5AD validators, checksums, and force semantics.
2. Pass `ecoda_runtime_export_csv stage4 0` to preparation, annotation, and merge arrays. Add it to all three OOM retry exports. The annotation worker's R process must inherit in-image `PYTHON_BIN`; `2.1.1_process_chunk.R` remains scientifically unchanged and its existing `system2(Sys.getenv("PYTHON_BIN"), ...)` call must execute directly inside the image.
3. Keep `1.0_stage_reference_maps.sh`, the one-time `2.0_create_scgate_db.R` host preflight, NAS sync, and all three watchdog processes host-side. Bind `HOME_REF_DIR` read-only for annotation workers and keep `SCGATE_DB_PATH` readable through the explicit project/aux bind; never download the cache from workers.
4. Stage 4 must stage/checksum reference maps before the stage4 bind-profile validation on first run, while the validation itself remains before any scheduler submission. Source-artifact metadata comparisons must use the same logical scratch-path namespace across container and host watchdogs.
5. The host Stage 4 preflight validates both the pinned scGate model DB and a local CellOntology dictionary cache. If the DB is valid but the dictionary is missing, `--validate-db-only` permits one host preflight to stage the dictionary; workers expose that validated repository in their R tempdir and never download it.

**Verification:** run `bash tests/test_annotation_stage_submitter.sh` and the annotation worker policy tests. Fake arrays must carry the runtime profile while scGate creation remains a single host-side job; invalid image/bind input must prevent all Stage 4 array submissions.

**Parallel work unit 4D — `Stage5Runtime` (exclusive Stage 5 ownership).**

Targets: `src/5_run_benchmark_methods/1_submit_hpc_array.sh`, `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh`, `run_r_sample_embedding_methods/1.1_run_worker.sh`, `run_transformation_zeroimp_analysis/1.1_run_worker.sh`, `matrix_watchdog.sh`, `matrix_gate.sh`.
   The canonical Stage 5 host-side `matrix_watchdog.sh` receives and re-emits the same runtime export for OOM retries; `matrix_gate.sh` remains host-side and recovers its repository path from SLURM spool execution.

1. Source/invoke the runtime helper before manifest/output handling in all three worker scripts. Use profile `stage5`; pass literal source paths so inner execution does not depend on `/var/spool/slurmd` paths. Preserve direct Rscript expansion and all force/analysis-pass flags.
2. In `method_spec`, initialize `METHOD_RUNTIME_NV=0`; standard `benchmark_analysis` MrVI/scPoli remain in the H200-pinned default class, while `--gpu-policy any`, non-default GPU methods, and every batch-effect pass use the flexible any-GPU partition list with no model constraint. Keep `pilotgm` in the CPU method branch (no GPU flag or `--nv`) because the authoritative containerization contract classifies PILOT-GM-VAE as CPU.
3. Construct `METHOD_TIME_LIMIT` with a 12-hour default-class limit and a three-hour flexible-class limit, both explicit overrides. Propagate the selected limit through worker arrays and OOM retries. Accept an optional `BENCHMARK_GPU_ANY_VRAM_PER_GPU` request for measured large datasets; record input H5AD size/dimensions and peak CUDA allocated/reserved memory in GPU worker logs.
4. In `submit_matrix`, construct the worker export with `ecoda_runtime_export_csv stage5 ${METHOD_RUNTIME_NV}`. Pass that exact export to the worker array, host-side matrix watchdog, and retries. Pass the same runtime export to Stage 5 H5AD and R environment preflight jobs. Unknown methods, invalid policy/time/VRAM configuration, invalid runtime mode, invalid `ECODA_APPTAINER_NV`, missing image, or missing manifest must fail before any `sbatch` call.
5. Isolate pilotpy's hard-coded `Results_PILOT/plots` side effect in a writable temporary working directory and restore the worker cwd before publishing its validated Feather output; no third-party scratch plot directory may target the read-only project bind.
6. Do not edit compatibility delegators under the three `run_*/*/1_submit_hpc_array.sh` paths; they already delegate to the canonical root submitter. Do not containerize `watchdog_main.sh` itself.

**Verification:** run `bash tests/test_benchmark_matrix_submitter.sh`, `bash tests/test_benchmark_dependencies.sh`, and `bash tests/test_benchmark_matrix_watchdog.sh`. Captured calls must show H200 constraint/time only for the default GPU class, any-GPU partition/no model constraint/short time for flexible jobs and batch passes, `--nv` only for GPU methods, no GPU flag for `pilotgm`/CPU methods, exact runtime exports on every worker/preflight/retry, unchanged pseudobulk dependencies, and zero submissions for unknown mode/method/image.
**Barrier B4:** all four stage owners pass their focused submitter/worker contract tests and `bash -n` checks before final cross-stage verification.

### Phase 5 — Preserve retries, artifacts, and host-mode behavior

**Parallel verification work units 5A–5C (read-only test ownership after all Phase 4 edits).**

1. `RetryOwner` runs `bash src/5_run_benchmark_methods/test_oom_retry.sh` and verifies runtime exports survive every OOM retry while preserving the existing `OUT_OF_MEMORY`-only classification, 128G→256G→500G clamp, reduced manifest, status file, scheduler IDs, and non-OOM fail-closed behavior. Add assertions only for runtime/image/GPU propagation; an image validation/startup failure must never become an OOM retry.
2. `WorkerRetryOwner` runs a focused fake-worker scenario from `tests/test_runtime_container.sh` proving that an Apptainer startup/image/bind failure exits terminally without `scontrol requeue`, while a post-entry scientific `No module named`/`ENOENT` error retains the existing `worker_retry.sh` behavior and cap.
3. `RegressionOwner` runs the existing artifact/run ownership checks that cover changed path boundaries: `bash tests/test_ecoda_run_common.sh`, `bash tests/test_atomic_artifact_writers.sh`, `bash tests/test_benchmark_sync.sh`, `bash tests/test_annotation_reuse_and_cutover.sh`, and `bash tests/test_pipeline_watchdog_ids.sh`. Do not relax checksum, run-owned path, absolute manifest, or NAS-sync barriers.

**Verification barrier B5:** run the complete focused union once, with no concurrent test processes touching shared temporary roots:

```bash
bash -n src/slurm_config.sh \
  src/utils/bash/ecoda_runtime.sh \
  src/utils/bash/build_ecoda_runtime.sh \
  src/utils/bash/h5ad_preflight_submit.sh \
  src/utils/bash/h5ad_preflight_worker.sh \
  src/utils/bash/r_environment_preflight_worker.sh \
  src/utils/bash/worker_retry.sh \
  src/utils/bash/ecoda_run_common.sh \
  src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  src/2_dataset_specific_preprocessing/stage2_watchdog.sh \
  src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  src/3_scrnaseq_preprocessing/1.1_run_worker.sh \
  src/3_scrnaseq_preprocessing/1.2_preprocess_watchdog.sh \
  src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  src/4_cell_type_annotation/1.2_prepare_chunks_worker.sh \
  src/4_cell_type_annotation/2.1_run_worker.sh \
  src/4_cell_type_annotation/3.2_merge_worker.sh \
  src/5_run_benchmark_methods/1_submit_hpc_array.sh \
  src/5_run_benchmark_methods/benchmark_submit_common.sh \
  src/5_run_benchmark_methods/watchdog_main.sh
bash tests/test_runtime_container.sh
bash tests/test_r_environment_preflight.sh
bash tests/test_stage2_submitter.sh
bash tests/test_preprocessing_stage_submitter.sh
bash tests/test_h5ad_preflight.sh
bash tests/test_annotation_stage_submitter.sh
bash tests/test_benchmark_matrix_submitter.sh
bash tests/test_benchmark_dependencies.sh
bash tests/test_benchmark_matrix_watchdog.sh
bash src/5_run_benchmark_methods/test_oom_retry.sh
bash tests/test_ecoda_run_common.sh
bash tests/test_atomic_artifact_writers.sh
bash tests/test_benchmark_sync.sh
bash tests/test_annotation_reuse_and_cutover.sh
bash tests/test_pipeline_watchdog_ids.sh
```

Expected result: all scripts print their existing `OK` markers; host mode remains the default; no source invokes `pixi run` for workers; every retry/gate retains runtime metadata; no checksum, ownership, dependency, or scientific-contract test is weakened.
The verification set also covers Stage 4 multi-view selection (`--datasets/--views`), container-only CombinedPBMC module suppression, and logical scratch-symlink source identity.

### Phase 6 — Bounded `_debug` rollout and production release gate

**Work unit 6A — `DebugRuntimeGate`.**

1. Set `ECODA_RUNTIME_MODE=apptainer`, `ECODA_RUNTIME_IMAGE`, and `ECODA_RUNTIME_MANIFEST` to the reviewed layout; leave host mode as the explicit rollback/default and do not mutate the host environment during worker execution.
2. Run the actual `_debug` Joanito five-sample path through Stage 3, Stage 4, and Stage 5 using the existing submitters with `--gpu-policy any` for the debug GPU methods so it can backfill on any compatible GPU; reserve the H200 requirement for standard default benchmark runs. Exercise at least one CPU R worker, one CPU Python worker, the annotation R-to-Python call, and separate MrVI/scPoli GPU smoke paths. Record input H5AD size/dimensions and peak CUDA memory, plus image SHA/runtime identity, explicit bind destinations, valid artifacts/checksums, complete run-owned manifests, and successful watchdog/retry status.
3. Do not launch full cohorts for this phase. If any image, bind, GPU, child-process, concurrency, or artifact check fails, restore host mode for diagnosis and stop production release; do not use a startup delay or ordinary retry as proof of correctness.

Run this gate with the checked-in `durable-hpc-gate-ecoda` profile, one launch and one unbounded wait, then one terminal inspect plus reviewer approval. A reviewed `COMPLETED` gate is required for B6.

**Work unit 6B — `ProductionRuntimeGate` (only after B6).**

1. Prepare a new durable gate for the exact approved production selection. The wrapper must source `src/slurm_config.sh`, set the reviewed `ECODA_RUNTIME_MODE=apptainer` image/manifest, and invoke the existing pipeline 2–5 submitters without changing datasets or method selections. Use the profile's remote clone `$HOME/ECODA_paper`, scratch `$HOME/scratch/ECODA_paper`, serialization group, and immutable command contract.
2. Preserve the handoff's prelaunch checks: Bamboo `HEAD == origin/master`, no active jobs, no environment lock/writer, valid selected manifest, reviewed immutable runtime A/B evidence, and compute-node R preflight first for Stage 5. Use the exact 30-row forced selection only when the approved artifact-repair gate is the intended production action; otherwise use the explicitly selected pipeline scope.
3. After the single durable wait, inspect every scheduler ID emitted by Stage 2–5 arrays, preflights, watchdogs, retries, and aggregate gates in one terminal audit. Require every worker artifact/schema/checksum, run-owned status, source identity, and NAS synchronization contract to pass. Reviewer approval is required before marking the runtime rollout complete.

**End-to-end verification:** the reviewed `_debug` gate proves the changed runtime surface and the reviewed production gate proves the selected pipelines complete with the same scientific/artifact contracts and no intermittent missing-library failures. A missing status, checksum, scheduler row, image manifest, or reviewer approval is a failure, not an optional warning.

## Critical files & anchors

- `src/slurm_config.sh` — host/container interpreter branches and runtime environment variables.
- `src/utils/bash/ecoda_runtime.sh` — image validation, explicit bind profiles, one worker boundary, and export contract.
- `src/utils/bash/build_ecoda_runtime.sh` — pinned Pixitainer/Rocky build, realized-environment copy, provenance manifest, and atomic image publication.
- `src/5_run_benchmark_methods/1_submit_hpc_array.sh` — method/resource dispatch, `mrvi`/`scpoli` GPU passthrough, `pilotgm` CPU correction, and Stage 5 runtime exports.
- `src/4_cell_type_annotation/2.0_create_scgate_db.R` and `2.1.1_process_chunk.R` — one-time host cache creation plus offline worker exposure of the HiTME CellOntology dictionary.
- `src/4_cell_type_annotation/2.1_run_worker.sh` — R worker plus direct R-to-Python child boundary; its scientific R source remains unchanged.

## Assumptions & contingencies

- Bamboo provides user-level Apptainer execution and a compute allocation suitable for SIF construction. If Apptainer or the pinned Pixitainer 0.8.3 helper is unavailable, fail closed at Phase 2 and preserve evidence for the cluster administrator; do not introduce a second packaging system in this change.
- `rockylinux:9` and Pixitainer `0.8.3` are fixed build inputs. The exact host Pixi and Apptainer versions are discovered on the build allocation, passed/recorded explicitly, and validated from the retained definition/manifest.
- The current documentation is authoritative for `--nv`: only `mrvi` and `scpoli` receive passthrough. The current code's `pilotgm` GPU grouping is corrected to CPU rather than silently giving it a CPU run on a GPU allocation or silently exposing a GPU to every Python method.
- The path-preserving fallback is intentionally user/path-specific and is used only when relocated-prefix probes fail. If it also fails, production routing is blocked; no unvalidated relocation, tar extraction, broad `/home`/`/srv` bind, or duplicate BeeGFS environment is accepted.
- Host-side `refresh_env.sh`, `setup_env_sbatch.sh`, reference-map staging, scGate cache creation, watchdog accounting, and NAS sync remain outside the read-only SIF. Their existing direct host runtime is intentional and does not weaken the worker isolation contract.
- Runtime construction observed a non-relocatable R prefix under `/opt/ecoda/py-cuda13`; the approved image uses the exact project-root path-preserving layout after adding `which`/`jq` to the Rocky base. Relocated images remain diagnostic evidence and are not routed.
- Existing unrelated working-tree changes are user-owned and must remain unstaged and unmodified throughout implementation and verification.

## Execution observations

- Host-side submitter validation was explicitly tested with `ECODA_RUNTIME_MODE=apptainer` and `ECODA_RUNTIME_IN_CONTAINER=0`; it retained the host `.pixi` interpreter paths.
- The path-preserving SIF passed the compute-node CPU Python/R/reticulate child probe and the H200 CUDA probe; the relocated SIF failed the R prefix probe and was not approved.
- The first 16-node image stress used the wrong default image path and failed closed; a corrected image-only 16-node durable gate passed with `IMAGE_FAILURES=0` and reviewer approval.
- The first forced scGate-cache debug gate downloaded and validated the model DB/dictionary on its host preflight, then failed closed at the immediate post-build cache check; the cache remained valid and a cache-reuse gate confirmed all Stage 3/4 arrays/watchdogs and annotation merge completed.
- The cache-reuse gate then failed closed before Stage 5 because NAS H5AD sidecars retained the scratch `PATH=`. `ecoda_compare_checksum_remote` now atomically canonicalizes and revalidates the remote sidecar path after digest/size equality.
- The sidecar-fixed debug gate reached all Stage 5 arrays/watchdogs, then failed closed at the aggregate gate because pilotpy attempted to create `Results_PILOT/plots` under the read-only project bind; `run_pilot` now isolates that third-party relative write in a writable temporary directory.
- The PILOT-isolated debug gate completed every worker/watchdog array, including MrVI/scPoli CUDA runs, but its aggregate gate failed closed on a SLURM spool-path source lookup; `matrix_gate.sh` now has the same repository recovery contract as the other host-side gates.
- The aggregate-gate-fixed debug run completed every Stage 3–5 array/watchdog and aggregate gate, then failed closed during final RDS validation because GloScope emitted alphabetically ordered sample IDs; `collapse_sample_metadata` now restores first-appearance H5AD order before result bundles are built.
- The H200-constrained sample-order debug gate was canceled after its scheduler estimate reached three days; the revised gate uses `--gpu-policy any` and the three-hour flexible worker limit so debug validation does not consume the H200 queue.
- The telemetry any-GPU gate `ecoda_runtime_debug_telemetry_20260830T125637Z` completed with one durable wait and a passing terminal accounting/artifact audit. MrVI and scPoli ran on an NVIDIA GeForce RTX 3090 (23.80 GiB); the input was 259,543,148 bytes, 6,000 cells, and 17,903 genes, with observed peak reserved CUDA memory 0.07 GiB (MrVI) and 0.06 GiB (scPoli). This is debug-scale evidence only; large batch-effect datasets require their own measured VRAM/time class.
- Gate-specific checks passed for both `_debug` benchmark and uncorrected batch-effect H5AD contracts (including annotation sidecars and 6,000-row schemas); the fresh Luna Max reviewer inspect approved the telemetry gate and set `release_eligible=true`.
- Metadata-only HDF5 sizing of all existing batch-effect outputs found 103,642–2,095,732 cells, 2,000-HVG dense fp32 proxies of 0.77–15.61 GiB, and full dense fp32 bounds of 10.15–180.93 GiB. These are planning bounds, not measured model VRAM; each production tier must use the new worker telemetry before selecting a VRAM/time class.
- Post-release advisory correction: `slurm_config.sh` now preserves inherited absolute `HPC_SCRATCH_DIR`, `HOME_REF_DIR`, and runtime image values inside `--containall --no-home` workers instead of deriving path contracts from container `HOME`; the runtime regression covers the nested source.
- Post-release advisory correction: the Python worker now passes `--device cuda` whenever `ECODA_APPTAINER_NV=1`, and the Python entry point rejects `auto`, `cpu`, or unavailable CUDA for MrVI/scPoli before reading a dataset. Focused dispatch and CUDA-guard regressions pass; a fresh durable debug gate is required before release approval.
- The first post-release CUDA gate captured the pre-sync Bamboo HEAD and was intentionally not released: after preserving the remote WIP in `stash@{0}` and fast-forwarding Bamboo to c293a89, its single terminal inspect recorded the repository-head fingerprint mismatch and `state=FAILED`.
- The path-preserving runtime rebuild gate `ecoda_runtime_build_retry_20260830T142233Z` produced image SHA-256 `d0949ae095c6eef4343e8aeb61409551bf1d36abaaa50a400dc80f04d48f1f69` with `GIT_REVISION=c293a89`; terminal audit, direct image-contract validation, and Luna Max review all passed.
- The source-consistent CUDA gate `ecoda_runtime_debug_cuda_clean_20260830T143643Z` completed all Stage 3–5 workers/watchdogs/aggregate with one passing terminal audit and Luna Max approval. MrVI reported CUDA true on an RTX 3090 (23.80 GiB), scPoli reported CUDA true on an NVIDIA RTX PRO 6000 Blackwell Server Edition (95.01 GiB), and both emitted GPU memory telemetry; release eligibility is true.
