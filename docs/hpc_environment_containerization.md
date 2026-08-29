# HPC Runtime Environment Containerization for ECODA

> Working technical findings and discussion record. This is not the final implementation plan or final summary. Implementation, workflow selection, and production validation remain separate later steps.

## 1. Scope and intended outcome

The intended outcome is to containerize the runtime environment used to spawn workers in **pipelines 2–5** while preserving the existing pipeline scripts and scientific behavior as far as possible.

The following are explicitly outside this runtime-containerization scope:

- Pipeline 1 NAS staging.
- Local macOS analysis workflows.
- Host-side Pixi environment installation or refresh.
- `src/utils/bash/refresh_env.sh`.
- `src/utils/bash/setup_env_sbatch.sh`.
- Creation of a second parallel set of pipeline scripts for containerized execution.

The same worker and scientific scripts should remain usable in both host and container execution modes. Differences should be confined to the runtime-launch boundary, configuration, and small shared helpers.

## 2. Problem motivating containerization

When many SLURM workers start concurrently on Bamboo, they all traverse the shared Pixi environment under the BeeGFS-backed project directory. Python imports and R package loading perform many filesystem metadata operations over thousands of files. The reported failure is intermittent missing-library behavior, including errors equivalent to `No such file or directory`, missing R package metadata, and Python module lookup failures.

The likely mechanism is a BeeGFS metadata and small-file I/O storm. The precise client-cache behavior remains a hypothesis, but the mitigation does not depend on proving the exact BeeGFS failure mechanism: reducing each worker's independent traversal of the shared environment is beneficial regardless of whether the immediate failure is a metadata timeout, stale negative lookup, or another transient filesystem error.

A SIF image stores the runtime as a compressed, read-only SquashFS filesystem inside one image file. BeeGFS then handles one image path instead of tens of thousands of independent environment paths. The image can still generate ordinary data reads from BeeGFS, especially on a cold node; containerization reduces metadata fan-out but does not eliminate all network I/O.

## 3. Current ECODA environment model

`src/slurm_config.sh` currently exports direct binaries from the HPC `py-cuda13` environment:

```text
PYTHON_BIN=${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/python
PIXI_RSCRIPT=${PROJECT_ROOT}/.pixi/envs/py-cuda13/bin/Rscript --vanilla
```

This direct R invocation is intentional. Worker processes must not activate Pixi on every task because the R activation hook can execute `R CMD javareconf` and write shared R configuration files.

Most pipeline workers and validators already call these variables rather than hard-coding `python` or `Rscript`. This provides a useful central integration point, but it does not mean that a global replacement is automatically safe.

The host environment setup scripts also use `PIXI_RSCRIPT` to run environment mutation and validation. In particular, `refresh_env.sh` runs the guarded R package installation and smoke checks through that variable. Therefore, the host mutation path and the container runtime path must remain distinguishable.

The realized HPC environment contains more than the packages represented directly in `pixi.lock`. `src/utils/setup_r_packages.R` installs additional SHA-pinned GitHub, Bioconductor, and CRAN packages, including the HiTME, STACAS, ProjecTILs, scATOMIC, EPIC, scECODA, MOFA2, GloScope, and scITD-related layers. Any packaging method that only reconstructs the lockfile environment must be checked for these additions.

Relevant repository files:

- [`src/slurm_config.sh`](../src/slurm_config.sh)
- [`src/utils/bash/refresh_env.sh`](../src/utils/bash/refresh_env.sh)
- [`src/utils/bash/setup_env_sbatch.sh`](../src/utils/bash/setup_env_sbatch.sh)
- [`src/utils/setup_r_packages.R`](../src/utils/setup_r_packages.R)

## 4. Candidate packaging approaches

### 4.1 Apptainer SIF from the realized environment

This is the primary candidate. Build one immutable SIF on Bamboo from the already-validated HPC `py-cuda13` environment. The image should preserve the installed custom R packages rather than re-solving a new environment.

Advantages:

- One runtime image file on BeeGFS.
- No per-worker extraction of tens of thousands of files.
- Captures the actual installed R and Python package state.
- Direct execution of internal Python and R binaries avoids Pixi activation at worker startup.
- Image hash and labels can provide a stable runtime identity.

Risks:

- Copying the environment to a new path can expose absolute-prefix references in scripts, R configuration, or compiled files.
- The SIF still requires reads from the image file on cold nodes.
- GPU execution depends on host NVIDIA drivers and Apptainer's `--nv` integration.
- The image is read-only; project outputs, logs, caches, and temporary files must be bound to writable host paths.
- Updating the realized host environment requires a deliberate image rebuild.

### 4.2 Pixi Pack or a generic archive

`pixi-pack` can package a lockfile-defined Pixi environment into an archive and `pixi-unpack` can restore it. This is useful for portable environment distribution, but it is not automatically equivalent to archiving ECODA's realized environment because custom R packages installed after Pixi are not necessarily represented in the lockfile package set.

A raw tar archive of the realized environment would preserve those files, but a worker would need to extract it to node-local storage. That adds extraction time, local disk usage, cleanup, and prefix-relocation concerns. It remains a reasonable fallback if Apptainer cannot be used.

### 4.3 Multiple environment directories on BeeGFS

This is not recommended. Replicating the environment across several directories does not remove the shared-filesystem metadata problem and increases storage and maintenance burden.

### 4.4 `cotainr` or a fresh Conda/YAML rebuild

The UniGe documentation recommends `cotainr` for creating a container from a Conda YAML file. That is appropriate for a clean image build, but a fresh YAML-based rebuild risks omitting ECODA's post-Pixi R packages or changing the realized dependency state. It should not be the first method for preserving the current validated environment.

## 5. Pixitainer assessment

[`pixitainer`](https://github.com/RaphaelRibes/pixitainer) is useful as a convenience tool, but its default behavior is not sufficient for ECODA.

The current Pixitainer SIF implementation:

1. Copies `pixi.toml` and `pixi.lock` into the image.
2. Installs Pixi inside the image.
3. Runs a fresh `pixi install ... --frozen` during the image build.
4. Uses a seamless runscript equivalent to:

   ```bash
   pixi run --locked --as-is -m /opt/conf/pixi.toml "$@"
   ```

That default path does not copy the existing `.pixi/envs/py-cuda13` directory and uses the Pixi runtime wrapper that ECODA deliberately avoids.

Pixitainer also supports the relevant escape hatches:

- `--manual` creates a direct shell-style entrypoint instead of the seamless `pixi run` entrypoint.
- `--no-install` skips installation of the project environment inside the image.
- `--add-file SRC:DEST` can copy an existing environment directory into the image.
- `--keep-def` retains the generated Apptainer definition file.
- `--dry-run` allows inspection of the generated definition before building.

Therefore, a candidate exact-copy build could conceptually use:

```bash
pixi containerize \
  --manual \
  --no-install \
  --env py-cuda13 \
  --base-image=rockylinux:9 \
  --add-file "$PROJECT_ROOT/.pixi/envs/py-cuda13:/opt/ecoda/py-cuda13" \
  --keep-def \
  --output "$PROJECT_ROOT/.containers/ecoda-py-cuda13.sif"
```

This is a feasibility-test command shape, not an approved production command. The generated definition must be inspected, and the copied environment must be tested for external symlinks and absolute-prefix behavior. The Pixitainer version should be pinned if it becomes part of the build process.

If this exact-copy mode works, Pixitainer can avoid maintaining much of the low-level Apptainer definition boilerplate. It should remain a build-time helper; pipeline workers should invoke direct internal Python/R binaries, not Pixitainer or Pixi. If the generated definition cannot preserve the environment correctly, a small project-specific `.def` file is safer and more transparent.

## 6. Image-build design

The image should be built from the validated Linux HPC environment on a Bamboo compute allocation, not from the macOS environment and not as a heavy build on the login node.

The build needs sufficient temporary space because Apptainer constructs an uncompressed temporary filesystem before producing the compressed SIF. `APPTAINER_TMPDIR` and, if needed, `APPTAINER_CACHEDIR` should point to suitable node-local or approved temporary storage.

The base image should be selected explicitly and tested against Bamboo's Rocky Linux host ABI. An arbitrary default base image should not be accepted for production. The image should carry provenance such as:

- Pixi lockfile hash.
- Git revision of the source tree used for the build.
- Installed environment identity.
- Pixitainer and Apptainer versions, if applicable.
- A checksum of the final SIF.

The SIF should live outside Git-tracked source files, for example in a dedicated ignored container directory or an approved HPC image directory.

### Prefix and relocation requirement

Copying the environment from the project path to `/opt/ecoda/py-cuda13` may require relocation handling. `PATH`, `R_HOME`, `RETICULATE_PYTHON`, and `LD_LIBRARY_PATH` are necessary runtime settings but are not proof that every package is relocatable.

The feasibility phase must check:

- Python console-script shebangs.
- R's compiled/default home path.
- Reticulate's selected Python.
- Shared-library dependencies.
- Symlinks escaping the copied environment.
- R package lazy-load databases and compiled extensions.

If relocation fails, alternatives are:

1. Use a path-preserving environment inside the SIF and mount project source at a different internal path such as `/workspace`, avoiding a host project bind over the embedded environment.
2. Use a tested Conda-compatible relocation mechanism.
3. Build a custom definition that recreates the environment and installs the pinned post-Pixi R packages during `%post`.

The first option is less portable between users but may be the most faithful to the current realized environment.

## 7. Runtime integration for pipelines 2–5

The goal is one interoperable runtime contract, not duplicate containerized and non-containerized script trees.

A useful conceptual mode is:

```text
ECODA_RUNTIME_MODE=host
ECODA_RUNTIME_MODE=apptainer
```

In host mode, the existing direct Pixi environment paths remain active. In Apptainer mode, the same worker scripts see direct paths to the binaries inside the image.

The scientific Python and R scripts should not need separate copies. The runtime layer should provide the correct executable paths and filesystem mounts for the selected mode.

### Preferred boundary: execute scheduled workers inside the image

The cleanest R/Python behavior is to cross the container boundary once per scheduled worker. The existing worker shell script runs inside the image, sources a container-aware runtime configuration, and invokes the internal Python/R binaries directly.

Then `2.1.1_process_chunk.R` remains unchanged:

```r
system2(Sys.getenv("PYTHON_BIN"), ...)
```

In host mode, `PYTHON_BIN` is the host Python executable. In container mode, `PYTHON_BIN` is the internal image Python executable. Both are ordinary executable paths, so R's `system2()` behaves identically.

This avoids launching a second nested Apptainer process when R calls Python.

### Alternative boundary: executable command wrappers

A host-side wrapper can preserve most existing call sites because they already use `PYTHON_BIN` and `PIXI_RSCRIPT`. The wrapper must be a real executable file, not a multiword string such as `apptainer exec image.sif python`.

If R itself is launched through a wrapper, the containerized R process must see the internal Python path. Otherwise, its `system2(PYTHON_BIN, ...)` call may invoke a host wrapper from inside the container and attempt an unsupported nested Apptainer launch.

A wrapper implementation therefore needs either:

- a container-detection branch that directly executes the internal Python when already inside Apptainer; or
- explicit in-container environment overrides for `PYTHON_BIN`, `PIXI_RSCRIPT`, `R_HOME`, and `RETICULATE_PYTHON`.

The scheduled-worker boundary is preferable because it makes this behavior natural and avoids nested-container edge cases. Either design can preserve one shared set of scientific scripts.

## 8. GPU mode selection

GPU selection is technically straightforward and does not require a second pipeline implementation.

Stage 5 already dispatches methods by method name and selects different worker/resource classes. The shared runtime launcher can consume an explicit per-job flag, for example:

```text
ECODA_APPTAINER_NV=1   # MrVI and scPoli
ECODA_APPTAINER_NV=0   # CPU methods
```

The Stage 5 submitter or worker-dispatch branch would set the flag only for `mrvi` and `scpoli`. The launcher then adds `--nv` only when the flag is enabled. CPU methods such as PILOT, QOT, PILOT-GM-VAE, and the R benchmark methods use the same launcher without GPU passthrough.

This is a small conditional runtime change, not a second set of worker scripts. The mapping should fail closed for unknown methods rather than silently assuming GPU or CPU behavior.

The GPU feasibility test must use the repository's actual Bamboo GPU partition, constraints, and resource configuration. A generic `public-gpu` example is not sufficient. Apptainer's `--nv` exposes host NVIDIA devices and driver libraries; it does not replace the host driver. CUDA/PyTorch compatibility must be tested on the actual H200 nodes.

## 9. Safe bind-mount design

Binding the entire `/home` and `/srv` trees is not the preferred production solution. It can expose unrelated data, hide image paths, fail when a source path is unavailable, and obscure which filesystem a worker is actually using.

The runtime should define a small bind profile with explicit source and destination paths:

| Resource | Access | Requirement |
| --- | --- | --- |
| Project source/configuration | Read-only where possible | Must contain scripts, `datasets.json`, and required auxiliary files |
| ECODA scratch data | Read/write | Must preserve the logical paths used by manifests and validators, or apply a complete path translation |
| NAS/reference maps | Read-only | Bind only when a worker needs them |
| Logs and run metadata | Read/write | Bind the exact output directories that workers update |
| Node-local temporary directory | Read/write | Bind the actual `$TMPDIR`/`/scratch` path used by the cluster |
| Image file | Read-only | Opened by Apptainer; not a data bind inside the image |

The host `$HOME/scratch` path is a symlink to a BeeGFS scratch target on UniGe clusters. The launcher should resolve the real source path before constructing the bind specification. The destination path must be chosen deliberately because run-owned manifests contain absolute paths.

The launcher should also avoid accidental host-runtime contamination:

- Set `PATH` to the internal environment first.
- Set `R_HOME`, `RETICULATE_PYTHON`, and `LD_LIBRARY_PATH` deliberately.
- Prevent host `PYTHONPATH` or user-site packages from overriding the image.
- Use `--cleanenv` only with a curated set of required SLURM, CUDA, threading, path, and user variables; otherwise test the inherited environment carefully.
- Keep image paths outside the host project bind if the image contains an environment at a path that the project bind would mask.

These bind rules solve the broad-bind problem without requiring separate scientific pipeline implementations.

## 10. Host-side environment mutation

The phrase “host-side environment mutation” refers only to the existing controlled process that creates or updates the shared Pixi/R environment:

- `refresh_env.sh`.
- `setup_env_sbatch.sh`.
- Their existing Pixi installation and guarded R-package setup.

Those scripts can remain as they are conceptually. They are not worker-runtime scripts and should not be executed inside the read-only SIF.

The intended lifecycle is:

1. Validate or refresh the host HPC Pixi environment outside Apptainer.
2. Build or rebuild the SIF manually from that validated environment.
3. Run pipelines 2–5 through the SIF.
4. If the host environment changes, deliberately rebuild and revalidate the SIF.

No automatic container rebuild is required for routine worker submission.

## 11. Feasibility and validation sequence

Before changing production callers:

1. Generate Pixitainer's definition with `--dry-run` and inspect it.
2. Verify that `--manual --no-install --add-file` copies the realized environment rather than reconstructing only the lockfile environment.
3. Build one candidate SIF on a compute allocation with sufficient temporary space.
4. Run a CPU smoke test covering Python imports, R imports, custom R packages, `R.home()`, Python `sys.prefix`, and reticulate's selected Python.
5. Verify that no runtime command resolves to the host `.pixi` directory.
6. Run a GPU smoke test with the actual Stage 5 GPU allocation and `--nv`.
7. Exercise the R-to-Python child call used by annotation workers.
8. Run a small concurrent `_debug` worker test for pipelines 3–5, followed by the standard `_debug` path through pipelines 3–5.
9. Verify artifact schemas, checksums, run-owned manifests, logs, temporary files, and retry behavior.
10. Only after these checks pass, route production pipeline 2–5 workers through the image.

The verification should distinguish:

- Runtime image correctness.
- Filesystem bind correctness.
- GPU passthrough correctness.
- R-to-Python process interoperability.
- Concurrent startup behavior.
- Scientific/artifact contract preservation.

## 12. Current findings and decisions

- Apptainer/SIF is a strong primary mitigation for the BeeGFS small-file startup problem.
- Multiple environment directories on BeeGFS are not a good solution.
- The default Pixitainer seamless mode is unsuitable because it uses `pixi run` at runtime.
- Pixitainer remains a viable build-time helper when used with manual mode, no environment installation, and explicit copying of the validated realized environment.
- The generated Apptainer definition should be retained for auditability.
- The exact current realized `py-cuda13` environment should be preserved where possible; rebuilding only from `pixi.lock` is insufficient without accounting for post-Pixi R packages.
- Pipelines 2–5 should share one interoperable runtime design rather than having separate containerized and non-containerized script trees.
- `refresh_env.sh` and `setup_env_sbatch.sh` remain outside the container and retain their host-side purpose.
- GPU mode can be selected by a method-based conditional in the existing Stage 5 dispatch path; only MrVI and scPoli should receive `--nv`.
- Broad `/home` and `/srv` binds should be replaced by explicit, path-aware bind profiles.
- R-to-Python interoperability is possible without duplicate scientific scripts when `PYTHON_BIN` is always an ordinary executable path in the current execution context. A worker-level container boundary is the cleanest way to ensure this.
- No production caller or pipeline script has been changed as part of this findings document.

## 13. Discussion log

### Entry 1 — Containerization findings and requested constraints

The requested documentation should explain why HPC environment containerization is being considered, how Apptainer/SIF and Pixitainer relate to the problem, what each approach does, and how a future implementation would be used. It should preserve traceability while implementation choices remain open.

The target is specifically pipelines 2–5, where many workers need the shared Python/R environment. Other repository components can remain non-containerized. The host Pixi environment is considered final or nearly final; any future host environment update can be followed by a deliberate SIF rebuild rather than an automatic runtime rebuild.

The desired architecture is one interoperable script set. Separate duplicate pipelines for containerized and non-containerized operation are not acceptable. Small shared launchers, runtime configuration changes, and method-dependent flags are acceptable if they keep the scientific workers common.

The feasibility tests should precede production caller changes. They must cover the actual realized environment, CPU and GPU runtime behavior, explicit binds, R-to-Python calls, concurrent `_debug` workers, and artifact contracts.

Future discussion entries should be appended below this section rather than rewriting earlier findings. The final summary and implementation plan are intentionally deferred until the workflow has been selected and validated.

### Entry 2 — Runtime integration clarifications

The containerized execution requirement applies to pipelines 2–5 because those stages spawn the workers that currently traverse the shared Pixi environment. Pipeline 1, local analysis, and the host environment-maintenance scripts do not need to be containerized.

The host Pixi environment is expected to remain the source of truth for image construction. `refresh_env.sh` and `setup_env_sbatch.sh` may continue using the existing host-side `PIXI_RSCRIPT`; they must not be routed into a read-only runtime image. If the environment changes later, the SIF can be rebuilt and revalidated manually.

GPU selection is expected to be a small method-dependent conditional in the existing Stage 5 dispatch path, not a second collection of pipeline scripts. MrVI and scPoli should receive Apptainer `--nv`; CPU methods should not. The launcher should reject unknown method/runtime combinations rather than silently selecting a mode.

Broad `/home` and `/srv` binds should be replaced with explicit path-aware binds. The project source, ECODA scratch, NAS/reference paths, logs, and node-local temporary storage have different access requirements. The `$HOME/scratch` symlink must be resolved deliberately, and bind destinations must remain compatible with absolute paths in run-owned manifests.

R-to-Python interoperability should remain identical in host and container modes by ensuring that `PYTHON_BIN` is always a real executable path in the current process namespace. A worker-level container boundary is preferred because R can then call the internal Python directly through `system2()` without nested Apptainer execution. If wrappers are used instead, they need an in-container bypass or explicit environment overrides.

The preferred Pixitainer role is a pinned build-time convenience helper using manual mode, no fresh environment installation, and an explicit `--add-file` copy of the validated realized `py-cuda13` environment. The generated definition remains part of the audit trail. Pixitainer itself should not be invoked by production workers.

## References

### Repository references

- [`src/slurm_config.sh`](../src/slurm_config.sh)
- [`src/utils/bash/refresh_env.sh`](../src/utils/bash/refresh_env.sh)
- [`src/utils/bash/setup_env_sbatch.sh`](../src/utils/bash/setup_env_sbatch.sh)
- [`src/utils/setup_r_packages.R`](../src/utils/setup_r_packages.R)
- [`src/4_cell_type_annotation/2.1.1_process_chunk.R`](../src/4_cell_type_annotation/2.1.1_process_chunk.R)
- [`docs/hpc_docs/applications_and_libraries.md`](hpc_docs/applications_and_libraries.md)
- [`docs/hpc_docs/storage_on_hpc.md`](hpc_docs/storage_on_hpc.md)

### External references

- [UniGe HPC applications and libraries](https://doc.eresearch.unige.ch/hpc/applications_and_libraries)
- [Apptainer definition files](https://apptainer.org/docs/user/latest/definition_files.html)
- [Apptainer build environment and temporary storage](https://apptainer.org/docs/user/latest/build_env.html)
- [Apptainer bind paths and mounts](https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html)
- [Apptainer GPU support](https://apptainer.org/docs/user/latest/gpu.html)
- [Pixitainer README](https://github.com/RaphaelRibes/pixitainer/blob/main/README.md)
- [Pixitainer SIF generation](https://raw.githubusercontent.com/RaphaelRibes/pixitainer/main/lib/sif.sh)
- [Pixitainer common installation logic](https://raw.githubusercontent.com/RaphaelRibes/pixitainer/main/lib/common.sh)
- [Pixi Pack documentation](https://github.com/prefix-dev/pixi/blob/main/docs/deployment/pixi_pack.md)
