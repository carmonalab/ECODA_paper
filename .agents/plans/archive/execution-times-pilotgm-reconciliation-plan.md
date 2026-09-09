# Execution-time reconciliation with PILOT-GM-VAE recovery

Plan slug: `execution-times-pilotgm-reconciliation`
Canonical plan artifact: `local://execution-times-pilotgm-reconciliation-plan.md`
Requested repository copy path: `.agents/plans/execution-times-pilotgm-reconciliation-plan.md`

## Context

The primary deliverable is to add the missing recorded PILOT-GM-VAE runtimes to `data/benchmark/embeddings/execution_times.feather` without rerunning benchmark methods, submitting HPC work, committing, pushing, pulling, or touching the authoritative runtime while the Bassez rolling pipeline can still write it. The secondary deliverable is to determine whether the conflicting Bassez rows came from the active Bassez rollout and replace them only when run-owned provenance proves that they are newer; otherwise leave every Bassez value unchanged.

The current local merged Feather has 1,376 rows, 12 dataset keys including `_debug`, and 118 method labels. It has exact four-column schema, unique `(dataset, method)` keys, finite runtimes, and current MD5 `4df267b7492cce9937db26a717d3aa24`. Retained per-task logs contain 1,478 rows and 1,388 unique keys. They provide 27 task-only PILOT-GM-VAE rows, while the merged file has 130 value conflicts overall: 118 for Bassez, 7 for GongSharma, and 5 for `_debug`.

The referenced Bassez rollout plan is a plausible source of updated Bassez times: its B5-BASSEZ command includes `pilotgm` and `trans,zeroimp` (`.agents/plans/20260830170000-bassez-rolling-pipeline-rollout-plan.md:88-98`), its prior B5 benchmark-rest checkpoint records a PILOT-GM-VAE watchdog timeout (`:449-475`), and its current-wave section treats B5-BASSEZ as a live promotion gate (`:630-645`). Local retained task logs alone do not establish that their Bassez values are newer than the merged rows.

## Approach

### 1. Freeze the active pipeline and establish an immutable local baseline

1. Do not read from or write to the active Bamboo scratch/NAS runtime paths while the Bassez rollout gate can still write or synchronize benchmark artifacts. Do not invoke `sbatch`, any submitter, any durable gate, SSH synchronization, `git pull`, `git push`, or `git commit`.
2. Work only from the local snapshot after the active gate reaches its terminal/reviewed state. Hash and preserve the current `data/benchmark/embeddings/execution_times.feather`, its `.md5` sidecar, and the `data/benchmark/checksums.md5` entry before constructing a candidate.
3. Preserve the canonical raw identifiers used by the Feather file: `Gongsharma_cmv_young_males` and method labels such as `PILOT-GM-VAE_hvg2000_highres`. Do not use notebook display normalization (`GongSharma`, stripped `_highres`) when writing the Feather.
4. Build an explicit reconciliation manifest outside the authoritative file. Each selected row records `dataset`, `method`, `time_secs`, `mem_GB`, source Feather path, source modification time, source run/gate identifier if available, and a status of `CURRENT_RUN_VERIFIED`, `RETAINED_TASK_LOG`, or `UNVERIFIED_RETAINED_TASK_LOG`.

### 2. Add the missing PILOT-GM-VAE rows first

Use the current reviewed run-owned log if one exists for a key after the active gate; otherwise use the retained task log named below. Preserve the source row's `mem_GB` and exact `time_secs`. Add the following 27 task-only rows; do not add a fabricated Gongsharma hvg3000 row because no such retained output/log exists.

| dataset | method | time_secs | retained source |
|---|---|---:|---|
| Adams | PILOT-GM-VAE_hvg1000_highres | 2658.040784 | `execution_times_pilotgm_Adams.feather` |
| Adams | PILOT-GM-VAE_hvg2000_highres | 2715.320060 | `execution_times_pilotgm_Adams.feather` |
| Adams | PILOT-GM-VAE_hvg2000_lowres | 1944.358890 | `execution_times_pilotgm_Adams.feather` |
| Adams | PILOT-GM-VAE_hvg3000_highres | 2705.844977 | `execution_times_pilotgm_Adams.feather` |
| Gongsharma_cmv_young_males | PILOT-GM-VAE_hvg1000_highres | 14725.618864 | `execution_times_pilotgm_Gongsharma_cmv_young_males.feather` |
| Gongsharma_cmv_young_males | PILOT-GM-VAE_hvg2000_highres | 13990.631948 | `execution_times_pilotgm_Gongsharma_cmv_young_males.feather` |
| Gongsharma_cmv_young_males | PILOT-GM-VAE_hvg2000_lowres | 11640.825297 | `execution_times_pilotgm_Gongsharma_cmv_young_males.feather` |
| Kfoury | PILOT-GM-VAE_hvg1000_highres | 347.068789 | `execution_times_pilotgm_Kfoury.feather` |
| Kfoury | PILOT-GM-VAE_hvg2000_highres | 351.430396 | `execution_times_pilotgm_Kfoury.feather` |
| Kfoury | PILOT-GM-VAE_hvg2000_lowres | 278.934396 | `execution_times_pilotgm_Kfoury.feather` |
| Kfoury | PILOT-GM-VAE_hvg3000_highres | 346.250744 | `execution_times_pilotgm_Kfoury.feather` |
| Kim | PILOT-GM-VAE_hvg1000_highres | 1487.952933 | `execution_times_pilotgm_Kim.feather` |
| Kim | PILOT-GM-VAE_hvg2000_highres | 1508.293093 | `execution_times_pilotgm_Kim.feather` |
| Kim | PILOT-GM-VAE_hvg2000_lowres | 910.281262 | `execution_times_pilotgm_Kim.feather` |
| Kim | PILOT-GM-VAE_hvg3000_highres | 1554.172378 | `execution_times_pilotgm_Kim.feather` |
| Pelka | PILOT-GM-VAE_hvg1000_highres | 7385.304739 | `execution_times_pilotgm_Pelka.feather` |
| Pelka | PILOT-GM-VAE_hvg2000_highres | 7377.742997 | `execution_times_pilotgm_Pelka.feather` |
| Pelka | PILOT-GM-VAE_hvg2000_lowres | 2304.164306 | `execution_times_pilotgm_Pelka.feather` |
| Pelka | PILOT-GM-VAE_hvg3000_highres | 7371.934789 | `execution_times_pilotgm_Pelka.feather` |
| Stephenson | PILOT-GM-VAE_hvg1000_highres | 2400.978251 | `execution_times_pilotgm_Stephenson.feather` |
| Stephenson | PILOT-GM-VAE_hvg2000_highres | 2405.345404 | `execution_times_pilotgm_Stephenson.feather` |
| Stephenson | PILOT-GM-VAE_hvg2000_lowres | 1866.583914 | `execution_times_pilotgm_Stephenson.feather` |
| Stephenson | PILOT-GM-VAE_hvg3000_highres | 2284.238830 | `execution_times_pilotgm_Stephenson.feather` |
| Wu | PILOT-GM-VAE_hvg1000_highres | 526.376027 | `execution_times_pilotgm_Wu.feather` |
| Wu | PILOT-GM-VAE_hvg2000_highres | 521.818132 | `execution_times_pilotgm_Wu.feather` |
| Wu | PILOT-GM-VAE_hvg2000_lowres | 358.570813 | `execution_times_pilotgm_Wu.feather` |
| Wu | PILOT-GM-VAE_hvg3000_highres | 519.055718 | `execution_times_pilotgm_Wu.feather` |

The required default-method consequence is 11 production `PILOT-GM-VAE_hvg2000_highres` rows: the four already merged rows plus these seven additions. If a post-gate run-owned log supplies a newer value for any of the 27 keys, use that row and record the run ID instead of the retained value. If a retained log lacks a checksum sidecar, do not invent a value or block the primary repair: preserve its exact finite row, mark it `UNVERIFIED_RETAINED_TASK_LOG`, and document the source path and timestamp in the reconciliation manifest.

### 3. Trace Bassez against the rolling-pipeline evidence, but default to preserving it

1. After the active Bassez gate is terminal and reviewed, inspect its immutable selection manifest, worker logs, execution-log Feather, terminal audit, reviewer evidence, and selected output checksums. Match Bassez rows by exact `(dataset, method)` and verify that the runtime source and result artifact belong to the same run/source identity.
2. Compare only Bassez rows with explicit current-run provenance. The local retained Bassez log snapshot is not sufficient: all 118 Bassez method rows differ from the merged file, and the merged Feather has no per-row run ID or timestamp.
3. If the reviewed B5-BASSEZ or corresponding reviewed benchmark-rest evidence identifies a newer Bassez timing row, replace only those exact keys and record the evidence path/run ID. If the evidence is missing, mixed, or only consists of old retained logs, leave all Bassez values unchanged and record `BASSEZ_UNRESOLVED_LEFT_UNCHANGED`.
4. Do not use filename order, `drop_duplicates(..., keep="last")`, or file modification time alone as a newest-runtime policy. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py:264-279` is not sufficient for this historical reconciliation because its ordering is task-log path order.

### 4. Install the candidate atomically without invalidating gate checks

1. Construct the candidate from the current merged table plus the approved additions/replacements; do not delete the 15 `_debug` rows that have no retained task log.
2. Validate exact columns `dataset`, `method`, `time_secs`, `mem_GB`, nonblank identifiers, finite `time_secs`, finite-or-NA `mem_GB`, and unique `(dataset, method)` keys. With Bassez unchanged, the expected row count is 1,403; Bassez replacements do not change that count.
3. Write the candidate Feather to a temporary same-filesystem path, generate a matching `.md5` sidecar, and atomically replace both under an exclusive local reconciliation lock. Preserve the existing canonical `PATH=` spelling in the sidecar.
4. Update only the `embeddings/execution_times.feather` MD5 line in `data/benchmark/checksums.md5`; preserve all other artifact hashes byte-for-byte. Verify the Feather sidecar and aggregate manifest after installation. The existing synchronization contract in `benchmark_submit_common.sh:980-1003` and `:1005-1063` is the reference for sidecar/manifest ordering, but it must not be invoked against the active pipeline.
5. Do not modify any H5AD, result RDS, embedding Feather, `datasets.json`, `pixi.toml`, or `pixi.lock`. Do not commit, push, pull, or synchronize the repaired local file to HPC in this change.

### 5. Add a narrowly scoped prevention fix after the data repair

The immediate data repair must not wait for source changes. After the Feather is repaired, make the timing recorder replay durable cached timings without recomputing methods:

1. Python benchmark cache hits currently validate a Feather and `continue` at `1.1.1_benchmark_methods_py.py:992-999`, returning before `log_execution_time()` at `:1122-1125`. Add an atomic runtime metadata sidecar tied to the output MD5. On a valid cache hit, validate the sidecar and re-emit its `time_secs`/`mem_GB`; if the sidecar is absent or mismatched, fail closed instead of silently skipping the timing row. Include the sidecar in the guarded benchmark sync/checksum selection.
2. Pipeline-B `trans`/`zeroimp` has the same gap: `benchmark_hpc_utils.R:713-716` exits on a valid RDS before the only `log_exec_row()` call at `:752-755`. Use the same sidecar contract for `<dataset>_trans.rds` and `<dataset>_zeroimp.rds`; on cache reuse, validate the output checksum and sidecar, then call `log_exec_row()` before exiting. Do not mutate the trans/zeroimp RDS result shape solely to store timing metadata.
3. Preserve the existing R benchmark-bundle replay behavior (`benchmark_pipeline.R` cached branches already call `log_exec_row()` from stored `exec_time`).
4. Add focused local tests for: Python cache-hit replay, Python missing/mismatched metadata fail-closed behavior, Pipeline-B trans/zeroimp cache-hit replay, and final merge rejection of duplicate/invalid runtime keys. These tests must use temporary fixtures and must not invoke Slurm, SSH, NAS synchronization, or pipeline submitters.

## Critical files & anchors

- `.agents/plans/20260830170000-bassez-rolling-pipeline-rollout-plan.md:88-98,449-475,630-645` — Bassez B5 method selection, PILOT-GM-VAE timeout evidence, and active promotion wave.
- `data/benchmark/embeddings/execution_times.feather` plus `.md5` and `data/benchmark/checksums.md5` — authoritative merged runtime table and gate-tracked hashes.
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:992-999,1104-1125` — Python output cache skip and runtime logging.
- `src/5_run_benchmark_methods/benchmark_hpc_utils.R:682-759` — Pipeline-B trans/zeroimp cache exit and runtime logging.
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py:264-283` — task-log merge, deduplication, and existing-log precedence.

## Verification

1. Before installation, read the baseline Feather and assert 1,376 rows, exact schema, no duplicate keys, and the recorded baseline MD5. Save a candidate diff listing only the 27 PILOT-GM-VAE additions plus any explicitly approved Bassez replacements.
2. After installation, run a Pixi Python validation that asserts:
   - 1,403 rows when Bassez is left unchanged;
   - all 27 listed task-only keys exist with the exact seconds above and source `mem_GB` values copied from their task logs;
   - 11 production `PILOT-GM-VAE_hvg2000_highres` rows exist;
   - `Gongsharma_cmv_young_males` is retained as the raw dataset key;
   - exact four-column schema, no blank identifiers, finite runtimes, and unique keys;
   - all non-runtime artifact hashes in `data/benchmark/checksums.md5` are unchanged;
   - the new Feather MD5 matches both its sidecar and the aggregate checksum entry.
3. If Bassez provenance is unresolved, assert all 118 pre-existing Bassez rows remain byte/value-identical and record the unresolved status in the reconciliation manifest. If provenance is resolved, assert only the explicitly selected Bassez keys changed and record their run evidence.
4. Run the focused cache-replay tests from Step 5 locally with temporary fixtures. Confirm no Slurm command, SSH command, submitter, durable gate, commit, push, or pull is invoked.
5. Confirm the working tree changes are limited to the requested runtime artifact/traceability plan and any explicitly approved recorder source/tests; preserve all unrelated pre-existing changes.

## Assumptions & contingencies

- The active Bassez rolling pipeline remains authoritative until its gate reaches terminal reviewed status. No local runtime replacement is installed while it can still synchronize the same benchmark namespace.
- If the active gate produces newer reviewed PILOT-GM-VAE logs for the seven missing datasets, those rows supersede the retained August task-log values. If not, the exact retained values in Step 2 are the approved fallback and are marked as retained-log evidence.
- If no reviewed Bassez run evidence ties a newer row to the current output artifact, Bassez remains unchanged; the primary PILOT-GM-VAE additions proceed independently.
- The repair is deliberately local and uncommitted. A later operator may synchronize it only after the active pipeline's gate and reviewer evidence are complete; that synchronization is not part of this plan.
