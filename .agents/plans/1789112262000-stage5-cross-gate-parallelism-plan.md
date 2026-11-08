# Stage 5 Cross-Gate Parallelism Plan (Draft)

Status: Draft; architecture and implementation contract require review before code changes.

## Problem

ECODA is intended to run independent datasets and methods concurrently. The existing Stage 5 submitter already dispatches independent rows through SLURM arrays and per-method watchdogs. The remaining limitation is concurrency across **separate durable gates** that target the same analysis root, especially:

```text
HPC_SCRATCH_DIR/batch_effect/uncorrected
NAS_TARGET_DIR/batch_effect/uncorrected
```

`src/5_run_benchmark_methods/benchmark_submit_common.sh` currently acquires one global Stage 5 synchronization owner at:

```text
stage5/sync/${LOCAL_ROOT}
```

(`analysis_merge_sync_cleanup`, current owner acquisition near line 1107). Two same-root gates therefore cannot safely reach synchronization concurrently. Using different durable serialization groups would bypass the coordination rather than solve it: workers could write the shared analysis tree, `checksums.md5`, execution-time merges, cleanup manifests, and NAS outputs concurrently.

Current behavior is fail-closed, but it serializes otherwise independent gates. Kidney Stage 5 is prepared and queued behind the active Lupus/Parkinson same-root recovery gate; it has not been launched.

## Findings and Evidence

- Independent datasets/methods in one Stage 5 selection manifest are already represented as array rows and method-specific watchdogs.
- Method dependencies are explicit: `prepare_pseudobulk` precedes composition/pseudobulk-related consumers through watchdog dependencies.
- Stage 5 output ownership is per method/dataset, but final synchronization currently converges on a shared root owner.
- `analysis_merge_sync_cleanup()` writes selected artifacts to the shared local/NAS roots and manages shared synchronization state.
- Durable gates use the `ecoda-benchmark` serialization group. This currently prevents two same-root gates from launching concurrently.
- A distinct serialization group is not a safe workaround; the shared owner and shared checksum/NAS merge remain the real synchronization boundary.
- Existing validated artifacts must remain outside every recovery selection. Any redesign must preserve targeted selection, fail-closed checksums, and run-owned provenance.

## Supported Path Without Refactor

Use one explicit Stage 5 selection manifest containing all independent dataset/method rows that share an analysis root. SLURM arrays and method watchdogs provide the intended parallelism; one gate owns the final shared synchronization boundary.

This is the current safe path for Kidney, Lupus, Parkinson, and independent methods if they are selected together in a new approved gate. Separate same-root gates remain serialized until this plan is implemented.

## Proposed Solution

Replace the global same-root synchronization owner with a two-level contract:

1. **Fine-grained dataset/method ownership**
   - Acquire ownership for the exact dataset/view/method artifact set, not the whole `ANALYSIS_ROOT`.
   - Keep existing artifact-level owners and run-owned records authoritative.
   - Reject overlap when two gates target the same dataset/view/method or the same mutable artifact.

2. **Serialized, transactional shared-output merge**
   - Workers write immutable per-run sync fragments/manifests under their run root.
   - A merge coordinator serializes only the shared mutable operations:
     - NAS copy/update for the selected artifact paths;
     - `checksums.md5` merge;
     - execution-time aggregate merge;
     - cleanup of obsolete files.
   - Merge inputs are explicit selected-row fragments; no broad filesystem scan.
   - Use an atomic lock/queue with durable ownership and stale-run recovery. A merge failure must leave validated existing outputs intact and fail closed.

3. **Durable-gate dependency model**
   - Independent compute gates may launch concurrently when their dataset/method ownership sets are disjoint.
   - Each gate performs its own terminal scheduler/accounting audit.
   - A final merge gate depends on all selected compute gates and is the sole owner of shared checksum/NAS consolidation for that analysis root.
   - Reviewer approval remains required before dependent pipeline stages proceed.

4. **Recovery semantics**
   - OOM retry remains method/manifest-row scoped.
   - Failed rows release only their exact dataset/method owners.
   - Successful rows and their sync fragments remain immutable and outside later repair selections.
   - Re-running a failed merge must reuse validated fragments and must not rerun workers.

## Implementation Plan

### Phase 1 — Contract and ownership design

- Define canonical owner keys for dataset/view/method and merge-root/pass.
- Define per-run sync fragment schema: artifact path, dataset, view, method, MD5, size, producer run, source/runtime identity.
- Define merge lock/queue state transitions and stale-owner recovery rules.
- Decide whether the final merge coordinator is a durable gate wrapper or a validator/merge-only worker with an explicit dependency manifest.
- Record the contract in code comments and this plan before implementation.

### Phase 2 — Refactor Stage 5 synchronization

- Update `benchmark_submit_common.sh` to stop acquiring one global owner for all same-root dataset rows.
- Emit and validate per-dataset/method sync fragments from watchdogs.
- Implement the serialized merge coordinator with atomic checksum and execution-time updates.
- Preserve NAS path validation, checksum parity, artifact records, and no-delete-on-failure behavior.
- Keep existing single-gate behavior compatible during migration where possible.

### Phase 3 — Durable gate integration

- Add explicit compute-gate and merge-gate dependency edges.
- Update emitted scheduler manifests and terminal audit contracts to include all compute, watchdog, and merge IDs.
- Ensure same-root disjoint compute gates can use distinct serialization groups only after ownership/merge contracts prove disjointness.

### Migration invariant

Distinct serialization groups must **not** invoke the current
`analysis_merge_sync_cleanup()` implementation. Before cross-gate compute
parallelism is enabled, the submitter must support a workers-only/deferred-sync
mode that emits immutable run-owned sync fragments and does not mutate shared
`execution_times.feather`, `checksums.md5`, or NAS state. Only the serialized
merge coordinator may perform those operations.

Until that mode exists, same-root gates remain in the `ecoda-benchmark`
serialization group and are queued rather than launched concurrently.

## Assumptions & contingencies

- `datasets.json`, `pixi.lock`, existing benchmark artifacts, and current gate
  evidence remain read-only for this refactor.
- The current Lupus/Parkinson remote run may finish after its local durable
  transport failure, but its existing gate is not releasable from that failed
  inspect. Recovery must use a new explicit validator/merge design; never
  reset or hand-edit its lifecycle state.
- Kidney Stage 5 remains queued until the current same-root synchronization
  boundary is terminally resolved. A smaller cohort does not remove the shared
  owner/checksum/NAS dependency.
- If the profile or submitter cannot represent workers-only/deferred-sync
  execution, stop before compute launch and generalize the wrapper contract
  rather than bypassing the durable gate.
- Update snapshot-backed wrapper examples and run reports.

### Phase 4 — Focused validation

Add deterministic tests for observable contracts:

- Two disjoint datasets in the same analysis root can compute concurrently without owner conflict.
- Same dataset/view/method overlap is rejected fail-closed.
- Concurrent fragments merge atomically into one checksum manifest.
- A failed merge preserves existing validated outputs and allows merge-only recovery.
- A stale merge owner is reclaimed only with terminal evidence.
- OOM retry affects only the failed manifest rows.
- Existing valid artifacts outside the selected rows remain unchanged.
- Durable terminal audit rejects missing scheduler IDs, checksum mismatch, malformed fragments, and incomplete merge state.

Use shell stubs and temporary roots; do not launch full cohorts for these tests.

### Phase 5 — Rollout and cleanup

- Run the focused contract tests and `_debug` path.
- Validate one explicitly scoped production recovery before enabling concurrent same-root gates broadly.
- Update operator documentation and the snapshot-backed workflow examples.
- Re-read `AGENTS.md` after implementation and update the Stage 5 concurrency boundary section to describe the implemented contract, supported launch topology, and any remaining limitation. The current `AGENTS.md` note is a temporary limitation record and must be checked against the final implementation rather than left stale.
- Preserve the current serialization behavior until the new merge contract is terminally validated.

## Acceptance Criteria

- Independent same-root dataset/method gates can run compute concurrently without a false global sync-owner conflict.
- Shared checksums/NAS state has one durable, atomic merge owner or queue.
- No worker reruns are needed for merge-only recovery.
- All scheduler IDs, artifacts, checksums, source/runtime identities, and reviewer approvals remain auditable.
- Existing validated outputs outside the explicit selection are never invalidated or overwritten.
- A documented fallback remains: one combined selection manifest/array for safe parallel execution.

## Current Run Boundary

- Lupus/Parkinson consumer recovery is the active same-root Stage 5 gate.
- Kidney Stage 5 has a corrected prepared gate and remains unlaunched until the active gate reaches terminal wait, inspection, and Luna Max review.
- The earlier Kidney prepared manifest with the incorrect runtime-manifest path is preserved as an unlaunched preparation error; do not launch it.
- Approved high-resolution recovery scope remains exactly six rows:
  `composition`, `pilot`, and `qot` for `Lupus_PBMC` and `Parkinson`.
  In uncorrected mode MrVI receives no cell-type column; `cell_type_high_res`
  is consumed by PILOT/QOT (and composition's metadata path), not MrVI.
  The submitted MrVI array `4399479` is extra evidence from an earlier
  over-broad recovery attempt and must not be counted as an approved
  affected-row recovery. Future recovery selections must exclude `mrvi`.
