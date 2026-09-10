# Pseudobulk calculation deduplication plan

## Context

The canonical Stage 5 pipeline currently has two different pseudobulk boundaries:

1. `src/utils/py/h5ad_pseudobulk.py` streams persisted CSR counts by `Sample`, creates a dense int64 sample-by-gene aggregate, and returns genes-by-samples data.
2. `src/5_run_benchmark_methods/benchmark_hpc_utils.R` transfers that aggregate through reticulate into a one-sample-per-column Seurat object, after which `src/utils/pseudobulk.R:get_pb()` calls `AggregateExpression()` again.

The legacy Seurat path instead materializes full cell-by-gene counts before calling `AggregateExpression()`. It is not the target production path.

Cell-type pseudobulks are separate: `process_pseudobulk_ct_fig()` currently loops over cell types, subsets the full Seurat object, applies the five-cell `Sample × cell_type` eligibility rule, aggregates and DESeq2-normalizes each subset, computes Euclidean distances, and averages each sample pair over successful cell types.

The goal is to remove redundant representations and repeated work without changing the benchmark method. Existing benchmark H5ADs, pseudobulk caches, result bundles, checksums, manifests, and gates are immutable and must not be overwritten or invalidated by this refactor.

## Decisions locked during design review

- Preserve current scientific and numerical semantics as closely as possible.
- Update only canonical maintained Stage 5 production paths. Do not spend the main refactor migrating historical notebook/Seurat code. Remove legacy pseudobulk code only after a maintained-caller audit proves it is unused.
- Remove the sample-level Seurat round trip from canonical preparation, fallback, and ordinary pseudobulk execution.
- Keep two explicit input contracts:
  - `Sample -> raw gene counts` for ordinary pseudobulk.
  - `(Sample, cell_type) -> raw gene counts + group cell counts` for cell-type pseudobulk.
- Preserve the CT contract: first-observation cell-type iteration order, sorted sample distance universe, at least five cells per `Sample × cell_type`, separate DESeq2 normalization per cell type, isolated failed-cell-type handling, and pairwise distance averaging.
- Keep the custom HDF5/CSR implementation as the production baseline. Evaluate Scanpy and decoupler only as candidates; adopt either only if a controlled benchmark shows a clear wall-time or peak-RSS win while preserving the contract.
- Preserve current raw-count, no-label-leakage, cache naming, published matrix orientation, and batch-only correction contracts.
- Use checked int64 accumulation, then reject any final count above the DESeq2/R integer ceiling (`2,147,483,647`) before reticulate/DESeq2. Never rely on unchecked int64 wraparound or arbitrary int64 conversion through the existing DESeq2 path.
- Record shared preparation timing separately from per-variant selection timing. New timing consumers must count shared aggregation/fit work once rather than repeat it for every variant.
- Do not modify `datasets.json`, `pixi.toml`, or `pixi.lock`.

## Exact output and ordering contracts

### Ordinary sample pseudobulk

- Source: `layers["counts"]`, never normalized/log-transformed `X`.
- Grouping: `Sample` only.
- Internal raw orientation: genes × samples, matching the current Python helper.
- Published `$pb` orientation: samples × genes, matching current RDS consumers.
- Sample IDs: first appearance in the H5AD aggregation, then canonical metadata order before publication.
- Metadata: first observation per sample for requested bookkeeping columns; biological labels remain evaluation-only and are never DESeq2 covariates.
- Raw values: finite, nonnegative, integer-valued counts.
- Ordinary DESeq2 settings: design `~1`, `blind=TRUE`, `correct_batch=FALSE`.
- Corrected batch settings: configured batch-only design, `blind=FALSE`, `correct_batch=TRUE`; biological labels remain absent.

### Cell-type pseudobulk

- Source: the same raw CSR count layer plus cell-level `Sample` and configured cell-type metadata.
- Grouping: present `(Sample, cell_type)` combinations only; do not fabricate the full Cartesian product.
- Eligibility: count raw cells in each `(Sample, cell_type)` group and retain groups with at least five cells. A cell type with fewer than two eligible samples is skipped.
- Cell-type order: first-observation order from the source metadata, excluding missing cell-type values.
- Final sample distance universe: `sort(unique(Sample))`, as in `process_pseudobulk_ct_fig()`.
- Per-cell-type sample order: canonical first-observation order before distance placement.
- Per-cell-type failure: skip that cell type, retain successful-cell-type counters, and fail closed only if no successful cell type contributes a sample pair.
- Distance result: preserve current Euclidean distance calculation and pairwise average denominator.

### Count and transport safety

- Validate negative values for every numeric CSR dtype, including signed integer data. The current integer branch does not reject negative values and must be corrected.
- Use checked int64 updates. Before every accumulator addition, verify `incoming <= INT64_MAX - current`; an observed maximum after unchecked `+=` is not an overflow check.
- The Stage 5 DESeq2 path must request an output ceiling of `INT_MAX = 2,147,483,647` from the Python reducer so unsafe values are rejected before crossing reticulate. The generic raw reducer may retain int64 output only when its caller does not cross into DESeq2.
- Keep int64 as the authoritative aggregation type. Conditional int32 is not part of the first refactor: post-aggregation downcasting cannot reduce peak accumulator memory, and int32 accumulation adds overflow branches without removing the DESeq2 ceiling.
- If an int64 cumulative sum would overflow, fail closed; never wrap, saturate, or continue with corrupted counts.

## Backend evaluation

Run a focused, local benchmark before selecting a library implementation. Do not launch full cohorts or change the pinned environment merely to test a candidate.

Compare:

1. The current custom HDF5/CSR reader.
2. A revised custom checked CSR reducer with no dense per-chunk grouped matrix.
3. `scanpy.get.aggregate()` with `layer="counts"` and sparse input.
4. `decoupler.pp.pseudobulk()` with `layer="counts"`, `mode="sum"`, `groups_col=None` and the cell-type grouping mode.

Use the `_debug` fixture plus deterministic synthetic sparse H5ADs that vary cell count, gene count, sparsity, unsorted sample order, and Sample × cell-type group occupancy. Measure wall time, peak RSS, output dtype, output density, and exact contract equality.

Evidence already reviewed:

- [Scanpy aggregate API](https://scanpy.readthedocs.io/en/stable/generated/scanpy.get.aggregate.html): supports sparse CSR/CSC input, multiple grouping columns, and experimental Dask arrays, but does not promise this repository's bounded H5AD behavior, first-observation metadata, or canonical ordering.
- [decoupler pseudobulk API](https://decoupler.readthedocs.io/en/latest/api/generated/decoupler.pp.pseudobulk.html): directly supports sample-only and Sample × group sums, layer selection, and `bsize`; its current implementation accepts sparse input but allocates dense floating-point profile/output arrays and can retain absent grouped combinations.

Neither library is a safe drop-in without adapters. Adopt one only if it is contract-equivalent and measurably faster or lower-RSS than the revised custom reducer. Otherwise keep the custom implementation and do not add a dependency.

## Implementation approach

### 1. Establish the oracle and maintained-caller inventory

- Add or extend a deterministic synthetic CSR H5AD oracle before changing the reducer.
- Cover unsorted sample order spanning chunk boundaries, unsorted cell-type order, a `Sample × cell_type` group below five cells, absent combinations, a failed/empty cell-type normalization path, and every ordinary pseudobulk variant.
- Compare old and new raw aggregates, normalized matrices, sample IDs, cell-type counters, distance placement, and cache payload fields. Existing stubs that only test forwarding are insufficient for numerical equivalence.
- Before the caller inventory, query the local CodeGraph resource. If no CodeGraph MCP resource is available, record that fact and use language-server references plus targeted source reads; deletion still requires an explicit maintained-caller audit.

### 2. Refactor the Python CSR reducer

Target: `src/utils/py/h5ad_pseudobulk.py`.

- Factor the shared H5AD/CSR validation and checked group-accumulation logic without changing the public sample-only result contract.
- Keep `aggregate_h5ad_counts_by_sample()` as the sample-only entry point, retaining genes × samples int64 output, first-seen sample IDs, first-observation metadata, missing-ID failures, and `chunk_size` behavior.
- Replace `grouped.toarray()` as the default accumulation mechanism. The selected reducer must avoid a dense `U × n_genes` per-chunk intermediate: prefer direct guarded CSR-row updates, or retain grouped chunks as sparse CSR and update only their stored nonzero positions.
- If sparse membership multiplication remains selected, apply a conservative per-chunk bound before multiplication. For every local group, prove `max_raw_value * n_cells_in_group <= limit` using division-based checks that cannot overflow, where `limit` is `INT64_MAX` for generic aggregation or the requested DESeq2 ceiling for Stage 5.
- Before every update to the persistent accumulator, verify `incoming <= limit - current` and only then add the values. Never rely on observing a maximum after an unchecked `+=`.
- Add a composite Sample × cell-type entry point backed by the same reducer. Return present groups, raw group cell counts, group metadata, exact group order, and a bounded representation that remains sparse until a per-cell-type DESeq2 call requires a dense matrix.
- Support an optional maximum output value so Stage 5 can reject values above `INT_MAX` before reticulate. Keep generic raw aggregation capable of returning int64 values when no DESeq2 boundary is requested.
- Add explicit negative-count rejection for signed integer CSR data and checked cumulative-overflow regression coverage.

### 3. Add a direct matrix/metadata DESeq2 boundary

Target: `src/utils/pseudobulk.R`.

- Keep `DESeq2.normalize()` as the scientific core, but separate fitting/transformation from final HVG selection so one fitted normalized/VST matrix can serve multiple variants.
- Add a direct raw matrix plus metadata path that accepts the genes × samples aggregate without constructing Seurat or calling `AggregateExpression()`.
- At this R boundary, independently validate finite, nonnegative, integer-valued counts and require every value to be `<= .Machine$integer.max` before constructing `DESeqDataSetFromMatrix`; then coerce through the supported R integer representation. This duplicate guard protects legacy adapters and future callers that bypass Python.
- Preserve the current VST fallback chain, batch-only correction, no-label-leakage rules, sample alignment, and published samples × genes orientation.
- For `hvg500`, `hvg1000`, `hvg2000`, and `hvg3000`, run the common full-gene DESeq2 size-factor/VST/variance calculation once and derive the requested prefixes from one variance ordering.
- Keep `schvg2000` separate because it restricts the raw gene universe before DESeq2. Preserve the current `hvg2000_bl` behavior during this optimization; do not combine a blacklist correction with the performance refactor.
- Enforce the `INT_MAX` gate before DESeq2. Do not pass arbitrary int64 totals through reticulate and rely on DESeq2 to coerce them.
- Retain `get_pb()`/`get_pb_deseq2()` only as an explicitly legacy Seurat boundary until the caller audit and cleanup phase. Canonical Stage 5 code must not call them.

### 4. Replace sample-level preparation and fallback integration

Targets:

- `src/5_run_benchmark_methods/benchmark_hpc_utils.R`
- `src/5_run_benchmark_methods/benchmark_pipeline.R`
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R`
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R`

- Add a raw aggregate/metadata loader that calls the Python reducer and hands the result directly to the matrix/metadata DESeq2 boundary.
- In the preparation worker, validate the H5AD contract, read metadata/HVG ranks, and compute the pending cache set before touching the counts layer. If every requested variant is checksum-valid and `--force` is absent, re-emit cached timings and perform no raw aggregation, sample-level object construction, or count-layer read. Aggregate only for missing variants or an explicit force.
- Remove `load_h5ad_pseudobulk_seurat()` from canonical preparation and missing-cache fallback. Do not create a one-sample-per-column Seurat object merely to invoke `AggregateExpression()`.
- Preserve cache names, cache stems, producer identities, checksum records, atomic publication, missing-only resume, `--force`, and batch-mode selection.
- Preserve the existing ordinary six-variant and batch-only `hvg2000` selection behavior.
- Make `prepare_pseudobulks_hpc()` perform one raw Sample aggregation and one shared full-gene DESeq2 fit where mathematically valid, then publish each variant from the shared fit. Keep `schvg2000` as a separate fit.
- Update missing-variant fallback in `load_pb_variants()` to use the H5AD raw aggregate path rather than requiring a Seurat object. Cache validation must precede aggregation there as well.
- Keep composition counts-free and ensure ordinary pseudobulk no longer loads full counts solely because CT pseudobulks used to require Seurat.
- Retain the full count-backed Seurat loader only for methods that genuinely need cell-level counts, such as scITD.

### 5. Define honest shared timing and memory records

- Keep the existing `$pb`, `$time_secs`, and `$mem_GB` fields for cache readers, and add a timing schema for newly produced records:
  - `aggregate_time_secs`
  - `shared_fit_time_secs`
  - `shared_time_secs`
  - `variant_time_secs`
  - `shared_mem_GB`
  - `timing_id`
  - `timing_schema`
- Define `shared_time_secs` exactly as `aggregate_time_secs + shared_fit_time_secs`. `variant_time_secs` includes only variant-specific work: separate `schvg2000` fitting when applicable, variant selection, and publication/cache work; common full-gene fitting is excluded.
- Define `timing_id` as the run-scoped identity of one shared aggregate/fit: `${ECODA_RUN_ID}:${cache_stem}:${view}:${analysis_pass_or_none}`. All variants produced from that aggregate/fit carry the same ID; different source snapshots, modes, or runs must not share an ID.
- For timing schema 2, set `$time_secs` to `variant_time_secs` for canonical consumers; legacy records without the schema retain their old inclusive interpretation. New consumers must use `shared_time_secs` once per `timing_id` and never sum it once per variant.
- Emit one `prepare_pseudobulk_shared` execution-log row keyed by `timing_id` and one variant-local row per requested variant. `run_mofa_hpc()` and `run_pseudobulk_hpc()` use variant-local time for combo results; method/report aggregation adds the shared row exactly once per timing ID.
- Before changing these fields, inventory every canonical consumer of `$time_secs`, `$mem_GB`, `log_exec_row()`, and `exec_time` (including preparation, fallback, method replay, report, merge, and validator paths), then migrate all of them to schema 2. No consumer may silently double-charge or omit shared time. Keep legacy cache records readable. If a cached record lacks the new timing fields, treat its existing `time_secs` as legacy-inclusive timing and do not fabricate a shared-time decomposition.
- Record `shared_mem_GB` as a shared-stage peak or upper bound, not an isolated component allocation: `peak_rss_gb()`/VmHWM is process-cumulative. If component-level memory is required, measure aggregation and DESeq2 phases in isolated subprocesses. Do not claim that sparse aggregation eliminates the unavoidable samples × genes dense DESeq2 allocation.

### 6. Replace repeated CT Seurat subsets

Target: `src/5_run_benchmark_methods/benchmark_methods_r.R` and the Stage 5 orchestration in `benchmark_pipeline.R`.
- Add a canonical CT routine that consumes the composite Sample × cell-type aggregate and raw group cell counts rather than a full Seurat object.
- During a metadata pass, assign stable first-seen IDs to present `(Sample, cell_type)` groups and determine five-cell eligibility before reading count values.
- During the single count pass, append sparse grouped contributions to run-owned temporary HDF5/CSR datasets keyed by composite group ID: contribution group IDs, CSR indptr, indices, and data plus group metadata. The store is append-only/chunked and bounded; do not retain an in-memory dictionary of all `Sample × cell_type × gene` values. Benchmark reducer/storage variants before implementation, but do not introduce a second-pass fallback that changes the one-pass contract. Create it below `<run-owned scratch>/pseudobulk_ct/<run_id>/<unique-token>`, acquire it exclusively, and write an ownership manifest containing run ID, PID/scheduler identity, source H5AD identity/checksum, stage, schema, and creation time. The manifest and path must make an interrupted store impossible to mistake for a canonical artifact.
- For each cell type, scan only its eligible contribution rows, merge them into one sparse per-CT matrix with checked int64 additions, materialize only that CT matrix at the DESeq2 boundary, normalize it with the existing per-CT settings, compute distances, and release it before the next CT where possible. Close/flush all HDF5 handles before cleanup, run cleanup from `finally`/`on.exit` on both success and ordinary failure, and require the temporary root to be absent before publishing a canonical result. Because `finally`/`on.exit` cannot run after SIGKILL/OOM, perform a no-compute stale-store audit/cleanup before reuse: remove only stores with a valid matching ownership manifest, a demonstrably dead owner, and an expired age/lock policy; never remove active or ambiguous stores, and fail closed on uncertainty. Cleanup failure is fatal; atomic rename is reserved for canonical artifact publication.
- Preserve `successful_cell_types`, `n_sample_pairs_contributed`, `n_ct_pair_contributions`, sorted final sample dimnames, first-seen CT order, per-CT error isolation, and the fail-closed all-failed guard.
- Change the ordinary pseudobulk worker to load only required metadata plus the H5AD path for CT aggregation. It must not materialize full counts after this change.
- Keep the full Seurat CT function only while the legacy caller audit requires it; remove it after canonical callers have migrated and focused tests cover the new path.

### 7. Legacy cleanup after caller audit

- Do not edit historical notebooks merely to make them use the new path.
- After canonical migration, confirm whether `run_benchmark_analysis()` and direct `get_pb()`/`get_pb_deseq2()` callers are maintained.
- If no maintained caller remains, remove the obsolete Seurat pseudobulk wrappers and their obsolete tests/docs. Preserve all existing artifacts and historical commit history.
- If a maintained caller remains, leave the legacy boundary isolated and documented rather than silently deleting it.

### 8. Documentation

Update `docs/ARCHITECTURE.md` and relevant source comments to document:

- raw CSR counts as the only pseudobulk source;
- separate Sample-only and Sample × cell-type contracts;
- direct matrix-to-DESeq2 processing;
- one shared DESeq2 fit for full-gene HVG variants;
- checked int64 accumulation plus the DESeq2 `INT_MAX` boundary;
- honest shared versus variant-local timing;
- the fact that Scanpy/decoupler are optional only if measured winners.

Do not document a library as adopted until the benchmark proves the required contract and performance win.

## Verification

No full-cohort HPC run is part of this refactor. Existing validated production artifacts remain untouched. Use focused local checks and `_debug` only.
1. **Python reducer contract — `tests/test_h5ad_pseudobulk.py`**
   - Exact sums on synthetic CSR H5ADs with samples crossing chunk boundaries.
   - First-seen sample order and first-observation metadata.
   - Rejection of negative signed-integer values, nonfinite values, noninteger floats, missing IDs, and cumulative int64 overflow.
   - Rejection before return/transport when a DESeq2-bound call exceeds `INT_MAX`.
   - Verification that the selected reducer proves the per-chunk bound before sparse multiplication, when applicable, and performs a checked guard before every persistent update.
   - No dense `U × n_genes` grouped chunk allocation in the selected implementation.

2. **Legacy-versus-new sample oracle — `tests/test_h5ad_pseudobulk.py`**
   - Build the same synthetic data through the old Seurat aggregation and the new direct matrix path.
   - Compare raw sums, normalized matrices within a documented numeric tolerance, gene selections, sample IDs/order, and all six variant names.
   - Exercise ordinary, uncorrected batch, and corrected batch parameter combinations.

3. **Shared-fit oracle — `tests/test_batch_effect_correction.R`**
   - Verify that full-gene `hvg500/1000/2000/3000` outputs match independent current fits within the declared tolerance.
   - Verify that `schvg2000` remains a separate pre-filtered fit.
   - Verify that valid existing cache records are reused without recomputation and that missing-only fallback computes only missing variants.
   - Call the direct R matrix/metadata helper with a value of `.Machine$integer.max + 1` and assert rejection before `DESeqDataSetFromMatrix`.
   - Verify timing schema 2 fields, exact shared-time sum, one shared timing identity, and no duplicate shared-time charge in new summaries.

4. **CT oracle — `tests/test_batch_effect_correction.R` or a new focused R test**
   - Use unsorted Sample × cell-type metadata with a group below five cells, absent combinations, and a deliberately failed CT normalization.
   - Compare eligible sample IDs, CT order, successful/failure counters, sorted final distance dimensions, pairwise denominators, and final distances against the existing routine.
   - Confirm that the ordinary pseudobulk worker no longer loads a full count-backed Seurat object solely for CT calculations.
   - Confirm the temporary HDF5/CSR group store is unique/run-owned, bounded, closed and removed on both success and failure, and never published as a canonical artifact.

5. **Backend benchmark**
   - Compare the custom revised reducer, Scanpy, and decoupler on `_debug` and deterministic synthetic sparse fixtures.
 - Use a fresh subprocess (or a resettable process with a newly reset RSS boundary) for every warm-up and measured repetition; do not collect multiple repetitions in one process because VmHWM is monotonic. Use one warm-up followed by at least five measured repetitions per fixture and backend, with the same environment and fixture order randomized between repetitions. Report median and p95 wall time plus median and maximum per-process peak RSS.
 - Require exact raw aggregate equality, matching dtype/order/metadata/empty-group behavior, and normalized-output equality within the declared tolerance. Define a clear win as at least 10% lower median wall time or peak RSS, with no more than a 5% regression in the other resource metric and no p95 contract/resource outlier.
   - Adopt a library only if it is contract-equivalent and meets that measured threshold; otherwise keep the custom reducer with no dependency change.

6. **Canonical `_debug` smoke path**
   - Exercise `_debug` separately for the exact configured views `benchmark_analysis`, `batch_effect_uncorrected`, and `batch_effect_corrected` from `datasets.json`.
   - Confirm cache files, checksums, result shapes, finite outputs, timing records, and unchanged input H5AD checksums.
   - Do not launch full cohorts or invoke Pipeline 1–5 recovery jobs for verification.

7. **Source/API audit**
 - Run focused language-server diagnostics/references on changed R/Python symbols where available.
 - Inventory and inspect every canonical `$time_secs`, `$mem_GB`, `log_exec_row()`, and `exec_time` consumer, including prep, fallback, method replay, report, merge, and validator paths; verify schema-2 shared timing is counted exactly once in each output.
 - Confirm no canonical Stage 5 caller still constructs the sample-level Seurat round trip or invokes `AggregateExpression()` for the already aggregated sample path.
 - Confirm the run-owned CT store has an ownership manifest and that stale-store cleanup is validator-only, conservative, and fail-closed; confirm legacy deletion, if performed, removed only demonstrably unused code and tests.


## Acceptance criteria

- canonical sample pseudobulk performs one raw Sample aggregation and passes its result directly to DESeq2;
- complete valid variant caches are reused lazily: metadata/HVG/cache validation occurs first, and no counts aggregation or count-backed object is created unless a variant is missing or `--force` is explicit;
- canonical CT pseudobulk performs one composite Sample × cell-type aggregation and preserves all existing CT result semantics;
- no unchecked dense per-chunk grouped matrix remains in the selected reducer;
- sparse multiplication is protected by a pre-multiplication bound and every persistent accumulation is checked before addition;
- the temporary CT group store is bounded, run-owned, cleaned up, and excluded from canonical artifacts;
- cumulative int64 overflow and the DESeq2 `INT_MAX` boundary fail closed;
- full-gene HVG variants share one DESeq2/VST fit where mathematically equivalent;
- published `$pb` orientation, cache names, sample ordering, batch behavior, and result contracts remain compatible;
- shared and variant-local timing are separately recorded and not double-counted in new reports;
- Scanpy/decoupler are adopted only if the controlled benchmark proves a contract-preserving performance win;
- all focused oracle, reducer, CT, backend, and `_debug` checks pass;
- existing production artifacts remain immutable and no unapproved full-cohort compute is launched.
