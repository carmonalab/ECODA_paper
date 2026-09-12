# Multi-batch corrected-mode support plan

## Status and approval boundary

Implementation is approved and complete in this local worktree. No HPC jobs, existing artifacts, gates, or scheduler state were launched or modified by this work.

The current worktree contains user-owned changes in `AGENTS.md`, `NOTES.md`, `datasets.json`, `notebooks/benchmark_analysis.rmd`, `src/utils/constants.R`, the existing derived-analysis plan, and `src/5_run_benchmark_methods/run_local_ecoda_derived.R`. The implementation must not overwrite or clean up those changes.

## Goal

Make future `batch_effect_corrected` processing accept the final scalar or list-valued `columns.batch` definitions in `datasets.json`, while preserving scalar behavior and leaving the completed benchmark and historical `batch_effect_uncorrected` contracts untouched.

The corrected batch method set is exactly:

- `ECODA_authors_HR`
- `ECODA_seuratres_2`
- `ECODA_authors_HR_NULL`
- `Pseudobulk`
- `GloScope`
- `PILOT`
- `MrVI`
- `QOT`

`pilot-gm-vae` is excluded from this work. It must not be added to the corrected suite or used in its verification. Existing ordinary-benchmark references are inspection-only in this scope; no ordinary benchmark artifact is changed or rerun.

## User-confirmed design decisions

1. A list in `columns.batch` represents independent additive technical factors, not a request for one interaction model everywhere.
2. Readers preserve the raw scalar/list/null JSON shape. Execution-boundary code resolves and validates the representation needed by each method.
3. Scalar-only APIs use an ephemeral deterministic composite category. Harmony receives the original ordered list natively.
4. MRVI receives an ephemeral composite column recreated in memory by its worker. It cannot receive a list: pinned scvi-tools `1.5.0.post1` defines `MRVI.setup_anndata(..., batch_key: str | None)` and registers one categorical obs column.
5. The HVG and MRVI paths use separate equivalent in-memory builders, as requested. Their equivalence must be tested against the same fixture and written encoding contract to control drift risk.
6. Every configured batch key must be constant within standardized `Sample`.
7. Missing, blank, exact batch sentinels, or non-estimable batch metadata fails closed. No `Unknown` imputation is allowed.
8. Corrected ECODA uses additive random-effect terms. Corrected Pseudobulk uses one composite end-to-end because the current DESeq2/limma path is scalar-oriented and `removeBatchEffect()` does not provide an arbitrary list interface.
9. `ECODA_authors_HR_NULL` uses the same corrected feature path as `ECODA_authors_HR`; only evaluation labels are shuffled with the existing deterministic null-control behavior.
10. Verification is limited to deterministic synthetic multi-key contracts, validator/no-launch checks, and the scalar `_debug` corrected smoke path. No benchmark or batch-uncorrected jobs are launched, and no old benchmark or batch-uncorrected pipeline is rerun.

## Evidence driving the change

- `datasets.json` currently contains scalar, `null`, and list-valued `columns.batch`. Lists include two-key and three-key cohorts such as Alzheimer, Breast_cancer, Kidney_KPMP_full, and Myocardial_infarction.
- `src/utils/datasets_io.R:1-53` and `src/utils/py/datasets_io.py:19-68` preserve those values unchanged as `batch_col`.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:104-153` passes one value directly to `scanpy.pp.highly_variable_genes(..., batch_key=...)`; the official API documents `batch_key` as `str | None`.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:303-321` calls `harmonypy.run_harmony()` with the current single value. The pinned Harmony source accepts `vars_use` as a string or list; the official Scanpy wrapper documents `key` as `str | Sequence[str]`.
- `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py:628-655` uses `Sample` for uncorrected views, one configured value for corrected views, and does not create a combined column.
- `src/5_run_benchmark_methods/benchmark_methods_r.R:103-165` fits `y ~ 1 + (1 | batch)` for one column and recenters CLR rows.
- `src/utils/pseudobulk.R:10-106` constructs a scalar `~ batch_col` design and calls `limma::removeBatchEffect(x, batch=...)`.
- `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py:1172-1178` forwards one `batch_key` to MRVI; the pinned source confirms one categorical batch field and one `_scvi_batch` tensor.
- `NOTES.md:306-444` defines the two-pass/no-leakage contract and the eight-method corrected suite, but still describes one confirmed technical column and mentions PILOT-GM-VAE in the Harmony group. The plan therefore updates §3.0/§3.2–3.3 explicitly; it does not inherit the old one-key wording or add pilot-gm-vae.
- `notebooks/batch_effect_analysis.rmd` does not exist. The present notebook files are `notebooks/batch_effect_analysis_uncorrected.rmd` and `notebooks/batch_effect_analysis_legacy.rmd`; the uncorrected notebook is explicitly read-only for this task.

### Consumer inventory

The implementation and tests must cover each scalar assumption independently:

| Consumer | Current scalar failure boundary | Planned representation |
|---|---|---|
| `datasets_io.R` / `datasets_io.py` | `batch_col` is exposed unchanged and must not silently change shape | Preserve raw scalar/list/null; validate at execution |
| Preprocessing membership/HVG | list membership and Scanpy `batch_key` are scalar-oriented | Validate all source columns; use ephemeral composite only for HVG |
| Harmony | current caller passes one value, but pinned `vars_use` supports a list | Pass the original ordered validated list |
| R CLR correction | `correct_clr_batch_lmm()` currently fits one `(1 | batch)` term | Fit additive random intercepts for every configured key |
| R Pseudobulk | DESeq2/limma path accepts one batch column | Use the composite scalar for batch removal; corrected DESeq2 design is exactly `~ 1` |
| MRVI loader/setup | loader stringifies a list and MRVI registers one categorical field | Load all source columns, build one in-memory composite, pass its scalar name |
| GloScope/PILOT/QOT | no local batch argument | Resolve only the exact corrected Harmony embedding |

No consumer may infer multi-key safety from another consumer's API.

Official API references:

- https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
- https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html
- https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.MRVI.html

## Contract to implement

### Configuration normalization

- Keep `datasets.json` unchanged.
- Treat a scalar string as one ordered batch key.
- Treat a nonempty list as an ordered batch-key list; a one-element list has the same execution semantics as its one scalar key.
- Reject `null` only for corrected mode, preserving current uncorrected behavior.
- Reject non-string, empty, duplicated, or whitespace-only key names.
- Resolve view-level column overrides before validating the batch list.
- Reject a batch key that equals the configured biological label or standardized sample identifier.
- Preserve the original configured columns in `adata.obs`; do not overwrite them.

### Per-cell/sample metadata validation

Run after view subsetting, sample-name standardization, and the 500-cell filter, before HVG/PCA/Harmony:

- Every declared key must exist in `adata.obs`.
- For configured batch columns only, treat actual missing values, blanks, and the case-insensitive sentinels `NA`, `nan`, `None`, `<NA>`, `n/a`, `null`, and `Unknown` as invalid. This rule must be identical in Python and R; `Unknown` in cell-type annotation columns is outside this rule.
- Every key must be constant within each standardized `Sample`; never let `read_h5ad_sample_metadata()`, `aggregate_h5ad_counts_by_sample()`, or `collapse_sample_metadata()` silently select the first cell.
- Each key must have at least two observed levels for corrected modeling; the combined key must also be estimable.
- Reject disconnected, rank-deficient, near-unique, or otherwise non-estimable correction designs rather than warning and continuing.
- Biological labels remain evaluation-only and never enter HVG selection, normalization, PCA, Harmony, CLR correction, DESeq2 design, or MRVI setup.

### Composite encoding

Define the exact cross-language encoding `ecoda_batch_composite_v1`:

- Normalize `columns.batch` to an ordered key vector. Composite construction is used only for two or more keys; a scalar or one-key list remains a direct column.
- Canonicalize each nonmissing scalar category value to typed Unicode text: character/factor values use `s:` followed by the exact label; booleans use `b:true` or `b:false`; signed integers use `i:` followed by base-10 decimal; finite IEEE-754 binary64 values use `f64:` followed by their 16-lowercase-hex-digit big-endian bit pattern. Other/list/date/object values are rejected. This removes R/Python differences such as `1` versus `1.0`.
- Encode the UTF-8 bytes of each key name and canonical value as lowercase hexadecimal. A token is exactly `ecoda_batch_composite_v1|<key-count>|<pair-1>;<pair-2>;...`, where each pair is `<key-byte-length>:<key-hex>,<value-byte-length>:<value-hex>`. Lengths count raw UTF-8 bytes, not characters; every hex field has exactly two characters per byte; delimiters never occur inside hex.
- Preserve configured key order in every token. Sort categorical levels by the raw UTF-8 token bytes only when constructing category metadata; never sort or reorder configured keys.
- Do not coerce missing values into a category. Record the encoding version, ordered key names, per-key levels, composite-level count, and a stable key-set/configuration fingerprint in run-owned metadata.
- Define the fingerprint input as exact bytes with no JSON/object serialization, locale dependence, or map iteration: `ecoda_batch_contract_v1\0` followed, in this fixed order, by `field("encoding","ecoda_batch_composite_v1")`, `field("keys", ordered_key_vector_v1)`, `field("scalarization", scalarization_id)`, `field("method", method_id)`, and `field("model", model_id)`. `field(name,value)` is `<name-byte-length>:<lowercase-hex-UTF-8(name)>,<value-byte-length>:<lowercase-hex-UTF-8(value)>;`, with lengths counting raw UTF-8 bytes. `ordered_key_vector_v1` is `key-count|<key-byte-length>:<lowercase-hex-UTF-8(key)>;...` in configured order. `scalarization_id` is exactly `direct_v1` for a scalar/one-key direct column or `composite_v1` for two or more keys. `method_id` is the exact canonical corrected method or `preprocess`; `model_id` is a fixed UTF-8 policy token (`hvg_composite_v1`, `harmony_native_list_v1`, `ecoda_additive_random_intercepts_v1`, `pseudobulk_composite_v1`, `mrvi_composite_v1`, or `embedding_consumer_harmony_v1`). Hash the resulting bytes with SHA-256 and record lowercase hexadecimal output. Source H5AD/config checksums remain separate.

Example: keys `site`, `tech` with values `A`, `x` encode as `ecoda_batch_composite_v1|2|4:73697465,3:733a41;4:74656368,3:733a78`.

Use the reserved in-memory obs name `__ecoda_batch_combined_v1`; assert that it is absent before construction and absent again before writing any corrected H5AD. The HVG, MRVI, and R pseudobulk builders may be separate implementations, but must produce identical golden-vector tokens and the same key-set fingerprint for the same scalar, two-key, and three-key fixtures.

### Method mapping

| Method/path | Corrected multi-key representation |
|---|---|
| Scanpy HVG | Build `__ecoda_batch_combined_v1` temporarily and pass that one string as `batch_key`; delete it immediately after HVG ranking. |
| Harmony | Pass the original ordered list of batch columns to the pinned `harmonypy.run_harmony(..., vars_use=list)`/equivalent supported API. Preserve the current orientation normalization and pass-qualified `obsm` keys. |
| ECODA authors, ECODA Leiden res-2, ECODA shuffled null | Use additive per-key random intercepts in the CLR correction model. The null variant shuffles only labels after the same corrected features are produced. |
| Pseudobulk | Build the composite at sample metadata level and pass one scalar composite column through the scalar DESeq2/limma batch-removal path. Corrected DESeq2 design is exactly `~ 1`; keep `blind=FALSE`, `correct_batch=TRUE`, and no biological design protection. |
| GloScope, PILOT, QOT | No method-local batch argument. Resolve the exact corrected Harmony embedding `X_pca_harmony_batch_effect_corrected_hvg2000`; missing keys remain hard errors. |
| MrVI | Load all original batch columns, recreate the equivalent composite in memory, and pass only its scalar name to `MRVI.setup_anndata(batch_key=...)`. Never pass a list and never persist the temporary column. |

### Artifact and compatibility identity

- Do not change benchmark or `batch_effect_uncorrected` H5ADs, RDS bundles, Feather files, checksums, manifests, or gates.
- Preserve existing pass-qualified corrected output filenames and required PCA/Harmony/Leiden keys; the new corrected H5AD contains the original batch columns but not the temporary composite obs column.
- New corrected run manifests and artifact metadata must carry the ordered source-key list, encoding version, key-set fingerprint, sample-constancy result, levels, and correction formula/mode. A checksum alone is not semantic proof of the key set.
- Since no corrected-mode data has been created, do not add a migration or broad invalidation path. Future reuse must require the new key-contract identity; a missing identity is invalid for newly produced corrected rows.
- Keep scalar execution behavior compatible: one scalar configured key remains a direct one-column correction, while list-valued configurations use the method-specific representations above.

## Implementation work packages

### 1. Preprocessing/HVG/Harmony

Target: `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py`

- Add corrected-mode batch normalization and strict metadata validation.
- Keep benchmark and uncorrected branches byte/behavior compatible: benchmark still uses `Sample` for HVG/Harmony, and uncorrected still uses `Sample` with Harmony disabled.
- Extend the view/process interface so corrected HVG receives the temporary scalar key while corrected Harmony receives the original list.
- Ensure the temporary key is removed on both success and exception paths, including the deterministic Seurat-v3 jitter retry.
- Pass list-valued Harmony keys through the pinned multi-key API without changing the existing output shape handling.
- Assert the temporary column is not persisted before `_write_h5ad_atomic()`.

### 2. Python method plumbing and MRVI

Target: `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`

- Resolve scalar/list batch configuration only for corrected mode.
- Load all configured batch columns in the minimal metadata path; do not stringify a list as a column name.
- Validate sample constancy and missing values against the corrected H5AD before any MRVI setup.
- Add the separate equivalent composite builder for MRVI and pass its scalar name to the pinned API.
- Preserve ordinary benchmark and uncorrected `technical_batch=None` behavior.
- Keep GloScope/PILOT/QOT on exact persisted Harmony keys; do not add pilot-gm-vae handling.

### 3. R configuration and corrected method plumbing

Targets:

- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R`
- `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_prepare_pseudobulk.R`
- `src/5_run_benchmark_methods/benchmark_hpc_utils.R`
- `src/5_run_benchmark_methods/benchmark_pipeline.R`
- `src/5_run_benchmark_methods/benchmark_methods_r.R`
- `src/utils/pseudobulk.R`
- `src/utils/py/h5ad_pseudobulk.py`

Changes:

- Normalize scalar/list `entry$batch_col` at the corrected execution boundary and load every required metadata column.
- Add strict R-side validation for names, values, sample constancy, levels, and sample order before any first-observation metadata collapse.
- Make the Python metadata/aggregation path validate all configured batch columns across the full selected cell metadata before `read_h5ad_sample_metadata()` or `aggregate_h5ad_counts_by_sample()` reduces rows to one record per `Sample`; never let those loaders' first-observation behavior decide a batch value.
- Extend `correct_clr_batch_lmm()` with dynamic safe aliases and one additive random-intercept term per configured key. Subtract the sum of fitted random effects, preserve dimensions/rownames, and recenter every row to exact zero sum.
- Keep biological labels out of the model formula and preserve the existing convergence/error policy.
- Build the R-side composite for corrected Pseudobulk at sample metadata level and pass it as the scalar `batch_col` through `get_pb_deseq2()`, `prepare_pseudobulks_hpc()`, and `DESeq2.normalize()`.
- In corrected Pseudobulk, force the exact batch-only DESeq2 design `~ 1`; use the composite only as the scalar removal batch. Do not carry the current contradictory `~ batch_col` behavior into the corrected contract.
- Apply additive CLR correction to all three corrected ECODA bundles, including the deterministic shuffled-label bundle; test the three-key Breast_cancer-shaped case.
- Leave GloScope's exact embedding resolution and all ordinary benchmark code unchanged.

### 4. Validator, submitter, and run metadata contracts

Targets:

- `src/5_run_benchmark_methods/1_submit_hpc_array.sh`
- `src/utils/bash/ecoda_run_common.sh`
- `src/utils/bash/ecoda_run_audit.sh`
- `src/utils/py/benchmark_h5ad_contract.py`
- `src/5_run_benchmark_methods/matrix_artifact_validator.py`
- `src/5_run_benchmark_methods/validate_benchmark_rds_contract.R`
- `NOTES.md`
- `notebooks/dataset_onboarding/README.md`

Changes:

- Add validator-only checks for scalar/list/null batch schema, duplicate/unknown key names, required obs columns, sample constancy, and corrected-only key-set identity.
- Ensure malformed multi-key configuration fails before any compute scheduler boundary and creates no scheduler IDs or worker manifests.
- Preserve the literal twelve-row `batch_effect_uncorrected` selection contract and all existing ordinary benchmark selection behavior.
- Extend corrected artifact metadata validation without changing old uncorrected artifact paths or method names.
- Ensure corrected H5AD/result/pseudobulk cache checks compare the ordered configured key identity, encoding version, scalarization policy, model policy, and fingerprint before reusing bytes; an audit file alone is insufficient.
- Ensure artifact checks distinguish source/config identity from content checksum and record the exact corrected method/key mapping.
- Update `NOTES.md` §3.0/§3.1/§3.2–3.3 and `notebooks/dataset_onboarding/README.md` to define list-valued batch semantics, the corrected-only eight-method matrix, per-method scalarization/native-list behavior, exact DESeq2 `~ 1` policy, sample-constancy validation, and the exclusion of pilot-gm-vae. Keep the historical uncorrected notebook and artifacts read-only.
- Do not add pilot-gm-vae to corrected validators or corrected defaults; preserve all existing ordinary benchmark method lists/references unchanged and never select pilotgm for corrected work.

### 5. Focused regression contracts

Targets:

- `tests/test_ecoda_run_common.sh`
- `tests/test_preprocessing_stage_submitter.sh`
- `tests/test_stage2_submitter.sh`
- `tests/test_artifact_contracts.py`
- `src/5_run_benchmark_methods/test_oom_retry.sh`
- `tests/test_multibatch_contracts.py`

Required observable coverage:

1. Scalar batch configuration remains accepted with unchanged uncorrected/benchmark selection semantics.
2. Readers preserve raw scalar/list/null values while execution-boundary normalization rejects malformed batch schema before compute.
3. Two- and three-key fixtures produce identical golden-vector composites in the separate HVG, MRVI, and R pseudobulk builders, including separator, escape, Unicode, sentinel, and fingerprint cases.
4. Composite encoding is collision-safe, deterministic, versioned, and rejects missing/sentinel/duplicate/unknown keys.
5. Within-Sample disagreement fails closed before first-observation metadata collapse.
6. HVG receives one temporary composite key and leaves no temporary obs column.
7. Harmony receives the original ordered list of keys and writes the existing exact corrected embedding key.
8. The MRVI loader requests each original batch column, never stringifies the list, and MRVI receives one composite key only.
9. ECODA additive correction preserves sample order, finite values, and exact CLR row zero sums; nonconvergence and rank deficiency fail closed.
10. Pseudobulk receives the scalar composite, uses corrected DESeq2 design `~ 1`, and applies batch-only removal with `correct_batch=TRUE`; the three-key case is covered.
11. GloScope/PILOT/QOT resolve only the exact corrected Harmony embedding; pilot-gm-vae is not selected.
12. A malformed multi-key selection/configuration triggers no `sbatch`, no scheduler manifest, and no run-owned compute state.
13. Existing uncorrected exact-selection, checksum, ownership, and validator-only NOOP tests continue to pass without launching jobs.
14. `NOTES.md` and `notebooks/dataset_onboarding/README.md` no longer assert scalar-only corrected batches or include pilot-gm-vae in the corrected method matrix; ordinary benchmark references remain unchanged.

## Dependency order

```text
Batch contract + synthetic fixtures
        |
        +--> preprocessing HVG/Harmony adapters
        +--> Python MRVI adapter
        +--> R ECODA/Pseudobulk adapters
        |
        +--> corrected artifact/run validators
        |
        +--> focused regression/no-launch checks
        |
        +--> scalar _debug corrected smoke (only after local contracts pass)
```

No step depends on or launches the old benchmark or batch-uncorrected pipelines.

## Verification after implementation approval

1. Read-only config/schema preflight over the current `datasets.json`; do not rewrite it.
2. Run focused synthetic Python/R contract checks with deterministic metadata and no H5AD production outputs.
3. Run shell syntax and no-launch submitter tests using scheduler stubs only.
4. Exercise only the `_debug` `batch_effect_corrected` scalar path in a fresh run-owned output root; verify finite required keys, unchanged source checksum/mtime, and no temporary composite column.
5. Validate existing benchmark and `batch_effect_uncorrected` contracts in validator-only mode; do not call submitters with production selections, `--force`, `prepare`, or `launch`.
6. Verify corrected cache reuse rejects missing/mismatched key-contract identity even when artifact bytes and checksums are valid.
7. If a future full-cohort corrected run is approved, select only explicitly named corrected dataset/view/method rows, record the expected row count and key contract, and route it through `durable-hpc-gate-ecoda`. That launch is outside this session and requires separate approval.

## Risks and mitigations

- **Separate builders can drift.** Mitigate with a shared written encoding contract, cross-builder fixture equality, and a recorded encoding version; do not silently accept mismatched fingerprints.
- **Composite level explosion.** Reject non-estimable or near-unique composite designs before correction and report the affected dataset/key set.
- **Confounding with biology.** Keep labels evaluation-only, reject direct label-key reuse, and report—not silently mask—confounding/estimation failures.
- **Temporary-column leakage.** Use cleanup guards and a pre-write assertion; run metadata, not `obs`, carries the key identity.
- **Old artifact invalidation.** Restrict code-path changes to corrected mode, preserve old uncorrected/benchmark validators and paths, and never rewrite existing files in this work.
