# DESeq2 batch-aware normalization + vst sparsity fix

## Context

Two coupled problems in `src/utils/pseudobulk.R`:

1. **Kfoury pseudobulk job failure (4295224):** `DESeq2::vst()` hard-stops with
   "less than 'nsub' rows with mean normalized count > 5" when a per-cell-type
   pseudobulk subset (Kfoury prostate-metastasis CT subsets, computed in
   `process_pseudobulk_ct_fig` → `get_pb_deseq2` → `DESeq2.normalize`) has
   < 1000 genes (default `nsub`) with mean normalized count > 5. The nsub check
   is **unconditional** in DESeq2's `vst()` (verified in package source); the
   first failing CT kills the whole task → array task FAILED → no NAS sync.
   The same failure will hit the local notebook CT combos on Kfoury.

2. **Batch correction exists but is dead code:** `DESeq2.normalize()` has a
   `batch_col` parameter (added 2026-07-15, commit `5b7aaca`) that builds
   `design = ~ batch_col` and calls `limma::removeBatchEffect(x, batch)`, but
   **no caller passes `batch_col`** — not `get_pb_deseq2()`, not
   `benchmark_pipeline.R`, not `notebooks/batch_effect_analysis.rmd`
   (which defines `batch_col` but never passes it). The batch-effect
   pseudobulk is therefore effectively uncorrected.

## Research conclusions (authoritative, with sources)

- VST/rlog **do not remove** batch variation: *"It uses the design formula to
  calculate the within-group variability (if `blind=FALSE`) or the
  across-all-samples variability (if `blind=TRUE`). It does **not** use the
  design to remove variation in the data."* (DESeq2 vignette, FAQ "Why after
  VST are there still batches in the PCA plot?"). The only removal step is
  `limma::removeBatchEffect()`.
- The design **does** affect transformed values via the dispersion trend:
  `blind=TRUE` treats group differences as noise → over-shrinkage; `blind=FALSE`
  uses within-group variability. Michael Love (support.bioconductor.org/p/9148734):
  *"I recommend blind=FALSE... It only looks at the design to know the global
  distribution of dispersion values."*
- `removeBatchEffect(mat, batch, design=model.matrix(~condition,...))` protects
  covariates via the `design` argument (vignette FAQ). **We deliberately do NOT
  use that here — see no-leakage decision below.**
- `design = ~ Sample` at pseudobulk level (one column per sample) is a saturated
  model → degenerate dispersion estimation → rejected for the benchmark
  (DESeq2 FAQ "dataset without replicates": No). The `preprocess.py`
  batch-aware-HVG analogy does not transfer to pseudobulk.
- No-leakage is the central premise of the project: bio labels are ground truth
  **only**, never inputs to any preprocessing/embedding step. Batch correction
  must therefore be batch-only (no bio protection). Known limitation (bio signal
  correlated with batch can be partially removed) is accepted and documented.

## Design decisions

### `DESeq2.normalize()` — new signature (src/utils/pseudobulk.R)

```r
DESeq2.normalize <- function(
  matrix, metadata, n_hvg = 2000,
  batch_col = NULL,      # batch column (batch-effect analysis only)
  blind = TRUE,          # explicit; benchmark = TRUE (legacy-equivalent)
  correct_batch = FALSE  # apply limma::removeBatchEffect(x, batch) — batch-only, NO design protection
)
```

- **Defaults preserve legacy behavior for every existing caller**
  (HPC benchmark, notebook benchmark): `batch_col=NULL, blind=TRUE,
  correct_batch=FALSE` → `design ~ 1`, `vst(blind=TRUE)` — identical to
  legacy imp 1 (`vst(matrix)` default) and equivalent to today's
  `vst(dds, blind=FALSE)` with `design ~ 1`.
- When `batch_col` given: `design = ~ batch_col`; `blind=FALSE`;
  `vst` uses the batch design for dispersion estimation (no removal yet).
- When `correct_batch=TRUE` (requires `batch_col`): after vst,
  `norm_matrix <- limma::removeBatchEffect(norm_matrix, batch = metadata[[batch_col]])`
  **without** a `design` argument — batch-only, no bio protection, no bio
  parameter anywhere in the function (enforces no-leakage by construction).
  If `batch_col` is NULL but `correct_batch=TRUE`: `warning()` + skip.
- **Sparsity fallback chain** (fixes Kfoury; applies regardless of blind/batch):
  1. `n_gt5 <- sum(MatrixGenerics::rowMeans(DESeq2::counts(dds, normalized = TRUE)) > 5)`
  2. `n_gt5 == 0` → `norm_matrix <- log2(counts(dds, normalized=TRUE) + 1)` + `message()`.
  3. else `tryCatch`:
     `DESeq2::vst(dds, blind = blind, nsub = min(1000, n_gt5))`
     → error → `DESeq2::varianceStabilizingTransformation(dds, blind = blind, fitType = "mean")`
     → error → `DESeq2::varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean")`
     → error → `log2(counts(dds, normalized=TRUE) + 1)` + `warning()`.
  - Normal path (n_gt5 ≥ 1000) is byte-identical to today.
  - `fitType = "mean"` avoids nls fragility on few genes; `blind=TRUE` fallback
    dodges rank-deficient batch designs on tiny subsets.

### `get_pb_deseq2()` — pass-through (src/utils/pseudobulk.R)

```r
get_pb_deseq2 <- function(seurat, sample_col = "Sample", hvg = NULL, n_hvg = 2000,
                          black_list = "none", batch_col = NULL, blind = TRUE, correct_batch = FALSE)
```

Passes the three new args to `DESeq2.normalize()`. All existing callers
unchanged (defaults).

### Label alignment — verified, NO changes

`create_result_bundle()` (benchmark_methods_r.R:422) re-aligns
`labels <- labels[rownames(feat_mat)]`; `get_labels()` returns labels named by
Sample. Imp 1's explicit line was moved there in the `3805890` refactor; the
"missing re-alignment" claim is not a bug. Do not re-add redundant code; a
comment in `process_pseudobulk_ct_fig` noting that alignment is centralized is
optional.

## Implementation tasks (ordered)

1. **AGENTS.md — no-leakage guardrail (highlight, central premise):**
   - Add to "Agent Guardrails & Domain Terms": bio labels (Status, sample.origin,
     cond, …) are **ground truth only**; they must NEVER be passed as a design
     covariate, batch key, or any other input to preprocessing, DESeq2
     normalization, batch correction, HVG selection, or embedding steps. Batch
     correction is batch-only (no `design` argument in `removeBatchEffect`).
     Violation = supervised analysis, invalidates the paper's premise.
   - Add `DESeq2.normalize()` semantics: `blind`/`batch_col`/`correct_batch`
     defaults; benchmark = `blind=TRUE` no correction; batch-effect analysis =
     `batch_col` + `blind=FALSE` + `correct_batch=TRUE`.

2. **src/utils/pseudobulk.R — rewrite `DESeq2.normalize()`** per Design decisions
   (new params, fallback chain, batch-only limma step). Keep the existing
   `black_list` no-op typo at line 84 intact (legacy-equivalent, documented).

3. **src/utils/pseudobulk.R — extend `get_pb_deseq2()`** with the three
   pass-through params.

4. **src/5_run_benchmark_methods/benchmark_methods_r.R — harden
   `process_pseudobulk_ct_fig()`:** wrap the per-CT `get_pb_deseq2()` call in
   `tryCatch`; on error, `warning()` naming the CT and `next` (skip it). One
   sparse CT must never kill the whole method again. (Belt-and-braces on top of
   the sparsity fix; the dist averaging already handles missing CTs via
   `count_mat`.)

5. **notebooks/batch_effect_analysis.rmd — wire batch correction (batch-only):**
   - Joanito (line 127): `get_pb_deseq2(seurat, sample_col="sample.ID",
     black_list="none", batch_col="seqtec", blind=FALSE, correct_batch=TRUE)`.
   - Stephenson (lines 218, 234): `DESeq2.normalize(pb, metadata,
     n_hvg=2000, batch_col="Site", blind=FALSE, correct_batch=TRUE)`.
   - Stephenson blacklist (line 374): `get_pb_deseq2(seurat_sub,
     sample_col="Sample", black_list="default", batch_col="Site",
     blind=FALSE, correct_batch=TRUE)`.
   - Combined PBMC (line 592): `DESeq2.normalize(combined_pseudobulk,
     combined_metadata, n_hvg=2000, batch_col="ds", blind=FALSE,
     correct_batch=TRUE)` — batch is the study column `ds`, NOT `Site`
     (Site is confounded with study).
   - Add short comments: bio labels used in `plot_mds` are evaluation-only;
     nothing bio-related is passed to the normalization (no leakage).
   - NOTE: `vst(blind=FALSE)` with `~batch` changes the transformed values for
     these analyses vs. the old `blind=TRUE` output — expected and intended.

6. **TODO.md — Phase 4 updates:** mark "verify `DESeq2.normalize()` `batch_col`
   is correctly implemented and does not get `'Sample'` as batch column" as
   done (it now is: batch-only, wired in the batch-effect notebook, and
   `"Sample"` is never used as batch column); note "Pseudobulk DESeq2+limma
   with `batch_col`" implemented for the notebook analyses.

7. **docs/ARCHITECTURE.md — minimal doc touch:** update the `pseudobulk.R`
   file-role entry (if listed) with the new params and the two usage modes.

8. **Validation (repo convention: no pipeline runs unless asked):**
   - Syntax: `Rscript -e 'invisible(lapply(c("src/utils/pseudobulk.R",
     "src/5_run_benchmark_methods/benchmark_methods_r.R"), function(f) parse(f)))'`
   - Grep all `DESeq2.normalize`/`get_pb_deseq2` callers → confirm defaults keep
     them behavior-identical (batch-free paths unchanged).
   - Notebook rendering + HPC runs are user-side; full pipeline validation
     happens later with the `_debug` dataset per AGENTS.md.

## Rollout — Kfoury resume

After implementation:
```
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Kfoury --methods pseudobulk
```
prepare_pseudobulk cache from job 4295221 is valid; per-combo RDS caches resume;
only the previously-failed CT combos recompute; sacct gate syncs to NAS on
success. Same command pattern for any other dataset/notebook re-runs.

## Risks / edge cases

- **Batch level with one sample** (rank-deficient design, `blind=FALSE`) →
  fallback chain degrades gracefully (fitType="mean" → blind=TRUE → log2) with
  warnings; results remain usable.
- **`correct_batch=TRUE` without `batch_col`** → warning + skip (defensive).
- **Combined PBMC**: `ds` column must exist in `combined_metadata` (it is
  created at lines 573-578); keep `batch="ds"`.
- Stephenson `Site` levels: verify ≥2 samples per site for a meaningful fit;
  fallbacks cover degenerate cases.
- No result-file/format changes for any benchmark bundle; only values for
  sparse-CT combos change (from "missing" to computed).

## Out of scope

- Full Phase 4 expansion (ECODA batch-associated CT removal, MrVI `batch_key`,
  GloScope/PILOT-GM-VAE on `X_pca_harmony`, HPC-izing the batch-effect
  notebook) — stays in TODO.md Phase 4.
- `preprocess.py` HVG batch-key logic (already correct: benchmark views use
  Sample, batch-effect views use batch_col).
- SVA/ComBat alternatives; bio-protected `removeBatchEffect(design=...)` mode.
