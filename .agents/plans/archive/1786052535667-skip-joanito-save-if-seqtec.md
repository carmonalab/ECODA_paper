# Skip in-place Joanito .rds save when `seqtec` already present

## Goal
In `src/2_dataset_specific_preprocessing/1.2.1_prepare_joanito.R`, avoid re-writing the large full Joanito `.rds` (`saveRDS(seurat, input)`) when the `seqtec` metadata column already exists on disk (i.e., a previous run already computed it). Saves time on re-runs; `_debug` subset build is unaffected (uses the in-memory object, where `seqtec` is already present).

## Changes

### 1. `src/2_dataset_specific_preprocessing/1.2.1_prepare_joanito.R`
Wrap the seqtec computation + save (lines 50–62) in a presence check:

```r
# -----------------------------------------------------------------------------
# 1. seqtec batch column (single source of truth — see header note)
# -----------------------------------------------------------------------------
if (!"seqtec" %in% colnames(seurat@meta.data)) {
  seurat$seqtec <- ifelse(
    seurat$dataset %in% c("CRC-SG1", "KUL5"),
    "5' seq",
    "3' seq"
  )
  # In-place save is idempotent. Must run AFTER 1_stage_data.sh and BEFORE
  # the preprocess array.
  saveRDS(seurat, input)
  message("Saved seqtec back to: ", input)
} else {
  message("seqtec already present — skipping recompute and in-place save.")
}
```

- Decision: when the column is present, skip **both** recompute and save (values were written by a previous run of this same script, so they are already per the mapping). Recompute is cheap, but the guard on presence is simpler and matches the intent.
- No changes to the `_debug` subset section (lines 64–153): it reads `seqtec` from the in-memory object, which is correct in both branches.

### 2. Header comment (lines 7–10)
Update step-1 description: "computes the `seqtec` batch column and saves the full object back in place (idempotent — recomputes seqtec each run)" → state that compute+save run only when `seqtec` is absent; otherwise the in-place save is skipped.

### 3. `docs/ARCHITECTURE.md` (lines 100 and 102)
Update the **Joanito prepare** bullet to note the fast path: the in-place `saveRDS` is skipped when `seqtec` already exists in the staged `.rds`. Line 102 ("after the `seqtec` computation and `saveRDS`") can stay or be reworded to "after the `seqtec` computation (and in-place save, on first run)" — keep consistent with line 100.

## Caveats (no action needed)
- If the `seqtec` mapping in this file ever changes, existing on-disk `.rds` files would keep stale values on re-run (the guard skips the save). Recovery = re-stage / delete the staged `.rds` and rerun. Accepted per user's request to prioritize speed on re-runs.
- Script must still run after `1_stage_data.sh` and before the preprocess array — unchanged.

## Validation
- Per repo rules: do not run pipeline scripts; validation happens when the full pipeline is tested with the `_debug` dataset.
- Static check only: on re-run with an existing staged `.rds`, the message "seqtec already present — skipping recompute and in-place save." appears and no `saveRDS` occurs; on a fresh staged `.rds` (no `seqtec` column), behavior is unchanged (compute + save + subset build).
