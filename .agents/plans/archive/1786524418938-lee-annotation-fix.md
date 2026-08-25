# Fix Lee stale-annotation + re-run pipeline (standardized scATOMIC/HiTME)

## Context (established facts)

- Lee's benchmark h5ad obs has scATOMIC columns (`layer_1..layer_6`, `scATOMIC_pred`, `S.Score`, `G2M.Score`, `Phase`, `classification_confidence`) and legacy columns (`scGate_multi`, `functional.cluster`) but **no** HiTME columns (`layer1/layer2/layer3`). `datasets.json` declares Lee `cell_type_high_res: "layer2"` → 4 failed tasks (trans 4311561_6, zeroimp 4311572_6, PILOT 4311499_6; pending scPoli 4311498_6 will fail too). MrVI passed (no ct col needed).
- Root cause: the original Lee rds carries pre-existing annotation columns. Preprocessing preserves all observation columns, and the union h5ad plus `get_seurat_obj_from_h5ad` (`seurat_utils.R:337`, `meta.data = meta_data`) carry them into the annotation worker's Seurat object, where:
  - `2.1.1_process_chunk.R:275` skips scATOMIC if `layer_1` already exists → scATOMIC silently skipped (stale values passed through).
  - `Run.HiTME` sees old scGate/ProjecTILs columns → produced nothing for Lee (no layer1/2/3).
- The current annotation run (4307176) shows NO scATOMIC/HiTME activity for Lee's tasks (323–336) → chunks were skipped as "already processed" (`2.1.1_process_chunk.R:211`) — Lee's feathers come from an older run.
- Smillie `prepare_pseudobulk` failure (4311518_8) was a separate transient (staged R library copy); re-run `--ds_name Smillie --methods prepare_pseudobulk` succeeded and its tail synced the whole benchmark dir.
- Kfoury's trans/zeroimp + pseudobulks were computed in an earlier session (skip-if-exists today) — user confirmed re-running trans/zeroimp --force.

## Sync-status answer (no action needed)

The Smillie re-run's tail synced `benchmark/` wholesale via `rsync -rlptDv` (`benchmark_submit_common.sh:265`) — **no `--delete`**, so nothing was removed from NAS. Exec-log merge (`1.1.2_merge_execution_times.py:124-127`) dedups on (dataset, method) keep="last" against the NAS `--existing-log` — the 93-row feather preserved all prior rows. The NAS snapshot is simply incomplete for Lee (pilot/scpoli/trans/zeroimp) and partial for scPoli (array still running); re-runs below complete it. `checksums.md5` is regenerated on every sync tail.

## Decisions (user-confirmed)

1. **Drop legacy annotation columns from final h5ad obs** during merge (standardized columns only).
2. **Re-annotate any dataset the audit flags** (not just Lee).
3. **Re-run Kfoury trans/zeroimp with --force**.

---

## Task 1 — Worker fix: drop known annotation columns before annotating

File: `src/4_cell_type_annotation/2.1.1_process_chunk.R`

- After `annot_cols` is built (lines 98–102), add a `LEGACY_ANNOT_COLS` constant: `c("scGate_multi", "functional.cluster")` plus any other known old-annotation names found in the audit (Task 3). Keep the whitelist semantics of `annot_cols` unchanged.
- In the per-sample loop, immediately after `seurat_obj <- get_seurat_obj_from_h5ad(...)` (line 250–254, before `NormalizeData`), drop every column in `union(annot_cols, LEGACY_ANNOT_COLS)` that exists in `seurat_obj@meta.data` (message with dropped names). This makes the line-275 scATOMIC guard always evaluate fresh (never sees stale `layer_1`) and hides old scGate/ProjecTILs columns from `Run.HiTME`.
- Guardrails: keep the existing SCRIPT_DIR recovery block, `PIXI_RSCRIPT`/`PYTHON_BIN` conventions untouched; no changes to `NormalizeData` or wall-clock budget logic.

## Task 2 — Merge fix: strip legacy annotation columns from final obs

File: `src/4_cell_type_annotation/3.1_merge_annotations.py`

- Add the same `LEGACY_ANNOT_COLS` list (single source: keep in sync with Task 1; comment pointing at both files).
- In `merge_annotations()`:
  - Extend the idempotency drop (lines 70–77) to include `LEGACY_ANNOT_COLS` (drop from `obs` before the join, like the whitelisted columns).
  - Extend the `orig_cols` computation (line 85) to exclude `LEGACY_ANNOT_COLS` so they never survive into the output obs.
- Result: final h5ad obs = original study metadata (Sample, condition, MGMT, IDH_status, …) + fresh standardized scATOMIC/HiTME columns only.

## Task 3 — Audit all benchmark datasets (identify affected)

Create a throwaway audit script (debug-cpu sbatch job, mirroring the earlier `~/inspect_lee.sh` pattern; `source ~/ECODA_paper/src/slurm_config.sh`, run `${PYTHON_BIN}` heredoc):

- For every dataset in the benchmark manifest (jq: `use_for_benchmark == true` and `views.benchmark_analysis != null`, skip `_*`): open `${HPC_SCRATCH_DIR}/<DS>/output/<output_file_name>` (backed) and print:
  - presence of `layer1/layer2/layer3` (HiTME) and `layer_1..layer_6`, `scATOMIC_pred`, `classification_confidence` (scATOMIC);
  - presence of legacy cols `scGate_multi`, `functional.cluster`, `*_UCell`, and any other column matching the whitelist names;
  - declared ct cols from datasets.json (`cell_type_low_res`/`cell_type_high_res`) vs actual obs columns.
- Also on the login node: grep each dataset's annotation logs (latest run, tasks from `chunks_manifest.txt` mapping) for `"scATOMIC attempt"` / `"HiTME attempt"` / `"Chunk already processed"` to flag datasets where annotation was skipped.
- Output: a table `<DS> | HiTME layer1/2/3? | scATOMIC layer_1..6? | legacy cols? | annotation actually ran?`. Any dataset with missing HiTME columns, missing scATOMIC activity, or legacy column collisions is **affected** and gets the Task 4 treatment (same as Lee).

## Task 4 — Clean + re-annotate affected datasets (Lee first)

Per affected dataset (Lee first), from `${PROJECT_ROOT}` on the login node:

```bash
./src/4_cell_type_annotation/1_prepare_chunks.sh production <DS> --force   # rebuilds union + chunks, deletes stale feathers
./src/4_cell_type_annotation/2_submit_hpc_array.sh <DS>                    # re-annotate with fixed worker
./src/4_cell_type_annotation/3_submit_merge.sh <DS>                        # merge feathers into view h5ad (+ NAS sync)
```

- The worker-side drop (Task 1) means no manual h5ad scrubbing is needed — stale columns are dropped from the Seurat metadata at sample start and legacy columns are stripped by the merge (Task 2).
- After the merge, verify via the audit-style check that the h5ad obs now has `layer1/layer2/layer3`, fresh `layer_1..6`, and no `scGate_multi`/`functional.cluster`.

## Task 5 — Re-run failed/pending benchmark tasks

After Lee (and any audited datasets) is re-annotated and merged:

```bash
# trans + zeroimp for Lee (and any other affected dataset)
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name Lee --analysis trans,zeroimp
# Kfoury consistency re-run (user-confirmed)
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name Kfoury --analysis trans,zeroimp --force
# PILOT for Lee
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods pilot
# scPoli for Lee — AFTER the currently pending array 4311498 finishes (its Lee task will fail; then re-run just Lee)
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods scpoli
# R benchmark methods (gloscope/mofa/pseudobulk/scitd) have NOT run yet for any dataset — full run (auto-prepends prepare_pseudobulk)
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh
```

- MrVI results for Lee stay valid (no ct column needed) — do not re-run.
- Each submitter's tail re-merges exec logs and re-syncs the whole `benchmark/` dir to NAS, completing the incomplete snapshot.
- Let the still-running arrays finish first (`4311497_11` MrVI/Zhang; `4311498` scPoli) — don't cancel them.

## Task 6 — Validation

- Lee h5ad obs: `layer1/layer2/layer3` present, legacy cols gone, scATOMIC cols fresh.
- `results/Lee_trans.rds`, `results/Lee_zeroimp.rds`, Lee PILOT/scPoli feathers exist in scratch and on NAS after syncs.
- Kfoury trans/zeroimp recomputed (--force) and re-synced.
- Audit table (Task 3) shows no remaining affected dataset.
- NAS completeness: spot-check `${NAS_TARGET_DIR}/benchmark/results/`, `embeddings/`, `pseudobulks/` for Lee + all datasets; `checksums.md5` regenerated.
- Exec-log continuity: `execution_times.feather` rows for Lee trans/zeroimp/pilot/scpoli present after final syncs.

## Contingency

- If HiTME still fails for Lee after the fix (e.g., brain tumor unsupported by HiTME/scGate models — will be visible in the new annotation `.err` as `HiTME error: ... - retrying...`), the fallback is a datasets.json change for Lee's `cell_type_low_res`/`cell_type_high_res` (e.g. scATOMIC `layer_1`/`layer_2` or the original `functional.cluster`). **datasets.json changes require explicit user approval** — stop and ask before touching it.

## Delivery

- Follow AGENTS.md Task Completion Workflow after implementation: move this plan to `.kilo/plans/archive/`, `git add .`, commit, push.
- Do NOT run pipeline scripts locally for validation (HPC-only); the audit/verify jobs run via sbatch debug-cpu per the earlier pattern.
