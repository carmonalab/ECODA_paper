# Clean-annotation wipe + pipeline-entry check (supersedes lee-annotation-fix)

Supersedes `.kilo/plans/1786524418938-lee-annotation-fix.md` (same problem, stronger design; Tasks 1–2 there are replaced by Tasks 1–2 here, HPC steps carry over). Archive both plans at delivery.

## Context (established facts)

- Preprocessing keeps ALL source-rds obs columns; the annotation union and the worker's Seurat object (`get_seurat_obj_from_h5ad` → `CreateSeuratObject(meta.data = meta_data)`, `src/utils/seurat_utils.R:335-338`) carry them all. Legacy annotation columns therefore leak into annotation and into the final view h5ad obs (Lee: `scGate_multi`, `functional.cluster`; Wu: `scGate_multi` + 22 `*_UCell` scores — verified by local rds/NAS inspection; all other datasets clean).
- The old worker kept `hitme_cols_keep`/`scatomic_cols` as a keep-whitelist after annotation — legacy remnants of a "minimal metadata" era; the whitelist must NOT apply to pre-existing columns.
- Old pipeline merges only ever wrote the `annot_cols` names into view h5ads; raw rds annotation columns are the only unknown source. Feathers stay in `${HPC_SCRATCH_DIR}/<DS>/output/` after merge (kept deliberately for the coverage gate, `3_submit_merge.sh:241`) — they are the record of exactly what THIS pipeline produced.
- Current worker/merge (commit `f60fa57`) drop a fixed name list + patterns; insufficient per user: no whitelist at all, wipe ALL pre-existing annotation columns, plus an entry check.

## Decisions (user-confirmed)

1. **Worker wipes ALL meta.data except the sample column** before annotation — no pattern lists needed; scATOMIC `layer_1` guard always fires fresh; `Run.HiTME` sees zero pre-existing columns.
2. **Merge drops pre-existing annotation columns by deny-patterns (Tier 1, unconditional)** and keeps ALL other metadata (QC cols, unlisted study cols — "minimal metadata is not an issue anymore").
3. **Tier-2 columns are dropped with a loud per-dataset warning** (exact-name match only, so author columns like a plain `cellCycle` or `Phase_variant` never collide).
4. **Entry check**: a dataset entering annotation is "clean" iff every pattern-matching obs column is also a column of its own feathers (i.e. produced by this pipeline); otherwise it is NOT clean → never skipped → chunks rebuilt → fresh re-annotation.
5. **Marker**: merge writes `adata.uns["ecoda_annotation_version"] = "1"` — "annotated here" provenance for validation/future logic.
6. **Re-annotation scope**: Lee + Wu (audit confirms no other dataset affected; any dataset the check flags gets the same treatment).
7. `hitme_cols_keep`/`scatomic_cols`/`annot_cols` are renamed to output-set constants (fresh-output extraction only, no whitelist semantics).

## Column classification (single source per file; mirror comments between R and Python — repo convention, no shared constant file)

- `LEGACY_ANNOT_TIER1` (patterns, unconditional drop): `^scGate`, `^functional\.cluster`, `_UCell$`, `^scATOMIC`, `^layer_?\d` (covers `layer1/2/3` + `layer_1..6`; a plain `layer` column survives).
- `LEGACY_ANNOT_TIER2` (exact names, drop + warning): `S.Score`, `G2M.Score`, `Phase`, `classification_confidence`, `cellCycle.G1S_UCell`, `cellCycle.G2M_UCell` (the two `cellCycle*` names are also matched by the `_UCell` suffix; exact-name listing is for the warning).
- `ANNOT_OUTPUT_COLS` (= old `annot_cols`, renamed): fresh-output extraction set only — the ONLY annotation columns allowed in final obs, and only from feathers.
- None of the Tier-2 exact names exist in the current corpus's raw files → Tier 2 currently drops zero author metadata; the warning documents it.

## Task 1 — Worker: wipe obs to Sample-only after the colname check (minimal Seurat construction) + safe perf items

Files: `src/4_cell_type_annotation/2.1.1_process_chunk.R` + `src/utils/seurat_utils.R`

- `2.1.1_process_chunk.R`: immediately after the `args$sample_colname %in% colnames(obs)` check (lines 202–204), collapse the once-per-chunk obs to the sample column only: `obs <- obs[, args$sample_colname, drop = FALSE]`, with a `message()` naming the dropped columns (count + names) — "wiping N pre-existing obs columns; building minimal Seurat objects (meta.data = Sample only)".
- `src/utils/seurat_utils.R`: `get_seurat_obj_from_h5ad` currently IGNORES the R-side `r_obs` for metadata and re-loads ALL columns from Python per sample (`meta_data <- as.data.frame(py_to_r(subset_py$obs))`), which would bypass the worker's wipe. Change that ONE line to build meta_data from the `r_obs` it already receives: `meta_data <- r_obs[r_indices, , drop = FALSE]` (keep `as.data.frame()` wrapper if needed for non-data.frame inputs).
  - NO signature change, NO new argument. `r_obs` stays the metadata source; `cell_names <- rownames(meta_data)` replaces the per-sample `subset_py$obs_names$values` read.
  - Behavior-identical for the benchmark caller (`load_benchmark_seurat`, `benchmark_hpc_utils.R:121`, passes the full obs → full meta; same rows/columns/rownames/factor levels — reticulate pandas-categorical factors keep all levels on subset, matching R factor subsetting).
  - Eliminates the per-sample full-obs `py_to_r` copies.
  - Add a comment: meta_data intentionally sourced from the R-side obs so annotation workers can pre-wipe metadata; callers needing full metadata pass the full obs.
- Perf item (a) — Python-side transpose in `extract_matrix` (user-approved, MUST validate): keep `py_mat$astype("float64")`, then `.T$tocsc()` instead of `.tocsc()` + R `t()`: `r_mat <- py_to_r(py_mat$astype("float64")$T$tocsc())`; set `dimnames(r_mat) <- list(gene_names, cell_names)` (drop the R-side `t()`). scipy `csr.T` returns a CSC view (O(1)); reticulate maps scipy CSC → R `dgCMatrix` directly. KEEP the `astype("float64")` cast (view `layers["counts"]` can be integer-dtype; scATOMIC needs doubles). Empty-counts fallback (`Matrix::sparseMatrix`) unchanged.
- Perf item (b) — hoist in the worker (user-approved): compute `layer_keys` (via the existing hardened `import_builtins(convert = FALSE)$list(adata$layers$keys())` — do NOT switch to `names()`; keys()-view crash history, commit `8aa3ba4`) and `counts_layer` ONCE per chunk, outside the per-sample loop (currently recomputed per sample ~lines 239–240). Pass `counts_layer` into the per-sample `get_seurat_obj_from_h5ad` call as today (no signature change).
- Perf item (c) — `cell_names <- rownames(meta_data)` (covered above); no `gene_names` param (rejected: signature change for negligible gain).
- REMOVE the drop block added in `f60fa57` (post-hoc `seurat_obj@meta.data` wipe after `get_seurat_obj_from_h5ad(...)`, before `NormalizeData`) — superseded: the object is minimal by construction.
- Edge: the line-202 stop fires first if the sample column is missing; `NormalizeData`, wall-clock budget, SCRIPT_DIR recovery block, `PIXI_RSCRIPT`/`PYTHON_BIN` untouched.

### Task 1 validation (user-requested — local, not HPC)

Throwaway unit test (pixi `default` env, which has reticulate + Seurat + anndata 0.12.10; run from `/tmp`, not committed):
1. Build a tiny h5ad (e.g. 20 cells × 10 genes, integer counts in `X` + `layers["counts"]`, obs with legacy columns `scGate_multi`/`functional.cluster`/`Phase` + `Sample`, obsm `X_pca_test`).
2. In R via reticulate: read backed, build `obs` df, call `get_seurat_obj_from_h5ad` TWICE per scenario — once with a `git stash`-style old-path copy (`.tocsc()` + `t()`, `py_to_r(subset_py$obs)`) and once with the new path — and assert byte-identical results:
   - scenario A (benchmark caller): full obs → meta.data has ALL columns, identical names/values/rownames; counts matrix values + dimnames identical; embeddings attached identically.
   - scenario B (worker): Sample-only obs → meta.data has ONLY `Sample`; counts identical.
3. Assert matrix equality exactly (`identical()` on `as.matrix` / `all.equal` on the sparse structure), and `storage.mode` is double.
4. Report old-vs-new equivalence in the delivery commit message.

If the local test cannot run (reticulate/python path mismatch), fall back to an HPC debug-cpu smoke job (same assertions, sbatch + heredoc pattern).

## Task 2 — Merge: tiered drop + warning + marker + invariant

File: `src/4_cell_type_annotation/3.1_merge_annotations.py`

- Rename `hitme_cols_keep`/`scatomic_cols`/`annot_cols` → `ANNOT_OUTPUT_COLS` (keep sub-lists). Comment: NOT a keep-whitelist — the only columns this pipeline emits; pre-existing versions are always dropped.
- `drop_cols = ANNOT_OUTPUT_COLS ∪ {Tier-1 pattern matches} ∪ {Tier-2 exact names present}` (patterns evaluated on `obs.columns`; the `legacy_ucell`/`legacy_prefix` blocks from `f60fa57` collapse into this).
- Idempotency drop (existing `existing = [c for c in drop_cols if c in obs.columns]`) unchanged in mechanics; extend its print to a clear warning: total count + names, Tier-2 names explicitly flagged as "possibly author metadata (exact-name match)".
- `orig_cols` computation unchanged (excludes `drop_cols`).
- Post-merge invariant (fail loudly): every final obs column matching Tier-1 patterns or Tier-2 exact names must be present in the feather columns (they came fresh from this run). If violated → `sys.exit(1)` with the offending names.
- Marker: after the join, set `adata.uns["ecoda_annotation_version"] = "1"` before `write_h5ad` (both the in-place and `--output-path` paths).
- Keep the no-unannotated-h5ad guard, match-rate guard, atomic write.

## Task 3 — Entry check ("clean or annotated here")

- `src/4_cell_type_annotation/1.1_prepare_chunks.py`: add `--check-clean` mode (no chunk/union writes): for each view h5ad in `output/*.h5ad` excluding `*_raw.h5ad` (h5py obs column-name read, mirroring `read_obs_keys_h5py` fail-closed style), collect columns matching Tier-1/Tier-2; read the feather column set from `output/annotations_chunk_*.feather` (schema-only read via pyarrow); `legacy = matches - feather_cols`. Print per view: matched names, counts, Tier-1 vs Tier-2 split, and whether it looks like previous annotation (Tier-1 names) or possibly author metadata (Tier-2 only). Exit codes: 0 = clean/no views, 2 = legacy found, 1 = error.
- `src/4_cell_type_annotation/1_prepare_chunks.sh`: in the per-dataset loop, when the ANNOTATED branch would skip (ANNOTATED==1 && FORCE_ARG==0), run the check first via the existing short-srun pattern (30 min/1 cpu/4G, `${PYTHON_BIN}`): if exit 2 → print the warning, set `ANNOTATED=0`, fall through to rebuild (chunks + feather deletion in production mode → re-annotation is then required and happens normally). If exit 0 → skip as today. Never use `--force` semantics to circumvent.
- Add the check output to the run summary (`SKIPPED_ANNOTATED` vs new `FLAGGED_LEGACY` list).

## Task 4 — Constant cleanup

- Both files: rename to `ANNOT_OUTPUT_COLS` (+ `HITME_OUTPUT_COLS`, `SCATOMIC_OUTPUT_COLS` sub-lists), `LEGACY_ANNOT_TIER1`/`LEGACY_ANNOT_TIER2`; cross-file mirror comments ("must mirror 2.1.1_process_chunk.R" / "…3.1_merge_annotations.py").
- Update the audit script `audit_annotation_state.sh` (transfer artifact) to the same tier lists.

## Task 5 — Audit (agent, NAS-side — no HPC needed)

- Re-run the NAS obs scan (`/Volumes/Shared/Projects/ECODA_paper/<DS>/output/*.h5ad`, ALL views incl. batch-effect views, excluding `*_raw.h5ad`) with the Tier-1/Tier-2 lists: table `<DS> | view | legacy cols (names+counts) | ct cols ok`. Expect: only Lee and Wu flagged. Report → confirms Task 6 scope.

## Task 6 — HPC execution (USER runs; agent has no SSH access)

From `${HOME}/ECODA_paper` on the bamboo login node, inside `tmux` (submit tails are login-bound and blocking):

```bash
git pull                                  # brings Tasks 1-4
# Task 6a — Lee + Wu: check flags them automatically → rebuild chunks, re-annotate, merge
./src/4_cell_type_annotation/1_prepare_chunks.sh production Lee
./src/4_cell_type_annotation/2_submit_hpc_array.sh Lee
./src/4_cell_type_annotation/3_submit_merge.sh Lee
./src/4_cell_type_annotation/1_prepare_chunks.sh production Wu
./src/4_cell_type_annotation/2_submit_hpc_array.sh Wu
./src/4_cell_type_annotation/3_submit_merge.sh Wu
# Task 6b — benchmark re-runs (wait for arrays 4311497/4311498 to finish first)
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name Lee --analysis trans,zeroimp
./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh --ds_name Kfoury --analysis trans,zeroimp --force
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods pilot
./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh --ds_name Lee --methods scpoli   # AFTER 4311498 finishes
./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh   # full R run (auto-prepends prepare_pseudobulk)
```

- MrVI not re-run (no ct col needed). Other datasets untouched (clean, no re-annotation — per user decision).
- No manual h5ad scrubbing needed: worker wipe + merge tiered drop cover both paths.

## Task 7 — Validation

- Agent (NAS-side, after syncs): Lee/Wu view h5ads have `layer1/2/3` + `layer_1..6` fresh, ZERO Tier-1/Tier-2 columns, `uns["ecoda_annotation_version"] == "1"`, declared ct cols present; other datasets unchanged.
- User (HPC logs): Lee/Wu annotation `.err` show `scATOMIC attempt`/`HiTME attempt` lines (annotation actually ran — no "Chunk already processed" skip), merge logs show the Tier drop warnings (or clean), `results/Lee_trans.rds`, `results/Lee_zeroimp.rds`, Lee PILOT/scPoli feathers exist; Kfoury trans/zeroimp recomputed; exec-log rows present after final syncs; NAS `benchmark/results/`, `embeddings/`, `pseudobulks/` complete; `checksums.md5` regenerated.
- Re-run the Task 5 audit table → no remaining affected dataset.

## Delivery

- AGENTS.md Task Completion Workflow: archive THIS plan + `1786524418938-lee-annotation-fix.md` into `.kilo/plans/archive/`, `git add .`, commit, push.

## Risks / contingency

- **HiTME may still fail for Lee** (brain tumor unsupported by scGate/HiTME models — visible in new `.err` as `HiTME error: … - retrying…` × 5). Fallback = datasets.json change for Lee's `cell_type_low_res`/`cell_type_high_res` — **requires explicit user approval; stop and ask** before touching datasets.json.
- **Value-level staleness** (stale values under pipeline column names, e.g. old pass-through `layer_1`): undetectable by column-name check. Mitigated by the worker wipe (no future run can produce it) + the corpus audit; entry-check subset test catches name-level legacy. Accepted limitation.
- Entry-check feather-subset test relies on feathers surviving in `output/` (true post-merge by design); a dataset whose feathers were deleted but is otherwise clean falls back to the normal rebuild path (worker wipe makes it safe).
- Do not run pipeline scripts locally (HPC-only); NAS-side reads are fine.
