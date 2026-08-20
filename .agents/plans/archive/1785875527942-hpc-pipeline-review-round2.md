# HPC Pipeline Review — Round 2 (src/2, src/3, src/4)

Scope: `src/2_dataset_specific_preprocessing/`, `src/3_scrnaseq_preprocessing/`,
`src/4_cell_type_annotation/`, `src/1_stage_data/`, `src/slurm_config.sh`,
`config_helper.R`, `docs/ARCHITECTURE.md`, `README.md`, `TODO.md`, `datasets.json`.

Prior review round (`.kilo/plans/1785866906173-hpc-pipeline-review-fixes.md`) is
fully implemented (commits 8202f92…6917bef) — its Phase-1 items were verified present.

---

## Findings (verified against code)

### A. Bug (factual error)

**A1. `src/1_stage_data/1_stage_data.sh` never stages Zhu (view-less dataset).**
The jq program iterates `views` entries only (`to_entries[] | .value.views | ...`).
Zhu has `"views": {}` in datasets.json, so `ZhuH_2023_37379396whole.rds` is never
rsynced to `${HPC_SCRATCH_DIR}/Zhu/data` — confirmed by running the exact jq locally
(Zhu absent from output; only 13 datasets emitted). But
`1.1.1_create_combinedpbmc_dataset.py` loads Zhu from the *dataset-level*
`file_names` in per-dataset layout (`${HPC_SCRATCH_DIR}/Zhu/data`) →
the CombinedPBMC step fails on HPC with a file-not-found inside `load_input()`.
→ Fix in Phase 1, step 1.

### B. Artefacts / dead code

**B1. Draft qmd files in `src/3_scrnaseq_preprocessing/`:**
`preprocess_gongsharma.qmd` and `TODO_STUMP_preprocess_sikkema.qmd` are
INTENTIONAL drafts for future implementation (Gongsharma: additional subsetting
conditions; Sikkema: future Lung dataset) — USER DECISION: keep in place, do not
delete. Nothing references them from scripts (grep: only TODO.md:163,367). They
must be explicitly marked as drafts in TODO.md and noted in AGENTS.md so they are
not mistaken for dead code. TODO.md:163 ("Cleanup preprocess_gongsharma.qmd if
necessary") is misleading and should be reworded.

**B2. `src/3_scrnaseq_preprocessing/__pycache__/1.1.1_preprocess.cpython-313.pyc`**
— regenerated local artifact (gitignored); delete again.

**B3. Dead variables:**
- `PIXI_R_LIB` (`slurm_config.sh:27`) — no consumers anywhere (R_LIBS_SITE
  wiring was removed when workers moved to `${PIXI_RSCRIPT}` = pixi env activation).
- `SLURM_ACCOUNT=""` (`slurm_config.sh:52`) — never used.
- `config_helper.R`: `force_overwrite` parameter and `path_output_chunks` list
  entry — no consumers (`1.1_prepare_chunks.py` and `2_submit_hpc_array.sh`
  construct the chunks dir themselves).

**B4. `README.md:146`** presents
`src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`
as an existing workflow step; the folder only contains `1.2_benchmark_methods_py.qmd`
(scripts are planned, not implemented).

### C. Inconsistencies across pipelines

**C1. Partition hardcoded in 5 places** (TODO Item B, still open):
`#SBATCH --partition=shared-cpu` in `1.1_submit_combinedpbmc.sh`,
`1.2_submit_joanito_batch_col.sh`, `1.1_run_worker.sh`, `2.1_run_worker.sh`;
`srun --partition=shared-cpu` in `1_prepare_chunks.sh` and
`2_submit_hpc_array.sh`. `SLURM_PARTITION` in slurm_config.sh stays unused.

**C2. `3_merge_annotations.py` is an orphan** (TODO Item C, still open):
no bash wrapper calls it; ARCHITECTURE.md:118/158 presents it as a pipeline stage
with a manual invocation.

**C3. Preprocess array has no single-dataset mode.** Chunk prep
(`1_prepare_chunks.sh <mode> <DS>`) and annotation submit (`2_submit_hpc_array.sh <DS>`)
support a dataset argument; the preprocess array `1_submit_hpc_array.sh` always
submits all datasets (its worker already supports `--ds_name`, but the task→dataset
index mapping is over all `jq keys[]`). Blocking for debug/test mode (Phase 2 Item A).

**C4. Annotation runs on every preprocessed view → double annotation.**
`1.1_prepare_chunks.py` chunks every h5ad in the output dir (e.g. Stephenson has
benchmark + batch-effect views; batch view is a cell superset of the benchmark
view). Overlapping samples are annotated twice, and `3_merge_annotations.py` then
merges BOTH views' feathers into EACH h5ad (safe via the `(sample, barcode)` join +
`drop_duplicates(keep="first")`, but wasteful and confusing). Needs a decision
(Phase 2, Item D).

**C5. Feather sample column naming.** `2.1.1_process_chunk.R` writes the sample
column as hardcoded lowercase `"sample"`; `3_merge_annotations.py` also hardcodes
`"sample"` for the feather side while reading `SAMPLE_COLNAME` ("Sample") for the
h5ad side. Works, but inconsistent — align both sides to `SAMPLE_COLNAME`.

**C6. `1.2.1_create_joanito_batch_col.R` hardcodes** `"Joanito"` + the input file
name instead of reading datasets.json. Acceptable for a dataset-specific step, but
a one-line `read_datasets_json()` lookup (or comment) would remove the duplicate
ground-truth source. Document-only, low priority.

### D. Robustness / verify items

**D1. `1.1.1_preprocess.py:213` skips views whose output file already exists** —
partial/corrupt/stale outputs (e.g. killed mid-write, or regenerated after a code
change) are silently never recomputed. Consider a `--force` flag (Phase 2, Item E).

**D2. Cluster verify items** (cannot be checked locally): CombinedPBMC step 64G
baseline; preprocess worker 16G for GongSharma (>100k cells); annotation worker
`--time=02:00:00` vs 5×2 retry timeouts; `aux/scGateDB.rds` is now committed to
git (8.8 KB) so the `2_submit_hpc_array.sh` "download once" `srun` path is mostly
bypassed — the comment/flow could mention the committed cache.

### E. Verified OK (factual correctness)

- Dispatcher glob `1.*_submit_*.sh` does not recurse into `1_submit_hpc.sh`.
- Array task↔dataset mapping consistent (`jq 'keys[]'` sorted + `sed -n` in worker).
- CSR-on-disk guarantee, `counts`-layer annotation input, feather naming from chunk
  file, `(sample, barcode)` merge key, sacct fail-closed gates before NAS sync,
  single-dataset sync in annotation submit — all match docs.
- Chunk counter is global across views within one dataset dir (`chunk_<N>.txt`
  continue across h5ad files) → no feather-name collision across views; chunk file
  line 1 pins the source h5ad.
- `apply_subset_vars` before sample-col standardization; empty-subset guard.
- `slurm_config.sh` vars match consumers (`CHUNKS_MANIFEST`, `SCGATE_DB_PATH`,
  `SAMPLE_COLNAME`, `PYTHON_BIN`, `PIXI_RSCRIPT`, `MAX_NUM_CHUNKS_PARALLEL`).

### F. Debug/test mode status (analysis for TODO.md — out of scope for implementation)

**NOT implemented.** No `_create_debug_dataset.R`, no `data/debug/`, no `_debug`
entry in datasets.json. TODO.md Step 2 is marked stale; TODO.md Item A (Phase 2)
exists and covers the plan at a high level.

Existing partial test support:
- `1_prepare_chunks.sh test <DS>` → `1.1_prepare_chunks.py --test` (1 sample/chunk;
  deliberately does NOT delete stale feathers).
- Preprocess worker `--ds_name` (but submit array covers all datasets — gap C3).
- Single-dataset args: chunk prep + annotation submit; `3_merge_annotations.py` is
  manual and works on any dir.

Gaps for a full debug run:
- `1_stage_data.sh` has no dataset filter (stages everything) and reads from NAS —
  debug files must be on NAS or staged manually (scp/rsync from `data/debug/`).
- `2_dataset_specific_preprocessing` steps run on FULL datasets (the Joanito
  seqtec step would process the full Joanito .rds, not the debug subset — debug
  subset must be created/staged as a standalone file).
- Preprocess array: no single-dataset mode (C3).
- `_debug` registration decision (datasets.json vs separate debug_datasets.json):
  adding `_debug` to datasets.json makes every `keys[]`-iterating script (staging,
  both arrays, NAS sync loops) see it — a separate `debug_datasets.json` or an
  array filter is needed.

---

## Phase 1 — Direct fixes (small; implement now)

### 1. Fix Zhu staging — `src/1_stage_data/1_stage_data.sh`
Replace the jq program so that view-less datasets emit their dataset-level
`file_names` (existing `sort -u` dedups overlaps with view inputs):

```bash
jq -r '
  to_entries[] |
  .key as $key |
  .value.folder_name as $folder |
  (if (.value.views | length) > 0 then .value.views | to_entries[] | .value.input_file_name else .value.file_names end) |
  if type == "array" then .[] else . end |
  "\($key) \($folder) \(.)"
' "${DATASETS_JSON_FILE}" | sort -u | ...
```

Sanity expectations after fix: Zhu now emitted; CombinedPBMC still skipped
(`folder_name == null` guard); Stephenson/Gongsharma rows unchanged (dedup).

Optional hardening in `1.1.1_create_combinedpbmc_dataset.py`
(`load_and_prepare_source`): before `load_input`, check
`(base_path / ds_name / "data" / name).exists()` (per-dataset layout) and raise a
clear "run 1_stage_data.sh first" error — improves DX, low risk.

### 2. Keep draft qmds, remove other artefacts, document drafts
- KEEP `src/3_scrnaseq_preprocessing/preprocess_gongsharma.qmd` and
  `TODO_STUMP_preprocess_sikkema.qmd` (USER DECISION: drafts for later
  implementation — do NOT delete).
- Delete `src/3_scrnaseq_preprocessing/__pycache__/` (regenerated artifact).
- TODO.md:163 — reword "Cleanup preprocess_gongsharma.qmd if necessary" to
  clearly mark the file as a KEPT DRAFT (to be implemented/cleaned when the
  GongSharma other-subsetting conditions are added, see "New Datasets to Be
  Added").
- TODO.md:367 — keep the Sikkema pointer, clarify draft status: "Possibly:
  Sikkema Lung (n = 165) — draft in
  `src/3_scrnaseq_preprocessing/TODO_STUMP_preprocess_sikkema.qmd` (to be
  implemented)".
- AGENTS.md — add a short note in the repo-structure section for
  `src/3_scrnaseq_preprocessing/`: both qmds are drafts for future
  implementation, not dead code.

### 3. Remove dead variables
- `slurm_config.sh`: delete `PIXI_R_LIB` (line 27) and `SLURM_ACCOUNT` (line 52;
  add back when an account is actually needed). Keep the comment on
  `SLURM_PARTITION` (used in step 4).
- `config_helper.R`: drop the `force_overwrite` parameter and the
  `path_output_chunks` list entry (verify no consumers — grep first).

### 4. Centralize SLURM partition (closes TODO Item B)
- Remove `#SBATCH --partition=shared-cpu` from `1.1_submit_combinedpbmc.sh`,
  `1.2_submit_joanito_batch_col.sh`, `1.1_run_worker.sh`, `2.1_run_worker.sh`.
- Pass `--partition="${SLURM_PARTITION}"` on the sbatch command line from:
  `1_submit_hpc.sh` (per step), `3_scrnaseq_preprocessing/1_submit_hpc_array.sh`,
  `4_cell_type_annotation/2_submit_hpc_array.sh` (array + the scGate DB `srun`),
  and the `srun` in `1_prepare_chunks.sh`.
- `SLURM_PARTITION` stays sourced (not exported) — all these scripts source
  slurm_config.sh in the same shell before invoking sbatch/srun.
- TODO.md Item B: mark done, keep the note that benchmark methods will later need
  their own partition (GPU) via `SLURM_PARTITION` override (TODO Step 6).

### 5. Feather sample column alignment (optional, do it — touches both sides together)
- `2.1.1_process_chunk.R`: `annot$sample <- target_sample` →
  `annot[[args$sample_colname]] <- target_sample` (already env-driven).
- `3_merge_annotations.py`: read the annotations-column name from
  `os.environ.get("SAMPLE_COLNAME", "Sample")` instead of hardcoded `"sample"`
  (lines 38, 42, 48).
- Keep the join logic unchanged.

### 6. README + docs updates
- README.md:146 — mark the python-methods `1_submit_hpc_array.sh` step as
  planned/not yet implemented (matches TODO.md).
- ARCHITECTURE.md — update only where behavior changed (partition centralization,
  staged files incl. Zhu, feather column naming); keep the `3_merge_annotations.py`
  manual-invocation presentation until Item C is decided.
- AGENTS.md — only if anything above contradicts it (no PIXI_R_LIB mentions found;
  verify after step 3).
- TODO.md — tick off completed items; add Phase 2 sections below.

---

## Phase 2 — TODO.md additions (large; out of scope for this implementation round)

### Item A — Debug/test execution mode (rewrite/expand existing Item A + stale Step 2)
Add a consolidated section with:
- **Status**: not implemented (no `_create_debug_dataset.R`, no `data/debug/`, no
  `_debug` in datasets.json); Step 2 marked stale.
- **Coverage matrix** (script → today → needed adaptation):
  | Script | Test support today | Gap for debug run |
  |---|---|---|
  | `1_stage_data/1_stage_data.sh` | none (stages all, NAS-only) | dataset filter; debug file must exist on NAS or be staged manually |
  | `2_dataset_specific_preprocessing/*` | none (full datasets) | debug subset is standalone; Joanito seqtec step still processes full Joanito — decide skip vs run |
  | `3_scrnaseq_preprocessing/1_submit_hpc_array.sh` | none (all datasets) | single-dataset/`--ds_name` array mode (maps DS→jq index, 1-task array) |
  | `3_scrnaseq_preprocessing/1.1_run_worker.sh` | `--ds_name` via preprocess.py | none |
  | `4_cell_type_annotation/1_prepare_chunks.sh` | `test <DS>` + `--test` (1 sample/chunk) | none |
  | `4_cell_type_annotation/2_submit_hpc_array.sh` | `<DS>` single-dataset | none |
  | `3_merge_annotations.py` | manual | none |
- **Wiring checklist**: `_create_debug_dataset.R` (5 Joanito samples, 500
  cells/sample, minimal obs cols incl. `seqtec`; outputs `.rds` + `.h5ad`),
  registration decision (`_debug` in datasets.json vs `debug_datasets.json` —
  NOTE: `_debug` in datasets.json pollutes every `keys[]` loop), preprocess array
  filter, validation (< 30 s, keys `X_pca`/`X_pca_harmony`, shape ~2500 × genes).

### Item C — Integrate `3_merge_annotations.py` (keep, add recommendation)
Recommend a `3.1_submit_merge.sh` wrapper (per-dataset, `--ds_name` arg, srun
compute node, `PYTHON_BIN`), run after the annotation array completes; alternative
is hooking it into `2_submit_hpc_array.sh` after the sync gate. Note the Item D
interaction (per-view feather sets).

### Item D (new) — Annotation scope decision: per-view vs per-dataset
Document C4: annotating every preprocessed view double-annotates overlapping
samples (Stephenson); `3_merge_annotations.py` merges all views' feathers into
each h5ad. Decide: keep per-view (simple, current) vs annotate once per dataset on
the union and reuse across views.

### Item E (new, optional) — Preprocess idempotency
`1.1.1_preprocess.py` "Already processed" skip (D1): add `--force` flag or write
via temp file + rename so partial outputs are never skipped.

### Verify items (cluster)
CombinedPBMC 64G baseline; preprocess worker 16G (GongSharma); annotation worker
2 h / 16 G; scGate DB committed-cache note (`2_submit_hpc_array.sh` comment).

---

## Validation (no pipeline runs — per AGENTS.md)

- `bash -n` on all touched `.sh`; `python -m py_compile` on touched `.py`;
  R parse (`Rscript -e 'parse(file=...)'`) on touched `.R`.
- jq smoke test for the new staging query: Zhu present, CombinedPBMC absent,
  output deduplicated (`sort -u`), count 14.
- Greps: no `PIXI_R_LIB`, `SLURM_ACCOUNT`, `--partition=shared-cpu` outside
  slurm_config.sh; both draft qmds still present; new TODO.md/AGENTS.md draft
  notes in place; `__pycache__` absent.
- Local debug-dataset sanity (optional, if `data/debug/` exists): preprocess
  `--ds_name` on the debug h5ad in `--force`-free clean dir.

## Open decisions (left for user; non-blocking for Phase 1)

1. (RESOLVED) The two draft qmds stay in place — drafts for later
   implementation; documented in TODO.md + AGENTS.md (Phase 1 step 2).
2. Item A: `_debug` in datasets.json vs separate `debug_datasets.json`.
3. Item C: dedicated merge wrapper vs hook into `2_submit_hpc_array.sh`.
4. Item D: annotation per-view vs per-dataset.
