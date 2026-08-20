# Fix Wu dataset preprocessing failure (job 4299048, task 12) — hidden `names` attribute on Assay5 features

## Context

Preprocess array 4299048 (src/3_scrnaseq_preprocessing): 12/14 tasks COMPLETED,
2 FAILED → gate failed closed → nothing synced to NAS.
- Task 6 (Kfoury): `FileNotFoundError: .../numpy/_pytesttester.py` — py-cuda13
  env corruption (documented class, AGENTS.md); other tasks imported numpy fine
  at the same time → partial file state, likely an env mutation racing with the
  array. SEPARATE from Wu; fix = env repair + re-run.
- Task 12 (Wu): `Error in private$.validate_aligned_mapping(value, "layers", ...)`
  "Expected column names ... Provided column names ..." (head/tail identical in
  the truncated cli output).

## Root cause (Wu) — fully diagnosed and reproduced locally

The Wu RDS (`/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/WuS_2021_34493872/output/WuS_2021_34493872.rds`,
385 MB, 98,294 cells, 20,828 genes, 26 samples, single Assay5 "RNA" + "counts" layer)
carries a hidden **`names` attribute on its Assay5 features**:
`attr(features_df, "dimnames")[[1]]` is a named character vector
(`c(LINC00115 = "LINC00115", ...)`). All expected data (Sample, subtype,
celltype_major, celltype_subset, counts) is present and correct — only this
attribute is the defect.

Failure chain (verified against anndataR@07612e4, the pinned SHA, and
SeuratObject 5.2.0):
1. `GetAssayData`/`LayerData` propagate the named features into the layer
   matrix rownames.
2. anndataR `from_Seurat` builds `var_names` from a data.frame rownames
   (plain, UNNAMED character vector) and transposes the layer
   (`to_py_matrix`), so layer colnames carry the names attribute.
3. `.validate_aligned_array` uses `identical()` (names-attribute sensitive):
   `identical(colnames(layer), var_names)` → FALSE although all elements are
   equal (`identical(unname(colnames(tp)), vn)` → TRUE). Element-wise equal,
   attribute differs → cli-truncated message shows identical head/tail → all
   previous fixes (layer reindexing/transposing, commit ec85cb2) targeted the
   wrong hypothesis.

Verified fix (direct R, pinned env): `dimnames(counts) <- list(unname(rownames(counts)), unname(colnames(counts)))`
before `CreateSeuratObject` → `write_h5ad` succeeds. `unname()` is a no-op for
clean objects (all other datasets unaffected).

Local env note: `import rpy2.robjects` segfaults in ALL local envs (default,
py-cpu) on macOS — pre-existing ABI issue, unrelated to Wu (HPC rpy2 works;
task 12 failed there with the R error, not a crash). Local validation of the R
side must use direct `Rscript`; the python-side pipeline can be validated with a
pre-built raw h5ad (bypassing rpy2) or on HPC.

## Code changes (implementation agent)

1. **`src/utils/seurat_utils.R` — `create_clean_seuratv5_object()`** (the fix;
   only callers are the rds→h5ad converter and the diagnostic script):
   - After the counts fallback chain (`counts` layer / `@layers[["X"]]` / `data`
     layer), immediately before `CreateSeuratObject`:
     ```r
     # Strip stray 'names' attributes from dimnames (Wu quirk: the source RDS
     # carries a names attribute on its Assay5 features; it propagates into
     # layer dimnames and makes anndataR's identical() layer validation fail
     # with a misleading validate_aligned_mapping error). unname() is a no-op
     # when the names attribute is absent.
     dimnames(counts) <- list(unname(rownames(counts)), unname(colnames(counts)))
     ```
   - After the list-column drop, right before `CreateSeuratObject`:
     ```r
     rownames(md) <- unname(rownames(md))
     ```
   - Note: `unname(NULL)` is `NULL`, so the fallback paths (dimnames-less
     matrices) are unaffected.

2. **`src/utils/py/preprocess_utils.py` — R block `align_assay_layers`**
   (defense in depth, no-op for clean objects): normalize the canonical
   feature vector once at function start (`features <- unname(rownames(a))`)
   and use it for all comparisons/reindexing/transposing.

3. **`src/3_scrnaseq_preprocessing/diagnose_layer_alignment.R`**: extend
   `describe_layer()`/`inspect()` to print whether `names(rownames(m))` /
   `names(colnames(m))` / `names(rownames(a))` are NULL — so this quirk is
   caught immediately in future diagnostics (the current script reported
   `identical=TRUE`, `setequal=TRUE` and missed it).

4. **`TODO.md`**: update Phase 6 — mark the Wu root cause (names attribute)
   and code fixes; keep the HPC re-run checkboxes open (below).

No changes to `datasets.json`, no NAS modifications (user confirmed: pipeline
fix only, NAS file stays pristine).

## Validation

### Local (Mac, pinned default env; direct Rscript — rpy2 is broken locally)

1. Full-flow check on the real RDS (proves the fix at the exact failure
   point):
   ```
   .pixi/envs/default/bin/Rscript -e '
     library(Seurat); library(anndataR)
     source("src/utils/seurat_utils.R")
     so <- readRDS("/Volumes/Shared/DataCollections/Standardized_SingleCell_Datasets/WuS_2021_34493872/output/WuS_2021_34493872.rds")
     so <- create_clean_seuratv5_object(so)
     write_h5ad(so, "/tmp/wu_fixed_raw.h5ad")
   '
   ```
   Expect: succeeds (before the fix: `validate_aligned_mapping` error).
2. Quirk-subset regression (2000 random cells; subsetting strips the quirk, so
   re-inject it via `attr(slot(a, "features"), "dimnames")[[1]] <- setNames(dn[[1]], dn[[1]])`
   before `saveRDS`) → same flow → succeeds; sanity-check that the UNFIXED
   converter still fails on it (faithful repro).
3. Clean-dataset sanity (unname is a no-op): run the same flow on one more
   RDS (e.g. Joanito `JoaI_2022_35773407`) → succeeds unchanged.
4. Optional python-side check: feed the fixed raw h5ad (step 1 output) into
   `1.1.1_preprocess.py`'s `process_view()` locally, bypassing the rpy2
   `load_input` path (load the h5ad directly), to confirm PCA/HVG/Harmony
   steps on Wu data. Not required if HPC validation passes.

### HPC (user steps — agents cannot run HPC)

1. Verify/repair the py-cuda13 env (Kfoury):
   - `ls .pixi/envs/py-cuda13/lib/python3.13/site-packages/numpy/_pytesttester.py`
     + `pixi run -e py-cuda13 python -c "import numpy; print(numpy.__version__)"`
   - If the file is missing or numpy is corrupt: env repair per AGENTS.md
     (`src/utils/bash/setup_env_sbatch.sh` first; definitive:
     `rm -rf .pixi/envs/py-cuda13 && pixi install -e py-cuda13 && pixi run -e py-cuda13 setup`).
     NEVER run env installs while an array is active.
2. Re-run Wu (validates the fix end-to-end, syncs its outputs):
   `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name Wu`
   → expect task COMPLETED + "synced to NAS" email; verify
   `WuS_2021_34493872_benchmark_analysis_ECODAprocessed.h5ad` on NAS, loads in
   scanpy with `X`, `layers["counts"]`, `X_pca_*` / `X_pca_harmony_*` keys.
3. Re-run Kfoury: `--ds_name Kfoury`.
4. Re-run the 12 COMPLETED-but-unsynced datasets (each run skips existing
   outputs and syncs): Adams, Bassez, CombinedPBMC, Gongsharma_cmv_young_males,
   Joanito, Kim, Lee, Pelka, Smillie, Stephenson, Zhang, Zhu
   (optionally: full array without `--ds_name`).

## Delivery (per AGENTS.md Task Completion Workflow)

After implementation + local validation: archive this plan to
`.kilo/plans/archive/`, `git add .`, commit (message style: `Preprocess Wu fix:
strip names attribute on Assay5 features (validate_aligned_mapping root cause)`),
push.

## Open questions / risks

- None blocking. Risks: (a) macOS rpy2 segfault blocks local python-path
  validation — mitigated via direct-R + HPC validation; (b) Kfoury numpy
  corruption may recur if env refreshes run during arrays — covered by the
  guardrail + repair steps.
