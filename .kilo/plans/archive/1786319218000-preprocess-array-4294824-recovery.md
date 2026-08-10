# Recover preprocessing array 4294824 (5 failure classes) + GongSharma memory fix

## Context

Array 4294824 (2026-08-10, `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`,
14 datasets, 128G/16c per task): 9 tasks COMPLETED, 5 failed. Gate correctly
failed closed → **nothing synced to NAS**. Task→dataset mapping (jq `keys[]`
order): 1 Adams, 2 Bassez, 3 CombinedPBMC, 4 Gongsharma_cmv_young_males,
5 Joanito, 6 Kfoury, 7 Kim, 8 Lee, 9 Pelka, 10 Smillie, 11 Stephenson,
12 Wu, 13 Zhang, 14 Zhu.

| Task | Dataset | State | Root cause (validated) |
|------|---------|-------|------------------------|
| 4 | Gongsharma_cmv_young_males | OUT_OF_MEMORY 0:125 | Two large SoundLife h5ad inputs fully loaded + concat exceed 128G |
| 5 | Joanito | FAILED | anndataR `write_h5ad` tags string cols `encoding='ascii'` even with non-ASCII bytes (NBSP `0xa0` in `Stage.TNM`) → `sc.read_h5ad` UnicodeDecodeError. **Reproduced locally** (anndataR roundtrip fails the same way). |
| 10 | Smillie | FAILED | scanpy `seurat_v3_paper` loess singular fit: quantized `log10(mean)` in small per-sample batches → `ValueError: b'reciprocal condition number 7.8464e-15'`. **Reproduced locally** with pinned scanpy 1.12.2 on synthetic integer counts. |
| 11 | Stephenson | RUNNING at gate | Not a code failure: array left squeue, sacct still reported RUNNING (stale/orphaned record). Verify + re-run. |
| 12 | Wu | FAILED | anndataR `validate_aligned_mapping`: layer colnames not `identical()` to assay features (mismatch in the middle of the gene list; head/tail identical). Not reproducible synthetically — Seurat normalizes layers through the public API, so the quirk lives inside the RDS. Needs defensive alignment + diagnostics. |
| 13 | Zhang | FAILED | Same loess singularity as Smillie (`There are other near singularities as well. 0.090619`). |
| 14 | Zhu | FAILED | R env corruption: `mime` lazyload DB missing (`…/py-cuda13/lib/R/library/mime/R/mime`); env has mixed R package builds (SeuratObject 4.5.1 / sp 4.5.3 vs R 4.5.2). |

## Fix A — loess jitter retry in `select_hvgs_ranked` [agent]

File: `src/3_scrnaseq_preprocessing/1.1.1_preprocess.py` (`select_hvgs_ranked`, ~line 19).

- Wrap the `sc.pp.highly_variable_genes(..., layer="counts", flavor=..., batch_key=..., check_values=True)` call in `try/except ValueError`.
- On `ValueError`: build a deterministic jittered copy of the counts layer:
  `rng = np.random.RandomState(42)`; `cnt = adata.layers["counts"].copy()`;
  `cnt.data = (cnt.data + rng.random_sample(cnt.data.shape).astype(np.float32) * 1e-6)`
  (counts is guaranteed CSR after `base_preprocessing`, so `.data` exists);
  store as `adata.layers["counts_jittered"]`.
- Retry once with `layer="counts_jittered"`, `check_values=False`; on a second
  `ValueError`, delete the jitter layer and re-raise the ORIGINAL exception.
- On success: assign `adata.var["hvg_rank"]`, print a clear warning
  ("seurat_v3 loess fit failed (…); retried with deterministic jittered counts"),
  and `del adata.layers["counts_jittered"]` so the jitter layer is never written
  to the output h5ad.
- Same behavior for `batch_key=None` and `batch_key=<col>` (both hit the same loess).
- Validation (already done during diagnosis): synthetic 100-cell zero-inflated
  integer counts fails with the pinned scanpy; the jittered retry succeeds and
  returns 3000 HVGs + `highly_variable_rank` (identical downstream contract).

## Fix B — meta.data ASCII sanitization in RDS→h5ad conversion [agent]

File: `src/utils/py/preprocess_utils.py` (`ro.r` block defining
`convert_rds_to_raw_h5ad`).

- Inside `convert_rds_to_raw_h5ad`, after `create_clean_seuratv5_object()` and
  before `write_h5ad()`, sanitize character `meta.data` columns to pure ASCII:
  `md[] <- lapply(md, function(x) if (is.character(x)) iconv(x, from = "latin1", to = "ASCII", sub = " ") else x)`.
  (anndataR records `encoding='ascii'` regardless of content; any non-ASCII byte
  then breaks `sc.read_h5ad`. NBSP → space; other non-ASCII → space. Clinical/
  sample metadata is ASCII-safe; lossy replacement is deliberate.)
- Cache note: the guard `if (!file.exists(output_path))` means a BROKEN cached
  `JoaI_…_raw.h5ad` on scratch would block regeneration. Hardening implemented
  in the fix commit (review follow-up): `load_single_input` catches
  `UnicodeDecodeError` on a cached `*_raw.h5ad`, deletes it and re-converts
  once — the manual HPC `rm` below is now belt-and-braces (keep it anyway).
  Sanitization covers character AND factor/ordered meta.data columns
  (factor-ness preserved), and the sanitize/align blocks only run when the
  cache is missing.

## Fix C — Wu: defensive layer alignment + diagnostics [agent]

File: `src/utils/py/preprocess_utils.py` (`convert_rds_to_raw_h5ad`).

- After `create_clean_seuratv5_object()`, if the RNA assay is `Assay5`, align
  every layer's rownames to the assay features before writing:
  - `features <- rownames(rna)`; for each `lyr` in `names(rna@layers)`: if
    `!identical(rownames(m), features)`, `message()` the count of differing
    positions (and the first ~20 differing gene names), then reindex
    `m <- m[features, , drop = FALSE]` and write back into `rna@layers[[lyr]]`.
  - Reindexing fails loudly if features ∉ rownames (duplicates/missing genes) —
    that is the desired clear failure (identifies the RDS inconsistency).
- No-op for the rebuilt path (fresh `CreateSeuratObject` layers are aligned by
  construction); protects the unchanged-original path (counts/data/X fallbacks
  all empty) where the mismatch originates.

## GongSharma memory fix — cap 5000 cells/sample via dataset-specific step [user implements]

Historical approach (git: `3a4711e`, `src/py/preprocess_gongsharma.qmd`):
`downsample_by_group(adata, group_key='specimen.specimenGuid', max_cells=5000, seed=42)`
per input file, then concat; `Sample = specimen.specimenGuid`.

New files (auto-discovered by `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`):
- `src/2_dataset_specific_preprocessing/1.4_submit_gongsharma.sh` — sbatch
  (1h / 16c / 128G), the standard SCRIPT_DIR-recovery block (scontrol `Command=`
  under sbatch, `BASH_SOURCE` fallback), `source slurm_config.sh`, `cd ${PROJECT_ROOT}`,
  then `${PYTHON_BIN} "${SCRIPT_DIR}/1.4.1_subset_gongsharma.py"`.
- `src/2_dataset_specific_preprocessing/1.4.1_subset_gongsharma.py` — read the two
  staged SoundLife h5ads from `${HPC_SCRATCH_DIR}/Gongsharma_cmv_young_males/data/`
  (names from `datasets.json` `file_names`), per-sample cap to 5000 cells
  (`np.random.RandomState(42)`), print before/after cell counts + n samples,
  write capped data.

Output decision (pick one):
- (a) **In-place overwrite of the two staged h5ads** — zero config, no
  `datasets.json` change; caveat: any re-run of `1_stage_data.sh` restores the
  raw files, so the subset step must be re-run after re-staging (the dispatcher
  order guarantees this for a normal run).
- (b) Write a single combined `Gongsharma_cmv_young_males_capped5000.h5ad` and
  update `datasets.json` `file_names`/view `input_file_name` — cleaner end
  state, but datasets.json is ground truth → requires explicit approval.

## HPC manual steps [REQUIRES USER — no agents on HPC]

1. Repair the py-cuda13 R env (login node; installs allowed):
   `~/.pixi/bin/pixi run -e py-cuda13 Rscript -e 'install.packages("mime", repos = "https://cloud.r-project.org")'`
   then `~/.pixi/bin/pixi run -e py-cuda13 setup` (idempotent). ABI warnings
   (SeuratObject/Matrix/sp built under other R versions) are benign. This also
   protects the annotation/benchmark workers that use the same env.
2. Delete the broken Joanito conversion cache:
   `rm "${HOME}/scratch/ECODA_paper/Joanito/output/JoaI_2022_35773407_Nofilt_whole_raw.h5ad"`
   (regenerates sanitized after Fix B).
3. Check Stephenson task 11: `sacct -j 4294824_11 --format=JobID,State,Elapsed,End -X`
   + `ls "${HOME}/scratch/ECODA_paper/Stephenson/output/"`; re-run
   `--ds_name Stephenson` if outputs missing.
4. After Fixes A–C are committed: re-run failed datasets individually
   (`./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name <DS>` for
   Joanito, Smillie, Wu, Zhang, Stephenson) — each run syncs its own outputs.
5. Sync the 9 COMPLETED-but-unsynced datasets by re-running their `--ds_name`
   (outputs are skipped as already-processed, then synced): Adams, Bassez,
   CombinedPBMC, Kfoury, Kim, Lee, Pelka.
6. After the GongSharma subset step is implemented: run
   `./src/2_dataset_specific_preprocessing/1_submit_hpc.sh`, then
   `./src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh --ds_name Gongsharma_cmv_young_males`.
