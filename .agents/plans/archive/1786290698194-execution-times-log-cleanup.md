# Execution-time log cleanup: deterministic naming + no residue accumulation

## Context (findings)

- Per-task exec logs `execution_times_task_<JOBID>_<TASKID>.feather` are current-pipeline
  artifacts, but older runs' logs (job ids 4291150/4291184/4291191/4294608/4294614/4294616)
  were never cleaned: `benchmark_merge_sync_cleanup` deletes only the **current** run's job
  ids (`benchmark_submit_common.sh:162-164`), so stale scratch logs get re-rsynced to the NAS
  on every run. The NAS copy of `benchmark/embeddings/` currently holds 10 stale task logs.
- Job-id naming existed because concurrent arrays share task ids 1..N. METHOD/ANALYSIS +
  DS_NAME are already unique per task and deterministic across runs → rename to
  `execution_times_<method>_<ds>.feather` so re-runs overwrite automatically.
- **Per-combo bundle filenames have a redundant method segment**:
  `results/<ds>_<method>_<combo>.rds` where every combo name is already method-prefixed
  (`GloScope_hvg2000_pcadims30`, `MOFA_hvg2000_factors2`, `Pseudobulk_hvg2000`,
  `scITD_hvg2000_factors2`) → `_debug_gloscope_GloScope_hvg2000_pcadims30.rds` duplicates
  the method. Fix: per-combo cache files become `<ds>_<combo>.rds` (six construction sites in
  `benchmark_pipeline.R`: lines 772, 850, 933, 961, 993, 1059).
- The **bundle keys** (combo names, e.g. `GloScope_hvg2000_pcadims30`) are a public API and
  MUST stay unchanged: the notebook accesses them literally in ~100 places (e.g.
  `result_list$bmark$Pelka$Pseudobulk_hvg2000$dist_mat`, benchmark_analysis.rmd:1231) and
  they double as the exec-log `method` column + plot labels. Only the FILE names change.
- NOT renamed (no duplication there): method-level `results/<ds>_<method>.rds` (the notebook
  loads exactly this via `load_hpc_benchmark_results`),
  `pseudobulks/<ds>_pseudobulk_<variant>.rds` (variant names like `hvg2000` have no method
  prefix; benchmark_hpc_utils.R:94,151,171 + 1.1.1_prepare_pseudobulk.R:76,87),
  `gloscope_dists/<ds>_gloscope_hvg<n>_pcadims<d>_dists.rds` (benchmark_pipeline.R:794 — raw
  GloScope distance cache, consumed only by the worker on cache hits; keep the lowercase
  infix, known cosmetic asymmetry vs the new bundle names).

## Design

- Worker per-task log: `execution_times_<METHOD>_<DS_NAME>.feather` (R/python workers) and
  `execution_times_<ANALYSIS>_<DS_NAME>.feather` (trans/zeroimp worker). Uniqueness across
  concurrent arrays holds (each array has a distinct METHOD/ANALYSIS); determinism means
  re-runs overwrite the same file.
- Merge scoped to the run's (method × dataset) cross product (preserves the "no stale logs
  from failed runs leak in" guarantee that `--job_ids` provided).
- Post-rsync cleanup deletes **this run's** `execution_times_<label>_<ds>.feather`
  per-task logs (label x dataset cross product loop; never matches the merged
  `execution_times.feather`), plus a separate legacy `execution_times_task_*`
  sweep (no current worker produces that naming, so it only hits stale files).
  Scoping means an overlapping submission's not-yet-merged logs are never
  deleted (review fix; the original delete-all glob was changed after review).
- Per-combo result caches: `results/<ds>_<combo>.rds` (method infix dropped; combo names are
  method-prefixed). Method-level skip-if-exists (`<ds>_<method>.rds`) and bundle keys unchanged.

## Tasks

1. **Rename per-task logs in the 3 workers**
   - `src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1_run_worker.sh:49`
   - `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1_run_worker.sh:49`
   - `src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1.1_run_worker.sh:53`
   - Change `LOG_FILE="${OUT_DIR}/execution_times_task_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.feather"`
     to `LOG_FILE="${OUT_DIR}/execution_times_${METHOD}_${DS_NAME}.feather"`
     (trans/zeroimp: `${ANALYSIS}` instead of `${METHOD}`).
   - Rewrite the explanatory comment: job-id no longer needed — concurrent arrays each have a
     distinct METHOD/ANALYSIS, and the deterministic name is overwritten on re-runs.

2. **Rework `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py`**
   - Replace `--job_ids` with `--labels <name>...` and `--datasets <ds>...`; glob the cross
     product `execution_times_<label>_<ds>.feather` instead of
     `execution_times_task_<jobid>_*.feather`. No back-compat needed (all 3 callers updated).
   - `--labels` without `--datasets` defaults to all datasets
     (`execution_times_<label>_*.feather`); `--datasets` without `--labels` errors
     (review fix: the original `[""]` fallback glob matched nothing, contradicting
     the help text).
   - Keep: `--existing-log` NAS continuity, dedup `(dataset, method)` keep="last",
     `--cleanup`/`--no-cleanup` (submitters keep passing `--no-cleanup`).
   - Update module docstring.

3. **Update `src/5_run_benchmark_methods/benchmark_submit_common.sh`**
   - Change `benchmark_merge_sync_cleanup <job_ids...>` → `benchmark_merge_sync_cleanup <labels...>`:
     build merge args `--labels "${@}" --datasets "${DATASET_NAMES[@]}"`; job ids are no longer
     needed by this function (still used by the wait/gate loop in submitters).
   - Replace the per-job-id cleanup loop (`:162-164`) with a scoped delete
     after successful rsync: loop the run's `--labels` x `DATASET_NAMES`
     cross product removing `execution_times_<label>_<ds>.feather`, plus a
     separate legacy sweep `rm -f .../execution_times_task_*.feather`
     (review fix: the original single `execution_times_*.feather` glob could
     delete another overlapping submission's unmerged logs).
   - Update the header comment (lines 23-31).

4. **Update the 3 submit scripts' merge call**
   - `run_r_sample_embedding_methods/1_submit_hpc_array.sh:233` → `benchmark_merge_sync_cleanup "${METHODS[@]}"`
     (METHODS includes auto-prepended `prepare_pseudobulk`).
   - `run_python_sample_embedding_methods/1_submit_hpc_array.sh:183` → `benchmark_merge_sync_cleanup "${METHODS[@]}"`.
   - `run_transformation_zeroimp_analysis/1_submit_hpc_array.sh:147` → `benchmark_merge_sync_cleanup "${ANALYSES[@]}"`.

5. **Deduplicate per-combo result bundle filenames** (`src/5_run_benchmark_methods/benchmark_pipeline.R`)
   - Change all six bundle_file constructions from `paste0(ds, "_<method>_", nm, ".rds")` to
     `paste0(ds, "_", nm, ".rds")`: lines 772 (gloscope), 850 (mofa), 933, 961, 993
     (pseudobulk ×3 loops), 1059 (scitd). Every `nm` is method-prefixed, so the infix is
     redundant.
   - Do NOT touch: bundle keys/`results[[nm]]` names, method-level `<ds>_<method>.rds`
     (1.1.1_run_benchmark_methods_r.R:65), `pseudobulks/` and `gloscope_dists/` names.
   - Update the driver header comment (benchmark_pipeline.R:672-674: "per-combo cache files
     (<ds>_<combo>.rds, skip-if-exists unless --force)") and the worker header comment
     (1.1.1_run_benchmark_methods_r.R:12-14).

6. **Docs**
   - `docs/ARCHITECTURE.md`: update lines 227, 239-242, 349, 363, 367, 466 (per-task log name
     `execution_times_<METHOD>_<DS>.feather`, merge args `--labels/--datasets`, cleanup-deletes-all
     semantics); grep for `_gloscope_`/`_mofa_`/`_pseudobulk_`/`_scitd_` per-combo naming
     references and update to `<ds>_<combo>.rds`.
   - Update header comments in the 3 `1_submit_hpc_array.sh` files (they describe the cleanup
     tail) and any per-combo naming mentions.

## One-time debris cleanup (manual, run by user on HPC login node after `git pull`)

```bash
# 1) HPC scratch stale exec logs (also auto-removed by the new cleanup on the next run)
rm -f "${HPC_SCRATCH_DIR}/benchmark/embeddings"/execution_times_task_*.feather

# 2) NAS stale exec logs (keep the merged execution_times.feather)
rm -f /srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper/benchmark/embeddings/execution_times_task_*.feather

# 3) Convert old per-combo bundle files <ds>_<method>_<Combo...>.rds -> <ds>_<Combo...>.rds
#    (mv preserves the expensive caches; the regex's [A-Za-z] guard leaves the method-level
#    <ds>_<method>.rds and _trans/_zeroimp files untouched — [A-Z] alone would MISS the
#    scITD combos, whose names start lowercase: scITD_*). Run on scratch AND on the NAS
#    results dir. For the _debug dataset, deleting instead of mv'ing is equivalent (fast recompute).
for D in "${HPC_SCRATCH_DIR}/benchmark/results" \
         "/srv/smednas515.unige.ch/carmona_smb/Projects/ECODA_paper/benchmark/results"; do
  for f in "${D}"/{,_}*_{gloscope,mofa,pseudobulk,scitd}_*.rds; do
    [[ -e "${f}" ]] || continue
    new="$(basename "${f}" | sed -E 's/^(.+)_(gloscope|mofa|pseudobulk|scitd)_([A-Za-z].*)$/\1_\3/')"
    mv "${f}" "${D}/${new}"
  done
done
# Note: the NAS checksums.md5 will still list old-style paths until the next run rewrites it —
# harmless (verify_md5_sidecar looks entries up per loaded file; method-level files are untouched).
```

The merged `execution_times.feather` (43 rows, deduped) is correct and stays; the notebook
reads it via the local `data/execution_times.feather` copy (refreshing that local copy is
optional).

## Validation

1. Rerun `./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh --ds_name _debug --methods gloscope,mofa,pseudobulk` (and a python-methods debug run).
2. Submit log shows: merge message with per-method/dataset logs, post-sync "Deleted
   per-task execution-time logs", and `embeddings/` on both scratch and NAS containing only
   `execution_times.feather`.
3. Run the same command twice back-to-back; confirm no `execution_times_*` per-task files
   accumulate and the merged feather row count stays stable (28 rows for the debug combo set
   + historical rows via `--existing-log`).
4. `results/` after the rerun contains `_debug_GloScope_hvg2000_pcadims30.rds`,
   `_debug_MOFA_hvg2000_factors2.rds`, `_debug_Pseudobulk_hvg2000.rds` etc. (no
   `_gloscope_`/`_mofa_`/`_pseudobulk_` infix), plus unchanged method-level
   `_debug_gloscope.rds`/`_debug_mofa.rds`/`_debug_pseudobulk.rds`; `pseudobulks/` and
   `gloscope_dists/` names unchanged.
5. Notebook smoke check: `load_hpc_benchmark_results()` (bundle keys unchanged) and the
   exec-time chunk still work.

## Risks / notes

- Bundle keys/combo names must NOT be renamed — notebook accesses them literally
  (benchmark_analysis.rmd:1231 and ~100 method-mapping sites); they also appear in the
  exec-log `method` column and figure labels. Only filenames change.
- Concurrent submissions of the two benchmark submitters are (and were) unsupported:
  the shared merged feather + full-dir rsync already interleave. The post-rsync
  cleanup is scoped to each run's own logs (review fix), so overlapping runs do
  not delete each other's per-task logs, but the merged-feather/rsync race remains. Document in ARCHITECTURE.md if not already.
- `--existing-log` continuity is unchanged, so partial `--ds_name _debug` runs still extend
  the full NAS log.
- DS_NAME values come from `datasets.json` keys (`[A-Za-z0-9_-]`, incl. `_debug`) —
  filename-safe, no sanitization needed.
- Known cosmetic asymmetry (accepted): `gloscope_dists/` keeps the lowercase
  `<ds>_gloscope_hvg<n>_pcadims<d>_dists.rds` naming (no duplication there); aligning it to
  `_GloScope_*` is possible later but out of scope.
