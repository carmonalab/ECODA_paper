# TODO.md update: `_debug` + Kfoury fully validated; full rollout running

## Goal

Update `TODO.md` to reflect: (1) `_debug` AND Kfoury ran end-to-end
successfully through all five pipelines (Stage 1 → Stage 5, incl. benchmark
methods), (2) Phase 3.1/3.2 HPC debug validation is now DONE, (3) Phase 6
Kfoury-related recovery items are resolved (env repair + Kfoury re-run),
(4) the remaining Phase 6 re-run/sync items are IN PROGRESS, covered by a new
"Phase 7 — Full dataset rollout" section currently running on HPC.

This is a pure documentation edit of `TODO.md`. No source changes, no commits.

## Edits (apply in order)

### 1. Phase 2 — replace the final unchecked item (line 30)

Old:
```
- [ ] After debug passes: run one real dataset (e.g. Kfoury) before full rollout. (partly done, src/5_run_benchmark_methods/run_r_sample_embedding_methods has now to be run on HPC after re-installing a clean environment (after fix))(started running the full pipeline on all datasets, see phase 6, still some debugging necessary, currenltly running src/3_scrnaseq_preprocessing and debugging)
```

New:
```
- [X] After debug passes: run one real dataset (e.g. Kfoury) before full rollout.
      (2026-08-11: `_debug` AND Kfoury ran end-to-end successfully through all
      five pipelines — Stage 1 staging, Stage 2 dataset-specific preprocessing,
      Stage 3 preprocess, Stage 4 annotation + merge, Stage 5 benchmark methods
      incl. R methods (GloScope/MOFA/Pseudobulk/scITD), transformation/
      zero-imputation analyses, and Python methods (MrVI/scPoli/PILOT). Full
      rollout on all remaining datasets is running — see Phase 7.)
```

### 2. Phase 3.1 — mark HPC debug validation DONE

Old (end of the 3.1 bullet):
```
  CODE COMPLETE (slurm_config.sh benchmark vars + ARCHITECTURE.md updated);
  HPC debug validation PENDING — smoke test:
  `./1_submit_hpc_array.sh --ds_name _debug --methods mrvi` (then scpoli, pilot);
  check `benchmark/embeddings/` feathers + `execution_times.feather`.
```

New:
```
  CODE COMPLETE (slurm_config.sh benchmark vars + ARCHITECTURE.md updated);
  HPC debug validation DONE (2026-08-11): `--ds_name _debug --methods mrvi`
  (then scpoli, pilot) + full Kfoury run — `benchmark/embeddings/` feathers +
  `execution_times.feather` confirmed.
```

Keep the "Combos run defaults-first" note (lines after) unchanged.

### 3. Phase 3.2 — mark HPC debug validation DONE

Old (end of the 3.2 bullet):
```
  target; see ARCHITECTURE.md). HPC debug validation PENDING — smoke test:
  Pipeline B first (`--ds_name _debug --analysis trans,zeroimp`), then Pipeline A
  `--ds_name _debug --methods prepare_pseudobulk,pseudobulk`; check
  `benchmark/results/`, `benchmark/pseudobulks/`, `execution_times.feather` on NAS.
```

New:
```
  target; see ARCHITECTURE.md). HPC debug validation DONE (2026-08-11):
  Pipeline B (`--ds_name _debug --analysis trans,zeroimp`) then Pipeline A
  (`--ds_name _debug --methods prepare_pseudobulk,pseudobulk`) + full Kfoury
  run — `benchmark/results/`, `benchmark/pseudobulks/`,
  `execution_times.feather` confirmed on NAS.
```

### 4. Phase 6 — item "(X) Joanito preprocessing needs to be run again" (line 24)

Append a note that the rerun is covered by the Phase 7 rollout, e.g.:
`([X] Joanito preprocessing needs to be run again — covered by the Phase 7 rollout, see below)`

### 5. Phase 6 — HPC manual steps: re-run/sync items → IN PROGRESS

Old:
```
- [ ] Re-run failed datasets individually (each run syncs its own outputs):
      `--ds_name` Joanito, Smillie, Wu, Zhang, Stephenson. (needs re-running after fixes)
- [ ] Sync the 9 COMPLETED-but-unsynced datasets by re-running their `--ds_name`
      (already-processed outputs are skipped, then synced): Adams, Bassez,
      CombinedPBMC, Kfoury, Kim, Lee, Pelka.
```

New:
```
- [ ] (IN PROGRESS — 2026-08-11, covered by the Phase 7 full rollout) Re-run
      failed datasets individually (each run syncs its own outputs):
      `--ds_name` Joanito, Smillie, Wu, Zhang, Stephenson.
- [ ] (IN PROGRESS — 2026-08-11, covered by the Phase 7 full rollout; Kfoury
      DONE) Sync the COMPLETED-but-unsynced datasets by re-running their
      `--ds_name` (already-processed outputs are skipped, then synced): Adams,
      Bassez, CombinedPBMC, Kim, Lee, Pelka.
```

### 6. Phase 6 — Array 4299048 recovery items

- Env repair (Kfoury numpy): `[ ]` → `[X]` — 2026-08-11: py-cuda13 env
  re-installed clean; Kfoury numpy corruption resolved (Kfoury ran end-to-end).
  Keep the "NEVER run env installs while an array is active" warning.
- "Re-run Wu" bullet: mark as covered by the Phase 7 rollout (keep the verify
  checklist: `WuS_2021_34493872_benchmark_analysis_ECODAprocessed.h5ad` on NAS).
- "Re-run Kfoury" bullet: `[ ]` → `[X]` — 2026-08-11: full end-to-end run
  COMPLETED, `Kfoury_2021_34719426_benchmark_analysis_ECODAprocessed.h5ad`
  synced to NAS.
- "Re-run the 12 COMPLETED-but-unsynced datasets" bullet: mark as covered by
  the Phase 7 rollout (keep the dataset list + optional full-array note).

### 7. GongSharma OOM fix note (end of Phase 6, line ~227)

Update "(HPC validation pending: cap log must show ...)" to state validation
is pending within the Phase 7 rollout (GongSharma is included in the current
run; check the cap log + NAS sync gate when its preprocess task completes).

### 8. New "Phase 7 — Full dataset rollout" section

Insert after the "Array 4299048 recovery" subsection (before
"## Human-managed tasks"):

```
## Phase 7 — Full dataset rollout [IN PROGRESS — USER RUNNING ON HPC]

(2026-08-11) `_debug` + Kfoury validated end-to-end (Stage 1 → Stage 5, incl.
all benchmark methods). Full rollout currently running for all remaining
datasets:

- Benchmark datasets: Adams, Bassez, Gongsharma_cmv_young_males, Kim, Lee,
  Pelka, Smillie, Stephenson (benchmark view), Wu, Zhang.
- Batch-effect datasets: Joanito, Stephenson (batch_effect view),
  CombinedPBMC.
- Zhu: no views (feeds only the CombinedPBMC stage-2 step 1.2) — confirm
  participation in the rollout (stage 2 only).

Notes:
- The rollout also covers the Phase 6 re-run/sync items (each run skips
  existing outputs and syncs its own; Kfoury already synced by its full run).
- GongSharma cap validation (cap log: 531,291 + 365,000 = 896,291 cells,
  max 5000 per sample) checked when its preprocess task reaches the NAS sync
  gate.
- After all datasets complete: verify NAS outputs (preprocessed h5ads +
  benchmark bundles), then resume Phase 3.4 (notebook adaptation), 3.3 (new
  methods: PILOT-GM-VAE, QOT, PULSAR), 3.6 (docs), 3.7 (SLURM config cleanup),
  and Phase 4 (batch effect analysis).
```

## Open items / out of scope (explicitly left unchanged)

- Phase 3.3, 3.4, 3.6, 3.7 and Phase 4 remain unchecked — still pending.
- Phase 5, Phase 6 code fixes: already `[X]`, untouched.
- `datasets.json`, README, ARCHITECTURE.md, AGENTS.md: untouched.

## Validation

- Re-read `TODO.md` after edits: no stale "PENDING" benchmark-validation text
  for 3.1/3.2; Phase 2 and Kfoury items all checked; Phase 7 section present.
- `git diff --stat` shows only `TODO.md` changed.
