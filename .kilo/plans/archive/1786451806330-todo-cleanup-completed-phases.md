# TODO.md cleanup: remove completed phases, fold rollout into Phase 3.1, commit + push

## Goal

Declutter `TODO.md` (pure documentation edit): delete all completed sections
(Phase 2, Phase 5, Phase 6 incl. the Array 4299048 subsection, Phase 7),
consolidate the benchmark-rollout work into a single new Phase 3.1 point, and
commit **only** `TODO.md` + push. Completed history stays in git.

## Decisions (user-confirmed)

- Include `run_transformation_zeroimp_analysis/1_submit_hpc_array.sh` as a 6th
  remaining command in the new 3.1 (trans + zeroimp still need to run for all
  remaining datasets).
- Renumber Phase 3 items sequentially: new 3.1 (run pipelines) → 3.2 New
  methods → 3.3 Notebook adaptation → 3.4 Docs → 3.5 SLURM config cleanup.
- Commit scope: `TODO.md` only. Do NOT stage `docs/ARCHITECTURE.md`,
  `src/4_cell_type_annotation/1.1_prepare_chunks.py`, or the untracked
  `.kilo/plans/*.md` files (separate in-flight item).
- Keep untouched: `## Human-managed tasks (not agent)`, `## Ideas for later`,
  `## Keep-draft notes`, and Phase 4.

## Edits to TODO.md (apply top-to-bottom; match by exact content, NOT line numbers — line numbers shift after each edit)

### 1. Header intro (lines 3–8)

Update the phase-ordering sentence: drop the `Phase 5` reference, mark Phase 2
done, drop the stale "changelog at the bottom" reference (no changelog exists
in the file). E.g.:

```
Implementation plan for the remaining pipeline work. Phases are ordered:
**Phase 1** (agent, done) → **Phase 2** (debug run on HPC, done) →
**Phase 3** (benchmark methods + full dataset rollout, agent + HPC) →
**Phase 4** (batch effect analysis) → human-managed tasks. Completed history
is preserved in git; see `git log`.
```

### 2. Delete the whole `## Phase 2 — Debug run on HPC [REQUIRES USER: connect NAS + log in to HPC]` section

All bullets `[X]` (incl. the "After debug passes: run one real dataset (e.g.
Kfoury)" bullet — its content is preserved in the new 3.1 context note).
Delete from the `## Phase 2` header through the blank line before
`## Phase 3 — src/5_run_benchmark_methods`.

### 3. Phase 3 — replace the 3.1 and 3.2 bullets with a single new 3.1 bullet

Old: the `- [x] **3.1 Python methods**: ...` bullet AND the
`- [ ] **3.2 R methods**: ...` bullet (both fully implemented + HPC-validated).

New (single bullet, unchecked — rollout in progress):

```
- [ ] **3.1 Run pipelines for all remaining datasets** [IN PROGRESS — USER
      RUNNING ON HPC]: (2026-08-11) `_debug` + Kfoury validated end-to-end
      (Stage 1 → Stage 5, incl. all benchmark methods); implementation +
      syntax/parse checks + HPC debug validation DONE (former 3.1/3.2).
      Remaining commands (preprocess rollout currently running; then):
      - `./src/4_cell_type_annotation/1_prepare_chunks.sh production`
      - `./src/4_cell_type_annotation/2_submit_hpc_array.sh`
      - `./src/4_cell_type_annotation/3_submit_merge.sh`
      - `./src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`
      - `./src/5_run_benchmark_methods/run_r_sample_embedding_methods/1_submit_hpc_array.sh`
      - `./src/5_run_benchmark_methods/run_transformation_zeroimp_analysis/1_submit_hpc_array.sh`
      Dataset coverage — benchmark: Adams, Bassez, Gongsharma_cmv_young_males,
      Kim, Lee, Pelka, Smillie, Stephenson (benchmark view), Wu, Zhang;
      batch-effect: Joanito, Stephenson (batch_effect view), CombinedPBMC.
      Zhu: no views (feeds only the CombinedPBMC stage-2 step 1.2) — stage 2
      only, confirm participation.
      GongSharma cap validation (cap log: 531,291 + 365,000 = 896,291 cells,
      max 5000 per sample) checked when its preprocess task reaches the NAS
      sync gate.
      After all datasets complete: verify NAS outputs (preprocessed h5ads +
      benchmark bundles), then resume 3.2 (new methods), 3.3 (notebook
      adaptation), 3.4 (docs), 3.5 (SLURM config cleanup), and Phase 4.
```

Keep the "Combos run defaults-first" note? No — it was part of the old 3.1
bullet; it is implementation detail now preserved in git history. Drop it.

### 4. Renumber the remaining Phase 3 bullets

- `- [ ] **3.3 New methods**:` → `- [ ] **3.2 New methods**:` (content unchanged)
- `- [ ] **3.4 Notebook adaptation**:` → `- [ ] **3.3 Notebook adaptation**:`
- `- [ ] **3.6 Docs**:` → `- [ ] **3.4 Docs**:`
- `- [ ] **3.7 SLURM config cleanup**:` → `- [ ] **3.5 SLURM config cleanup**:`

### 5. Delete the trailing Phase 3 bullet

`- [ ] Validation: bash -n/py_compile/R parse; debug-dataset run on HPC once
implemented.` — completed (folded into the new 3.1 note).

### 6. Phase 4 — unchanged.

### 7. Delete the whole `## Phase 5 — Annotation completeness guard [agent]` section

Both bullets `[x]`; section fully completed.

### 8. Delete the whole `## Phase 6 — Preprocess array 4294824 recovery [agent code fixes + user HPC steps]` section

Includes the `### Array 4299048 recovery [REQUIRES USER — agents cannot run HPC]`
subsection. All items completed (user confirmed; do not worry about any
unchecked box, e.g. "Re-run Wu"). Delete from the `## Phase 6` header through
the blank line before `## Phase 7`.

### 9. Delete the whole `## Phase 7 — Full dataset rollout [IN PROGRESS — USER RUNNING ON HPC]` section

Absorbed into the new Phase 3.1 bullet. Delete from the `## Phase 7` header
through the blank line before `## Human-managed tasks (not agent)`.

### 10. `## Human-managed tasks (not agent)`, `## Ideas for later`, `## Keep-draft notes` — untouched.

## Commit & push (implementation agent)

1. `git add TODO.md` — ONLY this file.
2. Commit (repo style: short imperative, no body needed), e.g.:
   `TODO cleanup: drop completed Phase 2/5/6/7, fold rollout into Phase 3.1`
3. `git push` (origin, current branch).
4. Do NOT touch the other modified/untracked files.

## Validation

- Re-read `TODO.md`: no `## Phase 2/5/6/7` headers remain; Phase 3 items are
  numbered 3.1–3.5 with no gaps; new 3.1 lists exactly the 6 commands above;
  the three keep-sections are intact.
- `git status --short`: shows `docs/ARCHITECTURE.md`,
  `src/4_cell_type_annotation/1.1_prepare_chunks.py` and the untracked
  `.kilo/plans/` files still uncommitted; `git show --stat HEAD` after push
  shows exactly `TODO.md` (1 file).
