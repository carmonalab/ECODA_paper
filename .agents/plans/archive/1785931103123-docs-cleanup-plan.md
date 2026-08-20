# Documentation Cleanup Plan (README / ARCHITECTURE / AGENTS / TODO)

## 1. Assessment of the external (Gemini) review

**Rating: mostly sound (7/10), needs three refinements.**

What is right:
- README is overloaded: the Workflow + HPC-execution sections (README.md:96–153) duplicate
  ARCHITECTURE.md's usage/role tables. Should be compressed to a quick-start summary + links.
- AGENTS.md is too big (193 lines) and duplicates ARCHITECTURE: "Pipeline Overview"
  (AGENTS.md:68–121) ≈ condensed copy of ARCHITECTURE's preprocessing/annotation sections;
  "HPC folder layout" (AGENTS.md:159–163) duplicates ARCHITECTURE.md:27–64; the Stage-3
  implementation bullets (AGENTS.md:91–101) duplicate TODO.md Phase 3.
- Single-source-of-truth per file type is the right framework.

What the review missed or got wrong:
1. **No mention of stale/factually wrong content** (user approved fixing these):
   - README.md:122–124 claims Python benchmark methods are "invoked automatically via
     rpy2" — WRONG. Verified: rpy2 exists only in `src/utils/preprocess_utils.py` (preprocess
     interop); the notebook reads Python method outputs from `.feather` files
     (benchmark_analysis.rmd:23, 1781). README's annotation quick-start also omits the
     `3.1_submit_merge.sh` step.
   - ARCHITECTURE.md:4 overview calls the repo "an R-based bioinformatics pipeline for
     benchmarking batch effect correction methods" — WRONG focus; it is a compositional
     data analysis / sample-representation benchmark (batch effect analysis is one stage).
   - Cohort counts: README says 11 cohorts (697 samples), AGENTS.md:71 says 12 QC
     notebooks; datasets.json has 14 non-underscore keys (11 `use_for_benchmark`, 3
     `use_for_batch_effect` incl. CombinedPBMC, Zhu view-less). Do NOT hardcode a number —
     phrase as "see datasets.json (ground truth)".
2. **"Remove task-tracking logs" is too blunt**: AGENTS.md:146–156 mixes stale changelog
   (renv removed, heredoc extraction, config_helper move — belongs in git history/TODO
   changelog) with LIVE invariants agents must keep (`PYTHON_BIN`/`PIXI_RSCRIPT`/
   `RETICULATE_PYTHON`, counts-layer input, feather naming, scGate DB cache). Keep the
   invariants, drop the history.
3. **No audit of ARCHITECTURE.md itself**: leftover `# TODO` markers (line 79 heading),
   empty "## Batch Effect Analysis # TODO" section (line 346), duplicated gene-reference
   provenance (line 151 vs README.md:64–70).

## 2. Target single-source-of-truth assignment

| Information type | Sole home |
|---|---|
| Paper summary, citation, key findings, install, quick-start commands | README.md |
| Pipeline call flow, file-role tables, HPC folder layout, env vars, memory baselines, invariants | docs/ARCHITECTURE.md |
| Agent guardrails, domain terminology, repo conventions agents must obey | AGENTS.md (short; pointer-heavy) |
| Task tracking (phases, method extensions, changelog) | TODO.md |

## 3. Ordered tasks

### 3.1 README.md (target ≈ 110–130 lines)
- **Workflow section (96–153)**: compress to one line per stage + a compact quick-start
  command block (stage → `1_submit_hpc.sh` → preprocess array → chunks → annotation array →
  `3.1_submit_merge.sh <DS>`); replace explicit per-script prose with links to
  ARCHITECTURE.md (preprocessing + annotation + usage anchors). Keep the "benchmark methods
  HPC pipeline PLANNED — see TODO.md" note.
- **Fix inaccuracies**:
  - Replace the rpy2 sentence with: Python methods run separately (currently
    `run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`; HPC pipeline
    planned, see TODO.md) and exchange data with the notebook via `.feather` files.
  - Add the missing merge step to the annotation commands.
  - Replace hardcoded "11 scRNA-seq cohorts (697 samples)" wording with a datasets.json
    reference (e.g. "11 benchmark + batch-effect datasets — see datasets.json for the
    authoritative list"). Verify 697 samples claim against the notebook if quick; otherwise
    soften it.
- Keep: title/badge, overview, key findings, repository contents, reference data
  (gene-file provenance stays README-only), installation, expected outputs, citation.
- Note under Stage 3 that rendering in RStudio is the current local flow, to be superseded
  by the HPC pipeline (Phase 3).

### 3.2 docs/ARCHITECTURE.md
- Fix overview sentence (line 4) to describe the actual project focus (sample-representation
  benchmark + batch effect analysis for scRNA-seq cohorts).
- Remove "#### Files # TODO" marker (line 79).
- Replace empty "## Batch Effect Analysis # TODO" (line 346) with 2–3 lines: notebook
  exists, methods to be added per TODO.md Phase 4 (ECODA CT removal, Pseudobulk DESeq2+limma,
  MrVI `batch_key`, GloScope/PILOT-GM-VAE on `X_pca_harmony`).
- Remove the gene-reference provenance dupe from the annotation NAS↔scratch bullet
  (line 151; provenance stays in README "Reference data").
- Add one line noting the benchmark call-flow (Layers 1–5) documents the current
  notebook-based pipeline and will be restructured by TODO Phase 3.
- All other content stays (it is the authority).

### 3.3 AGENTS.md (target ≈ 100–110 lines)
- **Delete** "Pipeline Overview" (68–121) → replace with a 6–8 line stage map (Stage 1–4,
  one line each) with anchors into ARCHITECTURE.md (preprocessing, annotation, benchmark).
  Keep only the agent-relevant one-liners: drafts-keep warning
  (`preprocess_gongsharma.qmd`, `TODO_STUMP_preprocess_sikkema.qmd`), CSR-on-disk invariant
  (pointer to ARCHITECTURE), `_debug` convention (pointer to ARCHITECTURE).
- **Move to TODO.md (dedupe, not duplicate)**: the Stage-3 implementation bullets
  (rename qmd→py, `1_submit_hpc_array.sh`/`1.1_run_worker.sh`, QOT/PULSAR feasibility) —
  already present in TODO Phase 3; leave a single pointer line "pending pipeline work: see
  TODO.md".
- **Prune "not finished yet" (146–156)**: delete stale changelog bullets (renv, heredoc,
  config_helper move — git history + TODO changelog); keep the live worker-environment
  invariants as one-liners: `PYTHON_BIN`/`PIXI_RSCRIPT`/`RETICULATE_PYTHON` from
  `slurm_config.sh` (never bare python/R), counts-layer input for annotation, feather
  naming, scGate DB cache path. Point to ARCHITECTURE for details.
- **Prune HPC section (127–169)**: keep login-node rule, ssh access, VPN note,
  stage/sync-back rules; replace the folder-layout block (159–163) with a pointer to
  ARCHITECTURE.md#hpc-folder-layout.
- **Compress**: "Major reviewer points" (8–24) → 3–4 lines pointing at TODO.md Phase 3/4;
  "Documentation" self-list (37–42) → one line; "R Modules" (123–125) → one line; batch
  effect dataset info (172–176) stays (unique, short, agent-relevant).
- **Fix**: AGENTS.md:73 references `1.2_submit_joanito.sh` → the step file is
  `1.2.1_prepare_joanito.R`; align cohort-count wording to "see datasets.json".
- Keep "General rules", datasets.json rules (compressed, incl. NAS paths), `data/` note,
  Domain Terminology (unchanged — unique to AGENTS.md).

### 3.4 TODO.md
- Compress the Phase-1 "DONE" section (9–102) into ~5 lines (what was done + "details in
  git log / changelog"); it duplicates the Changelog (177–189).
- Confirm Phase 3 items absorb everything removed from AGENTS.md (they already do —
  verify no wording drift, e.g. rename target `1.1.1_benchmark_methods_py.py`).
- Keep Phase 2 checkboxes, Phase 4, human-managed tasks, keep-draft notes, changelog.
- Add one changelog entry for this documentation cleanup.

## 4. Validation (no HPC runs, no pipeline scripts)
1. Re-read the four files; grep key terms (`_debug`, `CSR`, `X_pca_harmony`, rpy2,
   `1.1.1_benchmark_methods_py`, `3.1_submit_merge`) and confirm each appears in exactly
   one authoritative doc (dupes allowed only as pointers).
2. Check all internal anchors (`docs/ARCHITECTURE.md#...`) referenced from README/AGENTS
   exist (section ids are auto-generated from headings — keep heading wording when
   referencing).
3. Line-count sanity: AGENTS.md ≤ ~110, README ≤ ~130, ARCHITECTURE grows/shrinks
   minimally (moves, not rewrites).

## 5. Out of scope
- Restructuring ARCHITECTURE's benchmark Layers 1–5 content (still accurate for the
  notebook pipeline; Phase 3 will rewrite it).
- Any code/script changes; any datasets.json changes; HPC runs.
