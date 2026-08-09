# Plan: Simplest Tasks from TODO.md General Section (Lines 1–36)

## Context

- **Completed (prior session):** README.md partially updated (repo contents section); AGENTS.md extended with pipeline overview, R module table, domain terminology; ARCHITECTURE.md created.
- **Paper (now read):** "Cell type composition drives patient stratification in single-cell RNA-seq cohorts" — Halter, Andreatta, Carmona (2026). 11 cohorts, 697 samples. ECODA (CLR-transformed cell-type proportions) outperforms all state-of-the-art methods (MrVI, scPoli, PILOT, GloScope, GloProp, MOFA+, scITD) for unsupervised patient stratification. Performance driven by 5–18 HVCs (12–29% of cell types). Robust to batch effects and annotation strategies. Orders of magnitude faster. scECODA R package at github.com/carmonalab/scECODA.
- **Remaining scope:** Lines 2–36 — installation/usage docs, paper context, AGENTS refinement, partitioning strategy, naming conventions, script consolidation.

## Identified Simplest Tasks (ordered by complexity)

---

### Task 1: Add Installation & Usage Instructions to README.md

**Why simplest:** Pure documentation. All info already in `pixi.toml` + directory structure. Zero code risk.

Steps:
1. Add an **Installation** subsection to README.md (after "Repository Contents"):
   - Prerequisite: [Pixi](https://pixi.sh) package manager
   - Clone repo, then `pixi install`
   - Platform selection: default env (`py-cpu`) for macOS/development; `py-cuda13` for HPC with CUDA
   - `pixi run setup` — chains 7 R package install tasks (Seurat → anndataR → SignatuR → scATOMIC deps → scITD deps → HiTME deps → all R methods)
2. Add a **Usage / Workflow** subsection:
   - Stage 1 (QC): Open `QC_filtering/*.Rmd` in RStudio, render per dataset
   - Stage 2 (Preprocessing): `pixi run -e py-cpu python src/py/preprocess.py` (driven by `datasets.json`)
   - Stage 3 (Benchmark): Render `benchmark_analysis.rmd` in RStudio; Python methods auto-invoked via rpy2
   - Stage 4 (Batch Effect): Render `batch_effect_analysis.rmd`
   - HPC: `sbatch src/bash/cell_type_annotation/*.sh`; data copy via `src/bash/copy_data_from_nas_to_hpc_scratch.sh`
3. Add an **Expected Outputs** subsection:
   - `.feather` files from Python methods (distance matrices / embeddings)
   - Figures: MDS plots, PCA biplots, benchmark bar charts, separation metric tables
   - Execution time logs
4. Add a concise summary of the same to AGENTS.md (line 9 directive)

---

### Task 2: Add Paper-Derived Context to README.md and AGENTS.md (lines 10–21)

**Why simple:** README.md already has an "Overview" and "Key Findings" section from the prior session. AGENTS.md has domain terminology. This task fills gaps using the now-available paper content.

Steps:
1. **README.md** — verify the existing Overview/Key Findings sections are accurate vs the paper. If missing, add:
   - The 11 cohorts and 697 samples (can reference `datasets.json`)
   - The benchmarked methods list (MrVI, scPoli, PILOT, GloScope, GloProp, MOFA+, scITD, Pseudobulk, ECODA)
   - The batch effect robustness finding (Joanito dataset: ECODA batch ANOSIM=0.041 vs biological ANOSIM=0.640; pseudobulk: batch ANOSIM=0.706)
2. **AGENTS.md** — evaluate lines 17–21 ("as kilo code uses code base indexing…"):
   - Function signatures, folder structure, module dependencies: **skip** — kilo code indexes these automatically.
   - Core concepts (ECODA, CLR, HVCs, zero imputation, separation metrics): **already present** in Domain Terminology section.
   - Add a brief note: "Kilo Code indexes code structure and function signatures automatically. AGENTS.md focuses on domain concepts, pipeline logic, and project conventions that indexing cannot infer."
3. No additional paper-derived sections needed in AGENTS.md — the Domain Terminology and Pipeline Overview sections already cover the key concepts.

---

### Task 3: Optimize Information Partitioning Between README.md, AGENTS.md, and TODO.md (lines 25–31)

**Why simple:** Analysis + minor edits. No code changes.

Steps:
1. Compare current content across all three files:
   - **README.md** = human-facing: paper overview, installation, usage, outputs, citation. ✓ mostly done.
   - **AGENTS.md** = agent-facing: pipeline overview, R module table, domain terminology, HPC notes, datasets.json constraint, reviewer priorities. ✓ comprehensive.
   - **TODO.md** = actionable tasks only. Currently contains some descriptive text (lines 38–51: "Explain new cell type annotation pipeline", "Explain migration from R to Python"). These are more documentation tasks than TODO items — but they are self-contained sections, so leave them.
2. Check for redundancy:
   - Pipeline overview: present in both README.md ("Repository Contents") and AGENTS.md ("Pipeline Overview"). README has the high-level summary; AGENTS.md has the detailed 4-stage breakdown. This is the correct pattern — no change needed.
   - Domain terminology: only in AGENTS.md. Correct — agents need it, humans can reference the paper.
   - datasets.json: only in AGENTS.md. Correct — agents need the constraint, humans can read the file itself.
3. Result: current partitioning is sound. No moves needed. Document this finding as a note in the plan (no file edits for this task).

---

### Task 4: Audit & Document File Naming Conventions (line 24)

**Why simple:** Observational only. List files, check patterns, produce recommendations.

Steps:
1. List all key files by type:
   - `.rmd`: `benchmark_analysis.rmd`, `batch_effect_analysis.rmd` — lowercase_with_underscores ✓
   - `.qmd`: `benchmark_methods_py.qmd`, `preprocess_gongsharma.qmd`, `DRAFT_BARE_preprocess_sikkema.qmd` — lowercase ✓ (except DRAFT prefix)
   - `.py`: `preprocess.py` — lowercase ✓
   - `.R`: all in `src/R/` — lowercase_with_underscores ✓
   - `.sh`: `copy_data_from_nas_to_hpc_scratch.sh` — lowercase ✓
   - Directories: `QC_filtering/` (mixed case), `src/` (lowercase), `docs/` (lowercase)
2. Identify inconsistencies:
   - `QC_filtering/` uses uppercase — but this is a well-known acronym, acceptable.
   - `DRAFT_BARE_preprocess_sikkema.qmd` has uppercase prefix — intentional draft marker, acceptable.
   - `Process_data.ipynb` (referenced in TODO.md line 90) — if it still exists, check its naming.
   - `Preprocess_datasets.Rmd` (referenced in TODO.md lines 41, 47) — check if it still exists or was removed during migration.
3. Check if any old file names remain that should have been renamed:
   - `MAIN_Analysis.Rmd` → `benchmark_analysis.rmd` (done)
   - `Batch_effect.rmd` → `batch_effect_analysis.rmd` (done)
   - `Preprocess_datasets.Rmd` — check if it still exists
   - `Process_data.ipynb` — check if it still exists
4. Document findings. If inconsistencies are found, recommend specific renames; if naming is already consistent, mark line 24 as resolved.

---

### Task 5: Consolidate Duplicate or Redundant Scripts (line 35)

**Why last (highest risk):** May require understanding code semantics. Identification phase is safe; actual consolidation could be deferred.

Steps:
1. Check for duplicates:
   - `src/py/preprocess.py` vs `preprocess_gongsharma.qmd` — likely different (one is general, one is dataset-specific). Verify.
   - `benchmark_methods_py.qmd` vs any `.ipynb` — check if old notebooks remain.
   - R functions across `src/R/` files — grep for shared function names to detect duplication.
2. Check for orphaned/deprecated files:
   - `Preprocess_datasets.Rmd` — if it exists, is it superseded by `preprocess.py`?
   - `Process_data.ipynb` — is it superseded by `benchmark_methods_py.qmd`?
   - Any other files in the root that are not referenced in README.md or AGENTS.md?
3. For each finding: either propose a concrete consolidation (keep canonical, remove duplicate, update references) or document that the file serves a distinct purpose.
4. If no true duplicates: mark line 35 as resolved in TODO.md.

---

## Execution Order

1. **Task 1** — Installation/Usage in README.md (lowest risk, highest reader value)
2. **Task 2** — Paper context verification + AGENTS.md kilo code note
3. **Task 3** — Partitioning audit (quick check, likely no changes needed)
4. **Task 4** — File naming audit (observational, may trigger rename recommendations)
5. **Task 5** — Script consolidation audit (highest complexity, saved for last)
6. **Final pass** — Update TODO.md to mark resolved items

## Validation

- After Task 1: read final README.md, verify markdown renders correctly
- After Task 2: confirm AGENTS.md Domain Terminology and Pipeline Overview are consistent with paper
- After Task 5: if R files are touched, run `Rscript -e "source('src/R/load_all_functions.R')"` to verify loading
- Cross-check: no information in AGENTS.md contradicts README.md
