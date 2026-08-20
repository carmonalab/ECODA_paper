# De-bloat README.md + encode documentation organization in AGENTS.md

## Goal
Remove the SSH-disconnect/sync-email block from README.md (not relevant for human users) and keep a concise version of that rule set in AGENTS.md, plus add the documentation-organization table to AGENTS.md so future agents stop bloating README.md.

## Context (verified)
- README.md:118–135 is the only place in README.md mentioning `--sync-only`, `USER_EMAIL`, or disconnects.
- AGENTS.md:110 already has a concise SSH-disconnect bullet (drop-safe, `--sync-only` recovery, fail-closed, best-effort emails) — it is missing the `USER_EMAIL` export tip and the list of submitters that support `--sync-only`.
- ARCHITECTURE.md already documents all details (lines 82, 89, 149, 160, 243, 376, 381, 490) — no changes needed there.

## Tasks

1. **README.md — remove the bloated block**
   - Delete lines 118–135 (`**SSH disconnects are safe:** ... non-deliverable.`).
   - No replacement text: the pipeline code block is followed directly by the `- **Benchmark methods**` bullet.
   - Docs-only change; per AGENTS.md rules, do not run any pipeline scripts for validation.

2. **AGENTS.md — extend the existing SSH rule bullet (line 110)**
   - Keep it as one concise bullet. Fold in the two missing pieces:
     - `--sync-only` is supported by `1_submit_hpc.sh`, `3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, `4_cell_type_annotation/2_submit_hpc_array.sh`, and the three benchmark submitters (one comma-separated id per submitted array).
     - Sync-status emails: to receive them, `export USER_EMAIL="you@example.com"` in the login-node shell profile (`~/.bashrc`); otherwise Slurm falls back to `${USER}@unige.ch` which may be non-deliverable.

3. **AGENTS.md — add "Documentation organization" section**
   - Add a new short section (between "# General rules" and "# Repo structure", or directly under "# General rules") containing the user's 4-row table verbatim:
     - Project Summary & Citation → README.md (primary home)
     - Pipeline Call Graphs & HPC Layout → docs/ARCHITECTURE.md (keep here; remove detailed file-by-file logic trees from AGENTS.md)
     - Agent Guardrails & Domain Terms → AGENTS.md (high-level rules; reference ARCHITECTURE.md for step details)
     - Pending Tasks & Method Extensions → TODO.md (centralize planned conversions, script additions, reviewer extensions)
   - This is the guardrail preventing future README bloat. Note: AGENTS.md is currently lean (defers to ARCHITECTURE.md), so no immediate de-bloat work beyond the table.

## Validation
- `git diff` review: confirm README paragraph gone, list structure intact, AGENTS.md bullets render correctly.
- No pipeline scripts (.R/.py/.sh) executed — docs-only change.

## Completion workflow
- Move this plan to `.kilo/plans/archive/`.
- `git add .`; commit summarizing the README de-bloat + AGENTS.md doc-organization rules; push.

## Out of scope
- ARCHITECTURE.md content (already complete).
- TODO.md content.
- Any restructuring of AGENTS.md beyond the two edits above.
