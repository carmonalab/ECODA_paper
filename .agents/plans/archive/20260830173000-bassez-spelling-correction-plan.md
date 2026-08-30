# Correct the Bassez Stage 2 spelling typo before rollout

Plan slug: `20260830173000-bassez-spelling-correction`
Canonical plan artifact: `local://20260830173000-bassez-spelling-correction-plan.md`

## Context

The legacy Bassez step spelling is invalid; the user confirmed that only `Bassez` exists. Current occurrences are in the Stage 2 submitter, the Stage 2 watchdog, and one untracked review note. Correct every occurrence to the canonical internal step identifier `bassez_cellsubtype` before launching any Bassez or rolling-wave pipeline gate.

The correction must preserve the existing hook filename `1.6_submit_bassez.sh`, dataset key `Bassez`, raw output path, owner key behavior, dependency behavior, force behavior, and all unrelated working-tree files. Do not edit `datasets.json`.

## Approach

### 1. Correct the Stage 2 submitter identifiers

In `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`, replace only the legacy Bassez step spelling with `bassez_cellsubtype` at all five current locations:

- `usage()` step list.
- `step_script()` case arm mapping the step to `1.6_submit_bassez.sh`.
- `ALL_STEPS` ordered step list.
- Dataset-to-step filter mapping for `Bassez`.

Retain the existing canonical mapping:

```text
bassez_cellsubtype -> src/2_dataset_specific_preprocessing/1.6_submit_bassez.sh
```

Do not add the legacy Bassez spelling as a compatibility alias. The old spelling must be rejected as an unknown Stage 2 step before any run root, owner, or scheduler state is created.

### 2. Correct the Stage 2 watchdog identifiers

In `src/2_dataset_specific_preprocessing/stage2_watchdog.sh`, replace only the legacy Bassez step spelling with `bassez_cellsubtype` in both current case mappings:

- `stage2_step_script()` must resolve the canonical step to `1.6_submit_bassez.sh`.
- `stage2_step_outputs()` must resolve the canonical step to `${HPC_SCRATCH_DIR}/Bassez/data/BassezA_2021_33958794whole.rds`.

Leave all other Stage 2 step names, output paths, checksums, owner validation, and retry behavior unchanged.

### 3. Correct the explicitly requested review-note occurrence

In `.agents/plans/review_findings.md`, replace the single legacy Bassez step spelling occurrence with `bassez_cellsubtype`, preserving the surrounding finding text and every unrelated user-owned review note. This explicit typo-only correction is authorized by the user's request; do not otherwise rewrite or stage that review file unless the execution workflow explicitly includes it.

### 4. Add behavioral spelling regressions

Extend `tests/test_stage2_submitter.sh` using its existing fake `sbatch`/temporary-home harness:

1. Invoke the submitter with the legacy Bassez step spelling; require nonzero exit, an unknown-step diagnostic, no `sbatch` call, and no newly created run-owned scheduler state.
2. Invoke the submitter with `--datasets Bassez --steps baszez_cellsubtype` only if the test is intentionally checking another malformed spelling; otherwise use the canonical `--datasets Bassez --steps bas(s)ez_cellsubtype` literal exactly as `bassez_cellsubtype`; require successful test-mode submission and a `steps.tsv` row whose step column is exactly `bassez_cellsubtype`, whose script column is `1.6_submit_bassez.sh`, and whose output column is the existing Bassez RDS path.

Use the exact canonical literal `bassez_cellsubtype`; do not introduce a typo in the regression itself. Keep existing CombinedPBMC dependency and force assertions intact.

Extend `tests/test_stage2_watchdog.sh` with a minimal run-owned Bassez manifest fixture using the canonical `bassez_cellsubtype` step and the existing `1.6_submit_bassez.sh`/Bassez RDS output contract. Require the watchdog's manifest validation to accept the canonical row. Mutate only the step field to the legacy Bassez step spelling in a separate fixture and require fail-closed rejection before any retry submission. Preserve the existing Kfoury OOM-retry assertions.

### 5. Make this correction a rollout prerequisite

Do not launch the Bassez Pipeline 2 gate, the hook-backed Stage 2 rest wave, or any rolling Pipeline 3–5 gate until the spelling regressions and repository-wide literal check pass. After the correction is verified, use `bassez_cellsubtype` consistently in any new Stage 2 selection, manifest, owner, or gate evidence; never infer the step name from the script filename.

## Critical files & anchors

- `src/2_dataset_specific_preprocessing/1_submit_hpc.sh` — `usage()`, `step_script()`, `ALL_STEPS`, and dataset filter mapping.
- `src/2_dataset_specific_preprocessing/stage2_watchdog.sh` — `stage2_step_script()` and `stage2_step_outputs()`.
- `tests/test_stage2_submitter.sh` — existing fake scheduler harness for canonical acceptance and old-spelling rejection.
- `tests/test_stage2_watchdog.sh` — existing run-owned watchdog/OOM fixture for canonical and invalid step manifests.
- `.agents/plans/review_findings.md` — one explicitly authorized typo-only review-note correction.

## Verification

Run from the repository root after the edits:

```bash
bash -n src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  src/2_dataset_specific_preprocessing/stage2_watchdog.sh
bash tests/test_stage2_submitter.sh
bash tests/test_stage2_watchdog.sh
```

The submitter test must demonstrate both observable transitions: canonical `bassez_cellsubtype` produces the Bassez hook manifest row and the legacy Bassez step spelling fails before scheduler submission. The watchdog test must demonstrate canonical manifest acceptance and legacy-spelling fail-closed rejection.

Perform a repository-wide case-insensitive literal search for the legacy Bassez step spelling across tracked and untracked project files, including `.agents/plans/review_findings.md`; require zero matches. Separately search for `bassez_cellsubtype` and require matches in both Stage 2 scripts and the Stage 2 submitter regression. Do not treat a passing syntax check without these behavioral checks as sufficient.

Only after all checks pass may the approved Bassez-led rolling pipeline plan be executed. The first Bassez gate must use the corrected Stage 2 identifier and its fresh run-owned manifest; do not reuse a pre-correction run or gate evidence.

## Assumptions & contingencies

- `bassez_cellsubtype` is the canonical Stage 2 internal step name; `Bassez` remains the dataset key and `1.6_submit_bassez.sh` remains the worker script filename.
- If any occurrence of the legacy Bassez step spelling remains, stop before pipeline execution and correct that exact occurrence; do not suppress the search or add an alias.
- If the canonical Bassez fixture fails, inspect the existing `step_outputs()`/owner contract and correct only the fixture or spelling mapping; do not change `datasets.json` or broaden Stage 2 selection.
- Preserve unrelated user WIP and do not stage it as part of this correction. The review-note line is the sole explicitly authorized user-WIP edit.
