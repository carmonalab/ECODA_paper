# Reorder dataset-specific steps so numbering matches execution order (GongSharma 1.1, runs first)

## Goal

The GongSharma cap step is numbered `1.4` but must run before all other steps
(CombinedPBMC `1.1` in particular). Submission-order + dependency already exist
in `1_submit_hpc.sh` (commit 6e2ff70), but the **file numbering and docs order
still show CombinedPBMC first**. Renumber all step scripts with a full shift so
execution order = numeric order, fix a double-submission bug in the dispatcher,
and update all documentation.

**Decision (user-confirmed): full shift** — GongSharma becomes `1.1`, others shift up.

## Rename map (use `git mv` to preserve history)

| Old | New |
|---|---|
| `1.4_submit_gongsharma.sh` | `1.1_submit_gongsharma.sh` |
| `1.4.1_subset_gongsharma.py` | `1.1.1_subset_gongsharma.py` |
| `1.1_submit_combinedpbmc.sh` | `1.2_submit_combinedpbmc.sh` |
| `1.1.1_create_combinedpbmc_dataset.py` | `1.2.1_create_combinedpbmc_dataset.py` |
| `1.2_submit_joanito.sh` | `1.3_submit_joanito.sh` |
| `1.2.1_prepare_joanito.R` | `1.3.1_prepare_joanito.R` |
| `1.3_submit_kfoury_lowres_ct.sh` | `1.4_submit_kfoury_lowres_ct.sh` |
| `1.3.1_create_kfoury_lowres_ct.R` | `1.4.1_create_kfoury_lowres_ct.R` |

All in `src/2_dataset_specific_preprocessing/`. Rename order matters only for
clarity — do all `git mv` first, then edit contents (paths are resolved via
`SCRIPT_DIR` at runtime, so no collision issues).

## Tasks (ordered)

### 1. `git mv` all 8 files per the map above

### 2. Fix `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`

**Bug fix (double submission)**: the glob loop `for step_script in "${STEP_SCRIPTS[@]}"`
currently submits `1.4_submit_gongsharma.sh` AGAIN (the glob `1.*_submit_*.sh`
matches it) after it was already submitted as the cap prerequisite. Two
concurrent caps race on the deterministic temp path `*.capped_tmp.h5ad` +
`os.replace` in `1.1.1_subset_gongsharma.py` (`cap_file()`), and CombinedPBMC is
gated only on the first — re-opening the race the serialization was designed to
close. In the loop, skip the cap script by name:
`[[ "${step_name}" == "1.1_submit_gongsharma.sh" ]] && continue` (with a comment
explaining it is submitted separately above).

Name updates in the same file:
- Header comment (lines ~12-18): cap step → `1.1_submit_gongsharma.sh`; gated step → `1.2_submit_combinedpbmc.sh`.
- Body comment (lines ~77-83): same name updates.
- `CAP_STEP_SCRIPT="${SCRIPT_DIR}/1.4_submit_gongsharma.sh"` → `1.1_submit_gongsharma.sh`.
- Gate condition: `"1.1_submit_combinedpbmc.sh"` → `"1.2_submit_combinedpbmc.sh"`.
- Keep the `-f` existence guard and `--dependency=afterok` semantics unchanged; `--sync-only` unaffected.

### 3. Internal self-references in renamed files

- `1.1_submit_gongsharma.sh` (ex-1.4): header comment lines ~13/28-31
  (`1.4.1_subset_gongsharma.py` → `1.1.1_subset_gongsharma.py`;
  `1.1_submit_combinedpbmc.sh` → `1.2_submit_combinedpbmc.sh`); line 46 invocation
  → `"${SCRIPT_DIR}/1.1.1_subset_gongsharma.py"`.
- `1.1.1_subset_gongsharma.py` (ex-1.4.1): docstring lines ~30-37 —
  `1.4_submit_gongsharma.sh` → `1.1_submit_gongsharma.sh`,
  `1.1_submit_combinedpbmc.sh` → `1.2_submit_combinedpbmc.sh`,
  usage `1.4.1_subset_gongsharma.py` → `1.1.1_subset_gongsharma.py`.
- `1.2_submit_combinedpbmc.sh` (ex-1.1): line 24 →
  `"${SCRIPT_DIR}/1.2.1_create_combinedpbmc_dataset.py"`.
- `1.2.1_create_combinedpbmc_dataset.py` (ex-1.1.1): docstring usage lines ~5-6 →
  `1.2.1_create_combinedpbmc_dataset.py`.
- `1.3_submit_joanito.sh` (ex-1.2): line 22 →
  `"${SCRIPT_DIR}/1.3.1_prepare_joanito.R"`.
- `1.3.1_prepare_joanito.R` (ex-1.2.1): header lines 2 and 4 →
  `1.3.1_prepare_joanito.R` / `1.3_submit_joanito.sh`.
- `1.4_submit_kfoury_lowres_ct.sh` (ex-1.3): line 22 →
  `"${SCRIPT_DIR}/1.4.1_create_kfoury_lowres_ct.R"`.
- `1.4.1_create_kfoury_lowres_ct.R` (ex-1.3.1): header line 2 → `1.4.1_create_kfoury_lowres_ct.R`.

### 4. Documentation updates

- **`docs/ARCHITECTURE.md`**:
  - Line 82 (dispatcher row): enumerate steps in new order, GongSharma first —
    `1.1_submit_gongsharma.sh, 1.2_submit_combinedpbmc.sh, 1.3_submit_joanito.sh, 1.4_submit_kfoury_lowres_ct.sh`;
    "cap step (`1.4`) ... CombinedPBMC step (`1.1`)" → (`1.1`) / (`1.2`); `1.2.1_prepare_joanito.R` → `1.3.1_prepare_joanito.R`.
  - Line 101 (9 sbatch workers list): new names in new order.
  - Line 102: "GongSharma cap step `1.4` ... gates the CombinedPBMC step `1.1`" → `1.1` / `1.2`.
  - Line 103 (CombinedPBMC bullet): heading + `1.1.1_create_combinedpbmc_dataset.py` → `1.2`/`1.2.1`; "see the 1.4 bullet" → "see the 1.1 bullet".
  - Line 104 (Joanito bullet): `1.2_submit_joanito.sh`/`1.2.1_prepare_joanito.R` → `1.3`/`1.3.1`.
  - Line 105 (Kfoury bullet): `1.3_*` → `1.4_*`.
  - Line 106 (GongSharma bullet): `1.4_submit_gongsharma.sh`/`1.4.1_subset_gongsharma.py` → `1.1`/`1.1.1`; "gates the CombinedPBMC step (1.1)" → "(1.2)".
  - Line 107 (debug subset): `1.2.1_prepare_joanito.R` → `1.3.1_prepare_joanito.R`.
  - Line 723: `1.3.1_create_kfoury_lowres_ct.R` → `1.4.1_create_kfoury_lowres_ct.R`.
- **`AGENTS.md`**: lines 23-25 (9-worker list) new names/order, GongSharma first; line 72 `1.2.1_prepare_joanito.R` → `1.3.1_...`; line 87 both joanito (`1.3.1_...`) and kfoury (`1.4.1_...`) names.
- **`TODO.md`**: line 21 → `1.3_submit_joanito.sh`/`1.3.1_prepare_joanito.R`; lines 80-81 → `1.3.1_...` / `1.4.1_...`; lines 151-152 → `1.1_submit_gongsharma.sh` + `1.1.1_subset_gongsharma.py`.
- **`src/1_stage_data/1_stage_data.sh`** line 16 (comment): `1.2.1_prepare_joanito.R` → `1.3.1_prepare_joanito.R`.

## Do NOT touch

- `.kilo/plans/archive/*` — historical records keep old names (precedent: Kfoury renumber plan 1786047080874).
- `datasets.json` — no references to step scripts.
- Stale `src/2_dataset_specific_preprocessing/__pycache__/` (gitignored, harmless).
- HPC working copy: picks up renames via `git pull` (fileMode handling already configured).

## Validation (static only — no pipeline runs per AGENTS.md)

1. `ls src/2_dataset_specific_preprocessing/` shows contiguous `1.1_*` … `1.4_*`, GongSharma first.
2. `bash -n` on `1_submit_hpc.sh` and the four renamed `*_submit_*.sh` scripts.
3. Repo-wide grep for the 8 OLD names (`1.1_submit_combinedpbmc`, `1.1.1_create_combinedpbmc`, `1.2_submit_joanito`, `1.2.1_prepare_joanito`, `1.3_submit_kfoury`, `1.3.1_create_kfoury`, `1.4_submit_gongsharma`, `1.4.1_subset_gongsharma`) over `src/ docs/ README.md AGENTS.md TODO.md` → matches only in `.kilo/plans/archive/`.
4. Confirm `1_submit_hpc.sh` submits `1.1_submit_gongsharma.sh` once (cap block) and skips it in the loop; gate still on `1.2_submit_combinedpbmc.sh`.

## Risks

- Any missed doc reference is caught by validation step 3 (grep for old names).
- Runtime behavior unchanged apart from the double-submission removal (one cap job instead of two; CombinedPBMC still gated `afterok` on the cap).
- No commit/push unless the user asks (per repo rules); when committing, include this plan archived per the Task Completion Workflow.
