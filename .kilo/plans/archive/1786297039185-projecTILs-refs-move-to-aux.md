# Remove orphaned `ProjecTILs_references/` + add Figshare ref-map download fallback

## Verdict

**The local `ProjecTILs_references/` folder is orphaned — delete it.** Verified:
- No tracked code references the path (`grep -rn "ProjecTILs_references"` across R/sh/py/qmd/rmd/toml → zero matches).
- The HPC pipeline never touches the repo folder: it uses `NAS_REF_DIR`/`HOME_REF_DIR` (`src/slurm_config.sh:64-65`), staged by `src/4_cell_type_annotation/1_prepare_chunks.sh:65-71` into `$HOME/reference_atlases/sketched_200ct/` on the login node.
- The only `.gitignore` mention (`ProjectTILs_references/`, line 56) is a **dead typo'd rule** ("Project" vs "Projec") that never matched; the folder was ignored only via `*.rds`.

**Second request: pipeline must fetch ref maps NAS-first, Figshare-fallback.** Today `1_prepare_chunks.sh` only rsyncs from NAS (and hard-fails under `set -e` if NAS is unavailable). It needs a Figshare download fallback.

**Verified facts for the fallback** (Figshare API v2, article 26310994, "Light reference atlases in blood and tumors", CC BY 4.0, Garnica/Carmona/Andreatta 2024 — matches README "Reference data"):
- The article contains exactly the 4 files the pipeline needs (`2.1.1_process_chunk.R:95-98`), each with a stable download URL + MD5. **Local mirror MD5s match the API MD5s exactly** → URLs/checksums below are canonical:

| File | Figshare file id / URL | MD5 |
|---|---|---|
| `sketched_CD8T_human_ref_v1.rds` | `47714158` / `https://ndownloader.figshare.com/files/47714158` | `be86058ddafdd0154faf0485286b86e7` |
| `sketched_CD4T_human_ref_v2.rds` | `47714155` / `https://ndownloader.figshare.com/files/47714155` | `5540a0ee287e291528c96d476794b194` |
| `sketched_DC_human_ref_v2.rds` | `47714161` / `https://ndownloader.figshare.com/files/47714161` | `033d491ba7ca9bbf0badcae828e55b2c` |
| `sketched_MoMac_human_v1.rds` | `47714164` / `https://ndownloader.figshare.com/files/47714164` | `3043cd9058a8746d972c7be195b18e36` |

## Implementation steps

### 1. Delete the orphaned folder
`rm -rf ProjecTILs_references/` (untracked + ignored → no `git rm` needed, no git status noise).

### 2. `.gitignore` cleanup
Remove the dead typo'd line 56 (`ProjectTILs_references/`). No other change (`*.rds` stays as backstop; `!aux/*` untouched — no move into `aux/` anymore).

### 3. `src/4_cell_type_annotation/1_prepare_chunks.sh` — replace the staging block (lines 61-71)
Keep the NAS-first preference, add Figshare fallback, and fail closed. New logic (script already runs on the login node where NAS is mounted; note `set -euo pipefail` is active at line 12, so **rsync must be guarded with `if`/`||`**, not left unguarded):

1. Define the manifest as bash associative arrays (bash ≥4, standard on the cluster login node):
   - `REF_MAP_FILE_IDS["<filename>"]="<figshare file id>"`
   - `REF_MAP_MD5S["<filename>"]="<md5>"`
   - plain array `REF_MAP_NAMES=("sketched_CD8T_human_ref_v1.rds" "sketched_CD4T_human_ref_v2.rds" "sketched_DC_human_ref_v2.rds" "sketched_MoMac_human_v1.rds")`
2. **Skip if already staged**: if all 4 files exist in `HOME_REF_DIR` → message + skip (keeps current idempotence; also fixes the partial-state hole of the old non-empty-dir check).
3. **NAS first**: `if rsync -a "${NAS_REF_DIR}" "${HOME_REF_DIR}/"; then ...` — on failure print a warning and continue to fallback (do NOT abort under `set -e`).
4. **Figshare fallback**: for each file still missing in `HOME_REF_DIR`, download to a temp path inside `HOME_REF_DIR` (`curl -f -L --retry 3 -o "${HOME_REF_DIR}/.${f}.tmp.$$" "<url>"`), verify `md5sum` against the manifest, then `mv` into place; delete the temp file on any failure and print a clear error.
5. **Fail closed**: after both attempts, if any of the 4 files is still missing → `echo` error + `exit 1` (consistent with the script's existing fail-closed tail).

Reference comments: cite `https://doi.org/10.6084/m9.figshare.26310994` (light reference atlases, Garnica et al. 2024; carmonalab/Reference_maps), matching README provenance.

### 4. Docs
- `docs/ARCHITECTURE.md:66` (`NAS_REF_DIR` row): note Figshare fallback.
- `docs/ARCHITECTURE.md:143` (`1_prepare_chunks.sh` row): "stages ref maps from NAS → scratch" → "stages ref maps NAS → `$HOME` (skips if all 4 present; Figshare download fallback via DOI 10.6084/m9.figshare.26310994)".
- `README.md:171-178` ("Reference data"): add one bullet — the pipeline auto-stages the 4 sketched maps into `HOME_REF_DIR`, NAS first, else downloaded from the Figshare article (link already present); the local repo no longer carries a mirror.

## Validation (implementer)

- `rm -rf` done; `git status` shows no trace of `ProjecTILs_references/`; `git ls-files aux/` unchanged (3 files).
- `bash -n src/4_cell_type_annotation/1_prepare_chunks.sh` passes.
- Optional URL sanity check on the local Mac (do not run pipeline scripts per AGENTS.md): `curl -f -L -o /tmp/test.rds https://ndownloader.figshare.com/files/47714158 && md5 /tmp/test.rds` → `be86058ddafdd0154faf0485286b86e7`, then delete the temp file. (Use `md5sum` in the script itself — Linux HPC; macOS uses `md5`.)
- Full HPC validation of the updated script happens later with the `_debug` dataset, per repo convention.

## Out of scope / notes

- No changes to `2.1.1_process_chunk.R`, `config_helper.R`, `slurm_config.sh`, or `datasets.json`: workers keep reading `HOME_REF_DIR` via `path_ref`; staging remains a login-node concern in `1_prepare_chunks.sh`.
- No `aux/` move (previous plan superseded): `aux/` stays as-is, 100% committed (1.4 MB).
- The NAS copy (`NAS_REF_DIR`) is the in-house source and is left untouched; the deleted repo folder was only a local mirror (verified byte-identical to Figshare).
- No curl/wget precedent in `src/*.sh` — plain system `curl` on the login node is sufficient; no new pixi dependency.
