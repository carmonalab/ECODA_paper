# AGENTS.md: note on code search (prevent grep timeouts)

## Context

An agent ran `grep -rn "pattern" --exclude-dir=.git .` and timed out (>120 s). The repo root contains gitignored heavy directories that plain `grep -r .` does not skip:
- `data/` — 97 GB, ignored in `.gitignore`
- `.pixi/` — ignored in `.gitignore`
- also ignored: `plots/`, `ProjectTILs_references/`

`rg` (ripgrep) is NOT installed on this Mac (`command not found`), so recommending `rg` is not safe.

## Decision

Add a short "Code search" bullet under `# General rules` in `AGENTS.md` (after the existing two bullets, before the HPC working-dir bullet). This is a documentation-only change.

## Proposed text

```
- Search code with the built-in Grep/semantic-search tools or `git grep` (tracked files
  only). Never run raw `grep -rn "..." .` — it scans the gitignored 97 GB `data/` and
  `.pixi/` and will time out. If plain `grep` must be used, scope the path and exclude
  heavy dirs, e.g. `grep -rn "pattern" --exclude-dir={data,.pixi} src notebooks docs`.
```

Note: `--exclude-dir` only prunes directories found *during traversal* — the positional
roots (`src`, `notebooks`, `docs`) are always searched. So the example above never skips
the source folders; it only skips any `data`/`.pixi` subdirectories along the way.

## Validation

- Re-read AGENTS.md section for formatting consistency.
- Optional smoke test of the suggested command on a real pattern (e.g. `low_res_ct_col`) — plain grep is safe when scoped.
- No pipeline scripts run (per General rules).

## Out of scope

- Installing ripgrep or other tooling.
- Changing `.gitignore`.
