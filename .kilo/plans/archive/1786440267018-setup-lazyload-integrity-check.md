# HPC findings: transient NFS stale-view failures after env mutations (2026-08-11)

## Findings

Three consecutive preprocess arrays each failed exactly ~1 random task at the
rpy2/python import stage, with ENOENT/empty-read errors on files that **exist**:

| Array | Task | Node | Error |
|-------|------|------|-------|
| 4306707 (12:41) | 8 / Lee | cpu014 | ENOENT `abind/R/abind.rdb` (R ns load) |
| 4306765 (12:58) | 14 / Zhu | cpu027 | ENOENT `rlang/Meta/features.rds` (R ns load) |
| 4306830 (~13:05) | 14 / Zhu | cpu027 | `TypeError: 'NoneType' object is not iterable` at `import scanpy` — `importlib.metadata.version("anndata")` returned None because `anndata-0.12.19.dist-info/METADATA` read as empty (reproduced locally: empty METADATA → `version()` → None → `Version(None)` → TypeError) |

Verified:
- All failing files/dirs healthy from the login node; exactly one anndata
  dist-info; anndata==0.12.19 / scanpy==1.12.2 pinned; 13/14 tasks read the same
  env fine in every run (persistent corruption would fail ALL tasks identically).
- cpu027 COMPLETED the same task in run 1 → not a bad node (3 different nodes
  failed across runs; cpu022 clean control).
- Every failure sat within ~1h of an env mutation (~12:35 pixi prefix update via
  `pixi run`; 12:45–48 `install.packages("abind")` + `pixi run setup`).

Model: worker nodes serve **stale NFS client caches** (missing directory entries,
torn/empty file reads) after env mutations on the NFS home filesystem, until the
attribute-cache window passes. The login node always revalidates fresh, so the env
always looks healthy there.

Resolution: array 4306890 (~13:3x, after the cache window) passed all 14 tasks
and synced to NAS. **No code fix needed** — recovery is the idempotent re-run;
prevention is "wait up to ~1h (typically faster for most nodes) after env
mutations before arrays". The `PIXI_CACHE_DIR`
`~/.bashrc` workaround was explicitly rejected (user-local hack, not reproducible).

## AGENTS.md update (implementation task)

Append to the env-refresh paragraph (line 141), after the existing
"observed: ..." corruption list, one concise note:

> HPC (2026-08-11): after env mutations (pixi install/prefix update,
> install.packages, setup), worker nodes can serve stale NFS views — arrays may
> fail ~1 random task at the R/python import stage with ENOENT/empty reads on
> files that exist (files fine on the login node). Wait up to ~1h (typically
> faster for most nodes) after env changes before submitting arrays, or re-run
> (idempotent). Details:
> `.kilo/plans/1786440267018-setup-lazyload-integrity-check.md`.
