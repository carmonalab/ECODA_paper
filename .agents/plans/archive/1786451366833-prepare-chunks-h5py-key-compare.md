# Root fix: memory-lean multi-view key-set comparison in `1.1_prepare_chunks.py`

## Goal

Fix the Stephenson chunk-preparation OOM (`srun` job 4307030, oom_kill on gpu004,
`prepare_chunks_production_Stephenson.log`) at its root, so the `--mem=4G` srun
allocation in `1_prepare_chunks.sh` stays valid for every dataset and future large
multi-view datasets cannot regress.

## Root cause (context)

- `1_prepare_chunks.sh:275-282` runs each dataset's chunk build via
  `srun ... --mem=4G "${PYTHON_BIN}" 1.1_prepare_chunks.py`.
- Stephenson is the only large dataset with **2 views** (`benchmark_analysis` =
  Site=="Ncl" subset, `batch_effect_analysis` = full cohort). All other datasets
  have 1 view → `build_union()` takes the `len(h5ad_files) == 1` streaming h5py
  copy branch (`write_annotation_union`), which never opens the file in anndata.
- The multi-view branch (`1.1_prepare_chunks.py:149-157`) opens each view with
  `ad.read_h5ad(f, backed="r")` to compare `(sample, barcode)` key sets. Per the
  script's own docs (lines 18-27, verified on anndata 0.12.10/0.12.19 — the
  repo pins `==0.12.19`), backed open **eagerly materializes the whole `layers`
  group into RAM**. Preprocessed views keep the full raw-counts matrix in
  `layers["counts"]` (float64 CSR; full-cohort view ≈ 800k cells ≈ 6-10 GB) →
  exceeds the 4 GB cgroup → oom_kill.

## Scope

- **Changed**: `src/4_cell_type_annotation/1.1_prepare_chunks.py` (code + module
  docstring), `docs/ARCHITECTURE.md` (line 148: `1.1_prepare_chunks.py` row,
  "backed obs-only key-set reads" wording).
- **Untouched**: `1_prepare_chunks.sh` (`--mem=4G` stays), `datasets.json`,
  `pixi.toml`, all other scripts.
- **Out of scope**: the in-memory `sc.concat` path for truly partial-overlap views
  (`build_union()` step 2, lines 169-209) — still memory-heavy by design, keeps
  its warning to raise srun `--mem`; no current dataset reaches it (Stephenson
  union == batch-effect view → streaming path after the comparison).

## Implementation steps (apply in order)

### 1. Add module-level helpers in `1.1_prepare_chunks.py` (after `cell_keys()`, ~line 55)

- `_read_str_dataset(ds)` — read a 1-D string h5py dataset into a numpy `str`
  array:
  - fixed-length bytes (`dtype.kind == "S"`): `np.char.decode(arr, "utf-8")`
  - variable-length (h5py `string_dtype`, `dtype.kind == "O"`, object array of
    `bytes`): explicit per-element `x.decode("utf-8")` (do **not** rely on
    `astype(str)` on object arrays of bytes — may produce repr-style
    `"b'...'"`; if the element is not `bytes`, `str(x)` fallback)
  - numeric (kinds `iu f`): `arr.astype(str)` (mirrors the old
    `obs_names.astype(str)` semantics)
  - anything else: `RuntimeError` (fail closed, never silently mis-decode)
- `_read_obs_column(obs_grp, col)` — read one obs column as a str array:
  - missing column: let h5py raise `KeyError` (same semantics as
    `adata.obs[col]`)
  - `encoding-type` attr: read defensively (`str` or `bytes`, decode bytes)
  - `"categorical"`: `categories = _read_str_dataset(ds["categories"])`;
    `codes = ds["codes"][:].astype(np.int64)`; missing codes (`-1`) →
    `"nan"` via `np.where(codes >= 0, categories[np.clip(codes, 0, None)], "nan")`
    (matches the old pandas `astype(str)` NaN semantics; `-1` must never index
    `categories`)
  - plain string dataset (no categorical encoding): `_read_str_dataset(ds)`
  - other encodings: `RuntimeError` (fail closed)
- `read_obs_keys_h5py(path, sample_col)` — one `with h5py.File(path, "r")`:
  1. guard: `obs` group `encoding-type` must be `"dataframe"` (anndata 0.12
     layout; `"recarray"`/other → `RuntimeError` fail closed)
  2. `barcodes = _read_str_dataset(obs["_index"])`
  3. `samples = _read_obs_column(obs, sample_col)`
  4. row-count guard from the counts source (mirror `write_annotation_union`'s
     fallback at lines 103-109 and CSR checks at lines 115-120): `layers/counts`
     if present else `X`; its `encoding-type` must be `"csr_matrix"` else
     `RuntimeError`; read `indptr`; `n_rows = len(indptr) - 1` must equal
     `len(barcodes)` else `RuntimeError` (corrupt obs↔matrix mismatch, fail
     closed). `nnz = int(indptr[-1])` is available if wanted for logging.
  5. return `{f"{s}_{b}" for s, b in zip(samples, barcodes)}`

  Memory: two string columns + key set ≈ 100-200 MB for a full-cohort view
  (800k cells); indptr read is ~6 MB.

### 2. Replace the key-comparison loop in `build_union()` (lines ~149-157)

Old:

```python
    key_sets = []
    for f in h5ad_files:
        adata = ad.read_h5ad(f, backed="r")
        keys = set(cell_keys(adata, sample_col))
        key_sets.append(keys)
        adata.file.close()
        print(f"  {f.name}: {len(keys)} cells")
```

New:

```python
    key_sets = []
    for f in h5ad_files:
        keys = read_obs_keys_h5py(f, sample_col)
        key_sets.append(keys)
        print(f"  {f.name}: {len(keys)} cells")
```

- Update the leading comment (lines 149-151): explain h5py-only reads avoid the
  backed-open eager `layers` materialization (anndata 0.12.x, Stephenson
  full-cohort example) — keep the "if one view equals the union, stream it"
  logic description.
- `cell_keys()` stays (still used by the in-memory concat path, line ~177).
  `np` stays used (`np.argmax`, line ~160). No import changes (anndata/h5py/
  numpy/scipy/scanpy all already imported).

### 3. Update the module docstring (lines 29-35)

Replace "backed obs-only reads" with h5py-only obs reads (no anndata open; note
the eager-layers-backed-open reason and that this keeps chunk prep at ~100-200 MB
for full-cohort multi-view datasets).

### 4. Update `docs/ARCHITECTURE.md` (line 148, `1.1_prepare_chunks.py` row)

Change "(backed obs-only key-set reads; no in-memory concat)" to reflect the
h5py-only key-set read (obs `_index` + sample column with categorical decoding +
counts `indptr` row-count guard; no anndata open, no eager layers materialization
at 4G). Keep the sentence short (table style).

## Validation

1. **Local smoke test (implementer, non-pipeline — permitted)** — throwaway
   script in `/tmp` (not committed) run with the local pixi python env
   (`pixi run -e py-cpu python`, anndata `==0.12.19` pinned in pixi.toml;
   `pixi install -e py-cpu` first if missing):
   - write two tiny h5ads mimicking the preprocessed layout (one with the sample
     column as plain strings, one as categorical; both with barcodes as
     variable-length strings and `layers["counts"]` CSR)
   - assert `read_obs_keys_h5py(f, sample_col) == set(cell_keys(ad.read_h5ad(f), sample_col))`
     for both encodings
   - assert the row-count guard raises when obs has one extra row
   - assert a missing sample column raises `KeyError`
   - if the local env is unavailable, skip — the HPC rerun is authoritative
2. **HPC (user runs, after pulling the commit)**:
   - `./src/4_cell_type_annotation/1_prepare_chunks.sh production Stephenson` →
     expect: success, log shows `Stephenson: <N> cells` / `Union: <N> cells` /
     `Union == view ...; streaming copy` / `Union self-check OK`, srun still
     `--mem=4G` (no `1_prepare_chunks.sh` change)
   - full `./src/4_cell_type_annotation/1_prepare_chunks.sh production` → 0
     failed (single-view datasets unaffected; `_debug` still skipped as
     annotated)
   - then proceed with `2_submit_hpc_array.sh` as normal

## Risks / failure modes

- Unexpected string encoding on disk (fixed vs variable length, `<U`): helper
  raises loudly (fail closed) rather than silently producing wrong keys; if a
  real dataset trips it, inspect the actual `obs` encoding via h5py and extend
  `_read_str_dataset`/`_read_obs_column`.
- `-1` category codes (missing samples) mapped to `"nan"` to match old pandas
  semantics — keys stay deterministic and consistent across views.
- Single-view datasets: streaming path untouched, zero behavior change.
- Do not run HPC pipeline scripts from this repo for validation (AGENTS.md);
  validation is user-run on HPC.

## Commit (per AGENTS.md workflow)

1. Move this plan file to `.kilo/plans/archive/`
2. `git add .` (script + ARCHITECTURE.md + archived plan)
3. Commit (concise message, e.g. "fix: h5py-only key-set comparison in prepare_chunks (Stephenson OOM)")
4. Push; user pulls on the HPC clone and runs the HPC validation above.
