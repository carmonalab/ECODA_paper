# Plan: Replace bionty Gene.standardize with flat-file Ensembl 105 reference

## Context

`src/py/preprocess.py` calls `bt.Gene.standardize(adata.var_names, organism="human")` at line 144. This is **broken** — bionty v2.4.2 requires a configured laminDB instance, which was never set up.

**Decision**: Replace bionty with a flat-file approach using the existing `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` reference (same Ensembl 105 reference that STACAS uses). This is:
- Maximally reproducible (file in repo, pinned forever)
- Zero new dependencies (pandas only)
- 100% compatible with the existing STACAS pipeline
- Simple (~10 lines of Python)

## Implementation Tasks

### Task 1: Replace `bt.Gene.standardize()` in `src/py/preprocess.py`

**Remove** `import bionty as bt`
**Remove** the `adata.var_names = bt.Gene.standardize(...)` call in `base_preprocessing()`

**Add** a new function `standardize_gene_symbols(adata)` in `base_preprocessing()` that:

```
1. Load aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz with pd.read_csv(sep="\t")
2. Build two mapping dicts:
   - Gene name -> itself (identity, to standardize already-correct names)
   - Gene Synonym -> Gene name (for aliases)
   Merge them, keeping first occurrence per key (pd.Series.duplicated(keep="first"))
3. Apply: adata.var_names = [mapped.get(g, g) for g in adata.var_names]
4. adata.var_names_make_unique()  (already exists after the old call)
```

Edge cases handled:
- NaN/empty synonyms are skipped (filter with `.dropna()` and non-empty check)
- A synonym mapping to multiple gene names: `keep="first"` resolves this
- Genes not in the reference: keep original name unchanged

File path resolution: same `project_root` pattern already used in `main()`, or use `Path(__file__).resolve().parents[2] / "aux" / ...`

### Task 2: Clean up bionty dependency

If bionty is not used elsewhere in the Python codebase (confirmed: only in `preprocess.py`), remove `bionty` from `[dependencies]` in `pixi.toml`.

Check for any other Python files importing bionty first.

### Task 3: Update `TODO.md`

Replace the investigation TODO (lines 1-2) with a summary of the decision and implementation:
- bionty was replaced with a flat-file Ensembl 105 reference
- `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` is the single source of truth
- No compatibility issues with scATOMIC/HiTME since the same STACAS reference is used

## Validation

1. **Unit test**: Run a small Python snippet that loads the Ensembl 105 file, builds the alias map, and standardizes a test gene list including known aliases, symbols not in the reference, and edge cases.
2. **Run preprocess.py**: Execute with a small dataset to ensure no `CurrentInstanceNotConfigured` or other errors.
3. **Verify reproducibility**: The same gene name always produces the same output regardless of internet connection, date, or environment.

## Files Affected

| File | Change |
|------|--------|
| `src/py/preprocess.py` | Remove bionty import + call; add `standardize_gene_symbols()` |
| `pixi.toml` | Remove `bionty` dependency (if not used elsewhere) |
| `TODO.md` | Replace lines 1-2 with implementation summary |
