# Document provenance of the gene reference file (EnsemblGenes105)

## Context / goal

The download script for the gene reference file dropped out during earlier cleanup, and it should **not** be re-added. Instead, its provenance must be documented so it is clear where the file originally comes from.

Facts confirmed in the codebase:

- The file **still exists and is actively used**: `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` (committed to the repo) is loaded by `src/gene_utils.py:12` and used for gene name standardization in the preprocessing pipeline (`1.1.1_preprocess.py`).
- The **download logic is gone**: `GENE_REF_FILE`/`GENE_REF_URL` were removed from `slurm_config.sh` and the `GENE_REF_FILE` staging block from `1_prepare_chunks.sh` (recorded in `TODO.md:15,412`, `docs/ARCHITECTURE.md:137`). Nothing downloads the file anymore; it must already be present in `aux/`.
- **Original source**: Ensembl 105 human gene reference (GRCh38.p13), retrieved 14.02.2022 from the `aux/` folder of the carmonalab repo:
  `https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`
  (same source referenced in all QC-filtering notebooks, e.g. `notebooks/QC_filtering/scRNAseq_data_processing_SmillieC_2019_31348891.Rmd:267`)
- Path discrepancy note (historical, for accuracy): the old `GENE_REF_FILE` pointed to `${PROJECT_ROOT}/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` (project root), while the actual file lives at `aux/`. Document the canonical path as `aux/`.

## Changes (docs only — no code, no download script)

### 1. `README.md`

Add a short "Reference data" note under `### Repository Contents` (after the `src` bullet list, before `### Installation`), e.g.:

```
### Reference data

- `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`: Ensembl 105 human gene reference
  (GRCh38.p13) used for gene-name standardization in the preprocessing pipeline
  (`src/gene_utils.py`). Originally retrieved on 14.02.2022 from the `aux/` folder
  of the carmonalab `scRNAseq_data_processing` repository:
  https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz
  The file is committed to this repository — there is no download step.
```

### 2. `docs/ARCHITECTURE.md:137`

Extend the existing NAS ↔ Scratch data flow bullet instead of stating only that the staging block was removed. Replace:

```
`1_prepare_chunks.sh` only stages reference maps (gene standardization moved into `1.1.1_preprocess.py`; the `GENE_REF_FILE` staging block was removed).
```

with a version that also names the file's current location and origin, e.g.:

```
`1_prepare_chunks.sh` only stages reference maps (gene standardization moved into `1.1.1_preprocess.py`; the `GENE_REF_FILE` staging block was removed). The gene reference is committed to the repo at `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` (Ensembl 105, GRCh38.p13), originally downloaded 14.02.2022 from https://raw.githubusercontent.com/carmonalab/scRNAseq_data_processing/master/aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz; it is consumed by `src/gene_utils.py`.
```

## Validation

- No pipeline scripts are run (per AGENTS.md rules); this is a docs-only change.
- Verify the rendered markdown reads correctly (no broken links/tables) and that the documented path `aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz` matches the path in `src/gene_utils.py:12`.
- Confirm the file exists in `aux/` (already verified: present, committed Oct 2023).
