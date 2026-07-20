# preprocess.py
- hvgs: for non-batch views, make sure that sc.pp.highly_variable_genes(batch_key="Sample") is used
- for non-batch views
- Add blacklist as default for filtering genes
    - Load gene blacklist from file (e.g. from STACAS, see call to `default_black_list` in get_pb_deseq2 in src/R/pseudobulk.R)
        - Maybe save blacklist file to this repo for clarity and add explanation
    - Filter out blacklisted genes before HVG calculation
- if is_batch_view don't do clustering (not needed in this case)

# Batch effect
Should it be done once without batch correction, and once with? -> probably more important to only do WITH batch correction.

The final analysis for batch effect correction needs to be run on the following methods and batch effect mitigation strategies:
- ECODA: remove cell types that are significantly different across batches
    - Possibly re-use process_coda_fig() with added argument batch_col = NULL?
    - Use metadata column batch_col to test for batch-associated cell types (after clr-transformation), remove them and re-calculate clr-transformed cell type composition
        - Which statistical method to use and which threshold to use for significant difference between batches? -> Print warning naming cell types and Significance. -> Possibly depends on the number of batches:
            - If 2 batches: use t-test or Wilcoxon rank-sum test, p-value < 0.05
            - If >2 batches: use ANOVA or Kruskal-Wallis test, p-value < 0.05
- Pseudobulk (DESeq2 + limma)
    - Possibly re-use `process_pseudobulk_fig()` with added argument `batch_col`?
- MrVI (native batch effect correction)
- GloScope (run on harmony integrated space)
- PILOT-GM-VAE (run on harmony integrated space)


## Preprocessing
- handle hvg calculation using correct batch_key
    - datasets.json: add "columns" "batch"
        - Joanito batch was manually defined at the end of MAIN_analysis.rmd
    - update datasets with batch column mapping
    - handle case where multiple datasets are combined (e.g. in Batch_effect.rmd at "# Combined PBMC (Stephenson, GongSharma, Zhu)")

## Down-stream
- DESeq2.normalize():
    - check batch_col is correctly implemented
    - check that it does not get "Sample" as a batch column
- get_pb_deseq2(): add argument batch_col, implement batch_col

## Analysis
- Centralize analysis for batch effect in one script file, separate from MAIN_Analysis.rmd (currently fragmented, see Batch_effect.rmd and end of MAIN_Analysis.rmd)