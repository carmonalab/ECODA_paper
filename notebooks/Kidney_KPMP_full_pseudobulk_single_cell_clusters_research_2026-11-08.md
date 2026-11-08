# Kidney KPMP: two single-cell pseudobulk clusters — 2026-11-08

## Conclusion

The corrected pseudobulk single-cell specimens split into **28 negative-MDS1 specimens and 17 positive-MDS1 specimens**. The strongest observed distinction is an **upstream gene-availability/expression-coverage pattern**, already present in the downloaded source, not simply sex, disease, cell count, nominal assay or anatomy. Several leading MDS1 genes are zero in all diagnostically sampled cells from every left-group specimen but detected throughout the right group. This is consistent with a technical processing/reference/filtering difference; the exact upstream step is **not established**.

The split is **not first created by correction**. The same 28/17 groups also separate when MDS is fitted only to the 45 single-cell specimens using the existing uncorrected distances. In the original all-84 uncorrected figure, cell-versus-nucleus separation dominates MDS1 (86.5%), obscuring this within-single-cell structure.

Investigation only. No datasets/configuration, production matrices, benchmark results, existing plots or notebook implementation were changed. No normalization/correction was rerun. Only small sample metadata, distance/gene-summary artifacts and a bounded selected-gene source diagnostic were transferred; no full H5AD or expression matrix was copied locally.

## Exact result and cluster identity

- Original paired figures: `plots/batch_effect_corrected/Kidney_KPMP_full_mds_region.pdf` and `Kidney_KPMP_full_mds_batch_suspension_type.pdf`, Pseudobulk row.
- Corrected source inspected read-only on Bamboo: `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/batch_effect/corrected_final/results/Kidney_KPMP_full_batch_effect_corrected_final_pseudobulk.rds`, bundle `Pseudobulk_hvg2000`, MD5 **7c8eb867cfa27a54b3592da0da502f9c**. Actual MD5 matches the sidecar and plotting snapshot `data/batch_effect/corrected/20261105_corrected_heatmap/heatmap_scores.csv`.
- Stored features are 84 specimens × 2,000 genes. Their Euclidean distances reproduce the persisted distance matrix to tolerance **1e-8**. `cmdscale(..., k=2, eig=TRUE)` reproduces the plotted MDS percentages **21.91336% / 18.48031%**.
- Corrected sample Feather MD5 **48b33f4dd65055faf85de5fe804d2cdb**, validated against its source-path/MD5/size sidecar. All specimen IDs, ordering, biological labels, sex and source counts match the authoritative source metadata and previous validated composition audit.
- All 45 single-cell specimens have `tissue = kidney`, `sample_tissue_type = Biopsy`, `tissue_type = tissue`, `region.l1 = NA`, and nominal `assay = 10x 3' v3`. Missing region is not an anatomical region or proof of identical microscopic tissue content.
- Descriptive cluster definition: among `suspension_type == cell`, negative versus positive corrected MDS1. The visible gap is large: left values **−35.17 to −9.65**, right **+24.86 to +46.76**. This is a plot description, not a specimen exclusion rule.

### Membership and straightforward metadata

| Property | Left group | Right group |
|---|---:|---:|
| Specimens | 28 | 17 |
| AKI / CKD / Ref | 8 / 9 / 11 | 3 / 5 / 9 |
| Female / male | 15 / 13 | 11 / 6 |
| Median cells | 1,331.5 | 1,540 |
| Cell-count range | 519–4,690 | 517–9,027 |
| Median of per-specimen median detected genes/cell | 1,481.5 | 1,589 |
| Median of per-specimen median source `percent.mt` | 23.99% | 25.70% |
| All reference specimens | 11 `Sample115*_EO*` specimens | 9 `PRE*` specimens |
| Explicit historical-hg19 / CellXGene-rerun note in source supplement | 0 | 13 |

There is **no clean sex, AKI/CKD/Ref, low-cell, mitochondrial-fraction, acquisition-type or nominal chemistry split**. Ages, diabetes history and hypertension vary and overlap across groups. Both groups contain biopsy specimens with substantial cell counts; neither is simply the very low-count end of the cohort.

Reference naming families segregate perfectly:

- Left: `Sample1153_EO1/EO2/EO3`, `Sample1157_EO1/EO2/EO3`, `Sample1158_EO1/EO2/EO3`, `Sample1162_EO1/EO2`.
- Right: `PRE018_1`, `PRE019_4`, `PRE027`, `PRE038`, `PRE055_1`, `PRE062_1`, `PRE19_025`, `PRE19_05`, `PRE98sc`.

Right-group disease specimens are `S_1908_010118`, `S_1908_009836`, `S_1903_005471`, `S_1905_017555`, `S_1909_007154`, `S_2001_005409`, `S_1910_000048`, `S_1908_000898`. The remaining 17 single-cell disease specimens belong to the left group. Complete membership and available metadata are saved in the linked diagnostics CSV.

`experiment` is not a shared binary batch label: all 45 retained single-cell specimens have distinct PREMIERE experiment values in the sample Feather. Source metadata shows that `S-1905-017555` combines three libraries/experiments (PREMIERE15–17), so a sample Feather's single representative library does not capture all source libraries. The diagnostics CSV adds source-level unique-library counts while retaining the Feather's representative library values.

## Cell composition is not the clean dividing variable

Full cell-type counts from the previous authoritative `obs` audit were joined by exact original/plotted sample identity. Representative-cell `subclass.*` values in the per-sample Feather were **not** treated as sample composition.

| Broad author lineage, median specimen percentage | Left | Right |
|---|---:|---:|
| Proximal tubule | 22.98% | 18.41% |
| Thick ascending limb | 17.31% | 14.64% |
| Immune cells | 9.47% | 7.74% |
| Descending thin limb | 0% | 0% |

Ranges overlap substantially; both sides contain variable and sometimes unusual mixtures. This does not rule out a contribution from composition, but it does not supply a coherent binary anatomy/lineage explanation comparable to the gene-coverage pattern below.

## Leading expression drivers and direct source spot check

For Euclidean distances in the stored 2,000-gene feature space, classical MDS is equivalent to PCA up to axis sign. Gene coefficients were obtained by projecting centered stored features onto the exact plotted MDS1 coordinates. The largest coefficients are **MT-RNR2, MT-RNR1, SNHG29, MTATP6P1, GPX1, DANCR, CHASERR, SNHG6, OIP5-AS1 and GAS5**. MT-RNR2 and MT-RNR1 together account for **9.64% of the sum of squared MDS1 feature coefficients**, not 9.64% of total biological variance.

Corrected transformed-feature means differ markedly (these are VST/corrected units, **not log-fold changes**): MT-RNR2 left **1.91** versus right **16.97**; MT-RNR1 **2.28** versus **15.21**; SNHG29 **2.95** versus **12.04**; GPX1 **3.12** versus **11.23**.

To check whether the pattern predates ECODA correction, the Yggdrasil source H5AD `/home/users/h/halterc/scratch/ECODA_paper/Kidney_KPMP_full/data/Kidney_KPMP_full.h5ad` was opened read-only. Exactly **16 evenly spaced source observations per retained single-cell specimen** were read for eight named genes: **448 left-group and 272 right-group cells**. HDF5 CSR row slices were used; no complete expression matrix or whole-cohort aggregation was loaded or transferred. Both source `raw/X` and `X` were checked.

| Gene | Source `raw/X` detections, left / 448 | Source `raw/X` detections, right / 272 |
|---|---:|---:|
| MT-RNR1 | 0 | 271 |
| MT-RNR2 | 0 | 272 |
| SNHG29 | 0 | 212 |
| MTATP6P1 | 0 | 226 |
| GPX1 | 0 | 178 |
| DANCR | 0 | 120 |
| MALAT1, control | 387 | 243 |
| ACTB, control | 411 | 240 |

**Every one of the 28 left-group specimens** has zero sampled detections for the first six genes; **every one of the 17 right-group specimens** has detections for those genes. MT-RNR2 is detected in all 16 sampled cells of every right-group specimen. The left zeros are also present in source `X`, while ACTB/MALAT1 are detected on both sides. Thus this is not a general absence of usable source expression or simply lower cell counts.

These are **bounded diagnostic slices**, not full-specimen gene-detection frequencies or a proof that every cell of a left-group specimen lacks those genes. Nevertheless, coherent sample-wide group concordance across unrelated genes, together with the full-pseudobulk gene loadings, is strong evidence of upstream feature coverage/processing heterogeneity. **[INFERENCE]** Reference/feature filtering, gene-ID handling or upstream count processing is more plausible than a biological on/off program shared by all diseases and both sexes. The exact responsible step remains unverified.

Source `percent.mt` values overlap; absence of mitochondrial rRNA features should not be equated with a generally low mitochondrial-read fraction or blanket failure of the mitochondrial QC metric.

## A documented processing-history difference, with an important caveat

[Lake Supplementary Table 1, column AF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05769-3/MediaObjects/41586_2023_5769_MOESM6_ESM.zip) states for PREMIERE8, 15–17, 27, 28 and 30–38:

> Cell ranger was originally analyzed using Hg19 genome database. But for the CellXGene instance, the sample was re-run using GRch38.3.0

This covers all nine PRE references and four distinct disease specimens: `S-1908-009836`, `S-1905-017555` (three experiments), `S-1910-000048`, `S-1908-000898`. All 13 corresponding retained specimens lie in the right group. None of the 28 left-group specimens carries this note.

However, four right-group disease specimens are **not** flagged: `S_1908_010118` (PREMIERE6), `S_1903_005471` (PREMIERE12), `S_1909_007154` (PREMIERE20), `S_2001_005409` (PREMIERE25). The note is therefore **not a complete explanation** of the 28/17 division, and an absent note does not establish identical processing.

**Do not claim current mixed hg19/hg38 matrices.** The note explicitly says the CellXGene version was rerun with GRCh38. GEO records for both reference families currently state GRCh38 and otherwise report the same nominal 10x v3/NovaSeq/SoupX pipeline: [PRE exemplar PREMIERE30](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5554497&targ=self&form=text&view=full), [EO exemplar PREMIERE39](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5554506&targ=self&form=text&view=full).

The same supplement identifies both PRE and EO reference families as **CZI-HCA living-donor preimplant needle biopsies**, originally [GSE169285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169285). Original GEO records share the reported Michigan acquisition/preparation, Liberase TL dissociation and NovaSeq workflow. No inspected primary source establishes different PRE/EO centers, v2-versus-v3 chemistry, sequencing lanes or an EO suffix meaning. Naming-family concordance is provenance evidence, not proof of distinct acquisition methods or independent donors.

## Why correction does not remove it; why it becomes more visible

The actual plotted corrected result's batch contract records:

- Effective keys **`suspension_type` and `sex`**.
- Additive technical design **`~1 + batch_key_1 + batch_key_2`**, rank 3 and residual degrees of freedom 81.
- Intercept-only preservation; no condition labels, sample-series factor, experiment/run factor or gene-availability group in the correction design.

This model subtracts average suspension/sex effects. It cannot reconstruct gene values absent from one upstream source family or remove an unmodeled within-single-cell processing group merely because both groups are called `cell`.

The traced current producer path is raw-count aggregation, full-gene DESeq2 transformation, additive limma correction, **post-correction** variance ranking, selection of 2,000 genes, then Euclidean distances. There is no post-correction clipping or zero floor in the traced path. Relevant code: `src/5_run_benchmark_methods/benchmark_pipeline.R:805-915`, `src/5_run_benchmark_methods/benchmark_hpc_utils.R:815-960`, `src/utils/batch_contract.R:279-335`, `src/5_run_benchmark_methods/benchmark_methods_r.R:332-343`.

Corrected and legacy uncorrected selected-gene sets overlap in only **592 of 2,000 genes**. Therefore the paired figure is not a fixed-feature experiment isolating limma subtraction alone: the selected feature set also changes, and the uncorrected artifact is an explicitly retained historical producer.

Critically, using the **existing uncorrected single-cell distance submatrix**, the same groups already occupy disjoint within-single-cell MDS1 ranges (left −27.29 to −5.91; right +12.29 to +34.00; horizontal sign chosen for consistency with the original groups). MT-RNR1/2, SNHG29 and GPX1 are present in the uncorrected selected feature list too. The independent 45-sample diagnostic has MDS1 variance **30.5% uncorrected / 36.0% corrected**; these are not the original 84-sample axis percentages. See the saved before/after diagnostic.

**[INFERENCE]** Removing the dominant suspension contrast and reselecting genes exposes/reweights an already existing within-single-cell technical/coverage axis. It is not evidence that correction newly produced the underlying source zero pattern.

## Practical interpretation and unresolved provenance

Treat this as an **unmodeled upstream processing/feature-coverage concern** when interpreting pseudobulk biological separation. Do not automatically call either group biologically valid/invalid, exclude one by MDS side, add a post hoc correction factor or launch a rerun from these findings alone.

To establish the exact cause, the missing evidence is a source-by-source **reference feature manifest / alignment and gene-filtering/count-processing history** for all PREMIERE experiments, especially unflagged right-group PREMIERE6/12/20/25. Existing metadata, source slices and final results identify the pattern, but do not prove which upstream operation created it. The documented hg19/rerun history is a relevant clue, not a complete causal diagnosis. No count-level full-cohort reaggregation or raw-read reanalysis was performed.

## Saved evidence and verification

Run-owned directory: [`data/derived/kidney_pseudobulk_sc_clusters_20261108/`](../data/derived/kidney_pseudobulk_sc_clusters_20261108/).

- `sample_diagnostics.csv`: all 84 original plotted coordinates, single-cell group membership, source sample/clinical/QC metadata, reference family, historical-processing-note indicator, full-author broad compositions and uncorrected coordinates.
- `cluster_metadata_counts.csv`, `cluster_numeric_summary.csv`: categorical distributions and median/range summaries.
- `source_gene_spot_check.csv`: selected-gene detection counts and summed values from the 16-cell-per-specimen source slices, retaining matrix identity (`raw/X` versus `X`). Only retained single-cell specimens are included.
- `MDS1_gene_loadings.csv`: all 2,000 gene coefficients and group-mean summaries; no sample-by-gene matrix.
- `corrected_pseudobulk_clusters.png`: exact original 84-sample corrected coordinates, groups highlighted and sex shapes.
- `single_cell_only_uncorrected_vs_corrected.png`: independently fitted MDS from each existing 45-specimen distance submatrix, visually inspected; no expression/correction rerun.
- `manifest.json`: exact source checksums, correction contract, scope, diagnostic limitations and output MD5 values.

Observed checks: exact persisted distance/feature agreement; correct plotted eigenvalue percentages; metadata/sample identity and count agreement; 28/17 membership; source-zero/detection assertions for the bounded diagnostic; 592-gene overlap; disjoint groups in uncorrected single-cell MDS; both diagnostic images rendered and inspected. No permanent tests or production changes were needed for this investigation.
