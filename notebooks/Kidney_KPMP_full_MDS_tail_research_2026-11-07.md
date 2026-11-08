# Kidney KPMP ECODA-author MDS left tail — 2026-11-07

## Finding

The extreme left tail primarily separates **papilla/medulla from cortex**, not low-cell or incompletely annotated specimens. The tissue-region metadata, coherent medullary cell-type signatures, and same-donor cortex/papilla comparisons independently support this interpretation. There is no evidence in this audit that the tail was created by ECODA dropping cells or author annotations. This does not establish the absence of every upstream acquisition or processing bias.

Investigation only: no source/configuration changes, benchmark recomputation, scheduler submission, specimen exclusion, or modification of existing plots/results. Only bounded per-specimen metadata/count summaries were transferred from H5ADs; expression matrices were never read or transferred.

## Exact plotted inputs and checks

- Original figure: [`plots/batch_effect_corrected/Kidney_KPMP_full_mds.pdf`](../plots/batch_effect_corrected/Kidney_KPMP_full_mds.pdf). Both ECODA-author panels contain 84 specimens and 76 author `subclass.l3` features.
- Uncorrected input is selected as a legacy lane by `data/batch_effect/uncorrected_final/final_analysis_artifacts.tsv`: `data/batch_effect/uncorrected/results/Kidney_KPMP_full_batch_effect_uncorrected_composition.rds`, bundle `ECODA_authors_HR`, MD5 **a8c8d4a9e7e44cdcb81cf03a3737dccb**. Actual checksum matches its sidecar.
- Corrected input: read on the existing Bamboo source, `/srv/beegfs/scratch/users/h/halterc/ECODA_paper/batch_effect/corrected_final/results/Kidney_KPMP_full_batch_effect_corrected_final_composition.rds`, same bundle key, MD5 **bf364cf53e152ddc394888d08e15aa6d**. Actual checksum matches its sidecar and the local `20261105_corrected_heatmap/heatmap_scores.csv` plotting-provenance snapshot. No RDS was transferred.
- Reconstructed `stats::cmdscale(dist_mat, k=2, eig=TRUE)` coordinates reproduce the PDF's variance percentages: uncorrected **32.03937%, 16.34455%**; corrected **34.74680%, 12.22195%**. The diagnostic plot preserves these coordinates without alignment or refitting.
- Read only selected `obs` fields from the raw, uncorrected and corrected H5ADs under Yggdrasil `/home/users/h/halterc/scratch/ECODA_paper/Kidney_KPMP_full/{data,output}`. Raw source has **304,652 observations / 88 specimens**; both processed views have **303,360 / 84**.
- Resolved the historical specimen-ID normalization explicitly. Numeric-leading specimens retain hyphens after a `g` prefix in the persisted result; most other specimen IDs replace hyphens with underscores. The evidence CSV retains both original and plotted IDs. No fuzzy matching was used.
- **All 6,384 specimen × HR-type counts** in the uncorrected result, after subtracting its 0.5 pseudocount, match the authoritative uncorrected H5AD metadata exactly. Corrected and uncorrected H5AD HR counts also match exactly. All sample IDs, sample cell counts and biological labels agree.
- `log(counts + 0.5) - rowMeans(log(counts + 0.5))` reproduces the uncorrected feature matrix and its Euclidean distances to tolerance **1e-10**.

## Which specimens?

For a reproducible descriptive summary, define the extreme tail as **uncorrected MDS1 < -10**, yielding 11 specimens. This is not a discovered biological boundary, formal outlier test, or exclusion criterion; the broader left shoulder is a continuum.

| Rank along uncorrected MDS1 | Plotted specimen | Condition | Recorded region | Retained cells/nuclei | Uncorrected MDS1 | Corrected MDS1 |
|---:|---|---|---|---:|---:|---:|
| 1 | K2000094_2_R1 | Ref | Papilla | 6,488 | -19.04 | -16.80 |
| 2 | K1900387_4_R1 | Ref | Papilla | 5,815 | -18.29 | -17.58 |
| 3 | K1900175_1_R1 | Ref | Medulla | 3,730 | -16.84 | -14.96 |
| 4 | S_1908_000858_R1 | CKD | Medulla | 8,216 | -16.42 | -14.35 |
| 5 | KRP460P_1 | Ref | Papilla | 6,089 | -15.52 | -13.72 |
| 6 | K1800364_1 | Ref | Medulla | 2,256 | -13.81 | -10.70 |
| 7 | KRP462P_1 | Ref | Papilla | 1,678 | -13.63 | -11.51 |
| 8 | K1800430_3_R3 | Ref | Medulla | 1,291 | -12.26 | -8.55 |
| 9 | K1900387_1_R1 | Ref | Medulla | 5,652 | -11.49 | -8.47 |
| 10 | K1900019_1_R3 | Ref | Medulla | 2,468 | -11.31 | -7.16 |
| 11 | S_1910_000048 | CKD | Unrecorded | 1,724 | -11.03 | -14.09 |

The ten most negative uncorrected specimens are **all anatomically recorded papilla or medulla**, all nuclei; none are cortex. All four retained papilla specimens are among the seven most negative. The 11-specimen set contains 9 Ref and 2 CKD specimens.

The left shoulder includes additional medullary/mixed-region CKD specimens: `S_1908_000952_R1` (6,039 cells/nuclei; MDS1 -8.40), `S_1910_000149_R1` (13,763; -7.82), `S_1908_000905_R1` (5,113; -5.78), and `S_2001_000055_R1` (2,098; -5.69). Some have `region.l2 = C-M`, rather than pure medulla `M`.

Using the separate descriptive cutoff **corrected MDS1 < -10** gives 9 specimens: ranks **1, 2, 3, 4, 5, 6, 7, 11**, plus rank **14**, `S_1904_007298` (AKI, 2,617 cells, unrecorded region; corrected MDS1 -10.12). Axes are lane-specific; the same numeric cutoff is not an identical biological selection.

## Low counts or missing annotations?

| Metric | Extreme uncorrected tail (11) | Remaining specimens (73) |
|---|---:|---:|
| Median retained cells/nuclei | 3,730 | 2,891 |
| Range of retained cells/nuclei | 1,291–8,216 | 517–16,670 |
| Median detected HR author types | 56 | 52 |
| Nonmissing author HR annotation coverage | 100% | 100% |

The tail is not the low-count end of this cohort and does not have reduced author-annotation coverage. Across all 84 specimens, the Spearman correlation of uncorrected MDS1 with cell count is **-0.141**, inconsistent with a simple low-cell left-tail explanation.

For every retained specimen, the number of observations in the raw source equals the processed count exactly. The difference from 88 to 84 specimens is wholly explained by four source specimens below the pipeline's 500-cell specimen threshold: `S-1905-017562-R1` (352), `446_B1` (268), `S-1908-009883` (460), and `S-1905-018690` (212). Their counts sum to 1,292, exactly the source-to-processed cell difference. No additional per-cell loss occurred for the retained specimens in these files.

The smallest retained specimens include `PRE055_1` (517), `Sample1158_EO1` (519), and `PRE062_1`/`PRE19_05` (635 each); none belongs to the extreme left tail. A small count still warrants uncertainty about rare types, but it does not explain this particular geometry.

## How composition differs

Percentages below use actual observed counts, **without pseudocounts**, and are medians of specimen percentages, not pooled cell percentages. Broad `subclass.l1` PT includes altered PT states; it is not restricted to canonical PT-S1/2 and PT-S3.

| Cell type/group | Extreme tail median (%) | Remaining median (%) |
|---|---:|---:|
| All proximal-tubule states, broad PT | 2.67 | 27.21 |
| Descending thin limbs, broad DTL | 11.67 | 0.069 |
| Distal convoluted tubule, broad DCT | 0.044 | 5.48 |
| Principal cells, broad PC | 13.30 | 4.30 |
| Medullary thick ascending limb, M-TAL | 10.70 | 0.104 |
| Medullary fibroblasts, M-FIB | 7.88 | 0 |
| Ascending-vasa-recta endothelium, EC-AVR | 5.59 | 0.321 |
| Descending-vasa-recta endothelium, EC-DVR | 5.19 | 0.210 |

These are coherent anatomical shifts, not merely absent rare classes. Papilla and medulla are enriched for loops of Henle, medullary collecting ducts, medullary stroma and vasa recta; cortical proximal/distal-tubule and glomerular populations are depleted. The source atlas independently maps these types along the corticomedullary axis in [Figure 2](https://www.nature.com/articles/s41586-023-05769-3/figures/2).

Concrete examples:

- **K2000094_2_R1**, papilla: M-FIB **19.9%**, IMCD **18.0%**, ATL **8.3%**; canonical PT-S1/2 + PT-S3 only **0.062%**. IMCD = inner-medullary collecting duct; ATL = ascending thin limb.
- **K1900387_4_R1**, papilla: M-FIB **25.9%**, IMCD **13.7%**, ATL **10.4%**, DTL3 **8.2%**, degenerative IMCD **8.2%**.
- **KRP460P_1**, papilla: IMCD **31.1%**, degenerative IMCD **16.1%**, M-FIB **9.1%**.
- **S_1908_000858_R1**, CKD medulla: M-TAL **28.6%**, DTL1 **8.2%**, DTL2 **6.6%**, OMCD-PC **6.5%**, degenerative OMCD-PC **6.1%**. OMCD-PC = outer-medullary collecting-duct principal cells.
- **S_1910_000048**, CKD with unrecorded region: ATL **23.6%**, degenerative OMCD-PC **14.8%**, DTL1 **7.1%**, EC-DVR **5.6%**, IMCD **5.2%**. **[INFERENCE]** This resembles medullary sampling, but the missing region cannot be promoted to a known region.
- **S_1904_007298**, AKI with unrecorded region, farther left after correction: DTL1 **26.9%**, degenerative OMCD-PC **21.0%**, C-TAL **19.3%**. **[INFERENCE]** Also anatomically suggestive, not proof of region.

Among the 39 specimens with recorded anatomy, region as a three-level factor descriptively accounts for **82.5%** of uncorrected MDS1 variance and **83.0%** of corrected MDS1 variance. All 39 are nuclei, so this is not simply the whole-cohort cell-versus-nucleus contrast. These are in-sample descriptive R-squared values, not causal estimates or validated predictive performance. Region is unrecorded for the other 45 specimens.

### Same-donor evidence

| Donor | Cortex specimen / cells | Papilla specimen / cells | Broad PT: cortex → papilla | Uncorrected MDS1: cortex → papilla |
|---|---|---|---|---|
| KRP460 | KRP460C_1 / 4,718 | KRP460P_1 / 6,089 | 34.3% → 1.56% | +5.17 → -15.52 |
| KRP462 | KRP462C_1 / 3,473 | KRP462P_1 / 1,678 | 30.3% → 5.84% | +3.46 → -13.63 |

These pairs hold donor, condition, sex and suspension type constant. The cortex specimens lie in the main cloud while papilla specimens lie far left, including after correction. Acquisitions/libraries still differ, so pairing is strong anatomical evidence, not proof that every technical effect is absent.

## Why the separation survives correction

ECODA computes Euclidean distances between natural-log CLR compositions with 0.5 added to every count. A specimen with many medullary types present and many cortical types genuinely absent changes numerous log ratios together, producing a large compositional distance. Rare-type zeros/pseudocounts can affect distances, but the independently observed broad-lineage shifts show this is not only a rare-type sparsity artifact. MDS1's direction is arbitrary; the contrast and distances, not the word “left,” carry meaning.

**Important provenance distinction:** the actual plotted corrected bundle records technical fixed effects for **`suspension_type` + `sex`** in `batch_contract$ordered_source_keys`. Current `datasets.json` lists only `suspension_type`. This investigation describes the persisted plotted result, not a hypothetical rerun under today's configuration. Neither historical correction key explicitly represents anatomical region. Correction of technical means cannot make a papillary specimen compositionally equivalent to cortex or recreate cortical populations absent from the captured tissue.

The source atlas documents distinct scRNA/snRNA preparation and some modality-specific resolution/capture differences. It also notes that deeper medulla is less likely to be captured in needle biopsies. These mechanisms support anatomical and acquisition heterogeneity; they do not identify processing failures in these specimens. Fine author labels include adaptive/degenerative/cycling states, not only stable lineages. See [Lake et al., 2023](https://doi.org/10.1038/s41586-023-05769-3), especially Methods and Figure 2.

## Interpretation and next analysis

Do not discard these specimens merely because they are composition outliers: the strongest extremes look like valid samples from different kidney compartments. For AKI/CKD/Ref comparisons, **anatomical sampling is a confounder** because the tail is dominated by Ref papilla/medulla. A separately declared region-matched/cortex-only sensitivity analysis is more interpretable than post hoc exclusion by MDS position. Missing-region specimens would need source metadata resolution before claiming a full region-matched cohort. Any new selection must be explicit and separately owned; no selection was changed here.

Limits: no tissue histology, sequencing-level QC, raw read reprocessing, or sample-specific upstream protocol audit was performed. The metadata/count evidence supports anatomy as the dominant observed driver of MDS1; it cannot rule out every upstream capture bias, label error, injury effect, or acquisition difference.

## Saved evidence

All new outputs are in [`data/derived/kidney_mds_tail_20261107/`](../data/derived/kidney_mds_tail_20261107/):

- `sample_diagnostics.csv`: all 84 specimens, original/plotted IDs, coordinates in both views, counts, annotation coverage, anatomical/acquisition metadata, observed types and broad/fine percentages.
- `composition_percent.csv`: 84 × 76 observed fine-type percentages, excluding pseudocounts.
- `broad_composition_percent.csv`: broad author-lineage percentages.
- `celltype_contrasts.csv`: fine-type tail/rest median and mean percentages; descriptive only.
- `mds_by_anatomical_region.png`: visually inspected diagnostic using the original plotted coordinates; rank labels refer to the table/CSV.
- `manifest.json`: source identities, scope, checks and output MD5 values.

Primary source identification: [CELLxGENE collection API](https://api.cellxgene.cziscience.com/curation/v1/collections/bcb61471-2a44-4d00-a0af-ff085512674c) reports a 304,652-cell integrated adult kidney sc/sn dataset linked to Lake et al. This matches the repository's declared source cohort; it does not prove byte identity with the latest downloadable CELLxGENE version.
