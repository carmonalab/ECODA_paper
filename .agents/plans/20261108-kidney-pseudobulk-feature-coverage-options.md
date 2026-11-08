# Kidney KPMP residual pseudobulk split: discussion and options

Status: **brainstorming only; no implementation or compute approval**. No pipeline/configuration edits, new batch labels, filtering, exclusions, normalization/correction or cohort reruns were performed. This note integrates existing diagnostic evidence, additional checks of the saved metadata and primary-source research. It is not a launch plan.

## Recommended direction

1. **Name the uncertainty correctly:** an observed upstream feature-coverage signature, not a proven PRE-versus-EO library effect. Neither group is established as the defective group.
2. **Resolve original library-level provenance and feature availability first**, including the 39 nuclei, then assign evidence-backed groups.
3. If a defensible comparable gene universe exists, **prefer a common-quantifiable-feature sensitivity analysis before stronger correction**. Add a validated processing covariate only for remaining effects on genuinely comparable genes.
4. If references/counting definitions materially differ, prioritize restoring or uniformly requantifying the underlying measurements. If neither repair nor a comparable core is possible, stratify/restrict conclusions rather than forcing integration.
5. Keep all existing results immutable. Any future experiment uses separate run-owned outputs, fixed scope and explicit approval; a full-cohort count audit would require the repository's scoped durable-HPC workflow.

## 1. Which samples are currently implicated?

Existing evidence: [full investigation](../../notebooks/Kidney_KPMP_full_pseudobulk_single_cell_clusters_research_2026-11-08.md), [all-sample inventory](../../data/derived/kidney_pseudobulk_sc_clusters_20261108/sample_diagnostics.csv), [source spot-check summaries](../../data/derived/kidney_pseudobulk_sc_clusters_20261108/source_gene_spot_check.csv).

| Provisional description | Specimens | Reference specimens | Disease specimens | Meaning |
|---|---:|---:|---:|---|
| Source signature A; negative corrected pseudobulk MDS1 | 28 | 11 EO-named | 17 | Six sentinel genes were zero in all 16 sampled source cells per specimen; candidate measurement-support issue, not proven full-specimen absence |
| Source signature B; positive corrected pseudobulk MDS1 | 17 | 9 PRE-named | 8 | All six sentinels detected in each specimen's diagnostic slice; not automatically the more correct processing route |
| Nuclei | 39 | — | — | Not assessed by that single-cell source-gene spot check; do not assign an A/B processing identity by assumption |

Thus **all 45 single-cell specimens participate in the observed heterogeneity**. The 28 signature-A specimens are candidates for missing-feature support, but the 17 signature-B specimens are comparators—not a validated gold standard. The exact affected gene set and processing cause remain unresolved.

### Exact signature-A candidates

Reference specimens:

- `Sample1153_EO1`, `Sample1153_EO2`, `Sample1153_EO3`
- `Sample1157_EO1`, `Sample1157_EO2`, `Sample1157_EO3`
- `Sample1158_EO1`, `Sample1158_EO2`, `Sample1158_EO3`
- `Sample1162_EO1`, `Sample1162_EO2`

Disease specimens:

```
S_1904_007298
S_1904_008134
S_1908_000757
S_1908_000945
S_1908_009646
S_1908_009740
S_1908_009787
S_1908_010024
S_1908_010071
S_1910_000095
S_1910_000142
S_1910_000189
S_2001_000048
S_2001_000095
S_2001_000142
S_2001_019418
S_2002_007748
```

### Exact signature-B comparators

References: `PRE018_1`, `PRE019_4`, `PRE027`, `PRE038`, `PRE055_1`, `PRE062_1`, `PRE19_025`, `PRE19_05`, `PRE98sc`.

Disease: `S_1903_005471`, `S_1905_017555`, `S_1908_000898`, `S_1908_009836`, `S_1908_010118`, `S_1909_007154`, `S_1910_000048`, `S_2001_005409`.

PRE/EO describes the reference naming families only. EO also appears in disease-library names on **both** sides—for example left `1121-EO-2` and right `613-EO-1`. Parsing the substring `EO` is not a valid batch assignment rule.

### What other metadata already tells us

- Both sides contain AKI, CKD and Ref, both sexes, overlapping cell counts, detected-gene counts and mitochondrial percentages. All 45 are recorded as kidney biopsies, nominal 10x 3′ v3, with unrecorded anatomical region.
- Available age categories, diabetes history, hypertension and ethnicity do not give a clean binary explanation. Ethnicity is mostly unknown (21/28 left, 14/17 right). BMI, eGFR, percent cortex and percent medulla have **no usable numeric values for these 45 specimens** in the saved source metadata; their presence as column names does not make them usable explanatory covariates.
- The historical hg19-to-GRCh38-rerun comment applies to **13/17 signature-B and 0/28 signature-A** specimens. Four signature-B disease exceptions are PREMIERE6, 12, 20 and 25. It is an incomplete provenance clue, not an exact assignment rule or evidence of current mixed genome builds.
- Every retained sc specimen has a distinct representative PREMIERE experiment ID. Fitting one coefficient per experiment would saturate the sample-level model rather than discover a useful shared batch.
- The 84 retained specimens represent **89 source libraries**. `S_1905_017555` combines three libraries/experiments; nuclear `K1900174_2` combines four. Audit at library level before collapsing to specimens. A single representative Feather value must not silently decide a mixed-library specimen's correction label.
- Repeated donor units matter: the 17 signature-B specimens have 15 distinct recorded donor IDs; PRE019 appears in multiple specimens. Do not infer independence or replicate structure from filename prefixes/suffixes alone.

## 2. Proposed audit: establish measurement support, not just cluster membership

No new count audit was run for this discussion. Proposed scope is all **84 retained specimens / 89 contributing libraries**, not only the 2,000 selected pseudobulk genes and not a whole-repository scan.

### A. Build a provenance crosswalk

For each original library, record original and canonical specimen ID, donor, experiment, contributing cells, source accession, and the actual acquisition/processing fields available: collection center/date, preservation/dissociation, chemistry/kit, sequencing instrument/run/lane, alignment software/version/arguments, reference FASTA/GTF and hashes, feature manifest, exon/intron policy, ambient correction, gene filtering and gene-ID conversion.

Use source manifests/logs/reference metadata when possible; record unavailable fields as **unknown**, never infer them from GEO submitter institution, PRE/EO names, MDS position or a generic GRCh38 string. All parts of a specimen must agree before assigning one specimen-level provenance category. Otherwise retain ambiguity and resolve at library level; do not take the first row or majority category.

### B. Separate three meanings of zero

For every original library and gene, distinguish:

1. Feature existed and was quantified under a documented definition, with zero assigned counts.
2. Feature was absent/removed before matrices were combined: unavailable measurement, potentially encoded as zero later.
3. Original availability or quantification definition is unknown.

Audit stable gene IDs, original symbols/aliases, duplicate/merged/split IDs, version stripping, biotype filters and reference differences. A common `.var` table or matching gene names after concatenation is insufficient. AnnData explicitly documents that sparse outer joins zero-fill unavailable variables [R1]. A `.raw` slot does not recover the pre-merge feature manifest merely because it is named raw.

### C. Expand the count evidence

A future scoped audit should summarize gene × library and gene × specimen total counts plus detected-cell counts/fractions across the **full gene universe**, with library depth and cell counts. Include the six sentinels, positive controls and a systematic search for additional discontinuous features. Use per-lineage summaries only as evaluation diagnostics to ask whether composition/depth can explain a detection pattern—not to use author/condition labels to fit a gene mask or model.

Check original source matrices against processed `layers['counts']` for identical selected cells/features. Existing source zeros establish that ECODA correction did not create the six observed source patterns, but not which upstream step caused them. Include nuclei to identify shared or additional coverage patterns rather than assuming nuclei are unaffected.

### D. Freeze evidence-backed classifications

After reviewing the audit, freeze:

- Library/specimen provenance membership, including unresolved and mixed cases.
- A gene-availability/comparability mask, with a recorded reason per exclusion.
- Whether a gene is unavailable, ambiguously mapped, comparably measured, or unresolved.
- The analysis's intended sample universe and source revision.

Do not define the mask as “genes expressed in every sample,” “genes associated with condition,” or “genes that make the MDS lobes disappear.” True lineage-specific and low-expression zeros must remain possible. Do not turn provisional expression clusters into established technical metadata.

## 3. Remediation options and tradeoffs

| Option | When justified | Main benefit | Main limitation / risk | Position |
|---|---|---|---|---|
| Restore original counts / harmonize feature identity | Original equivalent measurements exist but were dropped, aliased or zero-filled during merging | Repairs bookkeeping at source; preserves potentially useful genes | Cannot recover truly uncounted molecules; symbol matching alone can merge distinct loci | Preferred if audit finds a reversible matrix-level problem |
| Uniform requantification from reads | References, annotation/counting definitions or mapping differ materially; common rows are not reliably comparable | Most direct route to consistent measurement | Broader authorized compute/data project; may change barcodes/genes; raw reads and full provenance required | Scientifically strongest when measurement definitions differ materially |
| Audited common-quantifiable gene universe | Comparability can be established for an informative core but some genes are structurally unsupported | Removes unsupported features without inventing expression | Loses information; may not remove assignment differences among retained genes | Preferred practical first sensitivity analysis if a defensible core exists |
| Add an independently validated processing covariate | Comparable genes still show supported processing-associated shifts | Fits the missing technical effect using existing additive machinery | Cannot repair structural unavailability; can remove real variation when provenance is confounded | Conditional second step, possibly combined with the mask |
| Stratify by processing provenance / analyze sc and sn separately | Comparable core or correction is uncertain, or as a robustness analysis | Avoids forcing unlike measurements into one corrected space | Smaller groups; changes the question; not a replacement for the full cohort | Required sensitivity option; possible fallback for restricted conclusions |
| Drop specimens | Independently demonstrated integrity/identity failure or no defensible comparison | Can exclude genuinely unusable measurements | Changes cohort and potentially biology/donor representation | Last resort, never based solely on a lobe or PRE/EO prefix |

Deleting only MT-RNR1/2 is an **ablation**, not a root-cause solution: the signature already involves non-MT genes, including GPX1. Blanket removal of all mitochondrial genes, lncRNAs or pseudogenes is also not a validated harmonization strategy. Conversely, stronger integration cannot reconstruct unmeasured expression; additive correction may center a structurally zero feature without recovering molecules [R3].

## 4. Where would gene exclusion belong?

### Pseudobulk-focused experiment — narrowest first

Proposed sequence:

**Preserved raw counts → raw sample aggregation → approved comparable-feature mask → DESeq2 transformation → joint technical correction → label-blind HVG selection → distances/MDS.**

The existing narrow seam is `src/5_run_benchmark_methods/benchmark_pipeline.R:815-860`, before `fit_pseudobulk_for_variants()` fits the full-gene aggregate; selection follows at `:899-915`. The mask is an analysis view, not a destructive modification of the original H5AD/raw counts.

This would test the pseudobulk issue without automatically rerunning canonical cell embeddings, changing author annotations, or rebuilding unrelated methods. A mask applied only after the current final 2,000 genes cannot assess effects on earlier VST/HVG selection; it remains a useful explicitly labeled diagnostic.

Normalization nuance: the current path uses default DESeq2 ratio size factors. Genes with zero geometric means do not supply usable ratios, so it is **not established that the observed group-zero genes currently bias size factors** [R4]. Distinguish genes used for normalization from the VST/HVG/distance feature universes, and evaluate any normalization change rather than assuming it.

### Dataset-wide feature harmonization — broader separate decision

If the audit establishes that other expression-based methods are affected, propose a shared comparable analysis gene universe **before Stage 3 normalization/HVG/PCA**, while retaining original complete source/count artifacts. This is a dataset-wide scientific input change with downstream consequences, not an automatic extension of a pseudobulk repair. Existing author-composition counts do not directly use expression genes; their remaining anatomy-related structure must not be “fixed” by this mask.

Cell-level Harmony correction alone is not a pseudobulk fix: this pseudobulk path reads preserved raw counts, not Harmony embeddings or a corrected cell-expression layer.

## 5. If adding a correction term, avoid circular or unidentifiable designs

The plotted historical result used suspension + sex; today's registry differs. Hold the chosen baseline technical specification fixed across comparisons, state whether it is the plotted historical baseline or a separately approved new target, and do not silently switch configurations.

A documented sc-specific processing contrast may eventually be estimable alongside suspension and sex. Alternatives include a full-rank three-level technical group (sc-A, sc-B, nuclei) **instead of** a separate suspension term, or suspension plus one explicitly within-sc contrast (zero for nuclei). The latter avoids claiming nuclei share sc-A's processing provenance. If the nucleus audit reveals additional groups, reassess rather than forcing three levels.

Temporary rank calculations on the existing 84-sample metadata (algebra only; no correction fitted):

| Illustrative design | Columns | Rank | Residual df |
|---|---:|---:|---:|
| Existing suspension + sex | 3 | 3 | 81 |
| Suspension + provisional three-level group + sex | 5 | 4 | 80 |
| Provisional three-level group + sex | 4 | 4 | 80 |
| Representative experiment ID + sex | 85 | 84 | 0 |

The second is redundant; the last is saturated. Provisional MDS labels were used solely to illustrate nesting/rank, **not adopted as correction labels**. Rank correctness is necessary but does not resolve biological confounding. Fit approved nuisance terms jointly on the same pre-correction transformed input, not as arbitrary sequential regressions of already corrected outputs.

No `condition.l1`, disease/Ref indicator or other ground-truth label may enter normalization, feature selection, correction design or protected design columns. Metadata-versus-condition tables are evaluation only. Standard condition-supervised nuisance estimation is not compatible with this label-blind benchmark contract.

## 6. Fair, staged comparison — proposed, not run

### Phase 1: frozen-transform diagnostic

After independently auditing a mask and any processing covariate, freeze samples, the pre-correction transformed input, gene list H, weights/scaling and distance choices. Compare:

| Features | Baseline technical design | Baseline + validated processing contrast |
|---|---|---|
| H | Reference | Correction-only diagnostic |
| H minus audited unsupported features | Mask-only diagnostic | Mask plus correction diagnostic |

Do not refill the HVG quota in this phase. Also compare arms on the same retained gene coordinates to avoid conflating correction with dimensionality. With fixed expression/design/weights, limma is gene-wise: deleting other genes does **not** alter retained-gene corrections [R3]. A feature-mask effect here is a change in which genes contribute to distances, not evidence that masking magically changed other genes' coefficients.

The currently published paired outputs cannot isolate correction alone: only 592/2,000 selected genes overlap and their historical transforms/producers differ. A controlled future input must be explicitly defined; simply reusing those two final matrices is not a fixed-feature correction experiment.

### Phase 2: intended end-to-end sensitivity

Only after Phase 1, separately authorize rebuilding the justified count-mask/normalization/HVG path. Keep sample unit/universe and target feature-selection rule fixed, report retained genes and cross-arm gene overlap, and distinguish normalization refitting from HVG reselection and from correction.

### Phase 3: robustness and acceptance

Evidence for an acceptable remedy should include:

- Documented feature support and stable ID handling, including unresolved limitations and per-library specimen membership.
- No structurally unsupported genes driving the candidate final distance space.
- Reduced **held-out** provenance/coverage signatures not used to define the mask/model; reduced association with a fitted batch label is expected by construction and insufficient alone.
- Stable within-provenance distances/neighborhoods and non-technical signals on an agreed common feature space; no collapse of sample variability merely to force source mixing.
- sc-only, nuclei-only and supported within-provenance analyses; account for repeated donors in interpretation.
- Condition-label summaries examined only after freezing choices; no optimization of ARI/ANOSIM or disease recovery to select the repair.
- Transparent loss of genes/samples and uncertainty. Perfect source mixing or a smaller MDS gap is **not** the success criterion.

A persisting split after a defensible mask may represent additional technical effects, composition, anatomical sampling or genuine biology. An additive-only collapse does not establish that the correction was appropriate. If measurement comparability cannot be established, restrict the inference rather than manufacturing equivalence [R5].

## Primary-method references

- **R1.** AnnData `concat`, sparse outer-join zero-padding warning: https://anndata.readthedocs.io/en/stable/generated/anndata.concat.html
- **R2.** 10x feature matrices and reference construction: https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/outputs/cr-outputs-mex-matrices ; https://www.10xgenomics.com/support/software/cell-ranger/9.0/analysis/inputs/cr-3p-references ; counting algorithm https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm . These establish possible mechanisms, not the specific KPMP version/cause.
- **R3.** limma official manual and correction implementation: https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf ; https://raw.githubusercontent.com/bioc/limma/RELEASE_3_22/R/removeBatchEffect.R
- **R4.** DESeq2 authors' size-factor documentation: https://rdrr.io/bioc/DESeq2/man/estimateSizeFactors.html
- **R5.** Luecken et al., primary integration benchmark and batch-removal/biological-conservation tradeoff: https://doi.org/10.1038/s41592-021-01336-8 . Its cell-level rankings do not prescribe a pseudobulk method here.
- Source-specific processing-history evidence and specimen/GEO links are documented in the existing research note linked above. The source supplement says the CellXGene subset was rerun with GRCh38; it does not establish present mixed assemblies.

## Work performed for this discussion

Read existing checked diagnostic artifacts and relevant pipeline boundaries; verified the existing sample-inventory checksum; printed exact 28/17 membership; checked additional saved clinical metadata and multi-library counts; calculated small model-matrix ranks. Researched primary methodological documentation. No new source count scan, pipeline implementation, gene mask, correction, exclusion, new batch metadata field, configuration edit, scheduler submission or cohort rerun.
