# ECODA batch effect correction method comparison

User-approved implementation, 2026-11-16. Root implements and executes the
scoped derived comparison; becm_plot_facts owns notebook plotting and plot
relocation; becm_input_facts owns focused contract tests and scientific review.

## Approved scope

- Maintained code and description: notebooks/batch_effect_ECODA_becm_comparison.rmd.
- Figures: plots/batch_effect_ECODA_becm_comparison/{boxplots,mds,summary}.
- Seven current cohorts excluding Lung: Alzheimer, Breast_cancer,
  Covid19_PBMC, Joanito, Kidney_KPMP_full, Lupus_PBMC, Stephenson.
- Reuse the reviewed local control with MD5 e75b1036df9d80d21cd633058957f55f
  and its native source receipts. Kidney uses condition.l1 only. Lupus retains
  the explicitly approved whole-sample majority_v1 policy; no cell filter changes.
- Screen original author and Harmony/Leiden CLR independently using Wilcoxon
  for two batches, Kruskal-Wallis otherwise; BH over every original cell type.
- Remove BH p<0.05 AND max-min batch median CLR >= threshold. Screen once,
  then recenter retained CLR columns; no repeated screen or second pseudocount.
- Thresholds 0:0.5:4 for ANOSIM; 1,2,3 for nine-row/four-column MDS.
- Additional author threshold2 variant removes the union with raw zero counts
  in >=50% of samples, irrespective of significance.
- Twenty-four variants per cohort: two uncorrected anchors, three fixed
  corrected methods, eighteen drop settings, one sparse-filter setting.
- ANOSIM statistic only; no permutation p-values or mean-separation figure.
  Biology is evaluation-only. No donor adjustment or biological covariates.
- Author rows pair with original author CLR. All Harmony rows pair with
  original non-Harmony Leiden CLR. Fewer than two retained types or identical
  samples are explicit UNSCORABLE variants, never fabricated zero scores.

## Execution and evidence

This is the user's explicitly requested initial noncanonical scoring test,
using existing named sample-composition results, outside canonical outputs.
No pipeline submitter, preprocessing, annotation, Harmony, limma fit, full H5AD,
expression data, pseudobulk or canonical result is created or replaced.
Run the small derived calculation on a Bamboo compute allocation from the
committed notebook code; return only compact sample distances, labels, scores,
cell-type tests and grouped plotting summaries. Preserve native inputs and
frozen historical execution snapshots/checksum inventories.

Validate existing artifacts and reviewed parent receipts before reading; keep
each dataset's numerical work separate and release its source matrices before
the next. Record exact input selection, code/source revision, runtime and output
checksums. A relocation map preserves historical plot MD5s and old/new paths.
Focused contract fixtures and visual PDF inspection precede completion.

## Progress

- Input selection prepared from the approved control: 35 native parent
  references, seven sample metadata tables, exact full sample universes.
- Notebook/core/test implementation and plot relocation in progress.
