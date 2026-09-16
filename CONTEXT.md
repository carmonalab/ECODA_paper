# ECODA Domain Glossary

Unsupervised patient stratification, cell-type composition analysis, and representation benchmarking across single-cell disease cohorts.

## Cohorts & Biological Specimens

**Benchmark Cohort**:
A curated single-cell RNA sequencing dataset evaluated across embedding and compositional methods against known biological groupings.
_Avoid_: Dataset (too generic), batch, study

**Sample**:
The atomic biological specimen and operational unit of cohort stratification, corresponding to a single profiled single-cell library or donor aliquot. In cohort metadata, heterogeneous upstream keys (such as `specimen` or `donor_id`) are mapped into this canonical unit.
_Avoid_: Unit (too abstract), Patient (ambiguous when one donor has multiple serial or regional biopsies), Observation, Donor

**Cohort View**:
A standardized analytical subset and processing mode of a benchmark cohort designed for a specific scientific evaluation, such as standard benchmarking or batch-effect comparison.
_Avoid_: Dataset slice, configuration view, subset

## Cell Annotations & Composition

**Cell-Type Identity**:
A validated biological cell category or curated annotation label assigned to individual single cells.
_Avoid_: Cluster (unless explicitly unsupervised), cell label, cell state

**Annotation Granularity**:
The hierarchical resolution or biological specificity of cell-type identities assigned to cells (e.g., broad lineages versus fine-grained subpopulations).
_Avoid_: Annotation resolution (confusing with numeric clustering parameters), annotation level, cell depth

**Operational Cluster**:
An unsupervised, data-driven cell grouping (such as a Leiden cluster) treated as a synthetic cell type for composition benchmarking without asserting biological identity.
_Avoid_: Cell type, biological cluster

**Deconvolution Component**:
An inferred continuous cell-composition proportion derived from mixture modeling without discrete per-cell classification.
_Avoid_: Cell type, topic, factor

## Sample Representations & Stratification

**Sample Representation**:
The overall mathematical characterization of a sample used for unsupervised stratification, spanning compositional vectors, latent embeddings, pseudobulks, or pairwise distance matrices.
_Avoid_: Patient profile, sample vector

**Compositional Representation**:
A sample profile defined solely by the relative abundances of cell types or operational components across that sample, optionally restricted to highly variable cell types.
_Avoid_: Expression profile, cell proportion matrix

**Sample Embedding**:
A continuous, low-dimensional coordinate representation of samples learned from single-cell gene expression, cell-type distributions, or optimal-transport models.
_Avoid_: Projection, latent code

**Pseudobulk Representation**:
A sample profile constructed by aggregating raw or normalized gene-expression counts across single cells within each sample.
_Avoid_: Bulk RNA-seq, pooled sample

**ECODA**:
An unsupervised analytical framework and representation that stratifies patient cohorts using centered log-ratio (CLR)-transformed cell-type composition vectors.
_Avoid_: scECODA (reserved for the standalone R software package), compositional clustering

**Sample Feature Space**:
An explicit coordinate representation where each sample is positioned as a point in a shared $P$-dimensional feature matrix.
_Avoid_: Coordinate space, embedding space (when distances are used instead)

**Sample Distance Space**:
A representation of samples defined strictly by a pairwise dissimilarity matrix ($N \times N$) without explicit coordinate embeddings.
_Avoid_: Feature matrix, metric space

**Result Bundle**:
The standardized output triad containing a method's sample feature matrix, pairwise distance matrix, and ground truth labels for downstream metric evaluation.
_Avoid_: Output bundle, embedding dictionary

## Evaluation & Technical Noise

**Ground Truth Biological Label**:
An external clinical, disease, or biological annotation used strictly for post-hoc evaluation of unsupervised stratification recovery.
_Avoid_: Target variable, dependent variable, covariate, class label, feature (never an input to models)

**Biological Group Recovery**:
The degree to which unsupervised sample representations or distances preserve and separate known ground truth biological groups.
_Avoid_: Clustering accuracy, classification score, supervised accuracy

**Batch Variable**:
A recorded technical factor (such as sequencing platform, dissociation protocol, or processing site) that should not govern biological sample separation.
_Avoid_: Confounder, experimental factor

**Batch Confounding**:
An experimental scenario in which technical batch variables correlate with ground truth biological labels, making technical noise difficult to distinguish from biological variation.
_Avoid_: Batch effect (too general; batch effect is technical noise, confounding is correlation with biology)

**Batch Invariance**:
The degree to which a sample representation is free from technical batch effects and mixes samples across batches within the same biological condition.
_Avoid_: Batch correction quality, batch removal score
