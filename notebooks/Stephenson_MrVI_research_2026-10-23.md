# Stephenson: what MrVI learns beyond whole-sample pseudobulk

Research date: 2026-10-23. Research-only: no training, HPC jobs, environment changes or runtime checks. This note distinguishes **published/source evidence**, **mathematical inference**, and **hypotheses about the local benchmark**. The user's strong CD14-only pseudobulk separation is an observation, not a result independently recomputed for this literature slice. Numerical comparisons and the exact active runner belong to the coordinating agent's existing-artifact analysis.

## Bottom line

**Yes: MrVI can learn sample-specific expression effects at inferred cell-state resolution without supplied cell-type annotations.** It does not first label a cell “CD14” and then compute a CD14 pseudobulk. It learns a continuous shared cell-state representation, then learns how the effect of a sample ID changes with that state. The resulting sample comparisons are local to reference cells, so an expression change confined to monocytes need not first survive averaging monocytes together with every other lineage. This is the intended mechanism, not proof of which cells drive this particular benchmark result. [V2, Results §2.1; Methods A.1–A.3][V2]; [official guide][guide].

The CD14-only result therefore **supports the presence of useful expression information in the selected compartment**; it does not contradict MrVI's motivation. Nor does it establish that MrVI's advantage is caused by CD14, by interferon genes, or by expression rather than within-CD14 state mixtures. Those are attribution questions requiring local distances or counterfactual outputs. [Inference from the mechanisms below; published myeloid evidence: V2 §2.3][V2].

**Directly relevant precedent:** V2's COVID application is the **same Stephenson study, Newcastle subset**, reporting exactly **418,768 cells and 55 patients**, matching the local cohort's reported totals. However, V2 used **10,000 HVGs** rather than the local 2,000-HVG MrVI run, and its additional myeloid stratification concerned **days since symptoms**, not just binary COVID/Healthy separation. Shared source/cohort size is not an exact reproduction of model settings or endpoints. [V2 §2.3 and D, COVID experiment][V2]; local measurements in §8.

## 1. Precisely what is annotation-free?

| Information | Required or optional? | What it does |
|---|---|---|
| Per-cell transcript counts | Input | Learns expression/state structure. |
| Per-cell **sample-of-origin ID** | Required for meaningful sample-specific modeling/comparisons | Tells the model which cells belong to each donor/specimen/perturbation and selects a learned sample embedding. IDs are not disease labels. |
| Technical/nuisance covariate | Supply when correction is needed; a single-batch case is allowed | E.g. site or library preparation; effects enter the expression decoder separately. |
| Cell-type/state annotation | **Not required**; optionally used | Can guide a mixture-of-Gaussians prior, or be used after fitting for interpretation/aggregation. |
| Disease/severity/age labels | Not required for exploratory sample stratification | Used afterward to interpret clusters or specify comparative DE/DA covariates. |

Sources: [V2 A.1, A.4 and B][V2]; [guide, Preliminaries][guide]; [actual `setup_anndata` source][setup]. **“Annotation-free” means no mandatory discrete cell partition, not no metadata, no assumptions, or no sample IDs.** Naming a discovered population “CD14 monocytes” still needs marker/annotation interpretation; discovering that it has a distinct sample-distance pattern does not. [V2 §2.3, A.3][V2].

V2 itself allows annotation-guided fitting: the IBD example explicitly used a cell-type prior, while the COVID example specifies a 20-component mixture prior and subsequently interprets states with original annotations. Thus “MrVI can run annotation-free” must not become “all published MrVI fits ignored all annotations.” [V2 B, D: COVID and IBD model parameters][V2].

## 2. Orientation and mechanism: shared `u`, sample-aware `z`

1. **`u` is sample-unaware/intrinsic cell state**: intended to represent broad cell types and finer continuous states after harmonizing sample and nuisance effects. **`z` is sample-aware**: it retains cell state plus target/sample effects, while remaining nuisance-corrected. The orientation is not reversed. These are learned modeling aims, not a guarantee of perfect statistical independence in any fitted dataset. [V2 §2.1, A.1][V2]; [guide latent-variable table][guide].
2. In V2's variational description, an encoder estimates a Gaussian `q(u | x)`. The sample-aware representation is a deterministic residual mapping, schematically `z(u,s) = u + f(u,s)` (with a learned projection when dimensions differ). The neural mapping combines the continuous state with a learned embedding of sample ID using cross-attention. Because `f` depends on **both** state and sample, sample A can differ from sample B in monocytes but barely differ in T cells. A mixture prior on `u` does not turn each cell into a supplied type label or make the latent states discrete. [V2 A.1–A.2, B][V2]; [guide equations][guide]; [actual `EncoderUZ` source][uz].
3. For a reference cell, keep its inferred `u` fixed and change the sample ID: obtain `z(u,s₁), z(u,s₂), …`. These are **model-predicted counterfactual versions of the same reference state**, not observed matched cells, nearest-neighbor pairing, or the raw average expression of cells in each sample. In current source, inference keeps the original sample index for encoding `u` and substitutes `cf_sample` specifically for the `u → z` mapping. [V2 A.3][V2]; [actual inference source][inference].
4. The negative-binomial count decoder maps `z` and nuisance covariates to gene-expression distributions, accounting for library size and overdispersion. This is conditional probabilistic expression modeling, with shared parameters learned from cells across samples. “Denoising” means using the model's fitted latent/decoded expectation rather than raw noisy counts; it does **not** establish recovery of biological truth. The exploratory distance is computed in **latent `z` space**, not by comparing decoded pseudobulk gene vectors; decoding is used to interpret gene-specific DE effects. [V2 A.1, A.4][V2]; [actual negative-binomial decoder][decoder]; [distance source][dist].

**Important version caveat:** V2 and the present guide write the encoder as `q(u | x)`. The inspected later PyTorch source actually passes the original sample index into `EncoderXU`, using sample-conditioned normalization and a sample embedding. Thus “sample-unaware” describes the intended *representation*, not a literal assertion that sample ID is never used computationally when inferring `u`. This implementation detail must not be backdated to V2. [V2 A.2][V2]; [guide, Inference][guide]; [actual later encoder][xu].

## 3. Local distances, averaging, and composition

For the basic mean-representation calculation, the sample distance at reference cell `n` is

`Dₙ(a,b) = || z(uₙ,a) − z(uₙ,b) ||₂`.

Each cell therefore defines a full sample-by-sample matrix. V2 clusters cells by their flattened local distance matrices (PCA, then Leiden), averages distances **within each discovered group**, and clusters samples from that group-specific matrix. This can discover populations according to *how samples differ*, rather than requiring predefined cell types. [V2 A.3; COVID analysis D][V2]; [actual distance implementation][dist].

**Averaging distances is not averaging expression first.** In an unnormalized arithmetic aggregate `D̄(a,b)=Σₙ wₙDₙ(a,b)` with nonnegative weights, opposite directions of sample shifts across cell states do not cancel as signed expression vectors can: a large monocyte-specific distance remains a positive contribution. But a rare affected state can still be diluted by many near-zero distances, and averaging destroys the information about *where* a difference occurred unless local matrices/grouped outputs are retained. This is a mathematical inference from the distance definition, not a guarantee of better separation. [V2 A.3][V2].

Weighting changes the question. With equal weight per observed reference cell, abundant states dominate, and donors yielding more cells contribute more reference states to the pooled average. Giving each donor equal reference mass or giving each cell type equal mass is a different aggregate. Current grouped reduction sums cell distances and divides by the number of cells in the group; it does not automatically equalize donors or cell types within that group. An all-identical dummy grouping therefore averages over the entire pooled reference set, not a type-balanced atlas. [Actual grouped reduction source][aggregation].

**Do not conflate the distance output with MrVI's DA module.** Local `z` distances compare predicted sample effects **at a fixed state**. DA instead compares sample/group densities of cells in `u` space using aggregated posteriors. V2's semi-synthetic example deliberately separates these: distance stratification is strong in the DE subset and small in DA-only subsets, while the DA analysis detects the altered abundances. Thus a benchmark exporting only averaged local distances is not automatically exporting MrVI's full differential-abundance information. Nevertheless reference-state weighting, estimation, and imperfect disentanglement can connect a fitted aggregate to composition; “purely expression, entirely composition-independent” would also be too strong. [V2 §2.2, A.4, Fig. S1][V2]; [guide, Differential abundance][guide].

## 4. What is and is not preserved relative to pseudobulk?

- **Whole-sample pseudobulk** collapses cells before comparing samples. Before transformations, its counts combine cell abundance, cell-specific RNA output and state-specific expression; unrelated compartments can dominate or obscure a localized program. A CD14-restricted pseudobulk removes much of that between-lineage mixing. This mixture explanation is mathematical intuition, not an exact equivalence to the benchmark's normalization/HVG/PCA pipeline. [Inference; local preprocessing must be checked separately.]
- **Per-cell-type pseudobulk** preserves expression separately for supplied types but still averages sub-states inside each type. MrVI instead evaluates a smooth state-conditioned sample response and can produce both latent distances and decoded gene effects. Calling it “annotation-free per-cell-type pseudobulk” is a helpful rough analogy for *retaining locality*, but mechanistically incorrect. [V2 A.2–A.4][V2].
- **Within-type expression and state abundance are not cleanly exclusive biological stories**: a good CD14 pseudobulk can reflect altered expression within a stable state, changed fractions of activated/IFN-responsive states inside the CD14 label, or both. The original Stephenson analysis resolves multiple CD14 states, so that possibility is not merely hypothetical at the dataset level. [Stephenson, Results “Mononuclear phagocytes and HSPC changes”, Fig. 2][Stephenson].

## 5. Defaults and version boundaries

The anchor paper is **bioRxiv v2, posted 10 May 2024**, not the earlier 2022 version or the later Nature Methods article. Its revision footnote explicitly says the generative model changed to cross-attention plus a mixture-of-Gaussians prior. V2 describes Euclidean distances and **optional** Monte Carlo uncertainty-based z-scoring. It does not justify assuming every analysis used normalization; its IBD Methods explicitly specifies normalized distances and 20 Monte Carlo samples. [V2 revision footnote; A.3; D, IBD Analysis][V2].

The current stable API/source inspected here cites the **2025 Nature Methods** publication and identifies its implementation as **PyTorch**, replacing the earlier JAX backend. Its `get_local_sample_distances` defaults are `use_mean=True`, `normalize_distances=False`, `norm='l2'`, `mc_samples=10`; Monte Carlo sampling matters when sampled/normalized representations are requested, not for the default mean-distance route. Source normalization subtracts a same-original-sample baseline mean and divides by its standard deviation, then averages Monte Carlo results. Therefore “normalized” here is uncertainty-based z-scoring, **not** normalization of cell-type composition or library-size normalization. [Current API][api]; [actual distance/default/normalization source][defaults].

The unmaintained standalone `YosefLab/mrvi` package also exposes mean, unnormalized L2 defaults and directs users to `scvi.external.MRVI`; its current `main` branch is not a timestamped copy of the code used for V2. Exact historical benchmark defaults require its actual scvi-tools version or saved model metadata, not today's documentation. No such historical-version identity is asserted in this note. [Standalone repository][standalone]; [standalone source][oldsource].

## 6. Support limitations: absent states and confounding

A counterfactual is not reliable simply because the network can return a vector. If sample B has no cells near the reference state, its predicted response there is extrapolation. V2 explicitly calls this insufficient information and filters inadmissible cell/sample pairs for **DE**, using sample-specific aggregated posterior support; Methods B describes a fifth-percentile density threshold. This does not establish that the benchmark's plain distance-export call applies the same filter. The inspected distance routine computes all sample counterfactuals without invoking an outlier filter. [V2 A.4, B “Out-of-distribution checks”][V2]; [actual local-statistics source][aggregation].

Similarly, “intrinsic” is a modeling decomposition, not a causal guarantee. Sample embeddings can capture biological differences, technical differences not modeled as nuisance, or confounding. The original GloScope paper's Stephenson example is particularly cautionary: uncorrected sample densities showed strong sequencing-site separation. [Inference from conditioning design; direct batch evidence: GloScope Fig. 7][GloScope].

This full-cohort batch example must not be transferred uncritically: the coordinating agent reports that the local cohort is restricted to `Site == 'Ncl'` and COVID/Healthy, excluding `BGCV10_CV0198` and `MH8919230`. The original GloScope Fig. 7 covers 143 samples across sites; the local matched comparison has 55. [GloScope Fig. 7][GloScope]; coordinating agent's configuration/artifact trace.

## 7. Primary-source biological context and GloScope contrast

**Direct MrVI V2 evidence:** the COVID analysis used **Newcastle only**, 418,768 cells from 55 patients, 10,000 HVGs, patient ID as target, no batch identifier, latent dimensions `u=5`, `z=30`. It discovered a monocyte/DC-associated sample stratification related to days since first symptoms; a lower-DSS patient group had an inferred IFN-related module including **GBP1, IFITM3 and IFI27**. This is evidence that the method can recover myeloid-local expression variation without disease labels as model targets. It is not a quantitative replication of the local all-sample ECODA setup or proof that local CD14 drives the benchmark. [V2 §2.3; D, COVID experiment][V2].

**Original Stephenson evidence:** the study reports type I/III interferon-response expression in monocytes, DCs and HSPCs across COVID severities; it resolves classical CD14, activated CD83+ CD14 and proliferating monocyte states, plus distinct CD16 states. Its strongest type-I-IFN response among myeloid subsets was in C1QA/B/C+ CD16 monocytes—not exclusively CD14. These observations make compartment-local expression/state signals biologically plausible while cautioning against attributing all IFN signal to CD14. The paper is **Nature Medicine (2021)**; MrVI V2's reference entry labels its journal incorrectly as Nature Methods. [Stephenson Figs. 1d, 2a,b,f and corresponding Results][Stephenson].

**GloScope uses a different route:** it embeds cells jointly (e.g. PCA), represents each sample by its distribution/density of cells in that space, and compares samples by the square root of symmetrized KL divergence. It does not align a fixed reference state counterfactually across samples. Because an entire distribution is retained, it can reflect both state proportions and expression/location/shape changes, without cell-type annotations. Its publication reports Gaussian-mixture and nearest-neighbor density strategies, and explicitly notes that composition-versus-DE distinctions can blur. These are general method properties, not proof of the local implementation's signal origin. [GloScope, Overview; simulations; Methods][GloScope].

The coordinating agent traces the local GloScope implementation to a nearest-neighbor estimator using 30 PCs; the publication frequently illustrates Gaussian-mixture densities with 10 PCs. Neither the paper's illustrative defaults nor its full-cohort site finding should be substituted for that local setup. [GloScope Figs. 2–7 and Methods][GloScope]; coordinating agent's runner trace.

## 8. Existing-artifact quantitative context (parent-measured)

The coordinating agent reports checking input size/MD5 and the same **55 sample IDs and labels (43 COVID, 12 Healthy)** across these representations, with **418,768 total cells** and **103,326 CD14 cells**. These are the coordinating agent's measurements, not new analyses performed by this research slice. The metrics were recomputed from the repository's full sample-distance matrices using ANOSIM with zero permutations, modularity at `k=3`, and mean Ward/PAM ARI:

| Existing representation | ANOSIM R | Modularity | Mean ARI |
|---|---:|---:|---:|
| ECODA author high-resolution | 0.39532 | 0.31537 | 0.31049 |
| Whole-sample pseudobulk, 2,000 HVGs | 0.28478 | 0.43297 | 0.43760 |
| GloScope, 2,000 HVGs / 30 PCs | 0.40633 | 0.52663 | 0.30548 |
| MrVI, 2,000 HVGs | 0.58047 | 0.80711 | 0.76817 |
| CD14 pseudobulk, 2,000 HVGs | 0.50146 | 0.43518 | 0.34640 |

Thus CD14 pseudobulk has substantially stronger ANOSIM than whole-sample pseudobulk, but **lower ARI**: “better separation” is not uniform across metrics. MrVI exceeds whole-sample pseudobulk on all three reported metrics; GloScope does not (its ARI is lower). Zero ANOSIM permutations here provide an R statistic, not a permutation P value. [Coordinating agent, existing-artifact measurements supplied 2026-10-23.]

The coordinating agent also reports a descriptive Spearman correlation across 1,485 unordered sample-pair distances of **0.81414 for MrVI versus CD14 pseudobulk**, compared with **0.50587 for MrVI versus whole-sample pseudobulk**. This is useful concordance, **not CD14 attribution**: pairs are non-independent and shared COVID/Healthy structure can contribute. The saved MrVI matrix is reported exactly symmetric with zero diagonal (Feather MD5 `1619edbd9427c09f94b68ede2c55007c`). [Same parent analysis.]

CD14 PC1 explains **29.52%** of variance; its ten largest absolute loadings are **IFI27, IFITM1, CLU, RNASE1, IFITM3, IFI6, OTOF, SIGLEC1, IFI44L, LGALS3BP**. This supports an IFN-related component in that particular pseudobulk representation, but does not identify the genes driving MrVI. Only **984 of the 2,000 HVGs** overlap between CD14 and whole-sample pseudobulk, so this is not a same-gene compartment restriction: subset-specific normalization and feature selection also change. [Same parent analysis.]

**Within-condition sensitivity (same parent calculation):** the closer MrVI/CD14 geometry persists within the **43 COVID samples**: Spearman **0.70295** for MrVI–CD14 versus **0.40096** for MrVI–whole-pseudobulk across 903 pairs. Within Healthy (66 pairs), the respective values are **0.44407** versus **−0.13700**; across COVID/Healthy (516 pairs), **0.83846** versus **0.45141**. Thus the overall concordance is not solely the binary between-group gap, but remains descriptive: sample pairs are dependent, no correlation P values are claimed, and this still does not isolate CD14's causal contribution to the trained model.

**Local evidence provenance supplied by the coordinating agent:** [composition RDS](../data/benchmark/results/Stephenson_composition.rds), [pseudobulk RDS](../data/benchmark/results/Stephenson_pseudobulk.rds), [GloScope RDS](../data/benchmark/results/Stephenson_gloscope.rds), [metadata RDS](../data/benchmark/results/Stephenson_metadata.rds), and [MrVI Feather](../data/benchmark/embeddings/Stephenson_hvg2000_mrvi_dists.feather). CD14 inputs are under `data/derived/stephenson_cd14_pseudobulk_20261022` (SHA256 manifest checked by parent); the existing loader is [`stephenson_ucell.R`, lines 250–309](../src/notebooks/stephenson_ucell.R), with existing [PCA score CSV](../plots/exploratory_Stephenson_CD14_LR_pseudobulk_hvg2000_PCA_scores.csv). Scoring follows [`scoring_metrics.R`](../src/utils/scoring_metrics.R): `calc_modularity(d,y,3)`, `clust_eval(d,y,return_mean=TRUE)`, and `vegan::anosim(d,y,permutations=0,parallel=1)`. Parent reconstructed and reordered both MrVI matrix axes to named sample metadata, then scored its supplied distance matrix directly. CD14 PC1 used `prcomp(center=TRUE, scale.=FALSE)`; distance correlations used lower triangles.

## 9. Strongest targeted attribution test without retraining

**Recommended analysis, not executed:** use the existing fitted MrVI model and/or saved per-cell distance tensor to compute matched sample-pair matrices for (i) CD14 references only, (ii) non-CD14 references only, and (iii) the original all-cell weighting. Use annotations only *post hoc*. Quantify separation with the same sample cohort, endpoint, distance convention and benchmark statistics; also report donor-balanced and cell-type-balanced reference weighting, and assess cell/sample support. For a simple unnormalized mean, verify the exact decomposition `D_all = p_CD14 D_CD14 + (1−p_CD14) D_nonCD14`. A corresponding reduction in the *disease-separating component* after leaving CD14 references out is much stronger attribution evidence than correlating whole MrVI distances with CD14 pseudobulk alone. [Proposed analysis based on local-matrix definition and aggregation][V2]; [source][aggregation].

Inspect finer CD14-state distance contributions and decoded counterfactual gene effects to test an IFN hypothesis, rather than identifying a gene program from geometry alone. Keep caveats explicit: leave-reference-group-out attribution isolates the contribution to the **distance aggregation**, not the effect of removing those cells from model training; all parameters were still learned jointly. If only the final averaged matrix was saved and no fitted model or per-cell outputs remain, CD14 attribution cannot be recovered exactly from that matrix alone. No new training is needed when the fitted model and compatible inputs remain available. [Inference from information discarded by aggregation; counterfactual/DE mechanisms in V2 A.3–A.4][V2].

**Availability in this repository:** the coordinating agent's runner trace reports that the fitted model was **not saved**, and `keep_cell=False` discarded the per-cell tensor. Consequently the saved pooled Feather cannot support the proposed CD14 leave-reference-out analysis on its own. Under the no-training constraint, the reachable result is the descriptive distance concordance above; exact attribution requires recovery of an already-fitted model or local outputs from elsewhere. This is a specific missing artifact, not a reason to substitute pseudobulk correlation as mechanistic proof. [Coordinating agent's current runner/output trace.]

## Sources and access

The V2 full HTML article, including main text, Methods A–D, supplementary captions and references, was accessible. Its displayed mathematical equations are partly image-based; the explicit formulas above were cross-checked against the official guide and actual model code. Later implementation details are pinned where possible and are not claimed to be the historical benchmark implementation.

[V2]: https://www.biorxiv.org/content/10.1101/2022.10.04.510898v2.full
[guide]: https://docs.scvi-tools.org/en/stable/user_guide/models/mrvi.html
[api]: https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.MRVI.html
[setup]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_model.py#L186-L222
[aggregation]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_model.py#L426-L663
[dist]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_model.py#L698-L755
[defaults]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_model.py#L665-L940
[uz]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_module.py#L177-L295
[xu]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_module.py#L297-L366
[inference]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_module.py#L536-L605
[decoder]: https://github.com/scverse/scvi-tools/blob/e48b259fadee5a30947d6c266edb64ca73da0091/src/scvi/external/mrvi/_module.py#L39-L175
[standalone]: https://github.com/YosefLab/mrvi
[oldsource]: https://github.com/YosefLab/mrvi/blob/main/src/mrvi/_model.py#L626-L731
[Stephenson]: https://www.nature.com/articles/s41591-021-01329-2
[GloScope]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11463121/
