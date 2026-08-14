# QOT script / pipeline hotfixes (traceability)

This document records every change made to the vendored QOT script and the
QOT pipeline glue, beyond the unmodified upstream code. It exists so that any
deviation from PennShenLab/QOT @ `28cd529880c1` (2024-12-04) is auditable.

## Vendored script: `src/5_run_benchmark_methods/run_python_sample_embedding_methods/qot_utils_re.py`

Vendored from https://github.com/PennShenLab/QOT (MIT license), commit
`28cd529880c1`. The file header documents the same two changes; this section
adds the rationale, reproduction and verification for each.

### Hotfix 1 — fall-through groups silently dropped (the 3-5-sample bug)

**Change**: in `Gaussian_Mixture_Representation`, the `elif num_samples == 2:`
tail became an `else` branch (identical single-Gaussian representation: group
mean, zero covariances, weight 1).

**Why**: the original chain was
`if num_samples >= num_components + min_samples_gmm: GMM; elif num_samples == 1: ...; elif num_samples == 2: ...`
with NO `else`. Any (sample, cell-type) group with
`3 <= num_samples < num_components + min_samples_gmm` fell through and was
silently dropped from `params`, leaving zero-distance rows/cols in
`QOT_Distance`.

**Status**: NOT triggered at the benchmark settings (`num_components=1`,
`min_samples_for_gmm=0` -> threshold 1), fixed for robustness per the user's
request. Verified by a synthetic unit test (2 samples x 2 cell types, 4/2/4/2
cells per group, `num_components=2, min_samples_gmm=4`): all four
(sample, ct) keys present in `params`, non-zero pairwise distances.

### Hotfix 2 — GMM crash on small groups (blocker at benchmark settings)

**Change**: in `Gaussian_Mixture_Representation`, the GMM branch condition
gained `and num_samples > data.shape[1]`:
`if num_samples >= num_components + min_samples_gmm and num_samples > data.shape[1]:`.

**Why**: a full-covariance `sklearn.mixture.GaussianMixture` requires a
positive-definite empirical covariance. The covariance of n points in d PCA
dims has rank <= n-1, so it is singular (or numerically indefinite) whenever
n <= d and sklearn raises
`ValueError: Fitting the mixture model failed because some components have
ill-defined empirical covariance` (Cholesky `potrf` failure). At the benchmark
settings (`num_components=1`, `min_samples_for_gmm=0` -> threshold 1) the GMM
branch fires for EVERY group: 1-cell groups are rejected outright
(`ensure_min_samples=2`), and on the `_debug` dataset (JoaI 5-sample subset,
d=50) 13 of 35 groups failed — down to 1-cell groups and up to 48-cell
groups. Full datasets will hit it too via rare cell types in few samples.
Groups with n <= d now fall through to the single-Gaussian branches
(`== 1`, `== 2`, and the `else` branch from Hotfix 1), which was the chain's
obvious intent. GMM behavior is unchanged for every group where the
covariance is estimable (n > d), i.e. `min_samples_for_gmm=0` semantics
preserved exactly. Verified: on `_debug`, all 16 groups with n > 50 fit
cleanly; all 13 failures had n <= 50.

**Status**: verified end-to-end on `_debug` (see below) and by synthetic unit
tests (1-cell group takes the `== 1` branch; 4-cell groups at
`num_components=2, min_samples_gmm=4, d=5` take the `else` branch).

## Pipeline glue (NOT vendored): workarounds in

`src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py`

### Wrapper workaround 1 — duplicate-key rename bug (upstream bug, fixed in OUR wrapper)

**Upstream bug**: `Extract_Info` renames `{type_cell, id, progession}` in ONE
`rename(columns=...)` dict literal. When the sample column equals the
status/progression column (both `"Sample"` in our call), the dict literal
collapses the duplicate key to `{'Sample': 'status'}`: BOTH `Sample` columns
get renamed to `status`, no `sampleID` column survives, and the GMM step
(`groupby(['sampleID', 'Cell_type'])` / `df.drop(['sampleID'])`) raises
KeyError. Verified by simulation with the repo's py-cpu pandas.

**Our fix (run_qot)**: create a distinct temp obs column
`adata.obs["_bench_prog"] = adata.obs["Sample"]` and pass
`progession="_bench_prog"`. Side benefit: no bio label is passed to the
script (no-leakage).

### Wrapper workaround 2 — NaN cell-type labels

`fill_unknown_ct()` fills NaN in the cell-type column with `"Unknown"`
(scPoli pattern). The script itself filters `Cell_type != 'Unknown'`; NaN
would otherwise leak into the GMM groupby keys.

## Verification

- Synthetic unit test (Hotfix 1): PASSED.
- End-to-end `_debug` run: `--method qot --hvg 2000` -> 5x5
  `_debug_hvg2000_highres_qot_dists.feather` with sample IDs as index,
  `execution_times_*` rows `QOT_hvg2000_highres` / `QOT_hvg2000_lowres`.

## Known latent upstream issues left unfixed (dead paths for this benchmark)

- `sns.barplot(ci=None)` deprecation warning (`plot_hor_vs_vert`, subgroup-DE
  path only).
- Hardcoded `Results_PILOT` paths (`compute_diff_expressions`).
- `compute_shapley_values` / `trajectory_analysis` not used by this benchmark.
