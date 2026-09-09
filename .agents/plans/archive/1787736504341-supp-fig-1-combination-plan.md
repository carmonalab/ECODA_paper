# Supp fig 1 combination plan

## Problem

`notebooks/benchmark_analysis.rmd` currently writes Supp fig 1A (the `GGally::ggpairs` pairwise metric matrix) and Supp fig 1B (the patchwork row of metric-versus-mean plots) as separate PDFs. The mean-row panels use a shared, overly broad limit and do not explicitly enforce square panels; the correlation matrix and mean plots also use inconsistent/default tick formatting.

## Goal

Produce one combined `Supp_fig_1.pdf` containing the pairwise correlation matrix and the metric-versus-mean row, while retaining the existing individual outputs only if they remain useful as intermediate artifacts. Make the metric-versus-mean panels square and give every axis four consistently formatted two-decimal ticks without dropping data.

## Decisions

- Reuse the notebook's established `patchwork::wrap_plots()` composition rather than adding `cowplot` or `gridExtra`.
- Keep the pairwise matrix square and wrap it as a patchwork element with `patchwork::wrap_elements()` because `ggpairs()` returns a `ggmatrix`.
- Derive each metric's range from finite observed values, add a small 5% span padding, and generate exactly four evenly spaced breaks. Apply the same per-metric range/breaks to the corresponding x/y scales so panels do not become visually flat across metrics.
- Use `scales::label_number(accuracy = 0.01)` (or an equivalent explicit formatter) for stable two-decimal labels. Keep `coord_cartesian()` semantics so smoothing/correlation inputs are not discarded.
- Use `coord_fixed()`/`theme(aspect.ratio = 1)` for each mean panel and choose a combined output size that gives the matrix and mean-row panels readable, equal-sided individual plots.
- Do not alter `datasets.json`, benchmark data, or pipeline outputs.

## Implementation

1. Update the Supp fig 1 chunk's range/break helper and `plot_mean_corr()` so each metric panel computes finite metric/mean bounds, pads degenerate ranges safely, creates four breaks, formats them to two decimals, and uses fixed aspect ratio.
2. Build `mean_row` with the existing `wrap_plots()` convention and wrap the `ggmatrix` with `wrap_elements(full = pairwise_matrix)`.
3. Compose a vertical `supp_fig_1` patchwork with tags `A` and `B`, save it as `Supp_fig_1.pdf`, and preserve the existing A/B PDFs as compatibility/intermediate outputs unless rendering proves they conflict.
4. Add a minimal, data-backed render/smoke path that evaluates only the Supp fig 1 chunk or an equivalent extracted R script, because the full notebook may fail when in-progress datasets/results are absent.

## Acceptance criteria

- The notebook source has one combined Supp fig 1 output using patchwork, with no dependency on new plotting packages.
- Mean-row panels have equal height and width and each metric panel uses four two-decimal ticks on both axes.
- Bounds are derived from available finite data with padding; negative silhouettes remain visible and no points are clipped.
- Pairwise and mean plots remain generated from the same `df_wide` data and retain their intended labels/smoothing.
- A focused render or smoke test succeeds against available benchmark data, and any inability to full-knit is reported with the exact missing/incomplete artifact rather than masking errors.
