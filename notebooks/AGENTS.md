# Interactive `benchmark_analysis.rmd` workflow

## User requirements

- Orchestrate subagents to read/explore/implement/plot when possible and it makes sense (e.g. for editing just a line of code it does not make sense to spawn a subagent).
- Keep one persistent interactive R session so benchmark objects remain in memory while figures are iterated.
- For figure changes, execute only the affected plotting code and save the resulting PDF/plot under `plots/` for inspection; do not rerun the whole notebook for every edit.
- Always check a figure update for:
  - Any text (e.g., title, axis labels) is legible (not too small, not too big, not overlapping, and not cut off).
  - Numbers do not have more than two significant digits and have consistent number of digits across a plot (including also sub plots within a figure).
  - Axis limits are appropriate for the data.
  - No points fall outside of the axis limits.
  - When possible, choose meaningful tick intervals (e.g. 0.00, 0.25, 0.50, 0.75, 1.00 etc. or -0.25, 0.00, 0.25, 0.50 etc. or 0.0, 0.5, 1.0 with additional subticks at 0.25 or 0.1 intervals). Usually 3-5 major ticks are enough (and minor ticks if needed, e.g. most points fall within 0-1 (so major ticks at 0.0, 0.5, 1.0) but some points might be around -0.1 to -0.2).
- the default scoring methods should be anosim, modularity (knn3), adjusted rand index (ari) (the other scoring methods are for comparison only and are shown in diagnostic supplemental figures such as Supp_fig_X2; Supp_fig_2 is the explicit five-metric screening exception).
- New additional figures will be highlighted with an "X" (e.g. Supp_fig_1 and a new figure after that should be added as Supp_fig_X2). Renumbering will be done in the end and will be initiated by the user.

## Findings

- The combined `Supp_fig_1.pdf` is currently written at 14 x 17.5 inches. The pairwise matrix was originally written at 7 x 7 inches, while the mean row was written at 15 x 3.5 inches. Composing the matrix into the wider combined figure enlarges its cells without enlarging fixed text grobs, so Figure 1A labels look smaller. This is primarily a grob-scaling/layout issue; the 1B dimensions and the `heights = c(4, 1)` composition also determine the final panel proportions.
- Legacy chunks around lines 1610--1711 use `eval=FALSE`. A blind `knitr::purl()` followed by `source()` does not preserve that execution policy and can run legacy code unintentionally.

## Launching the persistent R session

From the repository root, use the Pixi default environment:

```bash
pixi run -e default R
```

Do not use `R --vanilla` when relying on the repository `.Rprofile`; the profile points reticulate at the Pixi Python interpreter. If `--vanilla` is required, source `.Rprofile` explicitly.

Within the Oh My Pi harness, a long-lived REPL must be started with the process manager rather than a foreground shell command:

```text
hub start
  name=benchmark-r
  application=pixi
  args=["run", "-e", "default", "R"]
  cwd=<repository-root>
  pty=true
  ready.log=">"
```

Keep that process alive between figure edits. Send R expressions to it and inspect output with the corresponding process-manager operations; do not start a second R process for each plot.

## Initial or recovery execution

When a full environment rebuild is necessary (as specified by the user), execute the notebook faithfully so chunk options are honored:

```r
knitr::knit(
  input = "notebooks/benchmark_analysis.rmd",
  output = tempfile(fileext = ".md"),
  envir = .GlobalEnv
)
```

This full execution is for initial loading or recovery only. For normal figure iteration, evaluate the changed plotting chunk directly in an R session (either use an existing one or spin up a persistent one that should be closed upon task completion. ask user if task is completed and wait for confirmation before closing the r session. also remind the user to do so if you think the task is completed) and save its output with `ggsave()`.
