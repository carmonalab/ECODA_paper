# Run from the repository root: pixi run -e default Rscript tests/test_publication_zeroimp_panels.R
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("src/utils/constants.R")
source("src/utils/evaluation_scores.R")
source("src/notebooks/benchmark_analysis.R")

# Evaluate only the recipe selection and panel function from the maintained
# notebook. No artifact reads, benchmark evaluation, or figure writes occur.
lines <- readLines("notebooks/benchmark_analysis.rmd")
first <- which(startsWith(lines, "zeroimp_counts_order <-"))
last <- which(startsWith(lines, "p2 <- zeroimp_metric_panel("))
stopifnot(length(first) == 1L, length(last) == 1L, last > first)
setup <- parse(text = lines[first:(last - 1L)])

methods <- c(paste0("counts_all_", c("0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")),
  paste0("counts_zeros_", c("0.5", "2/3", "1", "5", "10", "20", "50", "100", "200")),
  paste0("multLN_", c("0.001", "0.01", "0.1", "1", "2", "5"), "%"))
df_zeroimp <- expand.grid(dataset = paste0("cohort", 1:5),
  method = c(methods, "multRepl_1%"), score = evaluation_metrics(), stringsAsFactors = FALSE)
df_zeroimp$value <- match(df_zeroimp$dataset, paste0("cohort", 1:5)) / 32 +
  match(df_zeroimp$method, c(methods, "multRepl_1%")) / 1024
original <- df_zeroimp
eval(setup)
stopifnot(identical(df_zeroimp, original), identical(zeroimp_order, methods),
  identical(levels(zeroimp_companion_scores$family), c("counts_all", "counts_zeros", "multLN")),
  !any(grepl("multRepl", zeroimp_companion_scores$method)))

panel <- zeroimp_metric_panel("anosim_score")
built <- ggplot_build(panel)
stopifnot(nrow(built$layout$layout) == 3L,
  identical(as.character(built$layout$layout$family), c("counts_all", "counts_zeros", "multLN")))
annotations <- panel$layers[[length(panel$layers)]]$data
stopifnot(nrow(annotations) == length(methods) - 1L,
  all(annotations$label == "*")) # Five exactly tied positive differences, normal approximation.

# Reordering source rows cannot silently change paired tests.
set.seed(31)
df_zeroimp <- df_zeroimp[sample(nrow(df_zeroimp)), ]
eval(setup)
reordered <- zeroimp_metric_panel("anosim_score")
stopifnot(identical(annotations, reordered$layers[[length(reordered$layers)]]$data))

# Missing one dataset/recipe/metric cell must fail before a partial panel appears.
df_zeroimp <- original[-1L, ]
stopifnot(inherits(try(eval(setup), silent = TRUE), "try-error"))
cat("Publication zero-replacement facet and pairing checks passed.\n")
