#!/usr/bin/env Rscript
# Regenerate Figure 3 (Figure 3A, 3B, 3C with Leiden_res_5, and composite Figure_3.pdf)
# and Supp_fig_18 companions directly from notebooks/benchmark_analysis.rmd chunks
# backed by data/derived/access_v2/analysis_release.rds.

if (basename(getwd()) == "notebooks") setwd("..")
root <- normalizePath(".")
Sys.setenv(
  PROJECT_ROOT = root,
  ECODA_BENCHMARK_RELEASE_RDS = file.path(root, "data/derived/access_v2/analysis_release.rds"),
  ECODA_BENCHMARK_RELEASE_MD5 = "fa146125e805593a70810229e47028ec"
)

rmd_lines <- readLines(file.path(root, "notebooks/benchmark_analysis.rmd"))
extract_chunks_between <- function(lines, start_line, end_line) {
  sub <- lines[start_line:end_line]
  out <- character()
  in_chunk <- FALSE
  for (ln in sub) {
    if (grepl("^```\\{r", ln)) {
      in_chunk <- TRUE
      next
    }
    if (in_chunk && grepl("^```\\s*$", ln)) {
      in_chunk <- FALSE
      next
    }
    if (in_chunk) out <- c(out, ln)
  }
  out
}

fig3_end <- grep('ggsave\\(file\\.path\\(path_plots, "Figure_3\\.pdf"\\)', rmd_lines)
stopifnot(length(fig3_end) == 1L)
# Include closing line + fence after Figure_3.pdf ggsave
fig3_end <- fig3_end + 3L
code_lines <- c(
  extract_chunks_between(rmd_lines, 44L, 275L),
  extract_chunks_between(rmd_lines, 1170L, fig3_end)
)
tmp_script <- tempfile(fileext = ".R")
writeLines(code_lines, tmp_script)
env <- new.env(parent = .GlobalEnv)
source(tmp_script, local = env)
unlink(tmp_script)
cat("DONE: Regenerated Figure_3_C_annotationmethods_pca_adams.pdf and Figure_3.pdf with Leiden_res_5.\n")
