# ============================================================
# LIBRARY LOADING — notebook attach list (20 packages)
# ============================================================
# Attach list for the two analysis notebooks
# (notebooks/benchmark_analysis.rmd, notebooks/batch_effect_analysis.rmd)
# via src/utils/load_all_functions.R. Repo-wide env verification (attach ∪
# namespaced-only ∪ worker/annotation packages) lives in
# src/utils/env_check.R — do NOT grow this list for namespaced-only or
# worker-only packages.
#
# Order encodes the former effective attach order (incl. tidyverse's member
# attach at its old position): dplyr/purrr/tibble sit right after tidyr,
# forcats after tibble — do NOT alphabetize (masking behavior depends on it).

my_packages <- c(
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggtext",
  "GGally",
  "scales",
  "Seurat",
  "stringr",
  "tidyr",
  "dplyr",
  "purrr",
  "tibble",
  "forcats",
  "funkyheatmap",
  "RColorBrewer",
  "gtools",
  "scECODA",
  "pheatmap",
  "patchwork",
  "writexl"
)

# Function to attempt loading all packages quietly
load_my_packages <- function(pkgs) {
  suppressPackageStartupMessages({
    # require() returns a logical vector (TRUE for success, FALSE for failure)
    loaded <- sapply(pkgs, require, character.only = TRUE, quietly = TRUE)
  })

  # Identify which packages failed to load
  missing_pkgs <- pkgs[!loaded]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are missing from the pixi environment: \n",
      paste(missing_pkgs, collapse = ", "),
      "\n\nFor the HPC py-cuda13 environment:\n",
      "  - conda-available packages are changed by the guarded environment entry points, which run `pixi install` only after their lock and no-active-job checks\n",
      "  - GitHub/Bioc/CRAN-pinned packages are installed by `src/utils/setup_r_packages.R` through `src/utils/bash/refresh_env.sh` or `setup_env_sbatch.sh`; do not run `pixi run setup` directly.\n"
    )
  }

  return(TRUE)
}

# 1. Attempt to load packages (will stop execution if any are missing)
invisible(load_my_packages(my_packages))

# Custom Negate operator
`%notin%` <- Negate(`%in%`)
