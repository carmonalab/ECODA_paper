# ============================================================
# LIBRARY LOADING
# ============================================================

my_packages <- c(
  "doParallel",
  "factoextra",
  "foreach",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "ggtext",
  "mclust",
  "plotly",
  "scales",
  "anndataR",
  "Seurat",
  "limma",
  "Matrix",
  "HiTME",
  "stringr",
  "tidyr",
  "tidyverse",
  "rstatix",
  "EPIC",
  "reticulate",
  "MOFA2",
  "scITD",
  "GloScope",
  "funkyheatmap",
  "RColorBrewer",
  "forcats",
  "gtools",
  "scECODA",
  "pheatmap",
  "ncdf4",
  "arrow",
  "pROC",
  "Rfast",
  "zoo",
  "parallelly",
  "patchwork",
  "dplyr"
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
      "\n\nPlease add them to your pixi.toml (e.g., as 'r-packagename') and run `pixi install`."
    )
  }

  return(TRUE)
}

# 1. Attempt to load packages (will stop execution if any are missing)
invisible(load_my_packages(my_packages))

# Custom Negate operator
`%notin%` <- Negate(`%in%`)
