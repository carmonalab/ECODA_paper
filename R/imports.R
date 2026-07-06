# ============================================================
# LIBRARY LOADING
# ============================================================

my_packages <- c(
  "pak",
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
  "Seurat",
  "limma",
  "SeuratData",
  "SeuratDisk",
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
  return(all(loaded)) # Returns TRUE only if ALL packages loaded successfully
}

# 1. Attempt to load packages
all_loaded <- load_my_packages(my_packages)

# 2. If any fail, run renv::restore() and try again
if (!all_loaded) {
  message("One or more packages are missing. Running renv::restore()...")

  # prompt = FALSE ensures scripts don't hang waiting for a user to type 'y'
  renv::restore(prompt = FALSE)

  # 3. Try loading again
  if (!load_my_packages(my_packages)) {
    stop(
      "renv::restore() finished, but some packages still could not be loaded. Please check for installation errors."
    )
  } else {
    message("Successfully restored and loaded all packages.")
  }
} else {
  # message("All packages successfully loaded.")
}

`%notin%` <- Negate(`%in%`)
