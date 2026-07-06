# ============================================================
# LIBRARY LOADING
# ============================================================

suppressPackageStartupMessages({
  # install.packages("devtools")
  library(compositions)
  library(doParallel)
  library(factoextra)
  library(foreach)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ggtext)
  library(mclust)
  library(plotly)
  library(scales)
  # remotes::install_version(package = "SeuratObject", version = "5.2.0")
  library(Seurat)
  library(limma)
  # remotes::install_github('satijalab/seurat-data')
  library(SeuratData)
  # install.packages("hdf5r")
  # remotes::install_github("mojaveazure/seurat-disk")
  library(SeuratDisk)
  library(HiTME) # for default_black_list
  library(stringr)
  library(tidyr)
  library(tidyverse)
  # install.packages("zCompositions")
  # library(zCompositions)
  # library(robCompositions)
  library(rstatix)
  # library(peakRAM)

  library(EPIC)
  library(reticulate)
  library(MOFA2)
  library(scITD)
  library(GloScope)
  library(funkyheatmap)
  library(RColorBrewer)
  library(forcats)
  library(gtools)
  library(scECODA)
  library(pheatmap)

  library(ncdf4)
  library(arrow)
  library(pROC)
  library(Rfast)
  library(zoo)

  library(parallelly)
  library(patchwork) # To combine the plots

  # Import at the end, otherwise the "select" function gets re-assigned by another package
  library(dplyr)
})

`%notin%` <- Negate(`%in%`)
