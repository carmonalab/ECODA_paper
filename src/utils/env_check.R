# ============================================================
# ENV VERIFICATION — repo-wide R package check (no attach)
# ============================================================
# Verifies that every R package the repo needs is INSTALLED, without
# attaching any of them (requireNamespace, not require/library). Covers:
#   - the 20-package notebook attach list (src/utils/imports.R),
#   - packages used namespaced-only (pkg::fn) by the notebooks / utils,
#   - worker/annotation-only packages (HPC workers, 4_cell_type_annotation).
# The loaders themselves (imports.R / imports_worker_*.R) only check their
# own attach lists; env_check.R is the full-repo safety net (catches e.g.
# DESeq2/vegan/robCompositions or an annotation-only package going missing
# before a job runs).
#
# Intended to run in the environment being validated. On macOS, the default
# environment is the local analysis surface; Linux-targeted packages provided
# only by the py-cuda13 feature are checked on HPC instead.

env_check_packages <- c(
  # --- notebook attach list (src/utils/imports.R, in attach order) ----------
  "ggplot2", "ggpubr", "ggrepel", "ggtext", "GGally", "scales", "Seurat",
  "stringr", "tidyr", "dplyr", "purrr", "tibble", "forcats", "funkyheatmap",
  "RColorBrewer", "gtools", "scECODA", "pheatmap", "patchwork", "writexl",
  # --- namespaced-only (pkg::fn in notebooks / utils; never attached) --------
  "arrow", "jsonlite", "zCompositions", "limma", "DESeq2", "mclust", "cluster",
  "igraph", "vegan", "Matrix", "proxy", "parallelly", "BiocParallel",
  "factoextra", "Hotelling",
  # --- worker / annotation-only (HPC workers + cell type annotation) ---------
  "anndataR", "reticulate", "doParallel", "foreach", "MOFA2", "scITD", "EPIC",
  "GloScope", "HiTME", "scATOMIC", "cutoff.scATOMIC", "scGate", "ProjecTILs",
  "SignatuR", "R.utils", "robCompositions"
)

if (identical(unname(Sys.info()[["sysname"]]), "Darwin")) {
  env_check_packages <- setdiff(env_check_packages, "robCompositions")
}

check_env_packages <- function(pkgs) {
  missing_pkgs <- character(0)
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are missing from the pixi environment: \n",
      paste(missing_pkgs, collapse = ", "),
      "\n\nFor the HPC py-cuda13 environment:\n",
      "  - conda-available packages are changed by the guarded environment entry points, which run `pixi install` only after their lock and no-active-job checks\n",
      "  - GitHub/Bioc/CRAN-pinned packages are installed by `src/utils/setup_r_packages.R` through `src/utils/bash/refresh_env.sh` or `setup_env_sbatch.sh`; do not run `pixi run setup` directly.\n"
    )
  }

  cat("env_check_packages OK (", length(pkgs), " packages present)\n", sep = "")
  return(invisible(TRUE))
}

invisible(check_env_packages(env_check_packages))
