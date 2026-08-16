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
# Intended to run on the HPC (login node, or sourced by the guarded
# env-refresh smoke checks in setup_env_sbatch.sh / refresh_env.sh). On the
# macOS `default` env HPC-only packages (py-cuda13 / linux-64-scoped or
# annotation-only) may legitimately be absent and WILL be reported — this is
# expected; do not "fix" it by installing them locally.

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
      "\n\nTwo install paths (HPC py-cuda13 env):\n",
      "  - conda-available packages: add 'r-<name>' to [dependencies] in pixi.toml, then run `pixi install`\n",
      "  - GitHub/Bioc/CRAN-pinned packages (Seurat, anndataR, MOFA2, scITD, HiTME, ...): run `pixi run -e py-cuda13 setup`\n"
    )
  }

  cat("env_check_packages OK (", length(pkgs), " packages present)\n", sep = "")
  return(invisible(TRUE))
}

invisible(check_env_packages(env_check_packages))
