# ============================================================
# LIBRARY LOADING — trans/zeroimp workers (subset of imports.R)
# ============================================================
# Slim import set for the ECODA transformation / zero-imputation workers
# (run_transformation_zeroimp_analysis/): obs-only backed reads, no counts
# matrix, NO Seurat. src/utils/imports.R remains the canonical
# env-verification list (guarded env-refresh smoke checks + notebook loader);
# namespace-qualified packages (zCompositions::, Hotelling::, vegan::,
# mclust::, cluster::, igraph::, Matrix::, arrow::, ...) load on demand and
# must NOT be attached here. scECODA is intentionally NOT attached: the
# transformation analysis calls the LOCAL datrans() (benchmark_pipeline.R:5),
# not an scECODA API — scECODA has no references in any worker code path.
# datrans() does bare-call makeCluster (parallel), registerDoParallel
# (doParallel) and foreach, so doParallel/foreach are attached (pulling
# parallel/iterators via Depends, exactly as the full imports.R does).
#
# Package order mirrors the relative order in imports.R (doParallel:6,
# foreach:8, reticulate:27, dplyr:46) to avoid name-masking regressions.

my_packages <- c(
  "doParallel",
  "foreach",
  "reticulate",
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

# Attempt to load packages (will stop execution if any are missing)
invisible(load_my_packages(my_packages))
