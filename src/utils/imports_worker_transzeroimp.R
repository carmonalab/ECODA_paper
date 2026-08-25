# ============================================================
# LIBRARY LOADING — trans/zeroimp workers
# ============================================================
# Slim import set for the ECODA transformation / zero-imputation workers
# (run_transformation_zeroimp_analysis/): obs-only backed reads, no counts
# matrix, NO Seurat. src/utils/imports.R is the notebook attach list (repo-
# wide env verification via src/utils/env_check.R);
# namespace-qualified packages (zCompositions::, Hotelling::, vegan::,
# mclust::, cluster::, igraph::, Matrix::, arrow::, ...) load on demand and
# must NOT be attached here. scECODA is intentionally NOT attached: the
# transformation analysis calls the LOCAL datrans() (benchmark_pipeline.R:5),
# not an scECODA API — scECODA has no references in any worker code path.
# datrans() does bare-call makeCluster (parallel), registerDoParallel
# (doParallel) and foreach, so doParallel/foreach are attached (pulling
# parallel/iterators via Depends).
#
# Attach order avoids masking regressions: doParallel/foreach before
# reticulate/dplyr.

my_packages <- c(
  "doParallel",
  "foreach",
  "reticulate",
  "dplyr",
  # zCompositions::cenGeoMean() calls survival::survreg() without a
  # namespace qualifier; attach survival so zero-imputation workers fail
  # closed at startup rather than halfway through the analysis.
  "survival"
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
