# ============================================================
# LIBRARY LOADING — R benchmark workers (subset of imports.R)
# ============================================================
# Slim import set for the R benchmark workers (prepare_pseudobulk and
# gloscope/mofa/pseudobulk/scitd via run_r_sample_embedding_methods/): only
# the packages these workers call BARE are attached. src/utils/imports.R
# remains the canonical env-verification list (guarded env-refresh smoke
# checks + notebook loader); namespace-qualified packages (DESeq2::, limma::,
# zCompositions::, vegan::, mclust::, cluster::, igraph::, Matrix::,
# GloScope::, arrow::, ...) load on demand and must NOT be attached here.
# MOFA2/scITD are attached per-method in 1.1.1_run_benchmark_methods_r.R
# (bare create_mofa / initialize_params + make_new_container) so gloscope/
# pseudobulk tasks never load them.
#
# Package order mirrors the relative order in imports.R (Seurat:18,
# reticulate:27, dplyr:46) to avoid name-masking regressions.

my_packages <- c(
  "Seurat",
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
