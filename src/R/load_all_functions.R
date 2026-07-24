# ============================================================
# Main Loader
# Sources all modular R files in dependency order
# ============================================================

# Load all modules in dependency order
source("src/R/imports.R")           # Package imports
source("src/R/constants.R")         # Constants & configurations
source("src/R/helpers.R")           # Utility helpers (exec_time, etc.)
source("src/R/math_utils.R")        # Math/statistics utilities (clr, cv, impute_zeros)
source("src/R/scoring_metrics.R")   # Separation scoring & clustering metrics
source("src/R/pseudobulk.R")        # Pseudobulk & DESeq2 functions
source("src/R/hvcs.R")              # Highly Variable Cell Type functions
source("src/R/seurat_utils.R")      # Seurat object utilities
source("src/R/plotting.R")          # PCA, MDS plotting functions
source("src/R/benchmark_methods_r.R") # Individual method processing functions
source("src/R/benchmark_pipeline.R") # Pipeline orchestration (run_benchmark, run_transformation_analysis, etc.)