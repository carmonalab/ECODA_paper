# ============================================================
# Main Loader
# Sources all modular R files in dependency order
# ============================================================

# Load all modules in dependency order
source("R/imports.R")           # Package imports
source("R/constants.R")         # Constants & configurations
source("R/helpers.R")           # Utility helpers (exec_time, etc.)
source("R/math_utils.R")        # Math/statistics utilities (clr, cv, impute_zeros)
source("R/scoring_metrics.R")   # Separation scoring & clustering metrics
source("R/pseudobulk.R")        # Pseudobulk & DESeq2 functions
source("R/hvcs.R")              # Highly Variable Cell Type functions
source("R/seurat_utils.R")      # Seurat object utilities
source("R/plotting.R")          # PCA, MDS plotting functions
source("R/benchmark_methods.R") # Individual method processing functions
source("R/pipeline.R")          # Pipeline orchestration (run_benchmark, run_transformation_analysis, etc.)