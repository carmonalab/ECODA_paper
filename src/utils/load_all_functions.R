# ============================================================
# Main Loader
# Sources all modular R files in dependency order
# ============================================================

# Load all modules in dependency order
source("src/utils/imports.R")
source("src/utils/datasets_io.R")
source("src/utils/constants.R")
source("src/utils/helpers.R")
source("src/utils/math_utils.R")
source("src/utils/scoring_metrics.R")
source("src/utils/pseudobulk.R")
source("src/utils/hvcs.R")
source("src/utils/seurat_utils.R")
source("src/utils/plotting.R")
source("src/benchmark/benchmark_methods_r.R")
source("src/benchmark/benchmark_pipeline.R")