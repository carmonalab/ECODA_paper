# ============================================================
# Worker Loader
# Sources the modular R files in dependency order, MINUS
# src/utils/imports.R (worker package attaches happen via
# imports_worker_core.R / imports_worker_transzeroimp.R) and
# src/utils/plotting.R (notebook-only).
# Notebooks keep using load_all_functions.R (full behavior).
# ============================================================

source("src/utils/datasets_io.R")
source("src/utils/constants.R")
source("src/utils/helpers.R")
source("src/utils/math_utils.R")
source("src/utils/scoring_metrics.R")
source("src/utils/pseudobulk.R")
source("src/utils/hvcs.R")
source("src/utils/seurat_utils.R")
source("src/5_run_benchmark_methods/benchmark_methods_r.R")
source("src/5_run_benchmark_methods/benchmark_pipeline.R")
