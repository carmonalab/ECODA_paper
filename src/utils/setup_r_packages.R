#!/usr/bin/env Rscript
# Install the R source/Bioconductor packages required by ECODA.
# This script is intentionally not a Pixi task: Pixi activation runs the
# r-base javareconf hook and writes shared R/etc files. Call it only through
# refresh_env.sh or setup_env_sbatch.sh after their environment lock and
# no-active-job checks have passed.

if (!identical(Sys.getenv("ECODA_ENV_MUTATION_GUARD"), "1")) {
  stop(
    "Refusing unguarded R environment mutation. Run ",
    "src/utils/bash/refresh_env.sh on the login node or submit ",
    "src/utils/bash/setup_env_sbatch.sh on Bamboo."
  )
}

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT is not set. Source src/slurm_config.sh first.")
}
setwd(project_root)

options(
  timeout = 9999999,
  repos = c(CRAN = "https://packagemanager.posit.co/cran/latest")
)

# Parallel source compilation: R CMD INSTALL builds with make; MAKEFLAGS=-jN
# compiles translation units in parallel. Without it every Rcpp-heavy package
# (Seurat, MOFA2, GloScope, ...) compiles single-threaded and setup takes
# 15-30+ min regardless of node cores. Cap at 8 on the shared login node;
# override with the R_SETUP_JOBS env var (e.g. 16 on a dedicated worker).
setup_jobs <- suppressWarnings(as.integer(Sys.getenv("R_SETUP_JOBS", unset = "")))
if (is.na(setup_jobs) || setup_jobs < 1) {
  setup_jobs <- min(8, max(2, parallel::detectCores() - 1))
}
Sys.setenv(MAKEFLAGS = paste0("-j", setup_jobs))
cat("Using", setup_jobs, "parallel compile jobs (MAKEFLAGS)")

# Keep pre-compiled conda-forge versions (r-igraph, r-leidenbase, ...) and only install
# missing packages. Without this, remotes upgrades the whole dependency closure from CRAN
# source (e.g. tries leidenbase 0.1.36 -> 0.1.37 even though Seurat has no version
# requirement on it), which fails on the HPC due to missing system dev libs.
Sys.setenv(R_REMOTES_UPGRADE = "never")
bioc_repos <- BiocManager::repositories()

# SHA-pinned GitHub packages (installed directly; remotes checks installed SHA and skips if matching)

# Basics
remotes::install_github("satijalab/seurat@e4cc89238b233e46e2cef4a2b866896cf0cb093c")
remotes::install_github("scverse/anndataR@07612e4f0e5c7efff05265eb916ae9ed94bf78b9")

# HiTME
remotes::install_github("carmonalab/SignatuR@0aded1db20d14f47b5b6074fcf5fc9f9c2e6279e")
remotes::install_github("carmonalab/STACAS@ea34e824e17df6316823db0a9b6322af0833e501")
remotes::install_github("carmonalab/ProjecTILs@1159e1778820180fd0233bbfc0d7f296e35bd25f")
remotes::install_github("carmonalab/HiTME@87f8d4147aaa468438df33f50cabf2005ac518bf")

# scATOMIC
# yikeshu0611/set does not contain a standard man/ (documentation) directory in its repository. When remotes::install_github attempts to run R CMD build on the source tree, R fails while building the documentation index (cannot open file 'man').
remotes::install_github("yikeshu0611/set@5360f0b85f0a618d9a05009b8553d397c4add828", build = FALSE)
if (!requireNamespace("dlm", quietly = TRUE)) remotes::install_version("dlm", version = "1.1-6.1", upgrade = "never")
if (!requireNamespace("Rmagic", quietly = TRUE)) remotes::install_version("Rmagic", version = "2.0.3", upgrade = "never")
remotes::install_github("inofechm/cutoff.scATOMIC@a43fd5e1e8f0e3b71ec446970e4316f305939a17")
remotes::install_github("abelson-lab/scATOMIC@d332cd5cf6a1ecef7c32d0adc4a862a4c47bcd95")

# scITD
if (!requireNamespace("edgeR", quietly = TRUE)) remotes::install_version("edgeR", version = "4.8.2", repos = bioc_repos, upgrade = "never")
if (!requireNamespace("sva", quietly = TRUE)) remotes::install_version("sva", version = "3.58.0", repos = bioc_repos, upgrade = "never")
if (!requireNamespace("scITD", quietly = TRUE)) remotes::install_version("scITD", version = "1.0.4", repos = bioc_repos, upgrade = "never")

if (!requireNamespace("limma", quietly = TRUE)) remotes::install_version("limma", version = "3.68.4", repos = bioc_repos, upgrade = "never")
if (!requireNamespace("MOFA2", quietly = TRUE)) remotes::install_version("MOFA2", version = "1.20.2", repos = bioc_repos, upgrade = "never")
if (!requireNamespace("GloScope", quietly = TRUE)) remotes::install_version("GloScope", version = "2.0.1", repos = bioc_repos, upgrade = "never")
if (!requireNamespace("thisutils", quietly = TRUE)) remotes::install_version("thisutils", version = "0.4.9", upgrade = "never")
remotes::install_github("GfellerLab/EPIC@50a4f404f96c2842b2891b517b4e3bfaa6c64b8f")
remotes::install_github("carmonalab/scECODA@055bc5efb51cac7c3b6e0ea7eae84420f0a8d21a")
remotes::install_github("funkyheatmap/funkyheatmap@d66dd6d65b4e29fb3f6a300627a4db367ad9db0c")

# ---------------------------------------------------------------------------
# Final integrity verification. pixi/conda do NOT verify installed package
# FILES (a corrupt/partial package dir — e.g. digest/Meta/package.rds missing
# after a concurrent or interrupted install — is invisible to a plain
# 'pixi install'), and conda (pixi install) + remotes (this script) both write
# into the same R library, so corruption can persist silently until a job
# hits it. Check every installed package has its Meta/package.rds and that the
# FULL loader package list loads; fail loudly with the wipe-and-reinstall repair.
# ---------------------------------------------------------------------------
check_lib_integrity <- function(lib) {
  bad <- character(0)
  for (p in list.dirs(lib, recursive = FALSE)) {
    if (!file.exists(file.path(p, "DESCRIPTION"))) next
    # R's own message-catalog component: copied at R build time, never
    # installed via R CMD INSTALL, so no Meta/package.rds by design
    # (false positive since the check was added).
    if (basename(p) == "translations") next
    if (!file.exists(file.path(p, "Meta", "package.rds"))) bad <- c(bad, p)
  }
  bad
}
bad_pkgs <- unique(unlist(lapply(.libPaths(), check_lib_integrity)))
if (length(bad_pkgs) > 0) {
  for (p in bad_pkgs) {
    pkg_line <- grep("^Package:", readLines(file.path(p, "DESCRIPTION"), n = 20), value = TRUE)
    cat("Missing Meta/package.rds:", p, paste(pkg_line, collapse = "; "))
  }
  stop(sprintf(
    "R library integrity check FAILED: %d package(s) missing Meta/package.rds (corrupt/partial install): %s. Repair by rerunning the guarded environment setup after verifying no jobs are active.",
    length(bad_pkgs), paste(head(bad_pkgs, 20), collapse = ", ")
  ))
}
# Full loader check: source src/utils/imports.R so this task validates exactly
# the package list that workers' load_all_functions.R requires (imports.R runs
# load_my_packages() at top level and stops with the same error a worker would
# produce). Catches failures that a plain install misses, e.g. a lockfile
# re-resolution pruning a conda R package that a remotes-installed package
# needs at load time (observed: 'there is no package called pracma' in jobs).
loader_error <- tryCatch({
  source("src/utils/imports.R")
  NULL
}, error = function(e) e)
if (!is.null(loader_error)) {
  stop(sprintf(
    "R library integrity check FAILED: full package load check (src/utils/imports.R) errored: %s. Repair by rerunning the guarded environment setup after verifying no jobs are active.",
    conditionMessage(loader_error)
  ))
}
cat("R library integrity check OK")
