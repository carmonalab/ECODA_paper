# ==============================================================================
# 2.0_create_scgate_db.R — One-time scGate model DB download / cache creation
# ==============================================================================
# Called by the canonical Stage 4 submitter on a compute node before the
# annotation array. The cache is never downloaded from annotation workers.

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

scgate_db_path <- Sys.getenv("SCGATE_DB_PATH", unset = file.path(project_root, "aux", "scGateDB.rds"))
scgate_db_branch <- Sys.getenv("SCGATE_DB_BRANCH", unset = "41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4")
args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
force <- "--force" %in% args

validate_scgate_model_group <- function(group, label) {
  if (!is.list(group) || is.null(names(group)) || !length(group) ||
      any(!nzchar(names(group)))) return(FALSE)
  required <- c("levels", "use_as", "name", "signature")
  all(vapply(group, function(model) {
    is.data.frame(model) && nrow(model) > 0L && all(required %in% colnames(model))
  }, logical(1)))
}

validate_scgate_cache <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  isTRUE(tryCatch({
    cached <- readRDS(path)
    is.list(cached) && is.list(cached$human) &&
      validate_scgate_model_group(cached$human$PBMC, "human$PBMC") &&
      validate_scgate_model_group(cached$human$HiTME, "human$HiTME")
  }, error = function(e) FALSE))
}

if (validate_only) {
  if (!validate_scgate_cache(scgate_db_path)) {
    stop("scGate DB cache failed validation (requires non-empty human$PBMC and human$HiTME): ", scgate_db_path)
  }
  message("scGate DB cache validation passed: ", scgate_db_path)
  quit(save = "no", status = 0)
}

if (file.exists(scgate_db_path) && !force) {
  if (validate_scgate_cache(scgate_db_path)) {
    message("scGate DB cache already exists and passed schema validation: ", scgate_db_path)
    quit(save = "no", status = 0)
  }
  stop("scGate DB cache exists but is invalid; rerun with --force to rebuild: ", scgate_db_path)
}

library(scGate)
message("Downloading scGate models DB (one-time, force_update=TRUE)...")
scGate_models_DB <- get_scGateDB(branch = scgate_db_branch, force_update = TRUE)
if (!is.list(scGate_models_DB) || is.null(scGate_models_DB$human) ||
    !validate_scgate_model_group(scGate_models_DB$human$PBMC, "human$PBMC") ||
    !validate_scgate_model_group(scGate_models_DB$human$HiTME, "human$HiTME")) {
  stop("Downloaded scGate DB has invalid human$PBMC or human$HiTME model entries")
}
tmp_path <- paste0(scgate_db_path, ".tmp.", Sys.getpid())
dir.create(dirname(scgate_db_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(scGate_models_DB, tmp_path)
if (!file.exists(tmp_path) || file.info(tmp_path)$size <= 0 || !validate_scgate_cache(tmp_path)) {
  unlink(tmp_path)
  stop("scGate DB cache write failed validation: ", tmp_path)
}
if (!file.rename(tmp_path, scgate_db_path)) {
  unlink(tmp_path)
  stop("Could not atomically install scGate DB cache: ", scgate_db_path)
}
message("scGate DB cache written and validated: ", scgate_db_path)
