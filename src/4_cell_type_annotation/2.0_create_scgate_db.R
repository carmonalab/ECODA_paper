# ==============================================================================
# 2.0_create_scgate_db.R — One-time scGate model DB download / cache creation
# ==============================================================================
# Called by the canonical Stage 4 submitter on a compute node before the
# annotation array. The cache is never downloaded from annotation workers.

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

scgate_db_path <- Sys.getenv("SCGATE_DB_PATH", unset = file.path(project_root, "aux", "scGateDB.rds"))
scgate_db_branch <- Sys.getenv("SCGATE_DB_BRANCH", unset = "41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4")
scgate_model_cache_dir <- Sys.getenv(
  "SCGATE_MODEL_CACHE_DIR",
  unset = file.path(Sys.getenv("HOME_REF_DIR", unset = dirname(scgate_db_path)), "scGate_models")
)
scgate_ontology_branch <- Sys.getenv("SCGATE_ONTOLOGY_BRANCH", unset = "master")
if (!nzchar(scgate_model_cache_dir) || !grepl("^/", scgate_model_cache_dir)) {
  stop("SCGATE_MODEL_CACHE_DIR must be an absolute path: ", scgate_model_cache_dir)
}
if (!nzchar(scgate_ontology_branch) ||
    grepl("[/\\r\\n]", scgate_ontology_branch, perl = TRUE)) {
  stop("SCGATE_ONTOLOGY_BRANCH is invalid: ", scgate_ontology_branch)
}
args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
validate_db_only <- "--validate-db-only" %in% args
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

ontology_repo_path <- function(cache_root) {
  file.path(cache_root, paste0("scGate_models-", scgate_ontology_branch))
}

validate_ontology_cache <- function(cache_root = scgate_model_cache_dir) {
  repo_path <- ontology_repo_path(cache_root)
  if (!dir.exists(repo_path)) return(FALSE)
  dictionary_paths <- list.files(
    repo_path,
    pattern = "dictionary\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(dictionary_paths) != 1L) return(FALSE)
  dictionary_info <- file.info(dictionary_paths[[1]])
  if (is.na(dictionary_info$size) || dictionary_info$size <= 0) return(FALSE)
  isTRUE(tryCatch({
    dictionary <- read.table(
      dictionary_paths[[1]],
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    is.data.frame(dictionary) &&
      nrow(dictionary) > 0L &&
      all(c("scGate_multi", "CellOntology_ID") %in% colnames(dictionary))
  }, error = function(e) FALSE))
}

validate_all_cache <- function() {
  validate_scgate_cache(scgate_db_path) && validate_ontology_cache()
}

if (validate_only) {
  if (!validate_all_cache()) {
    stop(
      "scGate annotation cache failed validation (requires model DB and local ",
      "CellOntology dictionary): ", scgate_db_path, " / ", scgate_model_cache_dir
    )
  }
  message(
    "scGate annotation cache validation passed: ",
    scgate_db_path, " / ", scgate_model_cache_dir
  )
  quit(save = "no", status = 0)
}
if (validate_db_only) {
  if (!validate_scgate_cache(scgate_db_path)) {
    stop("scGate DB cache failed validation: ", scgate_db_path)
  }
  message("scGate DB cache validation passed: ", scgate_db_path)
  quit(save = "no", status = 0)
}

if (file.exists(scgate_db_path) && !force) {
  if (validate_all_cache()) {
    message("scGate annotation cache already exists and passed schema validation.")
    quit(save = "no", status = 0)
  }
  if (validate_scgate_cache(scgate_db_path)) {
    library(scGate)
    message("scGate DB is valid; staging missing CellOntology dictionary cache...")
  } else {
    stop("scGate DB cache exists but is invalid; rerun with --force to rebuild: ", scgate_db_path)
  }
} else {
  library(scGate)
}

stage_ontology_cache <- function() {
  cache_parent <- dirname(scgate_model_cache_dir)
  dir.create(cache_parent, showWarnings = FALSE, recursive = TRUE)
  staging_root <- tempfile(
    pattern = ".scGate_models.",
    tmpdir = cache_parent
  )
  dir.create(staging_root, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(staging_root, recursive = TRUE, force = TRUE), add = TRUE)
  suppressMessages(scGate:::get_CellOntology_dictionary(
    destination = staging_root,
    force_update = TRUE,
    branch = scgate_ontology_branch,
    verbose = FALSE
  ))
  staged_repo <- ontology_repo_path(staging_root)
  if (!validate_ontology_cache(staging_root)) {
    stop("Downloaded CellOntology dictionary cache failed validation: ", staged_repo)
  }
  dir.create(scgate_model_cache_dir, showWarnings = FALSE, recursive = TRUE)
  target_repo <- ontology_repo_path(scgate_model_cache_dir)
  backup_repo <- paste0(target_repo, ".old.", Sys.getpid())
  if (dir.exists(target_repo) &&
      !file.rename(target_repo, backup_repo)) {
    stop("Could not stage old CellOntology dictionary cache aside: ", target_repo)
  }
  if (!file.rename(staged_repo, target_repo)) {
    if (dir.exists(backup_repo)) file.rename(backup_repo, target_repo)
    stop("Could not atomically install CellOntology dictionary cache: ", target_repo)
  }
  if (dir.exists(backup_repo)) {
    unlink(backup_repo, recursive = TRUE, force = TRUE)
  }
}

stage_ontology_cache()
if (!force && validate_scgate_cache(scgate_db_path)) {
  if (!validate_all_cache()) {
    stop("CellOntology dictionary was staged but combined scGate cache validation failed.")
  }
  message("scGate DB retained; local CellOntology dictionary cache written and validated.")
  quit(save = "no", status = 0)
}

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
if (!file.exists(tmp_path) || file.info(tmp_path)$size <= 0 ||
    !validate_scgate_cache(tmp_path)) {
  unlink(tmp_path)
  stop("scGate DB cache write failed validation: ", tmp_path)
}
if (!file.rename(tmp_path, scgate_db_path)) {
  unlink(tmp_path)
  stop("Could not atomically install scGate DB cache: ", scgate_db_path)
}
if (!validate_all_cache()) {
  stop("Installed scGate DB and dictionary cache failed combined validation.")
}
message("scGate DB and local CellOntology dictionary cache written and validated.")
