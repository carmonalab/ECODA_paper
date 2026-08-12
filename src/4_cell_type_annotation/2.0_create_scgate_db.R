# ==============================================================================
# 2.0_create_scgate_db.R — One-time scGate model DB download / cache creation
# ==============================================================================
# Called by 2_submit_hpc_array.sh (srun compute session, pixi run --as-is
# Rscript --vanilla) BEFORE the annotation array, so that 2.1.1_process_chunk.R
# workers load the scGate model DB from a local cache instead of downloading it
# in parallel (up to MAX_NUM_CHUNKS_PARALLEL concurrent downloads).
# Writes to ${SCGATE_DB_PATH} (default: ${PROJECT_ROOT}/aux/scGateDB.rds).
# Idempotent: exits early if the cache already exists.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

scgate_db_path <- Sys.getenv("SCGATE_DB_PATH", unset = file.path(project_root, "aux", "scGateDB.rds"))
scgate_db_branch <- Sys.getenv("SCGATE_DB_BRANCH", unset = "41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4")

if (scgate_db_path != "" && file.exists(scgate_db_path)) {
  message("scGate DB cache already exists at: ", scgate_db_path)
  quit(save = "no", status = 0)
}

library(scGate)

message("Downloading scGate models DB (one-time, force_update=TRUE)...")
scGate_models_DB <- get_scGateDB(branch = scgate_db_branch, force_update = TRUE)

dir.create(dirname(scgate_db_path), showWarnings = FALSE, recursive = TRUE)
tmp_path <- paste0(scgate_db_path, ".tmp.", Sys.getpid())
saveRDS(scGate_models_DB, tmp_path)
file.rename(tmp_path, scgate_db_path)
message("scGate DB cache written to: ", scgate_db_path)
