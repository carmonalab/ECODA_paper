# ==============================================================================
# Artifact records and atomic binary publication
#
# Deep module for immutable Stage 5 artifact persistence. Callers need only the
# publication, validation, and checked-load functions below; path hashing,
# checksum sidecars, run ownership, and atomic replacement stay behind this
# interface.
#
# Public interface:
#   artifact_record_path(), artifact_validate_record(), artifact_write_record()
#   artifact_checksum_ok(), read_rds_checked(), read_feather_checked()
#   save_rds_atomic()
#
# Dependencies are base R, tools, and optional digest (with shasum/sha256sum as
# the standard-library fallback). Arrow is loaded only by read_feather_checked.
# ===============================================================================

.artifact_canonical_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) stop("artifact path must be one non-empty string")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.artifact_sha256_text <- function(value) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(tolower(digest::digest(value, algo = "sha256", serialize = FALSE)))
  }
  commands <- c("shasum", "sha256sum")
  executable <- commands[nzchar(Sys.which(commands))][1L]
  if (is.na(executable)) stop("digest package or SHA-256 utility is required")
  temporary <- tempfile("ecoda_sha256_")
  on.exit(unlink(temporary), add = TRUE)
  writeBin(charToRaw(value), temporary)
  arguments <- if (identical(executable, "shasum")) {
    c("-a", "256", temporary)
  } else {
    temporary
  }
  output <- system2(executable, arguments, stdout = TRUE, stderr = TRUE)
  digest <- sub("[[:space:]].*$", "", output[grepl("^[[:xdigit:]]{64}", output)][1L])
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[0-9a-fA-F]{64}$", digest)) {
    stop("could not compute SHA-256 for artifact record path")
  }
  tolower(digest)
}

.artifact_context <- function(producer = NULL, run_id = NULL,
                              execution_log = FALSE) {
  if (is.null(run_id)) run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  if (is.null(producer)) {
    env_name <- if (execution_log) {
      "ECODA_EXECUTION_LOG_PRODUCER"
    } else {
      "ECODA_ARTIFACT_PRODUCER"
    }
    producer <- Sys.getenv(env_name, unset = "")
    if (!nzchar(producer) && !execution_log) {
      producer <- Sys.getenv("METHOD", unset = "")
      if (!nzchar(producer)) producer <- Sys.getenv("ANALYSIS", unset = "")
    }
    if (execution_log && !nzchar(producer)) producer <- "stage5_execution_log"
  }
  if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) return(NULL)
  if (!is.character(producer) || length(producer) != 1L || is.na(producer) ||
      !nzchar(producer) || grepl("[\r\n]", producer, perl = TRUE)) return(NULL)
  runs_root <- Sys.getenv("ECODA_RUNS_ROOT", unset = "")
  if (!nzchar(runs_root)) {
    scratch <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
    if (nzchar(scratch)) runs_root <- file.path(scratch, "_ecoda_runs")
  }
  if (!nzchar(runs_root) || !grepl("^/", runs_root)) return(NULL)
  list(
    run_id = run_id,
    producer = producer,
    runs_root = runs_root
  )
}

artifact_record_path <- function(path, run_id, runs_root = NULL) {
  canonical <- .artifact_canonical_path(path)
  if (!is.character(run_id) || length(run_id) != 1L ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) {
    stop("artifact record run ID is invalid")
  }
  if (is.null(runs_root)) {
    context <- .artifact_context(run_id = run_id, producer = "record")
    if (is.null(context)) stop("artifact record root is unavailable")
    runs_root <- context$runs_root
  } else {
    if (!is.character(runs_root) || length(runs_root) != 1L ||
        is.na(runs_root) || !grepl("^/", runs_root)) {
      stop("artifact record root must be absolute")
    }
    runs_root <- as.character(runs_root)
  }
  key <- substr(.artifact_sha256_text(canonical), 1L, 32L)
  file.path(runs_root, run_id, "manifests", "artifacts",
            paste0(key, ".record"))
}

.artifact_sidecar <- function(path, verify = FALSE, allow_canonical = FALSE) {
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) return(NULL)
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || !isTRUE(file.info(sidecar)$size > 0)) return(NULL)
  lines <- tryCatch(readLines(sidecar, warn = FALSE),
                    error = function(e) character())
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    return(NULL)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded <- sub("^PATH=", "", lines[[3L]])
  canonical <- .artifact_canonical_path(path)
  valid_path <- identical(recorded, as.character(path)) ||
    (allow_canonical && identical(recorded, canonical))
  if (!valid_path || !grepl("^[0-9a-fA-F]{32}$", md5) ||
      !grepl("^[1-9][0-9]*$", size) ||
      !identical(size, as.character(file.info(path)$size))) return(NULL)
  if (verify) {
    actual <- unname(tools::md5sum(path))
    if (length(actual) != 1L || is.na(actual) ||
        !identical(tolower(md5), tolower(actual))) return(NULL)
  }
  list(MD5 = tolower(md5), SIZE = size, PATH = recorded)
}

# Atomic RDS write plus a sidecar used by idempotency checks.
save_rds_atomic <- function(object, file, producer = NULL, run_id = NULL) {
  validate_rds_object <- function(value) {
    if (is.null(value) || (is.list(value) && length(value) == 0L)) {
      stop("RDS artifact is empty: ", file)
    }
    dimensions <- dim(value)
    if (length(dimensions) > 0L && any(dimensions <= 0L)) {
      stop("RDS artifact has empty dimensions: ", file)
    }
  }
  validate_rds_object(object)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp.", Sys.getpid())
  checksum_tmp <- paste0(file, ".md5.tmp.", Sys.getpid())
  sidecar <- paste0(file, ".md5")
  backup <- paste0(file, ".previous.", Sys.getpid())
  sidecar_backup <- paste0(sidecar, ".previous.", Sys.getpid())
  had_file <- file.exists(file)
  had_sidecar <- file.exists(sidecar)
  installed <- FALSE
  sidecar_installed <- FALSE
  digest <- NULL
  size <- NULL
  restore <- function() {
    if (installed && file.exists(file)) unlink(file)
    if (had_file && file.exists(backup)) file.rename(backup, file)
    if (sidecar_installed && file.exists(sidecar)) unlink(sidecar)
    if (had_sidecar && file.exists(sidecar_backup)) file.rename(sidecar_backup, sidecar)
    if (!had_file && file.exists(file)) unlink(file)
    if (!had_sidecar && file.exists(sidecar)) unlink(sidecar)
  }
  on.exit({
    for (temporary in c(tmp, checksum_tmp, backup, sidecar_backup)) {
      if (file.exists(temporary)) unlink(temporary)
    }
  }, add = TRUE)
  tryCatch({
    saveRDS(object, tmp)
    if (!file.exists(tmp) || file.info(tmp)$size <= 0) {
      stop("Empty RDS temporary file: ", tmp)
    }
    if (had_file && !isTRUE(file.link(file, backup))) {
      stop("Could not preserve existing RDS artifact: ", file)
    }
    if (had_sidecar && !isTRUE(file.link(sidecar, sidecar_backup))) {
      stop("Could not preserve existing RDS checksum: ", file)
    }
    if (!file.rename(tmp, file)) stop("Could not atomically install RDS: ", file)
    installed <- TRUE
    digest <- tolower(unname(tools::md5sum(file)))
    size <- as.character(file.info(file)$size)
    writeLines(c(
      paste0("MD5=", digest),
      paste0("SIZE=", size),
      paste0("PATH=", file)
    ), checksum_tmp)
    if (!file.rename(checksum_tmp, sidecar)) {
      stop("Could not atomically install RDS checksum: ", file)
    }
    sidecar_installed <- TRUE
  }, error = function(error) {
    restore()
    stop(error)
  })
  artifact_write_record(file, producer = producer, run_id = run_id,
                        md5 = digest, size = size)
  invisible(NULL)
}

.artifact_record_from_file <- function(path, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (is.null(context)) return(NULL)
  record_path <- artifact_record_path(
    path, context$run_id, runs_root = context$runs_root
  )
  if (!file.exists(record_path) || !isTRUE(file.info(record_path)$size > 0)) {
    return(NULL)
  }
  lines <- tryCatch(readLines(record_path, warn = FALSE),
                    error = function(e) character())
  keys <- c("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
  if (length(lines) != length(keys) ||
      !identical(sub("=.*$", "", lines), keys)) {
    stop("Artifact record has the wrong schema: ", record_path)
  }
  values <- sub("^[^=]*=", "", lines)
  record <- as.list(values)
  names(record) <- keys
  canonical <- .artifact_canonical_path(path)
  if (!identical(record$PATH, canonical) ||
      !identical(record$RUN_ID, context$run_id) ||
      !identical(record$PRODUCER, context$producer) ||
      !identical(record$STATE, "PUBLISHED") ||
      !grepl("^[0-9a-f]{32}$", record$MD5) ||
      !grepl("^[1-9][0-9]*$", record$SIZE)) {
    stop("Artifact record binding is invalid: ", record_path)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0 ||
      !identical(record$SIZE, as.character(info$size))) {
    stop("Artifact record SIZE mismatch: ", path)
  }
  sidecar <- .artifact_sidecar(path, allow_canonical = TRUE)
  if (is.null(sidecar) || !identical(sidecar$MD5, record$MD5) ||
      !identical(sidecar$SIZE, record$SIZE)) {
    stop("Artifact record checksum sidecar mismatch: ", path)
  }
  record
}

artifact_validate_record <- function(path, producer = NULL, run_id = NULL) {
  .artifact_record_from_file(path, producer = producer, run_id = run_id)
}

artifact_record_for_load <- function(path, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (!is.null(context)) {
    record_path <- artifact_record_path(
      path, context$run_id, runs_root = context$runs_root
    )
    if (file.exists(record_path)) {
      record <- .artifact_record_from_file(
        path, producer = context$producer, run_id = context$run_id
      )
      if (is.null(record)) {
        stop("Artifact checksum validation failed: ", path)
      }
      sidecar <- .artifact_sidecar(path, verify = TRUE)
      if (is.null(sidecar) ||
          !identical(sidecar$MD5, record$MD5) ||
          !identical(sidecar$SIZE, record$SIZE)) {
        stop("Artifact checksum validation failed: ", path)
      }
      return(record)
    }
  }
  sidecar <- .artifact_sidecar(path, verify = TRUE)
  if (is.null(sidecar)) {
    stop("Artifact checksum validation failed: ", path)
  }
  sidecar
}

artifact_write_record <- function(
  path, producer = NULL, run_id = NULL, md5 = NULL, size = NULL
) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (is.null(context)) return(invisible(NULL))
  canonical <- .artifact_canonical_path(path)
  sidecar <- .artifact_sidecar(path, allow_canonical = TRUE)
  if (is.null(sidecar)) {
    stop("Cannot publish artifact record without a valid sidecar: ", path)
  }
  if (is.null(md5)) md5 <- sidecar$MD5
  if (is.null(size)) size <- sidecar$SIZE
  if (!identical(tolower(as.character(md5)), sidecar$MD5) ||
      !identical(as.character(size), sidecar$SIZE)) {
    stop("Artifact record checksum does not match sidecar: ", path)
  }
  record_path <- artifact_record_path(
    canonical, context$run_id, runs_root = context$runs_root
  )
  dir.create(dirname(record_path), showWarnings = FALSE, recursive = TRUE)
  temporary <- paste0(record_path, ".tmp.", Sys.getpid())
  writeLines(c(
    paste0("PATH=", canonical),
    paste0("SIZE=", sidecar$SIZE),
    paste0("MD5=", sidecar$MD5),
    paste0("RUN_ID=", context$run_id),
    paste0("PRODUCER=", context$producer),
    "STATE=PUBLISHED"
  ), temporary)
  if (!file.rename(temporary, record_path)) {
    if (file.exists(temporary)) unlink(temporary)
    stop("Could not atomically install artifact record: ", record_path)
  }
  invisible(list(
    PATH = canonical, SIZE = sidecar$SIZE, MD5 = sidecar$MD5,
    RUN_ID = context$run_id, PRODUCER = context$producer,
    STATE = "PUBLISHED"
  ))
}

artifact_checksum_ok <- function(file, producer = NULL, run_id = NULL) {
  context <- .artifact_context(producer = producer, run_id = run_id)
  if (!is.null(context)) {
    record_path <- artifact_record_path(
      file, context$run_id, runs_root = context$runs_root
    )
    if (file.exists(record_path)) {
      return(isTRUE(tryCatch(
        !is.null(.artifact_record_from_file(
          file, producer = context$producer, run_id = context$run_id
        )),
        error = function(e) FALSE
      )))
    }
  }
  isTRUE(!is.null(.artifact_sidecar(file, verify = TRUE)))
}

read_rds_checked <- function(path, producer = NULL, run_id = NULL) {
  artifact_record_for_load(path, producer = producer, run_id = run_id)
  readRDS(path)
}

read_feather_checked <- function(path, producer = NULL, run_id = NULL) {
  artifact_record_for_load(path, producer = producer, run_id = run_id)
  arrow::read_feather(path)
}
