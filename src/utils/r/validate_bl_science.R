# Finite BL10 evidence; successful science does not repair failed orchestration.
bl_require <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
bl_path <- function(path, empty = FALSE) {
  bl_require(is.character(path) && length(path) == 1L && startsWith(path, "/") &&
    file.exists(path) && !dir.exists(path) && identical(normalizePath(path), path) &&
    (empty || file.info(path)$size > 0), paste("Unsafe/missing evidence:", path))
  path
}
bl_sha <- function(path) digest::digest(file = bl_path(path, TRUE), algo = "sha256", serialize = FALSE)
bl_fields <- function(path) {
  lines <- readLines(bl_path(path), warn = FALSE); keys <- sub("=.*$", "", lines)
  bl_require(all(grepl("=", lines, fixed = TRUE)) && !anyDuplicated(keys), "Malformed proof fields")
  stats::setNames(sub("^[^=]*=", "", lines), keys)
}
bl_json <- function(path) {
  object <- jsonlite::fromJSON(bl_path(path), simplifyVector = FALSE)
  check <- function(x) {
    if (is.list(x)) {
      bl_require(is.null(names(x)) || !anyDuplicated(names(x)), "Duplicate JSON evidence key")
      lapply(x, check)
    }
    invisible(TRUE)
  }
  check(object); object
}
bl_proof <- function(proof, empty = FALSE) {
  path <- bl_path(proof$path, empty)
  bl_require(identical(bl_sha(path), proof$sha256), paste("Changed evidence:", path))
  path
}
bl_python <- function(code, args) {
  python <- Sys.getenv("PYTHON_BIN", "")
  bl_require(startsWith(python, "/") && file.exists(python), "Configured Python required")
  output <- system2(python, c("-c", shQuote(code), vapply(args, shQuote, character(1))), stdout = TRUE)
  bl_require(is.null(attr(output, "status")), "Pinned metadata helper failed")
  jsonlite::fromJSON(paste(output, collapse = "\n"), simplifyVector = FALSE)
}
bl_option <- function(argv, name) {
  values <- c(vapply(which(argv == name), function(i) argv[i + 1L], character(1)),
              sub(paste0("^", name, "="), "", argv[startsWith(argv, paste0(name, "="))]))
  bl_require(length(values) == 1L && !is.na(values), paste("Wrapper lacks exact", name))
  values
}
bl_rows <- function(path, columns) {
  rows <- strsplit(readLines(bl_path(path), warn = FALSE), "\t", fixed = TRUE)
  bl_require(length(rows) > 0L && all(lengths(rows) == columns) &&
    all(vapply(rows, function(x) all(nzchar(x)), logical(1))), "Malformed finite TSV")
  do.call(rbind, rows)
}
bl_record <- function(path, original, producer) {
  path <- bl_path(path); size <- as.character(file.info(path)$size); md5 <- unname(tools::md5sum(path))
  expected <- c(PATH = path, SIZE = size, MD5 = md5)
  sidecar <- bl_fields(paste0(path, ".md5"))
  bl_require(setequal(names(sidecar), names(expected)) && identical(unname(sidecar[names(expected)]), unname(expected)), "Artifact checksum mismatch")
  key <- substr(digest::digest(path, algo = "sha256", serialize = FALSE), 1, 32)
  record_path <- file.path(original, "manifests/artifacts", paste0(key, ".record")); record <- bl_fields(record_path)
  bl_require(setequal(names(record), c(names(expected), "RUN_ID", "PRODUCER", "STATE")) &&
    identical(unname(record[names(expected)]), unname(expected)) &&
    identical(unname(record[c("RUN_ID", "PRODUCER", "STATE")]), c(basename(original), producer, "PUBLISHED")), "Original producer record mismatch")
  list(path = path, size = size, md5 = md5, record = as.list(record), record_path = record_path,
       record_sha256 = bl_sha(record_path), sidecar_sha256 = bl_sha(paste0(path, ".md5")))
}
bl_owner <- function(path, original) {
  sha <- substr(digest::digest(path, algo = "sha256", serialize = FALSE), 1, 32)
  key <- paste0(chartr("/:,\t |", "______", basename(path)), "_", sha)
  if (nchar(key, type = "bytes") > 100L) key <- paste0("artifact_", sha)
  file <- file.path(dirname(dirname(original)), "_ecoda_owners/artifact", key, "owner"); owner <- bl_fields(file)
  common <- file.path(Sys.getenv("PROJECT_ROOT"), "src/utils/bash/ecoda_run_common.sh")
  checked <- system2("bash", c("-c", shQuote('HPC_SCRATCH_DIR="$4"; source "$1"; _ecoda_artifact_owner_validate_dir "$2" "$3"'),
    "bl-owner", shQuote(common), shQuote(dirname(file)), shQuote(path),
    shQuote(dirname(dirname(original)))), stdout = TRUE, stderr = TRUE)
  bl_require(is.null(attr(checked, "status")), paste("Production owner validator rejected original evidence:", paste(checked, collapse = "\n")))
  bl_require(identical(unname(owner[c("PATH", "RUN_ID", "STAGE")]), c(path, basename(original), "stage5")) &&
    unname(owner["STATE"]) %in% c("OK", "FAIL", "ACTIVE"), "Foreign or malformed original ownership")
  list(path = file, sha256 = bl_sha(file), fields = as.list(owner), policy = "Original state preserved; completed worker/producer evidence validates science only")
}

bl_failed_gate <- function(control, evidence, original) {
  gate <- bl_json(evidence$gate_manifest); inspected <- bl_json(evidence$inspect); profile <- bl_json(evidence$profile)
  bl_require(identical(gate$state, "FAILED") && identical(inspected$state, "FAILED") &&
    identical(inspected$operation, "inspect") && identical(gate$gate_id, inspected$gate_id) &&
    identical(gate$audit, inspected$audit) && identical(gate$audit$passed, FALSE), "Original failed gate/inspect must remain failed")
  bl_require(identical(gate$profile, profile$name) && identical(gate$profile_digest, bl_sha(evidence$profile)) &&
    identical(gate$command_digest, digest::digest(charToRaw(enc2utf8(gate$exact_command)), algo = "sha256", serialize = FALSE)), "Gate/profile/wrapper digest differs")
  argv <- unlist(bl_python("import json,shlex,sys;print(json.dumps(shlex.split(sys.argv[1])))", gate$exact_command), use.names = FALSE)
  bl_require(bl_option(argv, "--run-id") == basename(original) && bl_option(argv, "--scratch-root") == dirname(dirname(original)), "Wrong original wrapper scope")
  ids <- c(unlist(gate$scheduler_array_ids, use.names = FALSE), unlist(gate$scheduler_watchdog_ids, use.names = FALSE))
  schedulers <- bl_rows(evidence$scheduler_ids, 2L)
  bl_require(length(ids) == 5L && !anyDuplicated(ids) && all(grepl("^[1-9][0-9]*$", ids)) &&
    nrow(schedulers) == 5L && !anyDuplicated(schedulers[, 2]) && setequal(ids, schedulers[, 2]), "Exact five original scheduler IDs required")
  accounting <- inspected$audit$accounting
  template <- "sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode"
  bl_require(identical(profile$policy$accounting_command, template) && identical(accounting$attempted, TRUE) &&
    identical(accounting$query_count, 1L) && identical(accounting$actual_exit, 0L) &&
    identical(unlist(accounting$scheduler_ids, use.names = FALSE), ids) &&
    identical(accounting$command, sub("{scheduler_ids}", paste(ids, collapse = ","), template, fixed = TRUE)), "Original accounting query/scope differs")
  raw <- rawToChar(readBin(evidence$accounting_raw, "raw", n = file.info(evidence$accounting_raw)$size))
  bl_require(identical(raw, accounting$stdout) && endsWith(raw, "\n") && !grepl("truncated", raw, fixed = TRUE), "Full original accounting bytes required")
  parts <- strsplit(strsplit(raw, "\n", fixed = TRUE)[[1L]], "|", fixed = TRUE)
  bl_require(length(parts) == 23L && all(lengths(parts) == 3L), "All 23 accounting rows required")
  rows <- do.call(rbind, parts)
  bl_require(!anyDuplicated(rows[, 1]) && all(grepl("^[1-9][0-9]*$", rows[, 1])) &&
    all(ids %in% rows[, 1]) && length(accounting$rows) == 23L, "Malformed/missing actual accounting IDs")
  for (i in seq_len(23L)) {
    parsed <- accounting$rows[[i]]
    bl_require(setequal(names(parsed), c("line", "job_id", "state", "exit_code")) && parsed$line == i &&
      identical(unlist(parsed[c("job_id", "state", "exit_code")], use.names = FALSE), rows[i, ]), "Parsed accounting differs from raw bytes")
  }
  failures <- control$expected_orchestration_failures
  bl_require(length(failures) == 2L && setequal(vapply(failures, `[[`, character(1), "role"), c("WATCHDOG", "AGGREGATE")), "Two explicit orchestration failures required")
  fail_ids <- vapply(failures, `[[`, character(1), "id")
  bl_require(!anyDuplicated(fail_ids) && all(fail_ids %in% ids), "Unknown orchestration failure ID")
  success <- !rows[, 1] %in% fail_ids
  bl_require(sum(success) == 21L && all(rows[success, 2] == "COMPLETED") && all(rows[success, 3] == "0:0") &&
    all(rows[!success, 2] == "FAILED") && all(rows[!success, 3] == "1:0"), "Unexpected worker failure or orchestration outcome")
  bl_require(control$worker_array_id %in% ids && rows[match(control$worker_array_id, rows[, 1]), 2] == "COMPLETED", "Worker array did not complete")
  bl_require(control$worker_array_id %in% unlist(gate$scheduler_array_ids, use.names = FALSE) &&
    identical(unname(schedulers[match(control$worker_array_id, schedulers[, 2]), 1]), "ARRAY"), "Worker array scheduler kind differs")
  for (failure in failures) {
    expected_kind <- if (failure$role == "AGGREGATE") "AGGREGATE_GATE" else "WATCHDOG"
    bl_require(identical(unname(schedulers[match(failure$id, schedulers[, 2]), 1]), expected_kind), "Failure is not the declared orchestration scheduler kind")
    bl_proof(failure$stdout, TRUE); bl_proof(failure$stderr, TRUE)
    bl_require(startsWith(failure$stdout$path, paste0(original, "/logs/")) &&
      endsWith(failure$stdout$path, paste0("_", failure$id, ".log")) &&
      identical(failure$stderr$path, sub("[.]log$", ".err", failure$stdout$path)), "Failure log ID differs")
  }
  list(gate = gate, wrapper = argv, rows = rows, raw = raw, accounting = accounting, ids = ids, failures = failures)
}

bl_validate_cli <- function(project, original, config, identity, paths, audit) {
  listing <- tempfile("bl-artifacts-", tmpdir = audit); on.exit(unlink(c(listing, paste0(listing, ".md5"))))
  writeLines(paths, listing)
  writeLines(c(paste0("MD5=", unname(tools::md5sum(listing))), paste0("SIZE=", file.info(listing)$size),
    paste0("PATH=", listing)), paste0(listing, ".md5"))
  saved <- Sys.getenv(c("ECODA_RUN_ROOT", "ECODA_RUN_ID", "ECODA_ARTIFACT_BINDING", "ANALYSIS_ROOT", "ANALYSIS_NAS_ROOT"), unset = NA_character_)
  on.exit(for (key in names(saved)) if (is.na(saved[[key]])) Sys.unsetenv(key) else do.call(Sys.setenv, setNames(list(saved[[key]]), key)), add = TRUE)
  binding <- bl_json(file.path(original, "manifests/artifact_bindings.json"))
  Sys.setenv(ECODA_RUN_ROOT = original, ECODA_RUN_ID = basename(original),
    ECODA_ARTIFACT_BINDING = file.path(original, "manifests/artifact_bindings.json"),
    ANALYSIS_ROOT = file.path(original, "artifacts/benchmark"), ANALYSIS_NAS_ROOT = file.path(binding$nas_root, "benchmark"))
  cli <- file.path(project, "src/5_run_benchmark_methods/validate_benchmark_rds_contract.R")
  args <- c("--vanilla", cli, "--artifact-list", listing, "--config", config,
    "--input-root", dirname(dirname(original)), "--source-identity", identity, "--source-identity-verified")
  output <- system2(file.path(R.home("bin"), "Rscript"), vapply(args, shQuote, character(1)), stdout = TRUE, stderr = TRUE)
  bl_require(is.null(attr(output, "status")) && any(output == "benchmark RDS artifact-list contract OK"), paste("Strict RDS CLI rejected science:", paste(output, collapse = "\n")))
  list(path = cli, sha256 = bl_sha(cli), output = output)
}

validate_bl_science <- function(control_path, control_md5, report, project) {
  bl_require(identical(unname(tools::md5sum(bl_path(control_path))), control_md5), "Control checksum differs")
  control <- bl_json(control_path)
  bl_require(identical(control$schema, "ecoda_bl_science_control_v1"), "Unknown BL science control")
  original <- normalizePath(control$original_run_root, mustWork = TRUE); audit <- normalizePath(control$audit_root, mustWork = TRUE)
  bl_require(identical(original, control$original_run_root) && identical(audit, control$audit_root) &&
    basename(dirname(original)) == "_ecoda_runs" && dirname(audit) == dirname(original) && audit != original &&
    identical(report, file.path(audit, "artifacts/evidence/bl_science.json")) &&
    startsWith(control_path, paste0(audit, "/manifests/")) && !file.exists(report) && !file.exists(paste0(report, ".md5")), "Fresh isolated audit/report scope required")
  evidence <- lapply(control$evidence, bl_proof)
  required <- c("gate_manifest", "inspect", "profile", "accounting_raw", "scheduler_ids", "requested_selection", "worker_manifest", "source_manifest", "runtime_identity", "source_identity", "input_artifacts")
  bl_require(setequal(names(evidence), required), "Exact evidence set required")
  fixed <- c(scheduler_ids = "scheduler_ids.tsv", requested_selection = "requested_selection.tsv", source_manifest = "source.manifest", runtime_identity = "runtime.identity", source_identity = "source_identity.json", input_artifacts = "input_artifacts.tsv")
  for (name in names(fixed)) bl_require(identical(evidence[[name]], file.path(original, "manifests", fixed[[name]])), "Original proof location differs")
  validation_manifest <- bl_fields(bl_proof(control$validation_source_manifest))
  bl_require(identical(unname(validation_manifest["SOURCE_ROOT"]), project), "Validation source mismatch")
  terminal <- bl_failed_gate(control, evidence, original)
  source <- bl_fields(evidence$source_manifest); runtime <- bl_fields(evidence$runtime_identity); source_root <- unname(source["SOURCE_ROOT"])
  bl_require(bl_option(terminal$wrapper, "--source-root") == source_root &&
    bl_sha(bl_option(terminal$wrapper, "--source-manifest")) == bl_sha(evidence$source_manifest) &&
    unname(source["AUX_ROOT"]) == file.path(source_root, "aux") &&
    bl_option(terminal$wrapper, "--runtime-image") == unname(runtime["RUNTIME_IMAGE"]) &&
    bl_option(terminal$wrapper, "--runtime-manifest") == unname(runtime["RUNTIME_MANIFEST"]), "Original source/runtime/aux wrapper differs")
  config <- file.path(source_root, "datasets.json")
  bl_require(bl_sha(config) == unname(source["DATASETS_SHA256"]) &&
    bl_sha(unname(runtime["RUNTIME_MANIFEST"])) == unname(runtime["RUNTIME_MANIFEST_SHA256"]) &&
    as.character(file.info(unname(runtime["RUNTIME_IMAGE"]))$size) == unname(runtime["RUNTIME_IMAGE_SIZE"]), "Original config/runtime identity differs")
  image_manifest <- bl_fields(unname(runtime["RUNTIME_MANIFEST"]))
  bl_require(grepl("^[a-f0-9]{40}$", unname(source["SOURCE_COMMIT"])) &&
    grepl("^[a-f0-9]{64}$", unname(runtime["RUNTIME_IMAGE_SHA256"])) &&
    identical(unname(runtime["RUNTIME_IMAGE_SHA256"]), unname(image_manifest["IMAGE_SHA256"])) &&
    identical(unname(runtime["RUNTIME_IMAGE"]), unname(image_manifest["IMAGE_PATH"])), "Original versioned runtime/source identity differs")
  request <- bl_rows(evidence$requested_selection, 4L); matrix_rows <- bl_rows(evidence$worker_manifest, 4L)
  bl_require(nrow(request) == 10L && !anyDuplicated(request[, 1]) && all(request[, 2] == "benchmark_analysis") &&
    all(request[, 3] == "pseudobulk") && all(request[, 4] == "hvg2000_bl") &&
    nrow(matrix_rows) == 10L && !anyDuplicated(matrix_rows[, 1]) &&
    setequal(apply(matrix_rows, 1, paste, collapse = "\t"), apply(request, 1, paste, collapse = "\t")) &&
    startsWith(evidence$worker_manifest, paste0(original, "/manifests/")), "Exact BL10 physical request differs")
  bl_require(bl_sha(bl_option(terminal$wrapper, "--selection-file")) == bl_sha(evidence$requested_selection) &&
    bl_sha(bl_option(terminal$wrapper, "--input-artifacts")) == bl_sha(evidence$input_artifacts), "Wrapper input/selection differs from frozen run")
  request <- matrix_rows # Actual physical row order binds task indices.
  bl_require(length(control$workers) == 10L && length(control$inputs) == 10L &&
    identical(vapply(control$workers, `[[`, character(1), "dataset"), request[, 1]) &&
    identical(vapply(control$inputs, `[[`, character(1), "dataset"), request[, 1]), "Worker/input index coverage differs")
  identity <- bl_json(evidence$source_identity); refs <- bl_rows(evidence$input_artifacts, 6L)
  identity_proof <- bl_fields(paste0(evidence$source_identity, ".md5"))
  bl_require(identical(unname(identity_proof[c("PATH", "MD5", "SIZE")]), c(evidence$source_identity,
    unname(tools::md5sum(evidence$source_identity)), as.character(file.info(evidence$source_identity)$size))), "Original source identity checksum differs")
  bl_require(nrow(refs) == 10L && !anyDuplicated(refs[, 1]) && setequal(refs[, 1], request[, 1]) &&
    all(refs[, 2] == "benchmark_analysis"), "Original input scope differs")
  stat_code <- "import json,os,sys;s=os.stat(sys.argv[1]);print(json.dumps({k:str(getattr(s,k)) for k in ['st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns']}))"
  input_stats <- sample_ids <- input_proofs <- list()
  for (i in seq_len(10L)) {
    ds <- request[i, 1]; ref <- refs[refs[, 1] == ds, ]; observed <- bl_python(stat_code, ref[3])
    bl_require(identical(observed, control$inputs[[i]]$stat) && observed$st_size == ref[5], "Current input observation changed")
    input_stats[[ds]] <- observed
    proof <- bl_fields(ref[6]); record_path <- file.path(original, "manifests/artifacts", paste0(substr(digest::digest(ref[3], algo = "sha256", serialize = FALSE), 1, 32), ".record")); record <- bl_fields(record_path)
    bl_require(identical(unname(proof[c("PATH", "MD5", "SIZE")]), unname(ref[c(3, 4, 5)])) &&
      identical(unname(record[c("PATH", "MD5", "SIZE", "PRODUCER", "RUN_ID", "STATE")]), c(unname(ref[c(3, 4, 5)]), "stage5_preflight", basename(original), "PUBLISHED")), "Original input preflight proof differs")
    status_path <- bl_proof(control$inputs[[i]]$preflight_status); status <- bl_fields(status_path)
    bl_require(startsWith(status_path, paste0(original, "/status/")) &&
      identical(unname(status[c("STATE", "RUN_ID", "DATASET", "VIEW")]), c("OK", basename(original), ds, "benchmark_analysis")), "Original input preflight status differs")
    entry <- Filter(function(e) identical(e$dataset, ds) && identical(e$view, "benchmark_analysis"), identity$entries)
    bl_require(length(entry) == 1L && entry[[1]]$path == ref[3] && entry[[1]]$md5 == ref[4] && as.character(entry[[1]]$size) == ref[5], "Input source identity differs")
    sample_ids[[ds]] <- unlist(entry[[1]]$sample_ids, use.names = FALSE)
    input_proofs[[ds]] <- list(reference = as.list(ref), record_path = record_path, record_sha256 = bl_sha(record_path), external_proof_sha256 = bl_sha(ref[6]), status = control$inputs[[i]]$preflight_status)
  }
  helper <- new.env(parent = globalenv())
  sys.source(file.path(project, "src/utils/pseudobulk.R"), helper)
  sys.source(file.path(project, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"), helper)
  old_aux <- Sys.getenv("ECODA_AUX_ROOT", unset = NA_character_)
  on.exit(if (is.na(old_aux)) Sys.unsetenv("ECODA_AUX_ROOT") else Sys.setenv(ECODA_AUX_ROOT = old_aux), add = TRUE)
  Sys.setenv(ECODA_AUX_ROOT = file.path(source_root, "aux"))
  policy <- helper$pseudobulk_blacklist()$policy
  artifacts <- timings <- list(); artifact_rows <- character()
  for (i in seq_len(10L)) {
    ds <- request[i, 1]; worker <- control$workers[[i]]
    bl_require(worker$task_index == i, "Worker index differs from physical manifest")
    stdout <- bl_proof(worker$stdout, TRUE); bl_proof(worker$stderr, TRUE)
    marker <- paste0("ECODA_WORKER_DISPATCH METHOD=pseudobulk R_SCRIPT=", source_root,
      "/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R ANALYSIS_PASS= FORCE_BENCHMARK=0")
    bl_require(startsWith(stdout, paste0(original, "/logs/")) &&
      endsWith(stdout, paste0("_", control$worker_array_id, "_", i, ".log")) &&
      identical(worker$stderr$path, sub("[.]log$", ".err", stdout)) &&
      any(readLines(stdout, warn = FALSE) == marker), "Worker dispatch/log index identity differs")
    cache_path <- file.path(original, "artifacts/benchmark/pseudobulks", paste0(ds, "_pseudobulk_hvg2000_bl.rds"))
    result_path <- file.path(original, "artifacts/benchmark/results", paste0(ds, "_Pseudobulk_hvg2000_bl.rds"))
    cache_producer <- unname(helper$PB_VARIANT_PRODUCERS[["hvg2000_bl"]])
    cache_record <- bl_record(cache_path, original, cache_producer); result_record <- bl_record(result_path, original, "stage5_pseudobulk")
    cache <- readRDS(cache_path); result <- readRDS(result_path)
    helper$validate_pseudobulk_timing_record(cache, "hvg2000_bl")
    matrix <- as.matrix(cache$pb); features <- as.matrix(result$feat_mat)
    bl_require(identical(cache$blacklist_policy, policy) && identical(result$blacklist_policy, policy) &&
      identical(cache$blacklist_audit, result$blacklist_audit) && identical(cache$timing_id, result$blacklist_source_timing_id) &&
      identical(cache$timing_id, result$timing_id) &&
      identical(dimnames(matrix), dimnames(features)) && isTRUE(all.equal(matrix, features, tolerance = 0)) &&
      setequal(rownames(matrix), sample_ids[[ds]]) && !anyDuplicated(sample_ids[[ds]]), "Coupled BL feature/sample/policy/timing identity differs")
    filtered <- helper$filter_pseudobulk_blacklist(t(matrix), "human_aux_five_groups")
    bl_require(identical(dimnames(filtered$counts), dimnames(t(matrix))), "Blocked features remain in BL output")
    timing <- bl_proof(worker$timing); timing_record <- bl_record(timing, original, "stage5_execution_log")
    bl_require(identical(timing, file.path(original, "logs", paste0("execution_times_pseudobulk_", ds, "_hvg2000_bl.feather"))), "Timing belongs to another worker")
    table <- as.data.frame(arrow::read_feather(timing))
    bl_require(setequal(names(table), c("dataset", "method", "time_secs", "mem_GB")) && nrow(table) == 2L &&
      all(table$dataset == ds) && setequal(table$method, c("prepare_pseudobulk_hvg2000_bl", "Pseudobulk_hvg2000_bl")) &&
      all(is.finite(table$time_secs)) && all(table$time_secs >= 0) &&
      table$time_secs[table$method == "prepare_pseudobulk_hvg2000_bl"] == cache$variant_time_secs &&
      table$time_secs[table$method == "Pseudobulk_hvg2000_bl"] == as.numeric(result$exec_time), "Scientific timing differs from BL artifacts")
    bl_require(is.numeric(table$mem_GB) && !any(is.nan(table$mem_GB)) &&
      all(is.na(table$mem_GB) | (is.finite(table$mem_GB) & table$mem_GB >= 0)), "Timing memory is malformed")
    artifacts[[ds]] <- list(cache = cache_record, result = result_record,
      owners = list(bl_owner(cache_path, original), bl_owner(result_path, original)),
      n_samples = nrow(matrix), feature_ids = colnames(matrix), timing_id = cache$timing_id,
      blacklist_policy = policy, blacklist_audit = cache$blacklist_audit)
    timings[[ds]] <- list(proof = timing_record, rows = table)
    artifact_rows <- c(artifact_rows, paste(cache_path, "pseudobulk", ds, "benchmark_analysis", "0", sep = "\t"),
      paste(result_path, "pseudobulk", ds, "benchmark_analysis", "0", sep = "\t"))
  }
  validation <- bl_validate_cli(project, original, config, evidence$source_identity, artifact_rows, audit)
  for (i in seq_len(10L)) {
    ds <- request[i, 1]; ref <- refs[refs[, 1] == ds, ]
    bl_require(identical(bl_python(stat_code, ref[3]), input_stats[[ds]]), "Input observation changed during audit")
    for (kind in c("cache", "result")) {
      recorded <- artifacts[[ds]][[kind]]
      bl_require(identical(bl_record(recorded$path, original, recorded$record$PRODUCER), recorded), "Artifact changed during audit")
    }
    bl_require(identical(bl_owner(artifacts[[ds]]$cache$path, original), artifacts[[ds]]$owners[[1]]) &&
      identical(bl_owner(artifacts[[ds]]$result$path, original), artifacts[[ds]]$owners[[2]]), "Owner evidence changed during audit")
    bl_proof(control$workers[[i]]$stdout, TRUE); bl_proof(control$workers[[i]]$stderr, TRUE)
    bl_proof(control$workers[[i]]$timing)
    bl_require(bl_sha(input_proofs[[ds]]$record_path) == input_proofs[[ds]]$record_sha256 &&
      bl_sha(ref[6]) == input_proofs[[ds]]$external_proof_sha256, "Input proof changed during audit")
  }
  lapply(control$evidence, bl_proof)
  for (failure in terminal$failures) { bl_proof(failure$stdout, TRUE); bl_proof(failure$stderr, TRUE) }
  bl_require(identical(unname(tools::md5sum(file.path(source_root, "aux/genes.blocklist.rds"))), policy$source_md5), "Frozen auxiliary policy changed during audit")
  bl_require(identical(unname(tools::md5sum(control_path)), control_md5), "Control changed during audit")
  receipt <- list(schema = "ecoda_bl_science_receipt_v1", status = "NOOP_VALIDATED", science_validated = TRUE,
    publication_completed = FALSE, original_gate_status = "FAILED", original_run_root = original, audit_root = audit,
    control = list(path = control_path, md5 = control_md5), validation_source = control$validation_source_manifest,
    original_evidence = control$evidence, original_source = as.list(source), original_runtime = as.list(runtime),
    original_accounting = list(raw_path = evidence$accounting_raw, sha256 = bl_sha(evidence$accounting_raw), raw_stdout = terminal$raw,
      rows = terminal$accounting$rows, query_count = 1L, scheduler_ids = terminal$ids, orchestration_failures = terminal$failures,
      worker_array_id = control$worker_array_id, mapping = "No JobIDRaw-to-array-index mapping inferred"),
    input_policy = "Original completed full preflight/expected MD5 at fit time reused; current references/stat observed before and after audit. No claim these stats prove immutability since preflight; no fresh whole-H5AD hash, H5AD consumption, or model fit.",
    input_observations = input_stats, input_proofs = input_proofs, worker_evidence = control$workers,
    artifacts = artifacts, timings = timings, production_validator = validation)
  dir.create(dirname(report), recursive = TRUE, showWarnings = FALSE)
  bl_require(identical(normalizePath(dirname(report)), dirname(report)) && !file.exists(report), "Unsafe/existing report path")
  temporary <- tempfile(".bl-science-", tmpdir = dirname(report)); on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(receipt, temporary, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 16)
  bl_require(file.link(temporary, report), "Cannot create exclusive science receipt")
  invisible(receipt)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); bl_require(length(args) == 4L, "Expected CONTROL MD5 REPORT PROJECT_ROOT")
  validate_bl_science(args[1], args[2], args[3], args[4])
}
