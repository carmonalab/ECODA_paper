# Actual finite audit/strict artifact CLI. Synthetic compact matrices, no fits.
root <- normalizePath(".")
source("src/utils/r/validate_bl_science.R")
Sys.setenv(PYTHON_BIN = file.path(root, ".pixi/envs/default/bin/python"), PROJECT_ROOT = root)
reject <- function(expr) stopifnot(inherits(tryCatch(suppressWarnings(force(expr)), error = base::identity), "error"))
tmp <- tempfile("bl-science-"); dir.create(tmp); tmp <- normalizePath(tmp)
original <- file.path(tmp, "scratch/_ecoda_runs/original"); audit <- file.path(dirname(original), "audit")
for (run in c(original, audit)) dir.create(file.path(run, "manifests"), recursive = TRUE)
dir.create(file.path(original, "logs")); dir.create(file.path(original, "status/preflight"), recursive = TRUE)
source_root <- file.path(tmp, "source/tree"); dir.create(file.path(source_root, "aux"), recursive = TRUE)
invisible(file.copy(file.path(root, "datasets.json"), file.path(source_root, "datasets.json")))
invisible(file.copy(file.path(root, "aux/genes.blocklist.rds"), file.path(source_root, "aux/genes.blocklist.rds")))
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = 16)
proof <- function(path) list(path = path, sha256 = bl_sha(path))
seal <- function(path, producer) {
  fields <- c(PATH = path, SIZE = as.character(file.info(path)$size), MD5 = unname(tools::md5sum(path)))
  writeLines(paste0(c("MD5", "SIZE", "PATH"), "=", fields[c("MD5", "SIZE", "PATH")]), paste0(path, ".md5"))
  rec <- file.path(original, "manifests/artifacts", paste0(substr(digest::digest(path, algo = "sha256", serialize = FALSE), 1, 32), ".record"))
  dir.create(dirname(rec), recursive = TRUE, showWarnings = FALSE)
  fields <- c(fields, RUN_ID = "original", PRODUCER = producer, STATE = "PUBLISHED")
  writeLines(paste0(names(fields), "=", fields), rec)
}
owner <- function(path, state) {
  key <- paste0(basename(path), "_", substr(digest::digest(path, algo = "sha256", serialize = FALSE), 1, 32))
  p <- file.path(tmp, "scratch/_ecoda_owners/artifact", key, "owner")
  dir.create(dirname(p), recursive = TRUE)
  writeLines(c("RUN_ID=original", "STAGE=stage5", paste0("PATH=", path), paste0("STATE=", state), "PID=123"), p)
}
config <- jsonlite::fromJSON(file.path(source_root, "datasets.json"), simplifyVector = FALSE)
datasets <- names(Filter(function(x) isTRUE(x$use_for_benchmark), config))
datasets <- setdiff(datasets[!startsWith(datasets, "_")], "Stephenson"); stopifnot(length(datasets) == 10L)
request <- paste(datasets, "benchmark_analysis", "pseudobulk", "hvg2000_bl", sep = "\t")
requested <- file.path(original, "manifests/requested_selection.tsv"); writeLines(request, requested)
matrix_path <- file.path(original, "manifests/matrix.tsv"); writeLines(request, matrix_path)
selection <- file.path(original, "manifests/selection.tsv")
writeLines(paste(datasets, "benchmark_analysis", "pseudobulk", sep = "\t"), selection)
image <- file.path(tmp, "image.sif"); writeLines("fixture image", image)
image_manifest <- paste0(image, ".manifest"); writeLines(c(paste0("IMAGE_PATH=", image), paste0("IMAGE_SHA256=", bl_sha(image))), image_manifest)
runtime <- file.path(original, "manifests/runtime.identity")
writeLines(c(paste0("RUNTIME_IMAGE=", image), paste0("RUNTIME_IMAGE_SHA256=", bl_sha(image)), paste0("RUNTIME_MANIFEST=", image_manifest),
  paste0("RUNTIME_MANIFEST_SHA256=", bl_sha(image_manifest)), paste0("RUNTIME_IMAGE_SIZE=", file.info(image)$size)), runtime)
source_manifest <- file.path(original, "manifests/source.manifest")
writeLines(c(paste0("SOURCE_ROOT=", source_root), paste0("SOURCE_COMMIT=", strrep("a", 40)), paste0("AUX_ROOT=", source_root, "/aux"),
  paste0("DATASETS_SHA256=", bl_sha(file.path(source_root, "datasets.json")))), source_manifest)
validation_manifest <- file.path(audit, "manifests/validation.manifest"); writeLines(paste0("SOURCE_ROOT=", root), validation_manifest)
invisible(file.copy(validation_manifest, file.path(audit, "manifests/source.manifest")))
invisible(file.copy(runtime, file.path(audit, "manifests/runtime.identity")))
helper <- new.env(parent = globalenv()); sys.source("src/utils/pseudobulk.R", helper)
Sys.setenv(ECODA_AUX_ROOT = file.path(source_root, "aux"))
policy <- helper$pseudobulk_blacklist()$policy
samples <- paste0("s", 1:6); labels <- setNames(factor(rep(c("a", "b"), each = 3)), samples)
features <- matrix(seq_len(18), 6, dimnames = list(samples, c("retained_A", "retained_B", "retained_C")))
blocked <- unname(helper$pseudobulk_blacklist()$groups$TCR[[1]])
blacklist_audit <- helper$filter_pseudobulk_blacklist(matrix(1, 4, 1,
  dimnames = list(c(blocked, colnames(features)), "fixture")), "human_aux_five_groups")$audit
refs <- entries <- workers <- inputs <- list()
for (i in seq_along(datasets)) {
  ds <- datasets[i]; input <- file.path(tmp, "inputs", ds, config[[ds]]$views$benchmark_analysis$output_file_name)
  dir.create(dirname(input), recursive = TRUE); writeLines("No expression or H5AD read in this fixture", input); seal(input, "stage5_preflight")
  md5 <- unname(tools::md5sum(input)); size <- file.info(input)$size
  refs[[i]] <- paste(ds, "benchmark_analysis", input, md5, size, paste0(input, ".md5"), sep = "\t")
  entries[[i]] <- list(dataset = ds, view = "benchmark_analysis", path = input, md5 = md5, size = size, sample_ids = samples)
  status <- file.path(original, "status/preflight", paste0(ds, ".status"))
  writeLines(c("STATE=OK", "RUN_ID=original", paste0("DATASET=", ds), "VIEW=benchmark_analysis"), status)
  stat <- bl_python("import json,os,sys;s=os.stat(sys.argv[1]);print(json.dumps({k:str(getattr(s,k)) for k in ['st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns']}))", input)
  inputs[[i]] <- list(dataset = ds, stat = stat, preflight_status = proof(status))
  cache <- list(pb = features, time_secs = 2, mem_GB = 0.25, aggregate_time_secs = 0,
    shared_fit_time_secs = 0, shared_time_secs = 0, variant_time_secs = 2, shared_mem_GB = NA_real_,
    timing_id = paste("original", ds, "benchmark_analysis", "hvg2000_bl", sep = ":"), timing_schema = 2L,
    blacklist_policy = policy, blacklist_audit = blacklist_audit, blacklist_aggregate_time_secs = 1)
  result <- list(feat_mat = features, labels = labels, dist_mat = dist(features),
    scores = list(sil_score = 0.1, mod_knn3_score = 0.2, cluster_score = 0.3, anosim_score = 0.4, lisi_score = 0.5),
    exec_time = 3, mem_GB = 0.25, variant_time_secs = 3, shared_time_secs = 0, timing_schema = 2L,
    timing_id = cache$timing_id, shared_mem_GB = NA_real_, blacklist_policy = policy,
    blacklist_audit = blacklist_audit, blacklist_source_timing_id = cache$timing_id)
  cache_path <- file.path(original, "artifacts/benchmark/pseudobulks", paste0(ds, "_pseudobulk_hvg2000_bl.rds"))
  result_path <- file.path(original, "artifacts/benchmark/results", paste0(ds, "_Pseudobulk_hvg2000_bl.rds"))
  for (p in c(cache_path, result_path)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, cache_path); seal(cache_path, "stage5_prepare_pseudobulk_hvg2000_bl"); owner(cache_path, "ACTIVE")
  saveRDS(result, result_path); seal(result_path, "stage5_pseudobulk"); owner(result_path, "FAIL")
  stdout <- file.path(original, "logs", paste0("5_matrix_300_", i, ".log")); stderr <- sub("log$", "err", stdout)
  writeLines(paste0("ECODA_WORKER_DISPATCH METHOD=pseudobulk R_SCRIPT=", source_root,
    "/src/5_run_benchmark_methods/run_r_sample_embedding_methods/1.1.1_run_benchmark_methods_r.R ANALYSIS_PASS= FORCE_BENCHMARK=0"), stdout)
  file.create(stderr)
  timing <- file.path(original, "logs", paste0("execution_times_pseudobulk_", ds, "_hvg2000_bl.feather"))
  arrow::write_feather(data.frame(dataset = ds, method = c("prepare_pseudobulk_hvg2000_bl", "Pseudobulk_hvg2000_bl"), time_secs = c(2, 3), mem_GB = 0.25), timing)
  seal(timing, "stage5_execution_log")
  workers[[i]] <- list(dataset = ds, task_index = i, stdout = proof(stdout), stderr = proof(stderr), timing = proof(timing))
}
refs_path <- file.path(original, "manifests/input_artifacts.tsv"); writeLines(unlist(refs), refs_path)
identity <- file.path(original, "manifests/source_identity.json"); write_json(list(schema = 1L, entries = entries), identity); seal(identity, "stage5_source_identity")
nas <- file.path(tmp, "nas"); dir.create(nas)
external_refs <- file.path(tmp, "input_artifacts.tsv"); stopifnot(file.rename(refs_path, external_refs))
invisible(bl_python(paste0("import json,pathlib,sys;sys.path.insert(0,sys.argv[1]);from derivative_artifacts import bind;",
  "print(json.dumps(str(bind(run_root=pathlib.Path(sys.argv[2]),release_id='fixture',nas_root=pathlib.Path(sys.argv[3]),",
  "stage='stage5',config_path=pathlib.Path(sys.argv[4]),selection=pathlib.Path(sys.argv[5]),input_artifacts=pathlib.Path(sys.argv[6])))))"),
  c(file.path(root, "src/utils/py"), original, nas, file.path(source_root, "datasets.json"), selection, external_refs)))
schedulers <- file.path(original, "manifests/scheduler_ids.tsv")
writeLines(c("PREFLIGHT\t100", "R_ENV_PREFLIGHT\t200", "ARRAY\t300", "WATCHDOG\t400", "AGGREGATE_GATE\t500"), schedulers)
raw_file <- file.path(audit, "manifests/accounting.raw.tsv")
success_ids <- as.character(c(100:109, 200, 300:309))
raw <- paste0(paste(c(paste(success_ids, "COMPLETED", "0:0", sep = "|"), "400|FAILED|1:0", "500|FAILED|1:0"), collapse = "\n"), "\n")
profile <- list(name = "fixture", policy = list(accounting_command = "sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode"))
profile_file <- file.path(audit, "manifests/profile.json"); write_json(profile, profile_file)
argv <- c("executor", "--run-id", "original", "--scratch-root", file.path(tmp, "scratch"),
 "--source-root", source_root, "--source-manifest", source_manifest, "--runtime-image", image,
 "--runtime-manifest", image_manifest, "--selection-file", requested, "--input-artifacts", refs_path)
command <- paste(vapply(argv, shQuote, character(1)), collapse = " ")
gate <- list(state = "FAILED", gate_id = "original-gate", profile = "fixture", profile_digest = bl_sha(profile_file),
 exact_command = command, command_digest = digest::digest(charToRaw(command), algo = "sha256", serialize = FALSE),
 scheduler_array_ids = c("100", "200", "300"), scheduler_watchdog_ids = c("400", "500"))
gate_file <- file.path(audit, "manifests/gate.json"); inspect_file <- file.path(audit, "manifests/inspect.json")
write_terminal <- function(raw) {
 writeChar(raw, raw_file, eos = NULL, useBytes = TRUE)
 parts <- strsplit(strsplit(raw, "\n", fixed = TRUE)[[1]], "|", fixed = TRUE)
 accounting <- list(attempted = TRUE, query_count = 1L, actual_exit = 0L, scheduler_ids = c("100", "200", "300", "400", "500"),
  command = "sacct -n -P -X -j 100,200,300,400,500 --format=JobIDRaw,State,ExitCode", stdout = raw,
  rows = lapply(seq_along(parts), function(i) list(line = i, job_id = parts[[i]][1], state = parts[[i]][2], exit_code = parts[[i]][3])))
 gate$audit <- list(passed = FALSE, accounting = accounting)
 write_json(gate, gate_file)
 write_json(list(state = "FAILED", operation = "inspect", gate_id = "original-gate", audit = gate$audit), inspect_file)
}
write_terminal(raw)
failures <- lapply(c("400", "500"), function(id) {
 stdout <- file.path(original, "logs", paste0("orchestration_", id, ".log")); stderr <- sub("log$", "err", stdout)
 writeLines("Original failed orchestration preserved", stdout); file.create(stderr)
 list(role = if (id == "400") "WATCHDOG" else "AGGREGATE", id = id, stdout = proof(stdout), stderr = proof(stderr))
})
control <- list(schema = "ecoda_bl_science_control_v1", audit_root = audit, original_run_root = original,
 validation_source_manifest = proof(validation_manifest), worker_array_id = "300", workers = workers, inputs = inputs,
 expected_orchestration_failures = failures)
evidence_paths <- list(gate_manifest = gate_file, inspect = inspect_file, profile = profile_file,
 accounting_raw = raw_file, scheduler_ids = schedulers, requested_selection = requested, worker_manifest = matrix_path,
 source_manifest = source_manifest, runtime_identity = runtime, source_identity = identity, input_artifacts = refs_path)
control$evidence <- lapply(evidence_paths, proof)
control_file <- file.path(audit, "manifests/control.json"); report <- file.path(audit, "artifacts/evidence/bl_science.json")
run <- function(ctl = control) {
 write_json(ctl, control_file)
 validate_bl_science(control_file, unname(tools::md5sum(control_file)), report, root)
}
Sys.setenv(RETICULATE_PYTHON = Sys.getenv("PYTHON_BIN"), RETICULATE_USE_MANAGED_VENV = "no")
science_paths <- unlist(lapply(datasets, function(ds) c(
 file.path(original, "artifacts/benchmark/pseudobulks", paste0(ds, "_pseudobulk_hvg2000_bl.rds")),
 file.path(original, "artifacts/benchmark/results", paste0(ds, "_Pseudobulk_hvg2000_bl.rds")))))
before <- tools::md5sum(c(unlist(evidence_paths), science_paths,
 vapply(workers, function(x) x$timing$path, character(1))))
before_modes <- file.info(names(before))$mode
receipt <- run()
stopifnot(receipt$status == "NOOP_VALIDATED", receipt$original_gate_status == "FAILED", !receipt$publication_completed,
 length(receipt$artifacts) == 10L, length(receipt$timings) == 10L, length(receipt$original_accounting$rows) == 23L,
 identical(receipt$original_accounting$raw_stdout, raw),
 receipt$artifacts[[1]]$owners[[1]]$fields$STATE == "ACTIVE", receipt$artifacts[[1]]$owners[[2]]$fields$STATE == "FAIL",
 identical(receipt$validation_source$path, file.path(audit, "manifests/source.manifest")),
 identical(receipt$validation_runtime$path, file.path(audit, "manifests/runtime.identity")),
 identical(tools::md5sum(names(before)), before), identical(file.info(names(before))$mode, before_modes))
reject(run()); unlink(report)
# An unmatched failed raw child is never filtered into a success projection.
write_terminal(sub("309|COMPLETED|0:0", "309|FAILED|1:0", raw, fixed = TRUE))
bad <- control; bad$evidence <- lapply(evidence_paths, proof); reject(run(bad)); write_terminal(raw)
control$evidence <- lapply(evidence_paths, proof)
bad <- control; bad$workers[[1]]$task_index <- 2L; reject(run(bad))
bad <- control; bad$workers <- bad$workers[-1]; reject(run(bad))
bad <- control; bad$evidence$accounting_raw$sha256 <- strrep("0", 64); reject(run(bad))
bad <- control; bad$inputs[[1]]$stat$st_size <- "1"; reject(run(bad))
bad_gate <- bl_json(gate_file); bad_gate$audit$accounting$query_count <- 2L
write_json(bad_gate, gate_file)
write_json(list(state = "FAILED", operation = "inspect", gate_id = "original-gate", audit = bad_gate$audit), inspect_file)
bad <- control; bad$evidence <- lapply(evidence_paths, proof); reject(run(bad)); write_terminal(raw)
control$evidence <- lapply(evidence_paths, proof)
owner_file <- receipt$artifacts[[1]]$owners[[1]]$path; owner_fields <- readLines(owner_file)
writeLines(sub("RUN_ID=original", "RUN_ID=foreign", owner_fields, fixed = TRUE), owner_file)
reject(run()); writeLines(owner_fields, owner_file)
first <- file.path(original, "artifacts/benchmark/results", paste0(datasets[1], "_Pseudobulk_hvg2000_bl.rds"))
saved <- readRDS(first); changed <- saved; changed$feat_mat[1, 1] <- changed$feat_mat[1, 1] + 1
saveRDS(changed, first); seal(first, "stage5_pseudobulk"); reject(run())
saveRDS(saved, first); seal(first, "stage5_pseudobulk")
saved_timing <- arrow::read_feather(workers[[1]]$timing$path)
changed <- saved_timing; changed$time_secs[1] <- 999
arrow::write_feather(changed, workers[[1]]$timing$path); seal(workers[[1]]$timing$path, "stage5_execution_log")
bad <- control; bad$workers[[1]]$timing <- proof(workers[[1]]$timing$path); reject(run(bad))
arrow::write_feather(saved_timing, workers[[1]]$timing$path); seal(workers[[1]]$timing$path, "stage5_execution_log")
stopifnot(!file.exists(report))
unlink(tmp, recursive = TRUE)
cat("PASS: finite BL10 validator-only receipt, actual20-artifact CLI, all23rows retained, oldfailure/owners untouched; no fits or scheduler calls\n")
