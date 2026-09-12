# ==============================================================================
# 1.1.1_run_benchmark_methods_r.R — Run one R benchmark method (gloscope,
# mofa, pseudobulk, scitd or composition) for one dataset (Pipeline A).
#
# Called by 1.1_run_worker.sh via ${PIXI_RSCRIPT} with:
#   --config_path --ds_name --view benchmark_analysis --method {gloscope,mofa,
#   pseudobulk,scitd,composition} --input_dir --results_dir --pseudobulk_dir
#   --gloscope_cache_dir --log_file [--force]
# Canonical ordinary and batch pseudobulk paths use the raw H5AD CSR Sample
# aggregate plus the direct matrix DESeq2 API.  They never materialize a
# sample-level Seurat object.  Count-backed Seurat remains reserved for
# genuine cell-level methods such as scITD; CT pseudobulk calls
# process_pseudobulk_ct_h5ad_fig() directly.
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") {
  stop("PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")
}

source(file.path(project_root, "src/utils/imports_worker_core.R"))
source(file.path(project_root, "src/utils/load_worker_functions.R"))
source(file.path(project_root, "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
# Cache records are run-owned and optional for legacy artifacts.  A record
# permits a size/sidecar comparison before the read boundary; the full MD5
# check below is still mandatory immediately before every readRDS().
ecoda_local_validate_run_id <- function(run_id, label = "artifact producer run ID") {
  if (!is.character(run_id) || length(run_id) != 1L ||
      is.na(run_id) || !nzchar(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id, perl = TRUE)) {
    stop(label, " is invalid")
  }
  run_id
}

ecoda_local_current_run_id <- function() {
  run_id <- Sys.getenv("ECODA_RUN_ID", unset = "")
  if (!nzchar(run_id)) {
    run_root <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
    if (nzchar(run_root)) {
      run_id <- basename(normalizePath(path.expand(run_root), mustWork = FALSE))
    }
  }
  if (!nzchar(run_id)) return(NULL)
  ecoda_local_validate_run_id(run_id, "ECODA_RUN_ID")
}

ecoda_local_record_run_id <- function(run_id = NULL) {
  if (!is.null(run_id)) {
    return(ecoda_local_validate_run_id(run_id))
  }
  explicit <- Sys.getenv("ECODA_ARTIFACT_PRODUCER_RUN_ID", unset = "")
  if (nzchar(explicit)) {
    return(ecoda_local_validate_run_id(
      explicit, "ECODA_ARTIFACT_PRODUCER_RUN_ID"
    ))
  }
  ecoda_local_current_run_id()
}

ecoda_local_runs_root <- function() {
  runs_root <- Sys.getenv("ECODA_RUNS_ROOT", unset = "")
  if (!nzchar(runs_root)) {
    run_root <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
    if (nzchar(run_root)) {
      runs_root <- dirname(path.expand(run_root))
    }
  }
  if (!nzchar(runs_root)) {
    scratch_root <- Sys.getenv("HPC_SCRATCH_DIR", unset = "")
    if (nzchar(scratch_root)) {
      runs_root <- file.path(path.expand(scratch_root), "_ecoda_runs")
    }
  }
  if (!nzchar(runs_root) || !grepl("^/", path.expand(runs_root))) {
    return(NULL)
  }
  path.expand(runs_root)
}

ecoda_local_canonical_path <- function(path) {
  if (!is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(path)) {
    stop("artifact path must be one non-empty string")
  }
  normalizePath(path.expand(path), mustWork = FALSE)
}

ecoda_local_record_path <- function(path, run_id) {
  run_id <- ecoda_local_validate_run_id(run_id)
  runs_root <- ecoda_local_runs_root()
  if (is.null(runs_root)) return(NULL)
  canonical <- ecoda_local_canonical_path(path)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package is required for artifact-record paths")
  }
  key <- digest::digest(canonical, algo = "sha256", serialize = FALSE)
  if (!grepl("^[0-9a-f]{64}$", key, perl = TRUE)) {
    stop("could not derive artifact-record key for ", canonical)
  }
  file.path(
    runs_root, run_id, "manifests", "artifacts",
    paste0(substr(key, 1L, 32L), ".record")
  )
}

ecoda_local_checksum_fields <- function(path) {
  sidecar <- paste0(path, ".md5")
  if (!file.exists(path) || !file.exists(sidecar)) return(NULL)
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) return(NULL)
  lines <- tryCatch(readLines(sidecar, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) != 3L ||
      any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    return(NULL)
  }
  fields <- setNames(sub("^[^=]*=", "", lines), c("MD5", "SIZE", "PATH"))
  if (!grepl("^[0-9a-f]{32}$", fields[["MD5"]], perl = TRUE) ||
      !grepl("^[0-9]+$", fields[["SIZE"]], perl = TRUE) ||
      !nzchar(fields[["PATH"]])) {
    return(NULL)
  }
  fields
}

ecoda_local_strict_checksum <- function(path, description = "RDS artifact") {
  fields <- ecoda_local_checksum_fields(path)
  if (is.null(fields) ||
      !identical(fields[["PATH"]], path) ||
      !identical(fields[["SIZE"]], as.character(file.info(path)$size))) {
    return(FALSE)
  }
  actual <- tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  isTRUE(length(actual) == 1L && !is.na(actual) &&
         identical(fields[["MD5"]], actual))
}

ecoda_local_read_record <- function(path, run_id = NULL, producer = NULL) {
  run_id <- ecoda_local_record_run_id(run_id)
  if (is.null(run_id)) {
    return(list(state = "unavailable", record_path = NULL))
  }
  record_path <- ecoda_local_record_path(path, run_id)
  if (is.null(record_path) || !file.exists(record_path)) {
    return(list(state = "absent", record_path = record_path, run_id = run_id))
  }
  lines <- tryCatch(readLines(record_path, warn = FALSE), error = function(e) NULL)
  expected_keys <- c("PATH", "SIZE", "MD5", "RUN_ID", "PRODUCER", "STATE")
  if (is.null(lines) || length(lines) != length(expected_keys) ||
      any(!nzchar(lines)) ||
      !identical(sub("=.*$", "", lines), expected_keys)) {
    return(list(state = "invalid", record_path = record_path, run_id = run_id))
  }
  record <- setNames(sub("^[^=]*=", "", lines), expected_keys)
  canonical <- ecoda_local_canonical_path(path)
  valid <- identical(record[["PATH"]], canonical) &&
    grepl("^[0-9]+$", record[["SIZE"]], perl = TRUE) &&
    grepl("^[0-9a-f]{32}$", record[["MD5"]], perl = TRUE) &&
    identical(record[["RUN_ID"]], run_id) &&
    nzchar(record[["PRODUCER"]]) &&
    identical(record[["STATE"]], "PUBLISHED")
  if (!is.null(producer) &&
      (!is.character(producer) || length(producer) != 1L ||
       is.na(producer) || !identical(record[["PRODUCER"]], producer))) {
    valid <- FALSE
  }
  if (!valid) {
    return(list(state = "invalid", record_path = record_path, run_id = run_id))
  }
  list(
    state = "valid", record_path = record_path, run_id = run_id,
    canonical_path = canonical, record = record
  )
}

ecoda_local_record_only_valid <- function(path, record_info) {
  if (!identical(record_info[["state"]], "valid")) return(FALSE)
  fields <- ecoda_local_checksum_fields(path)
  if (is.null(fields)) return(FALSE)
  info <- file.info(path)
  record <- record_info[["record"]]
  isTRUE(
    !is.na(info$size) && info$size > 0 &&
      identical(record[["SIZE"]], as.character(info$size)) &&
      identical(record[["MD5"]], fields[["MD5"]]) &&
      identical(record[["SIZE"]], fields[["SIZE"]]) &&
      identical(fields[["PATH"]], path)
  )
}

ecoda_local_cache_valid <- function(path, producer = NULL) {
  record_info <- ecoda_local_read_record(path, producer = producer)
  if (identical(record_info[["state"]], "invalid") ||
      identical(record_info[["state"]], "absent")) {
    # A run-bound cache without its expected record is not reusable.  Legacy
    # artifacts are handled only when no run/producer record context exists.
    return(FALSE)
  }
  if (identical(record_info[["state"]], "valid")) {
    return(isTRUE(ecoda_local_record_only_valid(path, record_info)))
  }
  # No record context is the explicit legacy path: strict sidecar/content
  # checking, never existence-only reuse.
  isTRUE(ecoda_local_strict_checksum(path))
}

ecoda_local_read_rds <- function(path, description = "RDS artifact",
                                 producer = NULL) {
  record_info <- ecoda_local_read_record(path, producer = producer)
  if (identical(record_info[["state"]], "invalid") ||
      identical(record_info[["state"]], "absent")) {
    stop(description, " artifact record is missing or malformed: ", path)
  }
  if (identical(record_info[["state"]], "valid") &&
      !ecoda_local_record_only_valid(path, record_info)) {
    stop(description, " artifact record does not match current file: ", path)
  }
  # This full checksum is intentionally the final check before deserialization.
  if (!ecoda_local_strict_checksum(path, description)) {
    stop(description, " checksum/content validation failed: ", path)
  }
  readRDS(path)
}


ecoda_local_write_record <- function(
  path, producer, run_id, checksum = NULL
) {
  run_id <- ecoda_local_validate_run_id(run_id)
  if (!is.character(producer) || length(producer) != 1L ||
      is.na(producer) || !nzchar(producer) ||
      grepl("[\r\n]", producer, perl = TRUE)) {
    stop("artifact producer must be one non-empty line")
  }
  if (is.null(checksum)) {
    if (!ecoda_local_strict_checksum(path)) {
      stop("cannot publish an artifact without a strict checksum: ", path)
    }
    checksum <- ecoda_local_checksum_fields(path)
  }
  if (is.null(checksum)) stop("artifact checksum fields are unavailable: ", path)
  record_path <- ecoda_local_record_path(path, run_id)
  if (is.null(record_path)) {
    stop("artifact-record root is unavailable for ", path)
  }
  dir.create(dirname(record_path), showWarnings = FALSE, recursive = TRUE)
  temporary <- paste0(record_path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(c(
    paste0("PATH=", ecoda_local_canonical_path(path)),
    paste0("SIZE=", checksum[["SIZE"]]),
    paste0("MD5=", checksum[["MD5"]]),
    paste0("RUN_ID=", run_id),
    paste0("PRODUCER=", producer),
    "STATE=PUBLISHED"
  ), temporary, useBytes = TRUE)
  if (!file.rename(temporary, record_path)) {
    stop("could not atomically publish artifact record: ", record_path)
  }
  invisible(record_path)
}

ecoda_local_publish_rds <- function(object, path, producer) {
  save_rds_atomic(object, path)
  run_id <- ecoda_local_current_run_id()
  if (!is.null(run_id)) {
    if (!ecoda_local_strict_checksum(path)) {
      stop("new RDS artifact failed strict publication checksum: ", path)
    }
    checksum <- ecoda_local_checksum_fields(path)
    ecoda_local_write_record(path, producer, run_id, checksum = checksum)
  }
  invisible(NULL)
}


raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_flags(raw_args)

for (req in c("config_path", "ds_name", "view", "method", "input_dir",
              "results_dir", "log_file")) {
  if (is.null(args[[req]]) || identical(args[[req]], TRUE)) {
    stop("Missing required --", req, " argument")
  }
}
force <- isTRUE(args[["force"]]) || identical(args[["force"]], "TRUE")

method <- args$method
script_identity_arg <- commandArgs(trailingOnly = FALSE)
script_identity <- sub(
  "^--file=",
  "",
  script_identity_arg[grepl("^--file=", script_identity_arg)][1L]
)
message(
  "ECODA_R_DISPATCH_PARSED method=", method,
  " source_root=", Sys.getenv("ECODA_SOURCE_ROOT", unset = ""),
  " script=", script_identity
)
if (!method %in% c("gloscope", "mofa", "pseudobulk", "scitd",
                   "composition")) {
  stop("Unknown method '", method,
       "' (expected gloscope, mofa, pseudobulk, scitd or composition)")
}
combo_supplied <- !is.null(args[["combo"]])
combo_token <- args[["combo"]]
if (combo_supplied && method != "gloscope") {
  stop("--combo is only supported for method gloscope")
}

# A combo shard is intentionally ordinary-only. Batch-effect GloScope keeps
# its existing single-combo behavior and pass-qualified artifact names.

# Method-specific attaches: MOFA2/scITD are needed only by their methods
# (bare create_mofa / initialize_params + make_new_container); gloscope needs
# only the installed namespace (GloScope::gloscope is called qualified).
# composition calls the EPIC::EPIC + GloScope::gloscopeProp drivers BARE
# (in benchmark_methods_r.R; EPIC/GloScope are not attached by any loader),
# so both must be attached for it.
if (method == "mofa") library(MOFA2)
if (method == "scitd") library(scITD)
if (method == "composition") {
  library(EPIC)
  library(GloScope)
}

config <- read_datasets_json(args$config_path, view = args$view)
ds <- args$ds_name
analysis_pass <- args[["analysis_pass"]]
if (!is.null(analysis_pass) && !analysis_pass %in% c("uncorrected", "corrected")) {
  stop("Unknown analysis pass: ", analysis_pass)
}
if (combo_supplied && !is.null(analysis_pass)) {
  stop("--combo is only supported for ordinary GloScope runs")
}
cache_stem <- if (is.null(analysis_pass)) {
  ds
} else {
  paste0(ds, "_batch_effect_", analysis_pass)
}
entry <- config[[ds]]
if (is.null(entry)) {
  stop("Dataset '", ds, "' not found in ", args$config_path)
}

h5ad_path <- get_h5ad_path(config, ds, args$view, args$input_dir)
if (!file.exists(h5ad_path)) {
  stop("Input h5ad not found: ", h5ad_path)
}
# Source identity is run-owned when available; the CT store manifest records
# this identity alongside its H5AD path and scheduler owner.
source_identity <- Sys.getenv("ECODA_SOURCE_IDENTITY", unset = "")
if (!nzchar(source_identity)) {
  run_root_for_identity <- Sys.getenv("ECODA_RUN_ROOT", unset = "")
  if (nzchar(run_root_for_identity)) {
    candidate_identity <- file.path(
      run_root_for_identity, "manifests", "source_identity.json"
    )
    if (file.exists(candidate_identity)) source_identity <- candidate_identity
  }
}
if (!nzchar(source_identity)) source_identity <- NULL
sample_col <- "Sample"
correct_batch_mode <- identical(analysis_pass, "corrected")
batch_keys <- NULL
batch_context <- NULL
python_batch_metadata <- NULL
batch_col <- NULL
h5ad_expected_batch_contract <- NULL
pseudobulk_batch_contract <- NULL
method_batch_contract <- NULL
if (correct_batch_mode) {
  if (is.null(entry$batch_col)) {
    stop("corrected batch-effect view requires a confirmed columns.batch")
  }
  batch_keys <- ecoda_batch_normalize_keys(
    entry$batch_col,
    sample_col = sample_col,
    biological_label = entry$label_col
  )
  batch_col <- if (length(batch_keys) >= 2L) {
    "__ecoda_batch_combined_v1"
  } else {
    batch_keys[[1L]]
  }
  h5ad_expected_batch_contract <- ecoda_hpc_batch_contract_identity(
    batch_keys,
    sample_col = sample_col,
    method_id = "preprocess",
    model_id = "hvg_composite_v1"
  )
  pseudobulk_batch_contract <- ecoda_hpc_batch_contract_identity(
    batch_keys,
    sample_col = sample_col,
    method_id = "Pseudobulk",
    model_id = "pseudobulk_composite_v1"
  )
  method_batch_contract <- switch(
    method,
    composition = ecoda_hpc_batch_contract_identity(
      batch_keys,
      sample_col = sample_col,
      method_id = "ECODA_authors_HR",
      model_id = "ecoda_additive_random_intercepts_v1"
    ),
    gloscope = ecoda_hpc_batch_contract_identity(
      batch_keys,
      sample_col = sample_col,
      method_id = "GloScope",
      model_id = "embedding_consumer_harmony_v1"
    ),
    pseudobulk = pseudobulk_batch_contract,
    NULL
  )
}
dir.create(args$results_dir, showWarnings = FALSE, recursive = TRUE)

method_rds_stem <- if (is.null(analysis_pass)) ds else cache_stem
method_rds <- file.path(
  args$results_dir,
  paste0(method_rds_stem, "_", method, ".rds")
)

# GloScope and composition use the embedding/obs-only loader.  MOFA and
# pseudobulk use a stricter metadata/HVG-only reader so complete cache paths
# never open or read count values.  scITD remains the sole canonical path
# below that materializes a count-backed Seurat object.
counts_free_method <- method %in% c("gloscope", "composition")
pseudobulk_metadata_method <- method %in% c("mofa", "pseudobulk")
dispatch_branch <- if (counts_free_method) {
  "counts_free"
} else if (pseudobulk_metadata_method) {
  "metadata"
} else if (method == "scitd") {
  "count_backed_seurat"
} else {
  "unknown"
}
message(
  "ECODA_R_DISPATCH_PRE method=", method,
  " branch=", dispatch_branch,
  " analysis_pass=", ifelse(is.null(analysis_pass), "", analysis_pass),
  " force=", force
)
embedding_key <- if (args$view == "batch_effect_corrected") {
  "X_pca_harmony_batch_effect_corrected_hvg2000"
} else if (args$view == "batch_effect_uncorrected") {
  "X_pca_batch_effect_uncorrected_hvg2000"
} else {
  "X_pca_benchmark_analysis_hvg2000"
}
  # Validate all cells before any loader is allowed to select a first row per
  # Sample.  The ordinary and uncorrected paths do not incur this pass.
  validation_method_id <- if (method == "composition") {
    "ECODA_authors_HR"
  } else if (method == "gloscope") {
    "GloScope"
  } else {
    "Pseudobulk"
  }
  validation_model_id <- if (method == "composition") {
    "ecoda_additive_random_intercepts_v1"
  } else if (method == "gloscope") {
    "embedding_consumer_harmony_v1"
  } else {
    "pseudobulk_composite_v1"
  }
  if (correct_batch_mode) {
    python_batch_metadata <- validate_h5ad_corrected_batch_metadata(
      h5ad_path = h5ad_path,
      batch_keys = as.list(unname(batch_keys)),
      sample_col = sample_col,
      biological_label = entry$label_col,
      method_id = validation_method_id,
      model_id = validation_model_id
    )
  }
required_hvg <- if (identical(args$view, "benchmark_analysis")) 3000L else 2000L
hvg_rank_genes <- NULL
embedding_matrices <- NULL
embedding_sample_ids <- NULL
if (pseudobulk_metadata_method) {
  ct_columns <- if (method == "pseudobulk" && is.null(analysis_pass)) {
    c(entry$cell_type_low_res, entry$cell_type_high_res)
  } else {
    character()
  }
  metadata_info <- load_h5ad_pseudobulk_metadata(
    h5ad_path,
    sample_col = sample_col,
    metadata_columns = unique(c(
      entry$label_col,
      if (correct_batch_mode) batch_keys else batch_col,
      ct_columns
    )),
    n_hvg = required_hvg,
    required_nonmissing_columns = unique(c(
      sample_col,
      entry$label_col,
      if (correct_batch_mode) batch_keys else batch_col
    )),
    expected_batch_contract = h5ad_expected_batch_contract,
    view = args$view,
    method = method
  )
  hvg_rank_genes <- metadata_info$hvg_rank_genes
} else if (method == "gloscope" || method == "composition") {
  # Both methods are counts-free, but GloScope has a deliberately minimal
  # metadata contract: it needs only Sample/label and its stored embeddings.
  # In particular, do not request the CT annotations used by composition and
  # never let GloScope fall through to the count-backed Seurat branch below.
  embedding_keys <- if (method == "gloscope" && is.null(analysis_pass)) {
    c(
      "X_pca_benchmark_analysis_hvg1000",
      "X_pca_benchmark_analysis_hvg2000",
      "X_pca_benchmark_analysis_hvg3000"
    )
  } else {
    embedding_key
  }
  composition_obs_columns <- if (
    method == "composition" &&
    is.null(analysis_pass) &&
    length(entry$not_suitable_for_auto_annotation) == 0
  ) {
    c("layer2", "scATOMIC_pred")
  } else {
    character()
  }
  obs_columns <- if (method == "gloscope") {
    c(
      sample_col,
      entry$label_col,
      if (correct_batch_mode) batch_keys else NULL
    )
  } else {
    c(
      sample_col,
      entry$label_col,
      entry$cell_type_low_res,
      entry$cell_type_high_res,
      if (correct_batch_mode) batch_keys else batch_col,
      composition_obs_columns
    )
  }
  adata <- load_h5ad_counts_free(
    h5ad_path,
    unique(obs_columns[!is.na(obs_columns) & nzchar(obs_columns)]),
    embedding_keys,
    obs_prefixes = if (method == "composition") "leiden_res_" else character(),
    view = args$view,
    method = method,
    expected_batch_contract = h5ad_expected_batch_contract
  )
  obs <- py_to_r(adata$obs)
  hvg_rank_genes <- get_hvg_rank_genes(adata)
  if (method == "gloscope") {
    embedding_names <- if (is.null(analysis_pass)) {
      c(
        hvg1000 = "X_pca_benchmark_analysis_hvg1000",
        hvg2000 = "X_pca_benchmark_analysis_hvg2000",
        hvg3000 = "X_pca_benchmark_analysis_hvg3000"
      )
    } else {
      c(hvg2000 = embedding_key)
    }
    embedding_matrices <- lapply(unname(embedding_names), function(key) {
      py_to_r(adata$obsm[[key]])
    })
    names(embedding_matrices) <- names(embedding_names)
    embedding_sample_ids <- as.character(obs[[sample_col]])
  }
} else {
  ad <- import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(h5ad_path, backed = "r")
  obs <- py_to_r(adata$obs)
  validate_benchmark_h5ad_contract(
    adata,
    obs = obs,
    view = args$view,
    method = method,
    expected_batch_contract = h5ad_expected_batch_contract
  )
  hvg_rank_genes <- get_hvg_rank_genes(adata)
}

if (!sample_col %in% colnames(obs)) {
  stop(sample_col, " not found in obs columns of ", h5ad_path)
}
blind_mode <- is.null(analysis_pass) || analysis_pass == "uncorrected"
if (correct_batch_mode) {
  missing_batch_keys <- setdiff(batch_keys, colnames(obs))
  if (length(missing_batch_keys) > 0L) {
    stop(
      "Confirmed batch column(s) missing from obs of ", h5ad_path, ": ",
      paste(missing_batch_keys, collapse = ", ")
    )
  }
  batch_context <- ecoda_hpc_batch_context(
    metadata = obs,
    batch_keys = as.list(unname(batch_keys)),
    sample_col = sample_col,
    biological_label = entry$label_col,
    python_metadata = python_batch_metadata
  )
  batch_col <- batch_context$scalar_batch_col
}
if (correct_batch_mode) {
  pseudobulk_batch_contract <- ecoda_hpc_augment_batch_contract(
    identity = pseudobulk_batch_contract,
    validation = batch_context$validation,
    method_id = "Pseudobulk",
    batch_keys = batch_keys,
    scalar_batch_col = batch_context$scalar_batch_col
  )
  method_batch_contract <- switch(
    method,
    composition = ecoda_hpc_augment_batch_contract(
      identity = method_batch_contract,
      validation = batch_context$validation,
      method_id = "ECODA_authors_HR",
      batch_keys = batch_keys
    ),
    gloscope = ecoda_hpc_augment_batch_contract(
      identity = method_batch_contract,
      validation = batch_context$validation,
      method_id = "GloScope",
      batch_keys = batch_keys,
      scalar_batch_col = batch_context$scalar_batch_col
    ),
    pseudobulk = pseudobulk_batch_contract,
    NULL
  )
}

if (!combo_supplied && ecoda_local_cache_valid(method_rds) && !force) {
  message("Method results already exist and passed checksum/record validation: ", method_rds)
  cached <- ecoda_local_read_rds(method_rds, "Method results")
  if (!is.list(cached)) {
    stop("Method results artifact is not a list: ", method_rds)
  }
  if (correct_batch_mode) {
    ecoda_hpc_validate_batch_contract(
      cached[["batch_contract"]],
      method_batch_contract,
      label = paste0("Method results ", ds, "/", method)
    )
    if (method == "composition") {
      required_composition_methods <- c(
        "ECODA_authors_HR",
        "ECODA_authors_HR_NULL",
        "ECODA_seuratres_2"
      )
      missing_composition_methods <- setdiff(
        required_composition_methods,
        names(cached)
      )
      if (length(missing_composition_methods) > 0L) {
        stop(
          "Method results ", ds, "/", method,
          " is missing corrected composition results: ",
          paste(missing_composition_methods, collapse = ", ")
        )
      }
      for (composition_method in required_composition_methods) {
        composition_contract <- ecoda_batch_augment_contract(
          identity = ecoda_hpc_batch_contract_identity(
            batch_keys,
            sample_col = sample_col,
            method_id = composition_method,
            model_id = "ecoda_additive_random_intercepts_v1"
          ),
          validation = batch_context$validation,
          correction_mode = "additive_random_intercepts",
          correction_formula = if (length(batch_keys) == 1L) {
            "y ~ 1 + (1 | batch)"
          } else {
            paste0(
              "y ~ 1 + ",
              paste0(
                "(1 | batch_key_", seq_along(batch_keys), ")",
                collapse = " + "
              )
            )
          }
        )
        ecoda_hpc_validate_batch_contract(
          cached[[composition_method]][["batch_contract"]],
          composition_contract,
          label = paste0(
            "Method result ", ds, "/", method, "/", composition_method
          )
        )
      }
    }
  }
  shared_replayed <- character()
  shared_rows <- list()
  for (nm in setdiff(
    names(cached),
    c("batch_contract", "batch_contract_identity")
  )) {
    value <- cached[[nm]]
    validate_hpc_timing_bundle(
      value, label = paste0("Method result ", ds, "/", method, "/", nm)
    )
    if ("timing_schema" %in% names(value)) {
      timing_id <- as.character(value[["timing_id"]])
      shared_method <- if ("shared_timing_method" %in% names(value)) {
        as.character(value[["shared_timing_method"]])
      } else if (grepl("^Pseudobulk_CT_", nm)) {
        ct_shared_timing_method_from_timing_id(
          timing_id, fallback_method = nm
        )
      } else {
        "prepare_pseudobulk_shared"
      }
      shared_key <- paste(shared_method, timing_id, sep = "\r")
      if (!shared_key %in% shared_replayed) {
        shared_rows[[shared_key]] <- list(
          method = shared_method,
          time = as.numeric(value[["shared_time_secs"]]),
          mem = value[["shared_mem_GB"]]
        )
        shared_replayed <- c(shared_replayed, shared_key)
      }
    }
    if (!is.null(value$exec_time)) {
      log_exec_row(ds, nm, value$exec_time, args$log_file,
                   mem_gb = value$mem_GB)
    }
  }
  for (shared in shared_rows) {
    log_exec_row(
      ds, shared[["method"]], shared[["time"]], args$log_file,
      mem_gb = shared[["mem"]]
    )
  }
  quit(save = "no", status = 0)
}
pb_variants <- NULL
seurat <- NULL
metadata <- NULL
labels <- NULL
if (method == "gloscope") {
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
}

if (method %in% c("mofa", "pseudobulk")) {
  # Both methods consume direct matrix pseudobulks.  Cache validation occurs
  # inside load_pb_variants before the H5AD raw counts pass; a complete cache
  # set therefore performs no aggregation and creates no Seurat object.
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method ", method)
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  requested_pb_variants <- if (method == "mofa" || is.null(analysis_pass)) {
    PB_VARIANT_NAMES
  } else {
    "hvg2000"
  }
  pb_variants <- load_pb_variants(
    seurat = NULL,
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir,
    ds = ds,
    force = force,
    log_file = args$log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = requested_pb_variants,
    h5ad_path = h5ad_path,
    view = args$view,
    analysis_pass = analysis_pass,
    run_id = ecoda_local_current_run_id(),
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_contract = if (correct_batch_mode) pseudobulk_batch_contract else NULL,
    expected_h5ad_batch_contract = h5ad_expected_batch_contract
  )
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
} else if (method == "composition") {
  # Obs-only path: no Seurat materialization. Consumes the backed h5ad obs
  # (cell-level metadata), the hvg2000 obsm PCA embedding (Avg_PCA_embedding)
  # and the precomputed hvg2000 pseudobulk variant (ECODA_deconv; the submit
  # script auto-prepends prepare_pseudobulk for composition).
  if (is.null(args$pseudobulk_dir) || identical(args$pseudobulk_dir, TRUE)) {
    stop("Missing required --pseudobulk_dir argument for method composition")
  }
  dir.create(args$pseudobulk_dir, showWarnings = FALSE, recursive = TRUE)
  # Keep a configured Leiden annotation source column alongside its legacy
  # RNA_snn_res.* alias. Parkinson declares leiden_res_5_* as its author
  # high-resolution column; renaming it in place would invalidate that config.
  obs <- rename_leiden_cols(
    obs,
    view = args$view,
    preserve_source = TRUE
  )
  metadata <- collapse_sample_metadata(obs, sample_col = sample_col)
  if (correct_batch_mode) {
    metadata <- ecoda_hpc_apply_batch_context(
      metadata,
      batch_context,
      sample_col = sample_col
    )
  }
  labels <- as.factor(metadata[[entry$label_col]])
  names(labels) <- metadata[[sample_col]]
  obsm_keys <- py_to_r(import_builtins(convert = FALSE)$list(
    adata$obsm$keys()
  ))
  emb_key <- embedding_key
  if (!emb_key %in% obsm_keys) {
    stop("Embedding '", emb_key, "' not found in adata.obsm of ", h5ad_path,
         ". Re-run preprocessing (1.1.1_preprocess.py) for this dataset.")
  }
  pca_emb <- py_to_r(adata$obsm[[emb_key]])
  if (is.null(rownames(pca_emb))) rownames(pca_emb) <- rownames(obs)
  colnames(pca_emb) <- paste0("PC_", seq_len(ncol(pca_emb)))
  pb_variants <- load_composition_pb_variants(
    sample_col = sample_col,
    hvg_rank_genes = hvg_rank_genes,
    pseudobulk_dir = args$pseudobulk_dir,
    ds = ds,
    log_file = args$log_file,
    cache_stem = cache_stem,
    batch_col = batch_col,
    blind = blind_mode,
    correct_batch = correct_batch_mode,
    variants = if (!is.null(analysis_pass)) "hvg2000" else PB_VARIANT_NAMES,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_contract = if (correct_batch_mode) pseudobulk_batch_contract else NULL,
    expected_h5ad_batch_contract = h5ad_expected_batch_contract
  )
} else {
  # scITD is the genuine cell-level count consumer.  Keep its existing
  # count-backed Seurat boundary; ordinary pseudobulk and CT never enter this
  # branch.
  seurat <- load_benchmark_seurat(
    adata, obs, sample_col = sample_col,
    fetch_embedding = NULL,
    counts_layer = "counts"
  )
  # Sample names are already standardized in the preprocessed obs
  # (1.1.1_preprocess.py): no standardize_sample_names() re-application
  # (kept only in the legacy Seurat path of run_benchmark_analysis).
  if (length(hvg_rank_genes) > 0) {
    VariableFeatures(seurat) <- hvg_rank_genes[
      seq_len(min(2000, length(hvg_rank_genes)))
    ]
  }
  seurat@misc$label_col <- entry$label_col
  seurat@misc$cell_type_low_res <- entry$cell_type_low_res
  seurat@misc$cell_type_high_res <- entry$cell_type_high_res
  metadata <- get_metadata(seurat)
  labels <- get_labels(seurat, entry$label_col)
}
message(
  "ECODA_R_DISPATCH_FINAL method=", method,
  " branch=", dispatch_branch,
  " seurat_is_null=", is.null(seurat),
  " embedding_matrices_is_null=", is.null(embedding_matrices)
)

results <- switch(
  method,
  gloscope = run_gloscope_hpc(
    seurat, metadata, label_col = entry$label_col,
    sample_col = sample_col,
    gloscope_cache_dir = args$gloscope_cache_dir,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file,
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    combo_token = combo_token,
    embedding_name = if (!is.null(analysis_pass)) {
      sub("^X_", "", embedding_key)
    } else {
      NULL
    },
    embedding_matrices = embedding_matrices,
    embedding_sample_ids = embedding_sample_ids,
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL
  ),
  mofa = run_mofa_hpc(
    metadata, labels, pb_variants,
    results_dir = args$results_dir,
    ds = ds,
    force = force,
    log_file = args$log_file
  ),
  pseudobulk = run_pseudobulk_hpc(
    seurat = NULL,
    labels = labels,
    pb_variants = pb_variants,
    sample_col = sample_col,
    results_dir = args$results_dir,
    ds = ds,
    force = force,
    log_file = args$log_file,
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    h5ad_path = h5ad_path,
    ct_col_low_res = entry$cell_type_low_res,
    ct_col_high_res = entry$cell_type_high_res,
    view = args$view,
    analysis_pass = analysis_pass,
    run_id = ecoda_local_current_run_id(),
    source_identity = source_identity,
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    batch_context = if (correct_batch_mode) batch_context else NULL,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL
  ),
  scitd = run_scitd_hpc(
    seurat, label_col = entry$label_col,
    hvg_sets = make_hvg_sets(hvg_rank_genes),
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file
  ),
  composition = run_composition_methods_hpc(
    labels, metadata, pca_emb, pb_variants[["hvg2000"]], obs,
    label_col = entry$label_col,
    ct_col_low_res = entry$cell_type_low_res,
    ct_col_high_res = entry$cell_type_high_res,
    sample_col = sample_col,
    results_dir = args$results_dir, ds = ds,
    force = force, log_file = args$log_file,
    seurat_res = if (!is.null(analysis_pass)) 2 else c(0.1, 0.4, 2, 5, 20),
    batch_mode = !is.null(analysis_pass),
    result_stem = cache_stem,
    batch_col = batch_col,
    corrected = correct_batch_mode,
    batch_keys = if (correct_batch_mode) as.list(unname(batch_keys)) else NULL,
    metadata_validation = if (correct_batch_mode) {
      batch_context$validation
    } else {
      NULL
    },
    batch_contract = if (correct_batch_mode) method_batch_contract else NULL,
    not_suitable_for_auto_annotation = if (
      is.null(entry$not_suitable_for_auto_annotation)
    ) {
      character(0)
    } else {
      entry$not_suitable_for_auto_annotation
    }
  )
)

if (combo_supplied) {
  message(
    "Saved GloScope combo ", combo_token,
    " per-combo bundle/cache; method-level RDS deferred to consolidation"
  )
} else {
  ecoda_local_publish_rds(
    results, method_rds, producer = paste0("stage5_", method)
  )
  message("Saved method results: ", method_rds, " (", length(results), " combos)")
}
message("--- ", method, " for ", ds, " complete ---")
