# ==============================================================================
# Standalone MOFAcellulaR derived-analysis wrapper
#
# This file is intentionally sourceable.  No package is loaded, no H5AD is
# opened, and no job is started while the file is sourced.  The command-line
# guard at the end is the only execution entry point.
#
# MOFAcellulaR is conditional in this repository.  It is not installed by this
# wrapper and it is not part of the ordinary Stage 5 method registry.  A
# benchmark_union run is accepted only with an independently reviewed,
# checksum-verified _debug pass record for the exact package commit supplied on
# the command line.
# ==============================================================================

.mofa_missing_sentinels <- c(
  "", "na", "nan", "none", "null", "<na>", "<null>", "unknown",
  "n/a", "unassigned", "not_assigned", "not_annotated", "unannotated"
)
.mofa_package_name <- "MOFAcellulaR"
.mofa_package_repository <- "saezlab/MOFAcellulaR"
.mofa_package_repository_url <- "https://github.com/saezlab/MOFAcellulaR"
.mofa_hvg_rank_n <- 2000L
.mofa_final_ngenes <- 15L
.mofa_seed <- 42L
.mofa_variants <- c("lowres", "highres")
.mofa_view_stages <- c(
  "filt_profiles", "filt_gex_byexpr", "filt_views_bysamples",
  "filt_views_bygenes", "filt_samples_bycov", "tmm_trns",
  "filt_gex_byhvg", "final_filt_views_bygenes"
)

.mofa_or <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}

.mofa_scalar <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(label, " must be one non-empty string")
  }
  value
}

.mofa_nonempty_character <- function(value, label) {
  if (!is.character(value) || !length(value) || anyNA(value) ||
      any(!nzchar(trimws(value)))) {
    stop(label, " must contain only non-empty strings")
  }
  as.character(value)
}

.mofa_is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path, perl = TRUE)
}

.mofa_path <- function(path, must_work = FALSE) {
  path <- .mofa_scalar(as.character(path), "path")
  path <- path.expand(path)
  if (must_work) return(normalizePath(path, mustWork = TRUE))
  # normalizePath(mustWork = FALSE) can retain a symlinked ancestor.  Resolve
  # the deepest existing ancestor so every source/output identity is stable.
  candidate <- path
  suffix <- character()
  repeat {
    if (file.exists(candidate) || dir.exists(candidate)) {
      current <- normalizePath(candidate, mustWork = TRUE)
      if (!length(suffix)) return(current)
      unresolved <- character()
      for (component in suffix) {
        if (!nzchar(component) || identical(component, ".")) next
        if (identical(component, "..")) {
          if (length(unresolved)) {
            unresolved <- unresolved[-length(unresolved)]
          } else {
            current <- dirname(current)
          }
          next
        }
        if (length(unresolved)) {
          unresolved <- c(unresolved, component)
          next
        }
        next_path <- file.path(current, component)
        if (file.exists(next_path) || dir.exists(next_path)) {
          current <- normalizePath(next_path, mustWork = TRUE)
        } else {
          unresolved <- component
        }
      }
      if (!length(unresolved)) return(current)
      return(normalizePath(do.call(file.path, c(list(current), as.list(unresolved))),
                           mustWork = FALSE))
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    suffix <- c(basename(candidate), suffix)
    candidate <- parent
  }
  normalizePath(path, mustWork = FALSE)
}

.mofa_path_within <- function(path, root) {
  path <- .mofa_path(path)
  root <- .mofa_path(root)
  boundary <- if (identical(root, .Platform$file.sep)) {
    root
  } else {
    paste0(root, .Platform$file.sep)
  }
  identical(path, root) || startsWith(path, boundary)
}

.mofa_is_symlink <- function(path) {
  link <- Sys.readlink(path)
  length(link) == 1L && !is.na(link) && nzchar(link)
}

.mofa_validate_component <- function(value, label) {
  value <- .mofa_scalar(value, label)
  if (value %in% c(".", "..") ||
      grepl("^[A-Za-z]:", value, perl = TRUE) ||
      grepl("[/\\\\]", value, perl = TRUE) ||
      grepl("[[:cntrl:]]", value, perl = TRUE)) {
    stop(label, " contains an unsafe path component")
  }
  value
}

.mofa_run_id <- function(value) {
  value <- .mofa_scalar(value, "run_id")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", value, perl = TRUE)) {
    stop("run_id contains unsafe characters")
  }
  value
}

.mofa_script_path <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  token <- full[startsWith(full, "--file=")]
  if (length(token)) {
    candidate <- sub("^--file=", "", token[[1L]])
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  }
  candidate <- file.path(getwd(), "src", "5_run_benchmark_methods",
                         "run_mofacellular.R")
  if (file.exists(candidate)) normalizePath(candidate, mustWork = TRUE) else NA_character_
}

.mofa_project_root <- function() {
  configured <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (nzchar(configured) && dir.exists(configured)) {
    return(.mofa_path(configured, must_work = TRUE))
  }
  script <- .mofa_script_path()
  if (!is.na(script)) dirname(dirname(dirname(script))) else .mofa_path(getwd(), TRUE)
}

.mofa_md5 <- function(path) {
  path <- .mofa_path(path, must_work = TRUE)
  if (dir.exists(path) || .mofa_is_symlink(path)) {
    stop("Cannot checksum a directory or symbolic link: ", path)
  }
  digest <- unname(tools::md5sum(path))
  if (length(digest) != 1L || is.na(digest) ||
      !grepl("^[0-9a-fA-F]{32}$", digest, perl = TRUE)) {
    stop("Could not compute MD5 for: ", path)
  }
  tolower(as.character(digest))
}

.mofa_atomic_text <- function(path, text, refuse_existing = FALSE) {
  path <- .mofa_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (refuse_existing && (file.exists(path) || .mofa_is_symlink(path))) {
    stop("Refusing to replace existing run-owned file: ", path)
  }
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(as.character(text), temporary, useBytes = TRUE)
  if (file.exists(path) || .mofa_is_symlink(path)) {
    stop("Refusing to replace existing run-owned file: ", path)
  }
  if (!file.rename(temporary, path)) stop("Could not atomically publish: ", path)
  invisible(path)
}

.mofa_validate_checksum_sidecar <- function(path, require_path = TRUE) {
  path <- .mofa_path(path, must_work = TRUE)
  if (dir.exists(path) || .mofa_is_symlink(path) || file.info(path)$size <= 0) {
    stop("Checksum target is missing, empty, or not a regular file: ", path)
  }
  sidecar <- paste0(path, ".md5")
  if (!file.exists(sidecar) || .mofa_is_symlink(sidecar) ||
      is.na(file.info(sidecar)$size) || file.info(sidecar)$size <= 0) {
    stop("Checksum sidecar is missing: ", sidecar)
  }
  lines <- readLines(sidecar, warn = FALSE)
  if (length(lines) != 3L ||
      !identical(sub("=.*$", "", lines), c("MD5", "SIZE", "PATH"))) {
    stop("Checksum sidecar has the wrong schema: ", sidecar)
  }
  md5 <- sub("^MD5=", "", lines[[1L]])
  size <- sub("^SIZE=", "", lines[[2L]])
  recorded_path <- sub("^PATH=", "", lines[[3L]])
  if (!grepl("^[0-9a-fA-F]{32}$", md5, perl = TRUE) ||
      !grepl("^[0-9]+$", size, perl = TRUE) || !nzchar(recorded_path) ||
      !identical(size, as.character(file.info(path)$size)) ||
      (isTRUE(require_path) && !identical(.mofa_path(recorded_path), path)) ||
      (!isTRUE(require_path) &&
       !identical(basename(recorded_path), basename(path)))) {
    stop("Checksum sidecar metadata mismatch: ", sidecar)
  }
  actual <- .mofa_md5(path)
  if (!identical(tolower(md5), actual)) {
    stop("Checksum mismatch: ", path)
  }
  list(MD5 = actual, SIZE = size, PATH = path)
}

.mofa_write_checksum_sidecar <- function(path, refuse_existing = TRUE) {
  path <- .mofa_path(path, must_work = TRUE)
  sidecar <- paste0(path, ".md5")
  text <- c(
    paste0("MD5=", .mofa_md5(path)),
    paste0("SIZE=", as.character(file.info(path)$size)),
    paste0("PATH=", path)
  )
  .mofa_atomic_text(sidecar, text, refuse_existing = refuse_existing)
  .mofa_validate_checksum_sidecar(path)
  invisible(sidecar)
}

.mofa_output_path <- function(path, output_dir) {
  path <- .mofa_path(path)
  output_dir <- .mofa_path(output_dir)
  if (!.mofa_path_within(path, output_dir) || identical(path, output_dir)) {
    stop("Run-owned output escapes output_dir: ", path)
  }
  path
}

.mofa_prepare_output_dir <- function(output_dir, run_id) {
  run_id <- .mofa_run_id(run_id)
  output_dir <- .mofa_path(output_dir)
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stop("output_dir must name a directory: ", output_dir)
  }
  if (!dir.exists(output_dir) &&
      !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output_dir: ", output_dir)
  }
  entries <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(entries)) {
    stop("output_dir must be a new empty run-owned directory: ", output_dir)
  }
  owner <- file.path(output_dir, ".mofacellular_run_owner")
  .mofa_atomic_text(owner, paste0("RUN_ID=", run_id), refuse_existing = TRUE)
  list(output_dir = output_dir, run_id = run_id, owner = owner)
}

.mofa_write_rds_atomic <- function(value, path, output_dir) {
  path <- .mofa_output_path(path, output_dir)
  if (file.exists(path) || .mofa_is_symlink(path) ||
      file.exists(paste0(path, ".md5")) ||
      .mofa_is_symlink(paste0(path, ".md5"))) {
    stop("Refusing to replace existing run-owned artifact: ", path)
  }
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"),
                        tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = TRUE)
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
      file.info(temporary)$size <= 0L) {
    stop("Atomic RDS write produced an empty file: ", path)
  }
  if (file.exists(path) || .mofa_is_symlink(path) ||
      !file.rename(temporary, path)) {
    stop("Could not atomically publish RDS: ", path)
  }
  .mofa_write_checksum_sidecar(path, refuse_existing = TRUE)
  .mofa_validate_checksum_sidecar(path)
  invisible(path)
}

.mofa_write_json_with_sidecar <- function(value, path, output_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to publish MOFAcellulaR manifests")
  }
  path <- .mofa_output_path(path, output_dir)
  text <- jsonlite::toJSON(
    value, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null",
    dataframe = "rows"
  )
  .mofa_atomic_text(path, paste0(text, "\n"), refuse_existing = TRUE)
  .mofa_write_checksum_sidecar(path, refuse_existing = TRUE)
  .mofa_validate_checksum_sidecar(path)
  invisible(path)
}

.mofa_parse_flags <- function(raw_args) {
  args <- list()
  i <- 1L
  while (i <= length(raw_args)) {
    flag <- raw_args[[i]]
    if (!startsWith(flag, "--")) stop("Unexpected positional argument: ", flag)
    name <- sub("^--", "", flag)
    if (grepl("=", name, fixed = TRUE)) {
      pieces <- strsplit(name, "=", fixed = TRUE)[[1L]]
      if (length(pieces) != 2L || !nzchar(pieces[[1L]]) ||
          pieces[[1L]] %in% names(args)) {
        stop("Malformed or repeated argument: ", flag)
      }
      args[[pieces[[1L]]]] <- pieces[[2L]]
      i <- i + 1L
    } else if (i < length(raw_args) && !startsWith(raw_args[[i + 1L]], "--")) {
      if (name %in% names(args)) stop("Repeated argument: --", name)
      args[[name]] <- raw_args[[i + 1L]]
      i <- i + 2L
    } else {
      if (name %in% names(args)) stop("Repeated argument: --", name)
      args[[name]] <- TRUE
      i <- i + 1L
    }
  }
  args
}

.mofa_usage <- function() {
  paste(
    "Usage:",
    "Rscript src/5_run_benchmark_methods/run_mofacellular.R",
    "--scope debug|benchmark_union --config_path <datasets.json>",
    "--input_dir <HPC scratch root> --output_dir <new run directory>",
    "--run_id <run id> --package_sha <40-hex Git SHA>",
    "[--debug_pass_record <reviewed debug record>]"
  )
}

.mofa_validate_sha <- function(value, label = "package_sha") {
  value <- .mofa_scalar(as.character(value), label)
  if (!grepl("^[0-9a-fA-F]{40}$", value, perl = TRUE)) {
    stop(label, " must be an exact 40-character Git commit SHA; mutable branches are forbidden")
  }
  tolower(value)
}

.mofa_get_package_sha <- function(args) {
  candidates <- list(
    args$package_sha,
    args$mofacellular_sha,
    args$mofacellular_commit,
    args$package_commit
  )
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (!length(candidates)) {
    stop("Missing required --package_sha (the exact verified MOFAcellulaR Git SHA)")
  }
  normalized <- vapply(candidates, .mofa_validate_sha, character(1),
                       label = "package_sha")
  if (length(unique(normalized)) != 1L) {
    stop("MOFAcellulaR provenance arguments disagree")
  }
  normalized[[1L]]
}

.mofa_read_config <- function(config_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read datasets.json")
  }
  config_path <- .mofa_path(config_path, must_work = TRUE)
  if (dir.exists(config_path) || .mofa_is_symlink(config_path)) {
    stop("config_path must name a regular JSON file")
  }
  config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
  if (!is.list(config) || is.null(names(config)) || !length(config)) {
    stop("datasets config must be a non-empty named object")
  }
  list(config = config, path = config_path, md5 = .mofa_md5(config_path))
}

.mofa_scope_datasets <- function(config, scope) {
  scope <- match.arg(scope, c("debug", "benchmark_union"))
  if (identical(scope, "debug")) {
    if (!"_debug" %in% names(config)) {
      stop("debug scope requires the explicit _debug dataset entry")
    }
    return("_debug")
  }
  eligible <- vapply(config, function(entry) {
    if (!is.list(entry)) return(FALSE)
    has_view <- is.list(entry$views) && "benchmark_analysis" %in% names(entry$views)
    isTRUE(entry$use_for_benchmark) || has_view
  }, logical(1))
  selected <- names(config)[eligible & !startsWith(names(config), "_")]
  if (!length(selected)) stop("benchmark_union contains no eligible datasets")
  selected
}

.mofa_annotation_column <- function(entry, name, dataset) {
  columns <- entry$columns
  value <- if (is.list(columns)) columns[[name]] else NULL
  if (is.null(value) && is.list(entry$views$benchmark_analysis)) {
    value <- entry$views$benchmark_analysis$columns[[name]]
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(dataset, " has no configured ", name, " annotation column")
  }
  value
}

.mofa_view_output_name <- function(entry, dataset) {
  view <- entry$views$benchmark_analysis
  if (!is.list(view)) stop(dataset, " has no benchmark_analysis view")
  output <- .mofa_or(view$output_file_name, view$output_file)
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(output) || output %in% c(".", "..") ||
      .mofa_is_absolute_path(output) || grepl("^[A-Za-z]:", output, perl = TRUE) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", output, perl = TRUE) ||
      grepl("[[:cntrl:]]", output, perl = TRUE)) {
    stop(dataset, " benchmark_analysis output_file_name is invalid")
  }
  output
}

.mofa_resolve_h5ad <- function(input_dir, dataset, entry) {
  input_dir <- .mofa_path(input_dir, must_work = TRUE)
  dataset <- .mofa_validate_component(dataset, "dataset")
  output <- .mofa_view_output_name(entry, dataset)
  candidates <- unique(c(
    file.path(input_dir, dataset, "output", output),
    file.path(input_dir, output),
    file.path(input_dir, dataset, output)
  ))
  candidates <- unique(vapply(candidates, .mofa_path, character(1)))
  if (any(!vapply(candidates, .mofa_path_within, logical(1), root = input_dir))) {
    stop(dataset, " benchmark_analysis source escapes input_dir")
  }
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop(dataset, " benchmark_analysis H5AD is missing")
  if (length(existing) > 1L) {
    identities <- lapply(existing, .mofa_validate_checksum_sidecar)
    hashes <- vapply(identities, `[[`, character(1), "MD5")
    if (length(unique(hashes)) != 1L) {
      stop(dataset, " has ambiguous benchmark_analysis H5AD paths")
    }
  }
  .mofa_path(existing[[1L]], must_work = TRUE)
}

.mofa_preflight_sources <- function(config_info, input_dir, scope) {
  selected <- .mofa_scope_datasets(config_info$config, scope)
  records <- list()
  blocked <- list()
  for (dataset in selected) {
    entry <- config_info$config[[dataset]]
    result <- tryCatch({
      if (!is.list(entry) || !is.list(entry$views) ||
          !"benchmark_analysis" %in% names(entry$views)) {
        stop("declared benchmark_analysis view is missing; batch-effect views are not substitutes")
      }
      low_col <- .mofa_annotation_column(entry, "cell_type_low_res", dataset)
      high_col <- .mofa_annotation_column(entry, "cell_type_high_res", dataset)
      label_col <- .mofa_annotation_column(entry, "label", dataset)
      source_path <- .mofa_resolve_h5ad(input_dir, dataset, entry)
      source_identity <- .mofa_validate_checksum_sidecar(source_path)
      debug_source <- NULL
      if (identical(scope, "debug")) {
        if (!identical(dataset, "_debug")) {
          stop("debug scope may preflight only the explicit _debug dataset")
        }
        debug_source <- .mofa_validate_debug_source(
          source_path, source_identity
        )
      }
      source_record <- list(
        dataset = dataset,
        view = "benchmark_analysis",
        h5ad_path = source_path,
        h5ad_md5 = source_identity$MD5,
        h5ad_size = as.numeric(source_identity$SIZE),
        sample_column = "Sample",
        label_column = label_col,
        cell_type_low_res = low_col,
        cell_type_high_res = high_col
      )
      if (!is.null(debug_source)) {
        source_record$debug_sample_ids <- debug_source$sample_ids
        source_record$debug_sample_count <- debug_source$sample_count
      }
      source_record
    }, error = function(error) {
      blocked[[dataset]] <<- conditionMessage(error)
      NULL
    })
    if (!is.null(result)) records[[dataset]] <- result
  }
  if (length(blocked)) {
    text <- paste(
      vapply(names(blocked), function(dataset) {
        paste0(dataset, ": ", blocked[[dataset]])
      }, character(1)),
      collapse = " | "
    )
    stop("MOFAcellulaR source preflight blocked datasets: ", text)
  }
  records
}

.mofa_validate_package_provenance <- function(package_sha) {
  package_sha <- .mofa_validate_sha(package_sha)
  if (!requireNamespace(.mofa_package_name, quietly = TRUE)) {
    stop(.mofa_package_name, " is not installed at the required verified Git SHA;",
         " this wrapper never installs it")
  }
  desc <- packageDescription(.mofa_package_name)
  reported <- unique(tolower(trimws(as.character(c(
    desc[["GithubSHA1"]], desc[["RemoteSha"]], desc[["GitSHA"]]
  )))))
  reported <- reported[!is.na(reported) & nzchar(reported) &
                         grepl("^[0-9a-f]{40}$", reported, perl = TRUE)]
  if (!length(reported) || !package_sha %in% reported) {
    stop(.mofa_package_name, " provenance does not expose the required exact Git SHA ",
         package_sha)
  }
  list(
    package = .mofa_package_name,
    repository = .mofa_package_repository,
    repository_url = .mofa_package_repository_url,
    commit = package_sha,
    version = as.character(.mofa_or(desc[["Version"]], NA_character_))
  )
}

.mofa_load_hvg_counts <- function(h5ad_path, annotation_column, label_column) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to stream selected H5AD counts")
  }
  module_dir <- normalizePath(
    file.path(.mofa_project_root(), "src", "utils", "py"), mustWork = TRUE
  )
  loader <- reticulate::import_from_path(
    "h5ad_counts_subset", path = module_dir, convert = FALSE
  )
  genes <- reticulate::py_to_r(loader$read_h5ad_hvg_genes(
    h5ad_path, as.integer(.mofa_hvg_rank_n)
  ))
  genes <- as.character(genes)
  if (length(genes) != .mofa_hvg_rank_n || anyNA(genes) ||
      any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("H5AD var[hvg_rank] did not yield exactly the top 2000 unique genes")
  }
  obs_columns <- unique(c("Sample", annotation_column, label_column))
  minimal <- loader$load_h5ad_counts_subset(
    h5ad_path,
    as.list(genes),
    as.list(obs_columns),
    as.integer(4096L)
  )
  obs <- reticulate::py_to_r(minimal$obs)
  counts <- reticulate::py_to_r(minimal$X)
  if (!is.data.frame(obs)) obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    counts <- tryCatch(as.matrix(counts), error = function(error) NULL)
  }
  if (is.null(counts) || length(dim(counts)) != 2L) {
    stop("Selected H5AD counts could not be represented as a two-dimensional matrix")
  }
  genes_from_var <- tryCatch(as.character(rownames(minimal$var)), error = function(error) NULL)
  if (is.null(genes_from_var) || length(genes_from_var) != length(genes)) {
    genes_from_var <- genes
  }
  rownames(counts) <- NULL
  if (ncol(counts) != length(genes)) {
    stop("Selected H5AD counts have the wrong gene dimension")
  }
  if (nrow(counts) != nrow(obs)) {
    stop("Selected H5AD counts and obs have inconsistent cell rows")
  }
  colnames(counts) <- genes_from_var
  list(counts = counts, obs = obs, genes = genes_from_var)
}

.mofa_load_h5ad_obs_free <- function(h5ad_path, obs_columns) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required to read H5AD metadata")
  }
  obs_columns <- .mofa_nonempty_character(
    as.character(obs_columns), "obs_columns"
  )
  module_dir <- normalizePath(
    file.path(.mofa_project_root(), "src", "utils", "py"), mustWork = TRUE
  )
  loader <- reticulate::import_from_path(
    "h5ad_obs_free", path = module_dir, convert = FALSE
  )
  obs <- reticulate::py_to_r(loader$load_h5ad_obs_free(
    h5ad_path, as.list(obs_columns)
  ))
  if (!is.data.frame(obs)) {
    obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  }
  for (column in intersect(obs_columns, colnames(obs))) {
    if (is.factor(obs[[column]])) obs[[column]] <- as.character(obs[[column]])
  }
  obs
}

.mofa_text_vector <- function(values, missing_as_sentinel = TRUE) {
  if (is.null(values)) return(character())
  out <- vapply(seq_along(values), function(index) {
    value <- values[[index]]
    if (is.null(value) || !length(value) || is.na(value[[1L]])) return(NA_character_)
    text <- as.character(value[[1L]])
    if (length(text) != 1L || is.na(text)) return(NA_character_)
    text <- trimws(text)
    if (!nzchar(text)) return(NA_character_)
    if (missing_as_sentinel && tolower(text) %in% .mofa_missing_sentinels) {
      return(NA_character_)
    }
    text
  }, character(1))
  unname(out)
}

.mofa_validate_debug_sample_record <- function(
  sample_ids, sample_count, label = "debug Sample provenance"
) {
  ids <- tryCatch(
    as.character(unlist(sample_ids, use.names = FALSE)),
    error = function(error) NULL
  )
  count_raw <- tryCatch(
    unlist(sample_count, use.names = FALSE),
    error = function(error) NULL
  )
  count <- suppressWarnings(as.numeric(count_raw))
  invalid_ids <- is.null(ids) || !length(ids) || anyNA(ids) ||
    any(!nzchar(trimws(ids))) || anyDuplicated(ids)
  invalid_count <- length(count) != 1L || is.na(count) ||
    !is.finite(count) || count != floor(count)
  if (invalid_ids || invalid_count) {
    stop(label, " is missing or invalid")
  }
  if (count != length(ids) || count < 2L) {
    stop(
      label, " must contain at least two unique Sample IDs; found ",
      length(ids), " IDs (", paste(ids, collapse = ", "), ")"
    )
  }
  list(
    sample_ids = sort(ids),
    sample_count = as.integer(count)
  )
}

.mofa_validate_debug_source <- function(h5ad_path, source_identity = NULL) {
  h5ad_path <- .mofa_path(h5ad_path, must_work = TRUE)
  identity <- .mofa_validate_checksum_sidecar(h5ad_path)
  if (!is.null(source_identity) &&
      (!is.list(source_identity) ||
       !identical(source_identity$PATH, identity$PATH) ||
       !identical(source_identity$MD5, identity$MD5) ||
       !identical(source_identity$SIZE, identity$SIZE))) {
    stop("_debug source checksum identity changed before metadata guard")
  }
  obs <- .mofa_load_h5ad_obs_free(h5ad_path, "Sample")
  if (!is.data.frame(obs) || !"Sample" %in% colnames(obs)) {
    stop("_debug benchmark_analysis H5AD has no Sample metadata column")
  }
  samples <- .mofa_text_vector(obs[["Sample"]], missing_as_sentinel = TRUE)
  if (!length(samples) || anyNA(samples)) {
    stop("_debug benchmark_analysis H5AD Sample IDs are missing or invalid")
  }
  sample_ids <- unique(samples)
  checked <- .mofa_validate_debug_sample_record(
    sample_ids = sample_ids,
    sample_count = length(sample_ids),
    label = "_debug benchmark_analysis Sample IDs"
  )
  list(
    path = identity$PATH,
    md5 = identity$MD5,
    size = as.numeric(identity$SIZE),
    sample_ids = checked$sample_ids,
    sample_count = checked$sample_count
  )
}
.mofa_validate_count_values <- function(counts) {
  values <- if (inherits(counts, "Matrix") &&
                "x" %in% methods::slotNames(counts)) {
    methods::slot(counts, "x")
  } else {
    as.vector(counts)
  }
  if (any(!is.finite(values)) || any(values < 0) ||
      any(values != floor(values))) {
    stop("Selected H5AD counts must be finite, nonnegative integer-valued counts")
  }
  invisible(TRUE)
}

.mofa_aggregate_profiles <- function(h5ad_path, annotation_column, label_column) {
  selected <- .mofa_load_hvg_counts(h5ad_path, annotation_column, label_column)
  obs <- selected$obs
  counts <- selected$counts
  required <- c("Sample", annotation_column, label_column)
  missing <- setdiff(required, colnames(obs))
  if (length(missing)) stop("H5AD obs is missing required columns: ", paste(missing, collapse = ", "))
  samples <- .mofa_text_vector(obs[["Sample"]], missing_as_sentinel = TRUE)
  cell_types <- .mofa_text_vector(obs[[annotation_column]], missing_as_sentinel = TRUE)
  score_labels <- .mofa_text_vector(obs[[label_column]], missing_as_sentinel = TRUE)
  if (length(samples) == 0L || nrow(counts) != length(samples)) {
    stop("H5AD has no usable Sample observations")
  }
  .mofa_validate_count_values(counts)
  source_sample_ids <- unique(samples[!is.na(samples)])
  valid_cell <- !is.na(samples) & !is.na(cell_types)
  # Missing/nonannotated cells are removed before any sample+cell-type grouping.
  # There is deliberately no Unassigned category and no fallback annotation.
  if (!any(valid_cell)) stop("No cells have both a valid Sample and annotation")
  valid_sample_ids <- unique(samples[valid_cell])
  sample_labels <- setNames(character(length(valid_sample_ids)), valid_sample_ids)
  for (sample_id in valid_sample_ids) {
    values <- score_labels[valid_cell & samples == sample_id]
    if (!length(values) || anyNA(values) || any(!nzchar(values))) {
      stop("Missing scoring label for retained sample: ", sample_id)
    }
    values <- unique(values)
    if (length(values) != 1L) {
      stop("Conflicting scoring labels for retained sample: ", sample_id)
    }
    sample_labels[[sample_id]] <- values[[1L]]
  }
  valid_indices <- which(valid_cell)
  group_key <- paste(samples[valid_cell], cell_types[valid_cell], sep = "\r")
  group_keys <- unique(group_key)
  group_index <- match(group_key, group_keys)
  profile_samples <- vapply(strsplit(group_keys, "\r", fixed = TRUE), `[[`,
                            character(1), 1L)
  profile_cell_types <- vapply(strsplit(group_keys, "\r", fixed = TRUE), `[[`,
                               character(1), 2L)
  profile_ids <- paste0("Sample=", profile_samples, ";cell_type=", profile_cell_types)
  if (anyDuplicated(profile_ids) || anyNA(profile_ids) || any(!nzchar(profile_ids))) {
    stop("Pseudobulk profile identifiers are not finite and unique")
  }
  pb_counts <- matrix(
    0, nrow = ncol(counts), ncol = length(group_keys),
    dimnames = list(colnames(counts), profile_ids)
  )
  cells_by_group <- split(valid_indices, group_index)
  for (group in seq_along(group_keys)) {
    cell_rows <- cells_by_group[[as.character(group)]]
    if (!length(cell_rows)) stop("Internal pseudobulk grouping produced an empty profile")
    pb_counts[, group] <- base::colSums(counts[cell_rows, , drop = FALSE])
  }
  coldata <- data.frame(
    Sample = profile_samples,
    cell_type = profile_cell_types,
    cell_counts = as.integer(vapply(cells_by_group, length, integer(1))),
    row.names = profile_ids,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    counts = pb_counts,
    coldata = coldata,
    source_sample_ids = source_sample_ids,
    retained_annotation_sample_ids = valid_sample_ids,
    sample_labels = sample_labels,
    genes = colnames(counts)
  )
}

.mofa_view_sample_ids <- function(view) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required to inspect MOFAcellulaR views")
  }
  metadata <- SummarizedExperiment::colData(view)
  if (!"Sample" %in% colnames(metadata)) {
    stop("MOFAcellulaR view has no Sample metadata")
  }
  ids <- as.character(metadata[["Sample"]])
  if (!length(ids)) return(character())
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("MOFAcellulaR view contains missing or duplicate Sample IDs")
  }
  ids
}

.mofa_view_ids <- function(views) {
  if (is.null(views) || !length(views)) return(list())
  if (is.null(names(views))) names(views) <- paste0("view", seq_along(views))
  lapply(views, .mofa_view_sample_ids)
}

.mofa_apply_filters <- function(profile_data, annotation_column, factor_count,
                                source_sample_ids) {
  if (!requireNamespace("MOFAcellulaR", quietly = TRUE)) {
    stop("MOFAcellulaR is unavailable; this wrapper never installs it")
  }
  cts <- sort(unique(as.character(profile_data$coldata$cell_type)))
  if (!length(cts) || anyNA(cts) || any(!nzchar(cts)) ||
      any(tolower(cts) %in% .mofa_missing_sentinels)) {
    stop("Pseudobulk views contain no valid annotated cell types")
  }
  pb_obj <- MOFAcellulaR::create_init_exp(
    counts = profile_data$counts,
    coldata = profile_data$coldata
  )
  views <- MOFAcellulaR::filt_profiles(
    pb_dat = pb_obj,
    cts = cts,
    ncells = 0,
    counts_col = "cell_counts",
    ct_col = "cell_type"
  )
  stage_ids <- list(filt_profiles = .mofa_view_ids(views))
  if (!length(views)) stop("filt_profiles produced no cell-type views")

  views <- MOFAcellulaR::filt_gex_byexpr(
    pb_dat_list = views,
    min.count = 5,
    min.prop = 0.25
  )
  stage_ids$filt_gex_byexpr <- .mofa_view_ids(views)

  views <- MOFAcellulaR::filt_views_bysamples(
    pb_dat_list = views,
    nsamples = 2
  )
  stage_ids$filt_views_bysamples <- .mofa_view_ids(views)

  views <- MOFAcellulaR::filt_views_bygenes(
    pb_dat_list = views,
    ngenes = 15
  )
  stage_ids$filt_views_bygenes <- .mofa_view_ids(views)

  views <- MOFAcellulaR::filt_samples_bycov(
    pb_dat_list = views,
    prop_coverage = 0.9
  )
  stage_ids$filt_samples_bycov <- .mofa_view_ids(views)

  views <- MOFAcellulaR::tmm_trns(
    pb_dat_list = views,
    scale_factor = 1000000
  )
  stage_ids$tmm_trns <- .mofa_view_ids(views)

  views <- MOFAcellulaR::filt_gex_byhvg(
    pb_dat_list = views,
    prior_hvg = NULL,
    var.threshold = 0
  )
  stage_ids$filt_gex_byhvg <- .mofa_view_ids(views)

  views <- MOFAcellulaR::filt_views_bygenes(
    pb_dat_list = views,
    ngenes = 15
  )
  stage_ids$final_filt_views_bygenes <- .mofa_view_ids(views)
  if (!length(views)) stop("Final 15-gene view filtering removed every view")
  final_view_samples <- stage_ids$final_filt_views_bygenes
  available_samples <- sort(unique(unlist(final_view_samples, use.names = FALSE)))
  if (length(available_samples) < 2L) {
    stop("MOFAcellulaR retained fewer than two samples after filtering")
  }
  if (identical(factor_count, 15L) && factor_count >= length(available_samples)) {
    stop("Production MOFAcellulaR factor count is not smaller than available samples")
  }
  view_names <- unique(unlist(lapply(stage_ids, names), use.names = FALSE))
  per_view_dropped <- setNames(lapply(view_names, function(view_name) {
    profile_ids <- .mofa_or(stage_ids$filt_profiles[[view_name]], character())
    list(
      filt_profiles = setdiff(source_sample_ids, profile_ids),
      filt_gex_byexpr = setdiff(profile_ids,
                                .mofa_or(stage_ids$filt_gex_byexpr[[view_name]],
                                         character())),
      filt_views_bysamples = setdiff(
        .mofa_or(stage_ids$filt_gex_byexpr[[view_name]], character()),
        .mofa_or(stage_ids$filt_views_bysamples[[view_name]], character())
      ),
      filt_views_bygenes = setdiff(
        .mofa_or(stage_ids$filt_views_bysamples[[view_name]], character()),
        .mofa_or(stage_ids$filt_views_bygenes[[view_name]], character())
      ),
      filt_samples_bycov = setdiff(
        .mofa_or(stage_ids$filt_views_bygenes[[view_name]], character()),
        .mofa_or(stage_ids$filt_samples_bycov[[view_name]], character())
      ),
      tmm_trns = setdiff(
        .mofa_or(stage_ids$filt_samples_bycov[[view_name]], character()),
        .mofa_or(stage_ids$tmm_trns[[view_name]], character())
      ),
      filt_gex_byhvg = setdiff(
        .mofa_or(stage_ids$tmm_trns[[view_name]], character()),
        .mofa_or(stage_ids$filt_gex_byhvg[[view_name]], character())
      ),
      final_filt_views_bygenes = setdiff(
        .mofa_or(stage_ids$filt_gex_byhvg[[view_name]], character()),
        .mofa_or(stage_ids$final_filt_views_bygenes[[view_name]],
                character())
      ),
      final = setdiff(
        source_sample_ids,
        .mofa_or(final_view_samples[[view_name]], character())
      )
    )
  }), view_names)
  list(
    views = views,
    stages = stage_ids,
    final_view_samples = final_view_samples,
    available_samples = available_samples,
    per_view_dropped_sample_ids = per_view_dropped,
    filtering = list(
      filt_profiles = list(ncells = 0L, counts_col = "cell_counts", ct_col = "cell_type"),
      filt_gex_byexpr = list(min.count = 5, min.prop = 0.25),
      filt_views_bysamples = list(nsamples = 2L),
      filt_views_bygenes = list(ngenes = 15L),
      filt_samples_bycov = list(prop_coverage = 0.9),
      tmm_trns = list(scale_factor = 1000000),
      filt_gex_byhvg = list(prior_hvg = NULL, var.threshold = 0),
      final_filt_views_bygenes = list(ngenes = 15L),
      hvg_rank_genes_in_memory = 2000L,
      biological_covariates = character(0)
    ),
    annotation_column = annotation_column
  )
}

.mofa_get_bundle_constructor <- function() {
  if (exists("create_result_bundle", envir = .GlobalEnv,
             inherits = FALSE) &&
      is.function(get("create_result_bundle", envir = .GlobalEnv))) {
    return(get("create_result_bundle", envir = .GlobalEnv))
  }
  root <- .mofa_project_root()
  methods_path <- file.path(root, "src", "5_run_benchmark_methods",
                            "benchmark_methods_r.R")
  scoring_path <- file.path(root, "src", "utils", "scoring_metrics.R")
  if (!file.exists(methods_path) || !file.exists(scoring_path)) {
    stop("Existing create_result_bundle source is unavailable")
  }
  contract_env <- new.env(parent = .GlobalEnv)
  if (requireNamespace("magrittr", quietly = TRUE)) {
    contract_env[["%>%"]] <- get("%>%", envir = asNamespace("magrittr"))
  }
  source(scoring_path, local = contract_env)
  source(methods_path, local = contract_env)
  if (!exists("create_result_bundle", envir = contract_env,
             inherits = FALSE) ||
      !is.function(contract_env$create_result_bundle)) {
    stop("Existing create_result_bundle source is invalid")
  }
  contract_env$create_result_bundle
}

.mofa_extract_factors <- function(model, expected_factors) {
  factors <- NULL
  if (exists("get_factors", envir = asNamespace("MOFA2"), inherits = FALSE)) {
    factor_groups <- tryCatch(
      MOFA2::get_factors(model, factors = "all"),
      error = function(error) NULL
    )
    if (is.list(factor_groups) && length(factor_groups) == 1L) {
      factors <- factor_groups[[1L]]
    }
  }
  # Keep a compatibility path for older MOFA2 builds whose get_factors()
  # generic is not exported.  This still reads the model's actual Z rows and
  # never creates or pads a latent row.
  if (is.null(factors) && methods::isS4(model) &&
      "expectations" %in% methods::slotNames(model)) {
    factors <- methods::slot(model, "expectations")[["Z"]]
  }
  if (is.null(factors)) stop("MOFA2 returned no actual factor rows")
  factors <- as.matrix(factors)
  if (length(dim(factors)) != 2L || nrow(factors) < 2L ||
      ncol(factors) != expected_factors || is.null(rownames(factors)) ||
      anyNA(rownames(factors)) || any(!nzchar(trimws(rownames(factors)))) ||
      anyDuplicated(rownames(factors)) || any(!is.finite(factors))) {
    stop("MOFA2 actual factor rows are missing, nonfinite, duplicated, or wrong-sized")
  }
  factors
}
.mofa_fit <- function(filtered, factor_count) {
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("MOFA2 is unavailable")
  }
  multiview_dat <- MOFAcellulaR::pb_dat2MOFA(
    pb_dat_list = filtered$views,
    sample_column = "Sample"
  )
  if (!is.data.frame(multiview_dat) || !all(c("view", "feature", "sample", "value") %in%
                                             colnames(multiview_dat)) ||
      !nrow(multiview_dat)) {
    stop("MOFAcellulaR did not produce a valid multi-view data frame")
  }
  # Only expression views are passed to MOFA.  Biological labels are kept in
  # profile_data$sample_labels and are never included in multiview_dat or
  # samples_metadata/model covariates.
  mofa_object <- MOFA2::create_mofa(multiview_dat)
  data_options <- MOFA2::get_default_data_options(mofa_object)
  model_options <- MOFA2::get_default_model_options(mofa_object)
  training_options <- MOFA2::get_default_training_options(mofa_object)
  model_options$num_factors <- as.integer(factor_count)
  model_options$spikeslab_weights <- FALSE
  training_options$convergence_mode <- "fast"
  training_options$seed <- as.integer(.mofa_seed)
  mofa_object <- MOFA2::prepare_mofa(
    object = mofa_object,
    data_options = data_options,
    model_options = model_options,
    training_options = training_options
  )
  model <- MOFA2::run_mofa(
    mofa_object, use_basilisk = FALSE, save_data = FALSE
  )
  factors <- .mofa_extract_factors(model, factor_count)
  list(
    model = model,
    factors = factors,
    multiview_rows = nrow(multiview_dat),
    options = list(
      num_factors = as.integer(factor_count),
      seed = as.integer(.mofa_seed),
      convergence_mode = "fast",
      spikeslab_weights = FALSE
    )
  )
}

.mofa_result_extra <- function(dataset, scope, run_id, source_identity,
                               package_info, annotation_column, label_column,
                               variant, factor_count, factors, profile_data,
                               filtered) {
  factor_ids <- rownames(factors)
  if (is.null(factor_ids) || any(!factor_ids %in% filtered$available_samples)) {
    stop("Actual MOFA factor rows contain IDs absent from retained views")
  }
  labels <- profile_data$sample_labels[factor_ids]
  if (length(labels) != length(factor_ids) || anyNA(labels) ||
      any(!nzchar(labels))) {
    stop("Actual MOFA factor rows cannot be matched to finite scoring labels")
  }
  dropped <- setdiff(profile_data$source_sample_ids, factor_ids)
  list(
    schema_version = 1L,
    dataset = dataset,
    method = paste0("MOFAcellulaR_hvg2000_factors", factor_count,
                    "_", variant),
    scope = scope,
    run_id = run_id,
    source_h5ad = source_identity$PATH,
    source_h5ad_md5 = source_identity$MD5,
    source_h5ad_size = as.numeric(source_identity$SIZE),
    sample_column = "Sample",
    annotation_column = annotation_column,
    label_column = label_column,
    labels_are_scoring_only = TRUE,
    model_covariates = character(0),
    package = package_info$package,
    package_repository = package_info$repository,
    package_repository_url = package_info$repository_url,
    package_commit = package_info$commit,
    package_version = package_info$version,
    factor_count = as.integer(factor_count),
    seed = as.integer(.mofa_seed),
    convergence_mode = "fast",
    spikeslab_weights = FALSE,
    hvg_rank_column = "hvg_rank",
    hvg_genes_requested = as.integer(.mofa_hvg_rank_n),
    hvg_genes_used = as.integer(length(profile_data$genes)),
    filtering = filtered$filtering,
    input_sample_ids = profile_data$source_sample_ids,
    factor_sample_ids = factor_ids,
    sample_ids = factor_ids,
    dropped_sample_ids = dropped,
    per_view_dropped_sample_ids = filtered$per_view_dropped_sample_ids,
    final_view_sample_ids = filtered$final_view_samples
  )
}

.mofa_make_bundle <- function(source_record, scope, run_id, package_info,
                              variant, factor_count) {
  annotation_column <- if (identical(variant, "lowres")) {
    source_record$cell_type_low_res
  } else {
    source_record$cell_type_high_res
  }
  source_before <- .mofa_validate_checksum_sidecar(source_record$h5ad_path)
  profile_data <- .mofa_aggregate_profiles(
    source_record$h5ad_path, annotation_column, source_record$label_column
  )
  filtered <- .mofa_apply_filters(
    profile_data = profile_data,
    annotation_column = annotation_column,
    factor_count = factor_count,
    source_sample_ids = profile_data$source_sample_ids
  )
  fit <- .mofa_fit(filtered, factor_count)
  factors <- fit$factors
  if (identical(scope, "benchmark_union") &&
      !identical(as.integer(factor_count), 15L)) {
    stop("benchmark_union requires exactly 15 factors")
  }
  if (identical(scope, "debug") && !identical(as.integer(factor_count), 2L)) {
    stop("debug requires exactly 2 factors")
  }
  factor_ids <- rownames(factors)
  labels <- profile_data$sample_labels[factor_ids]
  if (length(labels) != nrow(factors) || anyNA(labels)) {
    stop("Actual factor rows and scoring labels are not exactly aligned")
  }
  create_result_bundle <- .mofa_get_bundle_constructor()
  bundle_extra <- .mofa_result_extra(
    dataset = source_record$dataset,
    scope = scope,
    run_id = run_id,
    source_identity = source_before,
    package_info = package_info,
    annotation_column = annotation_column,
    label_column = source_record$label_column,
    variant = variant,
    factor_count = factor_count,
    factors = factors,
    profile_data = profile_data,
    filtered = filtered
  )
  # This is the actual existing repository constructor.  No factor rows are
  # padded and no synthetic score/features are added around it.
  bundle <- create_result_bundle(
    feat_mat = factors,
    labels = labels,
    dist_mat = stats::dist(factors),
    extra = bundle_extra
  )
  source_after <- .mofa_validate_checksum_sidecar(source_record$h5ad_path)
  if (!identical(source_before$MD5, source_after$MD5) ||
      !identical(source_before$SIZE, source_after$SIZE)) {
    stop("Source H5AD changed while building MOFAcellulaR result")
  }
  bundle
}

.mofa_filtering_matches <- function(filtering) {
  if (!is.list(filtering)) return(FALSE)
  required <- c(
    "filt_profiles", "filt_gex_byexpr", "filt_views_bysamples",
    "filt_views_bygenes", "filt_samples_bycov", "tmm_trns",
    "filt_gex_byhvg", "final_filt_views_bygenes"
  )
  if (!all(required %in% names(filtering))) return(FALSE)
  isTRUE(identical(as.numeric(filtering$filt_profiles$ncells), 0)) &&
    isTRUE(identical(as.numeric(filtering$filt_gex_byexpr$min.count), 5)) &&
    isTRUE(identical(as.numeric(filtering$filt_gex_byexpr$min.prop), 0.25)) &&
    isTRUE(identical(as.numeric(filtering$filt_views_bysamples$nsamples), 2)) &&
    isTRUE(identical(as.numeric(filtering$filt_views_bygenes$ngenes), 15)) &&
    isTRUE(identical(as.numeric(filtering$filt_samples_bycov$prop_coverage), 0.9)) &&
    isTRUE(identical(as.numeric(filtering$tmm_trns$scale_factor), 1000000)) &&
    is.null(filtering$filt_gex_byhvg$prior_hvg) &&
    isTRUE(identical(as.numeric(filtering$filt_gex_byhvg$var.threshold), 0)) &&
    isTRUE(identical(as.numeric(filtering$final_filt_views_bygenes$ngenes), 15)) &&
    isTRUE(identical(as.numeric(filtering$hvg_rank_genes_in_memory), 2000)) &&
    isTRUE(length(filtering$biological_covariates) == 0L)
}

# Dedicated validator for standalone MOFAcellulaR bundles.  `allow_dropped`
# defaults to FALSE for safety; callers must opt in explicitly for this method.
# The returned dropped_sample_ids is an auditable report, not a reconstructed
# row set.
validate_mofacellular_artifact <- function(
  path,
  expected_dataset = NULL,
  expected_method = NULL,
  expected_scope = NULL,
  expected_factor_count = NULL,
  expected_source_md5 = NULL,
  expected_input_samples = NULL,
  allow_dropped = FALSE
) {
  path <- .mofa_path(path, must_work = TRUE)
  sidecar <- .mofa_validate_checksum_sidecar(path)
  bundle <- readRDS(path)
  if (!is.list(bundle)) stop("MOFAcellulaR artifact is not an R list")
  required <- c(
    "scores", "feat_mat", "dist_mat", "labels", "dataset", "method", "scope",
    "source_h5ad_md5", "package_commit", "annotation_column", "factor_count",
    "sample_ids", "factor_sample_ids", "input_sample_ids", "dropped_sample_ids",
    "per_view_dropped_sample_ids", "filtering", "labels_are_scoring_only",
    "model_covariates"
  )
  missing <- setdiff(required, names(bundle))
  if (length(missing)) {
    stop("MOFAcellulaR artifact is missing: ", paste(missing, collapse = ", "))
  }
  method <- as.character(bundle$method)
  if (length(method) != 1L || !grepl(
    "^MOFAcellulaR_hvg2000_factors(2|15)_(lowres|highres)$",
    method,
    perl = TRUE
  )) {
    stop("MOFAcellulaR artifact method key is invalid")
  }
  if (!is.null(expected_method) && !identical(method, as.character(expected_method))) {
    stop("MOFAcellulaR artifact method key does not match expectation")
  }
  if (!is.character(bundle$dataset) || length(bundle$dataset) != 1L ||
      !nzchar(bundle$dataset) || (!is.null(expected_dataset) &&
                                  !identical(bundle$dataset, as.character(expected_dataset)))) {
    stop("MOFAcellulaR artifact dataset is invalid")
  }
  if (!identical(bundle$scope, "debug") && !identical(bundle$scope, "benchmark_union")) {
    stop("MOFAcellulaR artifact scope is invalid")
  }
  if (!is.null(expected_scope) && !identical(bundle$scope, as.character(expected_scope))) {
    stop("MOFAcellulaR artifact scope does not match expectation")
  }
  factor_count <- as.integer(bundle$factor_count)
  if (length(factor_count) != 1L || is.na(factor_count) ||
      !factor_count %in% c(2L, 15L) ||
      (!is.null(expected_factor_count) && factor_count != as.integer(expected_factor_count))) {
    stop("MOFAcellulaR artifact factor count is invalid")
  }
  if (!grepl("^[0-9a-fA-F]{32}$", as.character(bundle$source_h5ad_md5), perl = TRUE) ||
      !grepl("^[0-9a-fA-F]{40}$", as.character(bundle$package_commit), perl = TRUE) ||
      !isTRUE(bundle$labels_are_scoring_only) ||
      !is.character(bundle$model_covariates) || length(bundle$model_covariates) != 0L ||
      !.mofa_filtering_matches(bundle$filtering)) {
    stop("MOFAcellulaR provenance/filtering contract is invalid")
  }
  if (!is.character(bundle$annotation_column) || length(bundle$annotation_column) != 1L ||
      !nzchar(bundle$annotation_column) || identical(tolower(bundle$annotation_column), "unassigned")) {
    stop("MOFAcellulaR annotation column is invalid")
  }
  feat <- as.matrix(bundle$feat_mat)
  if (length(dim(feat)) != 2L || nrow(feat) < 2L || ncol(feat) != factor_count ||
      is.null(rownames(feat)) || anyNA(rownames(feat)) ||
      any(!nzchar(rownames(feat))) || anyDuplicated(rownames(feat)) ||
      any(!is.finite(feat))) {
    stop("MOFAcellulaR feature matrix has invalid actual factor rows")
  }
  factor_ids <- as.character(rownames(feat))
  if (!identical(as.character(bundle$sample_ids), factor_ids) ||
      !identical(as.character(bundle$factor_sample_ids), factor_ids) ||
      anyNA(bundle$labels) || is.null(names(bundle$labels)) ||
      !identical(sort(names(bundle$labels)), sort(factor_ids))) {
    stop("MOFAcellulaR labels are not aligned to actual factor rows")
  }
  if (!is.list(bundle$scores) || !length(bundle$scores)) {
    stop("MOFAcellulaR artifact scores are missing")
  }
  distance <- tryCatch(as.matrix(bundle$dist_mat), error = function(error) NULL)
  if (is.null(distance) || length(dim(distance)) != 2L ||
      nrow(distance) != nrow(feat) || ncol(distance) != nrow(feat) ||
      any(!is.finite(distance))) {
    stop("MOFAcellulaR distance matrix is invalid")
  }
  dropped <- as.character(bundle$dropped_sample_ids)
  if (length(dropped) && (anyNA(dropped) || any(!nzchar(dropped)) || anyDuplicated(dropped))) {
    stop("MOFAcellulaR dropped sample IDs are invalid")
  }
  input_ids <- as.character(bundle$input_sample_ids)
  if (!length(input_ids) || anyNA(input_ids) || any(!nzchar(input_ids)) ||
      anyDuplicated(input_ids) || any(!factor_ids %in% input_ids) ||
      !setequal(dropped, setdiff(input_ids, factor_ids))) {
    stop("MOFAcellulaR final dropped-sample accounting is invalid")
  }
  if (!is.null(expected_input_samples)) {
    expected_input_samples <- as.character(expected_input_samples)
    if (!setequal(input_ids, expected_input_samples)) {
      stop("MOFAcellulaR input sample universe does not match expectation")
    }
  }
  if (length(dropped) && (!isTRUE(allow_dropped) ||
                          !startsWith(method, "MOFAcellulaR_"))) {
    stop("Dropped samples are permitted only for explicitly opted-in MOFAcellulaR artifacts")
  }
  per_view <- bundle$per_view_dropped_sample_ids
  if (!is.list(per_view) || is.null(names(per_view))) {
    stop("MOFAcellulaR per-view dropped-sample report is missing")
  }
  for (view_name in names(per_view)) {
    report <- per_view[[view_name]]
    if (!is.list(report)) stop("MOFAcellulaR per-view drop report is malformed")
    for (ids in report) {
      ids <- as.character(ids)
      if (length(ids) && (anyNA(ids) || any(!nzchar(ids)) ||
                          anyDuplicated(ids) || any(!ids %in% input_ids))) {
        stop("MOFAcellulaR per-view dropped IDs are invalid")
      }
    }
  }
  if (!is.null(expected_source_md5) &&
      !identical(tolower(as.character(expected_source_md5)),
                 tolower(as.character(bundle$source_h5ad_md5)))) {
    stop("MOFAcellulaR source checksum does not match expectation")
  }
  list(
    valid = TRUE,
    path = path,
    md5 = sidecar$MD5,
    dataset = bundle$dataset,
    method = method,
    scope = bundle$scope,
    factor_count = factor_count,
    actual_sample_ids = factor_ids,
    dropped_sample_ids = dropped,
    source_h5ad_md5 = tolower(as.character(bundle$source_h5ad_md5))
  )
}

.mofa_output_filename <- function(dataset, variant, factor_count) {
  dataset <- .mofa_validate_component(dataset, "dataset")
  if (!variant %in% .mofa_variants) stop("unknown MOFAcellulaR variant: ", variant)
  paste0(dataset, "_MOFAcellulaR_hvg2000_factors", factor_count, "_",
         variant, ".rds")
}

.mofa_publish_bundle <- function(bundle, source_record, output_dir, scope,
                                 run_id, variant, factor_count) {
  filename <- .mofa_output_filename(source_record$dataset, variant, factor_count)
  path <- file.path(output_dir, filename)
  .mofa_write_rds_atomic(bundle, path, output_dir)
  expected_method <- paste0("MOFAcellulaR_hvg2000_factors", factor_count,
                            "_", variant)
  checked <- validate_mofacellular_artifact(
    path,
    expected_dataset = source_record$dataset,
    expected_method = expected_method,
    expected_scope = scope,
    expected_factor_count = factor_count,
    expected_source_md5 = source_record$h5ad_md5,
    expected_input_samples = bundle$input_sample_ids,
    allow_dropped = TRUE
  )
  if (!identical(checked$method, as.character(bundle$method))) {
    stop("Published MOFAcellulaR artifact method changed during validation")
  }
  list(
    dataset = source_record$dataset,
    variant = variant,
    path = path,
    md5 = checked$md5,
    size = as.numeric(file.info(path)$size),
    factor_count = as.integer(factor_count),
    sample_ids = checked$actual_sample_ids,
    dropped_sample_ids = checked$dropped_sample_ids
  )
}

.mofa_record_variant <- function(record) {
  list(
    path = record$path,
    md5 = record$md5,
    size = record$size,
    dataset = record$dataset,
    variant = record$variant,
    factor_count = record$factor_count,
    sample_ids = record$sample_ids,
    dropped_sample_ids = record$dropped_sample_ids
  )
}

.mofa_create_debug_record <- function(output_dir, run_id, package_info,
                                     config_info, source_record, artifacts) {
  if (!identical(source_record$dataset, "_debug")) {
    stop("debug-pass records may reference only _debug")
  }
  debug_provenance <- .mofa_validate_debug_sample_record(
    source_record$debug_sample_ids,
    source_record$debug_sample_count,
    label = "_debug source preflight Sample provenance"
  )
  variants <- setNames(
    lapply(.mofa_variants, function(variant) {
      candidate <- artifacts[vapply(artifacts, function(x) identical(x$variant, variant), logical(1))]
      if (length(candidate) != 1L) stop("debug pass did not produce both configured variants")
      .mofa_record_variant(candidate[[1L]])
    }),
    .mofa_variants
  )
  record <- list(
    schema_version = 1L,
    status = "PASS_PENDING_REVIEW",
    reviewed = FALSE,
    scope = "debug",
    run_id = run_id,
    package = package_info$package,
    package_repository = package_info$repository,
    package_repository_url = package_info$repository_url,
    package_sha = package_info$commit,
    package_version = package_info$version,
    config_path = config_info$path,
    config_md5 = config_info$md5,
    runner_path = .mofa_script_path(),
    runner_md5 = .mofa_md5(.mofa_script_path()),
    source = list(
      dataset = source_record$dataset,
      view = source_record$view,
      h5ad_path = source_record$h5ad_path,
      h5ad_md5 = source_record$h5ad_md5,
      h5ad_size = source_record$h5ad_size,
      sample_ids = debug_provenance$sample_ids,
      sample_count = debug_provenance$sample_count
    ),
    variants = variants,
    note = "Pass requires independent review before benchmark_union is enabled"
  )
  path <- file.path(output_dir, "mofacellular_debug_pass.json")
  .mofa_write_json_with_sidecar(record, path, output_dir)
  path
}

.mofa_load_debug_record <- function(path) {
  path <- .mofa_path(path, must_work = TRUE)
  identity <- .mofa_validate_checksum_sidecar(path)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read debug-pass records")
  }
  record <- if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
  if (!is.list(record)) stop("debug-pass record is not an object")
  list(record = record, identity = identity, path = path)
}

.mofa_validate_debug_pass_record <- function(path, package_sha, config_info,
                                             input_dir) {
  loaded <- .mofa_load_debug_record(path)
  record <- loaded$record
  required <- c("schema_version", "status", "reviewed", "scope", "package_sha",
                "config_path", "config_md5", "source", "variants")
  missing <- setdiff(required, names(record))
  if (length(missing)) stop("debug-pass record is missing: ", paste(missing, collapse = ", "))
  if (!identical(record$scope, "debug") ||
      !record$status %in% c("PASS", "DEBUG_PASS", "DEBUG_PASS_REVIEWED") ||
      !isTRUE(record$reviewed) ||
      !identical(.mofa_validate_sha(as.character(record$package_sha)), package_sha)) {
    stop("debug-pass record has not been reviewed for the requested exact package commit")
  }
  if (!identical(.mofa_path(as.character(record$config_path)), config_info$path) ||
      !identical(tolower(as.character(record$config_md5)), config_info$md5)) {
    stop("debug-pass record config identity does not match this run")
  }
  source <- record$source
  if (!is.list(source) || !identical(source$dataset, "_debug") ||
      !identical(source$view, "benchmark_analysis")) {
    stop("debug-pass record source is not the explicit _debug benchmark_analysis H5AD")
  }
  debug_entry <- config_info$config[["_debug"]]
  debug_path <- .mofa_resolve_h5ad(input_dir, "_debug", debug_entry)
  debug_identity <- .mofa_validate_checksum_sidecar(debug_path)
  if (!identical(.mofa_path(as.character(source$h5ad_path)), debug_path) ||
      !identical(tolower(as.character(source$h5ad_md5)), debug_identity$MD5) ||
      !identical(suppressWarnings(as.numeric(source$h5ad_size)),
                 as.numeric(debug_identity$SIZE))) {
    stop("debug-pass record source H5AD identity is stale or unverified")
  }
  debug_source <- .mofa_validate_debug_source(debug_path, debug_identity)
  debug_provenance <- .mofa_validate_debug_sample_record(
    source$sample_ids,
    source$sample_count,
    label = "debug-pass source Sample provenance"
  )
  if (!identical(debug_provenance$sample_ids, debug_source$sample_ids) ||
      !identical(debug_provenance$sample_count, debug_source$sample_count)) {
    stop("debug-pass record Sample provenance does not match the verified source")
  }
  if (!is.list(record$variants) ||
      !identical(sort(names(record$variants)), sort(.mofa_variants))) {
    stop("debug-pass record must contain exactly lowres and highres variants")
  }
  checked <- list()
  for (variant in .mofa_variants) {
    entry <- record$variants[[variant]]
    if (!is.list(entry) || is.null(entry$path) || is.null(entry$md5)) {
      stop("debug-pass record variant is malformed: ", variant)
    }
    artifact_path <- .mofa_path(as.character(entry$path), must_work = TRUE)
    artifact_identity <- .mofa_validate_checksum_sidecar(artifact_path)
    if (!identical(tolower(as.character(entry$md5)), artifact_identity$MD5)) {
      stop("debug-pass record variant checksum mismatch: ", variant)
    }
    checked[[variant]] <- validate_mofacellular_artifact(
      artifact_path,
      expected_dataset = "_debug",
      expected_scope = "debug",
      expected_factor_count = 2L,
      expected_source_md5 = debug_identity$MD5,
      allow_dropped = TRUE
    )
  }
  list(
    path = loaded$path,
    md5 = loaded$identity$MD5,
    package_sha = package_sha,
    source = debug_identity,
    variants = checked
  )
}

.mofa_main <- function(args) {
  allowed <- c(
    "scope", "config_path", "input_dir", "output_dir", "run_id",
    "package_sha", "mofacellular_sha", "mofacellular_commit", "package_commit",
    "debug_pass_record", "debug_pass", "help"
  )
  unknown <- setdiff(names(args), allowed)
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  if (isTRUE(args$help)) {
    cat(.mofa_usage(), "\n")
    return(invisible(NULL))
  }
  required <- c("scope", "config_path", "input_dir", "output_dir", "run_id")
  missing <- required[vapply(required, function(name) {
    is.null(args[[name]]) || identical(args[[name]], TRUE) ||
      !is.character(args[[name]]) || length(args[[name]]) != 1L ||
      is.na(args[[name]]) || !nzchar(args[[name]])
  }, logical(1))]
  if (length(missing)) stop("Missing required argument(s): ", paste0("--", missing, collapse = ", "))
  scope <- match.arg(as.character(args$scope), c("debug", "benchmark_union"))
  run_id <- .mofa_run_id(args$run_id)
  package_sha <- .mofa_get_package_sha(args)
  config_info <- .mofa_read_config(args$config_path)
  source_records <- .mofa_preflight_sources(config_info, args$input_dir, scope)
  package_info <- .mofa_validate_package_provenance(package_sha)
  debug_gate <- NULL
  debug_arg <- .mofa_or(args$debug_pass_record, args$debug_pass)
  if (identical(scope, "benchmark_union")) {
    if (is.null(debug_arg) || identical(debug_arg, TRUE) ||
        !is.character(debug_arg) || length(debug_arg) != 1L || !nzchar(debug_arg)) {
      stop("benchmark_union requires --debug_pass_record from a reviewed _debug pass")
    }
    debug_gate <- .mofa_validate_debug_pass_record(
      debug_arg, package_sha, config_info, args$input_dir
    )
  } else if (!is.null(debug_arg)) {
    stop("--debug_pass_record is accepted only for benchmark_union")
  }
  output <- .mofa_prepare_output_dir(args$output_dir, run_id)
  factor_count <- if (identical(scope, "debug")) 2L else 15L
  expected_artifact_count <- length(source_records) * length(.mofa_variants)
  if (!identical(expected_artifact_count, 2L) && identical(scope, "debug")) {
    stop("debug scope must have exactly two low/high artifact rows")
  }
  if (expected_artifact_count <= 0L) {
    stop("MOFAcellulaR scope produced no expected artifact rows")
  }
  artifacts <- list()
  for (dataset in names(source_records)) {
    source_record <- source_records[[dataset]]
    for (variant in .mofa_variants) {
      bundle <- .mofa_make_bundle(
        source_record = source_record,
        scope = scope,
        run_id = run_id,
        package_info = package_info,
        variant = variant,
        factor_count = factor_count
      )
      artifacts[[length(artifacts) + 1L]] <- .mofa_publish_bundle(
        bundle = bundle,
        source_record = source_record,
        output_dir = output$output_dir,
        scope = scope,
        run_id = run_id,
        variant = variant,
        factor_count = factor_count
      )
    }
  }
  debug_record_path <- NULL
  if (identical(scope, "debug")) {
    source_record <- source_records[["_debug"]]
    debug_record_path <- .mofa_create_debug_record(
      output_dir = output$output_dir,
      run_id = run_id,
      package_info = package_info,
      config_info = config_info,
      source_record = source_record,
      artifacts = artifacts
    )
  }
  manifest <- list(
    schema_version = 1L,
    status = "COMPLETED",
    analysis = "MOFAcellulaR",
    scope = scope,
    run_id = run_id,
    runner_path = .mofa_script_path(),
    runner_md5 = .mofa_md5(.mofa_script_path()),
    config_path = config_info$path,
    config_md5 = config_info$md5,
    package = package_info,
    debug_pass_record = debug_gate,
    debug_pass_record_path = debug_record_path,
    sources = unname(source_records),
    artifacts = lapply(artifacts, .mofa_record_variant)
  )
  manifest_path <- file.path(output$output_dir, "mofacellular_run_manifest.json")
  .mofa_write_json_with_sidecar(manifest, manifest_path, output$output_dir)
  message("Completed standalone MOFAcellulaR ", scope, " run ", run_id,
          " with ", length(artifacts), " validated artifacts")
  invisible(manifest)
}

.mofa_cli_invocation <- local({
  full <- commandArgs(trailingOnly = FALSE)
  token <- full[startsWith(full, "--file=")]
  if (!length(token)) {
    FALSE
  } else {
    candidate <- tryCatch(normalizePath(sub("^--file=", "", token[[1L]]),
                                        mustWork = TRUE),
                         error = function(error) NA_character_)
    target <- tryCatch(.mofa_script_path(), error = function(error) NA_character_)
    !is.na(candidate) && !is.na(target) && identical(candidate, target) &&
      identical(basename(candidate), "run_mofacellular.R")
  }
})

if (isTRUE(.mofa_cli_invocation)) {
  cli_result <- tryCatch(
    .mofa_main(.mofa_parse_flags(commandArgs(trailingOnly = TRUE))),
    error = function(error) {
      message("MOFAcellulaR wrapper FAILED: ", conditionMessage(error))
      quit(save = "no", status = 1L, runLast = FALSE)
    }
  )
}
