# ==============================================================================
# 2.1.1_process_chunk.R — Process one chunk of samples for cell type annotation
# ==============================================================================
# Called by 2.1_run_worker.sh (pixi run --as-is Rscript --vanilla)
# Expects a single argument: <path_to_chunk_txt>
# ==============================================================================

project_root <- Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("CRITICAL Error: PROJECT_ROOT not set. Source slurm_config.sh before calling this script.")

source(file.path(project_root, "config_helper.R"))
paths <- get_pipeline_config()

# Sourcing: seurat_utils.R (get_seurat_obj_from_h5ad) — lighter than
# load_all_functions.R (no package re-imports / conflicts with scGate & HiTME).
source(file.path(project_root, "src/utils/seurat_utils.R"))

library(scGate)
library(ProjecTILs)
library(SignatuR)
library(HiTME)
library(Seurat)
library(arrow)
library(scATOMIC)
library(cutoff.scATOMIC)

# Bounded cutoff.scATOMIC::em(): upstream bug — EM loop with default t=1e-64
# (below machine epsilon) and no iteration cap never terminates on
# ill-conditioned score mixtures (automatic_threshold jitters scores with
# runif, so the mle2 M-step never bit-stabilizes). See abelson-lab/scATOMIC
# issue <NUMBER>. t=1e-6 matches flexmix's tolerance; measured on real Pelka
# scores it yields thresholds identical to the bit-stable fit (max_iter=200
# bounds pathological calls; a smaller t would only cost a few iterations).
em_bounded <- function(data, D1, D2, t = 1e-6, max_iter = 200) {
  data_name <- unlist(strsplit(deparse(match.call()), "="))[2]
  data_name <- sub(",.*$", "", gsub(" ", "", data_name))
  start <- as.list(startval(data, D1, D2))
  D1b <- hash[[D1]]
  D2b <- hash[[D2]]
  lambda0 <- 0
  iter <- 0
  with(start, {
    while (abs(lambda0 - mean(lambda)) > t) {
      iter <- iter + 1
      if (iter >= max_iter) {
        message(sprintf("cutoff.scATOMIC::em: EM loop capped at %d iterations (t=%g)", max_iter, t))
        break
      }
      lambda <- mean(lambda)
      lambda0 <- lambda
      distr1 <- lambda * D1b(data, mu1, sigma1)
      distr2 <- (1 - lambda) * D2b(data, mu2, sigma2)
      lambda <- distr1 / (distr1 + distr2)
      mLL2 <- function(mu1, sigma1, mu2, sigma2)
        return(mLL(mu1, sigma1, mu2, sigma2, lambda, data, D1b, D2b))
      start <- as.list(log(c(mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2)))
      out <- bbmle::mle2(mLL2, start, "Nelder-Mead")
      coef <- out@coef
      coef_n <- names(coef)
      names(coef) <- NULL
      for (i in 1:4) assign(coef_n[i], exp(coef[i]))
    }
    out <- list(lambda = lambda, param = exp(out@coef), D1 = D1, D2 = D2,
                deviance = out@min, data = data, data_name = data_name, out = out, t = t)
    class(out) <- "em"
    return(out)
  })
}
em_ns <- asNamespace("cutoff.scATOMIC")
environment(em_bounded) <- em_ns
unlockBinding("em", em_ns)
assignInNamespace("em", em_bounded, em_ns)
lockBinding("em", em_ns)
message("Patched cutoff.scATOMIC::em with bounded EM loop (t=1e-6, max_iter=200)")

library(R.utils)

library(reticulate)

ad <- import("anndata", convert = FALSE)
set.seed(123)

# ==============================================================================
# HELPER FUNCTIONS & PARAMETER PARSING
# ==============================================================================

env_sample_col    <- Sys.getenv("SAMPLE_COLNAME")
env_tissue        <- Sys.getenv("TISSUE_TYPE")
env_normal_tissue <- Sys.getenv("NORMAL_TISSUE")

defaults <- list(
  chunk_file            = NULL,
  sample_colname        = if (env_sample_col != "") env_sample_col else "Sample",
  tissue_type           = if (env_tissue != "") env_tissue else "Tumor",
  normal_tissue         = if (env_normal_tissue != "") as.logical(env_normal_tissue) else TRUE
)

# Every checkpoint and final feather carries the complete dual-method schema.
# Missing upstream labels are represented as typed NA values, not omitted
# columns; the final dataset-level anchor check decides whether a method
# produced anything anywhere in the selected union.
HITME_REQUIRED_COLS <- c("layer1", "layer2", "layer3")
SCATOMIC_REQUIRED_COLS <- c(
  "layer_1", "layer_2", "layer_3", "layer_4", "layer_5", "layer_6",
  "scATOMIC_pred", "classification_confidence", "S.Score", "G2M.Score", "Phase"
)
HITME_OPTIONAL_COLS <- c(
  "IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell", "cellCycle.G2M_UCell"
)
HITME_OUTPUT_COLS <- c(HITME_OPTIONAL_COLS, HITME_REQUIRED_COLS)
SCATOMIC_OUTPUT_COLS <- SCATOMIC_REQUIRED_COLS
ANNOT_REQUIRED_COLS <- c(HITME_REQUIRED_COLS, SCATOMIC_REQUIRED_COLS)
ANNOT_OUTPUT_COLS <- c(HITME_OUTPUT_COLS, SCATOMIC_OUTPUT_COLS)
NUMERIC_ANNOT_COLS <- c(
  "IFN_UCell", "HeatShock_UCell", "cellCycle.G1S_UCell",
  "cellCycle.G2M_UCell", "classification_confidence", "S.Score", "G2M.Score"
)
# scATOMIC breast_mode is intentionally left at the upstream default FALSE for
# this uniform workflow. Do not pass breast_mode to run_scATOMIC() or
# create_summary_matrix(); tissue-specific mode would make cohorts incomparable.
raw_args <- commandArgs(trailingOnly = TRUE)
args <- defaults
if (length(raw_args) > 0) args$chunk_file <- raw_args[1]
# No breast_mode variable is derived here; upstream calls intentionally use FALSE.

# ============================= ANNOTATION SAFETY ============================
# Phase-5 T6: a method that annotated 0 cells, <2 unique cell types, or only
# NAs for a sample must NOT crash the worker (new tissues - brain/heart/
# kidney/pancreas - have no scGate/HiTME/scATOMIC models). Warn + keep the
# column NA/unclassified (the existing column-aligned chunk assembly handles
# the NA path), and record per-method stats (n cells, n types) into a
# per-chunk stats feather (output/annotation_stats_chunk_<N>.feather; never
# synced to NAS -- excluded from both rsyncs, mirroring annotation_tmp/). Keep
# the cutoff.scATOMIC::em patch and wall-time guards intact.
METHOD_COLS <- list(
  HiTME = HITME_REQUIRED_COLS,
  scATOMIC = SCATOMIC_REQUIRED_COLS
)

chunk_annot_stats <- list()

record_method_stats <- function(annot, sample_label) {
  for (method in names(METHOD_COLS)) {
    for (col in METHOD_COLS[[method]]) {
      if (!col %in% colnames(annot)) next
      vals <- annot[[col]]
      text_vals <- as.character(vals)
      valid_vals <- !is.na(vals) & nzchar(trimws(text_vals))
      n_annotated <- sum(valid_vals)
      n_types <- if (n_annotated > 0) length(unique(text_vals[valid_vals])) else 0
      chunk_annot_stats[[length(chunk_annot_stats) + 1]] <<- data.frame(
        sample = sample_label,
        method = method,
        column = col,
        n_cells = nrow(annot),
        n_annotated = n_annotated,
        n_types = n_types,
        stringsAsFactors = FALSE
      )
      if (n_annotated == 0) {
        warning(sprintf(
          "ANNOTATION SAFETY: %s column '%s' for sample %s is entirely NA (0/%d cells annotated); keeping NA (unclassified), no crash",
          method, col, sample_label, nrow(annot)
        ))
      } else if (n_types < 2) {
        warning(sprintf(
          "ANNOTATION SAFETY: %s column '%s' for sample %s has %d/%d cells annotated with only %d unique type(s) (< 2); annotation not informative for this tissue",
          method, col, sample_label, n_annotated, nrow(annot), n_types
        ))
      } else {
        message(sprintf(
          "ANNOTATION STATS: sample %s, %s (col '%s'): %d/%d cells, %d types",
          sample_label, method, col, n_annotated, nrow(annot), n_types
        ))
      }
    }
  }
}

typed_na <- function(column, n) {
  if (column %in% NUMERIC_ANNOT_COLS) return(rep(NA_real_, n))
  rep(NA_character_, n)
}

materialize_annotation <- function(meta, target_sample) {
  n <- nrow(meta)
  output_cols <- c(
    ANNOT_REQUIRED_COLS,
    intersect(HITME_OPTIONAL_COLS, colnames(meta))
  )
  annot <- list()
  for (column in output_cols) {
    annot[[column]] <- if (column %in% colnames(meta)) {
      meta[[column]]
    } else {
      typed_na(column, n)
    }
  }
  annot <- as.data.frame(annot, stringsAsFactors = FALSE, check.names = FALSE)
  annot$cell_barcode <- rownames(meta)
  annot$Sample <- rep(as.character(target_sample), n)
  annot[, c("Sample", "cell_barcode", output_cols), drop = FALSE]
}
call_scatomic <- function(run_fn, summary_fn, counts, normal_tissue) {
  predictions <- run_fn(counts)
  summary_fn(
    prediction_list = predictions,
    raw_counts = counts,
    normal_tissue = normal_tissue
  )
}

call_hitme <- function(run_fn, seurat_obj, scgate_model, ref_maps) {
  run_fn(
    object = seurat_obj,
    scGate.model = scgate_model,
    ref.maps = ref_maps,
    verbose = FALSE,
    ncores = 1
  )
}


checkpoint_valid <- function(path, target_sample) {
  tryCatch({
    if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
    checkpoint <- read_feather(path)
    needed <- c("Sample", "cell_barcode", ANNOT_REQUIRED_COLS)
    if (!all(needed %in% colnames(checkpoint)) || nrow(checkpoint) == 0) {
      return(FALSE)
    }
    samples <- as.character(checkpoint$Sample)
    barcodes <- as.character(checkpoint$cell_barcode)
    if (any(is.na(samples)) || any(is.na(barcodes)) ||
        any(!nzchar(trimws(samples))) || any(!nzchar(trimws(barcodes)))) return(FALSE)
    if (anyDuplicated(paste(samples, barcodes, sep = "\t"))) return(FALSE)
    for (column in intersect(NUMERIC_ANNOT_COLS, colnames(checkpoint))) {
      values <- checkpoint[[column]]
      numeric_values <- suppressWarnings(as.numeric(values))
      invalid <- !is.na(values) & is.na(numeric_values)
      if (any(invalid) || any(!is.finite(numeric_values[!is.na(numeric_values)]))) return(FALSE)
    }
    all(samples == as.character(target_sample))
  }, error = function(e) FALSE)
}

write_checksum_sidecar <- function(path) {
  digest <- unname(tools::md5sum(path))
  size <- file.info(path)$size
  if (is.na(digest) || is.na(size) || size <= 0) {
    stop("Cannot checksum annotation feather: ", path)
  }
  sidecar <- paste0(path, ".md5")
  sidecar_tmp <- paste0(sidecar, ".tmp.", Sys.getpid())
  writeLines(c(
    paste0("MD5=", digest),
    paste0("SIZE=", size),
    paste0("PATH=", path)
  ), sidecar_tmp)
  if (!file.rename(sidecar_tmp, sidecar)) {
    unlink(sidecar_tmp)
    stop("Could not atomically install annotation feather sidecar: ", sidecar)
  }
}

if (is.null(args$chunk_file) || !file.exists(args$chunk_file)) {
  stop("Valid 'chunk_file' parameter not parsed from execution context!")
}

# Wall-clock budget: SLURM_TIME_LIMIT is in MINUTES (Slurm sets it as the
# --time value, e.g. "120" for a 2h job) — multiply by 60 for the seconds
# scale. The 7200 s fallback matches the worker's #SBATCH --time=02:00:00.
# SLURM_JOB_START_TIME is not a standard Slurm env var, so proc.time()[3]
# (R launch, i.e. after staging) is the actual elapsed source; it under-counts
# wall time by the staging duration, which is a safe (conservative) direction.
wall_limit_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_TIME_LIMIT", "")))
if (is.na(wall_limit_s) || wall_limit_s <= 0) {
  wall_limit_s <- 7200
} else {
  wall_limit_s <- wall_limit_s * 60
}
job_start_s <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_START_TIME", "")))
wall_elapsed <- function() if (!is.na(job_start_s)) as.numeric(Sys.time()) - job_start_s else proc.time()[3]
wall_left <- function() wall_limit_s - wall_elapsed()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load ref maps ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Load scGate models ####
# The Stage 4 preflight creates and validates this cache exactly once. Workers
# must never download independently: parallel fallback downloads race the
# release contract and can produce different model versions.
scgate_db_path <- Sys.getenv("SCGATE_DB_PATH", unset = file.path(project_root, "aux", "scGateDB.rds"))
if (scgate_db_path == "" || !file.exists(scgate_db_path) ||
    file.info(scgate_db_path)$size <= 0) {
  stop("Validated scGate DB cache is required before annotation workers: ", scgate_db_path)
}
scGate_models_DB <- readRDS(scgate_db_path)
scGate_models_blood <- scGate_models_DB$human$PBMC
scGate_models_blood$MoMac <- scGate_models_blood$Monocyte
scGate_models_blood$Monocyte <- NULL
scGate_models_tumor <- scGate_models_DB$human$HiTME
scGate_models_blood <- c(scGate_models_blood, scGate_models_tumor[!names(scGate_models_tumor) %in% names(scGate_models_blood)])

### Load ProjecTILs ref maps ####
ref.maps_sketched <- list(
  CD8 = load.reference.map(file.path(paths$path_ref, "sketched_CD8T_human_ref_v1.rds")),
  CD4 = load.reference.map(file.path(paths$path_ref, "sketched_CD4T_human_ref_v2.rds")),
  DC = load.reference.map(file.path(paths$path_ref, "sketched_DC_human_ref_v2.rds")),
  MoMac = load.reference.map(file.path(paths$path_ref, "sketched_MoMac_human_v1.rds"))
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Process data ######
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunk_lines <- readLines(args$chunk_file)
h5ad_file <- chunk_lines[1]
samples_to_process <- chunk_lines[-1]

### Tap h5ad file ####
adata <- ad$read_h5ad(h5ad_file, backed = "r")
obs <- py_to_r(adata$obs)

# Hardening: warn (not stop) if the on-disk X format is not CSR. anndata only
# overrides selective row-indexing for backed CSR matrices; a backed CSC matrix
# falls back to a full in-memory read per-sample subset (silent OOM).
# Files produced by 1.1.1_preprocess.py are CSR by construction, so a non-CSR
# file here means preprocessing output was replaced/regenerated elsewhere.
x_format <- tryCatch(py_to_r(adata$X$format), error = function(e) NULL)
if (!is.null(x_format) && x_format != "csr") {
  warning("On-disk X format is '", x_format, "' (expected 'csr'); backed per-sample ",
          "subsetting will fully load the matrix into memory. Re-run preprocessing ",
          "(1.1.1_preprocess.py) to force CSR on-disk.")
}

if (!args$sample_colname %in% colnames(obs)) {
  stop(args$sample_colname, " not found in h5ad obs colnames!")
}

# Wipe ALL pre-existing obs columns except the sample column: annotation must
# start from a minimal Seurat object. Preprocessing keeps every source-rds obs
# column, so legacy annotation columns (scGate_*, functional.cluster, *_UCell,
# layer_1..6, ...) would otherwise leak into the annotation object via
# get_seurat_obj_from_h5ad's meta.data and silently skip scATOMIC (layer_1
# guard) or confuse Run.HiTME. No pattern lists needed — nothing but the
# sample column may enter annotation.
wiped_cols <- setdiff(colnames(obs), args$sample_colname)
if (length(wiped_cols) > 0) {
  message(sprintf(
    "Wiping %d pre-existing obs columns (%s); building minimal Seurat objects (meta.data = %s only)",
    length(wiped_cols), paste(wiped_cols, collapse = ", "), args$sample_colname
  ))
}
obs <- obs[, args$sample_colname, drop = FALSE]

# ==============================================================================
# PROCESSING LOOP
# ==============================================================================
message(paste("--- Starting processing for chunk file:", args$chunk_file, "---"))

# Feather name is derived from the chunk file (chunk_<N>.txt -> annotations_
# chunk_<N>.feather), not from the global Slurm task id.
annot_file <- file.path(
  paths$path_output,
  paste0("annotations_", sub("\\.txt$", ".feather", basename(args$chunk_file)))
)
final_feather_valid <- FALSE
if (file.exists(annot_file) && file.info(annot_file)$size > 0) {
  validator <- file.path(project_root, "src", "utils", "py", "annotation_contract.py")
  python_bin <- Sys.getenv("PYTHON_BIN")
  if (python_bin != "" && file.exists(validator)) {
    final_feather_valid <- identical(
      system2(
        python_bin,
        c(validator, "--path", annot_file, "--require-sidecar"),
        stdout = FALSE,
        stderr = FALSE
      ),
      0L
    )
  }
  if (final_feather_valid) {
    message(paste("Chunk already processed and validated:", annot_file))
  } else {
    warning("Existing annotation feather failed schema/checksum validation; rebuilding: ",
            annot_file)
    unlink(c(annot_file, paste0(annot_file, ".md5")))
  }
}
if (!final_feather_valid) {
  # Per-sample checkpoints are run-owned because paths$path_output is supplied
  # by the selected Stage 4 run.
  tmp_dir <- file.path(paths$path_output, "annotation_tmp",
                       sub("\\.txt$", "", basename(args$chunk_file)))
  # does NOT convert; materialize it to a Python list first so py_to_r returns
  # an R character vector usable with %in%. Computed once per chunk — the
  # layers group never changes between samples.
  layer_keys <- py_to_r(import_builtins(convert = FALSE)$list(adata$layers$keys()))
  counts_layer <- if ("counts" %in% layer_keys) "counts" else "X"
  if (counts_layer == "X") {
    # Annotation-union files carry the raw counts in X by design (minimal
    # layout: no layers group — see 1.1_prepare_chunks.py), so this is the
    # designed primary path for unions, not a warning case. Preprocessed
    # view files still carry layers["counts"] and take the "counts" branch.
    message("Union carries counts in X by design (no 'counts' layer in ", h5ad_file,
            "); using X as counts input for scATOMIC/HiTME.")
  }

  for (i in seq_along(samples_to_process)) {
    target_sample <- samples_to_process[i]
    sample_tmp <- file.path(tmp_dir, sprintf("sample_%02d.feather", i))
    if (file.exists(sample_tmp)) {
      if (checkpoint_valid(sample_tmp, target_sample)) {
        annot <- read_feather(sample_tmp)
        record_method_stats(annot, target_sample)
        message(paste("resume: sample", target_sample, "checkpoint validated; skipping"))
        next
      }
      warning("Existing sample checkpoint failed validation; rebuilding: ", sample_tmp)
      unlink(sample_tmp)
    }
    seurat_obj <- tryCatch({
      obj <- get_seurat_obj_from_h5ad(
        adata, obs, target_sample,
        sample_colname = args$sample_colname,
        counts_layer = counts_layer
      )
      message("  Adding log-normalized 'data' layer (NormalizeData) for HiTME...")
      NormalizeData(obj)
    }, error = function(e) {
      stop("structural annotation setup failed for sample ", target_sample,
           ": ", conditionMessage(e))
    })
    timeout <- max(60, ncol(seurat_obj) / 10000 * 60 * 10)
    attempt_timeout <- min(timeout, 1800)

    # Method calls are fault-tolerant: timeout/error leaves the method's
    # required columns absent so materialize_annotation fills typed NAs.
    # Setup and output materialization remain outside this boundary and
    # therefore fail the chunk on structural corruption.
    if (is.null(seurat_obj@meta.data[["layer_1"]])) {
      for (a in 1:5) {
        eff <- min(attempt_timeout, wall_left() - 300)
        if (eff < 60) {
          message(paste("  scATOMIC: insufficient wall time remaining (", round(wall_left()), "s); skipping attempt"))
          break
        }
        message(paste("  scATOMIC attempt", a, "with", round(eff), "s timeout"))
        result <- tryCatch({
          withTimeout({
            sca_results <- call_scatomic(
              run_scATOMIC,
              create_summary_matrix,
              seurat_obj@assays$RNA$counts,
              args$normal_tissue
            )
            "Complete"
          }, timeout = eff)
        }, TimeoutException = function(te) {
          message("  scATOMIC timeout, retrying...")
          NULL
        }, error = function(er) {
          message(paste("  scATOMIC error:", er$message, "- retrying..."))
          NULL
        })
        if (!is.null(result)) {
          added <- tryCatch({
            sca_cols <- intersect(SCATOMIC_OUTPUT_COLS, colnames(sca_results))
            if (length(sca_cols) > 0) {
              seurat_obj <- AddMetaData(seurat_obj, sca_results[, sca_cols, drop = FALSE])
            }
            TRUE
          }, error = function(er) {
            message(paste("  scATOMIC output materialization error:", er$message, "- retrying..."))
            FALSE
          })
          if (isTRUE(added)) break
        }
      }
    }

    ### HiTME annotation ####
    for (a in 1:5) {
      eff <- min(attempt_timeout, wall_left() - 300)
      if (eff < 60) {
        message(paste("  HiTME: insufficient wall time remaining (", round(wall_left()), "s); skipping attempt"))
        break
      }
      message(paste("  HiTME attempt", a, "with", round(eff), "s timeout"))
      result <- tryCatch({
        withTimeout({
          seurat_obj <- call_hitme(
            Run.HiTME,
            seurat_obj,
            if (args$tissue_type == "Blood") scGate_models_blood else scGate_models_tumor,
            ref.maps_sketched
          )
        }, timeout = eff)
      }, TimeoutException = function(te) {
        message("  HiTME timeout, retrying...")
        NULL
      }, error = function(er) {
        message(paste("  HiTME error:", er$message, "- retrying..."))
        NULL
      })
      if (!is.null(result)) break
    }

    meta <- seurat_obj@meta.data
    annot <- materialize_annotation(meta, target_sample)

    # Phase-5 T6 safety: 0-annotated/<2-types/all-NA method results warn +
    # record stats instead of propagating silently.
    record_method_stats(annot, target_sample)

    rm(seurat_obj)
    gc()
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    sample_tmp_write <- paste0(sample_tmp, ".tmp.", Sys.getpid())
    write_feather(annot, sample_tmp_write)
    if (!file.rename(sample_tmp_write, sample_tmp)) {
      unlink(sample_tmp_write)
      stop("Could not atomically install annotation checkpoint: ", sample_tmp)
    }
    message(paste("  Checkpoint written:", sample_tmp))

  }


  annotations_list <- list()
  for (i in seq_along(samples_to_process)) {
    annotations_list[[samples_to_process[i]]] <- read_feather(
      file.path(tmp_dir, sprintf("sample_%02d.feather", i))
    )
  }
  # Column-aligned assembly: a sample may legitimately lack annotation columns
  # (e.g., scATOMIC/HiTME produced nothing for it, or a stale checkpoint from an
  # older run); fill missing columns with NA instead of crashing the plain rbind
  # on heterogeneous checkpoint schemas (observed on resume: "numbers of columns
  # of arguments do not match", CombinedPBMC chunk_28 / Zhang chunk_13).
  all_cols <- unique(unlist(lapply(annotations_list, colnames)))
  annotations_list <- lapply(annotations_list, function(df) {
    missing <- setdiff(all_cols, colnames(df))
    if (length(missing) > 0) df[missing] <- NA
    df
  })
  annotations_df <- do.call(rbind, annotations_list)
  rownames(annotations_df) <- NULL
  annot_file_tmp <- paste0(annot_file, ".tmp.", Sys.getpid())
  write_feather(annotations_df, annot_file_tmp)
  if (!file.rename(annot_file_tmp, annot_file)) {
    unlink(annot_file_tmp)
    stop("Could not atomically install annotation feather: ", annot_file)
  }
  write_checksum_sidecar(annot_file)
  unlink(tmp_dir, recursive = TRUE)
  message(paste("Wrote annotations to:", annot_file))

  # Per-method annotation stats (Phase-5 T6): one row per sample; feeds the
  # post-run annotation-rate documentation step
  # (notebooks/dataset_onboarding/annotation_summary.json). The file name
  # must NOT match the annotations_chunk_*.feather globs (merge/coverage) and
  # is excluded from the canonical Stage 4 synchronization path.
  if (length(chunk_annot_stats) > 0) {
    stats_df <- do.call(rbind, chunk_annot_stats)
    rownames(stats_df) <- NULL
    stats_file <- sub(
      "^annotations_chunk_(.*)\\.feather$",
      "annotation_stats_chunk_\\1.feather",
      basename(annot_file)
    )
    stats_file <- file.path(dirname(annot_file), stats_file)
    stats_tmp <- paste0(stats_file, ".tmp.", Sys.getpid())
    write_feather(stats_df, stats_tmp)
    if (!file.rename(stats_tmp, stats_file)) {
      unlink(stats_tmp)
      stop("Could not atomically install annotation stats: ", stats_file)
    }
    message(paste("Wrote per-method annotation stats to:", stats_file))
  }
}

message(paste("---", args$chunk_file, "processing complete! ---"))
