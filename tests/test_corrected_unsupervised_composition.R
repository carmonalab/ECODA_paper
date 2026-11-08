# The two corrected composition arms share one real counts/CLR computation,
# while only the second applies sample-level technical limma.
root <- normalizePath(".")
Sys.setenv(PROJECT_ROOT = root,
           RETICULATE_PYTHON = file.path(root, ".pixi/envs/default/bin/python"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("src/utils/seurat_utils.R")
source("src/utils/math_utils.R")
source("src/utils/pseudobulk.R")
source("src/5_run_benchmark_methods/benchmark_hpc_utils.R")
source("src/5_run_benchmark_methods/benchmark_methods_r.R")
source("src/5_run_benchmark_methods/benchmark_pipeline.R")

assert_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
ids <- paste0("sample_", 1:12)
metadata <- data.frame(Sample = ids, batch = rep(c("A", "B"), each = 6),
                       biology = rep(c("case", "control"), 6))
key <- "leiden_res_2_batch_effect_corrected_hvg2000_harmony"
obs <- do.call(rbind, lapply(seq_along(ids), function(i) {
  data.frame(Sample = ids[[i]], batch = metadata$batch[[i]],
             biology = metadata$biology[[i]],
             author = rep(c("T", "B", "M"), c(3 + i %% 3, 5 + i %% 2, 7)),
             stringsAsFactors = FALSE)
}))
obs[[key]] <- as.character(match(obs$author, c("T", "B", "M")) - 1L)
obs[["RNA_snn_res.2"]] <- obs[[key]]
labels <- setNames(factor(metadata$biology), ids)
validation <- ecoda_batch_validate_metadata(obs, list("batch"), "Sample", "biology")

# Keep the shared scorer boundary observable without testing its algorithms.
calc_sep_score <- function(dist_mat, labels, ...) {
  stopifnot(identical(attr(dist_mat, "Labels"), names(labels)))
  as.list(setNames(rep(mean(dist_mat), 5),
    c("anosim_score", "mod_knn3_score", "cluster_score", "sil_score", "lisi_score")))
}
source_calls <- 0L
limma_inputs <- list()
process_original <- process_coda_fig
process_coda_fig <- function(..., ct_col) {
  if (identical(ct_col, key)) source_calls <<- source_calls + 1L
  process_original(..., ct_col = ct_col)
}
limma_original <- correct_clr_batch_limma
correct_clr_batch_limma <- function(feat_mat, ...) {
  limma_inputs[[length(limma_inputs) + 1L]] <<- feat_mat
  limma_original(feat_mat, ...)
}
exec_time <- function(expr) { force(expr); as.difftime(1, units = "secs") }
peak_rss_gb <- function() 0
log_exec_row <- function(...) invisible(NULL)
saved <- list()
save_rds_atomic <- function(object, file, ...) {
  saved[[basename(file)]] <<- object
}
artifact_checksum_ok <- function(...) FALSE
run <- function(data = obs) run_composition_methods_hpc(
  labels, metadata, NULL, NULL, data, "biology",
  ct_col_high_res = "author", results_dir = tempdir(), ds = "Fixture",
  batch_mode = TRUE, result_stem = "Fixture_corrected",
  corrected = TRUE, batch_keys = list("batch"), metadata_validation = validation
)
result <- run()
base <- result$ecoda_unsupervised_harmony
corrected <- result$ecoda_unsupervised_harmony_limma
stopifnot(
  identical(names(result), c("ECODA_authors_HR", "ECODA_authors_HR_NULL",
    "ecoda_unsupervised_harmony", "ecoda_unsupervised_harmony_limma", "batch_contract")),
  source_calls == 1L, length(limma_inputs) == 2L, # author HR and second arm
  identical(limma_inputs[[2L]], base$feat_mat),
  identical(base$counts, corrected$counts),
  identical(rownames(base$feat_mat), ids),
  identical(rownames(corrected$feat_mat), ids),
  identical(base$labels, corrected$labels),
  !isTRUE(all.equal(base$feat_mat, corrected$feat_mat)),
  identical(base$batch_contract$model_id, "harmony_leiden_clr_v1"),
  identical(corrected$batch_contract$model_id, "harmony_leiden_clr_limma_v1"),
  identical(base$batch_contract$fixed_effect_aliases, character()),
  !base$composition_history$sample_level_limma,
  corrected$composition_history$sample_level_limma,
  base$exec_time == 1, corrected$exec_time == 2,
  identical(base$composition_source_time_secs, corrected$composition_source_time_secs),
  identical(result$ECODA_authors_HR$feat_mat, result$ECODA_authors_HR_NULL$feat_mat),
  identical(result$ECODA_authors_HR_NULL$labels, shuffle_labels_deterministic(labels)),
  all(c("sil_score", "lisi_score") %in% names(base$scores)),
  all(c("sil_score", "lisi_score") %in% names(corrected$scores)),
  !("Fixture_corrected_ECODA_seuratres_2.rds" %in% names(saved))
)

# Exercise the actual publication validator with the newly produced bundle.
local({
  directory <- tempfile("harmony-composition-contract-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  config <- file.path(directory, "datasets.json")
  jsonlite::write_json(list(Fixture = list(
    columns = list(sample = "Sample", label = "biology", batch = list("batch")),
    views = list(batch_effect_corrected = list(output_file_name = "fixture.h5ad"))
  )), config, auto_unbox = TRUE)
  check <- function(value, method = "composition", final = FALSE, config_file = config) {
    output_dir <- if (final) file.path(directory, "batch_effect/corrected_final/results") else directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(output_dir, paste0("Fixture_batch_effect_corrected",
      if (final) "_final" else "", "_", method, ".rds"))
    saveRDS(value, path)
    writeLines(c(paste0("MD5=", unname(tools::md5sum(path))),
      paste0("SIZE=", file.info(path)$size), paste0("PATH=", path)), paste0(path, ".md5"))
    output <- system2(file.path(root, ".pixi/envs/default/bin/Rscript"),
      c("--vanilla", "src/5_run_benchmark_methods/validate_benchmark_rds_contract.R",
        "--artifact", shQuote(path), "--method", method, "--dataset", "Fixture",
        "--view", "batch_effect_corrected", "--batch-pass", "corrected",
        if (final) c("--analysis-variant", "corrected_final"),
        "--config", shQuote(config_file)), stdout = TRUE, stderr = TRUE)
    list(status = if (is.null(attr(output, "status"))) 0L else attr(output, "status"),
         output = paste(output, collapse = "\n"))
  }
  checked <- check(result)
  if (checked$status != 0L) stop(checked$output)
  wrong <- result
  wrong$ecoda_unsupervised_harmony$composition_history$sample_level_limma <- TRUE
  checked <- suppressWarnings(check(wrong))
  stopifnot(checked$status != 0L, grepl("mismatched source/correction history", checked$output))
  for (method in names(result)[names(result) != "batch_contract"]) {
    checked <- check(result[[method]], method, final = TRUE)
    if (checked$status != 0L) stop(checked$output)
  }
  wrong <- base
  wrong$composition_history$sample_level_limma <- TRUE
  checked <- suppressWarnings(check(wrong, "ecoda_unsupervised_harmony", final = TRUE))
  stopifnot(checked$status != 0L, grepl("mismatched source/correction history", checked$output))
  wrong <- base; wrong$batch_contract <- corrected$batch_contract
  checked <- suppressWarnings(check(wrong, "ecoda_unsupervised_harmony", final = TRUE))
  stopifnot(checked$status != 0L)
  wrong_config <- file.path(directory, "wrong-keys.json")
  jsonlite::write_json(list(Fixture = list(
    columns = list(sample = "Sample", label = "biology", batch = list("other_batch")),
    views = list(batch_effect_corrected = list(output_file_name = "fixture.h5ad"))
  )), wrong_config, auto_unbox = TRUE)
  checked <- suppressWarnings(check(base, "ecoda_unsupervised_harmony", final = TRUE,
                                    config_file = wrong_config))
  stopifnot(checked$status != 0L)
  checked <- suppressWarnings(check(result, "ecoda_unsupervised_harmony", final = TRUE))
  stopifnot(checked$status != 0L)
})

# All-cached and a second-arm-only miss reuse the verified shared CLR source.
artifact_checksum_ok <- function(path, ...) basename(path) %in% names(saved)
read_rds_checked <- function(path, ...) saved[[basename(path)]]
source_calls <- 0L; limma_inputs <- list()
cached <- run()
stopifnot(source_calls == 0L, length(limma_inputs) == 0L,
          identical(cached$ecoda_unsupervised_harmony$feat_mat, base$feat_mat))
saved$Fixture_corrected_ecoda_unsupervised_harmony_limma.rds <- NULL
repaired <- run()
stopifnot(source_calls == 0L, length(limma_inputs) == 1L,
          identical(limma_inputs[[1L]], base$feat_mat))
saved$Fixture_corrected_ecoda_unsupervised_harmony_limma.rds$counts[[1L]][[1L]] <-
  saved$Fixture_corrected_ecoda_unsupervised_harmony_limma.rds$counts[[1L]][[1L]] + 1
assert_error(run(), "disagree on shared counts/labels")
saved$Fixture_corrected_ecoda_unsupervised_harmony_limma.rds <-
  repaired$ecoda_unsupervised_harmony_limma
saved$Fixture_corrected_ecoda_unsupervised_harmony.rds$composition_history$sample_level_limma <- TRUE
assert_error(run(), "mismatched source/correction history")
saved$Fixture_corrected_ecoda_unsupervised_harmony.rds <- base
saved$Fixture_corrected_ecoda_unsupervised_harmony.rds$batch_contract <- corrected$batch_contract
assert_error(run(), "mismatched corrected batch contract")
bad <- obs; bad[[key]][[1L]] <- "wrong cluster"
assert_error(run(bad), "exact Harmony Leiden source")
bad <- obs; bad[[key]] <- NULL
assert_error(run(bad), "exact Harmony Leiden source")
bad <- obs; bad[[key]][[1L]] <- NA; bad[["RNA_snn_res.2"]][[1L]] <- NA
assert_error(run(bad), "cover every cell")
cat("PASS: shared Harmony counts/CLR, zero versus one sample limma, five-score boundary, strict histories/cache reuse\n")
