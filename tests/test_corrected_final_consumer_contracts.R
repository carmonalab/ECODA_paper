#!/usr/bin/env Rscript
# Focused no-compute barrier regression for corrected-final consumers.

script_arg <- commandArgs(trailingOnly = FALSE)[grepl(
  "^--file=", commandArgs(trailingOnly = FALSE)
)][1]
stopifnot(length(script_arg) == 1L)
root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg)), ".."),
  mustWork = TRUE
)

with_temporary <- function(code) {
  directory <- tempfile("ecoda-corrected-consumer-")
  dir.create(directory, recursive = TRUE)
  directory <- normalizePath(directory, mustWork = TRUE)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  eval(substitute(code), envir = environment())
}

with_temporary({
  config_path <- file.path(directory, "datasets.json")
  jsonlite::write_json(
    list(
      Synthetic = list(
        use_for_batch_effect = TRUE,
        columns = list(
          sample = "Sample",
          label = "label",
          batch = list("assay", "sex")
        ),
        views = list(
          batch_effect_corrected = list(
            input_file_name = "synthetic.h5ad",
            output_file_name = "synthetic-corrected.h5ad",
            subset_vars = list()
          )
        )
      )
    ),
    config_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  input_root <- file.path(directory, "input")
  h5ad_dir <- file.path(input_root, "Synthetic", "output")
  dir.create(h5ad_dir, recursive = TRUE)
  h5ad_path <- file.path(h5ad_dir, "synthetic-corrected.h5ad")
  analysis_root <- file.path(directory, "batch_effect", "corrected_final")
  metadata_dir <- file.path(analysis_root, "metadata")
  dir.create(metadata_dir, recursive = TRUE)
  metadata_path <- file.path(metadata_dir, "Synthetic_sample_metadata.feather")
  metadata_path <- normalizePath(metadata_path, mustWork = FALSE)
  selection_path <- file.path(directory, "selection.tsv")
  writeLines(
    "Synthetic\tbatch_effect_corrected\tbatch_effect_corrected",
    selection_path,
    useBytes = TRUE
  )
  report_path <- file.path(directory, "consumer-contract.json")

  fixture_python <- paste(
    "import hashlib, sys; from pathlib import Path;",
    "import anndata as ad, numpy as np, pandas as pd, pyarrow.feather as feather;",
    "sys.path.insert(0, sys.argv[2]);",
    "from src.utils.py.batch_contract import build_batch_contract_identity;",
    "h5ad = Path(sys.argv[1]); metadata = Path(sys.argv[3]);",
    "obs = pd.DataFrame({'Sample': ['s1','s2','s3','s4'],",
    "'label': ['A','A','B','B'], 'assay': ['a']*4,",
    "'sex': ['F','M','F','M']}, index=['c1','c2','c3','c4']);",
    "n = 2000; adata = ad.AnnData(X=np.ones((4,n), dtype=np.float32), obs=obs,",
    "var=pd.DataFrame({'hvg_rank': np.arange(1,n+1,dtype=float)},",
    "index=[f'g{i}' for i in range(n)]));",
    "adata.layers['counts'] = np.ones((4,n), dtype=np.int64);",
    "adata.obsm['X_pca_batch_effect_corrected_hvg2000'] = np.ones((4,2), dtype=np.float32);",
    "adata.obsm['X_pca_harmony_batch_effect_corrected_hvg2000'] = np.ones((4,2), dtype=np.float32);",
    "adata.uns['batch_contract'] = build_batch_contract_identity(['assay','sex'], sample_column='Sample', method_id='preprocess', model_id='hvg_composite_v1');",
    "adata.write_h5ad(str(h5ad)); feather.write_feather(obs, str(metadata));",
    "digest = hashlib.md5(metadata.read_bytes()).hexdigest();",
    "(metadata.parent / (metadata.name + '.md5')).write_text(f'MD5={digest}\\nSIZE={metadata.stat().st_size}\\nPATH={metadata}\\n')"
  )
  fixture_status <- system2(
    "pixi",
    c(
      "run", "-e", "default", "python", "-c",
      shQuote(fixture_python),
      shQuote(h5ad_path), shQuote(root), shQuote(metadata_path)
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  stopifnot(identical(fixture_status, 0L))

  validator <- file.path(
    root, "src", "5_run_benchmark_methods",
    "validate_corrected_final_consumer_contracts.R"
  )
  invoke <- function(output_path) {
    result <- system2(
      "pixi",
      c(
        "run", "-e", "default", "Rscript", "--vanilla", validator,
        "--config", config_path,
        "--selection", selection_path,
        "--analysis-root", analysis_root,
        "--input-root", input_root,
        "--output", output_path
      ),
      stdout = TRUE,
      stderr = TRUE,
      env = c(
        paste0("PROJECT_ROOT=", root),
        "ANALYSIS_VARIANT=corrected_final",
        "ANALYSIS_PASS=corrected",
        "ECODA_RUN_ID=consumer-contract-test"
      )
    )
    status <- attr(result, "status")
    if (is.null(status)) status <- 0L
    list(status = as.integer(status), output = result)
  }

  positive <- invoke(report_path)
  if (positive$status != 0L) {
    details <- if (file.exists(paste0(metadata_path, ".md5"))) {
      paste(readLines(paste0(metadata_path, ".md5")), collapse = "\\n")
    } else {
      "metadata sidecar is absent"
    }
    stop(
      "positive barrier failed: ", paste(positive$output, collapse = "\\n"),
      " ; ", details
    )
  }
  stopifnot(file.exists(report_path))
  report <- jsonlite::fromJSON(report_path, simplifyVector = FALSE)
  row <- report$rows[[1L]]
  stopifnot(
    identical(report$status, "CORRECTED_FINAL_CONSUMERS_VALIDATED"),
    identical(row$status, "OK"),
    identical(unname(unlist(row$estimable_batch_keys)), "sex"),
    identical(unname(unlist(row$non_estimable_batch_keys)), "assay"),
    identical(row$correction_state, "BATCH_CORRECTION"),
    grepl("batch_key_2", row$correction_formulas$composition, fixed = TRUE),
    file.exists(paste0(report_path, ".md5"))
  )

  constant_metadata <- data.frame(
    Sample = paste0("s", 1:4),
    label = c("A", "A", "B", "B"),
    assay = rep("a", 4),
    sex = rep("F", 4),
    stringsAsFactors = FALSE
  )
  arrow::write_feather(constant_metadata, metadata_path)
  digest <- tolower(unname(tools::md5sum(metadata_path)))
  writeLines(
    c(
      paste0("MD5=", digest),
      paste0("SIZE=", file.info(metadata_path)$size),
      paste0("PATH=", metadata_path)
    ),
    paste0(metadata_path, ".md5"),
    useBytes = TRUE
  )
  no_correction_path <- file.path(directory, "consumer-contract-no-correction.json")
  no_correction <- invoke(no_correction_path)
  stopifnot(no_correction$status == 0L)
  no_correction_report <- jsonlite::fromJSON(
    no_correction_path,
    simplifyVector = FALSE
  )
  stopifnot(
    identical(no_correction_report$rows[[1L]]$correction_state, "NO_CORRECTION"),
    identical(
      no_correction_report$rows[[1L]]$correction_formulas$composition,
      "NO_CORRECTION: no estimable technical batch key"
    )
  )
})

cat("corrected-final consumer barrier: OK\n")
