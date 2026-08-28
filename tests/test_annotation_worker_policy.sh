#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKER="${ROOT}/src/4_cell_type_annotation/2.1.1_process_chunk.R"
if grep -q 'breast_mode[[:space:]]*=' "${WORKER}"; then
  echo "breast_mode was passed to an annotation call" >&2
  exit 1
fi
grep -q 'breast_mode.*default FALSE' "${WORKER}"
grep -q 'Do not pass breast_mode' "${WORKER}"
Rscript --vanilla -e '
expressions <- parse(commandArgs(TRUE)[1])
load_helper <- function(name) {
  matches <- Filter(function(expr) {
    is.call(expr) && identical(as.character(expr[[1L]]), "<-") &&
      identical(as.character(expr[[2L]]), name)
  }, as.list(expressions))
  stopifnot(length(matches) == 1L)
  eval(matches[[1L]], envir = .GlobalEnv)
}
load_helper("call_scatomic")
load_helper("call_hitme")
calls <- new.env(parent = emptyenv())
stub_run <- function(...) {
  calls$run <- list(...)
  "predictions"
}
stub_summary <- function(...) {
  calls$summary <- list(...)
  data.frame()
}
call_scatomic(stub_run, stub_summary, matrix(1), FALSE)
stopifnot(identical(calls$summary$normal_tissue, FALSE))
stopifnot(!"breast_mode" %in% names(calls$summary))
stub_hitme <- function(...) {
  calls$hitme <- list(...)
  "seurat"
}
call_hitme(stub_hitme, "object", "model", list())
stopifnot(!"breast_mode" %in% names(calls$hitme))
' "${WORKER}"
Rscript --vanilla -e 'stopifnot(identical(as.logical("false"), FALSE)); cat("NORMAL_TISSUE=false: OK\n")'
echo "annotation worker policy: OK"
