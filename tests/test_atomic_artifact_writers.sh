#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-atomic-writers.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
pixi run Rscript --vanilla -e 'source(commandArgs(TRUE)[1]); d<-commandArgs(TRUE)[2]; dir.create(d,recursive=TRUE); r<-file.path(d,"result.rds"); save_rds_atomic(list(scores=1),r); stopifnot(artifact_checksum_ok(r)); old<-tools::md5sum(r); bad<-try(save_rds_atomic(list(),r),silent=TRUE); stopifnot(inherits(bad,"try-error"),artifact_checksum_ok(r),identical(tools::md5sum(r),old)); e<-file.path(d,"execution.feather"); write_feather_atomic(data.frame(dataset="Adams",method="m",time_secs=1,mem_GB=1),e); stopifnot(artifact_checksum_ok(e)); olde<-tools::md5sum(e); bad2<-try(write_feather_atomic(data.frame(),e),silent=TRUE); stopifnot(inherits(bad2,"try-error"),artifact_checksum_ok(e),identical(tools::md5sum(e),olde));' "${ROOT}/src/5_run_benchmark_methods/benchmark_hpc_utils.R" "${TMP_DIR}"
echo "atomic RDS/Feather writers: OK"
