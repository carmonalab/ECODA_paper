#!/bin/bash
# Finite, no-fit metadata projection of an already reviewed analysis control.
# One CPU job through the existing immutable controller/runtime worker pattern.
set -euo pipefail
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && "${ECODA_SOURCE_ROOT:-}" == /* ]] || exit 1
# Slurm executes its spool copy; recover the pinned source before configuration.
SCRIPT_DIR="${ECODA_SOURCE_ROOT%/}/src/utils/bash"
source "${SCRIPT_DIR}/../../slurm_config.sh"
source "${SCRIPT_DIR}/ecoda_run_common.sh"
source "${SCRIPT_DIR}/ecoda_runtime.sh"
MODE=controller
if [[ "${1:-}" == --worker || "${1:-}" == --validate ]]; then MODE="${1#--}"; shift; fi
SELF="${SCRIPT_DIR}/project_analysis_release.sh"
CONTROL=""; CONTROL_MD5=""; OPERATION=projection
while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || exit 1
  case "$1" in
    --analysis-release) CONTROL="$2" ;;
    --analysis-release-md5) CONTROL_MD5="$2" ;;
    --operation) OPERATION="$2" ;;
    *) echo 'Usage: project_analysis_release.sh --analysis-release ABS --analysis-release-md5 MD5 [--operation projection|render-pca]' >&2; exit 1 ;;
  esac
  shift 2
done
[[ "${ECODA_SOURCE_SNAPSHOT_REQUIRED:-0}" == 1 && -n "${ECODA_RUN_ID:-}" ]] || exit 1
ecoda_validate_run_id "${ECODA_RUN_ID}"
ecoda_require_source_script_path "${SELF}" "${ECODA_SOURCE_ROOT}" >/dev/null
SCRATCH="$(ecoda_realpath_existing "${HPC_SCRATCH_DIR}")"
RUN="${SCRATCH}/_ecoda_runs/${ECODA_RUN_ID}"
[[ -z "${ECODA_RUN_ROOT:-}" || "${ECODA_RUN_ROOT}" == "${RUN}" ]] || exit 1
if [[ "${MODE}" == controller ]]; then
  [[ ! -e "${RUN}" && ! -L "${RUN}" && "${CONTROL}" == "${SCRATCH}/_ecoda_inputs/"* &&
     "${CONTROL_MD5}" =~ ^[a-f0-9]{32}$ && -f "${CONTROL}" && ! -L "${CONTROL}" &&
     "$(ecoda_realpath_existing "${CONTROL}")" == "${CONTROL}" &&
     "$(ecoda_md5_file "${CONTROL}")" == "${CONTROL_MD5}" ]] || exit 1
  [[ "${OPERATION}" == projection || "${OPERATION}" == render-pca ]] || exit 1
  ecoda_init_run stage5 "${ECODA_RUN_ID}" >/dev/null
  ecoda_atomic_write "${RUN}/manifests/operation" "${OPERATION}\n"
else
  [[ -z "${CONTROL}${CONTROL_MD5}" && "${OPERATION}" == projection ]] || exit 1
  ecoda_open_run "${ECODA_RUN_ID}" stage5
  OPERATION="$(cat "${RUN}/manifests/operation")"
  [[ "${OPERATION}" == projection || "${OPERATION}" == render-pca ]] || exit 1
fi
export ECODA_RUN_ROOT="${RUN}" RETICULATE_PYTHON="${PYTHON_BIN}" ECODA_RUNTIME_PROFILE=stage5
OUTPUT="${RUN}/artifacts/analysis_projection"
PRODUCER=analysis_metadata_projection
if [[ "${OPERATION}" == render-pca ]]; then OUTPUT="${RUN}/artifacts/benchmark_pca"; PRODUCER=analysis_pca_render; fi
PURPOSE="${PRODUCER}"
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_write "${RUN}/metadata" "STAGE=stage5\nRUN_ID=${ECODA_RUN_ID}\nSTATE=ACTIVE\nPURPOSE=${PURPOSE}\n"
fi
STATUS="${RUN}/status/analysis_projection.status"
OWNED_FILES=()
failed() {
  local path owner
  ecoda_atomic_write "${STATUS}" "STATE=FAIL\nRUN_ID=${ECODA_RUN_ID}\n" || true
  for path in "${OWNED_FILES[@]:-}" "${OUTPUT}"; do
    [[ -n "${path}" ]] || continue
    owner="$(ecoda_artifact_owner_dir "${path}")" || continue
    if _ecoda_artifact_owner_validate_dir "${owner}" "${path}" &&
       [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]; then
      ecoda_artifact_owner_set_state "${path}" FAIL metadata_projection_failed >/dev/null || true
    fi
  done
  if [[ "${MODE}" == controller ]]; then ecoda_set_run_state FAIL metadata_projection_failed || true; fi
}
if [[ "${MODE}" != validate ]]; then trap failed EXIT; fi
if [[ "${MODE}" == controller ]]; then
  ecoda_atomic_copy "${CONTROL}" "${RUN}/manifests/analysis_release.rds"
  [[ "$(ecoda_md5_file "${RUN}/manifests/analysis_release.rds")" == "${CONTROL_MD5}" ]] || exit 1
  ecoda_atomic_copy "${ECODA_SOURCE_MANIFEST}" "${RUN}/manifests/source.manifest"
  ecoda_runtime_validate_submission apptainer
  for NAME in analysis_release.rds source.manifest runtime.identity operation; do
    ecoda_write_checksum "${RUN}/manifests/${NAME}"
    ecoda_write_artifact_record "${RUN}/manifests/${NAME}" analysis_projection_scope "${ECODA_RUN_ID}" >/dev/null
  done
  # Only control/receipt metadata are read on the controller; RNA-bearing RDS
  # payloads are opened by the single pinned-runtime compute worker below.
  ${PIXI_RSCRIPT} --vanilla - "${PROJECT_ROOT}" "${RUN}" "${CONTROL_MD5}" "${OPERATION}" <<'R'
args <- commandArgs(TRUE)
source(file.path(args[1], "src/utils/analysis_release.R"))
path <- file.path(args[2], "manifests/analysis_release.rds")
.analysis_checked_file(path, args[3])
control <- readRDS(path)
if (args[4] == "render-pca") {
  source(file.path(args[1], "src/utils/datasets_io.R"))
  source(file.path(args[1], "src/notebooks/benchmark_analysis.R"))
  datasets <- read_datasets_json(file.path(args[1], "datasets.json"), view="benchmark_analysis")
  datasets <- datasets[vapply(datasets,function(x)isTRUE(x$use_for_benchmark),logical(1)) & !grepl("_debug$",names(datasets))]
  control <- benchmark_pca_render_scope(control,datasets)
}
inventory <- analysis_release_input_inventory(control)
selection <- control$ledger[c("dataset", "view", "method_key")]
if (!nrow(selection) || anyDuplicated(selection)) stop("Projection selection must be finite and unique")
utils::write.table(selection, file.path(args[2], "manifests/selection.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
if (any(!startsWith(inventory$path, paste0(Sys.getenv("HPC_SCRATCH_DIR"), "/"))))
  stop("Every selected input, receipt and evidence file must be staged explicitly on HPC scratch")
utils::write.table(inventory, file.path(args[2], "manifests/inputs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
.analysis_checked_file(path, args[3])
R
  for NAME in inputs.tsv selection.tsv; do
    ecoda_write_checksum "${RUN}/manifests/${NAME}"
    ecoda_write_artifact_record "${RUN}/manifests/${NAME}" analysis_projection_scope "${ECODA_RUN_ID}" >/dev/null
  done
fi
validate_inputs() {
  local path md5 extra
  ecoda_validate_artifact_record "${RUN}/manifests/inputs.tsv" analysis_projection_scope "${ECODA_RUN_ID}" >/dev/null || return 1
  while IFS=$'\t' read -r path md5 extra; do
    [[ -n "${path}" && "${md5}" =~ ^[0-9a-f]{32}$ && -z "${extra}" ]] || return 1
    ecoda_require_input_ownership "${path}" "${ECODA_RUN_ID}" || return 1
    [[ "$(ecoda_md5_file "${path}")" == "${md5}" ]] || return 1
  done < "${RUN}/manifests/inputs.tsv"
}
validate_scope() {
  local name
  ecoda_runtime_validate_bound_run || return 1
  [[ "$(sed -n 's/^PURPOSE=//p' "${RUN}/metadata")" == "${PURPOSE}" ]] || return 1
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${RUN}/manifests/source.manifest" || return 1
  for name in analysis_release.rds source.manifest runtime.identity inputs.tsv operation selection.tsv; do
    ecoda_validate_artifact_record "${RUN}/manifests/${name}" analysis_projection_scope "${ECODA_RUN_ID}" >/dev/null || return 1
  done
  ${PIXI_RSCRIPT} --vanilla - "${PROJECT_ROOT}" "${RUN}" "${OPERATION}" <<'R' || return 1
a <- commandArgs(TRUE)
source(file.path(a[1], "src/utils/analysis_release.R"))
control <- readRDS(file.path(a[2], "manifests/analysis_release.rds"))
if (a[3] == "render-pca") {
  source(file.path(a[1], "src/utils/datasets_io.R"))
  source(file.path(a[1], "src/notebooks/benchmark_analysis.R"))
  datasets <- read_datasets_json(file.path(a[1], "datasets.json"), view="benchmark_analysis")
  datasets <- datasets[vapply(datasets,function(x)isTRUE(x$use_for_benchmark),logical(1)) & !grepl("_debug$",names(datasets))]
  control <- benchmark_pca_render_scope(control,datasets)
}
selection <- control$ledger[c("dataset", "view", "method_key")]
if (!nrow(selection) || anyDuplicated(selection)) stop("Invalid projection selection")
expected <- paste0(paste(apply(selection, 1L, paste, collapse="\t"), collapse="\n"), "\n")
path <- file.path(a[2], "manifests/selection.tsv")
if (!identical(rawToChar(readBin(path, "raw", n=file.info(path)$size)), expected))
  stop("Projection selection differs from the operation-scoped control")
R
  validate_inputs
}
validate_report() {
  local path hash extra owner
  validate_scope || return 1
  [[ -s "${STATUS}" && "$(sed -n 's/^STATE=//p' "${STATUS}")" == OK &&
     "$(sed -n 's/^RUN_ID=//p' "${STATUS}")" == "${ECODA_RUN_ID}" ]] || return 1
  ecoda_validate_artifact_record "${OUTPUT}/files.tsv" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
  for path in "${OUTPUT}" "${OUTPUT}/files.tsv"; do
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" &&
       "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done
  while IFS=$'\t' read -r path hash extra; do
    [[ "${path}" == "${OUTPUT}/"* && "${hash}" =~ ^[0-9a-f]{32}$ && -z "${extra}" ]] || return 1
    ecoda_validate_artifact_record "${path}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null || return 1
    [[ "${ECODA_ARTIFACT_RECORD_MD5}" == "${hash}" ]] || return 1
    owner="$(ecoda_artifact_owner_dir "${path}")" || return 1
    _ecoda_artifact_owner_validate_dir "${owner}" "${path}" || return 1
    [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" &&
       "${ECODA_ARTIFACT_OWNER_STAGE}" == stage5 && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]] || return 1
  done < "${OUTPUT}/files.tsv"
  if [[ "${OPERATION}" == render-pca ]]; then
    ${PIXI_RSCRIPT} --vanilla - "${PROJECT_ROOT}" "${OUTPUT}/render_manifest.rds" "$(ecoda_md5_file "${OUTPUT}/render_manifest.rds")" "$(ecoda_md5_file "${RUN}/manifests/analysis_release.rds")" <<'R'
a <- commandArgs(TRUE)
source(file.path(a[1],"src/utils/analysis_release.R"))
source(file.path(a[1],"src/notebooks/benchmark_analysis.R"))
report <- read_benchmark_pca_render(a[2],a[3],character())
if(!identical(report$control_md5,a[4]))stop("Rendered PCA control differs")
R
  fi
}
if [[ "${MODE}" == validate ]]; then
  validate_report
  [[ "$(sed -n 's/^STATE=//p' "${RUN}/status/terminal")" == OK &&
     "$(sed -n 's/^RUN_ID=//p' "${RUN}/status/terminal")" == "${ECODA_RUN_ID}" ]] || exit 1
  for NAME in scheduler_ids.tsv sbatch_submission.txt; do
    ecoda_validate_artifact_record "${RUN}/manifests/${NAME}" analysis_projection_submission "${ECODA_RUN_ID}" >/dev/null
  done
  trap - EXIT
  echo "ANALYSIS_PROJECTION_VALIDATED=${ECODA_RUN_ID}"
  exit 0
fi
if [[ "${MODE}" == worker ]]; then
  ecoda_runtime_reexec_worker stage5 "${SELF}" --worker
  cd "${PROJECT_ROOT}"
  validate_scope
  owner="$(ecoda_artifact_owner_dir "${OUTPUT}")"
  _ecoda_artifact_owner_validate_dir "${owner}" "${OUTPUT}"
  [[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE &&
     ! -e "${OUTPUT}" && ! -L "${OUTPUT}" ]]
  CONTROL_MD5="$(ecoda_md5_file "${RUN}/manifests/analysis_release.rds")"
${PIXI_RSCRIPT} --vanilla - "${PROJECT_ROOT}" "${RUN}" "${CONTROL_MD5}" "${OPERATION}" <<'R'
args <- commandArgs(TRUE)
source(file.path(args[1], "src/utils/analysis_release.R"))
source(file.path(args[1], "src/utils/evaluation_scores.R"))
source(file.path(args[1], "src/5_run_benchmark_methods/benchmark_pipeline.R"))
path <- file.path(args[2], "manifests/analysis_release.rds")
.analysis_checked_file(path, args[3])
control <- readRDS(path)
if (args[4] == "projection") {
  project_analysis_release(control, file.path(args[2], "artifacts/analysis_projection"),
    c(source = file.path(args[2], "manifests/source.manifest"), runtime = file.path(args[2], "manifests/runtime.identity"), control = path))
} else {
  suppressPackageStartupMessages({library(dplyr);library(ggplot2)})
  for(file in c("src/utils/datasets_io.R","src/utils/math_utils.R","src/utils/scoring_metrics.R",
      "src/utils/plotting.R","src/notebooks/benchmark_analysis.R"))source(file.path(args[1],file))
  datasets <- read_datasets_json(file.path(args[1],"datasets.json"),view="benchmark_analysis")
  datasets <- datasets[vapply(datasets,function(x)isTRUE(x$use_for_benchmark),logical(1)) & !grepl("_debug$",names(datasets))]
  render_benchmark_pca_release(control,file.path(args[2],"artifacts/benchmark_pca"),datasets,args[3])
}
.analysis_checked_file(path, args[3])
R
ecoda_runtime_validate_bound_run
validate_inputs
# Every path is from the projector's finite constructed output inventory.
while IFS=$'\t' read -r PATH_IN MD5_IN EXTRA; do
  [[ "${PATH_IN}" == "${OUTPUT}/"* && "${MD5_IN}" =~ ^[0-9a-f]{32}$ && -z "${EXTRA}" &&
     "$(ecoda_realpath_existing "${PATH_IN}")" == "${PATH_IN}" &&
     "$(ecoda_md5_file "${PATH_IN}")" == "${MD5_IN}" ]] || exit 1
  ecoda_artifact_owner_acquire "${PATH_IN}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
  OWNED_FILES+=("${PATH_IN}")
  ecoda_write_checksum "${PATH_IN}"
  ecoda_write_artifact_record "${PATH_IN}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
done < "${OUTPUT}/files.tsv"
ecoda_artifact_owner_acquire "${OUTPUT}/files.tsv" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
OWNED_FILES+=("${OUTPUT}/files.tsv")
ecoda_write_checksum "${OUTPUT}/files.tsv"
ecoda_write_artifact_record "${OUTPUT}/files.tsv" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
for PATH_IN in "${OWNED_FILES[@]}"; do
  ecoda_validate_artifact_record "${PATH_IN}" "${PRODUCER}" "${ECODA_RUN_ID}" >/dev/null
  ecoda_artifact_owner_set_state "${PATH_IN}" OK metadata_projection_complete >/dev/null
done
ecoda_artifact_owner_set_state "${OUTPUT}" OK metadata_projection_complete >/dev/null
  ecoda_atomic_write "${STATUS}" "STATE=OK\nRUN_ID=${ECODA_RUN_ID}\n"
  validate_report
  trap - EXIT
  exit 0
fi
validate_scope
mkdir -p "${OUTPUT%/*}"
ecoda_artifact_owner_acquire "${OUTPUT}" stage5 "${ECODA_RUN_ID}" 0 0 0 >/dev/null
RUNTIME_EXPORT="$(ecoda_runtime_export_csv stage5 0)"
set +e
JOB="$(sbatch --parsable --wait --partition="${SLURM_PARTITION}" --ntasks=1 --cpus-per-task=1 \
  --mem=64G --time=02:00:00 --output="${RUN}/logs/analysis_projection_%j.log" \
  --error="${RUN}/logs/analysis_projection_%j.err" --export="ALL,${RUNTIME_EXPORT},ECODA_RUN_ROOT=${RUN}" \
  "${SELF}" --worker)"
RC=$?
set -e
ecoda_atomic_write "${RUN}/manifests/sbatch_submission.txt" "${JOB}\n"
JOB="${JOB%%;*}"
[[ "${JOB}" =~ ^[0-9]+$ ]] || { echo 'ERROR: ambiguous projection submission; preserve run and inspect' >&2; exit 1; }
echo "ANALYSIS_PROJECTION_JOB_ID=${JOB}"
ecoda_atomic_write "${RUN}/manifests/scheduler_ids.tsv" "ARRAY\t${JOB}\n"
for NAME in scheduler_ids.tsv sbatch_submission.txt; do
  ecoda_write_checksum "${RUN}/manifests/${NAME}"
  ecoda_write_artifact_record "${RUN}/manifests/${NAME}" analysis_projection_submission "${ECODA_RUN_ID}" >/dev/null
done
[[ ${RC} -eq 0 ]] || exit "${RC}"
validate_report
ecoda_set_run_state OK metadata_projection_complete_scientific_origins_preserved
trap - EXIT
if [[ "${OPERATION}" == render-pca ]]; then
  printf 'BENCHMARK_PCA_RENDER_RDS=%s/render_manifest.rds\nNO_MODEL_FITS=1\n' "${OUTPUT}"
else
  printf 'ANALYSIS_PROJECTION_CONTROL=%s/analysis_release.rds\nNO_MODEL_FITS=1\n' "${OUTPUT}"
fi
