#!/bin/bash
# Real run/record/owner helpers; runtime binding and R audit are bounded stubs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-bl-wrapper.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'chmod -R u+w "${TMP}"; rm -rf "${TMP}"' EXIT
TREE="${TMP}/snapshot/tree"; SCRATCH="${TMP}/scratch"
mkdir -p "${TREE}/src/utils/bash" "${TREE}/src/utils/r" "${SCRATCH}/_ecoda_inputs" "${SCRATCH}/_ecoda_runs" "${TMP}/bin"
cp "${ROOT}/src/utils/bash/validate_bl_science.sh" "${ROOT}/src/utils/bash/ecoda_run_common.sh" "${TREE}/src/utils/bash/"
printf '# real scientific audit is covered independently\n' > "${TREE}/src/utils/r/validate_bl_science.R"
printf 'SOURCE_ROOT=%s\n' "${TREE}" > "${TREE}/source.manifest"
cat > "${TREE}/src/slurm_config.sh" <<'SH'
export PROJECT_ROOT="${TEST_TREE}" HPC_SCRATCH_DIR="${TEST_SCRATCH}"
export PYTHON_BIN="${TEST_PYTHON}" PIXI_RSCRIPT="${TEST_RSCRIPT}"
SH
cat > "${TREE}/src/utils/bash/ecoda_runtime.sh" <<'SH'
ecoda_runtime_validate_submission() {
  [[ "$1" == apptainer && "${ECODA_RUNTIME_PROFILE}" == stage5 ]]
  printf 'new validation runtime\n' > "${ECODA_RUN_ROOT}/manifests/runtime.identity"
  printf 'RUNTIME_BOUND\n' >> "${TEST_FLOW}"
}
ecoda_runtime_validate_bound_run() {
  [[ -s "${ECODA_RUN_ROOT}/manifests/runtime.identity" ]]
  cmp -s "${ECODA_SOURCE_MANIFEST}" "${ECODA_RUN_ROOT}/manifests/source.manifest"
}
SH
cat > "${TMP}/bin/Rscript" <<'SH'
#!/bin/bash
set -euo pipefail
[[ "$1" == --vanilla && "$2" == "${TEST_TREE}/src/utils/r/validate_bl_science.R" ]]
shift 2
source "${TEST_TREE}/src/utils/bash/ecoda_run_common.sh"
[[ "$1" == "${ECODA_RUN_ROOT}/manifests/control.json" && "$4" == "${TEST_TREE}" ]]
[[ "$(ecoda_md5_file "$1")" == "$2" ]]
for name in control.json source.manifest runtime.identity; do
  ecoda_validate_artifact_record "${ECODA_RUN_ROOT}/manifests/${name}" bl_science_scope "${ECODA_RUN_ID}" >/dev/null
done
owner="$(ecoda_artifact_owner_dir "$3")"
_ecoda_artifact_owner_validate_dir "${owner}" "$3"
[[ "${ECODA_ARTIFACT_OWNER_RUN}" == "${ECODA_RUN_ID}" && "${ECODA_ARTIFACT_OWNER_STATE}" == ACTIVE ]]
printf 'OWNER_ACTIVE_AND_SCOPE_SEALED\n' >> "${TEST_FLOW}"
if [[ "${TEST_BL_FAILURE:-0}" == 1 ]]; then
  printf 'preserved partial evidence\n' > "$3"
  exit 7
fi
printf '{"status":"NOOP_VALIDATED","publication_completed":false}\n' > "$3"
SH
chmod +x "${TMP}/bin/Rscript"
for command in sbatch sacct squeue; do
  printf '#!/bin/bash\nprintf forbidden >> "%s"\nexit 99\n' "${TMP}/scheduler.called" > "${TMP}/bin/${command}"
  chmod +x "${TMP}/bin/${command}"
done
export TEST_TREE="${TREE}" TEST_SCRATCH="${SCRATCH}" TEST_RSCRIPT="${TMP}/bin/Rscript"
export TEST_PYTHON="${ROOT}/.pixi/envs/default/bin/python" TEST_FLOW="${TMP}/flow"
export ECODA_SOURCE_ROOT="${TREE}" ECODA_SOURCE_MANIFEST="${TREE}/source.manifest" ECODA_SOURCE_SNAPSHOT_REQUIRED=1
export PATH="${TMP}/bin:${PATH}"
CONTROL="${SCRATCH}/_ecoda_inputs/control.json"
printf '{"schema":"wrapper-only-fixture"}\n' > "${CONTROL}"
source "${TREE}/src/slurm_config.sh"
source "${TREE}/src/utils/bash/ecoda_run_common.sh"
MD5="$(ecoda_md5_file "${CONTROL}")"
entry="${TREE}/src/utils/bash/validate_bl_science.sh"
invoke() {
  local id="$1" report="${SCRATCH}/_ecoda_runs/$1/artifacts/evidence/bl_science.json"
  ECODA_RUN_ID="${id}" ECODA_RUN_ROOT="${SCRATCH}/_ecoda_runs/${id}" \
    bash "${entry}" --control "${CONTROL}" --control-md5 "${MD5}" --report "${report}"
}
expect_fail() { if "$@" > "${TMP}/failure.log" 2>&1; then echo 'expected wrapper failure' >&2; exit 1; fi; }
invoke good > "${TMP}/success.log"
run="${SCRATCH}/_ecoda_runs/good"; report="${run}/artifacts/evidence/bl_science.json"
grep -q '^STATE=NOOP_VALIDATED$' "${run}/status/terminal"
ecoda_validate_artifact_record "${report}" bl_science_validation good >/dev/null
_ecoda_artifact_owner_validate_dir "$(ecoda_artifact_owner_dir "${report}")" "${report}"
[[ "${ECODA_ARTIFACT_OWNER_RUN}" == good && "${ECODA_ARTIFACT_OWNER_STATE}" == OK ]]
cmp "${CONTROL}" "${run}/manifests/control.json"
[[ "$(cat "${TEST_FLOW}")" == $'RUNTIME_BOUND\nOWNER_ACTIVE_AND_SCOPE_SEALED' ]]
before="$(ecoda_md5_file "${report}")"; expect_fail invoke good
[[ "$(ecoda_md5_file "${report}")" == "${before}" ]]
MD5=00000000000000000000000000000000; expect_fail invoke invalid_control
[[ ! -e "${SCRATCH}/_ecoda_runs/invalid_control" ]]
MD5="$(ecoda_md5_file "${CONTROL}")"
ECODA_RUN_ID=wrong_report expect_fail bash "${entry}" --control "${CONTROL}" --control-md5 "${MD5}" \
  --report "${SCRATCH}/_ecoda_runs/another/artifacts/evidence/bl_science.json"
[[ ! -e "${SCRATCH}/_ecoda_runs/wrong_report" ]]
TEST_BL_FAILURE=1 expect_fail invoke failed
run="${SCRATCH}/_ecoda_runs/failed"; report="${run}/artifacts/evidence/bl_science.json"
grep -q '^STATE=FAIL$' "${run}/status/terminal"
grep -q 'preserved partial evidence' "${report}"
_ecoda_artifact_owner_validate_dir "$(ecoda_artifact_owner_dir "${report}")" "${report}"
[[ "${ECODA_ARTIFACT_OWNER_RUN}" == failed && "${ECODA_ARTIFACT_OWNER_STATE}" == FAIL ]]
[[ ! -e "${TMP}/scheduler.called" ]]
echo 'PASS: fresh executor-bound validator run, sealed source/runtime/control, owned report, NOOP terminal and fail-closed partial preservation'
