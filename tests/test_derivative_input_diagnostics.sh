#!/bin/bash
# No runtime, scheduler or production data: preserve the actual wrapper API.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-resolver-diagnostic.XXXXXX")"
TMP="$(cd "${TMP}" && pwd -P)"
trap 'rm -rf "${TMP}"' EXIT
source "${ROOT}/src/utils/bash/ecoda_derivative_artifacts.sh"
export ECODA_RUN_ROOT="${TMP}/run" ECODA_RUN_ID=selected_run
export ECODA_ARTIFACT_BINDING="${ECODA_RUN_ROOT}/manifests/artifact_bindings.json"
export NAS_TARGET_DIR="${TMP}/nas" SLURM_JOB_ID=12345 SLURM_ARRAY_TASK_ID=2 SLURMD_NODENAME=node_fixture
mkdir -p "${ECODA_RUN_ROOT}/manifests"
printf 'BINDING_CONTENT_MUST_NOT_BE_PRINTED\n' > "${ECODA_ARTIFACT_BINDING}"
printf 'UNRELATED_SECRET_CONTENT\n' > "${ECODA_RUN_ROOT}/manifests/unrelated.txt"
PYTHON="${ROOT}/.pixi/envs/default/bin/python"
MODE=success
_ecoda_derivative_cli() {
  printf '%s\0' "$@" > "${TMP}/args"
  if [[ "${MODE}" == success ]]; then
    printf '/exact/resolved/path\n\n'
    return 0
  fi
  "${PYTHON}" - "${ECODA_RUN_ROOT}/manifests/input_artifacts.tsv" <<'PY'
import sys
try:
    open(sys.argv[1], 'rb')
except FileNotFoundError as error:
    print(f'FileNotFoundError: {error}',file=sys.stderr)
    raise SystemExit(37)
raise SystemExit('missing fixture must remain missing')
PY
}
ecoda_derivative_input Toy benchmark_analysis directory content > "${TMP}/success.out" 2> "${TMP}/success.err"
printf '/exact/resolved/path\n\n' > "${TMP}/expected.out"
cmp "${TMP}/success.out" "${TMP}/expected.out"
[[ ! -s "${TMP}/success.err" ]]
printf '%s\0' input --binding "${ECODA_ARTIFACT_BINDING}" --run-root "${ECODA_RUN_ROOT}" \
  --nas-root "${NAS_TARGET_DIR}" --dataset Toy --view benchmark_analysis --directory --verify-content > "${TMP}/expected.args"
cmp "${TMP}/args" "${TMP}/expected.args"
MODE=missing
rc=0
ecoda_derivative_input Toy benchmark_analysis directory content > "${TMP}/failure.out" 2> "${TMP}/failure.err" || rc=$?
[[ ${rc} -eq 37 && ! -s "${TMP}/failure.out" ]]
cmp "${TMP}/args" "${TMP}/expected.args"
grep -q 'FileNotFoundError:' "${TMP}/failure.err"
grep -q 'rc=37 dataset=Toy view=benchmark_analysis run_id=selected_run job_id=12345 task_id=2 node=node_fixture' "${TMP}/failure.err"
for path in "${ECODA_ARTIFACT_BINDING}" "${ECODA_RUN_ROOT}/manifests/input_artifacts.tsv" "${ECODA_RUN_ROOT}/manifests"; do
  grep -F -q "path=${path}" "${TMP}/failure.err"
done
if grep -E -q 'BINDING_CONTENT_MUST_NOT_BE_PRINTED|UNRELATED_SECRET_CONTENT|unrelated.txt' "${TMP}/failure.err"; then
  echo 'diagnostics exposed unrelated content or inventory' >&2; exit 1
fi
# Even unavailable native stat diagnostics must preserve the resolver status.
stat() { return 88; }
rc=0
ecoda_derivative_input Toy benchmark_analysis directory content > "${TMP}/stat.out" 2> "${TMP}/stat.err" || rc=$?
[[ ${rc} -eq 37 && ! -s "${TMP}/stat.out" ]]
grep -q 'DERIVATIVE_RESOLVER_STAT unavailable' "${TMP}/stat.err"
echo 'PASS: resolver success bytes/arguments unchanged; failure-only context and exact-path stat, original exit preserved, no contents or retries'
