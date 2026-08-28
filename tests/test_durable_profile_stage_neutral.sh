#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-profile-audit.XXXXXX")"
TMP_DIR="$(cd "${TMP_DIR}" && pwd)"
trap 'rm -rf "${TMP_DIR}"' EXIT
TMP_DIR="$(cd "${TMP_DIR}" && pwd -P)"
FAKE_HOME="${TMP_DIR}/home"
FAKE_REPO="${FAKE_HOME}/ECODA_paper"
FAKE_SCRATCH="${FAKE_HOME}/scratch/ECODA_paper"
FAKE_NAS="${FAKE_HOME}/nas"
mkdir -p "${FAKE_REPO}/src/2_dataset_specific_preprocessing" \
  "${FAKE_REPO}/src/3_scrnaseq_preprocessing" \
  "${FAKE_REPO}/src/4_cell_type_annotation" \
  "${FAKE_REPO}/src/5_run_benchmark_methods" \
  "${FAKE_REPO}/src/utils/bash" "${FAKE_SCRATCH}/_ecoda_runs" "${FAKE_NAS}"
for path in \
  src/2_dataset_specific_preprocessing/1_submit_hpc.sh \
  src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh \
  src/4_cell_type_annotation/1_submit_onboarding_stage.sh \
  src/5_run_benchmark_methods/1_submit_hpc_array.sh \
  src/utils/bash/ecoda_run_common.sh; do
  touch "${FAKE_REPO}/${path}"
done
cat > "${FAKE_REPO}/src/slurm_config.sh" <<'CONFIG'
#!/bin/bash
export PROJECT_ROOT="$HOME/ECODA_paper"
export HPC_SCRATCH_DIR="$HOME/scratch/ECODA_paper"
export NAS_TARGET_DIR="$HOME/nas"
CONFIG
printf '{"Bassez":{"columns":{"sample":"Sample","label":"Status"}}}\n' > "${FAKE_REPO}/datasets.json"
touch "${FAKE_REPO}/pixi.toml" "${FAKE_REPO}/pixi.lock" "${FAKE_REPO}/AGENTS.md"
for stage in stage2 stage3 stage4 stage5; do
  run="${FAKE_SCRATCH}/_ecoda_runs/${stage}_run"
  mkdir -p "${run}/manifests" "${run}/status"
  printf 'SCHEDULER_ID=%s\n' "${stage}" > "${run}/status/watchdog"
done
PROFILE="${ROOT}/.agents/skills/durable-hpc-gate-ecoda/references/profile.json"
[[ "$(jq -r '[.policy.audit_commands[].command, .policy.artifact_contracts[].command] | any(contains("ECODA_GATE_STAGE"))' "${PROFILE}")" == false ]]
run_entries() {
  local section="$1" name command
  while IFS=$'\t' read -r name command; do
    HOME="${FAKE_HOME}" bash -c "${command}" || {
      echo "profile command failed: ${name}" >&2
      return 1
    }
  done < <(jq -r --arg section "${section}" '.policy[$section][] | [.name,.command] | @tsv' "${PROFILE}")
}
run_entries audit_commands
run_entries artifact_contracts
BAD_ROOT="${FAKE_SCRATCH}/_ecoda_runs/bad_run"
mkdir -p "${BAD_ROOT}/status"
AUDIT_CMD="$(jq -r '.policy.audit_commands[] | select(.name == "terminal-run-owned-manifests-status") | .command' "${PROFILE}")"
set +e
HOME="${FAKE_HOME}" bash -c "${AUDIT_CMD}"
AUDIT_RC=$?
set -e
[[ ${AUDIT_RC} -ne 0 ]]
echo "durable profile stage-neutral audit: OK"
