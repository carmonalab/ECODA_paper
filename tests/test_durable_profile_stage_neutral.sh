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
PROFILE="${ROOT}/.agents/skills/durable-hpc-gate-ecoda/references/profile.json"
GLOBAL_GATE="${HOME}/.agents/skills/durable-hpc-gate/scripts/durable_hpc_gate.py"

# Snapshot-backed runs must not bind terminal identity to mutable canonical
# checkout fingerprints or broad run-root/status discovery.
[[ "$(jq -r '.policy.immutable_fingerprints | length' "${PROFILE}")" == "0" ]]
[[ "$(jq -r '[.policy.immutable_fingerprints[]?.name] | any(.[]; . == "datasets-pixi-lock-sha256" or . == "repository-head")' "${PROFILE}")" == "false" ]]
for broad_contract in \
  terminal-run-owned-manifests-status \
  terminal-emitted-scheduler-status-records \
  stage-neutral-run-root-contract \
  stage-neutral-scheduler-status-contract; do
  [[ "$(jq -r --arg name "${broad_contract}" '[.policy.audit_commands[]?.name, .policy.artifact_contracts[]?.name] | any(.[]; . == $name)' "${PROFILE}")" == "false" ]]
done
[[ "$(jq -r '[.policy.audit_commands[]?.command, .policy.artifact_contracts[]?.command] | any(.[]; contains("_ecoda_runs/*"))' "${PROFILE}")" == "false" ]]

# The retained profile contracts remain stage-neutral: canonical roots and
# wrapper presence are checked, while the durable gate owns exact-command
# identity, one accounting query, and reviewer release.
jq -e 'any(.policy.invariants[]; .name == "bamboo-repository-root")' "${PROFILE}" >/dev/null
jq -e 'any(.policy.audit_commands[]; .name == "terminal-repository-root")' "${PROFILE}" >/dev/null
jq -e 'any(.policy.audit_commands[]; .name == "terminal-canonical-artifact-roots")' "${PROFILE}" >/dev/null
jq -e 'any(.policy.audit_commands[]; .name == "terminal-canonical-wrapper-set")' "${PROFILE}" >/dev/null
jq -e 'any(.policy.artifact_contracts[]; .name == "scratch-and-nas-outputs-use-configured-roots")' "${PROFILE}" >/dev/null
jq -e '.policy.accounting_command == "sacct -n -P -X -j {scheduler_ids} --format=JobIDRaw,State,ExitCode" and .policy.require_scheduler_ids == true and .policy.reviewer_required == true' "${PROFILE}" >/dev/null

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

# Exercise the real synthetic profile parser with an explicit empty list.
SYNTHETIC_PROFILE="${TMP_DIR}/synthetic-empty-fingerprints.json"
cat > "${SYNTHETIC_PROFILE}" <<'JSON'
{
  "name": "synthetic-empty-fingerprints",
  "execution": {
    "mode": "synthetic",
    "runner_shell": "/bin/bash",
    "ssh_binary": null,
    "tmux_binary": null
  },
  "policy": {
    "manifest_constraints": {},
    "invariants": [],
    "audit_commands": [],
    "artifact_contracts": [],
    "immutable_fingerprints": [],
    "accounting_command": null,
    "require_scheduler_ids": false,
    "reviewer_required": true,
    "known_discrepancies": []
  }
}
JSON
EMPTY_MANIFEST="${TMP_DIR}/empty-fingerprint-manifest.json"
EMPTY_OUTPUT="${TMP_DIR}/empty-fingerprint-prepare.json"
EMPTY_WORKDIR="${TMP_DIR}/empty-fingerprint-workdir"
mkdir -p "${EMPTY_WORKDIR}"
python3 -B "${GLOBAL_GATE}" prepare \
  --manifest "${EMPTY_MANIFEST}" \
  --profile "${SYNTHETIC_PROFILE}" \
  --project synthetic \
  --gate-id empty-fingerprints \
  --remote-host synthetic-host \
  --remote-workdir "${EMPTY_WORKDIR}" \
  --exact-command true \
  --serialization-group empty-fingerprint-group \
  --tmux-session empty-fingerprint-tmux \
  --completion-channel "file:${TMP_DIR}/empty-fingerprint.done" \
  --remote-manifest "${TMP_DIR}/empty-fingerprint.remote.json" \
  --remote-runner "${TMP_DIR}/empty-fingerprint.runner.sh" \
  --remote-log "${TMP_DIR}/empty-fingerprint.log" \
  --remote-status "${TMP_DIR}/empty-fingerprint.status.json" \
  --output "${EMPTY_OUTPUT}" >/dev/null
jq -e '.ok == true and .state == "PREPARED"' "${EMPTY_OUTPUT}" >/dev/null
jq -e '.immutable_fingerprints == [] and .immutable_fingerprint_commands == [] and .exact_command == "true" and (.command_digest | test("^[0-9a-f]{64}$"))' "${EMPTY_MANIFEST}" >/dev/null

# The exact run-bound validator is scoped by explicit paths and cannot be
# satisfied by a separate run's evidence.
AUDIT_SCRIPT="${ROOT}/src/utils/bash/ecoda_run_audit.sh"
AUDIT_HELP="$(bash "${AUDIT_SCRIPT}" --help)"
[[ "${AUDIT_HELP}" == *"--run-root ABSOLUTE_RUN_ROOT"* ]]
[[ "${AUDIT_HELP}" == *"--stage STAGE"* ]]
[[ "${AUDIT_HELP}" == *"--selection ABSOLUTE_SELECTION"* ]]
[[ "${AUDIT_HELP}" == *"--source-manifest ABSOLUTE_SOURCE_MANIFEST"* ]]
[[ "${AUDIT_HELP}" == *"--runtime-identity ABSOLUTE_RUNTIME_IDENTITY"* ]]
[[ "${AUDIT_HELP}" != *"_ecoda_runs/*"* ]]
BAD_ROOT="${FAKE_SCRATCH}/_ecoda_runs/bad_run"
DECOY_ROOT="${FAKE_SCRATCH}/_ecoda_runs/decoy_run"
mkdir -p "${BAD_ROOT}" "${DECOY_ROOT}/status"
printf 'SCHEDULER_ID=decoy\n' > "${DECOY_ROOT}/status/watchdog"
set +e
bash "${AUDIT_SCRIPT}" \
  --run-root "${BAD_ROOT}" \
  --stage stage2 \
  --selection "${BAD_ROOT}/manifests/selection.tsv" \
  --source-manifest "${BAD_ROOT}/manifests/source.manifest" \
  --runtime-identity "${BAD_ROOT}/manifests/runtime.identity" \
  >/dev/null 2>&1
AUDIT_RC=$?
set -e
[[ ${AUDIT_RC} -ne 0 ]]

echo "durable profile stage-neutral audit: OK"
