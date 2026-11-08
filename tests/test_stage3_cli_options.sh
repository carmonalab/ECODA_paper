#!/bin/bash
# Exercise the real argument parser without runtime, cohort, or Slurm effects.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage3-cli.XXXXXX")"
trap 'rm -rf "${FIXTURE}"' EXIT
mkdir -p "${FIXTURE}/src/3_scrnaseq_preprocessing" "${FIXTURE}/src/utils/bash" "${FIXTURE}/bin"
cp "${ROOT}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh" \
  "${FIXTURE}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh"
cat > "${FIXTURE}/src/slurm_config.sh" <<'CONFIG'
PROJECT_ROOT="${SCRIPT_DIR}/../.."
SLURM_PARTITION=fixture
MAX_NUM_CHUNKS_PARALLEL=1
CONFIG
for helper in ecoda_runtime.sh ecoda_run_common.sh h5ad_preflight_submit.sh; do
  : > "${FIXTURE}/src/utils/bash/${helper}"
done
cat > "${FIXTURE}/bin/sbatch" <<'SBATCH'
#!/bin/bash
touch "${SBATCH_CALLED}"
exit 99
SBATCH
chmod +x "${FIXTURE}/bin/sbatch"
export PATH="${FIXTURE}/bin:${PATH}" SBATCH_CALLED="${FIXTURE}/sbatch.called"
SUBMITTER="${FIXTURE}/src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh"
check_help() {
  bash "${SUBMITTER}" "$@" --help > "${FIXTURE}/help.txt" 2>&1 || {
    cat "${FIXTURE}/help.txt" >&2
    return 1
  }
  [[ "$(head -n 1 "${FIXTURE}/help.txt")" == 'Usage: 1_submit_hpc_array.sh '* ]]
  [[ ! -e "${SBATCH_CALLED}" ]]
}
check_help --throttle 17
check_help --throttle=17
# These are the complete Stage 3 launch flags used for this derivative wave.
check_help --selection-file /unused/selection.tsv \
  --derivative-release fixture-release --input-artifacts /unused/input_artifacts.tsv \
  --throttle 17 --mem 128G --max-mem 500G
check_help --selection-file=/unused/selection.tsv \
  --derivative-release fixture-release --input-artifacts /unused/input_artifacts.tsv \
  --throttle=17 --mem=128G --max-mem=500G
echo 'PASS: Stage3 spaced/equal throttle and complete derivative launch flags parse without Slurm'
