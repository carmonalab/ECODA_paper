#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage2-watchdog.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status" "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/manifests" "${TMP_DIR}/home/scratch/ECODA_paper/Kfoury/data"
printf 'STAGE=stage2\nRUN_ID=run\nSTATE=ACTIVE\n' > "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/metadata"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sacct" <<'STUB'
#!/bin/bash
case "$*" in
  *1001*) printf 'OUT_OF_MEMORY\n' ;;
  *1002*) printf 'COMPLETED\n' ;;
  *) printf 'COMPLETED\n' ;;
esac
STUB
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
printf '%s\n' "$*" >> "${CAPTURE}"
printf '1002\n'
STUB
chmod +x "${TMP_DIR}/bin/sacct" "${TMP_DIR}/bin/sbatch"
printf 'rds payload\n' > "${TMP_DIR}/home/scratch/ECODA_paper/Kfoury/data/Kfoury_2021_34719426.rds"
OWNER_DIR="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners/stage2/kfoury_lowres_ct"
mkdir -p "${OWNER_DIR}"
printf 'RUN_ID=run\nSTATE=ACTIVE\nSTAGE=stage2\nKEY=kfoury_lowres_ct\n' > "${OWNER_DIR}/owner"
MANIFEST="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/manifests/steps.tsv"
JOB_FILE="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/manifests/jobs.tsv"
printf 'kfoury_lowres_ct\t%s\t%s\t-\t%s\n' "${ROOT}/src/2_dataset_specific_preprocessing/1.4_submit_kfoury_lowres_ct.sh" "${TMP_DIR}/home/scratch/ECODA_paper/Kfoury/data/Kfoury_2021_34719426.rds" "${OWNER_DIR}" > "${MANIFEST}"
printf 'kfoury_lowres_ct\t1001\n' > "${JOB_FILE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL=test@example.invalid STAGE2_FORCE=1 STAGE2_WATCHDOG_MAX_POLLS=1 \
  bash "${ROOT}/src/2_dataset_specific_preprocessing/stage2_watchdog.sh" run "${MANIFEST}" "${JOB_FILE}" 128G 256G shared-cpu 1000
[[ "$(grep '^STATE=' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == "STATE=OK" ]]
[[ -s "${TMP_DIR}/home/scratch/ECODA_paper/Kfoury/data/Kfoury_2021_34719426.rds.md5" ]]
[[ "$(grep '^STATE=' "${OWNER_DIR}/owner")" == "STATE=OK" ]]
[[ "$(grep -c '^SCHEDULER_ID=' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 2 ]]
[[ "$(grep -c '^SCHEDULER_ID=1001$' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 1 ]]
[[ "$(grep -c '^SCHEDULER_ID=1002$' "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/run/status/watchdog")" == 1 ]]
case "$(cat "${CAPTURE}")" in *"FORCE_PREPROCESS=1"*) ;; *) echo "OOM retry dropped FORCE_PREPROCESS=1" >&2; exit 1 ;; esac
echo "stage2 watchdog: OK"
