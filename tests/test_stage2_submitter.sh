#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/ecoda-stage2.XXXXXX")"
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/bin" "${TMP_DIR}/home"
CAPTURE="${TMP_DIR}/sbatch.calls"
export CAPTURE
cat > "${TMP_DIR}/bin/sbatch" <<'STUB'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${CAPTURE}"
N="$(wc -l < "${CAPTURE}" | tr -d '[:space:]')"
printf '71000%s\n' "${N}"
STUB
chmod +x "${TMP_DIR}/bin/sbatch"
OUTPUT="$(
  HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  STAGE2_SUBMITTER_TEST=1 bash "${ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
    --datasets CombinedPBMC,_debug --steps combinedpbmc,joanito --force
)"
RUN_ID="$(printf '%s\n' "${OUTPUT}" | sed -n 's/^STAGE2_RUN_ID=//p')"
[[ -n "${RUN_ID}" ]]
MANIFEST="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}/manifests/steps.tsv"
[[ "$(wc -l < "${MANIFEST}" | tr -d '[:space:]')" == 3 ]]
SCHEDULER_MANIFEST="${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_runs/${RUN_ID}/manifests/scheduler_ids.tsv"
[[ "$(wc -l < "${SCHEDULER_MANIFEST}" | tr -d '[:space:]')" == 4 ]]
CALLS="$(cat "${CAPTURE}")"
case "${CALLS}" in *"FORCE_PREPROCESS=1"*) ;; *) echo "force export missing" >&2; exit 1 ;; esac
case "${CALLS}" in *"--dependency=afterok:710001"*) ;; *) echo "CombinedPBMC cap dependency missing" >&2; exit 1 ;; esac
JOANITO_CALL="$(sed -n '3p' "${CAPTURE}")"
case "${JOANITO_CALL}" in *"--dependency="*) echo "Joanito was artificially serialized" >&2; exit 1 ;; esac
WATCHDOG_CALL="$(sed -n '4p' "${CAPTURE}")"
case "${WATCHDOG_CALL}" in *"--dependency=afterany:710001:710002:710003"*) ;; *) echo "aggregate watchdog dependency missing" >&2; exit 1 ;; esac

# Guarded CombinedPBMC legacy-raw migration: valid content and sidecar move to
# the canonical raw basename, with PATH rewritten and no duplicate left.
COMBINED_DIR="${TMP_DIR}/home/scratch/ECODA_paper/CombinedPBMC/data"
mkdir -p "${COMBINED_DIR}"
OLD="${COMBINED_DIR}/combined_pbmc_batch_effect_analysis.h5ad"
NEW="${COMBINED_DIR}/combined_pbmc.h5ad"
pixi run python -c 'import anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys; a=ad.AnnData(X=sp.csr_matrix([[1,0],[0,2]],dtype="float32"),obs=pd.DataFrame({"Sample":["s1","s2"],"cond":["Healthy","Healthy"],"batch":["A","B"]},index=["c1","c2"]),var=pd.DataFrame(index=["g1","g2"])); a.write_h5ad(sys.argv[1])' "${OLD}"
digest="$(md5sum "${OLD}" | cut -d' ' -f1)"
printf 'MD5=%s\nSIZE=%s\nPATH=%s\n' "${digest}" "$(wc -c < "${OLD}" | tr -d '[:space:]')" "${OLD}" > "${OLD}.md5"
rm -rf "${TMP_DIR}/home/scratch/ECODA_paper/_ecoda_owners"
: > "${CAPTURE}"
HOME="${TMP_DIR}/home" PATH="${TMP_DIR}/bin:${PATH}" USER_EMAIL="test@example.invalid" \
  STAGE2_SUBMITTER_TEST=1 bash "${ROOT}/src/2_dataset_specific_preprocessing/1_submit_hpc.sh" \
  --datasets CombinedPBMC --steps combinedpbmc >/dev/null
[[ -s "${NEW}" && -s "${NEW}.md5" ]]
[[ ! -e "${OLD}" && ! -e "${OLD}.md5" ]]
grep -q "^PATH=${NEW}$" "${NEW}.md5"

# Hook-level force propagation with a temporary slurm_config/python stub.
HOOK_ROOT="${TMP_DIR}/hook-project"
mkdir -p "${HOOK_ROOT}/src/2_dataset_specific_preprocessing" "${HOOK_ROOT}/bin"
cp "${ROOT}/src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh" "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/"
printf 'PROJECT_ROOT="%s"\nPYTHON_BIN="%s/bin/python"\n' "${HOOK_ROOT}" "${HOOK_ROOT}" > "${HOOK_ROOT}/src/slurm_config.sh"
printf '#!/bin/bash\nprintf "%%s\\n" "$*" > "%s/force.args"\n' "${HOOK_ROOT}" > "${HOOK_ROOT}/bin/python"
chmod +x "${HOOK_ROOT}/bin/python"
HOME="${TMP_DIR}/home" FORCE_PREPROCESS=1 bash "${HOOK_ROOT}/src/2_dataset_specific_preprocessing/1.5_submit_myocardial.sh"
case "$(cat "${HOOK_ROOT}/force.args")" in *"--force") ;; *) echo "myocardial hook dropped --force" >&2; exit 1 ;; esac
echo "stage2 submitter: OK"
