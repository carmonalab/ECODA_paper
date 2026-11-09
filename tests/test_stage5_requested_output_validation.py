"""Final produced-output reads are distinct from strict reusable-input reads."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind


def proof(path):
    Path(f"{path}.md5").write_text(
        f"MD5={hashlib.md5(path.read_bytes()).hexdigest()}\nSIZE={path.stat().st_size}\nPATH={path}\n")


def run():
    with tempfile.TemporaryDirectory(prefix="ecoda-produced-output-") as temporary:
        base = Path(temporary).resolve()
        scratch = base / "scratch"; nas = base / "nas"; nas.mkdir()
        run_root = scratch / "_ecoda_runs/same_run"
        manifests = run_root / "manifests"; manifests.mkdir(parents=True)
        (run_root / "metadata").write_text("STAGE=stage5\nRUN_ID=same_run\nSTATE=ACTIVE\n")
        config = base / "datasets.json"
        config.write_text(json.dumps({"Bassez": {"columns": {"sample": "Sample", "label": "biology",
            "cell_type_high_res": "HR"}, "views": {"benchmark_analysis": {"output_file_name": "source.h5ad"}}}}))
        source = base / "source.h5ad"; source.write_bytes(b"metadata identity fixture"); proof(source)
        (manifests / "source.manifest").write_text("source fixture\n")
        (manifests / "runtime.identity").write_text("runtime fixture\n")
        selection = manifests / "selection.tsv"
        selection.write_text("Bassez\tbenchmark_analysis\tcomposition\n")
        (manifests / "requested_selection.tsv").write_text(
            "Bassez\tbenchmark_analysis\tcomposition\tECODA_authors_HR\n")
        references = base / "inputs.tsv"
        references.write_text(f"Bassez\tbenchmark_analysis\t{source}\t{hashlib.md5(source.read_bytes()).hexdigest()}\t{source.stat().st_size}\t{source}.md5\n")
        binding = bind(run_root=run_root, release_id="test-release", nas_root=nas, stage="stage5",
            config_path=config, selection=selection, input_artifacts=references)
        identity = manifests / "source_identity.json"
        identity.write_text(json.dumps({"schema": 1, "entries": [{"dataset": "Bassez", "view": "benchmark_analysis",
            "path": str(source), "size": source.stat().st_size, "md5": hashlib.md5(source.read_bytes()).hexdigest(),
            "sample_ids": [f"s{i}" for i in range(1, 13)]}]})); proof(identity)
        lines = (ROOT / "src/5_run_benchmark_methods/stage5_dispatcher.sh").read_text().splitlines()
        functions = []
        for name in ("benchmark_artifacts_for", "stage5_recorded_producer", "stage5_allowed_producers_for",
                     "stage5_producer_allowed", "stage5_validate_path_component", "stage5_validate_reusable_artifact", "stage5_watchdog_status_value",
                     "stage5_watchdog_safe_label", "stage5_watchdog_state_for_method",
                     "stage5_validate_requested_output_context", "stage5_validate_requested_output_artifact",
                     "stage5_validate_requested_rds"):
            start = lines.index(f"{name}() {{")
            end = next(i for i in range(start + 1, len(lines)) if lines[i] == "}")
            functions.extend(lines[start:end + 1])
        (base / "functions.sh").write_text("\n".join(functions) + "\n")
        # Real selected composition orchestration, CLR math, scoring, atomic
        # persistence and producer records. No orchestration/validator stubs.
        (base / "produce.R").write_text(r'''
a <- commandArgs(TRUE)
setwd(a[1]); suppressPackageStartupMessages(library(dplyr))
source("src/utils/load_worker_functions.R")
source("src/5_run_benchmark_methods/benchmark_hpc_utils.R")
options(mc.cores=1)
obs <- do.call(rbind, lapply(1:12, function(i) data.frame(
  Sample=paste0("s",i), biology=if(i<=6) "A" else "B",
  HR=rep(c("T","B","M","NK"), c(4+i, 18-i, 5+i%%4, 3+i%%5)))))
labels <- setNames(factor(ifelse(1:12<=6,"A","B")),paste0("s",1:12))
invisible(run_composition_methods_hpc(labels,unique(obs[,c("Sample","biology")]),NULL,NULL,
  obs,"biology",ct_col_high_res="HR",results_dir=a[2],ds="Bassez",
  variants="ECODA_authors_HR"))
''')
        (base / "exercise.sh").write_text(r'''#!/bin/bash
set -euo pipefail
source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"
source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"
source "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"
source "$FIXTURE/functions.sh"
ecoda_open_run same_run stage5
RUN_ID=same_run; PASS_ARG=""; ANALYSIS_VARIANT=""; ANALYSIS_PASS=""
ANALYSIS_ROOT="$ECODA_RUN_ROOT/artifacts/benchmark"
ANALYSIS_NAS_ROOT="$NAS_TARGET_DIR/releases/test-release/$RUN_ID/benchmark"
export ANALYSIS_ROOT ANALYSIS_NAS_ROOT ANALYSIS_PASS ANALYSIS_VARIANT
SCRIPT_DIR="$PROJECT_ROOT/src/5_run_benchmark_methods"
SOURCE_IDENTITY="$ECODA_RUN_ROOT/manifests/source_identity.json"
DISPATCH_MANIFEST="$ECODA_RUN_ROOT/manifests/selection.tsv"
PENDING_SELECTION="$ECODA_RUN_ROOT/manifests/pending_selection.tsv"
SCHEDULER_FILE="$ECODA_RUN_ROOT/manifests/scheduler_ids.tsv"
cp "$DISPATCH_MANIFEST" "$PENDING_SELECTION"; ecoda_write_checksum "$PENDING_SELECTION"
mkdir -p "$ANALYSIS_ROOT/results" "$ECODA_RUN_ROOT/status/watchdogs"
benchmark_artifacts_for Bassez benchmark_analysis composition
for path in "${ARTIFACT_PATHS[@]}"; do
  ecoda_artifact_owner_acquire "$path" stage5 "$RUN_ID" 0 0 0 >/dev/null
done
export ECODA_ARTIFACT_PRODUCER=stage5_composition
"$PIXI_RSCRIPT" "$FIXTURE/produce.R" "$PROJECT_ROOT" "$ANALYSIS_ROOT/results" > "$FIXTURE/producer.log" 2>&1
path="$ANALYSIS_ROOT/results/Bassez_metadata.rds"
if stage5_validate_reusable_artifact "$path" composition > "$FIXTURE/reuse.log" 2>&1; then
  echo 'reuse input guard accepted ACTIVE output' >&2; exit 1
fi
grep -q 'input path has an active writer same_run' "$FIXTURE/reuse.log"
label=benchmark_analysis__composition__independent__parameters
matrix="$ECODA_RUN_ROOT/manifests/matrix_benchmark_analysis_composition_independent_parameters.tsv"
printf 'Bassez\tbenchmark_analysis\tcomposition\tECODA_authors_HR\n' > "$matrix"
printf 'LABEL=%s\nSTATE=OK\nSCHEDULER_ID=102\n' "$label" > "$ECODA_RUN_ROOT/status/watchdogs/$label.status"
printf 'ARRAY\t101\nWATCHDOG\t102\nAGGREGATE_GATE\t103\n' > "$SCHEDULER_FILE"
printf 'STATE=OK\nWATCHDOG_LABELS=%s\nSCHEDULER_ID=101\nSCHEDULER_ID=102\nSCHEDULER_ID=103\n' "$label" > "$ECODA_RUN_ROOT/status/aggregate"
stage5_validate_requested_rds > "$FIXTURE/semantic.log" 2>&1
# Real CLI validated metadata and the composition bundle; ownership remains ACTIVE.
grep -q 'artifact-list contract OK' "$FIXTURE/semantic.log"
for path in "${ARTIFACT_PATHS[@]}"; do
  ecoda_artifact_owner_validate "$path" "$RUN_ID" >/dev/null
  [[ "$ECODA_ARTIFACT_OWNER_STATE" == ACTIVE ]]
done
# A later same-run validation can encounter already-finalized OK owners.
# Preserve its existing reusable-artifact route without treating FAIL as reuse.
for path in "${ARTIFACT_PATHS[@]}"; do ecoda_artifact_owner_set_state "$path" OK fixture >/dev/null; done
stage5_validate_requested_rds > "$FIXTURE/semantic-terminal.log" 2>&1
record="$(ecoda_artifact_record_path "$ANALYSIS_ROOT/results/Bassez_metadata.rds" "$RUN_ID")"
mv "$record" "$record.absent"
if stage5_validate_requested_rds > "$FIXTURE/terminal-missing-record.log" 2>&1; then
  echo 'terminal pending output accepted without its producer record' >&2; exit 1
fi
mv "$record.absent" "$record"
for path in "${ARTIFACT_PATHS[@]}"; do ecoda_artifact_owner_set_state "$path" ACTIVE fixture >/dev/null; done
path="$ANALYSIS_ROOT/results/Bassez_metadata.rds"
reject() { if "$@" > "$FIXTURE/reject.log" 2>&1; then echo "unexpected success: $*" >&2; exit 1; fi; }
# Matching terminal compute is mandatory, not artifact presence alone.
cp "$ECODA_RUN_ROOT/status/aggregate" "$FIXTURE/aggregate.saved"
sed 's/STATE=OK/STATE=FAIL/' "$FIXTURE/aggregate.saved" > "$ECODA_RUN_ROOT/status/aggregate"
reject stage5_validate_requested_rds
cp "$FIXTURE/aggregate.saved" "$ECODA_RUN_ROOT/status/aggregate"
sed 's/SCHEDULER_ID=103/SCHEDULER_ID=999/' "$FIXTURE/aggregate.saved" > "$ECODA_RUN_ROOT/status/aggregate"
reject stage5_validate_requested_rds
cp "$FIXTURE/aggregate.saved" "$ECODA_RUN_ROOT/status/aggregate"
watchdog="$ECODA_RUN_ROOT/status/watchdogs/$label.status"
cp "$watchdog" "$FIXTURE/watchdog.saved"; sed 's/STATE=OK/STATE=FAIL/' "$FIXTURE/watchdog.saved" > "$watchdog"
reject stage5_validate_requested_rds
cp "$FIXTURE/watchdog.saved" "$watchdog"
owner="$(ecoda_artifact_owner_dir "$path")/owner"
cp "$owner" "$FIXTURE/owner.saved"
sed 's/RUN_ID=same_run/RUN_ID=foreign_run/' "$FIXTURE/owner.saved" > "$owner"
reject stage5_validate_requested_rds
cp "$FIXTURE/owner.saved" "$owner"
sed 's/STAGE=stage5/STAGE=stage3/' "$FIXTURE/owner.saved" > "$owner"
reject stage5_validate_requested_rds
cp "$FIXTURE/owner.saved" "$owner"
ecoda_artifact_owner_set_state "$path" FAIL fixture >/dev/null
reject stage5_validate_requested_rds
cp "$FIXTURE/owner.saved" "$owner"
# Even an exact-run owned artifact outside this requested expansion is forbidden.
reject stage5_validate_requested_output_artifact Bassez benchmark_analysis composition "$ANALYSIS_ROOT/results/unrequested.rds"
record="$(ecoda_artifact_record_path "$path" "$RUN_ID")"
cp "$record" "$FIXTURE/record.saved"
sed 's/PRODUCER=stage5_composition/PRODUCER=stage5_pseudobulk/' "$FIXTURE/record.saved" > "$record"
reject stage5_validate_requested_rds
cp "$FIXTURE/record.saved" "$record"
mv "$record" "$record.absent"; reject stage5_validate_requested_rds; mv "$record.absent" "$record"
cp "$PENDING_SELECTION" "$FIXTURE/pending.saved"
printf 'Bassez\tbenchmark_analysis\tpseudobulk\n' >> "$PENDING_SELECTION"
reject stage5_validate_requested_rds
cp "$FIXTURE/pending.saved" "$PENDING_SELECTION"
# Non-pending/reused rows continue through the unchanged strict input guard.
printf 'Other\tbenchmark_analysis\tcomposition\n' > "$PENDING_SELECTION"; ecoda_write_checksum "$PENDING_SELECTION"
reject stage5_validate_requested_rds
grep -q 'active writer same_run' "$FIXTURE/reject.log"
cp "$FIXTURE/pending.saved" "$PENDING_SELECTION"; ecoda_write_checksum "$PENDING_SELECTION"
cp "$path" "$FIXTURE/payload.saved"; chmod u+w "$path"; printf 'corrupt' >> "$path"
reject stage5_validate_requested_rds
cp "$FIXTURE/payload.saved" "$path"; chmod a-w "$path"
# Valid owner/producer/checksum cannot excuse malformed scientific metadata.
cp "$path.md5" "$FIXTURE/proof.saved"
chmod u+w "$path"
"$PIXI_RSCRIPT" -e 'a<-commandArgs(TRUE);x<-readRDS(a[1]);x$n_samples<-0;saveRDS(x,a[1])' "$path"
ecoda_write_checksum "$path"; rm "$record"
ecoda_write_artifact_record "$path" stage5_composition "$RUN_ID" >/dev/null
reject stage5_validate_requested_rds
grep -q 'metadata numeric fields are malformed' "$FIXTURE/reject.log"
chmod u+w "$path"; cp "$FIXTURE/payload.saved" "$path"; chmod a-w "$path"
cp "$FIXTURE/proof.saved" "$path.md5"; cp "$FIXTURE/record.saved" "$record"
stage5_validate_requested_rds > "$FIXTURE/semantic-final.log" 2>&1
# No validation promoted any owner or created any NAS output.
ecoda_artifact_owner_validate "$path" "$RUN_ID" >/dev/null
[[ "$ECODA_ARTIFACT_OWNER_STATE" == ACTIVE ]]
[[ -z "$(find "$NAS_TARGET_DIR" -type f -print)" ]]
printf 'PASS: actual producer + ACTIVE output owners + strict final RDS CLI; foreign/reused/tampered states rejected\n'
''')
        environment = dict(os.environ, PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(config),
            PYTHON_BIN=sys.executable, RETICULATE_PYTHON=sys.executable, RETICULATE_USE_MANAGED_VENV="no",
            HPC_SCRATCH_DIR=str(scratch), NAS_TARGET_DIR=str(nas), ECODA_RUNS_ROOT=str(scratch / "_ecoda_runs"),
            ECODA_OWNERS_ROOT=str(scratch / "_ecoda_owners"), ECODA_RUN_ROOT=str(run_root), ECODA_RUN_ID="same_run",
            ECODA_ARTIFACT_BINDING=str(binding), PIXI_RSCRIPT=str(ROOT / ".pixi/envs/default/bin/Rscript"),
            FIXTURE=str(base))
        completed = subprocess.run(["bash", str(base / "exercise.sh")], env=environment, capture_output=True, text=True)
        if completed.returncode:
            details = "\n".join(f"{p.name}:\n{p.read_text()}" for p in base.glob("*.log"))
            raise AssertionError(f"{completed.returncode}\n{completed.stdout}\n{completed.stderr}\n{details}")
        print(completed.stdout.strip())


if __name__ == "__main__":
    run()
