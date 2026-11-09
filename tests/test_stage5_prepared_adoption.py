"""Actual tiny cache producer, current RDS CLI and ownership/adoption boundaries."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prepared_metadata_fixture", ROOT / "tests/test_stage5_prepared_metadata.py")
metadata = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metadata)


def run(command, environment, valid=True):
    result = subprocess.run(command, env=environment, text=True, capture_output=True)
    assert (result.returncode == 0) == valid, (result.returncode, result.stdout, result.stderr)
    return result


def receipt(fixture):
    fixture["receipt"].write_text(json.dumps(fixture["value"], sort_keys=True) + "\n")
    return metadata.digest(fixture["receipt"])


def shell(fixture, context, body, valid=True):
    row = fixture["contexts"][context]
    environment = dict(row["env"], ECODA_AUX_ROOT=str(ROOT / "aux"),
        PIXI_RSCRIPT=shutil.which("Rscript"), RUN_ID=row["new_run"].name,
        SOURCE_IDENTITY=str(row["new_run"] / "manifests/source_identity.json"))
    code = 'set -euo pipefail\nsource "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n' \
        'source "$PROJECT_ROOT/src/utils/bash/ecoda_derivative_artifacts.sh"\n' \
        'source "$PROJECT_ROOT/src/5_run_benchmark_methods/stage5_prepared_adoption.sh"\n' + body
    return run(["bash", "-c", code, "fixture", str(fixture["receipt"]), receipt(fixture), context], environment, valid)


def real_caches(fixture):
    row = fixture["contexts"]["benchmark"]
    environment = dict(row["env"], ECODA_RUN_ROOT=str(row["old_run"]), ECODA_RUN_ID=row["old_run"].name,
        ECODA_AUX_ROOT=str(ROOT / "aux"), PROJECT_ROOT=str(ROOT))
    script = fixture["base"] / "actual_cache_producer.R"
    script.write_text(r'''
a<-commandArgs(TRUE);root<-a[1];old<-a[2]
e<-new.env(parent=globalenv())
for(p in c('src/utils/pseudobulk.R','src/utils/artifact_records.R','src/5_run_benchmark_methods/benchmark_pipeline.R','src/5_run_benchmark_methods/benchmark_hpc_utils.R'))sys.source(file.path(root,p),e)
genes<-c('XIST','SRY','retained_A','retained_B','retained_C')
counts<-matrix(seq_len(10),nrow=5,dimnames=list(genes,c('Stephenson_s1','Stephenson_s2')))
metadata<-data.frame(Sample=colnames(counts),row.names=colnames(counts))
e$load_h5ad_sample_aggregate<-function(...,metadata_columns)list(counts=counts,metadata=metadata,sample_ids=metadata$Sample)
e$fit_pseudobulk_deseq2<-function(counts,...)list(norm_matrix=counts,variance_order=rownames(counts))
e$exec_time<-function(expr){force(expr);as.difftime(1,units='secs')}
e$peak_rss_gb<-function()0.25
input<-read.delim(file.path(old,'manifests/input_artifacts.tsv'),header=FALSE)[1,3]
values<-e$prepare_pseudobulks_hpc(input,hvg_rank_genes=genes,variants=e$PB_VARIANT_NAMES,cache_stem='Stephenson',run_id=basename(old))
for(v in names(values))e$save_rds_atomic(values[[v]],file.path(old,'artifacts/benchmark/pseudobulks',paste0('Stephenson_pseudobulk_',v,'.rds')),producer=paste0('stage5_prepare_pseudobulk_',v))
log<-file.path(old,'logs/execution_times_prepare_pseudobulk_Stephenson.feather')
unlink(c(log,paste0(log,'.md5')))
e$emit_pseudobulk_timing_rows(values,'Stephenson',log_file=log)
''')
    run([shutil.which("Rscript"), str(script), str(ROOT), str(row["old_run"])], environment)
    for artifact in fixture["artifacts"]:
        if artifact["context"] != "benchmark":
            continue
        path = Path(artifact["source"])
        record = row["old_run"] / "manifests/artifacts" / (hashlib.sha256(str(path).encode()).hexdigest()[:32] + ".record")
        artifact.update(md5=hashlib.md5(path.read_bytes()).hexdigest(), size=path.stat().st_size,
            sidecar_sha256=metadata.digest(Path(str(path) + ".md5")), record_sha256=metadata.digest(record))
        if artifact["kind"] == "pseudobulk":
            # The actual failed controller may have marked valid preparation
            # owners FAIL. A reviewed copy never promotes those old owners.
            metadata.common(environment,
                'ecoda_artifact_owner_set_state "$1" FAIL original-controller-failed >/dev/null\n', path)
            owner = metadata.common(environment, 'ecoda_artifact_owner_dir "$1"\n', path).stdout.strip()
            artifact["owner_sha256"] = metadata.digest(Path(owner) / "owner")
    identity = row["new_run"] / "manifests/source_identity.json"
    source = Path(row["inputs"]).read_text().splitlines()[0].split("\t")
    identity.write_text(json.dumps(dict(schema=1, entries=[dict(dataset="Stephenson", view="benchmark_analysis",
        path=source[2], md5=source[3], size=source[4], sample_ids=["Stephenson_s1", "Stephenson_s2"])])))
    metadata.proof(identity)
    receipt(fixture)


def main():
    with tempfile.TemporaryDirectory(prefix="ecoda-prepared-cache-") as temporary:
        fixture = metadata.build_fixture(Path(temporary))
        real_caches(fixture)
        protected = {}
        for item in fixture["artifacts"]:
            for path in (Path(item["source"]), Path(item["source"] + ".md5")):
                protected[path] = (path.read_bytes(), path.stat().st_mode)
        for path in fixture["scratch"].joinpath("_ecoda_owners").glob("artifact/*/owner"):
            protected[path] = (path.read_bytes(), path.stat().st_mode)
        original = copy.deepcopy(fixture["value"])
        # Missing receipt proof and invalid selected scope fail before any copy.
        fixture["value"]["artifacts"].pop()
        shell(fixture, "benchmark", 'stage5_adopt_prepared_artifacts "$1" "$2" "$3"', False)
        # Failure froze only new small proof artifacts; use a fresh target run
        # in this fixture by removing those diagnostic controls, never old data.
        new = fixture["contexts"]["benchmark"]["new_run"]
        for name in ("prepared_adoption.receipt.json", "prepared_adoption.receipt.json.md5", "prepared_adoption.tsv.build"):
            for path in new.joinpath("manifests").glob(name + "*"):
                path.unlink()
        # Existing receipt record must also be removed only in this test reset.
        receipt_record = new / "manifests/artifacts" / (hashlib.sha256(str(new / "manifests/prepared_adoption.receipt.json").encode()).hexdigest()[:32] + ".record")
        receipt_record.unlink()
        fixture["value"] = original
        fixture["artifacts"] = fixture["value"]["artifacts"]
        shell(fixture, "benchmark", 'stage5_adopt_prepared_artifacts "$1" "$2" "$3"')
        report = json.loads((new / "manifests/prepared_adoption.json").read_text())
        assert report["operation"] == "COPY_WITHOUT_REFIT" and report["original_gate_state"] == "FAILED"
        assert len(report["artifacts"]) == 7
        # Exact current cache-reader reuse must succeed, preserving old timing
        # IDs. Its fitter/aggregate functions are replaced with hard failures.
        verify = fixture["base"] / "validate_copied_caches.R"
        verify.write_text(r'''
a<-commandArgs(TRUE);root<-a[1];new<-a[2]
e<-new.env(parent=globalenv())
for(p in c('src/utils/pseudobulk.R','src/utils/artifact_records.R','src/5_run_benchmark_methods/benchmark_hpc_utils.R'))sys.source(file.path(root,p),e)
e$fit_pseudobulk_deseq2<-function(...)stop('REFIT FORBIDDEN')
e$load_h5ad_sample_aggregate<-function(...)stop('COUNT READ FORBIDDEN')
for(v in e$PB_VARIANT_NAMES){p<-file.path(new,'artifacts/benchmark/pseudobulks',paste0('Stephenson_pseudobulk_',v,'.rds'));stopifnot(e$.pb_variant_cache_valid(p,v));x<-e$read_rds_checked(p,producer=e$PB_VARIANT_PRODUCERS[[v]]);stopifnot(startsWith(x$timing_id,'old_benchmark:'))}
''')
        environment = dict(fixture["contexts"]["benchmark"]["env"], ECODA_AUX_ROOT=str(ROOT / "aux"))
        run([shutil.which("Rscript"), str(verify), str(ROOT), str(new)], environment)
        # Exercise the real dispatcher reusable-preparation branch, rather
        # than replacing it with a cache-valid stub. It must need no watchdog.
        lines = (ROOT / "src/5_run_benchmark_methods/stage5_dispatcher.sh").read_text().splitlines()
        functions = []
        for name in ("stage5_validate_path_component", "benchmark_artifacts_for", "benchmark_rds_group_valid",
                     "stage5_validate_reusable_artifact", "benchmark_selected_artifacts_valid",
                     "stage5_recorded_producer", "stage5_allowed_producers_for", "stage5_producer_allowed",
                     "stage5_requested_preparation_dependency"):
            start = lines.index(name + "() {")
            end = next(index for index in range(start + 1, len(lines)) if lines[index] == "}")
            functions.extend(lines[start:end + 1])
        body = ('source "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"\n'
                'SCRIPT_DIR="$PROJECT_ROOT/src/5_run_benchmark_methods"\n'
                'PASS_ARG=""; ANALYSIS_PASS=""; ANALYSIS_VARIANT=""; ECODA_RECOVERY_ACTIVE=0\n'
                'METHOD_MATRIX_MODE=0; RDS_PREFLIGHT_DONE=""; RDS_PREFLIGHT_FAILED=""\n'
                'METHODS=(pseudobulk); ANALYSES=(_ecoda_none_); PREP_VIEWS=(); PREP_DATASETS=(); PREP_WATCHDOGS=()\n'
                'PENDING_SELECTION="$ECODA_RUN_ROOT/manifests/test_pending.tsv"\n'
                'printf "Stephenson\tbenchmark_analysis\tpseudobulk\\n" > "$PENDING_SELECTION"\n'
                + "\n".join(functions)
                + '\nbenchmark_selected_artifacts_valid Stephenson benchmark_analysis prepare_pseudobulk\n'
                'stage5_requested_preparation_dependency Stephenson benchmark_analysis\n'
                '[[ -z "$STAGE5_PREPARATION_DEPENDENCY" ]]\n')
        shell(fixture, "benchmark", body)
        shell(fixture, "benchmark", 'stage5_adopt_prepared_artifacts "$1" "$2" "$3"', False)
        for path, before in protected.items():
            assert (path.read_bytes(), path.stat().st_mode) == before, path
        print("PASS: real six-cache producer/schema/BL policy and copy-only reuse, FAIL old owners preserved, original timing retained, no fits/count reads, exact scope and create-only rejection")


if __name__ == "__main__":
    main()
