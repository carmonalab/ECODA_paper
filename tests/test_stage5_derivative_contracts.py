"""Exercise bound Stage5 roots and both language input resolvers with tiny files."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.utils.py.derivative_artifacts import bind, context_root, load
from src.utils.py.h5ad_source_identity import resolve_h5ad_path


def run():
    spec = importlib.util.spec_from_file_location("matrix_artifact_validator",
        ROOT / "src/5_run_benchmark_methods/matrix_artifact_validator.py")
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    for method, model in (("ecoda_unsupervised_harmony", "harmony_leiden_clr_v1"),
                          ("ecoda_unsupervised_harmony_limma", "harmony_leiden_clr_limma_v1")):
        identity = validator.build_batch_contract_identity(["batch"], sample_column="Sample",
            method_id=method, model_id=model)
        mode, formula = validator.batch_correction_spec_for_keys(method, ["batch"],
            effective_batch_keys=["batch"], non_estimable_batch_keys=[])
        harmony_only = method == "ecoda_unsupervised_harmony"
        recorded = dict(identity, effective_batch_keys=["batch"], non_estimable_batch_keys=[],
            correction_state="BATCH_CORRECTION", correction_mode=mode, correction_formula=formula,
            fixed_effect_aliases={} if harmony_only else {"batch": "batch_key_1"},
            correction_design_formula="~1" if harmony_only else "~1 + batch_key_1",
            design_rank=1 if harmony_only else 2, design_columns=1 if harmony_only else 2,
            design_residual_df=5 if harmony_only else 4)
        validator.validate_batch_contract_identity(identity, recorded, require_recorded=True,
            require_effective_metadata=True)
        wrong = dict(recorded, fixed_effect_aliases={"batch": "unexpected"})
        try:
            validator.validate_batch_contract_identity(identity, wrong, require_recorded=True,
                require_effective_metadata=True)
        except ValueError:
            pass
        else:
            raise AssertionError("wrong Harmony/limma aliases were accepted")
    with tempfile.TemporaryDirectory(prefix="ecoda-stage5-roots-") as temporary:
        base = Path(temporary).resolve()
        nas = base / "nas"; nas.mkdir()
        source_dir = base / "canonical"; source_dir.mkdir()
        cases = [("benchmark_analysis", "", "", "benchmark"),
            ("batch_effect_uncorrected", "uncorrected", "final", "batch_effect/uncorrected_final"),
            ("batch_effect_corrected", "corrected", "corrected_final", "batch_effect/corrected_final")]
        config = base / "datasets.json"
        config.write_text(json.dumps({"Toy": {"columns": {"sample": "sample", "label": "biology", "batch": "batch"},
            "views": {view: {"output_file_name": f"{view}.h5ad"} for view, *_ in cases}}}))
        for index, (view, analysis_pass, variant, context) in enumerate(cases):
            source = source_dir / f"{view}.h5ad"
            original = f"immutable {view} fixture".encode(); source.write_bytes(original)
            md5 = hashlib.md5(original).hexdigest(); proof = Path(f"{source}.md5")
            proof.write_text(f"MD5={md5}\nSIZE={len(original)}\nPATH={source}\n")
            inputs = base / f"inputs{index}.tsv"
            inputs.write_text(f"Toy\t{view}\t{source}\t{md5}\t{len(original)}\t{proof}\n")
            run_root = base / "scratch/_ecoda_runs" / f"stage5-{index}"
            manifests = run_root / "manifests"; manifests.mkdir(parents=True)
            (manifests / "source.manifest").write_text("source identity fixture\n")
            (manifests / "runtime.identity").write_text("runtime identity fixture\n")
            selection = manifests / "selection.tsv"; selection.write_text(f"Toy\t{view}\tpilot\n")
            binding = bind(run_root=run_root, release_id="test-release", nas_root=nas,
                stage="stage5", config_path=config, selection=selection, input_artifacts=inputs)
            value, _ = load(binding, config)
            root = context_root(value, context); nas_root = context_root(value, context, nas=True)
            env = dict(os.environ, PROJECT_ROOT=str(ROOT), DATASETS_JSON_FILE=str(config),
                PYTHON_BIN=sys.executable, HPC_SCRATCH_DIR=str(base / "scratch"), NAS_TARGET_DIR=str(nas),
                ECODA_RUN_ROOT=str(run_root), ECODA_ARTIFACT_BINDING=str(binding),
                ANALYSIS_ROOT=str(root), ANALYSIS_NAS_ROOT=str(nas_root),
                ANALYSIS_PASS=analysis_pass, ANALYSIS_VARIANT=variant,
                TEST_VIEW=view, TEST_SOURCE=str(source))
            shell = 'source "$PROJECT_ROOT/src/utils/bash/ecoda_run_common.sh"\n' \
                'source "$PROJECT_ROOT/src/utils/bash/ecoda_stage5_policy.sh"\n' \
                'ecoda_stage5_validate_identity "$ANALYSIS_PASS" "$ANALYSIS_VARIANT"\n' \
                'ecoda_stage5_validate_artifact_path "$ANALYSIS_ROOT/results/toy.rds"\n'
            result = subprocess.run(["bash", "-euo", "pipefail", "-c", shell], env=env,
                text=True, capture_output=True)
            assert result.returncode == 0, result.stderr
            audit_shell = shell + 'source "$PROJECT_ROOT/src/utils/bash/ecoda_run_audit_stage5.sh"\n' \
                '_audit_die() { echo "$*" >&2; return 1; }\n' \
                'STAGE_ARG=stage5\n' \
                'AUDIT_METADATA_PASS="$ANALYSIS_PASS"; AUDIT_METADATA_VARIANT="$ANALYSIS_VARIANT"\n' \
                'AUDIT_METADATA_ROOT="$ANALYSIS_ROOT"; AUDIT_METADATA_ANALYSIS_ROOT="$ANALYSIS_ROOT"\n' \
                'AUDIT_METADATA_ANALYSIS_NAS_ROOT="$ANALYSIS_NAS_ROOT"; AUDIT_METADATA_ANALYSIS_PASS="$ANALYSIS_PASS"\n' \
                'case "$ANALYSIS_VARIANT" in final) AUDIT_METADATA_ANALYSIS_LOG_PREFIX=execution_times_batch_effect_uncorrected_final_;; corrected_final) AUDIT_METADATA_ANALYSIS_LOG_PREFIX=execution_times_batch_effect_corrected_final_;; *) AUDIT_METADATA_ANALYSIS_LOG_PREFIX=execution_times_;; esac\n' \
                '_audit_stage5_identity\n'
            audit_result = subprocess.run(["bash", "-euo", "pipefail", "-c", audit_shell],
                env=env, text=True, capture_output=True)
            assert audit_result.returncode == 0, audit_result.stderr
            if analysis_pass == "corrected":
                composition_shell = shell + '_ecoda_stage5_artifacts_for Toy batch_effect_corrected composition\n' \
                    'for path in "${ECODA_BENCHMARK_ARTIFACTS[@]}"; do printf "%s\\t%s\\n" "$path" "$(ecoda_stage5_rds_artifact_method composition "$path")"; done\n'
                expansion = subprocess.run(["bash", "-euo", "pipefail", "-c", composition_shell],
                    env=env, text=True, capture_output=True)
                assert expansion.returncode == 0, expansion.stderr
                assert len(expansion.stdout.splitlines()) == 6
                assert [line.split("\t")[1] for line in expansion.stdout.splitlines()[2:]] == [
                    "ECODA_authors_HR", "ECODA_authors_HR_NULL", "ecoda_unsupervised_harmony", "ecoda_unsupervised_harmony_limma"]
            bad_env = dict(env, ANALYSIS_ROOT=str(base / "scratch/benchmark"))
            assert subprocess.run(["bash", "-euo", "pipefail", "-c", shell], env=bad_env,
                capture_output=True).returncode != 0
            with patch.dict(os.environ, env, clear=True):
                assert resolve_h5ad_path(base / "wrong", config, "Toy", view) == source
                if variant:
                    validator._validate_variant_root_binding(root, variant)
                    try:
                        validator._validate_variant_root_binding(base / "wrong", variant)
                    except ValueError:
                        pass
                    else:
                        raise AssertionError("unbound matrix artifact root accepted")
            r_code = 'source("src/5_run_benchmark_methods/benchmark_hpc_utils.R"); stopifnot(identical(get_h5ad_path(list(), "Toy", Sys.getenv("TEST_VIEW"), "/wrong"), Sys.getenv("TEST_SOURCE")))'
            result = subprocess.run([str(ROOT / ".pixi/envs/default/bin/Rscript"), "-e", r_code],
                cwd=ROOT, env=env, text=True, capture_output=True)
            assert result.returncode == 0, result.stderr
            artifact_name = f"Toy_batch_effect_{analysis_pass}_final_toy.rds" if variant else "toy.rds"
            artifact = root / "results" / artifact_name; artifact.parent.mkdir(parents=True)
            env["TEST_ARTIFACT"] = str(artifact)
            nas_artifact = nas_root / "results" / artifact_name; nas_artifact.parent.mkdir(parents=True)
            env["TEST_NAS_ARTIFACT"] = str(nas_artifact)
            r_boundary = '''
lines <- readLines("src/5_run_benchmark_methods/validate_benchmark_rds_contract.R")
start <- which(startsWith(lines, ".validate_variant_artifact_path <- function"))
end <- which(seq_along(lines) > start & lines == "}")[[1L]]
eval(parse(text=lines[start:end]))
validator_dir <- normalizePath("src/5_run_benchmark_methods")
config_path <- Sys.getenv("DATASETS_JSON_FILE")
batch_pass <- Sys.getenv("ANALYSIS_PASS"); batch <- nzchar(batch_pass); analysis_variant <- Sys.getenv("ANALYSIS_VARIANT")
analysis_root_env <- Sys.getenv("ANALYSIS_ROOT"); analysis_nas_root_env <- Sys.getenv("ANALYSIS_NAS_ROOT")
saveRDS(list(fixture="owned RDS root boundary"), Sys.getenv("TEST_ARTIFACT"))
stopifnot(file.copy(Sys.getenv("TEST_ARTIFACT"), Sys.getenv("TEST_NAS_ARTIFACT")))
.validate_variant_artifact_path(Sys.getenv("TEST_ARTIFACT"), "toy")
.validate_variant_artifact_path(Sys.getenv("TEST_NAS_ARTIFACT"), "toy")
analysis_root_env <- "/arbitrary/canonical/root"
stopifnot(inherits(try(.validate_variant_artifact_path(Sys.getenv("TEST_ARTIFACT"), "toy"), silent=TRUE), "try-error"))
stopifnot(inherits(try(.validate_variant_artifact_path(Sys.getenv("TEST_NAS_ARTIFACT"), "toy"), silent=TRUE), "try-error"))
'''
            result = subprocess.run([str(ROOT / ".pixi/envs/default/bin/Rscript"), "-e", r_boundary],
                cwd=ROOT, env=env, text=True, capture_output=True)
            assert result.returncode == 0, result.stderr
            assert source.read_bytes() == original
    print("PASS: three Stage5 derivative contexts, Python/R input references, and strict root rejection")


if __name__ == "__main__":
    run()
