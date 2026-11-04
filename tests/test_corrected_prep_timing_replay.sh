#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
"${ROOT}/.pixi/envs/default/bin/python" - "${ROOT}" <<'PY'
import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pandas as pd

root = Path(sys.argv[1])
entry = root / 'src/5_run_benchmark_methods/finalize_corrected_outputs.sh'
body = entry.read_text().split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
ns = {'__name__': 'timing_replay_fixture'}
exec(compile(body, str(entry), 'exec'), ns)
main = next(node for node in ast.parse(body).body if isinstance(node, ast.FunctionDef) and node.name == 'main')
for name in ('promotion_sources', 'execution_producer'):
    node = next(node for node in main.body if isinstance(node, ast.FunctionDef) and node.name == name)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), str(entry), 'exec'), ns)
# Execute the real post-validation shard-selection loop, not a copy of it.
start = next(i for i, node in enumerate(main.body) if isinstance(node, ast.Assign)
             and any(isinstance(target, ast.Name) and target.id == 'source_shards' for target in node.targets))
select_logs = compile(ast.fix_missing_locations(ast.Module(body=main.body[start:start + 2], type_ignores=[])), str(entry), 'exec')
rscript = str(root / '.pixi/envs/default/bin/Rscript')
os.environ.update(PROJECT_ROOT=str(root), PIXI_RSCRIPT=rscript + ' --vanilla')
with tempfile.TemporaryDirectory(prefix='ecoda-prep-timing-') as temporary:
    scratch = Path(temporary).resolve()
    runs = scratch / '_ecoda_runs'
    old = runs / 'original'
    report = runs / 'recovery_promotion'
    for path in (old / 'manifests', old / 'logs', report / 'manifests', report / 'logs'):
        path.mkdir(parents=True)
    (old / 'metadata').write_text('RUN_ID=original\nSTATE=ACTIVE\n')
    cache = scratch / 'Diabetes_pseudobulk_hvg2000.rds'
    os.environ.update(HPC_SCRATCH_DIR=str(scratch), ECODA_RUNS_ROOT=str(runs),
                      ECODA_RUN_ROOT=str(report), ECODA_RUN_ID=report.name)
    create = '''
source(file.path(Sys.getenv("PROJECT_ROOT"), "src/5_run_benchmark_methods/benchmark_hpc_utils.R"))
x <- list(pb=matrix(1, 2, 2), time_secs=0.25, mem_GB=1.5,
          aggregate_time_secs=2, shared_fit_time_secs=8, shared_time_secs=10,
          variant_time_secs=0.25, shared_mem_GB=1.5,
          timing_id="original:Diabetes:batch_effect_corrected:corrected", timing_schema=2L)
save_rds_atomic(x, commandArgs(TRUE)[[1]], producer="stage5_prepare_pseudobulk_hvg2000", run_id="original")
'''
    subprocess.run([rscript, '--vanilla', '-e', create, str(cache)], check=True)
    protected = {path: path.read_bytes() for path in scratch.rglob('*') if path.is_file()}
    log_module = ns['module']('timing_replay_logs', root / 'src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.2_merge_execution_times.py')
    ns.update(root=root, runs=runs, report=report, run_id=report.name,
              logical={('Diabetes', 'prepare_pseudobulk'): 'ADOPT'},
              decisions=[['Diabetes', ns['VIEW'], 'prepare_pseudobulk', 'ADOPT', str(cache), str(cache), 'original']],
              row_paths={('Diabetes', 'prepare_pseudobulk'): [cache]}, log_module=log_module, evidence=set())
    def shell(code, *args):
        return subprocess.run(['bash', '-c', 'source "$1/src/utils/bash/ecoda_run_common.sh"; shift; ' + code,
                               'fixture', str(root), *map(str, args)], check=True, text=True, stdout=subprocess.PIPE).stdout
    ns['shell'] = shell
    for _ in range(2):
        exec(select_logs, ns)
        target = report / 'logs/execution_times_batch_effect_corrected_final_prepare_pseudobulk_Diabetes.feather'
        measured = log_module._read_recorded_execution_log(target, producer='stage5_execution_log', run_id=report.name, require_record=True)
        actual = measured.set_index('method')[['dataset', 'time_secs', 'mem_GB']].sort_index()
        expected = pd.DataFrame({'dataset': ['Diabetes', 'Diabetes'], 'time_secs': [0.25, 10.0], 'mem_GB': [1.5, 1.5]},
                                index=['prepare_pseudobulk_hvg2000', 'prepare_pseudobulk_shared'])
        expected.index.name = 'method'
        pd.testing.assert_frame_equal(actual, expected)
        assert all(path.read_bytes() == value for path, value in protected.items()), 'cache or original evidence changed'
        assert not list((old / 'logs').iterdir()), 'historical log was fabricated'
    ns['logical'] = {('Diabetes', 'pseudobulk'): 'REUSE'}
    ns['decisions'] = [['Diabetes', ns['VIEW'], 'pseudobulk', 'REUSE', str(cache), str(cache), 'original']]
    try:
        exec(select_logs, ns)
    except ValueError as exc:
        assert 'missing original execution shard' in str(exc)
    else:
        raise AssertionError('missing result-method timing silently accepted')
print('PASS cache-derived measured preparation timing, repeat recovery, unchanged source and missing-result-log rejection')
PY
