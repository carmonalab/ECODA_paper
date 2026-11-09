"""No models/data: exercise actual resource functions and isolated device guards."""
import ast
import contextlib
import io
import re
import subprocess
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
DISPATCHER = ROOT / "src/5_run_benchmark_methods/stage5_dispatcher.sh"
WORKER = ROOT / "src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.1.1_benchmark_methods_py.py"


def allocation(view, policy, resource="base", partition=""):
    text = DISPATCHER.read_text()
    functions = "\n".join(re.search(r"^" + name + r"\(\) \{.*?^\}", text, re.M | re.S)[0]
                          for name in ("gpu_method_is_default", "method_spec", "method_spec_for_group"))
    script = r'''
set -eu
PROJECT_ROOT=/immutable/source
BENCHMARK_CPU_CONSTRAINT=EPYC-7742
BENCHMARK_GPU_CONSTRAINT=nvidia_h200_nvl
BENCHMARK_CPU_CPUS_PER_TASK=16
BENCHMARK_GPU_CPUS_PER_TASK=8
SLURM_PARTITION_BENCHMARK_CPU=cpu
MAX_NUM_CHUNKS_PARALLEL=40
BENCHMARK_CPU_TIME_LIMIT=12:00:00
BENCHMARK_GPU_DEFAULT_METHODS=mrvi,scpoli
BENCHMARK_GPU_COUNT=2
BENCHMARK_GPU_DEFAULT_PARTITION=gpu-default
BENCHMARK_GPU_ANY_PARTITION=gpu,private-gpu
BENCHMARK_GPU_ANY_VRAM_PER_GPU=80G
BENCHMARK_GPU_DEFAULT_ARRAY_THROTTLE=4
BENCHMARK_GPU_ANY_ARRAY_THROTTLE=3
BENCHMARK_GPU_DEFAULT_TIME_LIMIT=10:00:00
BENCHMARK_GPU_ANY_TIME_LIMIT=11:00:00
GPU_POLICY="$1"; PASS_ARG="$2"; PARTITION_ARG="$3"
''' + functions + r'''
method_spec_for_group mrvi "$4" "$5"
printf '%s\n' "$METHOD_GPU_POLICY" "$METHOD_RUNTIME_NV" "$METHOD_PARTITION" "$METHOD_TIME_LIMIT" "${METHOD_FLAGS[@]}"
'''
    result = subprocess.run(["/bin/bash", "-c", script, "policy", policy,
                             view.removeprefix("batch_effect_") if view.startswith("batch_effect_") else "",
                             partition, resource, view], text=True, capture_output=True)
    return result


for view in ("batch_effect_uncorrected", "batch_effect_corrected"):
    cpu = allocation(view, "auto"); assert cpu.returncode == 0, cpu.stderr
    assert cpu.stdout.splitlines()[:4] == ["cpu", "0", "cpu", "12:00:00"]
    assert cpu.stdout.splitlines()[4:] == ["--cpus-per-task=16"]
    any_gpu = allocation(view, "any"); assert any_gpu.returncode == 0, any_gpu.stderr
    assert any_gpu.stdout.splitlines() == ["any", "1", "gpu,private-gpu", "11:00:00",
                                           "--gres=gpu:2,VramPerGpu:80G", "--cpus-per-task=8"]
    default = allocation(view, "default"); assert default.returncode == 0, default.stderr
    assert default.stdout.splitlines() == ["default", "1", "gpu-default", "10:00:00",
                                          "--gpus=2", "--cpus-per-task=8"]
    explicit = allocation(view, "any", partition="custom-gpu")
    assert explicit.returncode == 0 and explicit.stdout.splitlines()[2] == "custom-gpu"
    assert "--gres=gpu:2,VramPerGpu:80G" in explicit.stdout
    assert allocation(view, "auto", partition="custom-gpu").returncode != 0

benchmark = allocation("benchmark_analysis", "auto")
assert benchmark.returncode == 0 and benchmark.stdout.splitlines()[:2] == ["default", "1"]
assert "--constraint=nvidia_h200_nvl" in benchmark.stdout
benchmark_cpu = allocation("benchmark_analysis", "auto", resource="cpu")
assert benchmark_cpu.returncode == 0 and benchmark_cpu.stdout.splitlines()[:2] == ["cpu", "0"]
assert "--constraint=EPYC-7742" in benchmark_cpu.stdout

tree = ast.parse(WORKER.read_text())
def function(name):
    return next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
available = True
cuda = SimpleNamespace(is_available=lambda: available, current_device=lambda: 0,
                       get_device_name=lambda index: "Mock GPU")
ns = {"torch": SimpleNamespace(cuda=cuda), "GPU_BACKED_METHODS": frozenset(("mrvi", "scpoli"))}
exec(compile(ast.Module(body=[function("validate_gpu_execution")], type_ignores=[]), str(WORKER), "exec"), ns)
validate = ns["validate_gpu_execution"]
def rejects(*args, **kwargs):
    try: validate(*args, **kwargs)
    except RuntimeError: return
    raise AssertionError("Invalid/implicit device was accepted")
for view in ("batch_effect_uncorrected", "batch_effect_corrected"):
    validate("mrvi", "cpu", view=view)
    validate("mrvi", "cuda", view=view)
    rejects("mrvi", "auto", view=view)
    available = False
    validate("mrvi", "cpu", view=view)
    rejects("mrvi", "cuda", view=view)
    available = True
validate("mrvi", "cuda", "hvg2000", view="benchmark_analysis")
validate("mrvi", "cpu", "hvg1000", view="benchmark_analysis")
rejects("mrvi", "cpu", "hvg2000", view="benchmark_analysis")

# The actual scientific method body changes only its accelerator argument.
calls = []
class FakeDistances:
    def __getitem__(self, key): assert key == "dummy_col"; return self
    def isel(self, **kwargs): assert kwargs == {"dummy_col_name": 0}; return self
    def to_pandas(self): return "stored_sample_distances"
class FakeMRVI:
    @staticmethod
    def setup_anndata(adata, **kwargs): calls.append(("setup", kwargs))
    def __init__(self, adata): calls.append(("model",))
    def train(self, **kwargs): calls.append(("train", kwargs))
    def get_local_sample_distances(self, **kwargs): calls.append(("distance", kwargs)); return FakeDistances()
ns.update(MRVI=FakeMRVI, np=SimpleNamespace(zeros=lambda size: [0] * size),
          _ordered_sample_ids=lambda adata: ["s1", "s2"],
          _align_square_frame=lambda value, ids, path: value,
          atomic_to_feather=lambda value, path, **kwargs: calls.append(("write", value, kwargs)))
exec(compile(ast.Module(body=[function("run_mrvi")], type_ignores=[]), str(WORKER), "exec"), ns)
results = []
for device in ("cpu", "cuda"):
    calls.clear()
    ns["run_mrvi"](SimpleNamespace(n_obs=2, obs={}), device, "not-written.feather", batch_key="technical_batch")
    assert calls[0] == ("setup", {"sample_key": "Sample", "batch_key": "technical_batch"})
    assert calls[2] == ("train", {"max_epochs": 50, "accelerator": device})
    results.append([call for call in calls if call[0] != "train"])
assert results[0] == results[1]

main = function("main")
log = next(n for n in main.body if isinstance(n, ast.If) and any(
    isinstance(value, ast.Constant) and isinstance(value.value, str) and value.value.startswith("MRVI_DEVICE")
    for value in ast.walk(n)))
ns["args"] = SimpleNamespace(method="mrvi", device="cuda", view="batch_effect_corrected")
out = io.StringIO()
with contextlib.redirect_stdout(out): exec(compile(ast.Module(body=[log], type_ignores=[]), str(WORKER), "exec"), ns)
assert "MRVI_DEVICE device=cuda name=Mock GPU view=batch_effect_corrected" in out.getvalue()
property_calls = []
cuda.current_device = lambda: property_calls.append("current_device") or 0
cuda.get_device_name = lambda index: property_calls.append("get_device_name") or "Mock GPU"
available = False
out = io.StringIO()
with contextlib.redirect_stdout(out): exec(compile(ast.Module(body=[log], type_ignores=[]), str(WORKER), "exec"), ns)
assert "name=unavailable" in out.getvalue() and property_calls == []
rejects("mrvi", "cuda", view="batch_effect_corrected")
available = True
ns["args"].view = "benchmark_analysis"
out = io.StringIO()
with contextlib.redirect_stdout(out): exec(compile(ast.Module(body=[log], type_ignores=[]), str(WORKER), "exec"), ns)
assert out.getvalue() == "" and property_calls == [], "Benchmark logger initialized CUDA outside fit timing"
print("PASS: actual batch CPU/default/any allocation and VRAM, unchanged benchmark policies, explicit CUDA/no fallback, unchanged technical-only method body and device provenance")
