# Phase 3.1 — Convert Python benchmark methods to HPC pipeline

Implements TODO.md:42-46 (Phase 3.1): rename + convert
`src/5_run_benchmark_methods/run_python_sample_embedding_methods/1.2_benchmark_methods_py.qmd`
into a CLI Python script, and add the SLURM array submit + worker scripts for
MrVI / scPoli / PILOT.

## Decisions (user-confirmed)

- **GPU methods**: MrVI + scPoli on `shared-gpu`, `--gpus=1`,
  `--constraint=nvidia_h200_nvl` (H200, gpu[005-006]; gpu006 has 4 GPUs),
  `--cpus-per-task=8`, `--mem=128G`. MrVI now runs in **pytorch** — drop the
  `jax` import from the qmd.
- **CPU method**: PILOT on `shared-cpu`, `--constraint=EPYC-7742` (V8,
  cpu[001-052]), `--cpus-per-task=16`, `--mem=128G`.
- Hardware is **pinned for runtime comparability**: same GPU model for all GPU
  tasks, same CPU model + core count + RAM for all CPU tasks.
- **Memory logging**: `execution_times.feather` gains a `mem_GB` column
  (peak RSS via `resource.getrusage().ru_maxrss`; Linux = KB, macOS = bytes —
  handle platform; convert to GB). Columns: `dataset`, `method`, `time_secs`
  (float seconds = R `as.numeric(difftime, units="secs")` equivalent),
  `mem_GB`. The notebook's `rbind(r_exec_times, p_exec_times)` tolerates the
  extra column (NA for R rows; only `time_secs` is plotted).
- **qmd kept as archive** — do NOT delete `1.2_benchmark_methods_py.qmd`.
- Feather naming, method-string format, and data semantics of the qmd are
  **preserved exactly** (R ingest `run_benchmark_analysis` file checks +
  `process_mrvi_fig`/`process_scpoli_fig`/`process_pilot_fig` +
  `constants.R` label map + notebook recodes depend on them).

## Files

Create (all in `src/5_run_benchmark_methods/run_python_sample_embedding_methods/`,
chmod +x the `.sh`):

| File | Role |
|---|---|
| `1.1.1_benchmark_methods_py.py` | CLI py script (replaces qmd logic) |
| `1.1.2_merge_execution_times.py` | concat per-task exec logs → `execution_times.feather` (run on login node after arrays) |
| `1.1_run_worker.sh` | `#SBATCH` worker: one (method, dataset) task |
| `1_submit_hpc_array.sh` | login-node submitter: one array per method + monitor + sacct gate + merge + NAS sync |

Modify:

- `src/slurm_config.sh` — benchmark hardware/config vars (below).
- `docs/ARCHITECTURE.md` — document the 4 new files + `benchmark/embeddings/`
  data flow (the `benchmark/{embeddings,plots}/` NAS target already exists at
  line 49).
- `TODO.md` — note 3.1 code complete, HPC debug validation pending.

## 1. `1.1.1_benchmark_methods_py.py`

`sys.path.insert(0, parents[2])` and reuse `src/datasets_io.read_datasets_json`
(import from repo root like `1.1.1_preprocess.py` does).

### CLI

```
--config_path <datasets.json>      # required
--ds_name <DS>                     # required
--view benchmark_analysis          # default benchmark_analysis
--method {mrvi,scpoli,pilot}       # required
--input_dir <dir>                  # dir holding the preprocessed view h5ad
--output_dir <dir>                 # feather output dir (created if missing)
--log_file <path>                  # per-task exec log feather (worker passes
                                   # execution_times_task_<TASKID>.feather;
                                   # default: <output_dir>/execution_times.feather
                                   # for local runs)
--hvg 1000 2000 3000               # default 1000 2000 3000 (nargs="+", int)
--force                            # recompute existing outputs
--device auto|cpu|cuda             # default auto (train accelerator; auto uses
                                   # the allocated GPU on shared-gpu nodes)
```

### Input resolution

- `entry = read_datasets_json(config_path, view=view)[ds_name]`; read
  `input_dir / entry["views"][view]["output_file"]` via `sc.read_h5ad(...)`,
  `.to_memory()` (mirrors qmd).
- CT annotation columns from datasets.json: lowres = `entry["low_res_ct_col"]`
  (MrVI ignores it; PILOT/scPoli use it), highres = `entry["hi_res_ct_col"]`.
- Sample column: `obs["Sample"]` (standardized by preprocess; `SAMPLE_COLNAME`).
- **HVG subset**: `var["hvg_rank"]` (stored by `1.1.1_preprocess.py`
  `select_hvgs_ranked`), sort ranks, take top-n.

### Method combos + outputs (legacy naming, ds = datasets.json key)

| Method | Combos | Output file |
|---|---|---|
| mrvi | lowres ct col; every n in `--hvg` | `{ds}_hvg{n}_mrvi_dists.feather` |
| scpoli | highres: n=2000 → dims [2,3,5,10,15]; n=1000/3000 → dims [15]; lowres: n=2000 → dims [15] | `{ds}_hvg{n}{_lowres|_highres}_scpoli_dims{d}_embs.feather` |
| pilot | highres: every n; lowres: n=2000 only | `{ds}_hvg{n}{_lowres|_highres}_pilot_dists.feather` |

Skip a combo if its output file exists and `--force` is not set (qmd skip
logic). Load the h5ad once per task, not per combo.

### Method bodies (preserve qmd semantics)

- **MrVI** (lowres only): `adata.obs["dummy_col"] = 0`;
  `MRVI.setup_anndata(adata, sample_key="Sample")`;
  `model = MRVI(adata); model.train(max_epochs=50, accelerator=args.device)`;
  `dists = model.get_local_sample_distances(keep_cell=False,
  groupby="dummy_col", batch_size=32)`; `dists["dummy_col"].isel(
  dummy_col_name=0).to_pandas().to_feather(...)`. (qmd used
  `accelerator="cpu"` — now `--device` default `auto` so GPU is used on
  shared-gpu nodes.)
- **scPoli** (dims loop): after HVG subset, `X = X.todense().astype("float32")`
  (dense required by scPoli, qmd does the same); per dim:
  `scPoli(adata, condition_keys="Sample", cell_type_keys=<ct col>,
  embedding_dims=dim, recon_loss="nb")`;
  `train(n_epochs=50, pretraining_epochs=40, eta=5,
  early_stopping_kwargs={<qmd values>})`;
  `get_conditional_embeddings()` → `pd.DataFrame(adata_emb.X,
  index=adata_emb.obs_names)` → `to_feather` (index preserved → last column =
  sample names, matching `process_scpoli_fig`'s `column_to_rownames(ncol)`).
- **PILOT**: consume the **preprocessed obsm** (this is the TODO's
  "consume preprocessed obsm" requirement — qmd recomputed PCA from scratch):
  `pl.tl.wasserstein_distance(adata, emb_matrix=f"X_pca_{view}_hvg{n}",
  clusters_col=<ct col>, sample_col="Sample", status="Sample")`;
  `adata.uns["EMD_df"].to_feather(...)`.

Feather writing: plain `DataFrame.to_feather(path)` with the pandas index
(sample names) — do NOT `reset_index()`. This reproduces the qmd layout the R
ingest functions expect.

### Exec-time/memory logging

- Per combo: `time.time()` around the method body only (excludes h5ad
  loading, like the qmd and R `exec_time()`), plus `ru_maxrss` after the run.
- Method strings (exact legacy format — `constants.R` + notebook recodes
  depend on them): `MrVI_hvg{n}`, `scPoli_hvg{n}_dims{d}{res}`,
  `PILOT_hvg{n}{res}`.
- Row = `{dataset: <ds key>, method: <method str>, time_secs: float,
  mem_GB: float}`. Append/read-modify-write into the per-task `--log_file`
  (single process per task → no concurrency issue). `scvi.settings.seed = 0`.
- Drop the `jax` import; keep `torch` (MrVI/scPoli pytorch). Keep
  `scvi.settings.seed = 0` and print `scvi.__version__` + `torch.cuda.is_available()`.

## 2. `1.1.2_merge_execution_times.py`

- CLI: `--output_dir <dir>` (default `benchmark/embeddings`), `--cleanup`
  (default on: delete per-task logs after merge).
- Read `execution_times_task_*.feather` from `output_dir`, concat,
  `drop_duplicates(subset=["dataset","method"], keep="last")` (matches the
  qmd's overwrite-on-rerun semantics), write
  `output_dir/execution_times.feather`.

## 3. `1.1_run_worker.sh`

`#SBATCH` defaults (overridden at submit time): `--job-name=benchmark_worker`,
`--time=12:00:00` (shared-* partition max), `--nodes=1`, `--ntasks=1`,
`--cpus-per-task=8`, `--mem=128G`, `--mail-type=END,FAIL`.

Body (follows `2.1_run_worker.sh` conventions):
- `set -euo pipefail`; `SCRIPT_DIR`; `source "${SCRIPT_DIR}/../slurm_config.sh"`;
  `cd "${PROJECT_ROOT}"`.
- `METHOD` from env (exported by submit script; error if unset).
- Read `DS_NAME` from `${BENCHMARK_MANIFEST}` via `sed -n
  "${SLURM_ARRAY_TASK_ID}p"` (submit script builds one manifest per method:
  `${HPC_SCRATCH_DIR}/benchmark_manifest_<method>.txt`, one DS per line).
- `OUT_DIR="${HPC_SCRATCH_DIR}/benchmark/embeddings"`;
  `mkdir -p "${OUT_DIR}"`;
  `LOG_FILE="${OUT_DIR}/execution_times_task_${SLURM_ARRAY_TASK_ID}.feather"`.
- Forward `--force` when `FORCE_BENCHMARK=1` (env from submit script).
- Run:
  `"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.1_benchmark_methods_py.py"
  --config_path "${DATASETS_JSON_FILE}" --ds_name "${DS_NAME}"
  --view benchmark_analysis --method "${METHOD}"
  --input_dir "${HPC_SCRATCH_DIR}/${DS_NAME}/output"
  --output_dir "${OUT_DIR}" --log_file "${LOG_FILE}" ${FORCE_FLAG}`
- Echo progress lines mirroring `2.1_run_worker.sh` style.

## 4. `1_submit_hpc_array.sh`

Usage: `./1_submit_hpc_array.sh [--ds_name <DS>] [--methods mrvi,scpoli,pilot]
[--partition <P>] [--force]`. Sources `slurm_config.sh`, `cd "${PROJECT_ROOT}"`,
`module load jq/1.6`.

1. Resolve datasets: all keys with `use_for_benchmark == true` AND a
   `benchmark_analysis` view (jq filter); skip `_*` keys unless
   `--ds_name` explicitly given (preprocess convention). `--ds_name`
   validates against datasets.json.
2. For each method in `--methods` (default all three), resolve resources:
   - mrvi, scpoli → partition `${SLURM_PARTITION_BENCHMARK_GPU}`,
     `--gpus="${BENCHMARK_GPU_COUNT}"`,
     `--constraint="${BENCHMARK_GPU_CONSTRAINT}"`,
     `--cpus-per-task="${BENCHMARK_GPU_CPUS_PER_TASK}"`
   - pilot → partition `${SLURM_PARTITION_BENCHMARK_CPU}`,
     `--constraint="${BENCHMARK_CPU_CONSTRAINT}"`,
     `--cpus-per-task="${BENCHMARK_CPU_CPUS_PER_TASK}"`
   - all: `--mem="${BENCHMARK_MEM}"`. `--partition <P>` overrides the
     per-method partition (for debug/private partitions).
   - Flags are passed on the **sbatch command line** (SLURM directives do not
     expand env vars — established convention).
3. Write `${HPC_SCRATCH_DIR}/benchmark_manifest_<method>.txt` (one DS per
   line, rebuilt every run); `export BENCHMARK_MANIFEST` per submission;
   `export FORCE_BENCHMARK=1` when `--force`.
4. `sbatch --array="1-${N}%${THROTTLE}"` — throttle: `BENCHMARK_GPU_ARRAY_THROTTLE`
   (4, matching the 4 H200s on gpu006) for GPU methods, `${MAX_NUM_CHUNKS_PARALLEL}`
   for CPU. Logs: `${LOGS_DIR}/5_benchmark_<method>_%A_%a.log/.err`.
   Collect each array job id.
5. Monitor: `while squeue -u "$USER" | grep -q <id>` per job (sleep 60),
   `sleep 30` after; then **fail-closed sacct gate** per job id (mirror
   `1_submit_hpc_array.sh`): every state row must be `COMPLETED`, else print
   `sacct --format=JobID,JobName,State,ExitCode` and exit 1 without syncing.
6. Merge exec logs: `"${PYTHON_BIN}" "${SCRIPT_DIR}/1.1.2_merge_execution_times.py"
   --output_dir "${HPC_SCRATCH_DIR}/benchmark/embeddings"` (login node).
7. NAS sync (login node; NAS unreachable → error exit): verify
   `ls "${NAS_TARGET_DIR}/.."`, `mkdir -p "${NAS_TARGET_DIR}/benchmark"`,
   `rsync -rlptDv "${HPC_SCRATCH_DIR}/benchmark/" "${NAS_TARGET_DIR}/benchmark/"`.

## 5. `slurm_config.sh` additions

```bash
# --- Benchmark methods (src/5_run_benchmark_methods) ---
# Hardware PINNED for runtime comparability (same model/cores/RAM within each
# resource class). GPU methods (MrVI, scPoli): shared-gpu, H200
# (nvidia_h200_nvl, gpu[005-006]), 8 cores, 128G. CPU method (PILOT):
# shared-cpu, EPYC-7742 (V8, cpu[001-052]), 16 cores, 128G. All env-overridable.
SLURM_PARTITION_BENCHMARK_GPU="shared-gpu"
SLURM_PARTITION_BENCHMARK_CPU="shared-cpu"
BENCHMARK_GPU_CONSTRAINT="nvidia_h200_nvl"
BENCHMARK_CPU_CONSTRAINT="EPYC-7742"
BENCHMARK_GPU_COUNT=1
BENCHMARK_GPU_CPUS_PER_TASK=8
BENCHMARK_CPU_CPUS_PER_TASK=16
BENCHMARK_MEM="128G"
BENCHMARK_GPU_ARRAY_THROTTLE=4   # 4 H200s on gpu006
```

## Data flow

```
1_submit_hpc_array.sh (login)
  └─ per method: manifest + sbatch array (partition/gpus/constraint/cpus/mem on CLI)
       └─ 1.1_run_worker.sh (worker, METHOD + manifest → DS)
            └─ 1.1.1_benchmark_methods_py.py
                 input:  ${HPC_SCRATCH_DIR}/<DS>/output/<view output_file> h5ad
                         (obsm X_pca_{view}_hvg{n} for PILOT, layers/var hvg_rank
                          for MrVI/scPoli; Sample obs col; ct annotation cols)
                 output: ${HPC_SCRATCH_DIR}/benchmark/embeddings/
                         {ds}_hvg{n}[_lowres|_highres]_{mrvi_dists,scpoli_dims<d>_embs,pilot_dists}.feather
                         execution_times_task_<TASKID>.feather
  └─ sacct gate → 1.1.2_merge_execution_times.py → execution_times.feather
  └─ rsync → ${NAS_TARGET_DIR}/benchmark/embeddings/ + execution_times.feather
```

## Validation

- `bash -n` on both `.sh`; `python -m py_compile` on both `.py`; ensure the
  new `.sh` are committed with the executable bit (repo convention `100755`).
- No HPC run in this task (AGENTS.md: validation happens on the HPC debug
  dataset after the full pipeline is implemented, user-run). Submit-script
  smoke test is part of the user's Phase 2/3 debug run:
  `./1_submit_hpc_array.sh --ds_name _debug --methods mrvi` (then scpoli, pilot)
  against the preprocessed `_debug` benchmark h5ad; check feathers + exec log;
  R-ingest compatibility is exercised by Phase 3.4.

## Out of scope (later TODO phases)

- 3.2 R methods on HPC; 3.3 new methods (PILOT-GM-VAE, QOT, PULSAR); 3.4
  notebook adaptation; 3.5 full partition strategy (this plan pins the 3.1
  hardware in slurm_config.sh, which 3.5 will formalize); 3.6 README/AGENTS
  docs; Phase 4 batch-effect view support (the CLI `--view` arg already
  supports it; feather naming for batch-effect views is a Phase 4 concern).
