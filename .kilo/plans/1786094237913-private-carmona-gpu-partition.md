# Add `private-carmona-gpu` node to SLURM partitions

## Goal & decisions (confirmed with user)

- **Stages 2–4** (`src/2_dataset_specific_preprocessing`, `src/3_scrnaseq_preprocessing`, `src/4_cell_type_annotation`): default `SLURM_PARTITION="shared-cpu,private-carmona-gpu"` (comma list — jobs may land on either partition, whichever frees resources first; order is not a strict preference).
- **Stage 5 (benchmark)**: private node must NOT be used for real benchmark runs (different CPU model would flaw cross-method runtime comparisons; its GPU ≠ H200 breaks the pinned GPU constraint). It stays usable for `_debug` testing via the existing `--partition <P>` override — which requires dropping the `--constraint` flags on override, otherwise jobs hang PENDING forever (constraint unsatisfiable on the private node).
- All partitions are passed on the sbatch/srun **command line** (no `#SBATCH --partition` directives exist anywhere), so `slurm_config.sh` is the single source of truth.

## Prerequisite (user-run on HPC login node — agents cannot SSH)

Confirm the partition actually exists as named and learn its requirements (account/QoS, GPU gres, CPU model, mem, time limit):

```bash
sinfo -a | grep -iE "private|carmona"
scontrol show partition private-carmona-gpu   # check AllowAccounts / AllowGroups / QOS
sinfo -N -p private-carmona-gpu -o "%n %c %m %G %C %l"
```

PS: I (the user) just quickly checked:
"""
(bamboo)-[halterc@login1 ~]$ sinfo -a | grep -iE "private|carmona"
private-burgi-gpu         up 7-00:00:00      0    n/a 
private-carmona-gpu       up 7-00:00:00      1   mix- gpu004
private-cui-cpu           up 7-00:00:00      4   mix- cpu[049-052]
private-cui-cpu           up 7-00:00:00      2  down* cpu[047-048]
private-drim-gpu          up 7-00:00:00      1   mix- gpu006
private-engelke-gpu       up 7-00:00:00      1   mix- gpu007
private-gapnl-cpu         up 7-00:00:00      3   mix- cpu[042-043,046]
private-hepia-gpu         up 7-00:00:00      1    mix gpu008
private-phylolab-cpu      up 7-00:00:00      1   mix- cpu041
private-kalousis-gpu      up 7-00:00:00      1   plnd gpu009
private-kalousis-gpu      up 7-00:00:00      1 drain* gpu010
private-mhn-cpu           up 7-00:00:00      1   mix- cpu040
(bamboo)-[halterc@login1 ~]$ scontrol show partition private-carmona-gpu   # check AllowAccounts / AllowGroups / QOS
PartitionName=private-carmona-gpu
   AllowGroups=private_carmona AllowAccounts=ALL AllowQos=ALL
   AllocNodes=ALL Default=NO QoS=N/A
   DefaultTime=00:01:00 DisableRootJobs=NO Exclusive=NO GraceTime=0 Hidden=NO
   MaxNodes=UNLIMITED MaxTime=7-00:00:00 MinNodes=0 LLN=NO MaxCPUsPerNode=UNLIMITED MaxCPUsPerSocket=UNLIMITED
   Nodes=gpu004
   PriorityJobFactor=4 PriorityTier=4 RootOnly=NO ReqResv=NO OverSubscribe=NO
   OverTimeLimit=NONE PreemptMode=OFF
   State=UP TotalCPUs=64 TotalNodes=1 SelectTypeParameters=NONE
   JobDefaults=(null)
   DefMemPerNode=UNLIMITED MaxMemPerNode=UNLIMITED
   TRES=cpu=62,mem=750G,node=1,billing=264,gres/gpu=1,gres/gpu:nvidia_h100_nvl=1
   TRESBillingWeights=CPU=1.0,Mem=0.25G,GRES/gpu=1,GRES/gpu:nvidia_a100-pcie-40gb=5,GRES/gpu:nvidia_a100_80gb_pcie=8,GRES/gpu:nvidia_geforce_rtx_2080_ti=2,GRES/gpu:nvidia_geforce_rtx_3080=3,GRES/gpu:nvidia_geforce_rtx_3090=5,GRES/gpu:nvidia_geforce_rtx_4090=8,GRES/gpu:nvidia_geforce_rtx_5090=10,GRES/gpu:nvidia_h100_nvl=14,GRES/gpu:nvidia_h200_nvl=17,GRES/gpu:nvidia_rtx_5000=9,GRES/gpu:nvidia_rtx_a5000=5,GRES/gpu:nvidia_rtx_a5500=5,GRES/gpu:nvidia_rtx_a6000=8,GRES/gpu:nvidia_rtx_pro_4500_blackwell=10,GRES/gpu:nvidia_rtx_pro_6000_blackwell=16,GRES/gpu:nvidia_titan_rtx=1,GRES/gpu:nvidia_titan_x=1,GRES/gpu:tesla_p100-pcie-12gb=1,GRES/gpu:tesla_v100-pcie-32gb=3

(bamboo)-[halterc@login1 ~]$ sinfo -N -p private-carmona-gpu -o "%n %c %m %G %C %l"
HOSTNAMES CPUS MEMORY GRES CPUS(A/I/O/T) TIMELIMIT
gpu004 64 768000 gpu:nvidia_h100_nvl:1(S:2),VramPerGpu:no_consume:94G 4/60/0/64 7-00:00:00
"""
-> It's currently not blocked fur us anymore due to previously low usage, but I will contact our IT department to blokc it again for us during office hours.

If the partition requires an account/QoS, carry out step 4 below as well.

## Steps

### 1. `src/slurm_config.sh`

- Line 84: `SLURM_PARTITION="shared-cpu,private-carmona-gpu"` — update the comment: comma list = jobs can land on either partition (used by stages 2–4 submit scripts); benchmark uses its own pinned vars below.
- Add after the benchmark block: `SLURM_PARTITION_PRIVATE="private-carmona-gpu"` (documented constant for `_debug` benchmark runs via `--partition ${SLURM_PARTITION_PRIVATE}`).
- Update the benchmark comment block (lines 66–79): note the private node is deliberately excluded from `SLURM_PARTITION_BENCHMARK_*` (CPU model differs → flawed runtime comparisons; GPU ≠ H200 → constraint mismatch). Keep the benchmark var values unchanged.
- Keep partition vars unexported plain assignments (existing style).

### 2. `src/5_run_benchmark_methods/run_python_sample_embedding_methods/1_submit_hpc_array.sh`

- When `--partition <P>` is passed (line 137–139), **drop the method's `--constraint` flag** from `EXTRA_FLAGS` (keep `--gpus`, `--cpus-per-task`, `--mem`): an explicit partition choice means the user accepts non-pinned hardware. This also fixes the pre-existing `--partition debug-cpu` usage (constraint `EPYC-7742` would otherwise never match).
- Update the header comment (lines 10–23) to document this.
- Print a notice when the override drops the constraint, e.g. `echo "  NOTE: --partition override drops the --constraint hardware pin"`.

### 3. Docs

- `docs/ARCHITECTURE.md`:
  - Benchmark workflow diagram (line 207) + "Hardware pinned for runtime comparability" bullet (lines 236–240): add that an explicit `--partition <P>` override drops the constraints (debug/private nodes; e.g. `_debug` smoke tests on the private node).
  - Benchmark submitter row (line 229): update the sentence "`--partition <P>` overrides the per-method partition (debug/private partitions)" to also mention the constraint drop.
  - `--device auto` bullet (line 258): note it also uses the GPU on the private node.
  - Usage block (line 275): add `./1_submit_hpc_array.sh --ds_name _debug --methods pilot --partition private-carmona-gpu`.
- Stage-2 dispatcher row (line 78) keeps its generic wording; no change needed beyond optionally noting the comma list.

### 4. Conditional: account/QoS required

Only if the prerequisite `scontrol show partition` output shows `AllowAccounts`/`AllowGroups`/`QOS` restrictions:

- Add `SLURM_ACCOUNT` (and `SLURM_QOS` if required) to `slurm_config.sh`.
- Pass `--account="${SLURM_ACCOUNT}"` (+ `--qos=...` if needed) on the sbatch/srun command lines in: `src/2_dataset_specific_preprocessing/1_submit_hpc.sh`, `src/3_scrnaseq_preprocessing/1_submit_hpc_array.sh`, `src/4_cell_type_annotation/1_prepare_chunks.sh`, `src/4_cell_type_annotation/2_submit_hpc_array.sh`, `src/4_cell_type_annotation/3_submit_merge.sh`.
- Add a warning if the var is unset at submit time.

## Validation (user-run on HPC login node)

1. Smoke test: `sbatch --partition="shared-cpu,private-carmona-gpu" --wrap "echo ok; scontrol show job \$SLURM_JOB_ID"` → confirm the job lands on a valid node and `sacct -j <id> -o Partition` shows one of the two partitions.
2. Preprocess: `./1_submit_hpc_array.sh --ds_name _debug` → tasks COMPLETED, NAS sync ok.
3. Annotation spot-check: `./1_prepare_chunks.sh test <DS>` and/or `2_submit_hpc_array.sh --ds_name _debug`.
4. Benchmark debug on private node: `./1_submit_hpc_array.sh --ds_name _debug --methods pilot --partition private-carmona-gpu` → job RUNNING (not PENDING on constraint), feathers + exec log produced.
5. Default benchmark unchanged: `./1_submit_hpc_array.sh --methods mrvi,scpoli,pilot` still pins H200 (`nvidia_h200_nvl`) / EPYC-7742 — verify via `scontrol show job` on a sample task.

## Risks / notes

- Comma-list semantics: jobs are not guaranteed to use the private node — whichever partition has free resources first wins. That is the requested behavior; if full private-node exclusivity is ever wanted, it requires a code edit (or optionally making `SLURM_PARTITION` env-overridable via `SLURM_PARTITION="${SLURM_PARTITION:-...}"` — out of scope for this change unless requested).
- Private-node workers get `PYTHON_BIN`/pixi env identically (sourced from `slurm_config.sh`); no interpreter changes. No NAS access on workers, same as other partitions — unaffected.
- Benchmark worker `#SBATCH --time=12:00:00` is the `shared-*` max; if the private partition has a lower time limit (check `%l` in the prerequisite), `_debug` GPU runs may need a time override — note only, do not change without user confirmation.
- `_debug` GPU runs on the private node use `--device auto` (MrVI) → uses whatever GPU the node has; acceptable for debugging only.
