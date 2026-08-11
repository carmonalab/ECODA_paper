# Important: new GPU types naming in Baobab

- Source: https://hpc-community.unige.ch/t/3772

- Created: 2024-12-19T13:02:10.121Z

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-12-19T13:02:10.205Z)

Dear users,
TLDR; we have changed the GPU type name, please update your script according to the new GRES value in the table[table](https://doc.eresearch.unige.ch/hpc/hpc_clusters#compute_nodes).
Longer explanation:
We have enabled a feature regarding GPUs on Baobab Yggdrasil and Bamboo clusters.
- Previously, we had a static definition of each GPU card on GPU servers. This was basic and Slurm didn’t know the affinity between the CPU socket and the GPU card. Nor did it know the affinities between GPUs. We had a limited “type” name to target GPUS.
- Now we’ve enabled NVML[NVML](https://slurm.schedmd.com/gres.conf.html) GPU auto-detection for Slurm. Slurm is now aware of the GPU topology and can make better resource allocations.
- The GPU type is also automatically inferred and we have a unique type name for each GPU model. You can look up the new GRES name for the requested GPU model in the table. As before, it’s not possible to request more than one GPU type. However, you can use a constraint to filter the GPUs you want to request.
This change should improve job performance as GPUs are now paired with the nearest CPU. This change was made because a user asked us [how to target two GPUs that are linked together] ([GPU][SLURM] How to request a pair of GPUs connected with an NVLINK?[[GPU][SLURM] How to request a pair of GPUs connected with an NVLINK?](https://hpc-community.unige.ch/t/gpu-slurm-how-to-request-a-pair-of-gpus-connected-with-an-nvlink/3741)).
Best regards


## Post 2 by @Malte.Algren (2025-01-14T12:00:46.615Z)

Hi, thanks for the update!
Before the new naming convention, we could request only ampere GPUs like:
` salloc ... --gres=gpu:ampere:1,VramPerGpu:2G`
What is the new convention if I want to restrict the architecture?
Cheers,
Malte


## Post 3 by @Yann.Sagon (2025-01-16T10:08:45.951Z)

Dear  @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren)
you can request a specific constraint instead of using the gpu type:
```
(bamboo)-[sagon@login1 ~]$ srun --constraint=COMPUTE_TYPE_AMPERE --partition=debug-gpu --gres=gpu:1,VramPerGpu:2G hostname
gpu001.bamboo
```
You can see which feature is available on the nodes:
```
(bamboo)-[sagon@login1 ~]$ sinfo -o "Node: %n | Gres: %Gres | Feature: %f" -p shared-gpu
Node: HOSTNAMES | Gres: GRESres | Feature: AVAIL_FEATURES
Node: gpu001 | Gres: gpu:nvidia_geforce_rtx_3090:8(S:0-1),VramPerGpu:no_consume:24Gres | Feature: EPYC-7742,V8,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPERE,SIMPLE_PRECISION_GPU,COMPUTE_MODEL_RTX_3090_25G
Node: gpu002 | Gres: gpu:nvidia_geforce_rtx_3090:8(S:0-1),VramPerGpu:no_consume:24Gres | Feature: EPYC-7742,V8,COMPUTE_CAPABILITY_8_6,COMPUTE_TYPE_AMPERE,SIMPLE_PRECISION_GPU,COMPUTE_MODEL_RTX_3090_25G
Node: gpu003 | Gres: gpu:nvidia_a100_80gb_pcie:4(S:0),VramPerGpu:no_consume:80Gres | Feature: EPYC-7302P,V8,DOUBLE_PRECISION_GPU,COMPUTE_CAPABILITY_8_0,COMPUTE_TYPE_AMPERE,COMPUTE_MODEL_A100_80G
```
Best
