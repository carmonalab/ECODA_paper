# Error when running GPU jobs on gpu002 (GPU #3 seems down)

- Source: https://hpc-community.unige.ch/t/3537

- Created: 2024-07-12T16:52:42.416Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2024-07-12T16:52:42.465Z)

## Primary informations
Username: rubinor
Cluster: baobab
Node: gpu002
## Description
Error when running GPU jobs on gpu002.
Error messages:
– Python
```
RuntimeError: No CUDA GPUs are available
```
– When the job starts, quick ssh to gpu002 and run nvidia-smi
```
No devices were found
```
Current usage of gpu002:
```
CfgTRES=cpu=10,mem=257000M,billing=10,gres/gpu=6,gres/gpu:titan=6
AllocTRES=cpu=8,mem=94G,gres/gpu=3,gres/gpu:titan=3
```
## Steps to Reproduce
The error only happens if gpu002 is already in use. To reproduce the error, run the command below n times until all GPUs are allocated.
Run some command on gpu002 requesting 1 GPU and 1 CPU, for instance:
```
sbatch run.sh
```
with “run.sh” containing:
```
#!/bin/bash
#SBATCH --partition=shared-gpu
#SBATCH --job-name=jobname
#SBATCH --time=00:10:00
#SBATCH --error=slurm.e%j 
#SBATCH --output=slurm.o%j
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu002

nvidia-smi
```
## Expected Result
When running the above steps, you should get an output file containing the usual “nvidia-smi” output
```
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 555.42.02              Driver Version: 555.42.02      CUDA Version: 12.5     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA TITAN X (Pascal)        On  |   00000000:83:00.0 Off |                  N/A |
| 23%   30C    P8              8W /  250W |       2MiB /  12288MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
```
## Actual Result
The content of the slurm output file is:
```
No devices were found
```
Edit it seems that GPU 3 is down on gpu002 (GPU ids start at 0)


## Post 2 by @Yann.Sagon (2024-08-20T14:12:18.885Z)

Dear @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) many thanks for the notification, it is a pity we didn’t saw it earlier, you know, vacations etc. I’ve drained the node and we’ll fix that.
Best
Yann
