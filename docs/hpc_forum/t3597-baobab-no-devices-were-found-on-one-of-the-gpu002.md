# Baobab: "No devices were found" on one of the GPU002

- Source: https://hpc-community.unige.ch/t/3597

- Created: 2024-08-19T14:48:21.368Z

- Posts: 4

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2024-08-19T14:48:21.405Z)

Hi HPC,
I allocated some resources on gpu002 and got a gpu. However when running code or nvidia-smi there is no device found (see picture).
Best,
Malte
image
image1888×425 22.8 KB
[image1888×425 22.8 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/53185932a23356cecf6b660f6508e69e661f3de2.png)


## Post 2 by @Raphael.Rubino (2024-08-20T07:16:32.423Z)

Hello,
Yes I think it is GPU#3 (index starts at 0) on baobab:gpu002 which is down.
I created a report a while ago:
Error when running GPU jobs on gpu002 (GPU #3 seems down)[Error when running GPU jobs on gpu002 (GPU #3 seems down)](https://hpc-community.unige.ch//hpc-community.unige.ch/t/error-when-running-gpu-jobs-on-gpu002-gpu-3-seems-down/3537) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Primary informations
Username: rubinor 
Cluster: baobab 
Node: gpu002 
Description
Error when running GPU jobs on gpu002. 
Error messages: 
– Python 
RuntimeError: No CUDA GPUs are available

– When the job starts, quick ssh to gpu002 and run nvidia-smi 
No devices were found

Current usage of gpu002: 
CfgTRES=cpu=10,mem=257000M,billing=10,gres/gpu=6,gres/gpu:titan=6
AllocTRES=cpu=8,mem=94G,gres/gpu=3,gres/gpu:titan=3

Steps to Reproduce
The error only happens if gpu002 is already in use. To rep…


## Post 3 by @Yann.Sagon (2024-08-20T14:13:14.559Z)

Indeed, @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) I’ve answered your previous post right now. And @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren) we’ll check the node, thanks.
Best


## Post 4 by @Yann.Sagon (2024-08-23T12:43:04.039Z)

This is fixed, the node as 6 GPUs again.
Best
