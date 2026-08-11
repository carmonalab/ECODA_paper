# Reducing wasted resources

- Source: https://hpc-community.unige.ch/t/3507

- Created: 2024-06-21T14:22:41.278Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Raphael.Rubino (2024-06-21T14:22:41.314Z)

It appears that a lot of GPU resources are wasted due to jobs using many CPUs for a single GPU. For instance on Yggdrasil, we can currently observe this for node gpu002:
```
CfgTRES=cpu=14,mem=385499M,billing=14,gres/gpu=8,gres/gpu:turing=8
AllocTRES=cpu=14,mem=150G,gres/gpu=3,gres/gpu:turing=3
```
which means that 5 GPUs are idle and the node is seen as fully allocated by Slurm.
Would it be possible to automatically route these jobs (more than 1 CPU per GPU) on nodes with more CPUs and allocate jobs with 1 GPU and 1 CPU to nodes with fewer CPUs?


## Post 2 by @Yann.Sagon (2024-07-05T06:55:57.897Z)

Dear @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) thanks for contacting us, we are aware of this issue and we already tried to tackle it, but unfortunately there is no such easy solution.
The good news is that all the newer GPU server have at least 64 CPUs so this is becoming less an issue.
What we wanted to apply is to enforce that you can’t request more than 2 CPUs per GPUs, but this isn’t in place yet.
What is already in place is a mechanism to prevent to use a GPU node without requesting for a GPU:
```
(yggdrasil)-[sagon@login1 ~]$ srun  --partition=debug-gpu hostname
srun: error: You are trying to submit on gpu partition without requesting gpu, do you really need to use a gpu node ?
srun: error: Unable to allocate resources: Invalid generic resource (gres) specification
```
