# Multi-gpu performance and number of tasks per node

- Source: https://hpc-community.unige.ch/t/3338

- Created: 2024-02-27T11:53:09.712Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Francesco.Marson (2024-02-27T11:53:09.767Z)

Hello HPC support,
## Primary information
Username: marson
Cluster: baobab
## Description
I am experimenting with the performance of palabos on gpu, and I noticed that the only way to get the expected performance from the code is to force having one mpi task per node and 1 gpu per task (`--ntasks=2`,`--ntasks-per-node=1`,`--exclusive`,`--gres=gpu:1`), while trying to use multiple gpus per node (`--ntasks=2`,`--ntasks-per-node=2`,`--exclusive`,`--gres=gpu:1`) leads to very bad performances (80MLUPS vs 8000MLUPS).
What am I doing wrong or misinterpreting in the `sbatch` options?
Thank you in advance for the help and support!
Cheers,
Francesco


## Post 2 by @Yann.Sagon (2024-02-27T13:10:03.656Z)

Hi @Francesco.Marson[@Francesco.Marson](https://hpc-community.unige.ch/u/francesco.marson)
how many GPUs do you expect to have per task?
In this variant:
```
srun --ntasks=2 --ntasks-per-node=1 --exclusive --gres=gpu:1
```
you have two task, one on each node, and one gpu on each node: it means each tasks (2) has a GPU. This works because you only have one task per node, but this is not a good solution  as this prevent other users to use the remaining resources of the compute node. Do not use the `--exclusive` flag unless necessary.
In this variant:
```
srun --ntasks=2 --ntasks-per-node=2 --exclusive --gres=gpu:1
```
You’re asking for two tasks on only one node and thus you get only one GPU. The two tasks have to share the single GPU, probably not what you want.
Tu summarize you thought you were using a single GPU in the first variant and multiple GPU on the second, but this is the inverse.
