# Srun error when runing job

- Source: https://hpc-community.unige.ch/t/3727

- Created: 2024-11-13T10:07:30.657Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-11-13T10:07:30.690Z)

## Primary informations
Username: mongin
Cluster:Baobab
## Description
When running a job with srun, I have an error:
```
srun: error: Couldn't find the specified plugin name for mpi/pmix_v3 looking at all files
srun: error: cannot find mpi plugin for mpi/pmix_v3
srun: error: MPI: Cannot create context for mpi/pmix_v3
srun: error: MPI: Unable to load any plugin
srun: error: Invalid MPI type 'pmix_v3', --mpi=list for acceptable types
```
## Steps to Reproduce
I am making `sbatch baobab_classify_SR.bash` in `[mongin@login1 classify_SR]`.
The batch file load the librarioes and call a pathon virtual env. I am able to run each part manually myself, so the problem is with `srun`, can’t figure out why.
The batch file:
```
#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --gpus=2
#SBATCH --partition=shared-gpu
#SBATCH --gres=VramPerGpu:25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mem=30000
#SBATCH --array=1,5,13

. ~/baobab_python_env_LLM3/bin/activate
ml  GCC/12.3.0  OpenMPI/4.1.5 PyTorch-bundle/2.1.2-CUDA-12.1.1

srun ~/baobab_python_env_LLM3/bin/python -u classify_SR.py ${SLURM_ARRAY_TASK_ID} > ./results/classify.out
```
I do : `sbatch baobab_classify_SR.bash`
## Expected Result
The file should launch the jobs, was working before (last week). Here it stops at launch, and I have the following errors in the slurm files:
```
srun: error: Couldn't find the specified plugin name for mpi/pmix_v3 looking at all files
srun: error: cannot find mpi plugin for mpi/pmix_v3
srun: error: MPI: Cannot create context for mpi/pmix_v3
srun: error: MPI: Unable to load any plugin
srun: error: Invalid MPI type 'pmix_v3', --mpi=list for acceptable types
```


## Post 2 by @Yann.Sagon (2024-11-13T11:05:42.361Z)

Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin) sorry about that.
A  quick working workaround is to force srun in your sbatch to use pmi2.
```
srun --mpi=pmi2 ~/baobab_python_env_LLM3/bin/python -u classify_SR.py ${SLURM_ARRAY_TASK_ID} > ./results/classify.out
```
We’ll investigate why it tries to use pmi3 and update this post.


## Post 3 by @Denis.Mongin (2024-11-13T12:18:25.102Z)

Perfect, this does the trick.
Thank you again for your reactivity.
Denis


## Post 4 by @Yann.Sagon (2024-11-14T09:34:39.845Z)

@Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
you can remove the workaround, this is fixed.
