# Multijob GPU partly failing: can't run 1 processes on 2 nodes

- Source: https://hpc-community.unige.ch/t/3728

- Created: 2024-11-13T12:44:37.915Z

- Posts: 8

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Denis.Mongin (2024-11-13T12:44:37.949Z)

## Primary informations
Username: mongin
Cluster: baobab
## Description
When lauching a multijob, part of the jobs are failing at start with the message:
```
srun: warning: can't run 1 processes on 2 nodes, setting nnodes to 1
srun: error: Unable to create step for job 13610226: Invalid generic resource (gres) specification
```
## Steps to Reproduce
I am lauching `sbatch baobab_classify_SR.bash`, which is the following:
```
#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --gpus=2
#SBATCH --partition=shared-gpu
#SBATCH --gres=VramPerGpu:25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mem=30000
#SBATCH --array=1-49

. ~/baobab_python_env_LLM3/bin/activate
ml  GCC/12.3.0  OpenMPI/4.1.5 PyTorch-bundle/2.1.2-CUDA-12.1.1

srun --mpi=pmi2 ~/baobab_python_env_LLM3/bin/python -u classify_SR.py ${SLURM_ARRAY_TASK_ID} > ./results/classify.out
```
Jobs 1 and 2 are running, but 3 to 11 failed. 12 worked, other are still pending.
The jobs failing report in the slurm files:
```
(baobab)-[mongin@login1 classify_SR]$ cat slurm-13609953_3.out
srun: warning: can't run 1 processes on 2 nodes, setting nnodes to 1
srun: error: Unable to create step for job 13610168: Invalid generic resource (gres) specification
```
## Expected Result
I would have expected all jobs to run the same way. I did run 7 jobs on this script without problem.
What am I doing wrong ?
Thank you for your help!


## Post 2 by @Adrien.Albert (2024-11-13T19:22:07.308Z)

Hi Denis
Denis.Mongin:
```
srun: warning: can't run 1 processes on 2 nodes, setting nnodes to 1
```
I’m not sure if this will solve the problem, but following the error message Could you try increasing the ntask to 2 :
```
#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --gpus=2
#SBATCH --partition=shared-gpu
#SBATCH --gres=VramPerGpu:25G
#SBATCH --ntasks=2       # <============= INCREASE HERE
#SBATCH --cpus-per-task 1
#SBATCH --mem=30000
#SBATCH --array=1-49
```
Best Regards,


## Post 3 by @Denis.Mongin (2024-11-14T07:36:01.131Z)

Hi Adrien
Ok, thanks, I will try.
Not sure why this behavior changed, this is the first time I have these error messages


## Post 4 by @Yann.Sagon (2024-11-14T10:57:30.782Z)

Dear @Denis.Mongin[@Denis.Mongin](https://hpc-community.unige.ch/u/denis.mongin)
Here I am with some news for your case. I checked your past jobs that failed and it seems that most of them were initially allocated on two compute nodes as you are requesting two GPU in total but not two gpu par task and converted to a single node job at run time as you had only one task, so not possible to spread that on more than one node.
The flat `--gres` will request a generic consumable, in your case GPUs, that will be allocated in each node.
In summary, you want two GPUs and one task. You should also request two CPUs per task, as it is better to have at least one CPU per GPU.
Using more than one task probably won’t work unless you set things up as explained here[here](https://pytorch.org/tutorials/beginner/dist_overview.html).
Find below an extract of what you should use:
```
#SBATCH --gpus-per-task=2
#SBATCH --partition=shared gpu
#SBATCH --gres=VramPerGpu:25G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G # Make sure you request at least as much memory as you get from the GPU card.
```
Best


## Post 5 by @Denis.Mongin (2024-11-14T11:44:40.265Z)

Ok, thank for the explanation @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon).
I will try and make some tests.
I had the experience that having two cpu on different node created problem of memory allocation when trying to allocate the VRAM memory on two GPU for the same model.
I will see if it the case with the config you propose. Weird that I did not come ccross this problem before though.
Thank you again


## Post 6 by @Yann.Sagon (2024-11-14T13:27:29.548Z)

Denis.Mongin:
> I had the experience that having two cpu on different node created problem of memory allocation when trying to allocate the VRAM memory on two GPU for the same model.
When you request one task with two CPUs, the resources are allocated on only node, this what you want.


## Post 7 by @Denis.Mongin (2024-11-14T13:34:40.062Z)

perfect then, it should work.
Thanks a bunch


## Post 8 by @Denis.Mongin (2024-11-14T14:43:31.296Z)

Works like a charm.
Thank you again for all the support and help.
Denis
