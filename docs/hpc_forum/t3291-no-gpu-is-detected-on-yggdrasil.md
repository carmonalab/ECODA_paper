# No GPU is detected on yggdrasil

- Source: https://hpc-community.unige.ch/t/3291

- Created: 2024-02-06T14:53:53.112Z

- Posts: 4

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-02-06T14:53:53.182Z)

Hi,
My GROMACS jobs on yggdrasil all failed last night. I resubmitted jobs this afternoon. One runs normally with gpu007, but another one only ran for a few sec on gpu008 and failed.
# Some lines of .out file:
```
Program:     gmx mdrun, version 2023.1
Source file: src/gromacs/taskassignment/findallgputasks.cpp (line 85)
MPI rank:    0 (out of 8)

Fatal error:
Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
detected.

[1707228579.325814] [gpu008:1539109:0]           ib_md.c:1234 UCX  WARN  IB: ibv_fork_init() was disabled or failed, yet a fork() has been issued.
[1707228579.325821] [gpu008:1539109:0]           ib_md.c:1235 UCX  WARN  IB: data corruption might occur when using registered memory.
```
# My batchfile:
```
#!/bin/bash
#SBATCH --job-name="*"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=*
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-gpu
#SBATCH --gpus=1
#========================================
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GROMACS/2023.1-CUDA-11.7.0
export OMP_NUM_THREADS=8

srun gmx_mpi mdrun -ntomp 8 -s *.tpr -deffnm * -nsteps -1 -v -pin on -nb gpu -noappend -cpi *.cpt
```
Best,
Jingze


## Post 2 by @Yann.Sagon (2024-02-07T07:45:43.047Z)

Dear @Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan)
I think I got it! Your job was running fine on `gpu007.yggdrasil` and no GPUs are seen when same job is running on `gpu008.yggdrasil`.
Gpu007 is equipped with 8 x Titan RTX cards and gpu008 is equipped with 8x V100 cards (see hpc:hpc_clusters [eResearch Doc][hpc:hpc_clusters [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/hpc_clusters#gpus_on_yggdrasil)). The V100 has a compute capability of 7.0. Until 06th of December 2023, all the software using GPUs we compiled were compiled for compute capability : 6.0,6.1,7.5,8.0,8.6. We added the 7.0 since then, but not all the software were recompiled.
We saw later the following statement from Nvidia:
> Each CUBIN file targets a specific compute capability version and is forward- compatible only with CUDA architectures of the same major version number; e.g., CUBIN files that target compute capability 1.0 are supported on all compute- capability 1.x (Tesla) devices but are not supported on compute-capability 2.0 (Fermi) devices.
It means we no GPU kernel were compatible for this card.
We are recompiling GROMACS right now.


## Post 3 by @Jingze.Duan (2024-02-07T09:13:36.036Z)

Thank you for your help! My jobs are all fine now.


## Post 4 by @Adrien.Albert (2024-02-08T13:55:40.484Z)

Dear @Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan)
Gromac has been recompiled.
