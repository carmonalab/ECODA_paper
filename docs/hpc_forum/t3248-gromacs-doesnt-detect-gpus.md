# GROMACS doesn't detect GPUs

- Source: https://hpc-community.unige.ch/t/3248

- Created: 2024-01-12T10:39:08.492Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Louis.Perrin (2024-01-12T10:39:08.581Z)

## Primary informations
Username: perrinlo
Cluster: Yggdrasil
Hello,
I have been running a few MD simulations using GROMACS for the past few days and got an error today: Fatal error: Cannot run short-ranged nonbonded interactions on a GPU because no GPU is detected. Hereunder is the batch file I used to run my job, which worked over the past few weeks but stopped working today. The command is supposed to run an MD simulation of a protein-membrane system and has been running without issue until today.
```
#!/bin/bash
#SBATCH --job-name="pore_16-"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=louis.perrin@unige.ch
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

srun gmx_mpi mdrun -v -deffnm step7_production -pin on -nb gpu -cpi step7_production.cpt
```
Best,
Louis
admin edit: job id 30122433


## Post 2 by @Yann.Sagon (2024-01-16T16:23:05.977Z)

@Maurice.Karrenbrock[@Maurice.Karrenbrock](https://hpc-community.unige.ch/u/maurice.karrenbrock) you had the same issue. You were using your own Gromacs version right?


## Post 3 by @Yann.Sagon (2024-01-16T16:24:59.180Z)

@Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan)   Related too? GROMACS efficiency is low on Baobab[GROMACS efficiency is low on Baobab](https://hpc-community.unige.ch/t/gromacs-efficiency-is-low-on-baobab/3254)


## Post 4 by @Yann.Sagon (2024-02-07T07:49:22.918Z)

Hi, for those who have the issue with Gromacs and GPUs. Can you let us know if it was maybe related to this issue? No GPU is detected on yggdrasil - #2 by Yann.Sagon[No GPU is detected on yggdrasil - #2 by Yann.Sagon](https://hpc-community.unige.ch/t/no-gpu-is-detected-on-yggdrasil/3291/2) . In summary, were you using gpu008.yggdrasil?
