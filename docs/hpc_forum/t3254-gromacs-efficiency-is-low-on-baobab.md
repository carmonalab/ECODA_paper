# GROMACS efficiency is low on Baobab

- Source: https://hpc-community.unige.ch/t/3254

- Created: 2024-01-16T09:03:55.887Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-01-16T09:03:55.978Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
what did you try:
I had running GROMACS on the baobab with the following script, but it had a low efficiency of 6.7 ns/day. My MD simulation system has about 600,000 atoms.
#!/bin/bash
#SBATCH --job-name=“*”
#SBATCH --mail-type=ALL
#SBATCH --mail-user= *
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-gpu
#SBATCH --gpus=1
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GROMACS/2023.1-CUDA-11.7.0
export OMP_NUM_THREADS=8
srun gmx_mpi mdrun -deffnm * -s *.tpr -v -cpi *.cpt -noappend -pin on -nb gpu
what was the expected result:
It should be an efficiency of about 15 ns/day, compared with the simulation results of my colleague who had the same setting of gpu on yggdrasil and similar size of system.


## Post 2 by @Gael.Rossignol (2024-04-04T09:26:29.150Z)

Hello,
I think this issue is solved, please read this post :
No GPU is detected on yggdrasil[No GPU is detected on yggdrasil](https://hpc-community.unige.ch/t/no-gpu-is-detected-on-yggdrasil/3291/4) HPC Technical[HPC Technical](https://hpc-community.unige.ch/c/hpc-technical/5)
> Dear @Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan) 
Gromac has been recompiled.
Best regards,
