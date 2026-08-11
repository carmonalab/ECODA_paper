# BadConstraints reason in squeue, job not launching

- Source: https://hpc-community.unige.ch/t/3773

- Created: 2024-12-19T14:04:13.499Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Louis.Perrin (2024-12-19T14:04:13.540Z)

Dear HPC staff,
First of all thank you for all your help this year.
I’m writing you regarding an issue that occurred today on the baobab cluster when I was trying to run a gromacs simulation: I submitted a job about two days ago that was queuing because no nods were available on the partition, but it appears like the reason changed to “BadConstraints” (job ID 14072175). When I tried to copy/paste the folder and launch another job to check what was going on, it appears that I couldn’t launch a new job with the exact same batch file (except for the job name). I got the following error when submitting my job, which didn’t appear this year with the same batch file:  sbatch: error: Batch job submission failed: Requested node configuration is not available
Is there something I am doing wrong in my batchfile by any chance? I join it hereunder:
```
#!/bin/bash
#SBATCH --job-name="SMDAGLC4"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=louis.perrin@unige.ch
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared-gpu
#SBATCH --gpus=turing:1
#========================================
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GROMACS/2023.1-CUDA-11.7.0
export OMP_NUM_THREADS=8

srun gmx_mpi grompp -f smd_conseq.mdp -o smd.tpr -c step6.6_equilibration.gro -r step5_input.pdb -p topol.top -n index.ndx -maxwarn 1
srun gmx_mpi mdrun -ntomp 8 -deffnm smd -pin on -nb gpu -cpi smd.cpt
```
Best regards,
Louis Perrin


## Post 2 by @Yann.Sagon (2024-12-19T15:58:35.535Z)




## Post 3 by @Yann.Sagon (2024-12-19T16:01:37.524Z)

Dear @Louis.Perrin[@Louis.Perrin](https://hpc-community.unige.ch/u/louis.perrin)
Louis.Perrin:
> First of all thank you for all your help this year.
Thanks!
Louis.Perrin:
> I submitted a job about two days ago that was queuing because no nods were available on the partition, but it appears like the reason changed to “BadConstraints”
The reason is because we changed the GPU type naming scheme[GPU type naming scheme](https://hpc-community.unige.ch/t/important-new-gpu-types-naming-in-baobab/3772).
Best
