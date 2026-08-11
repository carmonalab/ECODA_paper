# Failed job without .out

- Source: https://hpc-community.unige.ch/t/3340

- Created: 2024-02-28T10:03:52.794Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-02-28T10:03:52.886Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
I submitted a job on baobab and it failed because I received emails. But I didn’t get any .out file to find the error information. No new file was generated.
Job ID:7704708
My batchfile:
```
#!/bin/bash
#SBATCH --job-name="KAGdm4_100"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Jingze.Duan@etu.unige.ch
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

srun gmx_mpi grompp -f step6.1_equilibration.mdp -c em.gro -r em.gro -p topol.top -o eq1.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq1 -pin on -nb gpu
srun gmx_mpi grompp -f step6.2_equilibration.mdp -c eq1.gro -r eq1.gro -p topol.top -o eq2.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq2 -pin on -nb gpu
srun gmx_mpi grompp -f step6.3_equilibration.mdp -c eq2.gro -r eq2.gro -p topol.top -o eq3.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq3 -pin on -nb gpu
srun gmx_mpi grompp -f step6.4_equilibration.mdp -c eq3.gro -r eq3.gro -p topol.top -o eq4.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq4 -pin on -nb gpu
srun gmx_mpi grompp -f step6.5_equilibration.mdp -c eq4.gro -r eq4.gro -p topol.top -o eq5.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq5 -pin on -nb gpu
srun gmx_mpi grompp -f step6.6_equilibration.mdp -c eq5.gro -r eq5.gro -p topol.top -o eq6.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq6 -pin on -nb gpu
srun gmx_mpi grompp -f step6.7_equilibration.mdp -c eq6.gro -r eq6.gro -p topol.top -o eq7.tpr -n index.ndx -maxwarn 2
srun gmx_mpi mdrun -deffnm eq7 -pin on -nb gpu -cpi eq7.cpt
```


## Post 2 by @Yann.Sagon (2024-03-04T14:34:10.453Z)

Hi @Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan)
I have no idea what is causing the issue.
Why are you specifiying `OMP_NUM_THREADS=8` and at the same time specifying `--cpus-per-task=1`? both number should be the same I guess.
Where is the path where you launched the sbatch script?
Best
Yann
