# Srun error gpu006 task 6 Out Of Memory

- Source: https://hpc-community.unige.ch/t/3258

- Created: 2024-01-19T13:51:39.573Z

- Posts: 6

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jingze.Duan (2024-01-19T13:51:39.629Z)

If you are asking for help, try to provide information that can help us solve your issue, such as :
what did you try:
My job stopped after run only about 1h on baobab.
JobID: 6793261
what was the error message: (from .out file)
error: Detected 1 oom_kill event in StepId=6793261.0. Some of the step tasks have been OOM Killed.
srun: error: gpu006: task 6: Out Of Memory
the batchfile:
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
srun gmx_mpi mdrun -deffnm push -s push.tpr -v -pin on -nb gpu -cpi push.cpt -noappend
```


## Post 2 by @Adrien.Albert (2024-01-19T14:06:18.270Z)

@Jingze.Duan[@Jingze.Duan](https://hpc-community.unige.ch/u/jingze.duan)
It seems you forget to specify the memory you need for your job. By default slurm will assign 3 GB of memory.


## Post 3 by @Jingze.Duan (2024-01-19T14:09:16.605Z)

Thank you for your reply. But the output of my job was no more than 200 MB data for running 1h.
Btw, how can I specify the memory?


## Post 4 by @Adrien.Albert (2024-01-19T14:30:22.228Z)

Hi,
Is your log trace a real-time resource tracking?
From slurm documentation : Slurm Workload Manager - sbatch[Slurm Workload Manager - sbatch](https://slurm.schedmd.com/sbatch.html)
> –mem=<size>[units]
> Specify the real memory required per node. Default units are megabytes. Different units can be specified using the suffix [K|M|G|T]. Default value is DefMemPerNode and the maximum value is MaxMemPerNode. If configured, both parameters can be seen using the scontrol show config command. This parameter would generally be used if whole nodes are allocated to jobs (SelectType=select/linear). Also see –mem-per-cpu and –mem-per-gpu. The –mem, –mem-per-cpu and –mem-per-gpu options are mutually exclusive. If –mem, –mem-per-cpu or –mem-per-gpu are specified as command line arguments, then they will take precedence over the environment.
> NOTE: A memory size specification of zero is treated as a special case and grants the job access to all of the memory on each node.


## Post 5 by @Ludovic.Dumoulin (2024-01-29T15:49:44.968Z)

Hello,
I have similar issue since the update to 1.9.3 from 1.8.x.
It seems that it is coming from the GC of julia. This problem is solved with the new version 1.10 (from december 2023).
I’ll ask for the update of julia to 1.10.
Thank you,
Best regards
EDIT: you can manually call the GC using GC.gc()


## Post 6 by @Yann.Sagon (2024-01-30T14:04:49.900Z)

Hi, please try the new version we just installed: New software installed: Julia version 1.10.0[New software installed: Julia version 1.10.0](https://hpc-community.unige.ch/t/new-software-installed-julia-version-1-10-0/3278)
