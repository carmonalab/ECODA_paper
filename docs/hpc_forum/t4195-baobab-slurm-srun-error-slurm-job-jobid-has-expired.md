# [baobab][slurm] srun: error: Slurm job JOBID has expired

- Source: https://hpc-community.unige.ch/t/4195

- Created: 2026-01-16T14:03:41.041Z

- Tags: baobab, slurm

- Posts: 1

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-01-16T14:03:41.092Z)

Dear @hpc[@hpc](https://hpc-community.unige.ch/u/hpc) Team,
Running my sbtach job I got the following error in the out file:
```
(baobab)-[falkiewi@login1 ~]$ cat /home/users/f/falkiewi/fixing-cfg/apptainer/logs/6572146.out
srun: error: Slurm job 6572146 has expired
srun: Check SLURM_JOB_ID environment variable. Expired or invalid job 6572146
```
From the squeue perspective, the job is still “running”.
Could you please confirm that the problem is on my side (what could it be then? :thinking: ) not an SLURM problem?
Kind regards,
Maciej
