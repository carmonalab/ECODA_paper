# Jobs fail with error code 15/53 on baobab

- Source: https://hpc-community.unige.ch/t/3644

- Created: 2024-09-23T09:49:12.148Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-09-23T09:49:12.192Z)

Hi,
All my jobs, irrespective of the nature of the actual code being run, are quitting as soon as they are run, and no output or error files are produced. `seff` says the error code is 15
```
Job ID: 12695603
Cluster: baobab
User/Group: radhakrb/hpc_users
State: FAILED (exit code 15)
Cores: 1
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:00:00 core-walltime
Job Wall-clock time: 00:00:00
Memory Utilized: 1.29 MB
Memory Efficiency: 0.00% of 30.00 GB
```
while using `sacct` shows an exit code of 53
```
JobID           JobName  Partition    Account  AllocCPUS      State ExitCode 
------------ ---------- ---------- ---------- ---------- ---------- -------- 
12695603           mcmc private-d+     sonner          1     FAILED     0:53 
12695603.ba+      batch                sonner          1  CANCELLED     0:53 
12695603.ex+     extern                sonner          1  COMPLETED      0:0
```
Can the admins please have a look at the issue? Thanks!


## Post 2 by @Adrien.Albert (2024-09-24T09:42:34.114Z)

hi @Bharathkumar.Radhakrishnan[@Bharathkumar.Radhakrishnan](https://hpc-community.unige.ch/u/bharathkumar.radhakrishnan)
Please, could you share your sbatch and all relevant information about your job?
From Slurm Workload Manager - Job Exit Codes[Slurm Workload Manager - Job Exit Codes](https://slurm.schedmd.com/job_exit_code.html)
> # Job Exit Codes
> A job’s exit code (aka exit status, return code and completion code) is captured by Slurm and saved as part of the job record. For sbatch jobs, the exit code that is captured is the output of the batch script. For salloc jobs, the exit code will be the return value of the exit call that terminates the salloc session. For srun, the exit code will be the return value of the command that srun executes.
> Any non-zero exit code will be assumed to be a job failure and will result in a Job State of FAILED with a Reason of “NonZeroExitCode”.
This means that the exit code 15 originates from your script, indicating an issue within it rather than an error from Slurm itself.
