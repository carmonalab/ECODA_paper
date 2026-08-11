# Sbatch: I/O error writing script/environment to file

- Source: https://hpc-community.unige.ch/t/3464

- Created: 2024-06-01T07:57:45.667Z

- Posts: 6

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-06-01T07:57:45.706Z)

Hello,
There seems to be a problem with submitting batch jobs. When using sbatch with even a simple file like
```
#!/bin/bash
#SBATCH --job-name=trial
#SBATCH --partition=private-dpt-cpu
#SBATCH --time=00-00:01:00
#SBATCH -N 1 # total number of nodes
#SBATCH --mem=10G

echo Hello
```
I get an error saying
`sbatch: error: Batch job submission failed: I/O error writing script/environment to file`. Somehow running `srun echo hello` seems fine and works.
Can you please have a look?
Best
Bharath


## Post 2 by @Lorenzo.Bini (2024-06-01T08:17:06.904Z)

Hello, I’m getting the same error on shared-gpu partition, while everything seems fine on other —gpu partitions…


## Post 3 by @Bharathkumar.Radhakrishnan (2024-06-01T14:17:03.814Z)

This is the error I get now
```
Batch job submission failed: Unable to contact slurm controller (connect failure)
```


## Post 4 by @Stephen.Mulligan (2024-06-01T18:01:59.716Z)

also getting this error with sbatch, salloc, and squeue


## Post 5 by @Aravind.Srinivasan (2024-06-03T04:52:45.820Z)

I am also getting the same error with sbatch on shared and public-cpu.  `Batch job submission failed: Unable to contact slurm controller (connect failure)`


## Post 6 by @Yann.Sagon (2024-06-04T08:07:00.058Z)

Dear all, the issue was this one and it is now solved: [2024] Current issues on HPC Cluster - #11 by Yann.Sagon[[2024] Current issues on HPC Cluster - #11 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/11)
Best
Yann
