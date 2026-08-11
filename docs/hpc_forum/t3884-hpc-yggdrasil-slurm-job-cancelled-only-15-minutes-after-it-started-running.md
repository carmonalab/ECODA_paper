# [HPC][yggdrasil][slurm] Job cancelled only 15 minutes after it started running

- Source: https://hpc-community.unige.ch/t/3884

- Created: 2025-03-25T05:50:57.816Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Monika.Avila (2025-03-25T05:50:57.858Z)

Hello,
Username: avilamar
Cluster: yggdrasil
Subject: <slurm| job cancelled after 15 minutes running| R >
jobid: 39283257
My jobs are cancelled before they have finished even if I ask for 10 hours of wall time. The code runs perfectly in my Mac mini even for days.
- I submit a job with slurm using a batch file.  The code is in R.
- I am trying to estimate a model using R, I am using different libraries of machine learning.
- sbtach main.sh
- The path is yggdrasil/avilamar/P4_MonteCarloSimulations
I would appreciate your assistance on this issue
Best,
Monika
avilamar


## Post 2 by @Yann.Sagon (2025-03-25T09:29:58.685Z)

Dear @Monika.Avila[@Monika.Avila](https://hpc-community.unige.ch/u/monika.avila)
Monika.Avila:
> sbtach main.sh
How do you submit your job (full command line used to start your job) ?
