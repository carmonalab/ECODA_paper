# Salloc Srun Communication connection failure

- Source: https://hpc-community.unige.ch/t/3915

- Created: 2025-04-10T08:24:40.076Z

- Tags: yggdrasil, slurm

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yuzheng.Kang (2025-04-10T08:24:40.121Z)

## Primary informations
Username: Kang
Cluster: Yggdrasil
## Description
I am trying to run an interactive job via my ssh terminal. Turns out it shows me communication connection failure recently. It wasn’t like this before. Here is what I have
```
(cp_tf_env) (yggdrasil)-[kang@login1 emcee_cosmopower]$ salloc -n1 -c2 --partition=public-cpu --time=4:15:00 --cpus-per-task=2 --mem=8000
salloc: Pending job allocation 39992877
salloc: job 39992877 queued and waiting for resources
salloc: job 39992877 has been allocated resources
salloc: Granted job allocation 39992877
salloc: Nodes cpu051 are ready for job
srun: error: Task launch for StepId=39992877.interactive failed on node cpu051: Communication connection failure
srun: error: Application launch failed: Communication connection failure
srun: Job step aborted
salloc: Relinquishing job allocation 39992877
```
Any clue on how can I solve this issue?


## Post 2 by @Adrien.Albert (2025-04-10T09:50:25.229Z)

Hello
The issue have been fixed. Let me know everything is working fine


## Post 3 by @Yuzheng.Kang (2025-04-10T10:39:53.838Z)

Hi,
Thanks for your help! It’s working now.
