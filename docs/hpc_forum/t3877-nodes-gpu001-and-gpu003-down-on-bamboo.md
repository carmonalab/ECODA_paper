# Nodes gpu001 and gpu003 down on bamboo

- Source: https://hpc-community.unige.ch/t/3877

- Created: 2025-03-20T14:42:49.154Z

- Tags: bamboo, slurm

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Brian.Pulfer (2025-03-20T14:42:49.191Z)

## Primary informations
Username: pulfer
Cluster: bamboo
## Description
Nodes `gpu001` and `gpu003` are down / drained.
## Steps to Reproduce
`salloc -n 1 -c 1 --time 01:00:00 --nodelist gpu001`
or
`salloc -n 1 -c 1 --time 01:00:00 --nodelist gpu003`
## Expected Result
New shell session in gpu001
## Actual Result
The allocation just hangs.
```
salloc: Requested partition configuration not available now
salloc: Pending job allocation 293409
salloc: job 293409 queued and waiting for resources
```


## Post 2 by @Adrien.Albert (2025-03-21T12:28:59.188Z)

Hi Brian,
gpu001 was drained for:
[Slurm] Create public-interactive-gpu partition[[Slurm] Create public-interactive-gpu partition](https://hpc-community.unige.ch/t/slurm-create-public-interactive-gpu-partition/3879/1) HPC Announce[HPC Announce](https://hpc-community.unige.ch/c/hpc-announce/6)
> Dear Users, 
Cluster: Bamboo, Yggdrasil 
We would like to inform you about a recent change to the GPU partition on the system. The gpu001 partition has been modified from debug-gpu to public-interactive-gpu. This change introduces the following constraints for all users: 
Key Changes:

Partition Name: gpu001 has been moved to the public-interactive-gpu partition.
Constraints:

1 Job per User: Each user is limited to running a single job at a time.
1 GPU per Job: Each job will be allocated only o…
gpu003 is drained to prepare the next week maintenance.
