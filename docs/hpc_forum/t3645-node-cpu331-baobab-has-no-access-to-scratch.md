# Node cpu331.baobab has no access to scratch

- Source: https://hpc-community.unige.ch/t/3645

- Created: 2024-09-23T10:04:14.199Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Andrea.Serpolla (2024-09-23T10:04:14.251Z)

## Primary informations
Username: serpolla
Cluster: baobab
## Description
I noticed some of my jobs failing while trying to write files on scratch.
They were all running on `cpu331.baobab` and indeed, after opening a shell, I noticed scratch is unavailable on that computing node:
```
(baobab)-[serpolla@login1 ~]$ srun -p shared-cpu --nodelist=cpu331 --pty bash -i
srun: job 12696695 queued and waiting for resources
srun: job 12696695 has been allocated resources
(baobab)-[serpolla@cpu331 ~]$ ls /srv/beegfs/scratch/
<Empty dir>
```
Until now excluding just that one node seems sufficient, but I do not know if other nodes are experiencing the same issue.


## Post 2 by @Malte.Algren (2024-09-23T10:40:11.228Z)

I have the same issue. I have access to scratch on the login node but when using ‘salloc’ to go on a compute node, I do not have access anymore.
As far as I can see its only for some nodes: gpu017 I don’t have access but on cpu209 I do
Malte


## Post 3 by @Adrien.Albert (2024-09-23T11:32:49.579Z)

Hi all;
I am on it. I keep you informed as soon as possible


## Post 4 by @Andrea.Serpolla (2024-09-29T18:03:40.367Z)

It seems fine now :+1:


## Post 5 by @Yann.Sagon (2024-09-30T12:25:53.844Z)

Dear @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren) and @Andrea.Serpolla[@Andrea.Serpolla](https://hpc-community.unige.ch/u/andrea.serpolla) you can follow the issue there: [2024] Current issues on HPC Cluster - #25 by Yann.Sagon[[2024] Current issues on HPC Cluster - #25 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/25)
Every node was set in drain, it means every new jobs should be fine.
Best
