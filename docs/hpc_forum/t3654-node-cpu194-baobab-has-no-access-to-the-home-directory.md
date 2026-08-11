# Node cpu194.baobab has no access to the home directory

- Source: https://hpc-community.unige.ch/t/3654

- Created: 2024-09-29T18:07:53.034Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Andrea.Serpolla (2024-09-29T18:07:53.072Z)

## Primary informations
Username: `serpolla`
Cluster: Baobab
## Description
Computing node `cpu194` seems to have problems with the home directory:
```
(baobab)-[serpolla@login1 ~]$ srun -p shared-cpu --nodelist cpu194 --pty -t 10 bash -i
srun: job 12791225 queued and waiting for resources
srun: job 12791225 has been allocated resources
slurmstepd: error: couldn't chdir to `/home/users/s/serpolla': No such file or directory: going to /tmp instead
slurmstepd: error: couldn't chdir to `/home/users/s/serpolla': No such file or directory: going to /tmp instead
(baobab)-[serpolla@cpu194 tmp]$ ls /home
salt
```


## Post 2 by @Yann.Sagon (2024-09-30T12:32:51.785Z)

Dear @Andrea.Serpolla[@Andrea.Serpolla](https://hpc-community.unige.ch/u/andrea.serpolla) this is related to [2024] Current issues on HPC Cluster - #25 by Yann.Sagon[[2024] Current issues on HPC Cluster - #25 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/25)
You can not resubmit your job, this should be fixed.
Best
