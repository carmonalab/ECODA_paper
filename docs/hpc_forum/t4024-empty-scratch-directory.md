# Empty scratch directory

- Source: https://hpc-community.unige.ch/t/4024

- Created: 2025-07-31T16:23:13.783Z

- Tags: baobab

- Posts: 4

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Laure.Moinat (2025-07-31T16:23:13.832Z)

user: moinatl
Dear team,
After the problem with the login node earlier today, my scratch directory was completely empty. Is there an explication ?
Best,
Laure


## Post 2 by @Paul.Coppin (2025-07-31T16:32:39.497Z)

Same for me. It appears that the entire Baobab scratch is offline.


## Post 3 by @Adrien.Albert (2025-07-31T22:56:41.789Z)

Hello;
Sorry I forget to answer, scratch has been fixed and re-mount.
Let me know everything is working.


## Post 4 by @Adrien.Albert (2025-08-04T12:58:41.018Z)

[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/18) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: Baobab 
We experienced an issue on the scratch storage due to a timeout on one of the disks. This caused unresponsive data for a short period, which in turn led to a temporary crash of the BeeGFS service. 
Some users may have noticed latency on the scratch or error messages during this time. The concerned disk has been removed from the pool, and the issue is now resolved. 
Thank you for your understanding. 
– 
HPC Team 
Status : Resolved  green_circle
start: 2025-08-01T13:…


## Post 5 by @Adrien.Albert (2025-08-04T12:58:46.408Z)
