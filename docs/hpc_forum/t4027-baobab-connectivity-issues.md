# Baobab connectivity issues

- Source: https://hpc-community.unige.ch/t/4027

- Created: 2025-08-01T22:20:06.435Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @William.Ceva (2025-08-01T22:20:06.481Z)

Hello,
Whenever I connect to baobab (whether by x2go or direct ssh), the terminal freezes when I type “ls” or any similar command, in my default home directory.
If I instead cd into a directory under my home directory first, and then do ls, I do not have any problems.
Any help fixing this issue is appreciated,
Will Ceva


## Post 2 by @Malte.Algren (2025-08-04T07:39:17.594Z)

The issue still persists
For me, it is only when i ls/touch/cd etc into or near scratch
Malte


## Post 3 by @Adrien.Albert (2025-08-04T10:10:27.033Z)

hello
Here the descriptionof the issue. We apologize for the inconvenience caused during this long week-end :confused:
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
