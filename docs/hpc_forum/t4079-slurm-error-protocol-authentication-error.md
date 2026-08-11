# SLURM error : Protocol authentication error

- Source: https://hpc-community.unige.ch/t/4079

- Created: 2025-09-04T09:10:29.676Z

- Tags: baobab

- Posts: 6

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibaut.Chataing (2025-09-04T09:10:29.734Z)

Hello,
I have a problem with slurm in baobab.
When I try to run sbatch ou sinfo or squeue I get :
```
slurm_load_jobs error: Protocol authentication error
```
I tried a new terminal. It does not change anything.
Can you help, please ?
Best,
Thibaut


## Post 2 by @pablo.strasser1 (2025-09-04T09:22:12.469Z)

Same and when I disconnected, I’m unable to connect back in with ssh.


## Post 3 by @Yann.Sagon (2025-09-04T09:24:20.136Z)

Related to [2025] Current issues on HPC Cluster - #22 by Yann.Sagon[[2025] Current issues on HPC Cluster - #22 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/22)


## Post 4 by @Alexander.Froch (2025-09-04T09:27:31.061Z)

Not sure if this is related, but I think the login1 node of Baobab just crashed. I cannot login and my active sessions had timeouts


## Post 5 by @Yann.Sagon (2025-09-04T09:59:51.254Z)

Alexander.Froch:
> Not sure if this is related, but I think the login1 node of Baobab just crashed.
Yes, related. It is now fixed.
Best


## Post 6 by @Alexander.Froch (2025-09-04T12:30:22.860Z)

Thanks a lot Yann! Working again perfectly.
Cheers,
Alex
