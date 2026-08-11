# Yggdrasil - slurm not scheduling jobs

- Source: https://hpc-community.unige.ch/t/4158

- Created: 2025-12-05T14:04:21.369Z

- Tags: yggdrasil, slurm

- Posts: 4

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2025-12-05T14:04:21.421Z)

Dear @hpc[@hpc](https://hpc-community.unige.ch/u/hpc) Team,
There is a problem with slurm:
```
          43526839 public-cp                                                                 interactive falkiewi CG       0:04      1               cpu078 normal 2025-12-05T14:57:07
          43526838 public-cp                                                                 interactive falkiewi CG       2:47      1               cpu074 normal 2025-12-05T14:54:14
          43526841 public-cp                                                                 interactive falkiewi  R       6:12      1               cpu081 normal 2025-12-05T14:57:37
```
```
salloc: Pending job allocation 43526841
salloc: job 43526841 queued and waiting for resources
```
Kind regards,
Maciej Falkiewicz


## Post 2 by @Adrien.Albert (2025-12-05T14:45:34.083Z)

Hello,
I fixed something. Could you try again and let me know if it works fine ?


## Post 3 by @maciej.falkiewicz (2025-12-05T15:39:39.090Z)

Yes, my job has been granted allocation. Thank you for the quick fix!


## Post 4 by @maciej.falkiewicz (2025-12-05T17:16:07.644Z)

Dear @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/adrien.albert) ,
Now I experience the exact same problem on Baobab.
Kind regards,
Maciej Falkiewicz
