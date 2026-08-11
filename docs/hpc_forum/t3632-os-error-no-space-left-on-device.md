# OS Error, no space left on device

- Source: https://hpc-community.unige.ch/t/3632

- Created: 2024-09-14T10:10:57.557Z

- Tags: baobab

- Posts: 5

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Debajyoti.Sengupta (2024-09-14T10:10:57.600Z)

Hello HPC team,
I am currently on gpu002 on baobab, trying to run a script. It seems to keep failing with the error `OSError: [Errno 28] No space left on device`.
I have tried to clean up my scratch as much as possible. Is this something from my end still? If so, how may I go about sorting this issue out?


## Post 2 by @Brian.Pulfer (2024-09-15T10:50:36.626Z)

Hi all,
I am facing the same issue. The scratch partition is full (`49G` available when running `df -h`, which I guess is some spare space).
Best,
Brian


## Post 3 by @Adrien.Albert (2024-09-15T18:40:28.248Z)

Hi,
The BeeGFS scratch directory is almost full, with only 49GB remaining.
```
(baobab)-[root@gpu002 ~]$ df -t beegfs -h
Filesystem      Size  Used Avail Use% Mounted on
beegfs_dpnc     503T  3.6T  499T   1% /srv/beegfs/dpnc
beegfs_home     138T  119T   20T  86% /home
beegfs_scratch  1.5P  1.5P   49G 100% /srv/beegfs/scratch
```
We kindly ask everyone reading this message to clean up the scratch directory as soon as possible. A formal communication will be sent tomorrow, but immediate action is necessary to free up space.


## Post 4 by @Debajyoti.Sengupta (2024-09-17T06:48:10.257Z)

Hi Adrien,
Thanks, I am trying to clean up as much as feasible. I have also noticed my trainings running considerably slower, which involves i/o. I was wondering if this was related to the scratch being saturated?


## Post 5 by @Yann.Sagon (2024-09-17T08:04:48.803Z)

Dear @Debajyoti.Sengupta[@Debajyoti.Sengupta](https://hpc-community.unige.ch/u/debajyoti.sengupta)
probably yes. How many space do you need on the storage? I’m asking because you have several other alternatives:
- Use Bamboo, its home storage is much faster than other storage we provide, but limited to 1TB
- Request access to our fast storage on Baobab, very limited in size and no backup at all. Only for temporary data
- Use local storage on compute node. Only for ephemeral (job duration)
Check here[here](https://doc.eresearch.unige.ch/hpc/storage_on_hpc) for more details.
Best
