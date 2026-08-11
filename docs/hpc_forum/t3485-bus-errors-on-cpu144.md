# Bus errors on cpu144?

- Source: https://hpc-community.unige.ch/t/3485

- Created: 2024-06-12T08:45:41.488Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2024-06-12T08:45:41.545Z)

User: briel
Cluster: yggdrasil
I’ve been running large job array and have been receiving `bus error`’s for jobs that were ran on cpu144.
All other jobs seem to run fine.
These jobs fail within the first 10 minutes after starting.
I am reading files from `isilon` and writing to my home directory.
Is there a specific issue with this CPU?


## Post 2 by @Yann.Sagon (2024-06-14T09:34:44.256Z)

Dear @Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel) thanks for the notification. Indeed there was memory error on the compute nodes. We are checking.
Best
Yann


## Post 3 by @Matthias.Kruckow (2024-06-18T08:55:42.563Z)

The cpu 144 still produces bus errors from time to time. I got one a few hours ago.


## Post 4 by @Yann.Sagon (2024-06-28T09:23:37.641Z)

@Matthias.Kruckow[@Matthias.Kruckow](https://hpc-community.unige.ch/u/matthias.kruckow) thanks for the notification, we’ll replace memory. Node in drain until further notice.
