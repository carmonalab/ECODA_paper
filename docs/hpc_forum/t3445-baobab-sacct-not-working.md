# [baobab] sacct not working

- Source: https://hpc-community.unige.ch/t/3445

- Created: 2024-05-13T10:07:54.463Z

- Posts: 5

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-05-13T10:07:54.513Z)

what did you try: `$ sacct`
what was the error message:
```
sacct: error: slurm_persist_conn_open_without_init: failed to open persistent connection to host:slurm1:6819: Connection refused
sacct: error: Sending PersistInit msg: Connection refused
sacct: error: Problem talking to the database: Connection refused
```
Kind regards,
Maciej Falkiewicz


## Post 2 by @Adrien.Albert (2024-05-14T10:50:12.897Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
We are currently preparing the next maintenance and actions on slurm. There may be an interruption to slurm database queries (ex: sacct command), but this has no impact on jobs in progress. We are trying to minimize the interruption during lunch time  and it should not exceed 10-15 minutes.


## Post 4 by @Adrien.Albert (2024-05-15T12:31:18.550Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
It’s working for me with sinfo/sacct/squeue
Could you give me more information please


## Post 5 by @Genevieve.Savard (2024-06-13T12:34:31.892Z)

Hello,
this error just appeared on Yggdrasil today.
```
sacct: error: slurm_persist_conn_open_without_init: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Network is unreachable
sacct: error: Sending PersistInit msg: Network is unreachable
sacct: error: Problem talking to the database: Network is unreachable
```
cheers
Genevieve


## Post 6 by @Adrien.Albert (2024-06-13T12:48:26.533Z)

Dear @Genevieve.Savard[@Genevieve.Savard](https://hpc-community.unige.ch/u/genevieve.savard),
We are performing maintenance on Baobab; slurm utilities/functions may be interrupted due to an update of the server hosting the slurm database.
It should work now :slight_smile:
