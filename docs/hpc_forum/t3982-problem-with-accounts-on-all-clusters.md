# Problem with accounts on all clusters

- Source: https://hpc-community.unige.ch/t/3982

- Created: 2025-06-16T14:40:19.994Z

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-06-16T14:40:20.038Z)

It seems that the server in charge of the login/account association is down:
On On-Demand, the account field has no option.
On login, the command: $ sacctmgr show user where name=delislel withassoc
Gives:
sacctmgr: error: _open_persist_conn: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Connection refused
sacctmgr: error: Sending PersistInit msg: Connection refused
Thanks,
Lucille


## Post 2 by @Gael.Rossignol (2025-06-20T07:30:29.381Z)

Dear Lucille,
In order to update cluster Bamboo slurm version, we has to update global slurmdbd database. This is completely transparents jobside, but some informtions like sacctmgr are not available during this database update.
Sorry for inconvenience,
Best regards,


## Post 3 by @Malte.Algren (2025-06-22T19:26:59.461Z)

When will this update be finished?
I’m getting the same issue on Baobab when submitting Snakemake jobs.
Cheers,
Malte


## Post 4 by @Yann.Sagon (2025-06-25T08:11:55.660Z)

Hi, this was already finished the 16th of June, the duration of the interruption was only one hour.
What issue do you have @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren) ?


## Post 5 by @Malte.Algren (2025-06-25T08:50:00.709Z)

It has been fixed now.
Thanks!
Malte


## Post 6 by @Lucille.Delisle1 (2025-06-29T22:02:00.702Z)

Hi,
It seems that there is again an issue with the slurm database:
```
$ sacctmgr show user where name=delislel withassoc
sacctmgr: error: _open_persist_conn: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Connection timed out
sacctmgr: error: Sending PersistInit msg: Connection timed out
```
Thanks


## Post 7 by @Cedric.Renaud (2025-06-30T05:25:24.343Z)

Hello everyone, I have the same issue.
Thanks


## Post 8 by @Lucille.Delisle1 (2025-06-30T08:01:38.727Z)

Hi,
The problem is solved. Thank you.


## Post 9 by @Yann.Sagon (2025-06-30T08:25:29.630Z)

Hi, I’ve opened a ticket at schedmd (Slurm editor) as it seem it the service crashes for no reason. In the meantime, service restarted.


## Post 10 by @Gael.Rossignol (2025-09-22T13:20:41.144Z)

2 posts were split to a new topic: Slurm db unavailable[Slurm db unavailable](https://hpc-community.unige.ch/t/slurm-db-unavailable/4096)
