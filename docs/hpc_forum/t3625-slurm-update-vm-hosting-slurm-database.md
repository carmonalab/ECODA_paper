# [Slurm] Update VM hosting slurm database

- Source: https://hpc-community.unige.ch/t/3625

- Created: 2024-09-04T13:45:19.773Z

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-09-04T13:45:19.845Z)

Dear Users,
Please be advised that the VM hosting the SLURM database will undergo a reboot to finalize the his update. Access to the database may be impacted. (`sacct` - `squeue` - `sinfo` etc…)
The expected downtime is approximately less than 5 minutes.
Please note: Compute nodes and running jobs will not be impacted during this time.
Thank you for your understanding.
Best regards,


## Post 2 by @Adrien.Albert (2024-09-04T14:02:46.801Z)

The vm has been successfully updated in 43s, slurm command are back again  :sunglasses:
```
(admin)-[root@slurmdb ~]$ reboot
[root@localhost ~]$ ssh slurmdb
(admin)-[root@slurmdb ~]$ date
Wed Sep  4 15:55:43 CEST 2024
```
