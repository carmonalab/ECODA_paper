# Bamboo scheduled maintenance: 03-05 December 2024

- Source: https://hpc-community.unige.ch/t/3732

- Created: 2024-11-19T08:22:40.483Z

- Tags: bamboo

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-11-19T08:22:40.556Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we will perform software and hardware maintenance on the Bamboo HPC cluster from December 3 to 5, 2024, starting at 08:00 +0100.
You’ll receive an email when it ends. The cluster will be completely unavailable during this time, with no access at all (not even to retrieve files).
If you submit a job, ensure the expected wall time (duration) doesn’t overlap with the maintenance start, or it will be scheduled afterward. What will be done during this maintenance:
- Re-install all compute nodes with a new OS major version, Rocky 9. The compute nodes are in Rocky 8 currently.
- Upgrade Slurm to version 24.05.4
- Update all servers to the latest security and bugfix releases.
- Install new hard disks in home and scratch storage as spare
Thank you for your understanding. Best regards, The HPC Team


## Post 2 by @Adrien.Albert (2024-12-03T09:51:46.277Z)




## Post 3 by @Adrien.Albert (2024-12-05T11:29:06.284Z)
