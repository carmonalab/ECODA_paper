# Baobab scheduled maintenance: 05-07 November 2024

- Source: https://hpc-community.unige.ch/t/3691

- Created: 2024-10-18T14:45:36.196Z

- Tags: baobab

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-10-18T14:45:36.292Z)

Dear users,
As anounced by email we will perform software and hardware maintenance on the Baobab HPC cluster on 05-07 November 2024.
The maintenance will start at 08:00 +0100 and you will receive an email when the maintenance is finished.
The cluster will be completely unavailable during this time, with no access at all (not even to retrieve files).
If you submit a job in the meantime, make sure that the expected wall time (duration) does not overlap with the start of the maintenance, or your job will be scheduled after the maintenance.
What will be done during this maintenance:
- Upgrade OpenOnDemand to 3.1.7
- Update all servers to the latest security and bugfix releases (Rocky 8.10)
- Re-install all nodes with latest security and bugfix releases (Rocky 8.10)
- Install BeeGFS version 7.4.5
- Upgrade Slurm to version 24.05 (this is a major version upgrade Slurm Workload Manager - Release Notes 2)
- Upgrade infiniband card to latest firmware
Thank you for your understanding.


## Post 2 by @Adrien.Albert (2024-11-05T07:56:54.914Z)




## Post 3 by @Adrien.Albert (2024-11-07T14:04:13.553Z)




## Post 4 by @Adrien.Albert (2024-11-07T14:08:40.203Z)

Dear users, the maintenance is now finished, thank you for your patience.
What has been done:
- BeeGFS update to version 7.4.5
- Slurm updated to version 24.05.2
- Rocky updated to version 8.10
- New memory on two servers
- All compute nodes reinstalled, all servers updated
- Various minor fixes
Best regards
HPC Team
