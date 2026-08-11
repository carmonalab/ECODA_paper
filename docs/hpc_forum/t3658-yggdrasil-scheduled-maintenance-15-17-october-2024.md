# Yggdrasil scheduled maintenance: 15-17 October 2024

- Source: https://hpc-community.unige.ch/t/3658

- Created: 2024-09-30T14:55:19.671Z

- Tags: yggdrasil

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-09-30T14:55:19.751Z)

Dear users,
We will perform software and hardware maintenance on the Yggdrasil HPC cluster on 15-17 October 2024.
The maintenance will start at 08:00 +0100 and you will receive an email when the maintenance is finished.
The cluster will be completely unavailable during this time, with no access at all (not even to retrieve files).
If you submit a job in the meantime, make sure that the expected wall time (duration) does not overlap with the start of the maintenance, or your job will be scheduled after the maintenance.
What will be done during this maintenance:
- Update all servers to the latest security and bugfix releases (Rocky 8.10)
- Re-install all nodes with latest security and bugfix releases (Rocky 8.10)
- Install BeeGFS version 7.4.5
- Upgrade Slurm to version 24.05 (this is a major version upgrade Slurm Workload Manager[Slurm Workload Manager](https://slurm.schedmd.com/archive/slurm-24.05.0/release_notes.html))
- Install new memory on admin servers
Thank you for your understanding.
Best regards,
The HPC Team


## Post 2 by @Gael.Rossignol (2024-10-15T09:26:44.377Z)




## Post 3 by @Yann.Sagon (2024-10-17T07:35:02.703Z)




## Post 4 by @Yann.Sagon (2024-10-17T07:36:03.146Z)

Dear users, the maintenance is now finished, thank you for your patience.
What has been done:
- BeeGFS update to version 7.4.5
- Slurm updated to version 24.05.2
- Rocky updated to version 8.10
- New memory on two servers
- All compute nodes reinstalled, all servers updated
- Deploy CIFS Work Arround on login node
- Various minor fixes
Best regards
HPC Team
