# Yggdrasil scheduled maintenance: 28-31 July 2026

- Source: https://hpc-community.unige.ch/t/4351

- Created: 2026-07-17T14:33:01.933Z

- Tags: yggdrasil

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2026-07-17T14:33:02.104Z)

Dear users,
We will be performing software and hardware maintenance on the Yggdrasil HPC cluster from 28 to 31 July 2026, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
We remember you that you our two other HPC clusters are available during the maintenance!
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after.
What will be done during this maintenance
- All compute nodes will be reinstalled with latest Rocky 9
- Upgrade Slurm to version 26.05-1-1[26.05-1-1](https://github.com/SchedMD/slurm/releases/tag/slurm-26-05-1-1)
- Upgrade OpenOnDemand to 4.2.0
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The Baobab HPC Team


## Post 2 by @Adrien.Albert (2026-07-28T07:58:07.676Z)




## Post 3 by @Yann.Sagon (2026-07-30T15:00:27.877Z)




## Post 4 by @Yann.Sagon (2026-07-30T15:17:28.012Z)

Dear users,
The maintenance work is now complete. The following work was carried out:
- Slurm was updated to version 26.05.2
- Rocky was updated to version 9.8.
- OpenOnDemand was updated to version 4.2.
- All the servers were reinstalled with the latest security updates and patches.
- The Infiniband stack is now the one provided by Rocky Linux.
Best regards,
HPC Team
