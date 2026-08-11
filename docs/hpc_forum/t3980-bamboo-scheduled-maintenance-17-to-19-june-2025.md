# Bamboo scheduled maintenance: 17 to 19 June 2025

- Source: https://hpc-community.unige.ch/t/3980

- Created: 2025-06-11T15:16:25.207Z

- Tags: bamboo

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-06-11T15:16:25.308Z)

Dear users,
We will be performing software and hardware maintenance on the Bamboo HPC cluster from 17 to 19 June 2025, starting at 08:00 +0100. We apologize for such a short notice, next possible slot is too far in the future.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with Rocky 9.6
- Upgrade Slurm to version 24.11.5 or 25.05.0 if we have the time
- Update all servers to the latest security and bugfix releases.
- Various fixes such as RAID drivers change from legacy to nvme.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Yann.Sagon (2025-06-17T07:07:07.323Z)




## Post 3 by @Yann.Sagon (2025-06-19T07:50:32.825Z)

Dear users,
The maintenance work has now been completed and Bamboo is up and running, ready for your jobs!
We updated Slurm to version 25.05.0 and SlurmDB to this version for accounting on the three clusters.
The compute nodes are still using Rocky 9.5, as there is no driver available for the Infiniband network in version 9.6.
As always, please post feedback if you notice anything unusual.
Best regards,
The HPC Team


## Post 4 by @Yann.Sagon (2025-06-19T07:51:01.566Z)
