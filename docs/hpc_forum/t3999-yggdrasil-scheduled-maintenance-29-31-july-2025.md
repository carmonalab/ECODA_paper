# Yggdrasil scheduled maintenance: 29-31 July 2025

- Source: https://hpc-community.unige.ch/t/3999

- Created: 2025-07-07T13:53:59.920Z

- Tags: yggdrasil

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-07-07T13:54:00.024Z)

Dear users,
As just announced on the baobab-announce@ mailing list, we will be performing software and hardware maintenance on the Yggdrasil HPC cluster from 29 to 31 July 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with Rocky 9.6
- Upgrade Slurm to version 25.05.0.
- Update all servers to the latest security and bugfix releases.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Adrien.Albert (2025-07-29T06:54:54.274Z)




## Post 3 by @Adrien.Albert (2025-07-31T08:49:53.725Z)




## Post 4 by @Adrien.Albert (2025-07-31T08:55:28.603Z)

Dear users,
The scheduled maintenance has been successfully completed, and Yggdrasil is now up and fully operational, ready to run your jobs!
Here’s a summary of the key updates:
- Slurm has been upgraded to version 25.05.0.
- Compute nodes remain on Rocky Linux 9.5, as Infiniband drivers are not yet available for version 9.6.
- On the login node, we have introduced a limit on the number of concurrent tasks per user, to help prevent excessive load and ensure stability for all users.
As always, if you notice any issues or unusual behavior, please don’t hesitate to open a new post on hpc-community.
Best regards,
The HPC Team
