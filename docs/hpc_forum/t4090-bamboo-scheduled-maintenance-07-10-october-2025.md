# Bamboo scheduled maintenance: 07-10 October 2025

- Source: https://hpc-community.unige.ch/t/4090

- Created: 2025-09-17T09:24:19.190Z

- Tags: bamboo

- Posts: 4

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-09-17T09:24:19.305Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we will be performing software and hardware maintenance on the Bamboo HPC cluster from 07 to 10 October 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with Rocky 9.6
- Upgrade Slurm to version 25.05.3
- Update all servers to the latest security and bugfix releases.
- Change Slurm configuration to handle better newer AMD + GPU compute nodes (numa as socket)
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Gael.Rossignol (2025-10-07T06:13:11.488Z)




## Post 3 by @Yann.Sagon (2025-10-09T09:43:00.368Z)

Dear users,
Maintenance is now complete. Thank you for your patience. Have a happy computing day!
Best regards,
The HPC Team


## Post 4 by @Yann.Sagon (2025-10-09T09:43:31.452Z)




## Post 5 by @Raphael.Rubino (2025-10-13T08:02:55.732Z)

Hello,
Thanks for your work on Bamboo!
Not sure if it is since the maintenance, but GPU nodes on Bamboo are all in drain state.
Best regards


## Post 6 by @Yann.Sagon (2025-10-15T11:47:43.008Z)

@Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) this was the reason, solved now: [2025] Current issues on HPC Cluster - #25 by Yann.Sagon[[2025] Current issues on HPC Cluster - #25 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/25)
