# Baobab scheduled maintenance: 09-12 September 2025

- Source: https://hpc-community.unige.ch/t/4050

- Created: 2025-08-15T13:31:11.550Z

- Tags: baobab

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2025-08-15T13:31:11.653Z)

Dear users,
We will be performing software and hardware maintenance on the Baobab HPC cluster from 09 to 12 September 2025, starting at 08:00 +0100.
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap with the start of the maintenance, or it will be scheduled after. What will be done during this maintenance
- All compute nodes will be reinstalled with Rocky 9.6
- Upgrade Slurm to version 25.05.2.
- Update all servers to the latest security and bugfix releases.
- Use a private ip address for admin server: This is already the case for Bamboo. From the user POV: a compute node will be allowed to connect to the outside of UNIGE network only on port 80, 443 and 1049.
During the maintenance period we are very busy and provide limited user support.
Thank you for your understanding.
Best regards, The HPC Team


## Post 2 by @Adrien.Albert (2025-09-09T06:34:50.610Z)




## Post 3 by @Yann.Sagon (2025-09-10T07:20:37.575Z)

Dear users,
we’ll upgrade slurmdbd this morning. This means that commands such as sacct, sacctmgr won’t work on the three clusters.
Thanks for your understanding.
Best regards


## Post 4 by @Yann.Sagon (2025-09-10T16:29:07.802Z)

Dear users,
The maintenance work has now finished.
We still have a couple of problematic compute nodes that we’ll have to check on Friday, but we decided to put the cluster into production early so that you can take advantage of the “Jeûne Genevois” public holiday to carry out some heavy computations.
In the end, we didn’t modify the IP address of the admin server, so there are no changes for you.
Don’t forget to try out the new OpenOnDemand version. 4.0. on Baobab.   https://openondemand.baobab.hpc.unige.ch[https://openondemand.baobab.hpc.unige.ch](https://openondemand.baobab.hpc.unige.ch) and the Cluster status tab under “Clusters”
Best regards,
Your HPC team


## Post 5 by @Yann.Sagon (2025-09-15T07:11:06.503Z)
