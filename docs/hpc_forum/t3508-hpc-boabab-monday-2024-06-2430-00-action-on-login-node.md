# [HPC][Boabab] Monday 2024-06-24:12:30:00 Action on login node

- Source: https://hpc-community.unige.ch/t/3508

- Created: 2024-06-21T14:38:37.170Z

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2024-06-21T14:38:37.239Z)

Dear HPC Users,
You may have noticed issues with mounting CIFS/SMB shares, blocking some research project.
More info:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/15) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> All Clusters: GIO mount storage from NASAC
Dear User, 
We have recently discovered issues with mounting storage spaces using the CIFS/SMB protocol from NASAC. 
The problem comes from the update of Mellanox suite of packages, which manages our Infiniband (fast) network. Unfortunately, these packages disable the CIFS module, preventing CIFS/SMB storage mounts. 
Since the Mellanox suite is crucial for the cluster’s operation, we cannot remove it. Currently, Mellanox doesn’t provides a workaround. 
W…
As work-arround, we will apply a custom patch to the login1.baobab server. This will require a server restart, scheduled for June 24, 2024, at 12:30 PM. During this period, the server will be unavailable, and all active connections will be interrupted for approximately 25 minutes.
Please note that implementing this fix may affect the Infiniband status on the login node. Certain features may not work as usual. If these issues are significant, we may need to roll back the patch, requiring additional intervention.
This patch will not be deployed on compute nodes to guarantee optimal production.
We will keep you informed of our progress and notify you once the maintenance is complete.
Thank you for your understanding.
Best regards,
–
HPC Team


## Post 2 by @Adrien.Albert (2024-06-24T09:07:03.121Z)




## Post 3 by @Adrien.Albert (2024-06-24T09:07:42.514Z)




## Post 4 by @Adrien.Albert (2024-06-24T11:29:36.755Z)

Dear HPC Users,
Due to unexpected behavior, intervention on login1.baobab must be extended.
The login node is not accessible, but you can access Baobab and continue your work via OpenOnDemand (there’s also a terminal).
Thank you for your understanding.
–
HPC team


## Post 5 by @Adrien.Albert (2024-06-24T12:14:46.304Z)

Dear HPC Users,
While we successfully tested the update in our development environment, we encountered unexpected behavior in the production environment that prevented us from replicating the expected results. As a result, we have reverted the changes, and there are currently no modifications affecting login node.
Summary:
- Issue: We have observed differences in behavior between our development and production environments.
- Current status: No changes have been made, the problem is still ongoing.
We appreciate your understanding and patience as we work towards a resolution.
Best Regards,
–
HPC team


## Post 6 by @Adrien.Albert (2024-06-25T13:45:06.277Z)
