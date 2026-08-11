# [Yggdrasil] Could not resolve hostname gpu005: Name or service not known

- Source: https://hpc-community.unige.ch/t/3384

- Created: 2024-03-21T09:12:04.706Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Van-Khoa.Nguyen (2024-03-21T09:12:04.768Z)

Good morning,
Currently, I have a job allocated on gpu005. However, from the login node, I can not access the gpu005 compute node with ‘ssh gpu005’. Are there any problems happening on the cluster?
Bests


## Post 2 by @maciej.falkiewicz (2024-03-21T09:12:48.083Z)

Same here.
```
(base) (yggdrasil)-[falkiewi@login1 ksgan]$ salloc -n1 -c4 --mem=8000 --partition=public-cpu,shared-cpu,private-cui-cpu --time=1:00:00
salloc: Pending job allocation 32809996
salloc: job 32809996 queued and waiting for resources
salloc: job 32809996 has been allocated resources
salloc: Granted job allocation 32809996
salloc: Waiting for resource configuration
salloc: Nodes cpu117 are ready for job
srun: error: xgetaddrinfo: getaddrinfo(cpu117:6818) failed: Name or service not known
srun: error: slurm_set_addr: Unable to resolve "cpu117"
srun: error: _fwd_tree_get_addr: can't find address for host cpu117, check slurm.conf
srun: error: Task launch for StepId=32809996.interactive failed on node cpu117: Can't find an address, check slurm.conf
srun: error: Application launch failed: Can't find an address, check slurm.conf
srun: Job step aborted
salloc: Relinquishing job allocation 32809996
```


## Post 3 by @Yann.Sagon (2024-03-21T09:42:23.964Z)

@maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz) @Van-Khoa.Nguyen[@Van-Khoa.Nguyen](https://hpc-community.unige.ch/u/van-khoa.nguyen) thanks for the notification, we are fixing the issue. Please try again.
https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6[https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6)
