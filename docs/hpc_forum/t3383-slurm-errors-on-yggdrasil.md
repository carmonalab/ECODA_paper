# SLURM errors on Yggdrasil

- Source: https://hpc-community.unige.ch/t/3383

- Created: 2024-03-21T07:09:15.962Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Giuseppe.Chindemi (2024-03-21T07:09:16.020Z)

Hello,
Since yesterday I have troubles allocating resources on Yggdrasil, for example:
```
(yggdrasil)-[chindemi@login1 ~]$ salloc -N1 -n4 --mem=32G --gpus=1 --partition=shared-gpu --time=2:00:00
salloc: Pending job allocation 32808024
salloc: job 32808024 queued and waiting for resources
salloc: job 32808024 has been allocated resources
salloc: Granted job allocation 32808024
salloc: Nodes gpu007 are ready for job
srun: error: xgetaddrinfo: getaddrinfo(gpu007:6818) failed: Name or service not known
srun: error: slurm_set_addr: Unable to resolve "gpu007"
srun: error: _fwd_tree_get_addr: can't find address for host gpu007, check slurm.conf
srun: error: Task launch for StepId=32808024.interactive failed on node gpu007: Can't find an address, check slurm.conf
srun: error: Application launch failed: Can't find an address, check slurm.conf
srun: Job step aborted
salloc: Relinquishing job allocation 32808024
salloc: Job allocation 32808024 has been revoked.
```
Could you please have a look?
Cheers
Giuseppe


## Post 2 by @Yann.Sagon (2024-03-21T09:53:01.895Z)

Dear @Giuseppe.Chindemi[@Giuseppe.Chindemi](https://hpc-community.unige.ch/u/giuseppe.chindemi) this is solved.
https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6[https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/6)


## Post 3 by @Giuseppe.Chindemi (2024-03-21T12:28:51.790Z)

Great, thank you Yann!
