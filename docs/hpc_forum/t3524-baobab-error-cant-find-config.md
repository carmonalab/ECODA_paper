# Baobab error - can't find config

- Source: https://hpc-community.unige.ch/t/3524

- Created: 2024-07-03T17:19:35.146Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-07-03T17:19:35.193Z)

Over the last couple of hours, all my jobs have been failing with the error
```
srun: error: _fwd_tree_get_addr: can't find address for host cpu064, check slurm.conf
srun: error: Task launch for StepId=11154673.0 failed on node cpu064: Can't find an address, check slurm.conf
srun: error: Application launch failed: Can't find an address, check slurm.conf
srun: Job step aborted
```
Can the admins please have a look? Thanks!


## Post 2 by @Adrien.Albert (2024-07-03T17:42:47.771Z)

Hi,
I am on it:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/18) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Baobab: DNS issue
Dear users, 
We are currently experiencing DNS issue on Baobab, some of you have already encountered error messages as a result: 
srun: error: _fwd_tree_get_addr: can't find address for host cpu002, check slurm.conf
srun: error: Task launch for StepId=11157586.0 failed on node cpu002: Can't find an address, check slurm.conf
srun: error: Application launch failed: Can't find an address, check slurm.conf
srun: Job step aborted

We are working on it to resolve it as soon as possib…


## Post 3 by @Bharathkumar.Radhakrishnan (2024-07-03T17:44:58.534Z)

@Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert) thanks!
