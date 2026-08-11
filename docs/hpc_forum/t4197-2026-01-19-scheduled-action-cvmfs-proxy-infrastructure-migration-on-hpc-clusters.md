# [2026-01-19] Scheduled action: CVMFS Proxy Infrastructure Migration on HPC Clusters

- Source: https://hpc-community.unige.ch/t/4197

- Created: 2026-01-16T21:05:04.977Z

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2026-01-16T21:05:05.061Z)

Dear users,
We would like to inform you of an upcoming migration of the CVMFS Squid caching infrastructure across our HPC clusters. This update aims to improve performance, reliability, and cache consistency.
### Migration schedule
:date: Monday, 19 January 2026
:ten_thirty: Starting at 10:30
Order of clusters:
- Bamboo
- Yggdrasil
- Baobab
---
### Current setup
- A single Squid CVMFS server shared by the three clusters
- 30 GB proxy global cache
- 30 GB local cache on compute nodes
Example of current limitation:
> Reading multiple files totaling more than 20 GB forces both the compute node cache and the proxy cache to be flushed after each read. This significantly impacts performance for jobs relying on CVMFS.
> Combined with usage from multiple users, jobs, and three clusters sharing the same cache, the performance degradation can become substantial.
---
### What will change
- Each cluster will have its own dedicated Squid proxy for CVMFS
- The proxy global cache will be expanded to 100 GB
- Local cache on compute nodes will also be increased to 100 GB
---
### Expected benefits
These enhancements will:
- Improve CVMFS file read performance
- Greatly reduce the need for cache regeneration
- Provide cache autonomy and stability for each cluster independently
Our internal tests have confirmed that this migration can be performed live in production without any impact to running jobs.
---
We will keep you informed as the migration progresses.
If you have any questions, please feel free to reach out.


## Post 2 by @Adrien.Albert (2026-01-19T10:45:15.378Z)

Dear HPC users;
2026-01-19T10:44:00Z The CVMFS proxy  on Bamboo has been successfully migrated to the new one.
2026-01-19T10:55:00Z The CVMFS proxy  on Yggdrasil has been successfully migrated to the new one.
2026-01-19T11:10:00Z The CVMFS proxy  on Baobab has been successfully migrated to the new one.
