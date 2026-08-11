# Unable to run ug_slurm_usage_per_user.py script

- Source: https://hpc-community.unige.ch/t/3863

- Created: 2025-03-13T17:12:43.757Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Pablo.Strasser (2025-03-13T17:12:43.790Z)

## Primary informations
Username: strassep
Cluster: Baobab
## Description
```
(baobab)-[strassep@login1 ~]$ ug_slurm_usage_per_user.py --start 2025-01-01 --pi kalousis
Error: Failed to execute the sreport command. Command '['sreport', '--all_clusters', '-t', 'Hours', '--tres=billing', 'Cluster', 'UserUtilizationByAccount', 'users=strassep', 'Accounts=kalousis', 'start=2025-01-01', 'end=2025-03-13T18:08:44', 'Fo
rmat=Cluster,Login,Proper,Account,TresName,Used']' returned non-zero exit status 1.
stderr: sreport: error: slurm_persist_conn_open: No response to persist_init
sreport: error: Sending PersistInit msg: Resource temporarily unavailable
sreport: fatal: Problem connecting to the database: Resource temporarily unavailable
```
This seem to be a problem with the slurm controller.
Have a nice evening.


## Post 2 by @Adrien.Albert (2025-03-13T23:32:29.881Z)

Pablo.Strasser:
> sreport: fatal: Problem connecting to the database: Resource temporarily unavailable
Issue related to:
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/5) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
We encountered an issue with Slurm, the database is unreachable resulting errors executing slurm command (sinfo, sacct etc…) 
Resolution
Database service has been restarted. 
Status : Resolved green_circle
start: 2025-03-12T15:00:00Z (UTC) 
end: 2025-03-13T22:30:00Z (UTC)
