# Slurm db unavailable

- Source: https://hpc-community.unige.ch/t/4096

- Created: 2025-09-11T10:48:01.504Z

- Tags: slurm

- Posts: 3

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cyrus.Brueggimann (2025-09-11T10:48:01.504Z)

Good afternoon. I know a maintenance was started recently, but it has been reported to be done and I do happen to have the exact same issue as the OP. The account selector appears for most interactive session tabs in the baobab and yggdrasil cluster, but is empty with no options available.
image
image794×108 2.48 KB
[image794×108 2.48 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/c146fbcbed4be0fe20430afa7b9875278409790b.png)
Clicking on the “Launch” button greets me with a :
Failed to submit session with the following error:
```
sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
```
> - If this job failed to submit because of an invalid job name please ask your administrator to configure OnDemand to set the environment variable OOD_JOB_NAME_ILLEGAL_CHARS.
> - The Jupyter Lab session data for this session can be accessed under the staged root directory[staged root directory](https://openondemand.baobab.hpc.unige.ch/pun/sys/dashboard/files/fs/home/users/b/bruggim9/ondemand/data/sys/dashboard/batch_connect/sys/ug_JupyterLab/output/e9a1f202-9c67-443d-a3ed-5e930472f637).
The staged root directories for the failed sessions do not contain any output logs or any help regarding this issue.
I also happen to have the same issues as the OP regarding the login node :
`$ sacctmgr show user where name=bruggim9 withassoc`
sacctmgr: error: _open_persist_conn: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Connection refused
sacctmgr: error: Sending PersistInit msg: Connection refused
I still can create jobs manually and VS code sessions, but the Jupyter lab sessions, Advanced desktop sessions and R studio sessions are all not working for me. I never had this issue before (tried last week and everything was fine).
Many thanks for your guidance


## Post 2 by @Paul.Coppin (2025-09-12T06:52:03.574Z)

Hi,I also get errors trying to view the resource usage of jobs on Baobab:
```
[coppinp@login1 ~]$ seff 2428834_0
perl: error: _open_persist_conn: failed to open persistent connection to host:lunihpcslurm1.admin.unige.ch:6819: Connection refused
perl: error: Sending PersistInit msg: Connection refused
perl: error: Sending PersistInit msg: Connection refused
perl: error: DBD_GET_JOBS_COND failure: Unspecified error
Job not found.
Segmentation fault (core dumped)
```
If this is because the maintenance is not fully over yet as stated in the email, please ignore.


## Post 3 by @Gael.Rossignol (2025-09-22T13:23:05.322Z)

Dear Cyrus and Paul,
During last maintenance we had an issue, and the slurm db has been updated. So during update other cluster works fine but some requests was not possible.
Issue is now resolved. Sorry for inconvenience,
Best regards,
