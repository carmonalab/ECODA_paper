# Issue with RStudio on bamboo OOD

- Source: https://hpc-community.unige.ch/t/4235

- Created: 2026-03-03T10:13:09.500Z

- Tags: bamboo, openondemand

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2026-03-03T10:13:09.563Z)

Hi,
Since this morning I cannot use RStudio in OOD on bamboo:
The job is scheduled, the files are created:
```
$ ls -alh /home/users/d/delislel/ondemand/data/sys/dashboard/batch_connect/sys/ug_RStudio/output/79bf1c0c-2ccc-419a-918f-d31fb4edf3e1 
total 12K
drwxr-xr-x   3 delislel hpc_users    6 Mar  3 11:00 .
drwx------ 390 delislel unige      388 Mar  3 11:08 ..
-rwxr-xr-x   1 delislel hpc_users  216 Mar  3 11:00 before.sh
drwxr-xr-x   2 delislel hpc_users    2 Mar  3 11:00 etc
-rw-r--r--   1 delislel hpc_users 4.9K Mar  3 11:00 job_script_content.sh
-rw-r--r--   1 delislel hpc_users  707 Mar  3 11:00 job_script_options.json
-rwxr-xr-x   1 delislel hpc_users 3.2K Mar  3 11:00 script.sh
-rw-r--r--   1 delislel hpc_users  359 Mar  3 11:00 user_defined_context.json
```
But then it is requeued:
```
$ scontrol show jobid=3523703
JobId=3523703 JobName=sys/dashboard/sys/ug_RStudio
UserId=delislel(313457) GroupId=hpc_users(5000) MCS_label=N/A
Priority=0 Nice=0 Account=herrerap QOS=normal
JobState=PENDING Reason=user_env_retrieval_failed_requeued_held Dependency=(null)
```
Thanks,
Lucille


## Post 2 by @Lucille.Delisle1 (2026-03-03T12:00:45.797Z)

It seems to be solved. I don’t know if you did something.
