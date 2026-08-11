# No trace of submitted SLURM JOB

- Source: https://hpc-community.unige.ch/t/4241

- Created: 2026-03-10T15:00:21.958Z

- Posts: 10

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Marco.Froelich (2026-03-10T15:00:22.013Z)

## Primary informations
Username: mfroelich
Cluster: bamboo
## Description
I am submitting cpu-jobs to SLURM on Bamboo, but no output files are being created (.err and .out) and no email is being sent to my email. Squeue shows that the job is running. Nothing was happening (i.e. code did not seem to be running) so I cancel the job, and now squeue has been showing it with State CG for over 10 minutes.
I tested this again with a dummy sbatch script (see below), and the same thing happened.
I noticed as well that there other CG state jobs in squeue.
image
image382×453 9.65 KB
[image382×453 9.65 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/0326ced981fcd612412df8d7418b5abdaf889410.png)


## Post 2 by @Lucille.Delisle1 (2026-03-10T15:10:32.086Z)

Hi,
I got the same issue. If I look at the job it is in this type of state: `JobState=RUNNING Reason=Prolog`
See an example here, it can stay in Prolog for more than 20 minutes…
```
$ scontrol show jobid=3640369
JobId=3640369 JobName=Multiome_R_step6_compile
   UserId=delislel(313457) GroupId=hpc_users(5000) MCS_label=N/A
   Priority=1575015 Nice=0 Account=andreygu QOS=normal
   JobState=RUNNING Reason=Prolog Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:22:19 TimeLimit=12:00:00 TimeMin=N/A
   SubmitTime=2026-03-08T21:57:17 EligibleTime=2026-03-08T21:57:17
   AccrueTime=2026-03-08T21:57:17
   StartTime=2026-03-08T22:04:20 EndTime=2026-03-09T10:04:20 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2026-03-08T22:04:20 Scheduler=Main
   Partition=shared-cpu AllocNode:Sid=login1:1351239
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=cpu041
   BatchHost=cpu041
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   ReqTRES=cpu=1,mem=120G,node=1,billing=31
   AllocTRES=cpu=1,mem=120G,node=1,billing=31
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=120G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) LicensesAlloc=(null) Network=(null)
   Command=/home/users/d/delislel/scripts/limb-multiome-chudzik/cellranger_all_samples/06_compile_info_and_select_cells.sh
   SubmitLine=sbatch cellranger_all_samples/06_compile_info_and_select_cells.sh
   WorkDir=/home/users/d/delislel/scripts/limb-multiome-chudzik/
   StdErr=/home/users/d/delislel/scripts/limb-multiome-chudzik//cellranger_all_samples/logs_post_CRA/slurm-Multiome_R_step6_compile-3640369.err
   StdIn=/dev/null
   StdOut=/home/users/d/delislel/scripts/limb-multiome-chudzik//cellranger_all_samples/logs_post_CRA/slurm-Multiome_R_step6_compile-3640369.out
   TresPerTask=cpu=1
```


## Post 3 by @Lucille.Delisle1 (2026-03-10T15:11:59.229Z)

This was 2 days ago but today I still have a lot of issues like this.


## Post 4 by @Daniel.Gilfraga (2026-03-10T15:13:12.671Z)

Hi,
I am experiencing the same issue. When I submit CPU jobs to SLURM on Bamboo, they appear to start, but no `.out` or `.err` files are created.
When I check with `squeue --me`, the jobs appear to be running but nothing happens. When I cancel them, they remain listed in squeue –me with the same elapsed time.


## Post 5 by @Lucille.Delisle1 (2026-03-10T15:18:11.246Z)

If you use `squeue --me -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %.3C %.12r"` you will see the reason.


## Post 6 by @Marco.Froelich (2026-03-10T16:24:51.014Z)

UPDATE: I recieved emails of my jobs all at once about 2h after I cancelled them. Output files did not get written. My test script now runs and I get an email notification, but no file is written and it seems to crash immediately.


## Post 7 by @Marco.Froelich (2026-03-11T07:04:21.728Z)

UPDATE: Seems to be resolved. Jobs are running, are producing err and out files and email notifications are being sent.


## Post 8 by @Yann.Sagon (2026-03-11T11:58:01.986Z)

Hello,
this is probably related with this issue: [2026] Current issues on HPC Cluster - #8 by Yann.Sagon[[2026] Current issues on HPC Cluster - #8 by Yann.Sagon](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/8)
We have re opened the case at schedmd today.


## Post 9 by @Lucille.Delisle1 (2026-03-11T12:47:22.911Z)

I don’t think this is related to this issue because on my side it was not interactive jobs. I can give you the details if you need them. But I noticed that there were a lot of jobs running from the same user at this time probably launched at the same time so I think it is more linked to the number of request sent to the scheduler.
Best


## Post 10 by @Yann.Sagon (2026-03-11T13:45:14.248Z)

Still, it is probably related: in the background, the issue seems to be related with network communication issue between compute nodes and slurm controller. The trigger is indeed a high number of jobs running, but this isn’t a normal behavior. We tried to increase “some numbers” on Bamboo.
Please ping us here if you still have the issue on Bamboo. If the issue is solved, we’ll propagate the configuration to the other clusters.
