# Feature request : Allow the request for a tmp volume in SLURM

- Source: https://hpc-community.unige.ch/t/3591

- Created: 2024-08-12T10:54:59.383Z

- Posts: 4

- Category: 10

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maxime.Juventin (2024-08-12T10:54:59.420Z)

# Clusters : Yggdrasil, Bamboo
Hello,
For some intensive I/O computations, I have tried (on Yggdrasil and Bamboo) to copy my dataset at run time in the /tmp of my compute nodes, in order to avoid scratch I/O saturation[scratch I/O saturation](https://hpc-community.unige.ch/t/scratch-disk-i-o-saturation/3201). Despite the added copy time, the job is much faster by bypassing the reads and writes on the scratch. The only potential issue I foresee is the risk of saturating the /tmp if too many jobs are assigned to the same node.
For this reason, I would like to ask if you can add an option in the srun/sbatch to request an amount of tmp size ?
Thank you,
Cheers


## Post 2 by @Adrien.Albert (2024-08-12T13:26:45.436Z)

Hi @Maxime.Juventin[@Maxime.Juventin](https://hpc-community.unige.ch/u/maxime.juventin),
I implemented a local scratch share system about a year ago:
Local Share Directory Between Jobs on Compute[Local Share Directory Between Jobs on Compute](https://hpc-community.unige.ch/t/local-share-directory-beetween-jobs-on-compute/2893)
This shared space is cleaned after the last job on the node.
To find local storage information using Slurm, you can use the `scontrol` command:
```
(bamboo)-[alberta@login1 ~]$ scontrol show node cpu001
NodeName=cpu001 Arch=x86_64 CoresPerSocket=64 
   CPUAlloc=0 CPUEfctv=126 CPUTot=128 CPULoad=0.08
   AvailableFeatures=EPYC-7742,V8,TmpDisk800G
   ActiveFeatures=EPYC-7742,V8,TmpDisk800G
   Gres=(null)
   NodeAddr=cpu001 NodeHostName=cpu001 Version=23.11.7
   OS=Linux 4.18.0-513.24.1.el8_9.x86_64 #1 SMP Thu Apr 4 18:13:02 UTC 2024 
   RealMemory=512000 AllocMem=0 FreeMem=510129 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=63,127 
   State=IDLE ThreadsPerCore=1 TmpDisk=800000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=debug-cpu 
   BootTime=2024-08-08T15:13:22 SlurmdStartTime=2024-08-12T15:07:44
   LastBusyTime=2024-08-12T15:13:06 ResumeAfterTime=None
   CfgTRES=cpu=126,mem=500G,billing=126
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/a ExtSensorsWatts=0 ExtSensorsTemp=n/a
```
You can target nodes with a minimum of `XXXGB` of storage, but keep in mind that this does not guarantee that the node will have sufficient storage available, as the storage is shared among all running jobs.
On our end, it might be useful to configure a new constraint to specify node types with `XXXGB` of storage.
For example (I just test it), we can define a `TmpDiskXXXG` feature (5th line in the output):
```
(bamboo)-[alberta@login1 ~]$ scontrol show node cpu001
NodeName=cpu001 Arch=x86_64 CoresPerSocket=64 
   CPUAlloc=0 CPUEfctv=126 CPUTot=128 CPULoad=0.08
   AvailableFeatures=EPYC-7742,V8,TmpDisk800G
   ActiveFeatures=EPYC-7742,V8,TmpDisk800G
   Gres=(null)
   NodeAddr=cpu001 NodeHostName=cpu001 Version=23.11.7
   OS=Linux 4.18.0-513.24.1.el8_9.x86_64 #1 SMP Thu Apr 4 18:13:02 UTC 2024 
   RealMemory=512000 AllocMem=0 FreeMem=510129 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=63,127 
   State=IDLE ThreadsPerCore=1 TmpDisk=800000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=debug-cpu 
   BootTime=2024-08-08T15:13:22 SlurmdStartTime=2024-08-12T15:07:44
   LastBusyTime=2024-08-12T15:13:06 ResumeAfterTime=None
   CfgTRES=cpu=126,mem=500G,billing=126
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/a ExtSensorsWatts=0 ExtSensorsTemp=n/a
```
In this example, the constraint `ActiveFeatures=TmpDisk800G` indicates that the maximum local disk available is 800GB.
To start a job with this constraint:
```
(bamboo)-[alberta@login1 ~]$ srun --constraint=TmpDisk800G hostname
srun: job 53734 queued and waiting for resources
srun: job 53734 has been allocated resources
cpu001.bamboo
```
I’ll keep this idea in mind for our next team meeting and explore the potential benefits (or disadvantage) of it.


## Post 3 by @Maxime.Juventin (2024-08-12T16:14:45.013Z)

Adrien.Albert:
> You can target nodes with a minimum of `XXXGB` of storage, but keep in mind that this does not guarantee that the node will have sufficient storage available, as the storage is shared among all running jobs.
Ok, it is good to know that we can request a maximum local volume.
As SLURM can’t properly handle what the availability of the local storage, I understand that it is not possible to prevent several users to saturate the local scratch.
On the other, to avoid that my own jobs saturate the local scratch, I just adapted my batch to avoid having multiple jobs on the same node with `--ntask-per-node=4 --nodes=1` , it seems to work on Bamboo.


## Post 4 by @Adrien.Albert (2024-08-15T09:22:49.702Z)

Hi @Maxime.Juventin[@Maxime.Juventin](https://hpc-community.unige.ch/u/maxime.juventin)
My bad, the filter already exists:
The option `--tmp`[--tmp](https://slurm.schedmd.com/srun.html#OPT_tmp)  seems the most appropriate:
```
# Request Node with at least 500G storage
(baobab)-[alberta@login1 ~]$ srun --tmp 500G -p shared-cpu hostname
srun: job 11965317 queued and waiting for resources
srun: job 11965317 has been allocated resources
cpu329.baobab

# Check my job -> MinTmpDiskNode=500G
(baobab)-[alberta@login1 ~]$ scontrol show job 11965317
JobId=11965317 JobName=hostname
   UserId=alberta(401775) GroupId=hpc_users(5000) MCS_label=N/A
   Priority=1575002 Nice=0 Account=burgi QOS=normal
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=0 Reboot=0 ExitCode=0:0
   RunTime=00:00:01 TimeLimit=00:01:00 TimeMin=N/A
   SubmitTime=2024-08-15T11:05:45 EligibleTime=2024-08-15T11:05:45
   AccrueTime=2024-08-15T11:05:45
   StartTime=2024-08-15T11:05:57 EndTime=2024-08-15T11:05:58 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2024-08-15T11:05:57 Scheduler=Backfill
   Partition=shared-cpu AllocNode:Sid=login1:1013693
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=cpu329
   BatchHost=cpu329
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   ReqTRES=cpu=1,mem=3000M,node=1,billing=1
   AllocTRES=cpu=1,mem=3000M,node=1,billing=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=3000M MinTmpDiskNode=500G
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=hostname
   WorkDir=/home/users/a/alberta
   Power=
   
#  check the node -> TmpDisk=800000
(baobab)-[alberta@login1 ~]$ scontrol show node cpu329
NodeName=cpu329 Arch=x86_64 CoresPerSocket=64 
   CPUAlloc=122 CPUEfctv=126 CPUTot=128 CPULoad=20.11
   AvailableFeatures=EPYC-7742,V8
   ActiveFeatures=EPYC-7742,V8
   Gres=(null)
   NodeAddr=cpu329 NodeHostName=cpu329 Version=23.11.7
   OS=Linux 4.18.0-513.24.1.el8_9.x86_64 #1 SMP Thu Apr 4 18:13:02 UTC 2024 
   RealMemory=512000 AllocMem=490120 FreeMem=478682 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=63,127 
   State=MIXED ThreadsPerCore=1 TmpDisk=800000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=shared-cpu 
   BootTime=2024-07-19T11:10:48 SlurmdStartTime=2024-07-30T13:59:38
   LastBusyTime=2024-08-12T23:51:53 ResumeAfterTime=None
   CfgTRES=cpu=126,mem=500G,billing=126
   AllocTRES=cpu=122,mem=490120M
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/a ExtSensorsWatts=0 ExtSensorsTemp=n/a
```
But the following statement still applies:
Adrien.Albert:
> You can target nodes with a minimum of `XXXGB` of storage, but keep in mind that this does not guarantee that the node will have sufficient storage available, as the storage is shared among all running jobs.
