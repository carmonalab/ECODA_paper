# Nodes on ressources for private-gap-cpu

- Source: https://hpc-community.unige.ch/t/4037

- Created: 2025-08-08T13:47:10.477Z

- Tags: baobab

- Posts: 3

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Laure.Moinat (2025-08-08T13:47:10.527Z)

Username: moinatl
partition: private-gap-cpu
I don’t understand why the simulation are put on ‘Ressources’, because there is a lot of nodes available for 'private-gap-cpu.
What is the reason,
Best,
Laure


## Post 2 by @Adrien.Albert (2025-08-08T14:12:10.887Z)

Hi @Laure.Moinat[@Laure.Moinat](https://hpc-community.unige.ch/u/laure.moinat)
Could give more information:
Could you give the sbatch? The output of:
` scontrol show job <jobid>`
if the nodes is not in drain status, it requires times to slurm to recalculate the priority of your allocation. Even if a node is available it does not mean the requested resources are available on the node.
Requesting 32 cpu but only there are only 16 cpu  on the node. Slurm can not allocate the job and will wait other available resources.


## Post 3 by @Laure.Moinat (2025-08-08T14:22:05.147Z)

Here is the description of the situation :
```
  JobId=1965475 JobName=pres1.conf
   UserId=moinatl(401279) GroupId=hpc_users(5000) MCS_label=N/A
   Priority=1500050 Nice=0 Account=kasparia QOS=normal
   JobState=PENDING Reason=Resources Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=0 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=7-00:00:00 TimeMin=N/A
   SubmitTime=2025-08-08T15:58:06 EligibleTime=2025-08-08T15:58:06
   AccrueTime=2025-08-08T15:58:06
   StartTime=2025-08-09T01:43:39 EndTime=2025-08-16T01:43:39 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2025-08-08T16:14:58 Scheduler=Main
   Partition=private-gap-cpu AllocNode:Sid=login1:1982338
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList= SchedNodeList=cpu[222-223]
   NumNodes=1 NumCPUs=25 NumTasks=25 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   ReqTRES=cpu=25,mem=75000M,node=1,billing=43
   AllocTRES=(null)
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=3000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=pres1.conf
   WorkDir=/srv/beegfs/scratch/users/m/moinatl/Present1
```
Normally the partition has a lot a cpu available and there is only two job running so I am a bit confuse why the simulations are not running.
Thanks for your response.
