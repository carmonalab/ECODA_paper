# Private-dpnc-cpu partition reserved for maintenance

- Source: https://hpc-community.unige.ch/t/3802

- Created: 2025-01-29T10:32:03.684Z

- Tags: baobab, slurm

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2025-01-29T10:32:03.731Z)

Dear HPC admins,
I would like to better understand the reason my jobs have been queueing more than usual (some of them since yesterday).
The message I get is that the private-dpnc-cpu partition is reserved for maintenance and my jobs will only start on Saturday. The full log is attached below. Is there maintenance ongoing that I missed perhaps?
Thank you very much in advance.
Best wishes,
Vilius
> ➜ ~ scontrol show job 14474229
> JobId=14474229 JobName=97bc53fd-38f8-4e5a-8419-f7e0af9134ab
> UserId=cepaitis(404575) GroupId=private_dpnc(1014) MCS_label=N/A
> Priority=2703626 Nice=0 Account=schramm QOS=normal
> JobState=PENDING Reason=ReqNodeNotAvail,_Reserved_for_maintenance Dependency=(null)
> Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
> RunTime=00:00:00 TimeLimit=03:00:00 TimeMin=N/A
> SubmitTime=2025-01-28T15:05:25 EligibleTime=2025-01-28T15:05:25
> AccrueTime=2025-01-28T15:05:25
> StartTime=2025-02-01T01:43:20 EndTime=2025-02-01T04:43:20 Deadline=N/A
> SuspendTime=None SecsPreSuspend=0 LastSchedEval=2025-01-29T11:23:18 Scheduler=Main
> Partition=private-dpnc-cpu AllocNode:Sid=login1:3343805
> ReqNodeList=(null) ExcNodeList=(null)
> NodeList= SchedNodeList=cpu086
> NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0::
> ReqTRES=cpu=1,mem=5000M,node=1,billing=2
> AllocTRES=(null)
> Socks/Node=* NtasksPerN:B:S:C=0:0:: CoreSpec=*
> MinCPUsNode=1 MinMemoryNode=5000M MinTmpDiskNode=0
> Features=(null) DelayBoot=00:00:00
> OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
> Command=(null)
> WorkDir=/home/users/c/cepaitis/Pileup/weakly-supervised-search
> Comment=rule_prepare_data
> StdErr=/home/users/c/cepaitis/Pileup/weakly-supervised-search/.snakemake/slurm_logs/rule_prepare_data/14474229.log
> StdIn=/dev/null
> StdOut=/home/users/c/cepaitis/Pileup/weakly-supervised-search/.snakemake/slurm_logs/rule_prepare_data/14474229.log
> TresPerTask=cpu=1


## Post 2 by @Mario.AlvesCardoso (2025-01-29T10:37:13.649Z)

I’ve been experiencing the same issue! It also seems that most of the CPUs are being used by 1or 2 person/people.


## Post 3 by @Yann.Sagon (2025-01-30T10:57:04.075Z)

Dear members of the DPNC group, this is related to this post Unannounced maintenance on baobab? - #4 by Yann.Sagon[Unannounced maintenance on baobab? - #4 by Yann.Sagon](https://hpc-community.unige.ch/t/unannounced-maintenance-on-baobab/3797/4)
My advice is that DPNC members take some time to discuss internally why this reservation has been set up and why not all members are aware of it.
Best


## Post 4 by @Vilius.Cepaitis (2025-01-30T12:18:12.662Z)

Thank you Yann, I’ve raised this with the DPNC
