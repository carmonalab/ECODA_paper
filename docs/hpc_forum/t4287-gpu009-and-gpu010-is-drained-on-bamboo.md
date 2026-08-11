# Gpu009 and gpu010 is drained on bamboo

- Source: https://hpc-community.unige.ch/t/4287

- Created: 2026-04-29T11:03:20.785Z

- Tags: bamboo

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @pablo.strasser1 (2026-04-29T11:03:20.879Z)

Hi,
gpu009 and gpu010 are drained on bamboo
(bamboo)-[strassep@login1 compression-binary-tree-transformer]$ scontrol show node gpu010
NodeName=gpu010 Arch=x86_64 CoresPerSocket=24
CPUAlloc=0 CPUEfctv=94 CPUTot=96 CPULoad=0.01
AvailableFeatures=EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKWELL,nvidia_geforce_rtx_5090
ActiveFeatures=EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKWELL,nvidia_geforce_rtx_5090
Gres=gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G
NodeAddr=gpu010 NodeHostName=gpu010 Version=25.11.3
OS=Linux 5.14.0-570.58.1.el9_6.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Oct 31 13:55:05 UTC 2025
RealMemory=384000 AllocMem=0 FreeMem=155290 Sockets=4 Boards=1
CoreSpecCount=2 CPUSpecList=71,95
State=IDLE+DRAIN ThreadsPerCore=1 TmpDisk=7123000 Weight=10 Owner=N/A MCS_label=N/A
Partitions=private-kalousis-gpu,shared-gpu
BootTime=2026-02-26T09:43:05 SlurmdStartTime=2026-03-31T13:59:26
LastBusyTime=2026-04-29T12:48:47 ResumeAfterTime=None
CfgTRES=cpu=94,mem=375G,billing=231,gres/gpu=4,gres/gpu:nvidia_geforce_rtx_5090=4
AllocTRES=
CurrentWatts=n/a AveWatts=n/a
Reason=health_ps___blocked [root@2026-04-29T03:06:18]
(bamboo)-[strassep@login1 compression-binary-tree-transformer]$ scontrol show node gpu009
NodeName=gpu009 Arch=x86_64 CoresPerSocket=24
CPUAlloc=0 CPUEfctv=94 CPUTot=96 CPULoad=0.01
AvailableFeatures=EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKWELL,nvidia_geforce_rtx_5090
ActiveFeatures=EPYC-9654,V12,COMPUTE_CAPABILITY_12_0,COMPUTE_TYPE_BLACKWELL,nvidia_geforce_rtx_5090
Gres=gpu:nvidia_geforce_rtx_5090:4(S:0-3),VramPerGpu:no_consume:32G
NodeAddr=gpu009 NodeHostName=gpu009 Version=25.11.3
OS=Linux 5.14.0-570.58.1.el9_6.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Oct 31 13:55:05 UTC 2025
RealMemory=384000 AllocMem=0 FreeMem=185260 Sockets=4 Boards=1
CoreSpecCount=2 CPUSpecList=71,95
State=IDLE+DRAIN ThreadsPerCore=1 TmpDisk=7123000 Weight=10 Owner=N/A MCS_label=N/A
Partitions=private-kalousis-gpu,shared-gpu
BootTime=2026-02-26T09:43:06 SlurmdStartTime=2026-03-31T13:59:24
LastBusyTime=2026-04-29T12:48:47 ResumeAfterTime=None
CfgTRES=cpu=94,mem=375G,billing=231,gres/gpu=4,gres/gpu:nvidia_geforce_rtx_5090=4
AllocTRES=
CurrentWatts=n/a AveWatts=n/a
Reason=health_ps___blocked [root@2026-04-29T03:00:18]


## Post 2 by @Gael.Rossignol (2026-05-07T10:58:11.929Z)

Dear Pablo,
The servers have been back online for a while now. I hope everything is fine on your side after this issue.
Best regards,


## Post 3 by @pablo.strasser1 (2026-05-07T11:56:58.892Z)

Yes it is solved.I think it was solved the same day I repported it.
