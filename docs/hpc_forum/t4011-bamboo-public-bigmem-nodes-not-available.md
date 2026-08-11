# Bamboo public-bigmem nodes not available

- Source: https://hpc-community.unige.ch/t/4011

- Created: 2025-07-22T07:50:55.431Z

- Tags: bamboo

- Posts: 23

- Category: 8

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-07-22T07:50:55.474Z)

Hi,
I’ve been trying to launch a job on public-bigmem, but I get the following message: “Nodes required for job are DOWN, DRAINED or reserved for jobs in higher priority partitions” and my job is on hold, although when I write down `squeue -p public-bigmem` I see that nobody else is currently using that partition and more generally when typing `squeue` I don’t see so many nodes currently being used so I don’t understand why this happens and how long I will need to wait for my job to launch? Is that a normal behavior (like competition between public-bigmem and shared-bigmem with a higher priority by default to shared-bigmem possibly also explaining why sometimes the scripts on bigmem crash due to memory issues) or is there an issue currently with public-bigmem? I need to run a job lasting between 12h and 24h so cannot switch to shared-bigmem and need results by next week ideally so I hope the nodes are functional again soonish… :grimacing: :folded_hands:


## Post 2 by @nicolas.clairis (2025-07-22T14:30:45.912Z)

My scripts are still on hold since this morning… Any way to know how long this is going to last and/or which are the other partitions sharing the nodes with public-bigmem? Based on the squeue command there are like only 10 nodes currently engaged so I don’t understand why none is available for public-bigmem


## Post 3 by @Yann.Sagon (2025-07-22T14:56:10.786Z)

Dear @nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
I’ve checked and the issue is that a user launched a job array requesting 600 GB per items. The compute nodes on bigmem partition  has 1 TB of RAM, so I’m not sure why your job couldn’t be run, as this won’t prevent the other array items from being run — only one can be run at a time.
The good news is that most of the compute nodes have 512 GB or more, and I have updated the partition to use public-cpu. Your job is now running.
Best regards


## Post 4 by @Yann.Sagon (2025-07-22T14:58:01.425Z)

nicolas.clairis:
> Based on the squeue command there are like only 10 nodes currently engaged so I don’t understand why none is available for public-bigmem
The nodes associated with public-bigmem are statically defined: we have two nodes considered “bigmem” in Bamboo (1TB RAM).


## Post 5 by @nicolas.clairis (2025-07-22T15:13:29.704Z)

Thanks!! :folded_hands:
So you mean the “shared-bigmem” and “public-bigmem” use the same nodes? Because when I look at “public-bigmem” there is no other job running which is why I was a bit suprised that mine could not launch, but there is indeed a queue of like a 1000 jobs from 1 user on shared-bigmem which could explain why other jobs don’t manage to run.
Also, the fact that you swapped my job to public-cpu, won’t that make it crash if it requires high levels of RAM?


## Post 6 by @Yann.Sagon (2025-07-22T15:20:34.621Z)

nicolas.clairis:
> So you mean the “shared-bigmem” and “public-bigmem” use the same nodes?
shared-bigmem includes every “bigmem” nodes public and private if any.
nicolas.clairis:
> Because when I look at “public-bigmem” there is no other job running
Both nodes from public-bigmem had running jobs on them. How did you checked?
nicolas.clairis:
> Also, the fact that you swapped my job to public-cpu, won’t that make it crash if it requires high levels of RAM?
No, the memory you requested: 250GB will be reserved on public-cpu node without issue.


## Post 7 by @nicolas.clairis (2025-07-22T15:46:57.231Z)

Ok thanks for the clarifications!
Regarding “Both nodes from public-bigmem had running jobs on them. How did you checked?” I just typed `squeue -p public-bigmem` and when I do, only my jobs appear… while if I do the same on shared-bigmem I then see the list of ~1000 jobs which is currently using cpu044 and cpu045


## Post 8 by @nicolas.clairis (2025-07-23T07:38:37.568Z)

what’s the upper RAM limit for public-cpu? I see that my second job on public-bigmem launched yesterday is still on hold because the user using shared-bigmem is still only at job 600/1529 so I guess that public-bigmem will not be available until these jobs are done…
I couldn’t find the info about the RAM capacity for the different partitions in Bamboo on the doc hpc:slurm [eResearch Doc][hpc:slurm [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/slurm#cpu) and would be helpful to know it to know when it’s best to swap. The other job I have to launch would require around 500G of RAM so I don’t know if public-cpu is also appropriate or not


## Post 9 by @nicolas.clairis (2025-07-23T14:59:01.963Z)

Can you clarify the point regarding the amount of RAM on public-cpu please? Because I’m kind of stuck here. I tried to launch other jobs that may run in less than 12h on shared-bigmem but they are all on hold as well… I will try public-cpu but I’m afraid that it cannot go up to the amount of RAM I need for those jobs (at least 500Go)


## Post 10 by @Nicolas.Clairis1 (2025-07-24T10:15:03.504Z)

Ok thanks to chatgpt I managed to find that the public-cpu nodes are roughly limited to 500Go of RAM. One of my jobs crashed with the “out of memory” error message. I’m going to retry to launch it in public-cpu but ideally I would need more RAM and currently the bigmem partitions are still blocked by the same user who is still somewhere around job 640/1529 so I’d rather not wait 3 weeks until their jobs are done… Isn’t there a way to share resources so that when one job is done the same user does not always get priority on reusing that same node to be a bit more fair?
For now my jobs are running, but I’m afraid that they will probably crash again due to lack of RAM


## Post 11 by @Gael.Rossignol (2025-07-24T13:52:52.981Z)

Dear Nicolas,
You can use partition “shared-bigmem” on cluster Yggdrasil, nodes are available :
```
(yggdrasil)-[root@login1 ~]$ scontrol show nodes cpu[113-115,120-122]
NodeName=cpu113 Arch=x86_64 CoresPerSocket=8
   CPUAlloc=1 CPUEfctv=14 CPUTot=16 CPULoad=6.90
   AvailableFeatures=GOLD-6244,XEON_GOLD_6244,V9
   ActiveFeatures=GOLD-6244,XEON_GOLD_6244,V9
   Gres=(null)
   NodeAddr=cpu113 NodeHostName=cpu113 Version=24.11.1
   OS=Linux 5.14.0-503.14.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Nov 15 12:04:32 UTC 2024
   RealMemory=770000 AllocMem=3000 FreeMem=714574 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=7,15
   State=MIXED ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=public-bigmem,shared-bigmem
   BootTime=2025-05-09T11:39:03 SlurmdStartTime=2025-06-30T10:38:35
   LastBusyTime=2025-07-23T13:36:20 ResumeAfterTime=None
   CfgTRES=cpu=14,mem=770000M,billing=201
   AllocTRES=cpu=1,mem=3000M
   CurrentWatts=0 AveWatts=0

Node cpu114 not found
NodeName=cpu115 Arch=x86_64 CoresPerSocket=8
   CPUAlloc=8 CPUEfctv=14 CPUTot=16 CPULoad=8.00
   AvailableFeatures=GOLD-6244,XEON_GOLD_6244,V9
   ActiveFeatures=GOLD-6244,XEON_GOLD_6244,V9
   Gres=(null)
   NodeAddr=cpu115 NodeHostName=cpu115 Version=24.11.1
   OS=Linux 5.14.0-503.14.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Nov 15 12:04:32 UTC 2024
   RealMemory=770000 AllocMem=500000 FreeMem=508118 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=7,15
   State=MIXED ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=public-bigmem,shared-bigmem
   BootTime=2025-05-09T11:40:04 SlurmdStartTime=2025-06-30T10:38:35
   LastBusyTime=2025-07-23T12:29:07 ResumeAfterTime=None
   CfgTRES=cpu=14,mem=770000M,billing=201
   AllocTRES=cpu=8,mem=500000M
   CurrentWatts=0 AveWatts=0

NodeName=cpu120 Arch=x86_64 CoresPerSocket=18
   CPUAlloc=0 CPUEfctv=34 CPUTot=36 CPULoad=0.64
   AvailableFeatures=GOLD-6240,XEON_GOLD_6240,V9
   ActiveFeatures=GOLD-6240,XEON_GOLD_6240,V9
   Gres=(null)
   NodeAddr=cpu120 NodeHostName=cpu120 Version=24.11.1
   OS=Linux 5.14.0-503.14.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Nov 15 12:04:32 UTC 2024
   RealMemory=1546348 AllocMem=0 FreeMem=1539753 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=17,35
   State=IDLE ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=shared-bigmem,private-wesolowski-bigmem
   BootTime=2025-06-10T10:55:46 SlurmdStartTime=2025-06-30T10:38:42
   LastBusyTime=2025-07-24T15:43:00 ResumeAfterTime=None
   CfgTRES=cpu=34,mem=1546348M,billing=411
   AllocTRES=
   CurrentWatts=0 AveWatts=0

NodeName=cpu121 Arch=x86_64 CoresPerSocket=18
   CPUAlloc=0 CPUEfctv=34 CPUTot=36 CPULoad=0.00
   AvailableFeatures=GOLD-6240,XEON_GOLD_6240,V9
   ActiveFeatures=GOLD-6240,XEON_GOLD_6240,V9
   Gres=(null)
   NodeAddr=cpu121 NodeHostName=cpu121 Version=24.11.1
   OS=Linux 5.14.0-503.14.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Nov 15 12:04:32 UTC 2024
   RealMemory=1546348 AllocMem=0 FreeMem=1540581 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=17,35
   State=IDLE ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=shared-bigmem,private-wesolowski-bigmem
   BootTime=2025-06-10T10:55:47 SlurmdStartTime=2025-06-30T10:38:42
   LastBusyTime=2025-07-24T11:39:12 ResumeAfterTime=None
   CfgTRES=cpu=34,mem=1546348M,billing=411
   AllocTRES=
   CurrentWatts=0 AveWatts=0

NodeName=cpu122 Arch=x86_64 CoresPerSocket=18
   CPUAlloc=0 CPUEfctv=34 CPUTot=36 CPULoad=0.00
   AvailableFeatures=GOLD-6240,XEON_GOLD_6240,V9
   ActiveFeatures=GOLD-6240,XEON_GOLD_6240,V9
   Gres=(null)
   NodeAddr=cpu122 NodeHostName=cpu122 Version=24.11.1
   OS=Linux 5.14.0-503.14.1.el9_5.x86_64 #1 SMP PREEMPT_DYNAMIC Fri Nov 15 12:04:32 UTC 2024
   RealMemory=1546348 AllocMem=0 FreeMem=1540942 Sockets=2 Boards=1
   CoreSpecCount=2 CPUSpecList=17,35
   State=IDLE ThreadsPerCore=1 TmpDisk=150000 Weight=10 Owner=N/A MCS_label=N/A
   Partitions=shared-bigmem,private-wesolowski-bigmem
   BootTime=2025-06-10T10:55:57 SlurmdStartTime=2025-06-30T10:38:42
   LastBusyTime=2025-07-24T11:39:11 ResumeAfterTime=None
   CfgTRES=cpu=34,mem=1546348M,billing=411
   AllocTRES=
   CurrentWatts=0 AveWatts=0
```
cpu120,  cpu121, cpu122 have 1.5TB of ram and not used.
Best regards,


## Post 12 by @Nicolas.Clairis1 (2025-07-24T14:53:40.636Z)

sure thanks for the info but I have already migrated my data from Baobab to Bamboo for the same reason, and it would take multiple days to copy everything from Bamboo to Yggdrasil so I’d rather avoid doing another migration…


## Post 13 by @Yann.Sagon (2025-07-25T06:26:42.104Z)

Dear @Nicolas.Clairis1[@Nicolas.Clairis1](https://hpc-community.unige.ch/u/nicolas.clairis1)
Nicolas.Clairis1:
> it would take multiple days to copy everything from Bamboo to Yggdrasil
I don’t know which method you used to do the copy. To help you, we did a copy of all your scratch data from Bamboo to Yggdrasil this night and it is already done.
This is the command we used:
```
rsync -avH /srv/beegfs/scratch/users/c/clairis/ login1.yggdrasil:/srv/beegfs/scratch/users/c/clairis/copy_bamboo/
```


## Post 14 by @Yann.Sagon (2025-07-25T06:34:50.099Z)

Nicolas.Clairis1:
> Ok thanks to chatgpt I managed to find that the public-cpu nodes are roughly limited to 500Go of RAM.
Good to know chatgpt is more helpful than us:). So, what was the magical command? This may help other users.
Even if public-CPU compute nodes have a lot of RAM, they also have a lot of CPUs. The issue with using them instead of bigmem is that jobs requiring a lot of memory usually don’t need many CPUs. For example, if you request 500 GB and two CPU cores on a standard compute node, all the other CPUs on that node are rendered useless because all the RAM is in use. Bigmem compute nodes usually have faster CPUs with fewer cores. In summary, it is better to use a bigmem compute node if you need a lot of RAM. The issue is determining what “a lot of RAM” means. As many standard compute nodes are idle right now, this isn’t an issue if you use one.


## Post 15 by @nicolas.clairis (2025-07-25T15:32:50.705Z)

Ow thanks a lot for the copy! I got confused because I actually don’t need all the original data for these analyses (which is quite heavy and still stored on Baobab) but the preprocessed data which is on Bamboo and in fact, only part of it, but thanks a lot for the initiative anyway as even copying that still takes several hours in principle. I’ll try to move to Yggdrasil for these heavy analyses then, hoping to not get lost between the 3 partitions :sweat_smile:
Regarding the chatgpt command to know Bamboo cpu RAM limitation, let me know if it’s correct but chatgpt suggested:
`sinfo -p <partition> -N -o '%N %m'` which gives a number in MiB that can then be converted to Go by dividing it by 1024. From that, I get it that all public-cpu nodes are limited to ~500 Go of RAM, while bigmem cpus are roughly limited to 1000Go of RAM. This seemed consistent with empirical limitation when I tried to launch jobs with more than 500Go on public-cpu, but let me know if it’s correct or wrong?


## Post 16 by @Nicolas.Clairis1 (2025-07-28T09:21:15.953Z)

Btw by launching this command, I realize that Yggdrasil public-bigmem can not go higher than 750Go of RAM while Bamboo could go up to 1000 Go, but shared-bigmem has some nodes which can go up to 1500 Go of RAM… This is a very complicated space to navigate to ensure one finds the best configuration across all these parameters which are not homogeneous :sweat_smile: anyway, thanks for the transfer and let’s see if that works with one of these, I just hope we don’t have to pay everytime things crash


## Post 17 by @Yann.Sagon (2025-07-28T15:05:01.193Z)

Nicolas.Clairis1:
> I realize that Yggdrasil public-bigmem can not go higher than 750Go of RAM while Bamboo could go up to 1000 Go, but shared-bigmem has some nodes which can go up to 1500 Go of RAM… This is a very complicated space to navigate to ensure one finds the best configuration across all these parameters which are not homogeneous
Yes this is one of the issue we have with our clusters: they aren’t homogeneous and are upgraded during their whole life time.


## Post 18 by @Nicolas.Clairis1 (2025-08-08T15:15:08.096Z)

sure then I guess I just will patiently wait that the public-bigmem nodes are available in Bamboo as I really need as much RAM as possible (so the 750Go of Yggdrasil may not be enough) and as much time as possible (so shared-bigmem does not work)… I see that the user saturating shared-bigmem is around 1200/1500 jobs so let’s hope this liberates some nodes for public-bigmem at the end
Thanks for the help with Yggdrasil anyway and keep me informed if you change the policy or think of another potential solution please!


## Post 19 by @Yann.Sagon (2025-08-12T06:17:26.305Z)

Nicolas.Clairis1:
> I realize that Yggdrasil public-bigmem can not go higher than 750Go of RAM while Bamboo could go up to 1000 Go, but shared-bigmem has some nodes which can go up to 1500 Go of RAM…
You may be interested to know it is possible to list the partitions of the three clusters in one shot:
```
(baobab)-[root@login1 ~]$ sinfo -M all -p shared-bigmem,public-bigmem -N -o '%N %m %p %R'
CLUSTER: bamboo
NODELIST MEMORY PRIO_TIER PARTITION
cpu044 1024000 1 public-bigmem
cpu044 1024000 1 shared-bigmem
cpu045 1024000 1 public-bigmem
cpu045 1024000 1 shared-bigmem

CLUSTER: baobab
NODELIST MEMORY PRIO_TIER PARTITION
cpu203 512000 1 shared-bigmem
cpu218 512000 1 shared-bigmem
cpu219 512000 1 shared-bigmem
cpu245 256000 2 public-bigmem
cpu245 256000 1 shared-bigmem
cpu246 224000 2 public-bigmem
cpu246 224000 1 shared-bigmem
cpu312 1024000 1 shared-bigmem
cpu313 1024000 1 shared-bigmem

CLUSTER: yggdrasil
NODELIST MEMORY PRIO_TIER PARTITION
cpu113 770000 1 public-bigmem
cpu113 770000 1 shared-bigmem
cpu115 770000 1 public-bigmem
cpu115 770000 1 shared-bigmem
cpu120 1546348 1 shared-bigmem
cpu121 1546348 1 shared-bigmem
cpu122 1546348 1 shared-bigmem
```


## Post 20 by @nicolas.clairis (2025-08-18T12:47:24.213Z)

Hi thanks for the command.
I must say I’m a bit puzzled: All the batches on Bamboo public-bigmem are awaiting right now, and nobody is currently using shared-bigmem (see screenshot below).
image
image1338×257 8.92 KB
[image1338×257 8.92 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/cdc741af8779ef6f6d3e4bd6ab24e7cb8f5789de.png)
Could you tell us how we can check why that is the case and how I can track when they will be available? I thought that cpu044 and cpu045 were only used by public-bigmem and shared-bigmem but it seems they are also used by other partitions?
