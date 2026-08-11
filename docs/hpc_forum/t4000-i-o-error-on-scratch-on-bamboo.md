# I/O error on scratch on bamboo

- Source: https://hpc-community.unige.ch/t/4000

- Created: 2025-07-07T14:33:50.058Z

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-07-07T14:33:50.099Z)

Dear HPC,
I get a lot of `error reading 'whatever': Remote I/O error` on the scratch of bamboo.
Thanks,
Lucille


## Post 2 by @Adrien.Albert (2025-07-08T10:09:09.867Z)

Thank you @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)  for relaying. For more information about issue:
[2025] Current issues on HPC Cluster[[2025] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/16) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Description
Cluster: Bamboo 
We are currently experiencing an issue with the Bamboo scratch and Home storage. This share is temporarily unavailable. 
Our team is actively investigating the problem, and we will provide an update as soon as more information becomes available. 
We apologize for the inconvenience and appreciate your patience. 
Status : alert green_circle
start: 2025-07-07T14:50:00Z (UTC) 
end : 2025-07-07T15:10:00Z (UTC) 
Update
Filesystem impacted scratch: 
start: 2025-07-07T20:0…


## Post 3 by @nicolas.clairis (2025-07-10T09:40:03.700Z)

Hi
I thought that Bamboo was working again, but when I try to launch my Bamboo batches on public-cpu, I see that they appear as “(ReqNodeNotAvail, Reserved for maintenance)”. Is this related to the Bamboo issue or is it another issue with the public-cpu? I need to obtain some data soon, and these recurrent bugs of the cluster are worrying to be able to provide on time so I hope this gets fixed soon :sweat_smile:


## Post 4 by @Adrien.Albert (2025-07-10T10:37:58.353Z)

Hi Nicolas,
You’re really quick!
At the moment, there is no issue with the scratch filesystem. The message `(ReqNodeNotAvail, Reserved for maintenance)` means that the requested node is not available because it is reserved for maintenance.
In fact, we have created a reservation for a possible short maintenance tomorrow. If the maintenance is confirmed, an official communication will be sent out.
```
(bamboo)-[root@admin1 ~]$ scontrol show res
ReservationName=test_scratch StartTime=2025-07-11T09:00:00 EndTime=2025-07-11T13:00:00 Duration=04:00:00
   Nodes=cpu[001-052],gpu[001-003] NodeCnt=55 CoreCnt=6704 Features=(null) PartitionName=(null) Flags=IGNORE_JOBS,SPEC_NODES,ALL_NODES
   TRES=cpu=6704
   Users=alberta,root,sagon Groups=(null) Accounts=(null) Licenses=(null) State=INACTIVE BurstBuffer=(null)
   MaxStartDelay=(null)
```
I suspect you’re trying to use `srun` or `salloc` with a duration that exceeds the available time before the reservation starts.
For now, you can submit your job using `sbatch`.
Best regards,


## Post 5 by @nicolas.clairis (2025-07-10T12:06:59.747Z)

Hi
Yes it’s a long script (which is why I go via public-bigmem and not the other partitions, sorry for the typo and saying public-cpu in my previous post) so I request for 2 days but I guess this is why it doesn’t work if you guys have blocked it for tomorrow…
I am actually using sbatch so I don’t think that is the cause of the issue.
I guess I just need to wait for your feedback for the maintenance then.
Best,
Nicolas


## Post 6 by @nicolas.clairis (2025-07-11T07:28:52.820Z)

Any news as to when the maintenance will happen?


## Post 7 by @Adrien.Albert (2025-07-16T13:13:59.540Z)

Hello @nicolas.clairis[@nicolas.clairis](https://hpc-community.unige.ch/u/nicolas.clairis)
Sorry for the delay, we had the confirmation at the last minute:
[Bamboo] Short maintenance test on Bamboo – July 11, 2025 (09:00–13:00)[[Bamboo] Short maintenance test on Bamboo – July 11, 2025 (09:00–13:00)](https://hpc-community.unige.ch/t/bamboo-short-maintenance-test-on-bamboo-july-11-2025-09-00-13-00/4002) HPC Announce[HPC Announce](https://hpc-community.unige.ch/c/hpc-announce/6)
> Dear Bamboo users, 
Please be aware about a short maintenance on Bamboo for: 
three_o_clock Date:2025-07-11T07:00:00Z (UTC) 
nine_o_clock Duration: 09:00 to 13:00 (4 hours) 
wrench Purpose: 
We will be conducting scratch filesystem testing during this time to verify stability and performance after the recent issue. 
Bamboo login remains available but scratch will be umounted. All computes have been reserved during the mini maintenance. YOu can keep submitting jobs, once the reservation com…
