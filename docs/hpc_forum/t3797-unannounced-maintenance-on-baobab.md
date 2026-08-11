# Unannounced maintenance on baobab?

- Source: https://hpc-community.unige.ch/t/3797

- Created: 2025-01-27T08:03:09.725Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2025-01-27T08:03:09.770Z)

Hi,
When I request GPU resources on baobab i get the following error: ReqNodeNotAvail, Reserved for maintenance
Im using the following command:
`salloc -c4 --partition=private-dpnc-gpu, --time=00-12:00:00 --mem=32GB --gres=gpu:1,VramPerGpu:2G --exclude=`
Is there soonish maintenance on baobab?
Update: I think it is only on the private-dpnc-gpu
Cheers,
Malte
image
image1910×155 7.96 KB
[image1910×155 7.96 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/29dfec995ba8cc9f51dd23bacff0d8246b183dba.png)


## Post 2 by @Gael.Rossignol (2025-01-27T10:49:40.534Z)

Malte.Algren:
> private-dpnc-gpu
Dear Malte,
No maintenance are done on Baobab. All DPNC nodes have been reserved by Paul, as you are member of the dpnc group you can use the reservation.
Best regards,


## Post 3 by @Malte.Algren (2025-01-27T20:23:37.955Z)

I think there is an incorrect setting in your reservation implementation. If users want to
I know some of the DPNC GPUs are temporarily reserved for some group members, which I’m guessing is why when I use
`salloc -c4 --partition=private-dpnc-gpu, --time=00-12:00:00 --mem=32GB --gres=gpu:1,VramPerGpu:2G --exclude=`
I get the following error: ReqNodeNotAvail, Reserved for maintenance
However, i can allocate the resources reserved for them by doing:
`salloc -c4 --partition=private-dpnc-gpu, --time=00-12:00:00 --mem=32GB --gres=gpu:1,VramPerGpu:2G --exclude= --reservation=private_dpnc_gpu`
and I allocate one of the dpnc gpus.
Is this how the reservation of GPUs is supposed to work? I would imagine that the users who have reserved the GPUs are the only ones who can allocate them (based on username)?
Cheers,
Malte


## Post 4 by @Yann.Sagon (2025-01-28T08:52:03.051Z)

Dear @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren)
You can check the current reservations on our clusters like that:
```
(baobab)-[root@login1 ~]$ scontrol show reser
ReservationName=private_dpnc_cpu StartTime=2025-01-24T12:23:13 EndTime=2025-01-31T23:59:59 Duration=7-11:36:46
   Nodes=cpu[084-090,209-213,226-229,277] NodeCnt=17 CoreCnt=560 Features=(null) PartitionName=private-dpnc-cpu Flags=OVERLAP,IGNORE_JOBS,SPEC_NODES,PART_NODES
   TRES=cpu=560
   Users=(null) Groups=private_dpnc Accounts=(null) Licenses=(null) State=ACTIVE BurstBuffer=(null)
   MaxStartDelay=(null)

ReservationName=private_dpnc_gpu StartTime=2025-01-24T12:23:25 EndTime=2025-01-31T23:59:59 Duration=7-11:36:34
   Nodes=gpu[002,012,017,023-024,044,047,049] NodeCnt=8 CoreCnt=804 Features=(null) PartitionName=private-dpnc-gpu Flags=OVERLAP,IGNORE_JOBS,SPEC_NODES,PART_NODES
   TRES=cpu=804
   Users=(null) Groups=private_dpnc Accounts=(null) Licenses=(null) State=ACTIVE BurstBuffer=(null)
   MaxStartDelay=(null)
```
As you can see, the reserved nodes  are available to all members of the `private_dpnc` group. And you are a member of this group.
