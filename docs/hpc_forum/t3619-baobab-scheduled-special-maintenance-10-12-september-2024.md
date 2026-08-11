# Baobab scheduled special maintenance: 10-12 September 2024

- Source: https://hpc-community.unige.ch/t/3619

- Created: 2024-08-30T11:47:25.536Z

- Tags: baobab

- Posts: 3

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Yann.Sagon (2024-08-30T11:47:25.619Z)

Dear users,
as just announced on the baobab-announce@ mailing list, we will be performing special hardware maintenance on the Baobab HPC cluster from 10th to 12th September 2024 (inclusive).
The maintenance will start at 08:00 +0100 and you will receive an email when the maintenance is finished.
The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files).
If you submit a job in the meantime, make sure that the expected wall time (duration) does not overlap with the start of the maintenance, or your job will be scheduled after the maintenance.
What will be done during this maintenance:
We’ll be upgrading our Infiband fast networking stack from 40Gb QDR to 100Gb EDR. We already replaced ~150 Infiband network cards in our compute nodes and servers a few weeks ago, and now we’ll replace all our Infiband switches.
This will allow us to keep the cluster up to date and improve performance.
This isn’t replacing our routine maintenance which is planed from 2th to 3th October 2024  (an email with details will be sent)
Thank you for your understanding.
Best regards,
The HPC Team


## Post 2 by @Adrien.Albert (2024-09-10T06:42:33.240Z)




## Post 3 by @Yann.Sagon (2024-09-13T07:52:21.671Z)

Dear users,
The maintenance is over. This was the longest maintenance we have ever done: three full days in the data centre with Adrien, Gaël and myself.
We replaced 13 Infiniband switches by faster ones. We also took the opportunity to reroute all the Ethernet and power cables to ensure better cooling of the compute nodes.
The benefit for you: better bandwidth between nodes for access to storage and computation.
The benefit for us: easier to manage compute nodes, fewer network issues
Thank you for your patience, best regards
Yann for the HPC team


## Post 4 by @Yann.Sagon (2024-09-13T08:10:17.972Z)

Some random pictures of Baobab and us working in the DC
20240911_114151
20240911_1141513024×4032 2.97 MB
[20240911_1141513024×4032 2.97 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/c3654c2b0613c5f61c1ad4f70cc58043f3a86776.jpeg)
20240910_092645
20240910_0926453648×2736 2.86 MB
[20240910_0926453648×2736 2.86 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/fbfa577d1788faaa40a0c36b70c69ab0a811c936.jpeg)
20240910_142118
20240910_1421184032×3024 2.3 MB
[20240910_1421184032×3024 2.3 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/ada32da07b3547cedfc63f92048983e1cd00e894.jpeg)
20240911_114151
20240911_1141513024×4032 2.97 MB
[20240911_1141513024×4032 2.97 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/c3654c2b0613c5f61c1ad4f70cc58043f3a86776.jpeg)
20240911_165425
20240911_1654254032×3024 3.1 MB
[20240911_1654254032×3024 3.1 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/af03449db856541b625f6fb3c0213aa87d6d81aa.jpeg)
20240912_163630
20240912_1636302016×1512 1.44 MB
[20240912_1636302016×1512 1.44 MB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d162639275d2e060ad868cc5c0ff485cfcc28f53.jpeg)


## Post 5 by @Adrien.Albert (2024-09-17T09:02:15.900Z)
