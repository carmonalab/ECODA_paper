# CPU removal on Baobab

- Source: https://hpc-community.unige.ch/t/3600

- Created: 2024-08-20T08:54:15.429Z

- Tags: baobab

- Posts: 1

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Gael.Rossignol (2024-08-20T08:54:15.526Z)

Dear Users,
I would like to inform you that we will be removing the old CPUs from Baobab in a few weeks.
Generation
Modele
Freq
CPUcount
Architecture
Nodeset
V3
E5-2660V0
2.20GHz
16 cores
“Sandy Bridge-EP” (32 nm)
cpu[001-005,007-008,045-056,058]
V3
E5-2670V0
2.60GHz
16 cores
“Sandy Bridge-EP” (32 nm)
cpu[059,061-062]
V3
E5-4640V0
2.40GHz
32 cores
“Sandy Bridge-EP” (32 nm)
cpu[186]
These CPUs are currently used in public partitions, so the related partitions will be unavailable on this cluster once the servers are removed.
As Bamboo is about to be deployed in production and is currently running on newly installed hardware, I suggest beginning the migration process.
Public partitions have been deployed with the same names.
The reason for this removal is to gain space in the datacenter and reorganize the placement to reduce issues.
This will also help us provide more powerful EasyBuild software compiled for recent hardware.
Best regards,
