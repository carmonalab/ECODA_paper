# Mounting /cvmfs/atlas.cern.ch on baobab

- Source: https://hpc-community.unige.ch/t/3781

- Created: 2025-01-10T12:01:09.269Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Vilius.Cepaitis (2025-01-10T12:01:09.312Z)

Dear HPC admins,
Happy New Year!
I’m relying on `/cvmfs/` mount points for ATLAS-related work in the DPNC group. It seems that some of them are currently missing in baobab:
`source: no such file or directory: /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh`
The script above is the main entry point for setting up the analysis environment. This worked just before the break. Could you please take a look?
Thank you very much in advance.
Cheers,
Vilius


## Post 2 by @Yann.Sagon (2025-01-10T16:02:48.559Z)

Dear @Vilius.Cepaitis[@Vilius.Cepaitis](https://hpc-community.unige.ch/u/vilius.cepaitis)
Happy new year!
It was indeed crashed on login1 and only for atlas.cern.ch, weird. I’v restarted the service and it is working again.
Best
Yann
