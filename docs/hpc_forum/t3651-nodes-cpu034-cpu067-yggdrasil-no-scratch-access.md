# Nodes cpu034 & cpu067 yggdrasil no scratch access?

- Source: https://hpc-community.unige.ch/t/3651

- Created: 2024-09-27T06:46:40.700Z

- Tags: yggdrasil

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Max.Briel (2024-09-27T06:46:40.738Z)

I’m running job arrays on yggdrasil and jobs assigned to nodes cpu34 & cpu067 are failing with exit code 0:53.
I’m writing to `/srv/beegfs/scratch/`
Jobs assigned to these nodes crash immediately and no log files are created on scratch, indicating issues with reading/writing to scratch from these nodes.


## Post 2 by @Adrien.Albert (2024-09-27T07:32:43.818Z)

Hi @Max.Briel[@Max.Briel](https://hpc-community.unige.ch/u/max.briel)
Thank you for reporting this problem. The nodes have been drained from production for analysis.
