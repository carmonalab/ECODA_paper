# Changes to `apptainer` containers on HPC?

- Source: https://hpc-community.unige.ch/t/3457

- Created: 2024-05-22T15:25:38.418Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Imahn.Shekhzadeh1 (2024-05-22T15:25:38.461Z)

Hi,
Were there any recent changes to apptainer/singularity on the HPC in the past few days? I have a program installed in an apptainer/singularity container that ran a few days ago, but it is not running anymore even though I did not perform any changes on the container itself.
I was thinking of rebuilding it on the HPC, but before doing that, wanted to make there isn’t anything I need to adjust.
Best,
Imahn


## Post 2 by @Adrien.Albert (2024-05-24T09:49:00.322Z)

Hi @Imahn.Shekhzadeh1[@Imahn.Shekhzadeh1](https://hpc-community.unige.ch/u/imahn.shekhzadeh1)
Apptainer is not a module, it’s a system package. We haven’t updated it recently.
Are you sure you haven’t modified anything, a file in a bound directory, for example?
