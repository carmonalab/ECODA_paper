# /cvmfs/dampe.cern.ch missing on Baobab

- Source: https://hpc-community.unige.ch/t/4095

- Created: 2025-09-21T12:55:26.888Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Paul.Coppin (2025-09-21T12:55:26.944Z)

Dear HPC team,
It appears some cvmfs directories are not loaded on Baobab. Particularly:
/cvmfs/dampe.cern.ch
/cvmfs/geant4.cern.ch
/cvmfs/sft.cern.ch
```
(baobab)-[coppinp@login1 ~]$ cd /cvmfs
(baobab)-[coppinp@login1 cvmfs]$ cd dampe.cern.ch
-bash: cd: dampe.cern.ch: No such file or directory
(baobab)-[coppinp@login1 cvmfs]$ cd sft.cern.ch
-bash: cd: sft.cern.ch: No such file or directory
(baobab)-[coppinp@login1 cvmfs]$ cd geant4.cern.ch
-bash: cd: geant4.cern.ch: No such file or directory
```
Can you have a look?
(On Yggdrasil and Bamboo, these are accessible, and all is well)
All the best,
Paul


## Post 2 by @Alexander.Froch (2025-09-22T11:15:25.293Z)

Hi HPC team,
I have the same issue with `/cvmfs/atlas.cern.ch`
Cheers,
Alex


## Post 3 by @Gael.Rossignol (2025-09-22T13:18:24.618Z)

Hello,
We had an issue on the /tmp usage this morning but all seems to be ok now. Do you always have issue?
Best regards,


## Post 4 by @Paul.Coppin (2025-09-23T09:25:51.573Z)

Now they are all there. Thank you for resolving the problem!
