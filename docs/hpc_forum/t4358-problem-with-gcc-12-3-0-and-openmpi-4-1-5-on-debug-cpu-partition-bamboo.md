# Problem with GCC/12.3.0 and OpenMPI/4.1.5 on debug-cpu partition (bamboo)

- Source: https://hpc-community.unige.ch/t/4358

- Created: 2026-07-26T13:42:48.179Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Simon.Hug (2026-07-26T13:42:48.254Z)

## Primary informations
Username: hugs
Cluster: bamboo
## Description
modules GCC/12.3.0 and OpenMPI/4.1.5 fail to load on debug-cpu (which they did in the past)
## Steps to Reproduce
To test some code relying on the rjags odule on the debug-cpu  I encountered a problem when loading the GCC/12.3.0 and OpenMPI/4.1.5. More specifically
ml GCC/12.3.0  OpenMPI/4.1.5 (on debug-cpu partition)
fails
## Actual Result
ml GCC/12.3.0  OpenMPI/4.1.5
Module ERROR: Magic cookie ‘#%Module’ missing
In ‘/opt/ebmodules/all/Core/GCC/12.3.0.lua’
Please contact root@localhost[root@localhost](mailto:root@localhost)
ERROR: Unable to locate a modulefile for ‘OpenMPI/4.1.5’
Thanks for your help, simon
P.S. Actually there seems to be a larger problem, as also loading earlier versions fails:
module load      GCC/10.3.0  OpenMPI/4.1.1 R/4.1.0
Module ERROR: Magic cookie ‘#%Module’ missing
In ‘/opt/ebmodules/all/Core/GCC/10.3.0.lua’
Please contact root@localhost[root@localhost](mailto:root@localhost)
ERROR: Unable to locate a modulefile for ‘OpenMPI/4.1.1’
ERROR: Unable to locate a modulefile for ‘R/4.1.0’


## Post 2 by @Gael.Rossignol (2026-07-31T12:54:13.942Z)

Simon.Hug:
> ml GCC/12.3.0 OpenMPI/4.1.5
Dear simon,
Sorry for delay to answer.
I’m not sure whats is happening on this node so I drain this one. If you work on cpu002 all is working fine :
```
(base) (bamboo)-[rossigng@login1 ~]$ salloc
salloc: Granted job allocation 4174244
salloc: Nodes cpu002 are ready for job
(base) (bamboo)-[rossigng@cpu002 ~]$
(base) (bamboo)-[rossigng@cpu002 ~]$ ml GCC/12.3.0 OpenMPI/4.1.5
(base) (bamboo)-[rossigng@cpu002 ~]$ gcc --version
gcc (GCC) 12.3.0
Copyright (C) 2022 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

(base) (bamboo)-[rossigng@cpu002 ~]$
```
I will check to reinstall the node because I think something has been changed in the configuration.
Best regards


## Post 3 by @Simon.Hug (2026-07-31T13:11:27.258Z)

works like a charm, now. thanks


## Post 4 by @Gael.Rossignol (2026-07-31T13:34:14.909Z)

After node reinstallation all ios now working fine :
```
(base) (bamboo)-[rossigng@login1 ~]$ salloc
salloc: Granted job allocation 4174264
salloc: Nodes cpu001 are ready for job
(base) (bamboo)-[rossigng@cpu001 ~]$ module load GCC/10.3.0 OpenMPI/4.1.1 R/4.1.0
(base) (bamboo)-[rossigng@cpu001 ~]$ gcc --version
gcc (GCC) 10.3.0
Copyright (C) 2020 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```
Have a nice day.
