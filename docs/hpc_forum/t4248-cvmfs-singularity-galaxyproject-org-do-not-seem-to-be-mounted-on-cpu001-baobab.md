# Cvmfs singularity.galaxyproject.org do not seem to be mounted on cpu001 baobab

- Source: https://hpc-community.unige.ch/t/4248

- Created: 2026-03-12T10:44:49.663Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2026-03-12T10:44:49.716Z)

Hi,
From an interactive job:
```
(baobab)-[delislel@cpu001 A2B_all]$ ls /cvmfs/singularity.galaxyproject.org
ls: cannot access '/cvmfs/singularity.galaxyproject.org': No such file or directory
```
From the front end:
```
(baobab)-[delislel@login1 A2B_all]$ ls /cvmfs/singularity.galaxyproject.org
1  2  3  a  all  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
```
Thanks


## Post 2 by @Yann.Sagon (2026-03-13T15:56:12.135Z)

Hello,
I’ve tried to check what is going on without reboot (log, commands, clean cache etc) without luck. After a reboot it is working again. Still this isn’t very satisfying as workaround.
Best regards
