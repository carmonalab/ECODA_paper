# Mount issues on computing nodes of bamboo

- Source: https://hpc-community.unige.ch/t/3905

- Created: 2025-04-02T06:41:35.521Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lucille.Delisle1 (2025-04-02T06:41:35.565Z)

On the front-end of bamboo I can acess both CVMFS and acanas:
```
(bamboo)-[delislel@login1 ~]$ ls /acanas
biagetti  biosc  brunetti  celen  dpnc-test  ehret  hpc-exchange  m-BioinfoSupport  murrayry  sagon  seijo  test-astro  tonolla  wengery
(bamboo)-[delislel@login1 ~]$ ls /cvmfs/singularity.galaxyproject.org/
1  2  3  a  all  b  c  d  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
```
But in compute nodes I cannot:
```
(bamboo)-[delislel@login1 ~]$ salloc
salloc: Granted job allocation 294958
salloc: Nodes cpu001 are ready for job
(bamboo)-[delislel@cpu001 ~]$ ls /acanas
(bamboo)-[delislel@cpu001 ~]$ ls /acanas/m-BioinfoSupport
ls: cannot access '/acanas/m-BioinfoSupport': No such file or directory
(bamboo)-[delislel@cpu001 ~]$ ls /cvmfs/singularity.galaxyproject.org/
```
Thanks
PS: this is the cause of RStudio open on demand fail on bamboo (because it mounts /acanas/m-BioinfoSupport).


## Post 2 by @Yann.Sagon (2025-04-02T09:16:18.037Z)

Lucille.Delisle1:
> PS: this is the cause of RStudio open on demand fail on bamboo (because it mounts /acanas/m-BioinfoSupport).
Thanks for the debug!
We are awaiting for the change to be done by the storage team to solve the issue about NFS. For CVMFS, We’ll check.
Best


## Post 3 by @Yann.Sagon (2025-04-02T11:41:38.585Z)

Dear @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
Lucille.Delisle1:
> `ls /cvmfs/singularity.galaxyproject.org/`
This is solved!


## Post 4 by @Yann.Sagon (2025-04-02T13:58:35.528Z)

Dear @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
Lucille.Delisle1:
```
(bamboo)-[delislel@cpu001 ~]$ ls /acanas/m-BioinfoSupport
ls: cannot access '/acanas/m-BioinfoSupport': No such file or directory
```
This is solved as well:
```
(bamboo)-[root@cpu001 ~]$ ls -d /acanas/m-BioinfoSupport
/acanas/m-BioinfoSupport
```
