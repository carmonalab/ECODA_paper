# Rsync: [sender] read errors mapping "/srv/beegfs/scratch/users/s/strassep/dataset-offline-np/Test/0.npy": Remote I/O error (121) On Bamboo

- Source: https://hpc-community.unige.ch/t/4178

- Created: 2025-12-30T21:53:20.354Z

- Tags: bamboo

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @pablo.strasser1 (2025-12-30T21:53:20.396Z)

There is error on Bamboo when accessing files.


## Post 2 by @pablo.strasser1 (2025-12-30T22:14:24.424Z)

Looking more in detail at the problem, it only seem to affect some file. Doing a sha1sum on a directory I see some file can be read some other not. The failing are always the same. I believe some of the storage node are down or something like that. It also only affect scratch, not the home directory.


## Post 3 by @Nicolas.Clairis1 (2025-12-31T09:39:40.430Z)

Not sure if related, but since this morning, I just cannot access Bamboo in any way (tried both via the terminal and FileZilla and via different wifi but both give an error message suggesting Bamboo is completely off…
> “ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Operation timed out” this is the error message on the terminal


## Post 4 by @Nicolas.Clairis1 (2025-12-31T17:12:01.344Z)

the more general issue is there in case: Bamboo seems to have crashed - #6 by Nicolas.Clairis1[Bamboo seems to have crashed - #6 by Nicolas.Clairis1](https://hpc-community.unige.ch/t/bamboo-seems-to-have-crashed/4174/6)


## Post 5 by @pablo.strasser1 (2025-12-31T17:42:26.638Z)

Thanks I saw it.
Have a happy new year.


## Post 6 by @pablo.strasser1 (2026-01-05T08:16:55.422Z)

Baobab and Bamboo are reachable again, however, the storage problem for scratch is not yet solved. A sha1sum * in a directory with multiple file will show some file that are no readable.


## Post 7 by @Yann.Sagon (2026-01-05T09:23:25.262Z)

This is now solved. [2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185)
