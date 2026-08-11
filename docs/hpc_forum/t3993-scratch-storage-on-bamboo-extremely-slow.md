# Scratch storage on Bamboo extremely slow

- Source: https://hpc-community.unige.ch/t/3993

- Created: 2025-06-26T13:26:03.884Z

- Tags: bamboo

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Berk.Gercek (2025-06-26T13:26:03.925Z)

I’m currently running an interactive job on Bamboo (cpu023, shared-cpu) and the scratch partition is slow enough as to be unusable. Operations like `ls` take ~1-2s to perform, and a file transfer to local /tmp/ is proceeding at about 4-5 MB/s maximum speed.
Is there an issue with BeeGFS for that file system?
Thanks!
Berk


## Post 2 by @Lucille.Delisle1 (2025-06-26T13:34:48.642Z)

I was about to report the same. I did `time ls` and
```
real    0m38.679s
user    0m0.003s
sys     0m0.016s
```


## Post 3 by @Yann.Sagon (2025-06-27T07:20:04.099Z)

Hi, @Berk.Gercek[@Berk.Gercek](https://hpc-community.unige.ch/u/berk.gercek) and @Lucille.Delisle1[@Lucille.Delisle1](https://hpc-community.unige.ch/u/lucille.delisle1)
I think the issue was due to high load on the storage as for now it seems ok to me:
I copied a 12GB file from scratch to `/tmp` in 23 sec which is not exceptional but not bad:)
```
(bamboo)-[sagon@cpu001 scratch]$ time cp mybigfile /tmp

real    0m23.624s
user    0m0.013s
sys     0m6.796s
(bamboo)-[sagon@cpu001 scratch]$ ls -lah mybigfile
-rw-r--r-- 1 sagon unige 12G Jun 27 09:14 mybigfile
```
Can you confirm it is working for you as expected? Unfortunately it is hard to diagnose such performance issue because it is highly dependent on the users workload.


## Post 4 by @Lucille.Delisle1 (2025-06-27T08:01:05.649Z)

Hi @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon) ,
I confirm that this morning everything is fine.
Thanks
