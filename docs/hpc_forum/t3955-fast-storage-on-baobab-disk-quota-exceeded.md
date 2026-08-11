# Fast storage on Baobab - Disk quota exceeded

- Source: https://hpc-community.unige.ch/t/3955

- Created: 2025-05-20T10:45:23.319Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Malte.Algren (2025-05-20T10:45:23.362Z)

Hi,
I have moved some files out of the fast storage `/srv/fast/share/rodem` so I could move some other files in, but I’m getting hit by “Disk quota exceeded.”.
Is the fast disk lock cause I’m removing files to free up space?
Malte


## Post 2 by @Adrien.Albert (2025-05-21T10:05:24.787Z)

Hi @Malte.Algren[@Malte.Algren](https://hpc-community.unige.ch/u/malte.algren)
Indeed you have exceed your quota of 500 GB on `/srv/fast`.
```
(baobab)-[root@app8 ~]$ sudo xfs_quota -x -c 'quota -u algren -h' /srv/san
Disk quotas for User algren (404573)
Filesystem   Blocks  Quota  Limit Warn/Time    Mounted on
/dev/mapper/mpatha
             624.5G   500G   500G  00 [-none-] /srv/san
```
This means you can’t create any more files until your quota is restored.
Could you give the command you try to execute   ?
Where are you trying to move your data from and to?
