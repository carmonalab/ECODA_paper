# Job crashes w/o log

- Source: https://hpc-community.unige.ch/t/4008

- Created: 2025-07-17T13:55:28.984Z

- Tags: baobab, slurm

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Anton.Hanke (2025-07-17T13:55:29.031Z)

Dear HPC team,
Recently I have had problems with jobs on baobab running from the home storage.
My jobs run for some time without issues, but then randomly seem to stop, so apbruptly that both the log files of the job itself and the slurm out file end without any information as to why the job stopped.
This results in the job not writing required files, needed to continue the job.
Further dependency jobs within the same directory then fail within one second without generating a slurm out file at all.
Could you have some info/help with this please?
This has me stumped, as I can not even debug whether this is due to the content/commands within the job or a cluster side issue.
Example job would be: `slurm-1098047.out`
At: `/home/users/h/hankea/folding_bh3/oneopes`
Thank you for your help!
Best wishes,
Anton Hanke


## Post 2 by @Yann.Sagon (2025-07-17T14:41:05.627Z)

Dear @Anton.Hanke[@Anton.Hanke](https://hpc-community.unige.ch/u/anton.hanke)
The reason is that your home quote was exceeded:
```
/var/log/beegfs-meta-home_meta01.log:(3) Jul10 11:24:44 Worker4 [Quota Enforcement for create] >> User size quota exceeded. UID: 416089; GID: 5000
[...]
/var/log/beegfs-storage-home.log:(3) Jul17 10:50:44 Worker5-1 [WriteChunkFileMsg incoming] >> User size quota exceeded. UID: 416089; GID: 5000
```
You need to ensure that you have enough space in your home directory or better use the scratch space.
Best regards
Yann


## Post 3 by @Anton.Hanke (2025-07-17T14:56:20.322Z)

Ah thanks!
I was not aware.


## Post 4 by @Yann.Sagon (2025-07-17T15:06:11.142Z)

You can check here for more information and how to check your current quota/usage: hpc:storage_on_hpc [eResearch Doc][hpc:storage_on_hpc [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/storage_on_hpc#check_disk_usage_on_the_clusters)
