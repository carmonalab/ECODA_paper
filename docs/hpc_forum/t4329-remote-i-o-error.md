# Remote I/O error

- Source: https://hpc-community.unige.ch/t/4329

- Created: 2026-06-29T07:52:47.775Z

- Tags: baobab

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibaut.Chataing (2026-06-29T07:52:47.775Z)

Hi, It seems there is a new “Remote I/O error” on baobab :
> (baobab)-[chataint@login1 run_20260626_092030]$ cat run.log
> cat: run.log: Remote I/O error
Is there a problem with storage ?
Best,
Thibaut


## Post 2 by @Gael.Rossignol (2026-06-29T08:21:27.464Z)

Dear Thibaut,
No error was raised on the Baobab storage, could you please give me the full path to investigate?
Best regards


## Post 3 by @Thibaut.Chataing (2026-06-29T08:25:17.838Z)

Hi,
Thanks for helping.
> (baobab)-[chataint@login1 run_20260626_092030]$ echo $PWD
> /srv/beegfs/scratch/shares/schaerm/schaer2/video_sam2_pose/humanLISBET-paper/clusterAnalysis/results/run_20260626_092030


## Post 4 by @Gael.Rossignol (2026-06-29T09:29:33.176Z)

Dear Thibaut,
One disk bay was stuck and I unplug and replug the bay and all is now working fine.
Sorry for inconvenience,
