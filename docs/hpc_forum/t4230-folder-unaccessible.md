# Folder unaccessible

- Source: https://hpc-community.unige.ch/t/4230

- Created: 2026-02-23T13:53:43.277Z

- Tags: baobab

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Thibaut.Chataing (2026-02-23T13:53:43.330Z)

Hello,
I have an error trying to interact with some folders:
`baobab)-[chataint@login1 humanLISBET-paper]$ ls`
`ls: cannot access ‘dataset’: Communication error on send`
`dataset  quality_check.py  run_quality_check.sh`
`(baobab)-[chataint@login1 humanLISBET-paper]$ ll`
`ls: cannot access ‘dataset’: Communication error on send`
`total 47`
`d??? ? ?        ?             ?            ? dataset`
`-rw-r–r-- 1 chataint hpc_users 45095 Feb 23 14:43 quality_check.py`
`-rw-r–r-- 1 chataint hpc_users  2435 Feb 23 14:44 run_quality_check.sh`
`(baobab)-[chataint@login1 humanLISBET-paper]$ echo $PWD`
`/srv/beegfs/scratch/shares/schaerm/schaer2/video_sam2_pose/humanLISBET-paper`
What is happening ?
Best


## Post 2 by @Thibaut.Chataing (2026-02-23T13:54:24.633Z)

It’s back. I do not know what happened.


## Post 3 by @Adrien.Albert (2026-02-23T14:25:35.230Z)

Thibaut.Chataing:
> /srv/beegfs/scratch/shares/schaerm/schaer2/video_sam2_pose/humanLISBET-paper
Hello @Thibaut.Chataing[@Thibaut.Chataing](https://hpc-community.unige.ch/u/thibaut.chataing)
the issue has been fixed
