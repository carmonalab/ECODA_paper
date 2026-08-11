# Fast storage on Baobab seems to be empty

- Source: https://hpc-community.unige.ch/t/3808

- Created: 2025-02-03T13:50:44.210Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Matthew.Leigh (2025-02-03T13:50:44.257Z)

Hi HPC Team,
It seems that the /srv/fast storage is empty again on babab. I am not able to see or create and directories if I go onto:
```
(baobab)-[leighm@cpu331 fast]$ pwd
/srv/fast
(baobab)-[leighm@cpu331 fast]$ ls
(baobab)-[leighm@cpu331 fast]$ cd /srv/fast/share/rodem
bash: cd: /srv/fast/share/rodem: No such file or directory
```
Is this related to a previous issue around this time last year:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/2) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> BAOBAB: fast storage down
Dear users, 
the fast storage /srv/fast is currently down since the latest maintenance. 
We need to go on site to put it back on production, this should be done today. 
We apologize for the inconvenience. 
Thank you for your understanding. 
Best regards, 
Status : Solved green_circle :
start: 2024-01-10T02:26:00Z (UTC) 
end:2024-01-15T15:00:00Z (UTC)
All the best,
Matthew Leigh


## Post 2 by @Yann.Sagon (2025-02-03T16:04:00.987Z)

Dear @Matthew.Leigh[@Matthew.Leigh](https://hpc-community.unige.ch/u/matthew.leigh) it is solved as far as I can see. We wanted to convert the share to NFSv4 but we’ll do that during Baobab next maintenance Baobab scheduled maintenance: 18-21 February 2025[Baobab scheduled maintenance: 18-21 February 2025](https://hpc-community.unige.ch/t/baobab-scheduled-maintenance-18-21-february-2025/3809)
Best


## Post 3 by @Matthew.Leigh (2025-02-03T16:14:03.563Z)

Thank you. It is visible again
