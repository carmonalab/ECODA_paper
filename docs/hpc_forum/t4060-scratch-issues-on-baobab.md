# Scratch issues on baobab

- Source: https://hpc-community.unige.ch/t/4060

- Created: 2025-08-22T13:49:47.197Z

- Tags: baobab

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Christian.Scheulen (2025-08-22T13:49:47.256Z)

## Primary informations
Username: scheulen
Cluster: Baobab
## Description
Cannot access files on scratch space. The metadata for all files is still available, and folders are browsable, but trying to access any files (tested with text files and binary data) yields a `READ ERROR`. Applies to both personal and group scratch spaces.
Additionally, the problem seems to be confined to baobab. On bamboo, files on the scratch space are accessible without issues.
## Steps to Reproduce
Log in to baobab, switch to `/srv/beegfs/scratch/groups`, try to access some file.
## Expected Result
File opens and shows data
## Actual Result
`READ ERROR`


## Post 2 by @Adrien.Albert (2025-08-22T16:44:43.693Z)

Hello,
The issue is related to [2025] Current issues on HPC Cluster - #20 by Yann.Sagon[[2025] Current issues on HPC Cluster - #20 by Yann.Sagon](https://hpc-community.unige.ch/t/2025-current-issues-on-hpc-cluster/3788/20) and has been resolved.
