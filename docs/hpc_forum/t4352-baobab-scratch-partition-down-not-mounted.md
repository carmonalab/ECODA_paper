# Baobab scratch partition down - not mounted

- Source: https://hpc-community.unige.ch/t/4352

- Created: 2026-07-19T05:34:48.453Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Steven.Schramm (2026-07-19T05:34:48.517Z)

## Primary informations
Username: schramm
Cluster: Baobab
## Description
The scratch partition appears to be down. Using df does not list the partition as it normally does, and any access (such as ls) to the normal scratch location of /srv/beegfs/scratch hangs indefinitely
## Steps to Reproduce
ls /srv/beegfs/scratch
## Expected Result
That the contents of the directory are listed
## Actual Result
The command hangs until the terminal is forcefully closed. Everything works fine when on the home directory or other mounted partitions are accessed
Thank you in advance for looking into this!


## Post 2 by @Yann.Sagon (2026-07-20T07:02:56.893Z)

Hi, this is solved, see [2026] Current issues on HPC Cluster - #25 by Yann.Sagon[[2026] Current issues on HPC Cluster - #25 by Yann.Sagon](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/25). Best


## Post 3 by @Steven.Schramm (2026-07-20T07:14:11.488Z)

Thank you very much!
