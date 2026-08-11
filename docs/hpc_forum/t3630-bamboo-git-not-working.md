# [bamboo] git not working

- Source: https://hpc-community.unige.ch/t/3630

- Created: 2024-09-11T09:56:13.851Z

- Tags: bamboo

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-09-11T09:56:13.898Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support) ,
cluster: bamboo
user: falkiewi
what did you try:
```
cd scratch/
mkdir test
cd test
git init
```
what was the error message:
```
error: /srv/beegfs/scratch/users/f/falkiewi/test/.git/hooks/pre-merge-commit.sample: close error: Remote I/O error
fatal: cannot copy '/usr/share/git-core/templates/hooks/pre-merge-commit.sample' to '/srv/beegfs/scratch/users/f/falkiewi/test/.git/hooks/pre-merge-commit.sample': Remote I/O error
```
Thank you in advance,
Kind regards,
Maciej Falkiewicz


## Post 2 by @maciej.falkiewicz (2024-09-11T12:25:28.178Z)

@support[@support](https://hpc-community.unige.ch/groups/support) it is probably the whole filesystem that is broken.


## Post 3 by @Adrien.Albert (2024-09-12T09:30:35.269Z)

Hi @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/22) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Bamboo Scratch Storage Unavailable
Dear HPC Users, 
The scratch storage on Bamboo is currently unavailable due to an ongoing issue. Our team has already contacted the provider and we are actively working with them to resolve the situation as quickly as possible. 
Please note that the scratch storage have been unmounted on compute and login nodes and will remain unavailable until further notice. We will keep you updated as soon as we have more information on the situation. 
Thank you for your und…
