# Transferring large amount of files

- Source: https://hpc-community.unige.ch/t/3553

- Created: 2024-07-21T22:17:23.394Z

- Tags: baobab

- Posts: 3

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Omar.Darwish (2024-07-21T22:17:23.448Z)

Hi all,
I need to transfer from Baobab scratch folder 400 folders of around 5 GB each to an external cluster. I would use Globus as for other HPC clusters, but unfortunately I can not in this case.
I am using `rsync` to achieve my goal. I would normally open my terminal through ssh, and run it over the night. But for some reason, since yesterday, I can not login to Baobab (while I tried to Bamboo and I can safely login, with my ssh key for Unige).
Are there solutions for large scale files transfer? I need to share my files with some colleagues and store them in this other cluster, so that I can delete them on the Baobab scratch.
Thank you in advance.


## Post 2 by @Raphael.Rubino (2024-07-22T09:45:15.725Z)

Hello,
Personally, I would `tar czf` the 400 folders into one archive and `scp` the tarball to wherever I need it to be.
Best regards


## Post 3 by @Adrien.Albert (2024-07-22T11:10:44.373Z)

hi @Omar.Darwish[@Omar.Darwish](https://hpc-community.unige.ch/u/omar.darwish)
It should work again,
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/21) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Baobab: Login node down
Dear users, 
The login node on Baobab have crashed. The server  have been rebooted and is available again. 
We apologize for any inconveniance caused 
Thank you for your understanding. 
Status : Solved green_circle
start: 2024-07-21T18:42:00Z (UTC) 
end:Invalid date (UTC)
As I’m off today, I won’t investigate any further.
The HPC team have already talked about globus. We’ll be exploring the subject and analyzing the feasibility of this project.
For now, you can use rsync:
doc.eresearch.unige.ch[doc.eresearch.unige.ch](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=rsync#rsync)
### Rsync - hpc:best_practices [eResearch Doc][Rsync - hpc:best_practices [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=rsync#rsync)
This page gives best practices and tips on how to use the clusters Baobab and Yggdrasil.
Or the @Raphael.Rubino[@Raphael.Rubino](https://hpc-community.unige.ch/u/raphael.rubino) 's idea.
Best Regards
