# File transfer procedure for an external user

- Source: https://hpc-community.unige.ch/t/3666

- Created: 2024-10-02T09:55:14.326Z

- Tags: all

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Maria.Suveges (2024-10-02T09:55:14.362Z)

Hi all,
my collaborator has an external account on baobab. He attempted to transfer files to and from the cluster. As far as I know, there are two methods: FileZilla and ssh. When I tried to download and install FileZilla, my laptop warned me that it has a virus and did not install it. So I did not propose FileZilla to my collaborator (his laptop may not be so well protected as mine, I do not want to expose him to this danger). Unfortunately, ssh does not work straightforwardly for him. How can he go ahead to transfer files without FileZilla?
Thank you!
Maria Süveges


## Post 2 by @Yann.Sagon (2024-10-04T13:25:28.829Z)

Dear @Maria.Suveges[@Maria.Suveges](https://hpc-community.unige.ch/u/maria.suveges)
I have updated the doc[doc](https://doc.eresearch.unige.ch/hpc/access_the_hpc_clusters#from_windows1) with a warning and a new link.
Anyway, you are free to use another software than filezilla if willing to, as long as it supports sftp.
Best
Yann


## Post 3 by @Maria.Suveges (2024-10-07T08:20:29.029Z)

Dear Yann,
Thanks for your reply. Does your warning mean that non-sponsored versions exist and those might be safe?
Do you have any proposition what applications other than FileZilla can (or should) we use?
Thank you! Best,
Maria


## Post 4 by @Yann.Sagon (2024-10-07T08:34:25.779Z)

Dear @Maria.Suveges[@Maria.Suveges](https://hpc-community.unige.ch/u/maria.suveges) as stated in my previous answers, I’ve updated the doc with a new link: the version available from this new link is the non sponsored version, please give a try.
I don’t have any proposition for an an alternative solution as this one is working fine.
Best


## Post 5 by @Maria.Suveges (2024-10-07T09:58:35.634Z)

Thanks and apologies, I misread your reply, it indeed works. Best, Maria
