# Bamboo seems to have crashed

- Source: https://hpc-community.unige.ch/t/4174

- Created: 2025-12-17T14:12:49.734Z

- Tags: bamboo

- Posts: 14

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @nicolas.clairis (2025-12-17T14:12:49.780Z)

All my ongoing jobs suddenly crashed and I cannot access Bamboo via Filezilla anymore. Instead, I’m sent to a different folder which contains the following subfolders
image
image944×538 35.4 KB
[image944×538 35.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/7efc936cd2400ab29c48fec49145049c9a390d60.png)


## Post 2 by @Marco.Froelich (2025-12-17T14:15:44.921Z)

Second this: VSCode Remote - SSH login reports the following error:
Could not chdir to home directory <$HOME> Remote I/O error
-bash <$HOME>/.bash_profile: Remote I/O error
EDIT: Login possible again !


## Post 3 by @nicolas.clairis (2025-12-17T14:18:39.511Z)

Weirdly, I seem to have recovered my access both via Putty and FileZilla + possibility to launch jobs again and no way to know what happened nor why my scripts were abruptly stopped….
Side note: it’s not very fair to bill for the hours when we run the cluster when our jobs are stopped due to cluster issues. This has happened to me more than once and it would be nice to have a way to cancel these hours when it happens in the future…


## Post 4 by @Gael.Rossignol (2025-12-18T11:11:09.441Z)

Hello,
We didn’t get any faillure on Bamboo, I’m sorry to see your message. I check log and I wee 10 minutes downtime on storage servers.
image
image944×531 115 KB
[image944×531 115 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/eddee4b95ddcd42465c123137cd1e157c936703c.png)
We will investigate regarding this problem, sorry for inconvenience.
Best regards,


## Post 5 by @nicolas.clairis (2025-12-18T19:24:55.202Z)

thanks yes indeed it was super short but all jobs crashed abruptly and then it came back super quickly but in the meantime I still had to restart all the jobs that had crashed


## Post 6 by @Nicolas.Clairis1 (2025-12-31T09:37:55.798Z)

I think that it just happened again… This time I didn’t had any batch ongoing, but it’s totally impossible to connect to Bamboo since this morning


## Post 7 by @Hadrien.Cusin (2025-12-31T10:30:00.904Z)

Yes i have the same problem. It starts yesterday night, all of my job stop suddenly, and when i try some commands i get this answer “Remote I/O error” . Now i can’t even connect to my session, it says “ ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Connection timed out “


## Post 8 by @Cody.Cardenas (2025-12-31T10:58:48.357Z)

Just echoing the problem. I am also unable to ssh into bamboo; I get a port 22 time out error.


## Post 9 by @pablo.strasser1 (2025-12-31T14:15:56.226Z)

Got the same problem now.


## Post 10 by @Aurelie.Kamoun (2026-01-02T10:18:12.218Z)

Same here, i cannot connect to bamboo anymore (although i can still access other clusters): ssh: connect to host login1.bamboo.hpc.unige.ch port 22: Operation timed out


## Post 11 by @Nicolas.Clairis1 (2026-01-02T17:19:47.514Z)

the issue still seems ongoing so I guess more news next Monday when people are back from vacations


## Post 12 by @nicolas.clairis (2026-01-05T08:31:18.923Z)

it’s working again on my side!


## Post 13 by @nicolas.clairis (2026-01-05T08:40:55.215Z)

> Actually no, only the access is back to work (and the capacity to launch jobs) but then they crash immediately with the following error message:
> HDF5-DIAG: Error detected in HDF5 (1.8.12) thread 0:
> #000: H5Dio.c line 179 in H5Dread(): can’t read data
> major: Dataset
> minor: Read failed
> #001: H5Dio.c line 547 in H5D__read(): can’t read data
> major: Dataset
> minor: Read failed
> #002: H5Dchunk.c line 1836 in H5D__chunk_read(): unable to read raw data chunk
> major: Low-level I/O
> minor: Read failed
> #003: H5Dchunk.c line 2862 in H5D__chunk_lock(): unable to read raw data chunk
> major: Low-level I/O
> minor: Read failed
> #004: H5Fio.c line 113 in H5F_block_read(): read through metadata accumulator failed
> major: Low-level I/O
> minor: Read failed
> #005: H5Faccum.c line 258 in H5F_accum_read(): driver read request failed
> major: Low-level I/O
> minor: Read failed
> #006: H5FDint.c line 142 in H5FD_read(): driver read request failed
> major: Virtual File Layer
> minor: Read failed
> #007: H5FDsec2.c line 725 in H5FD_sec2_read(): file read failed: time = Mon Jan  5 09:30:30 2026
> , filename = ‘/srv/beegfs/scratch/users/c/clairis/fMRI_analysis/results/CAPS/CAPS/CAPS_SeedFree1c_K2__GM50.mat’, file descriptor = 787, errno = 121, error message = ‘Remote I/O error’, buf = 0x14d4e6778fa0, total read size = 500, bytes this sub-read = 500, bytes actually read = 18446744073709551615, offset = 523929
> major: Low-level I/O
> minor: Read failed
so I’m guessing it’s not fully fixed yet. The issue seems to be here [2026] Current issues on HPC Cluster - #3 by Yann.Sagon[[2026] Current issues on HPC Cluster - #3 by Yann.Sagon](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185/3) and is written as still ongoing so gonna wait for the announcement of the fix before trying again


## Post 14 by @Yann.Sagon (2026-01-05T09:21:18.691Z)

Dear all, we returned from vacation today, this is our first issue of the year :fireworks: !
It is now solved. [2026] Current issues on HPC Cluster[[2026] Current issues on HPC Cluster](https://hpc-community.unige.ch/t/2026-current-issues-on-hpc-cluster/4185)
