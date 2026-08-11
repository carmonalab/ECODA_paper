# [yggdrasil] Issues while downloading data using OOD interface

- Source: https://hpc-community.unige.ch/t/4359

- Created: 2026-07-27T07:27:33.010Z

- Tags: yggdrasil, openondemand

- Posts: 10

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Sebastien.Miche (2026-07-27T07:27:33.081Z)

Good morning everyone,
I’ve encountered issues this morning while trying to download data using the OOD interface (files tab).
I’ve tried downloading the folder containing my data, but most of the time, I get a message indicating a “Timeout while trying to determine directory size.” (folder size is 1.5 Gb).
I’ve tried downloading subfolders/files one by one, but again without success. I sometimes get a “Proxy error” or nothing is displayed at all and the download never starts.
Am I the only person encountering such issues or is it a well known problem? Please note that I’m working on an external network (Muséum) and I don’t have any other way (afaik) to download these files.
Thanks in advance for your answers,
Sébastien


## Post 2 by @Gael.Rossignol (2026-07-31T13:35:31.865Z)

Sebastien.Miche:
> I’ve tried downloading the folder containing my data, but most of the time, I get a message indicating a “Timeout while trying to determine directory size.” (folder size is 1.5 Gb).
Dear Sebastien,
Yggdrasil has been reinstalled this week, I just would like you test again to check if the problem persists?
Thanks for help,


## Post 3 by @Sebastien.Miche (2026-08-06T15:26:33.623Z)

Dear Gael,
Sorry for replying so late, I was caught up in other stuff.
Yes, the problem persists, I still have the same error message while trying to download files from yggdrasil using the OOD interface.
image
image1652×130 7.58 KB
[image1652×130 7.58 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/1a8604c1834da959d52947473946f93d53c5ca8f.png)


## Post 4 by @Adrien.Albert (2026-08-07T12:47:34.286Z)

Hello @Sebastien.Miche[@Sebastien.Miche](https://hpc-community.unige.ch/u/sebastien.miche)
I suspect that you are trying to download files from a directory containing a very large number of files, which may be causing the timeout.
In that case, you would likely experience the same issue from a standard terminal session.
Could you please provide the path to the directory you are trying to access so that we can investigate further?
Best,


## Post 5 by @Sebastien.Miche (2026-08-10T06:54:02.942Z)

Hi Adrien! Yes this is correct, I am trying to download an output folder containing output results from Orthofinder (which can generate large amount of subfolders and files especially when used on complete genomes).
Here is the link to the directory:
/srv/beegfs/scratch/users/m/miches/genomes_for_orthofinder/proteomes_orthofinder/OrthoFinder
Best


## Post 6 by @Adrien.Albert (2026-08-10T08:14:58.739Z)

Hello @Sebastien.Miche[@Sebastien.Miche](https://hpc-community.unige.ch/u/sebastien.miche),
The following issue has been reported here:
Open OnDemand – 10 May 24[Open OnDemand – 10 May 24](https://discourse.openondemand.org/t/proxy-error-when-downloading-folders-and-files-under-home-directory-of-file-menu-in-open-ondemand/3480)
### Proxy error when downloading folders and files under "Home Directory" of...[Proxy error when downloading folders and files under "Home Directory" of...](https://discourse.openondemand.org/t/proxy-error-when-downloading-folders-and-files-under-home-directory-of-file-menu-in-open-ondemand/3480)
Get Help
Hi Support,  When I download folders containing 5417 files and 70 sub-folders from “Home Directory” in “Files” menu, with a total size of 1.08GB, I occasionally encounter the following proxy timeout error. For example, when I attempt to download...
Indeed, the Open OnDemand timeout is set to 5 seconds by default, while the `du` command takes nearly 7 seconds to complete in your case:
```
(yggdrasil)-[root@app6 dashboard]$ time du -sh /srv/beegfs/scratch/users/m/miches/genomes_for_orthofinder/proteomes_orthofinder/OrthoFinder
1.5G	/srv/beegfs/scratch/users/m/miches/genomes_for_orthofinder/proteomes_orthofinder/OrthoFinder

real	0m6.823s
user	0m0.130s
sys	0m1.007s
```
For testing purposes, I have increased the timeout to 15 seconds.
Could you please restart your OOD session by clicking Developer or Help → Restart Web Server, then try again and let me know whether the issue is resolved?


## Post 7 by @Sebastien.Miche (2026-08-10T09:53:40.054Z)

I’ve restarted my session and tried again. This time I don’t get the timeout error message, but the download simply does not start, even after waiting for a while. Any idea of why this is happening?


## Post 8 by @Adrien.Albert (2026-08-10T15:31:55.563Z)

Hello Sebastien,
When using Open OnDemand directly to download a directory, the browser automatically creates and downloads it as a ZIP file. Depending on the size of the directory, it may take some time before the download starts.
However, for large directories, generating the ZIP file can exceed the 30-second timeout. I am not sure which parameter would need to be adjusted to increase this limit.
As a workaround, you can create a tar archive from a terminal by running:
```
tar cf genomes_for_orthofinder.tar genomes_for_orthofinder/
```
You can then download the resulting file through Open OnDemand.
Please let me know whether this works for you.


## Post 9 by @Sebastien.Miche (2026-08-11T08:38:40.487Z)

Adrien.Albert:
> `tar cf genomes_for_orthofinder.tar genomes_for_orthofinder/`
Hi Adrien,
I tried archiving the entire folder (700 MB) and downloading it but still had no success.
I managed to retrieve the files I needed but the rest of the data is unfortunately stuck for now.
Do you have any idea whether yggdrasil will be accessible again using ssh keys in the future?
Thanks for the help anyway!


## Post 10 by @Adrien.Albert (2026-08-11T12:36:57.255Z)

@Sebastien.Miche[@Sebastien.Miche](https://hpc-community.unige.ch/u/sebastien.miche)
Indeed the tar download does not end correctly. I will open an issue on OOD community.
Login1.yggdrasil.hpc.unige.ch is accessible via ssh, you can use rsync or filezilla to retrieve your data.
