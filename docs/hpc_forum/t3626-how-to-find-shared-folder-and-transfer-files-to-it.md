# How to find shared folder and transfer files to it

- Source: https://hpc-community.unige.ch/t/3626

- Created: 2024-09-04T15:29:10.820Z

- Posts: 2

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Michael.Silvestre (2024-09-04T15:29:10.861Z)

I would like to transfer files currently stored in OneDrive to a shared folder on the HPC cluster.
Adrien showed me the other day how to add it to my directory, but it was buried in sub-sub-sub folders and neither my colleague or I can find it anymore.
Additionnally, I would like to find it on FileZilla to transfer files to it. Can I even do it if the folder is not in my name ?
Thanks,
Michaël


## Post 2 by @Adrien.Albert (2024-09-10T07:28:13.867Z)

Hi @Michael.Silvestre[@Michael.Silvestre](https://hpc-community.unige.ch/u/michael.silvestre)
In my opinion, the easiest way is to connect your OneDrive to your laptop and use filezilla as if it were a local directory:
Here’s my OneDrive copy in Windows File Manager:
image
Then use FileZilla to transfer files from OneDrive on your laptop to the HPC cluster.
image
image1919×307 25 KB
[image1919×307 25 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/c297a5122861824f01ea209056521087a3e09490.png)
I tested and it’s working.
---
Michael.Silvestre:
> Adrien showed me the other day how to add it to my directory, but it was buried in sub-sub-sub folders and neither my colleague or I can find it anymore.
I sent you the path by email after creation:
`/home/share/pellizzari/priority`
`/srv/beegfs/scratch/shares/pellizzari/priority`
