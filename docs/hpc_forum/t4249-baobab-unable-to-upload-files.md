# Baobab - unable to upload files

- Source: https://hpc-community.unige.ch/t/4249

- Created: 2026-03-12T16:05:01.019Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Lilou.Dehondt (2026-03-12T16:05:01.084Z)

Hello, there seems to be with inodes on baobab, but I can’t figure out exactly what is happening. I am unable to upload scripts (no issues with scratch). Can you please help me understand what is happening? Thank you for your help!
image
image591×168 4.81 KB
[image591×168 4.81 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/bb1f89647aa78f31be35338ab6a7bf276bfcb734.png)


## Post 2 by @Yann.Sagon (2026-03-13T08:19:35.933Z)

dear @Lilou.Dehondt[@Lilou.Dehondt](https://hpc-community.unige.ch/u/lilou.dehondt)
Lilou.Dehondt:
> there seems to be with inodes on baobab
I don’t understand why you think there is an issue about inodes?
Lilou.Dehondt:
> I am unable to upload scripts (no issues with scratch)
What do you mean? You are uploading scripts from outside of the cluster and it isn’t working? With what software? What is the error?
I don’t understand your screenshot, please explain what bothers you? If it is the fact that you don’t see the inode count, this is normal, we use a distributed filesystem.
You are using 937GB of storage in  home (out of 1024), maybe you were blocked by the quota?


## Post 3 by @Lilou.Dehondt (2026-03-16T12:16:47.613Z)

Hi,thank you for the quick answer and sorry for the not so precise issue.
Last week I was unable to upload 2 scripts (bash and python) from my computer to HPC using Filezilla. When trying to find the issue, I looked at my storage and inodes, because I was unaware that I could not see the inode count.
Today it is working again. My guess is that it’s because last week I generated a lot of files using a script on the server before freeing up some space. So I probably went over the quota, and since I tried uplading my scripts right after freeing up the space (as seen on the screenshot), it might have not “registered” yet, so I was getting a storage error.
