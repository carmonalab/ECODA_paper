# Many CPU nodes drained on private-dpt-cpu, and a couple on private-dpt-gpu

- Source: https://hpc-community.unige.ch/t/3734

- Created: 2024-11-19T22:21:07.797Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-11-19T22:21:07.837Z)

Hi,
7 nodes on the private-dpt-cpu partition are in drain, and seem to have been in this state for a while (don’t know the exact amount of time). Because of this, there seems to be a lot of pile up of jobs. Can the admins please take a look?
In the same vein, gpu025, gpu026, and gpu027 also seem to constantly go into drain and then be back over the last week (currently only gpu025 is in drain). Can the admins please have a look at that as well?
Thanks a lot!
Screenshot 2024-11-19 at 14.17.09
Screenshot 2024-11-19 at 14.17.091762×752 153 KB
[Screenshot 2024-11-19 at 14.17.091762×752 153 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/35ef7576ac3823b321fd26166bb368d92eab685e.png)


## Post 2 by @Yann.Sagon (2024-11-20T10:02:22.736Z)

Dear @Bharathkumar.Radhakrishnan[@Bharathkumar.Radhakrishnan](https://hpc-community.unige.ch/u/bharathkumar.radhakrishnan)
We are upgrading the BIOS of every compute node (issue-grl) and you are right this may block some of your compute nodes. We’re restored them now.
Extra hint: in Baobab and Bamboo, we have a lot of idle compute node that you can use as well if 12h00 of compute time is enough for your job. In this case, specify a max time of 12h00 or less, and specify the partition as `---partition=private-dpt-cpu,shared-cpu`
The issue with gpu is a known issue, unfortunately there is not much we can do other than supervise them and restore them. They are used by many people for multiple application and it is hard to understand what trigger them to go in drain.


## Post 3 by @Yann.Sagon (2025-01-30T09:08:27.381Z)

A post was split to a new topic: How will we be charged if we use a shared partition instead of a private one?[How will we be charged if we use a shared partition instead of a private one?](https://hpc-community.unige.ch/t/how-will-we-be-charged-if-we-use-a-shared-partition-instead-of-a-private-one/3804)


## Post 4 by @Yann.Sagon (2025-01-30T09:35:39.012Z)
