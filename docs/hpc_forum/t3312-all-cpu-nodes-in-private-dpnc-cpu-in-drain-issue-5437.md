# All cpu nodes in private-dpnc-cpu in drain [issue-5437]

- Source: https://hpc-community.unige.ch/t/3312

- Created: 2024-02-15T08:54:04.737Z

- Posts: 4

- Category: 1

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @John.Raine (2024-02-15T08:54:04.793Z)

Hi,
Many cpu nodes, including all nodes on baobab that belong to private-dpnc-cpu, are in drain and have been unreliable for a few weeks.
What is “issue-5437” and how long is it expected for the nodes to be out of commission?
Cheers,
Johnny


## Post 2 by @Yann.Sagon (2024-02-16T15:38:21.000Z)

Dear @John.Raine[@John.Raine](https://hpc-community.unige.ch/u/john.raine), We are actually working for you!
We have many compute nodes on Baobab that aren’t new enough to upgrade Baobab Infiniband fabric. For this reason, we have ordered 150 Infiniband cards and are gradually replacing them. This takes some effort because we have to put the server in drain, unrack it, unscrew the IB card, replace it with a newer one, rerack the server, and hope that everything went well. This has to be done 150 times :slight_smile:


## Post 3 by @John.Raine (2024-02-19T10:35:09.334Z)

Thanks for the info Yann, is there somewhere in future we can look to see what the issues refer to?
One of our concerns was that with no nodes running for private-dpnc-cpu our Tier3 jobs we run on behalf of the WLCG were not able to start.
Cheers,
Johnny


## Post 4 by @Yann.Sagon (2024-02-22T15:42:08.687Z)

Hi,
unfortunately the issue named “issue-xxxx” are private issue we use them for our internal tasks. In this case we have already resumed some of your compute nodes and we plan to finish the replacement before the next cluster maintenance beginning of March. As this takes a lot of time to replace all that stuff, we preferred to drain the nodes instead of stopping the whole cluster as not all the nodes have to be upgraded.
