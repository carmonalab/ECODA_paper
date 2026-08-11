# Are there any changes to the cpu063 node?

- Source: https://hpc-community.unige.ch/t/3765

- Created: 2024-12-10T20:25:51.574Z

- Tags: all

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Bharathkumar.Radhakrishnan (2024-12-10T20:25:51.612Z)

Hi admins
From yesterday night (US time zone), I started getting SIGILL errors on cpu063. There was no change from my side on the executables, and they continue to run properly on other nodes.
Is cpu063 a new node with a different instruction set? Or were there any changes to cpu063? I notice that now the status of cpu063 says `plnd`
(I have not kept a close look on which nodes my program use, so I don’t know if I used cpu063 in the past or not)
Thanks a lot!


## Post 2 by @Yann.Sagon (2024-12-11T14:45:35.211Z)

Dear @Bharathkumar.Radhakrishnan[@Bharathkumar.Radhakrishnan](https://hpc-community.unige.ch/u/bharathkumar.radhakrishnan)
Bharathkumar.Radhakrishnan:
> From yesterday night (US time zone), I started getting SIGILL errors on cpu063. There was no change from my side on the executables, and they continue to run properly on other nodes.
You can check the characteristics[characteristics](https://doc.eresearch.unige.ch/hpc/hpc_clusters#cpus_on_baobab) of the node in our doc. This is a V4. The login node is a V10.
Please check[check](https://doc.eresearch.unige.ch/hpc/faq#illegal_instruction) here for workaround.
