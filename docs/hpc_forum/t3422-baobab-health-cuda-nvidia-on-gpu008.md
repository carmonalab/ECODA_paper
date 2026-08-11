# [baobab] health_cuda___nvidia on gpu008

- Source: https://hpc-community.unige.ch/t/3422

- Created: 2024-04-17T16:11:50.650Z

- Posts: 2

- Category: 13

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2024-04-17T16:11:50.700Z)

Dear @support[@support](https://hpc-community.unige.ch/groups/support) Team,
I noticed that `gpu008` node is in drained state because of `health_cuda___nvidia` since yesterday. When can we expect the issue to be resolved?
Kind regards,
Maciej Falkiewicz


## Post 2 by @Adrien.Albert (2024-04-17T16:27:48.468Z)

Dear @maciej.falkiewicz[@maciej.falkiewicz](https://hpc-community.unige.ch/u/maciej.falkiewicz)
We have a lot of other tasks running in parallel, dealing with user requests, production issues, and system administration tasks as quickly as possible. While waiting for this node to return to production (which should be processed tomorrow), you can use a node from the shared-gpu partition.
`shared-gpu up 12:00:00 2 idle gpu[025-026]`
We understand that the unavailability of the node may pose difficulties, and we are working to resolve the issue as soon as possible.
