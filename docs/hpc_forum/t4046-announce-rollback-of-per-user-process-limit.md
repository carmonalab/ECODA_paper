# [Announce] Rollback of Per-User Process Limit

- Source: https://hpc-community.unige.ch/t/4046

- Created: 2025-08-12T08:13:46.112Z

- Tags: yggdrasil

- Posts: 2

- Category: 6

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Adrien.Albert (2025-08-12T08:13:46.215Z)

Dear HPC Users,
Following the recent maintenance of Yggdrasil, we introduced a per-user process limit to maintain a healthy shared environment and ensure that one user’s activity does not unintentionally impact others.
Based on recent feedback, we have identified unexpected issues, particularly with VS Code, which can spawn a large number of processes with each login and does not clean them up after logout. This can quickly lead to reaching the process limit, preventing the creation of new processes (including SSH sessions) and blocking the affected user from reconnecting to the cluster.
For these reasons, we have decided to temporarily disable the process limit while we work on additional measures to resolve these issues.
Thank you for your understanding,


## Post 2 by @Nicolas.Clairis1 (2025-08-13T14:04:43.151Z)

thanks for taking our feedback into account! :folded_hands: I hope you manage to implement a similar per user limit on Bamboo!
